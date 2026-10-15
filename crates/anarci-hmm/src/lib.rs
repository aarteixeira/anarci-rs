//! In-process HMMER 3.4 scan engine.
//!
//! Reproduces what ANARCI's `hmmscan` subprocess produces — a `Vec<Hsp>` per
//! query sequence — by driving the HMMER 3.4 C pipeline directly through
//! `hmmer-sys`, with zero subprocess and zero text parsing.
//!
//! # Choreography (verified against HMMER 3.4 `src/hmmscan.c` serial loop and
//! Biopython's `Bio.SearchIO` hmmer3-text parser)
//!
//! On [`Engine::load`] we open `ALL.hmm`, read every profile, and build a fully
//! configured optimized profile (`P7_OPROFILE`) per model in `hmmscan` scan
//! mode (`p7_LOCAL`, `L_hint = 100`). `p7_FLogsumInit()` is called once.
//!
//! On [`Engine::scan`] we, per call: digitize the query into an `ESL_SQ`, take a
//! thread-local clone of every profile (so the shared profiles stay read-only),
//! and for each model run the same sequence of calls `hmmscan` runs —
//! `p7_pli_NewModel`, `p7_bg_SetLength`, `p7_oprofile_ReconfigLength`,
//! `p7_Pipeline` — accumulating into one `P7_TOPHITS`. After all loaded models we
//! `p7_tophits_SortBySortkey` + `p7_tophits_Threshold` and walk every reported
//! domain into an [`anarci_core::Hsp`].
//!
//! ## E-value
//! In scan mode `p7_pli_NewModel` sets `pli->Z = nmodels`, so after all models
//! `Z = nmodels`. The HSP `.evalue` field ANARCI consumes is Biopython's `hsp.evalue`
//! = the *independent* E-value = `exp(dom->lnP) * Z`. We compute that directly
//! with the number of loaded models (29 for exact, 7 for pan).
//! (`hmmscan`'s conditional E-value, `* domZ`, is `hsp.evalue_cond`, which ANARCI
//! does not use.)
//!
//! ## `impl_Init`
//! `hmmscan` calls `impl_Init()` once for processor-specific setup. On arm64
//! NEON (`impl_neon.h`) it only flushes subnormals when the SSE
//! `HAVE_FLUSH_ZERO_MODE` path is compiled in, which it is not on NEON — so it
//! is a no-op here and intentionally omitted. `p7_FLogsumInit()` is the only
//! mandatory one-time init and we make it.

use std::ffi::{CStr, CString};
use std::io::Write;
use std::os::raw::c_char;
use std::path::Path;
use std::ptr;
use std::sync::{Arc, Once};

use anarci_core::Hsp;
use anyhow::{anyhow, bail, Context, Result};
use hmmer_sys as ffi;

static LOGSUM_INIT: Once = Once::new();

/// `L_hint` HMMER uses to configure profiles in `hmmscan` (`src/hmmscan.c`).
const L_HINT: i32 = 100;

/// One loaded, fully configured optimized profile. The profile's name (e.g.
/// "mouse_H") lives inside `om->name` and is read back from the pipeline's hit
/// records during a scan, so we don't store it separately.
struct Model {
    om: *mut ffi::P7_OPROFILE,
}

/// A loaded set of HMMER profiles ready to scan query sequences in-process.
///
/// Built once with [`Engine::load`]; `scan` may be called concurrently from
/// many threads (see the `Sync` soundness note below).
pub struct Engine {
    profiles: Arc<Profiles>,
}

struct Profiles {
    abc: *mut ffi::ESL_ALPHABET,
    models: Vec<Model>,
}

// SAFETY (Send + Sync):
//
// After `load`, the `ESL_ALPHABET` and every `P7_OPROFILE` in `models` are
// treated as immutable for the rest of the `Engine`'s life. `scan` never
// mutates them: the only profile-mutating call in the pipeline,
// `p7_oprofile_ReconfigLength`, is applied exclusively to *thread-local clones*
// (`p7_oprofile_Clone`), never to the shared `om`. The shared profiles and
// alphabet are only ever read (by `p7_oprofile_Clone`, which copies, and by the
// pipeline through the clone). The thread cache retains an Arc to their owner,
// so cloned shared arrays outlive every scan. All per-scan mutable state — `ESL_SQ`,
// `P7_PIPELINE`, `P7_BG`, and `P7_TOPHITS` — is created,
// used, and destroyed within a single `scan` call on the calling thread and is
// never shared across threads. Therefore concurrent `&Engine` access from
// multiple threads touches the shared data read-only, which is sound.
unsafe impl Send for Profiles {}
unsafe impl Sync for Profiles {}

impl Engine {
    /// Load and configure every profile in an `ALL.hmm` (HMMER 3.4 plain-text or
    /// pressed) file.
    pub fn load(hmm_path: &Path) -> Result<Engine> {
        // HMMER's internal firsttime flag is not synchronized.
        LOGSUM_INIT.call_once(|| unsafe {
            ffi::p7_FLogsumInit();
        });

        let path_c = CString::new(
            hmm_path
                .to_str()
                .ok_or_else(|| anyhow!("hmm path is not valid UTF-8: {hmm_path:?}"))?,
        )?;

        let mut hfp: *mut ffi::P7_HMMFILE = ptr::null_mut();
        let mut errbuf = [0i8; 512];
        let status = unsafe {
            ffi::p7_hmmfile_Open(
                path_c.as_ptr(),
                ptr::null_mut(),
                &mut hfp,
                errbuf.as_mut_ptr() as *mut c_char,
            )
        };
        if status != ffi::eslOK as i32 || hfp.is_null() {
            let msg = unsafe { CStr::from_ptr(errbuf.as_ptr() as *const c_char) }
                .to_string_lossy()
                .into_owned();
            bail!("p7_hmmfile_Open({hmm_path:?}) failed (status {status}): {msg}");
        }

        let _file = Owned(hfp, ffi::p7_hmmfile_Close);
        let mut profiles = Profiles {
            abc: ptr::null_mut(),
            models: Vec::new(),
        };
        loop {
            let mut hmm = ptr::null_mut();
            let rs = unsafe { ffi::p7_hmmfile_Read(hfp, &mut profiles.abc, &mut hmm) };
            let hmm = Owned(hmm, ffi::p7_hmm_Destroy);
            if rs == ffi::eslEOF as i32 {
                break;
            }
            check_status(rs, "p7_hmmfile_Read")?;
            if hmm.0.is_null() {
                bail!("p7_hmmfile_Read returned NULL");
            }
            if profiles.abc.is_null() || unsafe { (*profiles.abc).type_ } != ffi::eslAMINO as i32 {
                bail!("profiles must use the HMMER amino-acid alphabet");
            }
            let m = unsafe { (*hmm.0).M };
            let name = unsafe { cstr_to_string((*hmm.0).name) }.context("profile has no NAME")?;
            let bg = Owned::checked(
                unsafe { ffi::p7_bg_Create(profiles.abc) },
                ffi::p7_bg_Destroy,
                "p7_bg_Create",
            )?;
            let gm = Owned::checked(
                unsafe { ffi::p7_profile_Create(m, profiles.abc) },
                ffi::p7_profile_Destroy,
                "p7_profile_Create",
            )?;
            let om = Owned::checked(
                unsafe { ffi::p7_oprofile_Create(m, profiles.abc) },
                ffi::p7_oprofile_Destroy,
                "p7_oprofile_Create",
            )?;
            check_status(
                unsafe { ffi::p7_ProfileConfig(hmm.0, bg.0, gm.0, L_HINT, ffi::p7_LOCAL as i32) },
                &format!("p7_ProfileConfig({name})"),
            )?;
            check_status(
                unsafe { ffi::p7_oprofile_Convert(gm.0, om.0) },
                &format!("p7_oprofile_Convert({name})"),
            )?;
            profiles.models.push(Model { om: om.into_raw() });
        }
        if profiles.models.is_empty() {
            bail!("no profiles read from {hmm_path:?}");
        }
        Ok(Engine {
            profiles: Arc::new(profiles),
        })
    }

    /// Load embedded profiles using a unique private temporary file. The file is
    /// removed after loading, including when parsing fails.
    pub fn load_embedded(bytes: &[u8]) -> Result<Engine> {
        let mut file = tempfile::NamedTempFile::new().context("create temporary HMM file")?;
        file.write_all(bytes).context("write temporary HMM file")?;
        file.flush().context("flush temporary HMM file")?;
        Self::load(file.path())
    }

    /// Number of loaded profiles.
    pub fn n_models(&self) -> usize {
        self.profiles.models.len()
    }

    /// Check that profiles satisfy the numbering layer's 128-state IMGT
    /// contract. Generic HMM scanning does not require this restriction.
    pub fn validate_numbering_profiles(&self) -> Result<()> {
        for model in &self.profiles.models {
            let name = unsafe { cstr_to_string((*model.om).name) }
                .context("numbering profile has no NAME")?;
            let valid_name = name.split_once('_').is_some_and(|(species, chain)| {
                !species.is_empty() && matches!(chain, "H" | "K" | "L" | "A" | "B" | "G" | "D")
            });
            if !valid_name {
                bail!("numbering profile NAME {name:?} must be species_chain with a nonempty species and chain H, K, L, A, B, G, or D");
            }
            let length = unsafe { (*model.om).M };
            if length != 128 {
                bail!("numbering profile {name:?} has {length} states; numbering requires 128 IMGT states");
            }
            // HMMER om->rf contains a sentinel then one annotation per model
            // state. The numbering layer interprets only `x` as a match state.
            let rf = unsafe { (*model.om).rf };
            let canonical_rf = !rf.is_null()
                && unsafe {
                    let annotation = CStr::from_ptr(rf).to_bytes();
                    annotation.len() == 129 && annotation[1..].iter().all(|&c| c == b'x')
                };
            if !canonical_rf {
                bail!("numbering profile {name:?} requires RF annotations with 128 x match states");
            }
        }
        Ok(())
    }

    /// Scan one sequence against every loaded profile, returning one [`Hsp`] per
    /// reported profile-domain — the same set ANARCI gets from `hmmscan`.
    ///
    /// `name` becomes the query name on the `ESL_SQ` (it does not affect scoring;
    /// HMMER only uses it for output labels we ignore). `seq` is the raw amino
    /// acid sequence (uppercase letters); it is digitized internally.
    pub fn scan(&self, name: &str, seq: &[u8]) -> Vec<Hsp> {
        self.scan_checked(name, seq)
            .expect("valid query and successful HMMER scan")
    }

    /// Scan a query with fallible validation and native error propagation.
    /// Empty sequences return no hits. HMMER's amino-acid letters A-Z (including
    /// ambiguous and unusual residues) are accepted in either case, as is `*`.
    pub fn scan_checked(&self, name: &str, seq: &[u8]) -> Result<Vec<Hsp>> {
        anarci_core::orchestrate::validate_scan_input(name, seq)?;
        if seq.is_empty() {
            return Ok(Vec::new());
        }
        THREAD_PROFILES.with(|cell| {
            let mut clones = cell.borrow_mut();
            clones.ensure(&self.profiles)?;
            self.scan_with(name, seq, &mut clones.profiles)
        })
    }

    fn scan_with(
        &self,
        name: &str,
        seq: &[u8],
        clones: &mut [*mut ffi::P7_OPROFILE],
    ) -> Result<Vec<Hsp>> {
        let name_c = CString::new(name)?;
        let seq_c = CString::new(seq)?;
        let sq = Owned::checked(
            unsafe {
                ffi::esl_sq_CreateFrom(
                    name_c.as_ptr(),
                    seq_c.as_ptr(),
                    ptr::null(),
                    ptr::null(),
                    ptr::null(),
                )
            },
            ffi::esl_sq_Destroy,
            "esl_sq_CreateFrom",
        )?;
        check_status(
            unsafe { ffi::esl_sq_Digitize(self.profiles.abc, sq.0) },
            "esl_sq_Digitize",
        )?;
        let n = seq.len() as i32; // validated before entering native code
        let pli = Owned::checked(
            unsafe {
                ffi::p7_pipeline_Create(
                    ptr::null(),
                    100,
                    100,
                    0,
                    ffi::p7_pipemodes_e_p7_SCAN_MODELS,
                )
            },
            ffi::p7_pipeline_Destroy,
            "p7_pipeline_Create",
        )?;
        let bg = Owned::checked(
            unsafe { ffi::p7_bg_Create(self.profiles.abc) },
            ffi::p7_bg_Destroy,
            "p7_bg_Create",
        )?;
        let th = Owned::checked(
            unsafe { ffi::p7_tophits_Create() },
            ffi::p7_tophits_Destroy,
            "p7_tophits_Create",
        )?;
        for om in clones.iter().copied() {
            unsafe {
                check_status(ffi::p7_pli_NewModel(pli.0, om, bg.0), "p7_pli_NewModel")?;
                check_status(ffi::p7_bg_SetLength(bg.0, n), "p7_bg_SetLength")?;
                check_status(
                    ffi::p7_oprofile_ReconfigLength(om, n),
                    "p7_oprofile_ReconfigLength",
                )?;
                check_status(
                    ffi::p7_Pipeline(pli.0, om, bg.0, sq.0, ptr::null_mut(), th.0),
                    "p7_Pipeline",
                )?;
                check_status(ffi::p7_pipeline_Reuse(pli.0), "p7_pipeline_Reuse")?;
            }
        }
        unsafe {
            check_status(
                ffi::p7_tophits_SortBySortkey(th.0),
                "p7_tophits_SortBySortkey",
            )?;
            check_status(
                ffi::p7_tophits_Threshold(th.0, pli.0),
                "p7_tophits_Threshold",
            )?;
            Ok(collect_hsps(th.0, self.n_models() as f64))
        }
    }
}

impl Drop for Profiles {
    fn drop(&mut self) {
        for m in &self.models {
            unsafe { ffi::p7_oprofile_Destroy(m.om) };
        }
        if !self.abc.is_null() {
            unsafe { ffi::esl_alphabet_Destroy(self.abc) };
        }
    }
}

/// Per-thread cache of profile clones, reconfigured in place each scan so we
/// never mutate the shared read-only profiles and never reallocate per scan.
struct ThreadProfiles {
    profiles: Vec<*mut ffi::P7_OPROFILE>,
    // HMMER clones share arrays and alphabet with their source profiles.
    owner: Option<Arc<Profiles>>,
}

impl ThreadProfiles {
    fn clear(&mut self) {
        for om in self.profiles.drain(..) {
            unsafe { ffi::p7_oprofile_Destroy(om) };
        }
        self.owner = None;
    }

    fn ensure(&mut self, profiles: &Arc<Profiles>) -> Result<()> {
        if self
            .owner
            .as_ref()
            .is_some_and(|owner| Arc::ptr_eq(owner, profiles))
        {
            return Ok(());
        }
        self.clear();
        for m in &profiles.models {
            let c = unsafe { ffi::p7_oprofile_Clone(m.om) };
            if c.is_null() {
                self.clear();
                bail!("p7_oprofile_Clone returned NULL");
            }
            self.profiles.push(c);
        }
        self.owner = Some(Arc::clone(profiles));
        Ok(())
    }
}

impl Drop for ThreadProfiles {
    fn drop(&mut self) {
        self.clear();
    }
}

thread_local! {
    static THREAD_PROFILES: std::cell::RefCell<ThreadProfiles> =
        std::cell::RefCell::new(ThreadProfiles { profiles: Vec::new(), owner: None });
}

/// Own a native allocation so every error path releases it.
struct Owned<T>(*mut T, unsafe extern "C" fn(*mut T));

impl<T> Owned<T> {
    fn checked(p: *mut T, destroy: unsafe extern "C" fn(*mut T), operation: &str) -> Result<Self> {
        if p.is_null() {
            bail!("{operation} returned NULL");
        }
        Ok(Self(p, destroy))
    }

    fn into_raw(mut self) -> *mut T {
        let p = self.0;
        self.0 = ptr::null_mut();
        p
    }
}

impl<T> Drop for Owned<T> {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { (self.1)(self.0) };
        }
    }
}

fn check_status(status: i32, operation: &str) -> Result<()> {
    if status != ffi::eslOK as i32 {
        bail!("{operation} failed with status {status}");
    }
    Ok(())
}

/// Walk a thresholded `P7_TOPHITS` into ANARCI `Hsp`s, in tophits order (hits
/// sorted by sortkey, domains in order) — the same order `hmmscan`'s output and
/// thus Biopython's HSP list use.
unsafe fn collect_hsps(th: *mut ffi::P7_TOPHITS, z: f64) -> Vec<Hsp> {
    let mut out: Vec<Hsp> = Vec::new();
    let nhits = (*th).N as isize;
    for i in 0..nhits {
        let hit = *(*th).hit.offset(i);
        // Only reported hits (mirrors hmmscan's output gate).
        if ((*hit).flags & ffi::p7_IS_REPORTED) == 0 {
            continue;
        }
        let name = cstr_to_string((*hit).name).unwrap_or_default();
        let ndom = (*hit).ndom as isize;
        for d in 0..ndom {
            let dom = (*hit).dcl.offset(d);
            if (*dom).is_reported == 0 {
                continue;
            }
            let ad = (*dom).ad;
            if ad.is_null() {
                continue;
            }
            let rf = cstr_to_string((*ad).rfline).unwrap_or_default();
            let pp = cstr_to_string((*ad).ppline).unwrap_or_default();
            // ANARCI consumes Biopython's parse of hmmscan's *text* output, so
            // every numeric field is the value printed there, at its print
            // precision — not the raw C float. We reproduce hmmscan's formats:
            //   bitscore  %6.1f  -> 1 decimal
            //   bias      %5.1f  on (dombias * LOG2R), NATS -> BITS -> 1 decimal
            //   i-Evalue  %9.2g  on exp(lnP) * Z   -> 2 significant figures
            let lnp = (*dom).lnP;
            let bias_bits = (*dom).dombias as f64 * ffi::eslCONST_LOG2R;
            out.push(Hsp {
                hit_id: name.clone(),
                hit_description: String::new(),
                evalue: independent_evalue(lnp, z),
                bitscore: round_1f((*dom).bitscore as f64),
                bias: round_1f(bias_bits),
                query_start: ((*ad).sqfrom - 1) as usize,
                query_end: (*ad).sqto as usize,
                hit_start: ((*ad).hmmfrom - 1) as usize,
                hit_end: (*ad).hmmto as usize,
                rf,
                pp,
                order: 0,
            });
        }
    }
    out
}

fn independent_evalue(lnp: f64, z: f64) -> f64 {
    round_2g(lnp.exp() * z)
}

/// Reproduce C `printf("%.1f", x)` then parse-back: the value at one-decimal
/// print precision, as Biopython reads it from hmmscan's text output.
fn round_1f(x: f64) -> f64 {
    format!("{x:.1}").parse().unwrap()
}

/// Reproduce C `printf("%.2g", x)`: two significant figures. Rust's `{:.1e}`
/// gives mantissa-with-one-decimal scientific (i.e. 2 sig figs); parsing it
/// back yields the same value Biopython parses from the `%9.2g` column.
fn round_2g(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    format!("{x:.1e}").parse().unwrap()
}

unsafe fn cstr_to_string(p: *const c_char) -> Option<String> {
    if p.is_null() {
        return None;
    }
    Some(CStr::from_ptr(p).to_string_lossy().into_owned())
}

#[cfg(test)]
mod robustness_tests {
    use super::*;

    // HMM resources are model parameters; queries in these tests are synthetic.
    fn model_bytes(name: &str, count: usize) -> Vec<u8> {
        let source = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../reference_data/dat/HMMs/ALL.hmm"
        ));
        let first = source.split_once("//\n").expect("profile record").0;
        let mut bytes = Vec::new();
        for index in 0..count {
            for line in first.lines() {
                if line.starts_with("NAME ") {
                    writeln!(bytes, "NAME  {name}_{index}").unwrap();
                } else {
                    writeln!(bytes, "{line}").unwrap();
                }
            }
            writeln!(bytes, "//").unwrap();
        }
        bytes
    }

    #[test]
    fn cache_tracks_identity_replaces_counts_and_keeps_source_alive() {
        let a = Engine::load_embedded(&model_bytes("synthetic_a", 1)).unwrap();
        let b = Engine::load_embedded(&model_bytes("synthetic_b", 1)).unwrap();
        let c = Engine::load_embedded(&model_bytes("synthetic_c", 2)).unwrap();
        let mut cache = ThreadProfiles {
            profiles: Vec::new(),
            owner: None,
        };
        cache.ensure(&a.profiles).unwrap();
        let first = cache.profiles[0];
        cache.ensure(&a.profiles).unwrap();
        assert_eq!(
            cache.profiles[0], first,
            "same-engine scan reuses its clone"
        );
        let weak = Arc::downgrade(&a.profiles);
        drop(a);
        assert!(
            weak.upgrade().is_some(),
            "shallow clone must retain its source"
        );
        cache.ensure(&b.profiles).unwrap();
        assert!(weak.upgrade().is_none());
        assert_eq!(cache.profiles.len(), 1);
        assert_eq!(
            unsafe { cstr_to_string((*cache.profiles[0]).name) }.unwrap(),
            "synthetic_b_0"
        );
        cache.ensure(&c.profiles).unwrap();
        assert_eq!(cache.profiles.len(), 2);
        cache.ensure(&b.profiles).unwrap();
        assert_eq!(
            cache.profiles.len(),
            1,
            "switching must replace rather than append"
        );
        assert_eq!(
            unsafe { cstr_to_string((*cache.profiles[0]).name) }.unwrap(),
            "synthetic_b_0"
        );
    }

    #[test]
    fn checked_scan_rejects_bad_input_and_accepts_ambiguous_residues() {
        let engine = Engine::load_embedded(&model_bytes("synthetic", 1)).unwrap();
        for seq in [b"A\0A".as_slice(), b"A A", b"A1A", b"A-A", &[0xff]] {
            assert!(engine.scan_checked("query", seq).is_err());
        }
        assert!(engine.scan_checked("bad\0name", b"AAA").is_err());
        assert!(engine.scan_checked("long", &vec![b'A'; 100_001]).is_err());
        assert!(engine.scan_checked("empty", b"").unwrap().is_empty());
        engine
            .scan_checked("ambiguous", b"ABJZOUXabjzoux*")
            .unwrap();
        engine.scan_checked("synthetic", &[b'A'; 128]).unwrap();
    }

    #[test]
    fn embedded_loads_and_scans_are_safe_in_parallel() {
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for i in 0..8 {
                handles.push(scope.spawn(move || {
                    let engine =
                        Engine::load_embedded(&model_bytes(&format!("synthetic_{i}"), 1)).unwrap();
                    for _ in 0..3 {
                        engine.scan_checked("query", &[b'A'; 96]).unwrap();
                    }
                    THREAD_PROFILES.with(|cell| {
                        let cache = cell.borrow();
                        assert!(Arc::ptr_eq(cache.owner.as_ref().unwrap(), &engine.profiles));
                        assert_eq!(
                            unsafe { cstr_to_string((*cache.profiles[0]).name) }.unwrap(),
                            format!("synthetic_{i}_0")
                        );
                    });
                }));
            }
            for handle in handles {
                handle.join().unwrap();
            }
        });
    }

    #[test]
    fn shared_engine_matches_across_worker_counts() {
        use rayon::prelude::*;
        let engine = Engine::load_embedded(&model_bytes("synthetic", 1)).unwrap();
        let queries = [vec![b'A'; 64], vec![b'X'; 96], vec![b'G'; 128]];
        let expected: Vec<_> = queries
            .iter()
            .map(|q| format!("{:?}", engine.scan_checked("query", q).unwrap()))
            .collect();
        for workers in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            let actual: Vec<_> = pool.install(|| {
                queries
                    .par_iter()
                    .map(|q| format!("{:?}", engine.scan_checked("query", q).unwrap()))
                    .collect()
            });
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn numbering_profiles_require_chain_names_and_128_states() {
        let resource = String::from_utf8(model_bytes("synthetic", 1)).unwrap();
        for name in [
            "synthetic",
            "_H",
            "synthetic_Q",
            "synthetic_extra_H",
            "synthetic_",
        ] {
            let bytes = resource.replace("NAME  synthetic_0", &format!("NAME  {name}"));
            let engine = Engine::load_embedded(bytes.as_bytes()).unwrap();
            assert!(engine.validate_numbering_profiles().is_err(), "{name}");
            engine.scan_checked("synthetic_query", &[b'A'; 64]).unwrap();
        }
        for chain in ["H", "K", "L", "A", "B", "G", "D"] {
            let bytes = resource.replace("NAME  synthetic_0", &format!("NAME  synthetic_{chain}"));
            let engine = Engine::load_embedded(bytes.as_bytes()).unwrap();
            engine.validate_numbering_profiles().unwrap();
        }
        let bytes = resource.replace("NAME  synthetic_0", "NAME  synthetic_H");
        let engine = Engine::load_embedded(bytes.as_bytes()).unwrap();
        // Alter model metadata only, with no biological query or scoring step.
        let om = engine.profiles.models[0].om;
        unsafe {
            (*om).M = 4;
        }
        let result = engine.validate_numbering_profiles();
        unsafe {
            (*om).M = 128;
        }
        assert!(result.unwrap_err().to_string().contains("requires 128"));
    }

    #[test]
    fn numbering_profiles_require_canonical_reference_annotations() {
        let bytes = String::from_utf8(model_bytes("synthetic", 1))
            .unwrap()
            .replace("NAME  synthetic_0", "NAME  synthetic_H");
        let engine = Engine::load_embedded(bytes.as_bytes()).unwrap();
        engine.validate_numbering_profiles().unwrap();
        let rf = unsafe { (*engine.profiles.models[0].om).rf };
        // HMMER uses a zero first byte to indicate absent RF annotations.
        let sentinel = unsafe { *rf };
        unsafe {
            *rf = 0;
        }
        let missing = engine.validate_numbering_profiles();
        unsafe {
            *rf = sentinel;
        }
        assert!(missing.unwrap_err().to_string().contains("requires RF"));
        unsafe {
            *rf.add(64) = b'-' as c_char;
        }
        let wrong = engine.validate_numbering_profiles();
        unsafe {
            *rf.add(64) = b'x' as c_char;
        }
        assert!(wrong.unwrap_err().to_string().contains("requires RF"));
        engine.validate_numbering_profiles().unwrap();
    }

    #[test]
    fn malformed_embedded_profiles_return_errors() {
        assert!(Engine::load_embedded(b"").is_err());
        assert!(Engine::load_embedded(b"not a profile\n").is_err());
        let mut truncated = model_bytes("synthetic", 1);
        truncated.truncate(truncated.len() / 2);
        assert!(Engine::load_embedded(&truncated).is_err());
    }

    #[test]
    fn independent_evalues_use_the_loaded_database_size() {
        unsafe {
            let mut ad: ffi::P7_ALIDISPLAY = std::mem::zeroed();
            ad.sqfrom = 1;
            ad.sqto = 1;
            ad.hmmfrom = 1;
            ad.hmmto = 1;
            let mut domain: ffi::P7_DOMAIN = std::mem::zeroed();
            domain.ad = &mut ad;
            domain.is_reported = 1;
            domain.lnP = 0.01_f64.ln();
            let mut hit: ffi::P7_HIT = std::mem::zeroed();
            hit.flags = ffi::p7_IS_REPORTED;
            hit.ndom = 1;
            hit.dcl = &mut domain;
            let mut hit_ptr = &mut hit as *mut _;
            let mut hits: ffi::P7_TOPHITS = std::mem::zeroed();
            hits.N = 1;
            hits.hit = &mut hit_ptr;
            assert_eq!(collect_hsps(&mut hits, 29.0)[0].evalue, 0.29);
            assert_eq!(collect_hsps(&mut hits, 7.0)[0].evalue, 0.07);
        }
    }
}
