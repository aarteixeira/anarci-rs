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
//! `p7_Pipeline` — accumulating into one `P7_TOPHITS`. After all 29 models we
//! `p7_tophits_SortBySortkey` + `p7_tophits_Threshold` and walk every reported
//! domain into an [`anarci_core::Hsp`].
//!
//! ## E-value
//! In scan mode `p7_pli_NewModel` sets `pli->Z = nmodels`, so after all models
//! `Z = 29`. The HSP `.evalue` field ANARCI consumes is Biopython's `hsp.evalue`
//! = the *independent* E-value = `exp(dom->lnP) * Z`. We compute that directly
//! with `Z = 29` rather than relying on `pli->Z`/`pli->domZ` state.
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
use std::os::raw::c_char;
use std::path::Path;
use std::ptr;

use anarci_core::Hsp;
use anyhow::{anyhow, bail, Context, Result};
use hmmer_sys as ffi;

/// Number of models in `ALL.hmm`; also the E-value search-space size `Z` that
/// `hmmscan` scan mode uses (`pli->Z = nmodels`).
const Z: f64 = 29.0;

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
// pipeline through the clone). All per-scan mutable state — `ESL_SQ`,
// `P7_PIPELINE`, `P7_BG`, `P7_TOPHITS`, and the cloned profiles — is created,
// used, and destroyed within a single `scan` call on the calling thread and is
// never shared across threads. Therefore concurrent `&Engine` access from
// multiple threads touches the shared data read-only, which is sound.
unsafe impl Send for Engine {}
unsafe impl Sync for Engine {}

impl Engine {
    /// Load and configure every profile in an `ALL.hmm` (HMMER 3.4 plain-text or
    /// pressed) file.
    pub fn load(hmm_path: &Path) -> Result<Engine> {
        // One-time global init (idempotent; HMMER guards with a `firsttime`).
        unsafe { ffi::p7_FLogsumInit() };

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

        let mut abc: *mut ffi::ESL_ALPHABET = ptr::null_mut();
        let mut models: Vec<Model> = Vec::new();

        loop {
            let mut hmm: *mut ffi::P7_HMM = ptr::null_mut();
            let rs = unsafe { ffi::p7_hmmfile_Read(hfp, &mut abc, &mut hmm) };
            if rs == ffi::eslEOF as i32 {
                break;
            }
            if rs != ffi::eslOK as i32 || hmm.is_null() {
                unsafe { ffi::p7_hmmfile_Close(hfp) };
                bail!("p7_hmmfile_Read failed with status {rs}");
            }

            let m = unsafe { (*hmm).M };
            let name = unsafe { cstr_to_string((*hmm).name) }
                .context("profile has no NAME")?;

            // Build and configure: gm in p7_LOCAL with L_hint=100, then optimize.
            let bg = unsafe { ffi::p7_bg_Create(abc) };
            let gm = unsafe { ffi::p7_profile_Create(m, abc) };
            let om = unsafe { ffi::p7_oprofile_Create(m, abc) };
            if bg.is_null() || gm.is_null() || om.is_null() {
                unsafe { ffi::p7_hmmfile_Close(hfp) };
                bail!("allocation failed configuring profile {name}");
            }
            let cs = unsafe { ffi::p7_ProfileConfig(hmm, bg, gm, L_HINT, ffi::p7_LOCAL as i32) };
            if cs != ffi::eslOK as i32 {
                bail!("p7_ProfileConfig({name}) failed with status {cs}");
            }
            let cv = unsafe { ffi::p7_oprofile_Convert(gm, om) };
            if cv != ffi::eslOK as i32 {
                bail!("p7_oprofile_Convert({name}) failed with status {cv}");
            }

            // gm, bg, and the raw hmm are no longer needed: om is self-contained.
            unsafe {
                ffi::p7_profile_Destroy(gm);
                ffi::p7_bg_Destroy(bg);
                ffi::p7_hmm_Destroy(hmm);
            }

            let _ = &name; // read above only to surface a clear error on a nameless profile
            models.push(Model { om });
        }

        unsafe { ffi::p7_hmmfile_Close(hfp) };

        if models.is_empty() {
            bail!("no profiles read from {hmm_path:?}");
        }
        Ok(Engine { abc, models })
    }

    /// Number of loaded profiles.
    pub fn n_models(&self) -> usize {
        self.models.len()
    }

    /// Scan one sequence against every loaded profile, returning one [`Hsp`] per
    /// reported profile-domain — the same set ANARCI gets from `hmmscan`.
    ///
    /// `name` becomes the query name on the `ESL_SQ` (it does not affect scoring;
    /// HMMER only uses it for output labels we ignore). `seq` is the raw amino
    /// acid sequence (uppercase letters); it is digitized internally.
    pub fn scan(&self, name: &str, seq: &[u8]) -> Vec<Hsp> {
        THREAD_PROFILES.with(|cell| {
            let mut clones = cell.borrow_mut();
            clones.ensure(&self.models);
            self.scan_with(name, seq, &mut clones.profiles)
        })
    }

    fn scan_with(&self, name: &str, seq: &[u8], clones: &mut [*mut ffi::P7_OPROFILE]) -> Vec<Hsp> {
        // Build a digital ESL_SQ for the query.
        let name_c = CString::new(name).unwrap_or_else(|_| CString::new("query").unwrap());
        let seq_c = CString::new(seq).expect("query sequence contains a NUL byte");
        let sq = unsafe {
            ffi::esl_sq_CreateFrom(
                name_c.as_ptr(),
                seq_c.as_ptr(),
                ptr::null(),
                ptr::null(),
                ptr::null(),
            )
        };
        assert!(!sq.is_null(), "esl_sq_CreateFrom returned NULL");
        let ds = unsafe { ffi::esl_sq_Digitize(self.abc, sq) };
        assert_eq!(ds, ffi::eslOK as i32, "esl_sq_Digitize failed");
        let n = unsafe { (*sq).n };

        // Per-scan, per-thread pipeline + background, exactly as hmmscan.
        let pli = unsafe {
            ffi::p7_pipeline_Create(
                ptr::null(),
                100,
                100,
                0, // long_targets = FALSE
                ffi::p7_pipemodes_e_p7_SCAN_MODELS,
            )
        };
        assert!(!pli.is_null(), "p7_pipeline_Create returned NULL");
        let bg = unsafe { ffi::p7_bg_Create(self.abc) };
        assert!(!bg.is_null(), "p7_bg_Create returned NULL");
        let th = unsafe { ffi::p7_tophits_Create() };
        assert!(!th.is_null(), "p7_tophits_Create returned NULL");

        for om in clones.iter().copied() {
            unsafe {
                ffi::p7_pli_NewModel(pli, om, bg);
                ffi::p7_bg_SetLength(bg, n as i32);
                ffi::p7_oprofile_ReconfigLength(om, n as i32);
                ffi::p7_Pipeline(pli, om, bg, sq, ptr::null_mut(), th);
                ffi::p7_pipeline_Reuse(pli);
            }
        }

        unsafe {
            ffi::p7_tophits_SortBySortkey(th);
            ffi::p7_tophits_Threshold(th, pli);
        }

        let hsps = unsafe { collect_hsps(th) };

        unsafe {
            ffi::p7_tophits_Destroy(th);
            ffi::p7_bg_Destroy(bg);
            ffi::p7_pipeline_Destroy(pli);
            ffi::esl_sq_Destroy(sq);
        }

        hsps
    }
}

impl Drop for Engine {
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
}

impl ThreadProfiles {
    /// Clone the engine's profiles into this thread the first time it scans.
    fn ensure(&mut self, models: &[Model]) {
        if self.profiles.len() == models.len() {
            return;
        }
        // First use on this thread.
        for m in models {
            let c = unsafe { ffi::p7_oprofile_Clone(m.om) };
            assert!(!c.is_null(), "p7_oprofile_Clone returned NULL");
            self.profiles.push(c);
        }
    }
}

impl Drop for ThreadProfiles {
    fn drop(&mut self) {
        for &om in &self.profiles {
            unsafe { ffi::p7_oprofile_Destroy(om) };
        }
    }
}

thread_local! {
    static THREAD_PROFILES: std::cell::RefCell<ThreadProfiles> =
        std::cell::RefCell::new(ThreadProfiles { profiles: Vec::new() });
}

/// Walk a thresholded `P7_TOPHITS` into ANARCI `Hsp`s, in tophits order (hits
/// sorted by sortkey, domains in order) — the same order `hmmscan`'s output and
/// thus Biopython's HSP list use.
unsafe fn collect_hsps(th: *mut ffi::P7_TOPHITS) -> Vec<Hsp> {
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
                evalue: round_2g(lnp.exp() * Z),
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
