//! EXPERIMENT (not a parity gate). Quantifies how much faster the in-process
//! HMMER scan can be made while keeping ANARCI output identical.
//!
//! Run explicitly (single test process, single thread for clean timing):
//!   cargo test -p anarci-hmm --release --test bench_stages -- --ignored --nocapture --test-threads=1
//!
//! Does NOT touch production crates. Replicates the engine's load + per-scan
//! choreography here with instrumentation, and reuses anarci_core::parse_hmmer_query
//! for parity checking against tests/fixtures/hsps.json.gz.

use anarci_core::{parse_hmmer_query, Hsp};
use flate2::read::GzDecoder;
use hmmer_sys as ffi;
use std::ffi::{c_char, CStr, CString};
use std::path::PathBuf;
use std::ptr;
use std::time::Instant;

const Z: f64 = 29.0;
const L_HINT: i32 = 100;
const ESL_SMALLX1: f64 = 5e-9;
const LOG2: f64 = 0.69314718055994529;

// evparam indices (hmmer.h enum p7_evparams_e)
const P7_MMU: usize = 0;
const P7_MLAMBDA: usize = 1;
const P7_VMU: usize = 2;
const P7_VLAMBDA: usize = 3;
const P7_FTAU: usize = 4;
const P7_FLAMBDA: usize = 5;

// pipeline default thresholds (p7_pipeline.c)
const F1: f64 = 0.02;
const F2: f64 = 1e-3;
const F3: f64 = 1e-5;

fn hmm_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference_data/dat/HMMs/ALL.hmm")
}
fn fx_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/hsps.json.gz")
}

// --- gumbel / exp survivor functions (exact ports of easel) ---
fn gumbel_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    let ey = -(-y).exp();
    if ey.abs() < ESL_SMALLX1 {
        -ey
    } else {
        1.0 - ey.exp()
    }
}
fn exp_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    if x < mu {
        1.0
    } else {
        (-lambda * (x - mu)).exp()
    }
}

fn load_seqs() -> Vec<(String, Vec<u8>)> {
    let f = std::fs::File::open(fx_path()).unwrap();
    let v: serde_json::Value = serde_json::from_reader(GzDecoder::new(f)).unwrap();
    v["sequences"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| {
            (
                s["id"].as_str().unwrap().to_string(),
                s["seq"].as_str().unwrap().as_bytes().to_vec(),
            )
        })
        .collect()
}

// ---------------------------------------------------------------------------
// A self-contained engine replica that exposes the FFI internals we need.
// Mirrors crates/anarci-hmm/src/lib.rs choreography exactly.
// ---------------------------------------------------------------------------
struct Eng {
    abc: *mut ffi::ESL_ALPHABET,
    oms: Vec<*mut ffi::P7_OPROFILE>, // one fully-configured optimized profile per model
}

unsafe fn cstr(p: *const c_char) -> String {
    if p.is_null() {
        String::new()
    } else {
        CStr::from_ptr(p).to_string_lossy().into_owned()
    }
}

impl Eng {
    fn load() -> Eng {
        unsafe { ffi::p7_FLogsumInit() };
        let path_c = CString::new(hmm_path().to_str().unwrap()).unwrap();
        let mut hfp: *mut ffi::P7_HMMFILE = ptr::null_mut();
        let mut errbuf = [0i8; 512];
        let st = unsafe {
            ffi::p7_hmmfile_Open(path_c.as_ptr(), ptr::null_mut(), &mut hfp, errbuf.as_mut_ptr())
        };
        assert_eq!(st, ffi::eslOK as i32, "open ALL.hmm");
        let mut abc: *mut ffi::ESL_ALPHABET = ptr::null_mut();
        let mut oms = Vec::new();
        loop {
            let mut hmm: *mut ffi::P7_HMM = ptr::null_mut();
            let rs = unsafe { ffi::p7_hmmfile_Read(hfp, &mut abc, &mut hmm) };
            if rs == ffi::eslEOF as i32 {
                break;
            }
            assert_eq!(rs, ffi::eslOK as i32);
            let m = unsafe { (*hmm).M };
            let bg = unsafe { ffi::p7_bg_Create(abc) };
            let gm = unsafe { ffi::p7_profile_Create(m, abc) };
            let om = unsafe { ffi::p7_oprofile_Create(m, abc) };
            unsafe { ffi::p7_ProfileConfig(hmm, bg, gm, L_HINT, ffi::p7_LOCAL as i32) };
            unsafe { ffi::p7_oprofile_Convert(gm, om) };
            unsafe {
                ffi::p7_profile_Destroy(gm);
                ffi::p7_bg_Destroy(bg);
                ffi::p7_hmm_Destroy(hmm);
            }
            oms.push(om);
        }
        unsafe { ffi::p7_hmmfile_Close(hfp) };
        assert_eq!(oms.len(), 29);
        Eng { abc, oms }
    }

    // Clone the shared profiles for per-thread use (single-threaded here).
    fn clones(&self) -> Vec<*mut ffi::P7_OPROFILE> {
        self.oms
            .iter()
            .map(|&om| unsafe { ffi::p7_oprofile_Clone(om) })
            .collect()
    }

    fn make_sq(&self, name: &str, seq: &[u8]) -> *mut ffi::ESL_SQ {
        let name_c = CString::new(name).unwrap_or_else(|_| CString::new("q").unwrap());
        let seq_c = CString::new(seq).unwrap();
        let sq = unsafe {
            ffi::esl_sq_CreateFrom(
                name_c.as_ptr(),
                seq_c.as_ptr(),
                ptr::null(),
                ptr::null(),
                ptr::null(),
            )
        };
        unsafe { ffi::esl_sq_Digitize(self.abc, sq) };
        sq
    }
}

impl Drop for Eng {
    fn drop(&mut self) {
        for &om in &self.oms {
            unsafe { ffi::p7_oprofile_Destroy(om) };
        }
        unsafe { ffi::esl_alphabet_Destroy(self.abc) };
    }
}

fn round_1f(x: f64) -> f64 {
    format!("{x:.1}").parse().unwrap()
}
fn round_2g(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    format!("{x:.1e}").parse().unwrap()
}

// Collect HSPs from a thresholded tophits (same as production collect_hsps).
unsafe fn collect_hsps(th: *mut ffi::P7_TOPHITS) -> Vec<Hsp> {
    let mut out = Vec::new();
    let nhits = (*th).N as isize;
    for i in 0..nhits {
        let hit = *(*th).hit.offset(i);
        if ((*hit).flags & ffi::p7_IS_REPORTED) == 0 {
            continue;
        }
        let name = cstr((*hit).name);
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
            let lnp = (*dom).lnP;
            let bias_bits = (*dom).dombias as f64 / LOG2;
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
                rf: cstr((*ad).rfline),
                pp: cstr((*ad).ppline),
                order: 0,
            });
        }
    }
    out
}

// Full production-equivalent scan: run p7_Pipeline over all 29 profiles.
unsafe fn scan_full(
    e: &Eng,
    clones: &mut [*mut ffi::P7_OPROFILE],
    name: &str,
    seq: &[u8],
) -> Vec<Hsp> {
    let sq = e.make_sq(name, seq);
    let n = (*sq).n as i32;
    let pli = ffi::p7_pipeline_Create(
        ptr::null(),
        100,
        100,
        0,
        ffi::p7_pipemodes_e_p7_SCAN_MODELS,
    );
    let bg = ffi::p7_bg_Create(e.abc);
    let th = ffi::p7_tophits_Create();
    for &om in clones.iter() {
        ffi::p7_pli_NewModel(pli, om, bg);
        ffi::p7_bg_SetLength(bg, n);
        ffi::p7_oprofile_ReconfigLength(om, n);
        ffi::p7_Pipeline(pli, om, bg, sq, ptr::null_mut(), th);
        ffi::p7_pipeline_Reuse(pli);
    }
    ffi::p7_tophits_SortBySortkey(th);
    ffi::p7_tophits_Threshold(th, pli);
    let hsps = collect_hsps(th);
    ffi::p7_tophits_Destroy(th);
    ffi::p7_bg_Destroy(bg);
    ffi::p7_pipeline_Destroy(pli);
    ffi::esl_sq_Destroy(sq);
    hsps
}

// ============================================================================
// TASK 1: stage breakdown
// ============================================================================
#[test]
#[ignore]
fn task1_stage_breakdown() {
    let e = Eng::load();
    let mut clones = e.clones();
    let seqs = load_seqs();

    // --- 1a: counters from full pipeline runs ---
    // Sum n_past_* across all profiles for all sequences via a per-seq pipeline
    // (so counters are not reset until we read them).
    let (mut tot_msv, mut tot_bias, mut tot_vit, mut tot_fwd) = (0u64, 0u64, 0u64, 0u64);
    let nseq = seqs.len() as u64;
    unsafe {
        for (name, seq) in &seqs {
            let sq = e.make_sq(name, seq);
            let n = (*sq).n as i32;
            let pli = ffi::p7_pipeline_Create(
                ptr::null(),
                100,
                100,
                0,
                ffi::p7_pipemodes_e_p7_SCAN_MODELS,
            );
            let bg = ffi::p7_bg_Create(e.abc);
            let th = ffi::p7_tophits_Create();
            for &om in clones.iter() {
                ffi::p7_pli_NewModel(pli, om, bg);
                ffi::p7_bg_SetLength(bg, n);
                ffi::p7_oprofile_ReconfigLength(om, n);
                ffi::p7_Pipeline(pli, om, bg, sq, ptr::null_mut(), th);
                ffi::p7_pipeline_Reuse(pli);
            }
            // p7_pipeline_Reuse does NOT reset the n_past_* counters, so after the
            // 29-profile loop they hold this query's totals across all 29 models.
            tot_msv += (*pli).n_past_msv;
            tot_bias += (*pli).n_past_bias;
            tot_vit += (*pli).n_past_vit;
            tot_fwd += (*pli).n_past_fwd;
            ffi::p7_tophits_Destroy(th);
            ffi::p7_bg_Destroy(bg);
            ffi::p7_pipeline_Destroy(pli);
            ffi::esl_sq_Destroy(sq);
        }
    }
    eprintln!("\n=== TASK 1a: avg # of 29 profiles reaching each stage (per query) ===");
    eprintln!("  reached MSV gate pass (n_past_msv):  {:.3}", tot_msv as f64 / nseq as f64);
    eprintln!("  reached bias gate pass (n_past_bias): {:.3}", tot_bias as f64 / nseq as f64);
    eprintln!("  reached Viterbi pass  (n_past_vit):  {:.3}", tot_vit as f64 / nseq as f64);
    eprintln!(
        "  reached Forward pass -> DOMAINDEF (n_past_fwd): {:.3}",
        tot_fwd as f64 / nseq as f64
    );

    // --- 1b: timing full pipeline vs filters-only ---
    let reps = 5usize;
    // warm
    for (name, seq) in &seqs {
        unsafe { scan_full(&e, &mut clones, name, seq) };
    }

    // full
    let t = Instant::now();
    for _ in 0..reps {
        for (name, seq) in &seqs {
            let _ = unsafe { scan_full(&e, &mut clones, name, seq) };
        }
    }
    let full_dt = t.elapsed().as_secs_f64() / reps as f64;

    // filters-only (MSV->bias->Vit->Fwd), no domaindef. Counts profiles passing F3.
    let mut fwd_pass_total = 0u64;
    let t = Instant::now();
    for _ in 0..reps {
        for (name, seq) in &seqs {
            fwd_pass_total += unsafe { filters_only(&e, &mut clones, name, seq) };
        }
    }
    let filt_dt = t.elapsed().as_secs_f64() / reps as f64;
    let fwd_pass_per_q = fwd_pass_total as f64 / (reps as u64 * nseq) as f64;

    let n = nseq as f64;
    eprintln!("\n=== TASK 1b: full pipeline vs filters-only (single thread, {reps} reps) ===");
    eprintln!("  full pipeline:   {:.3}s  ({:.1} seq/s)", full_dt, n / full_dt);
    eprintln!("  filters-only:    {:.3}s  ({:.1} seq/s)", filt_dt, n / filt_dt);
    eprintln!("  domaindef share of full scan time: {:.1}%", 100.0 * (full_dt - filt_dt) / full_dt);
    eprintln!("  filter-cascade share:              {:.1}%", 100.0 * filt_dt / full_dt);
    eprintln!("  avg profiles passing F3 (filters_only count): {:.3}", fwd_pass_per_q);
}

// Replicate the filter cascade exactly (MSV, bias, Viterbi, Forward), returning
// the number of profiles that pass F3 (i.e. would reach domaindef).
// Faithful port of p7_Pipeline lines up to n_past_fwd++ (SCAN mode).
unsafe fn filters_only(
    e: &Eng,
    clones: &mut [*mut ffi::P7_OPROFILE],
    name: &str,
    seq: &[u8],
) -> u64 {
    let sq = e.make_sq(name, seq);
    let n = (*sq).n as i32;
    let dsq = (*sq).dsq;
    let bg = ffi::p7_bg_Create(e.abc);
    // one-row omx, grow as needed per profile
    let oxf = ffi::p7_omx_Create((*clones[0]).M, 0, n);
    let mut npass = 0u64;
    for &om in clones.iter() {
        // p7_pli_NewModel: configure the bias-filter HMM for this model.
        ffi::p7_bg_SetFilter(bg, (*om).M, (*om).compo.as_ptr());
        ffi::p7_bg_SetLength(bg, n);
        ffi::p7_oprofile_ReconfigLength(om, n);
        ffi::p7_omx_GrowTo(oxf, (*om).M, 0, n);

        let mut nullsc = 0f32;
        ffi::p7_bg_NullOne(bg, dsq, n, &mut nullsc);

        // MSV
        let mut usc = 0f32;
        ffi::p7_MSVFilter(dsq, n, om, oxf, &mut usc);
        let mut seq_score = (usc as f64 - nullsc as f64) / LOG2;
        let mut p = gumbel_surv(
            seq_score,
            (*om).evparam[P7_MMU] as f64,
            (*om).evparam[P7_MLAMBDA] as f64,
        );
        if p > F1 {
            ffi::p7_omx_Reuse(oxf);
            continue;
        }
        // bias filter (do_biasfilter default TRUE)
        let mut filtersc = 0f32;
        ffi::p7_bg_FilterScore(bg, dsq, n, &mut filtersc);
        seq_score = (usc as f64 - filtersc as f64) / LOG2;
        p = gumbel_surv(
            seq_score,
            (*om).evparam[P7_MMU] as f64,
            (*om).evparam[P7_MLAMBDA] as f64,
        );
        if p > F1 {
            ffi::p7_omx_Reuse(oxf);
            continue;
        }
        // SCAN mode: reconfig rest length (profile already fully loaded; hfp NULL)
        ffi::p7_oprofile_ReconfigRestLength(om, n);
        // (NewModelThresholds only matters for use_bit_cutoffs; skip — not used here.)

        // Viterbi
        if p > F2 {
            let mut vfsc = 0f32;
            ffi::p7_ViterbiFilter(dsq, n, om, oxf, &mut vfsc);
            seq_score = (vfsc as f64 - filtersc as f64) / LOG2;
            p = gumbel_surv(
                seq_score,
                (*om).evparam[P7_VMU] as f64,
                (*om).evparam[P7_VLAMBDA] as f64,
            );
            if p > F2 {
                ffi::p7_omx_Reuse(oxf);
                continue;
            }
        }
        // Forward
        let mut fwdsc = 0f32;
        ffi::p7_ForwardParser(dsq, n, om, oxf, &mut fwdsc);
        seq_score = (fwdsc as f64 - filtersc as f64) / LOG2;
        p = exp_surv(
            seq_score,
            (*om).evparam[P7_FTAU] as f64,
            (*om).evparam[P7_FLAMBDA] as f64,
        );
        ffi::p7_omx_Reuse(oxf);
        if p > F3 {
            continue;
        }
        npass += 1;
    }
    ffi::p7_omx_Destroy(oxf);
    ffi::p7_bg_Destroy(bg);
    ffi::esl_sq_Destroy(sq);
    npass
}

// Filter cascade returning, per profile index that passes F3, its Forward
// seq bit-score (a cheap proxy for the final reported bitscore used to rank
// which profiles are worth aligning). Returns Vec<(profile_idx, fwd_bits)>.
unsafe fn filter_scores(
    e: &Eng,
    clones: &mut [*mut ffi::P7_OPROFILE],
    name: &str,
    seq: &[u8],
) -> Vec<(usize, f64)> {
    let sq = e.make_sq(name, seq);
    let n = (*sq).n as i32;
    let dsq = (*sq).dsq;
    let bg = ffi::p7_bg_Create(e.abc);
    let oxf = ffi::p7_omx_Create((*clones[0]).M, 0, n);
    let mut out = Vec::new();
    for (idx, &om) in clones.iter().enumerate() {
        ffi::p7_bg_SetFilter(bg, (*om).M, (*om).compo.as_ptr());
        ffi::p7_bg_SetLength(bg, n);
        ffi::p7_oprofile_ReconfigLength(om, n);
        ffi::p7_omx_GrowTo(oxf, (*om).M, 0, n);
        let mut nullsc = 0f32;
        ffi::p7_bg_NullOne(bg, dsq, n, &mut nullsc);
        let mut usc = 0f32;
        ffi::p7_MSVFilter(dsq, n, om, oxf, &mut usc);
        let mut seq_score = (usc as f64 - nullsc as f64) / LOG2;
        let mut p = gumbel_surv(seq_score, (*om).evparam[P7_MMU] as f64, (*om).evparam[P7_MLAMBDA] as f64);
        if p > F1 { ffi::p7_omx_Reuse(oxf); continue; }
        let mut filtersc = 0f32;
        ffi::p7_bg_FilterScore(bg, dsq, n, &mut filtersc);
        seq_score = (usc as f64 - filtersc as f64) / LOG2;
        p = gumbel_surv(seq_score, (*om).evparam[P7_MMU] as f64, (*om).evparam[P7_MLAMBDA] as f64);
        if p > F1 { ffi::p7_omx_Reuse(oxf); continue; }
        ffi::p7_oprofile_ReconfigRestLength(om, n);
        if p > F2 {
            let mut vfsc = 0f32;
            ffi::p7_ViterbiFilter(dsq, n, om, oxf, &mut vfsc);
            seq_score = (vfsc as f64 - filtersc as f64) / LOG2;
            p = gumbel_surv(seq_score, (*om).evparam[P7_VMU] as f64, (*om).evparam[P7_VLAMBDA] as f64);
            if p > F2 { ffi::p7_omx_Reuse(oxf); continue; }
        }
        let mut fwdsc = 0f32;
        ffi::p7_ForwardParser(dsq, n, om, oxf, &mut fwdsc);
        seq_score = (fwdsc as f64 - filtersc as f64) / LOG2;
        p = exp_surv(seq_score, (*om).evparam[P7_FTAU] as f64, (*om).evparam[P7_FLAMBDA] as f64);
        ffi::p7_omx_Reuse(oxf);
        if p > F3 { continue; }
        out.push((idx, seq_score));
    }
    ffi::p7_omx_Destroy(oxf);
    ffi::p7_bg_Destroy(bg);
    ffi::esl_sq_Destroy(sq);
    out
}

// Run the FULL p7_Pipeline (incl. domaindef) but only on the given subset of
// profile indices. Returns the collected, thresholded HSPs.
unsafe fn scan_subset(
    e: &Eng,
    clones: &mut [*mut ffi::P7_OPROFILE],
    subset: &[usize],
    name: &str,
    seq: &[u8],
) -> Vec<Hsp> {
    let sq = e.make_sq(name, seq);
    let n = (*sq).n as i32;
    let pli = ffi::p7_pipeline_Create(ptr::null(), 100, 100, 0, ffi::p7_pipemodes_e_p7_SCAN_MODELS);
    let bg = ffi::p7_bg_Create(e.abc);
    let th = ffi::p7_tophits_Create();
    for &idx in subset {
        let om = clones[idx];
        ffi::p7_pli_NewModel(pli, om, bg);
        ffi::p7_bg_SetLength(bg, n);
        ffi::p7_oprofile_ReconfigLength(om, n);
        ffi::p7_Pipeline(pli, om, bg, sq, ptr::null_mut(), th);
        ffi::p7_pipeline_Reuse(pli);
    }
    ffi::p7_tophits_SortBySortkey(th);
    ffi::p7_tophits_Threshold(th, pli);
    let hsps = collect_hsps(th);
    ffi::p7_tophits_Destroy(th);
    ffi::p7_bg_Destroy(bg);
    ffi::p7_pipeline_Destroy(pli);
    ffi::esl_sq_Destroy(sq);
    hsps
}

fn species() -> Vec<String> {
    ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"]
        .iter().map(|s| s.to_string()).collect()
}

// Compare state_vectors produced by `hsps` against the fixture for one sequence.
fn sv_matches(hsps: &[Hsp], seq_len: usize, exp_sv: &serde_json::Value) -> bool {
    let sp = species();
    let parsed = parse_hmmer_query(hsps, seq_len, 80.0, Some(&sp));
    let exp = exp_sv.as_array().unwrap();
    if parsed.state_vectors.len() != exp.len() {
        return false;
    }
    for (got, e) in parsed.state_vectors.iter().zip(exp.iter()) {
        let earr = e.as_array().unwrap();
        if got.len() != earr.len() {
            return false;
        }
        for (g, ev) in got.iter().zip(earr.iter()) {
            // ev = [[id, "m"/"i"/"d"], si_or_null]
            let pair = ev.as_array().unwrap();
            let idty = pair[0].as_array().unwrap();
            let eid = idty[0].as_u64().unwrap() as u8;
            let ety = idty[1].as_str().unwrap();
            let esi = if pair[1].is_null() { None } else { Some(pair[1].as_u64().unwrap() as usize) };
            let gty = match g.typ {
                anarci_core::StateType::M => "m",
                anarci_core::StateType::I => "i",
                anarci_core::StateType::D => "d",
            };
            if g.id != eid || gty != ety || g.si != esi {
                return false;
            }
        }
    }
    true
}

// ============================================================================
// TASK 2: score-all, align-winner
// ============================================================================
#[test]
#[ignore]
fn task2_score_all_align_winner() {
    let e = Eng::load();
    let mut clones = e.clones();

    // Load full fixture (need seq_len, state_vectors, and the per-profile hsps).
    let f = std::fs::File::open(fx_path()).unwrap();
    let v: serde_json::Value = serde_json::from_reader(GzDecoder::new(f)).unwrap();
    let seqs = v["sequences"].as_array().unwrap().clone();

    // --- Parity: state_vectors with "align top-K by Forward score" ---
    // For each K in {1,2,3}, run filter cascade, align only the K best by fwd
    // score (mapped via filter_scores -> subset), and check state_vector parity.
    for k in [1usize, 2, 3] {
        let mut sv_ok = 0usize;
        let mut empty_ok = 0usize;
        for s in &seqs {
            let id = s["id"].as_str().unwrap();
            let seq = s["seq"].as_str().unwrap().as_bytes();
            let seq_len = s["seq_len"].as_u64().unwrap() as usize;
            let mut fs = unsafe { filter_scores(&e, &mut clones, id, seq) };
            // rank by fwd bit score descending
            fs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let subset: Vec<usize> = fs.iter().take(k).map(|&(i, _)| i).collect();
            let hsps = unsafe { scan_subset(&e, &mut clones, &subset, id, seq) };
            if sv_matches(&hsps, seq_len, &s["state_vectors"]) {
                sv_ok += 1;
                if s["state_vectors"].as_array().unwrap().is_empty() {
                    empty_ok += 1;
                }
            }
        }
        eprintln!(
            "TASK 2 state_vector parity, align top-{k} by Forward score: {sv_ok}/{} match",
            seqs.len()
        );
        let _ = empty_ok;
    }

    // --- Parity: align ALL profiles passing F3 (== domaindef-only on candidates).
    // This is the "score-all, align-all-survivors" upper bound on state_vector
    // parity (it must equal full scan because the same profiles reach domaindef).
    let mut sv_ok_all = 0usize;
    for s in &seqs {
        let id = s["id"].as_str().unwrap();
        let seq = s["seq"].as_str().unwrap().as_bytes();
        let seq_len = s["seq_len"].as_u64().unwrap() as usize;
        let fs = unsafe { filter_scores(&e, &mut clones, id, seq) };
        let subset: Vec<usize> = fs.iter().map(|&(i, _)| i).collect();
        let hsps = unsafe { scan_subset(&e, &mut clones, &subset, id, seq) };
        if sv_matches(&hsps, seq_len, &s["state_vectors"]) {
            sv_ok_all += 1;
        }
    }
    eprintln!("TASK 2 state_vector parity, align ALL F3-survivors: {sv_ok_all}/{} match", seqs.len());

    // --- hit_table reachability: how many >=80 hits would be MISSING if we only
    // align top-K? For each seq, the fixture's >=80 hit_ids form the required
    // hit_table set. Count how many seqs' full >=80 set is covered by top-K subset.
    for k in [1usize, 2, 3] {
        let mut full_cover = 0usize;
        let mut tot_required = 0usize;
        let mut tot_covered = 0usize;
        for s in &seqs {
            let id = s["id"].as_str().unwrap();
            let seq = s["seq"].as_str().unwrap().as_bytes();
            // required >=80 hit_ids from the fixture's per-profile hsps
            let required: std::collections::HashSet<String> = s["hsps"].as_array().unwrap().iter()
                .filter(|h| h["bitscore"].as_f64().map_or(false, |b| b >= 80.0))
                .map(|h| h["hit_id"].as_str().unwrap().to_string())
                .collect();
            let mut fs = unsafe { filter_scores(&e, &mut clones, id, seq) };
            fs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let subset: Vec<usize> = fs.iter().take(k).map(|&(i, _)| i).collect();
            // map subset profile indices -> hit_ids via om name
            let got_ids: std::collections::HashSet<String> = subset.iter()
                .map(|&i| unsafe { cstr((*clones[i]).name) })
                .collect();
            tot_required += required.len();
            let covered = required.iter().filter(|r| got_ids.contains(*r)).count();
            tot_covered += covered;
            if covered == required.len() {
                full_cover += 1;
            }
        }
        eprintln!(
            "TASK 2 hit_table coverage, top-{k}: {full_cover}/{} seqs fully cover their >=80 set; \
             {tot_covered}/{tot_required} required >=80 rows present",
            seqs.len()
        );
    }

    // --- Timing: full scan vs score-all+align-top1 vs score-all+align-all-F3 ---
    let plain: Vec<(String, Vec<u8>)> = seqs.iter()
        .map(|s| (s["id"].as_str().unwrap().to_string(), s["seq"].as_str().unwrap().as_bytes().to_vec()))
        .collect();
    let reps = 3usize;
    // warm
    for (id, sq) in &plain { unsafe { scan_full(&e, &mut clones, id, sq) }; }

    let t = Instant::now();
    for _ in 0..reps { for (id, sq) in &plain { let _ = unsafe { scan_full(&e, &mut clones, id, sq) }; } }
    let full_dt = t.elapsed().as_secs_f64() / reps as f64;

    let t = Instant::now();
    for _ in 0..reps {
        for (id, sq) in &plain {
            let mut fs = unsafe { filter_scores(&e, &mut clones, id, sq) };
            fs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let subset: Vec<usize> = fs.iter().take(1).map(|&(i, _)| i).collect();
            let _ = unsafe { scan_subset(&e, &mut clones, &subset, id, sq) };
        }
    }
    let top1_dt = t.elapsed().as_secs_f64() / reps as f64;

    let t = Instant::now();
    for _ in 0..reps {
        for (id, sq) in &plain {
            let fs = unsafe { filter_scores(&e, &mut clones, id, sq) };
            let subset: Vec<usize> = fs.iter().map(|&(i, _)| i).collect();
            let _ = unsafe { scan_subset(&e, &mut clones, &subset, id, sq) };
        }
    }
    let allf3_dt = t.elapsed().as_secs_f64() / reps as f64;

    let n = plain.len() as f64;
    eprintln!("\n=== TASK 2 timing (single thread, {reps} reps) ===");
    eprintln!("  full scan (29 profiles, domaindef each): {:.3}s  ({:.1} seq/s)", full_dt, n / full_dt);
    eprintln!("  score-all + align top-1:                 {:.3}s  ({:.1} seq/s)  speedup {:.2}x",
        top1_dt, n / top1_dt, full_dt / top1_dt);
    eprintln!("  score-all + align all F3-survivors:      {:.3}s  ({:.1} seq/s)  speedup {:.2}x",
        allf3_dt, n / allf3_dt, full_dt / allf3_dt);
}

// ============================================================================
// TASK 3: other lossless wins
// ============================================================================
#[test]
#[ignore]
fn task3_lossless_wins() {
    let e = Eng::load();
    let mut clones = e.clones();
    let seqs = load_seqs();
    let n = seqs.len();

    // --- 3a: dedup ---
    let uniq: std::collections::HashMap<&[u8], &str> = seqs.iter()
        .map(|(id, s)| (s.as_slice(), id.as_str()))
        .collect();
    let n_uniq = uniq.len();
    eprintln!("\n=== TASK 3a: input dedup ===");
    eprintln!("  {n} sequences, {n_uniq} unique ({:.1}% redundant). Scans avoided: {} ({:.1}%).",
        100.0 * (n - n_uniq) as f64 / n as f64, n - n_uniq, 100.0 * (n - n_uniq) as f64 / n as f64);

    // timing: full scan all 996 vs scan-unique + replay
    let reps = 3usize;
    for (id, s) in &seqs { unsafe { scan_full(&e, &mut clones, id, s) }; }
    let t = Instant::now();
    for _ in 0..reps { for (id, s) in &seqs { let _ = unsafe { scan_full(&e, &mut clones, id, s) }; } }
    let all_dt = t.elapsed().as_secs_f64() / reps as f64;

    let unique_seqs: Vec<(&[u8], &str)> = uniq.iter().map(|(k, v)| (*k, *v)).collect();
    let t = Instant::now();
    for _ in 0..reps {
        for (s, id) in &unique_seqs { let _ = unsafe { scan_full(&e, &mut clones, id, s) }; }
    }
    let uniq_dt = t.elapsed().as_secs_f64() / reps as f64;
    eprintln!("  scan all 996:    {:.3}s  ({:.1} seq/s effective over 996)", all_dt, n as f64 / all_dt);
    eprintln!("  scan {n_uniq} unique: {:.3}s  -> effective {:.1} seq/s over 996  speedup {:.2}x",
        uniq_dt, n as f64 / uniq_dt, all_dt / uniq_dt);

    // --- 3b: pipeline reuse across scans with manual nmodels/Z reset ---
    // Reuse ONE pipeline + bg + tophits across all sequences, resetting the
    // counters that accumulate in SCAN mode (nmodels via Z), and verify HSP parity.
    test_pipeline_reuse(&e, &mut clones, &seqs);
}

// Reuse a single pipeline/bg across all scans. In SCAN mode p7_pli_NewModel does
// pli->nmodels++ and pli->Z = pli->nmodels, so Z grows past 29 across reused
// scans, corrupting E-values. We reset nmodels=0 (hence Z=29 after the 29-model
// loop) at the start of each sequence and check whether HSPs stay identical to a
// fresh-pipeline-per-scan baseline.
fn test_pipeline_reuse(e: &Eng, clones: &mut [*mut ffi::P7_OPROFILE], seqs: &[(String, Vec<u8>)]) {
    eprintln!("\n=== TASK 3b: per-thread pipeline reuse with manual nmodels reset ===");
    // baseline: fresh pipeline per scan (production behavior)
    let baseline: Vec<Vec<Hsp>> = seqs.iter()
        .map(|(id, s)| unsafe { scan_full(e, clones, id, s) })
        .collect();

    // reuse path
    let pli = unsafe {
        ffi::p7_pipeline_Create(ptr::null(), 100, 100, 0, ffi::p7_pipemodes_e_p7_SCAN_MODELS)
    };
    let bg = unsafe { ffi::p7_bg_Create(e.abc) };
    let mut mism = 0usize;
    let mut checked = 0usize;
    let reps = 3usize;
    let t = Instant::now();
    for rep in 0..reps {
        for (qi, (id, s)) in seqs.iter().enumerate() {
            let sq = e.make_sq(id, s);
            let n = unsafe { (*sq).n } as i32;
            let th = unsafe { ffi::p7_tophits_Create() };
            // manual reset of accumulating counters before the 29-model loop
            unsafe {
                (*pli).nmodels = 0;
                (*pli).nseqs = 0;
                (*pli).n_past_msv = 0;
                (*pli).n_past_bias = 0;
                (*pli).n_past_vit = 0;
                (*pli).n_past_fwd = 0;
            }
            for &om in clones.iter() {
                unsafe {
                    ffi::p7_pli_NewModel(pli, om, bg);
                    ffi::p7_bg_SetLength(bg, n);
                    ffi::p7_oprofile_ReconfigLength(om, n);
                    ffi::p7_Pipeline(pli, om, bg, sq, ptr::null_mut(), th);
                    ffi::p7_pipeline_Reuse(pli);
                }
            }
            unsafe {
                ffi::p7_tophits_SortBySortkey(th);
                ffi::p7_tophits_Threshold(th, pli);
            }
            let hsps = unsafe { collect_hsps(th) };
            if rep == 0 {
                checked += 1;
                if !hsps_eq(&hsps, &baseline[qi]) {
                    mism += 1;
                }
            }
            unsafe {
                ffi::p7_tophits_Destroy(th);
                ffi::esl_sq_Destroy(sq);
            }
        }
    }
    let reuse_dt = t.elapsed().as_secs_f64() / reps as f64;
    unsafe {
        ffi::p7_bg_Destroy(bg);
        ffi::p7_pipeline_Destroy(pli);
    }

    // baseline timing
    let t = Instant::now();
    for _ in 0..reps { for (id, s) in seqs { let _ = unsafe { scan_full(e, clones, id, s) }; } }
    let base_dt = t.elapsed().as_secs_f64() / reps as f64;

    let n = seqs.len() as f64;
    eprintln!("  HSP parity (reuse vs fresh-per-scan): {}/{} sequences identical", checked - mism, checked);
    eprintln!("  fresh pipeline per scan: {:.3}s  ({:.1} seq/s)", base_dt, n / base_dt);
    eprintln!("  reused pipeline (manual reset): {:.3}s  ({:.1} seq/s)  speedup {:.2}x",
        reuse_dt, n / reuse_dt, base_dt / reuse_dt);
}

fn hsps_eq(a: &[Hsp], b: &[Hsp]) -> bool {
    if a.len() != b.len() { return false; }
    for (x, y) in a.iter().zip(b.iter()) {
        if x.hit_id != y.hit_id
            || x.bitscore != y.bitscore
            || x.evalue != y.evalue
            || x.bias != y.bias
            || x.query_start != y.query_start
            || x.query_end != y.query_end
            || x.hit_start != y.hit_start
            || x.hit_end != y.hit_end
            || x.rf != y.rf
            || x.pp != y.pp
        {
            return false;
        }
    }
    true
}

// ============================================================================
// TASK 2b: Can a Forward-score cutoff losslessly pre-select a superset of the
// >=80 hits (so we run domaindef on fewer than 28.5 profiles but never drop a
// required row)? Measures the separation between Forward seq-score of profiles
// that END UP >=80 vs those that end <80, and what cutoff keeps 100% of >=80.
// ============================================================================
#[test]
#[ignore]
fn task2b_forward_cutoff_for_hit_table() {
    let e = Eng::load();
    let mut clones = e.clones();
    let f = std::fs::File::open(fx_path()).unwrap();
    let v: serde_json::Value = serde_json::from_reader(GzDecoder::new(f)).unwrap();
    let seqs = v["sequences"].as_array().unwrap();

    // For every (seq, profile-that-passes-F3): collect (fwd_score, final_bitscore).
    // final_bitscore from the fixture's hsps by hit_id (profiles producing no
    // reported hit have no fixture row -> treat as final bitscore = -inf).
    let mut min_fwd_for_ge80 = f64::INFINITY; // lowest fwd score among profiles that end >=80
    let mut max_fwd_for_lt80 = f64::NEG_INFINITY; // highest fwd score among profiles that end <80
    let mut n_ge80 = 0usize;
    let mut n_lt80 = 0usize;
    // count avg F3-survivors and avg >=80 per query, and the avg # survivors with
    // fwd_score >= (min_fwd_for_ge80) computed in a 2nd pass.
    let mut all: Vec<(f64, bool)> = Vec::new(); // (fwd_score, ends_ge80)
    let mut tot_survivors = 0usize;

    for s in seqs {
        let id = s["id"].as_str().unwrap();
        let seq = s["seq"].as_str().unwrap().as_bytes();
        // map hit_id -> max final bitscore in fixture
        let mut final_bs: std::collections::HashMap<String, f64> = std::collections::HashMap::new();
        for h in s["hsps"].as_array().unwrap() {
            let hid = h["hit_id"].as_str().unwrap().to_string();
            let bs = h["bitscore"].as_f64().unwrap_or(f64::NEG_INFINITY);
            let ent = final_bs.entry(hid).or_insert(f64::NEG_INFINITY);
            if bs > *ent { *ent = bs; }
        }
        let fs = unsafe { filter_scores(&e, &mut clones, id, seq) };
        tot_survivors += fs.len();
        for (idx, fwd) in fs {
            let pname = unsafe { cstr((*clones[idx]).name) };
            let fb = final_bs.get(&pname).copied().unwrap_or(f64::NEG_INFINITY);
            let ge80 = fb >= 80.0;
            all.push((fwd, ge80));
            if ge80 {
                n_ge80 += 1;
                if fwd < min_fwd_for_ge80 { min_fwd_for_ge80 = fwd; }
            } else {
                n_lt80 += 1;
                if fwd > max_fwd_for_lt80 { max_fwd_for_lt80 = fwd; }
            }
        }
    }

    let nseq = seqs.len() as f64;
    eprintln!("\n=== TASK 2b: Forward-score cutoff vs final >=80 ===");
    eprintln!("  avg F3-survivors/query (need full pipeline today): {:.2}", tot_survivors as f64 / nseq);
    eprintln!("  avg profiles ending >=80/query (hit_table rows):   {:.2}", n_ge80 as f64 / nseq);
    eprintln!("  Forward seq-score: min among >=80 profiles = {:.2} bits", min_fwd_for_ge80);
    eprintln!("  Forward seq-score: max among <80  profiles = {:.2} bits", max_fwd_for_lt80);
    eprintln!("  separable by a single global cutoff? {}",
        if min_fwd_for_ge80 > max_fwd_for_lt80 { "YES" } else { "NO (overlap)" });

    // If we set the cutoff at min_fwd_for_ge80 (keeps ALL >=80 losslessly), how
    // many profiles per query would we still domaindef?
    let cutoff = min_fwd_for_ge80;
    let kept: usize = all.iter().filter(|(f, _)| *f >= cutoff).count();
    eprintln!("  with lossless cutoff {:.2}: domaindef on {:.2} profiles/query (vs {:.2} F3-survivors)",
        cutoff, kept as f64 / nseq, tot_survivors as f64 / nseq);
    eprintln!("  (these are >=80-superset counts; <80 profiles above cutoff: {})",
        kept - n_ge80);
    let _ = (n_lt80,);
}
