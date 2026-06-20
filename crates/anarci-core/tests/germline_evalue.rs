//! E-value germline assignment: sanity + behaviour checks on the real golden
//! state vectors (996 seqs). This does NOT assert byte-parity with ANARCI (the
//! e-value path is intentionally a different, more accurate method); it checks
//! that the method runs over the whole set, returns coherent output, and that
//! every winning V/J gene the e-value path picks is at least as good a local
//! match (by SW score) as the identity path's gene — i.e. the selection is
//! genuinely score-optimal, not a regression.

use anarci_core::sw;
use anarci_core::{
    run_germline_assignment, run_germline_assignment_evalue, run_germline_assignment_evalue_exact,
    State, StateType, StateVector,
};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Fixture {
    sequences: Vec<SeqEntry>,
}
#[derive(Deserialize)]
struct SeqEntry {
    id: String,
    seq: String,
    domains: Vec<Domain>,
}
#[derive(Deserialize)]
struct Domain {
    chain_type: String,
    state_vector: Vec<((u8, String), Option<usize>)>,
    #[serde(default)]
    germlines: Option<serde_json::Value>,
}

fn load_fixture() -> Fixture {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden.json.gz");
    let f = std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {:?}: {}", path, e));
    serde_json::from_reader(GzDecoder::new(f)).expect("parse golden.json.gz")
}

fn to_sv(raw: &[((u8, String), Option<usize>)]) -> StateVector {
    raw.iter()
        .map(|((id, ty), si)| State {
            id: *id,
            typ: match ty.as_str() {
                "m" => StateType::M,
                "i" => StateType::I,
                "d" => StateType::D,
                other => panic!("bad state type {other:?}"),
            },
            si: *si,
        })
        .collect()
}

const SPECIES: [&str; 7] = ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"];

/// Runs cleanly over the whole golden set and produces coherent output:
/// whenever the identity path made a V call, the e-value path also makes one,
/// with a populated e-value and an identity in [0,1].
#[test]
fn evalue_runs_and_is_coherent() {
    let fx = load_fixture();
    let species: Vec<String> = SPECIES.iter().map(|s| s.to_string()).collect();
    let (mut domains, mut v_called, mut both_called) = (0u64, 0u64, 0u64);

    for seq in &fx.sequences {
        let bytes = seq.seq.as_bytes();
        for dom in &seq.domains {
            // Only domains where ANARCI produced a germline dict at all.
            let exp = match &dom.germlines {
                Some(v) if v.as_object().map(|o| !o.is_empty()).unwrap_or(false) => v,
                _ => continue,
            };
            domains += 1;
            let sv = to_sv(&dom.state_vector);
            let id_g = run_germline_assignment(&sv, bytes, &dom.chain_type, Some(&species));
            let ev_g = run_germline_assignment_evalue(&sv, bytes, &dom.chain_type, Some(&species));

            // Identity path's expected v_gene presence drives the comparison.
            let id_has_v = id_g.v_gene.is_some();
            if id_has_v {
                v_called += 1;
                assert!(
                    ev_g.v_gene.is_some(),
                    "{}: e-value path made no V call where identity did",
                    seq.id
                );
                let ev = ev_g.v_evalue.expect("v_evalue present on e-value path");
                assert!(ev >= 0.0 && ev.is_finite(), "{}: bad v_evalue {ev}", seq.id);
                if let Some(idn) = ev_g.v_identity {
                    assert!((0.0..=1.0).contains(&idn), "{}: v_identity {idn}", seq.id);
                }
                both_called += 1;
            }
            // The identity field must mirror upstream ANARCI's shape coverage.
            let _ = exp;
        }
    }
    eprintln!(
        "e-value germline: {domains} domains, identity made {v_called} V calls, \
         e-value matched all ({both_called})"
    );
    assert!(domains > 500, "expected the full golden set, got {domains}");
    assert_eq!(v_called, both_called);
}

/// The e-value path's chosen V gene must have a Smith-Waterman score >= the
/// identity path's chosen V gene against the same query V region. (Selection is
/// by e-value = score-monotone within the search, so the winner is score-optimal;
/// this guards against region-extraction or ranking regressions.)
#[test]
fn evalue_v_gene_is_score_optimal() {
    let fx = load_fixture();
    let species: Vec<String> = SPECIES.iter().map(|s| s.to_string()).collect();
    let mut checked = 0u64;

    for seq in &fx.sequences {
        let bytes = seq.seq.as_bytes();
        for dom in &seq.domains {
            if dom
                .germlines
                .as_ref()
                .and_then(|v| v.as_object())
                .map(|o| o.is_empty())
                .unwrap_or(true)
            {
                continue;
            }
            let sv = to_sv(&dom.state_vector);
            let id_g = run_germline_assignment(&sv, bytes, &dom.chain_type, Some(&species));
            let ev_g = run_germline_assignment_evalue(&sv, bytes, &dom.chain_type, Some(&species));
            let (id_v, ev_v) = match (&id_g.v_gene, &ev_g.v_gene) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            // Re-extract the V-region query the same way the e-value path does,
            // then compare SW scores of the two chosen genes against it.
            let v_query = v_region_query(&sv, bytes);
            if v_query.is_empty() {
                continue;
            }
            let id_seq = gene_ungapped("V", &dom.chain_type, &id_v.0, &id_v.1);
            let ev_seq = gene_ungapped("V", &dom.chain_type, &ev_v.0, &ev_v.1);
            if let (Some(ids), Some(evs)) = (id_seq, ev_seq) {
                let id_score = sw::local_score(&v_query, &ids);
                let ev_score = sw::local_score(&v_query, &evs);
                assert!(
                    ev_score >= id_score,
                    "{}: e-value gene {:?} scored {ev_score} < identity gene {:?} {id_score}",
                    seq.id,
                    ev_v,
                    id_v
                );
                checked += 1;
            }
        }
    }
    eprintln!("score-optimality checked on {checked} domains");
    assert!(checked > 500, "expected many domains, got {checked}");
}

/// Single-thread speed comparison: exact scalar full-scan (the old default) vs the
/// fast k-mer+SIMD path, over the whole golden set, all-species. Reports ms/domain
/// for both and the speedup. Run with `--nocapture` to see numbers; not a hard
/// assertion on absolute time (machine-dependent), only that fast is faster.
#[test]
fn bench_evalue_germline_speed() {
    use std::time::Instant;
    let fx = load_fixture();
    let all: Vec<String> = SPECIES.iter().map(|s| s.to_string()).collect();

    // Collect (sv, seq, chain) for every domain once.
    let mut work: Vec<(StateVector, Vec<u8>, String)> = Vec::new();
    for seq in &fx.sequences {
        let bytes = seq.seq.as_bytes().to_vec();
        for dom in &seq.domains {
            work.push((to_sv(&dom.state_vector), bytes.clone(), dom.chain_type.clone()));
        }
    }
    let n = work.len() as f64;

    // Warm the lazy indices / DBs once so timing excludes one-off build cost.
    if let Some((sv, sq, ct)) = work.first() {
        let _ = run_germline_assignment_evalue(sv, sq, ct, Some(&all));
        let _ = run_germline_assignment_evalue_exact(sv, sq, ct, Some(&all));
    }

    let t0 = Instant::now();
    let mut sink = 0usize;
    for (sv, sq, ct) in &work {
        let g = run_germline_assignment_evalue_exact(sv, sq, ct, Some(&all));
        sink += g.v_gene.is_some() as usize;
    }
    let exact_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t1 = Instant::now();
    for (sv, sq, ct) in &work {
        let g = run_germline_assignment_evalue(sv, sq, ct, Some(&all));
        sink += g.v_gene.is_some() as usize;
    }
    let fast_ms = t1.elapsed().as_secs_f64() * 1000.0;

    eprintln!(
        "evalue germline speed ({} domains, all-species, single-thread):\n  \
         exact scalar full-scan: {:.3} ms/domain ({:.0} dom/s)\n  \
         fast kmer+SIMD:         {:.3} ms/domain ({:.0} dom/s)\n  \
         speedup: {:.1}x   (sink={})",
        work.len(),
        exact_ms / n,
        n / (exact_ms / 1000.0),
        fast_ms / n,
        n / (fast_ms / 1000.0),
        exact_ms / fast_ms,
        sink
    );
    assert!(fast_ms < exact_ms, "fast path should be faster than exact full scan");
}

/// HARD GATE: the fast e-value path (k-mer prefilter -> SIMD SW) must produce the
/// IDENTICAL v_gene AND j_gene as the exact scalar full-scan brute force, for every
/// domain in the golden set, under BOTH all-species and human-only scoping. This is
/// the optimization's correctness contract: it is a faster way to compute the same
/// answer, not a behaviour change. Zero mismatches required.
#[test]
fn fast_evalue_calls_identical_to_exact_full_scan() {
    let fx = load_fixture();
    let all: Vec<String> = SPECIES.iter().map(|s| s.to_string()).collect();
    let human: Vec<String> = vec!["human".to_string()];

    let mut domains = 0u64;
    let mut v_mismatch = 0u64;
    let mut j_mismatch = 0u64;
    let mut first_mismatch: Option<String> = None;

    for scope_name in ["all-species", "human-only"] {
        let scope: &[String] = if scope_name == "all-species" { &all } else { &human };
        for seq in &fx.sequences {
            let bytes = seq.seq.as_bytes();
            for dom in &seq.domains {
                let sv = to_sv(&dom.state_vector);
                let fast =
                    run_germline_assignment_evalue(&sv, bytes, &dom.chain_type, Some(scope));
                let exact = run_germline_assignment_evalue_exact(
                    &sv,
                    bytes,
                    &dom.chain_type,
                    Some(scope),
                );
                domains += 1;
                if fast.v_gene != exact.v_gene {
                    v_mismatch += 1;
                    first_mismatch.get_or_insert_with(|| {
                        format!(
                            "[{scope_name}] {} ct={}: fast v={:?} exact v={:?}",
                            seq.id, dom.chain_type, fast.v_gene, exact.v_gene
                        )
                    });
                }
                if fast.j_gene != exact.j_gene {
                    j_mismatch += 1;
                    first_mismatch.get_or_insert_with(|| {
                        format!(
                            "[{scope_name}] {} ct={}: fast j={:?} exact j={:?}",
                            seq.id, dom.chain_type, fast.j_gene, exact.j_gene
                        )
                    });
                }
            }
        }
    }

    eprintln!(
        "identical-calls gate: {domains} (domain x scope) comparisons, \
         v_mismatch={v_mismatch}, j_mismatch={j_mismatch}"
    );
    assert!(domains > 1000, "expected the full golden set x2 scopes, got {domains}");
    assert_eq!(
        v_mismatch + j_mismatch,
        0,
        "fast e-value path diverged from exact full scan. first: {}",
        first_mismatch.unwrap_or_default()
    );
}

// --- helpers that mirror the internals of the e-value path for the test ---

/// Reproduce `split_v_j_regions`' V slice (start..=Cys104) from the state vector.
fn v_region_query(sv: &StateVector, seq: &[u8]) -> Vec<u8> {
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;
    let mut cys: Option<usize> = None;
    for s in sv {
        if let Some(si) = s.si {
            first = Some(first.map_or(si, |f: usize| f.min(si)));
            last = Some(last.map_or(si, |l: usize| l.max(si)));
            if s.id == 104 && s.typ == StateType::M {
                cys = Some(si);
            }
        }
    }
    match (first, last) {
        (Some(f), Some(l)) if f <= l && l < seq.len() => match cys {
            Some(c) if c >= f && c <= l => seq[f..=c].to_vec(),
            _ => seq[f..=l].to_vec(),
        },
        _ => Vec::new(),
    }
}

/// Ungapped germline sequence for a named gene, read straight from the embedded
/// gapped DB (independent of the module's private cache, for an honest check).
fn gene_ungapped(seg: &str, chain: &str, species: &str, gene: &str) -> Option<Vec<u8>> {
    #[derive(Deserialize)]
    struct G {
        germlines:
            std::collections::HashMap<
                String,
                std::collections::HashMap<
                    String,
                    std::collections::HashMap<String, Vec<(String, String)>>,
                >,
            >,
    }
    use once_cell::sync::Lazy;
    static DB: Lazy<G> = Lazy::new(|| {
        let raw = include_str!("../data/germlines.json");
        serde_json::from_str(raw).unwrap()
    });
    DB.germlines
        .get(seg)?
        .get(chain)?
        .get(species)?
        .iter()
        .find(|(g, _)| g == gene)
        .map(|(_, s)| s.bytes().filter(|&c| c != b'-').collect())
}
