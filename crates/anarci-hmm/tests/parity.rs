//! Primary gate: the in-process HMMER engine reproduces ANARCI's hmmscan.
//!
//! For each of the 996 fixture sequences we run `Engine::scan`, feed the
//! resulting `Vec<Hsp>` through `anarci_core::parse_hmmer_query` (bit-score
//! threshold 80, species human/mouse/rat/rabbit/rhesus/pig/alpaca), and assert
//! the produced state vectors equal the fixture's exactly (id, type, si per
//! element). We also report bitscore / evalue / query_start / query_end
//! agreement of the resulting `details` against the fixture.

use anarci_core::{parse_hmmer_query, State, StateType};
use anarci_hmm::Engine;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Fixture {
    sequences: Vec<SeqJson>,
}
#[derive(Deserialize)]
struct SeqJson {
    id: String,
    seq: String,
    seq_len: usize,
    state_vectors: Vec<Vec<((u8, String), Option<usize>)>>,
    details: Vec<DetailJson>,
}
#[derive(Deserialize)]
struct DetailJson {
    species: Option<String>,
    chain_type: Option<String>,
    bitscore: Option<f64>,
    evalue: Option<f64>,
    query_start: Option<i64>,
    query_end: Option<usize>,
}

fn species() -> Vec<String> {
    ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"]
        .iter()
        .map(|s| s.to_string())
        .collect()
}

fn tych(t: StateType) -> &'static str {
    match t {
        StateType::M => "m",
        StateType::I => "i",
        StateType::D => "d",
    }
}

fn hmm_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference_data/dat/HMMs/ALL.hmm")
}

/// Outcome of one sequence: ok flag, optional state-vector diff message, and
/// the per-domain details-agreement counters.
struct Outcome {
    id: String,
    sv_ok: bool,
    diff: Option<String>,
    n_details: usize,
    bitscore_ok: usize,
    evalue_ok: usize,
    qstart_ok: usize,
    qend_ok: usize,
}

#[test]
fn gate_engine_state_vector_parity() {
    let fx_path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/hsps.json.gz");
    let f = std::fs::File::open(&fx_path).unwrap_or_else(|e| panic!("open {:?}: {}", fx_path, e));
    let fx: Fixture = serde_json::from_reader(GzDecoder::new(f)).expect("parse hsps.json.gz");

    let engine = Engine::load(&hmm_path()).expect("load ALL.hmm");
    assert_eq!(engine.n_models(), 29, "expected 29 profiles in ALL.hmm");
    let sp = species();

    // Scan every sequence (in parallel; Engine is Send + Sync) and compare.
    let outcomes: Vec<Outcome> = fx
        .sequences
        .par_iter()
        .map(|s| {
            let hsps = engine.scan(&s.id, s.seq.as_bytes());
            let parsed = parse_hmmer_query(&hsps, s.seq_len, 80.0, Some(&sp));

            let mut diff: Option<String> = None;
            let mut bitscore_ok = 0;
            let mut evalue_ok = 0;
            let mut qstart_ok = 0;
            let mut qend_ok = 0;
            let mut n_details = 0;

            if parsed.state_vectors.len() != s.state_vectors.len() {
                diff = Some(format!(
                    "ndomains {} vs {}",
                    parsed.state_vectors.len(),
                    s.state_vectors.len()
                ));
            } else {
                'outer: for (di, (got_sv, exp_sv)) in parsed
                    .state_vectors
                    .iter()
                    .zip(s.state_vectors.iter())
                    .enumerate()
                {
                    if got_sv.len() != exp_sv.len() {
                        diff =
                            Some(format!("dom {di} sv len {} vs {}", got_sv.len(), exp_sv.len()));
                        break;
                    }
                    for (g, ((eid, ety), esi)) in got_sv.iter().zip(exp_sv.iter()) {
                        let &State { id, typ, si } = g;
                        if id != *eid || tych(typ) != ety.as_str() || si != *esi {
                            diff = Some(format!(
                                "dom {di}: got ({id},{},{:?}) want ({eid},{ety},{esi:?})",
                                tych(typ),
                                si
                            ));
                            break 'outer;
                        }
                    }

                    // details agreement (reported, not asserted individually
                    // except as parity counters)
                    let d = &parsed.details[di];
                    let e = &s.details[di];
                    n_details += 1;
                    // species/chain_type must match exactly (part of the gate).
                    assert_eq!(Some(d.species.as_str()), e.species.as_deref(), "{} species", s.id);
                    assert_eq!(
                        Some(d.chain_type.as_str()),
                        e.chain_type.as_deref(),
                        "{} chain_type",
                        s.id
                    );
                    if e.bitscore.map_or(false, |b| (d.bitscore - b).abs() < 1e-9) {
                        bitscore_ok += 1;
                    }
                    if e.evalue.map_or(false, |ev| {
                        (d.evalue - ev).abs() <= ev.abs() * 1e-12 + 1e-300
                    }) {
                        evalue_ok += 1;
                    }
                    if d.query_start.map(|v| v as i64) == e.query_start {
                        qstart_ok += 1;
                    }
                    if Some(d.query_end) == e.query_end {
                        qend_ok += 1;
                    }
                }
            }

            Outcome {
                id: s.id.clone(),
                sv_ok: diff.is_none(),
                diff,
                n_details,
                bitscore_ok,
                evalue_ok,
                qstart_ok,
                qend_ok,
            }
        })
        .collect();

    let mut ok = 0u64;
    let mut bad = 0u64;
    let mut samples: Vec<String> = Vec::new();
    let (mut tot_details, mut tot_bit, mut tot_ev, mut tot_qs, mut tot_qe) = (0, 0, 0, 0, 0);

    for o in &outcomes {
        if o.sv_ok {
            ok += 1;
        } else {
            bad += 1;
            if let Some(m) = &o.diff {
                if samples.len() < 25 {
                    samples.push(format!("{}: {m}", o.id));
                }
            }
        }
        tot_details += o.n_details;
        tot_bit += o.bitscore_ok;
        tot_ev += o.evalue_ok;
        tot_qs += o.qstart_ok;
        tot_qe += o.qend_ok;
    }

    eprintln!("\n=== Engine parity gate: ok={ok} bad={bad} (of {}) ===", outcomes.len());
    eprintln!(
        "details agreement over {tot_details} chosen domains: \
         bitscore {tot_bit}/{tot_details}  evalue {tot_ev}/{tot_details}  \
         query_start {tot_qs}/{tot_details}  query_end {tot_qe}/{tot_details}"
    );
    for s in &samples {
        eprintln!("  {s}");
    }
    assert_eq!(bad, 0, "{bad} sequences mismatch reference state vectors");
}
