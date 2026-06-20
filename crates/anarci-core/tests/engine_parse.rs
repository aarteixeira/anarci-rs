//! Gate A: the HMMER-output → state-vector transforms reproduce reference ANARCI
//! from the SAME HSPs. Oracle: tests/fixtures/hsps.json.gz (42k HSPs, 996 seqs),
//! captured by instrumenting `_parse_hmmer_query`. State vectors here are
//! pre-`check_for_j`.

use anarci_core::{parse_hmmer_query, Hsp, State, StateType};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Fixture {
    sequences: Vec<SeqJson>,
}
#[derive(Deserialize)]
struct SeqJson {
    id: String,
    seq_len: usize,
    hsps: Vec<HspJson>,
    state_vectors: Vec<Vec<((u8, String), Option<usize>)>>,
    details: Vec<DetailJson>,
}
#[derive(Deserialize)]
struct HspJson {
    hit_id: String,
    hit_description: String,
    evalue: f64,
    bitscore: f64,
    bias: f64,
    query_start: usize,
    query_end: usize,
    hit_start: usize,
    hit_end: usize,
    rf: String,
    pp: String,
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

#[test]
fn gate_a_parse_and_state_vectors() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/hsps.json.gz");
    let f = std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {:?}: {}", path, e));
    let fx: Fixture = serde_json::from_reader(GzDecoder::new(f)).expect("parse hsps.json.gz");
    let sp = species();

    let mut ok = 0u64;
    let mut bad = 0u64;
    let mut samples: Vec<String> = Vec::new();

    for s in &fx.sequences {
        let hsps: Vec<Hsp> = s
            .hsps
            .iter()
            .map(|h| Hsp {
                hit_id: h.hit_id.clone(),
                hit_description: h.hit_description.clone(),
                evalue: h.evalue,
                bitscore: h.bitscore,
                bias: h.bias,
                query_start: h.query_start,
                query_end: h.query_end,
                hit_start: h.hit_start,
                hit_end: h.hit_end,
                rf: h.rf.clone(),
                pp: h.pp.clone(),
                order: 0,
            })
            .collect();

        let parsed = parse_hmmer_query(&hsps, s.seq_len, 80.0, Some(&sp));

        let mut diff: Option<String> = None;
        if parsed.state_vectors.len() != s.state_vectors.len() {
            diff = Some(format!(
                "ndomains {} vs {}",
                parsed.state_vectors.len(),
                s.state_vectors.len()
            ));
        } else {
            'outer: for (di, (got_sv, exp_sv)) in
                parsed.state_vectors.iter().zip(s.state_vectors.iter()).enumerate()
            {
                if got_sv.len() != exp_sv.len() {
                    diff = Some(format!("dom {di} sv len {} vs {}", got_sv.len(), exp_sv.len()));
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
                // details
                let d = &parsed.details[di];
                let e = &s.details[di];
                let qs = d.query_start.map(|v| v as i64);
                if Some(d.species.as_str()) != e.species.as_deref()
                    || Some(d.chain_type.as_str()) != e.chain_type.as_deref()
                    || (d.bitscore - e.bitscore.unwrap()).abs() > 1e-9
                    || (d.evalue - e.evalue.unwrap()).abs() > (e.evalue.unwrap().abs() * 1e-12 + 1e-30)
                    || qs != e.query_start
                    || Some(d.query_end) != e.query_end
                {
                    diff = Some(format!(
                        "dom {di} details: got (sp={},ct={},bs={},qs={:?},qe={}) want (sp={:?},ct={:?},bs={:?},qs={:?},qe={:?})",
                        d.species, d.chain_type, d.bitscore, qs, d.query_end,
                        e.species, e.chain_type, e.bitscore, e.query_start, e.query_end
                    ));
                    break;
                }
            }
        }

        if let Some(m) = diff {
            bad += 1;
            if samples.len() < 25 {
                samples.push(format!("{}: {m}", s.id));
            }
        } else {
            ok += 1;
        }
    }

    eprintln!("\n=== Gate A (parse + state vectors): ok={ok} bad={bad} ===");
    for s in &samples {
        eprintln!("  {s}");
    }
    assert_eq!(bad, 0, "{bad} sequences mismatch reference parse");
}
