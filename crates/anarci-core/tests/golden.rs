//! Phase-1 correctness gate: reproduce reference ANARCI numbering + germline
//! assignment byte-for-byte from the same input state vectors.
//!
//! Oracle: tests/fixtures/golden.json.gz (996 seqs, 1013 domains, 6 schemes),
//! captured from conda `anarci 2024.05.21`.

use anarci_core::schemes::number_sequence_from_alignment;
use anarci_core::{run_germline_assignment, State, StateType, StateVector};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::BTreeMap;
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
    numbering: BTreeMap<String, NumberingJson>,
}
#[derive(Deserialize)]
struct NumberingJson {
    #[serde(default)]
    numbering: Option<Vec<((i32, String), String)>>,
    #[serde(default)]
    start: Option<i64>,
    #[serde(default)]
    end: Option<i64>,
    #[serde(default)]
    error: Option<String>,
}

fn load_fixture() -> Fixture {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/golden.json.gz");
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

const SCHEMES: [&str; 6] = ["imgt", "chothia", "kabat", "martin", "aho", "wolfguy"];

#[test]
fn golden_numbering_all_schemes() {
    let fx = load_fixture();
    // per-scheme (ok_matches, mismatches), with a few sample diffs
    let mut ok: BTreeMap<&str, u64> = SCHEMES.iter().map(|s| (*s, 0)).collect();
    let mut bad: BTreeMap<&str, u64> = SCHEMES.iter().map(|s| (*s, 0)).collect();
    let mut samples: Vec<String> = Vec::new();

    for seq in &fx.sequences {
        let bytes = seq.seq.as_bytes();
        for dom in &seq.domains {
            let sv = to_sv(&dom.state_vector);
            let ct = dom.chain_type.as_str();
            for scheme in SCHEMES {
                let expect = match dom.numbering.get(scheme) {
                    Some(n) => n,
                    None => continue,
                };
                let got = number_sequence_from_alignment(&sv, bytes, scheme, Some(ct));
                let mut mismatch: Option<String> = None;

                match (&expect.error, &got) {
                    (Some(err), Err(e)) => {
                        // fixture stores "AssertionError: <msg>"; ours is "<msg>".
                        let want = err.strip_prefix("AssertionError: ").unwrap_or(err);
                        if want != e.to_string() {
                            mismatch = Some(format!("err mismatch want {want:?} got {e}"));
                        }
                    }
                    (Some(err), Ok(_)) => {
                        mismatch = Some(format!("expected error {err:?} got Ok"));
                    }
                    (None, Err(e)) => {
                        mismatch = Some(format!("expected Ok got err {e}"));
                    }
                    (None, Ok(num)) => {
                        let exp_n = expect.numbering.as_ref().unwrap();
                        if num.residues.len() != exp_n.len() {
                            mismatch = Some(format!(
                                "len {} vs {}",
                                num.residues.len(),
                                exp_n.len()
                            ));
                        } else {
                            for (i, (r, ((pos, ins), aa))) in
                                num.residues.iter().zip(exp_n.iter()).enumerate()
                            {
                                if r.num != *pos
                                    || r.ins != ins.as_str()
                                    || (r.aa as char).to_string() != *aa
                                {
                                    mismatch = Some(format!(
                                        "pos {i}: got ({},{:?},{}) want ({},{:?},{})",
                                        r.num, r.ins, r.aa as char, pos, ins, aa
                                    ));
                                    break;
                                }
                            }
                        }
                        if mismatch.is_none() {
                            let gs = num.start.map(|v| v as i64);
                            let ge = num.end.map(|v| v as i64);
                            if gs != expect.start || ge != expect.end {
                                mismatch = Some(format!(
                                    "start/end got ({gs:?},{ge:?}) want ({:?},{:?})",
                                    expect.start, expect.end
                                ));
                            }
                        }
                    }
                }

                if let Some(m) = mismatch {
                    *bad.get_mut(scheme).unwrap() += 1;
                    if samples.len() < 25 {
                        samples.push(format!("[{scheme}] {} ct={ct}: {m}", seq.id));
                    }
                } else {
                    *ok.get_mut(scheme).unwrap() += 1;
                }
            }
        }
    }

    eprintln!("\n=== golden numbering per scheme ===");
    let mut total_bad = 0u64;
    for s in SCHEMES {
        eprintln!("  {:8} ok={:5} bad={:5}", s, ok[s], bad[s]);
        total_bad += bad[s];
    }
    if total_bad > 0 {
        eprintln!("\n--- sample mismatches ---");
        for s in &samples {
            eprintln!("  {s}");
        }
    }
    assert_eq!(total_bad, 0, "{total_bad} numbering mismatches vs reference");
}

#[test]
fn golden_germline_assignment() {
    let fx = load_fixture();
    let species: Vec<String> = ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut ok = 0u64;
    let mut bad = 0u64;
    let mut samples: Vec<String> = Vec::new();

    for seq in &fx.sequences {
        let bytes = seq.seq.as_bytes();
        for dom in &seq.domains {
            let exp = match &dom.germlines {
                Some(v) => v,
                None => continue,
            };
            let sv = to_sv(&dom.state_vector);
            let g = run_germline_assignment(&sv, bytes, &dom.chain_type, Some(&species));

            // Reference dict: {} (empty) or {'v_gene':[[sp,gene],id]|[null,null], 'j_gene':...}
            let exp_obj = exp.as_object().unwrap();
            let mut diff: Option<String> = None;
            if exp_obj.is_empty() {
                if !g.empty {
                    diff = Some("expected {} empty".into());
                }
            } else {
                // v_gene
                let (vg, vid) = parse_gene(exp_obj.get("v_gene"));
                let (jg, jid) = parse_gene(exp_obj.get("j_gene"));
                if g.v_gene != vg || !approx(g.v_identity, vid) {
                    diff = Some(format!("v got {:?},{:?} want {:?},{:?}", g.v_gene, g.v_identity, vg, vid));
                } else if g.j_gene != jg || !approx(g.j_identity, jid) {
                    diff = Some(format!("j got {:?},{:?} want {:?},{:?}", g.j_gene, g.j_identity, jg, jid));
                }
            }
            if let Some(d) = diff {
                bad += 1;
                if samples.len() < 25 {
                    samples.push(format!("{} ct={}: {d}", seq.id, dom.chain_type));
                }
            } else {
                ok += 1;
            }
        }
    }
    eprintln!("\n=== germline assignment ok={ok} bad={bad} ===");
    for s in &samples {
        eprintln!("  {s}");
    }
    assert_eq!(bad, 0, "{bad} germline mismatches vs reference");
}

fn parse_gene(v: Option<&serde_json::Value>) -> (Option<(String, String)>, Option<f64>) {
    // [[species, gene], identity] or [null, null]
    let arr = match v.and_then(|x| x.as_array()) {
        Some(a) => a,
        None => return (None, None),
    };
    let gene = arr.first().and_then(|g| g.as_array()).map(|g| {
        (
            g[0].as_str().unwrap().to_string(),
            g[1].as_str().unwrap().to_string(),
        )
    });
    let ident = arr.get(1).and_then(|x| x.as_f64());
    (gene, ident)
}

fn approx(a: Option<f64>, b: Option<f64>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => (x - y).abs() < 1e-9,
        _ => false,
    }
}
