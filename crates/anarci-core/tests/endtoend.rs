//! End-to-end orchestration gate (engine-independent): drive the FULL pipeline
//! (parse → check_for_j → number → germline → assemble) with a ReplayEngine that
//! serves the exact HSPs reference ANARCI saw, and compare anarci() output
//! (numbered, details incl. germlines, hit_table) to reference for all 996 seqs.

use anarci_core::{default_allow, run_anarci, Hsp, HmmEngine};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Deserialize, Clone)]
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

struct ReplayEngine {
    map: HashMap<String, Vec<HspJson>>,
}
impl HmmEngine for ReplayEngine {
    fn scan_one(&self, _name: &str, seq: &[u8]) -> Vec<Hsp> {
        let key = std::str::from_utf8(seq).unwrap();
        match self.map.get(key) {
            Some(hsps) => hsps
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
                .collect(),
            None => Vec::new(),
        }
    }
}

#[derive(Deserialize)]
struct RefFixture {
    sequences: Vec<RefSeq>,
}
#[derive(Deserialize)]
struct RefSeq {
    id: String,
    seq: String,
    numbered: Option<Vec<(Vec<((i32, String), String)>, Option<i64>, Option<i64>)>>,
    details: Option<Vec<RefDetail>>,
    hit_table: Vec<Vec<serde_json::Value>>,
}
#[derive(Deserialize)]
struct RefDetail {
    species: Option<String>,
    chain_type: Option<String>,
    evalue: Option<f64>,
    bitscore: Option<f64>,
    bias: Option<f64>,
    query_start: Option<i64>,
    query_end: Option<i64>,
    germlines: Option<serde_json::Value>,
}

fn load_gz<T: serde::de::DeserializeOwned>(rel: &str) -> T {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures").join(rel);
    let f = std::fs::File::open(&path).unwrap_or_else(|e| panic!("open {:?}: {}", path, e));
    serde_json::from_reader(GzDecoder::new(f)).expect("parse fixture")
}

fn species() -> Vec<String> {
    ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"]
        .iter()
        .map(|s| s.to_string())
        .collect()
}

fn parse_gene(v: Option<&serde_json::Value>) -> (Option<(String, String)>, Option<f64>) {
    let arr = match v.and_then(|x| x.as_array()) {
        Some(a) => a,
        None => return (None, None),
    };
    let gene = arr.first().and_then(|g| g.as_array()).map(|g| {
        (g[0].as_str().unwrap().to_string(), g[1].as_str().unwrap().to_string())
    });
    (gene, arr.get(1).and_then(|x| x.as_f64()))
}

#[test]
fn endtoend_imgt_parity() {
    let map: HashMap<String, Vec<HspJson>> = load_gz("replay_hsps.json.gz");
    let reff: RefFixture = load_gz("endtoend_imgt.json.gz");
    let engine = ReplayEngine { map };
    let sp = species();

    let seqs: Vec<(String, Vec<u8>)> =
        reff.sequences.iter().map(|s| (s.id.clone(), s.seq.clone().into_bytes())).collect();

    let results = run_anarci(&engine, &seqs, "imgt", &default_allow(), true, Some(&sp), 80.0)
        .expect("run_anarci");

    let mut ok = 0u64;
    let mut bad = 0u64;
    let mut samples: Vec<String> = Vec::new();

    for (r, exp) in results.iter().zip(reff.sequences.iter()) {
        let mut diff: Option<String> = None;

        // numbered
        match (&r.numbered, &exp.numbered) {
            (None, None) => {}
            (Some(_), None) => diff = Some("got Some numbered, want None".into()),
            (None, Some(_)) => diff = Some("got None numbered, want Some".into()),
            (Some(got), Some(want)) => {
                if got.len() != want.len() {
                    diff = Some(format!("ndomains {} vs {}", got.len(), want.len()));
                } else {
                    for (gd, (wnum, wstart, wend)) in got.iter().zip(want.iter()) {
                        if gd.residues.len() != wnum.len() {
                            diff = Some(format!("dom len {} vs {}", gd.residues.len(), wnum.len()));
                            break;
                        }
                        for (res, ((p, ins), aa)) in gd.residues.iter().zip(wnum.iter()) {
                            if res.num != *p
                                || res.ins != ins.as_str()
                                || (res.aa as char).to_string() != *aa
                            {
                                diff = Some(format!(
                                    "res got ({},{:?},{}) want ({},{:?},{})",
                                    res.num, res.ins, res.aa as char, p, ins, aa
                                ));
                                break;
                            }
                        }
                        if diff.is_none()
                            && (gd.start.map(|v| v as i64) != *wstart
                                || gd.end.map(|v| v as i64) != *wend)
                        {
                            diff = Some("start/end mismatch".into());
                        }
                        if diff.is_some() {
                            break;
                        }
                    }
                }
            }
        }

        // details (incl germlines)
        if diff.is_none() {
            match (&r.details, &exp.details) {
                (None, None) => {}
                (Some(gd), Some(wd)) if gd.len() == wd.len() => {
                    for (g, w) in gd.iter().zip(wd.iter()) {
                        if Some(g.species.as_str()) != w.species.as_deref()
                            || Some(g.chain_type.as_str()) != w.chain_type.as_deref()
                            || (g.bitscore - w.bitscore.unwrap()).abs() > 1e-9
                            || g.query_start.map(|v| v as i64) != w.query_start
                            || Some(g.query_end as i64) != w.query_end
                        {
                            diff = Some(format!(
                                "details got (sp={},ct={},bs={},qs={:?},qe={}) want (sp={:?},ct={:?},bs={:?},qs={:?},qe={:?})",
                                g.species, g.chain_type, g.bitscore, g.query_start, g.query_end,
                                w.species, w.chain_type, w.bitscore, w.query_start, w.query_end
                            ));
                            break;
                        }
                        // germlines
                        let gobj = w.germlines.as_ref().and_then(|v| v.as_object());
                        let germ = g.germlines.as_ref().unwrap();
                        match gobj {
                            Some(o) if o.is_empty() => {
                                if !germ.empty {
                                    diff = Some("germline want empty".into());
                                    break;
                                }
                            }
                            Some(o) => {
                                let (vg, vid) = parse_gene(o.get("v_gene"));
                                let (jg, jid) = parse_gene(o.get("j_gene"));
                                let close = |a: Option<f64>, b: Option<f64>| match (a, b) {
                                    (None, None) => true,
                                    (Some(x), Some(y)) => (x - y).abs() < 1e-9,
                                    _ => false,
                                };
                                if germ.v_gene != vg
                                    || !close(germ.v_identity, vid)
                                    || germ.j_gene != jg
                                    || !close(germ.j_identity, jid)
                                {
                                    diff = Some("germline mismatch".into());
                                    break;
                                }
                            }
                            None => {}
                        }
                    }
                }
                _ => diff = Some("details presence/len mismatch".into()),
            }
        }

        // hit_table (skip reference header row 0)
        if diff.is_none() {
            let want_rows = &exp.hit_table[1..];
            if r.hit_table.len() != want_rows.len() {
                diff = Some(format!("hit_table len {} vs {}", r.hit_table.len(), want_rows.len()));
            } else {
                for (g, w) in r.hit_table.iter().zip(want_rows.iter()) {
                    let id = w[0].as_str().unwrap();
                    let bs = w[3].as_f64().unwrap();
                    let qs = w[5].as_i64().unwrap();
                    let qe = w[6].as_i64().unwrap();
                    if g.id != id
                        || (g.bitscore - bs).abs() > 1e-9
                        || g.query_start as i64 != qs
                        || g.query_end as i64 != qe
                    {
                        diff = Some("hit_table row mismatch".into());
                        break;
                    }
                }
            }
        }

        if let Some(m) = diff {
            bad += 1;
            if samples.len() < 25 {
                samples.push(format!("{}: {m}", exp.id));
            }
        } else {
            ok += 1;
        }
    }

    eprintln!("\n=== end-to-end imgt parity: ok={ok} bad={bad} ===");
    for s in &samples {
        eprintln!("  {s}");
    }
    assert_eq!(bad, 0, "{bad} end-to-end mismatches vs reference anarci()");
}
