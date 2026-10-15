//! P4 gate: the dataset split is leakage-controlled. Two levels:
//!  * bookkeeping — no source cluster appears in more than one split; multi-V examples pair
//!    domains from a single split;
//!  * the real invariant — no near-duplicate (full-sequence or CDR3, length-tolerant) source
//!    domain crosses train↔test. This is the assertion an earlier per-domain / same-length
//!    version FAILED (4 identical full chains, 17 near-dup CDR3s across the split).

use anarci_core::{State, StateType};
use anarci_learned::data::{
    build_dataset, near_duplicate, splits_of, DomainRec, Split,
};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Golden {
    sequences: Vec<SeqRec>,
}
#[derive(Deserialize)]
struct SeqRec {
    seq: String,
    domains: Vec<DomRec>,
}
#[derive(Deserialize)]
struct DomRec {
    state_vector: Vec<((u8, String), Option<usize>)>,
}
fn to_sv(raw: &[((u8, String), Option<usize>)]) -> Vec<State> {
    raw.iter()
        .map(|((id, t), si)| State {
            id: *id,
            typ: match t.as_str() {
                "m" => StateType::M,
                "i" => StateType::I,
                _ => StateType::D,
            },
            si: *si,
        })
        .collect()
}
fn domains() -> Vec<DomainRec> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden.json.gz");
    let g: Golden = serde_json::from_reader(GzDecoder::new(std::fs::File::open(path).unwrap())).unwrap();
    let mut out = Vec::new();
    for s in &g.sequences {
        for d in &s.domains {
            out.push(DomainRec { seq: s.seq.as_bytes().to_vec(), sv: to_sv(&d.state_vector) });
        }
    }
    out
}

/// Bookkeeping: no cluster straddles splits; multi-V examples are single-split.
#[test]
fn no_cluster_leaks_across_splits() {
    let doms = domains();
    let examples = build_dataset(&doms, &[(66, 128), (39, 128)]);
    assert!(examples.len() > doms.len(), "dataset should include augmented examples");

    let mut cluster_split: HashMap<usize, Split> = HashMap::new();
    for ex in &examples {
        for &c in &ex.clusters {
            match cluster_split.get(&c) {
                Some(&prev) => assert_eq!(prev, ex.split, "cluster {c} in two splits"),
                None => {
                    cluster_split.insert(c, ex.split);
                }
            }
        }
    }
    for ex in examples.iter().filter(|e| e.kind == "multi_v") {
        assert_eq!(ex.clusters.len(), 2);
        assert_eq!(cluster_split[&ex.clusters[0]], cluster_split[&ex.clusters[1]], "multi_v crosses splits");
    }
    for sp in [Split::Train, Split::Val, Split::Test] {
        assert!(examples.iter().any(|e| e.split == sp), "split {sp:?} is empty");
        // junk must reach every split (E2 non-antibody rejection measurable on held-out)
        assert!(examples.iter().any(|e| e.split == sp && e.kind == "junk"), "no junk in {sp:?}");
    }
}

/// The real invariant: NO near-duplicate source domain crosses train↔test (full sequence or
/// CDR3, length-tolerant). Fails on the old per-domain / same-length clustering.
#[test]
fn no_near_duplicate_across_train_test() {
    let doms = domains();
    let splits = splits_of(&doms);
    let train: Vec<usize> = (0..doms.len()).filter(|&i| splits[i] == Split::Train).collect();
    let test: Vec<usize> = (0..doms.len()).filter(|&i| splits[i] == Split::Test).collect();
    assert!(!train.is_empty() && !test.is_empty());

    let mut leaks = 0;
    for &t in &test {
        for &r in &train {
            if near_duplicate(&doms[t], &doms[r]) {
                leaks += 1;
                break;
            }
        }
    }
    assert_eq!(leaks, 0, "{leaks} test domains have a near-duplicate in train (leakage)");
}

/// scFv/multi-domain records go to ONE split (both domains of a shared source sequence).
#[test]
fn multi_domain_records_stay_in_one_split() {
    let doms = domains();
    let splits = splits_of(&doms);
    // group domain indices by identical source sequence (= record)
    let mut by_seq: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (i, d) in doms.iter().enumerate() {
        by_seq.entry(d.seq.as_slice()).or_default().push(i);
    }
    let mut multi = 0;
    for idxs in by_seq.values().filter(|v| v.len() > 1) {
        multi += 1;
        let s0 = splits[idxs[0]];
        for &i in idxs {
            assert_eq!(splits[i], s0, "a multi-domain record's domains landed in different splits");
        }
    }
    assert!(multi > 0, "fixture should contain multi-domain records");
}
