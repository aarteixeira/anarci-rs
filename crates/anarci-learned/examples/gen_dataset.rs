//! P4 data pipeline (runnable): build a split, leakage-controlled training dataset from the
//! oracle golden domains + augmentation, and report its composition.
//!
//!   cargo run -p anarci-learned --example gen_dataset --release
//!
//! Prints per-split / per-kind example counts. The Test split is the held-out eval set for
//! *standard* numbering and *synthetic* capability suites (partials/multi-V) — derived only
//! from test-split domains, so no leakage. NOTE: genuine external eval data (real 5'-RACE
//! partial reads, real scFv/BiTE constructs, constant-region-bearing chains) is NOT bundled
//! and must be curated before the corresponding paper experiments (E5/E4/E6) are claimed —
//! see LEARNED_ENGINE_PLAN.md §5.

use anarci_core::{State, StateType};
use anarci_learned::data::{build_dataset, DomainRec, Split};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::BTreeMap;
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

fn main() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden.json.gz");
    let g: Golden = serde_json::from_reader(GzDecoder::new(std::fs::File::open(path).unwrap())).unwrap();
    let mut doms = Vec::new();
    for s in &g.sequences {
        for d in &s.domains {
            let sv = d
                .state_vector
                .iter()
                .map(|((id, t), si)| State {
                    id: *id,
                    typ: match t.as_str() {
                        "m" => StateType::M,
                        "i" => StateType::I,
                        _ => StateType::D,
                    },
                    si: *si,
                })
                .collect();
            doms.push(DomainRec { seq: s.seq.as_bytes().to_vec(), sv });
        }
    }

    let examples = build_dataset(&doms, &[(66, 128), (39, 128), (80, 128)]);

    // counts[(split, kind)] and residue totals
    let mut counts: BTreeMap<(String, &str), usize> = BTreeMap::new();
    let mut residues: BTreeMap<String, usize> = BTreeMap::new();
    for ex in &examples {
        let split = format!("{:?}", ex.split);
        *counts.entry((split.clone(), ex.kind)).or_default() += 1;
        *residues.entry(split).or_default() += ex.data.seq.len();
    }

    println!("source domains: {}", doms.len());
    println!("total examples: {}\n", examples.len());
    println!("{:<7} {:<10} {:>8}", "split", "kind", "count");
    for ((split, kind), n) in &counts {
        println!("{split:<7} {kind:<10} {n:>8}");
    }
    println!("\nresidues (labeled targets) per split:");
    for (split, n) in &residues {
        println!("  {split:<7} {n:>10}");
    }
    for sp in [Split::Train, Split::Val, Split::Test] {
        let n = examples.iter().filter(|e| e.split == sp).count();
        println!("split {:?}: {} examples", sp, n);
    }
}
