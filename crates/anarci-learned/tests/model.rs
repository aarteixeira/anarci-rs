//! P3 gate (feature = "model"): the model learns (overfits a tiny set), and its logits —
//! run through the constrained decoder — are always valid and numberable.
//! Run with: `cargo test -p anarci-learned --features model --release`.
#![cfg(feature = "model")]

use anarci_core::{number_sequence_from_alignment, State, StateType};
use anarci_learned::decode::label_to_index;
use anarci_learned::model::{logits_rows, overfit, tokenize};
use anarci_learned::{decode_sequence, encode_domain, viterbi_decode};
use candle_core::Device;
use flate2::read::GzDecoder;
use serde::Deserialize;
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
    chain_type: String,
}
fn load() -> Golden {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/golden.json.gz");
    serde_json::from_reader(GzDecoder::new(std::fs::File::open(path).unwrap())).unwrap()
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

#[test]
fn overfit_tiny_set_and_valid_decode() {
    let g = load();
    let mut examples = Vec::new();
    let mut probes = Vec::new(); // (residues, chain_type, oracle_sv) for the correctness check
    for s in g.sequences.iter().take(8) {
        let sv = to_sv(&s.domains[0].state_vector);
        let enc = encode_domain(&sv);
        let residues: Vec<u8> = enc.iter().map(|r| s.seq.as_bytes()[r.si]).collect();
        let targets: Vec<u32> = enc.iter().map(|r| label_to_index(r.label) as u32).collect();
        examples.push((tokenize(&residues), targets));
        probes.push((residues, s.domains[0].chain_type.clone(), sv));
    }

    let dev = Device::Cpu;
    // small model; 800 steps -> full memorization (loss ~0.016) so numbering matches the oracle.
    let (model, acc) = overfit(&examples, 128, 2, 256, 800, 3e-3, &dev).unwrap();
    assert!(acc > 0.95, "P3 learn-ability gate failed: overfit accuracy = {acc}");

    // The decoded output is (a) always valid + numberable AND (b) CORRECT vs the oracle
    // (not just "numbers OK" — a valid-but-wrong numbering must fail the gate).
    let mut correct_domains = 0usize;
    for (residues, ct, oracle_sv) in &probes {
        let rows = logits_rows(&model, residues, &dev).unwrap();
        let labels = viterbi_decode(&rows);
        assert_eq!(labels.len(), residues.len());
        let decoded = decode_sequence(&labels);
        for sv in &decoded {
            let ms: Vec<u8> = sv.iter().filter(|s| s.typ == StateType::M).map(|s| s.id).collect();
            assert!(ms.windows(2).all(|w| w[0] < w[1]), "invalid decoded numbering");
            let numbered =
                number_sequence_from_alignment(sv, residues, "imgt", Some(ct)).expect("must number");
            assert!(numbered.residues.iter().any(|r| r.aa != b'-'), "numbered to zero residues");
        }
        // correctness: the recovered domain equals the oracle (id, type) sequence.
        if decoded.len() == 1
            && decoded[0].iter().map(|s| (s.id, s.typ)).eq(oracle_sv.iter().map(|s| (s.id, s.typ)))
        {
            correct_domains += 1;
        }
    }
    assert!(
        correct_domains >= 7,
        "P3 correctness gate: only {correct_domains}/8 decoded numberings match the oracle"
    );
}
