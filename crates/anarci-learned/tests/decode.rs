//! P3 constrained-decode gate: the decoder's output is always valid (C2), recovers the
//! true labeling from confident logits, and resists noise-driven over-segmentation.

use anarci_core::{number_sequence_from_alignment, State, StateType};
use anarci_learned::decode::{index_to_label, label_to_index, viterbi_decode, NUM_LABELS};
use anarci_learned::{decode_sequence, encode_domain, Label};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
struct Golden {
    sequences: Vec<SeqRec>,
}
#[derive(Deserialize)]
struct SeqRec {
    domains: Vec<DomRec>,
}
#[derive(Deserialize)]
struct DomRec {
    state_vector: Vec<((u8, String), Option<usize>)>,
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

/// One-hot logits (score `hi` on `idx`, 0 elsewhere).
fn onehot(idxs: &[usize], hi: f32) -> Vec<Vec<f32>> {
    idxs.iter()
        .map(|&idx| {
            let mut row = vec![0.0f32; NUM_LABELS];
            row[idx] = hi;
            row
        })
        .collect()
}

/// A decoded state vector is valid: strictly-increasing match ids, no duplicate positions.
fn is_valid_sv(sv: &[State]) -> bool {
    let matches: Vec<u8> = sv.iter().filter(|s| s.typ == StateType::M).map(|s| s.id).collect();
    matches.windows(2).all(|w| w[0] < w[1])
}

#[test]
fn label_index_roundtrip() {
    for i in 0..NUM_LABELS {
        assert_eq!(label_to_index(index_to_label(i)), i, "index {i} did not round-trip");
    }
}

/// Confident (one-hot) logits on a real domain's true labels recover those labels, and the
/// decode reproduces the domain structure.
#[test]
fn viterbi_recovers_true_labels() {
    let g = load();
    let mut recovered = 0;
    for s in g.sequences.iter().take(50) {
        let sv = to_sv(&s.domains[0].state_vector);
        let true_labels: Vec<Label> = encode_domain(&sv).into_iter().map(|r| r.label).collect();
        let idxs: Vec<usize> = true_labels.iter().map(|&l| label_to_index(l)).collect();
        let got = viterbi_decode(&onehot(&idxs, 10.0));
        assert_eq!(got, true_labels, "viterbi failed to recover true labels for a real domain");
        recovered += 1;
    }
    assert_eq!(recovered, 50);
}

/// The core C2 guarantee: for ANY logits — random and adversarial (all favoring one match
/// position, which argmax would turn into dropped residues) — the decoded output is always
/// valid, every residue is accounted for (none dropped), and every domain numbers.
#[test]
fn viterbi_output_always_valid() {
    // adversarial: every residue's top choice is Match(50) (would collide under naive argmax)
    let dup = vec![label_to_index(Label::Match(50)); 40];
    // pseudo-random but deterministic logits (no rng in tests): vary top label by index
    let rand_like: Vec<Vec<f32>> = (0..60)
        .map(|i| {
            let mut row = vec![0.0f32; NUM_LABELS];
            row[(i * 37 + 11) % NUM_LABELS] = 5.0;
            row[(i * 13 + 3) % NUM_LABELS] = 3.0;
            row
        })
        .collect();

    for logits in [onehot(&dup, 8.0), rand_like] {
        let labels = viterbi_decode(&logits);
        assert_eq!(labels.len(), logits.len());
        let domains = decode_sequence(&labels);
        // validity: every domain strictly monotonic
        for sv in &domains {
            assert!(is_valid_sv(sv), "decoded domain is not strictly monotonic: {sv:?}");
            // and it numbers through the real pipeline without error
            let seq = vec![b'A'; labels.len()];
            number_sequence_from_alignment(sv, &seq, "imgt", Some("H")).expect("must number");
        }
        // no dropped residues: #residues placed in domains == #domain labels emitted
        let placed: usize = domains.iter().flat_map(|d| d.iter()).filter(|s| s.si.is_some()).count();
        let domain_labels = labels.iter().filter(|l| l.is_domain()).count();
        assert_eq!(placed, domain_labels, "residues were dropped during decode");
    }
}

/// Reset penalty resists over-segmentation: a single confident domain interrupted by one
/// weak decreasing blip stays ONE domain (the blip can't justify a spurious reset).
#[test]
fn viterbi_resists_oversegmentation() {
    // Strong ascending matches 1..=20, but residue 10 has a weak competing Match(5).
    let mut logits: Vec<Vec<f32>> = (1..=20u8)
        .map(|p| {
            let mut row = vec![0.0f32; NUM_LABELS];
            row[label_to_index(Label::Match(p))] = 10.0;
            row
        })
        .collect();
    // add a small decoy that would reset (5 < 10) if not penalized
    logits[9][label_to_index(Label::Match(5))] = 10.5; // slightly higher than the true 10.0

    let labels = viterbi_decode(&logits);
    let domains = decode_sequence(&labels);
    // Without the reset penalty the 10.5 decoy would win and fragment into 2 domains;
    // with it, the +0.5 gain can't cover the -6.0 reset, so it stays one domain.
    assert_eq!(domains.len(), 1, "reset penalty should prevent over-segmentation");
}

/// Non-finite logits (NaN / +Inf) must not silently drop residues to `Out`. A clean
/// ascending run with an injected +Inf and NaN must still place all residues in the domain.
#[test]
fn viterbi_robust_to_nonfinite_logits() {
    let mut logits: Vec<Vec<f32>> = (1..=10u8)
        .map(|p| {
            let mut row = vec![0.0f32; NUM_LABELS];
            row[label_to_index(Label::Match(p))] = 5.0;
            row
        })
        .collect();
    logits[5][label_to_index(Label::Match(6))] = f32::INFINITY; // true label, saturated (+Inf)
    logits[6][label_to_index(Label::Match(99))] = f32::NAN; // a NaN decoy on a non-true label

    let labels = viterbi_decode(&logits);
    // the +Inf keeps its label (definitely Match(6)); the NaN decoy is ignored; nothing is
    // dropped by a broken DP — the clean ascending numbering is fully recovered.
    let expected: Vec<Label> = (1..=10u8).map(Label::Match).collect();
    assert_eq!(labels, expected, "non-finite logits corrupted/dropped the decode");
    let placed: usize = decode_sequence(&labels)
        .iter()
        .flat_map(|d| d.iter())
        .filter(|s| s.si.is_some())
        .count();
    assert_eq!(placed, 10, "all 10 residues must be placed despite the non-finite logits");
}

#[test]
fn viterbi_empty() {
    assert!(viterbi_decode(&[]).is_empty());
}
