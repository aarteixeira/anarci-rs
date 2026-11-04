//! Bounded synthetic codec/decoder properties. No fixtures, model, or Candle required.
use anarci_core::{State, StateType};
use anarci_learned::decode::{index_to_label, label_to_index, viterbi_decode, NUM_LABELS};
use anarci_learned::{
    decode_domain, decode_sequence, encode_domain, strip_boundary_deletes, Label, LabeledRes,
};

fn state(id: u8, typ: StateType, si: Option<usize>) -> State {
    State { id, typ, si }
}

fn confident(labels: &[Label]) -> Vec<Vec<f32>> {
    labels
        .iter()
        .map(|&label| {
            let mut row = vec![-100.0; NUM_LABELS];
            row[label_to_index(label)] = 100.0;
            row
        })
        .collect()
}

#[test]
fn exhaustive_vocabulary_bijection_and_domain_partition() {
    let mut domain_count = 0;
    for index in 0..NUM_LABELS {
        let label = index_to_label(index);
        assert_eq!(label_to_index(label), index);
        assert_eq!(label.is_domain(), index < 256);
        domain_count += usize::from(label.is_domain());
    }
    assert_eq!(domain_count, 256);
}

#[test]
fn exhaustive_small_canonical_domains_roundtrip_with_gaps_and_owned_insertions() {
    // Eight independent internal columns, each M, D, or I+M: 3^8 shapes.
    for code in 0..3usize.pow(8) {
        let mut choices = code;
        let mut si = 7; // nonzero query origin must survive the codec
        let mut sv = vec![state(1, StateType::M, Some(si))];
        si += 1;
        for id in 2..=9 {
            match choices % 3 {
                0 => {
                    sv.push(state(id, StateType::M, Some(si)));
                    si += 1;
                }
                1 => sv.push(state(id, StateType::D, None)),
                _ => {
                    for _ in 0..2 {
                        sv.push(state(id, StateType::I, Some(si)));
                        si += 1;
                    }
                    sv.push(state(id, StateType::M, Some(si)));
                    si += 1;
                }
            }
            choices /= 3;
        }
        sv.push(state(10, StateType::M, Some(si)));
        let labels = encode_domain(&sv);
        assert!(labels.windows(2).all(|w| w[0].si + 1 == w[1].si));
        assert_eq!(labels.len(), sv.iter().filter(|s| s.si.is_some()).count());
        assert_eq!(decode_domain(&labels), sv, "shape {code}");
        let mut padded = vec![state(1, StateType::D, None)];
        padded.extend_from_slice(&sv);
        padded.push(state(11, StateType::D, None));
        assert_eq!(strip_boundary_deletes(&padded), sv);
        assert_eq!(strip_boundary_deletes(&strip_boundary_deletes(&padded)), sv);
    }
}

#[test]
fn every_separator_and_abutting_reset_preserve_absolute_indices() {
    for separator in [
        None,
        Some(Label::Out),
        Some(Label::Linker),
        Some(Label::Constant),
        Some(Label::Hinge),
        Some(Label::Other),
    ] {
        for first in [1, 23, 104, 127] {
            let mut labels = vec![
                Label::Match(first),
                Label::Ins(first + 1),
                Label::Match(first + 1),
            ];
            labels.extend(separator);
            let offset = labels.len();
            labels.extend([Label::Match(1), Label::Match(2)]);
            let domains = decode_sequence(&labels);
            assert_eq!(domains.len(), 2);
            assert_eq!(
                encode_domain(&domains[0]),
                vec![
                    LabeledRes {
                        si: 0,
                        label: Label::Match(first)
                    },
                    LabeledRes {
                        si: 1,
                        label: Label::Ins(first + 1)
                    },
                    LabeledRes {
                        si: 2,
                        label: Label::Match(first + 1)
                    },
                ]
            );
            assert_eq!(
                encode_domain(&domains[1]),
                vec![
                    LabeledRes {
                        si: offset,
                        label: Label::Match(1)
                    },
                    LabeledRes {
                        si: offset + 1,
                        label: Label::Match(2)
                    },
                ]
            );
        }
    }
}

#[test]
fn malformed_raw_codec_inputs_are_bounded_and_do_not_overflow() {
    // The raw Rust codec is permissive, unlike the validating Python boundary.
    // Do not call label_to_index on out-of-vocabulary labels (its input precondition).
    for id in 0..=u8::MAX {
        let malformed = [
            LabeledRes {
                si: usize::MAX,
                label: Label::Match(id),
            },
            LabeledRes {
                si: 0,
                label: Label::Match(id),
            }, // duplicate: last match wins
            LabeledRes {
                si: 1,
                label: Label::Ins(id),
            },
            LabeledRes {
                si: 2,
                label: Label::Out,
            },
        ];
        let decoded = decode_domain(&malformed);
        assert_eq!(
            decoded,
            vec![
                state(id, StateType::I, Some(1)),
                state(id, StateType::M, Some(0))
            ]
        );
        assert!(decode_domain(&[LabeledRes {
            si: 0,
            label: Label::Ins(id)
        }])
        .is_empty());
    }
    assert!(decode_domain(&[]).is_empty());
    assert!(strip_boundary_deletes(&[state(1, StateType::D, None)]).is_empty());
}

#[test]
fn confident_ascending_decoder_paths_roundtrip_exactly() {
    for start in [1, 23, 104, 120] {
        let mut labels = vec![Label::Out];
        for id in start..start + 7 {
            labels.push(Label::Match(id));
            labels.push(Label::Ins(id + 1));
        }
        labels.push(Label::Match(start + 7));
        labels.push(Label::Constant);
        assert_eq!(viterbi_decode(&confident(&labels)), labels);
        let domains = decode_sequence(&labels);
        assert_eq!(domains.len(), 1);
        let recovered = encode_domain(&domains[0]);
        assert_eq!(recovered.len(), labels.len() - 2);
        assert!(recovered.windows(2).all(|w| w[0].si + 1 == w[1].si));
    }
}

#[test]
fn seeded_decoder_outputs_preserve_every_domain_residue_in_query_order() {
    // Fixed LCG; short score matrices bound runtime and make failures reproducible.
    let mut seed = 0x5eed_u64;
    for case in 0..24 {
        let length = 1 + case % 12;
        let logits: Vec<Vec<f32>> = (0..length)
            .map(|_| {
                (0..NUM_LABELS)
                    .map(|_| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) % 101) as f32 - 50.0
                    })
                    .collect()
            })
            .collect();
        let labels = viterbi_decode(&logits);
        assert_eq!(
            labels,
            viterbi_decode(&logits),
            "case {case} is nondeterministic"
        );
        assert_eq!(labels.len(), length);
        let expected: Vec<_> = labels
            .iter()
            .enumerate()
            .filter(|(_, l)| l.is_domain())
            .map(|(si, &label)| LabeledRes { si, label })
            .collect();
        let recovered: Vec<_> = decode_sequence(&labels)
            .iter()
            .flat_map(|sv| encode_domain(sv))
            .collect();
        assert_eq!(
            recovered, expected,
            "case {case}: decoder reordered/dropped query residues"
        );
    }
}

#[test]
fn decoder_cannot_reorder_residues_when_insertion_owner_scores_decrease() {
    // Both inserts are legal relative to Match(1) in the current transition graph,
    // but grouping by owner during codec decode would reverse their query order.
    let preferred = [
        Label::Match(1),
        Label::Ins(3),
        Label::Ins(2),
        Label::Match(4),
    ];
    let labels = viterbi_decode(&confident(&preferred));
    let indices: Vec<_> = decode_sequence(&labels)
        .iter()
        .flat_map(|sv| sv.iter())
        .filter_map(|s| s.si)
        .collect();
    assert!(
        indices.windows(2).all(|w| w[0] < w[1]),
        "decoder reordered query residues: {labels:?} -> {indices:?}"
    );
    assert_eq!(
        indices.len(),
        labels.iter().filter(|l| l.is_domain()).count()
    );
}
