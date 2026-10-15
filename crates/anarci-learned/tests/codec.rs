//! C0 gate: the label codec is a lossless round-trip of real ANARCI state vectors,
//! and its decoded output, fed through the existing byte-exact `anarci-core` numbering,
//! reproduces the reference-ANARCI numbering.
//!
//! Oracle: `tests/fixtures/golden.json.gz` (996 seqs / 1013 domains, captured from
//! reference ANARCI conda 2024.05.21).

use anarci_core::{number_sequence_from_alignment, State, StateType};
use anarci_learned::{decode_domain, decode_sequence, encode_domain, Label, LabeledRes};
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
    numbering: SchemeNumberings,
}
#[derive(Deserialize)]
struct SchemeNumberings {
    imgt: ImgtNum,
}
#[derive(Deserialize)]
struct ImgtNum {
    #[serde(rename = "numbering")]
    residues: Vec<((i32, String), String)>,
}

fn load() -> Golden {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/golden.json.gz");
    let f = std::fs::File::open(&path).expect("open golden.json.gz");
    serde_json::from_reader(GzDecoder::new(f)).expect("parse golden.json.gz")
}

fn to_sv(raw: &[((u8, String), Option<usize>)]) -> Vec<State> {
    raw.iter()
        .map(|((id, t), si)| {
            let typ = match t.as_str() {
                "m" => StateType::M,
                "i" => StateType::I,
                "d" => StateType::D,
                other => panic!("bad state type {other:?}"),
            };
            State { id: *id, typ, si: *si }
        })
        .collect()
}

/// The make-or-break gate: `decode(encode(sv)) == sv` byte-for-byte, all 1013 domains.
#[test]
fn c0_roundtrip_all_golden_domains() {
    let g = load();
    let mut total = 0usize;
    let mut fails = Vec::new();
    for s in &g.sequences {
        for (di, dom) in s.domains.iter().enumerate() {
            total += 1;
            let sv = to_sv(&dom.state_vector);
            let rt = decode_domain(&encode_domain(&sv));
            if rt != sv {
                fails.push(format!("seq[{}] dom{di}", &s.seq[..8.min(s.seq.len())]));
            }
        }
    }
    assert_eq!(
        fails.len(),
        0,
        "codec round-trip failed on {}/{} domains: {:?}",
        fails.len(),
        total,
        &fails[..fails.len().min(5)]
    );
    assert_eq!(total, 1013, "expected 1013 golden domains, saw {total}");
}

/// End-to-end: codec output → existing `number_sequence_from_alignment` → oracle numbering.
/// Since the round-trip is byte-exact this must hold, but it proves the decoded state
/// vector actually plugs into the real numbering pipeline and reproduces ANARCI (imgt).
#[test]
fn decoded_sv_numbers_like_oracle_imgt() {
    let g = load();
    let mut checked = 0usize;
    let mut mismatches = 0usize;
    for s in &g.sequences {
        let seq = s.seq.as_bytes();
        for dom in &s.domains {
            let sv = decode_domain(&encode_domain(&to_sv(&dom.state_vector)));
            let got = match number_sequence_from_alignment(&sv, seq, "imgt", Some(&dom.chain_type)) {
                Ok(n) => n,
                Err(_) => continue, // known-errored numbering; skip
            };
            let want = &dom.numbering.imgt.residues;
            if got.residues.len() != want.len() {
                mismatches += 1;
                continue;
            }
            let ok = got.residues.iter().zip(want).all(|(r, ((pos, ins), aa))| {
                r.num == *pos && r.ins == ins.as_str() && (r.aa as char) == aa.chars().next().unwrap_or('-')
            });
            checked += 1;
            if !ok {
                mismatches += 1;
            }
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} imgt numbering mismatches via codec (checked {checked})");
    // all 1013 golden domains have non-empty imgt numbering, so the codec output must
    // number every one of them through the real pipeline (no silent mass-skip).
    assert_eq!(checked, 1013, "expected to check all 1013 domains, checked {checked}");
}

/// Whole-sequence segmentation recovers a single domain exactly (sanity for the
/// multi-domain path used by scFv/tandem/partial inputs).
#[test]
fn segmentation_recovers_single_domain() {
    let g = load();
    let dom = &g.sequences[0].domains[0];
    let sv = to_sv(&dom.state_vector);
    let seq_len = g.sequences[0].seq.len();

    // Build a whole-sequence label array: domain residues get their labels, the rest Out.
    let enc = encode_domain(&sv);
    let mut labels = vec![Label::Out; seq_len];
    for LabeledRes { si, label } in enc {
        labels[si] = label;
    }
    let domains = decode_sequence(&labels);
    assert_eq!(domains.len(), 1, "expected exactly one segmented domain");
    assert_eq!(domains[0], sv, "segmented+decoded domain != original state vector");
}

/// (id, type) projection of a state vector — structure without si (si round-trip is
/// already proven by `c0_roundtrip_all_golden_domains`).
fn idtyp(sv: &[State]) -> Vec<(u8, StateType)> {
    sv.iter().map(|s| (s.id, s.typ)).collect()
}

/// Labels of one domain's residues in sequence order (encode is already in si order).
fn domain_labels(sv: &[State]) -> Vec<Label> {
    encode_domain(sv).into_iter().map(|r| r.label).collect()
}

/// Two V domains with NO separator (abutting / zero-linker, e.g. tandem) must still
/// segment into two — via the match-id reset rule, not a separator label. This is the
/// hard multi-domain case the plan promises (§3a) and the one that was broken.
#[test]
fn segmentation_splits_abutting_domains() {
    let g = load();
    let sv_a = to_sv(&g.sequences[0].domains[0].state_vector);
    let sv_b = to_sv(&g.sequences[1].domains[0].state_vector);

    let mut labels = domain_labels(&sv_a);
    labels.extend(domain_labels(&sv_b)); // no Out between them

    let domains = decode_sequence(&labels);
    assert_eq!(domains.len(), 2, "abutting domains must split into two, got {}", domains.len());
    assert_eq!(idtyp(&domains[0]), idtyp(&sv_a), "first abutting domain structure wrong");
    assert_eq!(idtyp(&domains[1]), idtyp(&sv_b), "second abutting domain structure wrong");
}

/// Two V domains with a linker (Out run) between them — the easy multi-domain case.
#[test]
fn segmentation_splits_linked_domains() {
    let g = load();
    let sv_a = to_sv(&g.sequences[0].domains[0].state_vector);
    let sv_b = to_sv(&g.sequences[2].domains[0].state_vector);

    let mut labels = domain_labels(&sv_a);
    labels.extend([Label::Linker; 15]); // (G4S)3-style linker
    labels.extend(domain_labels(&sv_b));

    let domains = decode_sequence(&labels);
    assert_eq!(domains.len(), 2, "linked domains must split into two");
    assert_eq!(idtyp(&domains[0]), idtyp(&sv_a));
    assert_eq!(idtyp(&domains[1]), idtyp(&sv_b));
}

/// Synthetic/truncated SVs with a boundary delete only round-trip after normalization
/// (real ANARCI output never has one; window-cut partials in §P2 could).
#[test]
fn strip_boundary_deletes_enables_roundtrip() {
    let sv = vec![
        State { id: 1, typ: StateType::D, si: None }, // leading boundary delete
        State { id: 2, typ: StateType::M, si: Some(0) },
        State { id: 3, typ: StateType::M, si: Some(1) },
        State { id: 4, typ: StateType::D, si: None }, // trailing boundary delete
    ];
    // raw round-trip drops the boundary deletes (asymmetry)...
    assert_ne!(anarci_learned::decode_domain(&encode_domain(&sv)), sv);
    // ...but the canonical (stripped) form round-trips exactly.
    let canon = anarci_learned::strip_boundary_deletes(&sv);
    assert_eq!(anarci_learned::decode_domain(&encode_domain(&canon)), canon);
    assert_eq!(canon.len(), 2);
}

/// Out-of-contract high id must not panic (debug overflow on `hi + 1`).
#[test]
fn decode_domain_no_overflow_on_high_id() {
    let residues = [LabeledRes { si: 0, label: Label::Match(255) }];
    let sv = anarci_learned::decode_domain(&residues);
    assert_eq!(sv.len(), 1);
}

/// Segmentation edge cases: empty, all-outside, and match-less runs yield no domains.
#[test]
fn segmentation_edge_cases() {
    assert_eq!(decode_sequence(&[]).len(), 0, "empty input -> no domains");
    assert_eq!(decode_sequence(&[Label::Out; 20]).len(), 0, "all-Out -> no domains");
    assert_eq!(decode_sequence(&[Label::Constant; 10]).len(), 0, "all-Constant -> no domains");
    // a run of inserts with no match cannot form a domain
    assert_eq!(
        decode_sequence(&[Label::Ins(50), Label::Ins(50), Label::Out]).len(),
        0,
        "insert-only run -> no domain"
    );
}
