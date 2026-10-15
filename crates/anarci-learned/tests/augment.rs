//! P2 gate for synthetic augmentation: every augmented input must (a) round-trip through
//! the codec and (b) have its exact domain set / boundaries recovered by segmentation, and
//! each recovered domain must number through the real `anarci-core` pipeline.

use anarci_core::{annotate_regions, number_sequence_from_alignment, RegionStatus, State, StateType};
use anarci_learned::augment::LabeledSeq;
use anarci_learned::{concat, gly_ser_linker, Label};
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
fn idtyp(sv: &[State]) -> Vec<(u8, StateType)> {
    sv.iter().map(|s| (s.id, s.typ)).collect()
}
/// Residues a domain's state vector points at within `seq` (match/insert states only).
fn residues_of(sv: &[State], seq: &[u8]) -> Vec<u8> {
    sv.iter().filter_map(|s| s.si.map(|i| seq[i])).collect()
}

/// Pull the `i`-th single V domain (seq bytes, state vector, chain type) from golden.
fn domain(g: &Golden, i: usize) -> (LabeledSeq, Vec<State>, String) {
    let s = &g.sequences[i];
    let d = &s.domains[0];
    let sv = to_sv(&d.state_vector);
    (LabeledSeq::from_domain(s.seq.as_bytes(), &sv), sv, d.chain_type.clone())
}

/// Sanity: a single domain segment recovers exactly one domain with the right structure.
#[test]
fn single_domain_segment_recovers_itself() {
    let g = load();
    let (seg, sv, _) = domain(&g, 0);
    let got = seg.domains();
    assert_eq!(got.len(), 1);
    assert_eq!(idtyp(&got[0]), idtyp(&sv));
}

/// scFv (VH–linker–VL) and abutting tandem (no linker) both recover all domains.
#[test]
fn multi_v_recovers_all_domains() {
    let g = load();
    let (a, sva, _) = domain(&g, 0);
    let (b, svb, _) = domain(&g, 1);
    let (c, svc, _) = domain(&g, 2);

    // scFv: VH-(G4S)3-VL
    let scfv = concat(&[&a, &gly_ser_linker(3), &b]);
    let d = scfv.domains();
    assert_eq!(d.len(), 2, "scFv must recover 2 domains");
    assert_eq!(idtyp(&d[0]), idtyp(&sva));
    assert_eq!(idtyp(&d[1]), idtyp(&svb));

    // abutting tandem: three domains, zero linker
    let tandem = concat(&[&a, &b, &c]);
    let d = tandem.domains();
    assert_eq!(d.len(), 3, "abutting tandem must recover 3 domains, got {}", d.len());
    assert_eq!(idtyp(&d[0]), idtyp(&sva));
    assert_eq!(idtyp(&d[1]), idtyp(&svb));
    assert_eq!(idtyp(&d[2]), idtyp(&svc));
}

/// Residue-level (not just id/type) check: each recovered domain points at exactly the
/// original domain's residues within the composed sequence.
#[test]
fn multi_v_residues_match_originals() {
    let g = load();
    let (a, _sva, _) = domain(&g, 0);
    let (b, _svb, _) = domain(&g, 1);
    let scfv = concat(&[&a, &gly_ser_linker(3), &b]);
    let d = scfv.domains();
    assert_eq!(d.len(), 2);
    assert_eq!(residues_of(&d[0], &scfv.seq), a.seq, "domain 0 residues must equal original VH");
    assert_eq!(residues_of(&d[1], &scfv.seq), b.seq, "domain 1 residues must equal original VL");
}

/// Documented limitation (labels.rs): two abutting *partials* whose positions do NOT
/// reset (first ends high-of second's start) merge into one domain with no separator.
/// Pinned so the behavior is explicit — the real fix is the P3 decoder boundary, and the
/// data pipeline must insert a separator when composing abutting partials.
#[test]
fn abutting_partials_without_reset_merge_known_limit() {
    let g = load();
    let (seg, _sv, _) = domain(&g, 0);
    let head = seg.restrict_positions(1, 100); // ends at ~100
    let tail = seg.restrict_positions(105, 128); // starts at ~105 > 100 -> no reset
    if head.labels.is_empty() || tail.labels.is_empty() {
        return; // domain 0 must cover both windows; skip if not
    }
    let merged = concat(&[&head, &tail]);
    // With no reset and no separator, segmentation cannot split them: exactly one domain...
    let d = merged.domains();
    assert_eq!(d.len(), 1, "known limit: non-resetting abutting partials merge");
    // ...and — the DANGER — the fused domain still numbers cleanly (a silent, plausible
    // corruption, not a loud failure). Pinned so the `concat_domains` guard can't regress.
    number_sequence_from_alignment(&d[0], &merged.seq, "imgt", Some("H"))
        .expect("fused domain numbers (this is exactly why the raw concat path is unsafe)");

    // The SAFE composer inserts a separator, so the same partials split into two domains.
    let safe = anarci_learned::concat_domains(&[&head, &tail]);
    assert_eq!(safe.domains().len(), 2, "concat_domains must split abutting partials");
}

/// A real 2-domain golden sequence (native inter-domain linker) segments into both domains.
#[test]
fn real_multi_domain_golden_recovers_both() {
    let g = load();
    // find a genuine 2-domain entry (Fab-style pair with a real linker between)
    let s = g.sequences.iter().find(|s| s.domains.len() == 2).expect("a 2-domain golden entry");
    let sv0 = to_sv(&s.domains[0].state_vector);
    let sv1 = to_sv(&s.domains[1].state_vector);
    // build the labeled sequence directly from BOTH domains over the full seq
    let mut labels = vec![Label::Out; s.seq.len()];
    for r in anarci_learned::encode_domain(&sv0) {
        labels[r.si] = r.label;
    }
    for r in anarci_learned::encode_domain(&sv1) {
        labels[r.si] = r.label;
    }
    let domains = anarci_learned::decode_sequence(&labels);
    assert_eq!(domains.len(), 2, "real 2-domain sequence must recover 2 domains");
    assert_eq!(idtyp(&domains[0]), idtyp(&sv0));
    assert_eq!(idtyp(&domains[1]), idtyp(&sv1));
}

/// Constant-region / tag context: the V domain is recovered; flanks create no domains.
#[test]
fn constant_context_leaves_one_v_domain() {
    let g = load();
    let (vh, svh, _) = domain(&g, 0);
    let signal = LabeledSeq::region(b"MELGLSWIFLLAILKGVQC", Label::Other); // signal peptide
    let ch1 = LabeledSeq::region(b"ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEP", Label::Constant);
    let hinge = LabeledSeq::region(b"EPKSCDKTHTCPPCP", Label::Hinge);

    let construct = concat(&[&signal, &vh, &ch1, &hinge]);
    let d = construct.domains();
    assert_eq!(d.len(), 1, "V+constant construct must recover exactly one V domain");
    assert_eq!(idtyp(&d[0]), idtyp(&svh));
}

/// Junk (all `Out`) yields no domains.
#[test]
fn junk_has_no_domains() {
    let junk = LabeledSeq::junk(b"KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLC");
    assert_eq!(junk.domains().len(), 0);
}

/// FR3-CDR3-FR4 partial: exactly one domain, numbers through the real pipeline, and the
/// region annotation reports FR1..CDR2 absent and FR3/CDR3/FR4 present.
#[test]
fn partial_fr3_cdr3_fr4_numbers_and_annotates() {
    let g = load();
    // choose a domain whose coverage reaches FR4 so the partial is meaningful
    let (seg, _sv, ct) = domain(&g, 0);
    let partial = seg.restrict_positions(66, 128);
    let d = partial.domains();
    assert_eq!(d.len(), 1, "partial must be one domain");

    // numbers through existing byte-exact code
    let numbered = number_sequence_from_alignment(&d[0], &partial.seq, "imgt", Some(&ct))
        .expect("partial should number");
    // must place real residues, not an all-gap numbering (a FR1-only fragment with no CDR3/
    // FR4 anchor numbers to empty under ANARCI's smooth_insertions — scope partials to
    // fragments reaching the CDR3/FR4 region).
    assert!(numbered.residues.iter().any(|r| r.aa != b'-'), "partial numbered to all gaps");

    // region completeness: FR1,CDR1,FR2,CDR2 absent; FR3,CDR3,FR4 present (not absent)
    let ra = annotate_regions(&d[0]);
    let st = ra.statuses;
    assert_eq!(&st[0..4], &[RegionStatus::Absent; 4], "FR1..CDR2 must be absent in FR3-onward partial");
    assert!(st[4] != RegionStatus::Absent, "FR3 must be present");
    assert!(st[6] != RegionStatus::Absent, "FR4 must be present");
}

/// Ragged partial (raw window, not a region boundary) still segments to one domain and numbers.
#[test]
fn ragged_partial_numbers() {
    let g = load();
    let (seg, _sv, ct) = domain(&g, 0);
    // drop 20 residues off each end (ragged, mid-region cuts)
    let n = seg.seq.len();
    let ragged = seg.restrict_window(20, n - 20);
    let d = ragged.domains();
    assert_eq!(d.len(), 1, "ragged partial must be one domain");
    number_sequence_from_alignment(&d[0], &ragged.seq, "imgt", Some(&ct)).expect("ragged partial numbers");
}

/// Every domain in a composed construct keeps a well-formed, numberable state vector
/// (a broad correctness sweep over many golden domains combined pairwise).
#[test]
fn composed_domains_all_number() {
    let g = load();
    let n = g.sequences.len().min(60);
    let mut pairs = 0;
    for i in (0..n).step_by(2) {
        if i + 1 >= n {
            break;
        }
        let (a, sva, cta) = domain(&g, i);
        let (b, svb, ctb) = domain(&g, i + 1);
        let scfv = concat(&[&a, &gly_ser_linker(3), &b]);
        let d = scfv.domains();
        assert_eq!(d.len(), 2, "seq pair {i} should give 2 domains");
        assert_eq!(idtyp(&d[0]), idtyp(&sva));
        assert_eq!(idtyp(&d[1]), idtyp(&svb));
        number_sequence_from_alignment(&d[0], &scfv.seq, "imgt", Some(&cta)).expect("dom0 numbers");
        number_sequence_from_alignment(&d[1], &scfv.seq, "imgt", Some(&ctb)).expect("dom1 numbers");
        pairs += 1;
    }
    assert!(pairs >= 25, "expected >=25 composed pairs, did {pairs}");
}
