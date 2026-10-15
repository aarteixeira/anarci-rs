//! Synthetic training-data augmentation (P2), built entirely from codec operations.
//!
//! We own the oracle (the exact engine), so a labeled V domain — its residues plus their
//! [`Label`]s — is free. This module composes those into the *hard* inputs ANARCII refuses
//! or ignores and the HMM engine handles piecemeal, with **exact** labels:
//!
//! * **partials** — restrict a domain to an IMGT-position window (FR3-CDR3-FR4, …) or a
//!   ragged residue window; the retained residues keep their oracle labels.
//! * **multi-V** (scFv / tandem / bispecific) — concatenate 2–4 domains with a real or
//!   zero-length linker; each domain keeps its own numbering.
//! * **constant-context / topology** — flank a domain with constant/hinge/tag residues.
//! * **junk** — residues labeled entirely `Out` (a negative; no domain).
//!
//! Because a [`Label::Match`]/[`Label::Ins`] carries an *IMGT position*, not a sequence
//! index, concatenation needs no remapping — the labels compose directly. Every builder
//! returns a [`LabeledSeq`] whose `labels[i]` is the target for residue `seq[i]`.

use crate::labels::{decode_sequence, encode_domain, Label};
use anarci_core::State;

/// A residue sequence paired with its per-residue target labels (`labels.len() == seq.len()`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LabeledSeq {
    pub seq: Vec<u8>,
    pub labels: Vec<Label>,
}

impl LabeledSeq {
    /// A single V domain: its residues (the contiguous substring it covers) and their
    /// `Match`/`Ins` labels from the oracle state vector.
    pub fn from_domain(seq: &[u8], sv: &[State]) -> LabeledSeq {
        let enc = encode_domain(sv);
        let mut ls = LabeledSeq { seq: Vec::with_capacity(enc.len()), labels: Vec::with_capacity(enc.len()) };
        for r in enc {
            // oracle SVs always index within `seq`; guard so a malformed/hand-built SV
            // skips rather than panics.
            debug_assert!(r.si < seq.len(), "state index {} out of range (seq len {})", r.si, seq.len());
            if let Some(&aa) = seq.get(r.si) {
                ls.seq.push(aa);
                ls.labels.push(r.label);
            }
        }
        ls
    }

    /// A non-domain region: `residues` all tagged with one label (`Linker`, `Constant`,
    /// `Hinge`, `Other`, or `Out`). Real residues are supplied by the caller (P4 sources
    /// genuine constant/hinge sequences); the label is what matters for segmentation.
    pub fn region(residues: &[u8], label: Label) -> LabeledSeq {
        LabeledSeq { seq: residues.to_vec(), labels: vec![label; residues.len()] }
    }

    /// A junk negative: every residue is `Out` (no V domain).
    pub fn junk(residues: &[u8]) -> LabeledSeq {
        LabeledSeq::region(residues, Label::Out)
    }

    /// Restrict a domain segment to residues whose IMGT position lies in `lo..=hi`
    /// (a partial, e.g. `restrict_positions(66, 128)` = FR3-CDR3-FR4). Positions are
    /// monotonic, so the result is a contiguous fragment.
    pub fn restrict_positions(&self, lo: u8, hi: u8) -> LabeledSeq {
        self.filter(|label| match label {
            Label::Match(p) | Label::Ins(p) => (lo..=hi).contains(&p),
            _ => false,
        })
    }

    /// Restrict to a raw residue window `[start, end)` — a *ragged* cut that need not fall
    /// on a region boundary (models real 5'-RACE / missing-density fragments).
    pub fn restrict_window(&self, start: usize, end: usize) -> LabeledSeq {
        let end = end.min(self.seq.len());
        let start = start.min(end);
        LabeledSeq { seq: self.seq[start..end].to_vec(), labels: self.labels[start..end].to_vec() }
    }

    fn filter(&self, keep: impl Fn(Label) -> bool) -> LabeledSeq {
        let mut ls = LabeledSeq { seq: Vec::new(), labels: Vec::new() };
        for (&aa, &label) in self.seq.iter().zip(&self.labels) {
            if keep(label) {
                ls.seq.push(aa);
                ls.labels.push(label);
            }
        }
        ls
    }

    /// Decode this labeled sequence into its per-domain state vectors (N→C) via the codec.
    pub fn domains(&self) -> Vec<Vec<State>> {
        decode_sequence(&self.labels)
    }
}

/// Concatenate segments in order (residues and labels both), e.g.
/// `concat(&[vh, linker, vl])` for an scFv, or `concat(&[vh, vl])` for a zero-linker
/// abutting tandem. No remapping needed — labels carry IMGT positions, not indices.
pub fn concat(segments: &[&LabeledSeq]) -> LabeledSeq {
    let n: usize = segments.iter().map(|s| s.seq.len()).sum();
    let mut out = LabeledSeq { seq: Vec::with_capacity(n), labels: Vec::with_capacity(n) };
    for s in segments {
        out.seq.extend_from_slice(&s.seq);
        out.labels.extend_from_slice(&s.labels);
    }
    out
}

/// Safe multi-domain composer for the data pipeline: like [`concat`] but inserts a single
/// `Out` separator residue between consecutive domain segments so abutting domains — in
/// particular abutting **partials**, which `decode_sequence`'s reset rule cannot split —
/// cannot silently fuse into one plausible-but-wrong domain (a training-label poison).
///
/// Full-domain tandems (which already reset) are unaffected in meaning; the extra `Out` is
/// a one-residue inter-domain boundary the model learns anyway. Use this for **all**
/// composed training data; reserve raw [`concat`] for explicit zero-linker probes where
/// every segment is a full domain that provably resets.
pub fn concat_domains(segments: &[&LabeledSeq]) -> LabeledSeq {
    let mut out = LabeledSeq { seq: Vec::new(), labels: Vec::new() };
    for (i, s) in segments.iter().enumerate() {
        if i > 0 {
            out.seq.push(b'G'); // filler residue for the boundary
            out.labels.push(Label::Out); // separator: guarantees the next domain starts fresh
        }
        out.seq.extend_from_slice(&s.seq);
        out.labels.extend_from_slice(&s.labels);
    }
    out
}

/// A classic Gly-Ser scFv linker of the given repeat count: `(G4S)n`, labeled `Linker`.
pub fn gly_ser_linker(repeats: usize) -> LabeledSeq {
    let residues: Vec<u8> = "GGGGS".bytes().cycle().take(repeats * 5).collect();
    LabeledSeq::region(&residues, Label::Linker)
}
