//! Learned-engine label codec.
//!
//! The learned engine predicts, per residue, a [`Label`]; the existing byte-exact
//! numbering code then consumes the reconstructed state vector. This module is the
//! lossless map between an ANARCI state vector and its per-residue labels.
//!
//! Label design (verified lossless on all 1013 golden domains — see `tests/codec.rs`):
//! * `Match(p)` — residue occupies IMGT match state `p` (1..=128).
//! * `Ins(p)`   — residue is an insertion, **position-tagged** with its owning id.
//!   In ANARCI an insert takes the id of the *next* match column
//!   (`anarci-core/src/align.rs`), so the id is stored, not inferred: a bare "inherit
//!   the preceding match" rule is lossy on ~25% of domains (CDR3 111/112 split and the
//!   long-CDR3 id-116 block).
//! * `Out`/`Linker`/`Constant`/`Hinge`/`Other` — residue is outside any V/TCR-V domain
//!   (construct-topology classes; only appear in multi-domain/construct inputs).
//!
//! Delete states carry no residue: they are dropped by [`encode_domain`] and
//! reconstructed by [`decode_domain`] (a span position is a delete only if occupied by
//! neither a match nor an insert — the rule that avoids inventing spurious deletes).

use anarci_core::{State, StateType};
use std::collections::BTreeMap;
use std::ops::Bound;

/// Per-residue label predicted by the learned engine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Label {
    /// Residue occupies IMGT match state `p` (1..=128).
    Match(u8),
    /// Insertion, position-tagged with its owning id (the next match column).
    Ins(u8),
    /// Outside any domain — generic / unknown (default non-domain class).
    Out,
    /// Inter-domain linker (e.g. scFv (G4S)ₙ).
    Linker,
    /// Constant-region residue (CH1 / CL / CH2 / CH3).
    Constant,
    /// Hinge residue.
    Hinge,
    /// Tag / non-antibody / junk.
    Other,
}

impl Label {
    /// Whether this residue is part of a V/TCR-V domain (a `Match` or `Ins`).
    pub fn is_domain(self) -> bool {
        matches!(self, Label::Match(_) | Label::Ins(_))
    }
}

/// A residue tagged with its sequence index and label. `si` is the 0-based index into
/// the input sequence; the model does not predict it — it is the residue's position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LabeledRes {
    pub si: usize,
    pub label: Label,
}

/// Encode one domain's state vector into its per-residue labels.
///
/// Only match/insert states carry a residue and become labels, in state-vector order;
/// delete states are dropped (reconstructed by [`decode_domain`]).
pub fn encode_domain(sv: &[State]) -> Vec<LabeledRes> {
    sv.iter()
        .filter_map(|s| {
            let si = s.si?;
            let label = match s.typ {
                StateType::M => Label::Match(s.id),
                StateType::I => Label::Ins(s.id),
                StateType::D => return None,
            };
            Some(LabeledRes { si, label })
        })
        .collect()
}

/// Decode one domain's labeled residues back into a full state vector, reconstructing
/// delete states. Inverse of [`encode_domain`].
///
/// Non-domain labels (`Out`/`Linker`/…) are ignored here — segmentation into per-domain
/// residue runs happens in [`decode_sequence`] before this is called.
pub fn decode_domain(residues: &[LabeledRes]) -> Vec<State> {
    let mut m_si: BTreeMap<u8, usize> = BTreeMap::new();
    // inserts owned by each id, in stream order (preserves si order within a run).
    let mut ins_by_id: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for r in residues {
        match r.label {
            Label::Match(p) => {
                m_si.insert(p, r.si);
            }
            Label::Ins(p) => ins_by_id.entry(p).or_default().push(r.si),
            _ => {}
        }
    }
    if m_si.is_empty() {
        return Vec::new();
    }
    let lo = *m_si.keys().next().unwrap();
    let hi = *m_si.keys().next_back().unwrap();

    let mut sv = Vec::new();
    // inserts owned by ids below the match span (leading — e.g. synthetic partials).
    for (&id, sis) in ins_by_id.range(..lo) {
        for &si in sis {
            sv.push(State { id, typ: StateType::I, si: Some(si) });
        }
    }
    for id in lo..=hi {
        if let Some(sis) = ins_by_id.get(&id) {
            for &si in sis {
                sv.push(State { id, typ: StateType::I, si: Some(si) });
            }
        }
        if let Some(&si) = m_si.get(&id) {
            sv.push(State { id, typ: StateType::M, si: Some(si) });
        } else if !ins_by_id.contains_key(&id) {
            // span position occupied by neither a match nor an insert -> delete.
            sv.push(State { id, typ: StateType::D, si: None });
        }
    }
    // inserts owned by ids above the span (trailing). Bound::Excluded(hi) avoids a
    // `hi + 1` overflow if an out-of-contract id (u8::MAX) ever reaches here.
    for (&id, sis) in ins_by_id.range((Bound::Excluded(hi), Bound::Unbounded)) {
        for &si in sis {
            sv.push(State { id, typ: StateType::I, si: Some(si) });
        }
    }
    sv
}

/// Remove delete states before the first match / after the last match, putting a state
/// vector into the canonical form the codec round-trips (a domain span is defined by its
/// matches; a boundary delete is meaningless). Real ANARCI output is already canonical;
/// this is for **synthetic/truncated** state vectors (partial augmentation, §P2) where a
/// window cut could otherwise leave a leading/trailing delete that `encode`/`decode`
/// disagree on.
pub fn strip_boundary_deletes(sv: &[State]) -> Vec<State> {
    let first = sv.iter().position(|s| s.typ != StateType::D);
    let last = sv.iter().rposition(|s| s.typ != StateType::D);
    match (first, last) {
        (Some(f), Some(l)) => sv[f..=l].to_vec(),
        _ => Vec::new(),
    }
}

/// Segment a whole-sequence label array into per-domain state vectors, N→C.
///
/// A new domain begins at (a) the first domain label after a non-domain label
/// (linker/constant/…), **or** (b) a match-position **reset** — a `Match(q)` that does
/// not strictly exceed the previous match id in the current domain. Rule (b) is what
/// separates **abutting / zero-linker** domains (tandem scFv, bispecifics), where there
/// is no separator residue: the stream is `…Match(128) Match(1)…`, and within a real
/// domain match ids strictly increase, so a non-increase marks a domain boundary.
/// Without (b) abutting domains would merge into one run and collide in `decode_domain`.
///
/// **Known limit of rule (b):** it fires only on a position *decrease*. Two abutting
/// *partials* where the second starts at a higher position than the first ended (e.g. a
/// 1–100 fragment directly followed by a 105–128 fragment, no separator) do not reset and
/// silently merge. Full domains (the realistic tandem/scFv case) always reset because the
/// second starts near position 1; synthetic abutting *partials* must carry a separator
/// label. The general fix is the constrained decoder's explicit domain-boundary transition
/// (P3), not this segmenter.
///
/// Returns the decoded state vector for each domain in order of first residue.
pub fn decode_sequence(labels: &[Label]) -> Vec<Vec<State>> {
    let mut domains: Vec<Vec<State>> = Vec::new();
    let mut run: Vec<LabeledRes> = Vec::new();
    let mut last_match: Option<u8> = None;

    for (si, &label) in labels.iter().enumerate() {
        match label {
            Label::Match(p) => {
                // (b) position reset with no separator -> flush current (abutting) domain.
                if last_match.is_some_and(|lm| p <= lm) {
                    flush(&mut run, &mut domains);
                }
                run.push(LabeledRes { si, label });
                last_match = Some(p);
            }
            // ponytail: a leading Ins of an abutting next domain (id < last match) is a
            // non-real shape (domains start with FR1 matches); we don't reset on inserts.
            Label::Ins(_) => run.push(LabeledRes { si, label }),
            _ => {
                // (a) non-domain separator label.
                flush(&mut run, &mut domains);
                last_match = None;
            }
        }
    }
    flush(&mut run, &mut domains);
    domains
}

/// Decode and append the current run as a domain (if non-empty), then clear it.
fn flush(run: &mut Vec<LabeledRes>, domains: &mut Vec<Vec<State>>) {
    if !run.is_empty() {
        let d = decode_domain(run);
        if !d.is_empty() {
            domains.push(d);
        }
        run.clear();
    }
}
