//! Region-completeness annotation (F1a).
//!
//! Given a numbered domain's IMGT state vector, report which of the seven IMGT
//! regions (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4) are present and whether each is
//! fully covered or only partially covered — the signal needed to identify and
//! annotate partial chains (e.g. an FR3-CDR3-FR4-only fragment).
//!
//! Scheme-independent: coverage is computed from the IMGT match-state ids (always
//! `1..=128`, see [`crate::types::State`]), not from the output numbering, so the
//! annotation is identical regardless of the chosen output scheme.

use crate::types::{State, StateType};

/// Coverage status of one IMGT region within a numbered domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegionStatus {
    /// No residue of this region was sequenced.
    Absent,
    /// The region is covered, but the fragment starts and/or ends inside it.
    Partial,
    /// The coverage spans the region's full IMGT extent.
    Complete,
}

impl RegionStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            RegionStatus::Absent => "absent",
            RegionStatus::Partial => "partial",
            RegionStatus::Complete => "complete",
        }
    }
}

/// The seven IMGT regions in N→C order with their inclusive IMGT position ranges.
/// FR1 1–26, CDR1 27–38, FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, FR4 118–128.
/// These match the `region_string` used by `number_imgt` (see `schemes/mod.rs`).
pub const IMGT_REGIONS: [(&str, i32, i32); 7] = [
    ("fr1", 1, 26),
    ("cdr1", 27, 38),
    ("fr2", 39, 55),
    ("cdr2", 56, 65),
    ("fr3", 66, 104),
    ("cdr3", 105, 117),
    ("fr4", 118, 128),
];

/// Region-completeness annotation for one numbered domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionAnnotation {
    /// Status per region, in [`IMGT_REGIONS`] order.
    pub statuses: [RegionStatus; 7],
    /// Inclusive `[min, max]` IMGT match position carrying a residue, or `None`
    /// when the domain has no residue-bearing match state (degenerate).
    pub covered_imgt: Option<(i32, i32)>,
}

impl RegionAnnotation {
    /// `(region_name, status)` pairs in N→C order.
    pub fn pairs(&self) -> impl Iterator<Item = (&'static str, RegionStatus)> + '_ {
        IMGT_REGIONS
            .iter()
            .zip(self.statuses.iter())
            .map(|(&(name, _, _), &st)| (name, st))
    }
}

/// Annotate which IMGT regions a domain's state vector covers.
///
/// Coverage is taken from match states that carry a residue (`StateType::M` with
/// `si.is_some()`); deletions (no residue) and inserts (no backbone position) don't
/// define coverage. `validate_numbering` guarantees the numbered residues are a
/// contiguous segment, so the covered match positions form a single `[min, max]`
/// interval, and a region `[lo, hi]` is:
/// * `Absent`   if the interval doesn't overlap it (`hi < lo` or `min > hi`),
/// * `Complete` if the interval brackets it (`min <= lo && max >= hi`),
/// * `Partial`  otherwise (the fragment starts or ends inside the region).
pub fn annotate_regions(state_vector: &[State]) -> RegionAnnotation {
    let mut min_pos: Option<i32> = None;
    let mut max_pos: Option<i32> = None;
    for s in state_vector {
        if s.typ == StateType::M && s.si.is_some() {
            let p = s.id as i32;
            min_pos = Some(min_pos.map_or(p, |m| m.min(p)));
            max_pos = Some(max_pos.map_or(p, |m| m.max(p)));
        }
    }

    let mut statuses = [RegionStatus::Absent; 7];
    let covered_imgt = match (min_pos, max_pos) {
        (Some(lo), Some(hi)) => {
            for (i, &(_, r_lo, r_hi)) in IMGT_REGIONS.iter().enumerate() {
                statuses[i] = if hi < r_lo || lo > r_hi {
                    RegionStatus::Absent
                } else if lo <= r_lo && hi >= r_hi {
                    RegionStatus::Complete
                } else {
                    RegionStatus::Partial
                };
            }
            Some((lo, hi))
        }
        _ => None,
    };
    RegionAnnotation { statuses, covered_imgt }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{State, StateType};
    use RegionStatus::*;

    /// Build a state vector of residue-bearing match states for IMGT positions
    /// `lo..=hi` inclusive (one residue per position).
    fn match_span(lo: i32, hi: i32) -> Vec<State> {
        (lo..=hi)
            .enumerate()
            .map(|(i, p)| State { id: p as u8, typ: StateType::M, si: Some(i) })
            .collect()
    }

    fn statuses(lo: i32, hi: i32) -> [RegionStatus; 7] {
        annotate_regions(&match_span(lo, hi)).statuses
    }

    #[test]
    fn full_domain_all_complete() {
        let a = annotate_regions(&match_span(1, 128));
        assert_eq!(a.statuses, [Complete; 7]);
        assert_eq!(a.covered_imgt, Some((1, 128)));
    }

    #[test]
    fn fr3_cdr3_fr4_fragment() {
        // The canonical partial: sequenced only from FR3 (66) to the end.
        assert_eq!(statuses(66, 128), [Absent, Absent, Absent, Absent, Complete, Complete, Complete]);
    }

    #[test]
    fn fr2_onward_fragment() {
        assert_eq!(statuses(39, 128), [Absent, Absent, Complete, Complete, Complete, Complete, Complete]);
    }

    #[test]
    fn n_terminal_fragment_partial_fr2() {
        // 1..=53 covers FR1+CDR1 fully, ends inside FR2 (39–55) → FR2 partial.
        assert_eq!(statuses(1, 53), [Complete, Complete, Partial, Absent, Absent, Absent, Absent]);
        assert_eq!(annotate_regions(&match_span(1, 53)).covered_imgt, Some((1, 53)));
    }

    #[test]
    fn mid_fr3_start_partial() {
        // Starts inside FR3 (66–104) at 80 → FR3 partial; CDR3/FR4 complete.
        assert_eq!(statuses(80, 128), [Absent, Absent, Absent, Absent, Partial, Complete, Complete]);
    }

    #[test]
    fn cdr3_alone_only_cdr3() {
        // 105..=117 is exactly CDR3.
        assert_eq!(statuses(105, 117), [Absent, Absent, Absent, Absent, Absent, Complete, Absent]);
    }

    #[test]
    fn deletes_and_inserts_do_not_count_as_coverage() {
        // A delete state (si=None) at 200 and an insert state must not extend coverage.
        let mut sv = match_span(27, 38); // CDR1 only
        sv.push(State { id: 50, typ: StateType::M, si: None }); // delete in FR2 → not covered
        sv.push(State { id: 111, typ: StateType::I, si: Some(99) }); // insert → not a backbone pos
        let a = annotate_regions(&sv);
        assert_eq!(a.covered_imgt, Some((27, 38)));
        assert_eq!(a.statuses, [Absent, Complete, Absent, Absent, Absent, Absent, Absent]);
    }

    #[test]
    fn empty_is_all_absent() {
        let a = annotate_regions(&[]);
        assert_eq!(a.statuses, [Absent; 7]);
        assert_eq!(a.covered_imgt, None);
    }
}
