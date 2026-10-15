//! Constrained decoding — the "valid by construction" (C2) guarantee.
//!
//! The model emits, per residue, a score over [`NUM_LABELS`] labels. A plain argmax can
//! assign the *same* match position to two residues (silently dropping one in
//! [`crate::decode_domain`]) or place an insert before any match. [`viterbi_decode`] instead
//! finds the best-scoring label path that respects the numbering transition constraints, so
//! the decoded state vector is **always** valid: within a domain match positions strictly
//! increase (no duplicates ⇒ no dropped residues), inserts follow their match, and a
//! position reset (or a non-domain label) starts a new domain. Validity is a theorem here,
//! not a post-processing repair — the contrast with ANARCII's free autoregressive decode.

use crate::labels::Label;

/// Match(1..=128) → 0..128, Ins(1..=128) → 128..256, then Out/Linker/Constant/Hinge/Other.
pub const NUM_LABELS: usize = 261;
const OUT: usize = 256;
const LINKER: usize = 257;
const CONSTANT: usize = 258;
const HINGE: usize = 259;
const OTHER: usize = 260;

/// Canonical label ↔ index mapping (the model's output vocabulary).
pub fn label_to_index(l: Label) -> usize {
    match l {
        Label::Match(p) => (p - 1) as usize,
        Label::Ins(p) => 128 + (p - 1) as usize,
        Label::Out => OUT,
        Label::Linker => LINKER,
        Label::Constant => CONSTANT,
        Label::Hinge => HINGE,
        Label::Other => OTHER,
    }
}

/// Inverse of [`label_to_index`].
pub fn index_to_label(i: usize) -> Label {
    match i {
        0..=127 => Label::Match((i + 1) as u8),
        128..=255 => Label::Ins((i - 128 + 1) as u8),
        OUT => Label::Out,
        LINKER => Label::Linker,
        CONSTANT => Label::Constant,
        HINGE => Label::Hinge,
        _ => Label::Other,
    }
}

/// Penalty for a mid-stream domain reset (a match position that decreases with no
/// separator). Keeps noise from fragmenting one domain into many; a real new domain still
/// wins when its emission scores clear the penalty. Tunable; not a hard block.
const RESET_PENALTY: f32 = -6.0;

/// Hidden state = last match position in the current domain (0 = none / just reset).
/// Returns `(new_state, transition_score)` if emitting `label` from `last_match` is legal.
fn step(last_match: u8, label: Label) -> Option<(u8, f32)> {
    match label {
        Label::Match(q) => {
            if last_match == 0 || q > last_match {
                Some((q, 0.0)) // start a domain (any q, for partials) or continue it
            } else {
                Some((q, RESET_PENALTY)) // q<=last_match: abutting/new domain, discouraged
            }
        }
        // An insert must follow a match and tag a position beyond it; it leaves the last
        // match unchanged. Illegal before any match (can't start a domain on an insert).
        Label::Ins(p) => {
            if last_match > 0 && p > last_match {
                Some((last_match, 0.0))
            } else {
                None
            }
        }
        // A non-domain label ends the current domain (context resets to 0).
        _ => Some((0, 0.0)),
    }
}

const NEG_INF: f32 = f32::NEG_INFINITY;

/// Best-scoring valid label path (Viterbi) over per-residue `logits` (`L × NUM_LABELS`).
///
/// The result, fed to [`crate::decode_sequence`], always yields valid state vectors. Ties
/// break deterministically toward the lower label index (device-independent).
pub fn viterbi_decode(logits: &[Vec<f32>]) -> Vec<Label> {
    let l = logits.len();
    if l == 0 {
        return Vec::new();
    }
    const S: usize = 129; // hidden states 0..=128
    let mut dp = vec![NEG_INF; S];
    dp[0] = 0.0; // start: no domain open
    // back[i][s] = (label index emitted to reach state s at step i, previous state)
    let mut back: Vec<[(u16, u8); S]> = vec![[(0u16, 0u8); S]; l];

    for (i, row) in logits.iter().enumerate() {
        debug_assert_eq!(row.len(), NUM_LABELS);
        let mut ndp = vec![NEG_INF; S];
        for (m, &prev) in dp.iter().enumerate() {
            if prev == NEG_INF {
                continue;
            }
            for (yi, &emit) in row.iter().enumerate() {
                // Sanitize non-finite logits so the DP stays finite (a single NaN/+Inf
                // otherwise breaks the argmax chain and drops residues to `Out`). NaN is
                // garbage → strongly negative; ±Inf keeps its sign. The bound (1e4) dominates
                // any real logit yet is small enough that accumulation over a sequence stays
                // in f32's exact-integer range, so downstream residues keep their signal.
                const CLAMP: f32 = 1e4;
                let emit = if emit.is_nan() { -CLAMP } else { emit.clamp(-CLAMP, CLAMP) };
                if let Some((ns, trans)) = step(m as u8, index_to_label(yi)) {
                    let score = prev + emit + trans;
                    if score > ndp[ns as usize] {
                        ndp[ns as usize] = score;
                        back[i][ns as usize] = (yi as u16, m as u8);
                    }
                }
            }
        }
        dp = ndp;
    }

    // best final state
    let mut best_s = 0usize;
    let mut best = NEG_INF;
    for (s, &v) in dp.iter().enumerate() {
        if v > best {
            best = v;
            best_s = s;
        }
    }
    // backtrack
    let mut labels = vec![Label::Out; l];
    let mut s = best_s;
    for i in (0..l).rev() {
        let (yi, prev) = back[i][s];
        labels[i] = index_to_label(yi as usize);
        s = prev as usize;
    }
    labels
}
