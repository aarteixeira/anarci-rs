//! Constrained decoding of the per-residue label vocabulary.
//!
//! Within a domain, match positions increase and insertion owners never decrease.
//! Tracking the last emitted state, rather than only the last match, prevents the
//! codec from reordering query residues when it groups insertions by owner.
//! A position reset or non-domain label starts a new domain. These are structural
//! label constraints, not evidence of biological numbering accuracy.

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

/// State 0 is outside a domain; 1..=128 are matches and 129..=256 are inserts.
/// Keeping insertion ownership in the state prevents transitions that reorder residues.
fn step(previous: u16, label: Label) -> Option<(u16, f32)> {
    let was_insert = previous > 128;
    let position = if was_insert { previous - 128 } else { previous };
    match label {
        Label::Match(q) => {
            let q = u16::from(q);
            // Match(p) after Ins(p) completes that column; after Match(p) it resets.
            let reset = previous != 0 && if was_insert { q < position } else { q <= position };
            Some((q, if reset { RESET_PENALTY } else { 0.0 }))
        }
        Label::Ins(p) => {
            let p = u16::from(p);
            if previous != 0 && if was_insert { p >= position } else { p > position } {
                Some((128 + p, 0.0))
            } else {
                None
            }
        }
        _ => Some((0, 0.0)),
    }
}

const NEG_INF: f32 = f32::NEG_INFINITY;

/// Best-scoring valid label path (Viterbi) over per-residue `logits` (`L × NUM_LABELS`).
///
/// The result, fed to [`crate::decode_sequence`], always yields valid state vectors. Ties
/// break deterministically in state/label iteration order. Rows must contain NUM_LABELS scores.
pub fn viterbi_decode(logits: &[Vec<f32>]) -> Vec<Label> {
    let l = logits.len();
    if l == 0 {
        return Vec::new();
    }
    const S: usize = 257; // outside plus the last emitted match or insertion owner
    let mut dp = vec![NEG_INF; S];
    dp[0] = 0.0; // start: no domain open
    // back[i][s] = (label index emitted to reach state s at step i, previous state)
    let mut back: Vec<[(u16, u16); S]> = vec![[(0u16, 0u16); S]; l];

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
                if let Some((ns, trans)) = step(m as u16, index_to_label(yi)) {
                    let score = prev + emit + trans;
                    if score > ndp[ns as usize] {
                        ndp[ns as usize] = score;
                        back[i][ns as usize] = (yi as u16, m as u16);
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
