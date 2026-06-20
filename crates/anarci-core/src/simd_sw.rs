//! SIMD (inter-sequence batched) score-only Smith-Waterman, bit-exact with the
//! scalar [`crate::sw::local_score`].
//!
//! One query is aligned against up to 8 target genes at once: each i16 lane of an
//! `i16x8` holds the DP cell of a *different* target. This is the SWIPE/Opal
//! "inter-sequence" model — ideal here because the germline DB is many short
//! (~120 aa) references. The recurrence, gap penalties (`open=11, extend=1`) and
//! BLOSUM62 substitution scores are identical to the scalar Gotoh kernel, so the
//! returned score for every lane equals `local_score(query, target_lane)` exactly.
//!
//! Why i16 is safe: a local alignment over ~120-aa genes peaks near ~1300
//! (well under i16::MAX = 32767); the sentinel "negative infinity" is far from
//! i16::MIN even after repeated `-extend`, so no lane ever overflows or wraps.
//!
//! Targets in a batch are padded to the longest one. Padding columns are masked
//! per lane: their H is forced to 0 (a local alignment can always "restart" at 0)
//! and their E gap-state to the sentinel, so a shorter target's score is the same
//! as if it were aligned alone — no spurious carry from a neighbour's length.

use crate::sw::{GAP_EXTEND, GAP_OPEN};
use wide::{i16x8, CmpLt};

/// SIMD lane count (i16x8).
pub const LANES: usize = 8;

/// Sentinel "negative infinity" for the affine-gap states. Chosen so that
/// repeated `- GAP_EXTEND` across any realistic gene length stays well above
/// i16::MIN (no wrap), while still losing every `max` against a real score.
const NEG: i16 = -30_000;

/// BLOSUM62 row for residue byte `a` (uppercase letter), as i16, indexed by
/// `target_byte - b'A'` (0..26). Non-letters fall through to the -4 column at
/// index 26 of the padded row (matches the scalar `score()` fallback exactly).
#[inline]
fn blosum_row_i16(a: u8) -> [i16; 27] {
    let ai = (a.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
    let mut row = [-4i16; 27]; // index 26 = non-letter target -> -4
    if ai < 26 {
        for (bi, slot) in row.iter_mut().take(26).enumerate() {
            *slot = crate::sw::BLOSUM62[ai][bi] as i16;
        }
    } else {
        // Non-letter query: scalar `score` returns -4 for every target. Keep row at -4.
        for slot in row.iter_mut().take(26) {
            *slot = -4;
        }
    }
    row
}

/// Substitution-score index of a target residue byte (column into `blosum_row_i16`):
/// `b - b'A'` for letters, else 26 (the -4 fallback column).
#[inline]
fn sub_index(b: u8) -> usize {
    let bi = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
    if bi < 26 {
        bi
    } else {
        26
    }
}

/// Local-alignment scores of `query` against up to 8 `targets` (one per lane),
/// bit-exact with `local_score`. `targets.len()` must be in `1..=8`; lanes beyond
/// `targets.len()` are unused (their output is meaningless and not returned).
///
/// Returns an `[i16; 8]`; only the first `targets.len()` entries are valid.
///
/// Layout (SWIPE-style query profile): once per batch we build, for each possible
/// query-residue group (27 rows incl. the non-letter fallback) and each column j,
/// the 8-lane substitution vector `prof[group][j]`. Padding lanes get -4 there but
/// are also masked out per column by a precomputed `invalid[j]` vector (H forced to
/// 0, E to the sentinel). The inner DP cell is then pure SIMD: one profile load,
/// adds, maxes, and a blend — no per-cell scalar gather.
pub fn local_score_batch(query: &[u8], targets: &[&[u8]]) -> [i16; LANES] {
    let mut scratch = Scratch::default();
    score_one_batch(query, targets, &mut scratch)
}

/// Score MANY targets against one query with the SIMD kernel, reusing scratch
/// buffers across 8-gene batches (so a 1000-gene scan allocates O(n_max) once, not
/// per batch). Returns one i16 score per target, in order. Bit-exact with
/// `local_score`. This is the hot entry point for germline e-value scans.
pub fn local_score_many(query: &[u8], targets: &[&[u8]]) -> Vec<i16> {
    let mut out = Vec::with_capacity(targets.len());
    if targets.is_empty() {
        return out;
    }
    let mut scratch = Scratch::default();
    for chunk in targets.chunks(LANES) {
        let scores = score_one_batch(query, chunk, &mut scratch);
        out.extend_from_slice(&scores[..chunk.len()]);
    }
    out
}

/// Reusable per-scan scratch (grown to the largest batch seen, then reused).
#[derive(Default)]
struct Scratch {
    invalid: Vec<i16x8>,
    col_groups: Vec<[u8; LANES]>,
    prof: Vec<i16x8>,
    h_prev: Vec<i16x8>,
    h_cur: Vec<i16x8>,
    e: Vec<i16x8>,
}

/// Core: one 8-lane batch, using (and resizing) the caller's scratch buffers.
fn score_one_batch(query: &[u8], targets: &[&[u8]], s: &mut Scratch) -> [i16; LANES] {
    debug_assert!(!targets.is_empty() && targets.len() <= LANES);
    let m = query.len();
    let n_max = targets.iter().map(|t| t.len()).max().unwrap_or(0);
    if m == 0 || n_max == 0 {
        return [0; LANES];
    }

    let zero = i16x8::new([0; LANES]);
    let neg = i16x8::new([NEG; LANES]);
    let gap_open = GAP_OPEN as i16;
    let gap_ext = GAP_EXTEND as i16;

    // Per-lane real length and per-column padding mask.
    let mut len_lane = [0i16; LANES];
    for (l, t) in targets.iter().enumerate() {
        len_lane[l] = t.len() as i16;
    }
    let len_vec = i16x8::new(len_lane);
    s.invalid.clear();
    s.col_groups.clear();
    for j in 1..=n_max {
        let jvec = i16x8::new([j as i16; LANES]);
        s.invalid.push(len_vec.cmp_lt(jvec)); // ones where len < j (padding)
        let mut g = [26u8; LANES];
        for (l, t) in targets.iter().enumerate() {
            if j - 1 < t.len() {
                g[l] = sub_index(t[j - 1]) as u8;
            }
        }
        s.col_groups.push(g);
    }

    // Query profile: prof[group * n_max + (j-1)] = i16x8 of BLOSUM[group][target_lane[j]].
    // Only rows for groups occurring in `query` are filled (<=20 of 27 in practice).
    let mut present = [false; 27];
    for &b in query {
        present[sub_index(b)] = true;
    }
    s.prof.clear();
    s.prof.resize(27 * n_max, zero);
    for (group, &p) in present.iter().enumerate() {
        if !p {
            continue;
        }
        let row = blosum_row_i16(group_letter(group));
        let base = group * n_max;
        for (jc, col) in s.col_groups.iter().enumerate() {
            s.prof[base + jc] = i16x8::new([
                row[col[0] as usize],
                row[col[1] as usize],
                row[col[2] as usize],
                row[col[3] as usize],
                row[col[4] as usize],
                row[col[5] as usize],
                row[col[6] as usize],
                row[col[7] as usize],
            ]);
        }
    }

    s.h_prev.clear();
    s.h_prev.resize(n_max + 1, zero);
    s.h_cur.clear();
    s.h_cur.resize(n_max + 1, zero);
    s.e.clear();
    s.e.resize(n_max + 1, neg);
    let mut best = zero;

    for i in 1..=m {
        let qgroup = sub_index(query[i - 1]);
        let prof_row = &s.prof[qgroup * n_max..qgroup * n_max + n_max];
        s.h_cur[0] = zero;
        let mut f = neg;
        for j in 1..=n_max {
            let e_open = s.h_prev[j] - gap_open;
            let e_ext = s.e[j] - gap_ext;
            let mut e_j = e_open.max(e_ext);

            let f_open = s.h_cur[j - 1] - gap_open;
            let f_ext = f - gap_ext;
            f = f_open.max(f_ext);

            let diag = s.h_prev[j - 1] + prof_row[j - 1];
            let mut v = diag.max(e_j).max(f).max(zero);

            // Padding columns: force H to 0, E to the sentinel for those lanes.
            let inv = s.invalid[j - 1];
            v = inv.blend(zero, v);
            e_j = inv.blend(neg, e_j);

            s.h_cur[j] = v;
            s.e[j] = e_j;
            best = best.max(v);
        }
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }

    best.to_array()
}

/// Representative letter byte for a substitution group index (0..=25 -> 'A'+g,
/// 26 -> a non-letter sentinel). Used to fetch the correct BLOSUM row when building
/// the query profile: `blosum_row_i16` keys on the query byte, so we feed back a
/// byte whose `sub_index` is exactly `group`.
#[inline]
fn group_letter(group: usize) -> u8 {
    if group < 26 {
        b'A' + group as u8
    } else {
        b'*' // non-letter: blosum_row_i16 yields the all -4 row, matching scalar
    }
}

/// Single-pair SIMD score (places `target` in lane 0; rest of the batch is unused).
/// Bit-exact with `local_score`. Convenience for tests / one-off calls.
pub fn local_score_simd(query: &[u8], target: &[u8]) -> i16 {
    local_score_batch(query, &[target])[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sw::local_score;

    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

    /// Tiny deterministic xorshift PRNG (no rand dependency).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn aa_seq(&mut self, len: usize) -> Vec<u8> {
            (0..len).map(|_| AA[(self.next() % 20) as usize]).collect()
        }
        fn range(&mut self, lo: usize, hi: usize) -> usize {
            lo + (self.next() as usize) % (hi - lo + 1)
        }
    }

    /// SIMD single-pair score must EXACTLY equal the scalar score, over thousands
    /// of random amino-acid pairs of varied lengths.
    #[test]
    fn simd_equals_scalar_random_pairs() {
        let mut rng = Rng(0x9E3779B97F4A7C15);
        for _ in 0..20_000 {
            let (qlen, tlen) = (rng.range(1, 140), rng.range(1, 140));
            let q = rng.aa_seq(qlen);
            let t = rng.aa_seq(tlen);
            let scalar = local_score(&q, &t);
            let simd = local_score_simd(&q, &t) as i32;
            assert_eq!(simd, scalar, "q={:?} t={:?}", String::from_utf8_lossy(&q), String::from_utf8_lossy(&t));
        }
    }

    /// Batched (8-lane) score must equal the scalar score for every lane, even
    /// when targets have *different* lengths (padding-mask correctness).
    #[test]
    fn simd_batch_equals_scalar_mixed_lengths() {
        let mut rng = Rng(0xD1B54A32D192ED03);
        for _ in 0..3_000 {
            let qlen = rng.range(1, 140);
            let q = rng.aa_seq(qlen);
            let nt = rng.range(1, LANES);
            let mut targets: Vec<Vec<u8>> = Vec::with_capacity(nt);
            #[allow(clippy::needless_range_loop)]
            for _ in 0..nt {
                let tlen = rng.range(1, 140);
                targets.push(rng.aa_seq(tlen));
            }
            let refs: Vec<&[u8]> = targets.iter().map(|t| t.as_slice()).collect();
            let got = local_score_batch(&q, &refs);
            for (l, t) in targets.iter().enumerate() {
                let scalar = local_score(&q, t);
                assert_eq!(got[l] as i32, scalar, "lane {l}: q={:?} t={:?}",
                    String::from_utf8_lossy(&q), String::from_utf8_lossy(t));
            }
        }
    }

    /// Edge cases: identical sequences, single residues, highly similar pairs.
    #[test]
    fn simd_equals_scalar_edge_cases() {
        let cases: &[(&str, &str)] = &[
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHKLMNPQRSTVWY"),
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFAAGHIKLMNPQRSTVWY"),
            ("MKLVWQ", "MKLVWQ"),
            ("QVQLVQSGAEVKKPGAS", "QVQLVQSGAEVKKPGAT"),
            ("WWWWWW", "WWWWWW"),
            ("A", "A"),
            ("A", "WWWWWWWWWW"),
            ("WWWWWWWWWW", "A"),
        ];
        for (q, t) in cases {
            assert_eq!(
                local_score_simd(q.as_bytes(), t.as_bytes()) as i32,
                local_score(q.as_bytes(), t.as_bytes()),
                "q={q:?} t={t:?}"
            );
        }
    }
}
