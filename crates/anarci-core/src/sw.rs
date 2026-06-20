//! Smith-Waterman local alignment (Gotoh affine gaps) with BLOSUM62, plus the
//! Karlin-Altschul statistics RIOT uses to rank germline genes by e-value.
//!
//! This reproduces RIOT's amino-acid germline-assignment scoring
//! (`riot_na.alignment`, v5.0.3): a Striped-Smith-Waterman local alignment with
//! `gap_open=11`, `gap_extend=1`, the BLOSUM62 matrix, and the e-value
//!   E = m * n * 2^(-S)
//! where `S` is the raw alignment score, `m` the query length and `n` the gene
//! DB length. (RIOT passes the raw score, not the bit score, into its e-value
//! formula — see `skbio_alignment.align_aa`; we match that so rankings agree.)
//!
//! We compute only the optimal score (no traceback): germline ranking needs the
//! score, not the CIGAR. Score-only Gotoh is O(m·n) time, O(n) space.

/// AA Karlin-Altschul parameters for BLOSUM62 with gap-open 11 / extend 1, as
/// used by RIOT (`GumbellParams.AA`). lambda/K only affect the reported bit
/// score; gene *ranking* depends only on the raw score (see [`evalue`]).
/// Literals are copied verbatim from RIOT (more digits than f64 holds — the
/// extra digits round away harmlessly; we keep them for exact provenance).
#[allow(clippy::excessive_precision)]
pub const KA_LAMBDA_AA: f64 = 0.265_536_056_332_419_22;
#[allow(clippy::excessive_precision)]
pub const KA_K_AA: f64 = 0.043_186_874_595_437_463;

/// RIOT's SSW gap penalties for the amino-acid pipeline.
pub const GAP_OPEN: i32 = 11;
pub const GAP_EXTEND: i32 = 1;

/// BLOSUM62, indexed by `(byte - b'A')` for both residues (letters A..Z).
/// Generated from RIOT's `riot_na.data.constants.BLOSUM_62`. Letters with no
/// BLOSUM62 entry (J, O, U) get the matrix's stop/unknown penalty -4; queries
/// only contain the 20 standard residues (enforced by `validate_sequence`), so
/// those rows are never reached in practice.
#[rustfmt::skip]
pub(crate) const BLOSUM62: [[i32; 26]; 26] = [
    [4, -2, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -1, -2, -4, -1, -1, -1, 1, 0, -4, 0, -3, -1, -2, -1], // A
    [-2, 4, -3, 4, 1, -3, -1, 0, -3, -3, 0, -4, -3, 4, -4, -2, 0, -1, 0, -1, -4, -3, -4, -1, -3, 0], // B
    [0, -3, 9, -3, -4, -2, -3, -3, -1, -1, -3, -1, -1, -3, -4, -3, -3, -3, -1, -1, -4, -1, -2, -1, -2, -3], // C
    [-2, 4, -3, 6, 2, -3, -1, -1, -3, -3, -1, -4, -3, 1, -4, -1, 0, -2, 0, -1, -4, -3, -4, -1, -3, 1], // D
    [-1, 1, -4, 2, 5, -3, -2, 0, -3, -3, 1, -3, -2, 0, -4, -1, 2, 0, 0, -1, -4, -2, -3, -1, -2, 4], // E
    [-2, -3, -2, -3, -3, 6, -3, -1, 0, 0, -3, 0, 0, -3, -4, -4, -3, -3, -2, -2, -4, -1, 1, -1, 3, -3], // F
    [0, -1, -3, -1, -2, -3, 6, -2, -4, -4, -2, -4, -3, 0, -4, -2, -2, -2, 0, -2, -4, -3, -2, -1, -3, -2], // G
    [-2, 0, -3, -1, 0, -1, -2, 8, -3, -3, -1, -3, -2, 1, -4, -2, 0, 0, -1, -2, -4, -3, -2, -1, 2, 0], // H
    [-1, -3, -1, -3, -3, 0, -4, -3, 4, 3, -3, 2, 1, -3, -4, -3, -3, -3, -2, -1, -4, 3, -3, -1, -1, -3], // I
    [-1, -3, -1, -3, -3, 0, -4, -3, 3, 3, -3, 3, 2, -3, -4, -3, -2, -2, -2, -1, -4, 2, -2, -1, -1, -3], // J
    [-1, 0, -3, -1, 1, -3, -2, -1, -3, -3, 5, -2, -1, 0, -4, -1, 1, 2, 0, -1, -4, -2, -3, -1, -2, 1], // K
    [-1, -4, -1, -4, -3, 0, -4, -3, 2, 3, -2, 4, 2, -3, -4, -3, -2, -2, -2, -1, -4, 1, -2, -1, -1, -3], // L
    [-1, -3, -1, -3, -2, 0, -3, -2, 1, 2, -1, 2, 5, -2, -4, -2, 0, -1, -1, -1, -4, 1, -1, -1, -1, -1], // M
    [-2, 4, -3, 1, 0, -3, 0, 1, -3, -3, 0, -3, -2, 6, -4, -2, 0, 0, 1, 0, -4, -3, -4, -1, -2, 0], // N
    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4], // O
    [-1, -2, -3, -1, -1, -4, -2, -2, -3, -3, -1, -3, -2, -2, -4, 7, -1, -2, -1, -1, -4, -2, -4, -1, -3, -1], // P
    [-1, 0, -3, 0, 2, -3, -2, 0, -3, -2, 1, -2, 0, 0, -4, -1, 5, 1, 0, -1, -4, -2, -2, -1, -1, 4], // Q
    [-1, -1, -3, -2, 0, -3, -2, 0, -3, -2, 2, -2, -1, 0, -4, -2, 1, 5, -1, -1, -4, -3, -3, -1, -2, 0], // R
    [1, 0, -1, 0, 0, -2, 0, -1, -2, -2, 0, -2, -1, 1, -4, -1, 0, -1, 4, 1, -4, -2, -3, -1, -2, 0], // S
    [0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1, 0, -4, -1, -1, -1, 1, 5, -4, 0, -2, -1, -2, -1], // T
    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4], // U
    [0, -3, -1, -3, -2, -1, -3, -3, 3, 2, -2, 1, 1, -3, -4, -2, -2, -3, -2, 0, -4, 4, -3, -1, -1, -2], // V
    [-3, -4, -2, -4, -3, 1, -2, -2, -3, -2, -3, -2, -1, -4, -4, -4, -2, -3, -3, -2, -4, -3, 11, -1, 2, -2], // W
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -1], // X
    [-2, -3, -2, -3, -2, 3, -3, 2, -1, -1, -2, -1, -1, -2, -4, -3, -1, -2, -2, -2, -4, -1, 2, -1, 7, -2], // Y
    [-1, 0, -3, 1, 4, -3, -2, 0, -3, -3, 1, -3, -1, 0, -4, -1, 4, 0, 0, -1, -4, -2, -2, -1, -2, 4], // Z
];

#[inline]
fn score(a: u8, b: u8) -> i32 {
    let ai = (a.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
    let bi = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
    if ai < 26 && bi < 26 {
        BLOSUM62[ai][bi]
    } else {
        -4 // non-letter: treat as a strong mismatch (never hit on validated input)
    }
}

/// Optimal Smith-Waterman local-alignment **score** of `query` against `target`
/// with BLOSUM62 and affine gaps (open 11, extend 1, in RIOT's convention: a gap
/// of length L costs `open + (L-1)*extend`). Returns 0 if either is empty.
///
/// Gotoh recurrence, score only (no traceback), O(m·n) time / O(n) space.
pub fn local_score(query: &[u8], target: &[u8]) -> i32 {
    let (m, n) = (query.len(), target.len());
    if m == 0 || n == 0 {
        return 0;
    }
    // H: best local score ending at (i,j); E: gap in query (deletion in target row);
    // F: gap in target. Rolling rows over the target dimension (length n).
    let neg = i32::MIN / 2; // safe "negative infinity" without overflow on add
    let mut h_prev = vec![0i32; n + 1];
    let mut h_cur = vec![0i32; n + 1];
    let mut e = vec![neg; n + 1]; // E along the row, recomputed per i
    let mut best = 0i32;
    for i in 1..=m {
        let qi = query[i - 1];
        h_cur[0] = 0;
        let mut f = neg; // F resets at the start of each row
        for j in 1..=n {
            // E[j]: extend a gap in the query (consume target[j-1]); from H above or E above.
            let e_open = h_prev[j] - GAP_OPEN;
            let e_ext = e[j] - GAP_EXTEND;
            e[j] = e_open.max(e_ext);
            // F: extend a gap in the target (consume query[i-1]); from H left or F left.
            let f_open = h_cur[j - 1] - GAP_OPEN;
            let f_ext = f - GAP_EXTEND;
            f = f_open.max(f_ext);
            let diag = h_prev[j - 1] + score(qi, target[j - 1]);
            let v = diag.max(e[j]).max(f).max(0);
            h_cur[j] = v;
            if v > best {
                best = v;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }
    best
}

/// Sequence identity of the optimal local alignment of `query` vs `target`:
/// (identical aligned residue pairs) / (alignment columns, gaps included),
/// matching RIOT's `calculate_seq_identity` over the aligned span. Needs a
/// traceback, so this is heavier than [`local_score`]; call it only for the one
/// winning gene, not the whole DB. Returns `(identity, raw_score)`; identity is
/// `0.0` for empty inputs or a zero-length alignment.
pub fn local_identity(query: &[u8], target: &[u8]) -> (f64, i32) {
    let (m, n) = (query.len(), target.len());
    if m == 0 || n == 0 {
        return (0.0, 0);
    }
    let neg = i32::MIN / 2;
    // Full DP matrices for traceback. H/E/F as in Gotoh; ptr records the move
    // that produced H[i][j]: 0=stop(0), 1=diag, 2=up(gap in query/target row),
    // 3=left(gap in target/query col).
    let w = n + 1;
    let mut h = vec![0i32; (m + 1) * w];
    let mut e = vec![neg; (m + 1) * w];
    let mut f = vec![neg; (m + 1) * w];
    let mut ptr = vec![0u8; (m + 1) * w];
    let mut best = 0i32;
    let (mut bi, mut bj) = (0usize, 0usize);
    for i in 1..=m {
        let qi = query[i - 1];
        for j in 1..=n {
            let idx = i * w + j;
            e[idx] = (h[(i - 1) * w + j] - GAP_OPEN).max(e[(i - 1) * w + j] - GAP_EXTEND);
            f[idx] = (h[i * w + (j - 1)] - GAP_OPEN).max(f[i * w + (j - 1)] - GAP_EXTEND);
            let diag = h[(i - 1) * w + (j - 1)] + score(qi, target[j - 1]);
            let mut v = 0i32;
            let mut p = 0u8;
            if diag > v {
                v = diag;
                p = 1;
            }
            if e[idx] > v {
                v = e[idx];
                p = 2;
            }
            if f[idx] > v {
                v = f[idx];
                p = 3;
            }
            h[idx] = v;
            ptr[idx] = p;
            if v > best {
                best = v;
                bi = i;
                bj = j;
            }
        }
    }
    if best == 0 {
        return (0.0, 0);
    }
    // Traceback from the best cell to a stop (H==0), counting columns and matches.
    let (mut i, mut j) = (bi, bj);
    let (mut cols, mut matches) = (0usize, 0usize);
    while i > 0 && j > 0 {
        let idx = i * w + j;
        match ptr[idx] {
            1 => {
                cols += 1;
                if query[i - 1].eq_ignore_ascii_case(&target[j - 1]) {
                    matches += 1;
                }
                i -= 1;
                j -= 1;
            }
            2 => {
                cols += 1;
                i -= 1;
            }
            3 => {
                cols += 1;
                j -= 1;
            }
            _ => break, // stop (local start)
        }
    }
    if cols == 0 {
        (0.0, best)
    } else {
        (matches as f64 / cols as f64, best)
    }
}

/// E-value of a local alignment, RIOT-compatible: `E = m * n * 2^(-S)` where `S`
/// is the raw alignment score, `m` the query length, `n` the gene DB length.
/// Lower is a better hit. (This mirrors RIOT's `compute_evalue`, which — despite
/// the parameter name — exponentiates the *raw* score, so rankings match exactly.)
#[inline]
pub fn evalue(raw_score: i32, query_len: usize, db_len: usize) -> f64 {
    (query_len as f64) * (db_len as f64) * 2f64.powi(-raw_score)
}

/// Bit score `S' = (lambda*S - ln K) / ln 2` (BLOSUM62 AA params). Reported for
/// information; not used for ranking.
#[inline]
pub fn bit_score(raw_score: i32) -> f64 {
    (KA_LAMBDA_AA * raw_score as f64 - KA_K_AA.ln()) / std::f64::consts::LN_2
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Oracle values captured from RIOT's StripedSmithWaterman (scikit-bio) with
    /// AA_ALIGNER_PARAMS (gap_open=11, gap_extend=1, BLOSUM62). The score-only
    /// Gotoh here must reproduce the SSW optimal_alignment_score exactly.
    #[test]
    fn matches_ssw_oracle() {
        let cases: &[(&str, &str, i32)] = &[
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY", 116), // identical 20-mer
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHKLMNPQRSTVWY", 101),  // one deletion
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFAAGHIKLMNPQRSTVWY", 104), // two insertions
            ("MKLVWQ", "MKLVWQ", 34),
            ("QVQLVQSGAEVKKPGAS", "QVQLVQSGAEVKKPGAT", 78), // one mismatch at end
            ("WWWWWW", "WWWWWW", 66),
        ];
        for (q, t, want) in cases {
            let got = local_score(q.as_bytes(), t.as_bytes());
            assert_eq!(got, *want, "score q={q:?} t={t:?}");
        }
    }

    #[test]
    fn blosum_self_scores() {
        assert_eq!(score(b'Q', b'Q'), 5);
        assert_eq!(score(b'A', b'A'), 4);
        assert_eq!(score(b'W', b'W'), 11);
        assert_eq!(score(b'C', b'W'), -2);
        assert_eq!(score(b'w', b'c'), -2); // case-insensitive
    }

    #[test]
    fn empty_inputs() {
        assert_eq!(local_score(b"", b"ABC"), 0);
        assert_eq!(local_score(b"ABC", b""), 0);
        assert_eq!(local_identity(b"", b"ABC"), (0.0, 0));
    }

    /// Identity over the aligned span must match RIOT's `calculate_seq_identity`
    /// (matches / cigar-column-count), and the traceback score must equal
    /// `local_score`. Oracle from scikit-bio SSW + riot_na.
    #[test]
    fn identity_matches_riot_oracle() {
        let cases: &[(&str, &str, f64, i32)] = &[
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY", 1.0, 116),
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHKLMNPQRSTVWY", 0.95, 101),
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFAAGHIKLMNPQRSTVWY", 20.0 / 22.0, 104),
            ("QVQLVQSGAEVKKPGAS", "QVQLVQSGAEVKKPGAT", 16.0 / 17.0, 78),
            (
                "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYC",
                "DIQMTQSPSSLSASVGDKVTITCRASQGISNALAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQH",
                80.0 / 88.0,
                413,
            ),
        ];
        for (q, t, want_id, want_score) in cases {
            let (id, sc) = local_identity(q.as_bytes(), t.as_bytes());
            assert_eq!(sc, *want_score, "score q={q:?}");
            assert!((id - want_id).abs() < 1e-9, "identity q={q:?}: got {id} want {want_id}");
        }
    }

    #[test]
    fn lower_evalue_for_higher_score() {
        assert!(evalue(200, 100, 5000) < evalue(150, 100, 5000));
    }
}
