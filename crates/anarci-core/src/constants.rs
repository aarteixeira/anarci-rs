//! Numbering constants and small lookup helpers.

use once_cell::sync::Lazy;
use std::collections::HashMap;

// ALPHABET (53) and BLOSUM62_RAW (276) generated verbatim from reference schemes.py.
include!("generated_constants.rs");

/// Python `alphabet[i]`: `i == -1` (or any negative) -> `" "`; `i` in `0..=51` -> code.
#[inline]
pub fn alpha(i: i32) -> &'static str {
    if i < 0 {
        ALPHABET[52] // " "
    } else {
        ALPHABET[i as usize]
    }
}

/// `az = alphabet[:-1]` (52 codes A..ZZ); used by `get_imgt_cdr`.
#[inline]
pub fn az_imgt(i: usize) -> &'static str {
    ALPHABET[i]
}

/// `za = az[::-1]` (reverse of the 52-code az); used by `get_imgt_cdr`.
#[inline]
pub fn za_imgt(i: usize) -> &'static str {
    ALPHABET[51 - i]
}

/// CDR3-local 26-letter `az = "ABC...Z"` used by `get_cdr3_annotations`.
pub static AZ26: [&str; 26] = [
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
    "T", "U", "V", "W", "X", "Y", "Z",
];

/// CDR3-local 26-letter `za = "ZYX...A"` (reverse of AZ26).
pub static ZA26: [&str; 26] = [
    "Z", "Y", "X", "W", "V", "U", "T", "S", "R", "Q", "P", "O", "N", "M", "L", "K", "J", "I", "H",
    "G", "F", "E", "D", "C", "B", "A",
];

static BLOSUM62: Lazy<HashMap<(u8, u8), i32>> = Lazy::new(|| {
    let mut m = HashMap::with_capacity(BLOSUM62_RAW.len());
    for &(a, b, v) in BLOSUM62_RAW {
        m.insert((a, b), v);
    }
    m
});

/// BLOSUM62 score, trying `(a, b)` then `(b, a)` (the dict is stored asymmetrically).
/// Returns `None` if neither direction exists (Python would raise KeyError; the only
/// caller, `_get_wolfguy_L1`, catches that to try the reverse).
#[inline]
pub fn blosum62(a: u8, b: u8) -> Option<i32> {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    BLOSUM62
        .get(&(a, b))
        .or_else(|| BLOSUM62.get(&(b, a)))
        .copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alphabet_mapping() {
        assert_eq!(alpha(-1), " ");
        assert_eq!(alpha(0), "A");
        assert_eq!(alpha(25), "Z");
        assert_eq!(alpha(26), "AA");
        assert_eq!(alpha(51), "ZZ");
        assert_eq!(ALPHABET.len(), 53);
    }

    #[test]
    fn za_is_reverse_of_az() {
        assert_eq!(za_imgt(0), "ZZ");
        assert_eq!(za_imgt(51), "A");
        assert_eq!(AZ26[0], "A");
        assert_eq!(ZA26[0], "Z");
    }

    #[test]
    fn blosum_symmetry_via_fallback() {
        // ("W","L") is stored; ("L","W") is not -> fallback must find it.
        assert_eq!(blosum62(b'W', b'L'), Some(-2));
        assert_eq!(blosum62(b'L', b'W'), Some(-2));
        assert_eq!(blosum62(b'Y', b'Y'), Some(7));
    }
}
