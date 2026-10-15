//! DNA → protein: 6-frame translation and an adaptive frame-ordering counter.
//!
//! For NGS, reads from one library share a reading frame, so we keep a running
//! tally of which frames actually yield antibody domains and try the winners
//! first, short-circuiting once a frame numbers. Frame ids: 0,1,2 = forward at
//! offset 0/1/2; 3,4,5 = reverse-complement at offset 0/1/2.

use std::sync::atomic::{AtomicU64, Ordering};

pub const N_FRAMES: usize = 6;

/// Standard genetic code; stop → `*`, any codon containing a non-ACGT base → `X`.
fn codon(c: &[u8]) -> u8 {
    fn idx(b: u8) -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'T' | b'U' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    // TCAG order → the classic 64-entry table.
    const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    match (idx(c[0]), idx(c[1]), idx(c[2])) {
        (Some(a), Some(b), Some(d)) => AA[a * 16 + b * 4 + d],
        _ => b'X',
    }
}

fn revcomp(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Translate one frame (0..6) of `dna` into protein.
pub fn translate_frame(dna: &[u8], frame: usize) -> Vec<u8> {
    // Out-of-range frames produce no translation; public APIs reject these explicitly.
    if frame >= N_FRAMES {
        return Vec::new();
    }
    let (seq, off) = if frame < 3 {
        (dna.to_vec(), frame)
    } else {
        (revcomp(dna), frame - 3)
    };
    seq.get(off..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(codon)
        .collect()
}

/// Adaptive frame order shared across a batch: frames tried most-successful-first.
#[derive(Default)]
pub struct FrameOrder {
    hits: [AtomicU64; N_FRAMES],
}

impl FrameOrder {
    /// The 6 frame ids, ordered by past success (descending); ties keep 0..6 order,
    /// so a cold start tries the forward frames first.
    pub fn order(&self) -> [usize; N_FRAMES] {
        let mut idx = [0, 1, 2, 3, 4, 5];
        let h: [u64; N_FRAMES] = std::array::from_fn(|i| self.hits[i].load(Ordering::Relaxed));
        idx.sort_by(|&a, &b| h[b].cmp(&h[a]).then(a.cmp(&b)));
        idx
    }

    /// Record that `frame` produced a domain (updates the ordering for later reads).
    pub fn record(&self, frame: usize) {
        self.hits[frame].fetch_add(1, Ordering::Relaxed);
    }

    pub fn counts(&self) -> [u64; N_FRAMES] {
        std::array::from_fn(|i| self.hits[i].load(Ordering::Relaxed))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_known_codons() {
        // ATG GCT TAA -> M A *
        assert_eq!(translate_frame(b"ATGGCTTAA", 0), b"MA*");
        // frame 1 of the same drops the first base: TGG CTT AA -> W L
        assert_eq!(translate_frame(b"ATGGCTTAA", 1), b"WL");
        // reverse-complement frame 3: revcomp(ATGGCTTAA)=TTAAGCCAT -> L S H
        assert_eq!(translate_frame(b"ATGGCTTAA", 3), b"LSH");
        // ambiguity: a codon with N -> X
        assert_eq!(translate_frame(b"ATN", 0), b"X");
    }

    #[test]
    fn short_sequences_and_frame_bounds() {
        for seq in [b"".as_slice(), b"A", b"AT"] {
            for frame in 0..8 {
                assert!(translate_frame(seq, frame).is_empty());
            }
        }
    }

    #[test]
    fn frame_order_prefers_winners() {
        let fo = FrameOrder::default();
        assert_eq!(fo.order(), [0, 1, 2, 3, 4, 5]); // cold start
        fo.record(4);
        fo.record(4);
        fo.record(1);
        assert_eq!(fo.order()[0], 4); // most hits first
        assert_eq!(fo.order()[1], 1);
    }
}
