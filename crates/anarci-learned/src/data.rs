//! Training-data pipeline (P4): turn oracle domains + augmentation into a split,
//! leakage-controlled dataset.
//!
//! Antibody repertoires are highly redundant, so a random split lies — and so does a naive
//! one. The P4 red-team measured two real leaks in an earlier version and they are fixed
//! here:
//!  * **Split by source record, not per domain.** A multi-domain chain (scFv = VH+VL in one
//!    sequence) must go to ONE split, or the same physical molecule lands in both train and
//!    test. Domains are grouped into records by their full sequence.
//!  * **Length-tolerant similarity.** CDR3 (and full-sequence) near-duplicates differing by
//!    an indel must still cluster together; equal-length-only identity leaked 11% of test
//!    CDR3s. We use a gap-free best-overlap identity normalized by the shorter string.
//!
//! Records are clustered by **full-sequence OR CDR3** near-identity (union-find), whole
//! clusters go to train/val/test, and every synthetic derivative inherits its source
//! domain's split (composed multi-V examples pair domains within one split only). The
//! leakage invariant — no near-duplicate across splits — is asserted in `tests/data.rs`.

use crate::augment::{concat_domains, LabeledSeq};
use anarci_core::{State, StateType};
use std::collections::HashMap;

/// A source V domain: its residues (the full source sequence) and oracle state vector.
/// Domains sharing a `seq` belong to the same source record (e.g. the two halves of an scFv).
#[derive(Clone)]
pub struct DomainRec {
    pub seq: Vec<u8>,
    pub sv: Vec<State>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Split {
    Train,
    Val,
    Test,
}

/// One training example, tagged with its split and the source-domain clusters it derives
/// from (empty for junk). `kind` labels the augmentation type.
pub struct Example {
    pub data: LabeledSeq,
    pub split: Split,
    pub clusters: Vec<usize>,
    pub kind: &'static str,
}

/// Clustering thresholds (identity, length-tolerant). Full-sequence near-dups and CDR3
/// near-dups both trigger a cluster merge.
pub const FULL_THR: f32 = 0.95;
pub const CDR3_THR: f32 = 0.80;

/// CDR3 residues (IMGT 105–117) of a domain, in order.
pub fn cdr3_of(rec: &DomainRec) -> Vec<u8> {
    rec.sv
        .iter()
        .filter(|s| s.typ != StateType::D && (105..=117).contains(&(s.id as i32)))
        .filter_map(|s| s.si.map(|i| rec.seq[i]))
        .collect()
}

/// Length-tolerant identity: best gap-free overlap of the shorter string against the longer,
/// normalized by the shorter length. Equal-length inputs reduce to plain identity; an indel
/// near-duplicate (a shorter string that is ~a substring of a longer one) still scores high.
pub fn overlap_identity(a: &[u8], b: &[u8]) -> f32 {
    let (s, l) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if s.is_empty() {
        return 0.0;
    }
    let mut best = 0usize;
    for off in 0..=(l.len() - s.len()) {
        let m = s.iter().zip(&l[off..]).filter(|(x, y)| x == y).count();
        best = best.max(m);
    }
    best as f32 / s.len() as f32
}

/// Whether two full sequences are near-identical (length-gated so we don't scan wildly
/// different lengths — near-dups have similar length).
fn full_similar(a: &[u8], b: &[u8], thr: f32) -> bool {
    let (lo, hi) = (a.len().min(b.len()), a.len().max(b.len()));
    if hi == 0 || (lo as f32 / hi as f32) < thr {
        return false;
    }
    overlap_identity(a, b) >= thr
}

struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect() }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != c {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// Cluster domains leakage-safely: group into records by full sequence, union records that
/// share a near-identical full sequence OR any near-identical CDR3, and return the cluster
/// id (union-find root record) per domain.
pub fn cluster_records(domains: &[DomainRec], full_thr: f32, cdr3_thr: f32) -> Vec<usize> {
    // records: distinct full sequences
    let mut seq_to_rec: HashMap<&[u8], usize> = HashMap::new();
    let mut rec_of_domain = vec![0usize; domains.len()];
    let mut rec_seq: Vec<&[u8]> = Vec::new();
    let mut rec_cdr3: Vec<Vec<Vec<u8>>> = Vec::new();
    for (i, d) in domains.iter().enumerate() {
        let rid = *seq_to_rec.entry(d.seq.as_slice()).or_insert_with(|| {
            rec_seq.push(d.seq.as_slice());
            rec_cdr3.push(Vec::new());
            rec_seq.len() - 1
        });
        rec_of_domain[i] = rid;
        rec_cdr3[rid].push(cdr3_of(d));
    }

    let n = rec_seq.len();
    let mut uf = UnionFind::new(n);
    for i in 0..n {
        for j in (i + 1)..n {
            let merge = full_similar(rec_seq[i], rec_seq[j], full_thr)
                || rec_cdr3[i].iter().any(|ca| {
                    !ca.is_empty() && rec_cdr3[j].iter().any(|cb| overlap_identity(ca, cb) >= cdr3_thr)
                });
            if merge {
                uf.union(i, j);
            }
        }
    }
    (0..domains.len()).map(|i| uf.find(rec_of_domain[i])).collect()
}

/// Deterministic 70/15/15 split from a cluster's representative bytes (FNV-1a; stable across
/// runs — no RNG, so a cluster always lands in the same split).
pub fn split_for_rep(rep: &[u8]) -> Split {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in rep {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    match h % 100 {
        0..=69 => Split::Train,
        70..=84 => Split::Val,
        _ => Split::Test,
    }
}

/// Per-domain split via record clustering (every domain of a record shares one split).
pub fn splits_of(domains: &[DomainRec]) -> Vec<Split> {
    let cluster = cluster_records(domains, FULL_THR, CDR3_THR);
    let mut rep: HashMap<usize, &[u8]> = HashMap::new();
    for (i, d) in domains.iter().enumerate() {
        let e = rep.entry(cluster[i]).or_insert(d.seq.as_slice());
        if d.seq.as_slice() < *e {
            *e = d.seq.as_slice();
        }
    }
    (0..domains.len()).map(|i| split_for_rep(rep[&cluster[i]])).collect()
}

/// Whether two source domains are near-duplicates under the pipeline's thresholds (full
/// sequence or CDR3). Clustering unions every such pair, so no near-duplicate crosses splits.
pub fn near_duplicate(a: &DomainRec, b: &DomainRec) -> bool {
    if full_similar(&a.seq, &b.seq, FULL_THR) {
        return true;
    }
    let (ca, cb) = (cdr3_of(a), cdr3_of(b));
    !ca.is_empty() && !cb.is_empty() && overlap_identity(&ca, &cb) >= CDR3_THR
}

/// Build the dataset: leakage-safe cluster + split, then single / partial / multi-V / junk
/// examples, every derivative inheriting its source domain's split (composed examples pair
/// domains within one split only). `partial_windows` = IMGT (lo,hi) ranges emitted as partials.
pub fn build_dataset(domains: &[DomainRec], partial_windows: &[(u8, u8)]) -> Vec<Example> {
    let cluster = cluster_records(domains, FULL_THR, CDR3_THR);
    // stable representative per cluster: the lexicographically-smallest source sequence in it.
    let mut rep: HashMap<usize, &[u8]> = HashMap::new();
    for (i, d) in domains.iter().enumerate() {
        let e = rep.entry(cluster[i]).or_insert(d.seq.as_slice());
        if d.seq.as_slice() < *e {
            *e = d.seq.as_slice();
        }
    }
    let split_of = |c: usize| split_for_rep(rep[&c]);

    let mut out = Vec::new();
    // per-domain: single + partials
    for (i, d) in domains.iter().enumerate() {
        let c = cluster[i];
        let split = split_of(c);
        let seg = LabeledSeq::from_domain(&d.seq, &d.sv);
        out.push(Example { data: seg.clone(), split, clusters: vec![c], kind: "single" });
        for &(lo, hi) in partial_windows {
            let p = seg.restrict_positions(lo, hi);
            if !p.domains().is_empty() {
                out.push(Example { data: p, split, clusters: vec![c], kind: "partial" });
            }
        }
    }
    // multi-V: pair domains WITHIN the same split so a composed example never crosses splits.
    for target in [Split::Train, Split::Val, Split::Test] {
        let idxs: Vec<usize> = (0..domains.len()).filter(|&i| split_of(cluster[i]) == target).collect();
        for pair in idxs.chunks(2) {
            if let [a, b] = *pair {
                let sa = LabeledSeq::from_domain(&domains[a].seq, &domains[a].sv);
                let sb = LabeledSeq::from_domain(&domains[b].seq, &domains[b].sv);
                out.push(Example {
                    data: concat_domains(&[&sa, &sb]),
                    split: target,
                    clusters: vec![cluster[a], cluster[b]],
                    kind: "multi_v",
                });
            }
        }
    }
    // junk negatives (no source cluster -> cannot leak). Distribute round-robin across ALL
    // three splits so the held-out set can measure non-antibody rejection (E2).
    let splits = [Split::Train, Split::Val, Split::Test];
    for (k, junk_seq) in JUNK.iter().enumerate() {
        out.push(Example {
            data: LabeledSeq::junk(junk_seq.as_bytes()),
            split: splits[k % 3],
            clusters: vec![],
            kind: "junk",
        });
    }
    out
}

/// Non-antibody negative sequences (all `Out`): lysozyme + unrelated proteins. Extend from a
/// real non-antibody set in P6; kept small here (distributed one-per-split by `build_dataset`).
const JUNK: &[&str] = &[
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHM",
    "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQSRRPSSWRQEKI",
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGF",
    "MASNDYTQQATQSYGAYPTQPGQGYSQQSSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYG",
];
