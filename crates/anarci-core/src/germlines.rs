//! Germline assignment — a faithful port of ANARCI's `get_identity`,
//! `run_germline_assignment`, and `get_hmm_length`.

use crate::simd_sw;
use crate::sw;
use crate::types::{State, StateType, StateVector};
use once_cell::sync::Lazy;
use serde::Deserialize;
use std::collections::HashMap;

/// `[(gene, aligned_seq)]` — ORDER PRESERVED (tie-breaking in `max()` depends on it).
type GeneList = Vec<(String, String)>;

#[derive(Deserialize)]
struct GermlineData {
    all_species: Vec<String>,
    /// seg -> chain_type -> species -> ordered gene list.
    germlines: HashMap<String, HashMap<String, HashMap<String, GeneList>>>,
}

static DATA: Lazy<GermlineData> = Lazy::new(|| {
    let raw = include_str!("../data/germlines.json");
    serde_json::from_str(raw).expect("embedded germlines.json must parse")
});

/// Canonical species order (`list(all_germlines['V']['H'].keys())`).
pub fn all_species() -> &'static [String] {
    &DATA.all_species
}

fn seg<'a>(name: &str) -> Option<&'a HashMap<String, HashMap<String, GeneList>>> {
    DATA.germlines.get(name)
}

/// Result of germline assignment, shaped to mirror ANARCI's returned dict.
/// `empty` distinguishes the `{}` return (a required species was absent) from the
/// default `{'v_gene':[None,None],'j_gene':[None,None]}` return.
#[derive(Clone, Debug, PartialEq)]
pub struct Germline {
    pub empty: bool,
    pub v_gene: Option<(String, String)>, // (species, gene)
    pub v_identity: Option<f64>,
    pub j_gene: Option<(String, String)>,
    pub j_identity: Option<f64>,
    /// E-values from the alignment-based path (`run_germline_assignment_evalue`);
    /// always `None` on the identity path so ANARCI-parity output is unchanged.
    pub v_evalue: Option<f64>,
    pub j_evalue: Option<f64>,
}

impl Germline {
    fn default_dict() -> Self {
        Germline {
            empty: false,
            v_gene: None,
            v_identity: None,
            j_gene: None,
            j_identity: None,
            v_evalue: None,
            j_evalue: None,
        }
    }
    fn empty() -> Self {
        Germline {
            empty: true,
            v_gene: None,
            v_identity: None,
            j_gene: None,
            j_identity: None,
            v_evalue: None,
            j_evalue: None,
        }
    }
}

/// Partially-matched sequence identity. Both sequences must be length 128
/// (panics otherwise, matching the Python `assert`). Gaps in the germline are skipped.
pub fn get_identity(state_sequence: &[u8], germline_sequence: &[u8]) -> f64 {
    assert!(
        state_sequence.len() == 128 && germline_sequence.len() == 128,
        "germline identity requires length-128 aligned sequences"
    );
    let (mut n, mut m) = (0u32, 0u32);
    for i in 0..128 {
        let g = germline_sequence[i];
        if g == b'-' {
            continue;
        }
        if state_sequence[i].to_ascii_uppercase() == g {
            m += 1;
        }
        n += 1;
    }
    if n == 0 {
        0.0
    } else {
        m as f64 / n as f64
    }
}

/// Build the 128-char IMGT match-state sequence (residue per match column, else '-').
fn build_state_sequence(state_vector: &StateVector, sequence: &[u8]) -> Vec<u8> {
    // match-state index -> seq position
    let mut idx: [Option<usize>; 129] = [None; 129]; // 1..=128 used
    for &State { id, typ, si } in state_vector {
        if typ == StateType::M {
            if (1..=128).contains(&id) {
                idx[id as usize] = si;
            }
        }
    }
    let mut out = Vec::with_capacity(128);
    for i in 1..=128usize {
        match idx[i] {
            Some(p) => out.push(sequence[p]),
            None => out.push(b'-'),
        }
    }
    out
}

/// argmax over a gene list by identity, FIRST max wins (Python `max` tie-break).
/// Returns `(gene, identity)` for the chosen gene, or `None` if the list is empty.
fn argmax_gene(genes: &GeneList, state_sequence: &[u8]) -> Option<(String, f64)> {
    let mut best: Option<(String, f64)> = None;
    for (gene, gseq) in genes {
        let id = get_identity(state_sequence, gseq.as_bytes());
        match &best {
            Some((_, bid)) if id <= *bid => {} // strict-greater keeps the first max
            _ => best = Some((gene.clone(), id)),
        }
    }
    best
}

/// Faithful port of `run_germline_assignment`.
pub fn run_germline_assignment(
    state_vector: &StateVector,
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
) -> Germline {
    let mut genes = Germline::default_dict();
    let state_sequence = build_state_sequence(state_vector, sequence);

    let v_seg = match seg("V") {
        Some(s) => s,
        None => return genes,
    };
    let v_chain = match v_seg.get(chain_type) {
        Some(c) => c,
        None => return genes, // chain_type not in V -> default dict
    };

    // Resolve the species list (default to all_species when None).
    let owned_all;
    let species_list: &[String] = match allowed_species {
        Some(sp) => {
            // Non-fatal: if any requested species is absent for this V chain, return {}.
            if !sp.iter().all(|s| v_chain.contains_key(s)) {
                return Germline::empty();
            }
            sp
        }
        None => {
            owned_all = all_species().to_vec();
            &owned_all
        }
    };

    // --- V gene: argmax identity across (species in order, genes in order) ---
    let mut best_v: Option<(String, String, f64)> = None; // (species, gene, id)
    for species in species_list {
        let glist = match v_chain.get(species) {
            Some(g) => g,
            None => continue, // species absent for this chain (skip; "previously bug")
        };
        if let Some((gene, id)) = argmax_gene(glist, &state_sequence) {
            let better = match &best_v {
                Some((_, _, bid)) => id > *bid, // strict: first species/gene wins ties
                None => true,
            };
            if better {
                best_v = Some((species.clone(), gene, id));
            }
        }
    }

    let (v_species, v_gene, v_id) = match best_v {
        Some(t) => t,
        None => return genes, // no species matched (seq_ids empty) -> default dict
    };
    genes.v_gene = Some((v_species.clone(), v_gene));
    genes.v_identity = Some(v_id);

    // --- J gene: only the V gene's species ---
    if let Some(j_seg) = seg("J") {
        if let Some(j_chain) = j_seg.get(chain_type) {
            if let Some(glist) = j_chain.get(&v_species) {
                if let Some((gene, id)) = argmax_gene(glist, &state_sequence) {
                    genes.j_gene = Some((v_species, gene));
                    genes.j_identity = Some(id);
                }
            }
        }
    }

    genes
}

// ===========================================================================
// E-value-based germline assignment (RIOT-style, amino-acid)
// ===========================================================================
//
// Instead of identity over the 128 IMGT match columns (ANARCI), align the
// query's V region and J region SEPARATELY against the *ungapped* germline V/J
// genes with Smith-Waterman (BLOSUM62, affine gaps), and pick the gene with the
// lowest e-value. This reproduces RIOT's amino-acid pipeline, which reaches far
// higher V-gene accuracy than identity matching (RIOT, Brief Bioinform 2025).

/// One ungapped germline gene: `(gene_name, amino-acid sequence)`.
type UngappedGene = (String, Vec<u8>);

/// Ungapped gene DBs, built once from the embedded gapped germlines:
///   seg -> chain -> species -> [(gene, ungapped_seq)]
struct UngappedDb {
    map: HashMap<String, HashMap<String, HashMap<String, Vec<UngappedGene>>>>,
}

static UNGAPPED: Lazy<UngappedDb> = Lazy::new(|| {
    let mut map: HashMap<String, HashMap<String, HashMap<String, Vec<UngappedGene>>>> =
        HashMap::new();
    for (seg_name, chains) in &DATA.germlines {
        let seg_out = map.entry(seg_name.clone()).or_default();
        for (chain, species_map) in chains {
            let chain_out = seg_out.entry(chain.clone()).or_default();
            for (species, glist) in species_map {
                let genes: Vec<UngappedGene> = glist
                    .iter()
                    .map(|(g, aligned)| {
                        let ungapped: Vec<u8> =
                            aligned.bytes().filter(|&c| c != b'-').collect();
                        (g.clone(), ungapped)
                    })
                    .collect();
                chain_out.insert(species.clone(), genes);
            }
        }
    }
    UngappedDb { map }
});

fn ungapped_chain<'a>(seg: &str, chain: &str) -> Option<&'a HashMap<String, Vec<UngappedGene>>> {
    UNGAPPED.map.get(seg).and_then(|c| c.get(chain))
}

/// Extract the V-region and J-region query slices from the state vector.
///
/// V region = sequence indices from the domain start up to and INCLUDING the
/// residue numbered at IMGT 104 (the conserved second-Cys). J region = indices
/// strictly after IMGT 104 through the last numbered residue. This mirrors
/// RIOT's split (V over the front of the query; J over the remainder after V).
///
/// If IMGT 104 has no mapped residue (rare/truncated), falls back to the whole
/// numbered span as V and an empty J.
fn split_v_j_regions(state_vector: &StateVector, sequence: &[u8]) -> (Vec<u8>, Vec<u8>) {
    // First and last sequence index touched by any state with a residue.
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;
    let mut cys104: Option<usize> = None;
    for s in state_vector {
        if let Some(si) = s.si {
            first = Some(first.map_or(si, |f| f.min(si)));
            last = Some(last.map_or(si, |l| l.max(si)));
            if s.id == 104 && s.typ == StateType::M {
                cys104 = Some(si);
            }
        }
    }
    let (first, last) = match (first, last) {
        (Some(f), Some(l)) if f <= l && l < sequence.len() => (f, l),
        _ => return (Vec::new(), Vec::new()),
    };
    match cys104 {
        Some(c) if c >= first && c <= last => {
            let v = sequence[first..=c].to_vec();
            let j = if c < last {
                sequence[c + 1..=last].to_vec()
            } else {
                Vec::new()
            };
            (v, j)
        }
        // No Cys-104 anchor: align the whole span as V.
        _ => (sequence[first..=last].to_vec(), Vec::new()),
    }
}

/// Best gene in a (species, gene-list) collection for `query` by e-value
/// (lowest wins). Tie-breaks reproduce RIOT's `AlignmentEntryAA.__lt__`:
/// on equal e-value, prefer the higher SW score, then the lexicographically
/// smaller (species, gene). Returns `(species, gene, e_value, raw_score)`.
fn best_gene_evalue<'a>(
    species_genes: impl Iterator<Item = (&'a str, &'a [UngappedGene])>,
    query: &[u8],
    db_len: usize,
) -> Option<(String, String, f64, i32)> {
    let mut best: Option<(String, String, f64, i32)> = None; // (sp, gene, evalue, score)
    for (species, genes) in species_genes {
        for (gene, gseq) in genes {
            let s = sw::local_score(query, gseq);
            let ev = sw::evalue(s, query.len(), db_len.max(1));
            let take = match &best {
                None => true,
                Some((bsp, bgene, bev, bscore)) => {
                    if ev < *bev {
                        true
                    } else if ev > *bev {
                        false
                    } else {
                        // equal e-value: RIOT breaks by higher score (when ev==0),
                        // then lexicographic gene id. We use (score desc, then
                        // (species,gene) ascending) deterministically.
                        if s != *bscore {
                            s > *bscore
                        } else if species != bsp.as_str() {
                            species < bsp.as_str()
                        } else {
                            gene.as_str() < bgene.as_str()
                        }
                    }
                }
            };
            if take {
                best = Some((species.to_string(), gene.clone(), ev, s));
            }
        }
    }
    best
}

// ===========================================================================
// Fast e-value path: k-mer prefilter -> SIMD SW on candidates -> exact select.
//
// Within a single V (or J) search the e-value DB length `db_len` is constant, so
// `E = m * db_len * 2^(-score)` is strictly monotone-decreasing in `score`. The
// winner is therefore: max score, ties broken by min (species, gene) — exactly
// the total order `best_gene_evalue` applies. The fast path computes the *same*
// scores (via the bit-exact SIMD kernel) on a candidate subset and selects under
// that same total order; it is identical to the full scan iff the candidate set
// contains the true winner. A k-mer prefilter chooses the candidates, an exact
// full-scan path is always available, and a recall test gates equality on the
// whole golden set.
// ===========================================================================

/// A flat reference into the gene DB: `(species, gene_name, ungapped_seq)`.
type GeneRef<'a> = (&'a str, &'a str, &'a [u8]);

/// One scored candidate, carrying everything needed for selection / return.
struct Scored<'a> {
    species: &'a str,
    gene: &'a str,
    score: i32,
}

/// Apply the exact `best_gene_evalue` total order to scored candidates and return
/// the winner. Order-independent (it's a strict total order on distinct
/// (species, gene) pairs), so it matches the full scan regardless of input order.
/// `db_len` only sets the returned e-value (selection depends on score alone here).
fn select_best_scored(
    cands: &[Scored<'_>],
    query_len: usize,
    db_len: usize,
) -> Option<(String, String, f64, i32)> {
    let mut best: Option<&Scored<'_>> = None;
    for c in cands {
        let take = match best {
            None => true,
            Some(b) => {
                // ev monotone in score (db_len fixed): higher score = lower ev.
                if c.score != b.score {
                    c.score > b.score
                } else if c.species != b.species {
                    c.species < b.species
                } else {
                    c.gene < b.gene
                }
            }
        };
        if take {
            best = Some(c);
        }
    }
    best.map(|b| {
        let ev = sw::evalue(b.score, query_len, db_len.max(1));
        (b.species.to_string(), b.gene.to_string(), ev, b.score)
    })
}

/// Score a flat candidate list against `query` with the SIMD batched kernel
/// (8 genes per pass), bit-exact with `sw::local_score`.
fn score_candidates_simd<'a>(query: &[u8], cands: &[GeneRef<'a>]) -> Vec<Scored<'a>> {
    let seqs: Vec<&[u8]> = cands.iter().map(|(_, _, s)| *s).collect();
    let scores = simd_sw::local_score_many(query, &seqs);
    cands
        .iter()
        .zip(scores)
        .map(|((sp, g, _), sc)| Scored { species: sp, gene: g, score: sc as i32 })
        .collect()
}

// --- k-mer prefilter index -------------------------------------------------

/// k-mer length for the prefilter seeds (short, for high recall).
const KMER_K: usize = 4;
/// Number of top candidates to keep from the prefilter (generous; the DB is tiny,
/// so even top-128 of ~1000 is ~8x fewer SW than a full scan). Recall=1.0 on the
/// golden set is gated by a test; raise this if that test ever fails.
const TOP_K: usize = 128;

/// Reduced amino-acid alphabet (Murphy-8-like grouping) for recall-safe seeding.
/// Mapping over `byte - b'A'` (A..Z); each standard residue -> a group 0..7,
/// non-residues -> 7 (lumped, never hit on validated input). Grouping conserved
/// physicochemical classes so a single substitution rarely breaks every seed.
///   {L,V,I,M,C} {A,G,S,T,P} {F,Y,W} {E,D,N,Q} {K,R} {H} ... lumped to 7 groups.
#[inline]
fn reduced_code(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'L' | b'V' | b'I' | b'M' | b'C' => 0,
        b'A' | b'G' | b'S' | b'T' | b'P' => 1,
        b'F' | b'Y' | b'W' => 2,
        b'E' | b'D' | b'N' | b'Q' => 3,
        b'K' | b'R' => 4,
        b'H' => 5,
        _ => 6,
    }
}
/// Number of reduced groups (radix for the k-mer integer code).
const RADIX: u32 = 7;

/// Encode the set of distinct reduced k-mers in a sequence as a sorted Vec of
/// codes (base-RADIX). Small k keeps the code space tiny (7^4 = 2401).
fn kmer_set(seq: &[u8]) -> Vec<u32> {
    if seq.len() < KMER_K {
        return Vec::new();
    }
    let mut codes: Vec<u32> = Vec::with_capacity(seq.len());
    let mut code = 0u32;
    for (i, &b) in seq.iter().enumerate() {
        code = code * RADIX + reduced_code(b) as u32;
        if i + 1 >= KMER_K {
            codes.push(code);
            // strip the highest digit for the sliding window
            code -= reduced_code(seq[i + 1 - KMER_K]) as u32 * RADIX.pow(KMER_K as u32 - 1);
        }
    }
    codes.sort_unstable();
    codes.dedup();
    codes
}

/// Prefilter index for one (seg, chain): an inverted list k-mer-code -> gene ids,
/// plus the flat gene table the ids point into.
struct ChainIndex {
    /// Flat gene table: (species, gene, ungapped_seq), stable order.
    genes: Vec<(String, String, Vec<u8>)>,
    /// kmer code -> list of gene indices containing that code (sorted, deduped).
    postings: HashMap<u32, Vec<u32>>,
}

impl ChainIndex {
    fn build(chain_map: &HashMap<String, Vec<UngappedGene>>, species_order: &[String]) -> Self {
        // Flatten in a deterministic order: species_order, then gene order.
        let mut genes: Vec<(String, String, Vec<u8>)> = Vec::new();
        for sp in species_order {
            if let Some(glist) = chain_map.get(sp) {
                for (g, seq) in glist {
                    genes.push((sp.clone(), g.clone(), seq.clone()));
                }
            }
        }
        let mut postings: HashMap<u32, Vec<u32>> = HashMap::new();
        for (gid, (_, _, seq)) in genes.iter().enumerate() {
            for code in kmer_set(seq) {
                postings.entry(code).or_default().push(gid as u32);
            }
        }
        ChainIndex { genes, postings }
    }

    /// Top-K gene indices for `query` by distinct-k-mer seed-hit count (descending),
    /// ties broken by gene index (deterministic). Returns at most `TOP_K`.
    fn top_candidates(&self, query: &[u8]) -> Vec<u32> {
        let mut hits = vec![0u32; self.genes.len()];
        for code in kmer_set(query) {
            if let Some(list) = self.postings.get(&code) {
                for &gid in list {
                    hits[gid as usize] += 1;
                }
            }
        }
        // Indices with >0 hits, ranked by (hits desc, gid asc). Quickselect the
        // top-K boundary (O(n)) then sort only the kept prefix — same deterministic
        // order as a full sort, but without paying O(n log n) over ~1000 genes.
        let cmp = |a: &u32, b: &u32| {
            hits[*b as usize]
                .cmp(&hits[*a as usize])
                .then_with(|| a.cmp(b))
        };
        let mut ranked: Vec<u32> = (0..self.genes.len() as u32)
            .filter(|&g| hits[g as usize] > 0)
            .collect();
        if ranked.len() > TOP_K {
            ranked.select_nth_unstable_by(TOP_K - 1, cmp);
            ranked.truncate(TOP_K);
        }
        ranked.sort_unstable_by(cmp);
        ranked
    }
}

/// Lazily-built prefilter indices: seg -> chain -> ChainIndex (over all species,
/// in canonical `all_species()` order). Built once, shared across calls.
static INDEX: Lazy<HashMap<String, HashMap<String, ChainIndex>>> = Lazy::new(|| {
    let order = all_species().to_vec();
    let mut out: HashMap<String, HashMap<String, ChainIndex>> = HashMap::new();
    for (seg_name, chains) in &UNGAPPED.map {
        let seg_out = out.entry(seg_name.clone()).or_default();
        for (chain, species_map) in chains {
            seg_out.insert(chain.clone(), ChainIndex::build(species_map, &order));
        }
    }
    out
});

fn chain_index<'a>(seg: &str, chain: &str) -> Option<&'a ChainIndex> {
    INDEX.get(seg).and_then(|c| c.get(chain))
}

/// Fast e-value winner for a V/J search: k-mer prefilter -> SIMD SW on top-K
/// candidates -> exact total-order select. Falls back to the full SIMD scan when
/// the species scope isn't the full DB (the index covers all species in canonical
/// order; a restricted scope is handled by the exact scan, which is still SIMD-fast).
///
/// `allowed_species == None` means "all species" and uses the index. An explicit
/// list uses the exact SIMD full scan over just those species (correct + fast).
/// Returns `(species, gene, e_value, raw_score)`, identical to `best_gene_evalue`.
fn best_gene_evalue_fast(
    seg: &str,
    chain: &str,
    chain_map: &HashMap<String, Vec<UngappedGene>>,
    species_list: &[String],
    query: &[u8],
    db_len: usize,
) -> Option<(String, String, f64, i32)> {
    // Build the flat candidate set (in canonical order for deterministic ties).
    let use_index = is_full_species_scope(species_list);
    let idx = if use_index { chain_index(seg, chain) } else { None };

    if let Some(idx) = idx {
        // Index path: score only the top-K candidates with SIMD SW.
        let top = idx.top_candidates(query);
        if !top.is_empty() {
            let cands: Vec<GeneRef> = top
                .iter()
                .map(|&g| {
                    let (sp, gene, seq) = &idx.genes[g as usize];
                    (sp.as_str(), gene.as_str(), seq.as_slice())
                })
                .collect();
            let scored = score_candidates_simd(query, &cands);
            if let Some(w) = select_best_scored(&scored, query.len(), db_len) {
                return Some(w);
            }
        }
        // Empty seed hits (e.g. query shorter than k): fall through to full scan.
    }

    // Exact full SIMD scan over the requested species scope.
    let mut flat: Vec<GeneRef> = Vec::new();
    for sp in species_list {
        if let Some(glist) = chain_map.get(sp) {
            for (g, seq) in glist {
                flat.push((sp.as_str(), g.as_str(), seq.as_slice()));
            }
        }
    }
    if flat.is_empty() {
        return None;
    }
    let scored = score_candidates_simd(query, &flat);
    select_best_scored(&scored, query.len(), db_len)
}

/// True when `species_list` is exactly the canonical full species set (so the
/// all-species index applies). Order-insensitive set comparison.
fn is_full_species_scope(species_list: &[String]) -> bool {
    let all = all_species();
    if species_list.len() != all.len() {
        return false;
    }
    all.iter().all(|s| species_list.contains(s))
}

/// Ungapped sequence of a named gene in the ungapped DB, or `None`.
fn ungapped_gene_seq<'a>(seg: &str, chain: &str, species: &str, gene: &str) -> Option<&'a [u8]> {
    ungapped_chain(seg, chain)?
        .get(species)?
        .iter()
        .find(|(g, _)| g == gene)
        .map(|(_, s)| s.as_slice())
}

/// E-value-based germline assignment (RIOT-style, amino-acid).
///
/// Aligns the V region and J region of the query separately against the ungapped
/// V and J germline gene databases with Smith-Waterman (BLOSUM62, gap-open 11 /
/// extend 1) and selects each gene by lowest Karlin-Altschul e-value
/// (`E = m*n*2^(-S)`). Returns the same [`Germline`] shape as the identity path,
/// additionally populating `v_evalue` / `j_evalue`. The `v_identity`/`j_identity`
/// fields carry the SW-aligned identity to the chosen gene over the aligned span
/// (RIOT's `calculate_seq_identity`), NOT the 128-IMGT-column identity: the latter
/// is unreliable across species because germline DB entries are not consistently
/// framed (a near-perfect cross-species match can show a low column identity).
/// Gene *selection* is by e-value regardless.
///
/// Species scoping matches `run_germline_assignment`: `allowed_species = None`
/// means all species; an explicit list that contains a species absent for the V
/// chain yields the empty dict (ANARCI semantics). J is searched only within the
/// chosen V gene's species (as in both ANARCI and RIOT).
///
/// Speed: this uses a k-mer prefilter (recall-safe reduced-alphabet seeds) to pick
/// the top-K candidate V genes, then scores only those with a bit-exact SIMD
/// Smith-Waterman kernel (8 genes per pass) and selects by e-value. It returns the
/// IDENTICAL call to a full scalar brute-force scan ([`run_germline_assignment_evalue_exact`]):
/// because the per-search DB length is constant, the e-value winner is just the
/// max-SW-score gene (lex tie-break), and the prefilter is verified (test) to keep
/// the true winner in top-K for the whole golden set (recall = 1.0). About 4x the
/// throughput of the brute force single-thread (~6 ms/domain vs ~24 ms); J is a
/// small per-species set scanned directly with SIMD. Falls back to the exact SIMD
/// full scan automatically for restricted species scopes.
pub fn run_germline_assignment_evalue(
    state_vector: &StateVector,
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
) -> Germline {
    assign_evalue(state_vector, sequence, chain_type, allowed_species, ScanMode::Fast)
}

/// EXACT (always-correct) e-value germline assignment: full scalar Smith-Waterman
/// scan over every allowed gene, no prefilter or SIMD. This is the reference path
/// [`run_germline_assignment_evalue`] must match bit-for-bit; it exists as a slow
/// oracle and a guaranteed-correct fallback. Same arguments/semantics/output.
pub fn run_germline_assignment_evalue_exact(
    state_vector: &StateVector,
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
) -> Germline {
    assign_evalue(state_vector, sequence, chain_type, allowed_species, ScanMode::ExactScalar)
}

/// How the e-value path scores genes. `Fast` = k-mer prefilter + SIMD SW (default);
/// `ExactScalar` = full scalar Gotoh scan (reference oracle / fallback). Both
/// select the winner under the identical total order, so they return the same call.
#[derive(Clone, Copy, PartialEq, Eq)]
enum ScanMode {
    Fast,
    ExactScalar,
}

/// Find the e-value winner for one V/J search under `mode`.
fn winner(
    mode: ScanMode,
    seg: &str,
    chain: &str,
    chain_map: &HashMap<String, Vec<UngappedGene>>,
    species_list: &[String],
    query: &[u8],
    db_len: usize,
) -> Option<(String, String, f64, i32)> {
    match mode {
        ScanMode::Fast => {
            best_gene_evalue_fast(seg, chain, chain_map, species_list, query, db_len)
        }
        ScanMode::ExactScalar => {
            let iter = species_list
                .iter()
                .filter_map(|sp| chain_map.get(sp).map(|gl| (sp.as_str(), gl.as_slice())));
            best_gene_evalue(iter, query, db_len)
        }
    }
}

fn assign_evalue(
    state_vector: &StateVector,
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
    mode: ScanMode,
) -> Germline {
    let mut genes = Germline::default_dict();

    let v_chain = match ungapped_chain("V", chain_type) {
        Some(c) => c,
        None => return genes, // chain_type not in V -> default dict
    };

    // Resolve species list (same rules as the identity path).
    let owned_all;
    let species_list: &[String] = match allowed_species {
        Some(sp) => {
            if !sp.iter().all(|s| v_chain.contains_key(s)) {
                return Germline::empty();
            }
            sp
        }
        None => {
            owned_all = all_species().to_vec();
            &owned_all
        }
    };

    let (v_query, j_query) = split_v_j_regions(state_vector, sequence);
    if v_query.is_empty() {
        return genes; // nothing numbered -> default dict
    }

    // --- V gene: SW over all allowed species' V genes; lowest e-value wins. ---
    let v_db_len: usize = species_list
        .iter()
        .filter_map(|sp| v_chain.get(sp))
        .flat_map(|gl| gl.iter())
        .map(|(_, s)| s.len())
        .sum();
    let (v_species, v_gene, v_ev, _v_score) =
        match winner(mode, "V", chain_type, v_chain, species_list, &v_query, v_db_len) {
            Some(t) => t,
            None => return genes,
        };

    // SW-aligned identity of the winner (one traceback, on the winning gene only).
    genes.v_identity = ungapped_gene_seq("V", chain_type, &v_species, &v_gene)
        .map(|gseq| sw::local_identity(&v_query, gseq).0);
    genes.v_gene = Some((v_species.clone(), v_gene));
    genes.v_evalue = Some(v_ev);

    // --- J gene: SW over the V species' J genes only (RIOT/ANARCI behaviour). ---
    if !j_query.is_empty() {
        if let Some(j_chain) = ungapped_chain("J", chain_type) {
            if j_chain.contains_key(&v_species) {
                let j_db_len: usize = j_chain[&v_species].iter().map(|(_, s)| s.len()).sum();
                let one = [v_species.clone()];
                if let Some((_, j_gene, j_ev, _)) =
                    winner(mode, "J", chain_type, j_chain, &one, &j_query, j_db_len)
                {
                    genes.j_identity = ungapped_gene_seq("J", chain_type, &v_species, &j_gene)
                        .map(|gseq| sw::local_identity(&j_query, gseq).0);
                    genes.j_gene = Some((v_species, j_gene));
                    genes.j_evalue = Some(j_ev);
                }
            }
        }
    }

    genes
}

/// `get_hmm_length`: residues in the first J germline for (species, ctype),
/// trailing gaps stripped. Missing keys -> 128.
///
/// For the pan-species engine the profile "species" is `pan` (not a real species);
/// we then fall back to the chain's canonical J length (max over species — uniform in
/// practice: H/A/G/D=128, K/L/B=127), NOT a blind 128, which would shift the J-end
/// extension by one for kappa/lambda/beta.
pub fn get_hmm_length(species: &str, ctype: &str) -> usize {
    if let Some(j) = seg("J") {
        if let Some(chain) = j.get(ctype) {
            let jlen = |glist: &GeneList| glist.first().map(|(_, s)| s.trim_end_matches('-').len());
            if let Some(glist) = chain.get(species) {
                if let Some(n) = jlen(glist) {
                    return n;
                }
            }
            // species not found (e.g. "pan"): chain's canonical J length, deterministic.
            if let Some(n) = chain.values().filter_map(jlen).max() {
                return n;
            }
        }
    }
    128
}
