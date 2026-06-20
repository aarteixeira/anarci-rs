//! Germline assignment — a faithful port of ANARCI's `get_identity`,
//! `run_germline_assignment`, and `get_hmm_length`.

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
/// Cost: full SW against every allowed V gene (~1–2k genes), exactly (no lossy
/// prefilter — the e-value winner is the top SW score, kept exact). This is
/// ~10× the identity path single-threaded (~20 ms/domain across 8 species), but
/// `run_anarci` parallelises across sequences, so the default multi-core wall
/// time is comparable to the identity path. Restrict `allowed_species` to shrink
/// the DB when speed matters.
pub fn run_germline_assignment_evalue(
    state_vector: &StateVector,
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
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
    let v_iter = species_list
        .iter()
        .filter_map(|sp| v_chain.get(sp).map(|gl| (sp.as_str(), gl.as_slice())));
    let (v_species, v_gene, v_ev, _v_score) =
        match best_gene_evalue(v_iter, &v_query, v_db_len) {
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
            if let Some(jgenes) = j_chain.get(&v_species) {
                let j_db_len: usize = jgenes.iter().map(|(_, s)| s.len()).sum();
                let one = std::iter::once((v_species.as_str(), jgenes.as_slice()));
                if let Some((_, j_gene, j_ev, _)) = best_gene_evalue(one, &j_query, j_db_len) {
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
