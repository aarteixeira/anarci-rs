//! Germline assignment — a faithful port of ANARCI's `get_identity`,
//! `run_germline_assignment`, and `get_hmm_length`.

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
}

impl Germline {
    fn default_dict() -> Self {
        Germline {
            empty: false,
            v_gene: None,
            v_identity: None,
            j_gene: None,
            j_identity: None,
        }
    }
    fn empty() -> Self {
        Germline {
            empty: true,
            v_gene: None,
            v_identity: None,
            j_gene: None,
            j_identity: None,
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
