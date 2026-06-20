//! HMMER-output → state-vector transforms: faithful port of `_parse_hmmer_query`,
//! `_hmm_alignment_to_states`, `_domains_are_same`. Pure (no engine dependency):
//! the input is a list of HSPs (one per profile-domain) that the engine must supply.

use crate::germlines::get_hmm_length;
use crate::types::{State, StateType, StateVector};

/// One HMMER high-scoring pair (a profile-domain hit), mirroring the Biopython
/// `hsp` fields ANARCI consumes. Coordinates are 0-based Python slice indices.
#[derive(Clone, Debug)]
pub struct Hsp {
    pub hit_id: String, // "species_chain", e.g. "mouse_H"
    pub hit_description: String,
    pub evalue: f64,
    pub bitscore: f64,
    pub bias: f64,
    pub query_start: usize,
    pub query_end: usize,
    pub hit_start: usize,
    pub hit_end: usize,
    pub rf: String, // reference annotation ('x' at match cols)
    pub pp: String, // posterior-probability ('.' at delete cols)
    pub order: usize,
}

/// Details of one identified domain (ANARCI's `top_descriptions[i]`).
#[derive(Clone, Debug)]
pub struct DomainDetails {
    pub id: String,
    pub description: String,
    pub evalue: f64,
    pub bitscore: f64,
    pub bias: f64,
    pub query_start: Option<usize>, // overwritten to first state's seq index
    pub query_end: usize,
    pub species: String,
    pub chain_type: String,
}

/// One row of the hit table (excluding the header).
#[derive(Clone, Debug)]
pub struct HitRow {
    pub id: String,
    pub description: String,
    pub evalue: f64,
    pub bitscore: f64,
    pub bias: f64,
    pub query_start: usize,
    pub query_end: usize,
}

pub struct ParsedQuery {
    pub state_vectors: Vec<StateVector>,
    pub details: Vec<DomainDetails>,
    pub hit_table: Vec<HitRow>,
    /// True if species limiting was requested but no hit cleared the threshold,
    /// so ANARCI reverted to any species (ANARCI prints a message in this case).
    pub reverted_species: bool,
}

fn split_hit_id(hit_id: &str) -> (&str, &str) {
    hit_id.split_once('_').expect("hit_id must be species_chain")
}

/// `_domains_are_same`: do the two domains overlap on the query?
fn domains_are_same(a: &Hsp, b: &Hsp) -> bool {
    let (d1, d2) = if a.query_start <= b.query_start { (a, b) } else { (b, a) };
    d2.query_start < d1.query_end
}

/// `_hmm_alignment_to_states`: convert one HSP's alignment to a state vector.
pub fn hmm_alignment_to_states(hsp: &Hsp, n: usize, seq_length: usize) -> StateVector {
    let mut reference: Vec<u8> = hsp.rf.as_bytes().to_vec();
    let mut state: Vec<u8> = hsp.pp.as_bytes().to_vec();
    assert_eq!(
        reference.len(),
        state.len(),
        "Aligned reference and state strings had different lengths."
    );

    let mut hmm_start = hsp.hit_start;
    let mut hmm_end = hsp.hit_end;
    let mut seq_start = hsp.query_start;
    let mut seq_end = hsp.query_end;

    let (species, ctype) = split_hit_id(&hsp.hit_id);
    let hmm_length = get_hmm_length(species, ctype);

    // N-terminal extension (first domain only, small unmatched N-term).
    if hsp.order == 0 && hmm_start != 0 && hmm_start < 5 {
        let mut n_extend = hmm_start;
        if hmm_start > seq_start {
            n_extend = seq_start.min(hmm_start - seq_start);
        }
        let mut ns = vec![b'8'; n_extend];
        ns.extend_from_slice(&state);
        state = ns;
        let mut nr = vec![b'x'; n_extend];
        nr.extend_from_slice(&reference);
        reference = nr;
        seq_start -= n_extend;
        hmm_start -= n_extend;
    }

    // C-terminal extension to the J-element (single domain, half of FW4 seen).
    if n == 1 && seq_end < seq_length && (123 < hmm_end && hmm_end < hmm_length) {
        let n_extend = (hmm_length - hmm_end).min(seq_length - seq_end);
        state.extend(std::iter::repeat(b'8').take(n_extend));
        reference.extend(std::iter::repeat(b'x').take(n_extend));
        seq_end += n_extend;
        hmm_end += n_extend;
    }
    let _ = (hmm_end, seq_end); // updated only to mirror the source

    let mut h = 0usize; // index into hmm states (id = hmm_start + h + 1)
    let mut s = 0usize; // index into sequence positions (si = seq_start + s)
    let mut sv: StateVector = Vec::with_capacity(state.len());
    for i in 0..state.len() {
        let base = if reference[i] == b'x' { StateType::M } else { StateType::I };
        let (typ, si) = if state[i] == b'.' {
            (StateType::D, None)
        } else {
            (base, Some(seq_start + s))
        };
        let id = (hmm_start + h + 1) as u8;
        sv.push(State { id, typ, si });
        match typ {
            StateType::M => {
                h += 1;
                s += 1;
            }
            StateType::I => {
                s += 1;
            }
            StateType::D => {
                h += 1;
            }
        }
    }
    sv
}

/// `_parse_hmmer_query`: pick the best HSP per domain, threshold, order, and build
/// state vectors. `hsps` is every profile-domain hit for one query sequence.
pub fn parse_hmmer_query(
    hsps: &[Hsp],
    seq_len: usize,
    bit_score_threshold: f64,
    hmmer_species: Option<&[String]>,
) -> ParsedQuery {
    let mut hit_table: Vec<HitRow> = Vec::new();
    let mut domains: Vec<Hsp> = Vec::new();
    let mut details: Vec<DomainDetails> = Vec::new();
    let mut reverted = false;

    if !hsps.is_empty() {
        // Species limiting (with revert-to-any, matching ANARCI).
        let owned: Vec<Hsp>;
        let hsp_list: &[Hsp] = match hmmer_species {
            Some(species) if !species.is_empty() => {
                let mut correct: Vec<Hsp> = Vec::new();
                for hsp in hsps {
                    if hsp.bitscore >= bit_score_threshold {
                        for sp in species {
                            if hsp.hit_id.starts_with(sp.as_str()) {
                                correct.push(hsp.clone());
                            }
                        }
                    }
                }
                if !correct.is_empty() {
                    owned = correct;
                    &owned
                } else {
                    reverted = true;
                    owned = hsps.to_vec();
                    &owned
                }
            }
            _ => {
                owned = hsps.to_vec();
                &owned
            }
        };

        // Sort by e-value ascending (stable).
        let mut order: Vec<usize> = (0..hsp_list.len()).collect();
        order.sort_by(|&a, &b| {
            hsp_list[a].evalue.partial_cmp(&hsp_list[b].evalue).unwrap_or(std::cmp::Ordering::Equal)
        });

        for &idx in &order {
            let hsp = &hsp_list[idx];
            if hsp.bitscore >= bit_score_threshold {
                let mut new = true;
                for d in &domains {
                    if domains_are_same(d, hsp) {
                        new = false;
                        break;
                    }
                }
                hit_table.push(HitRow {
                    id: hsp.hit_id.clone(),
                    description: hsp.hit_description.clone(),
                    evalue: hsp.evalue,
                    bitscore: hsp.bitscore,
                    bias: hsp.bias,
                    query_start: hsp.query_start,
                    query_end: hsp.query_end,
                });
                if new {
                    domains.push(hsp.clone());
                    let (species, chain) = split_hit_id(&hsp.hit_id);
                    details.push(DomainDetails {
                        id: hsp.hit_id.clone(),
                        description: hsp.hit_description.clone(),
                        evalue: hsp.evalue,
                        bitscore: hsp.bitscore,
                        bias: hsp.bias,
                        query_start: Some(hsp.query_start),
                        query_end: hsp.query_end,
                        species: species.to_string(),
                        chain_type: chain.to_string(),
                    });
                }
            }
        }

        // Reorder domains by query_start (stable).
        let mut ord: Vec<usize> = (0..domains.len()).collect();
        ord.sort_by_key(|&i| domains[i].query_start);
        domains = ord.iter().map(|&i| domains[i].clone()).collect();
        details = ord.iter().map(|&i| details[i].clone()).collect();
    }

    // Build state vectors and finalise details.
    let ndomains = domains.len();
    let mut state_vectors: Vec<StateVector> = Vec::with_capacity(ndomains);
    for i in 0..ndomains {
        domains[i].order = i;
        let sv = hmm_alignment_to_states(&domains[i], ndomains, seq_len);
        // query_start <- first state's sequence index (post-extension).
        details[i].query_start = sv.first().and_then(|st| st.si);
        state_vectors.push(sv);
    }

    ParsedQuery { state_vectors, details, hit_table, reverted_species: reverted }
}
