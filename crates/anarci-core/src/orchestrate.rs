//! Top-level orchestration: faithful port of ANARCI's `run_hmmer` glue,
//! `check_for_j`, `number_sequences_from_alignment`, `anarci`, `run_anarci`,
//! `number`, `validate_sequence`, `validate_numbering`.
//!
//! The HMM engine is abstracted behind [`HmmEngine`] so this layer has no FFI
//! dependency. Batch processing parallelises across sequences with rayon.

use crate::align::{parse_hmmer_query, HitRow, Hsp};
use crate::germlines::{run_germline_assignment, run_germline_assignment_evalue, Germline};
use crate::regions::{annotate_regions, RegionAnnotation};
use crate::schemes::number_sequence_from_alignment;
use crate::types::{assertion, CResult, Numbered, State, StateType};
use rayon::prelude::*;
use std::collections::BTreeSet;

/// How germline (V/J gene) assignment is performed.
///
/// * `Identity` — ANARCI's identity over the 128 IMGT match columns (byte-for-byte
///   parity with reference ANARCI). The default for the exact/`ALL` engine.
/// * `Evalue` — RIOT-style: separate Smith-Waterman alignment of the V and J
///   regions against the ungapped germline genes, best gene chosen by e-value.
///   More accurate (RIOT, Brief Bioinform 2025). The default for the pan engine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GermlineMethod {
    Identity,
    Evalue,
}

impl GermlineMethod {
    /// Parse `"identity"` / `"evalue"` (case-insensitive). Unknown -> error.
    pub fn parse(s: &str) -> CResult<Self> {
        match s.to_lowercase().as_str() {
            "identity" => Ok(GermlineMethod::Identity),
            "evalue" | "e-value" | "e_value" => Ok(GermlineMethod::Evalue),
            other => Err(assertion!(
                "Unknown germline_method '{other}'; use 'identity' or 'evalue'."
            )),
        }
    }
}

/// An in-process HMM scan engine: scan one sequence, return its HSPs.
/// Implementations must be thread-safe (`scan_one` is called from rayon workers).
pub trait HmmEngine: Sync {
    fn scan_one(&self, name: &str, seq: &[u8]) -> Vec<Hsp>;
}

/// Alignment details for one numbered domain (ANARCI's per-domain dict).
#[derive(Clone, Debug)]
pub struct DomainInfo {
    pub id: String,
    pub description: String,
    pub evalue: f64,
    pub bitscore: f64,
    pub bias: f64,
    pub query_start: Option<usize>,
    pub query_end: usize,
    pub species: String,
    pub chain_type: String,
    pub scheme: String,
    pub query_name: String,
    pub germlines: Option<Germline>, // Some only if assign_germline
    /// IMGT region-completeness annotation (F1a), computed from the domain's
    /// state vector. Always present; the Python layer exposes it only when the
    /// caller opts in via `annotate_regions=True`, keeping the default dict
    /// byte-identical to reference ANARCI.
    pub regions: RegionAnnotation,
}

/// Per-sequence result (mirrors one slot of ANARCI's three output lists).
#[derive(Clone, Debug)]
pub struct SeqResult {
    pub numbered: Option<Vec<Numbered>>,
    pub details: Option<Vec<DomainInfo>>,
    pub hit_table: Vec<HitRow>,
    /// ANARCI prints a notice to stdout when species-limiting reverts to any species.
    pub reverted_species: bool,
}

const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// `validate_sequence`: reject overly long sequences or non-standard residues.
pub fn validate_sequence(sequence: &[u8]) -> CResult<()> {
    if sequence.len() >= 10000 {
        return Err(assertion!("Sequence too long."));
    }
    let mut unknown: Vec<u8> = sequence
        .iter()
        .map(|c| c.to_ascii_uppercase())
        .filter(|c| !AMINO_ACIDS.contains(c))
        .collect();
    if !unknown.is_empty() {
        unknown.sort_unstable();
        unknown.dedup();
        let letters: Vec<String> = unknown.iter().map(|&c| (c as char).to_string()).collect();
        return Err(assertion!(
            "Unknown amino acid letter found in sequence: {}",
            letters.join(", ")
        ));
    }
    Ok(())
}

/// `validate_numbering`: indices non-decreasing, and the numbered residues form a
/// contiguous substring of the (gap-stripped) input sequence.
fn validate_numbering(num: &Numbered, name: &str, seq: &[u8]) -> CResult<()> {
    let mut last: i32 = -1;
    let mut nseq: Vec<u8> = Vec::new();
    for r in &num.residues {
        if r.num < last {
            return Err(assertion!(
                "Numbering was found to decrease along the sequence {}. Please report.",
                name
            ));
        }
        last = r.num;
        if r.aa != b'-' {
            nseq.push(r.aa);
        }
    }
    let seq_nogap: Vec<u8> = seq.iter().copied().filter(|&c| c != b'-').collect();
    if !is_subslice(&seq_nogap, &nseq) {
        return Err(assertion!(
            "The algorithm did not number a contiguous segment for sequence {}. Please report",
            name
        ));
    }
    Ok(())
}

fn is_subslice(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    if needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// `check_for_j`: rescue very long CDR3s where the first pass missed the J region.
/// Mutates the single-domain state vector in place. Needs the engine for the rescan.
fn check_for_j(
    engine: &dyn HmmEngine,
    name: &str,
    seq: &[u8],
    state_vectors: &mut [Vec<State>],
    details: &mut [crate::align::DomainDetails],
) {
    if state_vectors.len() != 1 {
        return; // single-domain chains only
    }
    let ali = &state_vectors[0];
    let last = match ali.last() {
        Some(s) => s,
        None => return,
    };
    let last_state = last.id as i32;
    let last_si = match last.si {
        Some(v) => v,
        None => return, // (Python would error; never happens for a real FW4 end)
    };
    if last_state >= 120 {
        return;
    }
    if last_si + 30 >= seq.len() {
        return;
    }
    // Sequence index of conserved cysteine (IMGT 104), as a match state.
    let mut cys_si: Option<usize> = None;
    let mut cys_ai: Option<usize> = None;
    for (ai, s) in ali.iter().enumerate() {
        if s.id == 104 && s.typ == StateType::M {
            cys_si = s.si;
            cys_ai = Some(ai);
            break;
        }
    }
    let (cys_si, cys_ai) = match (cys_si, cys_ai) {
        (Some(a), Some(b)) => (a, b),
        _ => return,
    };

    // Rescan the remaining sequence after 104 with a low bit-score threshold.
    let sub = &seq[cys_si + 1..];
    let hsps = engine.scan_one(name, sub);
    let parsed = parse_hmmer_query(&hsps, sub.len(), 10.0, None);
    if parsed.state_vectors.is_empty() {
        return;
    }
    let re = &parsed.state_vectors[0];
    let re_first = re.first().map(|s| s.id as i32).unwrap_or(0);
    let re_last = re.last().map(|s| s.id as i32).unwrap_or(0);
    if !(re_last >= 126 && re_first <= 117) {
        return;
    }

    // V region up to and including 104.
    let v_region: Vec<State> = ali[..=cys_ai].to_vec();
    // J region: re-scan states with id >= 117 and a real residue, remapped.
    let j_region: Vec<State> = re
        .iter()
        .filter(|s| s.id as i32 >= 117 && s.si.is_some())
        .map(|s| State { id: s.id, typ: s.typ, si: Some(s.si.unwrap() + cys_si + 1) })
        .collect();
    if j_region.is_empty() {
        return;
    }
    let j_first_si = j_region[0].si.unwrap();
    // CDR region between V and J.
    let mut cdr_region: Vec<State> = Vec::new();
    let mut next = 105i32;
    for si in (cys_si + 1)..j_first_si {
        if next >= 116 {
            cdr_region.push(State { id: 116, typ: StateType::I, si: Some(si) });
        } else {
            cdr_region.push(State { id: next as u8, typ: StateType::M, si: Some(si) });
            next += 1;
        }
    }

    let mut new_sv = v_region;
    new_sv.extend(cdr_region);
    let j_last_si = j_region.last().unwrap().si.unwrap();
    new_sv.extend(j_region);
    state_vectors[0] = new_sv;
    details[0].query_end = j_last_si + 1;
}

#[allow(clippy::too_many_arguments)]
fn process_one(
    engine: &dyn HmmEngine,
    name: &str,
    seq: &[u8],
    scheme: &str,
    allow: &BTreeSet<String>,
    assign_germline: bool,
    allowed_species: Option<&[String]>,
    bit_score_threshold: f64,
    species_from_germline: bool,
    germline_method: GermlineMethod,
) -> CResult<SeqResult> {
    let hsps = engine.scan_one(name, seq);
    let mut parsed = parse_hmmer_query(&hsps, seq.len(), bit_score_threshold, allowed_species);

    // check_for_j mutates the alignment in place (long-CDR3 rescue).
    check_for_j(engine, name, seq, &mut parsed.state_vectors, &mut parsed.details);

    let mut hit_numbered: Vec<Numbered> = Vec::new();
    let mut hit_details: Vec<DomainInfo> = Vec::new();

    for di in 0..parsed.state_vectors.len() {
        let sv = &parsed.state_vectors[di];
        let det = &parsed.details[di];
        if !sv.is_empty() && allow.contains(&det.chain_type) {
            let numbered =
                number_sequence_from_alignment(sv, seq, scheme, Some(&det.chain_type))?;
            validate_numbering(&numbered, name, seq)?;
            // Pan engine: the HMM "species" is "pan", so derive species (and genes)
            // from germline assignment. Exact engine: keep ANARCI semantics.
            let want_germline = assign_germline || species_from_germline;
            let germ = if want_germline {
                Some(match germline_method {
                    GermlineMethod::Identity => {
                        run_germline_assignment(sv, seq, &det.chain_type, allowed_species)
                    }
                    GermlineMethod::Evalue => {
                        run_germline_assignment_evalue(sv, seq, &det.chain_type, allowed_species)
                    }
                })
            } else {
                None
            };
            let species = if species_from_germline {
                germ.as_ref()
                    .and_then(|g| g.v_gene.as_ref().map(|(sp, _)| sp.clone()))
                    .unwrap_or_else(|| det.species.clone())
            } else {
                det.species.clone()
            };
            hit_numbered.push(numbered);
            hit_details.push(DomainInfo {
                id: det.id.clone(),
                description: det.description.clone(),
                evalue: det.evalue,
                bitscore: det.bitscore,
                bias: det.bias,
                query_start: det.query_start,
                query_end: det.query_end,
                species,
                chain_type: det.chain_type.clone(),
                scheme: scheme.to_string(),
                query_name: name.to_string(),
                germlines: germ,
                regions: annotate_regions(sv),
            });
        }
    }

    let (numbered, details) = if !hit_numbered.is_empty() {
        (Some(hit_numbered), Some(hit_details))
    } else {
        (None, None)
    };
    Ok(SeqResult {
        numbered,
        details,
        hit_table: parsed.hit_table,
        reverted_species: parsed.reverted_species,
    })
}

/// Long form of a numbering scheme name, or an error matching ANARCI's message.
pub fn resolve_scheme(scheme: &str) -> CResult<&'static str> {
    match scheme.to_lowercase().as_str() {
        "m" | "martin" => Ok("martin"),
        "c" | "chothia" => Ok("chothia"),
        "k" | "kabat" => Ok("kabat"),
        "i" | "imgt" => Ok("imgt"),
        "a" | "aho" => Ok("aho"),
        "w" | "wolfguy" => Ok("wolfguy"),
        _ => Err(assertion!("Unrecognised or unimplemented scheme: {}", scheme)),
    }
}

/// Map of chain type to output class (kappa/lambda both "L").
pub fn chain_type_to_class(ct: &str) -> &str {
    match ct {
        "H" => "H",
        "K" | "L" => "L",
        "A" => "A",
        "B" => "B",
        "G" => "G",
        "D" => "D",
        other => other,
    }
}

/// The default allowed chain types: {H,K,L,A,B,G,D}.
pub fn default_allow() -> BTreeSet<String> {
    ["H", "K", "L", "A", "B", "G", "D"].iter().map(|s| s.to_string()).collect()
}

/// `anarci`: number a list of sequences (sequential).
#[allow(clippy::too_many_arguments)]
pub fn anarci(
    engine: &dyn HmmEngine,
    sequences: &[(String, Vec<u8>)],
    scheme: &str,
    allow: &BTreeSet<String>,
    assign_germline: bool,
    allowed_species: Option<&[String]>,
    bit_score_threshold: f64,
    species_from_germline: bool,
    germline_method: GermlineMethod,
) -> CResult<Vec<SeqResult>> {
    let scheme = resolve_scheme(scheme)?;
    sequences
        .iter()
        .map(|(name, seq)| {
            process_one(
                engine, name, seq, scheme, allow, assign_germline, allowed_species,
                bit_score_threshold, species_from_germline, germline_method,
            )
        })
        .collect()
}

/// `run_anarci`: number a batch in parallel across sequences (rayon).
#[allow(clippy::too_many_arguments)]
pub fn run_anarci(
    engine: &(dyn HmmEngine + Sync),
    sequences: &[(String, Vec<u8>)],
    scheme: &str,
    allow: &BTreeSet<String>,
    assign_germline: bool,
    allowed_species: Option<&[String]>,
    bit_score_threshold: f64,
    species_from_germline: bool,
    germline_method: GermlineMethod,
) -> CResult<Vec<SeqResult>> {
    let scheme = resolve_scheme(scheme)?;

    // Deduplicate identical sequences (lossless: the result depends only on the
    // sequence; the only per-input field, query_name, is restored per original below).
    // This is a free win on repetitive inputs (NGS sets often have many duplicates).
    let mut first_index: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
    let mut uniques: Vec<(&str, &[u8])> = Vec::new();
    let mut which: Vec<usize> = Vec::with_capacity(sequences.len());
    for (name, seq) in sequences {
        let u = *first_index.entry(seq.as_slice()).or_insert_with(|| {
            uniques.push((name.as_str(), seq.as_slice()));
            uniques.len() - 1
        });
        which.push(u);
    }

    let computed: Vec<SeqResult> = uniques
        .par_iter()
        .map(|(name, seq)| {
            process_one(
                engine, name, seq, scheme, allow, assign_germline, allowed_species,
                bit_score_threshold, species_from_germline, germline_method,
            )
        })
        .collect::<CResult<Vec<_>>>()?;

    // Replay per original, restoring query_name (the only name-dependent output field).
    Ok(sequences
        .iter()
        .zip(which.iter())
        .map(|((name, _), &u)| {
            let mut r = computed[u].clone();
            if let Some(details) = r.details.as_mut() {
                for d in details {
                    d.query_name = name.clone();
                }
            }
            r
        })
        .collect())
}

/// `number`: single sequence -> (numbering, chain class) or None.
/// Returns `Ok(None)` for sequences shorter than `min_length` or that don't number,
/// and (matching ANARCI) catches the scheme/chain `AssertionError` as `Ok(None)`.
/// `min_length` (default 70 at the API layer) and `bit_score_threshold` (default 80)
/// are exposed so callers can number partial fragments (F1b); lowering them lets a
/// short or marginal-scoring fragment through instead of being silently rejected.
pub fn number(
    engine: &dyn HmmEngine,
    sequence: &[u8],
    scheme: &str,
    allow: &BTreeSet<String>,
    allowed_species: Option<&[String]>,
    species_from_germline: bool,
    min_length: usize,
    bit_score_threshold: f64,
) -> CResult<Option<(Numbered, String)>> {
    validate_sequence(sequence)?;
    let scheme = resolve_scheme(scheme)?;
    if sequence.len() < min_length {
        return Ok(None);
    }
    let seqs = vec![("sequence_0".to_string(), sequence.to_vec())];
    // number() returns only (numbering, chain class); the germline call here exists
    // solely to derive species, which number() discards — so the cheaper identity
    // method is used regardless of engine, with no effect on the returned output.
    let res = match anarci(
        engine, &seqs, scheme, allow, false, allowed_species, bit_score_threshold,
        species_from_germline, GermlineMethod::Identity,
    ) {
        Ok(r) => r,
        // ANARCI catches AssertionError here (e.g. TCR with antibody scheme) -> (False, False).
        Err(_) => return Ok(None),
    };
    let r = &res[0];
    match (&r.numbered, &r.details) {
        (Some(nums), Some(dets)) if !nums.is_empty() => {
            let class = chain_type_to_class(&dets[0].chain_type).to_string();
            Ok(Some((nums[0].clone(), class)))
        }
        _ => Ok(None),
    }
}
