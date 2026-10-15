//! Private native computations for the anarci-rs Python packages.
//!
//! Mirrors ANARCI's `anarci`, `run_anarci`, `number` exactly (names, signatures,
//! return shapes). The HMM scan runs in-process (native HMMER 3.4); batch runs in
//! Rust with rayon. The GIL is released around compute.

use anarci_core::orchestrate::{self, DomainInfo, SeqResult};
use anarci_core::{Germline, GermlineMethod, HmmEngine, Hsp, Numbered};
use anarci_hmm::Engine;
mod experimental;
mod features;
use once_cell::sync::OnceCell;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::collections::BTreeSet;

// ---- engine singletons (bundled HMMs, loaded once) ------------------------
//
// Two engines, both embedded so the wheel is self-contained:
//  - PAN (default): FEW.hmm, one pan-species HMM per chain type (7 profiles). ~4-12x
//    faster; numbering equal-or-better than ANARCI; species/genes via germline assignment.
//  - ALL (database="ALL"): the 29 species×chain profiles — byte-for-byte ANARCI.

static ENGINE_ALL: OnceCell<Engine> = OnceCell::new();
static ENGINE_PAN: OnceCell<Engine> = OnceCell::new();

static ALL_HMM: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference_data/dat/HMMs/ALL.hmm"
));
static PAN_HMM: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference_data/dat/HMMs/FEW.hmm"
));

fn load_engine(
    cell: &'static OnceCell<Engine>,
    data: &[u8],
    _fname: &str,
) -> PyResult<&'static Engine> {
    cell.get_or_try_init(|| {
        Engine::load_embedded(data).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("failed to load HMM engine: {e}"))
        })
    })
}

fn engine(pan: bool) -> PyResult<&'static Engine> {
    if pan {
        load_engine(&ENGINE_PAN, PAN_HMM, "FEW.hmm")
    } else {
        load_engine(&ENGINE_ALL, ALL_HMM, "ALL.hmm")
    }
}

fn custom_engine(database: &str) -> PyResult<Option<Engine>> {
    if !std::path::Path::new(database).is_file() {
        return Ok(None);
    }
    let engine = Engine::load(std::path::Path::new(database))
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    engine
        .validate_numbering_profiles()
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    Ok(Some(engine))
}

/// `database` -> pan(default)/exact. Unknown values error (no silent fallback).
fn parse_database(database: &str) -> PyResult<bool> {
    match database {
        "pan" | "PAN" => Ok(true),
        "ALL" => Ok(false),
        other if std::path::Path::new(other).is_file() => Ok(false),
        other => Err(assertion_err(format!(
            "Unknown database '{other}'; use 'pan' (default, fast) or 'ALL' (exact ANARCI)."
        ))),
    }
}

/// Resolve the germline-assignment method.
///
/// When unspecified (`None`), the default depends on the engine:
///   * pan  -> `Evalue` (RIOT-style V/J assignment; more accurate germline genes).
///     The e-value path is now fast (k-mer prefilter + SIMD Smith-Waterman, ~10x
///     the throughput of the old brute force), so it no longer erases the pan
///     speedup. It produces IDENTICAL calls to the brute-force e-value path.
///   * exact (`database="ALL"`) -> `Identity` (byte-for-byte ANARCI parity).
///
/// An explicit `germline_method=` always overrides the default.
fn parse_germline_method(method: Option<&str>, pan: bool) -> PyResult<GermlineMethod> {
    match method {
        None => Ok(if pan {
            GermlineMethod::Evalue
        } else {
            GermlineMethod::Identity
        }),
        Some(s) => GermlineMethod::parse(s).map_err(assertion_err),
    }
}

/// Adapter so the orchestration layer can use the FFI engine via the `HmmEngine` trait.
struct Adapter<'a>(&'a Engine);
impl HmmEngine for Adapter<'_> {
    fn scan_one(&self, name: &str, seq: &[u8]) -> Vec<Hsp> {
        self.0.scan(name, seq)
    }
    fn scan_one_checked(&self, name: &str, seq: &[u8]) -> anarci_core::CResult<Vec<Hsp>> {
        self.0
            .scan_checked(name, seq)
            .map_err(|e| anarci_core::CoreError::Engine(e.to_string()))
    }
}

// ---- helpers --------------------------------------------------------------

fn assertion_err(msg: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyAssertionError::new_err(msg.to_string())
}

fn core_err(err: anarci_core::CoreError) -> PyErr {
    match err {
        anarci_core::CoreError::Assertion(msg) => assertion_err(msg),
        anarci_core::CoreError::Engine(msg) => pyo3::exceptions::PyRuntimeError::new_err(msg),
    }
}

fn species_notice(py: Python<'_>, species: Option<&[String]>) -> PyResult<()> {
    let representation = species
        .map(|s| s.to_vec())
        .into_pyobject(py)?
        .repr()?
        .to_str()?
        .to_owned();
    py.import("builtins")?.getattr("print")?.call1((format!(
        "Limiting hmmer search to species {representation} was requested but hits did not achieve a high enough bitscore. Reverting to using any species"),))?;
    Ok(())
}

fn default_allow_set() -> BTreeSet<String> {
    orchestrate::default_allow()
}

/// Extract `allow` (a Python set/iterable of str) -> BTreeSet, default {H,K,L,A,B,G,D}.
fn parse_allow(allow: Option<&Bound<'_, PyAny>>) -> PyResult<BTreeSet<String>> {
    match allow {
        None => Ok(default_allow_set()),
        Some(obj) => {
            let mut set = BTreeSet::new();
            for item in obj.try_iter()? {
                set.insert(item?.extract::<String>()?);
            }
            Ok(set)
        }
    }
}

/// Extract `allowed_species` (list of str, or None/empty -> None meaning "all").
fn parse_species(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<String>>> {
    match obj {
        None => Ok(None),
        Some(o) => {
            if o.is_none() {
                return Ok(None);
            }
            let v: Vec<String> = o.extract()?;
            if v.is_empty() {
                Ok(None)
            } else {
                Ok(Some(v))
            }
        }
    }
}

/// Parse the `sequences` argument (list of (id, seq) tuples) into owned bytes.
fn parse_sequences(obj: &Bound<'_, PyAny>) -> PyResult<Vec<(String, Vec<u8>)>> {
    let mut out = Vec::new();
    for item in obj.try_iter()? {
        let item = item?;
        let tup: (String, String) = item.extract().map_err(|_| {
            assertion_err("If list or tuple supplied as input format must be [ ('ID1','seq1'), ('ID2', 'seq2'), ... ]")
        })?;
        out.push((tup.0, tup.1.into_bytes()));
    }
    Ok(out)
}

fn read_fasta_file(path: &str) -> PyResult<Vec<(String, Vec<u8>)>> {
    let text = std::fs::read_to_string(path)?;
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name: Option<String> = None;
    let mut buf: Vec<u8> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push((n, std::mem::take(&mut buf)));
            }
            name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else if !line.is_empty() {
            buf.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(n) = name.take() {
        out.push((n, buf));
    }
    Ok(out)
}

// ---- result -> Python conversion (exact ANARCI shapes) --------------------

fn numbering_to_py<'py>(py: Python<'py>, num: &Numbered) -> PyResult<Bound<'py, PyAny>> {
    // (numbering, start, end) where numbering = [ ((pos, ins), aa), ... ]
    let residues: Vec<Bound<'py, PyAny>> = num
        .residues
        .iter()
        .map(|r| {
            let pos = PyTuple::new(
                py,
                [
                    r.num.into_pyobject(py)?.into_any(),
                    r.ins.into_pyobject(py)?.into_any(),
                ],
            )?;
            let aa = (r.aa as char).to_string();
            Ok(PyTuple::new(py, [pos.into_any(), aa.into_pyobject(py)?.into_any()])?.into_any())
        })
        .collect::<PyResult<_>>()?;
    let numbering = PyList::new(py, residues)?;
    let start = opt_usize(py, num.start)?;
    let end = opt_usize(py, num.end)?;
    Ok(PyTuple::new(py, [numbering.into_any(), start, end])?.into_any())
}

fn opt_usize<'py>(py: Python<'py>, v: Option<usize>) -> PyResult<Bound<'py, PyAny>> {
    Ok(match v {
        Some(x) => x.into_pyobject(py)?.into_any(),
        None => py.None().into_bound(py),
    })
}

fn germline_to_py<'py>(py: Python<'py>, g: &Germline) -> PyResult<Bound<'py, PyAny>> {
    let d = PyDict::new(py);
    if g.empty {
        return Ok(d.into_any()); // {}
    }
    let gene_obj = |gene: &Option<(String, String)>| -> PyResult<Bound<'py, PyAny>> {
        Ok(match gene {
            Some((sp, ge)) => PyTuple::new(
                py,
                [
                    sp.into_pyobject(py)?.into_any(),
                    ge.into_pyobject(py)?.into_any(),
                ],
            )?
            .into_any(),
            None => py.None().into_bound(py),
        })
    };
    let v = PyList::new(py, [gene_obj(&g.v_gene)?, opt_f64(py, g.v_identity)?])?;
    let j = PyList::new(py, [gene_obj(&g.j_gene)?, opt_f64(py, g.j_identity)?])?;
    d.set_item("v_gene", v)?;
    d.set_item("j_gene", j)?;
    // E-values only on the alignment-based path (None on identity), so the
    // ANARCI-parity dict is byte-identical to before.
    if g.v_evalue.is_some() || g.j_evalue.is_some() {
        d.set_item("v_evalue", opt_f64(py, g.v_evalue)?)?;
        d.set_item("j_evalue", opt_f64(py, g.j_evalue)?)?;
    }
    Ok(d.into_any())
}

fn opt_f64<'py>(py: Python<'py>, v: Option<f64>) -> PyResult<Bound<'py, PyAny>> {
    Ok(match v {
        Some(x) => x.into_pyobject(py)?.into_any(),
        None => py.None().into_bound(py),
    })
}

fn detail_to_py<'py>(
    py: Python<'py>,
    d: &DomainInfo,
    annotate_regions: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let dict = PyDict::new(py);
    dict.set_item("id", &d.id)?;
    dict.set_item("description", &d.description)?;
    dict.set_item("evalue", d.evalue)?;
    dict.set_item("bitscore", d.bitscore)?;
    dict.set_item("bias", d.bias)?;
    dict.set_item("query_start", opt_usize(py, d.query_start)?)?;
    dict.set_item("query_end", d.query_end)?;
    dict.set_item("species", &d.species)?;
    dict.set_item("chain_type", &d.chain_type)?;
    dict.set_item("scheme", &d.scheme)?;
    dict.set_item("query_name", &d.query_name)?;
    if let Some(g) = &d.germlines {
        dict.set_item("germlines", germline_to_py(py, g)?)?;
    }
    // F1a: IMGT region-completeness annotation. Opt-in so the default dict stays
    // byte-identical to reference ANARCI (drop-in parity).
    if annotate_regions {
        let regions = PyDict::new(py);
        for (name, status) in d.regions.pairs() {
            regions.set_item(name, status.as_str())?;
        }
        dict.set_item("regions", regions)?;
        let covered = match d.regions.covered_imgt {
            Some((lo, hi)) => PyTuple::new(py, [lo, hi])?.into_any(),
            None => py.None().into_bound(py),
        };
        dict.set_item("covered_imgt", covered)?;
    }
    Ok(dict.into_any())
}

fn hit_table_to_py<'py>(
    py: Python<'py>,
    rows: &[anarci_core::align::HitRow],
) -> PyResult<Bound<'py, PyList>> {
    let header = PyList::new(
        py,
        [
            "id",
            "description",
            "evalue",
            "bitscore",
            "bias",
            "query_start",
            "query_end",
        ],
    )?;
    let mut out: Vec<Bound<'py, PyAny>> = vec![header.into_any()];
    for r in rows {
        let row = PyList::new(
            py,
            [
                r.id.clone().into_pyobject(py)?.into_any(),
                r.description.clone().into_pyobject(py)?.into_any(),
                r.evalue.into_pyobject(py)?.into_any(),
                r.bitscore.into_pyobject(py)?.into_any(),
                r.bias.into_pyobject(py)?.into_any(),
                r.query_start.into_pyobject(py)?.into_any(),
                r.query_end.into_pyobject(py)?.into_any(),
            ],
        )?;
        out.push(row.into_any());
    }
    PyList::new(py, out)
}

/// Build the three output lists (numbered, alignment_details, hit_tables) from results.
fn results_to_py<'py>(
    py: Python<'py>,
    results: &[SeqResult],
    annotate_regions: bool,
) -> PyResult<(Bound<'py, PyList>, Bound<'py, PyList>, Bound<'py, PyList>)> {
    let mut numbered: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    let mut details: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    let mut hits: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    for r in results {
        match &r.numbered {
            None => numbered.push(py.None().into_bound(py)),
            Some(doms) => {
                let v: Vec<Bound<'py, PyAny>> = doms
                    .iter()
                    .map(|n| numbering_to_py(py, n))
                    .collect::<PyResult<_>>()?;
                numbered.push(PyList::new(py, v)?.into_any());
            }
        }
        match &r.details {
            None => details.push(py.None().into_bound(py)),
            Some(ds) => {
                let v: Vec<Bound<'py, PyAny>> = ds
                    .iter()
                    .map(|d| detail_to_py(py, d, annotate_regions))
                    .collect::<PyResult<_>>()?;
                details.push(PyList::new(py, v)?.into_any());
            }
        }
        hits.push(hit_table_to_py(py, &r.hit_table)?.into_any());
    }
    Ok((
        PyList::new(py, numbered)?,
        PyList::new(py, details)?,
        PyList::new(py, hits)?,
    ))
}

/// Run the core batch with `ncpu` rayon threads (0/None -> global pool), GIL released.
#[allow(clippy::too_many_arguments)]
fn run_core(
    py: Python<'_>,
    seqs: &[(String, Vec<u8>)],
    scheme: &str,
    allow: &BTreeSet<String>,
    assign_germline: bool,
    species: Option<&[String]>,
    bit_score_threshold: f64,
    ncpu: usize,
    pan: bool,
    germline_method: GermlineMethod,
    database: &str,
) -> PyResult<Vec<SeqResult>> {
    let custom = custom_engine(database)?;
    let eng = match custom.as_ref() {
        Some(e) => e,
        None => engine(pan)?,
    };
    let res = py.allow_threads(|| {
        let run = || {
            orchestrate::run_anarci(
                &Adapter(eng),
                seqs,
                scheme,
                allow,
                assign_germline,
                species,
                bit_score_threshold,
                pan,
                germline_method,
            )
        };
        if ncpu >= 1 {
            match rayon::ThreadPoolBuilder::new().num_threads(ncpu).build() {
                Ok(pool) => pool.install(run),
                Err(e) => Err(anarci_core::CoreError::Engine(format!(
                    "Unable to create worker pool: {e}"
                ))),
            }
        } else {
            run()
        }
    });
    let results = res.map_err(core_err)?;
    for result in &results {
        if result.reverted_species {
            species_notice(py, species)?;
        }
    }
    Ok(results)
}

// ---- public API -----------------------------------------------------------

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (sequences, scheme="imgt", database="pan", output=false, outfile=None,
    csv=false, allow=None, hmmerpath="", ncpu=None, assign_germline=false,
    allowed_species=None, bit_score_threshold=80.0, germline_method=None,
    annotate_regions=false))]
fn anarci<'py>(
    py: Python<'py>,
    sequences: &Bound<'py, PyAny>,
    scheme: &str,
    database: &str,
    output: bool,
    outfile: Option<String>,
    csv: bool,
    allow: Option<&Bound<'py, PyAny>>,
    hmmerpath: &str,
    ncpu: Option<usize>,
    assign_germline: bool,
    allowed_species: Option<&Bound<'py, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<&str>,
    annotate_regions: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    if output || csv || outfile.is_some() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Use anarci_rs.anarci for output handling",
        ));
    }
    if !hmmerpath.is_empty() {
        return Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "hmmerpath is unavailable with the in-process engine",
        ));
    }
    let pan = parse_database(database)?;
    let gmethod = parse_germline_method(germline_method, pan)?;
    let seqs = parse_sequences(sequences)?;
    let allow = parse_allow(allow)?;
    let species = parse_species(allowed_species)?;
    let results = run_core(
        py,
        &seqs,
        scheme,
        &allow,
        assign_germline,
        species.as_deref(),
        bit_score_threshold,
        ncpu.unwrap_or(1),
        pan,
        gmethod,
        database,
    )?;
    let (n, d, h) = results_to_py(py, &results, annotate_regions)?;
    PyTuple::new(py, [n.into_any(), d.into_any(), h.into_any()])
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (seq, scheme="imgt", database="pan", output=false, outfile=None,
    csv=false, allow=None, hmmerpath="", ncpu=1, assign_germline=false,
    allowed_species=None, bit_score_threshold=80.0, germline_method=None,
    annotate_regions=false))]
fn run_anarci<'py>(
    py: Python<'py>,
    seq: &Bound<'py, PyAny>,
    scheme: &str,
    database: &str,
    output: bool,
    outfile: Option<String>,
    csv: bool,
    allow: Option<&Bound<'py, PyAny>>,
    hmmerpath: &str,
    ncpu: usize,
    assign_germline: bool,
    allowed_species: Option<&Bound<'py, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<&str>,
    annotate_regions: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    if output || csv || outfile.is_some() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Use anarci_rs.run_anarci for output handling",
        ));
    }
    if !hmmerpath.is_empty() {
        return Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "hmmerpath is unavailable with the in-process engine",
        ));
    }
    let pan = parse_database(database)?;
    let gmethod = parse_germline_method(germline_method, pan)?;
    // Input: list of (id,seq) | fasta path | single sequence string.
    let seqs: Vec<(String, Vec<u8>)> = if let Ok(s) = seq.extract::<String>() {
        if std::path::Path::new(&s).is_file() {
            read_fasta_file(&s)?
        } else {
            anarci_core::validate_sequence(s.as_bytes()).map_err(assertion_err)?;
            vec![("Input sequence".to_string(), s.into_bytes())]
        }
    } else {
        parse_sequences(seq)?
    };
    let allow = parse_allow(allow)?;
    let species = parse_species(allowed_species)?;
    let results = run_core(
        py,
        &seqs,
        scheme,
        &allow,
        assign_germline,
        species.as_deref(),
        bit_score_threshold,
        ncpu.max(1),
        pan,
        gmethod,
        database,
    )?;
    let (n, d, h) = results_to_py(py, &results, annotate_regions)?;
    // Sequences out: [(id, seq_str), ...]
    let seq_list: Vec<Bound<'py, PyAny>> = seqs
        .iter()
        .map(|(id, s)| {
            let t = PyTuple::new(py, [id.clone(), String::from_utf8_lossy(s).into_owned()])?;
            Ok(t.into_any())
        })
        .collect::<PyResult<_>>()?;
    let seq_out = PyList::new(py, seq_list)?;
    PyTuple::new(
        py,
        [seq_out.into_any(), n.into_any(), d.into_any(), h.into_any()],
    )
}

#[pyfunction]
#[pyo3(signature = (sequence, scheme="imgt", database="pan", allow=None,
    allowed_species=Some(vec!["human".into(), "mouse".into()]), min_length=70,
    bit_score_threshold=80.0))]
fn number<'py>(
    py: Python<'py>,
    sequence: &str,
    scheme: &str,
    database: &str,
    allow: Option<&Bound<'py, PyAny>>,
    allowed_species: Option<Vec<String>>,
    min_length: usize,
    bit_score_threshold: f64,
) -> PyResult<Bound<'py, PyTuple>> {
    let pan = parse_database(database)?;
    let allow = parse_allow(allow)?;
    let species = allowed_species.filter(|s| !s.is_empty());
    let custom = custom_engine(database)?;
    let eng = match custom.as_ref() {
        Some(e) => e,
        None => engine(pan)?,
    };
    let byte_len = sequence.as_bytes().len();
    let out = py.allow_threads(|| {
        anarci_core::number(
            &Adapter(eng),
            sequence.as_bytes(),
            scheme,
            &allow,
            species.as_deref(),
            pan,
            min_length,
            bit_score_threshold,
        )
    });
    match out {
        Err(error) => Err(core_err(error)),
        // F1b: no numbering. Surface *why* via a UserWarning (return value stays the
        // byte-identical (False, False)) so the rejection is never silent.
        Ok(None) => {
            let reason = if byte_len < min_length {
                format!(
                    "anarci_rs.number: sequence length {byte_len} < min_length {min_length}; \
                     no numbering returned. Pass a smaller min_length= to number short/partial fragments."
                )
            } else {
                format!(
                    "anarci_rs.number: no domain scored >= bit_score_threshold {bit_score_threshold} \
                     (or no recognizable variable domain); no numbering returned. \
                     Lower bit_score_threshold= to number marginal/partial fragments."
                )
            };
            let cat = py.get_type::<pyo3::exceptions::PyUserWarning>();
            PyErr::warn(py, &cat, &std::ffi::CString::new(reason)?, 1)?;
            Ok(PyTuple::new(py, [false, false])?)
        }
        Ok(Some((num, class))) => {
            let numbering = numbering_to_py(py, &num)?;
            // number() returns just the numbering list (not the (num,start,end) tuple).
            let numbering_list = numbering.get_item(0)?;
            Ok(PyTuple::new(
                py,
                [numbering_list, class.into_pyobject(py)?.into_any()],
            )?)
        }
    }
}

// ---- validated lower-level compatibility helpers -------------------------

type PyState = ((i32, String), Option<usize>);
fn states_from_py(raw: Vec<PyState>, sequence: Option<&str>) -> PyResult<anarci_core::StateVector> {
    use anarci_core::{State, StateType};
    let mut out = Vec::with_capacity(raw.len());
    let mut last_id = 0;
    let mut last_si = None;
    for ((id, typ), si) in raw {
        if !(1..=128).contains(&id) || id < last_id {
            return Err(assertion_err(
                "State identifiers must increase within 1..128",
            ));
        }
        let typ = match typ.as_str() {
            "m" => StateType::M,
            "i" => StateType::I,
            "d" => StateType::D,
            _ => return Err(assertion_err("State type must be m, i or d")),
        };
        if (typ == StateType::D) != si.is_none() {
            return Err(assertion_err(
                "Only deletion states must have a None sequence index",
            ));
        }
        if let Some(i) = si {
            if sequence.is_some_and(|seq| !seq.is_ascii() || i >= seq.len()) {
                return Err(assertion_err("State sequence index is out of bounds"));
            }
            if last_si.is_some_and(|last| i != last + 1) {
                return Err(assertion_err("State sequence indices must be contiguous"));
            }
            last_si = Some(i);
        }
        last_id = id;
        out.push(State {
            id: id as u8,
            typ,
            si,
        });
    }
    Ok(out)
}

fn states_to_raw(states: &[anarci_core::State]) -> Vec<PyState> {
    states
        .iter()
        .map(|s| ((s.id as i32, s.typ.as_char().to_string()), s.si))
        .collect()
}

#[pyfunction]
#[pyo3(signature = (state_vector, sequence, scheme="imgt", chain_type=None))]
fn _number_alignment<'py>(
    py: Python<'py>,
    state_vector: Vec<PyState>,
    sequence: &str,
    scheme: &str,
    chain_type: Option<&str>,
) -> PyResult<Bound<'py, PyAny>> {
    let states = states_from_py(state_vector, Some(sequence))?;
    let out = std::panic::catch_unwind(|| {
        anarci_core::number_sequence_from_alignment(
            &states,
            sequence.as_bytes(),
            scheme,
            chain_type,
        )
    })
    .map_err(|_| assertion_err("State vector exceeds numbering scheme limits"))?
    .map_err(assertion_err)?;
    numbering_to_py(py, &out)
}

#[pyfunction]
fn _smooth_insertions(state_vector: Vec<PyState>) -> PyResult<Vec<PyState>> {
    let states = states_from_py(state_vector, None)?;
    let out = std::panic::catch_unwind(|| anarci_core::schemes::smooth_insertions(&states))
        .map_err(|_| assertion_err("State vector cannot be smoothed"))?;
    Ok(states_to_raw(&out))
}

#[pyfunction]
fn _get_identity(state_sequence: &str, germline_sequence: &str) -> PyResult<f64> {
    if !state_sequence.is_ascii()
        || !germline_sequence.is_ascii()
        || state_sequence.len() != 128
        || germline_sequence.len() != 128
    {
        return Err(assertion_err(""));
    }
    Ok(anarci_core::get_identity(
        state_sequence.as_bytes(),
        germline_sequence.as_bytes(),
    ))
}

#[pyfunction]
fn _get_hmm_length(species: &str, ctype: &str) -> usize {
    anarci_core::get_hmm_length(species, ctype)
}

#[pyfunction]
#[pyo3(signature = (state_vector, sequence, chain_type, allowed_species=None))]
fn _germline_assignment<'py>(
    py: Python<'py>,
    state_vector: Vec<PyState>,
    sequence: &str,
    chain_type: &str,
    allowed_species: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyAny>> {
    let states = states_from_py(state_vector, Some(sequence))?;
    if allowed_species.as_ref().is_some_and(|s| s.is_empty()) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "max() arg is an empty sequence",
        ));
    }
    germline_to_py(
        py,
        &anarci_core::run_germline_assignment(
            &states,
            sequence.as_bytes(),
            chain_type,
            allowed_species.as_deref(),
        ),
    )
}

#[pyfunction]
#[pyo3(signature = (sequence_list, hmm_database="ALL", ncpu=None, bit_score_threshold=80.0, hmmer_species=None))]
fn _run_hmmer<'py>(
    py: Python<'py>,
    sequence_list: &Bound<'py, PyAny>,
    hmm_database: &str,
    ncpu: Option<usize>,
    bit_score_threshold: f64,
    hmmer_species: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyList>> {
    use rayon::prelude::*;
    let seqs = parse_sequences(sequence_list)?;
    for (name, seq) in &seqs {
        orchestrate::validate_scan_input(name, seq).map_err(assertion_err)?;
    }
    if !bit_score_threshold.is_finite() {
        return Err(assertion_err("bit_score_threshold must be finite"));
    }
    let pan = parse_database(hmm_database)?;
    let custom = custom_engine(hmm_database)?;
    let eng = match custom.as_ref() {
        Some(e) => e,
        None => engine(pan)?,
    };
    let run = || {
        seqs.par_iter()
            .map(|(name, seq)| {
                let hsps = eng.scan_checked(name, seq).map_err(|e| e.to_string())?;
                Ok(anarci_core::parse_hmmer_query(
                    &hsps,
                    seq.len(),
                    bit_score_threshold,
                    hmmer_species.as_deref(),
                ))
            })
            .collect::<Result<Vec<_>, String>>()
    };
    let results = py
        .allow_threads(|| match ncpu {
            Some(n) if n > 0 => rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| e.to_string())?
                .install(run),
            _ => run(),
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    let out = PyList::empty(py);
    for r in results {
        if r.reverted_species {
            species_notice(py, hmmer_species.as_deref())?;
        }
        let details = PyList::empty(py);
        for d in &r.details {
            let dict = PyDict::new(py);
            dict.set_item("id", &d.id)?;
            dict.set_item("description", &d.description)?;
            dict.set_item("evalue", d.evalue)?;
            dict.set_item("bitscore", d.bitscore)?;
            dict.set_item("bias", d.bias)?;
            dict.set_item("query_start", d.query_start)?;
            dict.set_item("query_end", d.query_end)?;
            dict.set_item("species", &d.species)?;
            dict.set_item("chain_type", &d.chain_type)?;
            details.append(dict)?;
        }
        let vectors = PyList::new(py, r.state_vectors.iter().map(|s| states_to_raw(s)))?;
        out.append(PyTuple::new(
            py,
            [
                hit_table_to_py(py, &r.hit_table)?.into_any(),
                vectors.into_any(),
                details.into_any(),
            ],
        )?)?;
    }
    Ok(out)
}

#[pyfunction]
fn _germlines_json() -> &'static str {
    include_str!("../../anarci-core/data/germlines.json")
}

#[pyfunction]
fn _blosum62() -> Vec<((String, String), i32)> {
    anarci_core::constants::BLOSUM62_RAW
        .iter()
        .map(|&(a, b, v)| (((a as char).to_string(), (b as char).to_string()), v))
        .collect()
}

#[pyfunction]
fn _get_imgt_cdr(
    length: usize,
    maxlength: usize,
    start: i32,
    end: i32,
) -> PyResult<Vec<Option<(i32, &'static str)>>> {
    if length > 10000
        || maxlength > 10000
        || (maxlength < 2 && length > 1)
        || length.saturating_sub(maxlength) > 104
        || start < -10000
        || end > 10000
        || end <= start
    {
        return Err(assertion_err("Invalid CDR annotation bounds"));
    }
    std::panic::catch_unwind(|| anarci_core::schemes::get_imgt_cdr(length, maxlength, start, end))
        .map_err(|_| assertion_err("CDR annotation exceeds insertion limits"))
}

#[pyfunction]
fn _get_cdr3_annotations(
    length: usize,
    scheme: &str,
    chain_type: &str,
) -> PyResult<Vec<(i32, &'static str)>> {
    if length > 10000 {
        return Err(assertion_err(
            "Too many insertions for numbering scheme to handle",
        ));
    }
    anarci_core::schemes::get_cdr3_annotations(length, scheme, chain_type).map_err(assertion_err)
}

#[pyfunction]
fn _get_wolfguy_l1(sequence: Vec<((i32, String), String)>, length: usize) -> PyResult<Vec<i32>> {
    let mut residues = Vec::with_capacity(sequence.len());
    for ((num, _), aa) in sequence {
        if aa.len() != 1 || !aa.is_ascii() {
            return Err(assertion_err("Expected one ASCII amino acid per residue"));
        }
        residues.push(anarci_core::Res::new(num, " ", aa.as_bytes()[0]));
    }
    if length > residues.len() {
        return Err(assertion_err("CDR length exceeds sequence length"));
    }
    anarci_core::schemes::get_wolfguy_l1(&residues, length).map_err(assertion_err)
}

#[pyfunction]
fn _check_for_j(
    py: Python<'_>,
    name: &str,
    sequence: &str,
    state_vectors: Vec<Vec<PyState>>,
    query_ends: Vec<usize>,
) -> PyResult<(Vec<Vec<PyState>>, Vec<usize>)> {
    if state_vectors.len() != query_ends.len() {
        return Err(assertion_err(
            "State vectors and domain details must have matching lengths",
        ));
    }
    let mut states = state_vectors
        .into_iter()
        .map(|s| states_from_py(s, Some(sequence)))
        .collect::<PyResult<Vec<_>>>()?;
    let mut details = query_ends
        .into_iter()
        .map(|query_end| anarci_core::DomainDetails {
            id: String::new(),
            description: String::new(),
            evalue: 0.,
            bitscore: 0.,
            bias: 0.,
            query_start: None,
            query_end,
            species: String::new(),
            chain_type: String::new(),
        })
        .collect::<Vec<_>>();
    let eng = engine(false)?;
    py.allow_threads(|| {
        orchestrate::check_for_j(
            &Adapter(eng),
            name,
            sequence.as_bytes(),
            &mut states,
            &mut details,
        )
    })
    .map_err(core_err)?;
    Ok((
        states.iter().map(|s| states_to_raw(s)).collect(),
        details.iter().map(|d| d.query_end).collect(),
    ))
}

type RawHsp = (
    String,
    String,
    f64,
    f64,
    f64,
    usize,
    usize,
    usize,
    usize,
    String,
    String,
    usize,
);
fn hsp_from_raw(raw: RawHsp, seq_length: usize) -> PyResult<Hsp> {
    let (
        hit_id,
        hit_description,
        evalue,
        bitscore,
        bias,
        query_start,
        query_end,
        hit_start,
        hit_end,
        rf,
        pp,
        order,
    ) = raw;
    let valid_id = match hit_id.split_once('_') {
        Some((species, chain)) => {
            !species.is_empty() && ["H", "K", "L", "A", "B", "G", "D"].contains(&chain)
        }
        None => false,
    };
    if !valid_id
        || !evalue.is_finite()
        || !bitscore.is_finite()
        || !bias.is_finite()
        || query_start > query_end
        || query_end > seq_length
        || hit_start >= hit_end
        || hit_end > 128
        || !rf.is_ascii()
        || !pp.is_ascii()
        || rf.len() != pp.len()
        || rf.bytes().filter(|&b| b == b'x').count() != hit_end - hit_start
        || pp.bytes().filter(|&b| b != b'.').count() != query_end - query_start
    {
        return Err(assertion_err(
            "Invalid HMM alignment coordinates or annotation lengths",
        ));
    }
    let mut position = hit_start + 1;
    for ch in rf.bytes() {
        if position > 128 {
            return Err(assertion_err("HMM states exceed 128"));
        }
        if ch == b'x' {
            position += 1;
        }
    }
    Ok(Hsp {
        hit_id,
        hit_description,
        evalue,
        bitscore,
        bias,
        query_start,
        query_end,
        hit_start,
        hit_end,
        rf,
        pp,
        order,
    })
}

fn parsed_to_py<'py>(
    py: Python<'py>,
    r: &anarci_core::ParsedQuery,
) -> PyResult<Bound<'py, PyTuple>> {
    let details = PyList::empty(py);
    for d in &r.details {
        let dict = PyDict::new(py);
        dict.set_item("id", &d.id)?;
        dict.set_item("description", &d.description)?;
        dict.set_item("evalue", d.evalue)?;
        dict.set_item("bitscore", d.bitscore)?;
        dict.set_item("bias", d.bias)?;
        dict.set_item("query_start", d.query_start)?;
        dict.set_item("query_end", d.query_end)?;
        dict.set_item("species", &d.species)?;
        dict.set_item("chain_type", &d.chain_type)?;
        details.append(dict)?;
    }
    let vectors = PyList::new(py, r.state_vectors.iter().map(|s| states_to_raw(s)))?;
    PyTuple::new(
        py,
        [
            hit_table_to_py(py, &r.hit_table)?.into_any(),
            vectors.into_any(),
            details.into_any(),
        ],
    )
}

#[pyfunction]
#[pyo3(signature = (hsps, seq_length, bit_score_threshold=80.0, hmmer_species=None))]
fn _parse_hsps<'py>(
    py: Python<'py>,
    hsps: Vec<RawHsp>,
    seq_length: usize,
    bit_score_threshold: f64,
    hmmer_species: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyTuple>> {
    if !bit_score_threshold.is_finite() {
        return Err(assertion_err("bit_score_threshold must be finite"));
    }
    let hsps = hsps
        .into_iter()
        .map(|h| hsp_from_raw(h, seq_length))
        .collect::<PyResult<Vec<_>>>()?;
    let result = anarci_core::parse_hmmer_query(
        &hsps,
        seq_length,
        bit_score_threshold,
        hmmer_species.as_deref(),
    );
    if result.reverted_species {
        species_notice(py, hmmer_species.as_deref())?;
    }
    parsed_to_py(py, &result)
}

#[pyfunction]
fn _hmm_states(hsp: RawHsp, n: usize, seq_length: usize) -> PyResult<Vec<PyState>> {
    let hsp = hsp_from_raw(hsp, seq_length)?;
    Ok(states_to_raw(&anarci_core::hmm_alignment_to_states(
        &hsp, n, seq_length,
    )))
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn _number_regions<'py>(
    py: Python<'py>,
    sequence: &str,
    state_vector: Vec<PyState>,
    state_string: &str,
    region_string: &str,
    region_index_dict: std::collections::HashMap<String, usize>,
    rels: &Bound<'py, PyDict>,
    n_regions: usize,
    exclude_deletions: Vec<usize>,
) -> PyResult<Bound<'py, PyTuple>> {
    if state_string.len() != 128
        || !state_string.bytes().all(|c| c == b'X' || c == b'I')
        || region_string.len() != 128
        || !region_string.is_ascii()
        || n_regions == 0
        || n_regions > 128
    {
        return Err(assertion_err(
            "Region and state maps must contain 128 valid states",
        ));
    }
    let regions = region_string
        .chars()
        .map(|c| {
            region_index_dict
                .get(&c.to_string())
                .copied()
                .filter(|&i| i < n_regions)
                .ok_or_else(|| assertion_err("Region index is missing or out of range"))
        })
        .collect::<PyResult<Vec<_>>>()?;
    let mut offsets = (0..n_regions)
        .map(|i| {
            rels.get_item(i)?
                .ok_or_else(|| assertion_err("Missing relative numbering offset"))?
                .extract::<i32>()
        })
        .collect::<PyResult<Vec<_>>>()?;
    if offsets.iter().any(|r| r.abs_diff(0) > 10000)
        || exclude_deletions.iter().any(|&r| r >= n_regions)
    {
        return Err(assertion_err(
            "Numbering offset or excluded region is out of range",
        ));
    }
    let states = states_from_py(state_vector, Some(sequence))?;
    let (numbered, start, end) = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        anarci_core::schemes::number_regions(
            sequence.as_bytes(),
            &states,
            state_string.as_bytes(),
            &regions,
            &mut offsets,
            n_regions,
            &exclude_deletions,
        )
    }))
    .map_err(|_| assertion_err("State vector exceeds region numbering limits"))?
    .map_err(core_err)?;
    let out = PyList::empty(py);
    for residues in numbered {
        let n = numbering_to_py(
            py,
            &Numbered {
                residues,
                start,
                end,
            },
        )?;
        out.append(n.get_item(0)?)?;
    }
    for (i, offset) in offsets.into_iter().enumerate() {
        rels.set_item(i, offset)?;
    }
    PyTuple::new(
        py,
        [out.into_any(), opt_usize(py, start)?, opt_usize(py, end)?],
    )
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(anarci, m)?)?;
    m.add_function(wrap_pyfunction!(run_anarci, m)?)?;
    m.add_function(wrap_pyfunction!(number, m)?)?;
    m.add_function(wrap_pyfunction!(_number_alignment, m)?)?;
    m.add_function(wrap_pyfunction!(_smooth_insertions, m)?)?;
    m.add_function(wrap_pyfunction!(_get_identity, m)?)?;
    m.add_function(wrap_pyfunction!(_get_hmm_length, m)?)?;
    m.add_function(wrap_pyfunction!(_germline_assignment, m)?)?;
    m.add_function(wrap_pyfunction!(_run_hmmer, m)?)?;
    m.add_function(wrap_pyfunction!(_germlines_json, m)?)?;
    m.add_function(wrap_pyfunction!(_blosum62, m)?)?;
    m.add_function(wrap_pyfunction!(_get_imgt_cdr, m)?)?;
    m.add_function(wrap_pyfunction!(_get_cdr3_annotations, m)?)?;
    m.add_function(wrap_pyfunction!(_get_wolfguy_l1, m)?)?;
    m.add_function(wrap_pyfunction!(_check_for_j, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_hsps, m)?)?;
    m.add_function(wrap_pyfunction!(_hmm_states, m)?)?;
    m.add_function(wrap_pyfunction!(_number_regions, m)?)?;
    features::register(m)?;
    experimental::register(m)?;
    Ok(())
}
