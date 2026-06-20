//! `anarci_rs` — drop-in Python API for ANARCI with all-Rust internals.
//!
//! Mirrors ANARCI's `anarci`, `run_anarci`, `number` exactly (names, signatures,
//! return shapes). The HMM scan runs in-process (native HMMER 3.4); batch runs in
//! Rust with rayon. The GIL is released around compute.

use anarci_core::orchestrate::{self, DomainInfo, SeqResult};
use anarci_core::{Germline, GermlineMethod, HmmEngine, Hsp, Numbered};
use anarci_hmm::Engine;
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

static ALL_HMM: &[u8] =
    include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../reference_data/dat/HMMs/ALL.hmm"));
static PAN_HMM: &[u8] =
    include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../reference_data/dat/HMMs/FEW.hmm"));

fn load_engine(cell: &'static OnceCell<Engine>, data: &[u8], fname: &str) -> PyResult<&'static Engine> {
    cell.get_or_try_init(|| {
        let dir = std::env::temp_dir().join("anarci_rs_data");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join(fname);
        if std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0) != data.len() as u64 {
            std::fs::write(&path, data)?;
        }
        Engine::load(&path).map_err(|e| {
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

/// `database` -> pan(default)/exact. Unknown values error (no silent fallback).
fn parse_database(database: &str) -> PyResult<bool> {
    match database {
        "pan" | "PAN" => Ok(true),
        "ALL" => Ok(false),
        other => Err(assertion_err(format!(
            "Unknown database '{other}'; use 'pan' (default, fast) or 'ALL' (exact ANARCI)."
        ))),
    }
}

/// Resolve the germline-assignment method.
///
/// `None` defaults to identity germline assignment (ANARCI-compatible and fast — the
/// pan engine derives a species label from it on every call, so the ~10x-slower e-value
/// path must NOT be the default or it erases the pan speedup). Pass `germline_method="evalue"`
/// for the higher-accuracy RIOT-style V/J-gene assignment (worth it when you need accurate
/// germline genes, e.g. with `assign_germline=True`).
fn parse_germline_method(method: Option<&str>, _pan: bool) -> PyResult<GermlineMethod> {
    match method {
        None => Ok(GermlineMethod::Identity),
        Some(s) => GermlineMethod::parse(s).map_err(assertion_err),
    }
}

/// Adapter so the orchestration layer can use the FFI engine via the `HmmEngine` trait.
struct Adapter<'a>(&'a Engine);
impl HmmEngine for Adapter<'_> {
    fn scan_one(&self, name: &str, seq: &[u8]) -> Vec<Hsp> {
        self.0.scan(name, seq)
    }
}

// ---- helpers --------------------------------------------------------------

fn assertion_err(msg: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyAssertionError::new_err(msg.to_string())
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
            let pos = PyTuple::new(py, [r.num.into_pyobject(py)?.into_any(), r.ins.into_pyobject(py)?.into_any()])?;
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
            Some((sp, ge)) => PyTuple::new(py, [sp.into_pyobject(py)?.into_any(), ge.into_pyobject(py)?.into_any()])?.into_any(),
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

fn detail_to_py<'py>(py: Python<'py>, d: &DomainInfo) -> PyResult<Bound<'py, PyAny>> {
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
    Ok(dict.into_any())
}

fn hit_table_to_py<'py>(py: Python<'py>, rows: &[anarci_core::align::HitRow]) -> PyResult<Bound<'py, PyList>> {
    let header = PyList::new(
        py,
        ["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"],
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
) -> PyResult<(Bound<'py, PyList>, Bound<'py, PyList>, Bound<'py, PyList>)> {
    let mut numbered: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    let mut details: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    let mut hits: Vec<Bound<'py, PyAny>> = Vec::with_capacity(results.len());
    for r in results {
        match &r.numbered {
            None => numbered.push(py.None().into_bound(py)),
            Some(doms) => {
                let v: Vec<Bound<'py, PyAny>> =
                    doms.iter().map(|n| numbering_to_py(py, n)).collect::<PyResult<_>>()?;
                numbered.push(PyList::new(py, v)?.into_any());
            }
        }
        match &r.details {
            None => details.push(py.None().into_bound(py)),
            Some(ds) => {
                let v: Vec<Bound<'py, PyAny>> =
                    ds.iter().map(|d| detail_to_py(py, d)).collect::<PyResult<_>>()?;
                details.push(PyList::new(py, v)?.into_any());
            }
        }
        hits.push(hit_table_to_py(py, &r.hit_table)?.into_any());
    }
    Ok((PyList::new(py, numbered)?, PyList::new(py, details)?, PyList::new(py, hits)?))
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
) -> PyResult<Vec<SeqResult>> {
    let eng = engine(pan)?;
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
                Err(_) => run(),
            }
        } else {
            run()
        }
    });
    res.map_err(assertion_err)
}

// ---- public API -----------------------------------------------------------

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (sequences, scheme="imgt", database="pan", output=false, outfile=None,
    csv=false, allow=None, hmmerpath="", ncpu=None, assign_germline=false,
    allowed_species=None, bit_score_threshold=80.0, germline_method=None))]
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
) -> PyResult<Bound<'py, PyTuple>> {
    let _ = (output, outfile, csv, hmmerpath); // accepted for signature parity
    let pan = parse_database(database)?;
    let gmethod = parse_germline_method(germline_method, pan)?;
    let seqs = parse_sequences(sequences)?;
    let allow = parse_allow(allow)?;
    let species = parse_species(allowed_species)?;
    let results = run_core(
        py, &seqs, scheme, &allow, assign_germline, species.as_deref(),
        bit_score_threshold, ncpu.unwrap_or(1), pan, gmethod,
    )?;
    let (n, d, h) = results_to_py(py, &results)?;
    PyTuple::new(py, [n.into_any(), d.into_any(), h.into_any()])
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (seq, scheme="imgt", database="pan", output=false, outfile=None,
    csv=false, allow=None, hmmerpath="", ncpu=1, assign_germline=false,
    allowed_species=None, bit_score_threshold=80.0, germline_method=None))]
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
) -> PyResult<Bound<'py, PyTuple>> {
    let _ = (output, outfile, csv, hmmerpath);
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
        py, &seqs, scheme, &allow, assign_germline, species.as_deref(),
        bit_score_threshold, ncpu.max(1), pan, gmethod,
    )?;
    let (n, d, h) = results_to_py(py, &results)?;
    // Sequences out: [(id, seq_str), ...]
    let seq_list: Vec<Bound<'py, PyAny>> = seqs
        .iter()
        .map(|(id, s)| {
            let t = PyTuple::new(py, [id.clone(), String::from_utf8_lossy(s).into_owned()])?;
            Ok(t.into_any())
        })
        .collect::<PyResult<_>>()?;
    let seq_out = PyList::new(py, seq_list)?;
    PyTuple::new(py, [seq_out.into_any(), n.into_any(), d.into_any(), h.into_any()])
}

#[pyfunction]
#[pyo3(signature = (sequence, scheme="imgt", database="pan", allow=None,
    allowed_species=vec!["human".into(), "mouse".into()]))]
fn number<'py>(
    py: Python<'py>,
    sequence: &str,
    scheme: &str,
    database: &str,
    allow: Option<&Bound<'py, PyAny>>,
    allowed_species: Vec<String>,
) -> PyResult<Bound<'py, PyTuple>> {
    let pan = parse_database(database)?;
    let allow = parse_allow(allow)?;
    let species = if allowed_species.is_empty() { None } else { Some(allowed_species) };
    let eng = engine(pan)?;
    let out = py.allow_threads(|| {
        anarci_core::number(&Adapter(eng), sequence.as_bytes(), scheme, &allow, species.as_deref(), pan)
    });
    match out {
        Err(_) => Ok(PyTuple::new(py, [false, false])?), // scheme/chain error -> (False, False)
        Ok(None) => Ok(PyTuple::new(py, [false, false])?),
        Ok(Some((num, class))) => {
            let numbering = numbering_to_py(py, &num)?;
            // number() returns just the numbering list (not the (num,start,end) tuple).
            let numbering_list = numbering.get_item(0)?;
            Ok(PyTuple::new(py, [numbering_list, class.into_pyobject(py)?.into_any()])?)
        }
    }
}

#[pymodule]
fn anarci_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(anarci, m)?)?;
    m.add_function(wrap_pyfunction!(run_anarci, m)?)?;
    m.add_function(wrap_pyfunction!(number, m)?)?;
    Ok(())
}
