//! Python access to the shared CLI library. Native compute releases the GIL.
use anarci_cli::seqio::Record;
use anarci_cli::{Args, Eng, Fmt};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn error(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn python_json(py: Python<'_>, value: serde_json::Value) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("json")?
        .call_method1("loads", (value.to_string(),))?
        .unbind())
}
fn records(sequences: Vec<(String, String)>) -> Vec<Record> {
    sequences
        .into_iter()
        .map(|(id, seq)| Record {
            id,
            desc: String::new(),
            seq: seq.into_bytes(),
        })
        .collect()
}
#[allow(clippy::too_many_arguments)]
fn config(
    scheme: &str,
    database: &str,
    assign_germline: bool,
    allowed_species: Option<Vec<String>>,
    allow: Option<&Bound<'_, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<String>,
    ncpu: usize,
) -> PyResult<Args> {
    let engine = match database.to_ascii_lowercase().as_str() {
        "pan" => Eng::Pan,
        "all" => Eng::All,
        _ => return Err(error("database must be 'pan' or 'ALL'")),
    };
    let allow = allow
        .map(|value| {
            value
                .try_iter()?
                .map(|item| item?.extract::<String>())
                .collect::<PyResult<std::collections::BTreeSet<_>>>()
        })
        .transpose()?;
    let args = Args {
        scheme: scheme.into(),
        engine,
        assign_germline,
        species: allowed_species
            .filter(|s| !s.is_empty())
            .map(|s| s.join(","))
            .unwrap_or_else(|| "all".into()),
        allow,
        bit_score_threshold,
        germline_method,
        ncpu,
        ..Args::default()
    };
    args.validate().map_err(error)?;
    Ok(args)
}

#[pyfunction]
fn read_sequences(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let parsed = py
        .allow_threads(|| anarci_cli::parse_records(text.as_bytes()))
        .map_err(error)?;
    let values: Vec<_> = parsed.into_iter().map(|r| serde_json::json!({ "id":r.id, "description":r.desc, "sequence":String::from_utf8_lossy(&r.seq) })).collect();
    python_json(py, serde_json::json!(values))
}

#[pyfunction]
#[pyo3(signature=(dna, frame=0))]
fn translate_dna(py: Python<'_>, dna: &str, frame: usize) -> PyResult<String> {
    if frame >= 6 {
        return Err(error("frame must be between 0 and 5"));
    }
    if dna
        .bytes()
        .any(|b| !b"ACGTUNRYKMSWBDHV".contains(&b.to_ascii_uppercase()))
    {
        return Err(error("DNA contains an invalid nucleotide character"));
    }
    let translated =
        py.allow_threads(|| anarci_cli::translate::translate_frame(dna.as_bytes(), frame));
    String::from_utf8(translated).map_err(error)
}

#[pyfunction]
#[pyo3(signature=(text, format="pdb"))]
fn parse_structure(py: Python<'_>, text: &str, format: &str) -> PyResult<Py<PyAny>> {
    let chains = py
        .allow_threads(|| anarci_cli::parse_structure(text.as_bytes(), format))
        .map_err(error)?;
    let result: Vec<_> = chains.iter().map(|c| serde_json::json!({
        "id": c.id, "sequence": String::from_utf8_lossy(&c.seq()),
        "residues": c.residues.iter().map(|r| serde_json::json!({"aa": (r.aa as char).to_string(), "ca": r.ca})).collect::<Vec<_>>()
    })).collect();
    python_json(py, serde_json::json!(result))
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature=(sequences, scheme="imgt", database="pan", assign_germline=false, allowed_species=None, allow=None, bit_score_threshold=80.0, germline_method=None, ncpu=0))]
fn annotate_protein(
    py: Python<'_>,
    sequences: Vec<(String, String)>,
    scheme: &str,
    database: &str,
    assign_germline: bool,
    allowed_species: Option<Vec<String>>,
    allow: Option<&Bound<'_, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<String>,
    ncpu: usize,
) -> PyResult<Py<PyAny>> {
    let args = config(
        scheme,
        database,
        assign_germline,
        allowed_species,
        allow,
        bit_score_threshold,
        germline_method,
        ncpu,
    )?;
    let input = records(sequences);
    let value = py
        .allow_threads(|| {
            anarci_cli::annotate_protein(&input, &args)
                .and_then(|r| serde_json::to_value(r).map_err(|e| e.to_string()))
        })
        .map_err(error)?;
    python_json(py, value)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature=(sequences, scheme="imgt", database="pan", assign_germline=false, allowed_species=None, allow=None, bit_score_threshold=80.0, germline_method=None, ncpu=0))]
fn annotate_dna(
    py: Python<'_>,
    sequences: Vec<(String, String)>,
    scheme: &str,
    database: &str,
    assign_germline: bool,
    allowed_species: Option<Vec<String>>,
    allow: Option<&Bound<'_, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<String>,
    ncpu: usize,
) -> PyResult<Py<PyAny>> {
    let args = config(
        scheme,
        database,
        assign_germline,
        allowed_species,
        allow,
        bit_score_threshold,
        germline_method,
        ncpu,
    )?;
    let input = records(sequences);
    let value = py
        .allow_threads(|| {
            anarci_cli::annotate_dna(&input, &args)
                .and_then(|r| serde_json::to_value(r).map_err(|e| e.to_string()))
        })
        .map_err(error)?;
    python_json(py, value)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature=(text, format="pdb", scheme="imgt", database="pan", assign_germline=false, allowed_species=None, allow=None, bit_score_threshold=80.0, germline_method=None, ncpu=0))]
fn pair_structure(
    py: Python<'_>,
    text: &str,
    format: &str,
    scheme: &str,
    database: &str,
    assign_germline: bool,
    allowed_species: Option<Vec<String>>,
    allow: Option<&Bound<'_, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<String>,
    ncpu: usize,
) -> PyResult<Py<PyAny>> {
    let args = config(
        scheme,
        database,
        assign_germline,
        allowed_species,
        allow,
        bit_score_threshold,
        germline_method,
        ncpu,
    )?;
    let value = py
        .allow_threads(|| {
            anarci_cli::pair_structure(text.as_bytes(), format, &args)
                .and_then(|r| serde_json::to_value(r).map_err(|e| e.to_string()))
        })
        .map_err(error)?;
    python_json(py, value)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature=(sequences, output_format="json", input_type="protein", scheme="imgt", database="pan", assign_germline=false, allowed_species=None, allow=None, bit_score_threshold=80.0, germline_method=None, ncpu=0))]
fn format_annotations(
    py: Python<'_>,
    sequences: Vec<(String, String)>,
    output_format: &str,
    input_type: &str,
    scheme: &str,
    database: &str,
    assign_germline: bool,
    allowed_species: Option<Vec<String>>,
    allow: Option<&Bound<'_, PyAny>>,
    bit_score_threshold: f64,
    germline_method: Option<String>,
    ncpu: usize,
) -> PyResult<String> {
    let args = config(
        scheme,
        database,
        assign_germline,
        allowed_species,
        allow,
        bit_score_threshold,
        germline_method,
        ncpu,
    )?;
    let fmt = match output_format.to_ascii_lowercase().as_str() {
        "json" => Fmt::Json,
        "tsv" => Fmt::Tsv,
        "anarci" => Fmt::Anarci,
        _ => return Err(error("output_format must be json, tsv, or anarci")),
    };
    if !matches!(input_type, "protein" | "dna") {
        return Err(error("input_type must be protein or dna"));
    }
    let input = records(sequences);
    py.allow_threads(|| {
        let domains = if input_type == "dna" {
            anarci_cli::annotate_dna(&input, &args)?.domains
        } else {
            anarci_cli::annotate_protein(&input, &args)?
        };
        anarci_cli::render_output(&domains, fmt)
    })
    .map_err(error)
}

#[pyfunction]
fn cli(py: Python<'_>, args: Vec<String>) -> PyResult<i32> {
    py.allow_threads(|| anarci_cli::cli(args)).map_err(error)?;
    Ok(0)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(translate_dna, m)?)?;
    m.add_function(wrap_pyfunction!(annotate_protein, m)?)?;
    m.add_function(wrap_pyfunction!(annotate_dna, m)?)?;
    m.add_function(wrap_pyfunction!(parse_structure, m)?)?;
    m.add_function(wrap_pyfunction!(pair_structure, m)?)?;
    m.add_function(wrap_pyfunction!(format_annotations, m)?)?;
    m.add_function(wrap_pyfunction!(cli, m)?)?;
    Ok(())
}
