//! Validated Python access to the existing research codec and data utilities.
use anarci_core::{State, StateType};
use anarci_learned::{augment, data, labels, Label, LabeledRes, LabeledSeq};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeSet;

#[cfg(feature = "learned-model")]
#[path = "model.rs"]
mod model;

type StateRow = (u8, String, Option<usize>);
type Segment = (String, Vec<usize>);
type Record = (String, Vec<StateRow>);
type ExampleRow = (Segment, String, Vec<usize>, String);

fn invalid(message: &str) -> PyErr {
    PyValueError::new_err(message.to_owned())
}

fn sequence(s: &str) -> PyResult<&[u8]> {
    if !s.is_ascii() {
        return Err(invalid("sequences must contain ASCII characters"));
    }
    Ok(s.as_bytes())
}

fn label(i: usize) -> PyResult<Label> {
    if i >= anarci_learned::NUM_LABELS {
        return Err(invalid("label index must be in 0..261"));
    }
    Ok(anarci_learned::index_to_label(i))
}

fn parse_labels(values: Vec<usize>) -> PyResult<Vec<Label>> {
    values.into_iter().map(label).collect()
}

fn state_rows(states: Vec<State>) -> Vec<StateRow> {
    states
        .into_iter()
        .map(|s| (s.id, s.typ.as_char().to_ascii_uppercase().to_string(), s.si))
        .collect()
}

fn states(rows: Vec<StateRow>, length: Option<usize>) -> PyResult<Vec<State>> {
    let mut matches = BTreeSet::new();
    let mut previous = None;
    rows.into_iter()
        .map(|(id, kind, si)| {
            if !(1..=128).contains(&id) {
                return Err(invalid("state positions must be in 1..=128"));
            }
            let typ = match kind.as_str() {
                "M" => StateType::M,
                "I" => StateType::I,
                "D" => StateType::D,
                _ => return Err(invalid("state type must be M, I, or D")),
            };
            if (typ == StateType::D) != si.is_none() {
                return Err(invalid(
                    "M/I states require an index; D states require None",
                ));
            }
            if typ == StateType::M && !matches.insert(id) {
                return Err(invalid("a domain must have unique match positions"));
            }
            if let Some(index) = si {
                if length.is_some_and(|n| index >= n) {
                    return Err(invalid("state index is outside the sequence"));
                }
                if previous.is_some_and(|p| index <= p) {
                    return Err(invalid("state residue indices must strictly increase"));
                }
                previous = Some(index);
            }
            Ok(State { id, typ, si })
        })
        .collect()
}

fn segment((seq, ids): Segment) -> PyResult<LabeledSeq> {
    let bytes = sequence(&seq)?.to_vec();
    if bytes.len() != ids.len() {
        return Err(invalid("sequence and labels must have equal lengths"));
    }
    Ok(LabeledSeq {
        seq: bytes,
        labels: parse_labels(ids)?,
    })
}

fn segment_row(s: LabeledSeq) -> Segment {
    // All constructors accept ASCII sequences or construct ASCII linkers.
    (
        String::from_utf8(s.seq).expect("validated ASCII sequence"),
        s.labels
            .into_iter()
            .map(anarci_learned::label_to_index)
            .collect(),
    )
}

fn record((seq, sv): Record) -> PyResult<data::DomainRec> {
    let bytes = sequence(&seq)?.to_vec();
    let sv = states(sv, Some(bytes.len()))?;
    Ok(data::DomainRec { seq: bytes, sv })
}

fn records(rows: Vec<Record>) -> PyResult<Vec<data::DomainRec>> {
    rows.into_iter().map(record).collect()
}

fn split_name(split: data::Split) -> String {
    match split {
        data::Split::Train => "train",
        data::Split::Val => "val",
        data::Split::Test => "test",
    }
    .into()
}

fn window(lo: u8, hi: u8) -> PyResult<()> {
    if lo == 0 || hi > 128 || lo > hi {
        return Err(invalid("position windows require 1 <= lo <= hi <= 128"));
    }
    Ok(())
}

#[pyfunction]
fn exp_encode_domain(rows: Vec<StateRow>) -> PyResult<Vec<(usize, usize)>> {
    Ok(labels::encode_domain(&states(rows, None)?)
        .into_iter()
        .map(|r| (r.si, anarci_learned::label_to_index(r.label)))
        .collect())
}

#[pyfunction]
fn exp_decode_domain(rows: Vec<(usize, usize)>) -> PyResult<Vec<StateRow>> {
    let mut previous = None;
    let mut last_match = None;
    let residues = rows
        .into_iter()
        .map(|(si, index)| {
            if previous.is_some_and(|p| si <= p) {
                return Err(invalid("labeled residue indices must strictly increase"));
            }
            previous = Some(si);
            let label = label(index)?;
            if let Label::Match(p) = label {
                if last_match.is_some_and(|q| p <= q) {
                    return Err(invalid("domain match positions must strictly increase"));
                }
                last_match = Some(p);
            }
            Ok(LabeledRes { si, label })
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(state_rows(labels::decode_domain(&residues)))
}

#[pyfunction]
fn exp_strip_boundary_deletes(rows: Vec<StateRow>) -> PyResult<Vec<StateRow>> {
    Ok(state_rows(labels::strip_boundary_deletes(&states(
        rows, None,
    )?)))
}

#[pyfunction]
fn exp_decode_sequence(values: Vec<usize>) -> PyResult<Vec<Vec<StateRow>>> {
    Ok(labels::decode_sequence(&parse_labels(values)?)
        .into_iter()
        .map(state_rows)
        .collect())
}

#[pyfunction]
fn exp_label_to_index(index: usize) -> PyResult<usize> {
    Ok(anarci_learned::label_to_index(label(index)?))
}

#[pyfunction]
fn exp_index_to_label(index: usize) -> PyResult<usize> {
    exp_label_to_index(index)
}

#[pyfunction]
fn exp_viterbi_decode(py: Python<'_>, logits: Vec<Vec<f32>>) -> PyResult<Vec<usize>> {
    if logits
        .iter()
        .any(|row| row.len() != anarci_learned::NUM_LABELS)
    {
        return Err(invalid("each logits row must contain exactly 261 scores"));
    }
    if logits.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid("logits must contain finite scores"));
    }
    Ok(py
        .allow_threads(|| anarci_learned::viterbi_decode(&logits))
        .into_iter()
        .map(anarci_learned::label_to_index)
        .collect())
}

#[pyfunction]
fn exp_from_domain(seq: &str, rows: Vec<StateRow>) -> PyResult<Segment> {
    let bytes = sequence(seq)?;
    Ok(segment_row(LabeledSeq::from_domain(
        bytes,
        &states(rows, Some(bytes.len()))?,
    )))
}

#[pyfunction]
fn exp_region(seq: &str, index: usize) -> PyResult<Segment> {
    if index < 256 {
        return Err(invalid("region requires a non-domain label"));
    }
    Ok(segment_row(LabeledSeq::region(
        sequence(seq)?,
        label(index)?,
    )))
}

#[pyfunction]
fn exp_junk(seq: &str) -> PyResult<Segment> {
    Ok(segment_row(LabeledSeq::junk(sequence(seq)?)))
}

#[pyfunction]
fn exp_restrict_positions(value: Segment, lo: u8, hi: u8) -> PyResult<Segment> {
    window(lo, hi)?;
    Ok(segment_row(segment(value)?.restrict_positions(lo, hi)))
}

#[pyfunction]
fn exp_restrict_window(value: Segment, start: usize, end: usize) -> PyResult<Segment> {
    Ok(segment_row(segment(value)?.restrict_window(start, end)))
}

#[pyfunction]
fn exp_domains(value: Segment) -> PyResult<Vec<Vec<StateRow>>> {
    Ok(segment(value)?
        .domains()
        .into_iter()
        .map(state_rows)
        .collect())
}

fn concatenate(values: Vec<Segment>, separate: bool) -> PyResult<Segment> {
    let parts = values
        .into_iter()
        .map(segment)
        .collect::<PyResult<Vec<_>>>()?;
    let refs: Vec<_> = parts.iter().collect();
    Ok(segment_row(if separate {
        augment::concat_domains(&refs)
    } else {
        augment::concat(&refs)
    }))
}

#[pyfunction]
fn exp_concat(values: Vec<Segment>) -> PyResult<Segment> {
    concatenate(values, false)
}

#[pyfunction]
fn exp_concat_domains(values: Vec<Segment>) -> PyResult<Segment> {
    concatenate(values, true)
}

#[pyfunction]
fn exp_gly_ser_linker(repeats: usize) -> PyResult<Segment> {
    let size = repeats
        .checked_mul(5)
        .ok_or_else(|| invalid("linker size overflow"))?;
    // Check representable allocation size before the core's repeat multiplication.
    let mut probe = Vec::<u8>::new();
    probe
        .try_reserve_exact(size)
        .map_err(|_| invalid("linker is too large"))?;
    drop(probe);
    Ok(segment_row(augment::gly_ser_linker(repeats)))
}

#[pyfunction]
fn exp_cdr3_of(value: Record) -> PyResult<String> {
    Ok(String::from_utf8(data::cdr3_of(&record(value)?)).expect("validated ASCII sequence"))
}

#[pyfunction]
fn exp_overlap_identity(py: Python<'_>, a: &str, b: &str) -> PyResult<f32> {
    let (a, b) = (sequence(a)?, sequence(b)?);
    Ok(py.allow_threads(|| data::overlap_identity(a, b)))
}

#[pyfunction]
fn exp_cluster_records(
    py: Python<'_>,
    values: Vec<Record>,
    full_thr: f32,
    cdr3_thr: f32,
) -> PyResult<Vec<usize>> {
    if !full_thr.is_finite()
        || !cdr3_thr.is_finite()
        || !(0.0..=1.0).contains(&full_thr)
        || !(0.0..=1.0).contains(&cdr3_thr)
    {
        return Err(invalid("identity thresholds must be finite and in 0..=1"));
    }
    let values = records(values)?;
    Ok(py.allow_threads(|| data::cluster_records(&values, full_thr, cdr3_thr)))
}

#[pyfunction]
fn exp_split_for_rep(rep: &str) -> PyResult<String> {
    Ok(split_name(data::split_for_rep(sequence(rep)?)))
}

#[pyfunction]
fn exp_splits_of(py: Python<'_>, values: Vec<Record>) -> PyResult<Vec<String>> {
    let values = records(values)?;
    Ok(py
        .allow_threads(|| data::splits_of(&values))
        .into_iter()
        .map(split_name)
        .collect())
}

#[pyfunction]
fn exp_near_duplicate(py: Python<'_>, a: Record, b: Record) -> PyResult<bool> {
    let (a, b) = (record(a)?, record(b)?);
    Ok(py.allow_threads(|| data::near_duplicate(&a, &b)))
}

#[pyfunction]
fn exp_build_dataset(
    py: Python<'_>,
    values: Vec<Record>,
    partial_windows: Vec<(u8, u8)>,
) -> PyResult<Vec<ExampleRow>> {
    for &(lo, hi) in &partial_windows {
        window(lo, hi)?;
    }
    let values = records(values)?;
    Ok(py
        .allow_threads(|| data::build_dataset(&values, &partial_windows))
        .into_iter()
        .map(|e| {
            (
                segment_row(e.data),
                split_name(e.split),
                e.clusters,
                e.kind.into(),
            )
        })
        .collect())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    macro_rules! expose { ($($f:ident),+ $(,)?) => { $(m.add_function(wrap_pyfunction!($f, m)?)?;)+ }; }
    expose!(
        exp_encode_domain,
        exp_decode_domain,
        exp_strip_boundary_deletes,
        exp_decode_sequence,
        exp_label_to_index,
        exp_index_to_label,
        exp_viterbi_decode,
        exp_from_domain,
        exp_region,
        exp_junk,
        exp_restrict_positions,
        exp_restrict_window,
        exp_domains,
        exp_concat,
        exp_concat_domains,
        exp_gly_ser_linker,
        exp_cdr3_of,
        exp_overlap_identity,
        exp_cluster_records,
        exp_split_for_rep,
        exp_splits_of,
        exp_near_duplicate,
        exp_build_dataset
    );
    #[cfg(feature = "learned-model")]
    model::register(m)?;
    Ok(())
}
