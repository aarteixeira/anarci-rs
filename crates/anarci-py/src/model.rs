//! CPU bindings for the optional, untrained Candle research model.
use anarci_learned::{model, NUM_LABELS};
use candle_core::{DType, Device, Tensor};
use candle_nn::{VarBuilder, VarMap};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

const MAX_LEN: usize = 1024;

fn value_error(message: &str) -> PyErr {
    PyValueError::new_err(message.to_owned())
}

fn model_error(error: candle_core::Error) -> PyErr {
    PyRuntimeError::new_err(format!("research model failed: {error}"))
}

fn dimensions(dim: usize, layers: usize, ff: usize) -> PyResult<()> {
    if dim == 0 || ff == 0 {
        return Err(value_error("dim and ff must be positive"));
    }
    let count = dim
        .checked_mul(dim)
        .and_then(|n| n.checked_mul(4))
        .and_then(|n| {
            dim.checked_mul(ff)
                .and_then(|f| f.checked_mul(2))
                .and_then(|f| n.checked_add(f))
        })
        .and_then(|n| n.checked_mul(layers))
        .and_then(|n| {
            dim.checked_mul(MAX_LEN + model::VOCAB + NUM_LABELS)
                .and_then(|p| n.checked_add(p))
        });
    if count.map_or(true, |n| {
        n > isize::MAX as usize / std::mem::size_of::<f32>()
    }) {
        return Err(value_error(
            "model dimensions exceed the supported allocation size",
        ));
    }
    Ok(())
}

fn input_sequence(sequence: &str) -> PyResult<&[u8]> {
    if !sequence.is_ascii() || sequence.is_empty() || sequence.len() > MAX_LEN {
        return Err(value_error(
            "model sequences must contain 1..1024 ASCII characters",
        ));
    }
    Ok(sequence.as_bytes())
}

#[pyclass(name = "ExperimentalModel")]
pub struct PyModel {
    inner: model::Model,
}

#[pymethods]
impl PyModel {
    #[new]
    #[pyo3(signature = (dim=64, layers=2, ff=128))]
    fn new(py: Python<'_>, dim: usize, layers: usize, ff: usize) -> PyResult<Self> {
        dimensions(dim, layers, ff)?;
        let inner = py
            .allow_threads(|| {
                let variables = VarMap::new();
                let vb = VarBuilder::from_varmap(&variables, DType::F32, &Device::Cpu);
                model::Model::new(dim, layers, ff, vb)
            })
            .map_err(model_error)?;
        Ok(Self { inner })
    }

    fn logits(&self, py: Python<'_>, sequence: &str) -> PyResult<Vec<Vec<f32>>> {
        let sequence = input_sequence(sequence)?;
        py.allow_threads(|| model::logits_rows(&self.inner, sequence, &Device::Cpu))
            .map_err(model_error)
    }

    fn forward(&self, py: Python<'_>, tokens: Vec<Vec<u32>>) -> PyResult<Vec<Vec<Vec<f32>>>> {
        let length = tokens.first().map_or(0, Vec::len);
        if tokens.is_empty()
            || length == 0
            || length > MAX_LEN
            || tokens
                .iter()
                .any(|row| row.len() != length || row.iter().any(|&t| t >= model::VOCAB as u32))
        {
            return Err(value_error("tokens must be a nonempty rectangular batch with 1..1024 columns and values in 0..22"));
        }
        py.allow_threads(|| {
            let input = Tensor::new(tokens, &Device::Cpu)?;
            self.inner.forward(&input)?.to_vec3::<f32>()
        })
        .map_err(model_error)
    }
}

#[pyfunction]
fn exp_tokenize(sequence: &str) -> PyResult<Vec<u32>> {
    if !sequence.is_ascii() {
        return Err(value_error("sequences must contain ASCII characters"));
    }
    Ok(model::tokenize(sequence.as_bytes()))
}

#[pyfunction]
#[pyo3(signature = (examples, dim=64, layers=2, ff=128, steps=1, lr=0.001))]
fn exp_overfit(
    py: Python<'_>,
    examples: Vec<(String, Vec<usize>)>,
    dim: usize,
    layers: usize,
    ff: usize,
    steps: usize,
    lr: f64,
) -> PyResult<(PyModel, f32)> {
    dimensions(dim, layers, ff)?;
    if examples.is_empty() {
        return Err(value_error("overfit requires at least one labeled example"));
    }
    if !lr.is_finite() || lr <= 0.0 {
        return Err(value_error("learning rate must be finite and positive"));
    }
    let examples = examples
        .into_iter()
        .map(|(seq, targets)| {
            let sequence = input_sequence(&seq)?;
            if targets.len() != sequence.len() || targets.iter().any(|&i| i >= NUM_LABELS) {
                return Err(value_error(
                    "each example needs one label index in 0..261 per residue",
                ));
            }
            Ok((
                model::tokenize(sequence),
                targets.into_iter().map(|i| i as u32).collect(),
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    let (inner, accuracy) = py
        .allow_threads(|| model::overfit(&examples, dim, layers, ff, steps, lr, &Device::Cpu))
        .map_err(model_error)?;
    Ok((PyModel { inner }, accuracy))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyModel>()?;
    m.add_function(wrap_pyfunction!(exp_tokenize, m)?)?;
    m.add_function(wrap_pyfunction!(exp_overfit, m)?)?;
    Ok(())
}
