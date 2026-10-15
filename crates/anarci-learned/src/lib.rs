//! anarci-learned — learned numbering engine for anarci-rs.
//!
//! The learned engine predicts only the per-residue [`labels::Label`] stream; the
//! existing byte-exact `anarci-core` code turns the decoded state vectors into numbering
//! for every scheme. See `LEARNED_ENGINE_PLAN.md`.
//!
//! This crate is pure Rust and ML-runtime-free by default; the model (candle) lives
//! behind a Cargo feature so the default wheel stays lean.

pub mod augment;
pub mod data;
pub mod decode;
pub mod labels;
#[cfg(feature = "model")]
pub mod model;

pub use augment::{concat, concat_domains, gly_ser_linker, LabeledSeq};
pub use decode::{index_to_label, label_to_index, viterbi_decode, NUM_LABELS};
pub use labels::{
    decode_domain, decode_sequence, encode_domain, strip_boundary_deletes, Label, LabeledRes,
};
