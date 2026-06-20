//! anarci-core — pure-Rust core of anarci-rs.
//!
//! Faithful port of ANARCI's numbering schemes, germline assignment, and the
//! HMM-alignment → state-vector transform. No FFI, no Python — fully unit-testable
//! against fixtures captured from reference ANARCI (conda `anarci 2024.05.21`).

pub mod align;
pub mod constants;
pub mod germlines;
pub mod orchestrate;
pub mod schemes;
pub mod types;

pub use align::{hmm_alignment_to_states, parse_hmmer_query, DomainDetails, Hsp, ParsedQuery};
pub use germlines::{get_hmm_length, get_identity, run_germline_assignment, Germline};
pub use orchestrate::{
    anarci, chain_type_to_class, default_allow, number, resolve_scheme, run_anarci, validate_sequence,
    DomainInfo, HmmEngine, SeqResult,
};
pub use schemes::number_sequence_from_alignment;
pub use types::{CResult, CoreError, Numbered, Res, State, StateType, StateVector};
