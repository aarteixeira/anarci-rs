//! Embedded HMM engines (no runtime data dependency — a self-contained binary).
//! PAN (`FEW.hmm`, 7 pan-species profiles) is the fast default for NGS; ALL
//! (`ALL.hmm`, every species) is the byte-exact ANARCI-parity engine.

use anarci_core::align::Hsp;
use anarci_core::orchestrate::HmmEngine;
use anarci_hmm::Engine;

/// Bridges the HMMER [`Engine`] to the core's [`HmmEngine`] trait (same shape the
/// Python layer uses): `scan_one` delegates to `Engine::scan`.
pub struct Adapter<'a>(pub &'a Engine);

impl HmmEngine for Adapter<'_> {
    fn scan_one(&self, name: &str, seq: &[u8]) -> Vec<Hsp> {
        self.0.scan(name, seq)
    }
    fn scan_one_checked(&self, name: &str, seq: &[u8]) -> anarci_core::CResult<Vec<Hsp>> {
        self.0
            .scan_checked(name, seq)
            .map_err(|e| anarci_core::CoreError::Engine(format!("HMM scan failed: {e}")))
    }
}

const PAN_HMM: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference_data/dat/HMMs/FEW.hmm"
));
const ALL_HMM: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../reference_data/dat/HMMs/ALL.hmm"
));

/// Which bundled profile set to use.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum EngineKind {
    /// FEW.hmm — pan-species, ~4–12× faster. Default (right for NGS throughput).
    Pan,
    /// ALL.hmm — every species, byte-for-byte ANARCI parity.
    All,
}

/// Load a bundled profile set through the concurrency-safe native loader.
pub fn load(kind: EngineKind) -> Result<&'static Engine, String> {
    static PAN: once_cell::sync::OnceCell<Engine> = once_cell::sync::OnceCell::new();
    static ALL: once_cell::sync::OnceCell<Engine> = once_cell::sync::OnceCell::new();
    let (cell, data) = match kind {
        EngineKind::Pan => (&PAN, PAN_HMM),
        EngineKind::All => (&ALL, ALL_HMM),
    };
    cell.get_or_try_init(|| {
        Engine::load_embedded(data).map_err(|e| format!("load embedded HMM: {e}"))
    })
}
