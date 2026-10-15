//! Synthetic query tests; no biological reference inputs are used.
use anarci_core::orchestrate::validate_scan_input;
use anarci_core::{
    anarci, default_allow, number, run_anarci, CoreError, GermlineMethod, HmmEngine, Hsp,
};
use std::sync::atomic::{AtomicUsize, Ordering};

struct Engine {
    calls: AtomicUsize,
    fail: bool,
}
impl HmmEngine for Engine {
    fn scan_one(&self, _: &str, _: &[u8]) -> Vec<Hsp> {
        panic!("checked API must be used")
    }
    fn scan_one_checked(&self, _: &str, _: &[u8]) -> anarci_core::CResult<Vec<Hsp>> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        if self.fail {
            Err(CoreError::Engine("synthetic scan failure".into()))
        } else {
            Ok(vec![])
        }
    }
}

#[test]
fn batch_validation_precedes_deduplication_and_scanning() {
    let engine = Engine {
        calls: AtomicUsize::new(0),
        fail: false,
    };
    let rows = vec![
        ("ok".into(), b"AAA".to_vec()),
        ("bad\0name".into(), b"AAA".to_vec()),
    ];
    assert!(run_anarci(
        &engine,
        &rows,
        "imgt",
        &default_allow(),
        false,
        None,
        80.0,
        false,
        GermlineMethod::Identity
    )
    .is_err());
    assert!(anarci(
        &engine,
        &rows,
        "imgt",
        &default_allow(),
        false,
        None,
        80.0,
        false,
        GermlineMethod::Identity
    )
    .is_err());
    assert_eq!(engine.calls.load(Ordering::Relaxed), 0);
}

#[test]
fn hmm_alphabet_accepts_ambiguity_but_number_remains_strict() {
    let engine = Engine {
        calls: AtomicUsize::new(0),
        fail: false,
    };
    let rows = vec![("synthetic".into(), b"ABJZOUXabjzoux*".to_vec())];
    assert!(run_anarci(
        &engine,
        &rows,
        "imgt",
        &default_allow(),
        false,
        None,
        80.0,
        false,
        GermlineMethod::Identity
    )
    .is_ok());
    assert!(number(
        &engine,
        b"AXA",
        "imgt",
        &default_allow(),
        None,
        false,
        0,
        80.0
    )
    .is_err());
    assert!(validate_scan_input("query", &vec![b'A'; 100_000]).is_ok());
    assert!(validate_scan_input("query", &vec![b'A'; 100_001]).is_err());
    for seq in [b"A\0A".as_slice(), b"A A", b"A-A", &[0xff]] {
        assert!(validate_scan_input("query", seq).is_err());
    }
}

#[test]
fn invalid_numeric_settings_and_native_failures_propagate() {
    let engine = Engine {
        calls: AtomicUsize::new(0),
        fail: true,
    };
    for score in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(run_anarci(
            &engine,
            &[],
            "imgt",
            &default_allow(),
            false,
            None,
            score,
            false,
            GermlineMethod::Identity
        )
        .is_err());
        assert!(number(
            &engine,
            b"AAA",
            "imgt",
            &default_allow(),
            None,
            false,
            70,
            score
        )
        .is_err());
    }
    let rows = vec![("synthetic".into(), b"AAA".to_vec())];
    let error = run_anarci(
        &engine,
        &rows,
        "imgt",
        &default_allow(),
        false,
        None,
        80.0,
        false,
        GermlineMethod::Identity,
    )
    .unwrap_err();
    assert_eq!(error.to_string(), "synthetic scan failure");
}

#[test]
fn number_propagates_native_engine_failures() {
    let engine = Engine {
        calls: AtomicUsize::new(0),
        fail: true,
    };
    let error = number(
        &engine,
        b"AAA",
        "imgt",
        &default_allow(),
        None,
        false,
        0,
        80.0,
    )
    .unwrap_err();
    assert!(matches!(error, CoreError::Engine(_)));
}
