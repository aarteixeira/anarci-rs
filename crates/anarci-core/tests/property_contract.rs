//! Bounded, deterministic contracts using synthetic sequences and a checked fake scanner.
//! No native engine, historical fixtures, or property-test dependencies are required.
use anarci_core::orchestrate::{check_for_j, resolve_scheme, validate_scan_input};
use anarci_core::{
    anarci, default_allow, number, number_sequence_from_alignment, run_anarci, CResult, CoreError,
    GermlineMethod, HmmEngine, Hsp, SeqResult, State, StateType,
};
use std::sync::Mutex;

#[derive(Default)]
struct Scanner {
    calls: Mutex<Vec<(String, Vec<u8>)>>,
    fail: bool,
    hit: bool,
}

impl HmmEngine for Scanner {
    fn scan_one(&self, _: &str, _: &[u8]) -> Vec<Hsp> {
        panic!("orchestration must use the checked scanner entry point")
    }

    fn scan_one_checked(&self, name: &str, seq: &[u8]) -> CResult<Vec<Hsp>> {
        self.calls.lock().unwrap().push((name.into(), seq.into()));
        if self.fail {
            return Err(CoreError::Engine("synthetic checked failure".into()));
        }
        Ok(if self.hit {
            vec![Hsp {
                hit_id: "human_H".into(),
                hit_description: "synthetic full match".into(),
                evalue: 1e-20,
                bitscore: 100.0 + f64::from(seq[0]),
                bias: 0.0,
                query_start: 0,
                query_end: 128,
                hit_start: 0,
                hit_end: 128,
                rf: "x".repeat(128),
                pp: "*".repeat(128),
                order: 0,
            }]
        } else {
            vec![]
        })
    }
}

fn batch(
    scanner: &Scanner,
    rows: &[(String, Vec<u8>)],
    parallel: bool,
    scheme: &str,
    score: f64,
) -> CResult<Vec<SeqResult>> {
    if parallel {
        run_anarci(
            scanner,
            rows,
            scheme,
            &default_allow(),
            false,
            None,
            score,
            false,
            GermlineMethod::Identity,
        )
    } else {
        anarci(
            scanner,
            rows,
            scheme,
            &default_allow(),
            false,
            None,
            score,
            false,
            GermlineMethod::Identity,
        )
    }
}

#[test]
fn exhaustive_byte_alphabet_and_diagnostic_positions() {
    for byte in 0..=u8::MAX {
        for position in 0..3 {
            let mut seq = *b"AAA";
            seq[position] = byte;
            let result = validate_scan_input("synthetic", &seq);
            assert_eq!(
                result.is_ok(),
                byte.is_ascii_alphabetic() || byte == b'*',
                "byte {byte}, position {position}"
            );
            if let Err(error) = result {
                assert!(error
                    .to_string()
                    .contains(&format!("position {}", position + 1)));
            }
        }
    }
    for length in [0, 1, 99_999, 100_000, 100_001] {
        assert_eq!(
            validate_scan_input("", &vec![b'A'; length]).is_ok(),
            length <= 100_000
        );
    }
}

#[test]
fn every_invalid_batch_member_is_rejected_before_any_scan() {
    for parallel in [false, true] {
        for position in 0..4 {
            for (name, seq) in [
                ("bad\0name", b"AAA".to_vec()),
                ("bad-byte", b"AA-AA".to_vec()),
                ("too-long", vec![b'A'; 100_001]),
            ] {
                let scanner = Scanner::default();
                // Duplicate valid rows make validation-before-dedup observable too.
                let mut rows = vec![("valid".into(), b"AAA".to_vec()); 4];
                rows[position] = (name.into(), seq);
                assert!(matches!(
                    batch(&scanner, &rows, parallel, "imgt", 80.0),
                    Err(CoreError::Assertion(_))
                ));
                assert!(scanner.calls.lock().unwrap().is_empty());
            }
        }
        for score in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let scanner = Scanner::default();
            assert!(batch(&scanner, &[], parallel, "imgt", score).is_err());
            assert!(scanner.calls.lock().unwrap().is_empty());
        }
        let scanner = Scanner::default();
        assert!(batch(&scanner, &[], parallel, "unknown", 80.0).is_err());
        assert!(scanner.calls.lock().unwrap().is_empty());
    }
}

#[test]
fn checked_failures_are_distinct_from_successful_no_hit_in_all_entry_points() {
    let rows = vec![("synthetic".into(), b"AAA".to_vec())];
    for parallel in [false, true] {
        for fail in [false, true] {
            let scanner = Scanner {
                fail,
                ..Scanner::default()
            };
            match batch(&scanner, &rows, parallel, "imgt", 80.0) {
                Err(CoreError::Engine(message)) if fail => {
                    assert_eq!(message, "synthetic checked failure")
                }
                Ok(results) if !fail => {
                    assert_eq!(results.len(), 1);
                    assert!(results[0].numbered.is_none());
                    assert!(results[0].details.is_none());
                    assert!(results[0].hit_table.is_empty());
                }
                other => panic!("failure/no-hit conflated: {other:?}"),
            }
            assert_eq!(scanner.calls.lock().unwrap().len(), 1);
            let single = number(
                &scanner,
                b"AAA",
                "imgt",
                &default_allow(),
                None,
                false,
                0,
                80.0,
            );
            if fail {
                assert!(matches!(single, Err(CoreError::Engine(_))));
            } else {
                assert_eq!(single.unwrap(), None);
            }
        }
    }
}

#[test]
fn generated_duplicate_batches_restore_order_names_and_independent_results() {
    for size in [0, 1, 2, 7, 31] {
        let rows: Vec<_> = (0..size)
            .map(|i| (format!("query-{i}"), vec![b"ACD"[(i * 7 + 2) % 3]; 128]))
            .collect();
        let scanner = Scanner {
            hit: true,
            ..Scanner::default()
        };
        let mut got = batch(&scanner, &rows, true, "I", 80.0).unwrap();
        let unique: std::collections::BTreeSet<_> = rows.iter().map(|(_, s)| s).collect();
        assert_eq!(scanner.calls.lock().unwrap().len(), unique.len());
        let serial = batch(
            &Scanner {
                hit: true,
                ..Scanner::default()
            },
            &rows,
            false,
            "imgt",
            80.0,
        )
        .unwrap();
        assert_eq!(got.len(), rows.len());
        for (i, ((name, sequence), result)) in rows.iter().zip(&got).enumerate() {
            assert_eq!(result.numbered, serial[i].numbered);
            let details = &result.details.as_ref().unwrap()[0];
            assert_eq!(&details.query_name, name);
            assert_eq!(details.scheme, "imgt");
            assert_eq!(details.bitscore, 100.0 + f64::from(sequence[0]));
            assert_eq!((details.query_start, details.query_end), (Some(0), 128));
            assert!(details.germlines.is_none());
            assert_eq!(result.hit_table.len(), 1);
            let recovered: Vec<_> = result.numbered.as_ref().unwrap()[0]
                .residues
                .iter()
                .filter(|r| r.aa != b'-')
                .map(|r| r.aa)
                .collect();
            assert_eq!(&recovered, sequence);
        }
        if got.len() > 3 {
            got[0].details.as_mut().unwrap()[0].query_name.clear();
            assert_eq!(got[3].details.as_ref().unwrap()[0].query_name, "query-3");
        }
    }
}

#[test]
fn scheme_aliases_and_invalid_schemes_have_total_resolution() {
    for (short, long) in [
        ("i", "imgt"),
        ("k", "kabat"),
        ("c", "chothia"),
        ("m", "martin"),
        ("a", "aho"),
        ("w", "wolfguy"),
    ] {
        for spelling in [
            short.to_string(),
            short.to_uppercase(),
            long.to_string(),
            long.to_uppercase(),
        ] {
            assert_eq!(resolve_scheme(&spelling).unwrap(), long);
        }
    }
    for scheme in ["", " img t", "imgt\0", "snowman-☃"] {
        assert!(matches!(
            resolve_scheme(scheme),
            Err(CoreError::Assertion(_))
        ));
    }
}

#[test]
fn rescue_rejects_mismatched_or_out_of_bounds_states_without_scanning() {
    let scanner = Scanner::default();
    let detail = anarci_core::DomainDetails {
        id: "human_H".into(),
        description: "synthetic".into(),
        evalue: 1e-20,
        bitscore: 100.0,
        bias: 0.0,
        query_start: Some(0),
        query_end: 1,
        species: "human".into(),
        chain_type: "H".into(),
    };
    let mut states = vec![vec![]];
    assert!(check_for_j(&scanner, "synthetic", b"A", &mut states, &mut []).is_err());
    for si in [1, usize::MAX] {
        let mut states = vec![vec![State {
            id: 104,
            typ: StateType::M,
            si: Some(si),
        }]];
        assert!(check_for_j(
            &scanner,
            "synthetic",
            b"A",
            &mut states,
            &mut [detail.clone()]
        )
        .is_err());
    }
    assert!(scanner.calls.lock().unwrap().is_empty());
}

#[test]
fn malformed_public_alignment_states_return_errors_instead_of_panicking() {
    // These are public State values; catch_unwind makes a panic an explicit contract failure.
    let mut violations = Vec::new();
    for states in [
        vec![State {
            id: 0,
            typ: StateType::M,
            si: Some(0),
        }],
        vec![State {
            id: 129,
            typ: StateType::M,
            si: Some(0),
        }],
        vec![State {
            id: 23,
            typ: StateType::M,
            si: None,
        }],
        vec![State {
            id: 23,
            typ: StateType::M,
            si: Some(usize::MAX),
        }],
        // Individually in-bounds states must not let smoothing synthesize id zero.
        vec![
            State {
                id: 1,
                typ: StateType::D,
                si: None,
            },
            State {
                id: 1,
                typ: StateType::I,
                si: Some(0),
            },
            State {
                id: 1,
                typ: StateType::M,
                si: Some(1),
            },
            State {
                id: 23,
                typ: StateType::M,
                si: Some(2),
            },
        ],
    ] {
        let outcome = std::panic::catch_unwind(|| {
            number_sequence_from_alignment(&states, b"AAA", "imgt", Some("H"))
        });
        match outcome {
            Err(_) => violations.push(format!("panicked: {states:?}")),
            Ok(Ok(_)) => violations.push(format!("accepted: {states:?}")),
            Ok(Err(_)) => {}
        }
    }
    assert!(
        violations.is_empty(),
        "malformed public alignments: {violations:?}"
    );
}
