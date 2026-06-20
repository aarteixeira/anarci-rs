//! Coordinate/rf/pp sanity check on a single known sequence before the full
//! parity gate. Ground truth taken from tests/fixtures/hsps.json.gz (id
//! "p12e8_12E8_H", profile mouse_H, domain 0).

use anarci_hmm::Engine;
use std::path::PathBuf;

fn hmm_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../reference_data/dat/HMMs/ALL.hmm")
}

#[test]
fn smoke_p12e8_mouse_h() {
    let seq = b"EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD";

    let engine = Engine::load(&hmm_path()).expect("load ALL.hmm");
    assert_eq!(engine.n_models(), 29, "expected 29 profiles");

    let hsps = engine.scan("p12e8_12E8_H", seq);
    assert!(!hsps.is_empty(), "no HSPs produced");

    // Find the primary mouse_H domain (query_start == 0).
    let h = hsps
        .iter()
        .find(|h| h.hit_id == "mouse_H" && h.query_start == 0)
        .expect("no mouse_H domain at query_start 0");

    eprintln!(
        "mouse_H: bit={} eval={:.3e} bias={} qs={} qe={} hs={} he={}",
        h.bitscore, h.evalue, h.bias, h.query_start, h.query_end, h.hit_start, h.hit_end
    );
    eprintln!("rf: {}", h.rf);
    eprintln!("pp: {}", h.pp);

    assert!((h.bitscore - 173.4).abs() < 0.05, "bitscore {}", h.bitscore);
    assert!(
        (h.evalue - 3.3e-54).abs() < 3.3e-54 * 0.05,
        "evalue {:.3e}",
        h.evalue
    );
    assert!((h.bias - 0.3).abs() < 0.05, "bias {}", h.bias);
    assert_eq!(h.query_start, 0);
    assert_eq!(h.query_end, 120);
    assert_eq!(h.hit_start, 0);
    assert_eq!(h.hit_end, 128);
    assert_eq!(h.rf.len(), 128);
    assert_eq!(
        h.rf,
        "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    );
    assert_eq!(
        h.pp,
        "79*******.********************....9************************..***********.*********************************999999***************7"
    );
}
