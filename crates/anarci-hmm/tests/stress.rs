//! Stability/leak smoke: scan the full fixture set 30x (~30k scans) and confirm
//! RSS stays bounded — guards the per-scan alloc/free discipline.
use anarci_hmm::Engine;
use std::path::PathBuf;
#[test]
#[ignore] // run explicitly: cargo test -p anarci-hmm --release --test stress -- --ignored --nocapture
fn stress_no_leak() {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference_data/dat/HMMs/ALL.hmm");
    let e = Engine::load(&p).unwrap();
    let fx = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/hsps.json.gz");
    let f = std::fs::File::open(&fx).unwrap();
    let v: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(f)).unwrap();
    let seqs: Vec<(String,String)> = v["sequences"].as_array().unwrap().iter()
        .map(|s| (s["id"].as_str().unwrap().to_string(), s["seq"].as_str().unwrap().to_string())).collect();
    let mut total = 0u64;
    for _ in 0..30 {
        for (id, seq) in &seqs { total += e.scan(id, seq.as_bytes()).len() as u64; }
    }
    eprintln!("stress: {} scans, {} total HSPs", seqs.len()*30, total);
    let rss = rss_kb();
    eprintln!("RSS after {} scans: {} MB", seqs.len()*30, rss/1024);
    assert!(rss < 2_000_000, "RSS {} KB looks like a leak", rss);
}
fn rss_kb() -> u64 {
    let out = std::process::Command::new("ps").args(["-o","rss=","-p",&std::process::id().to_string()]).output().unwrap();
    String::from_utf8_lossy(&out.stdout).trim().parse().unwrap_or(0)
}
