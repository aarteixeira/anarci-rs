//! Fetch HMMER 3.4 + Easel at build time, build the static libraries (cached in
//! OUT_DIR), and generate FFI bindings with bindgen.
//!
//! Source acquisition (no silent acceptance of unverified bytes):
//!   1. If `$HMMER_TARBALL` is set, use that local `hmmer-3.4.tar.gz` (offline / CI).
//!   2. Otherwise download `HMMER_URL` with `curl`.
//! Either way the tarball's SHA-256 is verified against `HMMER_SHA256` before use.
//!
//! The extracted source + built `.a` files live under `OUT_DIR`, so they are cached
//! across Rust rebuilds and re-created (re-downloaded) after `cargo clean`.

use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::process::Command;

const HMMER_DIR: &str = "hmmer-3.4";
const HMMER_TARBALL_NAME: &str = "hmmer-3.4.tar.gz";
const HMMER_URL: &str = "http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz";
const HMMER_SHA256: &str = "ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3";

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let hmmer = out.join(HMMER_DIR);
    let easel = hmmer.join("easel");
    let libhmmer = hmmer.join("src/libhmmer.a");
    let libeasel = easel.join("libeasel.a");

    if !libhmmer.exists() || !libeasel.exists() {
        if !hmmer.join("configure").exists() {
            obtain_source(&out);
        }
        build_c(&hmmer);
        assert!(libhmmer.exists(), "configure/make did not produce {libhmmer:?}");
        assert!(libeasel.exists(), "configure/make did not produce {libeasel:?}");
    }

    // Link the static archives. Easel must come after hmmer (hmmer depends on
    // easel symbols); the linker resolves left-to-right for static archives.
    println!("cargo:rustc-link-search=native={}", hmmer.join("src").display());
    println!("cargo:rustc-link-search=native={}", easel.display());
    println!("cargo:rustc-link-lib=static=hmmer");
    println!("cargo:rustc-link-lib=static=easel");

    let bindings = bindgen::Builder::default()
        .header(manifest.join("wrapper.h").to_str().unwrap())
        .clang_arg(format!("-I{}", hmmer.join("src").display()))
        .clang_arg(format!("-I{}", easel.display()))
        .allowlist_function("p7_.*")
        .allowlist_type("P7_.*")
        .allowlist_var("p7_.*")
        .allowlist_function("esl_.*")
        .allowlist_type("ESL_.*")
        .allowlist_type("esl_.*")
        .allowlist_var("esl.*")
        .allowlist_var("eslOK")
        .allowlist_var("eslEOF")
        .allowlist_var("eslAMINO")
        .default_enum_style(bindgen::EnumVariation::Consts)
        .layout_tests(false)
        .generate()
        .expect("bindgen failed");

    bindings.write_to_file(out.join("bindings.rs")).expect("write bindings.rs");

    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=HMMER_TARBALL");
    println!("cargo:rerun-if-changed={}", libhmmer.display());
}

/// Obtain + verify the tarball, then extract `hmmer-3.4/` into `out`.
fn obtain_source(out: &Path) {
    let tarball = match std::env::var_os("HMMER_TARBALL") {
        Some(p) => {
            let p = PathBuf::from(p);
            assert!(p.exists(), "HMMER_TARBALL={p:?} does not exist");
            p
        }
        None => {
            let dest = out.join(HMMER_TARBALL_NAME);
            download(HMMER_URL, &dest);
            dest
        }
    };

    let actual = sha256_hex(&tarball);
    assert_eq!(
        actual, HMMER_SHA256,
        "HMMER tarball SHA-256 mismatch for {tarball:?}\n  expected {HMMER_SHA256}\n  got      {actual}\n\
         Refusing to build from unverified source. If you provided $HMMER_TARBALL, ensure it is the \
         official hmmer-3.4.tar.gz."
    );

    run(
        Command::new("tar").arg("xzf").arg(&tarball).arg("-C").arg(out),
        "tar extract",
    );
}

fn download(url: &str, dest: &Path) {
    let status = Command::new("curl")
        .args(["-fsSL", "--retry", "3", "-o"])
        .arg(dest)
        .arg(url)
        .status();
    match status {
        Ok(s) if s.success() => {}
        Ok(s) => panic!(
            "downloading {url} failed (curl exit {s}). For offline builds, set \
             $HMMER_TARBALL to a local hmmer-3.4.tar.gz."
        ),
        Err(e) => panic!(
            "could not run curl to download {url}: {e}. Install curl, or set \
             $HMMER_TARBALL to a local hmmer-3.4.tar.gz."
        ),
    }
}

fn sha256_hex(path: &Path) -> String {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {path:?}: {e}"));
    let digest = Sha256::digest(&bytes);
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

fn build_c(hmmer: &Path) {
    if !hmmer.join("Makefile").exists() {
        run(
            Command::new("./configure")
                .arg("--enable-static")
                .arg("--disable-shared")
                .current_dir(hmmer),
            "configure",
        );
    }
    let jobs = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    run(Command::new("make").arg(format!("-j{jobs}")).current_dir(hmmer), "make");
}

fn run(cmd: &mut Command, what: &str) {
    let status = cmd.status().unwrap_or_else(|e| panic!("failed to spawn {what}: {e}"));
    assert!(status.success(), "{what} failed with status {status}");
}
