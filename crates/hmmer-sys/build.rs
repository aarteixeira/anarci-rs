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
const HMMER_URL: &str = "https://eddylab.org/software/hmmer/hmmer-3.4.tar.gz";
const HMMER_SHA256: &str = "ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3";

fn main() {
    let host = std::env::var("HOST").expect("HOST");
    let target = std::env::var("TARGET").expect("TARGET");
    assert_eq!(
        host, target,
        "HMMER source builds require a native toolchain; cross compilation is not configured"
    );
    assert!(
        target.ends_with("-apple-darwin") || target.ends_with("-unknown-linux-gnu"),
        "Supported build targets are native macOS and glibc Linux"
    );
    assert!(
        matches!(
            std::env::var("CARGO_CFG_TARGET_ARCH").as_deref(),
            Ok("x86_64" | "aarch64")
        ),
        "Supported build architectures are x86_64 and aarch64"
    );
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let hmmer = out.join(HMMER_DIR);
    let easel = hmmer.join("easel");
    let libhmmer = hmmer.join("src/libhmmer.a");
    let libeasel = easel.join("libeasel.a");
    let config = format!(
        "{target}|{:?}|{:?}|{:?}|{:?}|pic-baseline-v1",
        std::env::var_os("CC"),
        std::env::var_os("CFLAGS"),
        std::env::var_os("CPPFLAGS"),
        std::env::var_os("LDFLAGS")
    );
    let stamp = out.join("native-build-config.txt");
    let changed = std::fs::read_to_string(&stamp).ok().as_deref() != Some(config.as_str());
    if changed && hmmer.join("Makefile").exists() {
        run(
            Command::new("make").arg("distclean").current_dir(&hmmer),
            "make distclean",
        );
    }

    if changed || !libhmmer.exists() || !libeasel.exists() {
        if !hmmer.join("configure").exists() {
            obtain_source(&out);
        }
        build_c(&hmmer);
        assert!(
            libhmmer.exists(),
            "configure/make did not produce {libhmmer:?}"
        );
        assert!(
            libeasel.exists(),
            "configure/make did not produce {libeasel:?}"
        );
        std::fs::write(stamp, config).expect("write native build configuration");
    }

    // Link the static archives. Easel must come after hmmer (hmmer depends on
    // easel symbols); the linker resolves left-to-right for static archives.
    println!(
        "cargo:rustc-link-search=native={}",
        hmmer.join("src").display()
    );
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

    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("write bindings.rs");

    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=HMMER_TARBALL");
    println!("cargo:rerun-if-env-changed=CC");
    println!("cargo:rerun-if-env-changed=CFLAGS");
    println!("cargo:rerun-if-env-changed=CPPFLAGS");
    println!("cargo:rerun-if-env-changed=LDFLAGS");
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
        Command::new("tar")
            .arg("xzf")
            .arg(&tarball)
            .arg("-C")
            .arg(out),
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
        let mut cflags = std::env::var("CFLAGS").unwrap_or_else(|_| "-O3".to_owned());
        cflags.push_str(" -fPIC");
        if std::env::var("CARGO_CFG_TARGET_ARCH").as_deref() == Ok("x86_64") {
            cflags.push_str(" -march=x86-64 -mtune=generic");
        }
        run(
            Command::new("./configure")
                .env("CFLAGS", cflags)
                .arg("--enable-static")
                .arg("--disable-shared")
                .current_dir(hmmer),
            "configure",
        );
    }
    let jobs = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    run(
        Command::new("make")
            .arg(format!("-j{jobs}"))
            .current_dir(hmmer),
        "make",
    );
}

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn {what}: {e}"));
    assert!(status.success(), "{what} failed with status {status}");
}
