//! Build the vendored HMMER 3.4 + Easel static libraries (cached) and generate
//! FFI bindings with bindgen.
//!
//! The C build (`./configure && make`) is the slow part, so we cache it: if
//! `src/libhmmer.a` and `easel/libeasel.a` already exist in the vendor tree we
//! skip straight to linking + bindgen. `cargo clean` does not touch the vendor
//! tree, so the C build survives across Rust rebuilds.

use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let hmmer = manifest.join("vendor/hmmer-3.4");
    let easel = hmmer.join("easel");

    let libhmmer = hmmer.join("src/libhmmer.a");
    let libeasel = easel.join("libeasel.a");

    if !libhmmer.exists() || !libeasel.exists() {
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

    // Generate bindings. Headers and the generated config headers both live in
    // the vendor tree; expose them on the include path.
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
        // Keep enums as plain integer constants for ergonomic FFI.
        .default_enum_style(bindgen::EnumVariation::Consts)
        .layout_tests(false)
        .generate()
        .expect("bindgen failed");

    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("write bindings.rs");

    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed=build.rs");
    // Trigger a rebuild if the C archives are removed.
    println!("cargo:rerun-if-changed={}", libhmmer.display());
}

fn build_c(hmmer: &Path) {
    // configure (autodetects SIMD: NEON on arm64, SSE on x86_64).
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
    run(
        Command::new("make").arg(format!("-j{jobs}")).current_dir(hmmer),
        "make",
    );
}

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn {what}: {e}"));
    assert!(status.success(), "{what} failed with status {status}");
}
