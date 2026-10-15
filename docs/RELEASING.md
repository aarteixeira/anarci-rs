# Release procedure

Build artifacts locally and in CI. Publication is a separate, explicitly authorized
step. The release workflow uploads build artifacts to its workflow run and has no
package-index credentials or publish job.

## Scope and support

The configured candidate matrix is CPython 3.9–3.13 on glibc Linux x86_64/ARM64 and
macOS x86_64/ARM64. Qualify each wheel with successful CI and a clean installation.
Record the actual results in the release evidence ledger. A local macOS result
qualifies only that interpreter, OS and architecture.

## Data rights and attribution

The official IMGT terms page inspected on 2026-10-15 states: “The IMGT® data are
provided to both public and private users under the CC BY-4.0 license.” It separates
tool use, for which it asks users to contact IMGT, from data/metadata licensing.
Source: https://www.imgt.org/about/termsofuse.php.

The current terms grant a general data/metadata license without a stated date
restriction on the data. CC BY 4.0 permits sharing and adaptation with attribution.
The relevant official paragraphs are retained in `licenses/IMGT-DATA-TERMS.txt`;
the complete license is in `licenses/CC-BY-4.0.txt`. The implementation is this project's Rust code; ANARCI provides compatibility
behavior and the provenance of the reference data. Retain the relevant BSD license
texts, source attribution and adaptation notices from `NOTICE.md`. The bundled
native components are HMMER and Easel.
No individual correspondence or separate distributor permission was obtained.

## Artifact checks

1. Run `python scripts/generate_rust_notices.py --check` to verify the retained
   dependency license texts match the lockfile. Complete safe unit/API checks and the independent implementation review. Use
   the selected synthetic test inputs; unrelated historical fixture sets are outside
   the release test suite.
2. Run `maturin build --release --interpreter python --out <wheel-directory>` and
   `maturin sdist --out <source-directory>` from the intended versioned checkout.
3. Inspect the wheel and sdist members. Require both Python packages, the native
   extension in the wheel, all license texts, and both HMM inputs in the sdist.
   Confirm the sdist contains only the Rust source, bundled inputs, Python package,
   selected Python tests, build metadata and required documentation. Keep workbench,
   private local files and arbitrary historical fixture sets out of both artifacts.
4. Install each wheel in a clean environment outside the checkout. Run
   `python -I -m unittest discover -s <extracted-sdist>/tests/python -v` and
   `anarci-rs --help`. Check `import anarci_rs` and `import anarci`.
5. Extract the sdist into a fresh directory, build a wheel there and repeat the clean
   installation checks. Do not use checkout-relative files for this build.
6. Run `python -m twine check <artifacts>`. Record SHA-256 checksums, versions,
   platform/interpreter tags, checks run, known limits and remote CI run links.
7. Check license/attribution completeness, resolve platform failures and review findings before proposing an
   upload. Obtain user authorization for the exact reviewed artifacts and index.

## Release evidence

Record only completed checks. Qualify remote platforms after their workflow runs
finish successfully. Keep the version, artifact hash and run link with each result.
