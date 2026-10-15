# PyPI release plan

## Outcome and scope

Release `anarci-rs` as an installable Python distribution. Keep the in-process Rust
HMMER engine, parallel batches, deduplication, exact and pan profiles, both germline
methods, partial-domain controls and region annotation. Expose the existing sequence
I/O, six-frame translation, structural parsing, pairing and contact reports through
Python as well as the CLI. Preserve the learned-engine research code and document its
experimental status; do not imply that a validated trained inference model exists.

No PyPI upload, push, deployment, paid service, or external data transmission is part
of implementation. Existing unrelated local edits must remain intact. Validation uses
synthetic inputs and a permitted antibody-only reference set; excluded biological
material and the unrelated workbench must not enter artifacts or tests.

## API contracts

1. `anarci_rs` is the extended API. Its existing pan defaults and extension arguments
   remain available. Native work runs in Rust with the GIL released.
2. `anarci` is the compatibility package. Its three high-level functions use the
   reference ANARCI signatures, defaults, return shapes, I/O and exception behavior.
   Existing callers can keep `import anarci`, including common helper/submodule
   imports. Both distributions must not be installed together; test isolation.
3. Inventory every existing callable/configuration. Record its Python entry point,
   tests, stability and any unsupported legacy behavior in a coverage matrix.
4. Reuse core and CLI library functions. Do not wrap subprocesses for computation or
   reimplement native algorithms in Python. Validate Python inputs before FFI calls.
5. Treat empty input, invalid arguments, no-hit results, multiple/partial domains,
   all numbering schemes and all chain types as explicit contracts.

## Work packages

1. **Plan review.** Independently challenge this plan and the API inventory. Resolve
   each actionable issue and record decisions before implementation.
2. **Python compatibility.** Build a mixed Maturin package with a private native
   extension. Supply compatibility wrappers, helper exports, FASTA/gzip input,
   ANARCI text/CSV output, exact argument order and reference error semantics.
   Preserve all extension options under `anarci_rs`.
3. **Native feature access.** Extract reusable CLI library entry points for sequence
   reading, translation, protein/DNA annotation, PDB/mmCIF parsing, structural
   pairing and serialized output. Make the CLI and Python use the same functions.
   Correct silent errors, engine/method routing and scheme-dependent CDR mapping.
4. **Native robustness.** Eliminate shared temporary-profile races, validate state
   vectors and numeric settings, and exercise concurrency and worker counts.
5. **Packaging.** Include embedded HMM inputs and all required license texts in
   source archives and wheels. Exclude private/unrelated material. Set a tested
   Python support range and truthful metadata. Validate independent source builds.
6. **Release automation.** Add clean-install and wheel CI for the declared Python/
   operating-system matrix. Start with Linux and macOS; Windows support requires
   successful native build evidence before being declared. Use pinned/checksummed
   upstream sources and secure artifact upload configuration with no automatic
   publication from ordinary commits.
7. **Documentation.** Provide quick starts, migration/compatibility notes, complete
   feature/configuration reference, CLI use, output schemas, performance method,
   installation/build instructions and a release checklist. Separate compatibility
   defaults from extended defaults, measured evidence from historical claims, and
   validated behavior from experimental research tools.
8. **Independent implementation review.** Challenge actual code and artifacts with
   concrete failure cases. Add meaningful regressions for confirmed defects. Repeat
   targeted review and verification until no unresolved release-blocking issue remains.

## Acceptance gates

- The native profile cache must retain its engine's buffers, key by a non-reused
  identity and replace clones on engine switches. Test alternating pan/exact engines,
  equal-profile-count engine reloads and engine destruction on a persistent worker.
- Protect HMMER's global initialization with `Once`; exercise concurrent first use.
- Every route must reject NUL, invalid alphabet, nonfinite thresholds and oversized
  inputs before C calls. Check native return codes and distinguish failures from
  valid no-hit results. Subprocess regressions must establish no interpreter abort.
- Derive HMM E-value scaling from the loaded database's model count; preserve exact
  29-profile results and test the seven-profile pan search-space convention.
- Build/install a fresh wheel and source archive without relying on the checkout.
- Check archive membership, metadata, notices, absence of unrelated data and runtime
  independence from an external HMMER installation.
- Compare the compatibility API against the pinned installed ANARCI reference with
  exact structures and explicit numeric tolerances where native floating-point
  representation requires them. Include signatures, errors, stdout/files, gzip,
  species/chain filters and helper imports; do not call numbering parity full API parity.
- Cover native features with synthetic tests, permitted reference proteins and
  synthetic structure geometry. Exercise all supported schemes, engines, germline
  methods, input formats and configuration branches. Test failing paths as well.
- Run relevant existing Rust tests using safe test selections. Do not load excluded
  fixture sets. Identify any validation gap this restriction leaves.
- Benchmark exact and pan Rust backends against reference with equal work/settings;
  record interpreter, OS, CPU, worker count, dataset identity, warmup and repetitions.
  Check the Python layer does not introduce a material regression in native batches.
- Run formatting/static checks and documentation examples, then install artifacts
  in a clean environment and rerun the smoke suite.
- Test the installed `anarci-rs` command: help, stdin, output file, pairing output,
  invalid flags and nonzero failure status. Test explicit input-type precedence.
- Declare a concrete initial build matrix: CPython 3.9 through 3.13, Linux x86_64
  and aarch64 (manylinux_2_28), and macOS arm64 / x86_64 (macOS 11+). Other platforms
  remain source-build candidates until verified. Inspect wheel linked libraries and
  avoid build-host CPU tuning. Distinguish local evidence from pending remote jobs.
- Pin the reference package version and source-file hashes in the validation record.
  Enumerate helper imports; test worker edge cases, FASTQ plus/quality checks and
  explicit unsupported `hmmerpath` use. Custom HMM files use native loading.
- Record independent review findings, resolutions, checks and platform evidence.
  A CI workflow definition is not evidence that remote CI passed.
- Record bundled-component licenses and data attribution terms. Do not declare
  release readiness while a required legal or platform gate remains unverified.

## Completion rule

Deliver local release artifacts, tests, documentation and an evidence ledger. Report
readiness only for demonstrated platforms and API contracts. List external gates
plainly. Publishing requires the user's authorization after artifact review.

## Plan red-team record

Review 1 found nine concrete gaps: engine-cache ownership, global initialization,
batch input/FFI validation, swallowed failures, CLI packaging, platform-specific
evidence, redistribution notices/terms, compatibility boundaries and pan E-value
scaling. The acceptance gates above now cover each issue. Implementation ownership
is split across native robustness, shared CLI features, Python compatibility,
experimental tools, packaging and independent validation. Review 2 checks these
revisions and implementation-specific risks. Open platform gates prevent an
unqualified claim of release readiness even if local software checks pass.

Review 2 approved the revised plan for implementation. Input validation must preserve
the documented ambiguity alphabet, including translated stop masking to `X`; the
legacy `number()` validator retains its own strict contract. No plan-level issue
remains open. Implementation and artifact evidence remain required.

The current official IMGT terms license data under CC BY 4.0, with no restriction
to a release date in that grant. The packaging work records the terms, attribution,
source provenance and full license separately from the project's BSD license.
This is our Rust implementation, with compatibility imports and bundled components;
the original ANARCI distribution is not the release artifact.
