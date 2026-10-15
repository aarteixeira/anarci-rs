# Release evidence

Version 0.1.0 is a release candidate. This file records completed source review
and qualification checks; the delivered manifest binds the final release commit
and artifact hashes to its CI run. No package-index upload has occurred.

## Implemented surface

The release preserves the Rust HMMER backend, pan/exact engines, native batches,
germline methods and partial-domain controls. It adds the ANARCI compatibility
namespace, native Python access to the shared CLI tools, and optional research
model bindings. See `FEATURE_COVERAGE.md` for the callable map and
`COMPATIBILITY.md` for explicit legacy differences.

## Independent review

The plan passed two review rounds. Implementation review found and resolved:

- Cached HMM profile ownership and engine-switch identity errors.
- Concurrent global initialization and unsafe shared temporary model files.
- Invalid batch inputs and swallowed native/CLI failures.
- Incorrect pan E-value search-space scaling.
- Custom HMM name, length and RF-annotation validation gaps.
- Explicit mmCIF/stdin routing, pairing output paths and scheme-independent CDR
  coordinate masks.
- TSV quoting and the meaning of region completeness flags.
- Legacy argument order, input containers, wildcard export names, output
  generation, species-filter edge cases and deprecated helper insertion rules.
- Missing source-build resources, license materials and installed CLI delivery.
- macOS compiler/libclang header alignment and Python 3.9-compatible source
  extraction, with tests for traversal, links, permissions and existing files.

The final independent source review found no remaining actionable issue within
its reviewed scope. That result is separate from artifact and platform tests.

## Completed native checks

| Check | Result |
|---|---|
| HMM engine robustness unit tests | 8 passed |
| Core native input contracts | 4 passed |
| Existing core source unit tests | 19 passed |
| Shared CLI library tests | 23 passed |
| Portable source extraction tests | 3 passed |
| Optional CPU model build | Passed |
| Experimental codec/data tests in optional build | 6 passed |
| Optional model shape/validation/training tests | 3 passed |

These checks use synthetic inputs and embedded immune-model resources. Historical
bulk scientific fixture suites were outside this release validation scope.

## Reference and integration scope

The reference record identifies the exact installed source and database hashes.
Tests cover all six numbering schemes, heavy/light reference outputs, seven
model-derived synthetic chain classes, both engines and germline methods,
multidomain/partial cases, helper outputs, malformed input, concurrent first use,
sequence/structure formats and installed command behavior. The saved-text parser
is checked with its optional Biopython dependency in a separate environment.

The public reference inputs are an identified mouse anti-lysozyme antibody pair;
the other chain-class inputs are deterministic model-derived synthetic consensuses.
These inputs establish tested software contracts rather than exhaustive biological
accuracy. HMM scores use explicitly recorded numeric tolerances.

## Performance

The completed matched-work benchmark and paired wrapper timing are described in
`PERFORMANCE.md`. Full durations and settings are delivered with the artifacts.

## Platform qualification

[CI run 38799456047](https://github.com/aarteixeira/anarci-rs/actions/runs/38799456047),
attempt 2, completed successfully on 2026-10-15.

| Platform | Architecture | CPython | Wheel baseline |
|---|---|---|---|
| Linux | x86_64 | 3.9 through 3.13 | manylinux_2_28 |
| Linux | aarch64 | 3.9 through 3.13 | manylinux_2_28 |
| macOS | x86_64 | 3.9 through 3.13 | macOS 11 |
| macOS | ARM64 | 3.9 through 3.13 | macOS 11 |

Each of the 20 wheels passed metadata checks, a clean installation, the 40-test
Python suite and the installed command check. The default builds skip the optional
saved-text parser test and three model tests; those features have separate passing
checks. The independent Linux source job compiled the extracted archive, ran the
selected Rust and packaging tests, and verified a clean installation. The optional
CPU-model job passed its feature build, notices check and model tests.

One Intel Python 3.10 runner stalled in repository checkout. Its cancelled attempt
was rerun on a fresh runner; the other 21 successful checks were retained. The
successful rerun reached the build and tests without source changes.

Local macOS ARM64 qualification also covers CPython 3.9 through 3.13. The actual
CI-built ARM64 Python 3.10 wheel was installed locally for the additional
performance measurements in `PERFORMANCE.md`.

Final artifact names, hashes, tags, GitHub artifact digests and job records are
delivered in the release manifest. The workflow contains no publication job.

Workflow checks also include actionlint 1.7.12 and a fresh Cargo-cache simulation.
A full locked fetch supplies the cross-platform metadata needed for the offline
notices check. Portable extraction validates member names/types before writing
regular files, and preserves existing destinations. The local release validator
uses POSIX symlink environments so standalone Python launchers retain their base
runtime's shared-library lookup.
