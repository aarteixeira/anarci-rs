# Release evidence

Version 0.1.0 is a local release candidate. The configured remote wheel matrix
has not run. No package-index upload has occurred.

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
- Legacy argument order, input containers, output generation, species-filter
  edge cases and deprecated helper insertion rules.
- Missing source-build resources, license materials and installed CLI delivery.

The final independent source review found no remaining actionable issue within
its reviewed scope. That result is separate from artifact and platform tests.

## Completed native checks

| Check | Result |
|---|---|
| HMM engine robustness unit tests | 8 passed |
| Core native input contracts | 4 passed |
| Existing core source unit tests | 19 passed |
| Shared CLI library tests | 23 passed |
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

Local macOS ARM64 wheels passed the 40-test installed-package suite on CPython
3.9, 3.10, 3.11, 3.12 and 3.13. The default builds skip the optional saved-text
parser test and three model tests; those features have separate passing checks.
Final clean wheel/source installation records and artifact hashes are delivered
beside the candidate files. Linux x86_64/aarch64 and macOS Intel require successful
remote CI before they qualify as tested release targets. The workflow includes
no publication job.

Local workflow checks include actionlint 1.7.12 and a fresh Cargo-cache simulation.
The latter reproduced a missing-platform-dependency failure in the offline notices
check; a full locked fetch before that check resolved it. The release validator
uses symlink environments on POSIX so standalone Python launchers retain their
base runtime's shared-library lookup.
