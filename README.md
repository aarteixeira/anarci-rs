# anarci-rs

Antibody and T-cell receptor variable-domain numbering with an in-process Rust
backend and Python interfaces. HMMER and the profile databases are bundled. Batch
scans run in Rust, release Python's global interpreter lock and deduplicate
identical sequences. Set `ncpu` to select multiple workers.

The distribution provides two interfaces:

- `anarci_rs` keeps the fast pan engine and the extended configuration options.
- `anarci` supplies compatibility entry points with reference ANARCI defaults.

## Install

Install a wheel built for your Python version and platform:

```bash
python -m pip install path/to/anarci_rs-0.1.0-<python>-<abi>-<platform>.whl
```

The distribution name is `anarci-rs`. Once a release is published on PyPI, install
it with `python -m pip install anarci-rs`. Supported build targets and source-build
requirements are in [BUILDING.md](docs/BUILDING.md). A wheel includes the native
engine; it does not need an external `hmmscan` executable.

Use a separate environment from the original ANARCI distribution. Both packages
provide the `anarci` import path, so installing them together can overwrite files.

## Number sequences from Python

```python
import anarci_rs

# records is a list of (identifier, amino_acid_sequence) pairs.
numbered, details, hits = anarci_rs.anarci(records, database="pan")

# File input, parallel batch, V/J calls and region coverage.
sequences, numbered, details, hits = anarci_rs.run_anarci(
    "antibodies.fasta",
    ncpu=4,
    assign_germline=True,
    germline_method="evalue",
    annotate_regions=True,
)

# Single sequence; lower these gates explicitly for partial domains.
positions, chain_class = anarci_rs.number(
    sequence, min_length=40, bit_score_threshold=65,
)
```

`anarci` returns three lists in input order. A numbered sequence contains a list of
domains, each `(positions, start, end)`. Each position is
`((position_number, insertion_code), amino_acid)`. Sequences with no recognized
domain have `None` entries. `run_anarci` also returns the input records.
`number` returns the first domain's positions and chain class, or `(False, False)`.
Its light-chain class is `L` for both kappa and lambda.

### Engines and configuration

| Option | Values and behavior |
|---|---|
| `database` | `"pan"` uses seven pan-species profiles; `"ALL"` uses the reference 29-profile set. Custom HMM files use the native loader. |
| `scheme` | `imgt`, `chothia`, `kabat`, `martin`, `aho`, `wolfguy`, plus ANARCI's short aliases. Some schemes apply only to antibody chains. |
| `allow` | A subset of `H`, `K`, `L`, `A`, `B`, `G`, `D`. |
| `allowed_species` | Species filter. Extended batch calls default to all supported species; the compatibility interface defaults to human/mouse. |
| `assign_germline` | Add V/J gene assignments to domain details. |
| `germline_method` | `identity` or `evalue`. The pan default is `evalue`; the exact default is `identity`. |
| `bit_score_threshold` | HMM domain acceptance cutoff; default 80. Profile sets have different score distributions. |
| `ncpu` | Number of Rust workers for batch calls. See the API reference for accepted edge values. |
| `annotate_regions` | Add IMGT region completeness and covered IMGT positions to domain details. |
| `min_length` | Single-sequence `number` length gate; default 70. |

Pan and exact profiles can place gaps differently and report different scores,
species and germlines. Choose `database="ALL"` and matching options when comparing
with the pinned reference ANARCI. The pan engine remains available explicitly even
when using the compatibility interface.

## Use existing ANARCI code

```python
from anarci import anarci, run_anarci, number

numbered, details, hits = anarci(records, assign_germline=True)
sequences, numbered, details, hits = run_anarci(records, 4)
positions, chain_class = number(sequence, scheme="imgt")
```

The compatibility functions use `database="ALL"` and human/mouse defaults. Their
computation still runs in Rust. The supported helper imports, output writers,
input behavior and intentional differences are listed in
[COMPATIBILITY.md](docs/COMPATIBILITY.md).

## Additional Python tools

The Python package also exposes native sequence-file reading, six-frame DNA
translation and annotation, PDB/mmCIF chain extraction, structure pairing and
output formatting. These call the same Rust library as the command line.
See [API.md](docs/API.md) for signatures, return fields and configuration.

Structural pairing and contact reports are geometry-based inferences. Contact
counts describe the supplied coordinates; they do not establish an epitope or
binding specificity independently of experimental evidence.

`anarci_rs.experimental` exposes the learned-engine label codec, constrained
decoder, augmentation and dataset utilities. The optional Candle model/training
code remains available in the Rust crate. The package supplies no validated
trained numbering checkpoint. See [EXPERIMENTAL.md](docs/EXPERIMENTAL.md).

## Command line

```bash
anarci-rs antibodies.fasta --engine pan --format tsv
anarci-rs reads.fastq --input-type dna --report-frames --format json
anarci-rs structure.pdb --pair --format json --output pairs.json
anarci-rs antibodies.fasta --engine all --scheme kabat --assign-germline
```

Run `anarci-rs --help` for the full option list. The Rust workspace also retains
the `anarci` binary (`cargo run --release -p anarci-cli -- ...`).

## Correctness and performance

Release checks compare the compatibility interface with a pinned reference
ANARCI installation and exercise the added APIs, failure paths and installed
artifacts. The release suite uses a small, identified reference set and synthetic
cases. Passing that suite establishes the tested contracts; it does not prove
identical behavior for every possible sequence or custom HMM.

Run the Python tests against an installed build:

```bash
python -m unittest discover -s tests/python -v
python scripts/validate_release.py --help
python scripts/benchmark_release.py --help
```

Performance depends on the profile set, batch size, worker count and duplicate
rate. Benchmark exact and pan modes separately and compare equal work. Current
release measurements and validation limits are in
[RELEASE_EVIDENCE.md](docs/RELEASE_EVIDENCE.md); historical measurements in
`BENCHMARKS.md` are tied to their original code and datasets.

## Development

```bash
maturin develop --release
maturin build --release
maturin sdist
```

See [BUILDING.md](docs/BUILDING.md) and [RELEASING.md](docs/RELEASING.md) for the
native toolchain, wheel targets, clean source builds and release checks.

| Path | Purpose |
|---|---|
| `crates/anarci-core` | Schemes, state vectors, orchestration, germlines and regions |
| `crates/anarci-hmm` | In-process HMMER engine |
| `crates/hmmer-sys` | Verified HMMER/Easel source build and native bindings |
| `crates/anarci-py` | Python native extension |
| `crates/anarci-cli` | Shared annotation/structure library and Rust command |
| `crates/anarci-learned` | Experimental codec, decoder, data tools and optional model |
| `python/` | Extended Python API and ANARCI compatibility package |

## Licenses and attribution

Project code is BSD-3-Clause. Bundled components and reference data have their
own notices and attribution requirements; see [NOTICE.md](NOTICE.md) and the
packaged license files. Review those terms before redistributing the artifacts.
