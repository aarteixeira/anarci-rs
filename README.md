# anarci-rs

A Rust reimplementation of [ANARCI](https://github.com/oxpig/ANARCI) (antibody & TCR
variable-domain numbering) with a Python wrapper. **Drop-in replacement**: change

```python
import anarci
```
to
```python
import anarci_rs as anarci
```
and the function names, signatures, and return values are identical.

All internals are in Rust. The batch path runs fully in Rust with rayon parallelism and
no Python in the hot loop. The HMM search runs in-process via a native HMMER 3.4 binding —
no `hmmscan` subprocess, no temp files, no text parsing.

## Two engines

| Engine | Select | Speed | Numbering |
|---|---|---|---|
| **pan** (default) | `database="pan"` | **~4× faster** (≈5.6× vs reference ANARCI; more on multicore + dedup) | 99.2% identical to ANARCI; ties on conserved IMGT anchors; *more robust* |
| **exact** | `database="ALL"` | baseline | **byte-for-byte identical** to stock ANARCI |

The **pan** engine uses one pan-species HMM per chain type (7 profiles vs ANARCI's 29);
species and germline genes are assigned by sequence identity. It reproduces ANARCI's IMGT
numbering on 99.2% of sequences and is actually *more robust* — it can't pick a wrong-species
profile (a real ANARCI failure mode: 7-species ANARCI mis-numbers ~5% by letting e.g.
`rhesus_K` win on a humanized kappa). The remaining ~0.8% are IMGT-legal gap-placement ties.
Use `database="ALL"` whenever you need output byte-identical to stock ANARCI.

## Why it's faster

Stock ANARCI spends ~98% of its time in the `hmmscan` subprocess and Biopython's text
parsing of its output (profiled: 67% subprocess, 31% parse, <2% in the numbering itself).
anarci-rs eliminates both by running HMMER in-process and doing the parse, numbering,
germline assignment, and batching in Rust. The **pan** engine adds a further ~4× by scanning
7 profiles instead of 29 (the HMM search is the real bottleneck), and the batch path
deduplicates identical sequences (a free win on repetitive data).

## Correctness

The **exact** engine (`database="ALL"`) is validated **byte-for-byte** against reference
ANARCI (conda `anarci 2024.05.21`, HMMER 3.4) on 996 sequences (1013 domains):

| Layer | Gate | Result |
|---|---|---|
| Numbering (imgt, chothia, kabat, martin, aho, wolfguy) | identical numbering from identical state vectors | 1013/1013 domains |
| Germline assignment | identical v/j gene + identity | all domains |
| HMMER-output → state vectors (`parse_hmmer_query`) | identical state vectors from identical HSPs | 996/996 |
| End-to-end `anarci()` (imgt) | identical numbered + details + germlines + hit_tables | 996/996 |
| Native HMM engine | identical state vectors vs stock hmmscan | 996/996 |

The **pan** engine (default) is validated by `scripts/validate_pan.py`: numbering identical
to reference ANARCI run against the same pan HMMs (996/996), 99.20% identical to stock ANARCI,
and tied on conserved-anchor correctness (99.21%).

The pinned reference HMM database and germlines are checked into `reference_data/`
(`ALL.hmm`, `FEW.hmm`, `germlines.py`) so results are reproducible.

**No silent failures or fallbacks.** Every error path raises explicitly; ANARCI's own
errors (e.g. `AssertionError` when numbering a TCR with an antibody-only scheme) are
reproduced exactly rather than worked around.

## Public API

Identical to ANARCI:

```python
import anarci_rs as anarci

# single sequence
numbering, chain_type = anarci.number("EVQLQ...SS", scheme="imgt")

# one or many; numbered / alignment_details / hit_tables out
numbered, details, hits = anarci.anarci([("id1", "EVQ..."), ("id2", "DIV...")],
                                        scheme="imgt", assign_germline=True)

# batch with native Rust parallelism (ncpu controls rayon threads)
seqs, numbered, details, hits = anarci.run_anarci("input.fasta", scheme="chothia", ncpu=8)

# byte-for-byte identical to stock ANARCI (slower exact engine)
numbering, chain_type = anarci.number("EVQLQ...SS", scheme="imgt", database="ALL")
```

All functions take `database="pan"` (default, fast) or `database="ALL"` (exact ANARCI).
Schemes: `imgt`, `chothia`, `kabat`, `martin`, `aho`, `wolfguy`. Chains: heavy, kappa,
lambda, TCR α/β/γ/δ. Species: human, mouse, rat, rabbit, rhesus, pig, alpaca, cow.

## Build

```bash
# Rust core tests (numbering, germline, parse, end-to-end)
cargo test -p anarci-core

# Build the Python extension (bundles HMMER + data; zero system deps)
maturin develop --release        # into the active venv
# or
maturin build --release          # produce a wheel
```

## Layout

```
crates/anarci-core   pure Rust: numbering schemes, germline, parse, orchestration
crates/hmmer-sys     FFI to HMMER 3.4 + Easel (fetched at build time, static-linked)
crates/anarci-hmm    safe in-process HMM scan engine -> HSPs
crates/anarci-py     PyO3 module `anarci_rs` (the drop-in API)
reference_data/      pinned ALL.hmm + germlines (the canonical reference data)
examples/            example FASTAs (from ANARCI, BSD-3) for the scripts/benchmark
tests/fixtures/      golden fixtures captured from reference ANARCI
scripts/             fixture generators + the accuracy/speed benchmark
```

### Build-time HMMER fetch
`hmmer-sys/build.rs` downloads `hmmer-3.4.tar.gz` (SHA-256-pinned, verified) and
compiles it once into `OUT_DIR`. The first build (and any build after `cargo clean`)
needs network + a C toolchain (`curl`, `tar`, `make`, a C compiler). For offline/CI
builds, point `HMMER_TARBALL` at a local copy of the official tarball:
`HMMER_TARBALL=/path/to/hmmer-3.4.tar.gz maturin build --release`.

## Benchmark

```bash
python scripts/benchmark.py --accuracy   # byte-for-byte vs reference ANARCI, all schemes
python scripts/benchmark.py --speed      # throughput vs reference, ncpu sweep
```

See `BENCHMARKS.md` for results. `PLAN.md` documents the full implementation plan.

## License

BSD-3-Clause (matching ANARCI, HMMER, and Easel, which are bundled).
