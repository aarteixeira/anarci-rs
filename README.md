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

## Why

Stock ANARCI spends ~98% of its time in the `hmmscan` subprocess and Biopython's text
parsing of its output (profiled: 67% subprocess, 31% parse, <2% in the numbering itself).
anarci-rs eliminates both by running HMMER in-process and doing the parse, numbering,
germline assignment, and batching in Rust.

## Correctness

anarci-rs is validated **byte-for-byte** against reference ANARCI (conda `anarci 2024.05.21`,
HMMER 3.4) on 996 sequences (1013 domains):

| Layer | Gate | Result |
|---|---|---|
| Numbering (imgt, chothia, kabat, martin, aho, wolfguy) | identical numbering from identical state vectors | 1013/1013 domains |
| Germline assignment | identical v/j gene + identity | all domains |
| HMMER-output → state vectors (`parse_hmmer_query`) | identical state vectors from identical HSPs | 996/996 |
| End-to-end `anarci()` (imgt) | identical numbered + details + germlines + hit_tables | 996/996 |
| Native HMM engine | identical state vectors vs stock hmmscan | see `crates/anarci-hmm/tests` |

The pinned reference HMM database and germlines are checked into `reference_data/`
(`ALL.hmm`, `germlines.py`) so both implementations use identical data.

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
```

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
crates/hmmer-sys     FFI to vendored HMMER 3.4 + Easel (static-linked)
crates/anarci-hmm    safe in-process HMM scan engine -> HSPs
crates/anarci-py     PyO3 module `anarci_rs` (the drop-in API)
reference_data/      pinned ALL.hmm + germlines (the canonical reference data)
tests/fixtures/      golden fixtures captured from reference ANARCI
scripts/             fixture generators + the accuracy/speed benchmark
```

## Benchmark

```bash
python scripts/benchmark.py --accuracy   # byte-for-byte vs reference ANARCI, all schemes
python scripts/benchmark.py --speed      # throughput vs reference, ncpu sweep
```

See `BENCHMARKS.md` for results. `PLAN.md` documents the full implementation plan.

## License

BSD-3-Clause (matching ANARCI, HMMER, and Easel, which are bundled).
