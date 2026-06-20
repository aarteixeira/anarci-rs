# Benchmarks

Reproduce with `python scripts/benchmark.py --accuracy` and `--speed`, run in an env
with both reference `anarci` (conda 2024.05.21) and `anarci_rs` installed.

Machine: Apple Silicon (arm64), 12 cores. Reference: ANARCI conda 2024.05.21 + HMMER 3.4.

## Accuracy — byte-for-byte vs reference ANARCI

Full `antibody_sequences.fasta` (**997 sequences**) × **6 schemes**, `assign_germline=True`,
all 7 species. Compared end-to-end: `numbered` (every `((pos,ins),aa)`), `alignment_details`
(species, chain_type, query_start/end, bitscore, germline v/j genes), via `run_anarci`.

| Scheme | Mismatches |
|---|---|
| imgt | 0 |
| chothia | 0 |
| kabat | 0 |
| martin | 0 |
| aho | 0 |
| wolfguy | 0 |
| **TOTAL** | **0 — PASS** |

Backed by Rust-side gates against fixtures captured from reference ANARCI (996 seqs /
1013 domains): numbering 1013/1013, germline assignment, parse→state-vectors 996/996,
native-engine state-vector parity 996/996, end-to-end `anarci()` 996/996.

## Speed — throughput vs reference ANARCI

`run_anarci`, scheme=imgt, `assign_germline=False`. seq/s (higher is better).

| N | ncpu | reference ANARCI | anarci-rs | speedup |
|---|---|---|---|---|
| 1,000 | 1 | 42.9 | 62.4 | 1.5× |
| 1,000 | 12 | 277.2 | 414.6 | 1.5× |
| 10,000 | 1 | 42.3 | 62.5 | 1.5× |
| 10,000 | 12 | 253.3 | **451.5** | **1.8×** |

### Honest reading of these numbers

The win is real but bounded, and it's worth being precise about why.

- **What anarci-rs removes:** the `hmmscan` subprocess spawn, temp-file I/O, and Biopython's
  text parsing of hmmscan output — and it replaces Python `multiprocessing` (pickling + IPC)
  with rayon (GIL released, shared memory).
- **What neither can remove:** the HMM search itself. Both implementations run the *same*
  HMMER 3.4 pipeline (MSV → Viterbi → Forward → domain definition) against **29 profiles per
  sequence**. That DP compute dominates: single-thread, anarci-rs reaches ~63 seq/s, which is
  essentially the raw HMMER compute ceiling (stock ANARCI's profiled "67% subprocess" was
  mostly `select.poll` *waiting on hmmscan to compute*, not pure overhead).
- So single-thread the gain is ~1.5× (overhead removed). The bigger structural win is
  **multicore scaling**: rayon scales near-linearly (62→451 seq/s, 7.2×) where ANARCI's
  multiprocessing scales sub-linearly (43→253, ~6×), widening the gap to **1.8×** at 12 cores.

We verified the engine is at the compute floor: per-thread reuse of the HMMER pipeline yields
0% (and is actually unsafe — `pli->Z` accumulates across reused scans in SCAN mode, corrupting
E-values), and the per-scan `bg`/`sq` allocations are negligible against the 29-profile DP.

**Going faster would require changing the algorithm** — e.g. an MSV/SSV pre-filter to skip
clearly-non-matching profiles, or scanning fewer profiles per query. Those change which hits
are reported and would break byte-for-byte parity, so they are out of scope for a drop-in
replacement (a future opt-in `--fast` mode could trade exactness for speed).

### The other wins (not in the table)
- **Zero runtime dependencies**: HMMER is statically linked and `ALL.hmm` is embedded in the
  extension. No `hmmscan` on PATH, no Biopython, no temp files. `pip install` and `import`.
- **No subprocess fragility**: no `multiprocessing` spawn/`__main__` re-import footguns, no
  temp-file cleanup, no parsing of a text format that can drift between HMMER versions.
