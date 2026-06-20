# Benchmarks

Reproduce with `python scripts/benchmark.py --accuracy` and `--speed`, run in an env
with both reference `anarci` (conda 2024.05.21) and `anarci_rs` installed.

Machine: Apple Silicon (arm64), 12 cores. Reference: ANARCI conda 2024.05.21 + HMMER 3.4.

## Engines

anarci-rs ships two engines. **pan** (default, 7 pan-species HMMs) is ~4× faster than
**exact** (29 species×chain HMMs, byte-identical to ANARCI), at equal accuracy.

| Engine | `database=` | seq/s (1 core) | numbering |
|---|---|---|---|
| exact | `"ALL"` | 61 | byte-for-byte ANARCI |
| **pan** (default) | `"pan"` | **240** | 99.2% = ANARCI; IMGT-anchor correctness tied (99.21% vs 99.21%) |

Pan is ~3.9× faster single-thread and **more robust** than 7-species ANARCI (it can't
pick a wrong-species profile — a failure mode that mis-numbers ~5% of sequences when ANARCI
is run with many species). The ~0.8% pan-vs-ANARCI numbering differences are IMGT-legal
gap-placement ties, not errors. Plus the batch path **deduplicates identical sequences**
(lossless; ~1.5× on the test set, more on repetitive NGS data). Validate with
`python scripts/validate_pan.py`.

## Accuracy — byte-for-byte vs reference ANARCI (exact engine)

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

| N | ncpu | ref ANARCI | rs exact | rs pan | pan/ref | pan/exact |
|---|---|---|---|---|---|---|
| 1,000 | 1 | 42.1 | 98.2 | 431.5 | 10.2× | 4.4× |
| 1,000 | 12 | 234.8 | 620.5 | 2334.1 | 9.9× | 3.8× |
| 10,000 | 1 | 38.0 | 855.7 | 3157.7 | 83.2× | 3.7× |
| 10,000 | 12 | 170.5 | 1873.7 | 5441.4 | 31.9× | 2.9× |

### Honest reading of these numbers

**Read `pan/exact` (~3.7–4.4×) as the clean engine speedup** — both run dedup, so it cancels,
leaving the pure 7-vs-29-profile win. The `pan/ref` column is inflated by two dataset-dependent
effects: (1) the benchmark replicates a base set, so N=10000 is ~94% duplicate sequences and
anarci-rs's **dedup** collapses them (reference re-scans every copy) — that's why N=10000 shows
30–83× and N=1000 (~37% duplicate) shows ~10×. On **fully diverse** data dedup does nothing,
and the realistic speedup over reference ANARCI is roughly **~4× (pan engine) × ~1.5× (in-process,
no subprocess/parse) ≈ 6× single-thread**, widening with cores (rayon scales better than ANARCI's
multiprocessing). Dedup is a real, free bonus on top — large on repetitive NGS data, zero on
unique data. (The earlier exact-only figures below — 1.5–1.8× — were measured before the pan
engine and dedup; they still describe the byte-exact `database="ALL"` path on diverse data.)

### (Historical) exact-engine-only reading

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
0% (it's safe with a manual `nmodels` reset, just pointless), and the per-scan `bg`/`sq`
allocations are negligible against the 29-profile DP. The hmmscan→hmmsearch (`SCAN_MODELS`→
`SEARCH_SEQS`) inversion the literature cites as 2–10× applies to *filter-dominated* workloads;
ours is domain-definition-dominated (88%, run on ~all 29 profiles regardless of loop order),
so the inversion's ceiling here is ~1.15× for substantial rework — not worthwhile. The real
lever is **scanning fewer profiles**, which is exactly what the pan engine does.

**Going faster would require changing the algorithm** — e.g. an MSV/SSV pre-filter to skip
clearly-non-matching profiles, or scanning fewer profiles per query. Those change which hits
are reported and would break byte-for-byte parity, so they are out of scope for a drop-in
replacement (a future opt-in `--fast` mode could trade exactness for speed).

## Germline assignment accuracy (`germline_method="evalue"`)

Optional RIOT-style e-value germline assignment is *more accurate* than ANARCI's identity
matching (the one numbering-adjacent axis where "better than ANARCI" is well-defined).
Truth set = the `riot_na` tool (itself an e-value SW assigner — a cross-tool agreement check,
not an independent gold standard), human-only V/J genes (where IMGT and OGRDB names overlap):

| method | V-gene agreement | J-gene agreement |
|---|---|---|
| identity (ANARCI) | 66.8% | 75.9% |
| **evalue** | **75.4%** | **88.5%** |

On V-gene disagreements, the e-value call matches `riot_na` where identity doesn't by ~7:1
(35 vs 5). It changes only `v_gene`/`j_gene` (+ `v_evalue`/`j_evalue`), never the numbering.
Honest limits: the RIOT-paper "97%" is not reproduced — the residual gap is the germline DB
(we use ANARCI's IMGT, riot_na uses OGRDB, whose mouse IDs can't name-match IMGT), and there
is no public curated residue/germline gold standard. Cost: ~10× the identity path per domain,
so it is **opt-in** (the pan engine's default species label uses the fast identity path).

### The other wins (not in the table)
- **Zero runtime dependencies**: HMMER is statically linked and `ALL.hmm` is embedded in the
  extension. No `hmmscan` on PATH, no Biopython, no temp files. `pip install` and `import`.
- **No subprocess fragility**: no `multiprocessing` spawn/`__main__` re-import footguns, no
  temp-file cleanup, no parsing of a text format that can drift between HMMER versions.
