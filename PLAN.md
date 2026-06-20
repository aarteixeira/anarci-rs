# anarci-rs — Implementation Plan

A Rust reimplementation of [ANARCI](https://github.com/oxpig/ANARCI) (antibody/TCR
variable-domain numbering) with a Python wrapper. Goal: **byte-for-byte drop-in
replacement** (`import anarci_rs as anarci`) with all internals in Rust, the batch
path fully parallel in Rust, and maximum performance — while reproducing reference
ANARCI's output exactly.

Reference pinned: **conda `anarci 2024.05.21`**, HMMER 3.4, against which all
correctness is measured.

---

## 0. Why this shape — the decisive measurement

Profiling stock ANARCI (1000 seqs, IMGT, single-thread, 27.8s total):

| Component | % of time |
|---|---|
| `hmmscan` subprocess (fork + `select.poll` wait) | **67%** |
| Biopython text parsing of hmmscan output (`SearchIO`) | **31%** |
| All numbering schemes (`number_*`) | 1.3% |
| State-vector build + germline assignment | <0.5% |

**The entire win is (1) running HMMER in-process and (2) deleting the text parse.**
The numbering logic is negligible to *run* but must still be ported exactly so the
batch loop is pure Rust with no Python per-sequence overhead and no GIL.

Confirmed enablers (verified this session):
- `ALL.hmm` has `RF yes`; `hmmscan` emits exact `RF` + `PP` alignment lines → a native
  HMMER-C FFI can read the precise alignment columns numbering depends on.
- `pyhmmer 0.12.1` runs the scan in-process and exposes `posterior_probabilities`
  **byte-identical** to the CLI, same hits/scores → perfect test-time oracle.
- HMMER 3.4 + Easel are **BSD-3** → can be vendored and statically linked into the wheel.

## Engine decision

**Native HMMER 3.4 FFI, statically bundled, driven in-process via `p7_Pipeline`,
parallelized across sequences with rayon.** Reached via a phased path that first ships
a pyhmmer-backed milestone (to prove the numbering port in isolation), then swaps in the
native engine. pyhmmer + stock `hmmscan` are **test-time oracles only — never a runtime
fallback.**

---

## 1. Repository layout

```
anarci-rs/
├── Cargo.toml                 # cargo workspace
├── pyproject.toml             # maturin build (module: anarci_rs)
├── crates/
│   ├── anarci-core/           # pure Rust: types, state-vector, ALL numbering schemes,
│   │                          #   germline assignment. No FFI, no Python. 100% unit-testable.
│   ├── hmmer-sys/             # thin -sys FFI: vendors HMMER3.4+Easel (pinned), build.rs
│   │                          #   static-links libhmmer.a + libeasel.a; bindgen exposes
│   │                          #   p7_* / P7_ALIDISPLAY incl. rfline/ppline.
│   ├── anarci-hmm/            # safe wrapper over hmmer-sys: scan pipeline, alignment
│   │                          #   extraction → state vectors, rayon batch.
│   └── anarci-py/             # PyO3 extension `anarci_rs`: mirrors ANARCI's public API.
├── python/anarci_rs/          # packaged data (HMMs, germlines) + thin __init__ re-exports.
├── data/                      # pinned ALL.hmm(+pressed) + germlines (generated Rust/bincode).
├── reference_data/            # (exists) ground-truth JSON + pinned reference artifacts.
├── tests/                     # golden numbering + end-to-end parity + differential.
├── benches/                   # criterion micro-bench + end-to-end throughput harness.
├── PLAN.md / README.md / BENCHMARKS.md
```

## 2. Public API parity (the drop-in surface)

Mirror exactly, with identical names, signatures, defaults, and return shapes:
- `anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None, csv=False,
  allow={...}, hmmerpath="", ncpu=None, assign_germline=False,
  allowed_species=['human','mouse'], bit_score_threshold=80)`
  → `(numbered, alignment_details, hit_tables)`
- `run_anarci(seq, ncpu=1, **kwargs)` → `(sequences, numbered, alignment_details, hit_tables)`
  (accepts list of tuples, a FASTA path, or a single string; batch path is **all-Rust+rayon**,
  not Python multiprocessing).
- `number(sequence, scheme="imgt", ...)` → `(numbering, chain_type)` or `(False, False)`.
- Plus the lower-level functions code may import: `run_hmmer`, `parse_hmmer_output`,
  `number_sequence_from_alignment`, `number_sequences_from_alignment`,
  `run_germline_assignment`, `get_identity`, `validate_sequence`, `read_fasta`, and a
  `schemes` submodule exposing every `number_*`.

Return values are Python tuples/lists matching ANARCI's exact nesting and `None`/`False`
sentinels (validated against the reference JSON).

## 3. Phased implementation, each phase with a hard verification gate

### Phase 0 — Scaffolding & data pinning
- Cargo workspace + maturin `pyproject.toml`; `import anarci_rs` builds (empty stubs).
- Vendor pinned data; assert sha256 (`ALL.hmm`=`cdb77c…`, `germlines.py`=`025ce6…`).
- Convert `germlines.py` (`all_germlines` dict) → a Rust-loadable form (codegen or bincode
  `include_bytes!`), with a round-trip test proving the Rust structure equals the Python dict.
- **Gate:** `cargo build` + `maturin develop` succeed; data checksums match; germline
  round-trip test passes.

### Phase 1 — Pure-Rust numbering core (`anarci-core`)
- Port verbatim from `schemes.py` per the extracted spec: constants (`alphabet` 53-entry,
  BLOSUM62 asymmetric, all 128-char `state_string`/`region_string`, `rels`, `exclude_deletions`,
  AHo/wolfguy/CDR3 tables), `smooth_insertions`, `_number_regions` (with in-place `rels`
  mutation + the `state_type!="i" or region is None` region guard), `get_imgt_cdr`,
  `gap_missing`, `get_cdr3_annotations`, `_get_wolfguy_L1`, and all ten `number_*` functions.
  Replicate known quirks (martin-heavy FW3 dead branch behavior; wolfguy un-gapped output;
  inclusive `end_index`).
- Port `_hmm_alignment_to_states` (incl. N/C-terminal extension heuristics), `_parse_hmmer_query`
  (domain dedup, bit-score threshold, species limiting), `get_hmm_length`, `check_for_j`,
  `run_germline_assignment`, `get_identity`.
- **Test harness:** instrument reference ANARCI to dump
  `(state_vector, sequence, scheme, chain_type) → expected_numbering` for the full test set;
  the Rust core must reproduce `expected_numbering` **byte-identical**.
- **Gate:** 100% match on golden numbering tests, all 6 schemes × all test sequences. Any
  mismatch is a hard failure to be root-caused (no tolerance).

### Phase 2 — In-process HMM engine
- **2a (de-risk the pipeline):** feed the Rust core alignments from pyhmmer (oracle) and
  produce full `anarci()` output; assert equality with reference `(numbered, alignment_details,
  hit_tables)` JSON byte-for-byte. Proves everything except the native engine.
- **2b (native engine):** `hmmer-sys` vendors HMMER 3.4 + Easel (git submodule pinned to tag
  `hmmer-3.4`); `build.rs` runs configure/make for static libs; bindgen allowlist incl.
  `p7_*`, `P7_ALIDISPLAY`, `esl_alidisplay_*`. Safe `anarci-hmm` wrapper: `impl_Init()` +
  `p7_FLogsumInit()` once; open pressed `ALL.hmm`; per-thread `P7_PIPELINE`+`P7_BG`; loop
  profiles `ReadMSV` → **set `pli->hfp`** → `p7_Pipeline` → extract
  `P7_ALIDISPLAY{aseq,model,mline,ppline,rfline,hmm/sq coords}` + bitscore/evalue → state vector.
- **Gate (engine exactness):** for every test sequence the native engine's per-domain
  RF+PP+aseq+coords and bitscore match stock `hmmscan`/pyhmmer **byte-identical**, and the set
  of hits over threshold is identical. Divergence = hard error, investigated, **never** routed
  around with a subprocess fallback.

### Phase 3 — Batch parallelism + PyO3 (`anarci-py`)
- rayon over sequences; `Python::detach` to drop the GIL around the batch.
- Wire the exact public API (Section 2) with exact error semantics: `AssertionError` where
  ANARCI raises (e.g. unimplemented scheme/chain — incl. the chothia/kabat/martin/wolfguy + TCR
  case, with the CLI's `allow`-set workaround replicated), `(False,False)` where `number()`
  catches, `<70aa` short-circuit, unknown-amino-acid validation. **No silent default ever.**
- **Gate:** end-to-end `anarci_rs` output == reference across all schemes + negative controls
  (lysozyme→None), multi-domain scFv, long CDR3, TCR α/β.

### Phase 4 — Benchmark (speed + accuracy) — final deliverable
- **Accuracy:** run the full `antibody_sequences.fasta` (997 seqs) × all 6 schemes with
  `assign_germline=True`; assert 100% identical to reference ANARCI; report mismatch count
  (target: 0). Add a larger external set if available.
- **Speed:** criterion micro-benches (numbering, germline) + end-to-end throughput at
  N=1k/10k/100k, ncpu sweep 1→12, vs reference baseline (~43 seq/s single, ~270 @12).
  Report speedup + scaling curve (expect near-linear, GIL-free).
- Deliver `BENCHMARKS.md` (tables + methodology).

### Phase 5 — Packaging & docs
- maturin wheel bundling HMMER static + data; `pip install` then `import anarci_rs as anarci`
  works with **zero system dependencies** (no installed `hmmscan` needed).
- `README.md`: drop-in usage, exact-parity statement, benchmark summary, build-from-source notes.

## 4. Correctness ladder (how exactness is guaranteed in all cases)

1. **Golden numbering tests** — Rust reproduces reference numbering from identical state
   vectors (isolates the scheme port).
2. **Engine exactness** — native FFI alignment == stock hmmscan/pyhmmer (RF/PP/coords/bitscore
   byte-identical).
3. **End-to-end parity** — full `anarci_rs` output == reference JSON across all schemes,
   edge cases, and negative controls.
4. **Differential testing** — large real + generated sequence sets cross-checked vs reference;
   any single mismatch fails the suite.
5. **Data pinning** — both implementations use the *same* `ALL.hmm` + germlines (checksummed),
   so drift is impossible.

## 5. No silent failures / no silent fallbacks (explicit policy)

- Every fallible operation returns `Result`/raises; **no defaulting, no best-effort coercion.**
- Replicate ANARCI's *own* explicit errors exactly (assertion messages preserved), including
  quirks — drop-in means matching behavior, not "fixing" it.
- If the native engine ever disagrees with the oracle, that is a **build-time hard error** to be
  fixed, not a runtime fallback to a subprocess.
- Unknown residues, too-short sequences, and no-hit sequences follow ANARCI's exact explicit
  behavior (raise / `None` / `False`), never silently substituted.

## 6. Risks & mitigations

| Risk | Mitigation |
|---|---|
| HMMER FFI build (autoconf, SSE, macOS arm64) | Start from `libhmmer-sys`'s proven `build.rs`; pin tag `hmmer-3.4`. Phase 2a (pyhmmer-backed) is a *delivery* milestone that ships value even if FFI slips — this is project sequencing, **not** a runtime fallback. |
| Exact alignment reproduction | HMMER C *is* the reference, so exact by construction; remaining risk is only extraction correctness, caught by the Phase-2b oracle gate. |
| `germlines.py` → Rust fidelity | Generate + checksum + round-trip equality test vs the Python dict. |
| Multi-domain / `check_for_j` corner cases | Explicit scFv + long-CDR3 golden cases. |
| PyO3 return-shape fidelity (tuples vs lists, sentinels) | Deep-equality assertions vs reference JSON. |

## 7. Deliverables

- `anarci-rs` workspace; pip-installable wheel; drop-in `import anarci_rs as anarci`,
  zero system deps.
- Full test suite (golden + end-to-end + differential) green.
- `BENCHMARKS.md` (speed + accuracy), `README.md`, this `PLAN.md`.
