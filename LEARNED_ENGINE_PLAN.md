# anarci-rs — Learned Numbering Engine (Design Proposal)

**Status: design only. Not implemented, not greenlit.** This is a spec for a *third*
engine — a learned (ML) numbering model — to sit alongside the existing `pan` and exact
`ALL` HMMER engines. It does **not** replace the byte-for-byte path; that stays the
correctness oracle.

Prompted by: "can we do better than ANARCII — MPS-compatible, more precise, faster?"
The grounded research behind this lives in the session memory `learned-numbering-engine.md`.

---

## 0. Why, and what we're beating

ANARCII (OPIG, BSD-3) replaces ANARCI's HMM with a **tiny autoregressive seq2seq
transformer** (~0.55M params speed / ~1.0M accuracy) that generates IMGT position tokens.
It's distilled from ANARCI labels on ~140–160M OAS chains, so its accuracy *ceiling is
ANARCI* (it claims parity, >99.99% conserved residues / >99.94% CDRs, not superiority).
Its real edge over alignment is **generalization**: VNARs, truncations, junk-flanked NGS
reads, rare species, and version stability. Documented limitations:

- **No Apple MPS path** — `configure_device` hardcodes cuda/cpu; on Apple it runs CPU-only.
- **Autoregressive decode loop** (~128 sequential steps/seq) is the dominant cost and the
  source of its failure modes (it can emit duplicate/non-monotonic/truncated numbering,
  caught only by Python post-processing heuristics).
- Headline **>90k seq/min is an A100 number**; no CPU throughput is reported.

The lever to beat it is a **simpler** model, not a bigger one. Numbering is **monotonic
per-residue labeling**, not translation — each output maps to exactly one input residue,
in order. Dropping the decoder gives three stacked wins:

1. **Speed** — one parallel forward pass instead of a 128-step loop. Measured on this Mac
   (architecture-shape micro-benchmark, random weights, batch 32, len 128): AR decoder
   (1.0M params) = 31 seq/s MPS / 12 CPU; encoder-only classifier (7.5M params) = **1,501
   seq/s MPS / 191 CPU** → 16–48× *despite being 7× larger*. (Note: our pan-HMMER engine
   at ~240–430 seq/s CPU already likely beats ANARCII's AR CPU path.)
2. **Precision** — a CRF/Viterbi head with illegal transitions masked to −∞ makes invalid
   numbering *representationally impossible*, eliminating ANARCII's failure modes by
   construction rather than patching them. (Gains land on the hard/rare tail; headline
   accuracy is already near ceiling.)
3. **MPS + CPU** — deployable in the maturin wheel with no PyTorch (see §6).

**Non-goals:** this engine is an approximation and will *not* be byte-for-byte ANARCI.
The exact `ALL` engine remains the parity oracle and the default for `database="ALL"`.

---

## 1. Where it fits

A new engine selected by `database="learned"`, parallel to `"pan"` and `"ALL"`, same
public API (`number`, `anarci`, `run_anarci`). Numbering output flows through the
*existing* scheme code (`schemes/`) so IMGT→Chothia/Kabat/Martin/Aho/Wolfguy conversion
is reused unchanged — the model only ever predicts the canonical IMGT state, exactly as
the HMM path produces state vectors today.

New components:
- `crates/anarci-ml` — pure-Rust inference (tract or candle), exposes the same
  `HmmEngine`-style interface (`scan(name, seq) -> state vector / DomainDetails`) so the
  orchestrator (`orchestrate.rs`) treats it as a drop-in engine.
- `training/` (Python, dev-only, not shipped in the wheel) — data distillation + training.
- Model weights bundled in the wheel (safetensors or ONNX), like `ALL.hmm` is today.

---

## 2. Problem framing → architecture

**Encoder-only, non-autoregressive, structured token classifier.**

```
amino-acid sequence  ──►  encoder (per-residue hidden states)
                              ├─► emission head ─► linear-chain CRF ─► Viterbi ─► per-residue IMGT label
                              ├─► pooled head ─► chain type (H/K/L/A/B/G/D)
                              └─► pooled head ─► species
```

- **Input / tokenization:** char-level amino acids (20 + X for non-standard), per-residue,
  matching ESM-2's tokenizer if warm-started. Long constructs (scFv, tandem, leader+Fv)
  handled by windowing (as ANARCII does) and/or a "outside-domain" label (see §3).
- **Backbone:** start from **ESM-2-8M** (7.8M params, 6 layers, dim 320, per-residue
  tokenized — drops straight into a token head). Rationale: evolutionary prior helps the
  rare-species / VNAR tail, which is the whole point of going learned, at near-zero size
  cost. Ablation: a from-scratch dim-256/4-layer encoder; escalate to ESM-2-35M only if
  8M underfits the rare tail. *Open: verify ESM-2 op coverage in tract/candle (§6).*
- **Emission head:** linear over the fixed label set (§3).
- **Structured head (the key piece):** linear-chain CRF. Transition matrix entries for
  label pairs that never occur in a valid numbering are set to −∞; decode with Viterbi.
  This guarantees the output is a locally-valid, monotonic, duplicate-free numbering. It
  also handles the IMGT CDR3 insertion *mirror* (111→111A→111B…→112B→112A→112, where the
  raw position number is non-monotonic) correctly, because the legal-transition mask is
  **derived empirically from the oracle's own output ordering**, not from a hand-coded
  "position must increase" rule (which would be wrong for the mirror).
- **Auxiliary heads:** chain-type and species as pooled classifiers off the same encoder
  (ANARCII spent a *separate* 0.5M-param model on chain classification; here it's free).

---

## 3. Label scheme (the part that needs care)

Each input residue gets exactly one label from a **finite, known** set:

- `OUTSIDE` — residue is not part of a variable domain (leader, tag, linker, trailing junk).
  The numbered domain is the maximal run of non-`OUTSIDE` labels → domain boundary
  detection falls out of the same head, no separate model.
- One label per **valid IMGT position-with-insertion** that actually occurs. We do **not**
  enumerate these from memory or from IMGT rules — we **derive the label set empirically**
  by running the exact engine over the training corpus and collecting every distinct
  `(position, insertion_letter)` pair observed. This is grounded (it's whatever ANARCI
  actually emits) and self-correcting for rare insertion sites. Expected ~150–200 classes.

Notes:
- **Deletions** need no label — a deleted IMGT position is simply one that no residue is
  labeled with (the position is absent from the output). Only *insertions* and *presence*
  are predicted.
- **Multi-domain** (scFv/tandem): `OUTSIDE`-separated runs. Whether one pass over a long
  construct suffices or we window per ANARCII is an **open question** to settle empirically
  — repeated position labels across two domains stress the CRF's local transition model.
  Start with windowing (lower risk), measure, then consider single-pass.
- **CRF transition mask** = the set of adjacent-label pairs observed in the oracle-numbered
  corpus, plus `OUTSIDE↔OUTSIDE`, `OUTSIDE→first-positions`, `last-positions→OUTSIDE`.
  Everything else → −∞. Validity is then a *structural property*, testable by fuzzing.

---

## 4. Training data — distill from our own exact engine

We are the oracle ANARCII had to borrow from ANARCI. The exact `ALL` engine is
byte-for-byte ANARCI, so it generates ground-truth labels directly.

- **Sequences:** OAS (as ANARCII used; check OAS terms), plus our existing example FASTAs,
  plus targeted rare-tail sets. Label every sequence with the exact engine (IMGT).
- **Coverage:** all chains (H/K/L + TCR A/B/G/D), all training species. Deliberately
  **enrich the hard tail** — rare insertions (e.g. 112-region), non-human/mouse species,
  truncations, N/C-terminal junk, VNARs — because that's where a learned model earns its
  keep and where uniform OAS sampling under-represents.
- **Augmentation:** synthetic truncations, random leader/tag/linker flanks, conservative
  point mutations — to teach the `OUTSIDE` label and robustness ANARCII gets from NGS noise.
- **Splits:** train/val/**held-out test by species and by germline family** (no leakage).
  A separate **rare-tail test** (VNAR, heavy truncation) measures the generalization claim.
- **Licensing:** OAS redistribution terms must be checked before bundling any derived data;
  the *model weights* are derived parameters, not the sequences. Flag for review.

---

## 5. Training procedure

- **Loss:** CRF negative log-likelihood (numbering) + cross-entropy (chain, species),
  weighted. Optionally a label-smoothing term on emissions.
- **Warm-start:** ESM-2-8M encoder weights; train head + fine-tune encoder.
- **Decode:** Viterbi (exact, deterministic). No sampling, no beam.
- **Stop criterion:** the §7 validation gates, not a fixed epoch count.
- **Hardware:** training on CUDA (or MPS for the small ablation); inference targets below.

---

## 6. Inference & deployment (MPS + CPU, no PyTorch in the wheel)

Ranked from the runtime research:

- **Primary: tract** (pure-Rust, ARM-tuned SIMD, deterministic, tiny runtime, int8).
  Treat as **CPU-only** (its Metal crate exists but is immature). Export path:
  PyTorch → ONNX → tract, optimized at build time. Verify ESM-2 + CRF op coverage on
  import (CRF Viterbi may be implemented in Rust outside the ONNX graph — likely cleaner).
- **Accelerator: candle** (pure-Rust, CPU + *opportunistic* Metal, loads safetensors
  directly, no conversion). Metal works but isn't fully mature → validate, don't trust blind.
- **Determinism policy:** **CPU is the canonical/reference path.** Float matmul order
  differs across devices, so an argmax/Viterbi tie could flip CPU↔MPS. GPU/MPS is a
  *validated-equal accelerator*, never the oracle. Document a tie-break rule.
- **Rejected:** PyTorch-MPS (the dep we removed), Core ML (Apple-only, breaks the
  cross-platform wheel), ggml/gguf (built for decoder-only LLMs, wrong architecture).
- **Quantization:** int8 dynamic on the CPU path (~2–4× typical, negligible accuracy loss
  at this scale); skip on MPS (patchier support).
- **Batching:** large batches dominate throughput on both backends — expose batch size in
  the engine, default conservative, document the GPU sweet spot.

---

## 7. Validation gates (go/no-go — these decide whether it ships)

Measured against the **exact engine as oracle** on held-out data:

| Gate | Threshold |
|---|---|
| Numbering agreement, conserved residues | ≥ 99.9% (match ANARCII parity) |
| Numbering agreement, CDR residues | ≥ 99.5% |
| Output structural validity (monotonic, no duplicates) | **100% by construction** — proven by CRF mask + fuzz |
| Rare-tail recall (VNAR / heavy truncation the exact engine *can't* number) | beats exact engine by a measured margin (the generalization win) |
| Speed, CPU single core, batched | ≥ 200 seq/s (≈ pan engine; ≫ ANARCII AR-CPU) |
| Speed, Apple MPS | ≥ 1,000 seq/s (target: beat ANARCII's A100 number on a laptop) |
| Determinism, same device | bit-identical across runs |
| Determinism, CPU vs MPS | ≥ 99.99% label agreement, documented tie policy |

**Kill criteria:** if conserved-residue agreement can't clear ~99.9% with ESM-2-35M, or if
structural validity can't be guaranteed, the learned engine doesn't ship — the HMMER
engines already cover the exact + fast cases, so a learned engine only earns its place by
winning the *generalization* tail without regressing the common case.

---

## 8. Phases & rough effort

| Phase | Work | Output |
|---|---|---|
| P0 | This design | ✅ (this doc) |
| P1 | Data pipeline: distill labels from exact engine, derive label set + transition mask, build splits | reproducible dataset + label spec |
| P2 | Model + training (PyTorch), feasibility metric on a slice (human H+L) | go/no-go on accuracy vs oracle |
| P3 | Export + Rust inference (`anarci-ml`, tract/candle), CRF Viterbi in Rust, determinism harness | working in-process engine |
| P4 | Integrate as `database="learned"`, run §7 gates, benchmark vs pan/exact/ANARCII | gated engine + BENCHMARKS update |
| P5 | Docs (README engine table, this plan → results) | shipped or documented-dead |

Realistic effort: **weeks**, dominated by P1–P2 (data + training + validation), not the
Rust integration. P2 is the gate — build the feasibility slice before committing to P3+.

---

## 9. Open questions / risks

- **Does a ~8–35M encoder actually match ANARCII's 99.99%?** Hypothesis from task structure,
  not a result. P2 feasibility answers it.
- **ESM-2 + CRF op coverage** in tract (ONNX) / candle — must verify before P3; CRF Viterbi
  likely lives in hand-written Rust, not the graph.
- **Multi-domain** (scFv/tandem) — windowing vs single-pass, unresolved (§3).
- **Cross-device float determinism** — inherent to any float model; mitigated by CPU-oracle.
- **OAS licensing** for any redistributed derived data (§4).
- **Label-set completeness** — empirically derived, so a rare insertion unseen in training
  would be unlabelable; mitigated by rare-tail enrichment + an explicit "unknown insertion"
  failure (no silent fallback).

