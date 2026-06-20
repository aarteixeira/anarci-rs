# anarci-rs — Partial-Chain & Multi-Domain Annotation (Design Proposal)

**Status: design only. Not implemented, not greenlit.** Plan for two requested features:
1. Identify & annotate **partial chains** (e.g. FR3-CDR3-FR4-only, or CDR2-onward fragments).
2. Automatically identify & annotate **multiple V regions** in one polypeptide (scFv,
   tandem scFv, and complex constructs like scFv-Fab-Fc-VHH).

The investigation behind this (code read + live tests in the reference env) is grounded
at the `file:line` references below. **Important framing: a large part of both features
already works today** — this plan is mostly about *annotation richness* and a few
robustness fixes, not building detection from scratch. See also [PLAN.md](PLAN.md) (the
HMMER engines) and [LEARNED_ENGINE_PLAN.md](LEARNED_ENGINE_PLAN.md) (the learned engine,
which handles truncations more naturally).

---

## 0. What already works (verified, do not rebuild)

**Multi-domain (Feature 2 core):** the engine already returns *all* non-overlapping V/TCR-V
domains, not just the best. The HMM scan emits one `Hsp` per reported profile-domain
(`crates/anarci-hmm/src/lib.rs:180`, `collect_hsps`); `parse_hmmer_query` sorts by e-value,
drops overlapping redundant hits, and orders survivors **N→C by query_start**
(`crates/anarci-core/src/align.rs:66, 187-233`); `process_one` numbers each domain whose
chain type is allowed and returns a **list of numberings + a parallel list of per-domain
detail dicts** per sequence (`crates/anarci-core/src/orchestrate.rs:246-289`; Python schema
`crates/anarci-py/src/lib.rs:270-297`). Non-V regions (constant domains, Fc, linkers) are
ignored because no V/J profile hits them. **Verified live**: scFv → 2 domains (H, K);
tandem VH-VL-VH → 3 (H, K, H); VH-constant-VL → 2 V domains, constant region correctly
skipped; linkers left unnumbered, never mis-numbered. So the V-domain enumeration of
scFv-Fab-Fc-VHH (VH, VL, VHH — VHH reported as chain_type H) **already happens.**

**Partial chains (Feature 1 core):** fragments are numbered with correct IMGT positions for
the covered region, gated only by bitscore (default 80). There is **no** conserved-residue
gate, no coverage gate. **Verified live** (full VH truncated):
FR3-CDR3-FR4 (IMGT 66→end) → numbered H, IMGT 66–128; FR2-onward (IMGT 39→end) → IMGT
39–128. HMMER local entry/exit is what allows the partial profile match
(`align.rs:90, 106` are only small terminal pads, not gates). Acceptance logic =
bitscore ≥ threshold (`align.rs:194`) + non-empty state vector + chain in `allow`
(`orchestrate.rs:249`) + `validate_numbering` contiguity check.

So the genuinely-new work is narrower than the asks sound. Below, each feature separates
**already works** / **gaps** / **implementation**.

---

## 1. Feature 1 — partial chains + region annotation

### Gaps
- **No explicit region annotation.** The output is a flat list of `(IMGT_position, ins), aa`.
  Which regions (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4) are present, and whether each is
  complete or partial, is left for the caller to infer from the numbers. The detail dict
  has no region field (`crates/anarci-py/src/lib.rs` — zero region references).
- **`number()` rejects fragments < 70 aa** (`orchestrate.rs:430`; parity with ANARCI
  `anarci.py:980`, comment: "ANARCI can number fragments of chains well"). `anarci()` /
  `run_anarci()` have *no* length gate, so partials already work there — only the
  single-sequence wrapper blocks short ones.
- **rs engine scores ~3–25 bits below reference** → borderline short partials that ANARCI
  accepts at threshold 80 are dropped by rs (verified: 50-aa FR1-CDR1-FR2 frag, ref 82.0
  accept vs rs 78.3 reject; at threshold 10 rs numbers it *identically* to ref). This is a
  calibration artifact, not a policy difference — but it's load-bearing for borderline
  partials.
- **Rejections are silent-ish** — a failed fragment returns bare `None`; the reason
  (bitscore X < threshold T) is only recoverable from the hit_table. Violates the
  no-silent-failure rule for a feature whose whole point is handling marginal inputs.

### Implementation
- **Region-completeness annotation (the main new piece) — engine-agnostic post-pass.**
  The IMGT region boundaries already exist internally as a `region_string`
  (`crates/anarci-core/src/schemes/mod.rs:432-433`; decoded: FR1 1–26, CDR1 27–38,
  FR2 39–55, CDR2 56–65, FR3 66–104, CDR3 105–117, FR4 118–128). After numbering, bucket
  each domain's residues into the 7 regions and report per region: **absent** / **partial**
  / **complete**, plus the covered IMGT span `[min,max]`.
  - FR spans are fixed-width → "complete" = full span present.
  - CDRs are variable-length → "complete" = the conserved anchors bounding the CDR are
    present (e.g. CDR3 complete iff Cys104 and Trp118 both present), since CDR insertion
    length is sequence-dependent. Do **not** fake completeness from a fixed width.
  - Surface as an additive detail-dict field, e.g.
    `regions: {FR1:"complete", CDR1:"complete", ..., FR3:"partial", CDR3:"absent"}` +
    `covered_imgt:[min,max]`. Additive → parity-safe (existing keys unchanged).
  - Works identically on HMM-engine or learned-engine output (it consumes `Numbered`).
- **Parameterize the `number()` length gate** — replace the hard `70` with a `min_length`
  arg (default 70 for parity). Don't silently lower it; make rejection a caller choice.
- **Surface rejection reasons** — when `numbered is None` but a sub-threshold hit exists,
  return/log "best bitscore X < threshold T (chain C)" so a partial that just missed is
  diagnosable. No new fallback, no auto-lowering — an explicit, reported decision.
- **Calibration note** — document (and ideally measure) an rs-appropriate default threshold
  or a per-engine offset so partials accepted by reference aren't dropped by rs at the same
  nominal value. Track separately from this feature; it's a pre-existing engine divergence.

---

## 2. Feature 2 — multiple V regions + full-construct annotation

### Gaps
- **No constant-domain / Fc / hinge identification.** Only V/J HMMs exist
  (`reference_data/dat/HMMs/`: FEW.hmm = 7 pan V/J profiles, ALL.hmm = 29 species×chain V/J).
  Constant regions are skipped, not labelled. ANARCI's own source flags constant-domain ID
  as explicit future work (`anarci.py:321`).
- **No VHH/nanobody flag** — a VHH is reported as `chain_type=H`, indistinguishable from a
  conventional VH.
- **No domain-class taxonomy or construct topology** — `details` carries `chain_type` ∈
  {H,K,L,A,B,G,D} only. Nothing says "VH vs CH1 vs Fc", no Fab/Fv grouping, no overall
  N→C construct map.
- **Coarse overlap resolution** — `domains_are_same` is a pure interval-overlap test
  (`align.rs:66`); any overlap merges. Fine for V-only, but a V and a (future) constant hit
  spanning a junction would wrongly drop one.
- **`number()` is single-domain** — returns `nums[0]` only (`orchestrate.rs:446-449`).
  Fine for the drop-in `number()` contract; note it for any new multi-domain convenience API.

### Implementation
Two tiers — ship the cheap one first; the constant-domain tier is the big lift.

**Tier A — annotate the V domains we already find (cheap, high value):**
- Add a `domain_class` field (`"V"` for now) and keep the N→C list = construct topology.
- **VHH detection**: post-numbering heuristic on H domains using the hallmark VHH framework-2
  substitutions (IMGT 37/44/45/47); set `domain_class:"VHH"` (or a `vhh:true` flag). This is
  a flag on an already-found/numbered domain — no HMM change. Validate against known
  VHH/VH sets; surface confidence, no silent guess.
- **`annotate()` convenience function** (optional, beyond the drop-in surface): returns a
  structured per-domain view (class, chain, species, covered span, regions from §1) +
  construct summary (e.g. `"VH→VL"` for an scFv). Additive; doesn't touch the ANARCI-shaped
  functions.

**Tier B — full construct topology incl. constant/Fc (large, optional):**
- **Add constant-domain HMM profiles** (CH1, CL-κ/λ, CH2, CH3, hinge). This is the
  load-bearing change. The scan/collect loop is already N-profile-generic; but:
  - update the hardcoded model count `Z` (`anarci-hmm/src/lib.rs:48`) — note it's already
    questionable for FEW.hmm (29 vs 7); fix as part of this.
  - constant profiles must follow the `species_chain` id convention or `split_hit_id`
    panics (`align.rs:61-63`).
  - make overlap resolution **class-aware** — a V and a C hit may legitimately coexist; only
    same-class hits should dedup.
  - in `orchestrate.rs`, do **not** run constant domains through the V-specific numbering
    (`number_sequence_from_alignment`); record them as annotated spans with `numbered:None`
    but populated `query_start/end`/class. `check_for_j` already no-ops for multi-domain.
  - recalibrate e-value/threshold logic when the profile set changes.
- **Output schema**: `domain_class ∈ {"V","VHH","C","Fc","hinge"}`; constant domains carry
  `numbered:None` + span + class. Construct topology = N→C list order.
- IMGT *does* define a constant-domain numbering scheme; full C-numbering is a further
  scope — MVP is identify + label spans, not number them. Flag as a decision point.

**Learned-engine path** (see LEARNED_ENGINE_PLAN.md): downstream schema (domain_class,
regions, spans) is identical. Detection differs — domain calls come from the model's
per-residue classification, and constant/Fc/hinge become additional learned label classes
(an extension of the `OUTSIDE` label), with no `Z`/e-value notion. Partial fragments are
handled naturally (no bitscore floor; CRF keeps numbering monotonic) — the learned engine
is the better long-term home for the marginal-partial and VHH-tail cases.

---

## 3. Shared output schema (additive, parity-safe)

Per-domain detail dict gains (existing keys unchanged):
- `domain_class`: `"V" | "VHH" | "C" | "Fc" | "hinge"` (V-only until Tier B).
- `regions`: `{FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4 → "absent"|"partial"|"complete"}`.
- `covered_imgt`: `[min_position, max_position]`.

Sequence-level (optional, via new `annotate()`): N→C ordered domain list + a topology
string (e.g. `"VH→VL"`, `"VHH"`, eventually `"VH-CH1 / VL-CL / Fc / VHH"`).

The ANARCI-shaped `number()`/`anarci()`/`run_anarci()` return values keep their current
structure; new fields are added inside the per-domain dicts only, so existing callers and
byte-parity checks are unaffected.

---

## 4. Validation gates
- **Multi-domain (already works)**: keep the existing scFv/tandem/constant-between-V tests
  green vs reference ANARCI (2/3/2 domains, correct chain types, N→C order).
- **Region annotation**: on full-length domains every region reports `complete`; on
  controlled truncations (FR3-CDR3-FR4, FR2-onward, N-term-only) the present/partial/absent
  map matches the truncation by construction; verify against the IMGT boundaries at
  `schemes/mod.rs:432`.
- **Partial robustness**: parameterized `min_length` and surfaced rejection reasons covered
  by tests; a rejection always carries a reason (no silent `None`).
- **VHH detection**: precision/recall on a labelled VHH vs VH set; report confidence.
- **Constant/Fc (Tier B)**: per-class identification accuracy on known Fab/Fc/full-mAb
  constructs; V-domain numbering must remain byte-identical to today (constant profiles must
  not perturb V-domain calls).
- **No silent failures anywhere** — every rejected domain/region/fragment surfaces a reason.

---

## 5. Phases & effort
| Phase | Work | Risk |
|---|---|---|
| F1a | Region-completeness annotation (post-pass) + schema field | low — boundaries already exist |
| F1b | Parameterize `min_length`, surface rejection reasons | low |
| F1c | rs-vs-ref bitscore calibration (separable, pre-existing) | medium |
| F2a | `domain_class` field + VHH heuristic + optional `annotate()` | low–medium |
| F2b | Constant/Fc/hinge HMMs, class-aware overlap, `Z` fix, recalibration, schema | **high** — needs constant-domain reference data + recalibration |

F1a + F1b + F2a are small and deliver most of the user-visible value on the V domains we
already find. F2b (full construct topology including constant regions) is a separate, larger
effort gated on obtaining/building constant-domain profiles — schedule it on its own.

---

## 6. Open questions / risks
- **Constant-domain reference data** — where do CH1/CL/CH2/CH3/hinge HMMs come from
  (build from IMGT C-domain alignments? reuse an existing tool's)? Licensing + provenance.
- **Number constant domains, or just label spans?** IMGT has a C-domain scheme; full
  C-numbering is extra scope. MVP = identify + label.
- **VHH heuristic reliability** — hallmark FR2 residues are indicative, not definitive;
  decide a confidence threshold and whether to ever assert vs report-uncertain.
- **Threshold recalibration** when adding profiles must not regress the byte-exact V-domain
  numbering of the `ALL` engine.
- **Overlap policy** for genuinely overlapping engineered domains — current pure-interval
  test may drop real domains; class-aware + score-based tie-break needed for Tier B.
