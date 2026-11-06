# Minimum-edit antibody design

Two experimental modes search an explicit substitution space:

- `humanize` finds the fewest edits that satisfy supplied human-reference targets
  and the declared protection/validation rules.
- `autonomy` finds the fewest edits that satisfy a reviewed VH-interface hypothesis.
  Its output is an engineered VH candidate, not a proven VHH or autonomous binder.

Both return a computational certificate with the search space, constraints and
software identities. Neither establishes biological humanization, immune tolerance,
solubility, stability or retained binding.

## A: minimum edits to a human-reference target

```python
from anarci_rs.antibody_design import humanize
from anarci_rs.species_reference import load_profile

result = humanize([{
    "id": "H",
    "sequence": mouse_vh_sequence,
    "domain_format": "VH",
    "profile": load_profile("outputs/species-reference/atlas-natural-v1/human-H-VH.json"),
    "threshold": {
        "kind": "percentile", "minimum": 50.0,
        "rationale": "Illustrative cohort-rank target; replace with a justified study criterion"
    },
    "template": {"v_gene": "IGHV3-23*01", "j_gene": "IGHJ4*02"},
    "mutable_positions": ["12", "13", "17", "19"],
    "core_positions": [],
    "antigen_positions": [],
    "interface_positions": []
}], max_mutations=4, max_seconds=30)
```

This example shows the schema. Its threshold, genes and mutable positions are not
universal biological recommendations, and a feasible result is not guaranteed.
An empty structural-protection list means no such evidence was supplied, not that
those residues were assessed as safe.

A domain requires `id`, `sequence`, `domain_format`, `profile`, `threshold`, and
exactly one of:

1. `allowed_substitutions`: an IMGT-position dictionary of explicit alternative
   residue lists, for example `{"12": ["L", "V"], "19": ["R"]}`; or
2. `template`: one explicitly selected human V gene and J gene from the bundled
   IMGT-aligned germlines. The V gene supplies framework positions before FR4;
   the J gene supplies FR4. Missing template residues are not invented.

The template compiler exposes `template_alphabet(domain, template)` for inspecting
potential substitutions before selecting a mutable subset. There is no automatic
nearest-gene choice or union across templates. A template restricts the edit
alphabet; it does not turn the pooled repertoire likelihood into a germline-family
or IGHV3-23-specific statistical model.

All references are supplied explicitly. Targets use either `raw_score` or
`percentile`, with an inclusive numeric `minimum` and a nonempty `rationale`.
Percentiles require a compatible calibration cohort and full modeled-coordinate
coverage. They use the existing empirical midrank rule exactly, including ties.
They are cohort ranks, not percentages of human identity or immune tolerance.
There is no built-in threshold meaning "humanized enough."

### Protections and paired input

Both modes preserve observed IMGT CDRs, insertions and cysteines by default. They
also protect existing glycine/proline residues and prohibit introducing Cys/Gly/Pro
without the task-specific exceptions described below. Vernier positions from
[Foote and Winter](https://pubmed.ncbi.nlm.nih.gov/1560463/), as tabulated in
[AbNumber](https://github.com/prihoda/AbNumber/blob/master/abnumber/common.py), are
mapped from native Kabat numbering through raw sequence indices. The implementation
pins the reference blob in its scope. It does not substitute Kabat position numbers
for IMGT numbers. Unresolved Kabat/IMGT residue correspondence fails closed.

Supply additional `protected_positions`, `core_positions`, `antigen_positions`
and `interface_positions` from independent evidence. An explicit substitution that
conflicts with a protection is an error, not a silently unlocked position.
Template-derived substitutions are filtered by the declared protection/mutable
policy, and every excluded template position is reported.

For paired humanization, pass a second domain with `domain_format: "VL"`, its own
K/L-compatible human reference and its own threshold. The solver minimizes the
**total** edits across both inputs. Every chain must satisfy its target, guards and
native validation. Missing light-chain references are unsupported; the current
OAS pilot atlas supplies VH references only. Separate chain scores do not constitute
a model of VH/VL coupling or prove that the paired antibody retains function.

Optional guards have the form `{"profile": model, "threshold": {...}}`. They use
the same exact threshold semantics. A guard is another computational constraint,
not an independent experimental assay.

## Search results and limits

| Status | Meaning |
|---|---|
| `parent_meets_target` | Zero edits satisfy every declared check. |
| `proven_minimum` | The returned witness passes all checks; every smaller edit cardinality was exhaustively rejected. |
| `proven_infeasible_within_scope` | Every candidate within the explicit alphabet and edit cap was decisively rejected. |
| `budget_exhausted` | Search was incomplete; no feasibility/minimality conclusion is substituted. |

The search includes zero edits and enumerates complete assignments in increasing
edit cardinality. It does not discard a mutation because its individual score or
guard effect is unfavorable; another mutation may compensate for it. Likewise,
intermediate motif creation does not discard a combination that repairs the motif.
No beam pruning or score-based branch pruning is used. The first accepted witness
is returned, not the highest-scoring witness among equal-size solutions or a list
of all alternatives.

At most 12 mutable positions are supported **in total**, including both chains.
Choose an explicit narrower set when a template produces more positions. Requests
are rejected rather than silently truncated. `max_mutations` is 0 through 12;
the default for humanization is 5. Defaults are 200,000 node visits, 100,000 complete
score screens, 10,000 individual native validators and a 30-second cooperative
search deadline. Integer budgets may be 0 through 500,000, and the deadline may
be at most 300 seconds. These caps deliberately limit resource use, not biology.

Each visited position subset and allele assignment costs a node; zero-edit
assessment therefore costs two nodes. Preflight numbering and model compilation
are separate from the search deadline and reported separately. Inputs are limited
to two sequences of at most 2,000 residues, with one numbered domain each. A native
call already in progress cannot be forcibly interrupted by the cooperative clock;
a possible overrun is reported. This is not a hard wall-clock sandbox.

Raw-sequence acceptance includes exact renumbering: residue indices, boundaries,
chain class and domain count must remain unchanged. Substitutions leave flanks
unchanged. New screened motif occurrences are checked at their raw sequence
locations rather than by total count alone. A successful no-hit or changed numbering
rejects a candidate. Operational exceptions or malformed native results raise
`EvaluationError`; they never become a proof of infeasibility.

The certificate exports `scope_manifest`, `scope_sha256`, per-layer candidate and
rejection counts, `completed_excluded_through`, the witness and validator evidence.
The scope binds parent maps, allele sets, thresholds, profiles, protections,
templates and implementation/native identities. Budget choices are recorded in a
separate run manifest/hash. Numeric decisions use inclusive binary64 comparisons
and authoritative full scores, without hidden epsilon or rounded-score comparisons.
The certificate proves an operational minimum under this evaluator and declared
space; it does not authenticate supplied provenance or prove biological efficacy.

## B: VH autonomy hypotheses

A conventional VH is not relabeled as a VHH. This mode tests a declared hypothesis
about its former VL interface. Its required proxy objective is a ceiling on the
count of `AVILMFWY` residues in a **frozen, explicit interface-position set**. This
simple count is not solvent-accessible surface area, solubility or stability. It
can reward a harmful buried substitution, so candidate selection still needs
structural review and experiments.

Natural-VHH likelihood is optional descriptive context, not an obligatory goal.
[Barthelemy et al.](https://pubmed.ncbi.nlm.nih.gov/18045863/) reported useful
human-VH autonomy substitutions uncommon in natural VH/VHH sequences. A candidate
can therefore pass the declared interface hypothesis even when its VHH score falls.
Every projection is explicitly labeled `query_format: VH`, `reference_format: VHH`,
and `percentile: null`. Missing modeled positions produce no partial-denominator
score. Ordinary species scoring still enforces its original format guard.

```python
from anarci_rs.vh_autonomy import engineer_vh_autonomy

result = engineer_vh_autonomy(
    domain={
        "id": "VH", "sequence": human_vh_sequence, "domain_format": "VH",
        "profile": human_reference,
        "threshold": {"kind": "raw_score", "minimum": selected_human_floor,
                      "rationale": "Predefined research retention floor, not immune-tolerance evidence"},
        "allowed_substitutions": reviewed_alphabet,
        "core_positions": protected_core_positions,
        "antigen_positions": protected_antigen_positions,
    },
    hypothesis="Describe the specific interface-engineering hypothesis",
    hydrophobic_maximum=selected_interface_count_ceiling,
    interface_evidence=reviewed_parent_evidence,
    interface_review={"source_id": "your-reviewed-source", "review_note": "Explain the editable-position rationale"},
    projection_profiles=[vhh_reference], max_mutations=3,
)
```

The input must be the isolated numbered VH sequence, not a full heavy chain with
constant regions. Output remains `engineered_VH_candidate`, with
`autonomous_function_established: false`, `structural_validation: false` and
`experimental_validation: false`, regardless of solver status. An unchanged parent
can meet a loose proxy target; that does not demonstrate existing autonomy.

### Parent structure evidence

`build_interface_evidence` in [autonomy_structure.py](../python/anarci_rs/autonomy_structure.py)
uses native PDB/mmCIF parsing. Supply the actual parsed VH and VL chain IDs and,
if available, explicit antigen-chain IDs. The helper requires exact unique
sequence correspondence and finite C-alpha coordinates. It uses only the VL
variable domain for interface proximity. No gapped sequence alignment, coordinate
imputation, inferred antigen role or automatic biological assembly selection occurs.

```python
from anarci_rs.autonomy_structure import build_interface_evidence

evidence = build_interface_evidence(
    human_vh_sequence, structure_text,
    vh_chain_id="YOUR_PARSED_VH_ID", vl_chain_id="YOUR_PARSED_VL_ID",
    antigen_chain_ids=[], structure_format="pdb",
    provenance={"source_id": "your-structure-source", "structure_kind": "experimental",
                "assembly_context": "asymmetric_unit", "reviewed": True},
)
```

C-alpha proximity describes the supplied parent coordinates; it does not establish
side-chain contacts, exposed surface, a biologically correct assembly or candidate
structural compatibility. Missing antigen chains remain unknown. The native parser
cannot reconstruct unobserved VL/antigen residues or return original author residue
IDs. Provenance and review flags remain caller assertions. Imported evidence is
checked for schema, parent identity and consistency: interface and antigen-protected
sets are rederived from the complete per-position distances, so dropping an antigen
lock alone is rejected. This is consistency checking, not source authentication.

Without structure, the default result is `needs_structure_evidence`. Explicit
`sequence_only=True` permits hypotheses only with a reviewed interface annotation,
`interface_positions`, an edit alphabet and a rationale. The result prominently
reports missing structural evidence. It does not claim a measured interface.

B may deliberately edit the former VL interface, unlike A. Position-specific
`protection_overrides`, each with a review reason, can relax Vernier or Gly/Pro rules
for a declared site. They cannot unlock CDRs, insertions, cysteines, antigen contacts,
core positions or caller immutable locks. No universal camelizing mutation recipe
is applied. Preserving CDR sequence does not preserve its conformation or replace
binding contributions from the original VL.

## JSON CLI

```sh
anarci-species humanize outputs/my-design/humanize-request.json \
  --output outputs/my-design/humanize-result.json
anarci-species autonomy outputs/my-design/autonomy-request.json \
  --output outputs/my-design/autonomy-result.json
```

The humanization request contains `domains` and optional search budgets. The autonomy
request uses the same keyword names as its Python API. Replace profile objects with
file paths relative to the request file; guards and `projection_profiles` use paths
as well. `interface_evidence` can be an evidence JSON path. Alternatively, supply a
`structure` object containing `path` and the keyword arguments of
`build_interface_evidence` to generate evidence before design. Unknown domain options
are rejected. Output files are never overwritten. Exit 0 means the request ran,
not that it found a design or established function; inspect the result status.
Operational evaluation faults exit 3 with `status: evaluation_error`.

## Evidence needed after candidate generation

For A, compare parent and candidates in the same paired antibody format. For B,
compare original paired Fv, isolated parent VH and candidates in matched constructs
and conditions. Predefine acceptance criteria, then measure expression/recovery,
concentration-dependent monomer/aggregation behavior, folding/stability, affinity,
specificity and the relevant functional assay. Independent structure prediction
is triage, not proof of autonomy. Immunogenicity needs a separate assessment.

Synthetic fixtures validate solver mathematics and failure handling only. Real OAS
mouse and IGHV3-23-derived sequences, plus public antibody structures, can validate
software integration; they do not establish function preservation. The current
reference limits remain those documented in [SPECIES_ADAPTATION.md](SPECIES_ADAPTATION.md).
Keep designs, structures and assay/run results in ignored output directories.
