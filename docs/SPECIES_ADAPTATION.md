# Species nativeness and adaptation

Score a numbered antibody against explicitly supplied repertoire references, then
propose conservative substitutions toward a selected reference. This is an
experimental independent-position baseline. It does not estimate affinity,
immunogenicity, organismal fitness or the probability that an antibody came from
a species. No trained species profile is silently bundled or downloaded.

## Score and propose

Install the current package, build the real-data atlas below, then use its profiles:

```sh
anarci-species audit outputs/species-reference/atlas/human-H-VH.json
anarci-species score "$SEQUENCE" --domain-format VH \
  --model outputs/species-reference/atlas/human-H-VH.json \
  --model outputs/species-reference/atlas/mouse-H-VH.json \
  --species human --species mouse --species dog
anarci-species propose "$SEQUENCE" --domain-format VH \
  --target-model outputs/species-reference/atlas/human-H-VH.json \
  --max-mutations 2 --max-candidates 6 --protect 42 49 \
  --output outputs/species-reference/proposals.json
```

The protected positions above illustrate syntax, not a universal protection set.
Supply positions appropriate to the antibody and its experimental/structural
context. All positions are IMGT, and raw sequence indices in results are zero-based.
The CLI refuses to overwrite output files. `python -m anarci_rs.species_cli` is an
alternative to the installed command.

Python exposes the same functions without changing legacy numbering signatures:

```python
from anarci_rs.species_reference import load_profile
from anarci_rs.species_adaptation import profile_sequence, propose_sequence

human = load_profile("outputs/species-reference/atlas/human-H-VH.json")
report = profile_sequence(sequence, [human], domain_format="VH")
proposals = propose_sequence(
    sequence, human, domain_format="VH", max_mutations=2,
    protected_positions=["42", "49"],
)
```

`domain_format` is required: `VH`, `VL` or `VHH`. An unpaired H-domain hit does not
establish VHH identity. `domain_index` selects a domain in a multidomain input;
other domains and unnumbered residues are left unchanged. Native IMGT numbering
uses bundled `ALL` by default, with explicit `pan` as an option. The input sequence
must contain canonical uppercase amino acids. Numbering must map exactly back to
the input residues and inclusive domain endpoints, or the operation fails.

## Build references from real OAS repertoires

The versioned [source manifest](../.agents/species/oas-sources.json) pins complete
public OAS/OAS-nano units, source URLs, study assignments, local SHA-256 checksums
and license evidence. It selects heavy-chain units with explicit study-level evidence of natural,
unselected repertoires. These are real repertoire records, not simulated species data.
The bounded selection is a reproducible pilot, not a representative sample of each
species' entire repertoire. Smaller download size was one selection criterion.

From a source checkout, with Python 3.11+ and the current package installed:

```sh
python3.12 scripts/build_oas_reference.py download \
  .agents/species/oas-sources.json \
  --directory outputs/species-reference/retrieved \
  --output outputs/species-reference/retrieved.json
python3.12 scripts/build_oas_reference.py prepare \
  outputs/species-reference/retrieved.json \
  --sample-size 2000 --seed 1729 \
  --output outputs/species-reference/cohorts.jsonl
python3.12 scripts/evaluate_species_reference.py \
  --records outputs/species-reference/cohorts.jsonl \
  --manifest outputs/species-reference/retrieved.json \
  --output-dir outputs/species-reference/atlas \
  --atlas-id oas-multispecies-pilot --version 2026-11-04-natural-v1
```

Downloads are explicit, bounded, complete compressed units; fetched shell scripts
are never executed. A checksum means the local complete download matches the
pinned bytes, not that the source published a cryptographic certificate. Changed
upstream units fail the checksum rather than silently redefining the reference.
The importer scans complete files and uses a seeded reservoir over eligible rows,
not the first records in each file. The manifest must declare reviewed natural
origin and an unselected library, with a citation and review note. `Disease=None`
and `Vaccine=None` are only metadata filters: they do not rule out synthetic or
panned libraries. The recipe explicitly excludes the Chen_2021 synthetic-library
study and the Rafique_2019/Rafique_2023 immunized or selected-library studies.
Their OAS metadata alone would have admitted them. The selected Jaffe_2022 donor
was also excluded because the primary donor table reports type 1 diabetes despite
OAS `Disease=None`. Source-specific review notes take precedence over a convenient
metadata label. Truncated domains, noncanonical residues,
nonproductive records, stops and inconsistent numbering are excluded with counts.
The observed numbered domain must reconstruct a substring of the source sequence;
terminal overhangs are removed before exact-domain deduplication.

OAS and OAS-nano use different numbering fields (`ANARCI_numbering` and
`ANARCII_numbering`). The importer records which one was used and preserves IMGT
insertion order, including reverse insertion ordering at positions 33, 61 and 112.
It does not infer VHH format solely from the donor species.

Studies are assigned to train, calibration or test in the manifest. Exact domain
duplicates and exact CDR3-plus-V-family proxy overlaps are removed with train-first
priority. This proxy is not a reconstructed clonal lineage, and unknown subject
identities remain explicitly unknown. It does not prove freedom from all sequence
ancestry or pretraining leakage.

### Current reference scope

| Species | Format | Scope and limits |
|---|---|---|
| Human | VH | Several training studies plus distinct calibration/test studies; selected donors and isotypes, not population validation |
| Mouse | VH | Strain labels retained in source metadata; unknown donors in some studies limit evidence |
| Alpaca | VHH | One reviewed natural-repertoire study; training-only, no calibration percentile or held-out performance claim |
| Llama | VHH | One reviewed natural-repertoire study; training-only, no calibration percentile or held-out performance claim |
| Dog | None | No canine entry in the audited OAS index; no score or fabricated surrogate |
| Light chains / paired Fv | None in this atlas | The code can fit separate K/L references, but this acquisition does not validate them or model VH/VL coupling |

OAS's official homepage states CC BY 4.0. The separate OAS-nano redistribution
license was not confirmed during acquisition; keep its raw and derived artifacts
local until that review is resolved. Cite OAS and each underlying study identified
in unit metadata. Do not treat an article license as automatic permission for every
associated dataset. The manifest documents this distinction.

## What the baseline computes

The initial panel uses framework positions 1-26, 39-55, 66-104 and 118-128.
Insertions and CDRs remain in the residue map but are not modeled by the default
framework scorer. In the atlas builder, a position must occur in at least 99% of
training records for each species within a chain/format group. The shared panel is
selected using training data only, before calibration or test scoring. Incomplete
training rows are excluded and counted; no amino acid is imputed. At least 60
positions must remain. Calibration/test records with incomplete coverage are not
assigned a percentile.

Training gives equal total weight to studies, then to donors within each study,
then to lineage proxies within a donor. Unique sequences share their lineage's
weight. The weights sum to the number of training studies. A positive pseudocount
(default 0.5 per amino acid) smooths each position's 20-residue distribution. The
raw score is the mean log probability over the fixed modeled positions. This is a
transparent, strongly regularized baseline; it does not learn interactions between
positions, structural compatibility or mutation epistasis.

The standalone `fit_profile` API instead uses the complete-training intersection
of the requested framework panel unless explicit positions are supplied. Its JSONL
record contract is `record_id`, `species`, `chain`, `domain_format`, `study`, `donor`,
`lineage`, `split`, `sequence` and a `numbering` dictionary. The sequence must equal
the canonical IMGT reconstruction. Study, qualified donor, lineage, raw sequence
and full-numbering hashes may not cross splits. Conflicting records fail instead
of silently passing an audit. Model files retain hashes/provenance, not raw training
sequences, and validate a deterministic content fingerprint on load.

Scores include per-position contributions, modeled coverage and an empirical
midrank percentile from held-out calibration sequences. The CDF is **unweighted
per unique sequence**, not donor-balanced and not an estimate of population immune
tolerance. Full-coverage calibration group counts accompany it. Percentiles from
different references are reference-relative ranks; they are not interchangeable
probabilities. Relative raw-score rankings are omitted unless the position panels,
fit settings and effective smoothing strength match. Scores never sum to 100%.

Evidence labels are conservative support summaries. `repertoire_baseline` requires
at least 60 positions, two training studies, a separate calibration study with at
least 20 full-coverage records, source units declared as reviewed natural/unselected
repertoires and resolved train/calibration group labels. Otherwise the result is `research_only`. Neither label means that provenance
was independently authenticated, a model was biologically validated or a mutation
was shown to preserve function. Sparse or missing species remain unsupported.

## Mutation planner contract

The planner searches a bounded set of substitutions that improve the target model's
raw score. It locks CDRs, insertion positions, cysteines, canonical disulfide positions,
and existing glycine/proline residues. It does not introduce cysteine, glycine or
proline. Caller protections add locks; an optional mutable-position list limits the
search further. These conservative rules do not identify buried residues, Vernier
roles, antigen contacts or chain interfaces: supply those from independent evidence.

Complete candidate sequences are rescored, and the planner retains alternatives
at different mutation budgets. Explicit same-format guard profiles may constrain
score decreases, but are not an independent functional assay. Beam width and evaluation
count bound the search; global optimality and epistasis modeling are explicitly false.
The planner will not silently score a VHH with a conventional human-VH model or
claim to convert a conventional VH into an autonomous VHH.

Raw-sequence proposals verify wild-type residues and renumber every candidate.
A chain-class change, changed domain count or changed residue map fails the numbering
preservation flag. New potential glycosylation/deamidation/isomerization motifs are
reported. `eligible_for_review` means only preserved numbering and no newly flagged
motif; it is not experimental approval. Low-level `propose_numbered` lacks the raw
sequence context and does not perform that renumbering/motif screen.

Every recommendation has evidence level `score_improvement_only` and
`experimental_validation: false`. Binding retention, expression, monomer content,
stability, specificity and immunogenicity require separate measurements. Keeping
CDR sequence unchanged does not guarantee unchanged CDR conformation or binding.

## Validation and next gates

[Reference tests](../tests/python/test_species_reference.py) use deliberately small
numerical fixtures to check arithmetic, leakage rejection and tamper detection.
[Adaptation tests](../tests/python/test_species_adaptation.py) check mutation and
coordinate safety. Those software fixtures are not biological reference data.
The [real-cohort evaluator](../scripts/evaluate_species_reference.py) separately
fits actual OAS records and writes study-held-out score/coverage distributions and
descriptive discrimination metrics. A single selected test study per species does
not establish broad generalization, and no confidence interval treating all
sequences as independent donors is reported.

Before a production humanization/caninization claim, expand donor/study coverage,
resolve dataset rights, compare against germline/grafting and established model
baselines, audit relatedness beyond the current proxy, and prospectively measure
functional retention on independent parental antibodies. Optional contextual models
and paired-chain scoring should earn inclusion by improving those outcomes.

Raw units, prepared cohorts, fitted profiles, evaluation results and mutation
examples belong in ignored output directories. Only the reusable code and pinned
retrieval recipe are distributed, under the [distribution policy](DISTRIBUTION.md).
