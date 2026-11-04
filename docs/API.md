# Python API

## Numbering

`anarci_rs.anarci(records, ...)` accepts `(identifier, protein)` pairs and returns
`(numbered, alignment_details, hit_tables)` in input order.

`anarci_rs.run_anarci(source, ...)` also accepts a FASTA file or a single protein
string. It returns `(records, numbered, alignment_details, hit_tables)`.

`anarci_rs.number(sequence, ...)` returns the first domain's positions and chain
class, or `(False, False)`. Its additional `min_length=70` and
`bit_score_threshold=80` parameters control rejection of short or weak matches.
The extended function warns when a valid sequence has no accepted domain.

Batch parameters are `scheme="imgt"`, `database="pan"`, `output=False`,
`outfile=None`, `csv=False`, `allow=None`, `hmmerpath=""`, `ncpu`,
`assign_germline=False`, `allowed_species=None`, `bit_score_threshold=80`,
`germline_method=None`, and `annotate_regions=False`. `anarci` defaults to one
worker when `ncpu=None`; the extended `run_anarci` defaults to `ncpu=1`.
`ncpu=0` uses the Rayon global worker pool. Worker-pool creation failures raise
an error. Negative worker counts are rejected.

`database="pan"` selects seven pan-species profiles, with evalue germline
assignment by default. `database="ALL"` selects 29 reference profiles, with
identity assignment by default. A custom HMM path uses native loading. Explicit
`germline_method="identity"` or `"evalue"` overrides the default.

Schemes are `imgt`, `chothia`, `kabat`, `martin`, `aho` and `wolfguy`. Short
aliases are `i`, `c`, `k`, `m`, `a` and `w`. Chain identifiers are `H`, `K`, `L`,
`A`, `B`, `G` and `D`; antibody-only schemes reject incompatible TCR domains.

`allowed_species` restricts profile/germline selection. Empty or `None` means
all species in extended batch calls. With pan profiles, reported species are
derived from germline matching. E-values use the actual loaded profile count.
`annotate_regions=True` adds a `regions` dictionary with `complete`, `partial`
or `absent` values and `covered_imgt` positions. Those regions use IMGT
boundaries independently of the selected display scheme.

`output=True` writes ANARCI text to stdout or `outfile`. With `csv=True`, an
output root is required and separate chain files are produced. See
[COMPATIBILITY.md](COMPATIBILITY.md) for legacy arguments and I/O details.

### Result coordinates

Each numbered domain is `(positions, start, end)`. Positions include gaps and
have the form `((number, insertion_code), residue)`. Numbered domain endpoints
are inclusive, as in ANARCI. Alignment-detail `query_start`/`query_end` describe
a zero-based, half-open HMM query span. These spans can differ at boundaries
when the numbering algorithm extends an alignment. Preserve the endpoint type
when slicing. Multiple domains remain in sequence order.

Batch input accepts the supported HMMER ambiguity alphabet. The legacy
single-sequence validation contract accepts the 20 standard amino acids and
rejects sequences with 10,000 or more residues. Native scans enforce their
own safe length bound before C calls. Invalid input raises a Python exception;
a valid sequence with no detected domain returns a no-hit result.

## Sequence and structure tools

These functions use the shared Rust annotation library. Computational work
releases the Python GIL. Input text is content, not a filename.

| Function | Input and result |
|---|---|
| `read_sequences(text)` | FASTA/FASTQ text; returns dictionaries with `id`, `description`, `sequence`. FASTQ checks headers and quality length. |
| `translate_dna(dna, frame=0)` | DNA string; returns protein. Frames 0, 1, 2 are forward; 3, 4, 5 use the reverse complement. Stops are `*`. |
| `annotate_protein(sequences, ...)` | `(id, protein)` pairs; returns one dictionary per numbered domain. |
| `annotate_dna(sequences, ...)` | `(id, DNA)` pairs; returns `domains` and `frame_counts`. Annotation masks translated stops to `X`. |
| `parse_structure(text, format="pdb")` | PDB or mmCIF text; returns chain dictionaries with sequence and residue C-alpha coordinates. |
| `pair_structure(text, format="pdb", ...)` | Structure text; returns inferred Fvs, unpaired domains, assemblies and contact reports. |
| `format_annotations(sequences, output_format="json", input_type="protein", ...)` | Sequence pairs; returns formatted `json`, `tsv` or `anarci` text. Input type is `protein` or `dna`. |

Annotation and pairing options are `scheme="imgt"`, `database="pan"`,
`assign_germline=False`, `allowed_species=None`, `allow=None`,
`bit_score_threshold=80`, `germline_method=None`, and `ncpu=0`.
These shared-library APIs select the bundled `pan` or `ALL` databases.
The high-level numbering API additionally accepts custom HMM files.

```python
import anarci_rs

with open("reads.fastq") as handle:
    reads = anarci_rs.read_sequences(handle.read())
result = anarci_rs.annotate_dna(
    [(row["id"], row["sequence"]) for row in reads],
    database="pan", ncpu=4, assign_germline=True,
)

with open("structure.cif") as handle:
    pairs = anarci_rs.pair_structure(handle.read(), format="cif")
```

Domain dictionaries contain `query`, `domain`, `chain_type`, `species`,
`scheme`, `start`, `end`, `evalue`, `bitscore`, `complete`, `regions`,
`numbering` and `numbered_seq`. DNA domains include `frame`. Requested
germlines add available V/J gene and V identity fields. Structured annotation
spans use zero-based, half-open alignment coordinates.

The `complete` field means that each of the seven IMGT regions is represented
at least partly. The TSV label `full` uses this same historical convention.
Inspect `regions` to identify partial terminal regions; this flag alone does
not establish that every terminal residue is present.

DNA annotation tries frames in an adaptive order and returns the first frame
with accepted domains. It does not enumerate all valid translations. With
parallel reads, a read that matches more than one frame can depend on the
adaptive worker history. Frame counts describe the frames selected by that run.

Pairing uses IMGT CDR coordinates even when another display scheme is selected.
Pairing and format labels are geometric heuristics. Missing coordinates are
excluded from contact calculations. Parsing returns sequence order and C-alpha
coordinates; original residue numbering and insertion codes are not returned.
The parser reads the first coordinate model and one selected alternate location.

## Core and research helpers

ANARCI-compatible state-vector, germline and scheme helpers are described in
[COMPATIBILITY.md](COMPATIBILITY.md). Experimental codec, augmentation, dataset
and optional model functions are described in [EXPERIMENTAL.md](EXPERIMENTAL.md).
HMMER pipeline failures raise an error rather than returning an empty result.

`anarci_rs.species_reference` fits, validates and scores supplied IMGT-numbered
repertoire cohorts. `anarci_rs.species_adaptation` exposes `profile_sequence`,
`profile_numbered`, `propose_sequence` and `propose_numbered`. These experimental
stdlib Python tools use the native engine for raw-sequence numbering; their
independent-position baseline and bounded mutation search do not predict binding
or immune tolerance. References and domain format are explicit inputs. See
[SPECIES_ADAPTATION.md](SPECIES_ADAPTATION.md) for contracts and real OAS preparation.
