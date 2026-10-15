# ANARCI compatibility

`anarci-rs` is a Rust implementation with Python compatibility entry points.
The distribution provides the `anarci` import namespace so existing numbering
code can use the Rust backend. It also provides `anarci_rs`, whose defaults and
extensions are separate from the compatibility contract.

```python
from anarci import anarci, run_anarci, number

numbered, details, hits = anarci(records, scheme="imgt")
records, numbered, details, hits = run_anarci(records, 4)
positions, chain_class = number(sequence)
```

Use a clean environment for this distribution. The original ANARCI package
uses the same `anarci` import namespace. Pip does not provide safe ownership
sharing for those files.

## High-level functions

The compatibility signatures are:

```python
anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None,
       csv=False, allow={"H", "K", "L", "A", "B", "G", "D"}, hmmerpath="",
       ncpu=None, assign_germline=False, allowed_species=["human", "mouse"],
       bit_score_threshold=80)

run_anarci(seq, ncpu=1, **kwargs)

number(sequence, scheme="imgt", database="ALL",
       allow={"H", "K", "L", "A", "B", "G", "D"},
       allowed_species=["human", "mouse"])
```

`run_anarci` accepts a FASTA filename, a single protein string, or a list/tuple
of identifier/sequence pairs. FASTA headers and gzipped FASTA input are
preserved. Output text and chain CSV files use the reference layout.
`number` validates sequence letters and the scheme before its short-sequence
gate, returns the first numbered domain, and silently returns `(False, False)`
for valid input shorter than 70 residues or with no accepted domain.

The compatibility defaults select the exact 29-profile database and identity
germline matching. Compute runs in Rust; `ncpu` selects Rust workers rather
than Python worker processes. No temporary `hmmscan` subprocess is used.

## Helper imports

| Import surface | Supported functions and data |
|---|---|
| `anarci` / `anarci.anarci` | High-level functions, `run_hmmer`, `parse_hmmer_output`, `number_sequence_from_alignment`, `number_sequences_from_alignment`, `run_germline_assignment`, `get_hmm_length`, `get_identity`, `check_for_j`, sequence validation and I/O helpers |
| `anarci.schemes` | IMGT/AHo and heavy/light Chothia, Kabat, Martin, Wolfguy functions; `smooth_insertions`, `get_imgt_cdr`, `get_cdr3_annotations`, `gap_missing`, `alphabet`, `blosum62` |
| `anarci.germlines` | `all_germlines` with the reference segment/chain/species/gene dictionary structure |
| Private legacy helpers | `_parse_hmmer_query`, `_hmm_alignment_to_states`, `_domains_are_same`, `_number_regions`, `_get_wolfguy_L1` in their original modules |

State-vector helpers accept `((position, "m" | "i" | "d"), sequence_index)`
tuples. Positions are 1 through 128; deletions use `None` as sequence index.
The native boundary validates indices and alignment metadata before numbering.
`run_hmmer` returns the reference alignment triples without starting an external
HMMER executable. The extended package exports core helpers too.

`parse_hmmer_output` reads previously saved HMMER text with Biopython. Install
the optional parser dependency when using it:

```bash
python -m pip install 'anarci-rs[legacy-parser]'
```

Normal numbering, native scans, file reading and output generation need no
Biopython dependency.

## Explicit differences

- A nonempty `hmmerpath` raises `NotImplementedError`. This argument selects an
  external executable in ANARCI; this implementation uses the bundled in-process
  engine. Remove that executable-path argument when migrating.
- `database` may select `pan`, `ALL` or a native-compatible HMM path. A custom
  numbering database must have 128-state profiles named `species_chain`, with a
  supported chain identifier and all reference positions marked `x` in RF
  annotations. Arbitrary HMMs require a different alignment mapping.
- Native failures raise errors. They are distinct from a successful no-hit scan.
  Invalid batch input is rejected before C calls, including NUL bytes and inputs
  above the safe native scan bound.
- Empty batches return empty lists. This is a supported extension to the
  reference wrapper's behavior on some empty inputs.
- Worker count 0 selects the global Rayon pool. Negative worker counts are
  rejected. Rust pools and the reference multiprocessing implementation have
  different scheduling and startup costs.
- Version metadata identifies `anarci-rs`; internal ANARCI file paths and
  executable/database installation layouts are not replicated.

These documented differences mean compatibility is a defined API contract,
not a claim that every use of the original implementation is interchangeable.
The exact profile set is the appropriate choice for reference comparisons;
pan is a separate performance option.

## Validation reference

The reference source hashes, installed distribution version, HMMER version,
profile hash, permitted test-input provenance and comparison tolerances are
recorded in [reference-version.json](reference-version.json). The local conda
artifact pin and source version can use different version labels; the hashes
identify the code actually compared. Release results and limits are in
[RELEASE_EVIDENCE.md](RELEASE_EVIDENCE.md).
