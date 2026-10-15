# Experimental Python tools

`anarci_rs.experimental` exposes the existing Rust label codec, constrained
decoder, augmentation and dataset helpers. These functions operate on supplied
labels and scores. Their correctness tests use synthetic inputs. They do not
establish the accuracy of a trained numbering model.

```python
from anarci_rs import experimental as exp

states = [(1, "M", 0), (2, "D", None), (3, "M", 1)]
encoded = exp.encode_domain(states)
decoded = exp.decode_domain(encoded)
fragment = exp.LabeledSeq.from_domain("AC", states)
joined = exp.concat_domains([fragment, fragment])
```

## Data contracts

An experimental state is `(position, "M" | "I" | "D", sequence_index | None)`.
Positions are 1 through 128. Sequence indices start at zero; delete states have
no sequence index. Compatibility helpers use ANARCI's different state shape,
`((position, "m" | "i" | "d"), sequence_index | None)`.

Labels are `("match", position)`, `("insert", position)`, or one of
`"out"`, `"linker"`, `"constant"`, `"hinge"`, `"other"`. The vocabulary
contains 261 labels. Label counts and state indices are checked at the Python
boundary. Invalid shapes, positions and nonfinite scores raise exceptions.

## Callable coverage

| Group | Functions |
|---|---|
| Codec | `encode_domain`, `decode_domain`, `strip_boundary_deletes`, `decode_sequence` |
| Decoder | `label_to_index`, `index_to_label`, `viterbi_decode` |
| Segments | `LabeledSeq.from_domain`, `.region`, `.junk`, `.restrict_positions`, `.restrict_window`, `.domains` |
| Composition | `concat`, `concat_domains`, `gly_ser_linker` |
| Dataset | `DomainRecord`, `cdr3_of`, `overlap_identity`, `cluster_records`, `split_for_rep`, `splits_of`, `near_duplicate`, `build_dataset` |

`viterbi_decode` accepts a finite `sequence_length x 261` score matrix.
The decoder enforces the Rust transition rules. Decoding valid state vectors
does not guarantee that their biological assignments are correct.

`concat` preserves labels directly. `concat_domains` inserts a `G` residue
labeled `out` between domains, so abutting partials remain separate. Use
`restrict_positions(lo, hi)` for an inclusive IMGT window, or
`restrict_window(start, end)` for a half-open residue window.

`DomainRecord(sequence, states)` holds a source sequence and one domain's states.
Clustering uses the Rust full-sequence and CDR3 identity rules, with defaults
`FULL_THR=0.95` and `CDR3_THR=0.80`. Splits are deterministic. These clusters
reflect sequence similarity, not shared target, epitope or specificity.
`build_dataset` constructs the Rust training examples, including its generated
negative examples and requested partial windows. Inspect its returned labels
and split assignments before using them for model evaluation.

## Optional model code

The Rust workspace retains the Candle model and training implementation.
The default wheel keeps this runtime optional. A source build with the
`learned-model` feature exposes the CPU model tools in Python:

```bash
maturin develop --release --features learned-model
```

An untrained model produces scores from initialized weights. Its training
helper measures accuracy on the supplied training examples. Neither output is
evidence of held-out numbering accuracy. There is no bundled validated trained
checkpoint or production learned-engine default.

CUDA and Metal remain Rust build features. GPU model execution needs a matching
toolchain and separate verification; ordinary wheels use the HMM backend on CPU.
