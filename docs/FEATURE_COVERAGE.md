# Python feature coverage

This table maps existing Rust functionality to the Python interfaces and
release checks. A supported input contract has tests for accepted and rejected
inputs. Test results and platform evidence are recorded separately in
[RELEASE_EVIDENCE.md](RELEASE_EVIDENCE.md).

| Rust functionality | Python access | Verification |
|---|---|---|
| Exact/pan profiles, native HMMER, parallel deduplicated batches | `anarci_rs.anarci`, `run_anarci`; compatibility equivalents | Safe reference comparisons, worker and engine-switch tests |
| Six numbering schemes and aliases | `scheme` and `schemes` helpers | Reference oracle for all schemes, synthetic chain-class tests |
| Species/chain filters | `allowed_species`, `allow` | Matching, excluded and empty filters |
| Identity/evalue germline assignment | `assign_germline`, `germline_method`; core helpers | Reference identity calls and both native methods |
| Partial/multiple domains and region status | Batch APIs; `number` gates; `annotate_regions` | Truncated and composed permitted inputs |
| State-vector numbering and HSP conversion | Core compatibility helpers | Synthetic state/HSP cases and malformed-input checks |
| FASTA/gzip, FASTQ, output files | Legacy I/O and `read_sequences` | Header preservation, malformed FASTQ, stdout/file/CSV |
| Six-frame translation and adaptive DNA annotation | `translate_dna`, `annotate_dna` | Synthetic codons, reverse frames and protein/DNA agreement |
| PDB/mmCIF chain and coordinate extraction | `parse_structure` | Synthetic coordinates, quoted/wrapped CIF, missing fields |
| Pairing, assemblies and contact reports | `pair_structure` | Synthetic chain geometry and IMGT contact invariance |
| JSON/TSV/ANARCI domain output | `format_annotations`; CLI | Output schema, quoted TSV identifiers and installed command |
| Command-line configuration and file routing | `anarci-rs`; shared Rust library | Help/stdin, explicit input formats, output paths, failure status |
| Learned labels and constrained decoder | `experimental` codec/decode functions | Vocabulary round trip and malformed/nonfinite score inputs |
| Augmentation and deterministic dataset tools | `experimental.LabeledSeq`, `DomainRecord` and helpers | Label preservation, separators, clustering/splits and windows |
| Optional CPU model construction/forward/training | `experimental.Model`, `tokenize`, `overfit`, source feature `learned-model` | Tiny synthetic shape/training tests in optional build |
| Optional CUDA/Metal model execution | Existing Rust Cargo features | Separate GPU/toolchain validation required |

The optional model interface operates on initialized or caller-trained models.
The release does not include a validated trained checkpoint. Geometry outputs
are contact-based inferences; sequence clusters express similarity rather
than biological equivalence.
