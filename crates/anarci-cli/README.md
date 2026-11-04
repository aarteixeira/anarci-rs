# Rust command and shared annotation library

The `anarci-cli` crate provides the `anarci` Rust command and the library used
by the Python feature bindings. The installed Python distribution exposes the
same command through `anarci-rs`.

```bash
cargo run --release -p anarci-cli -- antibodies.fasta --engine pan
anarci-rs reads.fastq --input-type dna --format json --report-frames
anarci-rs structure.cif --input-type cif --pair --output pairs.json --format json
```

Run `anarci-rs --help` for all options. Inputs include protein or DNA FASTA/FASTQ,
PDB and mmCIF; use `-` for standard input. Explicit `--input-type protein`,
`dna`, `pdb` or `cif` overrides filename detection. `mmcif` is an alias for `cif`.
Sequence input is uncompressed; decompress files before passing their content.

`--engine pan` selects the seven-profile database. `--engine all` selects the
reference 29-profile database. Configure the numbering scheme, HMM cutoff,
chain/species restrictions, worker count and identity/evalue germline method
with the documented flags. `--ncpu 0` uses the global Rayon pool.

JSON output retains numbered positions, alignment spans and IMGT region status.
TSV string fields are quoted when they contain tabs, newlines or quotes.
`complete=true` and TSV `full` mean all seven IMGT regions are represented at
least partly; inspect region status for terminal coverage.

DNA annotation masks translated stops to `X`, tries frames adaptively, and
returns the first successful frame per read. Structure parsing reads the first
model's protein C-alpha residues. Pairing uses IMGT CDR coordinates and geometric
contact heuristics. Contact reports describe the input coordinates and require
experimental context for biological interpretation.

See the root [README](../../README.md), [Python API](../../docs/API.md),
[build guide](../../docs/BUILDING.md) and
[release procedure](../../docs/RELEASING.md).
