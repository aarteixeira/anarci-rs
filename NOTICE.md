# Third-party components and attribution

anarci-rs is a reimplementation of ANARCI and bundles HMMER + Easel. All are
BSD-3-Clause, which requires retaining the notices below.

## ANARCI (Oxford Protein Informatics Group)
The numbering schemes, alignment-to-state-vector logic, germline assignment, and
the pipeline behaviour reproduced here derive from ANARCI. The reference data in
`reference_data/` (`germlines.py`, `dat/HMMs/ALL.hmm` and its pressed files) is the
data produced by the ANARCI build pipeline, pinned from the conda package
`anarci 2024.05.21` for exact reproducibility.
- Upstream: https://github.com/oxpig/ANARCI
- License: BSD-3-Clause. Dunbar J, Deane CM. "ANARCI: antigen receptor numbering
  and receptor classification." Bioinformatics (2016).

## HMMER 3.4 and Easel (HHMI / Sean R. Eddy, Rivas Lab)
Fetched at build time (`hmmer-3.4.tar.gz`, SHA-256-pinned and verified by
`crates/hmmer-sys/build.rs`) and statically linked into the Python extension.
The tarball bundles Easel and its own `LICENSE`.
- Upstream: http://hmmer.org
- License: BSD-3-Clause.

## IMGT
The antibody/TCR germline reference sequences underlying ANARCI's HMMs and
`germlines.py` ultimately derive from IMGT®, the international ImMunoGeneTics
information system (http://www.imgt.org). Users redistributing or relying on this
data should review IMGT's terms of use and cite IMGT accordingly.
