# Third-party components and attribution

The anarci-rs source code is BSD-3-Clause licensed. Bundled reference data has
additional attribution and data-license requirements, described below.

## ANARCI (Oxford Protein Informatics Group)

The numbering schemes, alignment-to-state-vector logic, germline assignment and
pipeline behavior derive from ANARCI. The reference data was pinned from conda
`anarci 2024.05.21`; it includes germline data and the ALL/FEW HMM profiles.

- Source: https://github.com/oxpig/ANARCI
- Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, Claire Marks.
- License: BSD-3-Clause; full retained text: `licenses/ANARCI-BSD-3-Clause.txt`.
- Citation: Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor
  classification. Bioinformatics (2016).

The data converted into JSON and the pan profile aggregation are adaptations.

## HMMER 3.4 and Easel

HMMER 3.4 and its bundled Easel library are fetched from the official distribution,
verified against SHA-256 and statically linked into the native extension.

- Source: http://hmmer.org
- Full copyright, contributor, license and disclaimer texts:
  `licenses/HMMER-LICENSE.txt` and `licenses/Easel-LICENSE.txt`.
- These files preserve the upstream attribution and applicable BSD terms.

## IMGT data

Copyright © 1995–2026 IMGT®, the international ImMunoGeneTics information system®.
The germline reference data underlying ANARCI ultimately derives from IMGT.

- Source and attribution: https://www.imgt.org
- Terms: https://www.imgt.org/about/termsofuse.php (inspected 2026-10-15).
- Current terms state that IMGT data are provided to public and private users under
  CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
- Full license: `licenses/CC-BY-4.0.txt`.
- Adaptations include ANARCI reference/profile processing, the germline JSON
  conversion and pan profile aggregation. Attribution does not imply endorsement.

The current terms grant a general license for data and metadata, including use by
private users. The relevant official paragraphs are retained in
`licenses/IMGT-DATA-TERMS.txt`. IMGT tool-use terms are separate from its data license.

## Rust dependencies

The native extension links Rust dependencies under their respective licenses.
`licenses/RUST-THIRD-PARTY.txt` retains their declared license identifiers and full
available license/copyright/notice texts for the locked default and optional CPU
model dependency graphs on the configured Linux and macOS architectures. Regenerate
it with `python scripts/generate_rust_notices.py` after dependency changes.
