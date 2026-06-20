#!/usr/bin/env python
"""Generate golden correctness fixtures from REFERENCE ANARCI (conda anarci 2024.05.21).

Run with a Python that has reference ANARCI installed (e.g. `conda install -c bioconda
anarci`):
  python scripts/gen_golden.py
Override the example-FASTA dir with $ANARCI_RS_EXAMPLES (defaults to ./examples).

Captures, per sequence, per identified domain:
  - state_vector (engine boundary) + alignment scalars  -> Phase 2 (engine) oracle
  - per-scheme numbering result OR the exact exception   -> Phase 1 (numbering) oracle
  - germline assignment genes dict                       -> Phase 1 (germline) oracle

Output: tests/fixtures/golden.json.gz   (gzipped to keep the repo light)

This is the single source of truth the Rust port must reproduce byte-for-byte.
No silent skips: any sequence that errors is recorded with its error, not dropped.
"""
import os, sys, gzip, json, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# Example FASTAs are bundled in ./examples (from ANARCI, BSD-3 — see NOTICE.md).
EX = os.environ.get("ANARCI_RS_EXAMPLES", os.path.join(ROOT, "examples"))

# All 7 supported species, least-restrictive (so every supported chain numbers).
SPECIES = ['human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca']
SCHEMES = ['imgt', 'chothia', 'kabat', 'martin', 'aho', 'wolfguy']
BIT_THRESHOLD = 80

import anarci
from anarci.anarci import (
    run_hmmer, check_for_j, number_sequence_from_alignment,
    run_germline_assignment, read_fasta,
)

sys.path.insert(0, os.path.join(ROOT, "reference_data"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_sequences():
    seqs = []
    seen = set()

    def add(name, s):
        s = s.strip().upper()
        if name in seen:
            return
        seen.add(name)
        seqs.append((name, s))

    # 1) full antibody set (997) — the core coverage
    for name, s in read_fasta(os.path.join(EX, "antibody_sequences.fasta")):
        add("ab_" + name.split()[0], s)
    # 2) 12e8 chains (heavy, light, light-dup, heavy-dup)
    for name, s in read_fasta(os.path.join(EX, "12e8.fasta")):
        add("p12e8_" + name.split('|')[0].replace(':', '_'), s)
    # 3) lysozyme negative control
    for name, s in read_fasta(os.path.join(EX, "lysozyme.fasta")):
        add("neg_" + name.split()[0], s)
    # 4) hand-picked edge cases (incl. TCR a/b, scFv multi-domain, long CDR3, too-short)
    sys.path.insert(0, HERE)
    from edge_cases import EDGE_CASES
    for name, s in EDGE_CASES:
        add(name, s)
    return seqs


def jsonable(x):
    """tuples->lists, recursively; leaves None/str/int/float as-is."""
    if isinstance(x, tuple):
        return [jsonable(v) for v in x]
    if isinstance(x, list):
        return [jsonable(v) for v in x]
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    return x


def main():
    sequences = load_sequences()
    print(f"loaded {len(sequences)} sequences", file=sys.stderr)

    # ---- run the reference HMM pipeline ONCE (scheme-independent state vectors) ----
    alignments = run_hmmer(sequences, hmm_database="ALL", hmmerpath="", ncpu=10,
                           bit_score_threshold=BIT_THRESHOLD, hmmer_species=SPECIES)
    # check_for_j mutates alignments in place; scheme arg is not materially used.
    check_for_j(sequences, alignments, "imgt")

    out_seqs = []
    n_dom = 0
    n_num = 0
    n_err = 0
    for i, (name, seq) in enumerate(sequences):
        hit_table, state_vectors, details = alignments[i]
        domains = []
        for d, (sv, det) in enumerate(zip(state_vectors, details)):
            ct = det["chain_type"]
            n_dom += 1
            dom = {
                "order": d,
                "species": det.get("species"),
                "chain_type": ct,
                "bitscore": det.get("bitscore"),
                "evalue": det.get("evalue"),
                "query_start": det.get("query_start"),
                "query_end": det.get("query_end"),
                "state_vector": jsonable(sv),
                "numbering": {},
            }
            # germline assignment (Phase 1 oracle)
            try:
                genes = run_germline_assignment(sv, seq, ct, allowed_species=SPECIES)
                dom["germlines"] = jsonable(genes)
            except Exception as e:  # record, never silently drop
                dom["germlines_error"] = "%s: %s" % (type(e).__name__, e)
            # per-scheme numbering (Phase 1 oracle)
            for scheme in SCHEMES:
                try:
                    numbering, start, end = number_sequence_from_alignment(sv, seq, scheme=scheme, chain_type=ct)
                    dom["numbering"][scheme] = {
                        "numbering": jsonable(numbering),
                        "start": start, "end": end,
                    }
                    n_num += 1
                except Exception as e:
                    dom["numbering"][scheme] = {"error": "%s: %s" % (type(e).__name__, e)}
                    n_err += 1
            domains.append(dom)
        out_seqs.append({"id": name, "seq": seq, "domains": domains,
                         "hit_table": jsonable(hit_table)})

    meta = {
        "anarci_version_str": anarci.__version__,
        "conda_version": "2024.05.21",
        "species": SPECIES,
        "schemes": SCHEMES,
        "bit_score_threshold": BIT_THRESHOLD,
        "n_sequences": len(sequences),
        "n_domains": n_dom,
        "n_numberings_ok": n_num,
        "n_numberings_err": n_err,
        "hmm_sha256": sha256(os.path.join(ROOT, "reference_data/dat/HMMs/ALL.hmm")),
        "germlines_sha256": sha256(os.path.join(ROOT, "reference_data/germlines.py")),
    }
    payload = {"meta": meta, "sequences": out_seqs}

    outdir = os.path.join(ROOT, "tests", "fixtures")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "golden.json.gz")
    with gzip.open(outpath, "wt", encoding="utf-8") as f:
        json.dump(payload, f, separators=(",", ":"), sort_keys=True)

    print(json.dumps(meta, indent=2), file=sys.stderr)
    print(f"wrote {outpath} ({os.path.getsize(outpath)} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
