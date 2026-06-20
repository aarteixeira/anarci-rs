#!/usr/bin/env python
"""Capture raw HMMER HSPs (RF/PP/coords/scores) + the resulting state vectors from
reference ANARCI, by instrumenting `_parse_hmmer_query`. This is the Gate-A oracle:
the Rust port of parse_hmmer_query + hmm_alignment_to_states must reproduce the
state_vectors and details from the SAME hsps. (State vectors here are BEFORE
check_for_j, which is validated later with the engine.)

  python scripts/capture_hsps.py
"""
import os, sys, gzip, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from gen_golden import load_sequences, jsonable, SPECIES, BIT_THRESHOLD  # reuse loader

import importlib
A = importlib.import_module("anarci.anarci")  # the module (package attr is shadowed by the fn)

captured = {}
_orig = A._parse_hmmer_query


def patched(query, bit_score_threshold=80, hmmer_species=None):
    hsps = []
    for hsp in query.hsps:
        hsps.append({
            "hit_id": hsp.hit_id,
            "hit_description": hsp.hit_description,
            "evalue": hsp.evalue,
            "bitscore": hsp.bitscore,
            "bias": hsp.bias,
            "query_start": hsp.query_start,
            "query_end": hsp.query_end,
            "hit_start": hsp.hit_start,
            "hit_end": hsp.hit_end,
            "rf": hsp.aln_annotation["RF"],
            "pp": hsp.aln_annotation["PP"],
        })
    result = _orig(query, bit_score_threshold, hmmer_species)
    hit_table, state_vectors, details = result
    captured[query.id] = {
        "seq_len": query.seq_len,
        "hsps": hsps,
        "state_vectors": jsonable(state_vectors),
        "details": [
            {k: d.get(k) for k in ("species", "chain_type", "bitscore", "evalue",
                                   "bias", "query_start", "query_end")}
            for d in details
        ],
    }
    return result


A._parse_hmmer_query = patched

sequences = load_sequences()
print(f"loaded {len(sequences)} sequences", file=sys.stderr)

alignments = A.run_hmmer(sequences, hmm_database="ALL", hmmerpath="", ncpu=10,
                         bit_score_threshold=BIT_THRESHOLD, hmmer_species=SPECIES)

# Map captured-by-query-id (which is the fasta index "0","1",...) back to our ids.
out = []
n_hsps = 0
for i, (name, seq) in enumerate(sequences):
    rec = captured.get(name)
    if rec is None:
        out.append({"id": name, "seq": seq, "seq_len": len(seq), "hsps": [],
                    "state_vectors": [], "details": []})
        continue
    n_hsps += len(rec["hsps"])
    out.append({"id": name, "seq": seq, **rec})

payload = {
    "meta": {"species": SPECIES, "bit_score_threshold": BIT_THRESHOLD,
             "n_sequences": len(sequences), "n_hsps": n_hsps},
    "sequences": out,
}
outdir = os.path.join(ROOT, "tests", "fixtures")
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, "hsps.json.gz")
with gzip.open(outpath, "wt", encoding="utf-8") as f:
    json.dump(payload, f, separators=(",", ":"))
print(f"captured {n_hsps} hsps across {len(out)} sequences -> {outpath} "
      f"({os.path.getsize(outpath)} bytes)", file=sys.stderr)
