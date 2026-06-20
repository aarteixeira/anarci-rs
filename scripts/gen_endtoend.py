#!/usr/bin/env python
"""End-to-end oracle for the orchestration layer.

Produces:
  tests/fixtures/replay_hsps.json.gz   : seq_string -> [raw hsp dicts] for EVERY
      run_hmmer call (main scan + check_for_j rescans), so a Rust ReplayEngine can
      drive the whole pipeline without HMMER.
  tests/fixtures/endtoend_imgt.json.gz : reference anarci() output (numbered,
      details, hit_table) per sequence for scheme=imgt, assign_germline=True.

  python scripts/gen_endtoend.py
"""
import os, sys, gzip, json, importlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from gen_golden import load_sequences, jsonable, SPECIES, BIT_THRESHOLD

A = importlib.import_module("anarci.anarci")

# Map name -> sequence for the CURRENT run_hmmer call (set by the run_hmmer patch),
# so the _parse_hmmer_query patch can key captured hsps by the real sequence string.
_current = {}
replay = {}

_orig_run_hmmer = A.run_hmmer
_orig_parse = A._parse_hmmer_query


def run_hmmer_patched(sequence_list, **kw):
    global _current
    _current = {name: seq for name, seq in sequence_list}
    return _orig_run_hmmer(sequence_list, **kw)


def parse_patched(query, bit_score_threshold=80, hmmer_species=None):
    seq = _current.get(query.id)
    if seq is not None:
        hsps = []
        for hsp in query.hsps:
            hsps.append({
                "hit_id": hsp.hit_id, "hit_description": hsp.hit_description,
                "evalue": hsp.evalue, "bitscore": hsp.bitscore, "bias": hsp.bias,
                "query_start": hsp.query_start, "query_end": hsp.query_end,
                "hit_start": hsp.hit_start, "hit_end": hsp.hit_end,
                "rf": hsp.aln_annotation["RF"], "pp": hsp.aln_annotation["PP"],
            })
        replay[seq] = hsps  # key by exact sequence string (full or sub)
    return _orig_parse(query, bit_score_threshold, hmmer_species)


A.run_hmmer = run_hmmer_patched
A._parse_hmmer_query = parse_patched

sequences = load_sequences()
print(f"loaded {len(sequences)} sequences", file=sys.stderr)

numbered, details, hit_tables = A.anarci(
    sequences, scheme="imgt", output=False, assign_germline=True,
    allowed_species=SPECIES, bit_score_threshold=BIT_THRESHOLD,
)

# Serialise reference output per sequence.
out = []
for i, (name, seq) in enumerate(sequences):
    nb = numbered[i]
    dt = details[i]
    ht = hit_tables[i]
    rec = {"id": name, "seq": seq}
    if nb is None:
        rec["numbered"] = None
        rec["details"] = None
    else:
        rec["numbered"] = [[jsonable(num), start, end] for (num, start, end) in nb]
        rec["details"] = [
            {k: d.get(k) for k in ("species", "chain_type", "evalue", "bitscore", "bias",
                                   "query_start", "query_end", "scheme", "query_name")}
            | {"germlines": jsonable(d.get("germlines"))}
            for d in dt
        ]
    rec["hit_table"] = jsonable(ht)  # includes header row
    out.append(rec)

outdir = os.path.join(ROOT, "tests", "fixtures")
os.makedirs(outdir, exist_ok=True)

with gzip.open(os.path.join(outdir, "endtoend_imgt.json.gz"), "wt", encoding="utf-8") as f:
    json.dump({"meta": {"scheme": "imgt", "species": SPECIES, "assign_germline": True,
                        "bit_score_threshold": BIT_THRESHOLD},
               "sequences": out}, f, separators=(",", ":"))

with gzip.open(os.path.join(outdir, "replay_hsps.json.gz"), "wt", encoding="utf-8") as f:
    json.dump(replay, f, separators=(",", ":"))

n_dom = sum(len(r["numbered"]) for r in out if r["numbered"])
n_none = sum(1 for r in out if r["numbered"] is None)
print(f"end-to-end: {len(out)} seqs, {n_dom} numbered domains, {n_none} None; "
      f"replay entries={len(replay)}", file=sys.stderr)
