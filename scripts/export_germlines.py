#!/usr/bin/env python
"""Serialize the pinned germlines.py `all_germlines` dict to JSON for Rust.

  python scripts/export_germlines.py

CRITICAL: gene insertion order is PRESERVED. `run_germline_assignment` picks the
germline via `max(...)`, which breaks ties by first-in-iteration-order, and many
J alleles are byte-identical -> ties are common. Genes are therefore emitted as an
ORDERED list of [gene, seq] pairs (not a sorted object). Species/segment order is
irrelevant (looked up by key), so those stay objects.

Structure: {seg: {chain_type: {species: [[gene, aligned_seq], ...]}}}
Output: crates/anarci-core/data/germlines.json
"""
import os, sys, json, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "reference_data"))
from germlines import all_germlines

out = {}
nseq = 0
for seg in all_germlines:                       # 'J','V'
    out[seg] = {}
    for ct in all_germlines[seg]:               # H,K,L,A,B,G,D
        out[seg][ct] = {}
        for sp in all_germlines[seg][ct]:        # species
            # PRESERVE insertion order of genes (dict order, Python 3.7+)
            genes = [[g, s] for g, s in all_germlines[seg][ct][sp].items()]
            out[seg][ct][sp] = genes
            nseq += len(genes)

# all_species = list(all_germlines['V']['H'].keys()) — order matters for the
# allowed_species=None path of run_germline_assignment.
all_species = list(all_germlines["V"]["H"].keys())
wrapper = {"all_species": all_species, "germlines": out}

outdir = os.path.join(ROOT, "crates", "anarci-core", "data")
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, "germlines.json")
# sort_keys=False everywhere: gene ORDER must survive.
blob = json.dumps(wrapper, separators=(",", ":")).encode("utf-8")
with open(outpath, "wb") as f:
    f.write(blob)

print("segments:", list(out.keys()))
print("total germline sequences:", nseq)
print("bytes:", len(blob))
print("sha256:", hashlib.sha256(blob).hexdigest())
print("wrote:", outpath)
