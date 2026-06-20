#!/usr/bin/env python
"""Build the pan-species HMM database FEW.hmm: ONE combined HMM per chain type
(H,K,L,A,B,G,D), pooling V+J IMGT-aligned germlines across ALL species, with the
same hmmbuild config ANARCI uses (--hand, RF reference, 128 match columns).

This is the default engine of anarci-rs: 7 profiles instead of ALL.hmm's 29, giving
a ~4x (up to ~12x with a chain-type prefilter) speedup. Species/gene are assigned
separately by germline identity (the HMM here only does chain-type + IMGT alignment).

  python scripts/build_pan_hmm.py
Requires hmmbuild + hmmpress on PATH (or set $HMMBUILD/$HMMPRESS).
"""
import os, sys, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "reference_data"))
from germlines import all_germlines  # noqa

HMMBUILD = os.environ.get("HMMBUILD", "hmmbuild")
HMMPRESS = os.environ.get("HMMPRESS", "hmmpress")
OUT = os.path.join(ROOT, "reference_data", "dat", "HMMs")
CHAINS = ["H", "K", "L", "A", "B", "G", "D"]
ALN_LEN = 128


def gather_rows(chain):
    rows = []
    for seg in ("V", "J"):
        if chain not in all_germlines[seg]:
            continue
        for sp in all_germlines[seg][chain]:
            for gene, seq in all_germlines[seg][chain][sp].items():
                assert len(seq) == ALN_LEN, (seg, chain, sp, gene, len(seq))
                rows.append((f"{seg}_{sp}_{gene}".replace(" ", "_"), seq))
    return rows


def write_stockholm(rows, path, name):
    with open(path, "w") as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {name}\n")
        w = max(len(n) for n, _ in rows) + 2
        for n, s in rows:
            f.write(f"{n.ljust(w)}{s}\n")
        f.write(f"{'#=GC RF'.ljust(w)}{'x' * ALN_LEN}\n")
        f.write("//\n")


def main():
    combined = os.path.join(OUT, "FEW.hmm")
    parts = []
    for ch in CHAINS:
        rows = gather_rows(ch)
        sto = os.path.join(OUT, f"pan_{ch}.sto")
        write_stockholm(rows, sto, f"pan_{ch}")
        hmm = os.path.join(OUT, f"pan_{ch}.hmm")
        r = subprocess.run([HMMBUILD, "--hand", "--amino", "-n", f"pan_{ch}", hmm, sto],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"hmmbuild FAILED for {ch}: {r.stderr[-500:]}")
        parts.append(hmm)
        # canonical J residue length for this chain (what get_hmm_length must return for "pan")
        jlen = None
        if ch in all_germlines["J"]:
            sp0 = next(iter(all_germlines["J"][ch]))
            jlen = len(next(iter(all_germlines["J"][ch][sp0].values())).rstrip("-"))
        print(f"pan_{ch}: rows={len(rows)} J_residue_len={jlen}")
    with open(combined, "w") as out:
        for p in parts:
            with open(p) as f:
                out.write(f.read())
        os.remove(p)  # tidy per-chain intermediates
    for ch in CHAINS:
        sto = os.path.join(OUT, f"pan_{ch}.sto")
        if os.path.exists(sto):
            os.remove(sto)
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        fp = combined + ext
        if os.path.exists(fp):
            os.remove(fp)
    r = subprocess.run([HMMPRESS, combined], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"hmmpress FAILED: {r.stderr}")
    print("built + pressed", combined, os.path.getsize(combined), "bytes")


if __name__ == "__main__":
    main()
