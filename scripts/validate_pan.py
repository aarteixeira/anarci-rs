#!/usr/bin/env python
"""Validate the pan-species engine of anarci-rs.

Gates (run with an env that has reference `anarci` + the anarci_rs wheel):
  python scripts/validate_pan.py

  1. PAN ENGINE CORRECT: anarci_rs(database='pan') numbering == reference ANARCI run
     against FEW.hmm (with get_hmm_length patched to the per-chain J length, which is
     what anarci-rs does). Proves the Rust pan engine reproduces HMMER-on-FEW.
  2. EXACT PRESERVED: anarci_rs(database='ALL') numbering == stock reference ANARCI.
  3. ACCURACY: anarci_rs(pan) vs stock ANARCI — % identical + IMGT-anchor correctness.
  4. SPEED: pan vs ALL throughput (anarci_rs).
"""
import os, sys, time, importlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from gen_golden import load_sequences  # noqa

import anarci  # reference
import anarci_rs
A = importlib.import_module("anarci.anarci")
from anarci.germlines import all_germlines as G

SP = ['human', 'mouse']  # default species
HMMDIR = os.path.join(ROOT, "reference_data", "dat", "HMMs")


def numbering_only(numbered):
    """Reduce anarci() 'numbered' to comparable form (drop None vs domains, keep tuples)."""
    if numbered is None:
        return None
    return [[[(tuple(p), aa) for (p, aa) in num], start, end] for (num, start, end) in numbered]


def cmp_numbering(a, b):
    """Count sequences whose numbering differs."""
    bad = 0
    examples = []
    for i, (x, y) in enumerate(zip(a, b)):
        nx, ny = numbering_only(x), numbering_only(y)
        if nx != ny:
            bad += 1
            if len(examples) < 5:
                examples.append(i)
    return bad, examples


def run_ref(seqs, hmm_path=None, patch_ghl=False):
    if patch_ghl:
        def ghl(species, ctype):
            try:
                return len(list(G['J'][ctype][species].values())[0].rstrip('-'))
            except KeyError:
                try:
                    return max(len(list(G['J'][ctype][sp].values())[0].rstrip('-'))
                               for sp in G['J'][ctype])
                except (KeyError, ValueError):
                    return 128
        A.get_hmm_length = ghl
    old = A.HMM_path
    if hmm_path:
        A.HMM_path = hmm_path
    try:
        numbered, details, _ = anarci.anarci(seqs, scheme="imgt", allowed_species=SP)
    finally:
        A.HMM_path = old
    return numbered, details


def setup_few_as_all():
    """run_hmmer looks for '<database>.hmm' (database='ALL') in HMM_path; make a dir
    where ALL.hmm IS FEW.hmm + pressed files, so reference ANARCI scans the pan DB."""
    import tempfile, shutil
    d = tempfile.mkdtemp(prefix="anarci_few_")
    shutil.copy(os.path.join(HMMDIR, "FEW.hmm"), os.path.join(d, "ALL.hmm"))
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        shutil.copy(os.path.join(HMMDIR, "FEW.hmm" + ext), os.path.join(d, "ALL.hmm" + ext))
    return d


ANCHORS = {"H": [(23, 'C'), (41, 'W'), (104, 'C')],
           "K": [(23, 'C'), (41, 'W'), (104, 'C')],
           "L": [(23, 'C'), (41, 'W'), (104, 'C')]}


def anchor_ok(numbered, details):
    """Count H/K/L domains whose conserved IMGT anchors C23/W41/C104 carry the right residue."""
    ok = tot = 0
    for nb, dt in zip(numbered, details):
        if not nb:
            continue
        for (num, _, _), d in zip(nb, dt):
            ct = d["chain_type"]
            if ct not in ANCHORS:
                continue
            pos = {p[0]: aa for (p, aa) in num}
            tot += 1
            if all(pos.get(n) == aa for n, aa in ANCHORS[ct]):
                ok += 1
    return ok, tot


def main():
    seqs = load_sequences()
    print(f"loaded {len(seqs)} sequences\n")

    # references
    ref_all, det_all = run_ref(seqs)                       # stock ANARCI
    few_dir = setup_few_as_all()
    ref_few, _ = run_ref(seqs, hmm_path=few_dir, patch_ghl=True)  # ANARCI on FEW.hmm

    # anarci_rs
    _, rs_pan, det_pan, _ = anarci_rs.run_anarci(seqs, scheme="imgt", allowed_species=SP)            # default pan
    _, rs_all, _, _ = anarci_rs.run_anarci(seqs, scheme="imgt", database="ALL", allowed_species=SP)  # exact

    b1, ex1 = cmp_numbering(rs_pan, ref_few)
    b2, ex2 = cmp_numbering(rs_all, ref_all)
    b3, ex3 = cmp_numbering(rs_pan, ref_all)

    print("=== GATE 1: pan engine == reference-on-FEW (numbering) ===")
    print(f"   mismatches: {b1}/{len(seqs)}  {'PASS' if b1==0 else 'FAIL '+str(ex1)}")
    print("=== GATE 2: exact mode == stock ANARCI (numbering) ===")
    print(f"   mismatches: {b2}/{len(seqs)}  {'PASS' if b2==0 else 'FAIL '+str(ex2)}")
    print("=== ACCURACY: pan vs stock ANARCI(human,mouse) ===")
    print(f"   identical: {len(seqs)-b3}/{len(seqs)} = {100*(len(seqs)-b3)/len(seqs):.2f}%")
    ok_pan, tot = anchor_ok(rs_pan, det_pan)
    ok_all, _ = anchor_ok(ref_all, det_all)
    print(f"   IMGT anchor correctness (C23/W41/C104): pan {ok_pan}/{tot}={100*ok_pan/tot:.2f}%  "
          f"stock {ok_all}/{tot}={100*ok_all/tot:.2f}%")

    # speed
    print("\n=== SPEED (anarci_rs, ncpu=1) ===")
    for db in ("pan", "ALL"):
        t0 = time.perf_counter()
        anarci_rs.run_anarci(seqs, scheme="imgt", database=db, ncpu=1)
        dt = time.perf_counter() - t0
        print(f"   {db:4}: {len(seqs)/dt:7.1f} seq/s")


if __name__ == "__main__":
    main()
