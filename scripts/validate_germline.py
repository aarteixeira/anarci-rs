#!/usr/bin/env python
"""Validate anarci-rs e-value germline assignment against the identity path and a
truth set.

Two comparisons are produced:

  1. e-value vs identity: how often the two methods (pan engine) disagree on the
     V-gene and J-gene call over the test set. This quantifies the *change*.

  2. accuracy vs a truth set: when `riot_na` is installed (the reference tool from
     the RIOT paper, NaturalAntibody), use its amino-acid V/J calls as TRUTH and
     report V-gene / J-gene agreement for BOTH our methods. riot_na is itself an
     e-value SW assigner, so high agreement of our e-value path with riot_na is a
     reproduction check, not an independent gold standard — see the printed caveat.

Allele handling: by default we compare at GENE level (strip the `*NN` allele
suffix), because germline DBs and allele calls drift between tools; pass
--allele to compare exact alleles.

Run after `maturin develop -m crates/anarci-py/Cargo.toml`:
    python scripts/validate_germline.py [--allele] [--limit N]
"""
import argparse
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EX = os.environ.get("ANARCI_RS_EXAMPLES", os.path.join(ROOT, "examples"))

import anarci_rs  # noqa: E402


def read_fasta(path):
    out, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def load_sequences():
    """Same set the golden fixture is built from (gen_golden.load_sequences)."""
    seqs, seen = [], set()

    def add(name, s):
        s = s.strip().upper()
        if name in seen:
            return
        seen.add(name)
        seqs.append((name, s))

    for name, s in read_fasta(os.path.join(EX, "antibody_sequences.fasta")):
        add("ab_" + name.split()[0], s)
    for name, s in read_fasta(os.path.join(EX, "12e8.fasta")):
        add("p12e8_" + name.split("|")[0].replace(":", "_"), s)
    for name, s in read_fasta(os.path.join(EX, "lysozyme.fasta")):
        add("neg_" + name.split()[0], s)
    sys.path.insert(0, HERE)
    try:
        from edge_cases import EDGE_CASES

        for name, s in EDGE_CASES:
            add(name, s)
    except Exception as e:  # noqa: BLE001
        print(f"warning: could not load edge_cases ({e})", file=sys.stderr)
    return seqs


def gene_only(call):
    return call.split("*")[0] if call else call


def first_domain_germlines(dets):
    """First numbered domain's germline dict for one sequence, or None.
    `dets` is the per-sequence list of domain dicts (None if nothing numbered)."""
    if not dets:
        return None
    return dets[0].get("germlines")


def call_anarci_rs(seqs, method, allele, allowed_species=None):
    """Map name -> (species, v_call, j_call) at gene or allele granularity."""
    _, _, details, _ = anarci_rs.run_anarci(
        seqs, scheme="imgt", database="pan", assign_germline=True, germline_method=method,
        allowed_species=allowed_species,
    )
    out = {}
    for (name, _), dets in zip(seqs, details):
        g = first_domain_germlines(dets)
        if not g:
            out[name] = (None, None, None)
            continue
        sp = g["v_gene"][0][0] if g["v_gene"][0] else None
        v = g["v_gene"][0][1] if g["v_gene"][0] else None
        j = g["j_gene"][0][1] if g["j_gene"][0] else None
        if not allele:
            v, j = gene_only(v), gene_only(j)
        out[name] = (sp, v, j)
    return out


def call_riot(seqs, allele):
    """Truth via riot_na (amino-acid). name -> (species, v_call, j_call)."""
    try:
        from riot_na import Scheme, create_riot_aa
    except Exception:  # noqa: BLE001
        return None
    riot = create_riot_aa()
    out = {}
    for name, seq in seqs:
        try:
            r = riot.run_on_sequence(name, seq, scheme=Scheme.IMGT)
            v, j = r.v_call, r.j_call
            sp = getattr(r, "locus_species", None) or None
            if not allele:
                v, j = gene_only(v), gene_only(j)
            out[name] = (sp, v, j)
        except Exception:  # noqa: BLE001
            out[name] = (None, None, None)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allele", action="store_true", help="compare exact alleles (default: gene level)")
    ap.add_argument("--limit", type=int, default=0, help="limit sequences (0 = all)")
    args = ap.parse_args()

    seqs = load_sequences()
    if args.limit:
        seqs = seqs[: args.limit]
    print(f"loaded {len(seqs)} sequences; comparing at "
          f"{'allele' if args.allele else 'gene'} level\n")

    # Tuple layout: (species, v, j). Helpers:
    V, J = 1, 2

    ident = call_anarci_rs(seqs, "identity", args.allele)
    evalue = call_anarci_rs(seqs, "evalue", args.allele)

    # ---- (1) e-value vs identity disagreement (same DB; the load-bearing result) ----
    both_v = [n for n in ident if ident[n][V] and evalue[n][V]]
    both_j = [n for n in ident if ident[n][J] and evalue[n][J]]
    v_diff = [n for n in both_v if ident[n][V] != evalue[n][V]]
    j_diff = [n for n in both_j if ident[n][J] != evalue[n][J]]
    print("=== (1) e-value vs identity (anarci-rs pan engine, all 8 species) ===")
    print(f"  V calls present in both: {len(both_v)};  disagree: {len(v_diff)} "
          f"({100*len(v_diff)/max(1,len(both_v)):.1f}%)")
    print(f"  J calls present in both: {len(both_j)};  disagree: {len(j_diff)} "
          f"({100*len(j_diff)/max(1,len(both_j)):.1f}%)")
    print("  sample V disagreements (identity -> evalue):")
    for n in v_diff[:10]:
        print(f"    {n}: {ident[n][V]} -> {evalue[n][V]}")
    print()

    # ---- (2) accuracy vs riot_na truth ----
    truth = call_riot(seqs, args.allele)
    if truth is None:
        print("=== (2) truth-set accuracy: SKIPPED (riot_na not installed) ===")
        return

    print("=== (2) V/J accuracy vs riot_na truth (RIOT paper reference) ===")
    print("  CAVEATS:")
    print("   - riot_na uses the OGRDB germline DB; we use ANARCI's IMGT DB. Gene")
    print("     naming/content differ, so absolute agreement is capped by the DB gap,")
    print("     not the method.")
    print("   - riot_na only supports human/mouse/alpaca; our pan search spans 8")
    print("     species, so unrestricted search can pick rat/rhesus/etc. genes riot_na")
    print("     cannot. The FAIR comparison restricts our search to riot_na's species.")
    print("   - riot_na is itself an e-value SW assigner, so this is a reproduction")
    print("     check, not an independent curated gold standard.\n")

    # 2a) unrestricted (all 8 species) — shows the species confound.
    for seg, label in ((V, "V"), (J, "J")):
        scored = [n for n in seqs_names(seqs) if truth[n][seg]]
        i_ok = sum(1 for n in scored if ident[n][seg] == truth[n][seg])
        e_ok = sum(1 for n in scored if evalue[n][seg] == truth[n][seg])
        print(f"  [all 8 species]  {label}-gene over {len(scored)} seqs: "
              f"identity {100*i_ok/max(1,len(scored)):.1f}%  e-value {100*e_ok/max(1,len(scored)):.1f}%")
    print()

    # 2b) HEADLINE: human-only. For mouse, OGRDB uses opaque IDs (IGHV-2UQD,...)
    # that cannot match IMGT names at all, so only human gives an interpretable
    # accuracy. Restrict OUR search to human and score against riot_na's human calls.
    human_truth = {n: truth[n] for n in truth if truth[n][0] == "human" and truth[n][V]}
    identH = call_anarci_rs(seqs, "identity", args.allele, allowed_species=["human"])
    evalueH = call_anarci_rs(seqs, "evalue", args.allele, allowed_species=["human"])
    print(f"  *** HEADLINE: HUMAN ONLY (IMGT/OGRDB names overlap) — {len(human_truth)} seqs ***")
    for seg, label in ((V, "V"), (J, "J")):
        scored = [n for n in human_truth if human_truth[n][seg]]
        i_ok = sum(1 for n in scored if identH[n][seg] == human_truth[n][seg])
        e_ok = sum(1 for n in scored if evalueH[n][seg] == human_truth[n][seg])
        print(f"      {label}-gene over {len(scored)} seqs: "
              f"identity {100*i_ok/max(1,len(scored)):.1f}%  ->  "
              f"e-value {100*e_ok/max(1,len(scored)):.1f}%")
    wins = [n for n in human_truth
            if evalueH[n][V] == human_truth[n][V] and identH[n][V] != human_truth[n][V]]
    losses = [n for n in human_truth
              if identH[n][V] == human_truth[n][V] and evalueH[n][V] != human_truth[n][V]]
    print(f"      V-gene: e-value RIGHT where identity WRONG: {len(wins)};  "
          f"identity RIGHT where e-value WRONG: {len(losses)}")
    for n in wins[:12]:
        print(f"        WIN  {n}: identity={identH[n][V]} -> evalue={evalueH[n][V]} (truth {human_truth[n][V]})")
    for n in losses[:8]:
        print(f"        LOSS {n}: identity={identH[n][V]} -> evalue={evalueH[n][V]} (truth {human_truth[n][V]})")


def seqs_names(seqs):
    return [n for n, _ in seqs]


if __name__ == "__main__":
    main()
