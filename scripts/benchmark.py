#!/usr/bin/env python
"""Benchmark anarci-rs vs reference ANARCI: ACCURACY (byte-for-byte) and SPEED.

Run in an env that has BOTH `anarci` (reference, conda 2024.05.21) and `anarci_rs`
(this project's wheel) installed, e.g. after `maturin develop` into the anarci-ref env:

  /Users/Andre.Teixeira/miniforge3/envs/anarci-ref/bin/python scripts/benchmark.py --accuracy
  /Users/Andre.Teixeira/miniforge3/envs/anarci-ref/bin/python scripts/benchmark.py --speed

The drop-in contract: anarci_rs exposes run_anarci / anarci / number with identical
signatures and return shapes; this script uses the SAME call for both modules.
"""
import argparse, sys, time, os, json

EX = "/tmp/ANARCI_src/Example_scripts_and_sequences"
SCHEMES = ["imgt", "chothia", "kabat", "martin", "aho", "wolfguy"]
SPECIES = ['human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca']
# CLI behaviour: antibody schemes ignore TCR chains (else AssertionError on A/B/G/D).
ALLOW_ALL = set(["H", "K", "L", "A", "B", "G", "D"])
ALLOW_IG = set(["H", "K", "L"])


def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def load_set(n=None):
    seqs = read_fasta(os.path.join(EX, "antibody_sequences.fasta"))
    if n is None:
        return seqs
    out = []
    i = 0
    while len(out) < n:
        name, s = seqs[i % len(seqs)]
        out.append((f"{name}_{i}", s))
        i += 1
    return out


# ---- deep comparison -------------------------------------------------------

def approx(a, b, rel=1e-6, abs_=1e-9):
    try:
        return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))
    except TypeError:
        return a == b


def cmp_numbering(a, b, path):
    # a,b: list of domains; each (numbering, start, end) or None
    if a is None or b is None:
        return [] if a == b else [f"{path}: numbered None vs not ({a is None},{b is None})"]
    if len(a) != len(b):
        return [f"{path}: ndomains {len(a)} vs {len(b)}"]
    diffs = []
    for di, (da, db) in enumerate(zip(a, b)):
        na, sa, ea = da
        nb, sb, eb = db
        if sa != sb or ea != eb:
            diffs.append(f"{path}.dom{di}: start/end ({sa},{ea}) vs ({sb},{eb})")
        if len(na) != len(nb):
            diffs.append(f"{path}.dom{di}: len {len(na)} vs {len(nb)}")
            continue
        for (pa, aaa), (pb, aab) in zip(na, nb):
            if tuple(pa) != tuple(pb) or aaa != aab:
                diffs.append(f"{path}.dom{di}: {pa}{aaa} vs {pb}{aab}")
                break
    return diffs


def cmp_details(a, b, path):
    if a is None or b is None:
        return [] if (a is None) == (b is None) else [f"{path}: details None mismatch"]
    if len(a) != len(b):
        return [f"{path}: details len {len(a)} vs {len(b)}"]
    diffs = []
    keys = ["species", "chain_type", "query_start", "query_end", "scheme"]
    for di, (da, db) in enumerate(zip(a, b)):
        for k in keys:
            if da.get(k) != db.get(k):
                diffs.append(f"{path}.dom{di}.{k}: {da.get(k)} vs {db.get(k)}")
        if not approx(da.get("bitscore", 0), db.get("bitscore", 0), rel=1e-3, abs_=0.2):
            diffs.append(f"{path}.dom{di}.bitscore: {da.get('bitscore')} vs {db.get('bitscore')}")
        ga, gb = da.get("germlines"), db.get("germlines")
        if (ga is None) != (gb is None):
            diffs.append(f"{path}.dom{di}.germlines presence")
        elif ga:
            for g in ("v_gene", "j_gene"):
                if ga.get(g, [None])[0] != gb.get(g, [None])[0]:
                    diffs.append(f"{path}.dom{di}.{g}: {ga.get(g)} vs {gb.get(g)}")
    return diffs


def run_accuracy(ref, rs, seqs, schemes):
    print(f"ACCURACY: {len(seqs)} sequences × {len(schemes)} schemes\n")
    total_bad = 0
    for scheme in schemes:
        allow = ALLOW_ALL if scheme in ("imgt", "aho") else ALLOW_IG
        kw = dict(scheme=scheme, allow=allow, assign_germline=True,
                  allowed_species=SPECIES, bit_score_threshold=80, ncpu=1)
        _, n_ref, d_ref, h_ref = ref.run_anarci(seqs, **kw)
        _, n_rs, d_rs, h_rs = rs.run_anarci(seqs, **kw)
        diffs = []
        for i in range(len(seqs)):
            diffs += cmp_numbering(n_ref[i], n_rs[i], f"seq{i}")
            diffs += cmp_details(d_ref[i], d_rs[i], f"seq{i}")
        bad = len(diffs)
        total_bad += bad
        print(f"  {scheme:8} mismatches={bad}")
        for d in diffs[:8]:
            print(f"      {d}")
    print(f"\nTOTAL mismatches: {total_bad}  ->  {'PASS' if total_bad == 0 else 'FAIL'}")
    return total_bad


def run_speed(ref, rs, sizes, ncpus):
    print("SPEED (seq/s, higher is better)\n")
    print(f"  {'N':>7} {'ncpu':>5} {'ref':>10} {'anarci_rs':>12} {'speedup':>8}")
    for n in sizes:
        seqs = load_set(n)
        for ncpu in ncpus:
            kw = dict(scheme="imgt", allow=ALLOW_ALL, assign_germline=False,
                      allowed_species=SPECIES, bit_score_threshold=80, ncpu=ncpu)
            t0 = time.perf_counter(); ref.run_anarci(seqs, **kw); tref = time.perf_counter() - t0
            t0 = time.perf_counter(); rs.run_anarci(seqs, **kw); trs = time.perf_counter() - t0
            sref, srs = n / tref, n / trs
            print(f"  {n:>7} {ncpu:>5} {sref:>10.1f} {srs:>12.1f} {srs/sref:>7.1f}x")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accuracy", action="store_true")
    ap.add_argument("--speed", action="store_true")
    ap.add_argument("--n", type=int, default=None, help="accuracy subset size")
    ap.add_argument("--schemes", default=",".join(SCHEMES))
    args = ap.parse_args()

    import anarci as ref
    import anarci_rs as rs

    if args.accuracy or not (args.accuracy or args.speed):
        seqs = load_set(args.n)
        rc = run_accuracy(ref, rs, seqs, args.schemes.split(","))
        if rc:
            sys.exit(1)
    if args.speed:
        run_speed(ref, rs, [1000, 10000], [1, os.cpu_count()])


if __name__ == "__main__":
    main()
