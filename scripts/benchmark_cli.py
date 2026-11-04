#!/usr/bin/env python3
"""Benchmark the `anarci` CLI (pan + exact engines) and ANARCII against ANARCI on the
golden set (996 seqs / 1013 domains). Golden *is* ANARCI's numbering.

Runs the available tools in the current env (all env-guarded):
  - the CLI (standalone binary) — always, both --engine pan and --engine all
  - reference ANARCI       — if `import anarci`   works (e.g. anarci-ref env)
  - ANARCII                — if `import anarcii`  works (e.g. anarcii env)

Build first:  cargo build --release -p anarci-cli
Run e.g.:     conda run -n anarcii   python scripts/benchmark_cli.py
              conda run -n anarci-ref python scripts/benchmark_cli.py
"""
import gzip, json, os, subprocess, sys, time, re
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI = os.path.join(REPO, "target", "release", "anarci")
FASTA = os.path.join(REPO, "target", "golden_bench.fasta")

g = json.load(gzip.open(os.path.join(REPO, "tests/fixtures/golden.json.gz")))
ids = [s["id"] for s in g["sequences"]]
seqs = [s["seq"] for s in g["sequences"]]
with open(FASTA, "w") as f:
    for i, s in zip(ids, seqs):
        f.write(f">{i}\n{s}\n")


def nonrgap(items):
    return {f"{n}{i}".strip(): aa for (n, i), aa in items if aa != "-"}


# golden reference, all domains and domain-0-only (for apples-to-apples with 1-domain tools)
golden = {s["id"]: [nonrgap(d["numbering"]["imgt"]["numbering"]) for d in s["domains"]] for s in g["sequences"]}
golden0 = {k: v[:1] for k, v in golden.items()}


def compare(ref, got, label, seq_s=None, limit=None):
    domok = ref_d = got_d = exact = rtot = rmatch = 0
    mispos = Counter()
    for sid in ids:
        e = ref.get(sid, [])
        o = got.get(sid, [])
        if limit:
            e, o = e[:limit], o[:limit]
        ref_d += len(e); got_d += len(o)
        domok += len(e) == len(o)
        for i in range(min(len(e), len(o))):
            rtot += len(e[i])
            rmatch += sum(1 for k, v in e[i].items() if o[i].get(k) == v)
            if e[i] == o[i]:
                exact += 1
            else:
                for k, v in e[i].items():
                    if o[i].get(k) != v:
                        mispos[int(re.sub(r"[A-Z]", "", k))] += 1
    n = min(ref_d, got_d)
    sp = f"{seq_s:6.0f} seq/s | " if seq_s else ""
    print(f"  {label:34} {sp}domains {got_d:4} | dom-exact {exact}/{n} ({100*exact/max(1,n):6.2f}%) "
          f"| per-res {100*rmatch/max(1,rtot):7.3f}%")
    if mispos and 100 * exact / max(1, n) < 99.5:
        print(f"       └ disagreements by IMGT pos: {mispos.most_common(6)}")


def run_cli(engine):
    t0 = time.time()
    out = subprocess.run([CLI, FASTA, "--engine", engine, "--scheme", "imgt", "--format", "json"],
                         capture_output=True, text=True)
    dt = time.time() - t0
    if out.returncode != 0:
        sys.exit(f"CLI ({engine}) failed: " + out.stderr[-400:])
    got = {}
    for r in sorted(json.loads(out.stdout), key=lambda r: (r["query"], r["domain"])):
        got.setdefault(r["query"], []).append({p["n"]: p["aa"] for p in r["numbering"]})
    return got, len(ids) / dt


print(f"golden: {len(ids)} seqs / {sum(len(v) for v in golden.values())} domains\n")
print("vs golden ANARCI numbering (IMGT):")
for eng in ("pan", "all"):
    got, sps = run_cli(eng)
    compare(golden, got, f"CLI --engine {eng} (all domains)", sps)
cli_all, _ = run_cli("all")

# --- ANARCII (single domain/seq by default) ---
try:
    from anarcii import Anarcii
    m = Anarcii(seq_type="antibody", mode="accuracy", batch_size=128, cpu=False, ncpu=1, verbose=False)
    m.number(seqs[:8])  # warm up (exclude model load)
    t0 = time.time(); out = m.number(seqs); dt = time.time() - t0
    aii = {}
    for sid, rec in zip(ids, out.values()):
        num = rec.get("numbering")
        aii[sid] = [nonrgap([(tuple(p), a) for p, a in num])] if num else []
    print("\nvs golden, DOMAIN-0 basis (apples-to-apples with ANARCII's 1-domain default):")
    compare(golden0, cli_all, "CLI --engine all", limit=1)
    compare(golden0, aii, f"ANARCII 2.0.8 (accuracy, {m.device})", len(ids) / dt, limit=1)
except ImportError:
    pass

# --- reference ANARCI (live) ---
try:
    from anarci import anarci as ref_anarci
    t0 = time.time(); numbering, _d, _ = ref_anarci(list(zip(ids, seqs)), scheme="imgt", output=False); dt = time.time() - t0
    live = {ids[i]: [nonrgap([(list(p), a) for p, a in num]) for num, _s, _e in (numbering[i] or [])] for i in range(len(ids))}
    print("\nvs live reference ANARCI:")
    compare(golden, live, "golden-vs-live (version drift)")
    compare(live, cli_all, "CLI --engine all vs live", len(ids) / dt)
except ImportError:
    pass
