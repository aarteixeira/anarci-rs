#!/usr/bin/env python3
"""C0 gate: prove the FIXED learned-engine label codec is a lossless round-trip
of real ANARCI state vectors.

Label set (per residue): M{p} (p in 1..128), INS{p} (position-tagged), OUT.
Deletions carry no residue -> reconstructed on decode.

decode rule (the fix the red-team forced):
  * INS is position-tagged with its actual id (the *next* match column id in ANARCI),
    NOT "inherit preceding match". Storing the id makes the CDR3 111/112 split and the
    check_for_j id-116 long-CDR3 block representable.
  * A match position p in the covered span [min_m, max_m] becomes a D-state ONLY if it
    is occupied by neither a match NOR an insert (insert-occupied columns are not deleted
    -> this is what kills the 4 spurious d116 inventions).
  * Reconstruction walks id lo..hi, emitting inserts-with-id==p (stream order) then the
    match/delete column for p. This reproduces the original si ordering.
"""
import gzip, json, sys
from collections import defaultdict

GOLDEN = "tests/fixtures/golden.json.gz"

def to_sv(raw):
    # raw entry: [[id, "m"|"i"|"d"], si_or_null]
    return [((e[0][0], e[0][1]), e[1]) for e in raw]

def encode(sv):
    """state vector -> list of per-residue labels (kind, id, si). Drops d-states."""
    labels = []
    for (pid, typ), si in sv:
        if typ == "m":
            labels.append(("M", pid, si))
        elif typ == "i":
            labels.append(("INS", pid, si))
        # 'd': no residue, dropped (reconstructed on decode)
    return labels

def decode(labels):
    """labels -> full state vector (with reconstructed d-states)."""
    m_ids = {pid for k, pid, si in labels if k == "M"}
    if not m_ids:
        return []
    lo, hi = min(m_ids), max(m_ids)
    ins_by_id = defaultdict(list)          # id -> [si, ...] in stream order
    m_si = {}                               # id -> si
    ins_below, ins_above = [], []           # trailing/leading inserts outside [lo,hi]
    for k, pid, si in labels:
        if k == "M":
            m_si[pid] = si
        else:  # INS
            if pid < lo:
                ins_below.append((pid, si))
            elif pid > hi:
                ins_above.append((pid, si))
            else:
                ins_by_id[pid].append(si)
    sv = []
    for pid, si in ins_below:              # inserts left of the first match (synthetic partials)
        sv.append(((pid, "i"), si))
    for pid in range(lo, hi + 1):
        for si in ins_by_id.get(pid, []):  # inserts belong *before* their (next-match) column
            sv.append(((pid, "i"), si))
        if pid in m_si:
            sv.append(((pid, "m"), m_si[pid]))
        elif pid not in ins_by_id:         # gap occupied by neither match nor insert -> delete
            sv.append(((pid, "d"), None))
        # else: insert-occupied column, no residue-less delete (kills the spurious d116)
    for pid, si in ins_above:
        sv.append(((pid, "i"), si))
    return sv

def main():
    d = json.load(gzip.open(GOLDEN))
    domains = [(s["id"], di, dom)
               for s in d["sequences"]
               for di, dom in enumerate(s["domains"])]
    n = len(domains)
    fails, ins_total, d_total = [], 0, 0
    for sid, di, dom in domains:
        sv = to_sv(dom["state_vector"])
        ins_total += sum(1 for (_, t), _ in sv if t == "i")
        d_total += sum(1 for (_, t), _ in sv if t == "d")
        rt = decode(encode(sv))
        if rt != sv:
            fails.append((sid, di, sv, rt))
    print(f"domains checked      : {n}")
    print(f"insert states        : {ins_total}")
    print(f"delete states        : {d_total}")
    print(f"round-trip PASS       : {n - len(fails)}/{n}")
    print(f"round-trip FAIL       : {len(fails)}/{n}")
    if fails:
        print("\nFIRST FAILURES:")
        for sid, di, sv, rt in fails[:3]:
            # show the first differing index
            for i in range(max(len(sv), len(rt))):
                a = sv[i] if i < len(sv) else None
                b = rt[i] if i < len(rt) else None
                if a != b:
                    print(f"  {sid}[{di}] idx {i}: orig={a} decoded={b}")
                    print(f"    orig slice   : {sv[max(0,i-2):i+3]}")
                    print(f"    decoded slice: {rt[max(0,i-2):i+3]}")
                    break
        sys.exit(1)
    print("\nC0 GATE: PASS — the position-tagged codec is lossless on all real domains.")

if __name__ == "__main__":
    main()
