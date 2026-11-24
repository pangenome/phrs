#!/usr/bin/env python3
"""Known-answer tests for conversion_sites.py on a synthetic target-centric PAF: 30 alleles (chr4_qarm) and
25 paralogs (chr10_qarm); paralogous sequence variants every 200 bp where the paralogs differ from the
target. Run: python3 scripts/conversion/test_conversion_sites.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conversion_sites import scan

L = 20000
T = "HG0#1#ctg:0-20000_chr4_qarm"
PSV = list(range(100, L, 200))


def cigar(mism, start=0, end=L, dels=()):
    mism = [p for p in sorted(set(mism)) if start <= p < end and not any(s <= p < e for s, e in dels)]
    ops, pos = [], start
    for p, op, n in sorted([(p, "X", 1) for p in mism] + [(s, "D", e - s) for s, e in dels]):
        if p > pos:
            ops.append(f"{p - pos}=")
        ops.append(f"{n}{op}")
        pos = p + n
    if end > pos:
        ops.append(f"{end - pos}=")
    return "".join(ops)


def rec(q, cg, ts=0, te=L, strand="+"):
    return f"{q}\t{L}\t{ts}\t{te}\t{strand}\t{T}\t{L}\t{ts}\t{te}\t0\t{te - ts}\t60\tcg:Z:{cg}\n"


def panel(plant=(), n_alleles=30, n_para=25, lone=False, dup_record=False, low_cov=False):
    """plant: PSV positions where the target carries the chr10q base (alleles differ, paralogs match)"""
    lines = [rec(T, cigar([]))]
    for i in range(n_alleles):
        dels = [(7000, 10000)] if low_cov and i >= 10 else []
        lines.append(rec(f"HG{i + 1}#1#c:0-{L}_chr4_qarm", cigar(list(plant) + [5 + i * 600], dels=dels)))
    for i in range(n_para):
        same = lone and i == 0
        m = [] if same else [p for p in PSV if p not in plant]
        lines.append(rec(f"HGp{i}#1#c:0-{L}_chr10_qarm", cigar(m), strand="-" if i % 2 else "+"))
        if dup_record:  # a second, shorter record of the same query over the planted region, all mismatches
            lines.append(rec(f"HGp{i}#1#c:0-{L}_chr10_qarm", "1000X", ts=8000, te=9000))
    return lines


def check(name, cond, detail):
    print(f"{'PASS' if cond else 'FAIL'} {name}: {detail}")
    return cond


ok = True
plant = (8100, 8300, 8500)
_, _, tr, pw, st = scan(panel(plant=plant))
ok &= check("planted 3-site tract", [(t["first_site"], t["last_site"], t["n_sites"], t["left_a_site"], t["right_a_site"],
            t["embedded"], t["donor_arm"]) for t in tr] == [(8100, 8500, 3, 7900, 8700, 1, "chr10_qarm")]
            and st["dropped_self"] == 1, tr)
ok &= check("power counts", pw == [("chr10_qarm", len(PSV) - 3, 3)], pw)
_, _, tr, _, _ = scan(panel())
ok &= check("no plant, no tract", tr == [], tr)
_, _, tr, _, _ = scan(panel(plant=(8100,)))
ok &= check("one site is not a tract", tr == [], tr)
_, _, tr, _, _ = scan(panel(plant=(100, 300)))
ok &= check("tract at the first PSV is not embedded", len(tr) == 1 and tr[0]["embedded"] == 0
            and tr[0]["left_a_site"] is None, tr)
_, _, tr, _, _ = scan(panel(plant=plant, low_cov=True))
ok &= check("fewer than 20 alleles over the tract, no call", tr == [], tr)
_, _, tr, _, _ = scan(panel(lone=True))
ok &= check("one paralog identical to the target, no Y-site", tr == [], tr)
_, _, tr, _, st = scan(panel(plant=plant, dup_record=True))
ok &= check("a second record of the same query does not count twice", len(tr) == 1 and tr[0]["n_sites"] == 3,
            [(t["first_site"], t["last_site"], t["n_sites"]) for t in tr])
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
