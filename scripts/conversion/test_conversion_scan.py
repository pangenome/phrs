#!/usr/bin/env python3
"""Known-answer tests for conversion_scan.py: planted tract, no tract, edge tract, low coverage,
reverse strand, ties, split records, and one lone identical paralog (the minimum-versus-median bias of
the first version, which called it). Run: python3 scripts/conversion/test_conversion_scan.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conversion_scan import scan

L, B = 20000, 2000
T = "HG0#1#ctg:0-20000_chr4_qarm"


def cigar(mism, start=0, end=L, dels=()):
    """'=' runs with one X at each mismatch position and a D run over each (s, e) in dels."""
    mism = [p for p in mism if not any(s <= p < e for s, e in dels)]
    ops, pos, events = [], start, sorted([(p, "X", 1) for p in mism] + [(s, "D", e - s) for s, e in dels])
    for p, op, n in events:
        if p > pos:
            ops.append(f"{p - pos}=")
        ops.append(f"{n}{op}")
        pos = p + n
    if end > pos:
        ops.append(f"{end - pos}=")
    return "".join(ops)


def mism_per_bin(counts):
    """counts[b] mismatches spread evenly in bin b"""
    out = []
    for b, n in enumerate(counts):
        out += [b * B + 1 + i * (B // (n + 1)) for i in range(n)] if n else []
    return out


def rec(q, cg, strand="+", ts=0, te=L):
    return f"{q}\t{L}\t{ts}\t{te}\t{strand}\t{T}\t{L}\t{ts}\t{te}\t0\t{te - ts}\t60\tcg:Z:{cg}\n"


def panel(plant=None, n_alleles=30, p1_strand="+", low_cov_bins=(), tie=False, lone=None, n_para=25):
    """plant: bins where the target carries the chr10q sequence (every chr10q copy identical to it);
    lone: bins where only one chr10q copy is identical to the target and the alleles differ by 1 site"""
    allele_counts = [0 if tie else 3] * (L // B)
    para_counts = [0 if tie else 20] * (L // B)
    p1_counts = list(para_counts)
    for b in plant or ():
        allele_counts[b], p1_counts[b], para_counts[b] = 20, 0, 0
    for b in lone or ():
        allele_counts[b], p1_counts[b] = 1, 0
    lines = [rec(T, cigar([]))]  # self record, must be dropped
    for i in range(n_alleles):
        dels = [(b * B, (b + 1) * B) for b in low_cov_bins] if i >= 10 else []
        lines.append(rec(f"HG{i + 1}#1#c:0-{L}_chr4_qarm", cigar(mism_per_bin(allele_counts), dels=dels)))
    lines.append(rec(f"HGp1#1#c:0-{L}_chr10_qarm", cigar(mism_per_bin(p1_counts)), strand=p1_strand))
    for i in range(2, n_para + 1):
        lines.append(rec(f"HGp{i}#1#c:0-{L}_chr10_qarm", cigar(mism_per_bin(para_counts))))
    return lines


def run(lines):
    return scan(lines, binsize=B, min_cov=0.75, min_alleles=20, min_arm_copies=20, min_delta=0.0025)


def check(name, cond, detail):
    print(f"{'PASS' if cond else 'FAIL'} {name}: {detail}")
    return cond


ok = True
_, _, bins, tr, st = run(panel(plant=(4, 5)))
ok &= check("planted embedded tract", len(tr) == 1 and (tr[0]["start"], tr[0]["end"], tr[0]["embedded"],
            tr[0]["paralog_arm"]) == (8000, 12000, 1, "chr10_qarm") and st["dropped_self"] == 1,
            [(t["start"], t["end"], t["embedded"], t["paralog_arm"]) for t in tr])
_, _, _, tr, _ = run(panel())
ok &= check("no plant, no tract", tr == [], tr)
_, _, _, tr, _ = run(panel(plant=(0, 1)))
ok &= check("edge tract is not embedded", len(tr) == 1 and tr[0]["embedded"] == 0 and tr[0]["start"] == 0,
            [(t["start"], t["embedded"]) for t in tr])
_, _, bins, tr, _ = run(panel(plant=(4, 5), low_cov_bins=(4, 5)))
ok &= check("low allele coverage, no call", tr == [] and bins[4]["n_alleles"] == 10, (tr, bins[4]["n_alleles"]))
_, _, _, tr, _ = run(panel(plant=(4, 5), p1_strand="-"))
ok &= check("reverse-strand donor", len(tr) == 1 and tr[0]["start"] == 8000, [(t["start"], t["end"]) for t in tr])
_, _, bins, tr, _ = run(panel(tie=True))
ok &= check("identical copies, no call", tr == [] and bins[3]["d_allele_median"] == bins[3]["d_arm_median"],
            (tr, bins[3]["d_allele_median"], bins[3]["d_arm_median"]))
_, _, bins, tr, _ = run(panel(lone=(4, 5)))
ok &= check("one identical paralog among many, no call", tr == [] and bins[4]["d_paralog_min"] == 0.0,
            (tr, bins[4]["d_paralog_min"]))
_, _, _, tr, _ = run(panel(plant=(4, 5), n_para=10))
ok &= check("fewer paralog copies than min_arm_copies, no call", tr == [], tr)
lines = panel(plant=(4, 5))
lines.append(rec("HG1#1#c:0-20000_chr4_qarm", "1000=", ts=8000, te=9000))  # shorter second record, same query
_, _, bins, tr, _ = run(lines)
ok &= check("split record keeps the longer one", bins[4]["n_alleles"] == 30 and len(tr) == 1, bins[4]["n_alleles"])
sys.exit(0 if ok else 1)
