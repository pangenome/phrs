#!/usr/bin/env python3
"""How close to each contig end the homolog alignment reaches, in a sweepGA PAF (child haplotype vs joint parent).

usage: zcat X.paf.gz | sweepga_end_reach.py LABEL > TSV
Homolog = target chromosome equal to the query chromosome (target names PARENT#joint#hH_chrN; chrX query also chrY), per parent
copy h1 and h2: the transmitted copy must reach the end (chr1p/chr5q: SELF.p chr5 vs DAD h2 stopped 846 kb short, h1 did not).
Per query contig and end: bases from the end to the first homolog-aligned base (uncovered at the end), homolog-covered
bp in the last 1 Mb, and the same for any target (all chromosomes). An end where the homolog stops early while other
chromosomes reach it is where a non-homologous winner can appear by default (chr1p/chr5q at FastGA -f32).
"""
import collections, sys
lab = sys.argv[1]; M = 1000000
L = {}; iv = collections.defaultdict(list); n = 0
for l in sys.stdin:
    f = l.split("\t", 9); q, ql, qs, qe, t = f[0], int(f[1]), int(f[2]), int(f[3]), f[5]
    n += 1; L[q] = ql; qc = q.split("#")[-1]; tc = t.split("_")[-1]
    hom = tc == qc or {qc, tc} == {"chrX", "chrY"}; th = t.split("#")[-1].split("_")[0]  # parent copy: h1 or h2
    for end, lo, hi in (("p", 0, min(M, ql)), ("q", max(0, ql - M), ql)):
        a, b = max(qs, lo), min(qe, hi)
        if a < b:
            iv[(q, end, "all")].append((a, b))
            if hom: iv[(q, end, "hom_" + th)].append((a, b))
def union(x):
    out = []
    for a, b in sorted(x):
        if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out
print("label\tquery\tlength\tend\th1_uncovered_at_end\th1_covered_bp_1mb\th2_uncovered_at_end\th2_covered_bp_1mb\tall_uncovered_at_end\tall_covered_bp_1mb")
for q in sorted(L):
    for end in ("p", "q"):
        r = []
        for k in ("hom_h1", "hom_h2", "all"):
            u = union(iv[(q, end, k)])
            if not u: r += [min(M, L[q]), 0]; continue
            gap = u[0][0] if end == "p" else L[q] - u[-1][1]
            r += [gap, sum(b - a for a, b in u)]
        print(f"{lab}\t{q}\t{L[q]}\t{end}\t" + "\t".join(map(str, r)))
print(f"{lab}: {n} PAF records", file=sys.stderr)
