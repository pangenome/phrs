#!/usr/bin/env python3
"""Drop PAF records that lie mostly in terminal telomere repeat (does (TTAGGG)n push long alignments out of sweepGA N:N?).

usage: zcat RAW.paf.gz | paf_drop_telomere.py SAMTOOLS QUERY_FA TARGET_FA OUT_BED | gzip > OUT.paf.gz
Telomere tract per contig end: from the end inward in 500 bp windows while (CCCTAA + TTAGGG) x 6 >= 50% of the window
(first 50 kb of the p end, last 50 kb of the q end, read with samtools faidx). A record is dropped when >= 50% of its
query span lies in a query tract, or >= 50% of its target span lies in a target tract.
out: filtered PAF on stdout, tracts in OUT_BED (contig, start, end, end label); counts on stderr.
"""
import subprocess, sys
S, qf, tf, bed = sys.argv[1:5]; W, E = 500, 50000
def tracts(fa):
    t = {}
    for l in open(fa + ".fai"):
        c, L = l.split("\t")[0], int(l.split("\t")[1]); n = min(E, L)
        for lab, a in (("p", 1), ("q", L - n + 1)):
            s = "".join(subprocess.run([S, "faidx", fa, f"{c}:{a}-{a + n - 1}"], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
            k = 0
            while k + W <= len(s):
                w = s[k:k + W] if lab == "p" else s[len(s) - k - W:len(s) - k]
                if (w.count("CCCTAA") + w.count("TTAGGG")) * 6 < 0.5 * W: break
                k += W
            if k: t.setdefault(c, []).append((0, k, lab) if lab == "p" else (L - k, L, lab))
    return t
tq, tt = tracts(qf), tracts(tf)
with open(bed, "w") as o:
    for t in (tq, tt):
        for c, v in t.items():
            for a, b, lab in v: o.write(f"{c}\t{a}\t{b}\t{lab}\n")
def frac(t, c, a, b): return sum(max(0, min(b, y) - max(a, x)) for x, y, _ in t.get(c, [])) / max(1, b - a)
n = d = 0
for l in sys.stdin:
    f = l.split("\t"); n += 1
    if frac(tq, f[0], int(f[2]), int(f[3])) >= 0.5 or frac(tt, f[5], int(f[7]), int(f[8])) >= 0.5: d += 1; continue
    sys.stdout.write(l)
print(f"tracts: query {sum(len(v) for v in tq.values())} ends, target {sum(len(v) for v in tt.values())} ends; records {n}, dropped {d} ({100 * d / max(1, n):.2f}%)", file=sys.stderr)
