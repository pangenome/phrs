#!/usr/bin/env python3
"""How far from the telomere do acrocentric p arms and the X PAR1 keep sharing sequence with other chromosomes, in the
genome-wide survey (Fig. 1A input: 100 kb CHM13 windows, IMPG queries with up to five hops at >= 98% identity)?
A window is SHARED when >= --min-haps haplotypes reach a chromosome other than its own (chrUn excluded). Per end, from
the telomere to the centromere start (data/chm13-annotations.bed): shared windows per 1 Mb bin, the last shared window
before the first run of --stop unshared windows, and the last shared window before the centromere."""
import argparse, gzip
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--windows", default="data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz")
ap.add_argument("--annot", default="data/chm13-annotations.bed")
ap.add_argument("--min-haps", type=int, default=47)
ap.add_argument("--stop", type=int, default=4)
a = ap.parse_args()
cen = {}
for l in open(a.annot):
    c, s, e, k = l.rstrip("\n").split("\t")[:4]
    if k == "Centromere":
        cen[c] = int(s)
sh = {}
for l in gzip.open(a.windows, "rt"):
    if l.startswith("chrom"):
        continue
    f = l.rstrip("\n").split("\t")
    c = f[0].split("#")[-1]
    parts = [] if len(f) < 8 or not f[7] else [p.rsplit("-", 1) for p in f[7].split(",")]
    sh[(c, int(f[1]))] = any(x not in (c, "chrUn") and int(n) >= a.min_haps for x, n in parts)
for c in ("chr13", "chr14", "chr15", "chr21", "chr22", "chrX"):
    n = cen[c] // 100000
    v = [sh.get((c, i * 100000), False) for i in range(n)]
    last_run, fail = 0, 0
    for i, x in enumerate(v):
        if x:
            last_run, fail = i + 1, 0
        else:
            fail += 1
            if fail >= a.stop:
                break
    last_any = max((i + 1 for i, x in enumerate(v) if x), default=0)
    bins = [sum(v[j:j + 10]) for j in range(0, n, 10)]
    print(f"{c}p: centromere at {cen[c] / 1e6:.2f} Mb; shared run from the telomere to {last_run / 10:.1f} Mb (stop after {a.stop} "
          f"unshared 100 kb windows); last shared window {last_any / 10:.1f} Mb; shared windows per Mb: {bins}")
