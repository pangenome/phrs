#!/usr/bin/env python3
"""Query gaps left by both parent haplotypes.

usage: g2_union.py HOMOLOG_DIR TOOL(mm|wf) OUT_TSV
Reads gaps/<tag>.<pair>.<chrom>.tsv of g2_homolog.sbatch (gaps >= 1 kb, outside N). A query interval is a
residual gap when no parent haplotype homolog covers it (intersection of the per-target query gaps).
out: one row per residual gap >= 1 kb; summary per tag, pair, chromosome on stdout.
"""
import collections, glob, os, sys

d, tool, out = sys.argv[1:4]
rows, summ = [], []
for f in sorted(glob.glob(f"{d}/gaps/*.tsv")):
    tag, pair, chrom = os.path.basename(f)[:-4].split(".")
    cov = [l.rstrip("\n").split("\t") for l in open(f"{d}/cov/{tag}.{pair}.{chrom}.tsv")][1:]
    tg = sorted({x[4] for x in cov if x[3] == tool and not x[4].endswith("chrY")})
    L = int(next(x[6] for x in cov if x[5] == "query"))
    g = collections.defaultdict(list)
    for l in list(open(f))[1:]:
        t = l.rstrip("\n").split("\t")
        if t[3] == tool and t[5] == "query" and t[4] in tg: g[t[4]].append((int(t[6]), int(t[7])))
    cur = g[tg[0]] if tg else []
    for t in tg[1:]:
        nxt = []
        for a, b in cur:
            for c, e in g[t]:
                s, z = max(a, c), min(b, e)
                if z - s >= 1000: nxt.append((s, z))
        cur = nxt
    best = min(tg, key=lambda t: sum(b - a for a, b in g[t])) if tg else "NA"
    for a, b in cur: rows.append((tag, pair, chrom, a, b, b - a, min(a, L - b), L))
    summ.append((tag, pair, chrom, L, len(tg), best, sum(b - a for a, b in g[best]) if tg else 0, sum(b - a for a, b in cur), sum(1 for a, b in cur if b - a >= 10000)))
with open(out, "w") as o:
    o.write("tag\tpair\tchrom\tstart\tend\tlen\tdist_to_end\tchrom_len\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("tag\tpair\tchrom\tlen\ttargets\tbest_target\tbest_uncovered_ge1kb\tboth_uncovered_ge1kb\tboth_gaps_ge10kb")
for s in summ: print("\t".join(map(str, s)))
tot = collections.Counter()
for s in summ: tot[(s[0], s[1])] += s[7]
for k, v in sorted(tot.items()): print(f"# {k[0]} {k[1]}: chromosomes {sum(1 for s in summ if s[:2] == k)}, residual bp {v}", file=sys.stderr)
