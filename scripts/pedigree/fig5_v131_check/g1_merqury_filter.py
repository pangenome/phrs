#!/usr/bin/env python3
"""Keep pieces that do not add merqury error k-mers.

usage: g1_merqury_filter.py WASHU SAMPLE BEFORE_TAG AFTER_TAG CORRECTIONS_TSV OUT_PREFIX
Per applied piece: merqury error k-mers (merqury/<tag>/<SAMPLE>/<SAMPLE>_hap{1,2}_only.bed.gz, k-mer start positions)
in [start - 1 kb, end + 1 kb) on the input assembly (BEFORE_TAG) and in [new_start - 1 kb, new_end + 1 kb) on the
corrected one (AFTER_TAG; polish moves positions by a few bp). keep = after <= before (rule: merqury must not rise).
out: OUT_PREFIX.tsv (locus, contig, start, end, mode, kmers_before, kmers_after, delta, keep), OUT_PREFIX.keep.txt
"""
import bisect, collections, gzip, sys

D, sm, bt, at, cf, out = sys.argv[1:7]
def load(t):
    d = collections.defaultdict(list)
    for h in (1, 2):
        for l in gzip.open(f"{D}/merqury/{t}/{sm}/{sm}_hap{h}_only.bed.gz", "rt"):
            c, s = l.split("\t")[:2]; d[c].append(int(s))
    for c in d: d[c].sort()
    return d
def cnt(d, c, a, z):
    p = d.get(c, []); return bisect.bisect_left(p, z) - bisect.bisect_left(p, a)
b, a = load(bt), load(at)
rows = [l.rstrip("\n").split("\t") for l in open(cf)]
rows = [dict(zip(rows[0], r)) for r in rows[1:] if r[0] == "applied"]
n = collections.Counter(); keep = []
with open(f"{out}.tsv", "w") as o:
    o.write("locus\tcontig\tstart\tend\tmode\tkmers_before\tkmers_after\tdelta\tkeep\n")
    for r in rows:
        c = r["contig"]; k0 = cnt(b, c, int(r["start"]) - 1000, int(r["end"]) + 1000)
        k1 = cnt(a, c, int(r["new_start"]) - 1000, int(r["new_end"]) + 1000); k = k1 <= k0
        n["keep" if k else "drop"] += 1; n["delta_kept"] += (k1 - k0) * k; n["delta_dropped"] += (k1 - k0) * (not k)
        if k: keep.append(r["locus"])
        o.write(f"{r['locus']}\t{c}\t{r['start']}\t{r['end']}\t{r['mode']}\t{k0}\t{k1}\t{k1 - k0}\t{'yes' if k else 'no'}\n")
open(f"{out}.keep.txt", "w").write("".join(f"{x}\n" for x in keep))
print(f"{sm}: {len(rows)} pieces, keep {n['keep']} (k-mer change {n['delta_kept']}), drop {n['drop']} (k-mer change {n['delta_dropped']})")
