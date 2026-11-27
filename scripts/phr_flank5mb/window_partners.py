#!/usr/bin/env python3
"""Partners of one 5 kb caller window, from the target's own all-vs-one PAF: records that overlap the window by
>= 3 kb with matches / block length >= 0.95 (an approximation of the impg query, which also projects through the
other PAFs). Per partner chromosome (from the _chrN_arm suffix): records, and how many come from contigs whose own
name holds 'chr' (CHM13, HG002), which the March caller counted under the contig name instead of the chromosome."""
import argparse, gzip, os
from collections import Counter

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--paf-dir", required=True)
ap.add_argument("--prefix", default="")
ap.add_argument("--seq", required=True)
ap.add_argument("--dist", type=int, required=True, help="window start, bp from the telomere")
ap.add_argument("--len", type=int, required=True, help="flank length")
a = ap.parse_args()
s, e = (a.dist, a.dist + 5000) if a.seq.endswith("_parm") else (a.len - a.dist - 5000, a.len - a.dist)
f = os.path.join(a.paf_dir, a.prefix + a.seq.translate(str.maketrans(":@#-", "____")) + ".paf.gz")
n, named = Counter(), Counter()
for l in gzip.open(f, "rt"):
    x = l.split("\t")
    ov = min(e, int(x[8])) - max(s, int(x[7]))
    if ov >= 3000 and int(x[9]) / int(x[10]) >= 0.95:
        ch = x[0].rsplit("_", 2)[-2]
        n[ch] += 1
        named[ch] += "chr" in x[0].split("#")[-1].split(":")[0]
print(f"{a.seq} window [{s}, {e}) ({a.dist} bp from the telomere): " +
      ", ".join(f"{c} {k} ({named[c]} from contigs named chr*)" for c, k in n.most_common()))
