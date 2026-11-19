#!/usr/bin/env python3
"""Per 2 kb window, same-chrom and interchrom winner of one class-winner table."""
import csv, gzip, sys
f, chrom, s0, s1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
w = {}
with gzip.open(f, "rt") as h:
    for r in csv.DictReader(h, delimiter="\t"):
        if r["chrom"] == chrom and s0 <= int(r["start"]) < s1:
            w.setdefault(int(r["start"]), {})[r["winner_class"]] = (r["other_seq"].split("#")[-1], float(r["estimated.identity"]))
print("start\tsame_chrom\tsame_id\tinterchrom\tinter_id\tsame_minus_inter")
for s in sorted(w):
    a, b = w[s].get("same_chrom", ("-", 0.0)), w[s].get("interchrom", ("-", 0.0))
    print(f"{s}\t{a[0]}\t{a[1]:.4f}\t{b[0]}\t{b[1]:.4f}\t{a[1] - b[1]:+.4f}")
