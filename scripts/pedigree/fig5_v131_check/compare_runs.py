#!/usr/bin/env python3
"""Compare two Fig5 drawn-run tables (whole_genome_ribbon_runs.tsv)."""
import csv, sys

def read(p):
    with open(p) as f:
        return [(r["query_chrom"], int(r["query_start"]), int(r["query_end"]), r["target_chrom"], r["donor_haplotype"],
                 float(r["mean_inter_identity"]), r["category"]) for r in csv.DictReader(f, delimiter="\t")]

def overlap(x, y):
    return x[0] == y[0] and x[1] < y[2] and y[1] < x[2]

a, b = read(sys.argv[1]), read(sys.argv[2])
print(f"A runs {len(a)}, B runs {len(b)}")
for tag, s, o in (("A_only", a, b), ("B_only", b, a)):
    n = [x for x in s if not any(overlap(x, y) for y in o)]
    print(f"{tag}: {len(n)}")
    for x in sorted(n):
        print(tag, x[0], x[1], x[2], x[2] - x[1], x[3], x[4], f"{x[5]:.4f}", x[6], sep="\t")
