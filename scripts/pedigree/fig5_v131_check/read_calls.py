#!/usr/bin/env python3
"""Per-read majority base at diagnostic sites, split by MAPQ.

usage: read_calls.py SAM SITES_TSV TRACT_START [MINQ] [PER_READ_TSV]
SAM: columns 1-6 and 10; SITES_TSV: pos, base A, base B
"""
import re, sys
sam, tsv, s0 = sys.argv[1], sys.argv[2], int(sys.argv[3]); minq = int(sys.argv[4]) if len(sys.argv) > 4 else 20
per = open(sys.argv[5], "w") if len(sys.argv) > 5 else None  # name, mapq, n_A, n_B
sites = {int(p): (a, b) for p, a, b in (l.split() for l in open(tsv) if l.strip())}

def bases_at(pos, cigar, seq):
    out, r, q = {}, pos, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            for p in sites:
                if r <= p < r + n: out[p] = seq[q + p - r]
            r += n; q += n
        elif op in "DN": r += n
        elif op in "IS": q += n
    return out

rows = {True: [0, 0, 0, 0], False: [0, 0, 0, 0]}  # A, B, tie, A with left anchor >= 5 kb
for l in open(sam):
    f = l.rstrip("\n").split("\t"); b = bases_at(int(f[3]), f[5], f[6])
    a = sum(b.get(p) == x for p, (x, y) in sites.items()); c = sum(b.get(p) == y for p, (x, y) in sites.items())
    if per: per.write(f"{f[0]}\t{f[4]}\t{a}\t{c}\n")
    hi = int(f[4]) >= minq; r = rows[hi]; r[0 if a > c else 1 if c > a else 2] += 1
    if a > c and s0 - int(f[3]) >= 5000: r[3] += 1
for hi in (True, False):
    r = rows[hi]
    print(f"{sam}\tMAPQ{'>=' if hi else '<'}{minq}\treads={sum(r[:3])}\tmostly_A={r[0]}\tmostly_B={r[1]}\ttie={r[2]}\tA_left_anchor_5kb={r[3]}")
