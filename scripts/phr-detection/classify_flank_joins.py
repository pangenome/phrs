#!/usr/bin/env python3
"""Classify the flanks flagged by scan_flank_joins.py (samtools faidx of each from the flank FASTA).
For each internal telomeric array (runs of 1 kb bins with >= 90% telomeric 6-mers, away from the telomere end): its
distance from the telomere, its length, its strand (C-rich CCCTAA or G-rich TTAGGG on the forward strand) relative to
the flank's own telomere (a p end reads CCCTAA, a q end TTAGGG on the forward strand: SAME means a second chromosome end
facing the same way, so the sequence between it and the flank end lies beyond an end), and whether an N run touches it
(within --touch bp). Classes: TEL+N (array touching an N run, the scaffolding-chimera signature), TEL (array, no N run
next to it), N (N runs only)."""
import argparse, re, subprocess, sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--fasta", required=True)
ap.add_argument("--flagged", required=True, help="scan_flank_joins.py TSV")
ap.add_argument("--samtools", required=True)
ap.add_argument("--touch", type=int, default=2000)
ap.add_argument("--out", required=True)
a = ap.parse_args()
names = [l.split("\t")[0] for l in open(a.flagged) if not l.startswith("flank\t")]
out = open(a.out, "w")
out.write("flank\tclass\tarray_dist_kb\tarray_kb\tarray_strand_vs_end\tN_runs(dist_from_telomere:len)\n")
cls = {}
for nm in names:
    r = subprocess.run([a.samtools, "faidx", a.fasta, nm], capture_output=True, text=True, check=True)
    s = "".join(r.stdout.split("\n")[1:]).upper(); L = len(s); p = nm.endswith("_parm")
    dist = lambda i: i if p else L - i  # bp from the telomere end
    nr = [(min(dist(m.start()), dist(m.end())), m.end() - m.start()) for m in re.finditer("N+", s)]
    ends = [d for nd, ln in nr for d in (nd, nd + ln)]  # both edges of each N run, bp from the telomere end
    bins = []
    for i in range(0, L - 999, 1000):
        b = s[i:i + 1000]; g, c = len(re.findall("TTAGGG", b)), len(re.findall("CCCTAA", b))
        if (g + c) * 6 >= 900 and min(dist(i), dist(i + 1000)) >= 2000:
            bins.append((i, "G" if g > c else "C"))
    arrays = []
    for i, st in bins:  # merge adjacent bins of the same strand
        if arrays and i - arrays[-1][1] <= 1000 and arrays[-1][2] == st:
            arrays[-1][1] = i + 1000
        else:
            arrays.append([i, i + 1000, st])
    own = "C" if p else "G"
    rows = []
    for s0, s1, st in arrays:
        d = min(dist(s0), dist(s1))
        touch = any(abs(e - d) <= a.touch or abs(e - max(dist(s0), dist(s1))) <= a.touch for e in ends)
        rows.append((d, s1 - s0, "SAME" if st == own else "OPPOSITE", touch))
    c = "TEL+N" if any(x[3] for x in rows) else ("TEL" if rows else "N")
    cls[nm] = c
    big = max(rows, key=lambda x: x[1]) if rows else None
    out.write(f"{nm}\t{c}\t{big[0] / 1000 if big else '.'}\t{big[1] / 1000 if big else '.'}\t{big[2] if big else '.'}\t"
              f"{','.join(f'{d}:{l}' for d, l in nr) or '.'}\n")
from collections import Counter
print(f"{len(names)} flagged flanks: {dict(Counter(cls.values()))}")
