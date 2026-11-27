#!/usr/bin/env python3
"""Contig coordinates of the candidate joins in the flagged flanks (scan_flank_joins.py TSV): every N run and every
internal telomeric array (runs of 1 kb bins with >= 90% telomeric 6-mers, >= 2 kb from the telomere end) of each flank,
from the flank sequence (samtools faidx of the flank FASTA). Flank names contig:a-b_chrN_arm are 1-based, so flank
index i is contig position a + i. Output: one row per feature with sample, haplotype, contig, 1-based start and end,
kind (N or TEL), its strand class for TEL (SAME = same orientation as the flank's own telomere), bp from the telomere end."""
import argparse, re, subprocess

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--fasta", required=True)
ap.add_argument("--flagged", required=True)
ap.add_argument("--samtools", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
NAME = re.compile(r"^(([^#]+)#([^#]+)#[^:]+):(\d+)-(\d+)_(chr[0-9XY]+)_(parm|qarm)$")
out = open(a.out, "w")
out.write("flank\tsample\thap\tcontig\tstart\tend\tkind\tstrand_vs_end\tdist_from_telomere\n")
n = 0
for l in open(a.flagged):
    if l.startswith("flank\t"):
        continue
    nm = l.split("\t")[0]
    m = NAME.match(nm)
    if not m:
        raise SystemExit(f"FATAL: cannot parse {nm}")
    ctg, smp, hap, fa, fb, arm = m.group(1), m.group(2), m.group(3), int(m.group(4)), int(m.group(5)), m.group(7)
    r = subprocess.run([a.samtools, "faidx", a.fasta, nm], capture_output=True, text=True, check=True)
    s = "".join(r.stdout.split("\n")[1:]).upper(); L = len(s); p = arm == "parm"
    if L != fb - fa + 1:
        raise SystemExit(f"FATAL: {nm} length {L} differs from its name")
    dist = lambda i0, i1: i0 if p else L - i1  # bp from the telomere end to the near edge of [i0, i1)
    feats = [(x.start(), x.end(), "N", ".") for x in re.finditer("N+", s)]
    arr = []
    for i in range(0, L - 999, 1000):
        b = s[i:i + 1000]; g, c = len(re.findall("TTAGGG", b)), len(re.findall("CCCTAA", b))
        if (g + c) * 6 >= 900 and dist(i, i + 1000) >= 2000:
            st = "G" if g > c else "C"
            if arr and i - arr[-1][1] <= 1000 and arr[-1][2] == st:
                arr[-1][1] = i + 1000
            else:
                arr.append([i, i + 1000, st])
    own = "C" if p else "G"
    feats += [(i0, i1, "TEL", "SAME" if st == own else "OPPOSITE") for i0, i1, st in arr]
    for i0, i1, k, sv in sorted(feats):
        out.write(f"{nm}\t{smp}\t{hap}\t{ctg}\t{fa + i0}\t{fa + i1 - 1}\t{k}\t{sv}\t{dist(i0, i1)}\n"); n += 1
print(f"{n} features written to {a.out}")
