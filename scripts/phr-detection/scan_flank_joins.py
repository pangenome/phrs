#!/usr/bin/env python3
"""Scaffolding joins inside the telomere-anchored flanks: N runs, and telomeric arrays away from the telomere end.

A chromosome end ends in its telomere. A flank that holds an internal array of telomeric repeat (>= --min-tel bp of 1 kb
bins with >= 90% TTAGGG/CCCTAA 6-mers) or an N run is a candidate join of two separate pieces, the signature of the
NA18982 chr18q / X PAR1 scaffolding chimera. Flank names: contig:start-end_chrN_parm|qarm (the telomere was trimmed,
so position 0 of a p flank and the last base of a q flank face it). Output TSV per flank with N runs or internal
telomeric arrays; summary on stdout."""
import argparse, gzip, re, sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("fasta")
ap.add_argument("--min-tel", type=int, default=2000)
ap.add_argument("--edge", type=int, default=2000, help="bins this close to the telomere end are the end itself")
ap.add_argument("--out", required=True)
a = ap.parse_args()
TEL = re.compile("TTAGGG|CCCTAA")
out = open(a.out, "w")
out.write("flank\tlength\tN_runs\tinternal_telomeric_kb\tfirst_internal_telomeric_bin_from_telomere\n")
n = hit = 0


def check(name, seq):
    global n, hit
    n += 1
    s = seq.upper(); L = len(s)
    nr = [(m.start(), m.end() - m.start()) for m in re.finditer("N+", s)]
    p = name.endswith("_parm")
    tel_bins = []
    for i in range(0, L - 999, 1000):
        d = i if p else L - i - 1000  # distance of the bin from the telomere end
        if d < a.edge:
            continue
        if len(TEL.findall(s[i:i + 1000])) * 6 >= 900:
            tel_bins.append(d)
    if nr or len(tel_bins) * 1000 >= a.min_tel:
        hit += 1
        out.write(f"{name}\t{L}\t{','.join(f'{x}:{y}' for x, y in nr) or '.'}\t{len(tel_bins)}\t{min(tel_bins) if tel_bins else '.'}\n")


name, chunks = None, []
for l in gzip.open(a.fasta, "rt"):
    if l.startswith(">"):
        if name:
            check(name, "".join(chunks))
        name, chunks = l[1:].split()[0], []
    else:
        chunks.append(l.strip())
if name:
    check(name, "".join(chunks))
print(f"{n} flanks scanned; {hit} with an N run or >= {a.min_tel} bp of internal telomeric array (written to {a.out})")
