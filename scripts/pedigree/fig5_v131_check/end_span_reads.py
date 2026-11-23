#!/usr/bin/env python3
"""Length from an anchor to the q-end telomere in reads that cross the anchor (array length check at a chromosome end).

usage: samtools view -F 0x904 BAM... CONTIG:A-A | end_span_reads.py A LABEL > TSV
For each primary read whose alignment covers reference position A (unique sequence before the end array): query position
at A (CIGAR walk on SEQ, which SAM stores in reference orientation), bases from there to the read end (toward the q
telomere), and whether the last 2 kb of SEQ is >= 50% TTAGGG (the read reaches the telomere).
out: label, read, mapq, bases_anchor_to_read_end, telomere_at_end (0/1), read_length
"""
import re, sys
A, lab = int(sys.argv[1]), sys.argv[2]
print("label\tread\tmapq\tanchor_to_end\ttelomere_at_end\tread_length")
n = 0
for l in sys.stdin:
    f = l.split("\t"); r, q = int(f[3]), 0  # r: 1-based reference position
    for k, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
        k = int(k)
        if op in "M=X":
            if r + k > A: q += max(0, A - r); break
            r += k; q += k
        elif op in "DN": r += k
        elif op in "IS": q += k
    else: continue
    s = f[9]; t = s[-2000:]
    tel = int(len(t) >= 2000 and t.count("TTAGGG") * 6 >= 0.5 * len(t))
    print(f"{lab}\t{f[0]}\t{f[4]}\t{len(s) - q}\t{tel}\t{len(s)}"); n += 1
print(f"{lab}: {n} reads cross {A}", file=sys.stderr)
