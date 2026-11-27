#!/usr/bin/env python3
"""Reads at a suspected assembly join (samtools view SAM text on stdin, reads aligned to the assembly).

For a join [--j0, --j1] (1-based, contig coordinates) we count primary alignments that SPAN it (start before j0 - anchor
and end after j1 + anchor) and alignments that STOP at it: end within --near bp of j0 with a right clip >= --min-clip, or
start within --near bp of j1 with a left clip >= --min-clip. For each stopping read we write its clipped sequence
(in contig orientation, the part beyond the join) to --fasta, with its length, the share of its bases in telomeric
6-mers (TTAGGG / CCCTAA), the same share in its last 2 kb, and the other alignments of the read (SA tag).
If the chromosome really continued across the join, long reads would span it; if the join is false, reads stop there
and their clipped part is the sequence that truly follows."""
import argparse, re, sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--j0", type=int, required=True)
ap.add_argument("--j1", type=int, required=True)
ap.add_argument("--anchor", type=int, default=5000)
ap.add_argument("--near", type=int, default=3000)
ap.add_argument("--min-clip", type=int, default=2000)
ap.add_argument("--min-mapq", type=int, default=0)
ap.add_argument("--fasta", required=True)
a = ap.parse_args()
OP = re.compile(r"(\d+)([MIDNSHP=X])")
tel = lambda s: (len(re.findall("TTAGGG", s)) + len(re.findall("CCCTAA", s))) * 6 / max(len(s), 1)
span, stopL, stopR, n = [], [], [], 0
fa = open(a.fasta, "w")
for l in sys.stdin:
    f = l.rstrip("\n").split("\t")
    flag, pos, mq, cig, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
    if flag & 0x904 or mq < a.min_mapq:  # primary, mapped
        continue
    n += 1
    ops = [(int(k), o) for k, o in OP.findall(cig)]
    end = pos + sum(k for k, o in ops if o in "MDN=X") - 1
    lc = ops[0][0] if ops[0][1] == "S" else 0
    rc = ops[-1][0] if ops[-1][1] == "S" else 0
    sa = next((t[5:] for t in f[11:] if t.startswith("SA:Z:")), "")
    others = ";".join(x.split(",")[0].split("#")[-1] + ":" + x.split(",")[1] + x.split(",")[2] + ",mq" + x.split(",")[4]
                      for x in sa.split(";") if x)
    if pos < a.j0 - a.anchor and end > a.j1 + a.anchor:
        span.append((f[0], pos, end, mq))
    elif abs(end - a.j0) <= a.near and rc >= a.min_clip:
        c = seq[len(seq) - rc:]
        stopL.append((f[0], pos, end, mq, rc, tel(c), tel(c[-2000:]), others))
        fa.write(f">{f[0]}|stops_before_join|end{end}\n{c}\n")
    elif abs(pos - a.j1) <= a.near and lc >= a.min_clip:
        c = seq[:lc]
        stopR.append((f[0], pos, end, mq, lc, tel(c), tel(c[:2000]), others))
        fa.write(f">{f[0]}|starts_after_join|start{pos}\n{c}\n")
print(f"primary alignments: {n}; spanning the join with {a.anchor} bp on both sides: {len(span)}")
for x in span[:10]:
    print("  span", x)
for lab, v in (("stop before the join (clipped part = what follows on the reads)", stopL),
               ("start after the join (clipped part = what precedes on the reads)", stopR)):
    print(f"{lab}: {len(v)}")
    for x in v:
        print(f"  {x[0][:12]} aln {x[1]}-{x[2]} mapq {x[3]} clip {x[4]} bp, telomeric share {x[5]:.2f} (far 2 kb {x[6]:.2f}), other alignments: {x[7] or '-'}")
