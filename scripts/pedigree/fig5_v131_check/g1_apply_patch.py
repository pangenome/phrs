#!/usr/bin/env python3
"""Replace one assembly interval with a contig piece.

usage: g1_apply_patch.py IN_FA.gz OUT_FA CONTIG T0 T1 PIECE_FA PIECE_NAME Q0 Q1 STRAND
Writes OUT_FA (plain FASTA, all records of IN_FA) with CONTIG[T0:T1] (0-based, end excluded) replaced by
PIECE_NAME[Q0:Q1] (reverse-complemented when STRAND is "-"). Prints old and new lengths and the md5 of the new piece.
"""
import gzip, hashlib, sys

inf, outf, ctg, t0, t1, pf, pn, q0, q1, st = sys.argv[1:11]
t0, t1, q0, q1 = int(t0), int(t1), int(q0), int(q1)
def read(f):
    op = gzip.open if f.endswith(".gz") else open
    name, seq = None, []
    with op(f, "rt") as h:
        for l in h:
            if l.startswith(">"):
                if name: yield name, "".join(seq)
                name, seq = l[1:].split()[0], []
            else: seq.append(l.strip())
    if name: yield name, "".join(seq)
piece = next(s for n, s in read(pf) if n == pn)[q0:q1]
if st == "-": piece = piece[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))
done = False
with open(outf, "w") as o:
    for n, s in read(inf):
        if n == ctg:
            new = s[:t0] + piece + s[t1:]; done = True
            print(f"{ctg}: replaced [{t0},{t1}) = {t1 - t0} bp by {len(piece)} bp; length {len(s)} -> {len(new)}; "
                  f"piece md5 {hashlib.md5(piece.upper().encode()).hexdigest()}", file=sys.stderr)
            s = new
        o.write(f">{n}\n")
        for i in range(0, len(s), 80): o.write(s[i:i + 80] + "\n")
if not done: sys.exit(f"ERROR contig {ctg} not in {inf}")
