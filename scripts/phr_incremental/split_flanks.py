#!/usr/bin/env python3
"""Split the telomere-anchored flank FASTA (bgzip or plain) into one single-line FASTA per flank, <sanitized>.fa
with its .fai, and write haplotypes.tsv (haplotype = name up to the second '#', flank name)."""
import gzip, os, sys

fa, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
san = lambda s: s.translate(str.maketrans(":@#-", "____"))
n = 0
hap = open(os.path.join(out, "haplotypes.tsv"), "w")


def flush(name, chunks):
    global n
    seq = "".join(chunks)
    p = os.path.join(out, san(name) + ".fa")
    with open(p, "w") as fh:
        fh.write(f">{name}\n{seq}\n")
    with open(p + ".fai", "w") as fh:
        fh.write(f"{name}\t{len(seq)}\t{len(name) + 2}\t{len(seq)}\t{len(seq) + 1}\n")
    hap.write("#".join(name.split("#")[:2]) + "\t" + name + "\n")
    n += 1


name, chunks = None, []
for l in (gzip.open if fa.endswith(".gz") else open)(fa, "rt"):
    if l.startswith(">"):
        if name:
            flush(name, chunks)
        name, chunks = l[1:].split()[0], []
    else:
        chunks.append(l.strip())
if name:
    flush(name, chunks)
hap.close()
sys.stderr.write(f"{n} flanks written to {out}\n")
