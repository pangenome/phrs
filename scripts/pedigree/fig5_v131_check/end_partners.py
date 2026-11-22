#!/usr/bin/env python3
"""Paralog partner ends of every contig end of one person (reads of a missing end sit on these ends).

usage: end_partners.py NORM_DIR SAMPLE MINIMAP2 THREADS WORK_DIR OUT_TSV (WORK_DIR: empty temp dir, removed by the caller)
Last (q) and first (p) 500 kb of every contig of both haplotypes, all against all (minimap2 asm20 -X, base level);
a partner = another end with >= 20 kb aligned (sum of records) at >= 0.95 identity (matches / block) to this end.
Same-chromosome same-arm ends of the other haplotype are kept and typed homolog; the rest are paralog.
out: end, partner, aligned_bp, identity (bp-weighted), partner_type (homolog | paralog)
"""
import collections, gzip, subprocess, sys
nd, sm, mm, th, T, out = sys.argv[1:7]; W = 500000
def read(f):
    name, seq = None, []
    for l in gzip.open(f, "rt"):
        if l[0] == ">":
            if name: yield name, "".join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    yield name, "".join(seq)
n = 0
with open(f"{T}/ends.fa", "w") as o:
    for h in (1, 2):
        for c, s in read(f"{nd}/{sm}#{h}.fa.gz"):
            o.write(f">{c}:p\n{s[:W]}\n>{c}:q\n{s[max(0, len(s) - W):]}\n"); n += 2
print(f"{sm} from {nd}: {n} ends of {W} bp", file=sys.stderr)
subprocess.run(f"{mm} -cx asm20 -X -t {th} -o {T}/a.paf {T}/ends.fa {T}/ends.fa", shell=True, check=True)
acc = collections.defaultdict(lambda: [0, 0])
for l in open(f"{T}/a.paf"):
    f = l.split("\t")
    if f[0] == f[5]: continue
    for a, b in ((f[0], f[5]), (f[5], f[0])): acc[(a, b)][0] += int(f[9]); acc[(a, b)][1] += int(f[10])
k = 0
with open(out, "w") as o:
    o.write("end\tpartner\taligned_bp\tidentity\tpartner_type\n")
    for (a, b), (m, bl) in sorted(acc.items()):
        if bl >= 20000 and m / bl >= 0.95:
            same = a.split(":")[0].split("#")[2] == b.split(":")[0].split("#")[2] and a.split(":")[1] == b.split(":")[1]
            o.write(f"{a}\t{b}\t{bl}\t{m / bl:.4f}\t{'homolog' if same else 'paralog'}\n"); k += 1
print(f"partner pairs kept: {k}", file=sys.stderr)
