#!/usr/bin/env python3
"""g1_locus task list for the end sweep: every contig end of one person, with reads of its paralog partner ends.

usage: endsweep_tasks.py NORM_DIR SAMPLE PARTNERS_TSV WINDOW_BP MAX_PARTNERS OUT_TSV
Per contig of both haplotypes, p end = first WINDOW_BP, q end = last WINDOW_BP (error = read error = that window, on
norm_V). Extra read regions = the 500 kb ends (end_partners.py) of its paralog partners, the MAX_PARTNERS with the most
aligned bp (the DAD chr9 control passed with 10 partners). Ends with no partner get no extra column value.
out: id, name, error, read error, extra (g1_locus_array.sbatch 5-column TSV)
"""
import collections, sys
nd, sm, pf, W, K, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
E = 500000
L = {}
for h in (1, 2):
    for l in open(f"{nd}/{sm}#{h}.fa.gz.fai"):
        c, n = l.split("\t")[:2]; L[c] = int(n)
par = collections.defaultdict(list)
for i, l in enumerate(open(pf)):
    t = l.rstrip("\n").split("\t")
    if i == 0: continue
    if t[4] == "paralog": par[t[0]].append((int(t[2]), t[1]))
def reg(end):
    c, a = end.rsplit(":", 1)
    if c not in L: sys.exit(f"ERROR partner contig {c} not in {nd}/{sm}#*.fa.gz.fai")
    return f"{c}:1-{min(E, L[c])}" if a == "p" else f"{c}:{max(1, L[c] - E + 1)}-{L[c]}"
n, capped, none = 0, 0, 0
with open(out, "w") as o:
    for c in L:
        if L[c] < 2 * W: sys.exit(f"ERROR contig {c} length {L[c]} < 2 x window {W}")
        for a in ("p", "q"):
            r = f"{c}:1-{W}" if a == "p" else f"{c}:{L[c] - W + 1}-{L[c]}"
            p = sorted(par[f"{c}:{a}"], reverse=True)
            capped += len(p) > K; none += not p
            n += 1; nm = f"{c.split('#')[2]}h{c.split('#')[1]}{a}_w{W // 1000}k"
            o.write(f"{n}\t{nm}\t{r}\t{r}\t{','.join(reg(e) for _, e in p[:K])}\n")
print(f"{sm}: {n} ends, window {W}; ends with no partner {none}, with > {K} partners (capped) {capped}", file=sys.stderr)
