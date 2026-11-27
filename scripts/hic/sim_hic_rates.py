#!/usr/bin/env python3
"""Hi-C null simulation, rate comparable with sim_porec_null.py evaluate: over all reads of the valid pairs (any true
window), the fraction placed in the PHR window of an end of another arm (homologous ends of the two haplotypes count
as the same arm), for pairs with MAPQ >= T on both reads (allValidPairs columns 11 and 12), T in THRESHOLDS. Joins
line i of ALL_VALID_PAIRS with rows 2i+1, 2i+2 of READS_TSV (sim_hic_null.py evaluate); the MAPQ of a row is the one
of the allValidPairs side at its assigned position.
usage: sim_hic_rates.py ALL_VALID_PAIRS READS_TSV OUT_TSV [--thresholds 0 1 2]"""
import argparse, re, sys
ap = argparse.ArgumentParser(); ap.add_argument("avp"); ap.add_argument("reads"); ap.add_argument("out")
ap.add_argument("--thresholds", type=int, nargs="+", default=[0, 1, 2]); a = ap.parse_args()
arm = lambda x: re.sub(r"_(MATERNAL|PATERNAL)", "", x)
tot = {t: 0 for t in a.thresholds}; mis = dict(tot); homo = dict(tot); pairs = dict(tot); n = 0
with open(a.avp) as f, open(a.reads) as r:
    if not r.readline().startswith("read\t"): sys.exit("FATAL reads.tsv header")
    for l in f:
        x = l.rstrip("\n").split("\t"); q = {(x[1], x[2]): int(x[10]), (x[4], x[5]): int(x[11])}; n += 1
        rows = [r.readline().rstrip("\n").split("\t") for _ in (0, 1)]
        if any(y[0] != x[0].split(":")[0] for y in rows): sys.exit(f"FATAL read order at pair {n}: {x[0]} vs {rows[0][0]}")
        for t in a.thresholds:
            if min(int(x[10]), int(x[11])) < t: continue
            pairs[t] += 1
            for y in rows:
                if (y[4], y[5]) not in q: sys.exit(f"FATAL side {y[4]}:{y[5]} not in pair {x[0]}")
                tot[t] += 1
                if y[7] == "phr" and y[6] != y[2]:
                    if arm(y[6]) == arm(y[2]): homo[t] += 1
                    else: mis[t] += 1
    if r.readline(): sys.exit("FATAL reads.tsv has more rows than 2 per valid pair")
with open(a.out, "w") as o:
    o.write("min_mapq\tpairs\treads\tin_phr_of_homologous_end\tin_phr_of_another_arm\trate_another_arm\n")
    for t in a.thresholds: o.write(f"{t}\t{pairs[t]}\t{tot[t]}\t{homo[t]}\t{mis[t]}\t{mis[t] / max(1, tot[t]):.6g}\n")
print(f"{n} valid pairs read\n" + open(a.out).read(), end="")
