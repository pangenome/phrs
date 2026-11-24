#!/usr/bin/env python3
"""Summarize geneconv_families.py output by copy (R1-M1, test 3a).

Per alignment family and level (units or arrays, over its haplotype batches): copies kept, cross-arm GI
fragments, those with simulation P < 0.05 / number of batches (Bonferroni over batches), and the copies
that recur in >= --min-frags significant cross-arm fragments of their batch. One converted copy appears
in one fragment per partner copy, so events are counted by copy, not by pair.
"""
import argparse, collections, csv, re, sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("outdirs", nargs="+", help="geneconv_families.py --out directories")
ap.add_argument("--min-frags", type=int, default=3)
a = ap.parse_args()
summ, frags = [], []
for d in a.outdirs:
    summ += list(csv.DictReader(open(f"{d}/summary.tsv"), delimiter="\t"))
    frags += list(csv.DictReader(open(f"{d}/fragments.tsv"), delimiter="\t"))
groups = collections.defaultdict(list)
for s in summ:
    groups[re.sub(r"_b\d+$", "", s["alignment"])].append(s)
print("alignment\tbatches\tcopies_kept\tfar_dropped\tcross_arm_GI\tcross_arm_GI_bonferroni\tcopies_recurrent\tcopies")
for g, S in sorted(groups.items()):
    nb = len(S)
    names = {s["alignment"] for s in S}
    F = [f for f in frags if f["alignment"] in names and f["pair_class"] == "cross_arm"]
    Fs = [f for f in F if float(f["sim_p"]) < 0.05 / nb]
    per = collections.defaultdict(collections.Counter)
    for f in Fs:
        for x in f["names"].split(";"):
            per[f["alignment"]][x] += 1
    rec = sorted(((n, x) for c in per.values() for x, n in c.items() if n >= a.min_frags), reverse=True)
    print(f"{g}\t{nb}\t{sum(int(s['n_kept']) for s in S if s['status'] == 'ok')}\t"
          f"{sum(int(s.get('n_far_dropped') or 0) for s in S if s['status'] == 'ok')}\t{len(F)}\t{len(Fs)}\t{len(rec)}\t"
          + ",".join(f"{x}:{n}" for n, x in rec))
    if any(s["status"] != "ok" for s in S):
        sys.stderr.write(f"{g}: {sum(s['status'] != 'ok' for s in S)} batches not run ({S[0]['status']})\n")
