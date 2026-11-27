#!/usr/bin/env python3
"""What the HPRC release-2 assembly QC says about the telomere-anchored flanks: HMM-Flagger v1.2.0 with HiFi and with ONT
(conservative, canonical, no_Hap BEDs: Err, Dup, Col, NNN) and NucFlag v0.3.3 with HiFi (MISJOIN, COLLAPSE, COLLAPSE_VAR,
COLLAPSE_OTHER, HET), one BED per haplotype (<dir>/<sample>_<hap>.bed). Per flank (names contig:a-b_chrN_arm, 1-based) we
sum the bp of each label in the flank and in its paper PHR (all-vs-all len TSV, region_start-region_end in flank
coordinates), and for the candidate joins (join_read_check.py output) whether a call lies within --pad bp of the join.
Output: per-flank TSV, per-join TSV, and a summary on stdout."""
import argparse, collections, csv, glob, os, re

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--fai", required=True, help="flank names (first column)")
ap.add_argument("--phr", required=True)
ap.add_argument("--tool", nargs=2, action="append", metavar=("NAME", "DIR"), required=True)
ap.add_argument("--joins", required=True)
ap.add_argument("--verdicts", required=True)
ap.add_argument("--pad", type=int, default=5000)
ap.add_argument("--out", required=True, help="prefix")
a = ap.parse_args()
NAME = re.compile(r"^(([^#]+)#([^#]+)#[^:]+):(\d+)-(\d+)_(chr[0-9XY]+)_(parm|qarm)$")
fl = {}
for l in open(a.fai):
    nm = l.split("\t")[0]; m = NAME.match(nm)
    if not m:
        raise SystemExit(f"FATAL: cannot parse {nm}")
    fl[nm] = (m.group(1), m.group(2), m.group(3), int(m.group(4)) - 1, int(m.group(5)))
phr = {r["seq"]: (0, 0) if r["region_start"] == r["region_end"] == "." else (int(r["region_start"]), int(r["region_end"]))
       for r in csv.DictReader(open(a.phr), delimiter="\t")}  # "." = no PHR
ctgs = {v[0] for v in fl.values()}
print(f"{len(fl)} flanks, {sum(1 for k in fl if k in phr and phr[k][1] > phr[k][0])} with a PHR, {len(ctgs)} contigs")
calls, have = {}, {}
for name, d in a.tool:
    c = collections.defaultdict(list); have[name] = set()
    for f in sorted(glob.glob(os.path.join(d, "*.bed"))):
        smp, hap = os.path.basename(f)[:-4].rsplit("_", 1); have[name].add((smp, hap))
        for l in open(f):
            if l.startswith(("track", "#")):
                continue
            x = l.split("\t")
            if x[0] in ctgs:
                c[x[0]].append((int(x[1]), int(x[2]), x[3]))
    calls[name] = c
    print(f"{name}: {len(have[name])} haplotype BEDs, calls on flank contigs {sum(len(v) for v in c.values())}")
labels = {n: sorted({k for v in calls[n].values() for *_, k in v}) for n in calls}
ov = lambda s, e, s0, e0: max(0, min(e, e0) - max(s, s0))
rows = []
for nm, (ctg, smp, hap, s0, e0) in fl.items():
    rs, re_ = phr.get(nm, (0, 0)); ps, pe = s0 + rs, s0 + re_
    r = {"flank": nm, "phr_bp": re_ - rs}
    for n in calls:
        r[n + ":covered"] = int((smp, hap) in have[n])
        for k in labels[n]:
            r[f"{n}:{k}"] = r[f"{n}:{k}:in_phr"] = 0
        for s, e, k in calls[n].get(ctg, []):
            r[f"{n}:{k}"] += ov(s, e, s0, e0); r[f"{n}:{k}:in_phr"] += ov(s, e, ps, pe)
    rows.append(r)
cols = list(rows[0])
with open(f"{a.out}.flanks.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
verdict = {r["flank"]: r["verdict"] for r in csv.DictReader(open(a.verdicts), delimiter="\t")}
for n in calls:
    cov = [r for r in rows if r[n + ":covered"]]
    print(f"\n{n}: {len(cov)} of {len(rows)} flanks have a QC BED (missing: {sorted({fl[r['flank']][1] for r in rows if not r[n + ':covered']})})")
    for k in labels[n]:
        b = [r[f"{n}:{k}"] for r in cov]; bp = [r[f"{n}:{k}:in_phr"] for r in cov if r["phr_bp"] > 0]
        print(f"  {k:15s} flanks with >= 1 bp {sum(x > 0 for x in b):6d}, >= 10 kb {sum(x >= 10000 for x in b):6d}, share of flank bp "
              f"{sum(b) / (500000 * len(cov)):.4f} | in PHR: flanks {sum(x > 0 for x in bp):6d} of {len(bp)}, share of PHR bp "
              f"{sum(bp) / max(1, sum(r['phr_bp'] for r in cov if r['phr_bp'] > 0)):.4f}")
    for vd in ("correct", "wrong", "wrong_tiny", "undecided"):
        v = [r for r in cov if verdict.get(r["flank"]) == vd]
        err = [k for k in labels[n] if k not in ("NNN", "HET")]
        print(f"  joins-table {vd:10s}: {len(v)} flanks, with any {'/'.join(err)} call: {sum(any(r[f'{n}:{k}'] > 0 for k in err) for r in v)}")
jr = []
for r in csv.DictReader(open(a.joins), delimiter="\t"):
    if r["tech"] != "ONT":
        continue
    ctg, s, e = r["contig"], int(r["start"]) - 1 - a.pad, int(r["end"]) + a.pad
    x = {"flank": r["flank"], "kind": r["kind"], "dist_from_telomere": r["dist_from_telomere"], "verdict": verdict[r["flank"]]}
    for n in calls:
        hits = collections.Counter(k for s1, e1, k in calls[n].get(ctg, []) if ov(s1, e1, s, e) > 0)
        x[n] = ",".join(f"{k}:{c}" for k, c in sorted(hits.items())) or "."
    jr.append(x)
with open(f"{a.out}.joins.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(jr[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jr)
print(f"\njoins (+-{a.pad} bp), by flank verdict: calls other than NNN/HET at the join")
for vd in ("correct", "wrong", "wrong_tiny", "undecided"):
    v = [x for x in jr if x["verdict"] == vd]
    print(f"  {vd:10s} {len(v):3d} joins | " + " | ".join(
        f"{n}: {sum(any(k.split(':')[0] not in ('NNN', 'HET') for k in x[n].split(',') if k != '.') for x in v)}" for n in calls))
