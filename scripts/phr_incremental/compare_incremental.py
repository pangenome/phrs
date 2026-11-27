#!/usr/bin/env python3
"""Compare incremental PHR calls (phr_incremental.py TSVs) with a reference all-vs-all run: region length per flank,
pass/fail per window (windows the reference scanned), rounds (haplotypes) the incremental run needed, and its work
(bp aligned x partner flanks) as a share of the full scan."""
import argparse, csv, glob, os, statistics as st

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--dir", required=True)
ap.add_argument("--calls", required=True)
ap.add_argument("--windows", required=True)
a = ap.parse_args()
new, summ = {}, {}
for f in sorted(glob.glob(os.path.join(a.dir, "*.tsv"))):
    for l in open(f):
        if l.startswith("#SUMMARY"):
            x = l.rstrip("\n").split("\t")
            summ[x[1]] = dict(kv.split("=") for kv in x[2:])
        elif not l.startswith("seq\t"):
            x = l.rstrip("\n").split("\t")  # 7 columns: haplotypes_checked; 6 (mode all, first version): every window checked
            if len(x) == 7 and x[5] == "0":
                continue  # window never scanned (beyond the stop)
            new.setdefault(x[0], {})[int(x[2])] = (x[3] == "True", x[4])
ref = {r["seq"]: r for r in csv.DictReader(open(a.calls), delimiter="\t") if r["seq"] in new}
rw = {}
for r in csv.DictReader(open(a.windows), delimiter="\t"):
    if r["seq"] in new:
        rw.setdefault(r["seq"], {})[int(r["dist_from_telomere"])] = r["passes_threshold"] == "True"
L = lambda r: 0 if r["region_start"] == "." else int(r["region_end"]) - int(r["region_start"])
same = 0
print("flank\tref_len\tnew_len\twindows_compared\tref_pass_new_fail\tref_fail_new_pass\tlast_round_needed\twork_share")
for s in sorted(new):
    if s not in summ:
        print(f"{s}\tNA\tNA\trun not finished"); continue
    nl = int(summ[s]["region_length"]); rl = L(ref[s]); same += nl == rl
    comp = [d for d in rw.get(s, {}) if d in new[s]]
    a1 = sum(rw[s][d] and not new[s][d][0] for d in comp); a2 = sum(new[s][d][0] and not rw[s][d] for d in comp)
    rounds = [int(v[1]) for v in new[s].values() if v[1] != "NA"]
    print(f"{s}\t{rl}\t{nl}\t{len(comp)}\t{a1}\t{a2}\t{max(rounds) if rounds else 'NA'}\t{int(summ[s]['work_bp']) / int(summ[s]['full_bp']):.3f}")
done = [s for s in new if s in summ]
print(f"SUMMARY: {same} of {len(done)} finished flanks have the reference region length; work share median "
      f"{st.median(int(summ[s]['work_bp']) / int(summ[s]['full_bp']) for s in done):.3f}")
