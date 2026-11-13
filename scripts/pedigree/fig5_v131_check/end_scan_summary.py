#!/usr/bin/env python3
"""Merge end_scan bins with >= 5 wrong positions into regions, per sample and version.
usage: end_scan_summary.py END_SCAN_DIR > regions.tsv
Region = adjacent 10 kb bins of one contig end with >= 5 wrong positions each. depth_ratio = region median MAPQ>=20 depth
over the median of all bins of that sample and version; low (< 0.5) suggests a false copy, high (> 2) a collapsed repeat."""
import csv, glob, os, statistics as st, sys

print("sample\tversion\tcontig\tend\tstart\tend_bp\tbins\twrong\tdepth_q20_median\tdepth_ratio\tflag")
for f in sorted(glob.glob(os.path.join(sys.argv[1], "end_scan.*.tsv"))):
    if f.endswith(".wrong.tsv"):
        continue
    _, sm, v, _ = os.path.basename(f).split(".")
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    med = st.median(float(r["depth_mapq20_median"]) for r in rows) or 1.0
    reg, cur = [], None
    for r in rows:
        w, s = int(r["wrong_positions"]), int(r["bin_start"])
        if w >= 5 and cur and cur["contig"] == r["contig"] and s == cur["last"] + 10000:
            cur["last"], cur["wrong"] = s, cur["wrong"] + w
            cur["d"].append(float(r["depth_mapq20_median"]))
        elif w >= 5:
            cur = {"contig": r["contig"], "end": r["end"], "first": s, "last": s, "wrong": w, "d": [float(r["depth_mapq20_median"])]}
            reg.append(cur)
        else:
            cur = None
    print(f"{sm} {v}: {len(rows)} bins, genome median depth {med}, regions {len(reg)}", file=sys.stderr)
    for g in reg:
        d = st.median(g["d"]); ratio = d / med
        flag = "low_depth" if ratio < 0.5 else ("high_depth" if ratio > 2 else "normal_depth")
        print(f"{sm}\t{v}\t{g['contig']}\t{g['end']}\t{g['first']}\t{g['last'] + 9999}\t{len(g['d'])}\t{g['wrong']}\t{d}\t{ratio:.2f}\t{flag}")
