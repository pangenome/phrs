#!/usr/bin/env python3
"""Plot KaryoScope annotate tracks (one row per sequence, one panel per feature set) for short sequences.

usage: plot_ks_tracks.py OUT_DIR COLORS_TSV FEATURE_SET[,FEATURE_SET...] ORDER(comma list of sequence names) OUT.png
Reads OUT_DIR/*.<feature_set>.smoothed.bed (seq, start, end, feature in columns 1-4); colours from the database colors.tsv
(feature_set, feature, color). x = bp from the start of each sequence.
"""
import collections, glob, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import Patch
od, cf, fss, order, out = sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4].split(","), sys.argv[5]
col = {(a, b): c for a, b, c in (l.rstrip("\n").split("\t") for l in list(open(cf))[1:])}
fig, axs = plt.subplots(len(fss), 1, figsize=(13, 1.0 + 0.42 * len(order) * len(fss)), squeeze=False)
for ax, fs in zip(axs[:, 0], fss):
    f = glob.glob(f"{od}/*.{fs}.smoothed.bed")
    if len(f) != 1: sys.exit(f"ERROR {len(f)} BEDs for {fs} in {od}")
    seen, L = {}, collections.Counter()
    for l in open(f[0]):
        if l.startswith(("#", "track")): continue
        t = l.rstrip("\n").split("\t"); s, a, b, ft = t[0], int(t[1]), int(t[2]), t[3]
        if s not in order: continue
        y = len(order) - 1 - order.index(s); c = col.get((fs, ft), "#DDDDDD")
        ax.barh(y, b - a, left=a, height=0.7, color=c, linewidth=0); seen[ft] = c; L[s] = max(L[s], b)
    missing = [s for s in order if s not in L]
    if missing: sys.exit(f"ERROR {fs}: no annotation for {missing}")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1], fontsize=8); ax.set_title(f"KaryoScope {fs}", fontsize=10, loc="left")
    ax.set_xlabel("bp from the start of each sequence", fontsize=8); ax.set_xlim(0, max(L.values()) * 1.01)
    ax.legend(handles=[Patch(color=c, label=k) for k, c in sorted(seen.items())], fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
plt.tight_layout(); plt.savefig(out, dpi=150); print(f"wrote {out}", file=sys.stderr)
