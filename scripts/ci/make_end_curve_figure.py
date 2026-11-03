#!/usr/bin/env python3
"""Spaghetti figure: paralog/allelic ratio vs distance from the telomere, one curve per chromosome end.

Left panel: every individual end as a faint curve (subsampled for legibility), with the cohort median
and per-community medians overlaid. The first ~6 kb is masked: the allelic denominator is invalid
there (allelic divergence 0.3-2.5% vs 0.09-0.16% in the interior), so those values are not
measurements and are not drawn.

Right panel: distribution of each end's own median ratio (bins >= 10 kb), coloured by community tier.
"""
import argparse, json, random, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

K15 = ("/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
       "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
MASK_KB = 6.0
TIERS = {"C7": "#c0392b", "C15": "#8e44ad", "C14": "#9b59b6",
         "C3": "#D95F02", "C11": "#1b6ca8", "C1": "#2a9d8f"}
TIER_LABEL = {"C7": "C7 acrocentric (ratio ~1)", "C15": "C15 PAR1", "C14": "C14 PAR2",
              "C3": "C3 (3q/9q/11p...)", "C11": "C11 (1p/5q/6q/8p)", "C1": "C1 (4q/10q)"}

ap = argparse.ArgumentParser(); ap.add_argument("--json", required=True)
ap.add_argument("--out", default="response/notes/telomere_end_curves"); a = ap.parse_args()
d = json.load(open(a.json))
comm = {l.split("\t")[0]: l.split("\t")[1] for l in open(K15).read().splitlines()[1:]}

def usable(rows):
    return [r for r in rows if r.get("ratio") and r["from_kb"] >= MASK_KB]

ends, arms = d["ends"], d["arms"]
end_curves = {k: usable(v) for k, v in ends.items()}
end_curves = {k: v for k, v in end_curves.items() if len(v) >= 6}
random.seed(3)
sample = random.sample(sorted(end_curves), min(400, len(end_curves)))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [2.1, 1]})

for k in sample:
    xs = [(r["from_kb"] + r["to_kb"]) / 2 for r in end_curves[k]]
    ys = [r["ratio"] for r in end_curves[k]]
    ax.plot(xs, ys, color="#b0b0b0", alpha=0.18, lw=0.45, zorder=1)

bins = sorted({(r["from_kb"], r["to_kb"]) for v in end_curves.values() for r in v})
def median_by_bin(subset):
    out = []
    for lo, hi in bins:
        vals = [r["ratio"] for v in subset for r in v if r["from_kb"] == lo and r["to_kb"] == hi]
        if len(vals) >= 20:
            out.append(((lo + hi) / 2, st.median(vals)))
    return zip(*out) if out else ((), ())

xs, ys = median_by_bin(end_curves.values())
ax.plot(xs, ys, color="black", lw=2.4, zorder=6, label="cohort median (n=%d ends)" % len(end_curves))
for c, col in TIERS.items():
    sub = {k: v for k, v in end_curves.items() if comm.get(k.split("|")[1]) == c}
    if len(sub) < 3:
        continue
    xs, ys = median_by_bin(sub.values())
    ax.plot(xs, ys, color=col, lw=1.6, zorder=5, label=TIER_LABEL[c] + " (n=%d)" % len(sub))
ax.axvspan(0, MASK_KB, color="#d9d9d9", hatch="///", edgecolor="#999999", alpha=0.55, zorder=0)
ax.text(MASK_KB * 1.08, ax.get_ylim()[1] if False else 60,
        "first ~6 kb: allelic denominator invalid\n(allelic div. 0.3-2.5% vs 0.09-0.16% interior)\nnot plotted",
        fontsize=6.5, va="top")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("distance from telomere-subtelomere junction (kb)")
ax.set_ylabel("paralog / allelic divergence ratio")
ax.set_title("A  one curve per chromosome end", fontsize=9, loc="left")
ax.legend(fontsize=5.6, loc="lower right", ncol=1, framealpha=0.9)
ax.set_ylim(0.5, 80)

# right: per-end median ratio, coloured by tier
data, labels, colors = [], [], []
for c in ["C7", "C15", "C14", "C1", "C11", "C3"]:
    vals = [st.median([r["ratio"] for r in v]) for k, v in end_curves.items()
            if comm.get(k.split("|")[1]) == c]
    if vals:
        data.append(vals); labels.append(c); colors.append(TIERS[c])
vals = [st.median([r["ratio"] for r in v]) for k, v in end_curves.items()
        if comm.get(k.split("|")[1]) not in TIERS]
if vals:
    data.append(vals); labels.append("other"); colors.append("#888888")
parts = ax2.boxplot(data, vert=True, patch_artist=True, widths=0.6, showfliers=False)
for p, c in zip(parts["boxes"], colors):
    p.set_facecolor(c); p.set_alpha(0.55)
for med in parts["medians"]:
    med.set_color("black"); med.set_linewidth(1.4)
ax2.set_xticklabels(labels, fontsize=7)
ax2.axhline(1, color="k", ls=":", lw=1)
ax2.set_yscale("log"); ax2.set_ylabel("per-end median ratio (bins >= 10 kb)")
ax2.set_title("B  tier structure, end by end", fontsize=9, loc="left")
plt.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(a.out + "." + ext, dpi=190, bbox_inches="tight")
print("wrote", a.out)
print("ends drawn:", len(sample), "of", len(end_curves))
