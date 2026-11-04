#!/usr/bin/env python3
"""ED Fig 6 - paralog/allelic divergence as a function of distance from the telomere.

Panel A: per-arm curves from histogram medians (robust to heavy tails), plus the cohort median.
Panel B: one curve per individual chromosome end (weighted ratios; texture, not magnitude), with
         community medians overlaid.
Panel C: per-arm median ratio grouped by community tier.

The first ~6 kb is masked throughout: allelic divergence there is 0.3-2.5% versus 0.09-0.16% in the
interior, which cannot be orthologous divergence, so the denominator is invalid and those values are
not plotted.

Inputs (see sources.tsv):
  /moosefs/erikg/phrs_tmp/telomere_end_curves.json   (scripts/ci/telomere_end_ratio_curves.py)
  arm-level Leipzig k=15 assignments (see sources.tsv)
Run from this directory: python3 figure_ed6_telomere_profile.py
"""
import json, random, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

J = "/moosefs/erikg/phrs_tmp/telomere_end_curves.json"
K15 = ("/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
       "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
MASK = 6.0
TIERS = {"C7": "#c0392b", "C15": "#8e44ad", "C14": "#9b59b6",
         "C3": "#D95F02", "C11": "#1b6ca8", "C1": "#2a9d8f"}
LAB = {"C7": "C7 acrocentric", "C15": "C15 PAR1", "C14": "C14 PAR2",
       "C3": "C3 (3q/9q/11p...)", "C11": "C11 (1p/5q/6q/8p)", "C1": "C1 (4q/10q)"}
d = json.load(open(J)); comm = {l.split("\t")[0]: l.split("\t")[1] for l in open(K15).read().splitlines()[1:]}
arms, ends = d["arms"], d["ends"]

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.2, 4.0), gridspec_kw={"width_ratios": [1.5, 1.5, 1]})

# ---- A: per-arm, histogram medians
xall, yall = [], []
for arm, rows in arms.items():
    r = [x for x in rows if x.get("ratio") and x["from_kb"] >= MASK]
    if len(r) < 5: continue
    xs = [(x["from_kb"] + x["to_kb"]) / 2 for x in r]; ys = [x["ratio"] for x in r]
    axA.plot(xs, ys, color="#9a9a9a", lw=0.6, alpha=0.45, zorder=2)
# cohort median per bin over arm medians
bins = sorted({(x["from_kb"], x["to_kb"]) for rows in arms.values() for x in rows if x["from_kb"] >= MASK})
mx, my = [], []
for lo, hi in bins:
    vals = [x["ratio"] for rows in arms.values() for x in rows
            if x.get("ratio") and x["from_kb"] == lo]
    if len(vals) >= 10:
        mx.append((lo + hi) / 2); my.append(st.median(vals))
axA.plot(mx, my, color="black", lw=2.4, zorder=6, label="cohort median (48 arms)")
axA.axvspan(0, MASK, color="#d9d9d9", hatch="///", edgecolor="#999999", alpha=0.55, zorder=0)
axA.set_xscale("log"); axA.set_yscale("log"); axA.set_ylim(0.5, 60)
axA.set_xlabel("distance from telomere (kb)"); axA.set_ylabel("paralog / allelic ratio")
axA.set_title("a  per chromosome arm (medians)", fontsize=9, loc="left"); axA.legend(fontsize=6.5, loc="lower right")

# ---- B: per-end spaghetti + community medians
usable = {k: [r for r in v if r.get("ratio") and r["from_kb"] >= MASK] for k, v in ends.items()}
usable = {k: v for k, v in usable.items() if len(v) >= 6}
random.seed(3)
for k in random.sample(sorted(usable), min(350, len(usable))):
    axB.plot([(r["from_kb"] + r["to_kb"]) / 2 for r in usable[k]], [r["ratio"] for r in usable[k]],
             color="#c0c0c0", lw=0.4, alpha=0.3, zorder=1)
for c, col in TIERS.items():
    sub = [v for k, v in usable.items() if comm.get(k.split("|")[1]) == c]
    if len(sub) < 3: continue
    bx, by = [], []
    for lo, hi in bins:
        vals = [r["ratio"] for v in sub for r in v if r["from_kb"] == lo]
        if len(vals) >= 20:
            bx.append((lo + hi) / 2); by.append(st.median(vals))
    axB.plot(bx, by, color=col, lw=1.7, zorder=5, label=LAB[c])
axB.axvspan(0, MASK, color="#d9d9d9", hatch="///", edgecolor="#999999", alpha=0.55, zorder=0)
axB.set_xscale("log"); axB.set_yscale("log"); axB.set_ylim(0.5, 60)
axB.set_xlabel("distance from telomere (kb)"); axB.set_ylabel("paralog / allelic ratio")
axB.set_title("b  per chromosome end (%d shown)" % len(usable), fontsize=9, loc="left")
axB.legend(fontsize=6, loc="lower right", ncol=1)

# ---- C: per-arm median ratio by tier
data, labs, cols = [], [], []
for c in ["C7", "C15", "C14", "C1", "C11", "C3"]:
    vals = []
    for arm, rows in arms.items():
        if comm.get(arm) != c: continue
        v = [x["ratio"] for x in rows if x.get("ratio") and x["from_kb"] >= 10]
        if len(v) >= 5: vals.append(st.median(v))
    if vals: data.append(vals); labs.append(c); cols.append(TIERS[c])
other = []
for arm, rows in arms.items():
    if comm.get(arm) in TIERS: continue
    v = [x["ratio"] for x in rows if x.get("ratio") and x["from_kb"] >= 10]
    if len(v) >= 5: other.append(st.median(v))
if other: data.append(other); labs.append("other"); cols.append("#777777")
bp = axC.boxplot(data, patch_artist=True, widths=0.62, showfliers=False)
for p, c in zip(bp["boxes"], cols): p.set_facecolor(c); p.set_alpha(0.55)
for m in bp["medians"]: m.set_color("black"); m.set_linewidth(1.3)
axC.set_xticklabels(labs, fontsize=7); axC.axhline(1, color="k", ls=":", lw=1)
axC.set_yscale("log"); axC.set_ylabel("per-arm median ratio (>=10 kb)")
axC.set_title("c  tiers, arm by arm", fontsize=9, loc="left")
plt.tight_layout()
fig.savefig("figure_ed6_telomere_profile.pdf", dpi=190, bbox_inches="tight")
fig.savefig("figure_ed6_telomere_profile.png", dpi=190, bbox_inches="tight")
print("wrote figure_ed6_telomere_profile.pdf/.png  arms=%d ends=%d" % (len(arms), len(usable)))
