#!/usr/bin/env python3
"""Figure 5 - concerted evolution of human subtelomeres.

a  Three-tier structure: for each Leiden community, divergence among its member arms divided by their
   allelic divergence. 1.0 = non-homologous copies as similar as alleles.
b  Allelic versus paralogous divergence, cumulative over aligned base pairs.
c  Paralogous/allelic ratio against distance from the telomere, one vote per chromosome arm; the first
   ~6 kb is excluded (TERRA promoter array, ratio undefined).
d  A non-homologous end inside the spread of the homologous ones: IMPG similarity of every candidate end
   to the daughter's paternal chr9q across the called tract, with all candidates in one graph and the
   father's chr9q rebuilt from his own reads.

Inputs are listed in sources.tsv. Run from this directory: python3 figure_fig5_concerted.py
"""
import json, csv, os, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Sizes for a 13 in wide figure printed at about 7 in: 13 pt here prints at about 7 pt
plt.rcParams.update({"font.size": 14, "axes.labelsize": 15, "xtick.labelsize": 13, "ytick.labelsize": 13,
                     "legend.fontsize": 12})
FS_TICK, FS_TITLE, FS_NOTE, FS_VAL = 13, 16, 12, 11

# FIG5_INPUT_DIR: a local copy of the three inputs below (same file names), when /moosefs is not mounted
IN = os.environ.get("FIG5_INPUT_DIR")
HET = os.path.join(IN, "hprc_het_final.json") if IN else "/moosefs/erikg/phrs_tmp/hprc_het_final.json"
CUR = os.path.join(IN, "telomere_end_curves.json") if IN else "/moosefs/erikg/phrs_tmp/telomere_end_curves.json"
K15 = (os.path.join(IN, "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv") if IN else
       "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
print(f"inputs: {HET} {CUR} {K15}")
TIER = {"C7": "#c0392b", "C15": "#8e44ad", "C14": "#9b59b6", "C1": "#2a9d8f",
        "C11": "#1b6ca8", "C3": "#D95F02", "C2": "#7f7f7f", "C5": "#7f7f7f",
        "C6": "#7f7f7f", "C12": "#7f7f7f", "C4": "#4d4d4d"}

d = json.load(open(HET)); cur = json.load(open(CUR))
comm = {l.split("\t")[0]: l.split("\t")[1] for l in open(K15).read().splitlines()[1:]}
allp = d["pairs"]
for r in allp:
    r["d_median_record"] = float(r["d_median_record"])
par = [r for r in allp if r["arm_a"] != r["arm_b"]]
ALE = {r["arm_a"]: r["d_median_record"] for r in allp if r["arm_a"] == r["arm_b"]}
ch = lambda a: a.rsplit("_", 1)[0]

fig = plt.figure(figsize=(13.0, 9.0))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.42)

# ---- a: tiers
ax = fig.add_subplot(gs[0, 0])
rows = []
for c in sorted(set(comm.values()), key=lambda z: int(z[1:])):
    arms = sorted(a for a, v in comm.items() if v == c)
    ia = [r["d_median_record"] for r in par
          if comm.get(r["arm_a"]) == c == comm.get(r["arm_b"]) and ch(r["arm_a"]) != ch(r["arm_b"])]
    if len(ia) < 1: continue
    al = [ALE[a] for a in arms if a in ALE]
    if not al: continue
    rows.append((c, st.median(ia) / st.median(al), len(ia), ",".join(arms)))
rows.sort(key=lambda r: r[1])
ys = range(len(rows))
ax.hlines(list(ys), 1.0, [r[1] for r in rows], color="#aaaaaa", lw=1)
ax.scatter([r[1] for r in rows], list(ys), s=70,
           c=[TIER.get(r[0], "#7f7f7f") for r in rows], zorder=3)
ax.axvline(1.0, color="k", lw=1.2)
ax.set_yticks(list(ys))
ax.set_yticklabels([f"{r[0]}  ({r[2]} pair{'s' if r[2]>1 else ''})" for r in rows], fontsize=FS_TICK)
ax.set_xscale("log"); ax.set_xlim(0.85, 120); ax.set_ylim(-0.6, len(rows) - 0.2)
ax.set_xlabel("divergence among member arms / allelic divergence")
ax.set_title("a  three exchange regimes", fontsize=FS_TITLE, loc="left")
ax.text(1.12, 0.5, "complete mixing: acrocentric short arms, PAR1", fontsize=FS_NOTE, color="#c0392b", va="center")
ax.text(5.0, 4.0, "cohesive", fontsize=FS_NOTE, color="#1b6ca8", va="center")
ax.text(16.0, 7.0, "facultative", fontsize=FS_NOTE, color="#555555", va="center")
_c4 = [i for i, r in enumerate(rows) if r[0] == "C4"]
if _c4:
    ax.text(rows[_c4[0]][1] * 1.15, _c4[0], "residual\n(forced pairing)", fontsize=FS_NOTE, color="#4d4d4d", va="center")

# ---- b: ECDF allelic vs paralogous
ax = fig.add_subplot(gs[0, 1])
cls = d["summary"]["classes"]
bpu = int(cls["allelic_between_indiv"].get("hist_bins_per_unit", 10000))
def pts(h): return [(i / bpu, w) for i, w in enumerate(h) if w]
def ecdf(p):
    tot = sum(w for _, w in p); xs, ys, run = [], [], 0.0
    for x, w in sorted(p):
        run += w; xs.append(x); ys.append(run / tot)
    return xs, ys
for key, lab, col in (("allelic_between_indiv", "allelic", "#4c72b0"),
                      ("paralog_inter_chrom", "paralogous", "#D95F02")):
    xs, ys = ecdf(pts(cls[key]["hist"]))
    ax.plot(xs, ys, color=col, lw=1.8,
            label=f"{lab}, median {cls[key]['d_sub_median_record']*100:.2f}%")
ax.set_xscale("log"); ax.set_xlim(1e-4, 0.2)
ax.set_xlabel("substitution divergence"); ax.set_ylabel("cumulative fraction of aligned bp")
ax.set_title("b  homologues versus non-homologues", fontsize=FS_TITLE, loc="left")
ax.legend(loc="lower right")

# ---- c: ratio vs distance
ax = fig.add_subplot(gs[1, 0])
MASK = 6.0
arms = cur["arms"]
bins = sorted({(x["from_kb"], x["to_kb"]) for rr in arms.values() for x in rr if x["from_kb"] >= MASK})
mx, my = [], []
for lo, hi in bins:
    vals = [x["ratio"] for rr in arms.values() for x in rr if x.get("ratio") and x["from_kb"] == lo]
    if len(vals) >= 10:
        mx.append((lo + hi) / 2); my.append(st.median(vals))
for rr in arms.values():
    xs = [(x["from_kb"] + x["to_kb"]) / 2 for x in rr if x.get("ratio") and x["from_kb"] >= MASK]
    ys = [x["ratio"] for x in rr if x.get("ratio") and x["from_kb"] >= MASK]
    if len(xs) >= 5: ax.plot(xs, ys, color="#9a9a9a", lw=0.5, alpha=0.4)
ax.plot(mx, my, color="k", lw=2.4, label="median, 48 arms")
ax.axvspan(1, MASK, facecolor="#eeeeee", hatch="///", edgecolor="#888888", lw=0)
ax.text(1.15, 0.68, "TERRA\npromoter array,\nexcluded", fontsize=FS_NOTE - 1, va="bottom",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.5))
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.6, 60); ax.set_xlim(1, 600)
ax.set_xlabel("distance from the telomere (kb)")
ax.set_ylabel("paralogous / allelic ratio")
ax.set_title("c  homogenization strongest nearest the telomere", fontsize=FS_TITLE, loc="left")
ax.legend(loc="lower right")

# ---- d: candidate competition at 9q, inside vs outside the tract
ax = fig.add_subplot(gs[1, 1])
cand = list(csv.DictReader(open("candidate_block_identity_9q.tsv"), delimiter="\t"))
cand = cand[::-1]
names = [c["candidate"] for c in cand]
num = lambda x: None if x == "NA" else float(x)  # NA: no alignment in that region
tin = [num(c["tract_dice"]) for c in cand]
tout = [num(c["outside_dice"]) for c in cand]
y = range(len(cand)); h = 0.38; X0 = 0.975
for vals, off, col, lab in ((tin, h/2, "#1b6ca8", "inside the tract"), (tout, -h/2, "#9ecae1", "outside the tract")):
    ax.barh([i + off for i, v in zip(y, vals) if v is not None], [v for v in vals if v is not None], height=h, color=col,
            label=lab, edgecolor="white")
    for i, v in zip(y, vals):  # value at the bar end; below the axis start or missing: text at the axis start
        txt = "no alignment" if v is None else f"{v:.4f}" + (" (off scale)" if v < X0 else "")
        ax.text(X0 + 0.0002 if v is None or v < X0 else v + 0.00015, i + off, txt, va="center", fontsize=FS_VAL)
ax.set_yticks(list(y)); ax.set_yticklabels(names, fontsize=FS_TICK)
ax.set_xlim(X0, 1.0025); ax.set_xticks([0.976, 0.982, 0.988, 0.994, 1.0])
ax.set_xlabel("IMPG dice similarity to the daughter's paternal chr9q")
ax.set_title("d  his chr3q, not his other chr9q, is the next match", fontsize=FS_TITLE, loc="left")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)

fig.suptitle("Concerted evolution of human subtelomeres: identical sequence at non-homologous ends,\n"
             "maintained without detectable exchange within a family", fontsize=15)
fig.savefig("figure_fig5_concerted.pdf", dpi=200, bbox_inches="tight")
fig.savefig("figure_fig5_concerted.png", dpi=200, bbox_inches="tight")
print("wrote figure_fig5_concerted.pdf/.png")
