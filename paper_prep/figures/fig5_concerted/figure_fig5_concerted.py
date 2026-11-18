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
import json, csv, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HET = "/moosefs/erikg/phrs_tmp/hprc_het_final.json"
CUR = "/moosefs/erikg/phrs_tmp/telomere_end_curves.json"
K15 = ("/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
       "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
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

fig = plt.figure(figsize=(13.0, 7.4))
gs = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.26)

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
ax.set_yticklabels([f"{r[0]}  ({r[2]} pair{'s' if r[2]>1 else ''})" for r in rows], fontsize=7)
ax.set_xscale("log"); ax.set_xlim(0.85, 60)
ax.set_xlabel("divergence among member arms / allelic divergence")
ax.set_title("a  three exchange regimes", fontsize=9, loc="left")
ax.text(1.03, 0.15, "obligate: acrocentric short arms, PAR1/PAR2", fontsize=6, color="#c0392b")
ax.text(2.0, 3.0, "facultative exchange", fontsize=6, color="#555555")
_c4 = [i for i, r in enumerate(rows) if r[0] == "C4"]
if _c4:
    ax.text(rows[_c4[0]][1] * 0.42, _c4[0] - 0.42, "residual (forced pairing)", fontsize=6, color="#4d4d4d")

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
                      ("paralog_inter_chrom", "paralogous (different chromosome)", "#D95F02")):
    xs, ys = ecdf(pts(cls[key]["hist"]))
    ax.plot(xs, ys, color=col, lw=1.8,
            label=f"{lab} (med {cls[key]['d_sub_median_record']*100:.2f}%)")
ax.set_xscale("log"); ax.set_xlim(1e-4, 0.2)
ax.set_xlabel("substitution divergence"); ax.set_ylabel("cumulative fraction of aligned bp")
ax.set_title("b  homologues versus non-homologues", fontsize=9, loc="left")
ax.legend(fontsize=6.5, loc="lower right")

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
ax.axvspan(0, MASK, color="#d9d9d9", hatch="///", edgecolor="#999", alpha=0.6)
ax.text(MASK * 1.1, 40, "TERRA promoter array\nexcluded (ratio undefined)", fontsize=6)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.6, 60)
ax.set_xlabel("distance from the telomere (kb)")
ax.set_ylabel("paralogous / allelic ratio")
ax.set_title("c  homogenization strongest nearest the telomere", fontsize=9, loc="left")
ax.legend(fontsize=6.5, loc="lower right")

# ---- d: candidate competition at 9q, inside vs outside the tract
ax = fig.add_subplot(gs[1, 1])
cand = list(csv.DictReader(open("candidate_block_identity_9q.tsv"), delimiter="\t"))
cand = cand[::-1]
names = [c["candidate"] for c in cand]
tin = [float(c["tract_block_identity"]) for c in cand]
tout = [float(c["outside_block_identity"]) for c in cand]
y = range(len(cand)); h = 0.38
ax.barh([i + h/2 for i in y], tin, height=h, color="#1b6ca8", label="inside the tract", edgecolor="white")
ax.barh([i - h/2 for i in y], tout, height=h, color="#9ecae1", label="outside the tract", edgecolor="white")
ax.set_yticks(list(y)); ax.set_yticklabels(names, fontsize=7)
ax.set_xlim(0.985, 1.0008)
for i, (a, b) in enumerate(zip(tin, tout)):
    ax.text(a + 0.00015, i + h/2, f"{a:.4f}", va="center", fontsize=6)
    ax.text(b + 0.00015, i - h/2, f"{b:.4f}", va="center", fontsize=6)
ax.set_xlabel("IMPG block identity (gap-compressed) to the daughter's paternal chr9q")
ax.set_title("d  a non-homologous end inside the spread of the homologues", fontsize=9, loc="left")
ax.legend(fontsize=6, loc="lower right")

fig.suptitle("Concerted evolution of human subtelomeres: identical sequence at non-homologous ends, "
             "maintained without detectable exchange within a family", fontsize=9.5)
fig.savefig("figure_fig5_concerted.pdf", dpi=200, bbox_inches="tight")
fig.savefig("figure_fig5_concerted.png", dpi=200, bbox_inches="tight")
print("wrote figure_fig5_concerted.pdf/.png")
