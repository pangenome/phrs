#!/usr/bin/env python3
"""Figures: the scope of subtelomeric divergence - allelic vs paralogous vs community.

Reads the cohort run of scripts/ci/hprc_heterolog_divergence.py (JSON with per-class
histograms) plus the k=15 community assignments, and draws four panels:

  A  ECDF of substitution divergence by class, log x: allelic / paralog in-community /
     paralog out-of-community. Medians marked; fold-change vs allelic annotated.
  B  Fold-change of each class median against the allelic baseline (the "many X over
     allelic" statement), with the ancient-duplication reference band (100-1000x).
  C  Per-community cohesion: median intra-community divergence divided by the same arms'
     out-of-community divergence, ordered, with pair counts.
  D  The two axes are separable: per arm pair, median divergence (magnitude) against the
     fraction of aligned bases at >= 99.5% identity (the near-identical tail), coloured by
     whether the pair is in-community. The tail separates the classes; the magnitude does not.

Usage:
  python3 scripts/ci/make_heterolog_figure.py --json /moosefs/erikg/phrs_tmp/hprc_het_hist.json \
      --out response/notes/heterolog_scope
"""
import argparse
import csv
import json
import math
import os
import statistics as st
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K15 = ("/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
       "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
C_ALLELE = "#4c72b0"
C_IN = "#D95F02"
C_OUT = "#757575"


def hist_to_points(hist, bins_per_unit=10000):
    """[(d, weight)] from a weighted histogram whose bin i covers d = i/bins_per_unit."""
    return [(i / bins_per_unit, w) for i, w in enumerate(hist) if w]


def ecdf(points):
    tot = sum(w for _, w in points)
    xs, ys, run = [], [], 0.0
    for d, w in sorted(points):
        run += w
        xs.append(d)
        ys.append(run / tot)
    return xs, ys


def median_of(points):
    tot = sum(w for _, w in points)
    run = 0.0
    for d, w in sorted(points):
        run += w
        if run >= 0.5 * tot:
            return d
    return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", default="response/notes/heterolog_scope")
    a = ap.parse_args()

    d = json.load(open(a.json))
    classes = d["summary"].get("classes", {})
    if not classes:
        sys.exit("no per-class histograms in this JSON - rerun the aggregator with hist dumping")
    pairs = [r for r in d["pairs"] if r["arm_a"] != r["arm_b"]]   # drop allelic self-pairs
    for r in pairs:
        for k in ("d_sub", "d_median_record", "frac_at_99.5pct", "aligned_bp"):
            try:
                r[k] = float(r[k])
            except (TypeError, ValueError):
                r[k] = float("nan")

    comm = {}
    with open(K15) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            comm[row["ChromArm"].strip()] = row["Community"].strip()

    alle = classes["allelic_between_indiv"]
    bpu = int(alle.get("hist_bins_per_unit", 10000))
    alle_pt = hist_to_points(alle["hist"], bpu) + hist_to_points(
        classes.get("allelic_within_indiv", alle)["hist"], bpu)
    intra_rows = [r for r in pairs if comm.get(r["arm_a"]) and comm.get(r["arm_a"]) == comm.get(r["arm_b"])]
    inter_rows = [r for r in pairs if comm.get(r["arm_a"]) and comm.get(r["arm_b"])
                  and comm[r["arm_a"]] != comm[r["arm_b"]]]
    # per-record medians are the fair comparison for pairs; rebuild a weighted point set
    # from each pair's per-record median weighted by its aligned bp (approximation, stated in caption)
    intra_pt = [(r["d_median_record"], r["aligned_bp"]) for r in intra_rows]
    inter_pt = [(r["d_median_record"], r["aligned_bp"]) for r in inter_rows]

    fig = plt.figure(figsize=(11.5, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)

    # ---- A: ECDF
    ax = fig.add_subplot(gs[0, 0])
    # record-level histograms only: exact, and consistent across curves
    inter_pt_rec = hist_to_points(classes["paralog_inter_chrom"]["hist"], bpu)
    intra_pt_rec = hist_to_points(classes["paralog_intra_chrom"]["hist"], bpu)
    for pt, lab, col in ((alle_pt, "allelic (same arm)", C_ALLELE),
                         (intra_pt_rec, "paralog, same chromosome", C_IN),
                         (inter_pt_rec, "paralog, different chromosome", C_OUT)):
        xs, ys = ecdf(pt)
        ax.plot(xs, ys, color=col, lw=1.8, label=f"{lab} (med {median_of(pt)*100:.2f}%)")
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 0.2)
    ax.set_xlabel("substitution divergence (d_sub)")
    ax.set_ylabel("cumulative fraction of aligned bp")
    ax.set_title("A  allelic vs paralogous divergence", fontsize=9, loc="left")
    ax.legend(fontsize=6.5, loc="lower right")
    ax.axvline(0.001, color=C_ALLELE, ls=":", lw=1)
    ax.text(0.00105, 0.04, "allelic median\n0.10%", color=C_ALLELE, fontsize=6)

    # ---- B: fold-change bars
    ax = fig.add_subplot(gs[0, 1])
    base = median_of(alle_pt) if alle_pt else 1e-3
    labs = ["allelic", "paralog\nall (record)", "paralog\nin-community\n(pair)", "paralog\nout-of-community\n(pair)",
            "tightest\npair", "ancient dup\nfamilies"]
    # record-level medians (exact) for allelic/all-paralog; pair-level medians for the
    # community split (labelled as such), which the aggregator does not accumulate per record
    vals = [1.0,
            median_of(inter_pt_rec) / base,
            median_of(intra_pt) / base,
            median_of(inter_pt) / base,
            min(r["d_sub"] for r in pairs if r["d_sub"] == r["d_sub"]) / base,
            100.0]
    cols = [C_ALLELE, C_OUT, C_IN, C_OUT, "#2a9d8f", "#cccccc"]
    bars = ax.bar(range(len(vals)), vals, color=cols, edgecolor="white")
    bars[-1].set_hatch("///")
    ax.set_yscale("log")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labs, fontsize=6.2)
    ax.set_ylabel("fold over allelic median")
    ax.set_title("B  how many X over allelic", fontsize=9, loc="left")
    for i, v in enumerate(vals):
        ax.text(i, v * 1.12, f"{v:.1f}x" if v < 100 else "100-1000x", ha="center", fontsize=6.2)
    ax.axhspan(100, 1000, color="#cccccc", alpha=0.35, zorder=0)
    ax.set_ylim(0.7, 2500)

    # ---- C: per-community cohesion
    ax = fig.add_subplot(gs[1, 0])
    rows = []
    for c in sorted(set(comm.values()), key=lambda z: int(z[1:])):
        memb = {x for x, v in comm.items() if v == c}
        ia = [r["d_sub"] for r in pairs if r["arm_a"] in memb and r["arm_b"] in memb]
        ob = [r["d_sub"] for r in pairs if (r["arm_a"] in memb) != (r["arm_b"] in memb)]
        if len(ia) >= 3 and len(ob) >= 10:
            rows.append((c, st.median(ia) / st.median(ob), len(ia), len(memb)))
    rows.sort(key=lambda r: r[1])
    ys = range(len(rows))
    ax.hlines(list(ys), 1.0, [r[1] for r in rows], color="#999999", lw=1)
    ax.plot([r[1] for r in rows], list(ys), "o", color=C_IN, ms=7)
    ax.axvline(1.0, color="k", lw=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{r[0]} ({r[3]} arms, {r[2]} pairs)" for r in rows], fontsize=7)
    ax.set_xlabel("within-community divergence / out-of-community divergence")
    ax.set_title("C  community cohesion is community-specific", fontsize=9, loc="left")
    ax.set_xlim(0.3, 1.15)
    ax.text(1.005, -0.45, "no cohesion", fontsize=6, ha="right", color="#666666")

    # ---- D: two axes separate
    ax = fig.add_subplot(gs[1, 1])
    for rs, lab, col in ((inter_rows, "out-of-community", C_OUT),
                         (intra_rows, "in-community", C_IN)):
        xs = [r["d_median_record"] for r in rs if r["d_median_record"] == r["d_median_record"]]
        ys2 = [r["frac_at_99.5pct"] for r in rs if r["frac_at_99.5pct"] == r["frac_at_99.5pct"]]
        ax.scatter(xs, ys2, s=9, color=col, alpha=0.55, label=lab, edgecolor="none")
    ax.set_xscale("log")
    ax.set_xlabel("per-record median divergence (magnitude)")
    ax.set_ylabel("fraction of bases at >= 99.5% identity (tail)")
    ax.set_title("D  the near-identical tail separates them, not the magnitude", fontsize=9, loc="left")
    ax.legend(fontsize=6.5, loc="upper right")
    ax.axvline(base, color=C_ALLELE, ls=":", lw=1)
    ax.text(base * 1.05, 0.95, "allelic", color=C_ALLELE, fontsize=6)

    # ---- E: best non-homolog partner per end
    bestp = {}
    for r in pairs:
        if r["arm_a"].rsplit("_", 1)[0] == r["arm_b"].rsplit("_", 1)[0]:
            continue                      # inter-chromosomal partners only
        for e1, e2 in ((r["arm_a"], r["arm_b"]), (r["arm_b"], r["arm_a"])):
            if r["d_median_record"] != r["d_median_record"]:
                continue
            if e1 not in bestp or r["d_median_record"] < bestp[e1][0]:
                bestp[e1] = (r["d_median_record"], e2)
    vals = sorted(v for v, _ in bestp.values())
    axv = fig.add_axes([0.09, -0.16, 0.86, 0.12])
    axv.scatter(vals, [0] * len(vals), s=26, color="#6a4c93", alpha=0.8, zorder=3)
    axv.set_ylim(-1, 1)
    axv.set_xlim(min(vals) * 0.8, max(vals) * 1.25)
    axv.axvline(base, color=C_ALLELE, ls=":", lw=1.2, zorder=2)
    axv.axvline(st.median(vals), color=C_IN, lw=1.5, zorder=2)
    axv.set_xscale("log")
    axv.set_yticks([])
    axv.set_xlabel("closest partner on another chromosome, median divergence per end",
                   fontsize=7.5)
    axv.tick_params(labelsize=7)
    axv.text(base * 1.03, 0.28, "allelic\n0.10%", color=C_ALLELE, fontsize=6)
    axv.text(st.median(vals) * 1.03, 0.28, "median end\n%.2f%%" % (st.median(vals) * 100),
             color=C_IN, fontsize=6)
    axv.set_title("E  every end's closest partner on ANOTHER chromosome", fontsize=9, loc="left")

    fig.suptitle("Subtelomeric divergence: allelic baseline, paralogous copies, and communities "
                 "(HPRCv2, 465 assemblies, %s window PAFs)" % d["summary"].get("n_files", "?"),
                 fontsize=9.5)
    for ext in ("pdf", "png"):
        fig.savefig(a.out + "." + ext, dpi=190, bbox_inches="tight")
    print("wrote", a.out + ".pdf/.png")

    print("\nsummary for caption:")
    print(f"  allelic median       {base*100:.3f}%")
    print(f"  paralog (record, all) {median_of(inter_pt_rec)*100:.3f}% ({median_of(inter_pt_rec)/base:.1f}x)")
    print(f"  paralog in-community {median_of(intra_pt)*100:.3f}% "
          f"({median_of(intra_pt)/base:.1f}x)")
    print(f"  paralog out-community{median_of(inter_pt)*100:.3f}% "
          f"({median_of(inter_pt)/base:.1f}x)")
    print(f"  tightest pair        {min(r['d_sub'] for r in pairs if r['d_sub']==r['d_sub'])*100:.3f}%")
    print(f"  in-community tail    {st.mean([r['frac_at_99.5pct'] for r in intra_rows]):.3f}")
    print(f"  out-community tail   {st.mean([r['frac_at_99.5pct'] for r in inter_rows]):.3f}")


if __name__ == "__main__":
    main()
