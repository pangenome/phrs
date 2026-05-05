#!/usr/bin/env python3
"""
Figure 2 — Within-community heterogeneity, the two-domain model, population history.

Panels:
  (a) Allele-vs-paralog Wilcoxon paired test, per multi-arm community.
  (b) Two-domain model: per-arm Spearman gradient + breakpoint composite (chr4p, chr4q, chr22q anchors).
  (c) Cross-arm superpop enrichment significance + Fst superpop matrix.
  (d) Out-of-Africa dendrogram from the Fst matrix (UPGMA).

Inputs are read directly from the upstream analysis directories under
/moosefs/guarracino/HPRCv2/PHR_III/. See sources.tsv for paths.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
BASE_HET = Path("/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity")
BASE_PLOT = Path("/moosefs/guarracino/HPRCv2/PHR_III/plots")
HERE = Path(__file__).resolve().parent

ALLELE_PARALOG_TSV = BASE_HET / "allele_vs_paralog_distance.tsv"
FST_MATRIX_TSV = BASE_HET / "fst_superpop_matrix.tsv"
CROSS_ARM_TSV = BASE_HET / "cross_arm_superpop_enrichment.tsv"
TWO_DOMAIN_TEST_TSV = BASE_PLOT / "two_domain_test.tsv"
TWO_DOMAIN_CHANGEPOINT_TSV = BASE_PLOT / "two_domain_changepoint.tsv"
TWO_DOMAIN_BINNED_TSV = BASE_PLOT / "two_domain_binned_means.tsv"
ITS_BREAKPOINT_TSV = BASE_PLOT / "its_breakpoint_coloc.tsv"

OUT_PDF = HERE / "figure_fig2.pdf"
OUT_PNG = HERE / "figure_fig2.png"


def read_tsv(path: Path):
    with open(path) as f:
        rdr = csv.DictReader((ln for ln in f if not ln.startswith("#")), delimiter="\t")
        return list(rdr)


def fnum(x):
    if x is None or x == "" or x == "NA":
        return float("nan")
    try:
        return float(x)
    except ValueError:
        return float("nan")


# ----------------------------------------------------------------------
# Panel 2a — allele-vs-paralog
# ----------------------------------------------------------------------
def panel_a(ax):
    rows = read_tsv(ALLELE_PARALOG_TSV)
    comm = [r for r in rows if r["community"] != "OVERALL"]
    overall = next(r for r in rows if r["community"] == "OVERALL")
    comm.sort(key=lambda r: r["community"])

    labels = [r["community"] for r in comm]
    pct_paralog = np.array([fnum(r["pct_paralog_closer"]) for r in comm])
    n_pairs = np.array([int(r["n_pairs"]) for r in comm])
    pvals = np.array([fnum(r["wilcoxon_p"]) for r in comm])
    direction = [r["direction"] for r in comm]

    y = np.arange(len(labels))
    colors = ["#bf3a2b" if d == "paralog_closer" else "#2c6fb1" for d in direction]
    ax.barh(y, pct_paralog, color=colors, edgecolor="black", linewidth=0.5)
    ax.axvline(50, color="grey", linestyle="--", lw=0.7, alpha=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of pairs where paralog is closer than allele", fontsize=8)
    ax.set_title("a  Allele vs. paralog distance",
                 fontsize=9, loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelsize=7)
    ax.invert_yaxis()

    for i, (pct, p, n) in enumerate(zip(pct_paralog, pvals, n_pairs)):
        # log10 p-value annotation
        if p > 0:
            log10p = math.log10(p)
            label = f"p=1e{int(round(log10p))}, n={n}"
        else:
            label = f"p<1e-300, n={n}"
        offset = pct + 2
        if offset > 75:
            offset = pct - 2
            ha = "right"
            tcol = "white"
        else:
            ha = "left"
            tcol = "black"
        ax.text(offset, i, label, va="center", ha=ha, fontsize=6, color=tcol)

    legend_handles = [
        mpatches.Patch(facecolor="#2c6fb1", edgecolor="black", label="allele closer (8/9)"),
        mpatches.Patch(facecolor="#bf3a2b", edgecolor="black", label="paralog closer (C7)"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="lower right",
              frameon=False, handlelength=1.0)
    overall_p = fnum(overall["wilcoxon_p"])
    overall_n = overall["n_pairs"]
    op = "p<1e-300" if overall_p == 0 else f"p={overall_p:.1e}"
    ax.text(0.98, -0.16, f"overall: {overall_n} pairs, {op} (allele closer)",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, style="italic")


# ----------------------------------------------------------------------
# Panel 2b — two-domain composite
# ----------------------------------------------------------------------
def panel_b(ax_top, ax_bottom):
    test_rows = read_tsv(TWO_DOMAIN_TEST_TSV)
    cp_rows = read_tsv(TWO_DOMAIN_CHANGEPOINT_TSV)
    binned = read_tsv(TWO_DOMAIN_BINNED_TSV)
    its = read_tsv(ITS_BREAKPOINT_TSV)

    # ----- top: per-arm Spearman ρ histogram with breakpoint kb scatter
    rho = []
    for r in test_rows:
        v = fnum(r["spearman_rho"])
        if not math.isnan(v):
            rho.append((r["arm"], v, r["supports_two_domain"] == "True"))
    rho.sort(key=lambda x: x[1])
    arms = [r[0] for r in rho]
    rho_vals = np.array([r[1] for r in rho])
    sup = np.array([r[2] for r in rho])

    n_total = len(test_rows)
    n_supports = sum(1 for r in test_rows if r["supports_two_domain"] == "True")
    n_cp = len(cp_rows)
    n_cp_better = sum(1 for r in cp_rows if r["two_segment_better"] == "True")

    x = np.arange(len(rho_vals))
    cols = np.where(sup, "#2c6fb1", "#cccccc")
    ax_top.bar(x, rho_vals, color=cols, edgecolor="none")
    ax_top.axhline(0, color="black", lw=0.6)
    ax_top.set_ylabel(r"Spearman $\rho$ (n_chrs vs distance)", fontsize=8)
    ax_top.set_xticks([])
    ax_top.set_xlabel(f"per-arm gradient — supports two-domain in {n_supports}/{n_total} arms",
                      fontsize=7)
    ax_top.set_ylim(-1.05, 0.6)
    ax_top.set_title("b  Two-domain model: gradient + breakpoint",
                     fontsize=9, loc="left", fontweight="bold")
    ax_top.tick_params(axis="y", labelsize=7)

    # Highlight focus arms (chr4_p, chr4_q, chr22_q)
    focus = {"chr4_p", "chr4_q", "chr22_q"}
    for i, a in enumerate(arms):
        if a in focus:
            ax_top.text(i, rho_vals[i] - 0.06, a.replace("_p", "p").replace("_q", "q"),
                        rotation=90, ha="center", va="top", fontsize=6,
                        color="#bf3a2b", fontweight="bold")
            ax_top.bar([i], [rho_vals[i]], color="#bf3a2b", edgecolor="black", linewidth=0.6)

    # Median rho annotation
    median_rho = float(np.median(rho_vals))
    ax_top.text(0.98, 0.95,
                f"median $\\rho$ = {median_rho:.2f}; piecewise > linear in {n_cp_better}/{n_cp} arms",
                transform=ax_top.transAxes, fontsize=7, color="black",
                ha="right", va="top", style="italic",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1))

    # ----- bottom: binned-mean curves for focus arms with breakpoint markers
    cp_by_arm = {r["arm"]: r for r in cp_rows}
    its_by_arm = {r["arm"]: r for r in its}
    focus_pairs = [
        ("chr4_p", "#1f77b4"),
        ("chr4_q", "#2ca02c"),
        ("chr22_q", "#d62728"),
    ]
    bp_lines = []
    for arm, color in focus_pairs:
        sub = [r for r in binned if r["chr_arm"] == arm]
        sub.sort(key=lambda r: fnum(r["dist_from_telomere"]))
        if not sub:
            continue
        d = np.array([fnum(r["dist_from_telomere"]) / 1000.0 for r in sub])  # kb
        n = np.array([fnum(r["mean_n_chrs"]) for r in sub])
        ax_bottom.plot(d, n, color=color, lw=1.4,
                       label=arm.replace("_p", "p").replace("_q", "q"))
        cp = cp_by_arm.get(arm)
        if cp:
            bp_kb = fnum(cp["breakpoint_kb"])
            ax_bottom.axvline(bp_kb, color=color, lw=0.8, linestyle="--", alpha=0.7)
            bp_lines.append((arm.replace("_p", "p").replace("_q", "q"), bp_kb, color))
    # Stack breakpoint annotations to avoid overlap
    if bp_lines:
        txt = "\n".join(f"{a}: bp={bp:.0f} kb" for a, bp, _ in bp_lines)
        ax_bottom.text(0.98, 0.95, txt, transform=ax_bottom.transAxes,
                       fontsize=6, ha="right", va="top",
                       bbox=dict(facecolor="white", edgecolor="grey", alpha=0.85, pad=2,
                                 linewidth=0.4))

    ax_bottom.set_xscale("symlog", linthresh=1)
    ax_bottom.set_xlim(left=1)
    ax_bottom.set_xlabel("distance from telomere (kb, symlog)", fontsize=8)
    ax_bottom.set_ylabel("mean # of chromosomes per window", fontsize=8)
    ax_bottom.tick_params(axis="both", labelsize=7)
    ax_bottom.legend(fontsize=7, loc="upper left", frameon=False, ncol=3, handlelength=1.5)
    ax_bottom.grid(axis="both", alpha=0.2, lw=0.4)

    # ITS coloc note
    near_bp = [(r["arm"], fnum(r["closest_its_to_bp_kb"])) for r in its]
    near_bp.sort(key=lambda x: x[1])
    n_within_50 = sum(1 for _, d in near_bp if d <= 50)
    ax_bottom.text(0.02, 0.02,
                   f"ITS-to-breakpoint: {n_within_50}/{len(its)} arms ≤ 50 kb",
                   transform=ax_bottom.transAxes, fontsize=7, ha="left", va="bottom",
                   style="italic")


# ----------------------------------------------------------------------
# Panel 2c — superpop composition + Fst matrix composite
# ----------------------------------------------------------------------
def panel_c(ax_left, ax_right):
    rows = read_tsv(CROSS_ARM_TSV)
    rows = [r for r in rows if r["arm"]]
    # Build a simple "significance heatmap": community x arm with -log10 BH-p
    pairs = [(r["community"], r["arm"], fnum(r["p_adjusted"])) for r in rows]

    # Sort communities and arms
    comms = sorted(set(p[0] for p in pairs))
    arms = sorted(set(p[1] for p in pairs))
    M = np.full((len(comms), len(arms)), np.nan)
    for c, a, p in pairs:
        i, j = comms.index(c), arms.index(a)
        M[i, j] = -math.log10(p) if p > 0 else 6
    sig_threshold = -math.log10(0.05)

    cmap = plt.cm.viridis
    im = ax_left.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=4)
    n_sig = int(np.nansum(M >= sig_threshold))
    n_total = int(np.sum(~np.isnan(M)))

    ax_left.set_xticks(range(len(arms)))
    ax_left.set_xticklabels([a.replace("_qarm", "q").replace("_parm", "p")
                             for a in arms], rotation=90, fontsize=7)
    ax_left.set_yticks(range(len(comms)))
    ax_left.set_yticklabels(comms, fontsize=7)

    # Mark significant cells
    for i in range(len(comms)):
        for j in range(len(arms)):
            v = M[i, j]
            if not np.isnan(v) and v >= sig_threshold:
                ax_left.text(j, i, "*", ha="center", va="center",
                             color="white", fontsize=8, fontweight="bold")

    ax_left.set_title("c  Cross-arm superpop enrichment",
                      fontsize=9, loc="left", fontweight="bold")
    ax_left.set_xlabel(f"arm — {n_sig}/{n_total} pairs BH-significant (q<0.05)", fontsize=7)
    cbar = plt.colorbar(im, ax=ax_left, fraction=0.046, pad=0.04)
    cbar.set_label(r"$-\log_{10}$ q", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # ----- right: Fst matrix
    fst_rows = read_tsv(FST_MATRIX_TSV)
    fst_rows = [r for r in fst_rows if r.get("superpop") and not r["superpop"].startswith("#")]
    pops = [r["superpop"] for r in fst_rows]
    F = np.array([[max(0.0, fnum(r[p])) for p in pops] for r in fst_rows])
    im2 = ax_right.imshow(F, aspect="auto", cmap="magma", vmin=0, vmax=0.16)
    ax_right.set_xticks(range(len(pops)))
    ax_right.set_xticklabels(pops, fontsize=7)
    ax_right.set_yticks(range(len(pops)))
    ax_right.set_yticklabels(pops, fontsize=7)
    for i in range(len(pops)):
        for j in range(len(pops)):
            v = F[i, j]
            tcol = "white" if v < 0.08 else "black"
            ax_right.text(j, i, f"{v:.3f}", ha="center", va="center",
                          color=tcol, fontsize=6)
    ax_right.set_title("Hudson Fst (superpop)", fontsize=8, loc="center")
    cbar2 = plt.colorbar(im2, ax=ax_right, fraction=0.046, pad=0.04)
    cbar2.set_label("Fst", fontsize=7)
    cbar2.ax.tick_params(labelsize=6)

    # AFR vs non-AFR Fst range annotation
    afr_idx = pops.index("AFR")
    non_afr = [F[afr_idx, j] for j in range(len(pops)) if j != afr_idx]
    ax_right.text(0.0, -0.18,
                  f"AFR vs non-AFR: Fst {min(non_afr):.2f}–{max(non_afr):.2f}",
                  transform=ax_right.transAxes, fontsize=7, ha="left", va="top",
                  style="italic")


# ----------------------------------------------------------------------
# Panel 2d — Out-of-Africa dendrogram from Fst (UPGMA, hand-rolled)
# ----------------------------------------------------------------------
def upgma(D, labels):
    """Return list of merge events (i, j, d, members) in agglomerative order."""
    D = D.copy().astype(float)
    n = D.shape[0]
    np.fill_diagonal(D, np.inf)
    clusters = {i: [labels[i]] for i in range(n)}
    sizes = {i: 1 for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    next_id = n
    events = []
    active = set(range(n))

    while len(active) > 1:
        # find min in active submatrix
        sub_idx = sorted(active)
        best = (None, None, np.inf)
        for ai in sub_idx:
            for bi in sub_idx:
                if ai >= bi:
                    continue
                if D[ai, bi] < best[2]:
                    best = (ai, bi, D[ai, bi])
        i, j, d = best
        new_id = next_id
        next_id += 1
        new_members = clusters[i] + clusters[j]
        new_size = sizes[i] + sizes[j]
        new_height = d / 2.0
        events.append({
            "id": new_id,
            "left": i,
            "right": j,
            "left_height": heights[i],
            "right_height": heights[j],
            "height": new_height,
            "members": new_members,
        })
        # Update D with weighted average distances
        new_row = np.full(D.shape[0] + 1, np.inf)
        for k in active:
            if k == i or k == j:
                continue
            new_row[k] = (sizes[i] * D[i, k] + sizes[j] * D[j, k]) / new_size
        # Add a new row/column
        D = np.pad(D, ((0, 1), (0, 1)), mode="constant", constant_values=np.inf)
        for k in active:
            if k == i or k == j:
                continue
            D[new_id, k] = new_row[k]
            D[k, new_id] = new_row[k]
        active.discard(i)
        active.discard(j)
        active.add(new_id)
        clusters[new_id] = new_members
        sizes[new_id] = new_size
        heights[new_id] = new_height

    return events, clusters


def panel_d(ax):
    fst_rows = read_tsv(FST_MATRIX_TSV)
    fst_rows = [r for r in fst_rows if r.get("superpop") and not r["superpop"].startswith("#")]
    pops = [r["superpop"] for r in fst_rows]
    F = np.array([[max(0.0, fnum(r[p])) for p in pops] for r in fst_rows])
    np.fill_diagonal(F, 0)

    events, clusters = upgma(F, pops)

    # Compute leaf x-positions in the order traversed by final tree
    last = events[-1]
    members_order = clusters[last["id"]]
    leaf_x = {lab: i for i, lab in enumerate(members_order)}

    # Build cluster x-positions recursively
    def cluster_x(node_id):
        if node_id < len(pops):
            return leaf_x[pops[node_id]]
        ev = next(e for e in events if e["id"] == node_id)
        return 0.5 * (cluster_x(ev["left"]) + cluster_x(ev["right"]))

    def cluster_h(node_id):
        if node_id < len(pops):
            return 0.0
        ev = next(e for e in events if e["id"] == node_id)
        return ev["height"]

    # Draw
    for ev in events:
        lx = cluster_x(ev["left"])
        rx = cluster_x(ev["right"])
        lh = cluster_h(ev["left"])
        rh = cluster_h(ev["right"])
        h = ev["height"]
        ax.plot([lx, lx], [lh, h], color="black", lw=1.2)
        ax.plot([rx, rx], [rh, h], color="black", lw=1.2)
        ax.plot([lx, rx], [h, h], color="black", lw=1.2)

    # Color leaves AFR in red, others in blue
    for lab, x in leaf_x.items():
        col = "#bf3a2b" if lab == "AFR" else "#2c6fb1"
        ax.text(x, -0.005, lab, ha="center", va="top",
                fontsize=8, fontweight="bold", color=col)

    # Cosmetic
    max_h = max(e["height"] for e in events)
    ax.set_ylim(-0.02, max_h * 1.15)
    ax.set_xlim(-0.5, len(pops) - 0.5)
    ax.set_xticks([])
    ax.set_ylabel("Fst / 2  (UPGMA height)", fontsize=8)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title("d  Out-of-Africa tree from cross-arm Fst",
                 fontsize=9, loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Annotation
    afr_idx = pops.index("AFR")
    non_afr = [F[afr_idx, j] for j in range(len(pops)) if j != afr_idx]
    ax.text(0.02, 0.95,
            f"AFR branch separates first; AFR–nonAFR Fst {min(non_afr):.2f}–{max(non_afr):.2f}",
            transform=ax.transAxes, fontsize=7, ha="left", va="top", style="italic")


# ----------------------------------------------------------------------
# Compose full figure
# ----------------------------------------------------------------------
def main():
    fig = plt.figure(figsize=(11.5, 9.2))
    gs_outer = gridspec.GridSpec(2, 2, figure=fig,
                                 hspace=0.55, wspace=0.35,
                                 left=0.07, right=0.97, bottom=0.07, top=0.95)
    # 2a top-left
    ax_a = fig.add_subplot(gs_outer[0, 0])
    panel_a(ax_a)

    # 2b top-right: split into two
    gs_b = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_outer[0, 1], hspace=0.55, height_ratios=[1, 1])
    ax_b_top = fig.add_subplot(gs_b[0])
    ax_b_bot = fig.add_subplot(gs_b[1])
    panel_b(ax_b_top, ax_b_bot)

    # 2c bottom-left: split into two columns
    gs_c = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_outer[1, 0], wspace=0.6, width_ratios=[1.2, 1.0])
    ax_c_l = fig.add_subplot(gs_c[0])
    ax_c_r = fig.add_subplot(gs_c[1])
    panel_c(ax_c_l, ax_c_r)

    # 2d bottom-right
    ax_d = fig.add_subplot(gs_outer[1, 1])
    panel_d(ax_d)

    fig.suptitle("Figure 2 — Within-community heterogeneity, the two-domain model, and population history",
                 fontsize=11, y=0.99, fontweight="bold")

    fig.savefig(OUT_PDF, bbox_inches="tight", dpi=300)
    fig.savefig(OUT_PNG, bbox_inches="tight", dpi=200)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
