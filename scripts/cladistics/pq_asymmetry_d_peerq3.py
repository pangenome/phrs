#!/usr/bin/env python3
"""
D-PeerQ3: P-vs-Q arm orientation asymmetry per Leiden k=15 community.

For each of the 15 arm-level Leiden communities, counts within-community
arm pairs by orientation type (P-P, Q-Q, P-Q) and tests whether same-
orientation pairs are enriched relative to a background of all pairwise
arm comparisons (Fisher exact test, BH-corrected).

Community assignments are sourced from:
  slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap/arm_order_community.tsv

These in turn derive from:
  /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv
  /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv
"""

import math
import sys
from itertools import combinations
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Arm-level community assignments (Leiden k=15, 41 signal-bearing arms)
# Sourced from arm_order_community.tsv (committed artifact)
# ---------------------------------------------------------------------------
ARM_COMMUNITIES = {
    "10q": "C1", "4q": "C1",
    "10p": "C2", "18p": "C2",
    "11p": "C3", "19p": "C3", "3q": "C3", "9q": "C3", "16q": "C3", "7p": "C3",
    "12q": "C4", "7q": "C4",
    "9p": "C5", "6p": "C5", "12p": "C5", "20q": "C5",
    "17q": "C6", "13q": "C6", "1q": "C6", "19q": "C6", "21q": "C6", "22q": "C6",
    "13p": "C7", "21p": "C7", "14p": "C7", "15p": "C7", "22p": "C7",
    "15q": "C8",
    "16p": "C9",
    "17p": "C10",
    "1p": "C11", "8p": "C11", "5q": "C11", "6q": "C11",
    "2q": "C12", "20p": "C12",
    "4p": "C13",
    "Xq": "C14", "Yq": "C14",
    "Xp": "C15", "Yp": "C15",
}

COMMUNITIES_ORDERED = [f"C{i}" for i in range(1, 16)]


def arm_side(arm: str) -> str:
    """Return 'p' or 'q' for an arm label like '10q' or 'Xp'."""
    return arm[-1]


def edge_type(arm1: str, arm2: str) -> str:
    s1, s2 = arm_side(arm1), arm_side(arm2)
    if s1 == s2:
        return "same"  # P-P or Q-Q
    return "cross"     # P-Q


def classify_pair(arm1: str, arm2: str) -> str:
    s1, s2 = arm_side(arm1), arm_side(arm2)
    if s1 == "p" and s2 == "p":
        return "PP"
    if s1 == "q" and s2 == "q":
        return "QQ"
    return "PQ"


def main():
    arms = list(ARM_COMMUNITIES.keys())
    n_total = len(arms)  # 41

    # Verify counts
    n_p = sum(1 for a in arms if arm_side(a) == "p")
    n_q = sum(1 for a in arms if arm_side(a) == "q")
    print(f"Total arms: {n_total}  P-arms: {n_p}  Q-arms: {n_q}", file=sys.stderr)
    assert n_total == 41, f"Expected 41 arms, got {n_total}"

    # Global pair counts (background for Fisher test)
    global_pp = math.comb(n_p, 2)   # C(21,2) = 210
    global_qq = math.comb(n_q, 2)   # C(20,2) = 190
    global_pq = n_p * n_q            # 21*20  = 420
    global_same = global_pp + global_qq  # 400
    global_total = math.comb(n_total, 2)  # 820
    print(f"Global: PP={global_pp} QQ={global_qq} PQ={global_pq} same={global_same} total={global_total}",
          file=sys.stderr)

    # Group arms by community
    community_arms: dict[str, list[str]] = {c: [] for c in COMMUNITIES_ORDERED}
    for arm, comm in ARM_COMMUNITIES.items():
        community_arms[comm].append(arm)

    # Per-community edge counts
    rows = []
    comm_labels = []
    p_vals = []

    for comm in COMMUNITIES_ORDERED:
        members = community_arms[comm]
        n = len(members)
        np_c = sum(1 for a in members if arm_side(a) == "p")
        nq_c = n - np_c

        obs_pp = math.comb(np_c, 2)
        obs_qq = math.comb(nq_c, 2)
        obs_pq = np_c * nq_c
        obs_same = obs_pp + obs_qq
        obs_total = math.comb(n, 2)

        # Fisher exact test: same-orientation enrichment within community vs background
        # 2x2 table:
        #            | same | cross |
        # within comm| obs_same | obs_pq |
        # outside    | global_same - obs_same | global_pq - obs_pq |
        out_same = global_same - obs_same
        out_pq = global_pq - obs_pq

        if obs_total == 0:
            # Single-arm community: no pairs, test undefined
            pval = float("nan")
            OR = float("nan")
        else:
            table = [[obs_same, obs_pq], [out_same, out_pq]]
            OR, pval = fisher_exact(table, alternative="greater")

        # Expected P-Q under null (global P/Q ratio)
        p_frac = n_p / n_total
        q_frac = n_q / n_total
        exp_pq_frac = 2 * p_frac * q_frac   # expected fraction of pairs that are P-Q
        exp_pq = exp_pq_frac * obs_total if obs_total > 0 else 0.0
        obs_pq_frac = obs_pq / obs_total if obs_total > 0 else float("nan")

        orientation_class = "pure-P" if nq_c == 0 else ("pure-Q" if np_c == 0 else "mixed")

        arms_str = ", ".join(sorted(members))
        rows.append({
            "community": comm,
            "members": arms_str,
            "n": n,
            "n_p": np_c,
            "n_q": nq_c,
            "orientation_class": orientation_class,
            "PP": obs_pp,
            "QQ": obs_qq,
            "PQ": obs_pq,
            "total_edges": obs_total,
            "PQ_frac_obs": obs_pq_frac,
            "PQ_frac_exp": exp_pq,
            "OR": OR,
            "p_fisher": pval,
        })
        comm_labels.append(comm)
        p_vals.append(pval)

    # BH correction (exclude NaN)
    valid_idx = [i for i, p in enumerate(p_vals) if not math.isnan(p)]
    valid_pvals = [p_vals[i] for i in valid_idx]
    if valid_pvals:
        reject, p_adj, _, _ = multipletests(valid_pvals, method="fdr_bh")
        adj_map = {valid_idx[i]: p_adj[i] for i in range(len(valid_idx))}
    else:
        adj_map = {}

    for i, row in enumerate(rows):
        row["p_adj_BH"] = adj_map.get(i, float("nan"))
        row["significant"] = (not math.isnan(row["p_adj_BH"])) and row["p_adj_BH"] < 0.05

    # Print results table
    print("\nD-PeerQ3: Per-community P-P / Q-Q / P-Q edge counts")
    print("=" * 100)
    hdr = f"{'Comm':<5} {'Members':<50} {'n':>2} {'nP':>3} {'nQ':>3} {'PP':>4} {'QQ':>4} {'PQ':>4} {'PQ%obs':>8} {'OR':>8} {'p_fisher':>10} {'p_BH':>10} {'sig':>4}"
    print(hdr)
    print("-" * 100)
    for row in rows:
        pq_pct = f"{row['PQ_frac_obs']*100:.1f}%" if not math.isnan(row["PQ_frac_obs"]) else "  N/A"
        or_str = f"{row['OR']:.2f}" if not math.isnan(row["OR"]) else "  N/A"
        pf_str = f"{row['p_fisher']:.4f}" if not math.isnan(row["p_fisher"]) else "   N/A"
        pb_str = f"{row['p_adj_BH']:.4f}" if not math.isnan(row["p_adj_BH"]) else "   N/A"
        sig_str = "*" if row["significant"] else ""
        print(f"{row['community']:<5} {row['members']:<50} {row['n']:>2} {row['n_p']:>3} {row['n_q']:>3} "
              f"{row['PP']:>4} {row['QQ']:>4} {row['PQ']:>4} {pq_pct:>8} {or_str:>8} {pf_str:>10} {pb_str:>10} {sig_str:>4}")

    print()
    print("Notes:")
    print("  OR = odds ratio for same-orientation enrichment (within-comm vs outside); alternative='greater'")
    print("  p_fisher = one-sided Fisher exact p-value; p_BH = BH-adjusted; * = p_BH < 0.05")
    print("  Singleton communities (C8, C9, C10, C13) have 0 edges: test undefined (NaN)")
    print()

    return rows


if __name__ == "__main__":
    main()
