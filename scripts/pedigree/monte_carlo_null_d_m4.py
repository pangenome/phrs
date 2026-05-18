#!/usr/bin/env python3
"""
D-M4 Monte Carlo null for the WashU 494/538 = 92% within-Leiden statistic.

Closes OPEN_REVIEWER_CONCERNS.md §D-M4: provides a permutation null distribution
for the pedigree within-community fraction, preserving per-arm patch-count
marginals, so reviewers can compare 92% against random-pairing expectation.

USAGE

  Default (in-repo, uses tabulated WashU within-community subset extracted from
  end-to-end-report/report/14_pedigree_recombination.md; falls back to the
  CEPH1463 11 cross-assembler features embedded from SURVEY_14 §1.6):

      python3 scripts/pedigree/monte_carlo_null_d_m4.py --reps 10000

  Full 538-patch WashU run (requires upstream TSV at /moosefs/.../pedigrees/):

      python3 scripts/pedigree/monte_carlo_null_d_m4.py \\
          --washu-tsv /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv \\
          --reps 10000

OUTPUTS

  paper_prep/figures/fig4/null_distribution_d_m4.png  (figure)
  stdout: textual summary block.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_PATH = REPO_ROOT / "end-to-end-report" / "report" / "14_pedigree_recombination.md"
FIG_OUT_DIR = REPO_ROOT / "paper_prep" / "figures" / "fig4"
FIG_OUT_DIR.mkdir(parents=True, exist_ok=True)

# 15-community arm-level Leiden partition (k=15) for HPRCv2 41 signal-bearing
# subtelomeric arms.  Recovered from
# /home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/enrichment/community_summary_table.tsv
# (Dropbox-vendored copy of the canonical
# hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv).
# Arm naming convention: chrN<p|q>.  41 arms total (acrocentric q-arms and most
# chr2/chr3/chr5/chr8/chr14/chr18 q-arms have no signal and are absent).
LEIDEN_K15 = {
    "C1":  ["chr4q",  "chr10q"],
    "C2":  ["chr10p", "chr18p"],
    "C3":  ["chr3q",  "chr11p", "chr15q", "chr19p"],
    "C4":  ["chr7q",  "chr12q"],
    "C5":  ["chr6p",  "chr9p",  "chr12p", "chr20q"],
    "C6":  ["chr1q",  "chr13q", "chr17q", "chr19q", "chr21q", "chr22q"],
    "C7":  ["chr13p", "chr14p", "chr15p", "chr21p", "chr22p"],
    "C8":  ["chr16p"],
    "C9":  ["chr7p",  "chr9q",  "chr16q"],
    "C10": ["chr17p"],
    "C11": ["chr1p",  "chr5q",  "chr6q",  "chr8p"],
    "C12": ["chr2q",  "chr20p"],
    "C13": ["chr4p"],
    "C14": ["chrXq", "chrYq"],
    "C15": ["chrXp", "chrYp"],
}
ARM_TO_COMM: dict[str, str] = {arm: c for c, arms in LEIDEN_K15.items() for arm in arms}
ALL_ARMS = sorted(ARM_TO_COMM.keys())
assert len(ALL_ARMS) == 41, f"expected 41 arms, got {len(ALL_ARMS)}"
assert len(set(ARM_TO_COMM.values())) == 15

# CEPH1463 cross-assembler validated parent features (n=11) — copied verbatim
# from SURVEY_14 §1.6 / end-to-end-report/report/14_pedigree_recombination.md
# §"CEPH1463 cross-assembler validated parent features (11 total)".  Each row
# = (parent, chrA, chrB, community).  Arm assignment within each chromosome
# (p vs q) is implied by the community (only one of {p,q} for each chromosome
# is a member of any given community in LEIDEN_K15).
CEPH_FEATURES = [
    ("NA12877", "chr1",  "chr19", "C6"),
    ("NA12877", "chr10", "chr18", "C2"),
    ("NA12877", "chr17", "chr19", "C6"),
    ("NA12877", "chr6",  "chr9",  "C5"),
    ("NA12878", "chr10", "chr18", "C2"),
    ("NA12878", "chr19", "chr22", "C6"),
    ("NA12878", "chr21", "chr22", "C7"),
    ("NA12878", "chr6",  "chr9",  "C5"),
    ("NA12889", "chr12", "chr9",  "C5"),
    ("NA12890", "chr12", "chr9",  "C5"),
    ("NA12892", "chr21", "chr22", "C6"),
]


# ---------------------------------------------------------------------------
# Arm normalisation
# ---------------------------------------------------------------------------

_ARM_PATTERNS = [
    (re.compile(r"^(chr\d+|chrX|chrY)_parm$"),   lambda m: m.group(1) + "p"),
    (re.compile(r"^(chr\d+|chrX|chrY)_qarm$"),   lambda m: m.group(1) + "q"),
    (re.compile(r"^(chr\d+|chrX|chrY)_p$"),      lambda m: m.group(1) + "p"),
    (re.compile(r"^(chr\d+|chrX|chrY)_q$"),      lambda m: m.group(1) + "q"),
    (re.compile(r"^(chr\d+|chrX|chrY)([pq])$"),  lambda m: m.group(1) + m.group(2)),
    (re.compile(r"^(chr\d+|chrX|chrY)([pq]):h[12]$"), lambda m: m.group(1) + m.group(2)),
]


def normalise_arm(s: str) -> str | None:
    """Return canonical chrNp / chrNq form or None if unrecognised."""
    if s is None:
        return None
    s = s.strip()
    for pat, fn in _ARM_PATTERNS:
        m = pat.match(s)
        if m:
            return fn(m)
    return None


# ---------------------------------------------------------------------------
# Report-table parsing: extracts WashU within-community patches from
# end-to-end-report/report/14_pedigree_recombination.md when the upstream
# all_pedigrees_patches.tsv is unavailable.
# ---------------------------------------------------------------------------

WASHU_TABLE_HEADERS = [
    "## WashU `gene_conversion_like`",
    "## WashU `crossover_like`",
    "## WashU `acros_like`",
]


def parse_report_washu_within_community(report_path: Path) -> list[tuple[str, str]]:
    """
    Walk the section headers WASHU_TABLE_HEADERS, parse the markdown table
    that immediately follows each, and return (source_arm, target_arm) pairs
    for every row.  All rows in these tables are within-community by
    construction (the report tabulates only the within-community subset).

    Columns of the in-report tables are:
      | # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
    target = Flank, source = Source (strip ':hX' suffix).
    """
    text = report_path.read_text()
    pairs: list[tuple[str, str]] = []

    for header in WASHU_TABLE_HEADERS:
        i = text.find(header)
        if i < 0:
            print(f"WARN: header not found in report: {header!r}", file=sys.stderr)
            continue
        block = text[i : i + 200_000]  # generous slice; tables are short
        rows_added = 0
        for line in block.splitlines():
            if not line.startswith("| ") or "---" in line:
                continue
            cols = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cols) < 12 or cols[0] in ("#", ""):
                continue
            try:
                int(cols[0])  # data row starts with an integer
            except ValueError:
                continue
            flank = cols[2]
            source_raw = cols[5]
            source_arm = normalise_arm(source_raw)
            target_arm = normalise_arm(flank)
            if source_arm is None or target_arm is None:
                print(f"WARN: skip unparsable row in {header!r}: src={source_raw!r} tgt={flank!r}", file=sys.stderr)
                continue
            if source_arm not in ARM_TO_COMM or target_arm not in ARM_TO_COMM:
                print(f"WARN: arm not in Leiden partition: {source_arm} or {target_arm}", file=sys.stderr)
                continue
            pairs.append((source_arm, target_arm))
            rows_added += 1
            # stop at next section break (subsequent table or '## ' line)
        # stop scanning at next '## ' header by clipping above (already 200k)
        # actually clip more carefully:
        next_hdr = block.find("\n## ", 1)
        if next_hdr > 0:
            # re-scan in clipped region
            clipped = block[:next_hdr]
            real_rows = 0
            tmp = []
            for line in clipped.splitlines():
                if not line.startswith("| ") or "---" in line:
                    continue
                cols = [c.strip() for c in line.strip().strip("|").split("|")]
                if len(cols) < 12 or cols[0] in ("#", ""):
                    continue
                try:
                    int(cols[0])
                except ValueError:
                    continue
                flank = cols[2]; source_raw = cols[5]
                sa = normalise_arm(source_raw); ta = normalise_arm(flank)
                if sa is None or ta is None: continue
                if sa not in ARM_TO_COMM or ta not in ARM_TO_COMM: continue
                tmp.append((sa, ta))
                real_rows += 1
            # replace
            pairs = pairs[: len(pairs) - rows_added] + tmp
            print(f"  parsed {real_rows} rows under {header!r}", file=sys.stderr)
        else:
            print(f"  parsed {rows_added} rows under {header!r}", file=sys.stderr)

    return pairs


# ---------------------------------------------------------------------------
# Upstream-TSV path (full 538 WashU patches) — used when --washu-tsv given
# ---------------------------------------------------------------------------

def load_washu_from_upstream_tsv(path: Path) -> list[tuple[str, str]]:
    """
    Load WashU inter-chr patches from
    /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv

    Subset to ds == 'WashU' and inter-chromosomal patches (which is implied by
    the HQ definition: is_interchr=True in patch construction; in this TSV all
    rows are HQ; just filter ds=='WashU').
    Returns every pair (within + cross + unknown community) — the MC null uses
    the full 538-patch marginal, and the observed within-community fraction is
    computed against community_status.
    """
    if not path.exists():
        raise FileNotFoundError(f"upstream TSV not found: {path}")
    pairs: list[tuple[str, str]] = []
    statuses: list[str] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols_needed = {"ds", "query_chr", "query_arm", "ref_chrarm", "community_status"}
        if not cols_needed.issubset(reader.fieldnames or set()):
            raise ValueError(f"expected columns {cols_needed} not in {reader.fieldnames}")
        for row in reader:
            if row["ds"] != "WashU":
                continue
            # query_full_arm = query_chr + '_' + query_arm.strip('arm')
            query_arm = normalise_arm(f"{row['query_chr']}_{row['query_arm']}")
            ref_arm = normalise_arm(row["ref_chrarm"])
            if query_arm is None or ref_arm is None:
                continue
            if query_arm not in ARM_TO_COMM or ref_arm not in ARM_TO_COMM:
                continue
            # inter-chr filter (drop intra-chromosome rows)
            if query_arm[:-1] == ref_arm[:-1]:
                continue
            pairs.append((ref_arm, query_arm))   # (source, target)
            statuses.append(row["community_status"])
    print(f"loaded {len(pairs)} WashU inter-chr patches from {path}", file=sys.stderr)
    obs_within = sum(1 for s in statuses if s == "within_community")
    print(f"  observed within-community per upstream TSV: {obs_within}/{len(statuses)} "
          f"= {obs_within/len(statuses):.3f}", file=sys.stderr)
    return pairs


# ---------------------------------------------------------------------------
# Monte Carlo core
# ---------------------------------------------------------------------------

def within_community_fraction(pairs: list[tuple[str, str]]) -> float:
    if not pairs:
        return float("nan")
    same = sum(1 for s, t in pairs if ARM_TO_COMM[s] == ARM_TO_COMM[t])
    return same / len(pairs)


def per_community_within_counts(pairs: list[tuple[str, str]]) -> dict[str, int]:
    """For each community c, count the patches whose BOTH arms sit in c
    (i.e. within-community patches assigned to c)."""
    out: Counter[str] = Counter()
    for s, t in pairs:
        cs, ct = ARM_TO_COMM[s], ARM_TO_COMM[t]
        if cs == ct:
            out[cs] += 1
    return dict(out)


def monte_carlo_null(
    pairs: list[tuple[str, str]],
    reps: int,
    rng: np.random.Generator,
    enforce_inter_chrom: bool = True,
) -> dict:
    """
    Null model: permute source and target arms independently, preserving the
    per-arm marginal (i.e. each arm appears as a source / target the same
    number of times as observed).  If enforce_inter_chrom, retry assignments
    that would put the same chromosome on both sides (since the input statistic
    is restricted to inter-chr patches).
    """
    n = len(pairs)
    if n == 0:
        raise ValueError("empty pairs")
    sources = np.array([s for s, _ in pairs])
    targets = np.array([t for _, t in pairs])

    source_arms = sources.copy()
    target_arms = targets.copy()

    null_fracs = np.empty(reps, dtype=float)
    null_per_comm = {c: np.zeros(reps, dtype=np.int32) for c in LEIDEN_K15}

    arm_to_idx = {a: i for i, a in enumerate(ALL_ARMS)}
    comm_idx_of_arm = np.array([list(LEIDEN_K15).index(ARM_TO_COMM[a]) for a in ALL_ARMS])
    comm_keys = list(LEIDEN_K15)

    src_idx = np.array([arm_to_idx[a] for a in source_arms])
    tgt_idx = np.array([arm_to_idx[a] for a in target_arms])
    chrom_of_arm = np.array([int(re.match(r"chr(\d+|X|Y)", a).group(1).replace("X","23").replace("Y","24")) for a in ALL_ARMS])

    src_chrom = chrom_of_arm[src_idx]
    tgt_chrom = chrom_of_arm[tgt_idx]

    src_comm = comm_idx_of_arm[src_idx]
    tgt_comm = comm_idx_of_arm[tgt_idx]

    t_start = time.time()
    for r in range(reps):
        # permute sources, permute targets independently — preserves both
        # marginals exactly
        s_perm = rng.permutation(src_idx)
        t_perm = rng.permutation(tgt_idx)

        if enforce_inter_chrom:
            # iteratively fix intra-chr collisions by swapping among s_perm
            # (small n, small loop ok)
            s_chrom = chrom_of_arm[s_perm]
            t_chrom = chrom_of_arm[t_perm]
            bad = np.where(s_chrom == t_chrom)[0]
            attempts = 0
            while len(bad) > 0 and attempts < 50:
                # pick another permutation of just the bad indices in s_perm
                shuffle = rng.permutation(bad)
                s_perm[bad] = s_perm[shuffle]
                s_chrom = chrom_of_arm[s_perm]
                bad = np.where(s_chrom == chrom_of_arm[t_perm])[0]
                attempts += 1
            # if still collisions, leave them — vanishingly rare for these inputs

        same_comm = comm_idx_of_arm[s_perm] == comm_idx_of_arm[t_perm]
        null_fracs[r] = same_comm.mean()
        # per-community counts: how many same-community pairs map to each comm?
        for ci, c in enumerate(comm_keys):
            null_per_comm[c][r] = int(((comm_idx_of_arm[s_perm] == ci) & (comm_idx_of_arm[t_perm] == ci)).sum())

    elapsed = time.time() - t_start
    print(f"  MC: {reps} reps on n={n} pairs in {elapsed:.2f}s ({reps/elapsed:.0f} reps/s)", file=sys.stderr)

    return {
        "n_pairs": n,
        "null_fracs": null_fracs,
        "null_per_comm": null_per_comm,
        "elapsed_s": elapsed,
    }


def summarise(name: str, pairs: list[tuple[str, str]], mc: dict, alpha: float = 0.05) -> dict:
    n = mc["n_pairs"]
    obs_frac = within_community_fraction(pairs)
    obs_count = int(round(obs_frac * n))
    null = mc["null_fracs"]
    mean = float(np.mean(null))
    ci_lo, ci_hi = (float(np.percentile(null, 100 * alpha / 2)),
                    float(np.percentile(null, 100 * (1 - alpha / 2))))
    # one-sided p (enrichment): fraction of null >= observed
    p_enrich = float((null >= obs_frac).mean())
    # depletion p (concern's wording): fraction of null <= observed
    p_deplete = float((null <= obs_frac).mean())
    pct = float((null < obs_frac).mean() * 100)

    # per-community within-community enrichment with BH-corrected q
    obs_per_comm = per_community_within_counts(pairs)
    rows: list[dict] = []
    for c in LEIDEN_K15:
        null_vec = mc["null_per_comm"][c]
        n_arms = len(LEIDEN_K15[c])
        obs_c = obs_per_comm.get(c, 0)
        nm = float(null_vec.mean())
        p_c = float((null_vec >= obs_c).mean()) if obs_c > 0 else 1.0
        rows.append({
            "community": c,
            "n_arms": n_arms,
            "arms": ",".join(LEIDEN_K15[c]),
            "obs_within": obs_c,
            "null_mean": nm,
            "p_enrich": p_c,
        })
    # BH q
    pvals = np.array([r["p_enrich"] for r in rows])
    order = np.argsort(pvals)
    n_tests = len(pvals)
    q = np.empty_like(pvals)
    prev = 1.0
    for rank_back, idx in enumerate(order[::-1]):
        i = n_tests - rank_back  # 1..n
        bh = pvals[idx] * n_tests / i
        prev = min(prev, bh)
        q[idx] = prev
    for r_, qv in zip(rows, q):
        r_["q_bh"] = float(qv)

    print(f"=== {name} ===", file=sys.stderr)
    print(f"  n pairs              : {n}", file=sys.stderr)
    print(f"  observed within-frac : {obs_frac:.4f}  ({obs_count}/{n})", file=sys.stderr)
    print(f"  null mean            : {mean:.4f}", file=sys.stderr)
    print(f"  null 95% CI          : [{ci_lo:.4f}, {ci_hi:.4f}]", file=sys.stderr)
    print(f"  observed percentile  : {pct:.2f}", file=sys.stderr)
    print(f"  p (enrichment)       : {p_enrich:.5f}   one-sided, obs vs null >=", file=sys.stderr)
    print(f"  p (depletion)        : {p_deplete:.5f}   one-sided, obs vs null <=", file=sys.stderr)
    return {
        "name": name,
        "n_pairs": n,
        "observed_within_frac": obs_frac,
        "observed_within_count": obs_count,
        "null_mean": mean,
        "null_ci_lo": ci_lo,
        "null_ci_hi": ci_hi,
        "observed_percentile": pct,
        "p_enrichment": p_enrich,
        "p_depletion": p_deplete,
        "per_community": rows,
    }


def plot_null(
    summaries: list[tuple[str, dict, np.ndarray]],
    out_png: Path,
):
    n_panels = len(summaries)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 4.0), squeeze=False)
    for ax, (name, s, null) in zip(axes[0], summaries):
        ax.hist(null, bins=40, color="#bbbbbb", edgecolor="white")
        ax.axvline(s["observed_within_frac"], color="#d62728", lw=2,
                   label=f"observed = {s['observed_within_frac']:.3f}")
        ax.axvline(s["null_mean"], color="#1f77b4", lw=1.5, ls="--",
                   label=f"null mean = {s['null_mean']:.3f}")
        ax.axvspan(s["null_ci_lo"], s["null_ci_hi"], color="#1f77b4", alpha=0.15,
                   label=f"null 95% CI [{s['null_ci_lo']:.3f}, {s['null_ci_hi']:.3f}]")
        ax.set_xlabel("within-Leiden-community fraction")
        ax.set_ylabel("permutation count")
        ax.set_title(f"{name} (n={s['n_pairs']})")
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.suptitle("D-M4 Monte Carlo null distribution for pedigree within-Leiden fraction",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    pdf = out_png.with_suffix(".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"  wrote {out_png}", file=sys.stderr)
    print(f"  wrote {pdf}", file=sys.stderr)


def ceph_features_to_pairs() -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for parent, ca, cb, comm in CEPH_FEATURES:
        # find the arm of ca within community comm, same for cb
        comm_arms = LEIDEN_K15[comm]
        arm_a = next((a for a in comm_arms if a.startswith(ca + "p") or a.startswith(ca + "q")), None)
        arm_b = next((a for a in comm_arms if a.startswith(cb + "p") or a.startswith(cb + "q")), None)
        if arm_a is None or arm_b is None:
            raise RuntimeError(f"could not map CEPH feature {parent} {ca}/{cb} {comm}")
        # canonicalise (source, target) order alphabetically
        s, t = sorted([arm_a, arm_b])
        pairs.append((s, t))
    return pairs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=10_000, help="Monte Carlo replicates [10000]")
    ap.add_argument("--seed", type=int, default=20260518)
    ap.add_argument("--washu-tsv", type=Path, default=None,
                    help="Optional path to upstream all_pedigrees_patches.tsv (full 538 WashU patches).")
    ap.add_argument("--ceph-denominator", type=int, default=None,
                    help="If provided, the candidate denominator for CEPH1463 before within-Leiden filter (records to report; not used in MC).")
    ap.add_argument("--out-png", type=Path, default=FIG_OUT_DIR / "null_distribution_d_m4.png")
    ap.add_argument("--out-json", type=Path, default=REPO_ROOT / "paper_prep" / "synthesis" / "ANALYSIS_D_M4.json")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # ---- WashU
    if args.washu_tsv is not None:
        print(f"[washu] loading from upstream TSV: {args.washu_tsv}", file=sys.stderr)
        washu_pairs = load_washu_from_upstream_tsv(args.washu_tsv)
        washu_source = f"full upstream all_pedigrees_patches.tsv (n={len(washu_pairs)})"
    else:
        if not REPORT_PATH.exists():
            print(f"FATAL: report not found at {REPORT_PATH}", file=sys.stderr)
            return 2
        print(f"[washu] upstream TSV not provided; reconstructing from {REPORT_PATH.name}", file=sys.stderr)
        washu_pairs = parse_report_washu_within_community(REPORT_PATH)
        washu_source = (f"in-tree report tables (n={len(washu_pairs)}; this is the within-community subset "
                        f"fully tabulated in §14: 133 gene_conv + 16 crossover + 30 acros_like-top30; "
                        f"the remaining 199 within-comm acros_like + 115 sandwich_same_hap + 1 complex + "
                        f"44 cross-comm patches are NOT in-tree and require /moosefs/ upstream TSV)")

    print(f"  WashU source: {washu_source}", file=sys.stderr)
    mc_washu = monte_carlo_null(washu_pairs, args.reps, rng, enforce_inter_chrom=True)
    s_washu = summarise("WashU", washu_pairs, mc_washu)
    s_washu["data_source"] = washu_source

    # ---- CEPH1463 11 cross-assembler features
    ceph_pairs = ceph_features_to_pairs()
    print(f"[ceph] loaded {len(ceph_pairs)} cross-assembler features (within-Leiden by design — all 11/11 = 100%)", file=sys.stderr)
    mc_ceph = monte_carlo_null(ceph_pairs, args.reps, rng, enforce_inter_chrom=True)
    s_ceph = summarise("CEPH1463 cross-assembler", ceph_pairs, mc_ceph)
    s_ceph["data_source"] = "embedded from end-to-end-report/report/14_pedigree_recombination.md §'CEPH1463 cross-assembler validated parent features (11 total)'"
    s_ceph["candidate_denominator_pre_within_leiden_filter"] = args.ceph_denominator  # None if not given

    # ---- output
    plot_null(
        [
            ("WashU", s_washu, mc_washu["null_fracs"]),
            ("CEPH1463 cross-assembler", s_ceph, mc_ceph["null_fracs"]),
        ],
        args.out_png,
    )

    summary_out = {
        "reps": args.reps,
        "seed": args.seed,
        "leiden_partition": LEIDEN_K15,
        "washu": {k: v for k, v in s_washu.items()},
        "ceph": {k: v for k, v in s_ceph.items()},
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary_out, indent=2, default=str))
    print(f"  wrote {args.out_json}", file=sys.stderr)

    print()
    print("=== FINAL SUMMARY ===")
    for s in (s_washu, s_ceph):
        print(f"[{s['name']}]")
        print(f"  observed: {s['observed_within_count']}/{s['n_pairs']} = {s['observed_within_frac']:.4f}")
        print(f"  null mean: {s['null_mean']:.4f}, 95% CI [{s['null_ci_lo']:.4f}, {s['null_ci_hi']:.4f}]")
        print(f"  observed percentile vs null: {s['observed_percentile']:.2f}")
        print(f"  p (enrichment, one-sided)  : {s['p_enrichment']:.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
