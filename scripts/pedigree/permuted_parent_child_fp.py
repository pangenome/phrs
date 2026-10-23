#!/usr/bin/env python3
"""
Permuted parent-child false-positive rate for the WashU pedigree subtelomere scan.

Referee #1, major point 2, asks for "a false-positive estimate from permuted
parent-child pairs".  The scan (POINT_BY_POINT R1-M2) is:

  for every 2 kb window of the PAN027 paternal query (= SELF.p, the child's
  paternal haplotype) against the father PAN011 (both haplotypes), a window is
  candidate-eligible only when the best NON-homologous parental match beats the
  best same-chromosome (homologous) parental match (positive identity margin),
  or when no homologous parental match exists at all in that window, AND the
  non-homologous identity is >= 0.95.  A tract is a run of adjacent (touching)
  eligible windows sharing one donor chromosome with total length >= 10 kb and
  mean identity >= 0.95.

WHAT IS PERMUTED, AND WHAT THAT MEANS

  The referee also suggests permuting which parental haplotype is treated as the
  transmitted copy.  That is NOT possible from this input: the class-winner
  table stores only the winner WITHIN the `same_chrom` class, which is already
  the maximum over the father's two haplotypes (both PAN011 h1 and h2 appear as
  `same_chrom` winners across windows; counts are reported in this JSON).  The
  per-haplotype identities required to swap the transmitted copy are not in
  this file.  We use instead the alternative the referee allows -- shuffle
  window labels preserving per-chromosome (per-arm) window counts:

  null A -- permute the complete per-window record (which window carries a
  positive-margin paralogous winner, and with which donor chromosome and
  identity) within each query chromosome.  Per chromosome this preserves the
  number of scanned windows, the number of windows in which the paralogue beats
  the homologue, and the joint (donor, identity) distribution of those windows.
  It destroys spatial order, hence both the clustering of candidate windows
  along a chromosome and the contiguity of a single donor.  Implemented exactly:
  a uniform permutation of the N_c records over the N_c window positions places
  the K_c eligible records on a uniformly random ordered K_c-subset of the
  positions, so positions are sampled without replacement and the shuffled
  eligible records are assigned to them.

  null B -- leave every window's own eligible/not status and identity in place
  and shuffle only the DONOR LABEL among the eligible windows of each
  chromosome.  This conditions on the observed pattern of "where the paralogous
  match wins" and tests the remaining claim: that a run shares one donor more
  often than chance given that pattern.

  Neither null is an orthogonal assay.  They bound the false-positive rate under
  a stated null; see the `limitations` block in the output.

STATISTICS REPORTED FOR THE OBSERVED DATA AND BOTH NULLS

  S1  total number of qualifying runs genome-wide
  S2  per candidate interval: number of qualifying runs overlapping it, and the
      longest qualifying run overlapping it
  S3  per candidate interval: number of eligible windows falling inside it
      (local density of paralogous winners), which the nulls can actually reach

OUTPUTS

  scripts/pedigree/permuted_parent_child_fp.json

USAGE

  python3 scripts/pedigree/permuted_parent_child_fp.py
  python3 scripts/pedigree/permuted_parent_child_fp.py --reps 20000 --seed 7
"""

import argparse
import json
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from identity_margin_distribution import (  # noqa: E402
    WINDOW, MIN_RUN_WINDOWS, MIN_RUN_MEAN_IDENTITY, MIN_WINDOW_IDENTITY,
    MARKER_CHROMS, TERMINAL_BP,
    load_windows, load_fai, is_run_eligible, state, margin,
    find_runs, qualifying, derive_candidates,
)


# --------------------------------------------------------------------------- #
# run detection over an explicit (position, donor, identity) sequence
# --------------------------------------------------------------------------- #
def runs_from_sequence(positions, donors, identities, stride=WINDOW):
    """Runs of adjacent positions (exactly `stride` bp apart) sharing a donor.

    Returns [(start, end, donor, n_windows, mean_identity)] for runs passing the
    >= MIN_RUN_WINDOWS and mean-identity filters.
    """
    out = []
    n = len(positions)
    i = 0
    while i < n:
        j = i + 1
        while (j < n and positions[j] == positions[j - 1] + stride
               and donors[j] == donors[i]):
            j += 1
        k = j - i
        if k >= MIN_RUN_WINDOWS:
            mean_id = sum(identities[i:j]) / k
            if mean_id >= MIN_RUN_MEAN_IDENTITY:
                out.append((positions[i], positions[j - 1] + stride, donors[i],
                            k, mean_id))
        i = j
    return out


def scan_runs(runs, intervals):
    """(total runs, {interval: n runs overlapping}, {interval: longest overlap bp})."""
    per_n = {name: 0 for name in intervals}
    per_max = {name: 0 for name in intervals}
    for (start, end, _donor, _k, _m) in runs:
        for name, (lo, hi) in intervals.items():
            if start < hi and end > lo:
                per_n[name] += 1
                per_max[name] = max(per_max[name], end - start)
    return len(runs), per_n, per_max


def summarise(values):
    vals = sorted(values)
    n = len(vals)
    return {
        "n": n,
        "mean": round(statistics.fmean(vals), 6),
        "sd": round(statistics.stdev(vals), 6) if n > 1 else 0.0,
        "min": vals[0],
        "max": vals[-1],
        "quantiles": {("%g" % q): vals[max(0, min(n - 1, int(round(q * (n - 1)))))]
                      for q in (0.0, 0.5, 0.9, 0.95, 0.99, 1.0)},
    }


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--similarity",
                    default="data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz")
    ap.add_argument("--fai",
                    default="data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai")
    ap.add_argument("--out", default="scripts/pedigree/permuted_parent_child_fp.json")
    ap.add_argument("--reps", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261016)
    args = ap.parse_args()

    windows = load_windows(args.similarity)
    chrom_lengths = load_fai(args.fai)

    # ---- audit: `same_chrom` winner is already the max over both haplotypes --
    same_chrom_hap = {}
    for rec in windows.values():
        if "same_chrom" in rec:
            seq = rec["same_chrom"]["other_seq"]
            hap = seq.split("_chr")[0]
            same_chrom_hap[hap] = same_chrom_hap.get(hap, 0) + 1

    # ---- observed ----
    obs_runs = qualifying(find_runs(windows))
    obs_runs.sort(key=lambda r: (-(r["end"] - r["start"]), r["chrom"], r["start"]))
    candidates = derive_candidates(windows, chrom_lengths, MARKER_CHROMS)
    intervals = {chrom: (candidates[chrom]["start"], candidates[chrom]["end"])
                 for chrom in MARKER_CHROMS}

    obs_runs_tuples = [(r["start"], r["end"], r["donor"], r["n"], r["mean_identity"])
                       for r in obs_runs]
    length_rank = {(r["chrom"], r["start"]): i + 1 for i, r in enumerate(obs_runs)}
    obs_total, obs_interval_n, obs_interval_max = scan_runs(obs_runs_tuples, intervals)

    # ---- per-chromosome machinery ----
    pos_by_chrom = {}
    for rec in windows.values():
        pos_by_chrom.setdefault(rec["chrom"], []).append(rec["start"])
    for c in pos_by_chrom:
        pos_by_chrom[c].sort()

    eligible = {}   # chrom -> sorted [(start, donor, identity)]
    for rec in windows.values():
        donor = is_run_eligible(rec)
        if donor is None:
            continue
        eligible.setdefault(rec["chrom"], []).append(
            (rec["start"], donor, rec["interchrom"]["identity"]))
    for c in eligible:
        eligible[c].sort()

    obs_interval_windows = {name: 0 for name in intervals}
    for c in intervals:
        lo, hi = intervals[c]
        obs_interval_windows[c] = sum(1 for rec in windows.values()
                                      if rec["chrom"] == c and lo <= rec["start"] < hi
                                      and is_run_eligible(rec) is not None)

    chroms = sorted(pos_by_chrom)
    rng = random.Random(args.seed)

    # ---- null A ----
    A_counts, A_win, A_maxlen = [], {k: [] for k in intervals}, {k: [] for k in intervals}
    A_nper = {k: [] for k in intervals}
    for _ in range(args.reps):
        total = 0
        per_n = {name: 0 for name in intervals}
        per_max = {name: 0 for name in intervals}
        inwin = {name: 0 for name in intervals}
        for c in chroms:
            el = eligible.get(c)
            if not el:
                continue
            plist = pos_by_chrom[c]
            k = len(el)
            idx = sorted(rng.sample(range(len(plist)), k))
            recs = el[:]
            rng.shuffle(recs)
            positions = [plist[i] for i in idx]
            donors = [r[1] for r in recs]
            idents = [r[2] for r in recs]
            runs = runs_from_sequence(positions, donors, idents)
            sub = {name: intervals[name] for name in intervals if name == c}
            if sub:
                total_i, n_i, max_i = scan_runs(runs, sub)
            else:
                total_i, n_i, max_i = len(runs), {}, {}
            total += total_i
            for name in n_i:
                per_n[name] += n_i[name]
                per_max[name] = max(per_max[name], max_i[name])
            if c in intervals:
                lo, hi = intervals[c]
                inwin[c] = sum(1 for p in positions if lo <= p < hi)
        A_counts.append(total)
        for name in intervals:
            A_win[name].append(per_n[name])
            A_maxlen[name].append(per_max[name])
            A_nper[name].append(inwin[name])

    rng = random.Random(args.seed + 1)
    # ---- null B: donor-label shuffle among eligible windows of a chromosome ----
    B_counts, B_win = [], {k: [] for k in intervals}
    B_maxlen = {k: [] for k in intervals}
    for _ in range(args.reps):
        total = 0
        per_n = {name: 0 for name in intervals}
        per_max = {name: 0 for name in intervals}
        for c in chroms:
            el = eligible.get(c)
            if not el:
                continue
            positions = [r[0] for r in el]
            donors = [r[1] for r in el]
            idents = [r[2] for r in el]
            rng.shuffle(donors)                # donors only; positions and ids fixed
            runs = runs_from_sequence(positions, donors, idents)
            total += len(runs)
            sub = {name: intervals[name] for name in intervals if name == c}
            if sub:
                _t, n_i, max_i = scan_runs(runs, sub)
                for name in n_i:
                    per_n[name] += n_i[name]
                    per_max[name] = max(per_max[name], max_i[name])
        B_counts.append(total)
        for name in intervals:
            B_win[name].append(per_n[name])
            B_maxlen[name].append(per_max[name])

    # ---- assemble ----
    def p_ge(obs, vals):
        return sum(1 for v in vals if v >= obs) / len(vals)

    report = {
        "generated_by": "scripts/pedigree/permuted_parent_child_fp.py",
        "inputs": {
            "similarity_table": args.similarity,
            "query_fai": args.fai,
            "comparison": "PAN027 paternal (PAN027#2) query vs PAN011 (father, both haplotypes) target",
        },
        "parameters": {
            "window_bp": WINDOW,
            "min_window_identity": MIN_WINDOW_IDENTITY,
            "min_run_windows": MIN_RUN_WINDOWS,
            "min_run_bp": MIN_RUN_WINDOWS * WINDOW,
            "min_run_mean_identity": MIN_RUN_MEAN_IDENTITY,
            "reps": args.reps,
            "seed_null_A": args.seed,
            "seed_null_B": args.seed + 1,
            "marker_chroms": MARKER_CHROMS,
            "terminal_bp_for_candidate_derivation": TERMINAL_BP,
            "p_value_resolution": "1/reps; a reported p_ge_observed of 0.0 means "
                                  "< 1/%d (no permutation out of %d reached the observed "
                                  "value)" % (args.reps, args.reps),
        },
        "scan_rule": {
            "eligible_window": "best interchromosomal parental identity >= min_window_identity "
                               "AND (no same-chromosome parental match OR interchromosomal identity "
                               "> same-chromosome identity)",
            "tract": "adjacent touching eligible windows, single donor chromosome, "
                     ">= min_run_bp, mean identity >= min_run_mean_identity",
        },
        "permutation_audit": {
            "same_chrom_winner_haplotype_counts": same_chrom_hap,
            "transmitted_copy_permutation_feasible": False,
            "why_not": "the class-winner table keeps only the per-window winner of the "
                       "same_chrom class, i.e. the maximum over the father's two haplotypes; "
                       "per-haplotype identities are not in this input, so 'which haplotype is "
                       "the transmitted copy' cannot be swapped. Label shuffling preserving "
                       "per-chromosome window counts is used instead (referee's stated alternative).",
        },
        "per_chromosome_marginals": {
            c: {"windows": len(pos_by_chrom[c]),
                "eligible_windows": len(eligible.get(c, []))}
            for c in chroms
        },
        "observed": {
            "n_qualifying_runs_genomewide": obs_total,
            "total_run_bp": sum(r["n"] for r in obs_runs) * WINDOW,
            "candidate_intervals": {
                chrom: {
                    "start": intervals[chrom][0],
                    "end": intervals[chrom][1],
                    "length_bp": intervals[chrom][1] - intervals[chrom][0],
                    "donor": candidates[chrom]["donor"],
                    "mean_identity": round(candidates[chrom]["mean_identity"], 6),
                    "n_windows": candidates[chrom]["n"],
                    "donor_haplotype_sequences": sorted(candidates[chrom]["seqs"]),
                    "homolog_absent_windows": candidates[chrom]["no_homolog_windows"],
                    "length_rank_among_qualifying_runs": length_rank[
                        (chrom, candidates[chrom]["start"])],
                    "eligible_windows_inside_interval": obs_interval_windows[chrom],
                    "n_qualifying_runs_overlapping": obs_interval_n[chrom],
                    "longest_overlapping_run_bp": obs_interval_max[chrom],
                }
                for chrom in MARKER_CHROMS
            },
            "top_runs_by_length": [
                {"chrom": r["chrom"], "start": r["start"], "end": r["end"],
                 "length_bp": r["end"] - r["start"], "donor": r["donor"],
                 "mean_identity": round(r["mean_identity"], 6),
                 "donor_haplotype_sequences": sorted(r["seqs"]),
                 "homolog_absent_windows": r["no_homolog_windows"]}
                for r in obs_runs[:20]
            ],
        },
        "null_A_permute_records_within_chromosome": {
            "description": "per-window record (eligible?, donor, identity) permuted within "
                           "each query chromosome; per-chromosome window counts, eligible-window "
                           "counts and donor/identity joint distribution preserved",
            "S1_total_qualifying_runs": {
                "observed": obs_total,
                "null": summarise(A_counts),
                "reps_reaching_observed": sum(1 for v in A_counts if v >= obs_total),
                "p_ge_observed": p_ge(obs_total, A_counts),
            },
            "S2_runs_overlapping_candidate": {
                name: {
                    "observed": obs_interval_n[name],
                    "null": summarise(A_win[name]),
                    "p_ge_observed": p_ge(obs_interval_n[name], A_win[name]),
                    "observed_longest_overlap_bp": obs_interval_max[name],
                    "null_longest_overlap_bp": summarise(A_maxlen[name]),
                }
                for name in intervals
            },
            "S3_eligible_windows_inside_candidate": {
                name: {
                    "observed": obs_interval_windows[name],
                    "null": summarise(A_nper[name]),
                    "p_ge_observed": p_ge(obs_interval_windows[name], A_nper[name]),
                }
                for name in intervals
            },
        },
        "null_B_shuffle_donor_labels_within_chromosome": {
            "description": "eligible/not status and identity fixed at each window; only the "
                           "donor chromosome label is shuffled among the eligible windows of "
                           "each query chromosome",
            "S1_total_qualifying_runs": {
                "observed": obs_total,
                "null": summarise(B_counts),
                "reps_reaching_observed": sum(1 for v in B_counts if v >= obs_total),
                "p_ge_observed": p_ge(obs_total, B_counts),
            },
            "S2_runs_overlapping_candidate": {
                name: {
                    "observed": obs_interval_n[name],
                    "null": summarise(B_win[name]),
                    "p_ge_observed": p_ge(obs_interval_n[name], B_win[name]),
                    "observed_longest_overlap_bp": obs_interval_max[name],
                    "null_longest_overlap_bp": summarise(B_maxlen[name]),
                }
                for name in intervals
            },
        },
        "limitations": [
            "A permutation null bounds the false-positive rate only under the null it "
            "enforces. It is not an orthogonal assay and cannot distinguish a real "
            "inter-chromosomal exchange from any process that produces the same "
            "window-level marginals (inherited paralogy, assembly artefact, mapping "
            "competition).",
            "Null A destroys spatial autocorrelation: the observed contiguity of eligible "
            "windows and of a single donor is exactly what it removes, so a count far above "
            "the null says the tract is not explained by window-level marginals alone, not "
            "that it is a real exchange.",
            "Null B conditions on the observed positions of eligible windows and only tests "
            "donor contiguity, so it inherits whatever produced those positions.",
            "Windows absent from the class-winner table (no alignment at all) are treated as "
            "absent windows and break runs; their absence is itself non-random and is "
            "preserved rather than permuted.",
            "estimated.identity in the class-winner table is the IMPG dice/coverage identity "
            "(equals dice.similarity to <=1e-6), not a base-level % identity. A 1% margin in "
            "this metric is not a 1% base-level divergence margin.",
            "All permutations are within one parent-child comparison (PAN027 paternal vs "
            "PAN011). Permuting which parental haplotype is the transmitted copy is not "
            "possible from this input, as documented under permutation_audit.",
        ],
    }

    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=False)
        fh.write("\n")

    # ---- stdout ----
    def fmt_p(k, n):
        return "<1/%d" % n if k == 0 else "%.4g" % (k / n)

    print("observed qualifying runs\t%d\t(%d bp)"
          % (obs_total, report["observed"]["total_run_bp"]))
    for name in intervals:
        c = report["observed"]["candidate_intervals"][name]
        print("candidate %s\t%d-%d\t%d bp\tdonor=%s\tmean_id=%.6f\truns_overlapping=%d"
              % (name, c["start"], c["end"], c["length_bp"], c["donor"],
                 c["mean_identity"], c["n_qualifying_runs_overlapping"]))
        k_run = sum(1 for v in A_win[name] if v >= obs_interval_n[name])
        k_win = sum(1 for v in A_nper[name] if v >= obs_interval_windows[name])
        k_runB = sum(1 for v in B_win[name] if v >= obs_interval_n[name])
        print("   null A runs_overlapping>=obs %d/%d p=%s | eligible_inside obs=%d null_mean=%.3f p=%s"
              % (k_run, args.reps, fmt_p(k_run, args.reps),
                 obs_interval_windows[name], statistics.fmean(A_nper[name]),
                 fmt_p(k_win, args.reps)))
        print("   null B runs_overlapping>=obs %d/%d p=%s (null_mean=%.3f null_max=%d)"
              % (k_runB, args.reps, fmt_p(k_runB, args.reps),
                 statistics.fmean(B_win[name]), max(B_win[name])))
    print("null A total runs: mean=%.4f sd=%.4f max=%d  reps>=obs=%d/%d p=%s"
          % (statistics.fmean(A_counts), statistics.stdev(A_counts), max(A_counts),
             sum(1 for v in A_counts if v >= obs_total), args.reps,
             fmt_p(sum(1 for v in A_counts if v >= obs_total), args.reps)))
    print("null B total runs: mean=%.4f sd=%.4f max=%d  reps>=obs=%d/%d p=%s"
          % (statistics.fmean(B_counts), statistics.stdev(B_counts), max(B_counts),
             sum(1 for v in B_counts if v >= obs_total), args.reps,
             fmt_p(sum(1 for v in B_counts if v >= obs_total), args.reps)))
    print("wrote\t%s" % args.out)


if __name__ == "__main__":
    sys.exit(main())
