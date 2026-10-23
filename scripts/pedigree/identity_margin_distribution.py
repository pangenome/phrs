#!/usr/bin/env python3
"""
Genome-wide identity-margin distribution for the WashU pedigree subtelomere scan.

Referee #1, major point 2, asked for a "genome-wide distribution of identity
margins locating the candidates".  The marginal logic the referee criticised is
(POINT_BY_POINT R1-M2): the scan is satisfied whenever a transmitted allele
resembles a paralogous chromosome end more than the father's own copy of the
same end.  We therefore define, for every 2 kb window of the PAN027 paternal
(PAN027#2 = SELF.p) query against the father (PAN011, both haplotypes):

    margin = identity(best non-homologous / interchromosomal parental match)
             - identity(best same-chromosome / homologous parental match)

Both terms are the per-window class winners already used by the Fig. 5 pipeline
(columns `winner_class` and `estimated.identity` of the IMPG class-winner
table).  `estimated.identity` in this table equals `dice.similarity` to <=1e-6
(verified below): it is the IMPG window-level estimated identity from
alignment-block coverage, NOT a base-level % identity recomputed from raw
alignment.  Every number printed here comes from the input table; nothing is
re-estimated.

OUTPUTS

  scripts/pedigree/identity_margin_distribution.tsv
      tidy 6-column TSV (see SCHEMA below); sections are distinguished by
      `record_type`.  No per-window dump of the whole genome is written.

SCHEMA (tab separated, one header line, no comment lines)

  record_type   label   value   value2   value3   note

  record_type in {summary, quantile, hist, candidate_run, candidate_window,
                  candidate_context}

USAGE

  python3 scripts/pedigree/identity_margin_distribution.py

  # explicit paths / different marker chromosomes
  python3 scripts/pedigree/identity_margin_distribution.py \
      --similarity data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz \
      --fai        data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai \
      --marker-chrom PAN027#2#chr9 --marker-chrom PAN027#2#chr5 \
      --out scripts/pedigree/identity_margin_distribution.tsv
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

WINDOW = 2000  # bp; the scan's window size (all rows are 2000 bp windows)

# Candidate (marker) definition -- fully data-derived, no hard-coded coordinates:
# for each marker query chromosome take the LONGEST run that satisfies the
# scan's own run filter (adjacent touching same-donor windows, >=10 kb, mean
# identity >=0.95, positive margin) and lies inside the terminal 100 kb.
MARKER_CHROMS = ["PAN027#2#chr9", "PAN027#2#chr5"]
TERMINAL_BP = 100_000
MIN_RUN_WINDOWS = 5          # >=10 kb
MIN_RUN_MEAN_IDENTITY = 0.95
MIN_WINDOW_IDENTITY = 0.95   # a window must reach this to be run-eligible

# margin histogram edges (half-open bins; exact zeros get their own bin)
NEG_EDGES = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2,
             -0.1, -0.05, -0.02, -0.01, -0.005, -0.002, -0.001]
POS_EDGES = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

QUANTILES = [0.0, 0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75,
             0.90, 0.95, 0.975, 0.99, 0.995, 0.999, 1.0]


# --------------------------------------------------------------------------- #
# input
# --------------------------------------------------------------------------- #
def load_windows(similarity_path):
    """Return {(chrom, start): {'chrom','start','end','interchrom','same_chrom'}}.

    Each class entry is a dict with identity / other_seq / other_chrom.  The
    class-winner table has at most one row per (window, winner_class), which is
    asserted rather than assumed.
    """
    windows = {}
    with gzip.open(similarity_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for required in ("chrom", "start", "end", "estimated.identity",
                         "winner_class", "other_seq", "other_chrom",
                         "dice.similarity", "intersection", "group.b.length"):
            if required not in ix:
                sys.exit("ERROR: input is missing column %r" % required)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            key = (p[ix["chrom"]], int(p[ix["start"]]))
            rec = windows.get(key)
            if rec is None:
                rec = windows[key] = {"chrom": key[0], "start": key[1],
                                      "end": int(p[ix["end"]])}
            cls = p[ix["winner_class"]]
            if cls in rec:
                sys.exit("ERROR: duplicate %s class row for window %s:%d"
                         % (cls, key[0], key[1]))
            rec[cls] = {"identity": float(p[ix["estimated.identity"]]),
                        "dice": float(p[ix["dice.similarity"]]),
                        "other_seq": p[ix["other_seq"]],
                        "other_chrom": p[ix["other_chrom"]],
                        "intersection": int(p[ix["intersection"]]),
                        "query_len": int(p[ix["group.b.length"]])}
    if not windows:
        sys.exit("ERROR: no windows read from %s" % similarity_path)
    return windows


def load_fai(fai_path):
    """Return {contig: length} in file order."""
    out = {}
    with open(fai_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            out[p[0]] = int(p[1])
    return out


# --------------------------------------------------------------------------- #
# per-window classification
# --------------------------------------------------------------------------- #
def margin(rec):
    """margin = best interchrom identity - best same_chrom identity.

    Returns None when either class is missing (see `state`).
    """
    i = rec.get("interchrom")
    s = rec.get("same_chrom")
    if i is None or s is None:
        return None
    return i["identity"] - s["identity"]


def state(rec):
    """'both' | 'interchrom_only' | 'same_chrom_only'."""
    if "interchrom" in rec and "same_chrom" in rec:
        return "both"
    if "interchrom" in rec:
        return "interchrom_only"
    return "same_chrom_only"


def is_run_eligible(rec):
    """A window is eligible for a run iff the interchromosomal match is the best
    parental match in that window (margin > 0), OR there is no homologous
    parental match at all in that window (`interchrom_only`: the homolog is
    absent, which is the strongest form of the same statement), AND the
    interchromosomal identity reaches MIN_WINDOW_IDENTITY.

    Returns the donor chromosome, or None.
    """
    i = rec.get("interchrom")
    if i is None or i["identity"] < MIN_WINDOW_IDENTITY:
        return None
    s = rec.get("same_chrom")
    if s is not None and i["identity"] <= s["identity"]:
        return None
    return i["other_chrom"]


def find_runs(windows):
    """Runs of adjacent (touching) eligible windows with a single donor
    chromosome.  Returns a list of dicts, unfiltered; `qualifying()` applies the
    >=10 kb / mean-identity >=0.95 filters."""
    by_chrom = defaultdict(list)
    for rec in windows.values():
        donor = is_run_eligible(rec)
        if donor is None:
            continue
        by_chrom[rec["chrom"]].append((rec["start"], rec["end"], donor,
                                       rec["interchrom"]["identity"],
                                       rec["interchrom"]["other_seq"],
                                       "same_chrom" not in rec))
    runs = []
    for chrom, rows in by_chrom.items():
        rows.sort()
        cur = None
        for start, end, donor, ident, seq, no_homolog in rows:
            if cur is not None and cur["end"] == start and cur["donor"] == donor:
                cur["end"] = end
                cur["n"] += 1
                cur["identities"].append(ident)
                cur["seqs"].add(seq)
                cur["no_homolog_windows"] += 1 if no_homolog else 0
            else:
                if cur is not None:
                    runs.append(cur)
                cur = {"chrom": chrom, "start": start, "end": end,
                       "donor": donor, "n": 1, "identities": [ident],
                       "seqs": {seq},
                       "no_homolog_windows": 1 if no_homolog else 0}
        if cur is not None:
            runs.append(cur)
    for r in runs:
        r["mean_identity"] = sum(r["identities"]) / r["n"]
        r["touching"] = r["n"] > 1
    return runs


def qualifying(runs):
    return [r for r in runs
            if r["n"] >= MIN_RUN_WINDOWS
            and r["mean_identity"] >= MIN_RUN_MEAN_IDENTITY]


# --------------------------------------------------------------------------- #
# candidate derivation
# --------------------------------------------------------------------------- #
def derive_candidates(windows, chrom_lengths, marker_chroms, terminal_bp=TERMINAL_BP):
    """The marker candidate on each marker chromosome: the longest qualifying
    run lying entirely inside the terminal `terminal_bp` of that chromosome."""
    runs = qualifying(find_runs(windows))
    out = {}
    for chrom in marker_chroms:
        if chrom not in chrom_lengths:
            sys.exit("ERROR: marker chromosome %r not in the .fai" % chrom)
        lo = chrom_lengths[chrom] - terminal_bp
        inside = [r for r in runs if r["chrom"] == chrom and r["end"] > lo
                  and r["start"] >= lo]
        if not inside:
            sys.exit("ERROR: no qualifying run in the terminal %d bp of %s"
                     % (terminal_bp, chrom))
        inside.sort(key=lambda r: (-(r["end"] - r["start"]), r["start"]))
        out[chrom] = inside[0]
    return out


# --------------------------------------------------------------------------- #
# histogram / quantiles
# --------------------------------------------------------------------------- #
def bin_label(v):
    if v == 0.0:
        return "=0"
    if v < 0:
        for lo, hi in zip(NEG_EDGES[:-1], NEG_EDGES[1:]):
            if lo <= v < hi:
                return "[%.4f,%.4f)" % (lo, hi)
        return "[%.4f,%.4f)" % (NEG_EDGES[0], NEG_EDGES[1])  # v == -1.0 edge
    for lo, hi in zip(POS_EDGES[:-1], POS_EDGES[1:]):
        if lo < v < hi:
            return "(%.4f,%.4f)" % (lo, hi)
    return "(%.4f,%.4f]" % (POS_EDGES[-2], POS_EDGES[-1])  # v == 1.0 edge


def quantile(sorted_vals, q):
    if not sorted_vals:
        sys.exit("ERROR: empty distribution")
    # nearest-rank (type 1): smallest value with CDF >= q
    n = len(sorted_vals)
    if q <= 0.0:
        return sorted_vals[0]
    if q >= 1.0:
        return sorted_vals[-1]
    idx = 0
    import math
    idx = max(0, min(n - 1, int(math.ceil(q * n)) - 1))
    return sorted_vals[idx]


def percentile_of(sorted_vals, v):
    """Percent of windows with a margin strictly below v (0..100)."""
    import bisect
    return 100.0 * bisect.bisect_left(sorted_vals, v) / len(sorted_vals)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--similarity",
                    default="data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz")
    ap.add_argument("--fai",
                    default="data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai")
    ap.add_argument("--out", default="scripts/pedigree/identity_margin_distribution.tsv")
    ap.add_argument("--marker-chrom", action="append", dest="markers",
                    default=None,
                    help="query chromosome to locate on the distribution "
                         "(repeatable; default chr9 + chr5)")
    args = ap.parse_args()

    windows = load_windows(args.similarity)
    chrom_lengths = load_fai(args.fai)
    markers = args.markers or MARKER_CHROMS

    # ---- audit the identity column semantics -----------------------------
    max_dice_delta = 0.0
    for rec in windows.values():
        for cls in ("interchrom", "same_chrom"):
            if cls in rec:
                d = abs(rec[cls]["identity"] - rec[cls]["dice"])
                if d > max_dice_delta:
                    max_dice_delta = d

    # ---- states and margins ----------------------------------------------
    counts = {"both": 0, "interchrom_only": 0, "same_chrom_only": 0}
    margins = []
    ties = 0
    for rec in windows.values():
        st = state(rec)
        counts[st] += 1
        if st == "both":
            m = margin(rec)
            margins.append(m)
            if m == 0.0:
                ties += 1
    margins.sort()
    n_margin = len(margins)
    mean_m = sum(margins) / n_margin
    var_m = sum((m - mean_m) ** 2 for m in margins) / (n_margin - 1)
    sd_m = var_m ** 0.5

    runs_all = find_runs(windows)
    runs_ok = qualifying(runs_all)
    runs_ok_sorted = sorted(runs_ok, key=lambda r: (-(r["end"] - r["start"]),
                                                    r["chrom"], r["start"]))
    rank = {}
    for i, r in enumerate(runs_ok_sorted, 1):
        rank[(r["chrom"], r["start"])] = i

    candidates = derive_candidates(windows, chrom_lengths, markers)

    rows = []

    def add(rec_type, label, value, value2="", value3="", note=""):
        rows.append([rec_type, str(label), str(value), str(value2),
                     str(value3), note])

    # ---------------- summary ----------------
    add("summary", "input_similarity_table", args.similarity, "", "",
        "IMPG class-winner table, PAN027 paternal (PAN027#2 = SELF.p) vs PAN011 (father)")
    add("summary", "input_query_fai", args.fai, "", "",
        "query (recipient) contig lengths")
    add("summary", "window_size_bp", WINDOW, "", "", "all scanned windows")
    add("summary", "windows_total", len(windows), "", "",
        "distinct (chrom,start) windows present in the table")
    add("summary", "windows_with_both_classes", counts["both"], "", "",
        "margin defined")
    add("summary", "windows_interchrom_only", counts["interchrom_only"], "", "",
        "no homologous parental match in the window (margin undefined; treated as "
        "the strongest positive statement in run detection)")
    add("summary", "windows_same_chrom_only", counts["same_chrom_only"], "", "",
        "no non-homologous parental match (margin undefined; run-ineligible)")
    add("summary", "max_abs_delta_estimated_identity_minus_dice", "%.3g" % max_dice_delta,
        "", "", "confirms estimated.identity is the IMPG dice/coverage identity, "
                "not a base-level % identity")
    add("summary", "margin_defined_n", n_margin, "", "", "")
    add("summary", "margin_mean", "%.6f" % mean_m, "", "", "")
    add("summary", "margin_sd", "%.6f" % sd_m, "", "", "sample sd (n-1)")
    add("summary", "margin_min", "%.6f" % margins[0], "", "", "")
    add("summary", "margin_max", "%.6f" % margins[-1], "", "", "")
    add("summary", "margin_exact_ties", ties, "", "",
        "windows where the paralogous and homologous best matches are exactly equal")
    for thr in (0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1):
        n = sum(1 for m in margins if m > thr)
        add("summary", "windows_margin_gt_%g" % thr, n, "",
            "%.4f%%" % (100.0 * n / n_margin),
            "count and percent of margin-defined windows")
    add("summary", "qualifying_runs_genomewide", len(runs_ok), "", "",
        "adjacent touching windows, same donor chromosome, >=%d kb, mean identity "
        ">=%.2f, every window with margin>0 or homolog-absent and identity>=%.2f"
        % (WINDOW * MIN_RUN_WINDOWS // 1000, MIN_RUN_MEAN_IDENTITY, MIN_WINDOW_IDENTITY))
    add("summary", "qualifying_runs_total_bp", sum(r["n"] for r in runs_ok) * WINDOW,
        "", "", "sum of run spans (windows x %d bp)" % WINDOW)

    # ---------------- quantiles ----------------
    for q in QUANTILES:
        add("quantile", "margin_p%g" % (q * 100.0), "%.6f" % quantile(margins, q),
            "", "", "nearest-rank quantile over the %d margin-defined windows" % n_margin)

    # ---------------- histogram ----------------
    hist = defaultdict(int)
    for m in margins:
        hist[bin_label(m)] += 1
    order = ["=0"]
    for lo, hi in zip(NEG_EDGES[:-1], NEG_EDGES[1:]):
        order.append("[%.4f,%.4f)" % (lo, hi))
    for lo, hi in zip(POS_EDGES[:-1], POS_EDGES[1:]):
        order.append("(%.4f,%.4f)" % (lo, hi))
    seen = set()
    for label in order:
        if label in seen:
            continue
        seen.add(label)
        add("hist", label, hist.get(label, 0), "",
            "%.4f%%" % (100.0 * hist.get(label, 0) / n_margin), "margin bin")
    for label in sorted(hist):
        if label not in seen:
            add("hist", label, hist[label], "",
                "%.4f%%" % (100.0 * hist[label] / n_margin), "margin bin")

    # ---------------- candidate markers ----------------
    for chrom in markers:
        r = candidates[chrom]
        add("candidate_run", chrom, r["start"], r["end"], r["n"] * WINDOW,
            "donor=%s; mean_identity=%.6f; donor_seqs=%s; no_homolog_windows=%d; "
            "length_rank_among_%d_qualifying_runs=%d"
            % (r["donor"], r["mean_identity"], ",".join(sorted(r["seqs"])),
               r["no_homolog_windows"], len(runs_ok),
               rank[(r["chrom"], r["start"])]))
        # per-window margins for this candidate
        n_pos = 0
        for start in range(r["start"], r["end"], WINDOW):
            rec = windows.get((chrom, start))
            if rec is None:
                sys.exit("ERROR: candidate window %s:%d missing" % (chrom, start))
            m = margin(rec)
            i = rec.get("interchrom")
            s = rec.get("same_chrom")
            if m is None:
                add("candidate_window", "%s:%d-%d" % (chrom, start, rec["end"]),
                    "NA", "%.6f" % i["identity"] if i else "NA",
                    "%.6f" % s["identity"] if s else "ABSENT",
                    "margin undefined; state=%s" % state(rec))
                continue
            if m > 0:
                n_pos += 1
            add("candidate_window", "%s:%d-%d" % (chrom, start, rec["end"]),
                "%.6f" % m, "%.6f" % i["identity"], "%.6f" % s["identity"],
                "percentile_of_margin=%.3f" % percentile_of(margins, m))
        add("candidate_context", chrom, r["start"], r["end"],
            "positive_margin_windows=%d/%d" % (n_pos, r["n"]),
            "fraction of candidate windows whose paralogous match beats the "
            "homologous match")

    # telomeric-60kb context: the block quoted in response/01_evidence_log_pedigree.txt
    for chrom in markers:
        if chrom not in chrom_lengths:
            continue
        lo = chrom_lengths[chrom] - 60_000
        pos = [(rec["start"], margin(rec)) for (c, s), rec in windows.items()
               if c == chrom and s >= lo and margin(rec) is not None
               and margin(rec) > 0]
        pos.sort()
        if not pos:
            add("candidate_context", chrom, lo, chrom_lengths[chrom], 0,
                "terminal 60 kb: no window with margin>0")
            continue
        pcts = [percentile_of(margins, m) for _s, m in pos]
        add("candidate_context", chrom, pos[0][0], pos[-1][0] + WINDOW, len(pos),
            "terminal 60 kb (%d-%d): windows_with_margin_gt_0=%d; "
            "margins=%s; percentile_range=%.3f..%.3f"
            % (lo, chrom_lengths[chrom], len(pos),
               ",".join("%.6f" % m for _s, m in pos),
               min(pcts), max(pcts)))

    with open(args.out, "w") as out:
        out.write("record_type\tlabel\tvalue\tvalue2\tvalue3\tnote\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    # ---- stdout summary (numbers only, no re-derivation) ----
    print("windows_total\t%d" % len(windows))
    print("margin_defined\t%d" % n_margin)
    print("interchrom_only\t%d" % counts["interchrom_only"])
    print("margin>0\t%d (%.4f%%)" % (sum(1 for m in margins if m > 0),
                                    100.0 * sum(1 for m in margins if m > 0) / n_margin))
    print("median\t%.6f\tp95\t%.6f\tp99\t%.6f"
          % (quantile(margins, 0.5), quantile(margins, 0.95), quantile(margins, 0.99)))
    print("qualifying_runs\t%d" % len(runs_ok))
    for chrom in markers:
        r = candidates[chrom]
        print("candidate\t%s\t%d\t%d\t%d bp\tdonor=%s\tmean_id=%.6f\trank=%d/%d"
              % (chrom, r["start"], r["end"], r["n"] * WINDOW, r["donor"],
                 r["mean_identity"], rank[(r["chrom"], r["start"])], len(runs_ok)))
    print("wrote\t%s" % args.out)


if __name__ == "__main__":
    main()
