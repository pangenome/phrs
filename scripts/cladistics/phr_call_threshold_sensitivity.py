#!/usr/bin/env python3
"""PHR-call threshold sensitivity -- Referee #3, major point 6.

The referee asks whether the claim that PHRs occur on 41 of 48 chromosome arms
depends on the PHR-detection thresholds (95% identity, >=5 alignments from >=2
chromosomes, 5-kb windows, 500-kb flank, 4 consecutive failed windows).

This script re-applies the exact window-walking logic of the calling code
(scripts/phr-detection/find-multichr-regions-incremental.py,
process_sequence()) to the per-window IMPG query table that produced the
published PHR call, varying every threshold that the per-window table actually
lets us vary, and reports the resulting arm count.

WHAT IS VARIABLE HERE (re-derived, no re-alignment needed)
  --min-diff-chrs   (distinct chromosomes required per window)  <- per-window n_chrs
  --min-consecutive-stop (consecutive failing windows)          <- window order
  --max             (inward scan limit from the telomere)        <- dist_from_telomere

Note on --min-consecutive-stop: the published scan stopped at the 4th consecutive
failing window, so the per-window table contains nothing beyond that point.
Stricter values (1, 2, 3) are computable; *looser* values (>4) are not, because
the windows they would need were never queried.

WHAT IS NOT VARIABLE HERE (baked into the per-window table; varying them needs a
re-run of the wfmash all-vs-all + a rebuild of the IMPG index. The wfmash and
impg binaries are available; the blocker is that the 18,827 all-vs-all PAFs are
88 GB and the run's paf.list points at /scratch paths that no longer exist, so a
full re-index is out of scope for this item)
  --min-identity 0.95   wfmash -p 95 produced no alignments below 95% identity,
                        and `impg query --min-identity 0.95` was applied before
                        the per-chromosome counts were taken.
  --min-count 5         the per-window table stores *chromosomes with >=5
                        alignments* (columns n_chrs / chrs_involved), not the raw
                        per-chromosome alignment counts, so the 5 can only be
                        raised/lowered by re-querying IMPG.
  --window 5000, --step 5000, --min-output-length 3000
                        window geometry and the 3-kb projection minimum are fixed
                        by the run that wrote the table.

Inputs (all read-only):
  --per-window   per-window IMPG query table of the published call
                 (seq, arm, self_chr, window_start, window_end,
                  dist_from_telomere, n_chrs, n_arms, chrs_involved,
                  arms_involved, passes_threshold)
                 Default: the moosefs copy, because the repo does not vendor it
                 (data/all-vs-all.1Mb.p95.id95.len.per_window.tsv does not exist).
  --flank-summary  data/all-vs-all.1Mb.p95.id95.len.tsv -- the published
                 per-flank call (one row per flank, region_start/region_end/
                 chrs_involved/arms_involved). Used as the replication gate and
                 as the PUBLISHED row of the output table.

Output: scripts/cladistics/phr_call_threshold_sensitivity.tsv

Gate (fail loud): re-walking the per-window table at the published settings must
reproduce, for every one of the 18,826 flanks, the published region_start,
region_end and the published arm set, and must reproduce 15,668 PHRs on 41 arms.
If it does not, the script prints the first few mismatches and exits non-zero.
"""
import argparse
import csv
import os
import statistics
import sys

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
ALL_ARMS = tuple(f"chr{c}{a}" for c in CHROMS for a in ("p", "q"))  # 48 arms

# Arms removed by the "no-acrocentric" variant of the call. This mirrors the
# actual content of data/chm13.phrs.no_acro.bed: relative to
# data/chm13.phrs.bed it drops every interval on chr13, chr14, chr15, chr21 and
# chr22 (both arms), keeping the PAR-bearing chrX intervals.
ACRO_ARMS = frozenset(
    f"chr{c}{a}" for c in ("13", "14", "15", "21", "22") for a in ("p", "q")
)
NON_ACRO_ARMS = tuple(a for a in ALL_ARMS if a not in ACRO_ARMS)  # 38 arms

# Published run parameters (scripts/phr-detection/PIPELINE.md, Stage 4).
PUB_MIN_CHRS = 2
PUB_CONSEC = 4
PUB_MAX_BP = 500_000
PUB_MIN_COUNT = 5
PUB_MIN_IDENTITY = 0.95
PUB_WINDOW = 5_000
PUB_MIN_OUTPUT_LEN = 3_000

# Flank dropped before PHR calling (scaffolding chimera, paper Methods and
# PIPELINE.md Stage 4); it is present in the per-window table but must be
# ignored so the replication matches the published 18,826-flank table.
# Filled in from the published table before the sweep; the reference arm set for
# the gained/lost columns.
PUBLISHED_ARMS = frozenset()

CHIMERA_SUBSTR = "NA18982#1#JBKABS010000018.1:83848607-84348606_chr18_qarm"

DEFAULT_PER_WINDOW = (
    "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.per_window.tsv"
)


def arm_of(seq_name):
    """chr10:134254995-134754994_chr10_qarm -> chr10q (same rule as the caller)."""
    parts = seq_name.rsplit("_", 2)
    if len(parts) >= 3 and parts[-1] in ("parm", "qarm") and parts[-2].startswith("chr"):
        return parts[-2] + ("p" if parts[-1] == "parm" else "q")
    raise SystemExit(f"cannot parse chromosome arm from sequence name: {seq_name}")


class Window:
    __slots__ = ("start", "end", "dist", "n_chrs", "arms")

    def __init__(self, start, end, dist, n_chrs, arms):
        self.start = start
        self.end = end
        self.dist = dist
        self.n_chrs = n_chrs
        self.arms = arms  # frozenset of chrNp/chrNq labels, empty if none


def read_per_window(path, window):
    """Return dict seq -> (arm, [Window, ...]) in scan order (telomere inward)."""
    if not os.path.exists(path):
        raise SystemExit(
            f"per-window table not found: {path}\n"
            "It is not vendored in the repo; re-fetch it with\n"
            "  scp HOST:/moosefs/guarracino/HPRCv2/PHR_III/"
            "all-vs-all.1Mb.p95.id95.len.per_window.tsv .\n"
            "and pass --per-window."
        )
    seqs = {}
    order = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        expect = ["seq", "arm", "self_chr", "window_start", "window_end",
                  "dist_from_telomere", "n_chrs", "n_arms", "chrs_involved",
                  "arms_involved", "passes_threshold"]
        if header != expect:
            raise SystemExit(f"unexpected per-window header:\n{header}\n!=\n{expect}")
        current = None
        for row in rd:
            seq = row[0]
            if seq != current:
                if seq in seqs:
                    raise SystemExit(
                        f"per-window rows for {seq} are not contiguous; "
                        "cannot re-walk the scan order"
                    )
                current = seq
                seqs[seq] = (arm_of(seq), [])
                order.append(seq)
            arms_field = row[9]
            arms = frozenset() if arms_field in ("", ".") else frozenset(arms_field.split(","))
            seqs[seq][1].append(
                Window(int(row[3]), int(row[4]), int(row[5]), int(row[6]), arms)
            )
    return seqs, order


def read_flank_summary(path):
    """seq -> (region_start, region_end, arms frozenset) from the published call."""
    out = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        expect = ["seq", "arm", "self_chr", "region_start", "region_end",
                  "chrs_involved", "arms_involved"]
        if header != expect:
            raise SystemExit(f"unexpected flank-summary header:\n{header}\n!=\n{expect}")
        for row in rd:
            if row[3] == ".":
                out[row[0]] = None
            else:
                arms = row[6]
                arms = frozenset() if arms in ("", ".") else frozenset(arms.split(","))
                out[row[0]] = (int(row[3]), int(row[4]), arms)
    return out


def walk(windows, min_chrs, consec_stop, max_bp, window):
    """Re-run process_sequence()'s region assembly on one flank.

    Returns (region_start, region_end, arms, n_windows_in_region) or None.
    Faithful to the calling code: windows are visited in scan order (telomere
    inward); a *passing* window flushes the pending failing windows into the
    region, and the scan stops at the consec_stop-th consecutive failing window
    (the failing run that triggers the stop is itself not part of any region).
    """
    region_idx = []
    pending = []
    consecutive_fail = 0
    for idx, w in enumerate(windows):
        if w.dist + window > max_bp:
            break  # beyond the scan limit: never queried for this setting
        if w.n_chrs >= min_chrs:
            region_idx.extend(pending)
            pending = []
            consecutive_fail = 0
            region_idx.append(idx)
        else:
            consecutive_fail += 1
            pending.append(idx)
            if consecutive_fail >= consec_stop:
                break
    if not region_idx:
        return None
    rw = [windows[i] for i in region_idx]
    arms = frozenset().union(*(w.arms for w in rw))
    return min(w.start for w in rw), max(w.end for w in rw), arms, len(rw)


def compute(seqs, order, min_chrs, consec_stop, max_bp, window):
    """Apply one setting to every flank. Returns dict seq -> (start, end, arms, nwin)."""
    out = {}
    for seq in order:
        arm, windows = seqs[seq]
        res = walk(windows, min_chrs, consec_stop, max_bp, window)
        out[seq] = (None if res is None else res, arm)
    return out


def setting_row(calls, flanks_total, label, param, value, min_chrs, consec_stop, max_bp):
    n_phrs = 0
    bearing = set()
    partner = set()
    lengths = []
    for seq, (res, own_arm) in calls.items():
        if res is None:
            continue
        n_phrs += 1
        bearing.add(own_arm)
        start, end, arms, _ = res
        lengths.append(end - start)
        partner.update(arms)
    non_acro = sorted(a for a in bearing if a not in ACRO_ARMS)
    return {
        "variant": label,
        "param_varied": param,
        "value": value,
        "min_diff_chrs": min_chrs,
        "min_consecutive_stop": consec_stop,
        "max_flank_bp": max_bp,
        "min_count_alignments_per_chr": PUB_MIN_COUNT,
        "min_identity": PUB_MIN_IDENTITY,
        "window_bp": PUB_WINDOW,
        "min_output_len_bp": PUB_MIN_OUTPUT_LEN,
        "n_flanks": flanks_total,
        "n_phrs": n_phrs,
        "pct_flanks_with_phr": round(100.0 * n_phrs / flanks_total, 2),
        "n_arms_bearing_phr_of_48": len(bearing),
        "arms_bearing_phr": ",".join(sorted(bearing, key=arm_sort_key)),
        "n_arms_bearing_phr_non_acro_of_38": len(non_acro),
        "arms_gained_vs_published": ",".join(
            sorted(bearing - PUBLISHED_ARMS, key=arm_sort_key)),
        "arms_lost_vs_published": ",".join(
            sorted(PUBLISHED_ARMS - bearing, key=arm_sort_key)),
        "n_arms_as_partner_of_48": len(partner),
        "median_phr_len_bp": int(statistics.median(lengths)) if lengths else 0,
        "mean_phr_len_bp": int(round(statistics.fmean(lengths))) if lengths else 0,
        "total_mb_bearing": round(sum(lengths) / 1e6, 2),
    }


def arm_sort_key(a):
    chr_part = a[3:-1]
    order = {c: i for i, c in enumerate(CHROMS)}
    return (order.get(chr_part, 99), a[-1])


COLUMNS = [
    "variant", "param_varied", "value", "min_diff_chrs", "min_consecutive_stop",
    "max_flank_bp", "min_count_alignments_per_chr", "min_identity", "window_bp",
    "min_output_len_bp", "n_flanks", "n_phrs", "pct_flanks_with_phr",
    "n_arms_bearing_phr_of_48", "arms_bearing_phr",
    "n_arms_bearing_phr_non_acro_of_38", "arms_gained_vs_published",
    "arms_lost_vs_published", "n_arms_as_partner_of_48",
    "median_phr_len_bp", "mean_phr_len_bp", "total_mb_bearing",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-window", default=DEFAULT_PER_WINDOW)
    ap.add_argument("--flank-summary",
                    default=os.path.join("data", "all-vs-all.1Mb.p95.id95.len.tsv"))
    ap.add_argument("--out", default=os.path.join(
        "scripts", "cladistics", "phr_call_threshold_sensitivity.tsv"))
    ap.add_argument("--window", type=int, default=PUB_WINDOW)
    args = ap.parse_args()

    print(f"reading per-window table {args.per_window}", file=sys.stderr)
    all_seqs, all_order = read_per_window(args.per_window, args.window)
    print(f"  {len(all_seqs)} flanks, {sum(len(w) for _, w in all_seqs.values())} windows",
          file=sys.stderr)

    published = read_flank_summary(args.flank_summary)
    print(f"reading published per-flank call {args.flank_summary}: "
          f"{len(published)} flanks", file=sys.stderr)

    # Drop the chimera so the walk covers exactly the published flank set.
    if CHIMERA_SUBSTR not in all_seqs:
        raise SystemExit(f"expected chimera flank not found in per-window table: "
                         f"{CHIMERA_SUBSTR}")
    seqs = {s: v for s, v in all_seqs.items() if s != CHIMERA_SUBSTR}
    order = [s for s in all_order if s != CHIMERA_SUBSTR]
    if set(seqs) != set(published):
        missing = sorted(set(published) - set(seqs))[:5]
        extra = sorted(set(seqs) - set(published))[:5]
        raise SystemExit(f"flank sets differ between per-window and summary table; "
                         f"only in summary: {missing}; only in per-window: {extra}")

    # ---------------- replication gate ----------------
    base = compute(seqs, order, PUB_MIN_CHRS, PUB_CONSEC, PUB_MAX_BP, args.window)
    bad = []
    for seq in order:
        (res, _), pub = base[seq], published[seq]
        if (res is None) != (pub is None):
            bad.append((seq, "presence", res is not None, pub is not None))
            continue
        if res is None:
            continue
        start, end, arms, _ = res
        if (start, end) != (pub[0], pub[1]):
            bad.append((seq, "coords", (start, end), (pub[0], pub[1])))
        elif arms != pub[2]:
            bad.append((seq, "arms", sorted(arms), sorted(pub[2])))
    if bad:
        for b in bad[:10]:
            print("MISMATCH", b, file=sys.stderr)
        raise SystemExit(
            f"replication FAILED: {len(bad)} of {len(order)} flanks disagree with "
            f"{args.flank_summary} at the published settings"
        )
    n_phrs = sum(1 for r, _ in base.values() if r is not None)
    arms = {arm_of(s) for s, (r, _) in base.items() if r is not None}
    print(f"replication gate PASSED: {n_phrs} PHRs on {len(arms)} arms, all "
          f"{len(order)} region_start/region_end/arm sets identical to the "
          f"published table", file=sys.stderr)

    # ---------------- sweep ----------------
    rows = []

    # Published row straight from the repo-local flank table (not re-derived).
    pub_calls = {}
    for seq, pub in published.items():
        own = arm_of(seq)
        if pub is None:
            pub_calls[seq] = (None, own)
        else:
            pub_calls[seq] = ((pub[0], pub[1], pub[2], 0), own)
    global PUBLISHED_ARMS
    PUBLISHED_ARMS = frozenset(own for (res, own) in pub_calls.values() if res is not None)
    rows.append(setting_row(pub_calls, len(published), "PUBLISHED_TABLE",
                            "none (data/all-vs-all.1Mb.p95.id95.len.tsv)", "",
                            PUB_MIN_CHRS, PUB_CONSEC, PUB_MAX_BP))

    rows.append(setting_row(base, len(seqs), "REPLICATION_FROM_PER_WINDOW",
                            "none (baseline re-walk)", "", PUB_MIN_CHRS,
                            PUB_CONSEC, PUB_MAX_BP))

    for k in (1, 3, 4, 5, 6, 8, 10, 12, 15):
        calls = compute(seqs, order, k, PUB_CONSEC, PUB_MAX_BP, args.window)
        rows.append(setting_row(calls, len(seqs), "SWEEP_MIN_DIFF_CHRS",
                                "min chromosomes with >=5 alignments per window",
                                k, k, PUB_CONSEC, PUB_MAX_BP))
        print(f"  min_diff_chrs={k}: done", file=sys.stderr)

    for c in (1, 2, 3):
        calls = compute(seqs, order, PUB_MIN_CHRS, c, PUB_MAX_BP, args.window)
        rows.append(setting_row(calls, len(seqs), "SWEEP_CONSECUTIVE_STOP",
                                "consecutive failing windows that end the scan",
                                c, PUB_MIN_CHRS, c, PUB_MAX_BP))
        print(f"  consecutive_stop={c}: done", file=sys.stderr)

    for m in (25_000, 50_000, 100_000, 150_000, 250_000, 350_000, 400_000):
        calls = compute(seqs, order, PUB_MIN_CHRS, PUB_CONSEC, m, args.window)
        rows.append(setting_row(calls, len(seqs), "SWEEP_FLANK_LIMIT",
                                "inward scan limit from the telomere (bp)",
                                m, PUB_MIN_CHRS, PUB_CONSEC, m))
        print(f"  max_flank={m}: done", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                            lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    print(f"wrote {args.out} ({len(rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
