#!/usr/bin/env python3
"""R1-minor-11: descriptive summary of PHRs and communities by HPRC superpopulation.

Verbatim ask (PEER_REVIEW.md, Referee #1, minor point 11): "The 232 individuals span
five HPRC superpopulations, yet no population-stratified analysis was performed. Given
the 'population-scale' framing, and that subtelomeric duplicon copy number varies
between populations, even a descriptive summary of PHR length and community membership
by superpopulation would strengthen the claim."

This script is DESCRIPTIVE. It does not test population differentiation; the formal
differentiation analysis is the Hudson F_ST work in scripts/popgen/matched_fst_d_m6.py.

WHAT IT COMPUTES, per superpopulation (AFR, AMR, EAS, EUR, SAS):
  * assemblies (individuals) and haplotype-resolved assemblies (haplotypes) observed
  * number of arms carrying at least one called PHR
  * number of called PHRs, and total PHR bp
  * PHR length distribution: median and IQR (linear-interpolation type-7 quantiles)
  * PHRs per haplotype: median and IQR
  * community membership distribution: PHR count and bp per arm-level Leiden community,
    as a fraction of that superpopulation's PHRs and PHR bp
  * PHR length distribution deciles p10 and p90, and the fraction of PHRs whose length
    equals the full 500 kb flank (i.e. right-censored at the flank boundary)
  * PHR abundance normalised for haplotype count: (this superpopulation's share of
    cohort PHRs) / (its share of cohort haplotypes).  1.0 = exactly proportional.
  * two label-permutation checks (superpopulation labels permuted among haplotypes, so
    within-haplotype clustering of community membership is preserved):
      - community composition (Cramer's V on the superpopulation x community bp table)
      - per-haplotype median PHR length (Kruskal-Wallis H)
    Both are reported as descriptive heterogeneity measures, not as population genetics.

NOTE ON THE TELOMERIC EDGE. The PHR caller walks inward from the telomeric window and
stops at the first failing window (scripts/phr-detection/
find-multichr-regions-incremental.py). A PHR is therefore only emitted when the
telomeric-most window passes, so every called region contains the telomeric boundary of
its flank: that is a property of the calling rule, not a result, and it is not tabulated.

INPUTS (all pre-existing; nothing is written except the output TSV)
  --phrs        PHR call table, one row per telomere-anchored 500 kb flank.
                data/all-vs-all.1Mb.p95.id95.len.tsv. Columns:
                seq, arm(parm|qarm), self_chr, region_start, region_end,
                chrs_involved, arms_involved. region_start/region_end are 0-based
                offsets into the flank sequence (0 = flank start; the telomere is at
                offset 0 for parm flanks and at the flank end for qarm flanks, see
                scripts/phr-detection/find-multichr-regions-incremental.py). A flank
                with region_start == "." is a flank where no PHR was called.
                Flank names are PanSN "<sample>#<haplotype>#<contig>:<start>-<end>_<chr>_<arm>".
  --communities arm-level Leiden community assignment (15 communities, 41 arms):
                data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv
  --labels      per-sample superpopulation labels. Repo-local default:
                slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/
                nearest_same_superpop_mds_distances.tsv (columns Sample, Superpopulation).
                That file covers exactly the 233 samples (232 HPRC individuals + CHM13)
                of the PHR call set.
  --labels      provenance: the repo-local table is a copy of the Superpopulation column
                of the upstream per-PHR assignment table (see --crosscheck-assign), taken
                through slides/v2-review-zoom/_revision_assets/v9/
                labels_superpop_stats_polish/make_labels_superpop_stats_polish.R.
  --crosscheck  optional second label source, formats accepted:
                  (a) per-PHR assignment table  Name, Sample, ..., Superpopulation
                  (b) sample-level metadata     ChildID, ..., Superpopulation
                Defaults: (a) /moosefs/guarracino/HPRCv2/PHR_III/similarity/
                hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv
                          (b) /moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv
                (b) lacks ~29 of the 232 individuals (see
                paper_prep/manuscript_revision/A1_sample_counts.md), so it can only
                confirm a subset; (a) covers the whole PHR set and its per-sample PHR
                counts are also compared with the call table. A missing or unreadable
                cross-check file is reported and skipped; a label CONFLICT is fatal.

OUTPUT
  scripts/cladistics/phr_by_superpopulation.tsv  (section, population, community, ...)
    section = superpopulation row  (community = "-")
    section = cohort row           (population = ALL_HPRC)
    section = community row        (population x community, fraction within that population)
    section = reference row        (CHM13 and CHM13 x community, reported separately)

CHM13v2.0 is one of the 465 assemblies but is not an HPRC individual; the upstream
heterogeneity tables carry it as EUR. It is reported as its own `reference` section and
is excluded from every superpopulation and cohort aggregate.

Fail-loud policy: the script aborts (exit 1, no partial TSV) if any haplotype lacks a
label, if a sample carries two different labels, if a label is outside the five
superpopulations, if a called PHR's arm has no community assignment, or if the optional
cross-check finds a label conflict.

Usage:
  python3 scripts/cladistics/phr_by_superpopulation.py
  python3 scripts/cladistics/phr_by_superpopulation.py --nperm 20000 --seed 20260201
"""
from __future__ import annotations

import argparse
import collections
import csv
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
REFERENCE_SAMPLE = "CHM13"
FLANK_LEN = 500000  # every flank in the call table is exactly 500 kb (verified below)

DEFAULT_PHRS = ROOT / "data/all-vs-all.1Mb.p95.id95.len.tsv"
DEFAULT_COMMUNITIES = ROOT / "data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv"
DEFAULT_LABELS = (ROOT / "slides/v2-review-zoom/_revision_assets/v9/"
                  "labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv")
DEFAULT_CROSSCHECK = Path("/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv")
DEFAULT_CROSSCHECK_ASSIGN = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
    "hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv")
OUT = ROOT / "scripts/cladistics/phr_by_superpopulation.tsv"

COLUMNS = [
    "section", "population", "community",
    "n_individuals", "n_haplotypes", "n_arms", "n_phrs", "phr_bp",
    "frac_phrs", "frac_bp", "phr_share_over_hap_share",
    "phrs_per_hap_median", "phrs_per_hap_q1", "phrs_per_hap_q3",
    "median_len_bp", "q1_len_bp", "q3_len_bp", "p10_len_bp", "p90_len_bp",
    "n_len_500kb", "pct_len_500kb",
]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def die(msg: str) -> "NoReturn":  # type: ignore[valid-type]
    sys.stderr.write("FATAL: " + msg + "\n")
    raise SystemExit(1)


def quantile(sorted_vals, p):
    """Type-7 (linear interpolation) quantile, same as R's default / numpy's 'linear'."""
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return float(sorted_vals[0])
    h = (n - 1) * p
    lo = int(h)
    hi = min(lo + 1, n - 1)
    return float(sorted_vals[lo]) + (h - lo) * float(sorted_vals[hi] - sorted_vals[lo])


def med_iqr(vals):
    s = sorted(vals)
    if not s:
        return float("nan"), float("nan"), float("nan")
    return quantile(s, 0.5), quantile(s, 0.25), quantile(s, 0.75)


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def load_communities(path: Path):
    """arm (e.g. 'chr10_parm') -> community ('C1'...'C15')"""
    comm = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            arm, c = r["arm"], r["community"]
            if arm in comm and comm[arm] != c:
                die(f"{path}: arm {arm} assigned to both {comm[arm]} and {c}")
            comm[arm] = c
    if not comm:
        die(f"{path}: no community assignments read")
    return comm


def load_labels(path: Path):
    """sample -> superpopulation from the per-sequence MDS label table.

    That table carries one row per called PHR, so it also yields the per-haplotype PHR
    counts used as an independent consistency gate on the call table.
    """
    need = {"Sample", "Superpopulation", "Haplotype"}
    lab = {}
    per_hap = collections.Counter()
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = need - set(rd.fieldnames or [])
        if missing:
            die(f"{path}: missing column(s) {sorted(missing)}")
        for r in rd:
            s, sp = r["Sample"].strip(), r["Superpopulation"].strip()
            if not s:
                continue
            per_hap[(s, r["Haplotype"].strip())] += 1
            if sp == "" or sp.upper() in {"NA", "NAN", "NONE"}:
                die(f"{path}: sample {s} has no superpopulation label ({sp!r})")
            if sp not in SUPERPOPS:
                die(f"{path}: sample {s} has label {sp!r}, not one of {SUPERPOPS}")
            if s in lab and lab[s] != sp:
                die(f"{path}: sample {s} labelled both {lab[s]} and {sp}")
            lab[s] = sp
    if not lab:
        die(f"{path}: no labels read")
    return lab, per_hap


def load_crosscheck(path: Path, labels, per_sample_call_counts, n_called):
    """Cross-check the primary labels against a second, named file.

    Accepts either schema seen in this project:
      * per-PHR assignment table with columns Name, Sample, ..., Superpopulation
        (upstream hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv). Because it holds
        one row per called PHR, its per-sample row counts are also compared with the
        call table's per-sample PHR counts.
      * sample-level production metadata with columns ChildID, Subpopulation,
        Superpopulation (upstream hprc-sequence-production.tsv). That table is known to
        omit ~29 of the 232 individuals, so it can only confirm a subset.

    Never modifies `labels`. Returns None (with a printed notice) when unavailable.
    """
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
    except OSError as e:
        print(f"NOTE: cross-check SKIPPED, cannot open {path}: {e}")
        return None
    per_phr = "Name" in header and "Sample" in header and "Superpopulation" in header
    per_sample = "ChildID" in header and "Superpopulation" in header
    if not (per_phr or per_sample):
        print(f"NOTE: cross-check SKIPPED, {path} has neither the per-PHR "
              f"(Name/Sample/Superpopulation) nor the sample-level "
              f"(ChildID/Superpopulation) schema")
        return None
    other = {}
    counts = collections.Counter()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if per_phr:
                smp = (r.get("Sample") or "").strip()
                sp = (r.get("Superpopulation") or "").strip()
                if smp:
                    counts[smp] += 1
            else:
                smp = (r.get("ChildID") or "").strip()
                sp = (r.get("Superpopulation") or "").strip()
            if smp and sp:
                if smp in other and other[smp] != sp:
                    die(f"{path}: sample {smp} carries both {other[smp]} and {sp}")
                other[smp] = sp

    n_cmp = n_conf = 0
    for smp, sp in sorted(other.items()):
        if smp == REFERENCE_SAMPLE or smp not in labels:
            continue
        n_cmp += 1
        if labels[smp] != sp:
            n_conf += 1
            print(f"CONFLICT {smp}: primary={labels[smp]} cross-check={sp}")
    if n_conf:
        die(f"{path}: {n_conf} superpopulation label conflict(s) with the primary source "
            f"{DEFAULT_LABELS.name}")
    schema = "per-PHR assignment" if per_phr else "sample-level metadata"
    missing = sum(1 for s in per_sample_call_counts if s not in other
                  and s != REFERENCE_SAMPLE)
    print(f"NOTE: label cross-check vs {path.name} ({schema}): {n_cmp} samples "
          f"labelled by both sources, 0 conflicts, {missing} of "
          f"{len(per_sample_call_counts) - (1 if REFERENCE_SAMPLE in per_sample_call_counts else 0)} "
          f"HPRC individuals absent from that file")
    if per_phr:
        if sum(counts.values()) != n_called:
            die(f"{path}: {sum(counts.values())} PHR rows but the call table has "
                f"{n_called} called PHRs")
        bad = {k: (per_sample_call_counts.get(k, 0), counts.get(k, 0))
               for k in set(counts) | set(per_sample_call_counts)
               if counts.get(k, 0) != per_sample_call_counts.get(k, 0)}
        if bad:
            die(f"{path}: per-sample PHR counts disagree for {len(bad)} sample(s), "
                f"e.g. {sorted(bad.items())[:5]}")
        print(f"NOTE: per-sample PHR counts identical between that file and the call "
              f"table across all {len(counts)} samples")
    return n_cmp, 0, missing


def load_phrs(path: Path, comm, labels):
    """Return list of dicts for called PHRs; abort on unlabelled sample or unassigned arm."""
    rows = []
    n_flanks = 0
    flank_lens = set()
    unlabelled = collections.Counter()
    no_comm = collections.Counter()
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"seq", "arm", "self_chr", "region_start", "region_end"}
        missing = need - set(rd.fieldnames or [])
        if missing:
            die(f"{path}: missing column(s) {sorted(missing)}")
        for r in rd:
            seq = r["seq"]
            n_flanks += 1
            base = seq
            if ":" in base:
                # contig:start-end_<chr>_<arm>; flank length from the inclusive coords
                try:
                    coords = base.split(":")[1].split("_")[0]
                    s, e = coords.split("-")
                    flank_lens.add(int(e) - int(s) + 1)
                except Exception:
                    pass
            parts = seq.split("#")
            if len(parts) < 3:
                die(f"{path}: cannot parse PanSN name {seq!r}")
            sample, hap = parts[0], parts[1]
            arm_bare = "parm" if r["arm"] == "parm" else "qarm"
            chrom = r["self_chr"].split(":")[0]
            arm = f"{chrom}_{arm_bare}"
            if r["region_start"] == "." or r["region_end"] == ".":
                continue  # flank with no called PHR
            rs, re = int(r["region_start"]), int(r["region_end"])
            if re <= rs:
                die(f"{path}: non-positive PHR extent for {seq}: {rs}-{re}")
            if sample not in labels:
                unlabelled[sample] += 1
                continue
            if arm not in comm:
                no_comm[arm] += 1
                continue
            rows.append({
                "seq": seq, "sample": sample, "hap": hap,
                "arm": arm, "community": comm[arm],
                "length": re - rs,
            })
    if unlabelled:
        die(f"{path}: {sum(unlabelled.values())} called PHRs from "
            f"{len(unlabelled)} unlabelled sample(s): "
            + ", ".join(sorted(unlabelled)[:20]))
    if no_comm:
        die(f"{path}: {sum(no_comm.values())} called PHRs on arms with no community "
            f"assignment: {sorted(no_comm)}")
    print(f"NOTE: flank table {path.name}: {n_flanks} flanks, "
          f"{len(rows)} called PHRs, flank lengths observed = {sorted(flank_lens)}")
    if flank_lens - {FLANK_LEN}:
        print("WARNING: flank lengths other than 500 kb are present; "
              "telomeric/distal edge flags still use the 500 kb offset convention")
    return rows


# ---------------------------------------------------------------------------
# descriptive statistics
# ---------------------------------------------------------------------------
def summarise(rows, n_ind, n_hap):
    """Descriptive summary of a set of called PHRs."""
    n = len(rows)
    bp = sum(r["length"] for r in rows)
    per_hap = collections.Counter()
    for r in rows:
        per_hap[(r["sample"], r["hap"])] += 1
    hp = list(per_hap.values())
    hp_med, hp_q1, hp_q3 = med_iqr(hp) if hp else (float("nan"),) * 3
    lens = sorted(r["length"] for r in rows)
    n_500 = sum(1 for v in lens if v == FLANK_LEN)
    return {
        "n_individuals": n_ind,
        "n_haplotypes": n_hap,
        "n_arms": len({r["arm"] for r in rows}),
        "n_phrs": n,
        "phr_bp": bp,
        "phrs_per_hap_median": hp_med, "phrs_per_hap_q1": hp_q1, "phrs_per_hap_q3": hp_q3,
        "median_len_bp": quantile(lens, 0.5) if n else float("nan"),
        "q1_len_bp": quantile(lens, 0.25) if n else float("nan"),
        "q3_len_bp": quantile(lens, 0.75) if n else float("nan"),
        "p10_len_bp": quantile(lens, 0.10) if n else float("nan"),
        "p90_len_bp": quantile(lens, 0.90) if n else float("nan"),
        "n_len_500kb": n_500,
        "pct_len_500kb": 100.0 * n_500 / n if n else float("nan"),
    }


def fmt(v):
    if isinstance(v, float):
        return "nan" if v != v else f"{v:.4f}".rstrip("0").rstrip(".")
    return str(v)


# ---------------------------------------------------------------------------
# permutation checks (haplotype-level label permutation)
# ---------------------------------------------------------------------------
def kruskal_h(groups):
    """Kruskal-Wallis H with tie correction. groups = list of lists of floats."""
    allv = [(v, i) for i, g in enumerate(groups) for v in g]
    n = len(allv)
    if n == 0:
        return float("nan")
    allv.sort()
    ranks = [0.0] * n
    i = 0
    tie_term = 0.0
    while i < n:
        j = i
        while j + 1 < n and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = r
        t = j - i + 1
        if t > 1:
            tie_term += t ** 3 - t
        i = j + 1
    rank_sum = [0.0] * len(groups)
    for (v, g), r in zip(allv, ranks):
        rank_sum[g] += r
    h = 0.0
    for g, rs in enumerate(rank_sum):
        ng = len(groups[g])
        if ng:
            h += rs * rs / ng
    h = 12.0 / (n * (n + 1)) * h - 3 * (n + 1)
    corr = 1.0 - tie_term / (n ** 3 - n) if n > 1 else 1.0
    return h / corr if corr > 0 else h


def permutation_checks(rows, nperm, seed):
    """Superpopulation labels are permuted across haplotypes while each haplotype keeps
    its own PHR set, so within-haplotype clustering of community membership is retained
    (PHRs are not independent observations; a PHR-level permutation would not be)."""
    pops = [p for p in SUPERPOPS if any(r["sample"] != REFERENCE_SAMPLE for r in rows)]
    hap_pop = {}
    hap_bp = collections.defaultdict(lambda: collections.Counter())
    hap_len = collections.defaultdict(list)
    for r in rows:
        if r["sample"] == REFERENCE_SAMPLE:
            continue
        key = (r["sample"], r["hap"])
        hap_pop[key] = r["population"]
        hap_bp[key][r["community"]] += r["length"]
        hap_len[key].append(r["length"])
    keys = sorted(hap_pop)
    communities = sorted({r["community"] for r in rows}, key=lambda c: int(c[1:]))
    if len(pops) < 2 or len(keys) < 10:
        print("NOTE: permutation checks skipped (too few superpopulations or haplotypes)")
        return None

    # per-haplotype median PHR length, precomputed
    hap_medlen = {k: quantile(sorted(v), 0.5) for k, v in hap_len.items()}
    vectors = {k: [hap_bp[k][c] for c in communities] for k in keys}
    # population totals under the observed assignment
    n_pop = collections.Counter(hap_pop.values())

    def stats(assign):
        """chi-square (community x superpopulation, bp) and Kruskal-Wallis H."""
        table = {p: [0] * len(communities) for p in pops}
        for k in keys:
            row = table[assign[k]]
            vec = vectors[k]
            for i, v in enumerate(vec):
                row[i] += v
        row_tot = {p: sum(table[p]) for p in pops}
        tot = sum(row_tot.values())
        chi2 = 0.0
        for p in pops:
            for i in range(len(communities)):
                col_tot = sum(table[q][i] for q in pops)
                e = row_tot[p] * col_tot / tot
                if e > 0:
                    chi2 += (table[p][i] - e) ** 2 / e
        groups = collections.defaultdict(list)
        for k in keys:
            groups[assign[k]].append(hap_medlen[k])
        return chi2, kruskal_h([groups[p] for p in pops])

    obs_chi2, obs_h = stats(hap_pop)
    rng = random.Random(seed)
    base = [hap_pop[k] for k in keys]
    ge_c = ge_l = 0
    for _ in range(nperm):
        perm = base[:]
        rng.shuffle(perm)
        c2, h = stats(dict(zip(keys, perm)))
        if c2 >= obs_chi2:
            ge_c += 1
        if h >= obs_h:
            ge_l += 1
    n_total_bp = sum(sum(v) for v in vectors.values())
    denom = n_total_bp * (min(len(pops), len(communities)) - 1)
    return {
        "n_hap": len(keys), "nperm": nperm, "seed": seed,
        "chi2": obs_chi2, "comm_p": (ge_c + 1) / (nperm + 1),
        "cramers_v": (obs_chi2 / denom) ** 0.5 if denom > 0 else float("nan"),
        "comm_df": (len(pops) - 1) * (len(communities) - 1),
        "h": obs_h, "len_p": (ge_l + 1) / (nperm + 1),
        "len_df": len(pops) - 1, "pops": pops, "n_comm": len(communities),
        "n_pop": dict(n_pop),
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--phrs", type=Path, default=DEFAULT_PHRS)
    ap.add_argument("--communities", type=Path, default=DEFAULT_COMMUNITIES)
    ap.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    ap.add_argument("--crosscheck", type=Path, default=DEFAULT_CROSSCHECK,
                    help="second label source, sample-level or per-PHR schema")
    ap.add_argument("--crosscheck-assign", type=Path, default=DEFAULT_CROSSCHECK_ASSIGN,
                    help="upstream per-PHR assignment table (Name/Sample/Superpopulation)")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260201)
    a = ap.parse_args()

    comm = load_communities(a.communities)
    labels, label_per_hap = load_labels(a.labels)
    rows = load_phrs(a.phrs, comm, labels)
    # attach sample labels once; refusals already raised in load_phrs
    for r in rows:
        r["population"] = labels[r["sample"]]

    n_total = len(rows)
    per_sample_call = collections.Counter(r["sample"] for r in rows)
    load_crosscheck(a.crosscheck, labels, per_sample_call, n_total)
    load_crosscheck(a.crosscheck_assign, labels, per_sample_call, n_total)

    # ---- consistency gates ------------------------------------------------
    n_total = len(rows)
    n_hprc = sum(1 for r in rows if r["sample"] != REFERENCE_SAMPLE)
    n_ref = n_total - n_hprc
    cohort_samples = {r["sample"] for r in rows if r["sample"] != REFERENCE_SAMPLE}
    cohort_haps = {(r["sample"], r["hap"]) for r in rows if r["sample"] != REFERENCE_SAMPLE}
    print(f"GATE total called PHRs = {n_total} (HPRC {n_hprc} + CHM13 {n_ref}); "
          f"HPRC individuals with >=1 PHR = {len(cohort_samples)}, "
          f"haplotypes = {len(cohort_haps)}")
    all_haps = {(r["sample"], r["hap"]) for r in rows}
    print(f"GATE haplotypes carrying >=1 called PHR (incl. CHM13) = {len(all_haps)}")

    # independent consistency gate: the label table has one row per called PHR, so its
    # per-haplotype row counts must equal the per-haplotype PHR counts of the call table
    call_per_hap = collections.Counter((r["sample"], r["hap"]) for r in rows)
    if sum(label_per_hap.values()) != n_total:
        die(f"label table {a.labels.name} has {sum(label_per_hap.values())} rows but the "
            f"call table has {n_total} called PHRs; they are not the same PHR set")
    diff = {k: (call_per_hap.get(k, 0), label_per_hap.get(k, 0))
            for k in set(call_per_hap) | set(label_per_hap)
            if call_per_hap.get(k, 0) != label_per_hap.get(k, 0)}
    if diff:
        die(f"per-haplotype PHR counts disagree between call table and label table for "
            f"{len(diff)} haplotype(s), e.g. {sorted(diff.items())[:5]}")
    print(f"GATE per-haplotype PHR counts identical between {a.phrs.name} and "
          f"{a.labels.name} across all {len(call_per_hap)} haplotypes")

    # ---- rows -------------------------------------------------------------
    out_rows = []

    def emit(section, pop, community, agg, frac_phrs, frac_bp, n_ind, n_hap,
             hap_share=None):
        rec = {"section": section, "population": pop, "community": community}
        rec.update({k: agg[k] for k in COLUMNS if k in agg})
        rec["n_individuals"] = n_ind
        rec["n_haplotypes"] = n_hap
        rec["frac_phrs"] = frac_phrs
        rec["frac_bp"] = frac_bp
        if hap_share:
            rec["phr_share_over_hap_share"] = frac_phrs / hap_share
        out_rows.append(rec)

    # cohort totals (232 HPRC individuals, CHM13 excluded)
    cohort_rows = [r for r in rows if r["sample"] != REFERENCE_SAMPLE]
    ind_by_pop = collections.defaultdict(set)
    hap_by_pop = collections.defaultdict(set)
    for r in cohort_rows:
        ind_by_pop[r["population"]].add(r["sample"])
        hap_by_pop[r["population"]].add((r["sample"], r["hap"]))
    cohort_bp = sum(r["length"] for r in cohort_rows)

    cohort_agg = summarise(cohort_rows, len(cohort_samples), len(cohort_haps))
    emit("cohort", "ALL_HPRC", "-", cohort_agg, 1.0, 1.0,
         len(cohort_samples), len(cohort_haps), hap_share=1.0)

    for p in SUPERPOPS:
        sub = [r for r in cohort_rows if r["population"] == p]
        if not sub:
            print(f"WARNING: no PHRs for superpopulation {p}")
            continue
        agg = summarise(sub, len(ind_by_pop[p]), len(hap_by_pop[p]))
        emit("superpopulation", p, "-", agg,
             len(sub) / len(cohort_rows), sum(r["length"] for r in sub) / cohort_bp,
             len(ind_by_pop[p]), len(hap_by_pop[p]),
             hap_share=len(hap_by_pop[p]) / len(cohort_haps))

    # community membership distribution, per population and cohort
    communities = sorted({r["community"] for r in rows}, key=lambda c: int(c[1:]))
    for p in ["ALL_HPRC"] + SUPERPOPS:
        base = cohort_rows if p == "ALL_HPRC" else [r for r in cohort_rows
                                                    if r["population"] == p]
        base_bp = sum(r["length"] for r in base)
        for c in communities:
            sub = [r for r in base if r["community"] == c]
            if not sub:
                continue
            agg = summarise(sub, len({r["sample"] for r in sub}),
                            len({(r["sample"], r["hap"]) for r in sub}))
            base_hap = (len(cohort_haps) if p == "ALL_HPRC"
                        else len({(r["sample"], r["hap"]) for r in base}))
            emit("community", p, c, agg,
                 len(sub) / len(base), sum(r["length"] for r in sub) / base_bp,
                 agg["n_individuals"], agg["n_haplotypes"],
                 hap_share=agg["n_haplotypes"] / base_hap)

    # CHM13 reference lines, reported separately and excluded from all of the above
    ref_rows = [r for r in rows if r["sample"] == REFERENCE_SAMPLE]
    if ref_rows:
        ref_agg = summarise(ref_rows, 1, 1)
        emit("reference", REFERENCE_SAMPLE, "-", ref_agg, float("nan"), float("nan"), 1, 1)
        ref_bp = sum(r["length"] for r in ref_rows)
        for c in communities:
            sub = [r for r in ref_rows if r["community"] == c]
            if not sub:
                continue
            agg = summarise(sub, 1, 1)
            emit("reference", REFERENCE_SAMPLE, c, agg,
                 len(sub) / len(ref_rows), sum(r["length"] for r in sub) / ref_bp, 1, 1)

    # ---- write ------------------------------------------------------------
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for rec in out_rows:
            w.writerow([fmt(rec.get(c, "")) for c in COLUMNS])
    print(f"WROTE {a.out} ({len(out_rows)} rows)")

    # ---- human-readable report -------------------------------------------
    print("\n== per superpopulation ==")
    hdr = ["pop", "ind", "hap", "arms", "phrs", "bp", "len_p10", "len_q1", "len_med",
           "len_q3", "len_p90", "phrs/hap_med", "pct500kb", "phr_share/hap_share"]
    print("\t".join(hdr))
    for rec in out_rows:
        if rec["section"] in ("cohort", "superpopulation"):
            print("\t".join(fmt(rec[c]) for c in [
                "population", "n_individuals", "n_haplotypes", "n_arms", "n_phrs",
                "phr_bp", "p10_len_bp", "q1_len_bp", "median_len_bp", "q3_len_bp",
                "p90_len_bp", "phrs_per_hap_median", "pct_len_500kb",
                "phr_share_over_hap_share"]))

    print("\n== fraction of PHR bp per community, by superpopulation ==")
    print("community\t" + "\t".join(["ALL_HPRC"] + SUPERPOPS))
    for c in communities:
        cells = []
        for p in ["ALL_HPRC"] + SUPERPOPS:
            hit = [r for r in out_rows
                   if r["section"] == "community" and r["population"] == p
                   and r["community"] == c]
            cells.append(f"{float(hit[0]['frac_bp']):.4f}" if hit else "-")
        print(c + "\t" + "\t".join(cells))

    print("\n== permutation checks (labels permuted across haplotypes) ==")
    chk = permutation_checks(rows, a.nperm, a.seed)
    if chk:
        print(f"n_haplotypes_tested={chk['n_hap']} nperm={chk['nperm']} seed={chk['seed']}")
        print(f"community x superpopulation (PHR bp, {chk['n_pop']}): "
              f"Cramers V={chk['cramers_v']:.4f} chi2={chk['chi2']:.1f} "
              f"df={chk['comm_df']} p_perm={chk['comm_p']:.4f}")
        print(f"per-haplotype median PHR length: H={chk['h']:.3f} "
              f"df={chk['len_df']} p_perm={chk['len_p']:.4f}")
        print("Both are DESCRIPTIVE heterogeneity measures on haplotype-level "
              "permutations; haplotypes of the same individual are not independent "
              "and no population-genetic inference should be drawn from them.")


if __name__ == "__main__":
    main()
