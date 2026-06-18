#!/usr/bin/env python3
"""
Summarise WashU pedigree interchromosomal patch tract lengths.

Run from the repository root:

    python3 scripts/pedigree/patch_tract_lengths.py

Default inputs are the moosefs pedigree patch tables used by the manuscript:

  * /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv
  * /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv

The primary analysis uses the direct WashU table and the direct
assembly-derived `patch_size` field.  The manuscript high-quality denominator
is reproduced by filtering to interchromosomal patches with min_score >= 0.8
and 500 bp <= patch_size <= 100 kb.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WASHU = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv"
)
DEFAULT_ALL = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv"
)
DEFAULT_TSV = REPO_ROOT / "scripts" / "pedigree" / "patch_tract_length_summary.tsv"
DEFAULT_LOWER_TSV = REPO_ROOT / "scripts" / "pedigree" / "patch_tract_lower_merge_summary.tsv"
DEFAULT_MD = REPO_ROOT / "paper_prep" / "_brainstorming" / "pedigree_patch_tract_lengths.md"
DEFAULT_LOWER_DIR = REPO_ROOT / "paper_prep" / "_brainstorming" / "pedigree_patch_tract_lower_untangle"
DEFAULT_GRAPH = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/"
    "washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og"
)
DEFAULT_TARGET_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle")

CONVERSION_RANGE = (22, 95)
CROSSOVER_RANGE = (318, 688)
NEAR_FACTOR = 2.0
GRAPH_WFMASH_SEGMENT_LENGTH_BP = 1000
ODGI_MERGE_DIST_BP = 1000

LOWER_PAIRINGS = {
    "PAN027_vs_PAN010": ("PAN027#1", "PAN027 maternal (hap1) vs PAN010 (mother)"),
    "PAN027_vs_PAN011": ("PAN027#2", "PAN027 paternal (hap2) vs PAN011 (father)"),
    "PAN028_vs_PAN027": ("PAN028#1", "PAN028 maternal (hap1) vs PAN027 (mother)"),
}


@dataclass(frozen=True)
class Patch:
    label: str
    pattern: str
    community_status: str
    overlaps_phr: str
    has_phr: str
    patch_size: int


@dataclass(frozen=True)
class LowerPatch:
    label: str
    pattern: str
    patch_size: int


def proportion(n: int, d: int) -> str:
    return f"{n / d:.6f}" if d else "NA"


def quantile(sorted_values: list[int], prob: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    h = (len(sorted_values) - 1) * prob
    lo = int(h)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = h - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


def fmt_num(value: float | None) -> str:
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.1f}"


def parse_path_info(path_name: str) -> tuple[str, str, str, str, str] | None:
    match = re.search(r"^(PAN\d+)#([12])#.*_(chr(?:[0-9]+|X|Y))_([pq]arm)$", path_name)
    if not match:
        return None
    sample, hap, chrom, arm = match.groups()
    chrarm = chrom + ("p" if arm == "parm" else "q")
    return sample, hap, chrom, arm, chrarm


def in_range(value: int, bounds: tuple[int, int]) -> bool:
    return bounds[0] <= value <= bounds[1]


def in_near_range(value: int, bounds: tuple[int, int]) -> bool:
    lo = bounds[0] / NEAR_FACTOR
    hi = bounds[1] * NEAR_FACTOR
    return lo <= value <= hi


def load_washu_hq(path: Path) -> list[Patch]:
    if not path.exists():
        raise FileNotFoundError(f"WashU patch table not found: {path}")
    patches: list[Patch] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "label",
            "patch_size",
            "min_score",
            "is_interchr",
            "pattern",
            "community_status",
            "overlaps_phr",
            "has_phr",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        for row in reader:
            size = int(row["patch_size"])
            if row["is_interchr"] != "True":
                continue
            if float(row["min_score"]) < 0.8:
                continue
            if not 500 <= size <= 100_000:
                continue
            patches.append(
                Patch(
                    label=row["label"],
                    pattern=row["pattern"],
                    community_status=row["community_status"],
                    overlaps_phr=row["overlaps_phr"],
                    has_phr=row["has_phr"],
                    patch_size=size,
                )
            )
    return patches


def load_aggregate_washu_count(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ds" not in (reader.fieldnames or []):
            return None
        return sum(1 for row in reader if row["ds"] == "WashU")


def make_summary_row(group_type: str, group_value: str, patches: Iterable[Patch]) -> dict[str, str]:
    rows = list(patches)
    lengths = sorted(p.patch_size for p in rows)
    n = len(lengths)
    conv_n = sum(in_range(x, CONVERSION_RANGE) for x in lengths)
    cross_n = sum(in_range(x, CROSSOVER_RANGE) for x in lengths)
    near_conv_n = sum(in_near_range(x, CONVERSION_RANGE) for x in lengths)
    near_cross_n = sum(in_near_range(x, CROSSOVER_RANGE) for x in lengths)
    min_len = lengths[0] if lengths else None
    return {
        "group_type": group_type,
        "group_value": group_value,
        "n": str(n),
        "length_field": "patch_size",
        "min_bp": fmt_num(float(min_len) if min_len is not None else None),
        "q25_bp": fmt_num(quantile(lengths, 0.25)),
        "median_bp": fmt_num(quantile(lengths, 0.50)),
        "q75_bp": fmt_num(quantile(lengths, 0.75)),
        "q90_bp": fmt_num(quantile(lengths, 0.90)),
        "q95_bp": fmt_num(quantile(lengths, 0.95)),
        "max_bp": fmt_num(float(lengths[-1]) if lengths else None),
        "conversion_22_95_n": str(conv_n),
        "conversion_22_95_prop": proportion(conv_n, n),
        "near_conversion_11_190_n": str(near_conv_n),
        "near_conversion_11_190_prop": proportion(near_conv_n, n),
        "crossover_318_688_n": str(cross_n),
        "crossover_318_688_prop": proportion(cross_n, n),
        "near_crossover_159_1376_n": str(near_cross_n),
        "near_crossover_159_1376_prop": proportion(near_cross_n, n),
        "conversion_range_observable": str(bool(min_len is not None and min_len <= CONVERSION_RANGE[1])),
        "crossover_range_observable": str(bool(min_len is not None and min_len <= CROSSOVER_RANGE[1])),
        "resolution_note": (
            f"current committed WashU patch table is from e50000.m{ODGI_MERGE_DIST_BP}; "
            f"wfmash segment-length={GRAPH_WFMASH_SEGMENT_LENGTH_BP} is not a tract-length bound; "
            "counts below observed min are not biological absence calls"
        ),
    }


def grouped_rows(patches: list[Patch]) -> list[dict[str, str]]:
    rows = [make_summary_row("overall", "all_high_quality_interchromosomal", patches)]
    for attr, group_type in [
        ("pattern", "pattern"),
        ("community_status", "community_status"),
        ("overlaps_phr", "overlaps_phr"),
        ("has_phr", "has_phr"),
        ("label", "transmission"),
    ]:
        grouped: dict[str, list[Patch]] = defaultdict(list)
        for patch in patches:
            grouped[getattr(patch, attr)].append(patch)
        for key in sorted(grouped, key=lambda k: (-len(grouped[k]), k)):
            rows.append(make_summary_row(group_type, key, grouped[key]))
    return rows


def run_lower_untangle(out_dir: Path, graph: Path, target_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for pair in LOWER_PAIRINGS:
        out_bed = out_dir / f"{pair}.e50000.m0.n1.bed.gz"
        err_path = out_dir / f"{pair}.e50000.m0.n1.err"
        targets = target_dir / f"targets_{pair}.txt"
        if out_bed.exists() and out_bed.stat().st_size > 0:
            continue
        cmd = [
            "odgi",
            "untangle",
            "-t",
            "48",
            "-i",
            str(graph),
            "-R",
            str(targets),
            "-e",
            "50000",
            "-m",
            "0",
            "-j",
            "0",
            "-n",
            "1",
        ]
        with out_bed.open("wb") as out_handle, err_path.open("wb") as err_handle:
            proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_handle)
            proc2 = subprocess.Popen(["gzip", "-c"], stdin=proc1.stdout, stdout=out_handle)
            assert proc1.stdout is not None
            proc1.stdout.close()
            rc2 = proc2.wait()
            rc1 = proc1.wait()
        if rc1 != 0 or rc2 != 0:
            raise RuntimeError(f"lower untangle failed for {pair}: odgi={rc1}, gzip={rc2}")


def load_lower_merge_patches(lower_dir: Path) -> list[LowerPatch]:
    patches: list[LowerPatch] = []
    for pair, (query_prefix, label) in LOWER_PAIRINGS.items():
        bed_path = lower_dir / f"{pair}.e50000.m0.n1.bed.gz"
        if not bed_path.exists():
            continue
        by_query: dict[str, list[tuple[int, int, str, str, str, float]]] = defaultdict(list)
        with gzip.open(bed_path, "rt") as handle:
            next(handle)
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or not fields[0].startswith(query_prefix):
                    continue
                query_info = parse_path_info(fields[0])
                ref_info = parse_path_info(fields[3])
                if query_info is None or ref_info is None:
                    continue
                _, ref_hap, ref_chr, _, ref_chrarm = ref_info
                by_query[fields[0]].append(
                    (int(fields[1]), int(fields[2]), ref_chr, ref_chrarm, ref_hap, float(fields[6]))
                )

        for query_name, segments in by_query.items():
            query_info = parse_path_info(query_name)
            if query_info is None:
                continue
            _, _, query_chr, _, query_chrarm = query_info
            segments.sort()
            merged: list[dict[str, object]] = []
            for start, end, ref_chr, ref_chrarm, ref_hap, score in segments:
                key = (ref_chrarm, ref_hap)
                if merged and merged[-1]["key"] == key and start <= int(merged[-1]["end"]):
                    merged[-1]["end"] = max(int(merged[-1]["end"]), end)
                    merged[-1]["min_score"] = min(float(merged[-1]["min_score"]), score)
                elif merged and merged[-1]["key"] == key and start == int(merged[-1]["end"]):
                    merged[-1]["end"] = end
                    merged[-1]["min_score"] = min(float(merged[-1]["min_score"]), score)
                else:
                    merged.append(
                        {
                            "start": start,
                            "end": end,
                            "ref_chr": ref_chr,
                            "ref_chrarm": ref_chrarm,
                            "ref_hap": ref_hap,
                            "key": key,
                            "min_score": score,
                        }
                    )

            for idx, item in enumerate(merged):
                length = int(item["end"]) - int(item["start"])
                if item["ref_chr"] == query_chr:
                    continue
                if float(item["min_score"]) < 0.8 or length <= 0 or length > 100_000:
                    continue
                pattern = "interchr_merged"
                left = merged[idx - 1] if idx > 0 else None
                right = merged[idx + 1] if idx + 1 < len(merged) else None
                if left and right and left["ref_chrarm"] == query_chrarm and right["ref_chrarm"] == query_chrarm:
                    if left["ref_hap"] != right["ref_hap"]:
                        pattern = "crossover_like"
                    elif left["ref_hap"] != item["ref_hap"]:
                        pattern = "gene_conversion_like"
                    else:
                        pattern = "sandwich_same_hap"
                patches.append(LowerPatch(label=label, pattern=pattern, patch_size=length))
    return patches


def lower_grouped_rows(patches: list[LowerPatch]) -> list[dict[str, str]]:
    rows = [make_summary_row("lower_merge_overall", "m0_n1_merged_interchromosomal", patches)]  # type: ignore[arg-type]
    for attr, group_type in [("pattern", "lower_merge_pattern"), ("label", "lower_merge_transmission")]:
        grouped: dict[str, list[LowerPatch]] = defaultdict(list)
        for patch in patches:
            grouped[getattr(patch, attr)].append(patch)
        for key in sorted(grouped, key=lambda k: (-len(grouped[k]), k)):
            rows.append(make_summary_row(group_type, key, grouped[key]))  # type: ignore[arg-type]
    for row in rows:
        row["resolution_note"] = (
            "lower-merge odgi untangle -e 50000 -m 0 -j 0 -n 1; "
            "adjacent best-hit intervals merged by donor chrarm/haplotype; "
            "not the m1000 high-quality patch denominator"
        )
    return rows


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_markdown(
    rows: list[dict[str, str]],
    lower_rows: list[dict[str, str]],
    path: Path,
    washu_path: Path,
    all_path: Path,
    aggregate_count: int | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    overall = rows[0]
    text = f"""# WashU Pedigree Patch Tract-Length Summary

Primary input table:
`{washu_path}`

Manuscript aggregate table cross-check:
`{all_path}`

Length field used: `patch_size` (assembly-derived child-flank interval length
from the WashU `odgi untangle` recombination patch table).

High-quality filter reproduced from `submission/paper.tex`: interchromosomal
patches with `min_score >= 0.8` and `500 <= patch_size <= 100000`.

Resolution audit:

- The committed WashU untangle BEDs are
  `PAN027_vs_PAN010.e50000.m1000.bed.gz`,
  `PAN027_vs_PAN011.e50000.m1000.bed.gz` and
  `PAN028_vs_PAN027.e50000.m1000.bed.gz`.
- The original `odgi untangle` command used `-e 50000 -m 1000`.
- The WashU PGGB graph was induced with wfmash `segment-length: 1000`; the
  PGGB log shows `wfmash -s 1000 -l 3000`. This is a graph-construction
  seed scale, not a hard lower bound on alignment or tract lengths.
- The observed minimum high-quality `patch_size` is {overall['min_bp']} bp.
  Therefore the current high-confidence patch table, which was produced from
  `odgi untangle -m 1000`, is effectively left-truncated at about 1 kb for
  tract-length comparison.
- A lower-merge rerun from the same graph is feasible and was run with
  `odgi untangle -e 50000 -m 0 -j 0 -n 1` for the three assayed transmissions.
  The companion parser merges adjacent best-hit intervals assigned to the same
  donor chrarm/haplotype and summarizes interchromosomal merged intervals with
  `min_score >= 0.8` and no 500 bp lower-size cutoff. This is a graph/untangle-
  resolved interval analysis, not the exact high-quality m1000 patch table used
  for the 494/538 within-community statistic.

Primary denominator: N = {overall['n']} high-quality WashU interchromosomal
candidate patches.
"""
    if aggregate_count is not None:
        text += (
            f"\nAggregate table cross-check: `ds == \"WashU\"` contains "
            f"{aggregate_count} rows.\n"
        )
    text += f"""
Overall length distribution: min {overall['min_bp']} bp; Q1 {overall['q25_bp']}
bp; median {overall['median_bp']} bp; Q3 {overall['q75_bp']} bp; Q90
{overall['q90_bp']} bp; Q95 {overall['q95_bp']} bp; max {overall['max_bp']} bp.

Comparison ranges:

- Short primate NCO/gene-conversion tract means cited in the manuscript:
  {CONVERSION_RANGE[0]}--{CONVERSION_RANGE[1]} bp.
- Longer primate CO-associated tract means cited in the manuscript:
  {CROSSOVER_RANGE[0]}--{CROSSOVER_RANGE[1]} bp.
- "Near" windows are descriptive two-fold windows around those cited ranges:
  11--190 bp and 159--1376 bp. They are not event validation criteria.

Current-table compatibility counts:

- In {CONVERSION_RANGE[0]}--{CONVERSION_RANGE[1]} bp conversion-like range:
  {overall['conversion_22_95_n']}/{overall['n']}
  ({float(overall['conversion_22_95_prop']) * 100:.1f}%; below current-table
  resolution).
- Near conversion-like two-fold window:
  {overall['near_conversion_11_190_n']}/{overall['n']}
  ({float(overall['near_conversion_11_190_prop']) * 100:.1f}%; below
  current-table resolution).
- In {CROSSOVER_RANGE[0]}--{CROSSOVER_RANGE[1]} bp CO-associated range:
  {overall['crossover_318_688_n']}/{overall['n']}
  ({float(overall['crossover_318_688_prop']) * 100:.1f}%; below current-table
  resolution).
- Near CO-associated two-fold window:
  {overall['near_crossover_159_1376_n']}/{overall['n']}
  ({float(overall['near_crossover_159_1376_prop']) * 100:.1f}%; overlaps the
  observed lower edge but is still descriptive only).

Interpretation: the high-quality WashU candidate patch intervals in the current
manuscript tables are not resolved in the 22--95 bp or 318--688 bp cited tract
ranges. The apparent zero counts in those bins are therefore a resolution limit
of the existing graph/untangle/patch table, not evidence that such tract lengths
are biologically absent. This analysis supports only a current-table
compatibility/proportion statement about assembly-derived patch sizes; it does
not validate event-level conversion or crossover mechanisms.

## Summary Table

Full tabular output: `scripts/pedigree/patch_tract_length_summary.tsv`.
Lower-merge tabular output, when generated:
`scripts/pedigree/patch_tract_lower_merge_summary.tsv`.

| group_type | group_value | n | min | q25 | median | q75 | q95 | max | conv 22-95 | CO 318-688 | near CO 159-1376 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
"""
    for row in rows:
        text += (
            f"| {row['group_type']} | {row['group_value']} | {row['n']} | "
            f"{row['min_bp']} | {row['q25_bp']} | {row['median_bp']} | "
            f"{row['q75_bp']} | {row['q95_bp']} | {row['max_bp']} | "
            f"{row['conversion_22_95_n']} | {row['crossover_318_688_n']} | "
            f"{row['near_crossover_159_1376_n']} |\n"
        )
    if lower_rows:
        lower_overall = lower_rows[0]
        text += f"""
## Lower-Merge Untangle Summary

Input lower-merge BEDs are generated by:
`bash scripts/pedigree/run_patch_tract_lower_merge.sh`

The generated intermediates are:
`paper_prep/_brainstorming/pedigree_patch_tract_lower_untangle/*.e50000.m0.n1.bed.gz`
and are intentionally not committed because they are large.

These rows are adjacent best-hit intervals merged by donor chrarm/haplotype
from `odgi untangle -e 50000 -m 0 -j 0 -n 1`, filtered to interchromosomal
merged intervals with `min_score >= 0.8` and length <= 100 kb. They are a
resolution check on the graph/untangle process, not the high-quality m1000 patch
denominator used for the manuscript's 494/538 within-community statistic.

Lower-merge overall: N = {lower_overall['n']}; min {lower_overall['min_bp']} bp;
Q1 {lower_overall['q25_bp']} bp; median {lower_overall['median_bp']} bp; Q3
{lower_overall['q75_bp']} bp; max {lower_overall['max_bp']} bp.

- In 22--95 bp conversion-like range:
  {lower_overall['conversion_22_95_n']}/{lower_overall['n']}
  ({float(lower_overall['conversion_22_95_prop']) * 100:.1f}%).
- In 318--688 bp CO-associated range:
  {lower_overall['crossover_318_688_n']}/{lower_overall['n']}
  ({float(lower_overall['crossover_318_688_prop']) * 100:.1f}%).
- Near CO-associated two-fold window:
  {lower_overall['near_crossover_159_1376_n']}/{lower_overall['n']}
  ({float(lower_overall['near_crossover_159_1376_prop']) * 100:.1f}%).

| group_type | group_value | n | min | q25 | median | q75 | q95 | max | conv 22-95 | CO 318-688 | near CO 159-1376 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
"""
        for row in lower_rows:
            text += (
                f"| {row['group_type']} | {row['group_value']} | {row['n']} | "
                f"{row['min_bp']} | {row['q25_bp']} | {row['median_bp']} | "
                f"{row['q75_bp']} | {row['q95_bp']} | {row['max_bp']} | "
                f"{row['conversion_22_95_n']} | {row['crossover_318_688_n']} | "
                f"{row['near_crossover_159_1376_n']} |\n"
            )
    path.write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--washu-patches", type=Path, default=DEFAULT_WASHU)
    parser.add_argument("--all-pedigrees", type=Path, default=DEFAULT_ALL)
    parser.add_argument("--out-tsv", type=Path, default=DEFAULT_TSV)
    parser.add_argument("--out-lower-tsv", type=Path, default=DEFAULT_LOWER_TSV)
    parser.add_argument("--out-md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--lower-untangle-dir", type=Path, default=DEFAULT_LOWER_DIR)
    parser.add_argument("--run-lower-untangle", action="store_true")
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--target-dir", type=Path, default=DEFAULT_TARGET_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    patches = load_washu_hq(args.washu_patches)
    if args.run_lower_untangle:
        run_lower_untangle(args.lower_untangle_dir, args.graph, args.target_dir)
    lower_patches = load_lower_merge_patches(args.lower_untangle_dir)
    rows = grouped_rows(patches)
    lower_rows = lower_grouped_rows(lower_patches) if lower_patches else []
    if not lower_rows:
        lower_rows = read_tsv_rows(args.out_lower_tsv)
    aggregate_count = load_aggregate_washu_count(args.all_pedigrees)
    write_tsv(rows, args.out_tsv)
    if lower_rows:
        write_tsv(lower_rows, args.out_lower_tsv)
    write_markdown(rows, lower_rows, args.out_md, args.washu_patches, args.all_pedigrees, aggregate_count)
    print(f"Wrote {args.out_tsv}")
    if lower_rows:
        print(f"Wrote {args.out_lower_tsv}")
    print(f"Wrote {args.out_md}")
    print(f"WashU high-quality interchromosomal patches: {len(patches)}")
    if lower_patches:
        print(f"Lower-merge m0/n1 merged interchromosomal intervals: {len(lower_patches)}")
    if aggregate_count is not None:
        print(f"Aggregate ds == WashU rows: {aggregate_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
