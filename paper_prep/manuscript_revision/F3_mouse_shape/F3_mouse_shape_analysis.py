#!/usr/bin/env python3
"""Mouse meiosis zygotene curve-shape contrast for manuscript revision F-3.

This script uses cached mouse PHR/contact TSVs already vendored in the repo.
It does not rebuild Hi-C matrices or touch the upstream /moosefs source tree.

Outputs are written next to this script:
  - F3_stage_series.tsv
  - F3_zygotene_contrasts.tsv
  - F3_input_inventory.tsv

The direct test is the paired zygotene-minus-comparator contrast in Spearman
rho, not the per-stage rho != 0 test. Resampling is clustered by mouse arm pair
for sequence-level rows and by mouse arm for the arm-pair/community-collapsed
rows, because pair entries that share arms are not independent.

The implementation intentionally uses only the Python standard library so the
revision artifact can be rerun without a conda/guix environment.
"""

from __future__ import annotations

import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent

STAGES = ("leptotene", "zygotene", "pachytene", "diplotene")
COMPARATORS = ("leptotene", "pachytene", "diplotene")
WINDOW = "1Mb"
RES = "50000bp"
SEED = 2463
N_BOOT = 2000

SEQLEVEL_TEMPLATE = (
    ROOT
    / "data/mouse_meiosis_sweep/seqlevel/1Mb/mouse_{stage}_phr_50000bp_seqlevel.tsv"
)
ARM_PAIR_TEMPLATE = (
    ROOT
    / "data/mouse_meiosis_sweep/zuo/1Mb/50000bp/zuo2021_{stage}_phr_pair_correlation.tsv"
)


def norm_pair(a: str, b: str) -> str:
    return "||".join(sorted((str(a), str(b))))


def parse_float(value: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out


def rank_average(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [float("nan")] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = avg_rank
        i = j
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 4:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    sx = math.sqrt(sum(v * v for v in dx))
    sy = math.sqrt(sum(v * v for v in dy))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum(a * b for a, b in zip(dx, dy)) / (sx * sy)


def spearman_from_rows(rows: list[dict[str, object]], x_col: str, y_col: str) -> float:
    x: list[float] = []
    y: list[float] = []
    for row in rows:
        xv = float(row[x_col])
        yv = float(row[y_col])
        if math.isfinite(xv) and math.isfinite(yv):
            x.append(xv)
            y.append(yv)
    if len(x) < 4:
        return float("nan")
    return pearson(rank_average(x), rank_average(y))


def load_seqlevel(stage: str) -> list[dict[str, object]]:
    path = Path(str(SEQLEVEL_TEMPLATE).format(stage=stage))
    rows: list[dict[str, object]] = []
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["stage"] = stage
            row["seq_pair"] = norm_pair(row["seq_a"], row["seq_b"])
            row["arm_pair"] = norm_pair(row["arm_a"], row["arm_b"])
            row["jaccard"] = parse_float(row["jaccard"])
            row["hic_contact_norm"] = parse_float(row["hic_contact_norm"])
            rows.append(row)
    return rows


def load_arm_pair(stage: str) -> list[dict[str, object]]:
    path = Path(str(ARM_PAIR_TEMPLATE).format(stage=stage))
    rows: list[dict[str, object]] = []
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["stage"] = stage
            row["arm_pair"] = norm_pair(row["arm_a"], row["arm_b"])
            row["mean_jaccard"] = parse_float(row["mean_jaccard"])
            row["hic_contact"] = parse_float(row["hic_contact"])
            rows.append(row)
    return rows


def count_unique(rows: list[dict[str, object]], col: str) -> int:
    return len({str(row[col]) for row in rows})


def bootstrap_clusters(
    base: list[dict[str, object]],
    cluster_col: str,
    stat_fn: Callable[[list[dict[str, object]]], float],
    n_boot: int,
    seed: int,
) -> list[float]:
    rng = random.Random(seed)
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in base:
        grouped[str(row[cluster_col])].append(row)
    clusters = sorted(grouped)
    out: list[float] = []
    for _ in range(n_boot):
        sample: list[dict[str, object]] = []
        for cluster in (rng.choice(clusters) for _ in clusters):
            sample.extend(grouped[cluster])
        value = stat_fn(sample)
        if math.isfinite(value):
            out.append(value)
    return out


def bootstrap_arms(
    base: list[dict[str, object]],
    stat_fn: Callable[[list[dict[str, object]]], float],
    n_boot: int,
    seed: int,
) -> list[float]:
    rng = random.Random(seed)
    arms = sorted({str(row["arm_a"]) for row in base} | {str(row["arm_b"]) for row in base})
    pair_lookup = {frozenset((str(row["arm_a"]), str(row["arm_b"]))): row for row in base}
    out: list[float] = []
    for _ in range(n_boot):
        sampled = [rng.choice(arms) for _ in arms]
        rows: list[dict[str, object]] = []
        for left in range(len(sampled)):
            for right in range(left + 1, len(sampled)):
                if sampled[left] == sampled[right]:
                    continue
                row = pair_lookup.get(frozenset((sampled[left], sampled[right])))
                if row is not None:
                    rows.append(row)
        value = stat_fn(rows)
        if math.isfinite(value):
            out.append(value)
    return out


def quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    pos = q * (len(sorted_values) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_values[lo]
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def summarize_bootstrap(dist: list[float]) -> tuple[float, float, float]:
    if not dist:
        return float("nan"), float("nan"), float("nan")
    ordered = sorted(dist)
    ci_lo = quantile(ordered, 0.025)
    ci_hi = quantile(ordered, 0.975)
    p_two = 2 * min(
        sum(1 for value in dist if value <= 0) / len(dist),
        sum(1 for value in dist if value >= 0) / len(dist),
    )
    return ci_lo, ci_hi, min(1.0, p_two)


def stage_series_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for stage in STAGES:
        seq = load_seqlevel(stage)
        arm = load_arm_pair(stage)
        rows.append(
            {
                "analysis": "sequence_level_exact_phr",
                "window": WINDOW,
                "resolution": RES,
                "stage": stage,
                "rho": spearman_from_rows(seq, "jaccard", "hic_contact_norm"),
                "n_rows": len(seq),
                "n_clusters": count_unique(seq, "arm_pair"),
                "cluster_unit": "arm_pair",
                "x": "per-PHR-pair Jaccard",
                "y": "Hi-C normalized contact at exact PHR coordinates",
                "source_path": str(SEQLEVEL_TEMPLATE).format(stage=stage),
            }
        )
        rows.append(
            {
                "analysis": "arm_pair_collapsed",
                "window": WINDOW,
                "resolution": RES,
                "stage": stage,
                "rho": spearman_from_rows(arm, "mean_jaccard", "hic_contact"),
                "n_rows": len(arm),
                "n_clusters": len({str(row["arm_a"]) for row in arm} | {str(row["arm_b"]) for row in arm}),
                "cluster_unit": "arm",
                "x": "mean Jaccard per unordered arm pair",
                "y": "community-analysis Hi-C contact per unordered arm pair",
                "source_path": str(ARM_PAIR_TEMPLATE).format(stage=stage),
            }
        )
    return rows


def seqlevel_contrast(comparator: str, seed_offset: int) -> dict[str, object]:
    z_lookup = {str(row["seq_pair"]): row for row in load_seqlevel("zygotene")}
    c_lookup = {str(row["seq_pair"]): row for row in load_seqlevel(comparator)}
    paired: list[dict[str, object]] = []
    for key in sorted(set(z_lookup) & set(c_lookup)):
        z = z_lookup[key]
        c = c_lookup[key]
        paired.append(
            {
                "seq_pair": key,
                "arm_pair": z["arm_pair"],
                "jaccard_zygotene": z["jaccard"],
                "hic_contact_norm_zygotene": z["hic_contact_norm"],
                f"jaccard_{comparator}": c["jaccard"],
                f"hic_contact_norm_{comparator}": c["hic_contact_norm"],
            }
        )

    def delta(rows: list[dict[str, object]]) -> float:
        rho_z = spearman_from_rows(rows, "jaccard_zygotene", "hic_contact_norm_zygotene")
        rho_c = spearman_from_rows(rows, f"jaccard_{comparator}", f"hic_contact_norm_{comparator}")
        return rho_z - rho_c

    rho_z = spearman_from_rows(paired, "jaccard_zygotene", "hic_contact_norm_zygotene")
    rho_c = spearman_from_rows(paired, f"jaccard_{comparator}", f"hic_contact_norm_{comparator}")
    boot = bootstrap_clusters(paired, "arm_pair", delta, N_BOOT, SEED + seed_offset)
    ci_lo, ci_hi, p_two = summarize_bootstrap(boot)
    return {
        "analysis": "sequence_level_exact_phr",
        "contrast": f"zygotene_minus_{comparator}",
        "rho_zygotene_common_rows": rho_z,
        "rho_comparator_common_rows": rho_c,
        "delta_rho": rho_z - rho_c,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_two_bootstrap": p_two,
        "n_rows_common": len(paired),
        "n_clusters": count_unique(paired, "arm_pair"),
        "cluster_unit": "arm_pair",
        "n_boot_valid": len(boot),
        "interpretation": "direct paired contrast; positive means zygotene has higher rho",
    }


def arm_pair_contrast(comparator: str, seed_offset: int) -> dict[str, object]:
    z_lookup = {str(row["arm_pair"]): row for row in load_arm_pair("zygotene")}
    c_lookup = {str(row["arm_pair"]): row for row in load_arm_pair(comparator)}
    paired: list[dict[str, object]] = []
    for key in sorted(set(z_lookup) & set(c_lookup)):
        z = z_lookup[key]
        c = c_lookup[key]
        paired.append(
            {
                "arm_pair": key,
                "arm_a": z["arm_a"],
                "arm_b": z["arm_b"],
                "mean_jaccard_zygotene": z["mean_jaccard"],
                "hic_contact_zygotene": z["hic_contact"],
                f"mean_jaccard_{comparator}": c["mean_jaccard"],
                f"hic_contact_{comparator}": c["hic_contact"],
            }
        )

    def delta(rows: list[dict[str, object]]) -> float:
        rho_z = spearman_from_rows(rows, "mean_jaccard_zygotene", "hic_contact_zygotene")
        rho_c = spearman_from_rows(rows, f"mean_jaccard_{comparator}", f"hic_contact_{comparator}")
        return rho_z - rho_c

    rho_z = spearman_from_rows(paired, "mean_jaccard_zygotene", "hic_contact_zygotene")
    rho_c = spearman_from_rows(paired, f"mean_jaccard_{comparator}", f"hic_contact_{comparator}")
    boot = bootstrap_arms(paired, delta, N_BOOT, SEED + seed_offset)
    ci_lo, ci_hi, p_two = summarize_bootstrap(boot)
    return {
        "analysis": "arm_pair_collapsed",
        "contrast": f"zygotene_minus_{comparator}",
        "rho_zygotene_common_rows": rho_z,
        "rho_comparator_common_rows": rho_c,
        "delta_rho": rho_z - rho_c,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "p_two_bootstrap": p_two,
        "n_rows_common": len(paired),
        "n_clusters": len({str(row["arm_a"]) for row in paired} | {str(row["arm_b"]) for row in paired}),
        "cluster_unit": "arm",
        "n_boot_valid": len(boot),
        "interpretation": "direct paired contrast; positive means zygotene has higher rho",
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def input_inventory() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for stage in STAGES:
        for analysis, template in (
            ("sequence_level_exact_phr", SEQLEVEL_TEMPLATE),
            ("arm_pair_collapsed", ARM_PAIR_TEMPLATE),
        ):
            path = Path(str(template).format(stage=stage))
            rows.append(
                {
                    "analysis": analysis,
                    "stage": stage,
                    "path": str(path),
                    "exists": path.exists(),
                    "bytes": path.stat().st_size if path.exists() else "NA",
                    "line_count": sum(1 for _ in path.open()) if path.exists() else "NA",
                }
            )
    return rows


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT_DIR / "F3_input_inventory.tsv", input_inventory())
    write_tsv(OUT_DIR / "F3_stage_series.tsv", stage_series_rows())

    contrast_rows: list[dict[str, object]] = []
    for i, comparator in enumerate(COMPARATORS):
        contrast_rows.append(seqlevel_contrast(comparator, seed_offset=100 * (i + 1)))
        contrast_rows.append(arm_pair_contrast(comparator, seed_offset=1000 + 100 * (i + 1)))
    write_tsv(OUT_DIR / "F3_zygotene_contrasts.tsv", contrast_rows)

    print(f"wrote {OUT_DIR / 'F3_input_inventory.tsv'}")
    print(f"wrote {OUT_DIR / 'F3_stage_series.tsv'}")
    print(f"wrote {OUT_DIR / 'F3_zygotene_contrasts.tsv'}")
    print("Slurm: not used; cached TSV inputs are small and resampling is local.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
