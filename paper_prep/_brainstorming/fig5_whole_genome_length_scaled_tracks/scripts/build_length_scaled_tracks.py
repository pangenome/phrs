#!/usr/bin/env python3
"""Build length-scaled whole-genome track tables for Fig5.

Inputs are already generated Fig5 products.  This script does not run aligners
or SweepGA; it reshapes strict untangle segments and the 2 kb query-grid 1:1
alignment support table into a common chromosome-row plotting format.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = Path(__file__).resolve().parents[1]

UNTANGLE = ROOT / "paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_segments.tsv"
BINNED = ROOT / "paper_prep/_brainstorming/fig5_whole_genome_alignment_overview/whole_genome_binned_support.tsv"
F32_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv"
WFMASH_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHR_RANK = {chrom: idx for idx, chrom in enumerate(CHR_ORDER)}

PAGES = {
    "PAN027_vs_PAN010": {
        "page_id": "PAN027_hap1_maternal_vs_PAN010",
        "page_order": 1,
        "comparison_id": "PAN027mat_vs_PAN010_joint",
        "query_sample": "PAN027",
        "query_haplotype": "PAN027#1",
        "page_label": "PAN027 maternal hap1 <- PAN010 mother",
    },
    "PAN027_vs_PAN011": {
        "page_id": "PAN027_hap2_paternal_vs_PAN011",
        "page_order": 2,
        "comparison_id": "PAN027pat_vs_PAN011_joint",
        "query_sample": "PAN027",
        "query_haplotype": "PAN027#2",
        "page_label": "PAN027 paternal hap2 <- PAN011 father",
    },
    "PAN028_vs_PAN027": {
        "page_id": "PAN028_hap1_maternal_vs_PAN027",
        "page_order": 3,
        "comparison_id": "PAN028mat_vs_PAN027_joint",
        "query_sample": "PAN028",
        "query_haplotype": "PAN028#1",
        "page_label": "PAN028 maternal hap1 <- PAN027 mother",
    },
}
COMPARISON_TO_PAGE = {value["comparison_id"]: value for value in PAGES.values()}

METHODS = {
    "untangle_strict_primary_path": {
        "method_order": 1,
        "method_family": "untangle",
        "method_label": "Untangle strict primary path",
        "source_detail": "strict primary-path untangle segments; no chopping",
    },
    "wfmash_p95_qgrid2kb_1to1_ani": {
        "method_order": 2,
        "method_family": "wfmash_p95_updated_bin",
        "method_label": "wfmash -p95 2 kb query-grid -> SweepGA 1:1 ANI",
        "source_detail": "wfmash -p95 updated binary; 2 kb query-grid chopped PAF; SweepGA --num-mappings 1:1 --scoring ani --scaffold-jump 0",
    },
    "sweepga_fastga_f32_qgrid2kb_1to1_ani": {
        "method_order": 3,
        "method_family": "sweepga_fastga_f32",
        "method_label": "SweepGA/FastGA -f32 2 kb query-grid -> 1:1 ANI",
        "source_detail": "SweepGA/FastGA --fastga-frequency 32 many:many -j 0; 2 kb query-grid chopped PAF; SweepGA --num-mappings 1:1 --scoring ani --scaffold-jump 0",
    },
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def chrom_key(chrom: str) -> tuple[int, str]:
    return (CHR_RANK.get(chrom, 999), chrom)


def as_int(value: str) -> int:
    if value == "":
        return 0
    return int(float(value))


def add_chrom_length(chrom_lengths: dict[tuple[str, str], int], page_id: str, chrom: str, length: int) -> None:
    if chrom in CHR_RANK and length > 0:
        key = (page_id, chrom)
        chrom_lengths[key] = max(chrom_lengths.get(key, 0), length)


def common_segment_fields() -> list[str]:
    return [
        "page_order",
        "page_id",
        "page_label",
        "query_sample",
        "query_haplotype",
        "comparison_id",
        "method_order",
        "method_id",
        "method_family",
        "method_label",
        "source_detail",
        "query_chrom",
        "query_chrom_length",
        "segment_start",
        "segment_end",
        "segment_length_bp",
        "display_state",
        "display_target_chrom",
        "display_support_bp",
        "support_fraction",
        "mean_identity",
        "source_row_type",
        "winner_target_chrom",
        "top_interchrom_target_chrom",
        "callout_event_id",
        "callout_label",
    ]


def page_values_from_pair(pair: str) -> dict[str, object]:
    if pair not in PAGES:
        raise KeyError(f"unknown untangle pair: {pair}")
    return PAGES[pair]


def page_values_from_comparison(comparison_id: str) -> dict[str, object] | None:
    return COMPARISON_TO_PAGE.get(comparison_id)


def add_untangle_segments(segments: list[dict[str, object]], chrom_lengths: dict[tuple[str, str], int]) -> None:
    method_id = "untangle_strict_primary_path"
    method = METHODS[method_id]
    for row in read_tsv(UNTANGLE):
        page = page_values_from_pair(row["pair"])
        qchrom = row["query_chrom"]
        qlen = max(as_int(row["query_source_end_0based_exclusive"]), as_int(row["native_query_end_0based_exclusive"]))
        add_chrom_length(chrom_lengths, page["page_id"], qchrom, qlen)
        start = as_int(row["native_query_start_0based"])
        end = as_int(row["native_query_end_0based_exclusive"])
        length = max(0, end - start)
        inter = row["interchromosomal"] == "1"
        target = row["target_chrom"] if inter else "same_chromosome"
        segments.append(
            {
                **page,
                "method_id": method_id,
                **method,
                "query_chrom": qchrom,
                "query_chrom_length": qlen,
                "segment_start": start,
                "segment_end": end,
                "segment_length_bp": length,
                "display_state": "interchromosomal" if inter else "same_chromosome",
                "display_target_chrom": target,
                "display_support_bp": length,
                "support_fraction": "1.000000",
                "mean_identity": row["identity"],
                "source_row_type": "untangle_segment",
                "winner_target_chrom": row["target_chrom"],
                "top_interchrom_target_chrom": row["target_chrom"] if inter else "none",
                "callout_event_id": row.get("event_id", ""),
                "callout_label": row.get("event_label", ""),
            }
        )


def add_binned_alignment_segments(segments: list[dict[str, object]], chrom_lengths: dict[tuple[str, str], int]) -> None:
    family_to_method = {
        "wfmash_p95_updated_bin": "wfmash_p95_qgrid2kb_1to1_ani",
        "sweepga_fastga_f32": "sweepga_fastga_f32_qgrid2kb_1to1_ani",
    }
    for row in read_tsv(BINNED):
        if row["method_family"] not in family_to_method or row["chop_length_bp"] != "2000":
            continue
        page = page_values_from_comparison(row["comparison_id"])
        if page is None:
            continue
        method_id = family_to_method[row["method_family"]]
        method = METHODS[method_id]
        qchrom = row["query_chrom"]
        qlen = as_int(row["query_chrom_length"])
        add_chrom_length(chrom_lengths, page["page_id"], qchrom, qlen)
        start = as_int(row["query_bin_start"])
        end = as_int(row["query_bin_end"])
        top_inter_bp = as_int(row["top_interchrom_support_bp"])
        winner_bp = as_int(row["winner_support_bp"])
        if row["no_support"] == "yes":
            state = "no_support"
            target = "no_support"
            support_bp = 0
        elif top_inter_bp > 0 and row["top_interchrom_target_chrom"] != "none":
            state = "interchromosomal"
            target = row["top_interchrom_target_chrom"]
            support_bp = top_inter_bp
        elif row["winner_target_family"] == "same_chromosome":
            state = "same_chromosome"
            target = "same_chromosome"
            support_bp = winner_bp
        else:
            state = "interchromosomal"
            target = row["winner_target_chrom"]
            support_bp = winner_bp
        segments.append(
            {
                **page,
                "method_id": method_id,
                **method,
                "query_chrom": qchrom,
                "query_chrom_length": qlen,
                "segment_start": start,
                "segment_end": min(end, qlen),
                "segment_length_bp": max(0, min(end, qlen) - start),
                "display_state": state,
                "display_target_chrom": target,
                "display_support_bp": support_bp,
                "support_fraction": row["support_fraction"],
                "mean_identity": row["winner_mean_identity"],
                "source_row_type": "one_mb_display_bin_from_2kb_query_grid_filtered_paf",
                "winner_target_chrom": row["winner_target_chrom"],
                "top_interchrom_target_chrom": row["top_interchrom_target_chrom"],
                "callout_event_id": row.get("callout_event_id", ""),
                "callout_label": row.get("callout_label", ""),
            }
        )


def build_chromosome_rows(chrom_lengths: dict[tuple[str, str], int]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    page_by_id = {value["page_id"]: value for value in PAGES.values()}
    for page_id, page in sorted(page_by_id.items(), key=lambda item: item[1]["page_order"]):
        for chrom in CHR_ORDER:
            length = chrom_lengths.get((page_id, chrom), 0)
            rows.append(
                {
                    "page_order": page["page_order"],
                    "page_id": page_id,
                    "page_label": page["page_label"],
                    "query_sample": page["query_sample"],
                    "query_haplotype": page["query_haplotype"],
                    "comparison_id": page["comparison_id"],
                    "query_chrom": chrom,
                    "query_chrom_order": CHR_RANK[chrom] + 1,
                    "query_chrom_length": length,
                }
            )
    return rows


def build_summary_rows(segments: list[dict[str, object]]) -> list[dict[str, object]]:
    totals: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
    row_counts: Counter[tuple[str, str, str, str]] = Counter()
    for row in segments:
        if row["display_state"] != "interchromosomal":
            continue
        key = (str(row["page_id"]), str(row["method_id"]), str(row["query_chrom"]), str(row["display_target_chrom"]))
        totals[key]["display_support_bp"] += int(row["display_support_bp"])
        totals[key]["segment_length_bp"] += int(row["segment_length_bp"])
        row_counts[key] += 1
    page_by_id = {value["page_id"]: value for value in PAGES.values()}
    out = []
    for key, counts in totals.items():
        page_id, method_id, qchrom, target = key
        page = page_by_id[page_id]
        method = METHODS[method_id]
        out.append(
            {
                "page_order": page["page_order"],
                "page_id": page_id,
                "page_label": page["page_label"],
                "method_order": method["method_order"],
                "method_id": method_id,
                "method_label": method["method_label"],
                "query_chrom": qchrom,
                "target_chrom": target,
                "display_support_bp": counts["display_support_bp"],
                "covered_bp": counts["segment_length_bp"],
                "row_count": row_counts[key],
            }
        )
    out.sort(key=lambda row: (int(row["page_order"]), int(row["method_order"]), chrom_key(str(row["query_chrom"])), -int(row["display_support_bp"]), str(row["target_chrom"])))
    return out


def build_manifest_rows() -> list[dict[str, object]]:
    return [
        {
            "artifact": "strict_untangle_segments",
            "path": UNTANGLE,
            "exists": UNTANGLE.exists(),
            "role": "untangle strict primary-path source",
            "settings": "no chopping; native query coordinates",
        },
        {
            "artifact": "whole_genome_binned_support",
            "path": BINNED,
            "exists": BINNED.exists(),
            "role": "common 1 Mb display-bin support table for filtered alignment PAFs",
            "settings": "filtered to method_family in {wfmash_p95_updated_bin,sweepga_fastga_f32} and chop_length_bp=2000",
        },
        {
            "artifact": "sweepga_f32_query_grid_manifest",
            "path": F32_MANIFEST,
            "exists": F32_MANIFEST.exists(),
            "role": "SweepGA/FastGA f32 2kb query-grid 1:1 source manifest",
            "settings": "raw whole-genome FASTA -> SweepGA/FastGA --fastga-frequency 32 many:many -j 0 -> pafchop query-grid 2kb -> SweepGA 1:1 ANI",
        },
        {
            "artifact": "wfmash_query_grid_manifest",
            "path": WFMASH_MANIFEST,
            "exists": WFMASH_MANIFEST.exists(),
            "role": "wfmash p95 2kb query-grid 1:1 source manifest",
            "settings": "whole-genome wfmash -p95 -> pafchop query-grid 2kb -> SweepGA 1:1 ANI",
        },
    ]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    segments: list[dict[str, object]] = []
    chrom_lengths: dict[tuple[str, str], int] = {}
    add_untangle_segments(segments, chrom_lengths)
    add_binned_alignment_segments(segments, chrom_lengths)
    segments.sort(
        key=lambda row: (
            int(row["page_order"]),
            int(row["method_order"]),
            chrom_key(str(row["query_chrom"])),
            int(row["segment_start"]),
            int(row["segment_end"]),
            str(row["display_target_chrom"]),
        )
    )
    chromosome_rows = build_chromosome_rows(chrom_lengths)
    summary_rows = build_summary_rows(segments)
    manifest_rows = build_manifest_rows()

    write_tsv(OUT_DIR / "length_scaled_track_segments.tsv", common_segment_fields(), segments)
    write_tsv(
        OUT_DIR / "length_scaled_track_chromosomes.tsv",
        [
            "page_order",
            "page_id",
            "page_label",
            "query_sample",
            "query_haplotype",
            "comparison_id",
            "query_chrom",
            "query_chrom_order",
            "query_chrom_length",
        ],
        chromosome_rows,
    )
    write_tsv(
        OUT_DIR / "length_scaled_track_summary.tsv",
        [
            "page_order",
            "page_id",
            "page_label",
            "method_order",
            "method_id",
            "method_label",
            "query_chrom",
            "target_chrom",
            "display_support_bp",
            "covered_bp",
            "row_count",
        ],
        summary_rows,
    )
    write_tsv(OUT_DIR / "length_scaled_track_manifest.tsv", ["artifact", "path", "exists", "role", "settings"], manifest_rows)
    print(f"wrote {len(segments)} segments, {len(chromosome_rows)} chromosome rows, {len(summary_rows)} interchromosomal summary rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
