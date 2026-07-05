#!/usr/bin/env python3
"""Build the compact audit table for pedigree PHR-exchange runs."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "submission/supplementary/Supplementary_Table_1_pedigree_runs.tsv"

PATERNAL_SUMMARY = (
    ROOT
    / "paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/"
    / "pan027_paternal_high_conf_interchrom_winning_tracts.tsv"
)
PATERNAL_DONORS = (
    ROOT
    / "paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/"
    / "donor_recipient_runs.tsv"
)
PAN027_MATERNAL_RUNS = (
    ROOT
    / "paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/"
    / "PAN027mat_vs_PAN010_joint/PAN027mat_vs_PAN010_joint.whole_genome_ribbon_runs.tsv"
)
PAN028_MATERNAL_RUNS = (
    ROOT
    / "paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/"
    / "PAN028mat_vs_PAN027_joint/PAN028mat_vs_PAN027_joint.whole_genome_ribbon_runs.tsv"
)
CHM13_SIZES = ROOT / "data/chm13.chrom.sizes"
CHM13_PHR_BED = ROOT / "data/chm13.phrs.bed"
ARM_AUDIT = ROOT / "submission/fig/MainFigures/arm_inclusion_audit.tsv"

FIELDS = [
    "run_id",
    "transmission",
    "call_class",
    "query_arm",
    "donor_arm",
    "query_interval_native_0based_halfopen",
    "donor_interval_native_0based_halfopen",
    "length_bp",
    "windows_2kb",
    "mean_nonhomolog_identity",
    "mean_samechr_identity",
    "delta_identity",
    "community_context",
    "chm13_query_phr_bed_interval",
    "chm13_donor_phr_bed_interval",
    "wfmash_status",
    "source_tables",
    "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_sizes() -> dict[str, int]:
    sizes: dict[str, int] = {}
    with CHM13_SIZES.open() as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def arm_for_position(chrom: str, start: int, sizes: dict[str, int]) -> str:
    if chrom == "chrY":
        return "p" if start < 10_000_000 else "q"
    return "p" if start < (sizes[chrom] / 2) else "q"


def arm_label(chrom: str, start: int, sizes: dict[str, int]) -> str:
    return f"{chrom}_{arm_for_position(chrom, start, sizes)}"


def short_arm(arm: str) -> str:
    chrom, side = arm.split("_")
    return f"{chrom}{side}"


def interval(chrom: str, start: str | int, end: str | int) -> str:
    return f"{chrom}:{start}-{end}"


def read_communities() -> dict[str, str]:
    rows = read_tsv(ARM_AUDIT)
    return {row["ChromArm"]: row["Community"] for row in rows}


def read_phr_bed(sizes: dict[str, int]) -> dict[str, str]:
    intervals: dict[str, list[str]] = {}
    with CHM13_PHR_BED.open() as handle:
        for line in handle:
            chrom, start_s, end_s, *_rest = line.rstrip("\n").split("\t")
            start = int(start_s)
            arm = arm_label(chrom, start, sizes)
            intervals.setdefault(arm, []).append(interval(chrom, start_s, end_s))
    phr = {arm: ";".join(vals) for arm, vals in intervals.items()}
    # The CHM13 PHR BED stores the PAR1 annotation on chrX.
    phr.setdefault("chrY_p", "PAR1 represented by chrX:1579-501579")
    return phr


def community_context(query_arm: str, donor_arm: str, communities: dict[str, str]) -> str:
    query_comm = communities.get(query_arm, "NA")
    donor_comm = communities.get(donor_arm, "NA")
    return f"{short_arm(query_arm)} {query_comm}; {short_arm(donor_arm)} {donor_comm}"


def phr_interval(arm: str, phr: dict[str, str]) -> str:
    return phr.get(arm, "no called interval in data/chm13.phrs.bed")


def paternal_summary_row(rows: list[dict[str, str]], query_chrom: str, target_chrom: str, start: int, end: int) -> dict[str, str]:
    for row in rows:
        if (
            row["query_chrom"] == query_chrom
            and row["target_chrom"] == target_chrom
            and int(row["query_start"]) == start
            and int(row["query_end"]) == end
        ):
            return row
    raise ValueError(f"missing paternal summary row: {query_chrom}:{start}-{end} -> {target_chrom}")


def paternal_donor_row(rows: list[dict[str, str]], query_chrom: str, target_chrom: str, start: int, end: int) -> dict[str, str] | None:
    for row in rows:
        if (
            row["query_chrom"] == query_chrom
            and row["target_chrom"] == target_chrom
            and int(row["query_start"]) == start
            and int(row["query_end"]) == end
        ):
            return row
    return None


def add_paternal(
    out_rows: list[dict[str, str]],
    *,
    run_id: str,
    call_class: str,
    query_chrom: str,
    target_chrom: str,
    start: int,
    end: int,
    wfmash_status: str,
    note: str,
    summary_rows: list[dict[str, str]],
    donor_rows: list[dict[str, str]],
    sizes: dict[str, int],
    communities: dict[str, str],
    phr: dict[str, str],
) -> None:
    summary = paternal_summary_row(summary_rows, query_chrom, target_chrom, start, end)
    donor = paternal_donor_row(donor_rows, query_chrom, target_chrom, start, end)
    query_arm = arm_label(query_chrom, start, sizes)
    if donor is None:
        donor_arm = arm_label(target_chrom, start, sizes) if target_chrom != "chrY" else "chrY_p"
        donor_interval = "not captured in donor-recipient zoom table"
    else:
        donor_start = int(donor["donor_start"])
        donor_end = int(donor["donor_end"])
        donor_arm = arm_label(target_chrom, donor_start, sizes) if target_chrom != "chrY" else "chrY_p"
        donor_interval = interval(donor["donor_seq"], donor_start, donor_end)
    out_rows.append(
        {
            "run_id": run_id,
            "transmission": "PAN027 paternal haplotype (child) vs PAN011 father joint haplotypes",
            "call_class": call_class,
            "query_arm": query_arm,
            "donor_arm": donor_arm,
            "query_interval_native_0based_halfopen": interval(query_chrom, start, end),
            "donor_interval_native_0based_halfopen": donor_interval,
            "length_bp": summary["bp"],
            "windows_2kb": summary["windows"],
            "mean_nonhomolog_identity": summary["mean_inter_identity"],
            "mean_samechr_identity": summary["mean_same_identity"],
            "delta_identity": summary["mean_delta"],
            "community_context": community_context(query_arm, donor_arm, communities),
            "chm13_query_phr_bed_interval": phr_interval(query_arm, phr),
            "chm13_donor_phr_bed_interval": phr_interval(donor_arm, phr),
            "wfmash_status": wfmash_status,
            "source_tables": f"{PATERNAL_SUMMARY.relative_to(ROOT)}; {PATERNAL_DONORS.relative_to(ROOT)}",
            "notes": note,
        }
    )


def add_maternal(
    out_rows: list[dict[str, str]],
    *,
    run_id: str,
    transmission: str,
    call_class: str,
    source_path: Path,
    selector: tuple[str, str] | None = None,
    query_chrom: str | None = None,
    target_chrom: str | None = None,
    note: str,
    sizes: dict[str, int],
    communities: dict[str, str],
    phr: dict[str, str],
) -> None:
    for row in read_tsv(source_path):
        if selector is not None:
            key, value = selector
            if row[key] != value:
                continue
        if query_chrom is not None and row["query_chrom"] != query_chrom:
            continue
        if target_chrom is not None and row["target_chrom"] != target_chrom:
            continue
        query_start = int(row["query_start"])
        donor_start = int(row["donor_start"])
        query_arm = arm_label(row["query_chrom"], query_start, sizes)
        donor_arm = arm_label(row["target_chrom"], donor_start, sizes)
        out_rows.append(
            {
                "run_id": run_id,
                "transmission": transmission,
                "call_class": call_class,
                "query_arm": query_arm,
                "donor_arm": donor_arm,
                "query_interval_native_0based_halfopen": interval(row["query_chrom"], row["query_start"], row["query_end"]),
                "donor_interval_native_0based_halfopen": interval(row["donor_seq"], row["donor_start"], row["donor_end"]),
                "length_bp": row["bp"],
                "windows_2kb": row["windows"],
                "mean_nonhomolog_identity": row["mean_inter_identity"],
                "mean_samechr_identity": row["mean_same_identity"],
                "delta_identity": row["delta_identity"],
                "community_context": community_context(query_arm, donor_arm, communities),
                "chm13_query_phr_bed_interval": phr_interval(query_arm, phr),
                "chm13_donor_phr_bed_interval": phr_interval(donor_arm, phr),
                "wfmash_status": "not separately evaluated with wfmash",
                "source_tables": str(source_path.relative_to(ROOT)),
                "notes": note,
            }
        )
        return
    raise ValueError(f"no maternal row matched {run_id}")


def main() -> None:
    sizes = read_sizes()
    communities = read_communities()
    phr = read_phr_bed(sizes)
    paternal_summary_rows = read_tsv(PATERNAL_SUMMARY)
    paternal_donor_rows = read_tsv(PATERNAL_DONORS)
    rows: list[dict[str, str]] = []

    add_paternal(
        rows,
        run_id="PAN027pat_PAR1_run1",
        call_class="PAR1 positive-control non-homologous winner",
        query_chrom="chrX",
        target_chrom="chrY",
        start=14000,
        end=106000,
        wfmash_status="not separately evaluated with wfmash",
        note="Obligate male PAR1 recombination control; all coordinates are assembly-native, not CHM13-lifted.",
        summary_rows=paternal_summary_rows,
        donor_rows=paternal_donor_rows,
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_paternal(
        rows,
        run_id="PAN027pat_PAR1_run2",
        call_class="PAR1 positive-control non-homologous winner",
        query_chrom="chrX",
        target_chrom="chrY",
        start=110000,
        end=156000,
        wfmash_status="not separately evaluated with wfmash",
        note="Second adjacent PAR1 positive-control run; grouped only by exact end-to-end donor/query continuity.",
        summary_rows=paternal_summary_rows,
        donor_rows=paternal_donor_rows,
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_paternal(
        rows,
        run_id="PAN027pat_chr5q_chr1p_candidate",
        call_class="putative autosomal PHR exchange",
        query_chrom="chr5",
        target_chrom="chr1",
        start=182052000,
        end=182080000,
        wfmash_status="not recovered by wfmash; wfmash places this window on chr5",
        note="Weaker aligner-dependent autosomal candidate; identity margin is narrow.",
        summary_rows=paternal_summary_rows,
        donor_rows=paternal_donor_rows,
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_paternal(
        rows,
        run_id="PAN027pat_chr9q_chr3q_candidate",
        call_class="putative autosomal PHR exchange",
        query_chrom="chr9",
        target_chrom="chr3",
        start=136168000,
        end=136188000,
        wfmash_status="concordant with wfmash -p 95 at >99.8% identity on chr3q",
        note="Autosomal candidate concordant across SweepGA/FastGA and wfmash; no ordinary chr3 crossover in manual WashU annotation.",
        summary_rows=paternal_summary_rows,
        donor_rows=paternal_donor_rows,
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_paternal(
        rows,
        run_id="PAN027pat_representative_acrocentric_run",
        call_class="representative acrocentric context run",
        query_chrom="chr13",
        target_chrom="chr21",
        start=7924000,
        end=7968000,
        wfmash_status="not separately evaluated with wfmash",
        note="Representative longest paternal acrocentric-context row; full acrocentric outputs remain in source run tables.",
        summary_rows=paternal_summary_rows,
        donor_rows=paternal_donor_rows,
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_maternal(
        rows,
        run_id="PAN027mat_representative_acrocentric_run",
        transmission="PAN027 maternal haplotype (child) vs PAN010 mother joint haplotypes",
        call_class="representative maternal acrocentric context run",
        source_path=PAN027_MATERNAL_RUNS,
        selector=("category", "acro_acro"),
        note="Longest displayed PAN027 maternal non-homologous run; maternal scan otherwise resolves acrocentric-acrocentric context.",
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_maternal(
        rows,
        run_id="PAN028mat_representative_acrocentric_run",
        transmission="PAN028 maternal haplotype (child) vs PAN027 mother joint haplotypes",
        call_class="representative maternal acrocentric context run",
        source_path=PAN028_MATERNAL_RUNS,
        selector=("category", "acro_acro"),
        note="Longest displayed PAN028 maternal acrocentric-acrocentric run.",
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_maternal(
        rows,
        run_id="PAN028mat_chr21p_chr4p_acro_other_run",
        transmission="PAN028 maternal haplotype (child) vs PAN027 mother joint haplotypes",
        call_class="maternal acrocentric-other context run",
        source_path=PAN028_MATERNAL_RUNS,
        selector=("category", "acro_other"),
        note="Longest displayed acrocentric-other maternal run.",
        sizes=sizes,
        communities=communities,
        phr=phr,
    )
    add_maternal(
        rows,
        run_id="PAN028mat_chr5q_chr1p_candidate",
        transmission="PAN028 maternal haplotype (child) vs PAN027 mother joint haplotypes",
        call_class="putative maternal chr5/chr1 PHR exchange context",
        source_path=PAN028_MATERNAL_RUNS,
        selector=("category", "chr5_chr1_candidate"),
        note="Displayed maternal chr5/chr1 candidate context; not treated as a primary paternal Fig. 5 claim.",
        sizes=sizes,
        communities=communities,
        phr=phr,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
