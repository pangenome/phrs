#!/usr/bin/env python3
"""Build the compact, paper-facing Extended Data Figure 2 tables.

The supported-term table is an exact field-preserving subset of the frozen V7
result ledger.  The arm table is the six-class post-inference display matrix in
the same community/tree order used by ED2.  No statistical inference occurs in
this script.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
V7 = REPO / "paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7"
HEATMAP = V7 / "community_attribution/functional_component_heatmap"
OUT = REPO / "submission/supporting_material"

TERM_COLUMNS = (
    "schema_version",
    "assignment",
    "hypothesis_index",
    "collection",
    "relation",
    "ontology",
    "namespace",
    "term_id",
    "term_name",
    "hypothesis_status",
    "a_T",
    "b_T",
    "c_T",
    "d_T",
    "M_T_genome_burden",
    "K_phr_total",
    "nonphr_total",
    "N_eligible_total",
    "phr_copy_fraction",
    "nonphr_copy_fraction",
    "fold_enrichment",
    "fold_enrichment_ha",
    "odds_ratio",
    "odds_ratio_ha",
    "expected_phr_copy_count",
    "copy_excess_over_expected",
    "p_exact_upper",
    "bh_q_within_collection",
    "holm_p_global",
    "primary_support",
    "support_status",
    "primary_support_rule",
    "finite_population_definition",
    "background_definition",
    "inference_status",
)

CLASS_ORDER = (
    "DUX4_ZGA_TRANSCRIPTION_NUCLEAR_ENVELOPE_CELL_CYCLE",
    "WASH_ENDOSOMAL_ACTIN_EXOCYST",
    "DDX11_HELICASE_CHROMOSOME",
    "SEPTIN14_SEPTIN_CYTOKINESIS",
    "WBP1L_CXCL12_SIGNALING",
    "RPL23A_RIBOSOMAL_NUCLEOLAR",
)

ARM_COLUMNS = (
    "class_id",
    "display_label",
    "full_display_label",
    "chrom_arm",
    "chrom_arm_label",
    "community",
    "community_order_position",
    "copy_burden",
)


def read_tsv(path: Path, *, gzipped: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if gzipped else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_supported_terms() -> list[dict[str, str]]:
    rows = read_tsv(V7 / "TERM_RESULTS.tsv.gz", gzipped=True)
    supported = [
        row
        for row in rows
        if row["assignment"] == "midpoint"
        and row["primary_support"] == "1"
        and row["support_status"] == "PRIMARY_SUPPORTED"
    ]
    supported.sort(key=lambda row: int(row["hypothesis_index"]))

    assert len(supported) == 209
    assert len(
        {(row["collection"], row["relation"], row["term_id"]) for row in supported}
    ) == 209
    assert all(row["K_phr_total"] == "187" for row in supported)
    assert all(row["nonphr_total"] == "31779" for row in supported)
    assert all(row["N_eligible_total"] == "31966" for row in supported)
    return supported


def build_arm_matrix() -> list[dict[str, str]]:
    matrix = read_tsv(HEATMAP / "functional_component_arm_matrix.tsv")
    ordering = read_tsv(HEATMAP / "functional_component_arm_orderings.tsv")
    positions = {
        row["chrom_arm"]: row["position"]
        for row in ordering
        if row["ordering"] == "community_order"
    }

    assert len(matrix) == 6 * 48
    assert set(row["class_id"] for row in matrix) == set(CLASS_ORDER)
    assert len(positions) == 48

    class_rank = {class_id: rank for rank, class_id in enumerate(CLASS_ORDER)}
    for row in matrix:
        row["community_order_position"] = positions[row["chrom_arm"]]
    matrix.sort(
        key=lambda row: (
            class_rank[row["class_id"]],
            int(row["community_order_position"]),
        )
    )

    for class_id in CLASS_ORDER:
        class_rows = [row for row in matrix if row["class_id"] == class_id]
        assert len(class_rows) == 48
        assert sorted(int(row["community_order_position"]) for row in class_rows) == list(
            range(1, 49)
        )
    return matrix


def main() -> None:
    write_tsv(
        OUT / "ED_Fig2_supported_ontology_terms.tsv",
        TERM_COLUMNS,
        build_supported_terms(),
    )
    write_tsv(
        OUT / "ED_Fig2_arm_functional_components.tsv",
        ARM_COLUMNS,
        build_arm_matrix(),
    )


if __name__ == "__main__":
    main()
