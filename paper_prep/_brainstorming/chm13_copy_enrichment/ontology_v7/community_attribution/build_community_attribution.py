#!/usr/bin/env python3
"""Build post-inference community attribution for supported V7 terms.

This program does not perform ontology inference. It starts from the frozen V7
term decisions and partitions each supported term's coordinate-copy contributors
across arm-level communities. Functional classes are display-only groupings of
shared functional-source signatures; they neither create p-values nor alter V7
support decisions.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence, Tuple


HERE = Path(__file__).resolve().parent
V7 = HERE.parent
REPO = HERE.parents[4]
TERM_RESULTS = V7 / "TERM_RESULTS.tsv.gz"
CONTRIBUTORS = V7 / "EXACT_TERM_CONTRIBUTORS.tsv.gz"
UPSTREAM_COMMUNITY_SUMMARY = V7 / "COMMUNITY_TERM_SUMMARY.tsv"
UPSTREAM_VALIDATION = V7 / "V7_VALIDATION.json"
ARM_ASSIGNMENTS = REPO / "data" / "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
COMMUNITY_BLOCKS = REPO / "submission" / "fig" / "MainFigures" / "community_blocks.tsv"
CYTOBANDS = REPO / "data" / "chm13v2.0_cytobands_allchrs.bed"

SCHEMA_VERSION = "chm13-physical-copy-ontology-v7-community-attribution-v1.0"
DOMINANCE_THRESHOLD = 0.5
POST_INFERENCE_ROLE = "POST_INFERENCE_DISPLAY_ONLY_NOT_A_HYPOTHESIS_OR_TEST"
REDUNDANCY_WARNING = (
    "Ontology contributor rows are redundant annotations, not independent biological systems; "
    "unique coordinate-copy burden is the nonredundant display unit."
)

TermKey = Tuple[str, str, str]


# Representative terms are prespecified display exemplars from exact primary-
# supported rows. Source membership, not term-name matching, defines each class.
CLASS_SPECS = (
    {
        "class_id": "DUX4_ZGA_TRANSCRIPTION_NUCLEAR_ENVELOPE_CELL_CYCLE",
        "display_label": "DUX4 / ZGA / transcription / nuclear envelope / cell cycle",
        "source_symbols": ("DUX4",),
        "representatives": (
            ("GO_BP", "direct", "GO:0070317"),
            ("Reactome", "direct", "R-HSA-9819196"),
            ("GO_MF", "direct", "GO:0000977"),
            ("GO_CC", "direct", "GO:0031965"),
        ),
        "interpretation": "DUX4-source physical copies carrying exact ZGA, transcription, nuclear-envelope, and cell-cycle terms.",
    },
    {
        "class_id": "WASH_ENDOSOMAL_ACTIN_EXOCYST",
        "display_label": "WASH / endosomal actin / exocyst",
        "source_symbols": ("WASHC1",),
        "representatives": (
            ("GO_CC", "direct", "GO:0071203"),
            ("GO_CC", "direct", "GO:0000145"),
            ("GO_BP", "direct", "GO:0034314"),
            ("GO_BP", "direct", "GO:0016197"),
        ),
        "interpretation": "WASHC1-source copies carrying exact WASH-complex, exocyst, actin-nucleation, and endosomal-transport terms.",
    },
    {
        "class_id": "DDX11_HELICASE_CHROMOSOME",
        "display_label": "DDX11 / helicase / chromosome",
        "source_symbols": ("DDX11",),
        "representatives": (
            ("GO_MF", "direct", "GO:0043139"),
            ("GO_MF", "direct", "GO:0160225"),
            ("GO_BP", "direct", "GO:0034085"),
            ("GO_BP", "direct", "GO:1990700"),
        ),
        "interpretation": "DDX11-source copies carrying exact helicase, G-quadruplex, cohesion, and chromatin terms.",
    },
    {
        "class_id": "RPL23A_RIBOSOMAL_NUCLEOLAR",
        "display_label": "RPL23A / ribosomal / nucleolar",
        "source_symbols": ("RPL23A",),
        "representatives": (
            ("GO_CC", "direct", "GO:0005730"),
            ("GO_MF", "direct", "GO:0019843"),
        ),
        "interpretation": "RPL23A-source copies carrying exact nucleolus and rRNA-binding terms.",
    },
    {
        "class_id": "SEPTIN14_SEPTIN_CYTOKINESIS",
        "display_label": "SEPTIN14 / septin / cytokinesis",
        "source_symbols": ("SEPTIN14", "SEPTIN14P20"),
        "representatives": (
            ("GO_CC", "direct", "GO:0031105"),
            ("GO_CC", "direct", "GO:0005940"),
            ("GO_BP", "direct", "GO:0061640"),
            ("GO_BP", "ancestor", "GO:0000910"),
            ("GO_BP", "direct", "GO:1905719"),
        ),
        "interpretation": "SEPTIN14-source copies carrying exact septin, cytokinesis, and perinuclear-localization terms.",
    },
    {
        "class_id": "WBP1L_CXCL12_SIGNALING",
        "display_label": "WBP1L / CXCL12 signaling",
        "source_symbols": ("WBP1L",),
        "representatives": (
            ("GO_BP", "direct", "GO:0038160"),
            ("GO_BP", "ancestor", "GO:0038146"),
            ("GO_BP", "ancestor", "GO:0038159"),
        ),
        "interpretation": (
            "WBP1L-source copies carrying exact CXCL12/CXCR4 signaling terms. The exact perinuclear-localization "
            "term has a SEPTIN14 contributor signature and is therefore retained in the SEPTIN14 class, not merged here."
        ),
    },
)


EXACT_ATTRIBUTION_FIELDS = (
    "schema_version", "analysis_layer", "assignment", "analysis_role", "collection", "relation",
    "ontology", "namespace", "term_id", "term_name", "support_status", "primary_support",
    "sensitivity_support", "p_exact_upper", "bh_q_within_collection", "by_q_within_collection",
    "holm_p_global", "bonferroni_p_global", "term_phr_burden_a_T", "community", "community_order",
    "community_contributor_copy_burden", "community_fraction_of_a_T", "dominant_community_flag",
    "dominance_threshold", "contributing_arms", "functional_source_ids_json",
    "functional_source_symbols_json", "gene_names_json", "coordinate_copy_count", "copy_ids_json",
    "summary_role",
)

SOURCE_SUMMARY_FIELDS = (
    "schema_version", "analysis_layer", "assignment", "community", "community_order", "community_arms",
    "community_ontology_eligible_phr_copy_burden", "functional_source_id", "functional_source_symbol",
    "supported_unique_coordinate_copy_burden", "fraction_of_community_eligible_burden",
    "supported_term_contributor_row_burden", "supported_exact_term_count", "gene_names_json",
    "copy_ids_json", "functional_class_ids_json", "ontology_redundancy_warning", "summary_role",
)

CLASS_DEFINITION_FIELDS = (
    "schema_version", "class_id", "display_label", "definition_status", "defining_source_symbols_json",
    "primary_supported_exact_term_count_for_signature", "primary_supported_unique_copy_burden",
    "representative_exact_terms_json", "interpretation", "statistical_status", "summary_role",
)

COMMUNITY_CLASS_FIELDS = (
    "schema_version", "analysis_layer", "row_type", "community", "community_order", "community_arms",
    "community_ontology_eligible_phr_copy_burden", "community_supported_term_contributor_row_burden",
    "community_unique_supported_coordinate_copy_burden", "class_id", "display_label", "class_rank_in_community",
    "dominant_functional_class_flag", "class_unique_coordinate_copy_burden",
    "class_fraction_of_community_eligible_burden", "class_fraction_of_all_class_copies",
    "contributing_arms", "functional_source_ids_json", "functional_source_symbols_json", "gene_names_json",
    "copy_ids_json", "representative_exact_terms_json", "unclassified_supported_coordinate_copy_burden",
    "ontology_redundancy_warning", "summary_role",
)

CLASS_COMMUNITY_FIELDS = (
    "schema_version", "analysis_layer", "class_id", "display_label", "total_unique_coordinate_copy_burden",
    "community_count", "communities_json", "physical_copy_contributors_per_community_json",
    "top_community", "top_community_copy_burden", "top_community_fraction", "distribution_status",
    "representative_exact_terms_json", "interpretation", "summary_role",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: stringify(row.get(field, "")) for field in fields})


def stringify(value: object) -> str:
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def json_list(values: Iterable[str]) -> str:
    return json.dumps(sorted(set(values)), ensure_ascii=False, separators=(",", ":"))


def term_key(row: Mapping[str, str]) -> TermKey:
    return row["collection"], row["relation"], row["term_id"]


def representative_json(keys: Iterable[TermKey], terms: Mapping[TermKey, Mapping[str, str]]) -> str:
    records = []
    for collection, relation, term_id in keys:
        row = terms[(collection, relation, term_id)]
        records.append({
            "collection": collection,
            "relation": relation,
            "term_id": term_id,
            "term_name": row["term_name"],
        })
    return json.dumps(records, ensure_ascii=False, separators=(",", ":"))


def load_community_metadata() -> tuple[dict[str, str], dict[str, int], dict[str, list[str]]]:
    arm_to_community: dict[str, str] = {}
    for row in read_tsv(ARM_ASSIGNMENTS):
        arm = row["ChromArm"].replace("_parm", "_p").replace("_qarm", "_q")
        if arm in arm_to_community and arm_to_community[arm] != row["Community"]:
            raise ValueError(f"conflicting community assignments for {arm}")
        arm_to_community[arm] = row["Community"]

    block_rows = read_tsv(COMMUNITY_BLOCKS)
    order = {row["Community"]: int(row["start_position"]) for row in block_rows}
    arms: dict[str, list[str]] = defaultdict(list)
    for arm, community in arm_to_community.items():
        arms[community].append(arm)
    return arm_to_community, order, {key: sorted(value) for key, value in arms.items()}


def q_arm_starts() -> dict[str, int]:
    starts: dict[str, list[int]] = defaultdict(list)
    with CYTOBANDS.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chromosome, start, _end, band = line.rstrip("\n").split("\t")[:4]
            if band.startswith("q"):
                starts[chromosome].append(int(start))
    return {chromosome: min(values) for chromosome, values in starts.items()}


def annotate_copy(row: Mapping[str, str], arm_to_community: Mapping[str, str],
                  q_starts: Mapping[str, int]) -> dict[str, str]:
    result = dict(row)
    chromosome = row["seqid"]
    midpoint = (int(row["start0"]) + int(row["end0"])) // 2
    arm = f"{chromosome}_{'p' if midpoint < q_starts[chromosome] else 'q'}"
    result["_arm"] = arm
    result["_community"] = arm_to_community.get(arm, "UNASSIGNED_NO_SIGNAL_COMMUNITY")
    return result


def unique_copy_records(rows: Iterable[Mapping[str, str]]) -> dict[str, Mapping[str, str]]:
    records: dict[str, Mapping[str, str]] = {}
    stable = ("seqid", "start0", "end0", "strand", "gene_name", "functional_source_id", "functional_source_symbol",
              "physical_copy_cn", "_arm", "_community")
    for row in rows:
        copy_id = row["copy_id"]
        if copy_id in records and any(records[copy_id][field] != row[field] for field in stable):
            raise ValueError(f"inconsistent metadata for coordinate copy {copy_id}")
        records[copy_id] = row
    return records


def supported_terms(term_rows: Sequence[Mapping[str, str]]) -> tuple[dict[TermKey, Mapping[str, str]], dict[TermKey, Mapping[str, str]]]:
    primary = {term_key(row): row for row in term_rows
               if row["assignment"] == "midpoint" and (row["support_status"] == "PRIMARY_SUPPORTED" or row["primary_support"] == "1")}
    sensitivity = {term_key(row): row for row in term_rows
                   if row["assignment"] == "overlap" and (row["support_status"] == "SENSITIVITY_SUPPORTED" or row["sensitivity_support"] == "1")}
    return primary, sensitivity


def build_exact_attribution(terms: Mapping[TermKey, Mapping[str, str]], contributors: Sequence[Mapping[str, str]],
                            assignment: str, analysis_layer: str, community_order: Mapping[str, int]) -> list[dict[str, object]]:
    partition_field = "v7_midpoint_partition" if assignment == "midpoint" else "v7_overlap_partition"
    grouped: dict[tuple[TermKey, str], list[Mapping[str, str]]] = defaultdict(list)
    seen_edges: set[tuple[TermKey, str]] = set()
    for row in contributors:
        key = term_key(row)
        if key not in terms or row[partition_field] != "PHR":
            continue
        edge = (key, row["copy_id"])
        if edge in seen_edges:
            raise ValueError(f"duplicate exact term-copy contributor edge: {edge}")
        seen_edges.add(edge)
        grouped[(key, row["_community"])].append(row)

    output = []
    for (key, community), rows in grouped.items():
        term = terms[key]
        burden = sum(int(row["physical_copy_cn"]) for row in rows)
        a_t = int(term["a_T"])
        fraction = burden / a_t
        output.append({
            "schema_version": SCHEMA_VERSION,
            "analysis_layer": analysis_layer,
            "assignment": assignment,
            "analysis_role": term["analysis_role"],
            "collection": term["collection"],
            "relation": term["relation"],
            "ontology": term["ontology"],
            "namespace": term["namespace"],
            "term_id": term["term_id"],
            "term_name": term["term_name"],
            "support_status": term["support_status"],
            "primary_support": term["primary_support"],
            "sensitivity_support": term["sensitivity_support"],
            "p_exact_upper": term["p_exact_upper"],
            "bh_q_within_collection": term["bh_q_within_collection"],
            "by_q_within_collection": term["by_q_within_collection"],
            "holm_p_global": term["holm_p_global"],
            "bonferroni_p_global": term["bonferroni_p_global"],
            "term_phr_burden_a_T": a_t,
            "community": community,
            "community_order": community_order.get(community, 9999),
            "community_contributor_copy_burden": burden,
            "community_fraction_of_a_T": fraction,
            "dominant_community_flag": int(fraction >= DOMINANCE_THRESHOLD),
            "dominance_threshold": DOMINANCE_THRESHOLD,
            "contributing_arms": "|".join(sorted({row["_arm"] for row in rows})),
            "functional_source_ids_json": json_list(row["functional_source_id"] for row in rows),
            "functional_source_symbols_json": json_list(row["functional_source_symbol"] for row in rows),
            "gene_names_json": json_list(row["gene_name"] for row in rows),
            "coordinate_copy_count": len({row["copy_id"] for row in rows}),
            "copy_ids_json": json_list(row["copy_id"] for row in rows),
            "summary_role": POST_INFERENCE_ROLE,
        })
    return sorted(output, key=lambda row: (row["analysis_layer"], row["collection"], row["relation"],
                                           row["term_id"], row["community_order"], row["community"]))


def source_class_map() -> dict[str, list[str]]:
    result: dict[str, list[str]] = defaultdict(list)
    for spec in CLASS_SPECS:
        for source in spec["source_symbols"]:
            result[source].append(spec["class_id"])
    return result


def build_source_summary(primary_terms: Mapping[TermKey, Mapping[str, str]], contributors: Sequence[Mapping[str, str]],
                         community_order: Mapping[str, int], community_arms: Mapping[str, list[str]]) -> list[dict[str, object]]:
    primary_rows = [row for row in contributors if row["v7_midpoint_partition"] == "PHR" and term_key(row) in primary_terms]
    all_midpoint_copies = unique_copy_records(row for row in contributors if row["v7_midpoint_partition"] == "PHR")
    total_by_community = Counter(row["_community"] for row in all_midpoint_copies.values())
    grouped: dict[tuple[str, str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in primary_rows:
        grouped[(row["_community"], row["functional_source_id"], row["functional_source_symbol"])].append(row)
    classes = source_class_map()
    output = []
    for (community, source_id, source_symbol), rows in grouped.items():
        copies = unique_copy_records(rows)
        unique_burden = sum(int(row["physical_copy_cn"]) for row in copies.values())
        output.append({
            "schema_version": SCHEMA_VERSION,
            "analysis_layer": "HEADLINE_PRIMARY_MIDPOINT",
            "assignment": "midpoint",
            "community": community,
            "community_order": community_order.get(community, 9999),
            "community_arms": "|".join(community_arms.get(community, [])),
            "community_ontology_eligible_phr_copy_burden": total_by_community[community],
            "functional_source_id": source_id,
            "functional_source_symbol": source_symbol,
            "supported_unique_coordinate_copy_burden": unique_burden,
            "fraction_of_community_eligible_burden": unique_burden / total_by_community[community],
            "supported_term_contributor_row_burden": sum(int(row["physical_copy_cn"]) for row in rows),
            "supported_exact_term_count": len({term_key(row) for row in rows}),
            "gene_names_json": json_list(row["gene_name"] for row in copies.values()),
            "copy_ids_json": json_list(copies),
            "functional_class_ids_json": json_list(classes.get(source_symbol, [])),
            "ontology_redundancy_warning": REDUNDANCY_WARNING,
            "summary_role": POST_INFERENCE_ROLE,
        })
    return sorted(output, key=lambda row: (row["community_order"], row["community"],
                                           -row["supported_unique_coordinate_copy_burden"], row["functional_source_symbol"],
                                           row["functional_source_id"]))


def build_class_outputs(primary_terms: Mapping[TermKey, Mapping[str, str]], contributors: Sequence[Mapping[str, str]],
                        community_order: Mapping[str, int], community_arms: Mapping[str, list[str]]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    primary_rows = [row for row in contributors if row["v7_midpoint_partition"] == "PHR" and term_key(row) in primary_terms]
    all_midpoint = unique_copy_records(row for row in contributors if row["v7_midpoint_partition"] == "PHR")
    supported_copies = unique_copy_records(primary_rows)
    eligible_by_community = Counter(row["_community"] for row in all_midpoint.values())
    supported_by_community = Counter(row["_community"] for row in supported_copies.values())
    redundant_by_community = Counter()
    for row in primary_rows:
        redundant_by_community[row["_community"]] += int(row["physical_copy_cn"])

    definitions = []
    class_data: dict[str, dict[str, object]] = {}
    for spec in CLASS_SPECS:
        sources = set(spec["source_symbols"])
        rows = [row for row in primary_rows if row["functional_source_symbol"] in sources]
        copies = unique_copy_records(rows)
        reps = tuple(spec["representatives"])
        missing = [key for key in reps if key not in primary_terms]
        if missing:
            raise ValueError(f"display class {spec['class_id']} has non-supported representative terms: {missing}")
        for key in reps:
            if not any(term_key(row) == key for row in rows):
                raise ValueError(f"representative {key} lacks the class contributor signature {spec['source_symbols']}")
        all_signature_terms = {term_key(row) for row in rows}
        definition_status = "PRIMARY_SUPPORTED" if copies else "EVALUATED_NOT_PRIMARY_SUPPORTED"
        rep_json = representative_json(reps, primary_terms)
        definitions.append({
            "schema_version": SCHEMA_VERSION,
            "class_id": spec["class_id"],
            "display_label": spec["display_label"],
            "definition_status": definition_status,
            "defining_source_symbols_json": json_list(sources),
            "primary_supported_exact_term_count_for_signature": len(all_signature_terms),
            "primary_supported_unique_copy_burden": sum(int(row["physical_copy_cn"]) for row in copies.values()),
            "representative_exact_terms_json": rep_json,
            "interpretation": spec["interpretation"],
            "statistical_status": "DISPLAY_ONLY_NO_NEW_P_VALUE_V7_DECISIONS_UNCHANGED",
            "summary_role": POST_INFERENCE_ROLE,
        })
        class_data[spec["class_id"]] = {"spec": spec, "rows": rows, "copies": copies, "rep_json": rep_json}

    class_rows = []
    class_copy_ids_by_community: dict[str, set[str]] = defaultdict(set)
    for class_id, data in class_data.items():
        spec = data["spec"]
        copies = data["copies"]
        grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
        for row in copies.values():
            grouped[row["_community"]].append(row)
        class_total = len(copies)
        provisional = []
        for community, rows in grouped.items():
            burden = sum(int(row["physical_copy_cn"]) for row in rows)
            class_copy_ids_by_community[community].update(row["copy_id"] for row in rows)
            provisional.append({
                "schema_version": SCHEMA_VERSION,
                "analysis_layer": "HEADLINE_PRIMARY_MIDPOINT",
                "row_type": "FUNCTIONAL_CLASS",
                "community": community,
                "community_order": community_order.get(community, 9999),
                "community_arms": "|".join(community_arms.get(community, [])),
                "community_ontology_eligible_phr_copy_burden": eligible_by_community[community],
                "community_supported_term_contributor_row_burden": redundant_by_community[community],
                "community_unique_supported_coordinate_copy_burden": supported_by_community[community],
                "class_id": class_id,
                "display_label": spec["display_label"],
                "class_unique_coordinate_copy_burden": burden,
                "class_fraction_of_community_eligible_burden": burden / eligible_by_community[community],
                "class_fraction_of_all_class_copies": burden / class_total,
                "contributing_arms": "|".join(sorted({row["_arm"] for row in rows})),
                "functional_source_ids_json": json_list(row["functional_source_id"] for row in rows),
                "functional_source_symbols_json": json_list(row["functional_source_symbol"] for row in rows),
                "gene_names_json": json_list(row["gene_name"] for row in rows),
                "copy_ids_json": json_list(row["copy_id"] for row in rows),
                "representative_exact_terms_json": data["rep_json"],
                "ontology_redundancy_warning": REDUNDANCY_WARNING,
                "summary_role": POST_INFERENCE_ROLE,
            })
        class_rows.extend(provisional)

    # Rank the display classes mechanically within each community. Tied maxima
    # are all marked dominant; these labels have no inferential meaning.
    rows_by_community: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in class_rows:
        rows_by_community[str(row["community"])].append(row)
    final_class_rows = []
    # Include the full k=15 assignment catalog. A community with zero eligible
    # midpoint copies is still an informative zero, not a missing row.
    all_communities = sorted(community_arms, key=lambda community: (community_order.get(community, 9999), community))
    for community in all_communities:
        rows = rows_by_community.get(community, [])
        ranked_burdens = sorted({int(row["class_unique_coordinate_copy_burden"]) for row in rows}, reverse=True)
        for row in rows:
            burden = int(row["class_unique_coordinate_copy_burden"])
            row["class_rank_in_community"] = ranked_burdens.index(burden) + 1
            row["dominant_functional_class_flag"] = int(burden == ranked_burdens[0])
            row["unclassified_supported_coordinate_copy_burden"] = (
                supported_by_community[community] - len(class_copy_ids_by_community[community])
            )
            final_class_rows.append(row)
        if not rows:
            final_class_rows.append({
                "schema_version": SCHEMA_VERSION,
                "analysis_layer": "HEADLINE_PRIMARY_MIDPOINT",
                "row_type": "COMMUNITY_TOTAL_NO_DEFINED_CLASS",
                "community": community,
                "community_order": community_order.get(community, 9999),
                "community_arms": "|".join(community_arms.get(community, [])),
                "community_ontology_eligible_phr_copy_burden": eligible_by_community[community],
                "community_supported_term_contributor_row_burden": redundant_by_community[community],
                "community_unique_supported_coordinate_copy_burden": supported_by_community[community],
                "unclassified_supported_coordinate_copy_burden": supported_by_community[community],
                "ontology_redundancy_warning": REDUNDANCY_WARNING,
                "summary_role": POST_INFERENCE_ROLE,
            })

    summaries = []
    for class_id, data in class_data.items():
        copies = data["copies"]
        counts = Counter(row["_community"] for row in copies.values())
        total = sum(counts.values())
        top_community, top_burden = sorted(counts.items(), key=lambda item: (-item[1], community_order.get(item[0], 9999), item[0]))[0]
        top_fraction = top_burden / total
        ordered_counts = {community: counts[community] for community in sorted(counts, key=lambda value: (community_order.get(value, 9999), value))}
        summaries.append({
            "schema_version": SCHEMA_VERSION,
            "analysis_layer": "HEADLINE_PRIMARY_MIDPOINT",
            "class_id": class_id,
            "display_label": data["spec"]["display_label"],
            "total_unique_coordinate_copy_burden": total,
            "community_count": len(counts),
            "communities_json": json_list(counts),
            "physical_copy_contributors_per_community_json": json.dumps(ordered_counts, separators=(",", ":")),
            "top_community": top_community,
            "top_community_copy_burden": top_burden,
            "top_community_fraction": top_fraction,
            "distribution_status": "ONE_COMMUNITY_DOMINATES" if top_fraction >= DOMINANCE_THRESHOLD else "DISTRIBUTED_NO_COMMUNITY_AT_50_PERCENT",
            "representative_exact_terms_json": data["rep_json"],
            "interpretation": data["spec"]["interpretation"],
            "summary_role": POST_INFERENCE_ROLE,
        })

    return definitions, sorted(final_class_rows, key=lambda row: (row["community_order"], row["community"],
                                                                    row.get("class_rank_in_community", 9999), row.get("class_id", ""))), summaries


def format_term_names(rep_json: str, limit: int = 3) -> str:
    terms = json.loads(rep_json)
    return "; ".join(f"{row['term_name']} ({row['term_id']}, {row['relation']})" for row in terms[:limit])


def build_report(primary_terms: Mapping[TermKey, Mapping[str, str]], source_rows: Sequence[Mapping[str, object]],
                 class_rows: Sequence[Mapping[str, object]], class_summaries: Sequence[Mapping[str, object]],
                 community_order: Mapping[str, int]) -> str:
    source_by_community: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in source_rows:
        source_by_community[str(row["community"])].append(row)
    classes_by_community: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in class_rows:
        if row["row_type"] == "FUNCTIONAL_CLASS":
            classes_by_community[str(row["community"])].append(row)

    totals_by_community = {str(row["community"]): row for row in class_rows}
    communities = sorted(totals_by_community, key=lambda value: (community_order.get(value, 9999), value))
    lines = [
        "# V7 community functional attribution", "",
        "## Interpretation", "",
        "This is a **post-inference attribution** of the 209 exact midpoint `PRIMARY_SUPPORTED` V7 term rows. "
        "It introduces no ontology test, no functional-class p-value, and no change to any V7 support decision. "
        "The inferential object remains each exact `(collection, relation, term_id)` row.", "",
        "Physical-copy burden is counted at coordinate resolution: N distinct coordinate copies contribute N units, "
        "even when they share a gene name, source gene, family, locus family, arm community, or display class. "
        "Ontology rows are redundant annotations and must not be interpreted as independent biological systems. "
        "The class labels below are post-inference display summaries only.", "",
        "The primary universe contains 187 ontology-eligible midpoint PHR copies; 186 distinct copies contribute to "
        "at least one supported exact term. Any-overlap-supported rows are retained only as a labeled sensitivity layer "
        "in `COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv` and are not used to define the headline classes.", "",
        "## Clean community–function mappings", "",
        "| Display class | Copy burden by community | Pattern | Representative exact supported terms |", 
        "|---|---|---|---|",
    ]
    for row in class_summaries:
        counts = json.loads(str(row["physical_copy_contributors_per_community_json"]))
        mapping = ", ".join(f"{community}: {count}" for community, count in counts.items())
        pattern = (f"{row['top_community']} carries {row['top_community_copy_burden']}/{row['total_unique_coordinate_copy_burden']}"
                   if row["distribution_status"] == "ONE_COMMUNITY_DOMINATES" else
                   f"distributed across {row['community_count']} communities")
        lines.append(f"| {row['display_label']} | {mapping} | {pattern} | {format_term_names(str(row['representative_exact_terms_json']))} |")

    lines.extend([
        "", "C1 is the cleanest mapping: all DUX4-source contributors fall in C1, tying that community to exact ZGA, "
        "transcription-regulatory, nuclear-envelope, and cell-cycle signals. WASHC1 and DDX11 signatures are "
        "distributed across multiple communities rather than constituting a single-community system. RPL23A and "
        "SEPTIN14 likewise span communities, with their exact nucleolar/rRNA-binding and septin/cytokinesis terms "
        "retained as separate source-defined displays.", "",
        "WBP1L-source copies robustly carry exact CXCL12/CXCR4 signaling rows. The supported exact term `protein "
        "localization to perinuclear region of cytoplasm` has a SEPTIN14 contributor signature, so it remains in the "
        "SEPTIN14 class rather than being conflated with WBP1L/CXCL12 signaling.", "",
        "## Community synopsis", "",
        "The repeated contributor-row burden below is supplied only as an audit of annotation incidence; because one "
        "copy can contribute to many related direct and ancestor terms, it is ontology-redundant. Unique-copy burden "
        "is the interpretable physical-copy summary.", "",
        "| Community | Eligible PHR copies | Unique supported copies | Redundant supported term-copy rows | Dominant display class(es) | Contributing source symbols |",
        "|---|---:|---:|---:|---|---|",
    ])
    for community in communities:
        sources = source_by_community.get(community, [])
        total_row = totals_by_community[community]
        eligible = int(total_row["community_ontology_eligible_phr_copy_burden"])
        copy_ids = set()
        redundant = int(total_row["community_supported_term_contributor_row_burden"])
        symbols = set()
        for row in sources:
            copy_ids.update(json.loads(str(row["copy_ids_json"])))
            symbols.add(str(row["functional_source_symbol"]))
        dominant = [str(row["display_label"]) for row in classes_by_community.get(community, [])
                    if int(row["dominant_functional_class_flag"]) == 1]
        lines.append(f"| {community} | {eligible} | {len(copy_ids)} | {redundant} | {'; '.join(dominant) or 'no defined display class'} | {', '.join(sorted(symbols)) or 'none'} |")

    immune_or_olfactory = [row for row in primary_terms.values()
                           if any(token in row["term_name"].lower() for token in ("immune", "immun", "olfact", "odorant"))]
    lines.extend([
        "", "## Boundary and negative findings", "",
        ("No broad immune or olfactory term is a headline-supported exact midpoint enrichment in V7. "
         if not immune_or_olfactory else
         "Broad immune/olfactory labels require term-by-term review because at least one matching exact row is supported. "),
        "Source and class labels describe annotation-bearing CHM13 physical copies, not expression, protein activity, "
        "dosage, retained pseudogene function, biological independence of neighboring copies, or population prevalence. "
        "The any-overlap layer is a sensitivity analysis and should not replace the midpoint headline.", "",
        "## Reproducibility", "",
        "Run `python3 build_community_attribution.py` from any directory, then "
        "`python3 -m unittest -v test_community_attribution.py`. JSON arrays are used for copy IDs because the atomic "
        "V7 copy identifiers themselves contain pipe characters. `OUTPUT_MANIFEST.sha256.tsv` records the exact files "
        "covered by release hashes.", "",
    ])
    return "\n".join(lines)


def validate_outputs(output: Path, upstream_validation: Mapping[str, object],
                     primary_terms: Mapping[TermKey, Mapping[str, str]], sensitivity_terms: Mapping[TermKey, Mapping[str, str]],
                     contributors: Sequence[Mapping[str, str]], exact_rows: Sequence[Mapping[str, object]],
                     source_rows: Sequence[Mapping[str, object]], definitions: Sequence[Mapping[str, object]],
                     class_rows: Sequence[Mapping[str, object]], upstream_community_rows: Sequence[Mapping[str, str]],
                     expected_communities: set[str], report: str) -> dict[str, object]:
    checks = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    check("v7_validation_input_status", upstream_validation.get("status") == "PASS",
          {"status": upstream_validation.get("status"), "checks": upstream_validation.get("checks_total")})

    burdens = Counter()
    for row in exact_rows:
        burdens[(row["assignment"], row["collection"], row["relation"], row["term_id"])] += int(row["community_contributor_copy_burden"])
    primary_ok = all(burdens[("midpoint",) + key] == int(term["a_T"]) for key, term in primary_terms.items())
    sensitivity_ok = all(burdens[("overlap",) + key] == int(term["a_T"]) for key, term in sensitivity_terms.items())
    check("every_primary_term_community_burden_sums_to_a_T", primary_ok,
          {"primary_supported_terms": len(primary_terms)})
    check("every_sensitivity_term_community_burden_sums_to_a_T", sensitivity_ok,
          {"overlap_supported_terms": len(sensitivity_terms)})
    regenerated = {
        (str(row["assignment"]), str(row["collection"]), str(row["relation"]), str(row["term_id"]), str(row["community"])):
        int(row["community_contributor_copy_burden"]) for row in exact_rows
    }
    upstream = {
        (row["assignment"], row["collection"], row["relation"], row["term_id"], row["community"]):
        int(row["community_contributor_copy_cn"]) for row in upstream_community_rows
    }
    check("exact_attribution_reconciles_to_upstream_v7_community_summary", regenerated == upstream,
          {"regenerated_rows": len(regenerated), "upstream_rows": len(upstream)})
    check("headline_contains_only_primary_supported_midpoint_rows", all(
        row["analysis_layer"] != "HEADLINE_PRIMARY_MIDPOINT" or
        (row["assignment"] == "midpoint" and row["support_status"] == "PRIMARY_SUPPORTED" and row["primary_support"] == "1")
        for row in exact_rows), {"headline_rows": sum(row["analysis_layer"] == "HEADLINE_PRIMARY_MIDPOINT" for row in exact_rows)})
    check("overlap_rows_explicitly_labeled_sensitivity", all(
        row["assignment"] != "overlap" or
        (row["analysis_layer"] == "SENSITIVITY_ANY_OVERLAP" and row["support_status"] == "SENSITIVITY_SUPPORTED")
        for row in exact_rows), {"sensitivity_rows": sum(row["assignment"] == "overlap" for row in exact_rows)})
    check("dominant_term_community_flags_mechanical", all(
        int(row["dominant_community_flag"]) == int(float(row["community_fraction_of_a_T"]) >= DOMINANCE_THRESHOLD)
        for row in exact_rows), {"threshold": DOMINANCE_THRESHOLD})
    check("exact_attribution_retains_atomic_coordinate_copies", all(
        int(row["coordinate_copy_count"]) == len(json.loads(str(row["copy_ids_json"]))) == int(row["community_contributor_copy_burden"])
        for row in exact_rows), "coordinate_copy_count == JSON copy list length == physical-copy burden for every row")

    primary_contributors = [row for row in contributors if row["v7_midpoint_partition"] == "PHR" and term_key(row) in primary_terms]
    expected_source: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    expected_redundant = Counter()
    for row in primary_contributors:
        key = row["_community"], row["functional_source_id"], row["functional_source_symbol"]
        expected_source[key].add(row["copy_id"])
        expected_redundant[key] += int(row["physical_copy_cn"])
    source_ok = len(source_rows) == len(expected_source) and all(
        set(json.loads(str(row["copy_ids_json"]))) == expected_source[(row["community"], row["functional_source_id"], row["functional_source_symbol"])] and
        int(row["supported_term_contributor_row_burden"]) == expected_redundant[(row["community"], row["functional_source_id"], row["functional_source_symbol"])]
        for row in source_rows)
    check("community_source_summaries_reproducible_from_exact_contributors", source_ok,
          {"community_source_rows": len(source_rows)})
    prohibited = {"p_value", "p_exact_upper", "bh_q_within_collection", "by_q_within_collection", "holm_p_global", "bonferroni_p_global"}
    class_fields = set(CLASS_DEFINITION_FIELDS) | set(COMMUNITY_CLASS_FIELDS) | set(CLASS_COMMUNITY_FIELDS)
    check("functional_classes_create_no_p_values_or_support_decisions", not (prohibited & class_fields) and all(
        row["statistical_status"] == "DISPLAY_ONLY_NO_NEW_P_VALUE_V7_DECISIONS_UNCHANGED" for row in definitions),
        {"prohibited_fields_present": sorted(prohibited & class_fields), "classes": len(definitions)})
    check("all_required_functional_classes_grounded_in_primary_terms_and_copies", len(definitions) == len(CLASS_SPECS) and all(
        row["definition_status"] == "PRIMARY_SUPPORTED" and int(row["primary_supported_unique_copy_burden"]) > 0 and
        len(json.loads(str(row["representative_exact_terms_json"]))) > 0 for row in definitions),
        {"classes": [row["class_id"] for row in definitions]})
    covered_communities = {row["community"] for row in class_rows}
    check("every_k15_community_has_a_summary_row_including_zero_burdens", covered_communities == expected_communities,
          {"communities": sorted(covered_communities)})
    required_statements = ("post-inference attribution", "Ontology rows are redundant", "post-inference display summaries only")
    check("report_states_inference_and_redundancy_boundaries", all(statement in report for statement in required_statements),
          {"required_statements": list(required_statements)})
    check("broad_immune_and_olfactory_not_headline_supported", not any(
        any(token in row["term_name"].lower() for token in ("immune", "immun", "olfact", "odorant"))
        for row in primary_terms.values()), "no matching exact primary-supported term names")

    errors = [row["check"] for row in checks if row["status"] != "PASS"]
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "checks_passed": len(checks) - len(errors),
        "checks_total": len(checks),
        "checks": checks,
        "primary_supported_exact_terms": len(primary_terms),
        "overlap_sensitivity_supported_exact_terms": len(sensitivity_terms),
        "dominant_community_threshold": DOMINANCE_THRESHOLD,
        "statistical_boundary": "POST_INFERENCE_ATTRIBUTION_ONLY_NO_NEW_TESTS",
        "ontology_redundancy_warning": REDUNDANCY_WARNING,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(output: Path) -> None:
    names = (
        "COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv",
        "COMMUNITY_SOURCE_COPY_SUMMARY.tsv",
        "FUNCTIONAL_CLASS_DEFINITIONS.tsv",
        "COMMUNITY_FUNCTIONAL_CLASSES.tsv",
        "FUNCTIONAL_CLASS_COMMUNITY_SUMMARY.tsv",
        "COMMUNITY_FUNCTIONAL_ATTRIBUTION_REPORT.md",
        "COMMUNITY_FUNCTIONAL_ATTRIBUTION_VALIDATION.json",
        "build_community_attribution.py",
        "test_community_attribution.py",
    )
    paths = {name: (HERE / name if name in {"build_community_attribution.py", "test_community_attribution.py"}
                    else output / name) for name in names}
    write_tsv(output / "OUTPUT_MANIFEST.sha256.tsv", ("sha256", "bytes", "path"), (
        {"sha256": sha256(paths[name]), "bytes": paths[name].stat().st_size, "path": name} for name in names
    ))


def generate(output: Path = HERE) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    upstream_validation = json.loads(UPSTREAM_VALIDATION.read_text(encoding="utf-8"))
    upstream_community_rows = read_tsv(UPSTREAM_COMMUNITY_SUMMARY)
    term_rows = read_tsv(TERM_RESULTS)
    raw_contributors = read_tsv(CONTRIBUTORS)
    arm_to_community, community_order, community_arms = load_community_metadata()
    q_starts = q_arm_starts()
    contributors = [annotate_copy(row, arm_to_community, q_starts) for row in raw_contributors]
    primary_terms, sensitivity_terms = supported_terms(term_rows)

    exact_rows = (
        build_exact_attribution(primary_terms, contributors, "midpoint", "HEADLINE_PRIMARY_MIDPOINT", community_order) +
        build_exact_attribution(sensitivity_terms, contributors, "overlap", "SENSITIVITY_ANY_OVERLAP", community_order)
    )
    exact_rows.sort(key=lambda row: (0 if row["analysis_layer"] == "HEADLINE_PRIMARY_MIDPOINT" else 1,
                                     row["collection"], row["relation"], row["term_id"], row["community_order"], row["community"]))
    source_rows = build_source_summary(primary_terms, contributors, community_order, community_arms)
    definitions, class_rows, class_summaries = build_class_outputs(primary_terms, contributors, community_order, community_arms)
    report = build_report(primary_terms, source_rows, class_rows, class_summaries, community_order)

    write_tsv(output / "COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv", EXACT_ATTRIBUTION_FIELDS, exact_rows)
    write_tsv(output / "COMMUNITY_SOURCE_COPY_SUMMARY.tsv", SOURCE_SUMMARY_FIELDS, source_rows)
    write_tsv(output / "FUNCTIONAL_CLASS_DEFINITIONS.tsv", CLASS_DEFINITION_FIELDS, definitions)
    write_tsv(output / "COMMUNITY_FUNCTIONAL_CLASSES.tsv", COMMUNITY_CLASS_FIELDS, class_rows)
    write_tsv(output / "FUNCTIONAL_CLASS_COMMUNITY_SUMMARY.tsv", CLASS_COMMUNITY_FIELDS, class_summaries)
    (output / "COMMUNITY_FUNCTIONAL_ATTRIBUTION_REPORT.md").write_text(report, encoding="utf-8")
    validation = validate_outputs(output, upstream_validation, primary_terms, sensitivity_terms, contributors,
                                  exact_rows, source_rows, definitions, class_rows, upstream_community_rows,
                                  set(community_arms), report)
    (output / "COMMUNITY_FUNCTIONAL_ATTRIBUTION_VALIDATION.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_manifest(output)
    if validation["status"] != "PASS":
        raise RuntimeError(f"community attribution validation failed: {validation['errors']}")
    return validation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE, help="output directory (default: script directory)")
    args = parser.parse_args()
    validation = generate(args.output.resolve())
    print(f"Community attribution validation: {validation['status']} ({validation['checks_passed']}/{validation['checks_total']})")


if __name__ == "__main__":
    main()
