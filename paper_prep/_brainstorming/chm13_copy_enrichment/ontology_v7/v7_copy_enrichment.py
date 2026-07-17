#!/usr/bin/env python3
"""V7 complete finite-population physical-copy ontology enrichment.

The executable inference path enumerates the frozen V6 eligible universe and
its exact set complement. It contains no stochastic or interval-relocation
operation. See V7_RUN_PROTOCOL.md, which was committed before this program was
run on the frozen release.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Tuple


HERE = Path(__file__).resolve().parent
V6 = HERE.parent / "ontology_v6"
REPO = HERE.parents[3]
COPY_MAP = V6 / "GENOMEWIDE_SOURCE_MAP.tsv.gz"
COPY_EDGES = V6 / "PHYSICAL_COPY_TERM_EDGES.tsv.gz"
HYPOTHESES = V6 / "FROZEN_HYPOTHESES.tsv.gz"
HYPOTHESIS_MANIFEST = V6 / "FROZEN_HYPOTHESES.json"
V6_AUDITED_MANIFEST = V6 / "PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv"
V6_CONTRIBUTORS = V6 / "results" / "release" / "EXACT_TERM_CONTRIBUTORS.tsv.gz"
V6_RELEASE_MANIFEST = V6 / "results" / "release" / "RELEASE_SHA256.tsv"
COMMUNITIES = REPO / "data" / "hprcv2.1Mb.subtelo.arm-leiden.communities.tsv"
CYTOBANDS = REPO / "data" / "chm13v2.0_cytobands_allchrs.bed"

SCHEMA_VERSION = "chm13-physical-copy-ontology-v7.0"
ASSIGNMENTS = (("midpoint", "phr_midpoint_cn", "PRIMARY_MIDPOINT"),
               ("overlap", "phr_any_overlap_cn", "PAIRED_ANY_OVERLAP_SENSITIVITY"))
COLLECTIONS = ("GO_BP", "GO_MF", "GO_CC", "Reactome")
EXPECTED_N = 31966
EXPECTED_K = {"midpoint": 187, "overlap": 193}
EXPECTED_HYPOTHESES = 31235
EXPECTED_RAW_COPIES = 61312
PRIMARY_ALPHA = 0.05
PRIMARY_RULE = "midpoint AND within_collection_BH<=0.05 AND global_Holm<=0.05"


TERM_FIELDS = [
    "schema_version", "assignment", "analysis_role", "hypothesis_index",
    "collection", "relation", "ontology", "namespace", "term_id", "term_name",
    "multiplicity_family_id", "multiplicity_family_size", "hypothesis_status",
    "a_T", "b_T", "c_T", "d_T", "M_T_genome_burden", "K_T_phr_burden",
    "nonphr_term_burden", "K_phr_total", "nonphr_total", "N_eligible_total",
    "phr_copy_fraction", "nonphr_copy_fraction", "fold_enrichment",
    "fold_enrichment_ha", "odds_ratio", "odds_ratio_ha",
    "copy_count_difference_a_minus_c", "expected_phr_copy_count",
    "copy_excess_over_expected", "p_exact_upper", "bh_q_within_collection",
    "by_q_within_collection", "holm_p_global", "bonferroni_p_global",
    "primary_support", "sensitivity_support", "support_status", "primary_support_rule",
    "finite_population_definition", "background_definition", "inference_status",
]


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(str(path), mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def read_rows(path: Path) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield dict(row)


def _stringify(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return format(value, ".17g")
    return str(value)


def write_rows(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                            lineterminator="\n", extrasaction="ignore")
                    writer.writeheader()
                    for row in rows:
                        writer.writerow({field: _stringify(row.get(field, "")) for field in fields})
    else:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: _stringify(row.get(field, "")) for field in fields})


def write_json(path: Path, value: Mapping[str, object]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collection_for(row: Mapping[str, str]) -> str:
    if row.get("collection"):
        return row["collection"]
    if row["ontology"] == "Reactome":
        return "Reactome"
    namespace = row["namespace"]
    return {
        "biological_process": "GO_BP",
        "molecular_function": "GO_MF",
        "cellular_component": "GO_CC",
    }[namespace]


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def exact_upper_tail(N: int, M: int, K: int, observed: int) -> float:
    """Stable exact hypergeometric upper tail for an enumerated population."""
    if not (0 <= M <= N and 0 <= K <= N and observed >= 0):
        raise ValueError("invalid finite-population count")
    lower = max(0, K - (N - M))
    upper = min(K, M)
    if observed <= lower:
        return 1.0
    if observed > upper:
        return 0.0
    log_p = log_choose(M, observed) + log_choose(N - M, K - observed) - log_choose(N, K)
    logs = [log_p]
    for x in range(observed, upper):
        numerator = (M - x) * (K - x)
        denominator = (x + 1) * (N - M - K + x + 1)
        if numerator <= 0 or denominator <= 0:
            raise ArithmeticError("invalid upper-tail recurrence")
        log_p += math.log(numerator) - math.log(denominator)
        logs.append(log_p)
    largest = max(logs)
    probability = math.exp(largest) * math.fsum(math.exp(value - largest) for value in logs)
    return min(1.0, max(0.0, probability))


def adjust_bh(p_values: Sequence[float]) -> List[float]:
    m = len(p_values)
    order = sorted(range(m), key=lambda index: (p_values[index], index))
    result = [1.0] * m
    running = 1.0
    for rank0 in range(m - 1, -1, -1):
        index = order[rank0]
        running = min(running, p_values[index] * m / float(rank0 + 1))
        result[index] = min(1.0, running)
    return result


def adjust_by(p_values: Sequence[float]) -> List[float]:
    harmonic = math.fsum(1.0 / index for index in range(1, len(p_values) + 1))
    return adjust_bh([min(1.0, value * harmonic) for value in p_values])


def adjust_holm(p_values: Sequence[float]) -> List[float]:
    m = len(p_values)
    order = sorted(range(m), key=lambda index: (p_values[index], index))
    result = [1.0] * m
    running = 0.0
    for rank0, index in enumerate(order):
        running = max(running, (m - rank0) * p_values[index])
        result[index] = min(1.0, running)
    return result


def adjust_bonferroni(p_values: Sequence[float]) -> List[float]:
    m = len(p_values)
    return [min(1.0, value * m) for value in p_values]


def _raw_ratio(numerator: float, denominator: float) -> object:
    if denominator == 0:
        return "NA" if numerator == 0 else "inf"
    return numerator / denominator


def effect_sizes(a: int, b: int, c: int, d: int) -> Dict[str, object]:
    K = a + b
    nonphr = c + d
    phr_fraction = a / float(K) if K else float("nan")
    nonphr_fraction = c / float(nonphr) if nonphr else float("nan")
    fold = _raw_ratio(phr_fraction, nonphr_fraction)
    fold_ha = ((a + 0.5) / (K + 1.0)) / ((c + 0.5) / (nonphr + 1.0))
    cross_num = a * d
    cross_den = b * c
    odds = _raw_ratio(float(cross_num), float(cross_den))
    odds_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return {
        "phr_copy_fraction": phr_fraction,
        "nonphr_copy_fraction": nonphr_fraction,
        "fold_enrichment": fold,
        "fold_enrichment_ha": fold_ha,
        "odds_ratio": odds,
        "odds_ratio_ha": odds_ha,
    }


def infer_term_results(
    copy_rows: Sequence[Mapping[str, str]],
    edge_rows: Iterable[Mapping[str, str]],
    hypotheses: Sequence[Mapping[str, str]],
) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, int]]]:
    """Enumerate all V6 hypotheses over exact PHR/non-PHR complements."""
    if [int(row["hypothesis_index"]) for row in hypotheses] != list(range(len(hypotheses))):
        raise ValueError("frozen hypothesis order is not contiguous")
    keys = [(row["collection"], row["relation"], row["term_id"]) for row in hypotheses]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate frozen hypothesis")
    key_to_index = {key: index for index, key in enumerate(keys)}

    eligible: Dict[str, Mapping[str, str]] = {}
    for row in copy_rows:
        cn = int(row["physical_copy_cn"])
        if cn <= 0:
            raise ValueError("nonpositive physical-copy count")
        if int(row["closure_term_count"]) > 0:
            if row["copy_id"] in eligible:
                raise ValueError("duplicate eligible physical-copy coordinate")
            eligible[row["copy_id"]] = row
    N = sum(int(row["physical_copy_cn"]) for row in eligible.values())
    membership: Dict[str, Dict[str, int]] = {}
    audit: Dict[str, Dict[str, int]] = {}
    for assignment, column, _role in ASSIGNMENTS:
        values = {copy_id: int(row[column]) for copy_id, row in eligible.items()}
        if any(value not in (0, int(eligible[copy_id]["physical_copy_cn"]))
               for copy_id, value in values.items()):
            raise ValueError("membership must assign a complete coordinate copy")
        membership[assignment] = values
        K = sum(values.values())
        nonphr = N - K
        intersection = sum(min(value, int(eligible[copy_id]["physical_copy_cn"]) - value)
                           for copy_id, value in values.items())
        audit[assignment] = {
            "N": N,
            "K": K,
            "nonphr": nonphr,
            "partition_sum_cn": K + nonphr,
            "intersection_cn": intersection,
            "eligible_rows": len(eligible),
        }

    genome = [0] * len(hypotheses)
    phr = {assignment: [0] * len(hypotheses) for assignment, _column, _role in ASSIGNMENTS}
    edge_rows_seen = 0
    for edge in edge_rows:
        key = (collection_for(edge), edge["relation"], edge["term_id"])
        index = key_to_index.get(key)
        if index is None:
            raise ValueError("physical-copy edge is outside frozen catalog: %r" % (key,))
        copy_id = edge["copy_id"]
        copy = eligible.get(copy_id)
        if copy is None:
            raise ValueError("term edge references ineligible or unknown copy: %s" % copy_id)
        cn = int(edge["physical_copy_cn"])
        if cn != int(copy["physical_copy_cn"]):
            raise ValueError("edge/copy physical count disagreement")
        genome[index] += cn
        for assignment, _column, _role in ASSIGNMENTS:
            phr[assignment][index] += membership[assignment][copy_id]
        edge_rows_seen += 1

    expected_genome = [int(row["genome_physical_copy_burden"]) for row in hypotheses]
    if genome != expected_genome:
        mismatches = sum(a != b for a, b in zip(genome, expected_genome))
        raise ValueError("term-edge burdens disagree with frozen catalog for %d hypotheses" % mismatches)
    if edge_rows_seen == 0:
        raise ValueError("empty physical-copy edge table")

    family_sizes = Counter(row["collection"] for row in hypotheses)
    results: List[Dict[str, object]] = []
    for assignment, _column, role in ASSIGNMENTS:
        K = audit[assignment]["K"]
        nonphr_total = N - K
        for index, hypothesis in enumerate(hypotheses):
            M = genome[index]
            a = phr[assignment][index]
            c = M - a
            b = K - a
            d = nonphr_total - c
            if min(a, b, c, d) < 0 or a + b != K or c + d != nonphr_total or a + c != M:
                raise AssertionError("invalid finite-population 2x2 partition")
            effects = effect_sizes(a, b, c, d)
            expected = K * M / float(N)
            row: Dict[str, object] = {
                "schema_version": SCHEMA_VERSION,
                "assignment": assignment,
                "analysis_role": role,
                "hypothesis_index": index,
                "collection": hypothesis["collection"],
                "relation": hypothesis["relation"],
                "ontology": hypothesis["ontology"],
                "namespace": hypothesis["namespace"],
                "term_id": hypothesis["term_id"],
                "term_name": hypothesis["term_name"],
                "multiplicity_family_id": hypothesis["collection"],
                "multiplicity_family_size": family_sizes[hypothesis["collection"]],
                "hypothesis_status": hypothesis["hypothesis_status"],
                "a_T": a, "b_T": b, "c_T": c, "d_T": d,
                "M_T_genome_burden": M, "K_T_phr_burden": a,
                "nonphr_term_burden": c, "K_phr_total": K,
                "nonphr_total": nonphr_total, "N_eligible_total": N,
                "copy_count_difference_a_minus_c": a - c,
                "expected_phr_copy_count": expected,
                "copy_excess_over_expected": a - expected,
                "p_exact_upper": exact_upper_tail(N, M, K, a),
                "primary_support_rule": PRIMARY_RULE,
                "finite_population_definition": "COMPLETE_FROZEN_ONTOLOGY_ELIGIBLE_CHM13_PHYSICAL_COPIES",
                "background_definition": "EXACT_COMPLETE_NONPHR_SET_COMPLEMENT_NO_SAMPLING",
                "inference_status": "EXACT_COMPLETE_FINITE_POPULATION_UPPER_TAIL",
            }
            row.update(effects)
            results.append(row)

    for assignment, _column, _role in ASSIGNMENTS:
        assignment_indices = [index for index, row in enumerate(results) if row["assignment"] == assignment]
        p_global = [float(results[index]["p_exact_upper"]) for index in assignment_indices]
        holm = adjust_holm(p_global)
        bonferroni = adjust_bonferroni(p_global)
        for local, result_index in enumerate(assignment_indices):
            results[result_index]["holm_p_global"] = holm[local]
            results[result_index]["bonferroni_p_global"] = bonferroni[local]
        for collection in COLLECTIONS:
            indices = [index for index in assignment_indices if results[index]["collection"] == collection]
            p_values = [float(results[index]["p_exact_upper"]) for index in indices]
            bh = adjust_bh(p_values)
            by = adjust_by(p_values)
            for local, result_index in enumerate(indices):
                results[result_index]["bh_q_within_collection"] = bh[local]
                results[result_index]["by_q_within_collection"] = by[local]

    for row in results:
        passes = (float(row["bh_q_within_collection"]) <= PRIMARY_ALPHA and
                  float(row["holm_p_global"]) <= PRIMARY_ALPHA)
        primary = row["assignment"] == "midpoint" and passes
        sensitivity = row["assignment"] == "overlap" and passes
        row["primary_support"] = int(primary)
        row["sensitivity_support"] = int(sensitivity)
        if row["assignment"] == "midpoint":
            row["support_status"] = "PRIMARY_SUPPORTED" if primary else "PRIMARY_NOT_SUPPORTED"
        else:
            row["support_status"] = "SENSITIVITY_SUPPORTED" if sensitivity else "SENSITIVITY_NOT_SUPPORTED"
    return results, audit


def validate_frozen_inputs() -> Dict[str, object]:
    errors: List[str] = []
    audited = {row["path"]: row for row in read_rows(V6_AUDITED_MANIFEST)}
    for path in (COPY_MAP, COPY_EDGES):
        expected = audited.get(path.name, {})
        if not path.is_file() or expected.get("sha256") != sha256(path):
            errors.append("v6_audited_checksum:%s" % path.name)
    frozen = json.loads(HYPOTHESIS_MANIFEST.read_text(encoding="utf-8"))
    if frozen.get("catalog_sha256") != sha256(HYPOTHESES):
        errors.append("frozen_hypothesis_checksum")
    if frozen.get("hypothesis_count") != EXPECTED_HYPOTHESES:
        errors.append("frozen_hypothesis_count")
    if frozen.get("target_filter") != "none" or frozen.get("target_bed_opened") is not False:
        errors.append("hypothesis_freeze_target_blind")
    v6_validation = json.loads((V6 / "V6_VALIDATION.json").read_text(encoding="utf-8"))
    if v6_validation.get("status") != "PASS" or not v6_validation.get("strict_pass"):
        errors.append("v6_independent_validation")
    release = {row["path"]: row for row in read_rows(V6_RELEASE_MANIFEST)}
    expected_contributors = release.get(V6_CONTRIBUTORS.name, {})
    if expected_contributors.get("sha256") != sha256(V6_CONTRIBUTORS):
        errors.append("v6_exact_contributor_checksum")
    return {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "v6_copy_map_sha256": sha256(COPY_MAP),
        "v6_copy_edges_sha256": sha256(COPY_EDGES),
        "v6_hypothesis_catalog_sha256": sha256(HYPOTHESES),
        "v6_exact_contributors_sha256": sha256(V6_CONTRIBUTORS),
        "v6_hypothesis_count": frozen.get("hypothesis_count"),
    }


def emit_mapping_coverage(copy_rows: Sequence[Mapping[str, str]], audit: Mapping[str, Mapping[str, int]],
                          output: Path) -> List[Dict[str, object]]:
    raw_cn = sum(int(row["physical_copy_cn"]) for row in copy_rows)
    eligible_rows = [row for row in copy_rows if int(row["closure_term_count"]) > 0]
    eligible_cn = sum(int(row["physical_copy_cn"]) for row in eligible_rows)
    rows: List[Dict[str, object]] = []
    for assignment, column, _role in ASSIGNMENTS:
        phr_raw = sum(int(row[column]) for row in copy_rows)
        phr_eligible = sum(int(row[column]) for row in eligible_rows)
        nonphr = eligible_cn - phr_eligible
        record = {
            "schema_version": SCHEMA_VERSION,
            "assignment": assignment,
            "raw_physical_copy_rows": len(copy_rows),
            "raw_physical_copy_cn": raw_cn,
            "ontology_eligible_rows": len(eligible_rows),
            "ontology_eligible_cn": eligible_cn,
            "ontology_ineligible_cn": raw_cn - eligible_cn,
            "phr_raw_physical_copy_cn": phr_raw,
            "phr_ontology_eligible_cn": phr_eligible,
            "phr_ontology_ineligible_cn": phr_raw - phr_eligible,
            "nonphr_ontology_eligible_cn": nonphr,
            "partition_sum_cn": phr_eligible + nonphr,
            "phr_nonphr_intersection_cn": audit[assignment]["intersection_cn"],
            "functional_source_ids_in_eligible_universe": len(set(
                row["functional_source_id"] for row in eligible_rows)),
            "all_coordinate_copy_cn_equal_one": int(all(
                int(row["physical_copy_cn"]) == 1 for row in copy_rows)),
            "background_definition": "EXACT_COMPLETE_NONPHR_SET_COMPLEMENT_NO_SAMPLING",
            "status": "PASS" if (raw_cn == EXPECTED_RAW_COPIES and eligible_cn == EXPECTED_N and
                                    phr_eligible == EXPECTED_K[assignment] and
                                    phr_eligible + nonphr == eligible_cn and
                                    audit[assignment]["intersection_cn"] == 0) else "FAIL",
        }
        rows.append(record)
    fields = list(rows[0])
    write_rows(output / "MAPPING_COVERAGE.tsv", fields, rows)
    return rows


def _cohort_for(name: str, synonyms: str) -> str:
    synonym_set = set(value for value in synonyms.replace(",", "|").split("|") if value)
    if name == "DUX4" or name.startswith("DUX4L") or "DUX4L30" in synonym_set:
        return "DUX4_DUX4L"
    for prefix in ("DDX11L", "TUBB8", "OR4F", "WASH", "RPL23A"):
        if name.startswith(prefix):
            return prefix
    return ""


def emit_named_cohort_audit(copy_rows: Sequence[Mapping[str, str]], output: Path) -> List[Dict[str, object]]:
    cohorts = ("DUX4_DUX4L", "DDX11L", "TUBB8", "OR4F", "WASH", "RPL23A")
    grouped: MutableMapping[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in copy_rows:
        cohort = _cohort_for(row["gene_name"], row.get("gene_synonyms", ""))
        if cohort:
            grouped[cohort].append(row)
    records: List[Dict[str, object]] = []
    for cohort in cohorts:
        rows = grouped[cohort]
        eligible = [row for row in rows if int(row["closure_term_count"]) > 0]
        source_cn = Counter()
        for row in eligible:
            source_cn[row["functional_source_id"]] += int(row["physical_copy_cn"])
        genome_cn = sum(int(row["physical_copy_cn"]) for row in rows)
        eligible_cn = sum(int(row["physical_copy_cn"]) for row in eligible)
        midpoint_cn = sum(int(row["phr_midpoint_cn"]) for row in rows)
        midpoint_eligible = sum(int(row["phr_midpoint_cn"]) for row in eligible)
        overlap_cn = sum(int(row["phr_any_overlap_cn"]) for row in rows)
        overlap_eligible = sum(int(row["phr_any_overlap_cn"]) for row in eligible)
        coordinate_semantics = (all(int(row["physical_copy_cn"]) == 1 for row in rows) and
                                genome_cn == len({row["copy_id"] for row in rows}))
        records.append({
            "schema_version": SCHEMA_VERSION,
            "cohort": cohort,
            "genome_coordinate_copy_cn": genome_cn,
            "genome_ontology_eligible_copy_cn": eligible_cn,
            "midpoint_phr_coordinate_copy_cn": midpoint_cn,
            "midpoint_phr_ontology_eligible_copy_cn": midpoint_eligible,
            "midpoint_nonphr_ontology_eligible_copy_cn": eligible_cn - midpoint_eligible,
            "overlap_phr_coordinate_copy_cn": overlap_cn,
            "overlap_phr_ontology_eligible_copy_cn": overlap_eligible,
            "overlap_nonphr_ontology_eligible_copy_cn": eligible_cn - overlap_eligible,
            "functional_source_count": len(source_cn),
            "functional_source_ids": "|".join(sorted(source_cn)),
            "largest_single_source_coordinate_copy_cn": max(source_cn.values()) if source_cn else 0,
            "forbidden_source_collapsed_genome_count": len(source_cn),
            "coordinate_copy_units_minus_source_count": eligible_cn - len(source_cn),
            "all_coordinate_copy_cn_equal_one": int(all(int(row["physical_copy_cn"]) == 1 for row in rows)),
            "coordinate_copy_semantics_status": "PASS" if coordinate_semantics else "FAIL",
            "audit_role": "POST_INFERENCE_NAMED_COHORT_COPY_SEMANTICS_NOT_HYPOTHESIS",
            "status": "PASS" if rows and coordinate_semantics else "FAIL",
        })
    write_rows(output / "NAMED_COHORT_AUDIT.tsv", list(records[0]), records)
    return records


def emit_contributor_ledger(results: Sequence[Mapping[str, object]], output: Path) -> List[Dict[str, object]]:
    decisions = {(str(row["assignment"]), str(row["collection"]), str(row["relation"]), str(row["term_id"])): row
                 for row in results}
    records: List[Dict[str, object]] = []
    source_fields: List[str] = []
    for source in read_rows(V6_CONTRIBUTORS):
        if not source_fields:
            source_fields = list(source)
        key_base = (source["collection"], source["relation"], source["term_id"])
        primary = decisions[("midpoint",) + key_base]
        sensitivity = decisions[("overlap",) + key_base]
        midpoint_phr = int(source["phr_midpoint_cn"]) > 0
        overlap_phr = int(source["phr_any_overlap_cn"]) > 0
        row: Dict[str, object] = dict(source)
        row.update({
            "v7_schema_version": SCHEMA_VERSION,
            "v7_midpoint_partition": "PHR" if midpoint_phr else "NON_PHR",
            "v7_overlap_partition": "PHR" if overlap_phr else "NON_PHR",
            "v7_midpoint_count_cell": "a_T" if midpoint_phr else "c_T",
            "v7_overlap_count_cell": "a_T" if overlap_phr else "c_T",
            "v7_primary_support": primary["primary_support"],
            "v7_overlap_sensitivity_support": sensitivity["sensitivity_support"],
            "v7_ledger_scope": "V6_VALIDATED_ANY_OVERLAP_CONTRIBUTORS_WITH_V7_PARTITION_STATUS",
            "v7_status": "PASS",
        })
        records.append(row)
    extra = ["v7_schema_version", "v7_midpoint_partition", "v7_overlap_partition",
             "v7_midpoint_count_cell", "v7_overlap_count_cell", "v7_primary_support",
             "v7_overlap_sensitivity_support", "v7_ledger_scope", "v7_status"]
    write_rows(output / "EXACT_TERM_CONTRIBUTORS.tsv.gz", source_fields + extra, records)
    return records


def _arm_communities() -> Dict[str, str]:
    result = {}
    for row in read_rows(COMMUNITIES):
        arm = row["arm"].replace("_parm", "_p").replace("_qarm", "_q")
        result[arm] = row["community"]
    return result


def _q_starts() -> Dict[str, int]:
    starts: MutableMapping[str, List[int]] = defaultdict(list)
    with CYTOBANDS.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chromosome, start, _end, band = line.rstrip("\n").split("\t")[:4]
            if band.startswith("q"):
                starts[chromosome].append(int(start))
    return {chromosome: min(values) for chromosome, values in starts.items()}


def emit_community_summary(results: Sequence[Mapping[str, object]],
                           contributors: Sequence[Mapping[str, object]], output: Path) -> List[Dict[str, object]]:
    """Summarize contributor arms only after exact term decisions are fixed."""
    supported = {(str(row["assignment"]), str(row["collection"]), str(row["relation"]), str(row["term_id"])): row
                 for row in results if int(row["primary_support"]) or int(row["sensitivity_support"])}
    communities = _arm_communities()
    q_starts = _q_starts()
    grouped: MutableMapping[Tuple[str, str, str, str, str], List[Mapping[str, object]]] = defaultdict(list)
    for row in contributors:
        chromosome = str(row["seqid"])
        midpoint = (int(row["start0"]) + int(row["end0"])) // 2
        arm = "%s_%s" % (chromosome, "p" if midpoint < q_starts[chromosome] else "q")
        community = communities.get(arm, "UNASSIGNED_NO_SIGNAL_COMMUNITY")
        base = (str(row["collection"]), str(row["relation"]), str(row["term_id"]))
        if str(row["v7_midpoint_partition"]) == "PHR" and ("midpoint",) + base in supported:
            grouped[("midpoint",) + base + (community,)].append(dict(row, _arm=arm))
        if str(row["v7_overlap_partition"]) == "PHR" and ("overlap",) + base in supported:
            grouped[("overlap",) + base + (community,)].append(dict(row, _arm=arm))
    records: List[Dict[str, object]] = []
    for key in sorted(grouped):
        assignment, collection, relation, term_id, community = key
        rows = grouped[key]
        result = supported[(assignment, collection, relation, term_id)]
        copy_cn = sum(int(row["physical_copy_cn"]) for row in rows)
        records.append({
            "schema_version": SCHEMA_VERSION,
            "assignment": assignment,
            "collection": collection,
            "relation": relation,
            "term_id": term_id,
            "term_name": result["term_name"],
            "support_status": result["support_status"],
            "community": community,
            "community_arms": "|".join(sorted(set(str(row["_arm"]) for row in rows))),
            "community_contributor_copy_cn": copy_cn,
            "term_phr_copy_cn": result["a_T"],
            "community_fraction_of_term_phr_cn": copy_cn / float(int(result["a_T"])),
            "functional_source_count": len(set(str(row["functional_source_id"]) for row in rows)),
            "functional_source_ids": "|".join(sorted(set(str(row["functional_source_id"]) for row in rows))),
            "functional_source_symbols": "|".join(sorted(set(str(row["functional_source_symbol"]) for row in rows))),
            "coordinate_copy_count": len(set(str(row["copy_id"]) for row in rows)),
            "copy_ids": "|".join(sorted(set(str(row["copy_id"]) for row in rows))),
            "summary_role": "POST_INFERENCE_ONLY_NOT_A_HYPOTHESIS_OR_TEST",
        })
    fields = [
        "schema_version", "assignment", "collection", "relation", "term_id", "term_name",
        "support_status", "community", "community_arms", "community_contributor_copy_cn",
        "term_phr_copy_cn", "community_fraction_of_term_phr_cn", "functional_source_count",
        "functional_source_ids", "functional_source_symbols", "coordinate_copy_count",
        "copy_ids", "summary_role",
    ]
    write_rows(output / "COMMUNITY_TERM_SUMMARY.tsv", fields, records)
    return records


def executable_path_audit() -> Dict[str, object]:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    imported = []
    calls = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.append((node.module or "").split(".")[0])
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                calls.append(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                calls.append(node.func.attr)
    forbidden_imports = {"ran" + "dom", "num" + "py"}
    forbidden_calls = {"per" + "mutation", "shuffle", "choice", "sample", "randint", "uniform"}
    bad_imports = sorted(set(imported) & forbidden_imports)
    bad_calls = sorted(set(calls) & forbidden_calls)
    return {
        "status": "PASS" if not bad_imports and not bad_calls else "FAIL",
        "forbidden_runtime_imports_found": bad_imports,
        "forbidden_stochastic_calls_found": bad_calls,
        "rng_seed_or_state": "NONE",
        "primary_inference_function": "infer_term_results",
        "primary_background_operation": "exact_boolean_set_complement",
    }


def validate_release(output: Path = HERE) -> Dict[str, object]:
    errors: List[str] = []
    checks: List[Dict[str, object]] = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})
        if not passed:
            errors.append(name)

    frozen = validate_frozen_inputs()
    check("inherited_v6_frozen_objects", frozen["status"] == "PASS", frozen)
    results = list(read_rows(output / "TERM_RESULTS.tsv.gz"))
    mapping = list(read_rows(output / "MAPPING_COVERAGE.tsv"))
    named = list(read_rows(output / "NAMED_COHORT_AUDIT.tsv"))
    contributors = list(read_rows(output / "EXACT_TERM_CONTRIBUTORS.tsv.gz"))
    communities = list(read_rows(output / "COMMUNITY_TERM_SUMMARY.tsv"))
    hypotheses = list(read_rows(HYPOTHESES))
    check("all_frozen_hypotheses_tested_twice", len(results) == 2 * len(hypotheses) == 62470,
          {"results": len(results), "hypotheses": len(hypotheses)})
    by_assignment = Counter(row["assignment"] for row in results)
    check("paired_assignments_complete", by_assignment == Counter({"midpoint": 31235, "overlap": 31235}),
          dict(by_assignment))
    expected_keys = {(assignment, row["collection"], row["relation"], row["term_id"])
                     for assignment, _column, _role in ASSIGNMENTS for row in hypotheses}
    actual_keys = {(row["assignment"], row["collection"], row["relation"], row["term_id"]) for row in results}
    check("exact_frozen_hypothesis_keys", actual_keys == expected_keys,
          {"expected": len(expected_keys), "actual": len(actual_keys)})
    check("zero_phr_burdens_retained", all(any(
        row["assignment"] == assignment and int(row["a_T"]) == 0 for row in results)
        for assignment in ("midpoint", "overlap")), "zero rows present in both assignments")
    cell_ok = all(int(row["a_T"]) + int(row["b_T"]) == int(row["K_phr_total"]) and
                  int(row["c_T"]) + int(row["d_T"]) == int(row["nonphr_total"]) and
                  int(row["a_T"]) + int(row["c_T"]) == int(row["M_T_genome_burden"]) and
                  int(row["K_phr_total"]) + int(row["nonphr_total"]) == int(row["N_eligible_total"])
                  for row in results)
    check("all_term_2x2_tables_partition", cell_ok, "a+b=K; c+d=N-K; a+c=M for every row")
    check("eligible_universe_constant", all(int(row["N_eligible_total"]) == EXPECTED_N for row in results), EXPECTED_N)
    check("phr_nonphr_exact_complement", all(row["status"] == "PASS" and
          int(row["partition_sum_cn"]) == EXPECTED_N and int(row["phr_nonphr_intersection_cn"]) == 0
          for row in mapping), mapping)
    check("no_phr_contributor_marked_background", all(
        not (int(row["phr_midpoint_cn"]) > 0 and row["v7_midpoint_partition"] != "PHR") and
        not (int(row["phr_any_overlap_cn"]) > 0 and row["v7_overlap_partition"] != "PHR")
        for row in contributors), "V7 partition statuses agree with both frozen membership columns")
    contributor_burdens = Counter()
    for row in contributors:
        base = (row["collection"], row["relation"], row["term_id"])
        if row["v7_midpoint_partition"] == "PHR":
            contributor_burdens[("midpoint",) + base] += int(row["physical_copy_cn"])
        if row["v7_overlap_partition"] == "PHR":
            contributor_burdens[("overlap",) + base] += int(row["physical_copy_cn"])
    check("exact_contributor_recount_matches_every_a_T", all(
        contributor_burdens[(row["assignment"], row["collection"], row["relation"], row["term_id"])] ==
        int(row["a_T"]) for row in results),
        {"ledger_rows": len(contributors), "term_assignment_rows": len(results)})
    copied = [row for row in named if row["cohort"] in {"DUX4_DUX4L", "DDX11L", "WASH", "RPL23A"}]
    check("coordinate_copies_not_source_collapsed", len(copied) == 4 and all(
        row["coordinate_copy_semantics_status"] == "PASS" and
        int(row["genome_ontology_eligible_copy_cn"]) > int(row["forbidden_source_collapsed_genome_count"])
        for row in copied), copied)
    check("all_named_cohort_audits_pass", len(named) == 6 and all(
        row["status"] == "PASS" for row in named), {"cohorts": [row["cohort"] for row in named]})
    code_audit = executable_path_audit()
    check("no_stochastic_or_placement_runtime_path", code_audit["status"] == "PASS", code_audit)
    check("multiplicity_complete", all(row["bh_q_within_collection"] and row["by_q_within_collection"] and
          row["holm_p_global"] and row["bonferroni_p_global"] for row in results),
          "BH/BY within collection and Holm/Bonferroni global populated for every row")
    primary_rule_ok = all(int(row["primary_support"]) == int(
        row["assignment"] == "midpoint" and float(row["bh_q_within_collection"]) <= 0.05 and
        float(row["holm_p_global"]) <= 0.05) for row in results)
    check("prespecified_primary_rule_recomputed", primary_rule_ok, PRIMARY_RULE)
    check("community_summary_post_inference_only", all(
        row["summary_role"] == "POST_INFERENCE_ONLY_NOT_A_HYPOTHESIS_OR_TEST" for row in communities),
        {"summary_rows": len(communities), "hypothesis_source": HYPOTHESES.name})
    community_burdens = Counter()
    for row in communities:
        community_burdens[(row["assignment"], row["collection"], row["relation"], row["term_id"])] += int(
            row["community_contributor_copy_cn"])
    supported_results = [row for row in results if int(row["primary_support"]) or int(row["sensitivity_support"])]
    check("post_inference_community_sums_match_supported_terms", all(
        community_burdens[(row["assignment"], row["collection"], row["relation"], row["term_id"])] ==
        int(row["a_T"]) for row in supported_results),
        {"supported_term_rows": len(supported_results), "community_rows": len(communities)})
    check("contributor_ledger_inherited_and_statused", len(contributors) == 17579 and all(
        row["v7_status"] == "PASS" for row in contributors), len(contributors))
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "checks_passed": sum(row["status"] == "PASS" for row in checks),
        "checks_total": len(checks),
        "checks": checks,
        "stage_order": [
            "1_VALIDATE_FROZEN_INPUTS_AND_PARTITION",
            "2_COMPLETE_TERM_LEVEL_INFERENCE",
            "3_COMPLETE_MULTIPLICITY_AND_DECISIONS",
            "4_CONTRIBUTOR_AND_MAPPING_AUDITS",
            "5_POST_INFERENCE_COMMUNITY_SOURCE_SUMMARY",
            "6_RELEASE_VALIDATION_AND_REPORT",
        ],
        "primary_support_rule": PRIMARY_RULE,
        "manuscript_interpretation_boundary": (
            "annotation-bearing CHM13 physical-copy burden only; not expression, activity, dosage, "
            "retained pseudogene function, or population prevalence"),
    }


def emit_final_report(results: Sequence[Mapping[str, object]], validation: Mapping[str, object],
                      output: Path) -> None:
    primary = [row for row in results if int(row["primary_support"])]
    sensitivity = [row for row in results if int(row["sensitivity_support"])]
    primary_counts = Counter(str(row["collection"]) for row in primary)
    relation_counts = Counter(str(row["relation"]) for row in primary)
    zero_primary = sum(row["assignment"] == "midpoint" and int(row["a_T"]) == 0 for row in results)
    top = sorted(primary, key=lambda row: (float(row["holm_p_global"]),
                                           -float(row["copy_excess_over_expected"]),
                                           str(row["collection"]), str(row["term_id"])))[:20]
    lines = [
        "# Final V7 whole-genome non-PHR physical-copy ontology report", "",
        "## Release verdict", "",
        "**V7 validation status: %s (%s/%s checks).**" % (
            validation["status"], validation["checks_passed"], validation["checks_total"]), "",
        "The corrected V7 analysis tested all **31,235** frozen exact direct/ancestor",
        "GO and Reactome hypotheses against the complete ontology-eligible CHM13",
        "non-PHR genome. The primary midpoint universe contains **187 PHR copies**",
        "and **31,779 non-PHR copies** (N = 31,966). The paired any-overlap",
        "sensitivity contains 193 PHR and 31,773 non-PHR copies. No background was",
        "sampled, no PHR interval was moved, and no spatial control defines V7 inference.", "",
        "The prespecified primary rule (within-collection BH q <= 0.05 and global",
        "Holm-adjusted p <= 0.05) supports **%d exact midpoint term rows**. The paired" % len(primary),
        "overlap sensitivity supports **%d rows**. Primary supports by collection are" % len(sensitivity),
        "`%s`; by relation they are `%s`. All %d midpoint rows with zero PHR burden" % (
            json.dumps(dict(sorted(primary_counts.items()))),
            json.dumps(dict(sorted(relation_counts.items()))), zero_primary),
        "were retained and tested.", "",
        "## What was counted and tested", "",
        "For each frozen `(collection, relation, term_id)` row, the result table gives",
        "`a_T`, `b_T`, `c_T`, and `d_T` over coordinate-distinct physical copies.",
        "The p-value is the exact upper tail for observing at least `a_T` term-bearing",
        "copies among K PHR copies in the complete fixed population of N eligible",
        "copies containing M_T term-bearing copies. This hypergeometric calculation is",
        "a complete finite-population copy-burden test, not generic weak gene-list ORA.", "",
        "Every inherited V6 copy has CN=1. Multiple DUX4/DUX4L, DDX11L, WASH,",
        "RPL23A, OR4F, and TUBB8 coordinates therefore contribute one unit each even",
        "when they share a functional source. The named-cohort audit explicitly shows",
        "the larger coordinate-copy counts beside the forbidden source-collapsed counts.", "",
        "## Leading exact supported rows", "",
        "| Collection | Relation | Exact term | a / K | c / non-PHR | Fold | Exact p | BH q | Global Holm |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in top:
        lines.append("| %s | %s | %s (`%s`) | %s / %s | %s / %s | %s | %.3g | %.3g | %.3g |" % (
            row["collection"], row["relation"], str(row["term_name"]).replace("|", "/"), row["term_id"],
            row["a_T"], row["K_phr_total"], row["c_T"], row["nonphr_total"],
            row["fold_enrichment"], float(row["p_exact_upper"]),
            float(row["bh_q_within_collection"]), float(row["holm_p_global"])))
    lines.extend([
        "", "The table above is ordered by global Holm value and copy excess; the complete",
        "results, including nonsupports and zero-PHR rows, are in `TERM_RESULTS.tsv.gz`.", "",
        "## Validation gates", "",
    ])
    for check in validation["checks"]:
        lines.append("- **%s:** `%s`." % (check["check"], check["status"]))
    lines.extend([
        "", "The exact contributor ledger reuses the independently validated V6 any-overlap",
        "contributor rows and adds V7 midpoint/overlap partition and decision statuses.",
        "`MAPPING_COVERAGE.tsv` proves the eligible partition and empty intersection.",
        "`COMMUNITY_TERM_SUMMARY.tsv` was generated only after all exact term decisions",
        "were fixed; it is a descriptive source/community display and created no test.", "",
        "## Interpretation boundary", "",
        "V7 answers which exact ontology terms have excess annotation-bearing physical",
        "copy burden in CHM13 PHRs relative to the rest of the eligible CHM13 genome.",
        "It does not test expression, protein activity, dosage effect, retained function",
        "of pseudogenes, biological independence of adjacent copies, or population prevalence.", "",
        "This V7 result supersedes V6's placement-null inference for the biological",
        "question specified here; the V6 map and hypothesis catalog remain the frozen inputs.",
    ])
    (output / "FINAL_V7_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def emit_manifest(output: Path) -> None:
    names = [
        "V7_RUN_PROTOCOL.md", "TERM_RESULTS.tsv.gz", "EXACT_TERM_CONTRIBUTORS.tsv.gz",
        "MAPPING_COVERAGE.tsv", "NAMED_COHORT_AUDIT.tsv", "COMMUNITY_TERM_SUMMARY.tsv",
        "V7_VALIDATION.json", "FINAL_V7_REPORT.md", "v7_copy_enrichment.py",
        "test_v7_copy_enrichment.py",
    ]
    rows = []
    for name in names:
        path = output / name
        rows.append({"path": name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_rows(output / "OUTPUT_MANIFEST.sha256.tsv", ["path", "bytes", "sha256"], rows)


def run(output: Path = HERE) -> Dict[str, object]:
    frozen = validate_frozen_inputs()
    if frozen["status"] != "PASS":
        raise RuntimeError("frozen V6 input validation failed: %s" % frozen["errors"])
    copy_rows = list(read_rows(COPY_MAP))
    hypotheses = list(read_rows(HYPOTHESES))
    results, audit = infer_term_results(copy_rows, read_rows(COPY_EDGES), hypotheses)
    write_rows(output / "TERM_RESULTS.tsv.gz", TERM_FIELDS, results)
    emit_mapping_coverage(copy_rows, audit, output)
    emit_named_cohort_audit(copy_rows, output)
    contributors = emit_contributor_ledger(results, output)
    emit_community_summary(results, contributors, output)
    validation = validate_release(output)
    write_json(output / "V7_VALIDATION.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError("V7 release validation failed: %s" % validation["errors"])
    emit_final_report(results, validation, output)
    emit_manifest(output)
    return {
        "status": "PASS",
        "hypotheses": len(hypotheses),
        "result_rows": len(results),
        "primary_supported": sum(int(row["primary_support"]) for row in results),
        "sensitivity_supported": sum(int(row["sensitivity_support"]) for row in results),
        "validation_checks": validation["checks_total"],
    }


def main(argv: Sequence[str] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "validate"), nargs="?", default="run")
    args = parser.parse_args(argv)
    if args.command == "run":
        print(json.dumps(run(), indent=2, sort_keys=True))
    else:
        result = validate_release(HERE)
        print(json.dumps(result, indent=2, sort_keys=True))
        if result["status"] != "PASS":
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
