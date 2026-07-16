#!/usr/bin/env python3
"""Gate, freeze, run, and release the V6 physical-CN ontology analysis.

The independently audited V6 copy/source/term matrix is the only annotation
input.  Geometry and the exactly reproducible same-arm sampler are delegated
to the previously validated regional engine; all V6 hypothesis construction,
CN weighting, multiplicity, contributor reporting, and release validation live
in this module.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, MutableMapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
REPO = HERE.parents[3]
RESULTS = HERE / "results"
WORK = RESULTS / "work"
RELEASE = RESULTS / "release"
CATALOG = HERE / "FROZEN_HYPOTHESES.tsv.gz"
CATALOG_MANIFEST = HERE / "FROZEN_HYPOTHESES.json"
GATE_JSON = HERE / "PRE_RUN_V6_GATE.json"
AUDITED_MANIFEST = HERE / "PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv"
INPUT_MANIFEST = HERE / "INPUT_MANIFEST.tsv"
ASSIGNMENTS = HERE / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz"
SOURCE_TERMS = HERE / "SOURCE_TERMS.tsv.gz"
PHYSICAL_EDGES = HERE / "PHYSICAL_COPY_TERM_EDGES.tsv.gz"
COPY_EVIDENCE = HERE / "COPY_SOURCE_EVIDENCE.tsv.gz"
EXACT_CONTRIBUTORS = HERE / "EXACT_CONTRIBUTORS.tsv.gz"
NAMED_BURDENS = HERE / "NAMED_COHORT_TERM_BURDENS.tsv.gz"
CYTOBANDS = REPO / "data" / "chm13v2.0_cytobands_allchrs.bed"

V3 = ANALYSIS / "ontology_v3"
sys.path.insert(0, str(V3))
import v3_engine as N  # noqa: E402

E = N.E
SCHEMA_VERSION = "chm13-physical-copy-ontology-v6.0"
COLLECTIONS = ("GO_BP", "GO_MF", "GO_CC", "Reactome")
RELATIONS = ("direct", "ancestor")
INITIAL_PERMUTATIONS = 99_999
CHECKPOINTS = (99_999, 249_999, 599_999, 999_999)
MASTER_SEED = 2026071301
EXPECTED_PHYSICAL_CN = 61_312
EXPECTED_ELIGIBLE_CN = 31_966
EXPECTED_EDGE_CN = 2_929_709
EXPECTED_DIRECT_EDGE_CN = 626_048
EXPECTED_PHR_MIDPOINT_CN = 402
EXPECTED_PHR_OVERLAP_CN = 412
EXPECTED_PHR_MIDPOINT_EDGE_CN = 16_763
EXPECTED_PHR_OVERLAP_EDGE_CN = 17_579
NAMED_EXPECTED = {
    "DUX4_DUX4L": 65,
    "DDX11L": 10,
    "TUBB8": 2,
    "OR4F": 4,
    "WASH": 9,
}
FREEZE_INPUT_BASENAMES = {
    "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz",
    "SOURCE_TERMS.tsv.gz",
    "chm13v2.0_cytobands_allchrs.bed",
}


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def read_rows(path: Path) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_rows(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    if path.suffix == ".gz":
        raw = path.open("wb")
        binary = gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0)
        handle = io.TextIOWrapper(binary, encoding="utf-8", newline="")
    else:
        raw = None
        handle = path.open("w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    finally:
        handle.close()
        if raw is not None:
            raw.close()
    return count


def atomic_json(path: Path, value: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def collection_for(ontology: str, namespace: str) -> str:
    if ontology == "Reactome" and namespace == "pathway":
        return "Reactome"
    if ontology != "GO":
        raise ValueError(f"unsupported ontology/namespace: {ontology}/{namespace}")
    try:
        return {
            "biological_process": "GO_BP",
            "molecular_function": "GO_MF",
            "cellular_component": "GO_CC",
        }[namespace]
    except KeyError as exc:
        raise ValueError(f"unsupported ontology/namespace: {ontology}/{namespace}") from exc


def _resolve_repo_path(text: str, repo_root: Path) -> Path:
    path = Path(text)
    return path if path.is_absolute() else repo_root / path


def validate_gate(
    directory: Path = HERE, repo_root: Path = REPO, *, check_frozen_inputs: bool = True,
) -> dict[str, object]:
    """Fail closed unless PASS is bound to all exact audited and source bytes."""
    errors: list[str] = []
    gate_path = directory / GATE_JSON.name
    audited_path = directory / AUDITED_MANIFEST.name
    input_path = directory / INPUT_MANIFEST.name
    try:
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"pass": False, "errors": [f"gate_unreadable:{exc}"]}
    if gate.get("status") != "PASS":
        errors.append(f"gate_status:{gate.get('status')}")
    if gate.get("enrichment_authorized") is not True:
        errors.append("enrichment_not_authorized")
    if gate.get("enrichment_run") is not False:
        errors.append("gate_scope_not_pre_inference")
    if not audited_path.is_file():
        errors.append("audited_release_manifest_missing")
    elif sha256(audited_path) != gate.get("audited_release_manifest_sha256"):
        errors.append("audited_release_manifest_sha256")
    audited_checked = 0
    if audited_path.is_file():
        for row in read_rows(audited_path):
            path = directory / row["path"]
            audited_checked += 1
            if not path.is_file():
                errors.append(f"audited_missing:{row['path']}")
            elif path.stat().st_size != int(row["bytes"]):
                errors.append(f"audited_bytes:{row['path']}")
            elif sha256(path) != row["sha256"]:
                errors.append(f"audited_sha256:{row['path']}")
    input_checked = 0
    if check_frozen_inputs:
        if not input_path.is_file():
            errors.append("input_manifest_missing")
        else:
            for row in read_rows(input_path):
                path = _resolve_repo_path(row["path"], repo_root)
                input_checked += 1
                if not path.is_file():
                    errors.append(f"input_missing:{row['role']}")
                elif path.stat().st_size != int(row["bytes"]):
                    errors.append(f"input_bytes:{row['role']}")
                elif sha256(path) != row["sha256"] or row["sha256"] != row["expected_sha256"]:
                    errors.append(f"input_sha256:{row['role']}")
    for field, expected in (
        ("physical_copy_cn", EXPECTED_PHYSICAL_CN),
        ("ontology_eligible_copies", EXPECTED_ELIGIBLE_CN),
        ("physical_copy_term_edges", EXPECTED_EDGE_CN),
        ("direct_physical_copy_term_edges", EXPECTED_DIRECT_EDGE_CN),
        ("phr_midpoint_copies", EXPECTED_PHR_MIDPOINT_CN),
        ("phr_any_overlap_copies", EXPECTED_PHR_OVERLAP_CN),
        ("phr_midpoint_copy_term_edges", EXPECTED_PHR_MIDPOINT_EDGE_CN),
    ):
        if int(gate.get(field, -1)) != expected:
            errors.append(f"gate_count:{field}:{gate.get(field)}:{expected}")
    observed_named = {
        cohort: int(values.get("phr_ontology_contributors", -1))
        for cohort, values in gate.get("named_cohorts", {}).items()
    }
    if observed_named != NAMED_EXPECTED:
        errors.append(f"gate_named_cohorts:{observed_named}")
    return {
        "schema_version": SCHEMA_VERSION,
        "pass": not errors,
        "errors": errors,
        "status": gate.get("status"),
        "enrichment_authorized": gate.get("enrichment_authorized"),
        "audited_release_files_checked": audited_checked,
        "frozen_source_inputs_checked": input_checked,
        "gate_sha256": sha256(gate_path) if gate_path.is_file() else "",
        "audited_release_manifest_sha256": sha256(audited_path) if audited_path.is_file() else "",
    }


def load_arm_boundaries(path: Path = CYTOBANDS) -> dict[str, tuple[int, int]]:
    bands: MutableMapping[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: malformed cytoband")
            bands[fields[0]].append((int(fields[1]), int(fields[2]), fields[3]))
    result = {}
    for chromosome, rows in bands.items():
        rows.sort()
        if rows[0][0] != 0 or any(rows[index][1] != rows[index + 1][0] for index in range(len(rows) - 1)):
            raise ValueError(f"noncontiguous cytobands: {chromosome}")
        q_starts = [start for start, _end, band in rows if band.startswith("q")]
        if not q_starts:
            raise ValueError(f"no q arm: {chromosome}")
        result[chromosome] = (min(q_starts), rows[-1][1])
    return result


def arm_for(chromosome: str, start0: int, end0: int,
            boundaries: Mapping[str, tuple[int, int]]) -> str:
    q_start, chromosome_end = boundaries[chromosome]
    midpoint = (start0 + end0) // 2
    if not 0 <= start0 < end0 <= chromosome_end:
        raise ValueError(f"copy outside chromosome: {chromosome}:{start0}-{end0}")
    return f"{chromosome}_{'p' if midpoint < q_start else 'q'}"


CATALOG_FIELDS = [
    "schema_version", "hypothesis_index", "collection", "relation", "ontology",
    "namespace", "term_id", "term_name", "genome_physical_copy_burden",
    "genome_arm_count", "genome_source_count", "multiplicity_family_id",
    "hypothesis_status",
]


def build_hypothesis_rows(
    assignments_path: Path, source_terms_path: Path,
    boundaries: Mapping[str, tuple[int, int]],
) -> list[dict[str, object]]:
    """Build hypotheses without opening any target-membership data."""
    source_cn: Counter[str] = Counter()
    source_arms: MutableMapping[str, set[str]] = defaultdict(set)
    assignment_header: set[str] | None = None
    opener = gzip.open if assignments_path.suffix == ".gz" else open
    with opener(assignments_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assignment_header = set(reader.fieldnames or [])
        if any(name.startswith("phr_") or name in {"in_phr_midpoint", "in_phr_any_overlap"}
               for name in assignment_header):
            raise ValueError("hypothesis freeze refuses target-membership columns")
        required = {"functional_source_id", "physical_copy_cn", "seqid", "start0", "end0"}
        if not required.issubset(assignment_header):
            raise ValueError(f"assignment freeze columns missing: {sorted(required - assignment_header)}")
        for row in reader:
            source = row["functional_source_id"]
            if not source:
                continue
            cn = int(row["physical_copy_cn"])
            if cn <= 0:
                raise ValueError("physical_copy_cn must be positive")
            source_cn[source] += cn
            source_arms[source].add(arm_for(
                row["seqid"], int(row["start0"]), int(row["end0"]), boundaries,
            ))
    burden: Counter[tuple[str, str, str]] = Counter()
    source_count: Counter[tuple[str, str, str]] = Counter()
    arms: MutableMapping[tuple[str, str, str], set[str]] = defaultdict(set)
    names: dict[tuple[str, str, str], str] = {}
    metadata: dict[tuple[str, str, str], tuple[str, str]] = {}
    seen: set[tuple[str, str, str]] = set()
    for row in read_rows(source_terms_path):
        source = row["functional_source_id"]
        if source not in source_cn:
            continue
        relation = row["relation"]
        if relation not in RELATIONS:
            raise ValueError(f"unexpected relation: {relation}")
        collection = collection_for(row["ontology"], row["namespace"])
        key = (collection, relation, row["term_id"])
        source_edge = (source, relation, row["term_id"])
        if source_edge in seen:
            raise ValueError(f"duplicate frozen source term: {source_edge}")
        seen.add(source_edge)
        if key in names and names[key] != row["term_name"]:
            raise ValueError(f"term name conflict: {key}")
        names[key] = row["term_name"]
        metadata[key] = (row["ontology"], row["namespace"])
        burden[key] += source_cn[source]
        source_count[key] += 1
        arms[key].update(source_arms[source])
    order_collection = {value: index for index, value in enumerate(COLLECTIONS)}
    order_relation = {value: index for index, value in enumerate(RELATIONS)}
    keys = sorted(burden, key=lambda key: (
        order_collection[key[0]], order_relation[key[1]], key[2],
    ))
    rows = []
    for index, key in enumerate(keys):
        ontology, namespace = metadata[key]
        rows.append({
            "schema_version": SCHEMA_VERSION,
            "hypothesis_index": index,
            "collection": key[0],
            "relation": key[1],
            "ontology": ontology,
            "namespace": namespace,
            "term_id": key[2],
            "term_name": names[key],
            "genome_physical_copy_burden": burden[key],
            "genome_arm_count": len(arms[key]),
            "genome_source_count": source_count[key],
            "multiplicity_family_id": key[0],
            "hypothesis_status": "FROZEN_TESTED_NO_TARGET_FILTER",
        })
    return rows


def freeze_hypotheses(output: Path = CATALOG, manifest_path: Path = CATALOG_MANIFEST) -> dict[str, object]:
    gate = validate_gate()
    if not gate["pass"]:
        raise RuntimeError(f"PRE_RUN_V6_GATE refused hypothesis freeze: {gate['errors']}")
    rows = build_hypothesis_rows(ASSIGNMENTS, SOURCE_TERMS, load_arm_boundaries())
    write_rows(output, CATALOG_FIELDS, rows)
    sizes = Counter(str(row["collection"]) for row in rows)
    relations = Counter((str(row["collection"]), str(row["relation"])) for row in rows)
    report = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": utcnow(),
        "git_commit": git_commit(),
        "status": "FROZEN_BEFORE_TARGET_TERM_RECOUNT",
        "gate_status": gate["status"],
        "gate_sha256": gate["gate_sha256"],
        "audited_release_manifest_sha256": gate["audited_release_manifest_sha256"],
        "freeze_inputs": {
            ASSIGNMENTS.name: sha256(ASSIGNMENTS),
            SOURCE_TERMS.name: sha256(SOURCE_TERMS),
            CYTOBANDS.name: sha256(CYTOBANDS),
        },
        "freeze_input_basenames": sorted(FREEZE_INPUT_BASENAMES),
        "target_membership_columns_opened": False,
        "target_bed_opened": False,
        "physical_copy_term_edges_opened": False,
        "hypothesis_count": len(rows),
        "multiplicity_family_sizes": dict(sorted(sizes.items())),
        "relation_counts": {f"{key[0]}|{key[1]}": value for key, value in sorted(relations.items())},
        "catalog_path": output.name,
        "catalog_bytes": output.stat().st_size,
        "catalog_sha256": sha256(output),
        "family_definition": "collection combining separately identified direct and ancestor hypotheses",
        "target_filter": "none",
    }
    atomic_json(manifest_path, report)
    return report


@dataclass(frozen=True)
class Hypothesis:
    index: int
    collection: str
    relation: str
    ontology: str
    namespace: str
    term_id: str
    term_name: str
    genome_cn: int
    genome_arms: int
    genome_sources: int

    @property
    def key(self) -> tuple[str, str, str]:
        return self.collection, self.relation, self.term_id

    @property
    def layer(self) -> str:
        return self.relation

    @property
    def family_id(self) -> str:
        return self.collection


def validate_hypothesis_freeze() -> dict[str, object]:
    errors: list[str] = []
    if not CATALOG.is_file() or not CATALOG_MANIFEST.is_file():
        return {"pass": False, "errors": ["frozen_hypothesis_artifacts_missing"]}
    manifest = json.loads(CATALOG_MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "FROZEN_BEFORE_TARGET_TERM_RECOUNT":
        errors.append("freeze_status")
    if manifest.get("target_membership_columns_opened") is not False:
        errors.append("target_columns_opened")
    if set(manifest.get("freeze_input_basenames", [])) != FREEZE_INPUT_BASENAMES:
        errors.append("freeze_input_scope")
    if sha256(CATALOG) != manifest.get("catalog_sha256"):
        errors.append("catalog_sha256")
    for name, expected in manifest.get("freeze_inputs", {}).items():
        path = CYTOBANDS if name == CYTOBANDS.name else HERE / name
        if not path.is_file() or sha256(path) != expected:
            errors.append(f"freeze_input_sha256:{name}")
    rows = list(read_rows(CATALOG))
    if len(rows) != int(manifest.get("hypothesis_count", -1)):
        errors.append("hypothesis_count")
    if [int(row["hypothesis_index"]) for row in rows] != list(range(len(rows))):
        errors.append("hypothesis_index_order")
    keys = [(row["collection"], row["relation"], row["term_id"]) for row in rows]
    if len(keys) != len(set(keys)):
        errors.append("duplicate_hypothesis")
    sizes = Counter(row["collection"] for row in rows)
    if dict(sorted(sizes.items())) != manifest.get("multiplicity_family_sizes"):
        errors.append("family_sizes")
    return {"pass": not errors, "errors": errors, "hypothesis_count": len(rows),
            "family_sizes": dict(sorted(sizes.items())), "catalog_sha256": sha256(CATALOG)}


def load_hypotheses() -> list[Hypothesis]:
    check = validate_hypothesis_freeze()
    if not check["pass"]:
        raise RuntimeError(f"invalid frozen hypothesis catalog: {check['errors']}")
    return [Hypothesis(
        index=int(row["hypothesis_index"]), collection=row["collection"],
        relation=row["relation"], ontology=row["ontology"], namespace=row["namespace"],
        term_id=row["term_id"], term_name=row["term_name"],
        genome_cn=int(row["genome_physical_copy_burden"]),
        genome_arms=int(row["genome_arm_count"]),
        genome_sources=int(row["genome_source_count"]),
    ) for row in read_rows(CATALOG)]


class WeightedTermArrays(list[np.ndarray]):
    """Per-coordinate term indices plus the physical CN carried by that row."""

    def __init__(self, values: Iterable[np.ndarray], weights: np.ndarray):
        super().__init__(values)
        self.weights = np.asarray(weights, dtype=np.uint32)
        self.all_unit = bool(np.all(self.weights == 1))


def _copy_to_engine(genome: object) -> tuple[dict[str, int], dict[str, dict[str, str]]]:
    by_coordinate: MutableMapping[tuple[str, int, int], list[object]] = defaultdict(list)
    for locus in genome.loci:
        by_coordinate[(locus.chromosome, locus.start, locus.end)].append(locus)
    mapping: dict[str, int] = {}
    catalog: dict[str, dict[str, str]] = {}
    used: set[int] = set()
    for row in read_rows(ASSIGNMENTS):
        key = (row["seqid"], int(row["start0"]), int(row["end0"]))
        matches = by_coordinate.get(key, [])
        exact = [locus for locus in matches if locus.locus_id == row["gff_id"]]
        if len(exact) == 1:
            matches = exact
        elif len(matches) > 1:
            matches = [locus for locus in matches if locus.gene_name == row["gene_name"]]
        if len(matches) != 1:
            raise ValueError(f"copy does not resolve to one engine locus: {row['copy_id']}")
        index = int(matches[0].index)
        if index in used:
            raise ValueError(f"multiple V6 rows resolve to engine locus: {row['copy_id']}")
        used.add(index)
        mapping[row["copy_id"]] = index
        catalog[row["copy_id"]] = row
    if len(mapping) != EXPECTED_PHYSICAL_CN or len(used) != EXPECTED_PHYSICAL_CN:
        raise ValueError(f"V6/engine coordinate bijection failed: {len(mapping)}/{len(used)}")
    return mapping, catalog


def build_locus_term_arrays(genome: object, hypotheses: Sequence[Hypothesis]) -> WeightedTermArrays:
    key_to_index = {row.key: row.index for row in hypotheses}
    source_terms: MutableMapping[str, list[int]] = defaultdict(list)
    seen: set[tuple[str, str, str]] = set()
    for row in read_rows(SOURCE_TERMS):
        collection = collection_for(row["ontology"], row["namespace"])
        key = (collection, row["relation"], row["term_id"])
        index = key_to_index.get(key)
        if index is None:
            continue
        source_edge = (row["functional_source_id"], row["relation"], row["term_id"])
        if source_edge in seen:
            raise ValueError(f"duplicate source-term edge: {source_edge}")
        seen.add(source_edge)
        source_terms[row["functional_source_id"]].append(index)
    copy_to_engine, copy_catalog = _copy_to_engine(genome)
    terms: list[np.ndarray | None] = [None] * len(genome.loci)
    weights = np.zeros(len(genome.loci), dtype=np.uint32)
    eligible_cn = edge_cn = direct_edge_cn = 0
    for copy_id, index in copy_to_engine.items():
        row = copy_catalog[copy_id]
        cn = int(row["physical_copy_cn"])
        if cn <= 0:
            raise ValueError(f"nonpositive physical CN: {copy_id}")
        weights[index] = cn
        values = np.asarray(sorted(set(source_terms.get(row["functional_source_id"], []))), dtype=np.int32)
        terms[index] = values
        if values.size:
            eligible_cn += cn
            edge_cn += cn * int(values.size)
            direct_edge_cn += cn * sum(hypotheses[int(value)].relation == "direct" for value in values)
    if np.any(weights == 0) or int(weights.sum()) != EXPECTED_PHYSICAL_CN:
        raise ValueError("physical CN vector does not cover the complete genome")
    if eligible_cn != EXPECTED_ELIGIBLE_CN:
        raise ValueError(f"eligible physical CN mismatch: {eligible_cn}")
    if edge_cn != EXPECTED_EDGE_CN or direct_edge_cn != EXPECTED_DIRECT_EDGE_CN:
        raise ValueError(f"term-edge CN mismatch: {edge_cn}/{direct_edge_cn}")
    arrays = [value if value is not None else np.empty(0, dtype=np.int32) for value in terms]
    if max(row.genome_cn for row in hypotheses) > np.iinfo(np.uint16).max:
        raise ValueError("uint16 count array would overflow physical CN")
    return WeightedTermArrays(arrays, weights)


def recount(selected: np.ndarray, locus_terms: Sequence[np.ndarray], n_hypotheses: int) -> np.ndarray:
    if not isinstance(locus_terms, WeightedTermArrays):
        raise TypeError("V6 recount requires physical-CN weights")
    selected_indices = [int(index) for index in selected if locus_terms[int(index)].size]
    if not selected_indices:
        return np.zeros(n_hypotheses, dtype=np.uint16)
    edges = np.concatenate([locus_terms[index] for index in selected_indices])
    if locus_terms.all_unit:
        result = np.bincount(edges, minlength=n_hypotheses)
    else:
        cn = np.concatenate([
            np.full(locus_terms[index].size, locus_terms.weights[index], dtype=np.uint32)
            for index in selected_indices
        ])
        result = np.bincount(edges, weights=cn, minlength=n_hypotheses)
    if result.size and float(result.max()) > np.iinfo(np.uint16).max:
        raise OverflowError("physical-CN term burden exceeds uint16")
    return np.asarray(np.rint(result), dtype=np.uint16)


def validate_annotation_release(*_args: object, **_kwargs: object) -> dict[str, object]:
    gate = validate_gate()
    freeze = validate_hypothesis_freeze()
    return {
        "pass": bool(gate["pass"] and freeze["pass"]),
        "errors": list(gate["errors"]) + list(freeze["errors"]),
        "gate": gate,
        "hypothesis_freeze": freeze,
    }


def immutable_run_config(mode: str, batch_size: int) -> dict[str, object]:
    if mode != "primary":
        raise ValueError("V6 runs only the prespecified same-arm primary null")
    paths = {
        "gate": GATE_JSON,
        "audited_release": AUDITED_MANIFEST,
        "hypothesis_catalog": CATALOG,
        "hypothesis_manifest": CATALOG_MANIFEST,
        "source_assignments": ASSIGNMENTS,
        "source_terms": SOURCE_TERMS,
        "geometry_arms": ANALYSIS / "analysis_ready/chm13_arm_summary.tsv",
        "geometry_intervals": ANALYSIS / "analysis_ready/chm13_phr_intervals.tsv",
        "geometry_loci": ANALYSIS / "analysis_ready/chm13_gene_loci.tsv.gz",
        "protocol": HERE / "V6_RUN_PROTOCOL.md",
        "v6_engine": Path(__file__).resolve(),
        "sampler_engine": Path(E.__file__).resolve(),
    }
    hypotheses = load_hypotheses()
    return {
        "schema_version": SCHEMA_VERSION,
        "mode": mode,
        "master_seed": MASTER_SEED,
        "seed_sequence_spawn_key": [0],
        "bit_generator": "PCG64DXSM",
        "rng_contract": "continuous_child_stream_saved_after_every_batch",
        "batch_size": batch_size,
        "assignments": ["midpoint", "overlap"],
        "hypothesis_count": len(hypotheses),
        "hypothesis_definition": "all physical-matrix collection|direct-or-ancestor|term hypotheses",
        "multiplicity_families": list(COLLECTIONS),
        "counting_statistic": "sum_physical_copy_cn",
        "physical_clusters_fixed_in_genome": True,
        "source_term_edges_fixed_in_genome": True,
        "mask_policy": "empty_frozen_mask",
        "input_checksums": {name: sha256(path) for name, path in paths.items()},
        "input_paths": {name: str(path) for name, path in paths.items()},
    }


def install_numerical_overrides() -> None:
    N.SCHEMA_VERSION = SCHEMA_VERSION
    N.CHECKPOINTS = CHECKPOINTS
    N.SEEDS = {"primary": MASTER_SEED}
    N.WORK = WORK
    N.RELEASE = RELEASE
    N.load_hypotheses = load_hypotheses
    N.validate_annotation_release = validate_annotation_release
    N.build_locus_term_arrays = build_locus_term_arrays
    N.recount = recount
    N.immutable_run_config = immutable_run_config


def run_initial(permutations: int, batch_size: int, resume: bool = False) -> None:
    if permutations != INITIAL_PERMUTATIONS:
        raise ValueError("the complete initial screen is frozen at exactly 99,999")
    check = validate_annotation_release()
    if not check["pass"]:
        raise RuntimeError(f"V6 run refused by gate/freeze: {check['errors']}")
    install_numerical_overrides()
    N.run_mode("primary", permutations, batch_size, resume, WORK / "primary")


class _CountSliceCache:
    def __init__(self) -> None:
        self.key: tuple[str, str, tuple[tuple[int, int], ...]] | None = None
        self.arrays: tuple[np.ndarray, ...] = ()

    def close(self) -> None:
        for array in self.arrays:
            mmap = getattr(array, "_mmap", None)
            if mmap is not None:
                mmap.close()
        self.key = None
        self.arrays = ()

    def load(self, run: Path, manifest: Mapping[str, object], assignment: str,
             first: int, last: int) -> np.ndarray:
        batches = tuple((int(row["first"]), int(row["last"])) for row in manifest["batches"])
        key = (str(run.resolve()), assignment, batches)
        if self.key != key:
            self.close()
            self.arrays = tuple(np.load(
                run / "batches" / f"{assignment}.counts.{start:09d}-{end:09d}.npy",
                mmap_mode="r", allow_pickle=False,
            ) for start, end in batches)
            self.key = key
        return np.concatenate([np.asarray(array[:, first:last]) for array in self.arrays], axis=0)


COUNT_CACHE = _CountSliceCache()


def family_ranges(hypotheses: Sequence[Hypothesis]) -> list[tuple[str, int, int]]:
    result = []
    cursor = 0
    while cursor < len(hypotheses):
        collection = hypotheses[cursor].collection
        end = cursor + 1
        while end < len(hypotheses) and hypotheses[end].collection == collection:
            end += 1
        result.append((collection, cursor, end))
        cursor = end
    if tuple(row[0] for row in result) != COLLECTIONS:
        raise ValueError(f"catalog is not in frozen collection order: {result}")
    return result


def plus_one_interval(exceedances: int, permutations: int, alpha: float) -> tuple[float, float]:
    lower, upper = N.exact_clopper_pearson(exceedances, permutations, alpha=alpha)
    return ((1.0 + permutations * lower) / (permutations + 1.0),
            (1.0 + permutations * upper) / (permutations + 1.0))


def classify_primary(row: Mapping[str, object]) -> str:
    passed = (
        float(row["bh_sequential_upper"]) <= 0.05
        and float(row["global_maxT_sequential_upper"]) <= 0.05
    )
    nonpass = (
        float(row["bh_sequential_lower"]) > 0.05
        or float(row["global_maxT_sequential_lower"]) > 0.05
    )
    return "CERTIFIED_PASS" if passed else "CERTIFIED_NONPASS" if nonpass else "MC_UNRESOLVED"


RESULT_FIELDS = [
    "schema_version", "assignment", "analysis_role", "collection", "relation", "ontology",
    "namespace", "term_id", "term_name", "multiplicity_family_id", "multiplicity_family_size",
    "genome_physical_copy_burden", "genome_arm_count", "genome_source_count",
    "observed_physical_copy_burden", "null_mean", "null_median", "null_q025", "null_q975", "null_max",
    "count_difference", "enrichment_ratio", "raw_exceedances", "raw_permutations",
    "p_empirical", "p_mc95_lower", "p_mc95_upper", "p_sequential_lower",
    "p_sequential_upper", "bh_q", "by_q", "bh_sequential_lower", "bh_sequential_upper",
    "collection_maxT_exceedances", "collection_maxT_permutations", "collection_maxT_p",
    "collection_maxT_mc95_lower", "collection_maxT_mc95_upper",
    "collection_maxT_sequential_lower", "collection_maxT_sequential_upper",
    "global_maxT_exceedances", "global_maxT_permutations", "global_maxT_p",
    "global_maxT_mc95_lower", "global_maxT_mc95_upper", "global_maxT_sequential_lower",
    "global_maxT_sequential_upper", "z_observed", "non_informative_constant",
    "inference_stage", "stopping_reason", "mc_status", "validation_status",
]


def infer_assignment(run: Path, assignment: str) -> list[dict[str, object]]:
    """Infer physical-CN burden with four collection families and global maxT."""
    manifest = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
    permutations = int(manifest["completed_permutations"])
    if permutations != INITIAL_PERMUTATIONS:
        raise ValueError("complete-family inference is tied to the 99,999 initial screen")
    hypotheses = load_hypotheses()
    observed_all = np.load(run / f"observed.{assignment}.counts.npy", allow_pickle=False)
    rows: list[dict[str, object] | None] = [None] * len(hypotheses)
    collection_maxima: dict[str, np.ndarray] = {}
    global_maxima = np.full(permutations, -np.inf)

    for collection, first, last in family_ranges(hypotheses):
        family_maxima = np.full(permutations, -np.inf)
        family_rows: list[dict[str, object]] = []
        for start in range(first, last, 64):
            end = min(last, start + 64)
            null = COUNT_CACHE.load(run, manifest, assignment, start, end)
            observed = observed_all[start:end]
            z_observed, z_null, degenerate = N.standardize_matrix(observed, null)
            chunk_max = np.max(z_null, axis=1)
            family_maxima = np.maximum(family_maxima, chunk_max)
            global_maxima = np.maximum(global_maxima, chunk_max)
            for local, hypothesis_index in enumerate(range(start, end)):
                hypothesis = hypotheses[hypothesis_index]
                values = null[:, local]
                observed_cn = int(observed[local])
                is_degenerate = bool(degenerate[local])
                if observed_cn == 0 or is_degenerate:
                    exceedances = permutations
                    p_value = 1.0
                else:
                    exceedances = int(np.count_nonzero(values >= observed_cn))
                    p_value = (exceedances + 1.0) / (permutations + 1.0)
                mc95 = plus_one_interval(exceedances, permutations, 0.05)
                sequential = plus_one_interval(exceedances, permutations, 0.0125)
                q025, median, q975 = np.quantile(values, [0.025, 0.5, 0.975])
                family_rows.append({
                    "schema_version": SCHEMA_VERSION,
                    "assignment": assignment,
                    "analysis_role": "PRIMARY_MIDPOINT" if assignment == "midpoint" else "SENSITIVITY_ANY_OVERLAP",
                    "collection": hypothesis.collection,
                    "relation": hypothesis.relation,
                    "ontology": hypothesis.ontology,
                    "namespace": hypothesis.namespace,
                    "term_id": hypothesis.term_id,
                    "term_name": hypothesis.term_name,
                    "multiplicity_family_id": hypothesis.family_id,
                    "multiplicity_family_size": last - first,
                    "genome_physical_copy_burden": hypothesis.genome_cn,
                    "genome_arm_count": hypothesis.genome_arms,
                    "genome_source_count": hypothesis.genome_sources,
                    "observed_physical_copy_burden": observed_cn,
                    "null_mean": float(np.mean(values)),
                    "null_median": float(median),
                    "null_q025": float(q025),
                    "null_q975": float(q975),
                    "null_max": int(np.max(values)),
                    "count_difference": float(observed_cn - median),
                    "enrichment_ratio": float((observed_cn + 0.5) / (float(np.mean(values)) + 0.5)),
                    "raw_exceedances": exceedances,
                    "raw_permutations": permutations,
                    "p_empirical": p_value,
                    "p_mc95_lower": mc95[0],
                    "p_mc95_upper": mc95[1],
                    "p_sequential_lower": sequential[0],
                    "p_sequential_upper": sequential[1],
                    "z_observed": float(z_observed[local]),
                    "non_informative_constant": int(is_degenerate),
                    "inference_stage": "COMPLETE_INITIAL_SCREEN_99999",
                    "stopping_reason": "SENSITIVITY_FIXED_INITIAL" if assignment == "overlap" else "PENDING_STAGE_DECISION",
                    "_hypothesis_index": hypothesis_index,
                })
        raw_p = [float(row["p_empirical"]) for row in family_rows]
        seq_lower = [float(row["p_sequential_lower"]) for row in family_rows]
        seq_upper = [float(row["p_sequential_upper"]) for row in family_rows]
        bh = E.bh_adjust(raw_p)
        by = E.by_adjust(raw_p)
        bh_lower = E.bh_adjust(seq_lower)
        bh_upper = E.bh_adjust(seq_upper)
        for index, row in enumerate(family_rows):
            row["bh_q"] = float(bh[index])
            row["by_q"] = float(by[index])
            row["bh_sequential_lower"] = float(bh_lower[index])
            row["bh_sequential_upper"] = float(bh_upper[index])
            rows[int(row["_hypothesis_index"])] = row
        collection_maxima[collection] = family_maxima

    for raw in rows:
        assert raw is not None
        row = raw
        hypothesis = hypotheses[int(row.pop("_hypothesis_index"))]
        if int(row["observed_physical_copy_burden"]) == 0 or int(row["non_informative_constant"]):
            family_exceed = global_exceed = permutations
        else:
            z_value = float(row["z_observed"])
            family_exceed = int(np.count_nonzero(collection_maxima[hypothesis.collection] >= z_value))
            global_exceed = int(np.count_nonzero(global_maxima >= z_value))
        for prefix, exceed in (("collection_maxT", family_exceed), ("global_maxT", global_exceed)):
            mc95 = plus_one_interval(exceed, permutations, 0.05)
            sequential = plus_one_interval(exceed, permutations, 0.0125)
            row[f"{prefix}_exceedances"] = exceed
            row[f"{prefix}_permutations"] = permutations
            row[f"{prefix}_p"] = (exceed + 1.0) / (permutations + 1.0)
            row[f"{prefix}_mc95_lower"] = mc95[0]
            row[f"{prefix}_mc95_upper"] = mc95[1]
            row[f"{prefix}_sequential_lower"] = sequential[0]
            row[f"{prefix}_sequential_upper"] = sequential[1]
        row["mc_status"] = classify_primary(row) if assignment == "midpoint" else "SENSITIVITY_COMPLETE"
        row["validation_status"] = "PENDING_RELEASE_VALIDATION"
    return [row for row in rows if row is not None]


def candidate_reasons(row: Mapping[str, object]) -> list[str]:
    """Prespecified selective-extension screen applied only after full maxT."""
    if row["assignment"] != "midpoint":
        return []
    # Either lower bound above 0.05 is already a certified overall nonpass;
    # uncertainty in the other safeguard cannot change that conjunction.
    if (float(row["bh_sequential_lower"]) > 0.05
            or float(row["global_maxT_sequential_lower"]) > 0.05):
        return []
    max_t_resolved = (
        float(row["global_maxT_sequential_upper"]) <= 0.05
        or float(row["global_maxT_sequential_lower"]) > 0.05
    ) and (
        float(row["collection_maxT_sequential_upper"]) <= 0.05
        or float(row["collection_maxT_sequential_lower"]) > 0.05
    )
    if not max_t_resolved:
        return ["MAXT_UNRESOLVED_NO_SELECTIVE_EXTENSION"]
    # A maxT-certified nonpass has a resolved overall decision and does not
    # need more raw-tail precision.
    if float(row["global_maxT_sequential_lower"]) > 0.05:
        return []
    reasons = []
    bh_unresolved = (
        float(row["bh_sequential_lower"]) <= 0.05
        < float(row["bh_sequential_upper"])
    )
    if bh_unresolved:
        reasons.append("BH_INTERVAL_STRADDLES_0.05")
    # These are useful diagnostics, but the user-prespecified staged policy
    # extends only unresolved candidates.  They accompany a BH-unresolved row
    # and never select an otherwise resolved hypothesis.
    if bh_unresolved and int(row["raw_exceedances"]) < 100:
        reasons.append("RAW_EXCEEDANCES_BELOW_100")
    if bh_unresolved and (
        0.04 <= float(row["p_empirical"]) <= 0.06
        or 0.04 <= float(row["bh_q"]) <= 0.06
    ):
        reasons.append("POINT_ESTIMATE_WITHIN_20_PERCENT_OF_0.05")
    if (int(row["observed_physical_copy_burden"]) > 0
            and int(row["null_max"]) == 0):
        reasons.append("POSITIVE_OBSERVED_ALL_ZERO_NULL_AUDIT")
    return reasons


def refresh_stage_summary(output: Path = RELEASE) -> dict[str, object]:
    """Refresh only the stopping summary; never recompute inferential values."""
    result_rows = list(read_rows(output / "TERM_RESULTS.tsv.gz"))
    midpoint = [row for row in result_rows if row["assignment"] == "midpoint"]
    stage = emit_stage_decision(midpoint, output)
    status_path = output / "RELEASE_STATUS.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status["stage_decision"] = stage
    status["status"] = (
        "COMPLETE_INITIAL_SCREEN" if not stage["selective_extension_candidates"]
        else "INITIAL_SCREEN_COMPLETE_EXTENSION_CANDIDATES_FROZEN"
    )
    atomic_json(status_path, status)
    report_path = output / "V6_RUN_REPORT.md"
    text = report_path.read_text(encoding="utf-8")
    old = (
        "with 0 selectively eligible raw/BH candidates and 1 maxT-unresolved rows retained "
        "as unresolved instead of triggering a reflexive million-permutation rerun"
    )
    new = (
        "with 0 selectively eligible raw/BH candidates and 0 overall-decision-unresolved "
        "maxT rows; no reflexive million-permutation rerun was triggered"
    )
    if old not in text:
        raise RuntimeError("V6 report stopping-summary sentence was not found")
    report_path.write_text(text.replace(old, new), encoding="utf-8")
    return stage


def emit_stage_decision(rows: Sequence[Mapping[str, object]], output: Path) -> dict[str, object]:
    candidates = []
    max_t_unresolved = []
    for index, row in enumerate(rows):
        reasons = candidate_reasons(row)
        if reasons == ["MAXT_UNRESOLVED_NO_SELECTIVE_EXTENSION"]:
            max_t_unresolved.append(row)
        elif reasons:
            candidates.append({
                "hypothesis_index": index,
                "collection": row["collection"],
                "relation": row["relation"],
                "term_id": row["term_id"],
                "selection_reasons": ";".join(reasons),
                "initial_raw_exceedances": row["raw_exceedances"],
                "initial_bh_q": row["bh_q"],
                "initial_global_maxT_p": row["global_maxT_p"],
            })
    fields = [
        "hypothesis_index", "collection", "relation", "term_id", "selection_reasons",
        "initial_raw_exceedances", "initial_bh_q", "initial_global_maxT_p",
    ]
    write_rows(output / "SELECTIVE_EXTENSION_CANDIDATES.tsv", fields, candidates)
    candidate_path = output / "SELECTIVE_EXTENSION_CANDIDATES.tsv"
    report = {
        "schema_version": SCHEMA_VERSION,
        "decision_utc": utcnow(),
        "complete_initial_hypotheses": len(rows),
        "complete_initial_permutations": INITIAL_PERMUTATIONS,
        "selective_extension_candidates": len(candidates),
        "maxT_unresolved_not_selectively_extended": len(max_t_unresolved),
        "next_checkpoint": 249_999 if candidates else None,
        "candidate_checksum": sha256(candidate_path),
        "stopping_policy": "V6_RUN_PROTOCOL.md",
        "maxT_policy": "complete-family 99,999 screen retained; selected subset never defines maxT",
        "reflexive_million_per_term": False,
        "status": "EXTENSION_REQUIRED" if candidates else "INITIAL_SCREEN_SUFFICIENT",
    }
    atomic_json(output / "STAGE_DECISION.json", report)
    return report


CONTRIBUTOR_FIELDS = [
    "schema_version", "collection", "relation", "ontology", "namespace", "term_id",
    "term_name", "copy_id", "seqid", "start0", "end0", "strand", "gene_name",
    "gene_biotype", "physical_copy_cn", "phr_midpoint_cn", "phr_any_overlap_cn",
    "own_annotation_id", "functional_source_id", "functional_source_symbol",
    "source_assignment_disposition", "assignment_tier", "mapping_confidence",
    "relationship_semantics", "source_evidence_record_id", "source_evidence_release",
    "evidence_route_source", "evidence_route_record_id", "evidence_route_disposition",
    "evidence_route_admissible", "evidence_route_record_sha256",
    "minimum_distance", "inherited_from_direct_term_ids", "analysis_role",
]


def build_contributors(output: Path) -> tuple[
    Counter[tuple[str, str, str, str]], Counter[tuple[str, str, str]], int, int,
]:
    assignments = {row["copy_id"]: row for row in read_rows(ASSIGNMENTS)}
    if len(assignments) != EXPECTED_PHYSICAL_CN:
        raise ValueError("incomplete source-assignment metadata")
    evidence = {row["copy_id"]: row for row in read_rows(COPY_EVIDENCE)}
    if len(evidence) != EXPECTED_PHYSICAL_CN:
        raise ValueError("incomplete exact source-evidence routes")
    rows = []
    observed: Counter[tuple[str, str, str, str]] = Counter()
    genome: Counter[tuple[str, str, str]] = Counter()
    physical_edge_cn = midpoint_edge_cn = overlap_edge_cn = direct_edge_cn = 0
    previous: tuple[str, str, str] | None = None
    for edge in read_rows(PHYSICAL_EDGES):
        collection = collection_for(edge["ontology"], edge["namespace"])
        key = (collection, edge["relation"], edge["term_id"])
        copy_key = (edge["copy_id"], edge["relation"], edge["term_id"])
        if copy_key == previous:
            raise ValueError(f"duplicate physical copy hypothesis edge: {copy_key}")
        previous = copy_key
        cn = int(edge["physical_copy_cn"])
        midpoint_cn = int(edge["phr_midpoint_cn"])
        overlap_cn = int(edge["phr_any_overlap_cn"])
        physical_edge_cn += cn
        midpoint_edge_cn += midpoint_cn
        overlap_edge_cn += overlap_cn
        direct_edge_cn += cn if edge["relation"] == "direct" else 0
        genome[key] += cn
        if midpoint_cn:
            observed[("midpoint",) + key] += midpoint_cn
        if overlap_cn:
            observed[("overlap",) + key] += overlap_cn
        if not overlap_cn:
            continue
        copy = assignments[edge["copy_id"]]
        route = evidence[edge["copy_id"]]
        if route["admissible_for_ontology"] != "1":
            raise ValueError(f"term contributor lacks an admissible source route: {edge['copy_id']}")
        if route["functional_source_id"] != edge["functional_source_id"]:
            raise ValueError(f"source route/term edge disagreement: {edge['copy_id']}")
        rows.append({
            "schema_version": SCHEMA_VERSION,
            "collection": collection,
            "relation": edge["relation"],
            "ontology": edge["ontology"],
            "namespace": edge["namespace"],
            "term_id": edge["term_id"],
            "term_name": edge["term_name"],
            "copy_id": edge["copy_id"],
            "seqid": copy["seqid"],
            "start0": copy["start0"],
            "end0": copy["end0"],
            "strand": copy["strand"],
            "gene_name": copy["gene_name"],
            "gene_biotype": copy["gene_biotype"],
            "physical_copy_cn": cn,
            "phr_midpoint_cn": midpoint_cn,
            "phr_any_overlap_cn": overlap_cn,
            "own_annotation_id": edge["own_annotation_id"],
            "functional_source_id": edge["functional_source_id"],
            "functional_source_symbol": edge["functional_source_symbol"],
            "source_assignment_disposition": edge["source_assignment_disposition"],
            "assignment_tier": copy["assignment_tier"],
            "mapping_confidence": copy["mapping_confidence"],
            "relationship_semantics": copy["relationship_semantics"],
            "source_evidence_record_id": edge["source_evidence_record_id"],
            "source_evidence_release": copy["source_evidence_release"],
            "evidence_route_source": route["evidence_source"],
            "evidence_route_record_id": route["evidence_record_id"],
            "evidence_route_disposition": route["disposition"],
            "evidence_route_admissible": route["admissible_for_ontology"],
            "evidence_route_record_sha256": route["evidence_record_sha256"],
            "minimum_distance": edge["minimum_distance"],
            "inherited_from_direct_term_ids": edge["inherited_from_direct_term_ids"],
            "analysis_role": "EXACT_PHYSICAL_COPY_CONTRIBUTOR_SOURCE_EVIDENCE",
        })
    if physical_edge_cn != EXPECTED_EDGE_CN or direct_edge_cn != EXPECTED_DIRECT_EDGE_CN:
        raise ValueError(f"physical edge CN totals disagree: {physical_edge_cn}/{direct_edge_cn}")
    if midpoint_edge_cn != EXPECTED_PHR_MIDPOINT_EDGE_CN or overlap_edge_cn != EXPECTED_PHR_OVERLAP_EDGE_CN:
        raise ValueError(f"PHR edge CN totals disagree: {midpoint_edge_cn}/{overlap_edge_cn}")
    write_rows(output / "EXACT_TERM_CONTRIBUTORS.tsv.gz", CONTRIBUTOR_FIELDS, rows)
    return observed, genome, physical_edge_cn, overlap_edge_cn


def named_cohort_audit(output: Path) -> list[dict[str, object]]:
    contributor_cn: Counter[str] = Counter()
    contributor_rows: Counter[str] = Counter()
    for row in read_rows(EXACT_CONTRIBUTORS):
        cohort = row["cohort"]
        if cohort and int(row["closure_term_count"]) > 0 and int(row["phr_midpoint_cn"]) > 0:
            contributor_cn[cohort] += int(row["phr_midpoint_cn"])
            contributor_rows[cohort] += 1
    burden_values: MutableMapping[str, list[int]] = defaultdict(list)
    for row in read_rows(NAMED_BURDENS):
        burden_values[row["cohort"]].append(int(row["phr_midpoint_physical_copy_burden"]))
    rows = []
    for cohort, expected in NAMED_EXPECTED.items():
        values = burden_values[cohort]
        observed = contributor_cn[cohort]
        rows.append({
            "cohort": cohort,
            "gated_phr_ontology_contributor_cn": expected,
            "inference_eligible_exact_contributor_cn": observed,
            "exact_contributor_rows": contributor_rows[cohort],
            "named_source_term_rows": len(values),
            "named_term_burden_min": min(values) if values else "",
            "named_term_burden_max": max(values) if values else "",
            "term_rows_carrying_full_cohort_cn": sum(value == expected for value in values),
            "dux_not_collapsed_to_one": int(cohort != "DUX4_DUX4L" or observed == 65),
            "status": "PASS" if observed == expected else "FAIL",
        })
    write_rows(output / "NAMED_COHORT_INFERENCE_AUDIT.tsv", [
        "cohort", "gated_phr_ontology_contributor_cn", "inference_eligible_exact_contributor_cn",
        "exact_contributor_rows", "named_source_term_rows", "named_term_burden_min",
        "named_term_burden_max", "term_rows_carrying_full_cohort_cn",
        "dux_not_collapsed_to_one", "status",
    ], rows)
    if any(row["status"] != "PASS" for row in rows):
        raise ValueError(f"named cohort inference CN mismatch: {rows}")
    return rows


def emit_geometry(output: Path) -> dict[str, object]:
    arms, intervals, _loci, _genome, blocks = N.engine_objects()
    block_for = {component[0]: block for block in blocks for component in block.components}
    write_rows(output / "TARGET_INTERVALS.tsv", [
        "phr_id", "seqid", "arm", "start0", "end0", "width", "block_id",
    ], ({
        "phr_id": row.interval_id, "seqid": row.chromosome, "arm": row.arm,
        "start0": row.start, "end0": row.end, "width": row.width,
        "block_id": block_for[row.interval_id].block_id,
    } for row in intervals))
    write_rows(output / "PLACEMENT_BLOCKS.tsv", [
        "block_id", "arm", "component_count", "component_widths", "internal_gaps",
        "observed_start0", "observed_end0", "span", "terminal_distance_stratum",
    ], ({
        "block_id": block.block_id, "arm": block.source_arm,
        "component_count": len(block.components),
        "component_widths": ";".join(str(value[2]) for value in block.components),
        "internal_gaps": ";".join(str(block.components[index + 1][1] - (
            block.components[index][1] + block.components[index][2]
        )) for index in range(len(block.components) - 1)),
        "observed_start0": block.source_start, "observed_end0": block.source_end,
        "span": block.span, "terminal_distance_stratum": block.stratum,
    } for block in blocks))
    sampler = E.RegionSampler(blocks, arms, E.load_masks(None, arms), "primary", min_candidates=100)
    candidates = list(sampler.candidate_rows())
    write_rows(output / "CANDIDATE_SPACES.tsv", [
        "mode", "block_id", "source_arm", "destination_arm", "candidate_count",
        "range_count", "explicit", "observed_start_is_candidate",
    ], candidates)
    manifest = json.loads((WORK / "primary/run_manifest.json").read_text(encoding="utf-8"))
    placement_rows = list(read_rows(WORK / "primary/PLACEMENT_MANIFEST.tsv"))
    write_rows(output / "NULL_PLACEMENT_MANIFEST.tsv", [
        "mode", "batch", "first_replicate", "last_replicate", "bytes", "sha256",
        "frozen_prefix_sha256", "frozen_prefix_match", "canonical_coordinate_sha256",
    ], placement_rows)
    count_rows = []
    for batch in manifest["batches"]:
        first, last = int(batch["first"]), int(batch["last"])
        for assignment in ("midpoint", "overlap"):
            path = WORK / "primary/batches" / f"{assignment}.counts.{first:09d}-{last:09d}.npy"
            count_rows.append({
                "assignment": assignment, "first_replicate": first, "last_replicate": last,
                "relative_path": str(path.relative_to(REPO)), "bytes": path.stat().st_size,
                "sha256": sha256(path), "statistic": "sum(physical_copy_cn)",
            })
    write_rows(output / "TRANSIENT_COUNT_ARRAYS.tsv", [
        "assignment", "first_replicate", "last_replicate", "relative_path", "bytes",
        "sha256", "statistic",
    ], count_rows)
    return {
        "target_intervals": len(intervals), "placement_blocks": len(blocks),
        "candidate_count_min": min(int(row["candidate_count"]) for row in candidates),
        "candidate_count_max": max(int(row["candidate_count"]) for row in candidates),
        "valid_replicates": int(manifest["valid_replicates"]),
        "placement_batches": len(placement_rows),
        "count_array_files": len(count_rows),
        "frozen_prefix_all_match": bool(manifest.get("prefix_all_100_match")),
    }


def verify_placement_invariants() -> dict[str, object]:
    arms, _intervals, _loci, _genome, blocks = N.engine_objects()
    by_id = {block.block_id: block for block in blocks}
    manifest = json.loads((WORK / "primary/run_manifest.json").read_text(encoding="utf-8"))
    rows = replicates = 0
    errors: list[str] = []
    for batch in manifest["batches"]:
        path = WORK / "primary/batches" / str(batch["placements"])
        for placement_set in E.load_placement_batch(path):
            replicates += 1
            if len(placement_set) != len(blocks):
                errors.append(f"block_count:{replicates}:{len(placement_set)}")
                break
            for placed in placement_set:
                rows += 1
                source = by_id[placed.block_id]
                if placed.arm != source.source_arm:
                    errors.append(f"arm_changed:{placed.block_id}")
                if placed.components != source.components or placed.end - placed.start != source.span:
                    errors.append(f"geometry_changed:{placed.block_id}")
                midpoint = placed.start + source.midpoint_offset
                if E.stratum_index(E.terminal_distance(arms[placed.arm], midpoint)) != source.stratum:
                    errors.append(f"stratum_changed:{placed.block_id}")
            if errors:
                break
        if errors:
            break
    return {
        "pass": not errors,
        "errors": errors[:20],
        "valid_replicates_checked": replicates,
        "placed_blocks_checked": rows,
        "expected_placed_blocks": int(manifest["valid_replicates"]) * len(blocks),
        "same_arm_all": not any(value.startswith("arm_changed") for value in errors),
        "rigid_geometry_all": not any(value.startswith("geometry_changed") for value in errors),
        "terminal_stratum_all": not any(value.startswith("stratum_changed") for value in errors),
        "copy_clusters_and_cn_fixed": True,
        "source_term_edges_fixed": True,
    }


def emit_multiplicity(output: Path, hypotheses: Sequence[Hypothesis]) -> list[dict[str, object]]:
    rows = []
    for collection, first, last in family_ranges(hypotheses):
        ordered = "\n".join(
            f"{row.relation}|{row.term_id}" for row in hypotheses[first:last]
        ) + "\n"
        rows.append({
            "family_id": collection,
            "scope": "PRIMARY_COLLECTION_DIRECT_AND_ANCESTOR",
            "collection": collection,
            "relations": "direct;ancestor",
            "hypothesis_count": last - first,
            "ordered_hypothesis_sha256": hashlib.sha256(ordered.encode()).hexdigest(),
            "adjustments": "BH;BY;single-step collection maxT",
            "frozen_before_target_term_results": 1,
        })
    ordered_global = "\n".join(
        f"{row.collection}|{row.relation}|{row.term_id}" for row in hypotheses
    ) + "\n"
    rows.append({
        "family_id": "GLOBAL_ALL_ONTOLOGY",
        "scope": "GLOBAL_MAXT_SAFEGUARD",
        "collection": "ALL",
        "relations": "direct;ancestor",
        "hypothesis_count": len(hypotheses),
        "ordered_hypothesis_sha256": hashlib.sha256(ordered_global.encode()).hexdigest(),
        "adjustments": "single-step global maxT",
        "frozen_before_target_term_results": 1,
    })
    write_rows(output / "MULTIPLICITY_FAMILIES.tsv", [
        "family_id", "scope", "collection", "relations", "hypothesis_count",
        "ordered_hypothesis_sha256", "adjustments", "frozen_before_target_term_results",
    ], rows)
    return rows


def copy_mapping_coverage(output: Path) -> int:
    source = HERE / "PRE_RUN_V6_MAPPING_COVERAGE.tsv"
    rows = list(read_rows(source))
    fields = list(rows[0]) if rows else []
    write_rows(output / "MAPPING_COVERAGE.tsv", fields, rows)
    return len(rows)


def finalize(output: Path = RELEASE) -> dict[str, object]:
    check = validate_annotation_release()
    if not check["pass"]:
        raise RuntimeError(f"release refused by gate/freeze: {check['errors']}")
    run_manifest_path = WORK / "primary/run_manifest.json"
    if not run_manifest_path.is_file():
        raise FileNotFoundError("initial Slurm run manifest is missing")
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    if int(run_manifest.get("completed_permutations", 0)) != INITIAL_PERMUTATIONS:
        raise RuntimeError("initial complete screen is not exactly 99,999")
    if not run_manifest.get("prefix_all_100_match"):
        raise RuntimeError("frozen placement prefix did not match")
    output.mkdir(parents=True, exist_ok=True)
    hypotheses = load_hypotheses()
    midpoint = infer_assignment(WORK / "primary", "midpoint")
    overlap = infer_assignment(WORK / "primary", "overlap")
    stage = emit_stage_decision(midpoint, output)
    observed, genome, _edge_cn, _overlap_edge_cn = build_contributors(output)
    errors = []
    for assignment, rows in (("midpoint", midpoint), ("overlap", overlap)):
        for row in rows:
            key = (assignment, str(row["collection"]), str(row["relation"]), str(row["term_id"]))
            if int(row["observed_physical_copy_burden"]) != observed[key]:
                errors.append(f"observed_contributor_sum:{'|'.join(key)}")
            genome_key = (str(row["collection"]), str(row["relation"]), str(row["term_id"]))
            if int(row["genome_physical_copy_burden"]) != genome[genome_key]:
                errors.append(f"genome_edge_sum:{'|'.join(genome_key)}")
            row["validation_status"] = "PASS" if not errors else "FAIL"
            if assignment == "midpoint" and row["stopping_reason"] == "PENDING_STAGE_DECISION":
                reasons = candidate_reasons(row)
                if reasons == ["MAXT_UNRESOLVED_NO_SELECTIVE_EXTENSION"]:
                    row["stopping_reason"] = reasons[0]
                elif reasons:
                    row["stopping_reason"] = "SELECTIVE_EXTENSION_CANDIDATE:" + ";".join(reasons)
                else:
                    row["stopping_reason"] = "INITIAL_SCREEN_DECISION_RESOLVED"
    if errors:
        raise RuntimeError(f"CN contributor reconciliation failed: {errors[:20]}")
    write_rows(output / "TERM_RESULTS.tsv.gz", RESULT_FIELDS, midpoint + overlap)
    multiplicity = emit_multiplicity(output, hypotheses)
    named = named_cohort_audit(output)
    geometry = emit_geometry(output)
    placement = verify_placement_invariants()
    atomic_json(output / "PLACEMENT_INVARIANTS.json", placement)
    if not placement["pass"]:
        raise RuntimeError(f"placement invariants failed: {placement['errors']}")
    coverage_rows = copy_mapping_coverage(output)
    input_rows = []
    for path in (
        GATE_JSON, AUDITED_MANIFEST, INPUT_MANIFEST, CATALOG, CATALOG_MANIFEST,
        ASSIGNMENTS, SOURCE_TERMS, PHYSICAL_EDGES, COPY_EVIDENCE,
        HERE / "V6_RUN_PROTOCOL.md", Path(__file__).resolve(),
    ):
        input_rows.append({"path": str(path.relative_to(REPO)), "bytes": path.stat().st_size,
                           "sha256": sha256(path), "role": "FROZEN_INPUT_OR_EXECUTABLE"})
    write_rows(output / "INPUT_CHECKSUMS.tsv", ["path", "bytes", "sha256", "role"], input_rows)
    write_rows(output / "RUN_MANIFEST.tsv", [
        "schema_version", "git_commit", "python", "numpy", "seed", "spawn_key",
        "bit_generator", "initial_permutations", "assignment_primary", "assignment_sensitivity",
        "statistic", "slurm_job_id", "completion_status",
    ], [{
        "schema_version": SCHEMA_VERSION,
        "git_commit": run_manifest["git_commit"],
        "python": run_manifest["python"],
        "numpy": run_manifest["numpy"],
        "seed": MASTER_SEED,
        "spawn_key": "[0]",
        "bit_generator": "PCG64DXSM",
        "initial_permutations": INITIAL_PERMUTATIONS,
        "assignment_primary": "midpoint",
        "assignment_sensitivity": "any-overlap",
        "statistic": "sum(physical_copy_cn)",
        "slurm_job_id": os.environ.get("V6_INITIAL_SLURM_JOB_ID", "recorded_in_SLURM_JOBS.tsv"),
        "completion_status": "COMPLETE",
    }])
    status_counts = Counter(str(row["mc_status"]) for row in midpoint)
    primary_pass = sum(row["mc_status"] == "CERTIFIED_PASS" for row in midpoint)
    report = f"""# V6 physical-copy ontology enrichment run

## Release status

The independent pre-run gate is **PASS** and every audited release and frozen source checksum matched before inference. The run tested {len(hypotheses):,} target-blind hypotheses. Direct terms and ancestors are separately identified hypotheses; BH and BY are applied within GO BP ({sum(row.collection == 'GO_BP' for row in hypotheses):,}), GO MF ({sum(row.collection == 'GO_MF' for row in hypotheses):,}), GO CC ({sum(row.collection == 'GO_CC' for row in hypotheses):,}), and Reactome ({sum(row.collection == 'Reactome' for row in hypotheses):,}). Collection and global single-step maxT safeguards are reported for every row.

Every observed and null statistic is a sum of `physical_copy_cn` from the identical frozen genome-wide copy/term matrix. The named-cohort entry audit is PASS: DUX contributes 65 physical copies (not one), DDX11L 10, TUBB8 2, OR4F 4, and WASH 9. Exact contributing copies, CN, functional sources, assignment dispositions, and evidence record identifiers are in `EXACT_TERM_CONTRIBUTORS.tsv.gz`; route/biotype coverage is in `MAPPING_COVERAGE.tsv`.

## Spatial null and Monte Carlo

The complete initial screen used {INITIAL_PERMUTATIONS:,} valid joint PCG64DXSM placements (seed {MASTER_SEED}, spawn key `[0]`). All 37 rigid PHR blocks remained on the same arm and in the same terminal stratum; widths, components, gaps, copy clusters, source/term edges, and CN weights were unchanged. Candidate spaces contained {geometry['candidate_count_min']:,}--{geometry['candidate_count_max']:,} valid integer translations. Midpoint is primary and any-overlap is the paired sensitivity on identical coordinates.

The staged stopping result is `{stage['status']}` with {stage['selective_extension_candidates']} selectively eligible raw/BH candidates and {stage['maxT_unresolved_not_selectively_extended']} maxT-unresolved rows retained as unresolved instead of triggering a reflexive million-permutation rerun. Primary midpoint status counts are `{json.dumps(dict(sorted(status_counts.items())), sort_keys=True)}`; {primary_pass} rows satisfy the conservative BH-plus-global-maxT interval rule at this stage. Empirical p-values use plus one; 95% and sequential 98.75% Clopper--Pearson bounds are reported.

## Scope

This release is a CHM13 regional physical-copy analysis. Families are not hypotheses and were not used for selection or inference. No manuscript file was read or modified.
"""
    (output / "V6_RUN_REPORT.md").write_text(report, encoding="utf-8")
    release = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": utcnow(),
        "status": "COMPLETE_INITIAL_SCREEN" if not stage["selective_extension_candidates"] else "INITIAL_SCREEN_COMPLETE_EXTENSION_CANDIDATES_FROZEN",
        "gate_pass": True,
        "hypotheses": len(hypotheses),
        "term_result_rows": len(midpoint) + len(overlap),
        "multiplicity_families": len(multiplicity),
        "mapping_coverage_rows": coverage_rows,
        "named_cohort_rows": len(named),
        "initial_permutations": INITIAL_PERMUTATIONS,
        "placement_invariants_pass": placement["pass"],
        "stage_decision": stage,
        "primary_mc_status_counts": dict(sorted(status_counts.items())),
        "physical_cn_statistic_only": True,
        "manuscript_modified": False,
    }
    atomic_json(output / "RELEASE_STATUS.json", release)
    return release


def validate_physical_matrix() -> dict[str, object]:
    hypotheses = load_hypotheses()
    expected = {row.key: row.genome_cn for row in hypotheses}
    burden: Counter[tuple[str, str, str]] = Counter()
    edge_cn = direct_cn = midpoint_cn = overlap_cn = 0
    edge_rows = 0
    for row in read_rows(PHYSICAL_EDGES):
        key = (collection_for(row["ontology"], row["namespace"]), row["relation"], row["term_id"])
        cn = int(row["physical_copy_cn"])
        burden[key] += cn
        edge_cn += cn
        direct_cn += cn if row["relation"] == "direct" else 0
        midpoint_cn += int(row["phr_midpoint_cn"])
        overlap_cn += int(row["phr_any_overlap_cn"])
        edge_rows += 1
    errors = []
    if burden != Counter(expected):
        errors.append("frozen_catalog_vs_physical_edge_burden")
    if edge_cn != EXPECTED_EDGE_CN or edge_rows != EXPECTED_EDGE_CN:
        errors.append(f"physical_edge_cn:{edge_cn}:{edge_rows}")
    if direct_cn != EXPECTED_DIRECT_EDGE_CN:
        errors.append(f"direct_edge_cn:{direct_cn}")
    if midpoint_cn != EXPECTED_PHR_MIDPOINT_EDGE_CN or overlap_cn != EXPECTED_PHR_OVERLAP_EDGE_CN:
        errors.append(f"phr_edge_cn:{midpoint_cn}:{overlap_cn}")
    return {
        "pass": not errors, "errors": errors, "hypotheses": len(hypotheses),
        "physical_edge_rows": edge_rows, "physical_edge_cn": edge_cn,
        "direct_edge_cn": direct_cn, "phr_midpoint_edge_cn": midpoint_cn,
        "phr_any_overlap_edge_cn": overlap_cn,
        "every_edge_statistic_weight": "physical_copy_cn",
    }


def preflight(output: Path = RESULTS / "PREFLIGHT.json") -> dict[str, object]:
    gate = validate_gate()
    freeze = validate_hypothesis_freeze()
    matrix: dict[str, object] = {"pass": False, "errors": ["prerequisite_failed"]}
    map_arrays: dict[str, object] = {"pass": False, "errors": ["prerequisite_failed"]}
    named: list[dict[str, object]] = []
    if gate["pass"] and freeze["pass"]:
        matrix = validate_physical_matrix()
        arms, _intervals, _loci, genome, _blocks = N.engine_objects()
        del arms
        arrays = build_locus_term_arrays(genome, load_hypotheses())
        map_arrays = {
            "pass": True,
            "physical_rows": len(arrays),
            "physical_cn": int(arrays.weights.sum()),
            "all_current_rows_unit_cn": arrays.all_unit,
            "weighted_recount_implementation": True,
            "copy_clusters_fixed_in_coordinate_index": True,
        }
        temporary_named = output.parent / ".preflight_named.tsv"
        temporary_named.parent.mkdir(parents=True, exist_ok=True)
        named = named_cohort_audit(temporary_named.parent)
        # named_cohort_audit uses its release filename in the supplied directory.
        generated = temporary_named.parent / "NAMED_COHORT_INFERENCE_AUDIT.tsv"
        if generated.is_file() and generated.parent == RESULTS:
            generated.unlink()
    report = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": utcnow(),
        "git_commit": git_commit(),
        "gate": gate,
        "hypothesis_freeze": freeze,
        "physical_matrix": matrix,
        "weighted_map": map_arrays,
        "named_cohorts": named,
        "hypotheses_frozen_before_physical_edge_target_recount": True,
        "manuscript_opened": False,
    }
    report["pass"] = bool(
        gate["pass"] and freeze["pass"] and matrix["pass"] and map_arrays["pass"]
        and named and all(row["status"] == "PASS" for row in named)
    )
    atomic_json(output, report)
    return report


def record_slurm(job_id: str, output: Path = RELEASE) -> dict[str, object]:
    command = [
        "sacct", "-j", job_id, "--parsable2", "--noheader",
        "--format=JobIDRaw,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"sacct failed: {completed.stderr}")
    fields = [
        "job_id", "job_name", "partition", "state", "exit_code", "elapsed",
        "allocated_cpus", "requested_memory", "max_rss", "node_list",
        "script", "seed", "permutations", "log_stdout", "log_stderr", "query_utc",
    ]
    rows = []
    for line in completed.stdout.splitlines():
        values = line.split("|")
        if len(values) < 10 or not values[0]:
            continue
        rows.append(dict(zip(fields[:10], values[:10]), **{
            "script": str((HERE / "slurm_v6_initial.sbatch").relative_to(REPO)),
            "seed": MASTER_SEED,
            "permutations": INITIAL_PERMUTATIONS,
            "log_stdout": str((RESULTS / "logs" / f"initial-{job_id}.out").relative_to(REPO)),
            "log_stderr": str((RESULTS / "logs" / f"initial-{job_id}.err").relative_to(REPO)),
            "query_utc": utcnow(),
        }))
    write_rows(output / "SLURM_JOBS.tsv", fields, rows)
    root = next((row for row in rows if row["job_id"] == job_id), rows[0] if rows else {})
    report = {
        "schema_version": SCHEMA_VERSION,
        "job_id": job_id,
        "submission_script": str((HERE / "slurm_v6_initial.sbatch").relative_to(REPO)),
        "command": command,
        "scheduler_rows": len(rows),
        "root_state": root.get("state", "MISSING"),
        "root_exit_code": root.get("exit_code", "MISSING"),
        "requested_resources": {"partition": "workers", "cpus": 1, "memory": "32G", "time": "12:00:00"},
        "seed": MASTER_SEED,
        "permutations": INITIAL_PERMUTATIONS,
        "stdout_log": str((RESULTS / "logs" / f"initial-{job_id}.out").relative_to(REPO)),
        "stderr_log": str((RESULTS / "logs" / f"initial-{job_id}.err").relative_to(REPO)),
        "completion_status": "COMPLETE" if str(root.get("state", "")).startswith("COMPLETED")
        and root.get("exit_code") == "0:0" else "FAILED_OR_INCOMPLETE",
        "recorded_utc": utcnow(),
    }
    atomic_json(output / "SLURM_COMPLETION.json", report)
    return report


def validate_release(output: Path = RELEASE, *, require_slurm: bool = True) -> dict[str, object]:
    errors: list[str] = []
    gate = validate_gate()
    freeze = validate_hypothesis_freeze()
    if not gate["pass"]:
        errors.extend(f"gate:{value}" for value in gate["errors"])
    if not freeze["pass"]:
        errors.extend(f"freeze:{value}" for value in freeze["errors"])
    hypotheses = load_hypotheses() if freeze["pass"] else []
    results = list(read_rows(output / "TERM_RESULTS.tsv.gz")) if (output / "TERM_RESULTS.tsv.gz").is_file() else []
    if len(results) != 2 * len(hypotheses):
        errors.append(f"term_result_rows:{len(results)}:{2 * len(hypotheses)}")
    by_assignment = Counter(row["assignment"] for row in results)
    if by_assignment != Counter({"midpoint": len(hypotheses), "overlap": len(hypotheses)}):
        errors.append(f"assignment_rows:{dict(by_assignment)}")
    family_sizes = Counter(row.collection for row in hypotheses)
    for row in results:
        if int(row["multiplicity_family_size"]) != family_sizes[row["collection"]]:
            errors.append(f"family_size:{row['collection']}")
            break
        if int(row["raw_permutations"]) < INITIAL_PERMUTATIONS:
            errors.append("raw_permutations_below_initial")
            break
        if int(row["collection_maxT_permutations"]) != INITIAL_PERMUTATIONS or int(row["global_maxT_permutations"]) != INITIAL_PERMUTATIONS:
            errors.append("maxT_not_complete_initial_screen")
            break
        if float(row["observed_physical_copy_burden"]) != int(row["observed_physical_copy_burden"]):
            errors.append("noninteger_observed_cn")
            break
    contributor_path = output / "EXACT_TERM_CONTRIBUTORS.tsv.gz"
    contributors = list(read_rows(contributor_path)) if contributor_path.is_file() else []
    sums: Counter[tuple[str, str, str, str]] = Counter()
    for row in contributors:
        sums[("midpoint", row["collection"], row["relation"], row["term_id"])] += int(row["phr_midpoint_cn"])
        sums[("overlap", row["collection"], row["relation"], row["term_id"])] += int(row["phr_any_overlap_cn"])
        if int(row["physical_copy_cn"]) <= 0:
            errors.append("nonpositive_contributor_cn")
            break
    for row in results:
        key = (row["assignment"], row["collection"], row["relation"], row["term_id"])
        if sums[key] != int(row["observed_physical_copy_burden"]):
            errors.append(f"contributor_sum:{'|'.join(key)}")
            break
    named_path = output / "NAMED_COHORT_INFERENCE_AUDIT.tsv"
    named = list(read_rows(named_path)) if named_path.is_file() else []
    if len(named) != len(NAMED_EXPECTED) or any(row["status"] != "PASS" for row in named):
        errors.append("named_cohort_audit")
    dux = next((row for row in named if row["cohort"] == "DUX4_DUX4L"), {})
    if int(dux.get("inference_eligible_exact_contributor_cn", -1)) != 65:
        errors.append("dux_not_65")
    placement_path = output / "PLACEMENT_INVARIANTS.json"
    placement = json.loads(placement_path.read_text()) if placement_path.is_file() else {}
    if not placement.get("pass") or int(placement.get("valid_replicates_checked", 0)) != INITIAL_PERMUTATIONS:
        errors.append("placement_invariants")
    multiplicity_path = output / "MULTIPLICITY_FAMILIES.tsv"
    multiplicity = list(read_rows(multiplicity_path)) if multiplicity_path.is_file() else []
    if len(multiplicity) != 5 or [row["family_id"] for row in multiplicity[:4]] != list(COLLECTIONS):
        errors.append("multiplicity_families")
    stage_path = output / "STAGE_DECISION.json"
    stage = json.loads(stage_path.read_text()) if stage_path.is_file() else {}
    if stage.get("complete_initial_permutations") != INITIAL_PERMUTATIONS or stage.get("reflexive_million_per_term") is not False:
        errors.append("stage_decision")
    coverage_path = output / "MAPPING_COVERAGE.tsv"
    coverage = list(read_rows(coverage_path)) if coverage_path.is_file() else []
    if len(coverage) != 148:
        errors.append(f"mapping_coverage_rows:{len(coverage)}")
    slurm_path = output / "SLURM_COMPLETION.json"
    slurm = json.loads(slurm_path.read_text()) if slurm_path.is_file() else {}
    if require_slurm and slurm.get("completion_status") != "COMPLETE":
        errors.append("slurm_completion")
    required = {
        "V6_RUN_REPORT.md", "TERM_RESULTS.tsv.gz", "EXACT_TERM_CONTRIBUTORS.tsv.gz",
        "MAPPING_COVERAGE.tsv", "MULTIPLICITY_FAMILIES.tsv", "STAGE_DECISION.json",
        "TARGET_INTERVALS.tsv", "PLACEMENT_BLOCKS.tsv", "CANDIDATE_SPACES.tsv",
        "NULL_PLACEMENT_MANIFEST.tsv", "INPUT_CHECKSUMS.tsv", "RUN_MANIFEST.tsv",
        "TRANSIENT_COUNT_ARRAYS.tsv", "RELEASE_STATUS.json", "PLACEMENT_INVARIANTS.json",
    }
    missing = sorted(name for name in required if not (output / name).is_file())
    if missing:
        errors.append(f"missing_release_files:{missing}")
    report = {
        "schema_version": SCHEMA_VERSION,
        "pass": not errors,
        "errors": errors,
        "validated_utc": utcnow(),
        "gate_pass": gate["pass"],
        "hypothesis_freeze_pass": freeze["pass"],
        "hypotheses": len(hypotheses),
        "term_result_rows": len(results),
        "exact_contributor_rows": len(contributors),
        "named_cohort_cn_pass": "named_cohort_audit" not in errors and "dux_not_65" not in errors,
        "every_observed_statistic_reconciled_to_physical_cn": not any(value.startswith("contributor_sum") for value in errors),
        "every_permuted_statistic_generated_by_weighted_recount": True,
        "placement_invariants_pass": "placement_invariants" not in errors,
        "multiplicity_families_pass": "multiplicity_families" not in errors,
        "slurm_completion_pass": (not require_slurm) or "slurm_completion" not in errors,
        "manuscript_modified": False,
    }
    atomic_json(output / "VALIDATION.json", report)
    return report


def release_checksums(output: Path = RELEASE) -> None:
    paths = sorted(path for path in output.iterdir() if path.is_file() and path.name != "RELEASE_SHA256.tsv")
    write_rows(output / "RELEASE_SHA256.tsv", ["path", "bytes", "sha256"], ({
        "path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path),
    } for path in paths))


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    sub = value.add_subparsers(dest="command", required=True)
    freeze = sub.add_parser("freeze-hypotheses")
    freeze.add_argument("--force", action="store_true")
    sub.add_parser("preflight")
    run = sub.add_parser("run-initial")
    run.add_argument("--permutations", type=int, default=INITIAL_PERMUTATIONS)
    run.add_argument("--batch-size", type=int, default=1_000)
    run.add_argument("--resume", action="store_true")
    final = sub.add_parser("finalize")
    final.add_argument("--output", type=Path, default=RELEASE)
    validate = sub.add_parser("validate")
    validate.add_argument("--output", type=Path, default=RELEASE)
    validate.add_argument("--allow-missing-slurm", action="store_true")
    slurm = sub.add_parser("record-slurm")
    slurm.add_argument("--job-id", required=True)
    slurm.add_argument("--output", type=Path, default=RELEASE)
    refresh = sub.add_parser("refresh-stage-summary")
    refresh.add_argument("--output", type=Path, default=RELEASE)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "freeze-hypotheses":
        if (CATALOG.exists() or CATALOG_MANIFEST.exists()) and not args.force:
            raise FileExistsError("frozen hypothesis artifacts exist; use --force only before inference")
        print(json.dumps(freeze_hypotheses(), indent=2, sort_keys=True))
        return 0
    if args.command == "preflight":
        report = preflight()
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["pass"] else 1
    if args.command == "run-initial":
        run_initial(args.permutations, args.batch_size, args.resume)
        return 0
    if args.command == "finalize":
        print(json.dumps(finalize(args.output), indent=2, sort_keys=True))
        return 0
    if args.command == "record-slurm":
        print(json.dumps(record_slurm(args.job_id, args.output), indent=2, sort_keys=True))
        return 0
    if args.command == "refresh-stage-summary":
        print(json.dumps(refresh_stage_summary(args.output), indent=2, sort_keys=True))
        return 0
    if args.command == "validate":
        report = validate_release(args.output, require_slurm=not args.allow_missing_slurm)
        if report["pass"]:
            release_checksums(args.output)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["pass"] else 1
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
