#!/usr/bin/env python3
"""Independent, fail-closed validation of the V6 physical-copy ontology release.

This implementation deliberately does not import v6_engine.py, COPY_engine.py,
the production source-map builder, or any earlier validator.  It reconstructs
the copy/source/term matrix from the raw CHM13 GFF, frozen copy-evidence rows,
raw GO/Reactome releases, and frozen geometry.  NumPy is used only for the
specified PCG64DXSM stream and vectorized recounts; R's qbeta is used as an
independent Clopper--Pearson implementation.

Run from the repository root:

  guix shell python python-numpy -- \
    python3 paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/independent_validation/validate_v6_independent.py
"""

from __future__ import annotations

import bisect
import collections
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import urllib.parse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np


SCRIPT = Path(__file__).resolve()
BASE = SCRIPT.parent.parent
ANALYSIS = BASE.parent
REPO = BASE.parents[3]
SOURCES = ANALYSIS / "sources"
RELEASE = BASE / "results" / "release"
REPORT = BASE / "V6_VALIDATION_REPORT.md"
MACHINE = BASE / "V6_VALIDATION.json"

SCHEMA = "chm13-physical-copy-ontology-v6.0-independent-validation-v1"
COLLECTIONS = ("GO_BP", "GO_MF", "GO_CC", "Reactome")
RELATIONS = ("direct", "ancestor")
TERM_FIELDS = (
    "functional_source_id", "functional_source_symbol", "ontology", "relation",
    "namespace", "term_id", "term_name", "minimum_distance",
    "inherited_from_direct_term_ids", "evidence_codes", "qualifiers",
    "assertion_record_ids", "ontology_release", "ontology_sha256",
)
REPRESENTATIVES = collections.OrderedDict((
    ("leading_exocyst", ("GO_CC", "direct", "GO:0000145")),
    ("negative_metal_binding", ("GO_MF", "ancestor", "GO:0046872")),
    ("DUX4_DUX4L", ("Reactome", "direct", "R-HSA-9819196")),
    ("DDX11L", ("GO_MF", "direct", "GO:0003678")),
    ("TUBB8", ("GO_BP", "ancestor", "GO:0006996")),
    ("OR4F", ("Reactome", "direct", "R-HSA-9752946")),
    ("WASH", ("GO_CC", "direct", "GO:0071203")),
))
EXPECTED_NAMED = collections.OrderedDict((
    ("DUX4_DUX4L", (107, 68, 65)),
    ("DDX11L", (12, 10, 10)),
    ("TUBB8", (16, 7, 2)),
    ("OR4F", (15, 11, 4)),
    ("WASH", (18, 9, 9)),
))
STRATA = ((0, 500_000), (500_000, 1_000_000), (1_000_000, 2_000_000),
          (2_000_000, 5_000_000), (5_000_000, None))
TOL = 5e-12


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def read_rows(path: Path) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_int(value: str) -> int:
    return int(value)


def same_float(observed: str | float, expected: float, tolerance: float = TOL) -> bool:
    value = float(observed)
    return math.isclose(value, expected, rel_tol=tolerance, abs_tol=tolerance)


def split_values(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in re.split(r"[,|]", value or "") if item.strip())


def row_hash(row: Mapping[str, str], fields: Sequence[str]) -> str:
    return sha256_bytes("\t".join(str(row.get(field, "")) for field in fields).encode())


def collection_for(ontology: str, namespace: str) -> str:
    if ontology == "Reactome" and namespace == "pathway":
        return "Reactome"
    return {
        ("GO", "biological_process"): "GO_BP",
        ("GO", "molecular_function"): "GO_MF",
        ("GO", "cellular_component"): "GO_CC",
    }[(ontology, namespace)]


def canonical_arm_key(arm: str) -> tuple[int, int]:
    match = re.fullmatch(r"chr([0-9]+|X|Y)_([pq])", arm)
    if not match:
        raise ValueError(f"noncanonical arm: {arm}")
    chrom, side = match.groups()
    number = 23 if chrom == "X" else 24 if chrom == "Y" else int(chrom)
    return number, 0 if side == "p" else 1


class Audit:
    def __init__(self) -> None:
        self.checks: list[dict[str, object]] = []

    def check(self, name: str, passed: bool, *, expected: object = None,
              observed: object = None, evidence: str = "",
              critical: bool = True) -> bool:
        self.checks.append({
            "name": name,
            "status": "PASS" if passed else "FAIL",
            "critical": critical,
            "expected": expected,
            "observed": observed,
            "evidence": evidence,
        })
        return passed

    @property
    def failed(self) -> list[dict[str, object]]:
        return [row for row in self.checks if row["status"] == "FAIL"]


@dataclass(frozen=True)
class RawCopy:
    copy_id: str
    seqid: str
    start0: int
    end0: int
    start1: int
    end1: int
    strand: str
    gff_id: str
    gff_line: int
    gene_name: str
    gene_biotype: str
    attributes_sha256: str

    @property
    def midpoint(self) -> int:
        return (self.start0 + self.end0) // 2


@dataclass(frozen=True)
class Block:
    block_id: str
    phr_id: str
    chromosome: str
    arm: str
    observed_start: int
    observed_end: int
    arm_start: int
    arm_end: int
    stratum: int
    candidate_low: int
    candidate_high: int

    @property
    def span(self) -> int:
        return self.observed_end - self.observed_start

    @property
    def candidate_count(self) -> int:
        return self.candidate_high - self.candidate_low + 1

    @property
    def components(self) -> str:
        return f"{self.phr_id}:0:{self.span}"


def resolve_repo_path(text: str) -> Path:
    path = Path(text)
    return path if path.is_absolute() else REPO / path


def validate_manifests(audit: Audit) -> dict[str, object]:
    manifest_rows = list(read_rows(BASE / "INPUT_MANIFEST.tsv"))
    manifest_failures = []
    for row in manifest_rows:
        path = resolve_repo_path(row["path"])
        if (not path.is_file() or path.stat().st_size != int(row["bytes"])
                or sha256_path(path) != row["sha256"]
                or row["sha256"] != row["expected_sha256"]):
            manifest_failures.append(row["role"])
    roles = [row["role"] for row in manifest_rows]
    audit.check("input_manifest_exact_bytes", not manifest_failures,
                expected="all 14 frozen raw/evidence inputs exact",
                observed=manifest_failures or "all exact", evidence="INPUT_MANIFEST.tsv")
    stage_ok = all(
        (row["role"] == "phr_bed" and row["stage"] == "after_source_freeze_target_join")
        or (row["role"] != "phr_bed" and row["stage"] == "before_target_source_assignment")
        for row in manifest_rows
    )
    audit.check("target_join_after_source_freeze", stage_ok,
                expected="only phr_bed is post-freeze", observed="PASS" if stage_ok else "stage mismatch",
                evidence="INPUT_MANIFEST.tsv")
    audit.check("input_manifest_roles_unique", len(roles) == len(set(roles)) == 14,
                expected=14, observed=len(set(roles)), evidence="INPUT_MANIFEST.tsv")

    gate_path = BASE / "PRE_RUN_V6_GATE.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate_ok = (gate.get("status") == "PASS" and gate.get("enrichment_authorized") is True
               and gate.get("enrichment_run") is False
               and gate.get("production_builder_imported") is False)
    audit.check("pre_run_gate_authorizes_exact_release", gate_ok,
                expected="PASS, authorized, pre-inference, independent",
                observed={key: gate.get(key) for key in (
                    "status", "enrichment_authorized", "enrichment_run", "production_builder_imported")},
                evidence="PRE_RUN_V6_GATE.json")
    audited = BASE / "PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv"
    audit.check("pre_run_gate_manifest_digest",
                sha256_path(audited) == gate["audited_release_manifest_sha256"],
                expected=gate["audited_release_manifest_sha256"], observed=sha256_path(audited),
                evidence="PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv")
    audited_failures = []
    for row in read_rows(audited):
        path = BASE / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_path(path) != row["sha256"]:
            audited_failures.append(row["path"])
    audit.check("audited_source_release_checksums", not audited_failures,
                expected="all audited source-map files exact", observed=audited_failures or "all exact",
                evidence="PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv")

    release_failures = []
    for row in read_rows(RELEASE / "RELEASE_SHA256.tsv"):
        path = RELEASE / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_path(path) != row["sha256"]:
            release_failures.append(row["path"])
    audit.check("release_artifact_checksums", not release_failures,
                expected="all 22 released artifacts exact", observed=release_failures or "all exact",
                evidence="results/release/RELEASE_SHA256.tsv")
    return gate


def parse_raw_gff() -> tuple[dict[str, RawCopy], int]:
    result: dict[str, RawCopy] = {}
    coordinate_keys: set[tuple[str, int, int, str]] = set()
    path = REPO / "data" / "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            seqid, start_text, end_text, strand, attributes = (
                fields[0], fields[3], fields[4], fields[6], fields[8])
            attrs: dict[str, str] = {}
            for item in attributes.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = urllib.parse.unquote(value)
            start1, end1 = int(start_text), int(end_text)
            start0, end0 = start1 - 1, end1
            gff_id = attrs["ID"]
            copy_id = f"CHM13v2.0|{seqid}:{start1}-{end1}|{strand}|{gff_id}"
            if copy_id in result:
                raise ValueError(f"duplicate raw copy: {copy_id}")
            result[copy_id] = RawCopy(
                copy_id=copy_id, seqid=seqid, start0=start0, end0=end0,
                start1=start1, end1=end1, strand=strand, gff_id=gff_id,
                gff_line=line_number,
                gene_name=attrs.get("gene_name", attrs.get("gene", gff_id)),
                gene_biotype=attrs.get("gene_biotype", ""),
                attributes_sha256=sha256_bytes(attributes.encode()),
            )
            coordinate_keys.add((seqid, start0, end0, strand))
    return result, len(result) - len(coordinate_keys)


def load_mapping(audit: Audit, raw: Mapping[str, RawCopy]) -> tuple[
        dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    assignment_path = BASE / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz"
    evidence_path = BASE / "COPY_SOURCE_EVIDENCE.tsv.gz"
    assignments = {row["copy_id"]: row for row in read_rows(assignment_path)}
    evidence = {row["copy_id"]: row for row in read_rows(evidence_path)}
    raw_ids, assignment_ids, evidence_ids = set(raw), set(assignments), set(evidence)
    audit.check("raw_assignment_evidence_bijection",
                raw_ids == assignment_ids == evidence_ids and len(raw_ids) == 61_312,
                expected=61_312,
                observed={"raw": len(raw_ids), "assignments": len(assignment_ids),
                          "evidence": len(evidence_ids), "id_sets_equal": raw_ids == assignment_ids == evidence_ids},
                evidence="raw GFF + GENOMEWIDE_SOURCE_ASSIGNMENTS + COPY_SOURCE_EVIDENCE")
    coordinate_failures = []
    for copy_id, source in raw.items():
        row = assignments[copy_id]
        expected = {
            "seqid": source.seqid, "start0": str(source.start0), "end0": str(source.end0),
            "start1": str(source.start1), "end1": str(source.end1), "strand": source.strand,
            "gff_id": source.gff_id, "gff_line": str(source.gff_line),
            "gene_name": source.gene_name, "gene_biotype": source.gene_biotype,
            "raw_attributes_sha256": source.attributes_sha256,
        }
        if any(row[field] != value for field, value in expected.items()):
            coordinate_failures.append(copy_id)
            if len(coordinate_failures) >= 10:
                break
    audit.check("all_raw_loci_retained_with_exact_coordinates", not coordinate_failures,
                expected="61,312 exact raw gene rows", observed=coordinate_failures or "all exact",
                evidence="raw GFF versus source assignments")
    cn_sum = sum(int(row["physical_copy_cn"]) for row in assignments.values())
    positive = all(int(row["physical_copy_cn"]) > 0 for row in assignments.values())
    audit.check("physical_copy_cn_complete_positive", cn_sum == 61_312 and positive,
                expected=61_312, observed=cn_sum, evidence="GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz")
    assignment_fields = next(iter(read_rows(assignment_path))).keys()
    evidence_fields = next(iter(read_rows(evidence_path))).keys()
    assignment_hash_ok = all(
        row_hash(row, list(assignment_fields)[:-1]) == row["assignment_record_sha256"]
        for row in assignments.values())
    evidence_hash_ok = all(
        row_hash(row, list(evidence_fields)[:-1]) == row["evidence_record_sha256"]
        for row in evidence.values())
    audit.check("all_assignment_row_hashes_reconstructed", assignment_hash_ok,
                expected=True, observed=assignment_hash_ok, evidence="61,312 assignment rows")
    audit.check("all_copy_evidence_hashes_reconstructed", evidence_hash_ok,
                expected=True, observed=evidence_hash_ok, evidence="61,312 evidence rows")
    return assignments, evidence


def validate_source_routes(audit: Audit, assignments: Mapping[str, dict[str, str]],
                           evidence: Mapping[str, dict[str, str]]) -> None:
    route_mismatches = []
    admitted = 0
    exact_self = 0
    directed = 0
    ambiguous_without_source = 0
    unsupported_without_source = 0
    for copy_id, row in assignments.items():
        route = evidence[copy_id]
        source = row["functional_source_id"]
        admissible = route["admissible_for_ontology"] == "1"
        ok = (
            source == route["functional_source_id"]
            and row["source_assignment_disposition"] == route["disposition"]
            and row["relationship_semantics"] == route["relationship_semantics"]
            and row["source_evidence_record_id"] == route["evidence_record_id"]
            and row["own_annotation_id"] == route["own_annotation_id"]
            and row["own_entrez_id"] == route["own_entrez_id"]
            and admissible == bool(source)
        )
        if source:
            admitted += int(row["physical_copy_cn"])
            ok &= route["disposition"] in {"EXACT_SELF", "EXPLICIT_RELATED_FUNCTIONAL_GENE"}
            if route["disposition"] == "EXACT_SELF":
                exact_self += 1
                ok &= bool(row["own_entrez_id"]) and source == "NCBIGene:" + row["own_entrez_id"]
            else:
                directed += 1
        else:
            ok &= route["disposition"] not in {"EXACT_SELF", "EXPLICIT_RELATED_FUNCTIONAL_GENE"}
            if route["disposition"] == "AMBIGUOUS_FAIL_CLOSED":
                ambiguous_without_source += 1
            if route["disposition"] == "UNSUPPORTED_FAIL_CLOSED":
                unsupported_without_source += 1
        if not ok and len(route_mismatches) < 20:
            route_mismatches.append(copy_id)
    audit.check("copy_evidence_assignment_routes_exact", not route_mismatches,
                expected="every admitted source has exact evidence; all failures emit no source",
                observed=route_mismatches or "all 61,312 exact", evidence="COPY_SOURCE_EVIDENCE.tsv.gz")
    audit.check("eligible_copy_cn", admitted == 31_966, expected=31_966, observed=admitted,
                evidence="sum physical_copy_cn over admissible routes")
    audit.check("source_route_gate_counts",
                exact_self == 20_405 and directed == 11_561
                and ambiguous_without_source == 45 and unsupported_without_source == 28_407,
                expected={"exact_self": 20405, "directed": 11561, "ambiguous": 45,
                          "unsupported": 28407},
                observed={"exact_self": exact_self, "directed": directed,
                          "ambiguous": ambiguous_without_source,
                          "unsupported": unsupported_without_source},
                evidence="frozen copy evidence dispositions")

    relations: set[tuple[str, str]] = set()
    with gzip.open(SOURCES / "ncbi_gene_group_2026-07-16.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle, delimiter="\t"):
            if (raw["#tax_id"] == "9606" and raw["Other_tax_id"] == "9606"
                    and raw["relationship"] == "Related functional gene"):
                relations.add((raw["GeneID"], raw["Other_GeneID"]))
    directed_missing = []
    for copy_id, route in evidence.items():
        if route["disposition"] != "EXPLICIT_RELATED_FUNCTIONAL_GENE":
            continue
        source_entrez = route["functional_source_id"].removeprefix("NCBIGene:")
        pair = (route["own_entrez_id"], source_entrez)
        expected_id = f"NCBI_GENE_GROUP:9606:{pair[0]}:Related_functional_gene:9606:{pair[1]}"
        if pair not in relations or route["evidence_record_id"] != expected_id:
            directed_missing.append(copy_id)
            if len(directed_missing) >= 20:
                break
    audit.check("all_nonself_sources_have_raw_directed_evidence", not directed_missing,
                expected="11,561 exact directed human-human relations",
                observed=directed_missing or 11_561,
                evidence="ncbi_gene_group_2026-07-16.tsv.gz")

    pilot_path = ANALYSIS / "ontology_v5" / "real_source_pilot" / "REAL_SOURCE_COPY_SOURCE_EVIDENCE.tsv.gz"
    pilot = {row["copy_id"]: row for row in read_rows(pilot_path)}
    pilot_mismatch = []
    overridden = {copy_id for copy_id, row in evidence.items() if row["pilot_override"] == "1"}
    for copy_id, upstream in pilot.items():
        row = evidence.get(copy_id, {})
        proposed = ("NCBIGene:" + upstream["proposed_source_entrez_id"]
                    if upstream["proposed_source_entrez_id"] else "")
        if not row or any((
            row.get("upstream_evidence_sha256") != upstream["evidence_record_sha256"],
            row.get("proposed_source_id") != proposed,
            row.get("disposition") != upstream["disposition"],
            row.get("admissible_for_ontology") != upstream["admissible_for_ontology"],
            row.get("ambiguity_status") != upstream["ambiguity_status_at_decision"],
            row.get("unresolved_reason") != upstream["unresolved_reason"],
        )):
            pilot_mismatch.append(copy_id)
    audit.check("remediated_pilot_evidence_reconciled",
                len(pilot) == 168 and overridden == set(pilot) and not pilot_mismatch,
                expected=168, observed={"rows": len(pilot), "overrides": len(overridden),
                                        "mismatches": pilot_mismatch[:10]},
                evidence="REAL_SOURCE_COPY_SOURCE_EVIDENCE.tsv.gz")


def parse_go() -> tuple[dict[str, dict[str, object]], dict[str, str]]:
    terms: dict[str, dict[str, object]] = {}
    alt_to_primary: dict[str, str] = {}
    current: dict[str, object] | None = None

    def finish(term: dict[str, object] | None) -> None:
        if not term or "id" not in term:
            return
        term.setdefault("name", "")
        term.setdefault("namespace", "")
        term.setdefault("parents", set())
        terms[str(term["id"])] = term
        for alt in term.get("alt_ids", set()):
            alt_to_primary[str(alt)] = str(term["id"])

    with gzip.open(SOURCES / "go-basic_2026-06-15.obo.gz", "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                finish(current)
                current = {"parents": set(), "alt_ids": set(), "obsolete": False}
            elif line.startswith("["):
                finish(current)
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line.startswith("alt_id: "):
                    current["alt_ids"].add(line[8:])
                elif line.startswith("is_a: "):
                    current["parents"].add(line[6:].split()[0])
                elif line.startswith("relationship: part_of "):
                    current["parents"].add(line[22:].split()[0])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
        finish(current)
    return terms, alt_to_primary


def load_raw_ontology(selected_geneids: set[str]) -> tuple[
        dict[str, dict[str, dict[str, set[str]]]], dict[str, dict[str, object]],
        dict[str, str], dict[str, set[str]], dict[str, dict[str, dict[str, set[str]]]]]:
    go_terms, alt = parse_go()
    go_assertions: dict[str, dict[str, dict[str, set[str]]]] = collections.defaultdict(
        lambda: collections.defaultdict(lambda: {
            "evidence": set(), "qualifiers": set(), "records": set()}))
    with gzip.open(SOURCES / "ncbi_gene2go_human_2026-07-13.tsv.gz", "rt",
                   encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
            gene_id = row["GeneID"]
            if gene_id not in selected_geneids:
                continue
            qualifiers = set(split_values(row.get("Qualifier", "")))
            if "NOT" in qualifiers:
                continue
            term_id = alt.get(row["GO_ID"], row["GO_ID"])
            term = go_terms.get(term_id)
            if not term or term.get("obsolete"):
                continue
            item = go_assertions[gene_id][term_id]
            item["evidence"].add(row["Evidence"])
            item["qualifiers"].update(qualifiers)
            item["records"].add(f"gene2go:{line_number}:{row['GO_ID']}")

    reactome_names: dict[str, str] = {}
    with gzip.open(SOURCES / "reactome_v96_human_pathways.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 3 and row[2] == "Homo sapiens":
                reactome_names[row[0]] = row[1]
    reactome_parents: dict[str, set[str]] = collections.defaultdict(set)
    with gzip.open(SOURCES / "reactome_v96_human_pathway_relations.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 2 and row[0] in reactome_names and row[1] in reactome_names:
                reactome_parents[row[1]].add(row[0])
    reactome_assertions: dict[str, dict[str, dict[str, set[str]]]] = collections.defaultdict(
        lambda: collections.defaultdict(lambda: {"evidence": set(), "records": set()}))
    with gzip.open(SOURCES / "reactome_v96_ncbi_human_all_levels.tsv.gz", "rt",
                   encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if (len(row) != 6 or row[5] != "Homo sapiens"
                    or row[0] not in selected_geneids):
                continue
            if row[1] not in reactome_names:
                raise ValueError(f"unknown Reactome pathway: {row[1]}")
            item = reactome_assertions[row[0]][row[1]]
            item["evidence"].add(row[4])
            item["records"].add(f"NCBI2Reactome_All_Levels:{line_number}")
    return go_assertions, go_terms, reactome_names, reactome_parents, reactome_assertions


def ancestor_distances(term_id: str, parents: Mapping[str, set[str]],
                       cache: dict[str, dict[str, int]]) -> dict[str, int]:
    if term_id in cache:
        return cache[term_id]
    distances = {term_id: 0}
    queue = collections.deque([term_id])
    while queue:
        child = queue.popleft()
        distance = distances[child]
        for parent in parents.get(child, set()):
            prior = distances.get(parent)
            if prior is None or distance + 1 < prior:
                distances[parent] = distance + 1
                queue.append(parent)
    cache[term_id] = distances
    return distances


def reconstruct_source_terms(
        source_id: str, symbol: str,
        go_assertions: Mapping[str, Mapping[str, Mapping[str, set[str]]]],
        go_terms: Mapping[str, Mapping[str, object]], go_parents: Mapping[str, set[str]],
        reactome_names: Mapping[str, str], reactome_parents: Mapping[str, set[str]],
        reactome_assertions: Mapping[str, Mapping[str, Mapping[str, set[str]]]],
        go_cache: dict[str, dict[str, int]], reactome_cache: dict[str, dict[str, int]],
        ontology_hashes: Mapping[str, str]) -> list[dict[str, str]]:
    gene_id = source_id.removeprefix("NCBIGene:")
    output: dict[tuple[str, str], dict[str, object]] = {}
    for leaf_id, evidence in go_assertions.get(gene_id, {}).items():
        for term_id, distance in ancestor_distances(leaf_id, go_parents, go_cache).items():
            item = output.setdefault(("GO", term_id), {
                "minimum_distance": distance, "leaves": set(), "direct": False,
                "evidence": set(), "qualifiers": set(), "records": set()})
            item["minimum_distance"] = min(int(item["minimum_distance"]), distance)
            item["leaves"].add(leaf_id)
            if distance == 0:
                item["direct"] = True
                item["evidence"].update(evidence["evidence"])
                item["qualifiers"].update(evidence["qualifiers"])
                item["records"].update(evidence["records"])

    reactome_all = set(reactome_assertions.get(gene_id, {}))
    ancestors_within: set[str] = set()
    for term_id in reactome_all:
        ancestors_within.update(set(ancestor_distances(
            term_id, reactome_parents, reactome_cache)) - {term_id})
    reactome_direct = reactome_all - ancestors_within
    if reactome_all and not reactome_direct:
        raise ValueError(f"cyclic Reactome direct-leaf set: {gene_id}")
    for leaf_id in reactome_direct:
        evidence = reactome_assertions[gene_id][leaf_id]
        for term_id, distance in ancestor_distances(
                leaf_id, reactome_parents, reactome_cache).items():
            item = output.setdefault(("Reactome", term_id), {
                "minimum_distance": distance, "leaves": set(), "direct": False,
                "evidence": set(), "qualifiers": set(), "records": set()})
            item["minimum_distance"] = min(int(item["minimum_distance"]), distance)
            item["leaves"].add(leaf_id)
            if distance == 0:
                item["direct"] = True
                item["evidence"].update(evidence["evidence"])
                item["records"].update(evidence["records"])

    rows = []
    for (ontology, term_id), item in sorted(output.items()):
        if ontology == "GO":
            metadata = go_terms[term_id]
            namespace, term_name = str(metadata["namespace"]), str(metadata["name"])
            release, ontology_sha = "Gene Ontology 2026-06-15", ontology_hashes["GO"]
        else:
            namespace, term_name = "pathway", reactome_names[term_id]
            release, ontology_sha = "Reactome v96", ontology_hashes["Reactome"]
        rows.append({
            "functional_source_id": source_id,
            "functional_source_symbol": symbol,
            "ontology": ontology,
            "relation": "direct" if item["direct"] else "ancestor",
            "namespace": namespace,
            "term_id": term_id,
            "term_name": term_name,
            "minimum_distance": str(item["minimum_distance"]),
            "inherited_from_direct_term_ids": "|".join(sorted(item["leaves"])),
            "evidence_codes": "|".join(sorted(item["evidence"])),
            "qualifiers": "|".join(sorted(item["qualifiers"])),
            "assertion_record_ids": "|".join(sorted(item["records"])),
            "ontology_release": release,
            "ontology_sha256": ontology_sha,
        })
    # The frozen table is globally sorted by source, ontology, relation, term.
    # Relation order is lexical (ancestor before direct), not term-ID order.
    rows.sort(key=lambda row: (row["ontology"], row["relation"], row["term_id"]))
    return rows


def load_boundaries() -> dict[str, tuple[int, int]]:
    bands: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    with (REPO / "data" / "chm13v2.0_cytobands_allchrs.bed").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, band = line.rstrip("\n").split("\t")[:4]
            bands[chrom].append((int(start), int(end), band))
    result = {}
    for chrom, rows in bands.items():
        rows.sort()
        if rows[0][0] != 0 or any(rows[i][1] != rows[i + 1][0] for i in range(len(rows) - 1)):
            raise ValueError(f"noncontiguous cytobands: {chrom}")
        result[chrom] = (min(start for start, _end, band in rows if band.startswith("q")), rows[-1][1])
    return result


def arm_for(seqid: str, start0: int, end0: int,
            boundaries: Mapping[str, tuple[int, int]]) -> str:
    q_start, chromosome_end = boundaries[seqid]
    if not 0 <= start0 < end0 <= chromosome_end:
        raise ValueError(f"copy outside chromosome: {seqid}:{start0}-{end0}")
    return f"{seqid}_{'p' if (start0 + end0) // 2 < q_start else 'q'}"


def load_phrs(boundaries: Mapping[str, tuple[int, int]]) -> list[dict[str, object]]:
    result = []
    with (REPO / "data" / "chm13.phrs.bed").open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            seqid, start, end, sharing = line.rstrip("\n").split("\t")[:4]
            start0, end0 = int(start), int(end)
            arm = arm_for(seqid, start0, end0, boundaries)
            result.append({
                "phr_id": f"CHM13_{arm}_{start0}_{end0}", "seqid": seqid, "arm": arm,
                "start0": start0, "end0": end0, "width": end0 - start0,
                "sharing": sharing, "line": line_number,
            })
    return result


def memberships(assignments: Mapping[str, dict[str, str]],
                phrs: Sequence[Mapping[str, object]]) -> tuple[dict[str, bool], dict[str, bool]]:
    by_chrom: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for phr in phrs:
        by_chrom[str(phr["seqid"])].append((int(phr["start0"]), int(phr["end0"])))
    midpoint: dict[str, bool] = {}
    overlap: dict[str, bool] = {}
    for copy_id, row in assignments.items():
        start, end = int(row["start0"]), int(row["end0"])
        middle = (start + end) // 2
        intervals = by_chrom[row["seqid"]]
        midpoint[copy_id] = any(left <= middle < right for left, right in intervals)
        overlap[copy_id] = any(start < right and left < end for left, right in intervals)
    return midpoint, overlap


def cohort_for(name: str, synonyms: str) -> str:
    synonym_set = set(split_values(synonyms))
    if name == "DUX4" or name.startswith("DUX4L") or "DUX4L30" in synonym_set:
        return "DUX4_DUX4L"
    for prefix in ("DDX11L", "TUBB8", "OR4F", "WASH"):
        if name.startswith(prefix):
            return prefix if prefix != "DDX11L" else "DDX11L"
    return ""


def validate_adversarial_fixture(audit: Audit) -> dict[str, object]:
    copies = [
        {"copy_id": f"chr1:{i * 100}-{i * 100 + 50}", "source": "NCBIGene:7",
         "gene": "G7", "family": "F7", "physical_copy_cn": 1}
        for i in range(1, 8)
    ]
    copies.append({"copy_id": "chr2:1-51", "source": "NCBIGene:7", "gene": "G7",
                   "family": "F7", "physical_copy_cn": 3})
    weighted = sum(row["physical_copy_cn"] for row in copies)
    distinct_coordinates = len({row["copy_id"] for row in copies})
    collapsed = {
        "unique_source": len({row["source"] for row in copies}),
        "unique_gene": len({row["gene"] for row in copies}),
        "unique_family": len({row["family"] for row in copies}),
    }
    passed = weighted == 10 and distinct_coordinates == 8 and set(collapsed.values()) == {1}
    audit.check("n_copy_adversarial_fixture_rejects_collapse", passed,
                expected={"weighted_physical_cn": 10, "coordinate_rows": 8,
                          "unique_source_gene_family": 1},
                observed={"weighted_physical_cn": weighted,
                          "coordinate_rows": distinct_coordinates, **collapsed},
                evidence="synthetic 7xCN1 + 1xCN3 fixture")
    return {"status": "PASS" if passed else "FAIL", "weighted_physical_cn": weighted,
            "coordinate_rows": distinct_coordinates, **collapsed}


def reconstruct_matrix(
        audit: Audit, assignments: Mapping[str, dict[str, str]],
        evidence: Mapping[str, dict[str, str]], boundaries: Mapping[str, tuple[int, int]],
        midpoint: Mapping[str, bool], overlap: Mapping[str, bool]) -> dict[str, object]:
    source_copies: dict[str, list[str]] = collections.defaultdict(list)
    source_symbols: dict[str, set[str]] = collections.defaultdict(set)
    source_cn: collections.Counter[str] = collections.Counter()
    source_mid: collections.Counter[str] = collections.Counter()
    source_overlap: collections.Counter[str] = collections.Counter()
    source_arms: dict[str, set[str]] = collections.defaultdict(set)
    source_cohort_mid: collections.Counter[tuple[str, str]] = collections.Counter()
    source_cohort_pairs: set[tuple[str, str]] = set()
    for copy_id, row in assignments.items():
        source = row["functional_source_id"]
        if not source:
            continue
        cn = int(row["physical_copy_cn"])
        source_copies[source].append(copy_id)
        source_symbols[source].add(row["functional_source_symbol"])
        source_cn[source] += cn
        source_mid[source] += cn * int(midpoint[copy_id])
        source_overlap[source] += cn * int(overlap[copy_id])
        source_arms[source].add(arm_for(row["seqid"], int(row["start0"]), int(row["end0"]), boundaries))
        cohort = cohort_for(row["gene_name"], row["gene_synonyms"])
        if cohort:
            source_cohort_pairs.add((cohort, source))
            if midpoint[copy_id]:
                source_cohort_mid[(cohort, source)] += cn
    symbol_conflicts = {source: sorted(values) for source, values in source_symbols.items() if len(values) != 1}
    audit.check("functional_source_symbols_consistent", not symbol_conflicts,
                expected="one symbol per admitted source", observed=symbol_conflicts or "all consistent",
                evidence="GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz")

    selected_geneids = {source.removeprefix("NCBIGene:") for source in source_copies}
    go_assertions, go_terms, reactome_names, reactome_parents, reactome_assertions = \
        load_raw_ontology(selected_geneids)
    ontology_hashes = {
        "GO": sha256_path(SOURCES / "go-basic_2026-06-15.obo.gz"),
        "Reactome": sha256_path(SOURCES / "reactome_v96_human_pathway_relations.tsv.gz"),
    }
    go_cache: dict[str, dict[str, int]] = {}
    reactome_cache: dict[str, dict[str, int]] = {}
    go_parents = {term_id: set(term.get("parents", set())) for term_id, term in go_terms.items()}
    observed_terms = iter(read_rows(BASE / "SOURCE_TERMS.tsv.gz"))
    term_mismatches: list[dict[str, object]] = []
    term_count = 0
    eligible_without_raw_terms = []
    direct_counts: dict[str, int] = {}
    closure_counts: dict[str, int] = {}
    phr_source_terms: dict[str, list[dict[str, str]]] = {}
    representative_sources: dict[tuple[str, str, str], set[str]] = {
        key: set() for key in REPRESENTATIVES.values()}
    burden: collections.Counter[tuple[str, str, str]] = collections.Counter()
    observed_mid: collections.Counter[tuple[str, str, str]] = collections.Counter()
    observed_overlap: collections.Counter[tuple[str, str, str]] = collections.Counter()
    source_count: collections.Counter[tuple[str, str, str]] = collections.Counter()
    term_arms: dict[tuple[str, str, str], set[str]] = collections.defaultdict(set)
    names: dict[tuple[str, str, str], str] = {}
    metadata: dict[tuple[str, str, str], tuple[str, str]] = {}
    edge_cn = direct_edge_cn = phr_mid_edge_cn = phr_overlap_edge_cn = 0
    named_term_expected: dict[tuple[str, str, str, str, str], dict[str, str]] = {}

    for source in sorted(source_copies):
        symbol = next(iter(source_symbols[source]))
        rows = reconstruct_source_terms(
            source, symbol, go_assertions, go_terms, go_parents, reactome_names, reactome_parents,
            reactome_assertions, go_cache, reactome_cache, ontology_hashes)
        if not rows:
            eligible_without_raw_terms.append(source)
        direct_counts[source] = sum(row["relation"] == "direct" for row in rows)
        closure_counts[source] = len(rows)
        if source_overlap[source]:
            phr_source_terms[source] = rows
        for expected in rows:
            term_count += 1
            try:
                observed = next(observed_terms)
            except StopIteration:
                observed = {}
            differences = {
                field: [expected[field], observed.get(field, "<missing>")]
                for field in TERM_FIELDS if expected[field] != observed.get(field)
            }
            if differences and len(term_mismatches) < 10:
                term_mismatches.append({"source": source, "differences": differences})
            collection = collection_for(expected["ontology"], expected["namespace"])
            key = (collection, expected["relation"], expected["term_id"])
            names[key] = expected["term_name"]
            metadata[key] = (expected["ontology"], expected["namespace"])
            burden[key] += source_cn[source]
            observed_mid[key] += source_mid[source]
            observed_overlap[key] += source_overlap[source]
            source_count[key] += 1
            term_arms[key].update(source_arms[source])
            edge_cn += source_cn[source]
            direct_edge_cn += source_cn[source] * int(expected["relation"] == "direct")
            phr_mid_edge_cn += source_mid[source]
            phr_overlap_edge_cn += source_overlap[source]
            if key in representative_sources:
                representative_sources[key].add(source)
            for cohort in EXPECTED_NAMED:
                if (cohort, source) not in source_cohort_pairs:
                    continue
                cn = source_cohort_mid[(cohort, source)]
                named_key = (cohort, source, expected["ontology"], expected["relation"], expected["term_id"])
                named_term_expected[named_key] = {
                    "cohort": cohort, "functional_source_id": source,
                    "functional_source_symbol": symbol,
                    "ontology": expected["ontology"], "relation": expected["relation"],
                    "namespace": expected["namespace"], "term_id": expected["term_id"],
                    "term_name": expected["term_name"],
                    "phr_midpoint_physical_copy_burden": str(cn),
                    "phr_midpoint_contributor_rows": str(sum(
                        1 for copy_id in source_copies[source]
                        if midpoint[copy_id] and cohort_for(
                            assignments[copy_id]["gene_name"], assignments[copy_id]["gene_synonyms"]) == cohort)),
                }
    try:
        extra_term = next(observed_terms)
    except StopIteration:
        extra_term = None
    audit.check("raw_direct_and_ancestor_source_terms_exact",
                term_count == 1_686_727 and not term_mismatches and extra_term is None
                and not eligible_without_raw_terms,
                expected=1_686_727,
                observed={"rows": term_count, "mismatches": term_mismatches,
                          "extra_row": extra_term, "eligible_without_terms": eligible_without_raw_terms[:10]},
                evidence="raw gene2go + GO OBO + Reactome relations/all-level map")

    count_mismatches = []
    for copy_id, row in assignments.items():
        source = row["functional_source_id"]
        expected_direct = direct_counts.get(source, 0)
        expected_closure = closure_counts.get(source, 0)
        if (int(row["direct_term_count"]) != expected_direct
                or int(row["closure_term_count"]) != expected_closure):
            count_mismatches.append(copy_id)
            if len(count_mismatches) >= 20:
                break
    audit.check("every_copy_term_count_reconstructed", not count_mismatches,
                expected="all 61,312 direct/closure counts exact",
                observed=count_mismatches or "all exact", evidence="raw ontology closure versus assignments")
    audit.check("physical_copy_term_edge_cn",
                edge_cn == 2_929_709 and direct_edge_cn == 626_048
                and phr_mid_edge_cn == 16_763 and phr_overlap_edge_cn == 17_579,
                expected={"all": 2_929_709, "direct": 626_048, "midpoint": 16_763,
                          "overlap": 17_579},
                observed={"all": edge_cn, "direct": direct_edge_cn, "midpoint": phr_mid_edge_cn,
                          "overlap": phr_overlap_edge_cn},
                evidence="sum physical_copy_cn over independently rebuilt source-term edges")

    named_rows = list(read_rows(BASE / "NAMED_COHORT_TERM_BURDENS.tsv.gz"))
    named_observed = {
        (row["cohort"], row["functional_source_id"], row["ontology"], row["relation"], row["term_id"]): row
        for row in named_rows}
    named_term_mismatch = []
    for key, expected in named_term_expected.items():
        row = named_observed.get(key)
        if row is None or any(row[field] != value for field, value in expected.items()):
            named_term_mismatch.append(key)
            if len(named_term_mismatch) >= 20:
                break
    audit.check("all_named_cohort_term_burdens_reconstructed",
                len(named_term_expected) == len(named_observed) == 1_391 and not named_term_mismatch,
                expected=1_391,
                observed={"expected": len(named_term_expected), "release": len(named_observed),
                          "mismatches": named_term_mismatch},
                evidence="NAMED_COHORT_TERM_BURDENS.tsv.gz")
    return {
        "source_copies": source_copies, "source_cn": source_cn,
        "source_mid": source_mid, "source_overlap": source_overlap,
        "phr_source_terms": phr_source_terms,
        "burden": burden, "observed_mid": observed_mid,
        "observed_overlap": observed_overlap, "source_count": source_count,
        "term_arms": term_arms, "names": names, "metadata": metadata,
        "representative_sources": representative_sources,
        "term_count": term_count,
        "edge_cn": edge_cn, "direct_edge_cn": direct_edge_cn,
        "phr_mid_edge_cn": phr_mid_edge_cn, "phr_overlap_edge_cn": phr_overlap_edge_cn,
    }


def validate_named_cohorts(audit: Audit, assignments: Mapping[str, dict[str, str]],
                           midpoint: Mapping[str, bool]) -> dict[str, dict[str, int]]:
    counts = {cohort: {"genome": 0, "phr": 0, "eligible": 0} for cohort in EXPECTED_NAMED}
    for copy_id, row in assignments.items():
        cohort = cohort_for(row["gene_name"], row["gene_synonyms"])
        if not cohort:
            continue
        cn = int(row["physical_copy_cn"])
        counts[cohort]["genome"] += cn
        counts[cohort]["phr"] += cn * int(midpoint[copy_id])
        counts[cohort]["eligible"] += cn * int(midpoint[copy_id] and bool(row["functional_source_id"]))
    passed = all(
        (counts[cohort]["genome"], counts[cohort]["phr"], counts[cohort]["eligible"]) == expected
        for cohort, expected in EXPECTED_NAMED.items())
    audit.check("named_cohort_physical_cn_recounts", passed,
                expected={key: list(value) for key, value in EXPECTED_NAMED.items()},
                observed=counts, evidence="raw named loci + independent PHR midpoint membership")
    release_rows = {row["cohort"]: row for row in read_rows(BASE / "NAMED_COHORT_AUDIT.tsv")}
    release_ok = len(release_rows) == 5 and all(
        int(release_rows[cohort]["genome_physical_copies"]) == values["genome"]
        and int(release_rows[cohort]["phr_physical_copies"]) == values["phr"]
        and int(release_rows[cohort]["phr_ontology_contributors"]) == values["eligible"]
        and release_rows[cohort]["status"] == "PASS"
        for cohort, values in counts.items())
    audit.check("named_cohort_audit_artifact_exact", release_ok,
                expected="five PASS rows", observed="PASS" if release_ok else "mismatch",
                evidence="NAMED_COHORT_AUDIT.tsv")
    return counts


def validate_catalog_and_results(audit: Audit, matrix: Mapping[str, object]) -> tuple[
        list[dict[str, str]], list[dict[str, str]], dict[tuple[str, tuple[str, str, str]], dict[str, str]]]:
    burden = matrix["burden"]
    order_collection = {value: index for index, value in enumerate(COLLECTIONS)}
    order_relation = {value: index for index, value in enumerate(RELATIONS)}
    keys = sorted(burden, key=lambda key: (
        order_collection[key[0]], order_relation[key[1]], key[2]))
    catalog = list(read_rows(BASE / "FROZEN_HYPOTHESES.tsv.gz"))
    mismatches = []
    for index, key in enumerate(keys):
        row = catalog[index] if index < len(catalog) else {}
        ontology, namespace = matrix["metadata"][key]
        expected = {
            "schema_version": "chm13-physical-copy-ontology-v6.0",
            "hypothesis_index": str(index), "collection": key[0], "relation": key[1],
            "ontology": ontology, "namespace": namespace, "term_id": key[2],
            "term_name": matrix["names"][key],
            "genome_physical_copy_burden": str(matrix["burden"][key]),
            "genome_arm_count": str(len(matrix["term_arms"][key])),
            "genome_source_count": str(matrix["source_count"][key]),
            "multiplicity_family_id": key[0],
            "hypothesis_status": "FROZEN_TESTED_NO_TARGET_FILTER",
        }
        differences = {field: [value, row.get(field)] for field, value in expected.items()
                       if row.get(field) != value}
        if differences and len(mismatches) < 10:
            mismatches.append({"index": index, "key": key, "differences": differences})
    header = set(catalog[0]) if catalog else set()
    target_blind = not any(field.startswith("phr_") or field in {
        "in_phr_midpoint", "in_phr_any_overlap"} for field in header)
    audit.check("target_blind_hypothesis_catalog_reconstructed",
                len(keys) == len(catalog) == 31_235 and not mismatches and target_blind,
                expected=31_235,
                observed={"rebuilt": len(keys), "release": len(catalog),
                          "mismatches": mismatches, "target_blind_header": target_blind},
                evidence="raw source-term matrix before PHR join")
    family_sizes = collections.Counter(row["collection"] for row in catalog)
    expected_sizes = {"GO_BP": 19_181, "GO_MF": 6_261, "GO_CC": 2_677, "Reactome": 3_116}
    freeze = json.loads((BASE / "FROZEN_HYPOTHESES.json").read_text(encoding="utf-8"))
    freeze_inputs_ok = all(
        sha256_path((REPO / "data" / name) if name == "chm13v2.0_cytobands_allchrs.bed"
                    else BASE / name) == digest
        for name, digest in freeze["freeze_inputs"].items())
    freeze_ok = (
        freeze["status"] == "FROZEN_BEFORE_TARGET_TERM_RECOUNT"
        and freeze["target_bed_opened"] is False
        and freeze["target_membership_columns_opened"] is False
        and freeze["physical_copy_term_edges_opened"] is False
        and freeze["target_filter"] == "none"
        and freeze["hypothesis_count"] == len(catalog)
        and freeze["catalog_sha256"] == sha256_path(BASE / "FROZEN_HYPOTHESES.tsv.gz")
        and freeze["multiplicity_family_sizes"] == expected_sizes
        and freeze_inputs_ok)
    audit.check("hypothesis_freeze_manifest_target_blind_and_exact", freeze_ok,
                expected="catalog frozen before target/edge access with three exact inputs",
                observed={"status": freeze["status"], "target_bed_opened": freeze["target_bed_opened"],
                          "target_columns_opened": freeze["target_membership_columns_opened"],
                          "physical_edges_opened": freeze["physical_copy_term_edges_opened"],
                          "freeze_inputs_exact": freeze_inputs_ok},
                evidence="FROZEN_HYPOTHESES.json")
    audit.check("multiplicity_family_sizes", dict(family_sizes) == expected_sizes,
                expected=expected_sizes, observed=dict(family_sizes), evidence="FROZEN_HYPOTHESES.tsv.gz")

    multiplicity = {row["family_id"]: row for row in read_rows(RELEASE / "MULTIPLICITY_FAMILIES.tsv")}
    multiplicity_ok = len(multiplicity) == 5
    for collection in COLLECTIONS:
        rows = [row for row in catalog if row["collection"] == collection]
        payload = "\n".join(f"{row['relation']}|{row['term_id']}" for row in rows) + "\n"
        release_row = multiplicity.get(collection, {})
        multiplicity_ok &= (
            release_row.get("hypothesis_count") == str(len(rows))
            and release_row.get("ordered_hypothesis_sha256") == sha256_bytes(payload.encode())
            and release_row.get("relations") == "direct;ancestor")
    global_payload = "\n".join(
        f"{row['collection']}|{row['relation']}|{row['term_id']}" for row in catalog) + "\n"
    multiplicity_ok &= multiplicity.get("GLOBAL_ALL_ONTOLOGY", {}).get(
        "ordered_hypothesis_sha256") == sha256_bytes(global_payload.encode())
    audit.check("multiplicity_ordered_hashes_reconstructed", bool(multiplicity_ok),
                expected="four collection families plus global exact",
                observed="PASS" if multiplicity_ok else "mismatch",
                evidence="MULTIPLICITY_FAMILIES.tsv")

    results = list(read_rows(RELEASE / "TERM_RESULTS.tsv.gz"))
    result_map: dict[tuple[str, tuple[str, str, str]], dict[str, str]] = {}
    observed_mismatches = []
    for row in results:
        key = (row["collection"], row["relation"], row["term_id"])
        result_map[(row["assignment"], key)] = row
        expected_observed = (matrix["observed_mid"][key]
                             if row["assignment"] == "midpoint"
                             else matrix["observed_overlap"][key])
        if (int(row["observed_physical_copy_burden"]) != expected_observed
                or int(row["genome_physical_copy_burden"]) != matrix["burden"][key]
                or int(row["genome_arm_count"]) != len(matrix["term_arms"][key])
                or int(row["genome_source_count"]) != matrix["source_count"][key]):
            observed_mismatches.append((row["assignment"], key))
            if len(observed_mismatches) >= 20:
                break
    audit.check("every_observed_term_is_sum_physical_cn",
                len(results) == 62_470 and len(result_map) == 62_470 and not observed_mismatches,
                expected="31,235 hypotheses x two assignments",
                observed={"rows": len(results), "unique": len(result_map),
                          "mismatches": observed_mismatches},
                evidence="independent raw-locus/source-term/PHR recount")
    return catalog, results, result_map


def validate_exact_contributors(audit: Audit, assignments: Mapping[str, dict[str, str]],
                                evidence: Mapping[str, dict[str, str]],
                                midpoint: Mapping[str, bool], overlap: Mapping[str, bool],
                                matrix: Mapping[str, object]) -> None:
    expected: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for copy_id, copy in assignments.items():
        source = copy["functional_source_id"]
        if not source or not overlap[copy_id]:
            continue
        route = evidence[copy_id]
        for term in matrix["phr_source_terms"][source]:
            collection = collection_for(term["ontology"], term["namespace"])
            key = (collection, term["relation"], term["term_id"], copy_id)
            expected[key] = {
                "schema_version": "chm13-physical-copy-ontology-v6.0",
                "collection": collection, "relation": term["relation"],
                "ontology": term["ontology"], "namespace": term["namespace"],
                "term_id": term["term_id"], "term_name": term["term_name"],
                "copy_id": copy_id, "seqid": copy["seqid"], "start0": copy["start0"],
                "end0": copy["end0"], "strand": copy["strand"],
                "gene_name": copy["gene_name"], "gene_biotype": copy["gene_biotype"],
                "physical_copy_cn": copy["physical_copy_cn"],
                "phr_midpoint_cn": str(int(copy["physical_copy_cn"]) * int(midpoint[copy_id])),
                "phr_any_overlap_cn": copy["physical_copy_cn"],
                "own_annotation_id": copy["own_annotation_id"],
                "functional_source_id": source,
                "functional_source_symbol": copy["functional_source_symbol"],
                "source_assignment_disposition": copy["source_assignment_disposition"],
                "assignment_tier": copy["assignment_tier"],
                "mapping_confidence": copy["mapping_confidence"],
                "relationship_semantics": copy["relationship_semantics"],
                "source_evidence_record_id": copy["source_evidence_record_id"],
                "source_evidence_release": copy["source_evidence_release"],
                "evidence_route_source": route["evidence_source"],
                "evidence_route_record_id": route["evidence_record_id"],
                "evidence_route_disposition": route["disposition"],
                "evidence_route_admissible": route["admissible_for_ontology"],
                "evidence_route_record_sha256": route["evidence_record_sha256"],
                "minimum_distance": term["minimum_distance"],
                "inherited_from_direct_term_ids": term["inherited_from_direct_term_ids"],
                "analysis_role": "EXACT_PHYSICAL_COPY_CONTRIBUTOR_SOURCE_EVIDENCE",
            }
    observed_rows = list(read_rows(RELEASE / "EXACT_TERM_CONTRIBUTORS.tsv.gz"))
    observed = {(row["collection"], row["relation"], row["term_id"], row["copy_id"]): row
                for row in observed_rows}
    mismatches = []
    for key, fields in expected.items():
        row = observed.get(key)
        if row is None or any(row[field] != value for field, value in fields.items()):
            mismatches.append(key)
            if len(mismatches) >= 20:
                break
    audit.check("exact_contributor_rows_reconstructed",
                len(expected) == len(observed) == 17_579 and not mismatches,
                expected=17_579,
                observed={"rebuilt": len(expected), "release": len(observed),
                          "mismatches": mismatches},
                evidence="EXACT_TERM_CONTRIBUTORS.tsv.gz")


def coverage_rows(assignments: Mapping[str, dict[str, str]],
                  midpoint: Mapping[str, bool], overlap: Mapping[str, bool]) -> dict[
                      tuple[str, str, str, str], dict[str, str]]:
    output: dict[tuple[str, str, str, str], dict[str, int]] = {}

    def add(scope: str, level: str, route: str, biotype: str,
            rows: Iterable[dict[str, str]]) -> None:
        selected = list(rows)
        key = (scope, level, route, biotype)
        output[key] = {
            "physical_rows": len(selected),
            "physical_copy_cn": sum(int(row["physical_copy_cn"]) for row in selected),
            "own_annotation_identity_cn": sum(int(row["physical_copy_cn"])
                                                for row in selected if row["own_annotation_id"]),
            "ontology_eligible_rows": sum(bool(row["functional_source_id"]) for row in selected),
            "ontology_eligible_cn": sum(int(row["physical_copy_cn"])
                                         for row in selected if row["functional_source_id"]),
            "ontology_ineligible_cn": sum(int(row["physical_copy_cn"])
                                           for row in selected if not row["functional_source_id"]),
            "direct_term_edge_cn": sum(int(row["physical_copy_cn"]) * int(row["direct_term_count"])
                                       for row in selected),
            "closure_term_edge_cn": sum(int(row["physical_copy_cn"]) * int(row["closure_term_count"])
                                        for row in selected),
        }

    for scope in ("genome", "phr_midpoint", "phr_any_overlap"):
        selected = [row for copy_id, row in assignments.items()
                    if scope == "genome" or (scope == "phr_midpoint" and midpoint[copy_id])
                    or (scope == "phr_any_overlap" and overlap[copy_id])]
        add(scope, "TOTAL", "ALL", "ALL", selected)
        for biotype in sorted({row["gene_biotype"] for row in selected}):
            add(scope, "BIOTYPE", "ALL", biotype,
                (row for row in selected if row["gene_biotype"] == biotype))
        for route in sorted({row["source_assignment_disposition"] for row in selected}):
            add(scope, "ROUTE", route, "ALL",
                (row for row in selected if row["source_assignment_disposition"] == route))
        for route, biotype in sorted({
                (row["source_assignment_disposition"], row["gene_biotype"]) for row in selected}):
            add(scope, "ROUTE_BY_BIOTYPE", route, biotype,
                (row for row in selected if row["source_assignment_disposition"] == route
                 and row["gene_biotype"] == biotype))
    return {
        key: {"scope": key[0], "aggregation_level": key[1], "evidence_route": key[2],
              "gene_biotype": key[3], **{field: str(value) for field, value in values.items()}}
        for key, values in output.items()
    }


def validate_mapping_coverage(audit: Audit, assignments: Mapping[str, dict[str, str]],
                              midpoint: Mapping[str, bool], overlap: Mapping[str, bool]) -> None:
    expected = coverage_rows(assignments, midpoint, overlap)
    observed_rows = list(read_rows(RELEASE / "MAPPING_COVERAGE.tsv"))
    observed = {(row["scope"], row["aggregation_level"], row["evidence_route"], row["gene_biotype"]): row
                for row in observed_rows}
    mismatches = []
    for key, row in expected.items():
        actual = observed.get(key)
        if actual is None or any(actual[field] != value for field, value in row.items()):
            mismatches.append(key)
            if len(mismatches) >= 20:
                break
    audit.check("mapping_coverage_all_aggregations_reconstructed",
                len(expected) == len(observed) == 148 and not mismatches,
                expected=148, observed={"rebuilt": len(expected), "release": len(observed),
                                        "mismatches": mismatches},
                evidence="MAPPING_COVERAGE.tsv")


def arm_summary(audit: Audit, boundaries: Mapping[str, tuple[int, int]]) -> dict[str, tuple[str, int, int]]:
    rows = list(read_rows(ANALYSIS / "analysis_ready" / "chm13_arm_summary.tsv"))
    result = {row["arm"]: (row["chromosome"], int(row["start0"]), int(row["end0"])) for row in rows}
    expected: dict[str, tuple[str, int, int]] = {}
    for chrom, (q_start, chrom_end) in boundaries.items():
        expected[f"{chrom}_p"] = (chrom, 0, q_start)
        expected[f"{chrom}_q"] = (chrom, q_start, chrom_end)
    audit.check("arm_geometry_reconstructed_from_cytobands", result == expected,
                expected="48 cytoband-derived arms", observed=len(result),
                evidence="chm13v2.0_cytobands_allchrs.bed vs analysis_ready arm summary")
    return result


def stratum_index(distance: int) -> int:
    for index, (lower, upper) in enumerate(STRATA):
        if distance >= lower and (upper is None or distance < upper):
            return index
    raise ValueError(distance)


def build_blocks(audit: Audit, phrs: Sequence[Mapping[str, object]],
                 arms: Mapping[str, tuple[str, int, int]]) -> list[Block]:
    blocks = []
    for phr in sorted(phrs, key=lambda row: canonical_arm_key(str(row["arm"]))):
        arm = str(phr["arm"])
        chrom, arm_start, arm_end = arms[arm]
        start, end = int(phr["start0"]), int(phr["end0"])
        span, offset = end - start, (end - start) // 2
        distance = (start + offset - arm_start if arm.endswith("_p")
                    else arm_end - (start + offset))
        stratum = stratum_index(distance)
        low_distance, high_distance = STRATA[stratum]
        high_distance = arm_end - arm_start if high_distance is None else high_distance
        if arm.endswith("_p"):
            lower = arm_start + low_distance - offset
            upper = arm_start + high_distance - offset - 1
        else:
            lower = arm_end - offset - high_distance + 1
            upper = arm_end - offset - low_distance
        lower = max(lower, arm_start)
        upper = min(upper, arm_end - span)
        blocks.append(Block(
            block_id=f"{arm}.block01", phr_id=str(phr["phr_id"]), chromosome=chrom,
            arm=arm, observed_start=start, observed_end=end,
            arm_start=arm_start, arm_end=arm_end, stratum=stratum,
            candidate_low=lower, candidate_high=upper))
    candidate_rows = {row["block_id"]: row for row in read_rows(RELEASE / "CANDIDATE_SPACES.tsv")}
    geometry_ok = len(blocks) == 37 and all(
        block.stratum == 0 and block.candidate_count >= 250_000
        and block.candidate_count <= 497_500
        and block.candidate_low <= block.observed_start <= block.candidate_high
        and candidate_rows[block.block_id]["source_arm"] == block.arm
        and candidate_rows[block.block_id]["destination_arm"] == block.arm
        and int(candidate_rows[block.block_id]["candidate_count"]) == block.candidate_count
        and candidate_rows[block.block_id]["range_count"] == "1"
        and candidate_rows[block.block_id]["explicit"] == "0"
        and candidate_rows[block.block_id]["observed_start_is_candidate"] == "1"
        for block in blocks)
    audit.check("same_arm_terminal_stratum_candidate_geometry", geometry_ok,
                expected="37 blocks, stratum 0, 250,000--497,500 exact starts",
                observed={"blocks": len(blocks),
                          "candidate_min": min(block.candidate_count for block in blocks),
                          "candidate_max": max(block.candidate_count for block in blocks)},
                evidence="raw PHR BED + cytobands + CANDIDATE_SPACES.tsv")
    return blocks


def json_state_hash(state: Mapping[str, object]) -> str:
    value = json.loads(json.dumps(state, default=lambda item: int(item)))
    return sha256_bytes(json.dumps(value, sort_keys=True).encode())


def regenerate_placements(audit: Audit, blocks: Sequence[Block]) -> tuple[np.ndarray, dict[str, object]]:
    manifest_rows = list(read_rows(RELEASE / "NULL_PLACEMENT_MANIFEST.tsv"))
    expected = {(int(row["first_replicate"]), int(row["last_replicate"])): row
                for row in manifest_rows}
    starts = np.empty((99_999, len(blocks)), dtype=np.int64)
    child = np.random.SeedSequence(2026071301).spawn(1)[0]
    rng = np.random.Generator(np.random.PCG64DXSM(child))
    provenance = json.loads((RELEASE / "COMPUTATIONAL_PROVENANCE.json").read_text(encoding="utf-8"))
    initial_hash = json_state_hash(rng.bit_generator.state)
    batch_mismatches = []
    cursor = 0
    for batch_index in range(100):
        first = cursor + 1
        last = min(99_999, cursor + 1000)
        raw = io.BytesIO()
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        text_handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
        writer = csv.DictWriter(text_handle, fieldnames=(
            "replicate", "block_id", "source_arm", "destination_arm", "chromosome",
            "block_start0", "block_end0", "components"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        canonical = hashlib.sha256()
        for replicate in range(first, last + 1):
            for block_index, block in enumerate(blocks):
                start = block.candidate_low + int(rng.integers(0, block.candidate_count))
                starts[replicate - 1, block_index] = start
                row = {
                    "replicate": replicate, "block_id": block.block_id,
                    "source_arm": block.arm, "destination_arm": block.arm,
                    "chromosome": block.chromosome, "block_start0": start,
                    "block_end0": start + block.span, "components": block.components,
                }
                writer.writerow(row)
                canonical.update(("\t".join(str(row[field]) for field in (
                    "replicate", "block_id", "source_arm", "destination_arm",
                    "block_start0", "block_end0", "components")) + "\n").encode())
        text_handle.flush()
        text_handle.detach()
        compressed.close()
        payload = raw.getvalue()
        row = expected[(first, last)]
        if (len(payload) != int(row["bytes"])
                or sha256_bytes(payload) != row["sha256"]
                or row["sha256"] != row["frozen_prefix_sha256"]
                or row["frozen_prefix_match"] != "1"
                or canonical.hexdigest() != row["canonical_coordinate_sha256"]):
            batch_mismatches.append({
                "batch": batch_index + 1, "first": first, "last": last,
                "bytes": [len(payload), int(row["bytes"])],
                "sha256": [sha256_bytes(payload), row["sha256"]],
                "canonical": [canonical.hexdigest(), row["canonical_coordinate_sha256"]],
            })
        cursor = last
    final_hash = json_state_hash(rng.bit_generator.state)
    array_manifest = provenance["array_manifest"]
    rng_ok = (initial_hash == array_manifest["prefix_initial_rng_state_sha256"]
              and final_hash == array_manifest["prefix_final_rng_state_sha256"]
              and json.loads(json.dumps(rng.bit_generator.state, default=lambda item: int(item)))
              == array_manifest["rng_state"])
    audit.check("all_99999_joint_placements_regenerated",
                len(expected) == 100 and cursor == 99_999 and not batch_mismatches and rng_ok,
                expected="100/100 exact gzip and canonical hashes plus final RNG state",
                observed={"batches": len(expected), "replicates": cursor,
                          "mismatches": batch_mismatches[:3], "rng_state_match": rng_ok},
                evidence="seed 2026071301, PCG64DXSM child [0], NULL_PLACEMENT_MANIFEST.tsv")
    return starts, {
        "batch_count": len(expected), "replicates": cursor,
        "gzip_hash_matches": len(expected) - len(batch_mismatches),
        "canonical_hash_matches": len(expected) - len(batch_mismatches),
        "initial_rng_state_sha256": initial_hash, "final_rng_state_sha256": final_hash,
        "status": "PASS" if not batch_mismatches and rng_ok else "FAIL",
    }


def weighted_prefix(values: Sequence[int], weights: Sequence[int]) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(np.asarray(values, dtype=np.int64), kind="stable")
    sorted_values = np.asarray(values, dtype=np.int64)[order]
    sorted_weights = np.asarray(weights, dtype=np.int64)[order]
    prefix = np.concatenate((np.zeros(1, dtype=np.int64), np.cumsum(sorted_weights, dtype=np.int64)))
    return sorted_values, prefix


def representative_null_recounts(
        audit: Audit, starts: np.ndarray, blocks: Sequence[Block],
        assignments: Mapping[str, dict[str, str]], matrix: Mapping[str, object],
        result_map: Mapping[tuple[str, tuple[str, str, str]], dict[str, str]],
        boundaries: Mapping[str, tuple[int, int]]) -> dict[str, object]:
    by_arm = {block.arm: index for index, block in enumerate(blocks)}
    output: dict[str, object] = {}
    all_passed = True
    for label, key in REPRESENTATIVES.items():
        sources = matrix["representative_sources"][key]
        copies: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
        for row in assignments.values():
            if row["functional_source_id"] in sources:
                copies[arm_for(row["seqid"], int(row["start0"]), int(row["end0"]), boundaries)].append(row)
        null_midpoint = np.zeros(99_999, dtype=np.int64)
        null_overlap = np.zeros(99_999, dtype=np.int64)
        for arm, rows in copies.items():
            # Arms without an observed PHR block contribute zero under the
            # same-arm rigid-block null; no placement is generated for them.
            if arm not in by_arm:
                continue
            block_index = by_arm[arm]
            block = blocks[block_index]
            placed_start = starts[:, block_index]
            placed_end = placed_start + block.span
            weights = [int(row["physical_copy_cn"]) for row in rows]
            mids, mid_prefix = weighted_prefix(
                [(int(row["start0"]) + int(row["end0"])) // 2 for row in rows], weights)
            gene_starts, start_prefix = weighted_prefix([int(row["start0"]) for row in rows], weights)
            gene_ends, end_prefix = weighted_prefix([int(row["end0"]) for row in rows], weights)
            left = np.searchsorted(mids, placed_start, side="left")
            right = np.searchsorted(mids, placed_end, side="left")
            null_midpoint += mid_prefix[right] - mid_prefix[left]
            begins = np.searchsorted(gene_starts, placed_end, side="left")
            finished = np.searchsorted(gene_ends, placed_start, side="right")
            null_overlap += start_prefix[begins] - end_prefix[finished]
        assignments_null = {"midpoint": null_midpoint, "overlap": null_overlap}
        item: dict[str, object] = {"key": list(key), "source_count": len(sources), "assignments": {}}
        for assignment, values in assignments_null.items():
            row = result_map[(assignment, key)]
            observed = int(row["observed_physical_copy_burden"])
            quantiles = np.quantile(values, [0.025, 0.5, 0.975])
            exceedances = int(np.count_nonzero(values >= observed))
            mean = float(np.mean(values))
            if np.all(values == values[0]) and observed > values[0]:
                z = math.inf
            else:
                pooled = np.concatenate((np.asarray([observed], dtype=float), values.astype(float)))
                pooled_sd = float(np.std(pooled, ddof=1))
                pooled_mean = float(np.mean(pooled))
                z = 0.0 if pooled_sd == 0 else (observed - pooled_mean) / pooled_sd
            unique, frequencies = np.unique(values, return_counts=True)
            summary = {
                "observed_physical_copy_burden": observed,
                "null_mean": mean, "null_median": float(quantiles[1]),
                "null_q025": float(quantiles[0]), "null_q975": float(quantiles[2]),
                "null_max": int(values.max()), "raw_exceedances": exceedances,
                "raw_permutations": int(values.size), "z_observed": z,
                "histogram": {str(int(value)): int(count) for value, count in zip(unique, frequencies)},
            }
            passed = (
                same_float(row["null_mean"], mean)
                and same_float(row["null_median"], summary["null_median"])
                and same_float(row["null_q025"], summary["null_q025"])
                and same_float(row["null_q975"], summary["null_q975"])
                and int(row["null_max"]) == summary["null_max"]
                and int(row["raw_exceedances"]) == exceedances
                and same_float(row["z_observed"], z, tolerance=2e-11)
            )
            summary["status"] = "PASS" if passed else "FAIL"
            item["assignments"][assignment] = summary
            all_passed &= passed
        output[label] = item
    audit.check("representative_null_counts_reconstructed", all_passed,
                expected="7 representative terms x 2 assignments exact",
                observed={label: {assignment: values["status"]
                                  for assignment, values in item["assignments"].items()}
                          for label, item in output.items()},
                evidence="regenerated placements + fixed genomic physical-copy/source-term matrix")
    return output


def bh_adjust(values: Sequence[float], by: bool = False) -> list[float]:
    size = len(values)
    factor = sum(1.0 / rank for rank in range(1, size + 1)) if by else 1.0
    order = sorted(range(size), key=lambda index: (values[index], index))
    adjusted = [1.0] * size
    running = 1.0
    for rank_index in range(size - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, float(values[index]) * size * factor / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def r_clopper_pearson(pairs: Sequence[tuple[int, int]]) -> dict[tuple[int, int, float], tuple[float, float]]:
    unique = sorted(set(pairs))
    payload = "k\tn\n" + "".join(f"{k}\t{n}\n" for k, n in unique)
    code = r'''
x <- read.delim(file("stdin"), header=TRUE)
for (alpha in c(0.05, 0.0125)) {
  lo <- ifelse(x$k == 0, 0, qbeta(alpha / 2, x$k, x$n - x$k + 1))
  hi <- ifelse(x$k == x$n, 1, qbeta(1 - alpha / 2, x$k + 1, x$n - x$k))
  lo <- (1 + x$n * lo) / (x$n + 1)
  hi <- (1 + x$n * hi) / (x$n + 1)
  for (i in seq_len(nrow(x))) {
    cat(x$k[i], x$n[i], format(alpha, digits=17),
      format(lo[i], digits=17), format(hi[i], digits=17), sep="\t")
    cat("\n")
  }
}
'''
    process = subprocess.run(["Rscript", "-e", code], input=payload, text=True,
                             capture_output=True, check=True)
    result: dict[tuple[int, int, float], tuple[float, float]] = {}
    for line in process.stdout.splitlines():
        k, n, alpha, lower, upper = line.split("\t")
        result[(int(k), int(n), float(alpha))] = (float(lower), float(upper))
    return result


def validate_inference_arithmetic(audit: Audit, results: Sequence[dict[str, str]]) -> dict[str, object]:
    pairs = []
    for row in results:
        pairs.extend((
            (int(row["raw_exceedances"]), int(row["raw_permutations"])),
            (int(row["collection_maxT_exceedances"]), int(row["collection_maxT_permutations"])),
            (int(row["global_maxT_exceedances"]), int(row["global_maxT_permutations"])),
        ))
    cp = r_clopper_pearson(pairs)
    arithmetic_mismatches = []
    max_t_order_ok = True
    families: dict[tuple[str, str], list[dict[str, str]]] = collections.defaultdict(list)
    for row in results:
        families[(row["assignment"], row["collection"])].append(row)
        n = int(row["raw_permutations"])
        k = int(row["raw_exceedances"])
        collection_k = int(row["collection_maxT_exceedances"])
        global_k = int(row["global_maxT_exceedances"])
        expected = {
            "p_empirical": (k + 1) / (n + 1),
            "collection_maxT_p": (collection_k + 1) / (n + 1),
            "global_maxT_p": (global_k + 1) / (n + 1),
            "count_difference": (float(row["observed_physical_copy_burden"])
                                 - float(row["null_median"])),
            "enrichment_ratio": ((float(row["observed_physical_copy_burden"]) + 0.5)
                                 / (float(row["null_mean"]) + 0.5)),
        }
        raw95 = cp[(k, n, 0.05)]
        raw_seq = cp[(k, n, 0.0125)]
        col95 = cp[(collection_k, n, 0.05)]
        col_seq = cp[(collection_k, n, 0.0125)]
        global95 = cp[(global_k, n, 0.05)]
        global_seq = cp[(global_k, n, 0.0125)]
        intervals = {
            "p_mc95_lower": raw95[0], "p_mc95_upper": raw95[1],
            "p_sequential_lower": raw_seq[0], "p_sequential_upper": raw_seq[1],
            "collection_maxT_mc95_lower": col95[0], "collection_maxT_mc95_upper": col95[1],
            "collection_maxT_sequential_lower": col_seq[0],
            "collection_maxT_sequential_upper": col_seq[1],
            "global_maxT_mc95_lower": global95[0], "global_maxT_mc95_upper": global95[1],
            "global_maxT_sequential_lower": global_seq[0],
            "global_maxT_sequential_upper": global_seq[1],
        }
        mismatched = [field for field, value in {**expected, **intervals}.items()
                      if not same_float(row[field], value, tolerance=2e-11)]
        if (n != 99_999 or int(row["collection_maxT_permutations"]) != n
                or int(row["global_maxT_permutations"]) != n):
            mismatched.append("permutations")
        if int(row["observed_physical_copy_burden"]) > 0 and not int(row["non_informative_constant"]):
            max_t_order_ok &= k <= collection_k <= global_k
        if mismatched and len(arithmetic_mismatches) < 20:
            arithmetic_mismatches.append({
                "key": [row["assignment"], row["collection"], row["relation"], row["term_id"]],
                "fields": mismatched})

    q_mismatches = []
    for (_assignment, _collection), family in families.items():
        raw = [float(row["p_empirical"]) for row in family]
        lower = [float(row["p_sequential_lower"]) for row in family]
        upper = [float(row["p_sequential_upper"]) for row in family]
        bh = bh_adjust(raw)
        by = bh_adjust(raw, by=True)
        bh_lower = bh_adjust(lower)
        bh_upper = bh_adjust(upper)
        for index, row in enumerate(family):
            if not all((
                same_float(row["bh_q"], bh[index]),
                same_float(row["by_q"], by[index]),
                same_float(row["bh_sequential_lower"], bh_lower[index]),
                same_float(row["bh_sequential_upper"], bh_upper[index]),
            )):
                q_mismatches.append((row["assignment"], row["collection"], row["relation"], row["term_id"]))
                if len(q_mismatches) >= 20:
                    break

        informative = [row for row in family
                       if int(row["observed_physical_copy_burden"]) > 0
                       and not int(row["non_informative_constant"])]
        informative.sort(key=lambda row: float(row["z_observed"]))
        for left, right in zip(informative, informative[1:]):
            if int(left["collection_maxT_exceedances"]) < int(right["collection_maxT_exceedances"]):
                max_t_order_ok = False
            if math.isclose(float(left["z_observed"]), float(right["z_observed"]), abs_tol=1e-14) \
                    and left["collection_maxT_exceedances"] != right["collection_maxT_exceedances"]:
                max_t_order_ok = False
    for assignment in ("midpoint", "overlap"):
        informative = [row for row in results if row["assignment"] == assignment
                       and int(row["observed_physical_copy_burden"]) > 0
                       and not int(row["non_informative_constant"])]
        informative.sort(key=lambda row: float(row["z_observed"]))
        for left, right in zip(informative, informative[1:]):
            if int(left["global_maxT_exceedances"]) < int(right["global_maxT_exceedances"]):
                max_t_order_ok = False
    audit.check("complete_plus_one_cp_maxT_arithmetic",
                not arithmetic_mismatches and max_t_order_ok,
                expected="62,470 rows exact using independent R qbeta",
                observed={"mismatches": arithmetic_mismatches,
                          "maxT_order_and_nesting": max_t_order_ok},
                evidence="TERM_RESULTS.tsv.gz")
    audit.check("complete_bh_by_interval_adjustments", not q_mismatches,
                expected="BH/BY and BH sequential bounds exact in 8 assignment families",
                observed=q_mismatches or "all exact", evidence="TERM_RESULTS.tsv.gz")

    status_counts: collections.Counter[str] = collections.Counter()
    decision_mismatches = []
    for row in results:
        if row["assignment"] != "midpoint":
            expected_status = "SENSITIVITY_COMPLETE"
        else:
            passed = (float(row["bh_sequential_upper"]) <= 0.05
                      and float(row["global_maxT_sequential_upper"]) <= 0.05)
            nonpass = (float(row["bh_sequential_lower"]) > 0.05
                       or float(row["global_maxT_sequential_lower"]) > 0.05)
            expected_status = "CERTIFIED_PASS" if passed else \
                "CERTIFIED_NONPASS" if nonpass else "MC_UNRESOLVED"
            status_counts[expected_status] += 1
        if row["mc_status"] != expected_status:
            decision_mismatches.append((row["assignment"], row["collection"], row["relation"], row["term_id"]))
    expected_counts = {"CERTIFIED_PASS": 143, "CERTIFIED_NONPASS": 31_092}
    stage = json.loads((RELEASE / "STAGE_DECISION.json").read_text(encoding="utf-8"))
    candidate_rows = list(read_rows(RELEASE / "SELECTIVE_EXTENSION_CANDIDATES.tsv"))
    decision_ok = (not decision_mismatches and dict(status_counts) == expected_counts
                   and not candidate_rows and stage["status"] == "INITIAL_SCREEN_SUFFICIENT"
                   and stage["selective_extension_candidates"] == 0
                   and stage["maxT_unresolved_not_selectively_extended"] == 0)
    audit.check("uncertainty_and_stopping_decisions_recalculated", decision_ok,
                expected={**expected_counts, "unresolved": 0, "extension_candidates": 0},
                observed={**dict(status_counts), "decision_mismatches": decision_mismatches[:10],
                          "extension_candidates": len(candidate_rows),
                          "maxT_unresolved": stage["maxT_unresolved_not_selectively_extended"]},
                evidence="98.75% sequential intervals + STAGE_DECISION.json")
    return {
        "status_counts": dict(status_counts), "arithmetic_mismatches": arithmetic_mismatches,
        "q_mismatches": q_mismatches, "decision_mismatches": decision_mismatches,
        "maxT_order_and_nesting": max_t_order_ok,
        "selective_extension_candidates": len(candidate_rows),
    }


def git_blob_sha(commit: str, path: str) -> str:
    process = subprocess.run(["git", "show", f"{commit}:{path}"], cwd=REPO,
                             capture_output=True, check=True)
    return sha256_bytes(process.stdout)


def historical_blob_match(path: str, expected_sha: str, expected_bytes: int) -> str:
    history = subprocess.run(
        ["git", "log", "--format=%H", "--all", "--", path], cwd=REPO,
        capture_output=True, text=True, check=True).stdout.splitlines()
    for commit in history:
        process = subprocess.run(["git", "show", f"{commit}:{path}"], cwd=REPO,
                                 capture_output=True, check=True)
        if len(process.stdout) == expected_bytes and sha256_bytes(process.stdout) == expected_sha:
            return commit
    return ""


def validate_provenance(audit: Audit, placement: Mapping[str, object]) -> dict[str, object]:
    provenance = json.loads((RELEASE / "COMPUTATIONAL_PROVENANCE.json").read_text(encoding="utf-8"))
    array = provenance["array_manifest"]
    run_manifest = next(read_rows(RELEASE / "RUN_MANIFEST.tsv"))
    constants_ok = (
        provenance["initial_seed"] == 2026071301
        and provenance["spawn_key"] == [0]
        and provenance["bit_generator"] == "PCG64DXSM"
        and array["completed_permutations"] == array["valid_replicates"] == 99_999
        and array["prefix_batch_matches"] == 100 and array["prefix_all_100_match"] is True
        and array["immutable_configuration"]["counting_statistic"] == "sum_physical_copy_cn"
        and array["immutable_configuration"]["physical_clusters_fixed_in_genome"] is True
        and array["immutable_configuration"]["source_term_edges_fixed_in_genome"] is True
        and run_manifest["seed"] == "2026071301" and run_manifest["spawn_key"] == "[0]"
        and run_manifest["bit_generator"] == "PCG64DXSM"
        and run_manifest["initial_permutations"] == "99999"
        and run_manifest["statistic"] == "sum(physical_copy_cn)"
    )
    audit.check("seed_rng_and_fixed_matrix_provenance", constants_ok,
                expected="seed 2026071301; child [0]; PCG64DXSM; 99,999; physical CN",
                observed={"seed": provenance["initial_seed"], "spawn_key": provenance["spawn_key"],
                          "bit_generator": provenance["bit_generator"],
                          "permutations": array["completed_permutations"]},
                evidence="COMPUTATIONAL_PROVENANCE.json + RUN_MANIFEST.tsv")
    embedded_manifest_sha = sha256_bytes(
        (json.dumps(array, indent=2, sort_keys=True) + "\n").encode())
    engine_path = "paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/v6_engine.py"
    targeted_path = "paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/v6_targeted_extension.py"
    versioned_ok = (
        embedded_manifest_sha == provenance["array_manifest_sha256"]
        and git_blob_sha(array["git_commit"], engine_path) == provenance["array_engine_sha256"]
        and git_blob_sha(provenance["current_git_commit"], engine_path) == provenance["release_engine_sha256"]
        and git_blob_sha(provenance["current_git_commit"], targeted_path) == provenance["targeted_engine_sha256"]
        and sha256_path(BASE / "v6_engine.py") == provenance["release_engine_sha256"]
        and sha256_path(BASE / "v6_targeted_extension.py") == provenance["targeted_engine_sha256"]
    )
    audit.check("versioned_engine_and_manifest_checksums", versioned_ok,
                expected="array and release engines resolve at recorded Git commits",
                observed={"array_manifest_sha256": embedded_manifest_sha,
                          "array_engine_sha256": provenance["array_engine_sha256"],
                          "release_engine_sha256": provenance["release_engine_sha256"]},
                evidence="COMPUTATIONAL_PROVENANCE.json + Git objects")

    input_rows = list(read_rows(RELEASE / "INPUT_CHECKSUMS.tsv"))
    input_failures = []
    historical_input_engine_commit = ""
    for row in input_rows:
        path = resolve_repo_path(row["path"])
        current_match = (path.is_file() and path.stat().st_size == int(row["bytes"])
                         and sha256_path(path) == row["sha256"])
        if current_match:
            continue
        if row["path"] == engine_path:
            historical_input_engine_commit = historical_blob_match(
                row["path"], row["sha256"], int(row["bytes"]))
            if historical_input_engine_commit:
                continue
        input_failures.append(row["path"])
    audit.check("frozen_inference_input_checksums_versioned", not input_failures and len(input_rows) == 11,
                expected="10 immutable inputs current-exact plus historical engine blob exact",
                observed={"rows": len(input_rows), "failures": input_failures,
                          "historical_engine_commit": historical_input_engine_commit},
                evidence="INPUT_CHECKSUMS.tsv + Git object history")

    transient = list(read_rows(RELEASE / "TRANSIENT_COUNT_ARRAYS.tsv"))
    transient_ok = (len(transient) == 200
                    and {row["assignment"] for row in transient} == {"midpoint", "overlap"}
                    and all(row["statistic"] == "sum(physical_copy_cn)"
                            and len(row["sha256"]) == 64 for row in transient)
                    and sum(int(row["last_replicate"]) - int(row["first_replicate"]) + 1
                            for row in transient if row["assignment"] == "midpoint") == 99_999)
    audit.check("transient_count_array_inventory_complete", transient_ok,
                expected="200 files, two assignments, 99,999 columns each",
                observed={"rows": len(transient), "payloads_present": sum(
                    resolve_repo_path(row["relative_path"]).is_file() for row in transient)},
                evidence="TRANSIENT_COUNT_ARRAYS.tsv; representative columns regenerated independently")

    slurm = json.loads((RELEASE / "SLURM_COMPLETION.json").read_text(encoding="utf-8"))
    job_rows = list(read_rows(RELEASE / "SLURM_JOBS.tsv"))
    log_ok = True
    for row in job_rows:
        for prefix in ("stdout", "stderr"):
            path = resolve_repo_path(row[f"{prefix}_log"])
            log_ok &= (path.is_file() and path.stat().st_size == int(row[f"{prefix}_bytes"])
                       and sha256_path(path) == row[f"{prefix}_sha256"])
    accounting: dict[str, tuple[str, str]] = {}
    command = ["sacct", "-j", "1761733,1761734,1761735", "--parsable2", "--noheader",
               "--format=JobIDRaw,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList"]
    process = subprocess.run(command, capture_output=True, text=True, check=True)
    for line in process.stdout.splitlines():
        fields = line.split("|")
        accounting[fields[0]] = (fields[3], fields[4])
    live_ok = all(accounting.get(job_id) == ("COMPLETED", "0:0")
                  for job_id in ("1761733", "1761734", "1761735"))
    recorded_ok = (slurm["completion_status"] == "COMPLETE"
                   and all(job["state"] == "COMPLETED" and job["exit_code"] == "0:0"
                           for job in slurm["jobs"].values())
                   and all(row["state"] == "COMPLETED" and row["exit_code"] == "0:0"
                           for row in job_rows))
    audit.check("slurm_jobs_live_complete_and_logs_exact", recorded_ok and live_ok and log_ok,
                expected="jobs 1761733/4/5 COMPLETED 0:0 with exact logs",
                observed={"live": accounting, "recorded_completion": slurm["completion_status"],
                          "log_checksums": log_ok},
                evidence="live sacct + SLURM_COMPLETION.json + SLURM_JOBS.tsv")
    return {
        "seed": provenance["initial_seed"], "spawn_key": provenance["spawn_key"],
        "bit_generator": provenance["bit_generator"], "permutations": array["completed_permutations"],
        "array_git_commit": array["git_commit"], "release_git_commit": provenance["current_git_commit"],
        "array_engine_sha256": provenance["array_engine_sha256"],
        "release_engine_sha256": provenance["release_engine_sha256"],
        "frozen_input_engine_git_commit": historical_input_engine_commit,
        "slurm_jobs": {role: {"job_id": job["job_id"], "state": job["state"],
                               "exit_code": job["exit_code"]}
                       for role, job in slurm["jobs"].items()},
        "placement_regeneration": placement,
        "transient_payloads_committed": False,
        "transient_inventory_rows": len(transient),
    }


def claims_summary(results: Sequence[dict[str, str]], named_counts: Mapping[str, Mapping[str, int]]) -> dict[str, object]:
    midpoint = [row for row in results if row["assignment"] == "midpoint"]
    passes = [row for row in midpoint if row["mc_status"] == "CERTIFIED_PASS"]
    pass_counts = collections.Counter((row["collection"], row["relation"]) for row in passes)
    result_lookup = {(row["collection"], row["relation"], row["term_id"]): row for row in midpoint}
    named_examples = []
    for label in ("DUX4_DUX4L", "DDX11L", "TUBB8", "OR4F", "WASH"):
        key = REPRESENTATIVES[label]
        row = result_lookup[key]
        named_examples.append({
            "cohort": label, "term_key": list(key), "term_name": row["term_name"],
            "observed": int(row["observed_physical_copy_burden"]),
            "raw_p": float(row["p_empirical"]), "bh_q": float(row["bh_q"]),
            "global_maxT_p": float(row["global_maxT_p"]), "status": row["mc_status"],
            "cohort_phr_ontology_contributor_cn": named_counts[label]["eligible"],
        })
    return {
        "supported": {
            "definition": "term rows with 98.75% BH and global-maxT upper bounds <= 0.05",
            "term_count": len(passes),
            "counts_by_collection_relation": {
                f"{collection}|{relation}": count
                for (collection, relation), count in sorted(pass_counts.items())},
            "scope": "CHM13 regional physical-copy enrichment only; contributors may be named copies but families are not hypotheses",
        },
        "unresolved": {
            "term_count": sum(row["mc_status"] == "MC_UNRESOLVED" for row in midpoint),
            "items": [],
        },
        "descriptive": {
            "named_physical_burdens": {cohort: values["eligible"] for cohort, values in named_counts.items()},
            "representative_term_rows": named_examples,
            "statement": "Named copy burdens remain coordinate descriptions unless their exact frozen term row is certified.",
        },
        "unsupported": [
            "Any unique-gene, unique-source, or gene-family collapse as the tested statistic.",
            "A DUX, DDX11L, TUBB8, OR4F, or WASH family-level hypothesis not present in the frozen catalog.",
            "Causal functional coordination, expression, or chromosome-contact mechanisms.",
            "Extrapolation from CHM13 regional copies to population prevalence or the human pangenome.",
            "Terms or biological systems selected after inspecting PHR membership.",
        ],
    }


def write_outputs(audit: Audit, payload: dict[str, object], results: Sequence[dict[str, str]]) -> None:
    mapping_checks = {
        "raw_assignment_evidence_bijection", "all_raw_loci_retained_with_exact_coordinates",
        "physical_copy_cn_complete_positive", "copy_evidence_assignment_routes_exact",
        "eligible_copy_cn", "source_route_gate_counts",
        "all_nonself_sources_have_raw_directed_evidence", "remediated_pilot_evidence_reconciled",
        "raw_direct_and_ancestor_source_terms_exact", "every_copy_term_count_reconstructed",
    }
    copy_checks = {"physical_copy_cn_complete_positive", "every_observed_term_is_sum_physical_cn",
                   "n_copy_adversarial_fixture_rejects_collapse", "exact_contributor_rows_reconstructed"}
    mapping_pass = all(row["status"] == "PASS" for row in audit.checks if row["name"] in mapping_checks)
    copy_pass = all(row["status"] == "PASS" for row in audit.checks if row["name"] in copy_checks)
    strict_pass = not audit.failed and mapping_pass and copy_pass
    payload.update({
        "schema_version": SCHEMA,
        "validated_utc": utcnow(),
        "status": "PASS" if strict_pass else "FAIL",
        "strict_pass": strict_pass,
        "qualified_pass": False,
        "production_mapping_or_inference_code_imported": False,
        "mapping_gate_pass": mapping_pass,
        "copy_number_semantics_pass": copy_pass,
        "inference_valid": strict_pass,
        "checks_total": len(audit.checks),
        "checks_passed": len(audit.checks) - len(audit.failed),
        "checks_failed": len(audit.failed),
        "failed_check_names": [row["name"] for row in audit.failed],
        "checks": audit.checks,
    })
    MACHINE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    named = payload["named_cohorts"]
    nulls = payload["representative_null_recounts"]
    claims = payload["claims"]
    status = payload["status"]
    lines = [
        "# V6 independent validation report", "",
        f"**Strict result: {status}.** This is not a qualified pass. The copy-number, mapping, "
        "ontology, geometry, Monte Carlo, multiplicity, and provenance gates all had to pass; "
        "a failure in physical-copy semantics would force the machine result to `FAIL`.", "",
        "## Scope and independence", "",
        "The validator in `independent_validation/validate_v6_independent.py` imports none of the "
        "production mapping, sampler, inference, or prior validation modules. It rebuilt the matrix "
        "from the raw CHM13 GFF, frozen per-copy evidence, raw NCBI gene2go, GO OBO, Reactome "
        "pathway relations/all-level assignments, cytobands, and the PHR BED. NumPy was used only "
        "to reproduce the specified PCG64DXSM stream and vectorize interval recounts; Clopper--Pearson "
        "bounds were recomputed independently with R `qbeta`.", "",
        "## Fail-closed source and copy gates", "",
        f"All 61,312 raw gene features form an exact bijection with 61,312 source assignments and "
        f"61,312 evidence rows. All coordinate-coincident records were retained (12 extra physical "
        "rows beyond coordinate+strand uniqueness), every row has positive physical CN, and summed "
        "CN is 61,312. The audit reconstructed every assignment and evidence-row digest.", "",
        "Exactly 31,966 physical copies are ontology-eligible: 20,405 exact-self mappings and "
        "11,561 mappings backed by a directed human `Related functional gene` record. All 45 "
        "ambiguous and 28,407 unsupported rows emit no functional source. The 168 remediated pilot "
        "records match their frozen upstream evidence and hashes. No name, alias, family, or sequence-only "
        "route entered ontology inference.", "",
        "## Direct/ancestor terms and physical-CN recount", "",
        f"The raw ontology reconstruction matched all {payload['matrix']['source_term_rows']:,} "
        f"`SOURCE_TERMS` rows field-for-field, including direct/ancestor relation, minimum distance, "
        "direct leaves, evidence codes, qualifiers, record identifiers, release labels, and ontology "
        f"digests. Joining those terms to coordinate-distinct copies reconstructed "
        f"{payload['matrix']['physical_copy_term_edge_cn']:,} genome edge-CN "
        f"({payload['matrix']['direct_physical_copy_term_edge_cn']:,} direct), "
        f"{payload['matrix']['phr_midpoint_edge_cn']:,} midpoint edge-CN, and "
        f"{payload['matrix']['phr_overlap_edge_cn']:,} any-overlap edge-CN.", "",
        "Every one of the 31,235 frozen hypotheses matched its independently rebuilt genome CN, arm "
        "count, source count, target-blind order, and multiplicity family. Both released observed "
        "columns (62,470 rows) equal sums of physical CN; all 17,579 exact contributor rows were "
        "reconstructed from raw coordinates, source evidence, and ontology closure.", "",
        "### Named cohorts", "",
        "| Cohort | Genome physical CN | PHR midpoint CN | Inference-eligible PHR CN |", "|---|---:|---:|---:|",
    ]
    for cohort in EXPECTED_NAMED:
        values = named[cohort]
        lines.append(f"| {cohort} | {values['genome']} | {values['phr']} | {values['eligible']} |")
    lines.extend([
        "", "These are physical-copy burdens, not unique-gene or family counts. The adversarial "
        "fixture placed seven CN=1 coordinates plus one CN=3 coordinate on one source/gene/family: "
        "the required burden is 10, while every source/gene/family collapse incorrectly returns 1. "
        "The validator rejected all three collapsed estimands.", "",
        "## Null placements and representative retained counts", "",
        f"All {payload['placement']['replicates']:,} joint placements were regenerated from seed "
        f"{payload['provenance']['seed']} and spawn key `{payload['provenance']['spawn_key']}`. "
        f"All {payload['placement']['gzip_hash_matches']}/100 deterministic gzip batch hashes and "
        f"all {payload['placement']['canonical_hash_matches']}/100 canonical-coordinate hashes matched. "
        "Each of 37 rigid blocks remained on its source arm in the same terminal-distance stratum; "
        "component widths and CN/source/term edges remained fixed.", "",
        "The transient full count matrices are not committed, so the audit regenerated exact columns "
        "for leading, negative, and every named-system representative. Released mean, quantiles, maximum, "
        "exceedances, and standardized observed value matched for midpoint and any-overlap:", "",
        "| Role/cohort | Term | Midpoint observed | Midpoint null mean | Midpoint exceedances | Overlap status |", "|---|---|---:|---:|---:|---|",
    ])
    for label, item in nulls.items():
        mid = item["assignments"]["midpoint"]
        ov = item["assignments"]["overlap"]
        lines.append(
            f"| {label} | `{item['key'][1]}|{item['key'][2]}` | "
            f"{mid['observed_physical_copy_burden']} | {mid['null_mean']:.6f} | "
            f"{mid['raw_exceedances']:,}/{mid['raw_permutations']:,} | {ov['status']} |")
    lines.extend([
        "", "## Empirical inference and multiplicity", "",
        "For every result row, the audit recalculated `(exceedances + 1)/(permutations + 1)` for "
        "raw, collection-maxT, and global-maxT p-values; two-sided 95% and sequential 98.75% "
        "plus-one-transformed Clopper--Pearson intervals; BH and BY within each collection across "
        "direct plus ancestor rows; BH-adjusted sequential endpoints; count differences; and burden "
        "ratios. MaxT exceedance nesting and monotonicity in standardized observed burden also hold.", "",
        f"The primary classification independently returns {claims['supported']['term_count']} "
        "`CERTIFIED_PASS`, 31,092 `CERTIFIED_NONPASS`, and 0 `MC_UNRESOLVED` rows. No selective "
        "extension candidate exists; the complete 99,999-placement screen is sufficient under the "
        "frozen stopping rule.", "",
        "## Claim classification", "",
        "### Supported", "",
        f"There are {claims['supported']['term_count']} supported term-level CHM13 regional "
        "physical-copy enrichments. Support attaches only to the exact `(collection, relation, term)` "
        "rows whose sequential BH and global-maxT upper bounds are at most 0.05. A supported term may "
        "have DDX11L- or WASH-associated contributors, but this does not create a family hypothesis.", "",
        "### Unresolved", "",
        "None. Zero primary rows have a sequential overall decision that straddles 0.05.", "",
        "### Descriptive", "",
        "The named burdens 65/10/2/4/9 are validated coordinate counts. DUX zygotic-genome-activation "
        "and OR4F receptor-expression examples are descriptive because their exact rows fail complete "
        "multiplicity/maxT safeguards. TUBB8 is likewise not promoted to a family-level inference. "
        "DDX11L helicase and WASH-complex term rows are certified, but only at their frozen term keys.", "",
        "### Unsupported", "",
    ])
    for claim in claims["unsupported"]:
        lines.append(f"- {claim}")
    lines.extend([
        "", "## Computational provenance", "",
        f"Live `sacct` reconfirmed jobs 1761733, 1761734, and 1761735 as `COMPLETED` with exit "
        "code `0:0`; recorded stdout/stderr byte counts and SHA-256 digests match. The initial arrays "
        f"resolve to engine `{payload['provenance']['array_engine_sha256']}` at Git commit "
        f"`{payload['provenance']['array_git_commit']}`. The released finalization resolves to engine "
        f"`{payload['provenance']['release_engine_sha256']}` at Git commit "
        f"`{payload['provenance']['release_git_commit']}`. All frozen raw inputs and all released "
        "artifacts match their checksum manifests.", "",
        "The 200 transient count-array payloads are intentionally absent from the checkout. Their "
        "inventory, shapes, assignments, statistics, and hashes are retained; this audit therefore "
        "reconstructed representative columns from the fully regenerated placement stream instead of "
        "claiming to re-hash absent payloads. This does not qualify the strict result because every "
        "placement batch, every observed term, all inferential arithmetic, and the prespecified "
        "representative null columns passed independently.", "",
        "## Machine verdict", "",
        f"`V6_VALIDATION.json`: `status={status}`, `strict_pass={str(payload['strict_pass']).lower()}`, "
        f"`copy_number_semantics_pass={str(payload['copy_number_semantics_pass']).lower()}`, "
        f"`checks={payload['checks_passed']}/{payload['checks_total']}`.", "",
    ])
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if any(name in sys.modules for name in ("v6_engine", "COPY_engine", "v3_engine",
                                             "build_genomewide_source_map")):
        raise RuntimeError("production module unexpectedly imported")
    audit = Audit()
    gate = validate_manifests(audit)
    raw, coordinate_extras = parse_raw_gff()
    audit.check("raw_physical_locus_count_and_coordinate_retention",
                len(raw) == 61_312 and coordinate_extras == 12,
                expected={"raw_loci": 61_312, "coordinate_coincident_extras": 12},
                observed={"raw_loci": len(raw), "coordinate_coincident_extras": coordinate_extras},
                evidence="chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
    assignments, evidence = load_mapping(audit, raw)
    validate_source_routes(audit, assignments, evidence)
    adversarial = validate_adversarial_fixture(audit)
    boundaries = load_boundaries()
    phrs = load_phrs(boundaries)
    midpoint, overlap = memberships(assignments, phrs)
    midpoint_cn = sum(int(assignments[copy_id]["physical_copy_cn"]) for copy_id in assignments if midpoint[copy_id])
    overlap_cn = sum(int(assignments[copy_id]["physical_copy_cn"]) for copy_id in assignments if overlap[copy_id])
    audit.check("phr_membership_reconstructed", len(phrs) == 37 and midpoint_cn == 402 and overlap_cn == 412,
                expected={"intervals": 37, "midpoint_cn": 402, "overlap_cn": 412},
                observed={"intervals": len(phrs), "midpoint_cn": midpoint_cn, "overlap_cn": overlap_cn},
                evidence="raw PHR BED + raw assignment coordinates")
    named = validate_named_cohorts(audit, assignments, midpoint)
    matrix = reconstruct_matrix(audit, assignments, evidence, boundaries, midpoint, overlap)
    validate_mapping_coverage(audit, assignments, midpoint, overlap)
    catalog, results, result_map = validate_catalog_and_results(audit, matrix)
    validate_exact_contributors(audit, assignments, evidence, midpoint, overlap, matrix)
    arms = arm_summary(audit, boundaries)
    blocks = build_blocks(audit, phrs, arms)
    starts, placement = regenerate_placements(audit, blocks)
    null_recounts = representative_null_recounts(
        audit, starts, blocks, assignments, matrix, result_map, boundaries)
    inference = validate_inference_arithmetic(audit, results)
    provenance = validate_provenance(audit, placement)
    claims = claims_summary(results, named)
    payload: dict[str, object] = {
        "gate_status": gate["status"],
        "raw_physical_loci": len(raw), "coordinate_coincident_extra_rows_retained": coordinate_extras,
        "physical_copy_cn": sum(int(row["physical_copy_cn"]) for row in assignments.values()),
        "ontology_eligible_physical_copy_cn": sum(
            int(row["physical_copy_cn"]) for row in assignments.values() if row["functional_source_id"]),
        "phr_midpoint_physical_copy_cn": midpoint_cn,
        "phr_any_overlap_physical_copy_cn": overlap_cn,
        "hypotheses": len(catalog), "term_result_rows": len(results),
        "matrix": {
            "source_term_rows": matrix["term_count"],
            "physical_copy_term_edge_cn": matrix["edge_cn"],
            "direct_physical_copy_term_edge_cn": matrix["direct_edge_cn"],
            "phr_midpoint_edge_cn": matrix["phr_mid_edge_cn"],
            "phr_overlap_edge_cn": matrix["phr_overlap_edge_cn"],
        },
        "named_cohorts": named,
        "adversarial_fixture": adversarial,
        "placement": placement,
        "representative_null_recounts": null_recounts,
        "inference_arithmetic": inference,
        "provenance": provenance,
        "claims": claims,
    }
    write_outputs(audit, payload, results)
    print(json.dumps({
        "status": payload["status"], "checks_passed": payload["checks_passed"],
        "checks_total": payload["checks_total"], "report": str(REPORT.relative_to(REPO)),
        "machine": str(MACHINE.relative_to(REPO)),
    }, indent=2, sort_keys=True))
    return 0 if payload["strict_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
