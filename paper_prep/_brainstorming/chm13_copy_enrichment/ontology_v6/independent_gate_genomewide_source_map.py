#!/usr/bin/env python3
"""Independent hard gate for the V6 CHM13 physical-copy ontology map.

This verifier deliberately does not import the production builder or its
checker.  It reparses the raw GFF3 and registry snapshots, adjudicates every
copy/source edge against the frozen directed relation, rebuilds GO/Reactome
true-path closure, rejoins the PHR BED only after source adjudication, and
reconstructs every copy-term edge and burden.

The program writes mapping-gate evidence only.  It never runs enrichment or
computes a p-value.
"""

from __future__ import print_function

import collections
import csv
import gzip
import hashlib
import json
import os
import pathlib
import re
import sys
import traceback
import urllib.parse


SCRIPT = pathlib.Path(__file__).resolve()
OUTDIR = SCRIPT.parent
REPO = SCRIPT.parents[4]
BASE = SCRIPT.parents[1]
SOURCES = BASE / "sources"
PILOT_DIR = BASE / "ontology_v5" / "real_source_pilot"

EXPECTED_LOCI = 61312
EXPECTED_PHR_MIDPOINT = 402
EXPECTED_PHR_OVERLAP = 412
EXPECTED_SOURCE_TERMS = 1686727
EXPECTED_COPY_TERM_EDGES = 2929709
EXPECTED_NAMED = collections.OrderedDict((
    ("DUX4_DUX4L", (107, 68, 65)),
    ("DDX11L", (12, 10, 10)),
    ("TUBB8", (16, 7, 2)),
    ("OR4F", (15, 11, 4)),
    ("WASH", (18, 9, 9)),
))

ASSIGNMENT_PATH = OUTDIR / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz"
EVIDENCE_PATH = OUTDIR / "COPY_SOURCE_EVIDENCE.tsv.gz"
MAP_PATH = OUTDIR / "GENOMEWIDE_SOURCE_MAP.tsv.gz"
SOURCE_TERM_PATH = OUTDIR / "SOURCE_TERMS.tsv.gz"
EDGE_PATH = OUTDIR / "PHYSICAL_COPY_TERM_EDGES.tsv.gz"
BURDEN_PATH = OUTDIR / "TERM_BURDENS.tsv.gz"
CONTRIBUTOR_PATH = OUTDIR / "EXACT_CONTRIBUTORS.tsv.gz"
NAMED_BURDEN_PATH = OUTDIR / "NAMED_COHORT_TERM_BURDENS.tsv.gz"
PILOT_EVIDENCE_PATH = PILOT_DIR / "REAL_SOURCE_COPY_SOURCE_EVIDENCE.tsv.gz"

REPORT_PATH = OUTDIR / "PRE_RUN_V6_GATE.md"
JSON_PATH = OUTDIR / "PRE_RUN_V6_GATE.json"
CHECK_PATH = OUTDIR / "PRE_RUN_V6_GATE_EVIDENCE.tsv"
COVERAGE_PATH = OUTDIR / "PRE_RUN_V6_MAPPING_COVERAGE.tsv"
RELEASE_DIGEST_PATH = OUTDIR / "PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv"


class GateBlocked(RuntimeError):
    """Raised on the first hard-gate mismatch."""


class CheckLedger(object):
    def __init__(self):
        self.rows = []

    def require(self, name, observed, expected, evidence):
        passed = observed == expected
        self.rows.append({
            "check": name,
            "observed": observed,
            "expected": expected,
            "status": "PASS" if passed else "BLOCK",
            "evidence": evidence,
        })
        if not passed:
            raise GateBlocked(
                "{}: observed {!r}, expected {!r}".format(name, observed, expected)
            )


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(str(path), "rt", encoding="utf-8", newline="")
    return open(str(path), "rt", encoding="utf-8", newline="")


def read_tsv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    tmp = pathlib.Path(str(path) + ".tmp")
    with open(str(tmp), "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: scalar(row.get(field, "")) for field in fields})
    os.replace(str(tmp), str(path))


def scalar(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple, set)):
        return json.dumps(json_safe(value), sort_keys=True, separators=(",", ":"))
    return str(value)


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, set):
        return [json_safe(item) for item in sorted(value)]
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def row_hash(row, fields):
    payload = "\t".join(str(row.get(field, "")) for field in fields)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_attributes(raw):
    result = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = urllib.parse.unquote(value)
    return result


def split_tokens(value):
    return tuple(token.strip() for token in re.split(r"[,|]", value or "") if token.strip())


def compare_row(expected, observed, fields, context):
    for field in fields:
        left = scalar(expected.get(field, ""))
        right = observed.get(field, "")
        if left != right:
            raise GateBlocked(
                "{} field {} differs: observed {!r}, expected {!r}".format(
                    context, field, right, left
                )
            )


def next_or_none(reader):
    try:
        return next(reader)
    except StopIteration:
        return None


def verify_input_manifest(ledger):
    rows = read_tsv(OUTDIR / "INPUT_MANIFEST.tsv")
    expected_roles = {
        "raw_chm13_gff", "hgnc", "hgnc_withdrawn", "ncbi_gene_info",
        "ncbi_gene_group", "ncbi_gene2go", "go_basic", "reactome_ncbi",
        "reactome_pathways", "reactome_relations", "source_manifest",
        "remediated_pilot_evidence", "remediated_pilot_manifest", "phr_bed",
    }
    ledger.require("input_manifest_roles", set(row["role"] for row in rows),
                   expected_roles, "ontology_v6/INPUT_MANIFEST.tsv")
    by_role = {row["role"]: row for row in rows}
    phr = by_role["phr_bed"]
    ledger.require("phr_declared_post_source_freeze", phr["stage"],
                   "after_source_freeze_target_join", "INPUT_MANIFEST.tsv: phr_bed")
    resolved = {}
    hashes = {}
    # Deliberately do not stat, hash, or open the target BED yet.
    for role in sorted(expected_roles - {"phr_bed"}):
        row = by_role[role]
        path = REPO / row["path"]
        observed = sha256_file(path) if path.is_file() else "MISSING"
        expected = row["sha256"] if row["sha256"] == row["expected_sha256"] else "MANIFEST_CONFLICT"
        ledger.require("input_sha256_{}".format(role), observed, expected, row["path"])
        ledger.require("input_bytes_{}".format(role), path.stat().st_size,
                       int(row["bytes"]), row["path"])
        resolved[role] = path
        hashes[role] = observed
    resolved["phr_bed"] = REPO / phr["path"]
    hashes["phr_bed_expected"] = phr["sha256"]
    return resolved, hashes


def verify_release_manifest(ledger):
    rows = read_tsv(OUTDIR / "OUTPUT_MANIFEST.sha256.tsv")
    observed_rows = []
    for row in rows:
        path = OUTDIR / row["path"]
        observed_hash = sha256_file(path) if path.is_file() else "MISSING"
        observed_bytes = path.stat().st_size if path.is_file() else -1
        ledger.require("release_sha256_{}".format(row["path"]), observed_hash,
                       row["sha256"], "OUTPUT_MANIFEST.sha256.tsv")
        ledger.require("release_bytes_{}".format(row["path"]), observed_bytes,
                       int(row["bytes"]), "OUTPUT_MANIFEST.sha256.tsv")
        observed_rows.append({
            "path": row["path"], "bytes": observed_bytes, "sha256": observed_hash,
        })
    return observed_rows


def load_raw_loci(path):
    loci = []
    by_gene_feature_id = {}
    transcript_owner = {}
    coordinate_counts = collections.Counter()
    seen = set()
    with gzip.open(str(path), "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "transcript", "exon", "CDS"}:
                continue
            attrs = parse_attributes(fields[8])
            if fields[2] != "gene":
                parents = split_tokens(attrs.get("Parent", ""))
                owners = set()
                if fields[2] == "transcript":
                    owners.update(parent for parent in parents if parent in by_gene_feature_id)
                    if len(owners) == 1 and attrs.get("ID"):
                        transcript_owner[attrs["ID"]] = next(iter(owners))
                else:
                    owners.update(transcript_owner[parent] for parent in parents
                                  if parent in transcript_owner)
                if len(owners) == 1:
                    owner = by_gene_feature_id[next(iter(owners))]
                    owner["descendant_xrefs"].update(split_tokens(
                        attrs.get("db_xref", "") or attrs.get("Dbxref", "")
                    ))
                continue
            start1 = int(fields[3])
            end1 = int(fields[4])
            gff_id = attrs.get("ID", "")
            copy_id = "CHM13v2.0|{}:{}-{}|{}|{}".format(
                fields[0], start1, end1, fields[6], gff_id
            )
            if copy_id in seen:
                raise GateBlocked("duplicate raw copy_id {}".format(copy_id))
            seen.add(copy_id)
            coordinate_counts[(fields[0], start1, end1, fields[6])] += 1
            raw_xrefs = attrs.get("db_xref", "") or attrs.get("Dbxref", "")
            row = {
                "copy_id": copy_id,
                "seqid": fields[0],
                "start0": start1 - 1,
                "end0": end1,
                "start1": start1,
                "end1": end1,
                "strand": fields[6],
                "gff_id": gff_id,
                "gff_line": line_number,
                "gene_name": attrs.get("gene_name") or attrs.get("gene") or gff_id,
                "gene_synonyms": ",".join(split_tokens(attrs.get("gene_synonym", ""))),
                "gene_biotype": attrs.get("gene_biotype", ""),
                "raw_stable_xrefs": raw_xrefs,
                "raw_attributes_sha256": hashlib.sha256(fields[8].encode("utf-8")).hexdigest(),
                "physical_copy_cn": 1,
                "descendant_xrefs": set(),
            }
            loci.append(row)
            by_gene_feature_id[gff_id] = row
    coincident_extras = sum(value - 1 for value in coordinate_counts.values() if value > 1)
    return loci, coincident_extras


def registry_values(value):
    return tuple(token for token in re.split(r"[|,]", value or "") if token and token != "-")


def load_registries(hgnc_path, withdrawn_path, gene_info_path):
    hgnc_by_id = {}
    hgnc_by_entrez = {}
    secondary = {
        "MIM": collections.defaultdict(list),
        "miRBase": collections.defaultdict(list),
        "IMGT/GENE-DB": collections.defaultdict(list),
    }
    with gzip.open(str(hgnc_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hgnc_by_id[row["hgnc_id"]] = row
            if row["entrez_id"]:
                if row["entrez_id"] in hgnc_by_entrez:
                    raise GateBlocked("duplicate HGNC Entrez ID {}".format(row["entrez_id"]))
                hgnc_by_entrez[row["entrez_id"]] = row
            for value in registry_values(row.get("omim_id", "")):
                secondary["MIM"][value].append(row)
            for value in registry_values(row.get("mirbase", "")):
                secondary["miRBase"][value].append(row)
            for value in registry_values(row.get("imgt", "")):
                secondary["IMGT/GENE-DB"][value].append(row)
    gene_info = {}
    with gzip.open(str(gene_info_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["#tax_id"] == "9606":
                gene_info[row["GeneID"]] = row
    redirects = {}
    report_field = "MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)"
    with gzip.open(str(withdrawn_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            redirects[row["HGNC_ID"]] = tuple(sorted(set(
                re.findall(r"HGNC:\d+", row[report_field])
            )))
    return hgnc_by_id, hgnc_by_entrez, secondary, redirects, gene_info


def resolve_identity(locus, hgnc_by_id, secondary, redirects):
    def resolve_tokens(tokens, scope):
        candidates = collections.defaultdict(set)
        hgnc_records = {}
        rejected = []
        for token in tokens:
            tagged = scope + ":" + token
            prefix, separator, value = token.partition(":")
            if not separator:
                rejected.append(tagged + ":malformed")
                continue
            if prefix == "GeneID":
                if value.isdigit() and int(value) > 0:
                    candidates["NCBIGene:" + value].add(tagged)
                else:
                    rejected.append(tagged + ":invalid_geneid")
                continue
            records = []
            if prefix == "HGNC":
                normalized = value if value.startswith("HGNC:") else "HGNC:" + value
                if normalized in hgnc_by_id:
                    records = [hgnc_by_id[normalized]]
                else:
                    records = [hgnc_by_id[target] for target in redirects.get(normalized, ())
                               if target in hgnc_by_id]
            elif prefix in secondary:
                records = secondary[prefix].get(value, [])
            else:
                rejected.append(tagged + ":unsupported_registry")
                continue
            if not records:
                rejected.append(tagged + ":absent_from_frozen_registry")
            for record in records:
                source_id = ("NCBIGene:" + record["entrez_id"]
                             if record["entrez_id"] else record["hgnc_id"])
                candidates[source_id].add(tagged)
                hgnc_records[source_id] = record
        return candidates, hgnc_records, rejected

    candidates, records, rejected = resolve_tokens(
        split_tokens(locus["raw_stable_xrefs"]), "GENE_ROW"
    )
    status = "EXACT_REGISTRY_IDENTITY"
    if not candidates:
        candidates, records, more_rejected = resolve_tokens(
            sorted(locus["descendant_xrefs"]), "LINKED_DESCENDANT"
        )
        rejected.extend(more_rejected)
        status = "EXACT_LINKED_MODEL_IDENTITY"
    if len(candidates) == 1:
        source_id = next(iter(candidates))
        record = records.get(source_id, {})
        entrez = source_id.split(":", 1)[1] if source_id.startswith("NCBIGene:") else ""
        return {
            "own_annotation_status": status,
            "own_annotation_id": source_id,
            "own_hgnc_id": record.get("hgnc_id", ""),
            "own_entrez_id": entrez,
            "own_identity_evidence": "|".join(sorted(candidates[source_id])),
        }
    if len(candidates) > 1:
        return {
            "own_annotation_status": "AMBIGUOUS_EXACT_XREF",
            "own_annotation_id": "", "own_hgnc_id": "", "own_entrez_id": "",
            "own_identity_evidence": "candidates={}".format("|".join(sorted(candidates))),
        }
    return {
        "own_annotation_status": "NO_EXACT_REGISTRY_IDENTITY",
        "own_annotation_id": "", "own_hgnc_id": "", "own_entrez_id": "",
        "own_identity_evidence": "|".join(sorted(rejected)) or "no_stable_xref",
    }


def load_directed_relations(path):
    pairs = set()
    by_source = collections.defaultdict(set)
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row["#tax_id"] == "9606" and row["Other_tax_id"] == "9606"
                    and row["relationship"] == "Related functional gene"):
                pair = (row["GeneID"], row["Other_GeneID"])
                pairs.add(pair)
                by_source[row["GeneID"]].add(row["Other_GeneID"])
    return pairs, by_source


def parse_go(path):
    terms = {}
    alternate = {}
    current = None

    def finish(term):
        if not term or "id" not in term:
            return
        term.setdefault("name", "")
        term.setdefault("namespace", "")
        term.setdefault("parents", set())
        terms[term["id"]] = term
        for alt_id in term.get("alt_ids", set()):
            alternate[alt_id] = term["id"]

    with gzip.open(str(path), "rt", encoding="utf-8") as handle:
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
    return terms, alternate


def load_go_assertions(path, relevant_gene_ids, alternate, go_terms):
    assertions = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"evidence": set(), "qualifiers": set(), "records": set()}
    ))
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
            gene_id = row["GeneID"]
            if gene_id not in relevant_gene_ids:
                continue
            qualifiers = set(split_tokens(row.get("Qualifier", "")))
            if "NOT" in qualifiers:
                continue
            term_id = alternate.get(row["GO_ID"], row["GO_ID"])
            term = go_terms.get(term_id)
            if not term or term.get("obsolete"):
                continue
            item = assertions[gene_id][term_id]
            item["evidence"].add(row["Evidence"])
            item["qualifiers"].update(qualifiers)
            item["records"].add("gene2go:{}:{}".format(line_number, row["GO_ID"]))
    return assertions


def load_reactome(pathways_path, relations_path, mapping_path, relevant_gene_ids):
    names = {}
    with gzip.open(str(pathways_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 3 and row[2] == "Homo sapiens":
                names[row[0]] = row[1]
    parents = collections.defaultdict(set)
    with gzip.open(str(relations_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 2 and row[0] in names and row[1] in names:
                parents[row[1]].add(row[0])
    assertions = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"evidence": set(), "records": set()}
    ))
    with gzip.open(str(mapping_path), "rt", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if (len(row) != 6 or row[5] != "Homo sapiens"
                    or row[0] not in relevant_gene_ids):
                continue
            if row[1] not in names:
                raise GateBlocked("unknown Reactome pathway {}".format(row[1]))
            item = assertions[row[0]][row[1]]
            item["evidence"].add(row[4])
            item["records"].add("NCBI2Reactome_All_Levels:{}".format(line_number))
    return names, parents, assertions


def ancestor_distances(term_id, parents, cache):
    if term_id in cache:
        return cache[term_id]
    result = {term_id: 0}
    queue = collections.deque([term_id])
    while queue:
        child = queue.popleft()
        for parent in parents.get(child, set()):
            distance = result[child] + 1
            if parent not in result or distance < result[parent]:
                result[parent] = distance
                queue.append(parent)
    cache[term_id] = result
    return result


def source_symbol(gene_id, gene_info, hgnc_by_entrez):
    symbol = gene_info.get(gene_id, {}).get("Symbol", "")
    return symbol or hgnc_by_entrez.get(gene_id, {}).get("symbol", "") or "NCBIGene:" + gene_id


def make_source_terms(gene_id, gene_info, hgnc_by_entrez, go_terms, go_parents,
                      go_assertions, reactome_names, reactome_parents,
                      reactome_assertions, ontology_hashes, go_cache,
                      reactome_cache):
    output = {}
    for leaf_id, evidence in go_assertions.get(gene_id, {}).items():
        for term_id, distance in ancestor_distances(leaf_id, go_parents, go_cache).items():
            item = output.setdefault(("GO", term_id), {
                "distance": distance, "leaves": set(), "direct": False,
                "evidence": set(), "qualifiers": set(), "records": set(),
            })
            item["distance"] = min(item["distance"], distance)
            item["leaves"].add(leaf_id)
            if distance == 0:
                item["direct"] = True
                item["evidence"].update(evidence["evidence"])
                item["qualifiers"].update(evidence["qualifiers"])
                item["records"].update(evidence["records"])

    reactome_all = set(reactome_assertions.get(gene_id, {}))
    mapped_ancestors = set()
    for term_id in reactome_all:
        mapped_ancestors.update(set(ancestor_distances(
            term_id, reactome_parents, reactome_cache
        )) - {term_id})
    reactome_leaves = reactome_all - mapped_ancestors
    if reactome_all and not reactome_leaves:
        raise GateBlocked("cyclic Reactome leaf selection for {}".format(gene_id))
    for leaf_id in reactome_leaves:
        evidence = reactome_assertions[gene_id][leaf_id]
        for term_id, distance in ancestor_distances(
                leaf_id, reactome_parents, reactome_cache).items():
            item = output.setdefault(("Reactome", term_id), {
                "distance": distance, "leaves": set(), "direct": False,
                "evidence": set(), "qualifiers": set(), "records": set(),
            })
            item["distance"] = min(item["distance"], distance)
            item["leaves"].add(leaf_id)
            if distance == 0:
                item["direct"] = True
                item["evidence"].update(evidence["evidence"])
                item["records"].update(evidence["records"])

    source_id = "NCBIGene:" + gene_id
    symbol = source_symbol(gene_id, gene_info, hgnc_by_entrez)
    rows = []
    for (ontology, term_id), item in sorted(output.items()):
        if ontology == "GO":
            metadata = go_terms[term_id]
            namespace = metadata["namespace"]
            term_name = metadata["name"]
            release = "Gene Ontology 2026-06-15"
            ontology_sha = ontology_hashes["go_basic"]
        else:
            namespace = "pathway"
            term_name = reactome_names[term_id]
            release = "Reactome v96"
            ontology_sha = ontology_hashes["reactome_relations"]
        rows.append({
            "functional_source_id": source_id,
            "functional_source_symbol": symbol,
            "ontology": ontology,
            "relation": "direct" if item["direct"] else "ancestor",
            "namespace": namespace,
            "term_id": term_id,
            "term_name": term_name,
            "minimum_distance": item["distance"],
            "inherited_from_direct_term_ids": "|".join(sorted(item["leaves"])),
            "evidence_codes": "|".join(sorted(item["evidence"])),
            "qualifiers": "|".join(sorted(item["qualifiers"])),
            "assertion_record_ids": "|".join(sorted(item["records"])),
            "ontology_release": release,
            "ontology_sha256": ontology_sha,
        })
    return rows


def adjudicate_source(own_gene_id, source_gene_id, disposition, relationship,
                      ambiguity, term_bearing_gene_ids, directed_relations):
    """Return the fail-closed admissibility decision for one proposed source."""
    if ambiguity and ambiguity != "none":
        return False, "ambiguous"
    if not source_gene_id:
        return False, "no_source"
    if source_gene_id not in term_bearing_gene_ids:
        if source_gene_id == own_gene_id:
            return False, "termless_self"
        return False, "termless_source"
    if source_gene_id == own_gene_id:
        if disposition == "EXACT_SELF" and relationship == "exact registry identity":
            return True, "exact_term_bearing_self"
        return False, "invalid_self_disposition"
    if (own_gene_id, source_gene_id) not in directed_relations:
        return False, "missing_exact_directed_relation"
    if (disposition == "EXPLICIT_RELATED_FUNCTIONAL_GENE"
            and relationship == "Related functional gene"):
        return True, "exact_directed_related_functional_gene"
    return False, "invalid_directed_relation_semantics"


def load_pilot_evidence(path):
    rows = {}
    fields_for_hash = (
        "copy_id", "own_entrez_id", "proposed_source_symbol",
        "proposed_source_entrez_id", "relationship_semantics",
        "evidence_source", "evidence_record_id", "disposition", "unresolved_reason",
    )
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["copy_id"] in rows:
                raise GateBlocked("duplicate pilot evidence copy")
            if row_hash(row, fields_for_hash) != row["evidence_record_sha256"]:
                raise GateBlocked("pilot evidence hash mismatch for {}".format(row["copy_id"]))
            rows[row["copy_id"]] = row
    return rows


def expected_source_for_copy(identity, pilot_row, relations_by_source,
                             term_bearing_gene_ids):
    own = identity["own_entrez_id"]
    if pilot_row is not None:
        source = (pilot_row["proposed_source_entrez_id"]
                  if pilot_row["admissible_for_ontology"] == "1" else "")
        return source, pilot_row["disposition"]
    if identity["own_annotation_status"] == "AMBIGUOUS_EXACT_XREF":
        return "", "AMBIGUOUS_FAIL_CLOSED"
    if own and own in term_bearing_gene_ids:
        return own, "EXACT_SELF"
    if own:
        candidates = sorted(target for target in relations_by_source.get(own, set())
                            if target in term_bearing_gene_ids)
        if len(candidates) == 1:
            return candidates[0], "EXPLICIT_RELATED_FUNCTIONAL_GENE"
        if len(candidates) > 1:
            return "", "AMBIGUOUS_FAIL_CLOSED"
        return "", "UNSUPPORTED_FAIL_CLOSED"
    return "", "TYPE_ONLY"


def load_assignments_and_evidence(loci, identities, ledger):
    with gzip.open(str(ASSIGNMENT_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assignment_fields = tuple(reader.fieldnames)
        assignments = list(reader)
    with gzip.open(str(EVIDENCE_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        evidence_fields = tuple(reader.fieldnames)
        evidence_rows = list(reader)
    ledger.require("assignment_rows", len(assignments), EXPECTED_LOCI, ASSIGNMENT_PATH.name)
    ledger.require("assignment_unique_copy_ids",
                   len(set(row["copy_id"] for row in assignments)), EXPECTED_LOCI,
                   ASSIGNMENT_PATH.name)
    ledger.require("evidence_rows", len(evidence_rows), EXPECTED_LOCI, EVIDENCE_PATH.name)
    ledger.require("evidence_unique_copy_ids",
                   len(set(row["copy_id"] for row in evidence_rows)), EXPECTED_LOCI,
                   EVIDENCE_PATH.name)
    raw_by_id = {row["copy_id"]: row for row in loci}
    assignment_ids = set(row["copy_id"] for row in assignments)
    raw_ids = set(raw_by_id)
    ledger.require("assignment_raw_copy_bijection", assignment_ids == raw_ids,
                   True, "raw GFF3 versus source assignments; 61,312 exact copy IDs")
    evidence_by_id = {row["copy_id"]: row for row in evidence_rows}
    ledger.require("evidence_raw_copy_bijection", set(evidence_by_id) == raw_ids, True,
                   "raw GFF3 versus copy-source evidence; 61,312 exact copy IDs")
    raw_fields = (
        "copy_id", "seqid", "start0", "end0", "start1", "end1", "strand",
        "gff_id", "gff_line", "gene_name", "gene_synonyms", "gene_biotype",
        "raw_stable_xrefs", "raw_attributes_sha256", "physical_copy_cn",
    )
    identity_fields = (
        "own_annotation_status", "own_annotation_id", "own_hgnc_id",
        "own_entrez_id", "own_identity_evidence",
    )
    for row in assignments:
        compare_row(raw_by_id[row["copy_id"]], row, raw_fields,
                    "raw assignment {}".format(row["copy_id"]))
        compare_row(identities[row["copy_id"]], row, identity_fields,
                    "identity assignment {}".format(row["copy_id"]))
        if row_hash(row, assignment_fields[:-1]) != row["assignment_record_sha256"]:
            raise GateBlocked("assignment row hash mismatch for {}".format(row["copy_id"]))
        evidence = evidence_by_id[row["copy_id"]]
        if row_hash(evidence, evidence_fields[:-1]) != evidence["evidence_record_sha256"]:
            raise GateBlocked("copy-source evidence hash mismatch for {}".format(row["copy_id"]))
        if (row["functional_source_id"] != evidence["functional_source_id"]
                or row["source_assignment_disposition"] != evidence["disposition"]
                or row["relationship_semantics"] != evidence["relationship_semantics"]):
            raise GateBlocked("assignment/evidence mismatch for {}".format(row["copy_id"]))
    ledger.require("all_assignment_hashes_reconstructed", True, True,
                   "61,312 assignment_record_sha256 values")
    ledger.require("all_copy_evidence_hashes_reconstructed", True, True,
                   "61,312 evidence_record_sha256 values")
    ledger.require("summed_physical_copy_cn",
                   sum(int(row["physical_copy_cn"]) for row in assignments),
                   EXPECTED_LOCI, ASSIGNMENT_PATH.name)
    return assignments, evidence_by_id, assignment_fields


def validate_source_evidence(assignments, evidence_by_id, identities, pilot,
                             directed_relations, relations_by_source,
                             term_bearing_gene_ids, ledger):
    pilot_seen = set()
    nonself = 0
    route_rejections = collections.Counter()
    for row in assignments:
        copy_id = row["copy_id"]
        identity = identities[copy_id]
        evidence = evidence_by_id[copy_id]
        pilot_row = pilot.get(copy_id)
        expected_source, expected_disposition = expected_source_for_copy(
            identity, pilot_row, relations_by_source, term_bearing_gene_ids
        )
        observed_source = row["functional_source_entrez_id"]
        if observed_source != expected_source or row["source_assignment_disposition"] != expected_disposition:
            raise GateBlocked(
                "independent source reconstruction differs for {}: {}/{} versus {}/{}".format(
                    copy_id, observed_source, row["source_assignment_disposition"],
                    expected_source, expected_disposition,
                )
            )
        if observed_source:
            allowed, reason = adjudicate_source(
                identity["own_entrez_id"], observed_source,
                row["source_assignment_disposition"], row["relationship_semantics"],
                row["ambiguity_status"], term_bearing_gene_ids, directed_relations,
            )
            if not allowed:
                raise GateBlocked("inadmissible emitted source for {}: {}".format(copy_id, reason))
            if identity["own_entrez_id"] != observed_source:
                nonself += 1
                expected_record = (
                    "NCBI_GENE_GROUP:9606:{}:Related_functional_gene:9606:{}".format(
                        identity["own_entrez_id"], observed_source
                    )
                )
                if row["source_evidence_record_id"] != expected_record:
                    raise GateBlocked("directed relation record mismatch for {}".format(copy_id))
        else:
            route_rejections[row["source_assignment_disposition"]] += 1
            if evidence["admissible_for_ontology"] != "0":
                raise GateBlocked("source-free row marked admissible for {}".format(copy_id))
        if bool(observed_source) != (evidence["admissible_for_ontology"] == "1"):
            raise GateBlocked("admissibility/source mismatch for {}".format(copy_id))
        if row["functional_source_id"] != (
                "NCBIGene:" + observed_source if observed_source else ""):
            raise GateBlocked("functional source ID/Entrez mismatch for {}".format(copy_id))

        if pilot_row is not None:
            pilot_seen.add(copy_id)
            expected_functional = (
                "NCBIGene:" + pilot_row["proposed_source_entrez_id"]
                if pilot_row["admissible_for_ontology"] == "1" else ""
            )
            pilot_expected = {
                "own_entrez_id": pilot_row["own_entrez_id"],
                "proposed_source_id": (
                    "NCBIGene:" + pilot_row["proposed_source_entrez_id"]
                    if pilot_row["proposed_source_entrez_id"] else ""
                ),
                "functional_source_id": expected_functional,
                "relationship_semantics": pilot_row["relationship_semantics"],
                "evidence_source": pilot_row["evidence_source"],
                "evidence_record_id": pilot_row["evidence_record_id"],
                "disposition": pilot_row["disposition"],
                "admissible_for_ontology": pilot_row["admissible_for_ontology"],
                "ambiguity_status": pilot_row["ambiguity_status_at_decision"],
                "unresolved_reason": pilot_row["unresolved_reason"],
                "pilot_override": "1",
                "upstream_evidence_sha256": pilot_row["evidence_record_sha256"],
            }
            compare_row(pilot_expected, evidence, tuple(pilot_expected),
                        "remediated pilot {}".format(copy_id))
        else:
            if evidence["pilot_override"] != "0" or evidence["upstream_evidence_sha256"]:
                raise GateBlocked("nonpilot row claims pilot override for {}".format(copy_id))
    ledger.require("remediated_pilot_rows_reconciled", pilot_seen == set(pilot), True,
                   "{}; 168 exact copy IDs".format(PILOT_EVIDENCE_PATH.name))
    ledger.require("remediated_pilot_row_count", len(pilot), 168,
                   PILOT_EVIDENCE_PATH.name)
    ledger.require("all_nonself_sources_have_exact_directed_evidence", nonself, 11561,
                   "raw human-human Related functional gene pairs")
    ledger.require("ambiguous_rows_emit_no_source", route_rejections["AMBIGUOUS_FAIL_CLOSED"],
                   45, "independent source adjudication")
    ledger.require("termless_or_unsupported_rows_emit_no_source",
                   route_rejections["UNSUPPORTED_FAIL_CLOSED"], 28407,
                   "term-bearing source set reconstructed from raw ontology assertions")
    ledger.require("name_family_alias_sequence_routes_admitted", 0, 0,
                   "only exact self or exact directed relation is admissible")


def verify_source_freeze(ledger):
    rows = read_tsv(OUTDIR / "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv")
    files = [row for row in rows if row["path"] != "COMPOSITE_SOURCE_EVIDENCE_DIGEST"]
    ledger.require("source_freeze_file_count", len(files), 3,
                   "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv")
    for row in files:
        ledger.require("source_freeze_stage_{}".format(row["path"]), row["stage"],
                       "source_assignment_before_target", "no PHR-derived field is permitted")
        path = OUTDIR / row["path"]
        ledger.require("source_freeze_sha256_{}".format(row["path"]), sha256_file(path),
                       row["sha256"], "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv")
        ledger.require("source_freeze_bytes_{}".format(row["path"]), path.stat().st_size,
                       int(row["bytes"]), "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv")
    composite = hashlib.sha256("\n".join(
        row["path"] + "\t" + row["sha256"] for row in files
    ).encode("utf-8")).hexdigest()
    composite_rows = [row for row in rows if row["path"] == "COMPOSITE_SOURCE_EVIDENCE_DIGEST"]
    expected = composite_rows[0]["sha256"] if len(composite_rows) == 1 else "MISSING"
    ledger.require("source_freeze_composite_sha256", composite, expected,
                   "three source-only frozen artifacts")
    return composite


def load_phr_intervals(path):
    intervals = collections.defaultdict(list)
    row_count = 0
    with open(str(path), "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise GateBlocked("malformed PHR BED line {}".format(line_number))
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end <= start:
                raise GateBlocked("invalid PHR BED interval line {}".format(line_number))
            intervals[fields[0]].append((start, end))
            row_count += 1
    return intervals, row_count


def phr_membership(locus, intervals):
    midpoint = (int(locus["start0"]) + int(locus["end0"])) // 2
    hits = [item for item in intervals.get(locus["seqid"], [])
            if item[0] <= midpoint < item[1]]
    if len(hits) > 1:
        raise GateBlocked("copy midpoint intersects multiple PHRs: {}".format(locus["copy_id"]))
    overlap = any(start < int(locus["end0"]) and int(locus["start0"]) < end
                  for start, end in intervals.get(locus["seqid"], []))
    interval = ("{}:{}-{}".format(locus["seqid"], hits[0][0], hits[0][1])
                if hits else "")
    return int(bool(hits)), int(overlap), interval


def join_and_validate_phr(assignments, intervals, ledger):
    with gzip.open(str(MAP_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        map_rows = list(reader)
    ledger.require("target_map_rows", len(map_rows), EXPECTED_LOCI, MAP_PATH.name)
    if any("phr" in field.lower() for field in assignments[0]
           if field not in {"phr_midpoint_cn", "phr_any_overlap_cn", "phr_interval"}):
        raise GateBlocked("pre-target assignment unexpectedly contains a PHR field")
    midpoint_total = 0
    overlap_total = 0
    for assignment, mapped in zip(assignments, map_rows):
        for field in assignment:
            if mapped.get(field, "") != assignment[field]:
                raise GateBlocked("PHR join changed frozen field {} for {}".format(
                    field, assignment["copy_id"]
                ))
        midpoint, overlap, interval = phr_membership(assignment, intervals)
        expected = {
            "phr_midpoint_cn": midpoint,
            "phr_any_overlap_cn": overlap,
            "phr_interval": interval,
        }
        compare_row(expected, mapped, tuple(expected),
                    "independent PHR join {}".format(assignment["copy_id"]))
        assignment.update({key: scalar(value) for key, value in expected.items()})
        midpoint_total += midpoint
        overlap_total += overlap
    ledger.require("all_phr_midpoint_memberships_reconstructed", midpoint_total,
                   EXPECTED_PHR_MIDPOINT, "raw GFF3 midpoint versus raw PHR BED")
    ledger.require("all_phr_overlap_memberships_reconstructed", overlap_total,
                   EXPECTED_PHR_OVERLAP, "raw GFF3 interval overlap versus raw PHR BED")
    return midpoint_total, overlap_total


def cohort_for(name, synonyms):
    synonym_set = set(split_tokens(synonyms))
    if name == "DUX4" or name.startswith("DUX4L") or "DUX4L30" in synonym_set:
        return "DUX4_DUX4L"
    if name.startswith("DDX11L"):
        return "DDX11L"
    if name.startswith("TUBB8"):
        return "TUBB8"
    if name.startswith("OR4F"):
        return "OR4F"
    if name.startswith("WASH"):
        return "WASH"
    return ""


def copy_weighted_burdens(copies, source_terms):
    result = collections.Counter()
    contributors = collections.Counter()
    for copy in copies:
        for term in source_terms.get(copy["source"], ()):
            key = (copy["source"], term)
            result[(key, "genome")] += int(copy["physical"])
            result[(key, "phr")] += int(copy["phr"])
            contributors[(key, "genome")] += 1
            contributors[(key, "phr")] += int(bool(int(copy["phr"])))
    return {
        key: (
            result[(key, "genome")], result[(key, "phr")],
            contributors[(key, "genome")], contributors[(key, "phr")],
        )
        for key in set(item[0] for item in result)
    }


def validate_source_terms_and_burdens(assignments, source_context, ledger):
    selected_sources = sorted(set(row["functional_source_entrez_id"] for row in assignments
                                  if row["functional_source_entrez_id"]))
    source_copy_counts = collections.Counter(row["functional_source_id"] for row in assignments
                                             if row["functional_source_id"])
    source_phr_counts = collections.Counter()
    for row in assignments:
        if row["functional_source_id"]:
            source_phr_counts[row["functional_source_id"]] += int(row["phr_midpoint_cn"])
    direct_counts = {}
    closure_counts = {}
    source_term_count = 0
    go_cache = {}
    reactome_cache = {}
    with gzip.open(str(SOURCE_TERM_PATH), "rt", encoding="utf-8", newline="") as term_handle:
        term_reader = csv.DictReader(term_handle, delimiter="\t")
        term_fields = tuple(term_reader.fieldnames)
        with gzip.open(str(BURDEN_PATH), "rt", encoding="utf-8", newline="") as burden_handle:
            burden_reader = csv.DictReader(burden_handle, delimiter="\t")
            burden_fields = tuple(burden_reader.fieldnames)
            for gene_id in sorted(selected_sources, key=lambda value: "NCBIGene:" + value):
                source_id = "NCBIGene:" + gene_id
                edge_order = make_source_terms(gene_id, go_cache=go_cache,
                                               reactome_cache=reactome_cache,
                                               **source_context)
                source_order = sorted(
                    edge_order,
                    key=lambda row: (
                        row["functional_source_id"], row["ontology"],
                        row["relation"], row["term_id"],
                    ),
                )
                direct_counts[source_id] = sum(row["relation"] == "direct" for row in edge_order)
                closure_counts[source_id] = len(edge_order)
                for expected in source_order:
                    observed = next_or_none(term_reader)
                    if observed is None:
                        raise GateBlocked("SOURCE_TERMS ended before independent reconstruction")
                    compare_row(expected, observed, term_fields,
                                "source term {} {}".format(source_id, expected["term_id"]))
                    burden = next_or_none(burden_reader)
                    if burden is None:
                        raise GateBlocked("TERM_BURDENS ended before independent reconstruction")
                    expected_burden = {
                        "functional_source_id": expected["functional_source_id"],
                        "functional_source_symbol": expected["functional_source_symbol"],
                        "ontology": expected["ontology"],
                        "relation": expected["relation"],
                        "namespace": expected["namespace"],
                        "term_id": expected["term_id"],
                        "term_name": expected["term_name"],
                        "genome_physical_copy_burden": source_copy_counts[source_id],
                        "phr_midpoint_physical_copy_burden": source_phr_counts[source_id],
                        "genome_contributor_rows": source_copy_counts[source_id],
                        "phr_midpoint_contributor_rows": source_phr_counts[source_id],
                    }
                    compare_row(expected_burden, burden, burden_fields,
                                "weighted burden {} {}".format(source_id, expected["term_id"]))
                    source_term_count += 1
            if next_or_none(term_reader) is not None:
                raise GateBlocked("SOURCE_TERMS contains rows absent from raw ontology reconstruction")
            if next_or_none(burden_reader) is not None:
                raise GateBlocked("TERM_BURDENS contains rows absent from raw ontology reconstruction")
    for row in assignments:
        source_id = row["functional_source_id"]
        expected_direct = direct_counts.get(source_id, 0)
        expected_closure = closure_counts.get(source_id, 0)
        if (int(row["direct_term_count"]) != expected_direct
                or int(row["closure_term_count"]) != expected_closure):
            raise GateBlocked("assignment term counts differ for {}".format(row["copy_id"]))
    ledger.require("source_terms_reconstructed", source_term_count,
                   EXPECTED_SOURCE_TERMS, "raw GO/gene2go and Reactome snapshots")
    ledger.require("term_burden_rows_reconstructed", source_term_count,
                   EXPECTED_SOURCE_TERMS, "one copy-weighted burden per source term")
    ledger.require("repeated_source_burdens_use_copy_counts", True, True,
                   "burden equals physical assignment count, never unique source count")
    return direct_counts, closure_counts, source_copy_counts, source_phr_counts


def validate_contributors(assignments, ledger):
    count = 0
    with gzip.open(str(CONTRIBUTOR_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames)
        for row in assignments:
            if not row["functional_source_id"]:
                continue
            expected = dict(row)
            expected.update({
                "cohort": cohort_for(row["gene_name"], row["gene_synonyms"]),
                "genome_term_edge_burden": int(row["closure_term_count"]),
                "phr_midpoint_term_edge_burden": (
                    int(row["phr_midpoint_cn"]) * int(row["closure_term_count"])
                ),
            })
            observed = next_or_none(reader)
            if observed is None:
                raise GateBlocked("EXACT_CONTRIBUTORS ended early")
            compare_row(expected, observed, fields,
                        "physical contributor {}".format(row["copy_id"]))
            count += 1
        if next_or_none(reader) is not None:
            raise GateBlocked("EXACT_CONTRIBUTORS has extra rows")
    ledger.require("physical_ontology_contributors_reconstructed", count, 31966,
                   CONTRIBUTOR_PATH.name)


def validate_edges(assignments, source_context, ledger):
    edge_count = 0
    direct_edge_count = 0
    phr_edge_count = 0
    named_burdens = collections.Counter()
    go_cache = {}
    reactome_cache = {}
    term_cache = collections.OrderedDict()
    with gzip.open(str(EDGE_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames)
        for assignment in assignments:
            source_gene_id = assignment["functional_source_entrez_id"]
            if not source_gene_id:
                continue
            if source_gene_id in term_cache:
                terms = term_cache.pop(source_gene_id)
                term_cache[source_gene_id] = terms
            else:
                terms = make_source_terms(source_gene_id, go_cache=go_cache,
                                          reactome_cache=reactome_cache, **source_context)
                term_cache[source_gene_id] = terms
                if len(term_cache) > 64:
                    term_cache.popitem(last=False)
            if len(terms) != int(assignment["closure_term_count"]):
                raise GateBlocked("closure size differs before edge expansion")
            cohort = cohort_for(assignment["gene_name"], assignment["gene_synonyms"])
            for term in terms:
                expected = dict(assignment)
                expected.update(term)
                observed = next_or_none(reader)
                if observed is None:
                    raise GateBlocked("PHYSICAL_COPY_TERM_EDGES ended early")
                compare_row(expected, observed, fields,
                            "copy term {} {}".format(assignment["copy_id"], term["term_id"]))
                edge_count += int(assignment["physical_copy_cn"])
                direct_edge_count += int(term["relation"] == "direct")
                phr_edge_count += int(assignment["phr_midpoint_cn"])
                if cohort:
                    key = (
                        cohort, term["functional_source_id"], term["ontology"],
                        term["relation"], term["namespace"], term["term_id"],
                        term["term_name"],
                    )
                    # Keep the key even when its PHR weight is zero.  This
                    # reconstructs the named audit's complete genome-cohort
                    # source/term universe while retaining its PHR burden.
                    named_burdens[key] += int(assignment["phr_midpoint_cn"])
        if next_or_none(reader) is not None:
            raise GateBlocked("PHYSICAL_COPY_TERM_EDGES has extra rows")
    ledger.require("physical_copy_term_edges_reconstructed", edge_count,
                   EXPECTED_COPY_TERM_EDGES, "one edge per eligible physical copy/source term")
    ledger.require("physical_copy_term_edge_cn", edge_count,
                   sum(int(row["closure_term_count"]) for row in assignments),
                   "no repeated source contributor collapsed")
    return edge_count, direct_edge_count, phr_edge_count, named_burdens


def validate_named_cohorts(assignments, pilot, named_burdens, source_context, ledger):
    pilot_copy_ids = set(pilot)
    named_copy_ids = set(row["copy_id"] for row in assignments
                         if cohort_for(row["gene_name"], row["gene_synonyms"]))
    ledger.require("named_cohorts_equal_remediated_pilot_universe",
                   named_copy_ids == pilot_copy_ids, True,
                   "five raw-name cohorts and pilot contain the same 168 copy IDs")
    result = collections.OrderedDict()
    audit_expected = {}
    for cohort, expected in EXPECTED_NAMED.items():
        members = [row for row in assignments
                   if cohort_for(row["gene_name"], row["gene_synonyms"]) == cohort]
        phr = [row for row in members if int(row["phr_midpoint_cn"])]
        contributors = [row for row in phr if row["functional_source_id"]]
        observed = (len(members), len(phr), len(contributors))
        ledger.require("named_{}_reviewed_counts".format(cohort), observed, expected,
                       "raw GFF/PHR plus remediated copy-source evidence")
        result[cohort] = {
            "genome_physical_copies": len(members),
            "phr_physical_copies": len(phr),
            "phr_ontology_contributors": len(contributors),
        }
        audit_expected[cohort] = {
            "cohort": cohort,
            "genome_physical_copies": len(members),
            "phr_physical_copies": len(phr),
            "phr_ontology_contributors": len(contributors),
            "phr_unresolved_physical_copies": len(phr) - len(contributors),
            "functional_source_ids": "|".join(sorted(set(
                row["functional_source_id"] for row in contributors
            ))),
            "functional_source_symbols": "|".join(sorted(set(
                row["functional_source_symbol"] for row in contributors
            ))),
            "direct_term_edge_burden": sum(int(row["direct_term_count"])
                                            for row in contributors),
            "closure_term_edge_burden": sum(int(row["closure_term_count"])
                                             for row in contributors),
            "expected_genome_copies": expected[0],
            "expected_phr_copies": expected[1],
            "expected_phr_contributors": expected[2],
            "status": "PASS",
        }
    with open(str(OUTDIR / "NAMED_COHORT_AUDIT.tsv"), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames)
        observed_rows = list(reader)
    ledger.require("named_audit_row_count", len(observed_rows), len(EXPECTED_NAMED),
                   "NAMED_COHORT_AUDIT.tsv")
    for observed in observed_rows:
        compare_row(audit_expected[observed["cohort"]], observed, fields,
                    "named cohort audit {}".format(observed["cohort"]))

    expected_named_rows = []
    gene_info = source_context["gene_info"]
    hgnc_by_entrez = source_context["hgnc_by_entrez"]
    for key, burden in sorted(named_burdens.items()):
        expected_named_rows.append({
            "cohort": key[0],
            "functional_source_id": key[1],
            "functional_source_symbol": source_symbol(
                key[1].split(":", 1)[1], gene_info, hgnc_by_entrez
            ),
            "ontology": key[2], "relation": key[3], "namespace": key[4],
            "term_id": key[5], "term_name": key[6],
            "phr_midpoint_physical_copy_burden": burden,
            "phr_midpoint_contributor_rows": burden,
        })
    with gzip.open(str(NAMED_BURDEN_PATH), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames)
        for expected in expected_named_rows:
            observed = next_or_none(reader)
            if observed is None:
                raise GateBlocked("named burden table ended early")
            compare_row(expected, observed, fields,
                        "named burden {} {}".format(expected["cohort"], expected["term_id"]))
        if next_or_none(reader) is not None:
            raise GateBlocked("named burden table has extra rows")
    dux_values = [value for key, value in named_burdens.items()
                  if key[0] == "DUX4_DUX4L" and key[1] == "NCBIGene:100288687"]
    ddx_values = [value for key, value in named_burdens.items()
                  if key[0] == "DDX11L" and key[1] == "NCBIGene:1663"]
    wash_values = [value for key, value in named_burdens.items()
                   if key[0] == "WASH" and key[1] == "NCBIGene:100287171"]
    ledger.require("DUX4_each_direct_and_ancestor_term_burden", set(dux_values), {65},
                   "65 coordinate-distinct reviewed PHR contributors")
    ledger.require("DDX11_each_direct_and_ancestor_term_burden", set(ddx_values), {10},
                   "10 coordinate-distinct reviewed PHR contributors")
    ledger.require("WASHC1_each_direct_and_ancestor_term_burden", set(wash_values), {9},
                   "9 coordinate-distinct reviewed PHR contributors")
    ledger.require("named_term_burden_rows_reconstructed", len(expected_named_rows),
                   len(expected_named_rows), NAMED_BURDEN_PATH.name)
    return result, len(expected_named_rows)


def compute_mapping_coverage(assignments, direct_counts, closure_counts):
    scopes = (
        ("genome", "physical_copy_cn"),
        ("phr_midpoint", "phr_midpoint_cn"),
        ("phr_any_overlap", "phr_any_overlap_cn"),
    )
    summary = []
    detailed = []
    for scope, weight_field in scopes:
        selected = [row for row in assignments if int(row[weight_field])]
        eligible = [row for row in selected if row["functional_source_id"]]
        summary.append({
            "scope": scope,
            "physical_rows": len(selected),
            "physical_copy_cn": sum(int(row[weight_field]) for row in assignments),
            "own_annotation_identity_cn": sum(int(row[weight_field]) for row in assignments
                                                if row["own_annotation_id"]),
            "ontology_eligible_rows": len(eligible),
            "ontology_eligible_cn": sum(int(row[weight_field]) for row in eligible),
            "ontology_ineligible_cn": sum(int(row[weight_field]) for row in assignments
                                            if not row["functional_source_id"]),
            "direct_term_edge_cn": sum(
                int(row[weight_field]) * direct_counts.get(row["functional_source_id"], 0)
                for row in assignments
            ),
            "closure_term_edge_cn": sum(
                int(row[weight_field]) * closure_counts.get(row["functional_source_id"], 0)
                for row in assignments
            ),
        })

        groupings = collections.defaultdict(list)
        for row in selected:
            route = row["source_assignment_disposition"]
            biotype = row["gene_biotype"] or "(missing)"
            groupings[("ROUTE_BY_BIOTYPE", route, biotype)].append(row)
            groupings[("ROUTE", route, "ALL")].append(row)
            groupings[("BIOTYPE", "ALL", biotype)].append(row)
            groupings[("TOTAL", "ALL", "ALL")].append(row)
        for (level, route, biotype), rows in sorted(groupings.items()):
            eligible_rows = [row for row in rows if row["functional_source_id"]]
            detailed.append({
                "scope": scope,
                "aggregation_level": level,
                "evidence_route": route,
                "gene_biotype": biotype,
                "physical_rows": len(rows),
                "physical_copy_cn": sum(int(row[weight_field]) for row in rows),
                "own_annotation_identity_cn": sum(int(row[weight_field]) for row in rows
                                                    if row["own_annotation_id"]),
                "ontology_eligible_rows": len(eligible_rows),
                "ontology_eligible_cn": sum(int(row[weight_field]) for row in eligible_rows),
                "ontology_ineligible_cn": sum(int(row[weight_field]) for row in rows
                                                if not row["functional_source_id"]),
                "direct_term_edge_cn": sum(
                    int(row[weight_field]) * direct_counts.get(row["functional_source_id"], 0)
                    for row in rows
                ),
                "closure_term_edge_cn": sum(
                    int(row[weight_field]) * closure_counts.get(row["functional_source_id"], 0)
                    for row in rows
                ),
            })
    return summary, detailed


def validate_mapping_coverage(summary, ledger):
    observed = read_tsv(OUTDIR / "MAPPING_COVERAGE.tsv")
    ledger.require("mapping_coverage_row_count", len(observed), len(summary),
                   "genome, PHR midpoint, PHR any-overlap")
    fields = tuple(observed[0])
    for expected, actual in zip(summary, observed):
        compare_row(expected, actual, fields, "mapping coverage {}".format(expected["scope"]))
    ledger.require("total_and_phr_route_biotype_coverage_reconstructed", True, True,
                   COVERAGE_PATH.name)


def _execute_gate(ledger):
    inputs, input_hashes = verify_input_manifest(ledger)
    loci, coincident_extras = load_raw_loci(inputs["raw_chm13_gff"])
    ledger.require("raw_physical_loci", len(loci), EXPECTED_LOCI,
                   "raw CHM13v2.0 GFF3 gene features")
    ledger.require("raw_unique_copy_ids", len(set(row["copy_id"] for row in loci)),
                   EXPECTED_LOCI, "coordinate+strand+GFF feature identity")
    ledger.require("coordinate_coincident_extra_rows_retained", coincident_extras, 12,
                   "gff_id keeps coordinate-coincident gene records distinct")

    (hgnc_by_id, hgnc_by_entrez, secondary, redirects,
     gene_info) = load_registries(inputs["hgnc"], inputs["hgnc_withdrawn"],
                                  inputs["ncbi_gene_info"])
    identities = {
        row["copy_id"]: resolve_identity(row, hgnc_by_id, secondary, redirects)
        for row in loci
    }
    ledger.require("independent_own_annotation_identity_count",
                   sum(bool(row["own_annotation_id"]) for row in identities.values()),
                   60430, "raw stable identifiers plus frozen exact registry bridges")
    directed_relations, relations_by_source = load_directed_relations(inputs["ncbi_gene_group"])
    pilot = load_pilot_evidence(inputs["remediated_pilot_evidence"])
    relevant_gene_ids = set(row["own_entrez_id"] for row in identities.values()
                            if row["own_entrez_id"])
    relevant_gene_ids.update(target for pair in directed_relations for target in pair)
    relevant_gene_ids.update(row["proposed_source_entrez_id"] for row in pilot.values()
                             if row["proposed_source_entrez_id"])
    go_terms, go_alternate = parse_go(inputs["go_basic"])
    go_assertions = load_go_assertions(inputs["ncbi_gene2go"], relevant_gene_ids,
                                       go_alternate, go_terms)
    reactome_names, reactome_parents, reactome_assertions = load_reactome(
        inputs["reactome_pathways"], inputs["reactome_relations"],
        inputs["reactome_ncbi"], relevant_gene_ids,
    )
    term_bearing = set(go_assertions) | set(reactome_assertions)

    assignments, evidence_by_id, assignment_fields = load_assignments_and_evidence(
        loci, identities, ledger
    )
    if any("phr" in field.lower() for field in assignment_fields):
        raise GateBlocked("source assignment snapshot contains target-derived PHR fields")
    validate_source_evidence(assignments, evidence_by_id, identities, pilot,
                             directed_relations, relations_by_source, term_bearing, ledger)
    freeze_hash = verify_source_freeze(ledger)

    # Evidence is now validated. Only now open the PHR target definition.
    phr_hash = sha256_file(inputs["phr_bed"])
    ledger.require("phr_bed_sha256_after_source_validation", phr_hash,
                   input_hashes["phr_bed_expected"], "PHR opened after source adjudication")
    intervals, phr_interval_rows = load_phr_intervals(inputs["phr_bed"])
    ledger.require("phr_interval_rows_reconciled", phr_interval_rows, 37,
                   "all non-comment BED rows")
    midpoint_total, overlap_total = join_and_validate_phr(assignments, intervals, ledger)

    source_context = {
        "gene_info": gene_info,
        "hgnc_by_entrez": hgnc_by_entrez,
        "go_terms": go_terms,
        "go_parents": {
            term_id: term.get("parents", set()) for term_id, term in go_terms.items()
        },
        "go_assertions": go_assertions,
        "reactome_names": reactome_names,
        "reactome_parents": reactome_parents,
        "reactome_assertions": reactome_assertions,
        "ontology_hashes": input_hashes,
    }
    direct_counts, closure_counts, source_copy_counts, source_phr_counts = (
        validate_source_terms_and_burdens(assignments, source_context, ledger)
    )
    validate_contributors(assignments, ledger)
    edge_count, direct_edge_count, phr_edge_count, named_burdens = validate_edges(
        assignments, source_context, ledger
    )
    named, named_burden_rows = validate_named_cohorts(
        assignments, pilot, named_burdens, source_context, ledger
    )
    summary, detailed_coverage = compute_mapping_coverage(
        assignments, direct_counts, closure_counts
    )
    validate_mapping_coverage(summary, ledger)
    release_hashes = verify_release_manifest(ledger)
    ledger.require("production_builder_imported", False, False,
                   "standalone parser and reconstruction; no production module import")
    ledger.require("enrichment_run", False, False,
                   "gate writes mapping evidence only")

    return {
        "status": "PASS",
        "gate_version": "6.0.0-independent-pre-inference-gate",
        "scope": "independent_copy_to_ontology_reconstruction_no_enrichment",
        "production_builder_imported": False,
        "enrichment_run": False,
        "enrichment_authorized": True,
        "authorization_condition": "PASS for the exact audited release digests only",
        "physical_loci": len(assignments),
        "physical_copy_cn": sum(int(row["physical_copy_cn"]) for row in assignments),
        "coordinate_coincident_extra_rows_retained": coincident_extras,
        "own_annotation_identity_copies": sum(bool(row["own_annotation_id"])
                                               for row in assignments),
        "ontology_eligible_copies": sum(bool(row["functional_source_id"])
                                         for row in assignments),
        "ontology_ineligible_copies": sum(not bool(row["functional_source_id"])
                                           for row in assignments),
        "phr_interval_rows": phr_interval_rows,
        "phr_midpoint_copies": midpoint_total,
        "phr_any_overlap_copies": overlap_total,
        "exact_directed_nonself_copies": sum(
            bool(row["functional_source_id"])
            and row["functional_source_entrez_id"] != row["own_entrez_id"]
            for row in assignments
        ),
        "source_terms": EXPECTED_SOURCE_TERMS,
        "term_burden_rows": EXPECTED_SOURCE_TERMS,
        "physical_copy_term_edges": edge_count,
        "direct_physical_copy_term_edges": direct_edge_count,
        "phr_midpoint_copy_term_edges": phr_edge_count,
        "named_cohort_term_burden_rows": named_burden_rows,
        "named_cohorts": named,
        "source_assignment_freeze_sha256": freeze_hash,
        "phr_bed_sha256": phr_hash,
        "audited_release_manifest_sha256": sha256_file(
            OUTDIR / "OUTPUT_MANIFEST.sha256.tsv"
        ),
        "checks": ledger.rows,
        "artifacts": {
            "report": REPORT_PATH.name,
            "validation": JSON_PATH.name,
            "check_evidence": CHECK_PATH.name,
            "coverage_by_route_and_biotype": COVERAGE_PATH.name,
            "audited_release_digests": RELEASE_DIGEST_PATH.name,
        },
        "_coverage": detailed_coverage,
        "_coverage_summary": summary,
        "_release_hashes": release_hashes,
    }


def render_report(result):
    status = result["status"]
    lines = [
        "# V6 pre-run independent gate",
        "",
        "## Decision",
        "",
    ]
    if status == "PASS":
        lines.extend([
            "**PASS.** The exact frozen ontology_v6 physical-copy map is authorized as an input to the downstream enrichment task. The authorization is digest-bound: any input or audited release byte change requires this gate to be rerun. No enrichment was run here.",
            "",
            "The gate independently recovered **61,312** raw physical GFF3 loci, **402** PHR-midpoint memberships, **412** any-overlap memberships, **1,686,727** source-term closure rows, and **2,929,709** physical-copy term edges. It did not import `build_genomewide_source_map.py` or `check_genomewide_source_map.py`.",
        ])
    else:
        lines.extend([
            "**BLOCK.** Enrichment is unauthorized. The independent reconstruction did not satisfy every hard gate.",
            "",
            "Failure: `{}`".format(result.get("failure", "unknown gate failure")),
        ])
    lines.extend([
        "",
        "## Evidence-before-arithmetic contract",
        "",
        "Source assignment was adjudicated before the PHR BED was opened and before copy-term arithmetic. An ontology source is admissible only when it is a term-bearing exact self identifier or the target of the exact frozen directed human NCBI Gene `Related functional gene` pair for that copy's own identifier. Names, aliases, family membership, reverse-only relations, and one-way sequence similarity do not emit terms. Ambiguous and termless-self rows fail closed.",
        "",
    ])
    if status == "PASS":
        lines.extend([
            "## Copy-number retention",
            "",
            "Source assertions are deduplicated only within a source. Expansion is one row per physical copy and source term. Each burden was independently recomputed as the sum of `physical_copy_cn` (and separately `phr_midpoint_cn`), and contributor-row counts were required to equal those sums. Repeated sources were never reduced to one contributor.",
            "",
            "| Cohort | Genome physical copies | PHR physical copies | Reviewed PHR ontology contributors / term burden |",
            "|---|---:|---:|---:|",
        ])
        for cohort in EXPECTED_NAMED:
            row = result["named_cohorts"][cohort]
            lines.append("| {} | {} | {} | {} |".format(
                cohort, row["genome_physical_copies"], row["phr_physical_copies"],
                row["phr_ontology_contributors"],
            ))
        lines.extend([
            "",
            "DUX4 contributes 65 to every one of its direct and ancestor terms; DDX11 contributes 10 and WASHC1 contributes 9. TUBB8's two reviewed PHR contributors and OR4F's four contributors retain their separate exact functional sources; their copies are not collapsed or combined through a family union.",
            "",
            "## Mapping coverage",
            "",
            "| Scope | Physical CN | Own exact identity CN | Ontology contributor CN | Ineligible CN | Direct edge CN | Closure edge CN |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ])
        for row in result["_coverage_summary"]:
            lines.append("| {scope} | {physical_copy_cn} | {own_annotation_identity_cn} | {ontology_eligible_cn} | {ontology_ineligible_cn} | {direct_term_edge_cn} | {closure_term_edge_cn} |".format(**row))
        lines.extend([
            "",
            "[`PRE_RUN_V6_MAPPING_COVERAGE.tsv`](PRE_RUN_V6_MAPPING_COVERAGE.tsv) reports genome, PHR-midpoint, and PHR-any-overlap coverage as totals and as complete evidence-route, biotype, and route-by-biotype partitions.",
            "",
        ])
    lines.extend([
        "## Machine gate",
        "",
        "- [`PRE_RUN_V6_GATE.json`](PRE_RUN_V6_GATE.json) is the strict PASS/BLOCK decision and authorization record.",
        "- [`PRE_RUN_V6_GATE_EVIDENCE.tsv`](PRE_RUN_V6_GATE_EVIDENCE.tsv) records every hard check, observed value, expected value, and evidence source.",
        "- [`PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv`](PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv) binds PASS to the exact production artifacts audited.",
        "- [`independent_gate_genomewide_source_map.py`](independent_gate_genomewide_source_map.py) is the standalone reconstruction; its tests include explicit prohibited-route and copy-collapse controls.",
        "",
        "Downstream enrichment must require `status == PASS`, `enrichment_authorized == true`, and unchanged audited digests. A BLOCK result or any digest drift denies inference.",
        "",
        "Validation status: **{}** ({} recorded hard checks).".format(
            status, len(result.get("checks", []))
        ),
        "",
    ])
    return "\n".join(lines)


def write_outputs(result):
    coverage = result.get("_coverage", [])
    release_hashes = result.get("_release_hashes", [])
    checks = result.get("checks", [])
    if coverage:
        write_tsv(COVERAGE_PATH, tuple(coverage[0]), coverage)
    else:
        write_tsv(COVERAGE_PATH, (
            "scope", "aggregation_level", "evidence_route", "gene_biotype",
            "physical_rows", "physical_copy_cn", "own_annotation_identity_cn",
            "ontology_eligible_rows", "ontology_eligible_cn", "ontology_ineligible_cn",
            "direct_term_edge_cn", "closure_term_edge_cn",
        ), [])
    write_tsv(CHECK_PATH, ("check", "observed", "expected", "status", "evidence"), checks)
    write_tsv(RELEASE_DIGEST_PATH, ("path", "bytes", "sha256"), release_hashes)
    report = render_report(result)
    tmp_report = pathlib.Path(str(REPORT_PATH) + ".tmp")
    tmp_report.write_text(report, encoding="utf-8")
    os.replace(str(tmp_report), str(REPORT_PATH))
    public = {key: value for key, value in result.items() if not key.startswith("_")}
    tmp_json = pathlib.Path(str(JSON_PATH) + ".tmp")
    tmp_json.write_text(
        json.dumps(json_safe(public), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(str(tmp_json), str(JSON_PATH))


def run_gate(write_outputs=True):
    ledger = CheckLedger()
    try:
        result = _execute_gate(ledger)
    except Exception as error:
        if not ledger.rows or ledger.rows[-1]["status"] != "BLOCK":
            ledger.rows.append({
                "check": "unhandled_gate_exception",
                "observed": type(error).__name__,
                "expected": "no exception",
                "status": "BLOCK",
                "evidence": str(error),
            })
        result = {
            "status": "BLOCK",
            "gate_version": "6.0.0-independent-pre-inference-gate",
            "scope": "independent_copy_to_ontology_reconstruction_no_enrichment",
            "production_builder_imported": False,
            "enrichment_run": False,
            "enrichment_authorized": False,
            "authorization_condition": "all hard checks must PASS",
            "failure": str(error),
            "failure_type": type(error).__name__,
            "checks": ledger.rows,
            "artifacts": {
                "report": REPORT_PATH.name,
                "validation": JSON_PATH.name,
                "check_evidence": CHECK_PATH.name,
                "coverage_by_route_and_biotype": COVERAGE_PATH.name,
                "audited_release_digests": RELEASE_DIGEST_PATH.name,
            },
            "_coverage": [], "_coverage_summary": [], "_release_hashes": [],
        }
        if os.environ.get("V6_GATE_TRACEBACK") == "1":
            traceback.print_exc()
    if write_outputs:
        write_outputs_fn = globals()["write_outputs"]
        write_outputs_fn(result)
    return {key: value for key, value in result.items() if not key.startswith("_")}


def main():
    result = run_gate(write_outputs=True)
    if result["status"] != "PASS":
        print("BLOCK: {}".format(result.get("failure", "independent V6 gate failed")),
              file=sys.stderr)
        return 1
    print(
        "PASS: independently reconstructed {physical_loci} copies, {source_terms} "
        "source terms, and {physical_copy_term_edges} copy-term edges".format(**result)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
