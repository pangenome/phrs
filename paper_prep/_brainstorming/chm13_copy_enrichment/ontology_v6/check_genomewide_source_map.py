#!/usr/bin/env python3
"""Independently check the V6 genome-wide copy/source mapping release."""

import collections
import csv
import gzip
import hashlib
import json
import pathlib
import re
import urllib.parse


SCRIPT = pathlib.Path(__file__).resolve()
OUTDIR = SCRIPT.parent
REPO = SCRIPT.parents[4]
GFF = REPO / "data" / "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
PHR = REPO / "data" / "chm13.phrs.bed"
GENE_GROUP = (
    REPO / "paper_prep/_brainstorming/chm13_copy_enrichment/sources/"
    "ncbi_gene_group_2026-07-16.tsv.gz"
)
EXPECTED_COHORTS = {
    "DUX4_DUX4L": (107, 68, 65),
    "DDX11L": (12, 10, 10),
    "TUBB8": (16, 7, 2),
    "OR4F": (15, 11, 4),
    "WASH": (18, 9, 9),
}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_attributes(raw):
    attrs = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = urllib.parse.unquote(value)
    return attrs


def raw_gff_copies():
    result = {}
    with gzip.open(str(GFF), "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            copy_id = "CHM13v2.0|{}:{}-{}|{}|{}".format(
                fields[0], fields[3], fields[4], fields[6], attrs["ID"]
            )
            if copy_id in result:
                raise AssertionError("duplicate raw copy ID")
            result[copy_id] = (
                fields[0], int(fields[3]) - 1, int(fields[4]), fields[6],
                attrs["ID"], line_number,
            )
    return result


def load_phrs():
    result = collections.defaultdict(list)
    with open(str(PHR), "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                result[fields[0]].append((int(fields[1]), int(fields[2])))
    return result


def exact_relations():
    result = set()
    with gzip.open(str(GENE_GROUP), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row["#tax_id"] == "9606" and row["Other_tax_id"] == "9606"
                    and row["relationship"] == "Related functional gene"):
                result.add((row["GeneID"], row["Other_GeneID"]))
    return result


def verify_manifests():
    inputs = read_tsv(OUTDIR / "INPUT_MANIFEST.tsv")
    if len(inputs) != 14:
        raise AssertionError("input manifest does not pin 14 inputs")
    for row in inputs:
        path = REPO / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise AssertionError("input manifest byte mismatch: {}".format(path))
        observed = sha256_file(path)
        if observed != row["sha256"] or observed != row["expected_sha256"]:
            raise AssertionError("input manifest hash mismatch: {}".format(path))
    phr_rows = [row for row in inputs if row["role"] == "phr_bed"]
    if len(phr_rows) != 1 or phr_rows[0]["stage"] != "after_source_freeze_target_join":
        raise AssertionError("PHR input stage is not explicitly post-freeze")

    outputs = read_tsv(OUTDIR / "OUTPUT_MANIFEST.sha256.tsv")
    if len(outputs) < 17:
        raise AssertionError("output manifest is incomplete")
    for row in outputs:
        path = OUTDIR / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise AssertionError("output manifest byte mismatch: {}".format(path))
        if sha256_file(path) != row["sha256"]:
            raise AssertionError("output manifest hash mismatch: {}".format(path))


def verify_freeze():
    freeze = read_tsv(OUTDIR / "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv")
    if len(freeze) != 4:
        raise AssertionError("source freeze must contain three files and a composite")
    file_rows = [row for row in freeze if row["path"] != "COMPOSITE_SOURCE_EVIDENCE_DIGEST"]
    for row in file_rows:
        if row["stage"] != "source_assignment_before_target" or "phr" in row["path"].lower():
            raise AssertionError("target data leaked into the source freeze")
        path = OUTDIR / row["path"]
        if path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise AssertionError("source freeze file mismatch")
    observed = hashlib.sha256("\n".join(
        row["path"] + "\t" + row["sha256"] for row in file_rows
    ).encode("utf-8")).hexdigest()
    composite = [row for row in freeze if row["path"] == "COMPOSITE_SOURCE_EVIDENCE_DIGEST"]
    if len(composite) != 1 or composite[0]["sha256"] != observed:
        raise AssertionError("source freeze composite mismatch")
    with gzip.open(str(OUTDIR / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz"),
                   "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if any("phr" in field.lower() for field in header):
        raise AssertionError("target-blind assignment schema contains a PHR field")
    return observed


def verify_source_terms_and_burdens(source_copy_counts, source_phr_counts):
    source_path = OUTDIR / "SOURCE_TERMS.tsv.gz"
    burden_path = OUTDIR / "TERM_BURDENS.tsv.gz"
    source_counts = collections.Counter()
    direct_counts = collections.Counter()
    row_count = 0
    with gzip.open(str(source_path), "rt", encoding="utf-8", newline="") as source_handle:
        with gzip.open(str(burden_path), "rt", encoding="utf-8", newline="") as burden_handle:
            source_reader = csv.DictReader(source_handle, delimiter="\t")
            burden_reader = csv.DictReader(burden_handle, delimiter="\t")
            while True:
                try:
                    source = next(source_reader)
                except StopIteration:
                    source = None
                try:
                    burden = next(burden_reader)
                except StopIteration:
                    burden = None
                if source is None or burden is None:
                    if source is not None or burden is not None:
                        raise AssertionError("source-term and burden tables differ in length")
                    break
                keys = (
                    "functional_source_id", "functional_source_symbol", "ontology",
                    "relation", "namespace", "term_id", "term_name",
                )
                if any(source[key] != burden[key] for key in keys):
                    raise AssertionError("source-term/burden key mismatch")
                source_id = source["functional_source_id"]
                expected_genome = source_copy_counts[source_id]
                expected_phr = source_phr_counts[source_id]
                if (int(burden["genome_physical_copy_burden"]) != expected_genome
                        or int(burden["genome_contributor_rows"]) != expected_genome
                        or int(burden["phr_midpoint_physical_copy_burden"]) != expected_phr
                        or int(burden["phr_midpoint_contributor_rows"]) != expected_phr):
                    raise AssertionError("N-copy burden invariant failed for {}".format(source_id))
                if source["relation"] not in {"direct", "ancestor"}:
                    raise AssertionError("invalid term relation")
                if source["relation"] == "ancestor" and int(source["minimum_distance"]) <= 0:
                    raise AssertionError("ancestor has nonpositive distance")
                if source["relation"] == "direct" and int(source["minimum_distance"]) != 0:
                    raise AssertionError("direct term has nonzero distance")
                source_counts[source_id] += 1
                direct_counts[source_id] += source["relation"] == "direct"
                row_count += 1
    return source_counts, direct_counts, row_count


def verify_edges(assignments, source_counts):
    observed_per_copy = collections.Counter()
    edge_count = 0
    current_copy = None
    current_keys = set()
    with gzip.open(str(OUTDIR / "PHYSICAL_COPY_TERM_EDGES.tsv.gz"),
                   "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            copy_id = row["copy_id"]
            assignment = assignments.get(copy_id)
            if assignment is None or not assignment["functional_source_id"]:
                raise AssertionError("term edge names absent/ineligible copy")
            if copy_id != current_copy:
                current_copy = copy_id
                current_keys = set()
            key = (row["ontology"], row["relation"], row["namespace"], row["term_id"])
            if key in current_keys:
                raise AssertionError("duplicate copy-term edge")
            current_keys.add(key)
            if (row["functional_source_id"] != assignment["functional_source_id"]
                    or row["own_annotation_id"] != assignment["own_annotation_id"]
                    or row["physical_copy_cn"] != "1"
                    or row["phr_midpoint_cn"] != assignment["phr_midpoint_cn"]
                    or row["phr_any_overlap_cn"] != assignment["phr_any_overlap_cn"]):
                raise AssertionError("copy edge lost identity or CN")
            observed_per_copy[copy_id] += 1
            edge_count += 1
    for copy_id, assignment in assignments.items():
        expected = int(assignment["closure_term_count"])
        if observed_per_copy[copy_id] != expected:
            raise AssertionError("copy closure edge count mismatch: {}".format(copy_id))
        if assignment["functional_source_id"] and expected != source_counts[assignment["functional_source_id"]]:
            raise AssertionError("copy did not inherit complete source closure")
    return edge_count


def verify_named_audits(assignments):
    audit = {row["cohort"]: row for row in read_tsv(OUTDIR / "NAMED_COHORT_AUDIT.tsv")}
    if set(audit) != set(EXPECTED_COHORTS):
        raise AssertionError("named audit cohort set mismatch")
    for cohort, expected in EXPECTED_COHORTS.items():
        row = audit[cohort]
        observed = (
            int(row["genome_physical_copies"]),
            int(row["phr_physical_copies"]),
            int(row["phr_ontology_contributors"]),
        )
        if observed != expected or row["status"] != "PASS":
            raise AssertionError("named cohort mismatch: {}".format(cohort))
    dux = [row for row in assignments.values()
           if row["phr_midpoint_cn"] == "1"
           and row["functional_source_symbol"] == "DUX4"
           and (row["gene_name"] == "DUX4" or row["gene_name"].startswith("DUX4L")
                or "DUX4L30" in set(re.split(r"[,|]", row["gene_synonyms"])))]
    if len(dux) != 65:
        raise AssertionError("DUX4-source contributor recount is not 65")
    burden_rows = read_tsv(OUTDIR / "NAMED_COHORT_TERM_BURDENS.tsv.gz")
    dux_burdens = [row for row in burden_rows
                   if row["cohort"] == "DUX4_DUX4L"
                   and row["functional_source_symbol"] == "DUX4"]
    if not dux_burdens or any(
            int(row["phr_midpoint_physical_copy_burden"]) != 65
            or int(row["phr_midpoint_contributor_rows"]) != 65
            for row in dux_burdens):
        raise AssertionError("DUX4 direct/ancestor burden is not uniformly 65")


def check_release():
    verify_manifests()
    freeze_hash = verify_freeze()
    raw = raw_gff_copies()
    if len(raw) != 61312:
        raise AssertionError("raw GFF recount is not 61,312")

    assignments_list = read_tsv(OUTDIR / "GENOMEWIDE_SOURCE_MAP.tsv.gz")
    assignments = {row["copy_id"]: row for row in assignments_list}
    source_assignments = read_tsv(OUTDIR / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz")
    evidence = read_tsv(OUTDIR / "COPY_SOURCE_EVIDENCE.tsv.gz")
    if (len(assignments_list) != 61312 or len(assignments) != 61312
            or len(source_assignments) != 61312 or len(evidence) != 61312):
        raise AssertionError("map/assignment/evidence is not a 61,312-row bijection")
    if set(assignments) != set(raw) or {row["copy_id"] for row in evidence} != set(raw):
        raise AssertionError("output copy IDs do not equal raw GFF copy IDs")
    for target, source in zip(assignments_list, source_assignments):
        if target["copy_id"] != source["copy_id"]:
            raise AssertionError("target/source map ordering differs")
        for key in source:
            if target[key] != source[key]:
                raise AssertionError("PHR join changed a frozen assignment field")
        raw_row = raw[target["copy_id"]]
        if (target["seqid"], int(target["start0"]), int(target["end0"]),
                target["strand"], target["gff_id"], int(target["gff_line"])) != raw_row:
            raise AssertionError("raw coordinate identity mismatch")
        if target["physical_copy_cn"] != "1":
            raise AssertionError("physical copy CN is not one")

    intervals = load_phrs()
    midpoint_total = 0
    overlap_total = 0
    for row in assignments_list:
        midpoint = (int(row["start0"]) + int(row["end0"])) // 2
        midpoint_hit = any(start <= midpoint < end for start, end in intervals.get(row["seqid"], []))
        overlap_hit = any(start < int(row["end0"]) and int(row["start0"]) < end
                          for start, end in intervals.get(row["seqid"], []))
        if int(row["phr_midpoint_cn"]) != int(midpoint_hit):
            raise AssertionError("PHR midpoint flag mismatch")
        if int(row["phr_any_overlap_cn"]) != int(overlap_hit):
            raise AssertionError("PHR overlap flag mismatch")
        midpoint_total += midpoint_hit
        overlap_total += overlap_hit
    if (midpoint_total, overlap_total) != (402, 412):
        raise AssertionError("independent PHR recount mismatch")

    relations = exact_relations()
    evidence_by_copy = {row["copy_id"]: row for row in evidence}
    for copy_id, assignment in assignments.items():
        record = evidence_by_copy[copy_id]
        if assignment["functional_source_id"] != record["functional_source_id"]:
            raise AssertionError("assignment/evidence source mismatch")
        if not assignment["functional_source_id"]:
            if record["admissible_for_ontology"] != "0":
                raise AssertionError("ineligible row marked admissible")
            continue
        if record["admissible_for_ontology"] != "1":
            raise AssertionError("eligible row lacks admissible evidence")
        own_id = assignment["own_entrez_id"]
        source_id = assignment["functional_source_entrez_id"]
        if own_id == source_id:
            if assignment["source_assignment_disposition"] != "EXACT_SELF":
                raise AssertionError("self source lacks exact-self disposition")
        elif ((own_id, source_id) not in relations
              or assignment["relationship_semantics"] != "Related functional gene"):
            raise AssertionError("non-self source lacks exact directed evidence")
    if any(row["functional_source_id"] for row in assignments.values()
           if row["source_assignment_disposition"] in {
               "AMBIGUOUS_FAIL_CLOSED", "UNSUPPORTED_FAIL_CLOSED", "TYPE_ONLY", "UNRESOLVED"
           }):
        raise AssertionError("fail-closed row emitted a functional source")

    source_copy_counts = collections.Counter(
        row["functional_source_id"] for row in assignments.values()
        if row["functional_source_id"]
    )
    source_phr_counts = collections.Counter()
    for row in assignments.values():
        if row["functional_source_id"]:
            source_phr_counts[row["functional_source_id"]] += int(row["phr_midpoint_cn"])
    source_counts, direct_counts, source_term_rows = verify_source_terms_and_burdens(
        source_copy_counts, source_phr_counts
    )
    if any(source_counts[source] <= 0 or direct_counts[source] <= 0 for source in source_copy_counts):
        raise AssertionError("eligible source lacks direct/closure terms")
    edge_count = verify_edges(assignments, source_counts)
    verify_named_audits(assignments)

    validation = json.loads((OUTDIR / "GENOMEWIDE_SOURCE_MAP_VALIDATION.json").read_text())
    if (validation["status"] != "PASS" or validation["enrichment_run"] is not False
            or validation["source_assignment_freeze_sha256"] != freeze_hash
            or validation["physical_loci"] != 61312
            or validation["physical_copy_term_edges"] != edge_count
            or validation["source_terms"] != source_term_rows):
        raise AssertionError("validation JSON summary mismatch")
    if not all(row["passed"] for row in validation["checks"]):
        raise AssertionError("builder validation contains a failed check")
    return {
        "status": "PASS",
        "physical_loci": len(assignments),
        "phr_midpoint_copies": midpoint_total,
        "source_term_rows": source_term_rows,
        "physical_copy_term_edges": edge_count,
        "source_assignment_freeze_sha256": freeze_hash,
        "enrichment_run": False,
    }


def main():
    result = check_release()
    print("PASS: independently checked {physical_loci} copies, "
          "{source_term_rows} source terms, and {physical_copy_term_edges} edges".format(**result))


if __name__ == "__main__":
    main()
