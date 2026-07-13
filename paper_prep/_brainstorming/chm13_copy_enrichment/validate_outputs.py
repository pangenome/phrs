#!/usr/bin/env python3
"""Independently validate frozen CHM13 copy-to-term output invariants."""

import collections
import csv
import gzip
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "outputs"
REPO = HERE.parents[2]


def read_gzip_tsv(name):
    with gzip.open(str(OUTPUT / name), "rt", encoding="utf-8", newline="") as inp:
        return list(csv.DictReader(inp, delimiter="\t"))


def read_tsv(name):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as inp:
        return list(csv.DictReader(inp, delimiter="\t"))


def load_builder():
    spec = importlib.util.spec_from_file_location("term_builder", HERE / "build_term_maps.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    builder = load_builder()
    source_manifest = builder.read_source_manifest_and_verify()
    require(len(source_manifest) == 7, "expected seven frozen external source objects")

    copies = read_gzip_tsv("copy_universe.tsv.gz")
    copy_ids = [row["copy_id"] for row in copies]
    copy_by_id = {row["copy_id"]: row for row in copies}
    require(len(copy_ids) == len(set(copy_ids)), "duplicate copy_id in universe")

    gff_loci = builder.read_gff(builder.DEFAULT_GFF)
    gff_by_id = {row.copy_id: row for row in gff_loci}
    require(set(copy_ids) == set(gff_by_id), "copy universe is not exactly the audited GFF")
    for copy_id, row in copy_by_id.items():
        locus = gff_by_id[copy_id]
        require(row["seqid"] == locus.seqid, "seqid mismatch for " + copy_id)
        require(int(row["start_1based"]) == locus.start, "start mismatch for " + copy_id)
        require(int(row["end_1based_inclusive"]) == locus.end, "end mismatch for " + copy_id)
        require(row["gene_biotype"] == locus.gene_biotype, "biotype mismatch for " + copy_id)

    terms = read_gzip_tsv("copy_to_term.tsv.gz")
    term_keys = [(x["copy_id"], x["source"], x["term_id"]) for x in terms]
    require(len(term_keys) == len(set(term_keys)), "duplicate copy/source/term inflation")
    require(set(x["copy_id"] for x in terms).issubset(copy_by_id), "term copy outside universe")
    require(
        all("symbol" not in x["mapped_via"].lower() for x in terms),
        "symbol-based annotation route detected",
    )

    metadata = read_gzip_tsv("term_metadata.tsv.gz")
    metadata_keys = [(x["source"], x["term_id"]) for x in metadata]
    require(len(metadata_keys) == len(set(metadata_keys)), "duplicate term metadata key")
    metadata_set = set(metadata_keys)
    require(
        all((x["source"], x["term_id"]) in metadata_set for x in terms),
        "copy term without metadata",
    )
    go_namespaces = {
        x["namespace"] for x in metadata if x["source"] == "GO" and x["is_obsolete"] == "false"
    }
    require(
        {"biological_process", "molecular_function", "cellular_component"}.issubset(go_namespaces),
        "GO namespaces incomplete",
    )

    hierarchy = read_gzip_tsv("term_hierarchy.tsv.gz")
    hierarchy_keys = [
        (x["source"], x["child_term_id"], x["parent_term_id"], x["relation"])
        for x in hierarchy
    ]
    require(len(hierarchy_keys) == len(set(hierarchy_keys)), "duplicate hierarchy edge")
    for row in hierarchy:
        require((row["source"], row["child_term_id"]) in metadata_set, "missing child metadata")
        require((row["source"], row["parent_term_id"]) in metadata_set, "missing parent metadata")

    diagnostics = read_gzip_tsv("copy_mapping_diagnostics.tsv.gz")
    require(
        {x["copy_id"] for x in diagnostics} == set(copy_ids) and len(diagnostics) == len(copy_ids),
        "diagnostics are not one-to-one with copy universe",
    )
    rejected = {
        x["copy_id"]
        for x in diagnostics
        if x["stable_mapping_status"].startswith("rejected_")
    }
    require(
        not any(x["copy_id"] in rejected and x["source"] != "biotype" for x in terms),
        "rejected stable mapping received a functional term",
    )

    # Every copy has exactly one explicit biotype term and it matches the GFF.
    biotype_terms = collections.defaultdict(list)
    source_mapped = collections.defaultdict(set)
    for row in terms:
        source_mapped[row["source"]].add(row["copy_id"])
        if row["source"] == "biotype":
            biotype_terms[row["copy_id"]].append(row["term_id"])
    for copy_id, copy in copy_by_id.items():
        require(
            biotype_terms[copy_id] == ["BIOTYPE:" + copy["gene_biotype"]],
            "missing or inflated biotype term for " + copy_id,
        )

    coverage = read_tsv("coverage_by_biotype.tsv")
    for row in coverage:
        biotype = row["gene_biotype"]
        selected = copies if biotype == "ALL" else [x for x in copies if x["gene_biotype"] == biotype]
        if row["source"] == "stable_id":
            mapped = {
                x["copy_id"]
                for x in diagnostics
                if x["stable_mapping_status"].startswith("mapped_")
            }
        else:
            mapped = source_mapped[row["source"]]
        expected = sum(x["copy_id"] in mapped for x in selected)
        require(int(row["n_loci"]) == len(selected), "coverage denominator mismatch")
        require(int(row["n_mapped_loci"]) == expected, "coverage numerator mismatch")

    phr_ids = {x["copy_id"] for x in copies if x["in_phr_audited_target"] == "1"}
    expected_phr = {x.copy_id for x in builder.read_gff(builder.DEFAULT_PHR_GFF)}
    require(phr_ids == expected_phr, "PHR target flag differs from audited PHR GFF")

    audits = read_tsv("hand_audit_examples.tsv")
    require(bool(audits), "hand audit examples absent")
    require(len(audits) == len({x["copy_id"] for x in audits}), "duplicate hand-audit copy")
    require(set(x["copy_id"] for x in audits).issubset(copy_by_id), "audit copy outside universe")
    require(
        any("pseudogene" in x["gene_biotype"] for x in audits),
        "hand audit lacks pseudogenes",
    )
    require(any(x["gene_biotype"] == "rRNA" for x in audits), "hand audit lacks rDNA control")

    summary = json.loads((OUTPUT / "validation_summary.json").read_text(encoding="utf-8"))
    require(summary["physical_loci"] == len(copies), "summary physical-locus count mismatch")
    require(summary["term_rows"] == len(terms), "summary term-row count mismatch")
    require(summary["term_metadata_rows"] == len(metadata), "summary metadata count mismatch")
    require(summary["hierarchy_edges"] == len(hierarchy), "summary hierarchy count mismatch")
    require(summary["validation"]["symbol_annotation_joins"] == 0, "symbol joins reported")
    require(
        summary["validation"]["pseudogene_parent_inheritance_operations"] == 0,
        "pseudogene inheritance reported",
    )
    print(
        "validated {:,} physical loci, {:,} unique copy-term rows, {:,} metadata terms, "
        "and {:,} hierarchy edges".format(len(copies), len(terms), len(metadata), len(hierarchy))
    )


if __name__ == "__main__":
    main()
