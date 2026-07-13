#!/usr/bin/env python3
"""Build copy-aware CHM13 functional term maps from frozen sources.

Only exact, stable identifier joins are permitted.  Gene symbols are retained
for display and hand-audit selection, but are never indexed for annotation.
"""

import argparse
import collections
import csv
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_GFF = REPO / "data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
DEFAULT_PHR_GFF = REPO / "data/phrs.genes.gff3"
SOURCE_DIR = HERE / "sources"
OUTPUT_DIR = HERE / "outputs"

TERM_COLUMNS = [
    "copy_id",
    "source",
    "term_id",
    "annotation_relation",
    "evidence_codes",
    "support_count",
    "mapped_via",
    "stable_gene_ids",
]

HAND_AUDIT = {
    "OR4F subtelomeric receptor family": (
        "OR4F4",
        "OR4F5",
        "OR4F13P",
        "OR4F14P",
        "OR4F15",
    ),
    "TUBB8/TUBB8B family with pseudogenes": (
        "TUBB8",
        "TUBB8B",
        "TUBB8P4",
        "TUBB8P8",
        "TUBB8P11",
    ),
    "DUX4-like multicopy family": ("DUX4", "DUX4L6", "DUX4L10", "DUX4L24"),
    "45S rDNA arrays without GFF stable xrefs": ("RNA45S1", "RNA45S3", "RNA45S4"),
}


@dataclass(frozen=True)
class StableXrefs:
    gene_ids: Tuple[str, ...] = ()
    hgnc_ids: Tuple[str, ...] = ()
    omim_ids: Tuple[str, ...] = ()
    mirbase_ids: Tuple[str, ...] = ()
    imgt_ids: Tuple[str, ...] = ()


@dataclass(frozen=True)
class Locus:
    copy_id: str
    seqid: str
    start: int
    end: int
    strand: str
    gff_id: str
    liftoff_copy_num_id: str
    gene_name: str
    gene_biotype: str
    gene_ids: Tuple[str, ...]
    hgnc_ids: Tuple[str, ...]
    omim_ids: Tuple[str, ...]
    mirbase_ids: Tuple[str, ...]
    imgt_ids: Tuple[str, ...]

    @property
    def xrefs(self):
        return StableXrefs(
            self.gene_ids,
            self.hgnc_ids,
            self.omim_ids,
            self.mirbase_ids,
            self.imgt_ids,
        )


@dataclass(frozen=True)
class MappingResult:
    status: str
    hgnc_id: str = ""
    gene_ids: Tuple[str, ...] = ()
    routes: Tuple[str, ...] = ()
    detail: str = ""


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as inp:
        for chunk in iter(lambda: inp.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(text):
    attributes = collections.defaultdict(list)
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key].append(unquote(value))
    return attributes


def normalize_hgnc(value):
    value = value.strip()
    while value.startswith("HGNC:"):
        value = value[5:]
    return "HGNC:" + value if value else ""


def stable_xrefs(attributes):
    values = []
    for key in ("db_xref", "Dbxref"):
        for field in attributes.get(key, []):
            values.extend(x for x in field.split(",") if x)
    found = collections.defaultdict(set)
    for value in values:
        if ":" not in value:
            continue
        namespace, identifier = value.split(":", 1)
        if namespace == "GeneID":
            found["gene"].add(identifier)
        elif namespace == "HGNC":
            found["hgnc"].add(normalize_hgnc(identifier))
        elif namespace == "MIM":
            found["omim"].add(identifier)
        elif namespace == "miRBase":
            found["mirbase"].add(identifier)
        elif namespace == "IMGT/GENE-DB":
            found["imgt"].add(identifier)
    return StableXrefs(
        tuple(sorted(found["gene"])),
        tuple(sorted(found["hgnc"])),
        tuple(sorted(found["omim"])),
        tuple(sorted(found["mirbase"])),
        tuple(sorted(found["imgt"])),
    )


def parse_gff_lines(lines):
    seen_copy_ids = set()
    for line_number, line in enumerate(lines, 1):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "gene":
            continue
        attrs = parse_attributes(fields[8])
        gff_id = attrs.get("ID", [""])[0]
        if not gff_id:
            raise ValueError("gene row {} has no ID".format(line_number))
        seqid, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
        # Coordinate + feature ID is the physical key.  copy_num_ID is not a
        # safe primary key because the 45S array features omit it.
        copy_id = "CHM13v2.0|{}:{}-{}|{}|{}".format(seqid, start, end, strand, gff_id)
        if copy_id in seen_copy_ids:
            raise ValueError("duplicate physical copy key: {}".format(copy_id))
        seen_copy_ids.add(copy_id)
        xrefs = stable_xrefs(attrs)
        yield Locus(
            copy_id=copy_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            gff_id=gff_id,
            liftoff_copy_num_id=attrs.get("copy_num_ID", [""])[0],
            gene_name=attrs.get("gene_name", attrs.get("gene", [""]))[0],
            gene_biotype=attrs.get("gene_biotype", ["unknown"])[0],
            gene_ids=xrefs.gene_ids,
            hgnc_ids=xrefs.hgnc_ids,
            omim_ids=xrefs.omim_ids,
            mirbase_ids=xrefs.mirbase_ids,
            imgt_ids=xrefs.imgt_ids,
        )


def read_gff(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt", encoding="utf-8") as inp:
        return list(parse_gff_lines(inp))


def split_multi(value):
    return tuple(x for x in (value or "").split("|") if x)


class HGNCResolver:
    """Resolve stable GFF cross-references to at most one HGNC record."""

    INDEX_FIELDS = {
        "hgnc_id": "hgnc",
        "entrez_id": "gene",
        "omim_id": "omim",
        "mirbase": "mirbase",
        "imgt": "imgt",
    }

    def __init__(self, records, withdrawn_records=None):
        self.records = {row.get("hgnc_id", ""): row for row in records if row.get("hgnc_id")}
        self.indices = {name: collections.defaultdict(set) for name in self.INDEX_FIELDS.values()}
        for hgnc_id, row in self.records.items():
            for field, index_name in self.INDEX_FIELDS.items():
                for value in split_multi(row.get(field, "")):
                    if field == "hgnc_id":
                        value = normalize_hgnc(value)
                    self.indices[index_name][value].add(hgnc_id)
        # HGNC's withdrawn report is an exact identifier-lifecycle bridge.
        # Splits remain multi-valued here and are therefore rejected by
        # resolve(); entries without a replacement remain unmapped.
        for row in withdrawn_records or []:
            withdrawn_id = row.get("HGNC_ID", "")
            report = row.get("MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)", "")
            replacements = set(re.findall(r"HGNC:\d+", report))
            replacements.intersection_update(self.records)
            if withdrawn_id and replacements:
                self.indices["hgnc"][withdrawn_id].update(replacements)

    def resolve(self, xrefs):
        provided = []
        for index_name, values in (
            ("hgnc", xrefs.hgnc_ids),
            ("gene", xrefs.gene_ids),
            ("omim", xrefs.omim_ids),
            ("mirbase", xrefs.mirbase_ids),
            ("imgt", xrefs.imgt_ids),
        ):
            for value in values:
                candidates = self.indices[index_name].get(value, set())
                if len(candidates) > 1:
                    return MappingResult(
                        "rejected_ambiguous_stable_xref",
                        gene_ids=tuple(sorted(xrefs.gene_ids)),
                        detail="{}:{} maps to {}".format(
                            index_name, value, "|".join(sorted(candidates))
                        ),
                    )
                if candidates:
                    provided.append((index_name, value, next(iter(candidates))))
        hgnc_candidates = {x[2] for x in provided}
        if len(hgnc_candidates) > 1:
            return MappingResult(
                "rejected_conflicting_stable_xrefs",
                gene_ids=tuple(sorted(xrefs.gene_ids)),
                detail=";".join("{}:{}={}".format(*x) for x in provided),
            )
        hgnc_id = next(iter(hgnc_candidates)) if hgnc_candidates else ""
        gene_ids = set(xrefs.gene_ids)
        routes = set()
        if hgnc_id:
            record_gene_ids = set(split_multi(self.records[hgnc_id].get("entrez_id", "")))
            if gene_ids and record_gene_ids and not gene_ids.intersection(record_gene_ids):
                return MappingResult(
                    "rejected_conflicting_stable_xrefs",
                    gene_ids=tuple(sorted(gene_ids)),
                    detail="GFF GeneID {} conflicts with {} Entrez {}".format(
                        "|".join(sorted(gene_ids)), hgnc_id, "|".join(sorted(record_gene_ids))
                    ),
                )
            gene_ids.update(record_gene_ids)
            for index_name, value, resolved_id in provided:
                if index_name == "hgnc" and value != resolved_id:
                    routes.add("exact_withdrawn_hgnc_redirect")
                else:
                    routes.add("exact_{}_to_hgnc".format(index_name))
            return MappingResult(
                "mapped_unique_hgnc",
                hgnc_id,
                tuple(sorted(gene_ids)),
                tuple(sorted(routes)),
            )
        if gene_ids:
            return MappingResult(
                "mapped_direct_geneid_only",
                gene_ids=tuple(sorted(gene_ids)),
                routes=("direct_geneid",),
                detail="GeneID absent from current HGNC complete set",
            )
        any_stable = any(
            (xrefs.hgnc_ids, xrefs.omim_ids, xrefs.mirbase_ids, xrefs.imgt_ids)
        )
        if any_stable:
            detail = []
            for namespace, values in (
                ("HGNC", xrefs.hgnc_ids),
                ("MIM", xrefs.omim_ids),
                ("miRBase", xrefs.mirbase_ids),
                ("IMGT", xrefs.imgt_ids),
            ):
                detail.extend("{}:{}".format(namespace, value) for value in values)
            return MappingResult(
                "unmapped_stable_xref_not_in_hgnc",
                detail="no approved or uniquely redirected HGNC record for " + "|".join(detail),
            )
        return MappingResult(
            "unmapped_no_stable_xref",
            detail="GFF locus supplies no GeneID, HGNC, MIM, miRBase, or IMGT identifier",
        )


def map_locus(locus, resolver):
    return resolver.resolve(locus.xrefs)


class TermAccumulator:
    def __init__(self):
        self.data = {}

    def add(
        self,
        copy_id,
        source,
        term_id,
        relation,
        evidence="",
        mapped_via="",
        stable_gene_id="",
    ):
        key = (copy_id, source, term_id)
        value = self.data.setdefault(
            key,
            {"relations": set(), "evidence": set(), "routes": set(), "gene_ids": set()},
        )
        if relation:
            value["relations"].add(relation)
        if evidence:
            value["evidence"].add(evidence)
        if mapped_via:
            value["routes"].add(mapped_via)
        if stable_gene_id:
            value["gene_ids"].add(stable_gene_id)

    def rows(self):
        result = []
        for (copy_id, source, term_id), value in sorted(self.data.items()):
            support = max(1, len(value["evidence"]))
            result.append(
                {
                    "copy_id": copy_id,
                    "source": source,
                    "term_id": term_id,
                    "annotation_relation": "|".join(sorted(value["relations"])),
                    "evidence_codes": "|".join(sorted(value["evidence"])),
                    "support_count": str(support),
                    "mapped_via": "|".join(sorted(value["routes"])),
                    "stable_gene_ids": "|".join(sorted(value["gene_ids"])),
                }
            )
        return result


def open_snapshot(name):
    return gzip.open(str(SOURCE_DIR / name), "rt", encoding="utf-8")


def load_hgnc():
    path = SOURCE_DIR / "hgnc_complete_set_2026-07-10.tsv.gz"
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as inp:
        return list(csv.DictReader(inp, delimiter="\t"))


def load_hgnc_withdrawn():
    path = SOURCE_DIR / "hgnc_withdrawn_2026-07-10.tsv.gz"
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as inp:
        return list(csv.DictReader(inp, delimiter="\t"))


def parse_obo(path):
    terms = {}
    alt_ids = {}
    edges = set()
    current = None

    def finish(term):
        if not term or not term.get("id"):
            return
        term_id = term["id"][0]
        terms[term_id] = {
            "source": "GO",
            "term_id": term_id,
            "term_name": term.get("name", [""])[0],
            "namespace": term.get("namespace", [""])[0],
            "is_obsolete": term.get("is_obsolete", ["false"])[0],
            "description": term.get("def", [""])[0],
        }
        for alt_id in term.get("alt_id", []):
            alt_ids[alt_id] = term_id
        for parent in term.get("is_a", []):
            edges.add(("GO", term_id, parent.split()[0], "is_a"))
        for relation in term.get("relationship", []):
            parts = relation.split()
            if len(parts) >= 2 and parts[1].startswith("GO:"):
                edges.add(("GO", term_id, parts[1], parts[0]))

    with gzip.open(str(path), "rt", encoding="utf-8") as inp:
        for raw in inp:
            line = raw.rstrip("\n")
            if line == "[Term]":
                finish(current)
                current = collections.defaultdict(list)
            elif line.startswith("["):
                finish(current)
                current = None
            elif current is not None and ": " in line:
                key, value = line.split(": ", 1)
                current[key].append(value)
        finish(current)
    return terms, alt_ids, edges


def load_gene2go(alt_ids, go_terms):
    annotations = collections.defaultdict(lambda: collections.defaultdict(set))
    rejected = collections.Counter()
    path = SOURCE_DIR / "ncbi_gene2go_human_2026-07-13.tsv.gz"
    with gzip.open(str(path), "rt", encoding="utf-8") as inp:
        reader = csv.DictReader((line.lstrip("#") for line in inp), delimiter="\t")
        for row in reader:
            qualifier = row["Qualifier"]
            if "NOT" in qualifier.split("|"):
                rejected["negative_qualifier"] += 1
                continue
            original = row["GO_ID"]
            term_id = alt_ids.get(original, original)
            if term_id not in go_terms:
                rejected["term_absent_from_go_basic"] += 1
                continue
            annotations[row["GeneID"]][term_id].add(row["Evidence"] or "unspecified")
    return annotations, rejected


def load_reactome():
    annotations = collections.defaultdict(lambda: collections.defaultdict(set))
    mapping_path = SOURCE_DIR / "reactome_v96_ncbi_human_all_levels.tsv.gz"
    with gzip.open(str(mapping_path), "rt", encoding="utf-8") as inp:
        for line in inp:
            gene_id, term_id, _url, _name, evidence, species = line.rstrip("\n").split("\t")
            if species != "Homo sapiens":
                raise ValueError("non-human Reactome row in frozen snapshot")
            annotations[gene_id][term_id].add(evidence or "unspecified")

    terms = {}
    pathway_path = SOURCE_DIR / "reactome_v96_human_pathways.tsv.gz"
    with gzip.open(str(pathway_path), "rt", encoding="utf-8") as inp:
        for line in inp:
            term_id, name, species = line.rstrip("\n").split("\t")
            terms[term_id] = {
                "source": "Reactome",
                "term_id": term_id,
                "term_name": name,
                "namespace": "pathway",
                "is_obsolete": "false",
                "description": "Homo sapiens Reactome pathway (release 96)",
            }
    edges = set()
    relation_path = SOURCE_DIR / "reactome_v96_human_pathway_relations.tsv.gz"
    with gzip.open(str(relation_path), "rt", encoding="utf-8") as inp:
        for line in inp:
            parent, child = line.rstrip("\n").split("\t")[:2]
            if parent not in terms or child not in terms:
                raise ValueError("Reactome hierarchy endpoint absent from pathway metadata")
            edges.add(("Reactome", child, parent, "is_a_subpathway_of"))
    return annotations, terms, edges


def hgnc_groups(record):
    ids = split_multi(record.get("gene_group_id", ""))
    names = split_multi(record.get("gene_group", ""))
    if len(ids) != len(names):
        raise ValueError(
            "HGNC group ID/name length mismatch for {}: {} vs {}".format(
                record.get("hgnc_id"), len(ids), len(names)
            )
        )
    return tuple(zip(ids, names))


def deterministic_gzip_dict_writer(path, fieldnames, rows):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with open_text_wrapper(gz) as text:
                writer = csv.DictWriter(
                    text, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


class open_text_wrapper:
    """Close a text wrapper without closing its caller-owned gzip stream."""

    def __init__(self, binary):
        import io

        self.wrapper = io.TextIOWrapper(binary, encoding="utf-8", newline="")

    def __enter__(self):
        return self.wrapper

    def __exit__(self, exc_type, exc, tb):
        self.wrapper.flush()
        self.wrapper.detach()


def write_tsv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_source_manifest_and_verify():
    path = SOURCE_DIR / "SOURCE_MANIFEST.tsv"
    records = []
    with path.open(encoding="utf-8") as inp:
        for row in csv.DictReader(inp, delimiter="\t"):
            snapshot = HERE / row["snapshot"]
            if not snapshot.exists() or sha256(snapshot) != row["snapshot_sha256"]:
                raise ValueError("frozen source missing or checksum mismatch: {}".format(snapshot))
            records.append(row)
    return records


def coverage_rows(loci, term_rows, mappings):
    by_source_copy = collections.defaultdict(set)
    for row in term_rows:
        by_source_copy[row["source"]].add(row["copy_id"])
    stable_mapped = {
        copy_id
        for copy_id, result in mappings.items()
        if result.status.startswith("mapped_")
    }
    all_biotypes = sorted({x.gene_biotype for x in loci})
    rows = []
    for source in ("stable_id", "biotype", "GO", "HGNC_group", "Reactome"):
        mapped_set = stable_mapped if source == "stable_id" else by_source_copy[source]
        for biotype in ["ALL"] + all_biotypes:
            selected = loci if biotype == "ALL" else [x for x in loci if x.gene_biotype == biotype]
            mapped = sum(x.copy_id in mapped_set for x in selected)
            rows.append(
                {
                    "source": source,
                    "gene_biotype": biotype,
                    "n_loci": str(len(selected)),
                    "n_mapped_loci": str(mapped),
                    "coverage_percent": "{:.3f}".format(100.0 * mapped / len(selected)),
                }
            )
    return rows


def audit_rows(loci, mappings, per_copy_counts):
    result = []
    for family, names in HAND_AUDIT.items():
        for name in names:
            selected = sorted(
                (x for x in loci if x.gene_name == name),
                key=lambda x: (x.seqid, x.start, x.end, x.gff_id),
            )[:4]
            for locus in selected:
                mapping = mappings[locus.copy_id]
                result.append(
                    {
                        "audit_family": family,
                        "selection_rule": "exact gene_name in frozen hand-audit list; display only",
                        "gene_name": name,
                        "copy_id": locus.copy_id,
                        "locus": "{}:{}-{}({})".format(
                            locus.seqid, locus.start, locus.end, locus.strand
                        ),
                        "gene_biotype": locus.gene_biotype,
                        "gff_stable_xrefs": stable_xref_display(locus),
                        "mapping_status": mapping.status,
                        "resolved_hgnc_id": mapping.hgnc_id,
                        "resolved_gene_ids": "|".join(mapping.gene_ids),
                        "go_terms": str(per_copy_counts[(locus.copy_id, "GO")]),
                        "hgnc_group_terms": str(
                            per_copy_counts[(locus.copy_id, "HGNC_group")]
                        ),
                        "reactome_terms": str(per_copy_counts[(locus.copy_id, "Reactome")]),
                    }
                )
    return result


def stable_xref_display(locus):
    fields = []
    for prefix, values in (
        ("GeneID", locus.gene_ids),
        ("HGNC", locus.hgnc_ids),
        ("MIM", locus.omim_ids),
        ("miRBase", locus.mirbase_ids),
        ("IMGT", locus.imgt_ids),
    ):
        for value in values:
            fields.append(value if prefix == "HGNC" else "{}:{}".format(prefix, value))
    return "|".join(fields)


def write_report(
    loci,
    phr_count,
    coverage,
    mappings,
    term_rows,
    metadata,
    edges,
    audit,
    source_records,
    go_rejected,
):
    coverage_lookup = {
        (x["source"], x["gene_biotype"]): x for x in coverage
    }
    statuses = collections.Counter(x.status for x in mappings.values())
    routes = collections.Counter(
        "|".join(x.routes) if x.routes else "none" for x in mappings.values()
    )
    pseudogene_biotypes = sorted(
        {x.gene_biotype for x in loci if "pseudogene" in x.gene_biotype}
    )
    ncrna_biotypes = sorted(
        {
            x.gene_biotype
            for x in loci
            if x.gene_biotype != "protein_coding"
            and "pseudogene" not in x.gene_biotype
        }
    )
    lines = [
        "# Frozen CHM13 functional annotation term maps",
        "",
        "## Outcome",
        "",
        "This directory freezes a stable-identifier-only mapping from **{:,} physical CHM13v2.0 gene loci** to GO, HGNC gene groups, Reactome pathways, and the {} explicit GFF3 biotypes. The audited PHR target contains {:,} loci and is an exact subset of that universe. Every output annotation row uses the coordinate-anchored `copy_id` from `copy_universe.tsv.gz`; no gene symbol, alias, edit distance, case folding, or pseudogene-parent relation is used to assign a term.".format(
            len(loci), len({x.gene_biotype for x in loci}), phr_count
        ),
        "",
        "## Identity model and joins",
        "",
        "The physical key is `CHM13v2.0|seqid:start-end|strand|GFF_ID`. It remains unique for all loci, including 45S rDNA features that omit Liftoff `copy_num_ID`; the original value is retained separately. GFF3 `db_xref` and `Dbxref` are both parsed. Exact `GeneID` and `HGNC` joins are preferred. Exact OMIM, miRBase, and IMGT identifiers may bridge to one HGNC record through the frozen HGNC fields `omim_id`, `mirbase`, and `imgt`. Withdrawn HGNC IDs are redirected only through an exact, single-target HGNC merged report; withdrawn splits are ambiguous and rejected. A bridge is accepted only when it is one-to-one and concordant with every other supplied stable identifier. Symbols are display-only.",
        "",
        "Pseudogenes are ordinary physical loci. A pseudogene gets only its own stable-ID annotations. The pipeline contains no parent-gene parser or inheritance operation. Loci without a usable stable identifier still receive their own biotype term and remain explicit in diagnostics.",
        "",
        "Pre-existing enrichment spreadsheets, reports, and symbol-based draft notes under `_brainstorming/` are not inputs to this build. Biological claims from those drafts were neither copied nor used to fill mapping gaps.",
        "",
        "Mapping status | Loci",
        "--- | ---:",
    ]
    lines.extend("{} | {:,}".format(k, v) for k, v in sorted(statuses.items()))
    lines += [
        "",
        "Accepted exact route | Loci",
        "--- | ---:",
    ]
    lines.extend("{} | {:,}".format(k, v) for k, v in sorted(routes.items()))
    lines += [
        "",
        "## Coverage",
        "",
        "Coverage counts physical loci, not unique symbols or genes. `HGNC_group` means at least one HGNC gene-group term; stable-ID coverage can therefore exceed functional-source coverage.",
        "",
        "Source | All loci mapped | Coverage",
        "--- | ---: | ---:",
    ]
    for source in ("stable_id", "biotype", "GO", "HGNC_group", "Reactome"):
        row = coverage_lookup[(source, "ALL")]
        lines.append(
            "{} | {:,}/{:,} | {}%".format(
                source,
                int(row["n_mapped_loci"]),
                int(row["n_loci"]),
                row["coverage_percent"],
            )
        )
    lines += [
        "",
        "Full source-by-biotype coverage is in `outputs/coverage_by_biotype.tsv`. Required edge categories are summarized below.",
        "",
        "Biotype | Loci | Stable ID | GO | HGNC group | Reactome",
        "--- | ---: | ---: | ---: | ---: | ---:",
    ]
    focus = pseudogene_biotypes + ncrna_biotypes + ["protein_coding"]
    for biotype in focus:
        n = int(coverage_lookup[("biotype", biotype)]["n_loci"])
        cells = []
        for source in ("stable_id", "GO", "HGNC_group", "Reactome"):
            row = coverage_lookup[(source, biotype)]
            cells.append("{}/{} ({}%)".format(row["n_mapped_loci"], n, row["coverage_percent"]))
        lines.append("{} | {:,} | {}".format(biotype, n, " | ".join(cells)))
    lines += [
        "",
        "## Term products",
        "",
        "`copy_to_term.tsv.gz` has {:,} unique `(copy_id, source, term_id)` rows. Duplicate source assertions are collapsed into sorted evidence sets; `support_count` counts distinct evidence codes, never join multiplicity. `term_metadata.tsv.gz` has {:,} terms. `term_hierarchy.tsv.gz` has {:,} GO/Reactome child-parent edges. GO metadata includes `biological_process`, `molecular_function`, or `cellular_component` namespace and obsolete status. Only direct GO annotations are assigned to copies; hierarchy edges are frozen separately so a downstream method must opt in explicitly to ancestor propagation.".format(
            len(term_rows), len(metadata), len(edges)
        ),
        "",
        "GO assertions with a `NOT` qualifier are excluded. GO identifiers absent from the matching `go-basic` snapshot are also rejected: `{}`.".format(
            ", ".join("{}={}".format(k, v) for k, v in sorted(go_rejected.items()))
            or "none"
        ),
        "",
        "## Unmapped-copy and inflation diagnostics",
        "",
        "`copy_mapping_diagnostics.tsv.gz` contains exactly one row per physical locus, with stable mapping route and per-source term counts. `unmapped_copies.tsv.gz` emits one reason per missing functional source. The build fails on duplicate physical keys, a PHR copy outside the universe, ambiguous stable-ID indices, conflicting cross-references, HGNC group ID/name mismatches, missing term metadata, non-human Reactome rows, hierarchy endpoints absent from metadata, duplicate long-table keys, or output copy IDs outside the universe.",
        "",
        "## Hand-audited duplicated families",
        "",
        "Selection uses the exact, frozen display-name lists in `HAND_AUDIT`; it is not an annotation join. Up to four coordinate-distinct copies per listed name are written to `outputs/hand_audit_examples.tsv`. This deliberately includes pseudogenes and the unmapped rDNA control.",
        "",
        "Family | Rows | Distinct copies | Stable-mapped | Notes",
        "--- | ---: | ---: | ---: | ---",
    ]
    for family in HAND_AUDIT:
        rows = [x for x in audit if x["audit_family"] == family]
        stable = sum(x["mapping_status"].startswith("mapped_") for x in rows)
        note = (
            "45S loci lack stable GFF xrefs and correctly receive biotype only"
            if family.startswith("45S")
            else "coordinate-distinct loci retain independent copies; terms follow each locus's own ID"
        )
        lines.append("{} | {} | {} | {} | {}".format(family, len(rows), len({x['copy_id'] for x in rows}), stable, note))
    lines += [
        "",
        "## Frozen provenance",
        "",
        "The checked-in `sources/*.gz` objects are the analysis inputs; `sources/SOURCE_MANIFEST.tsv` records complete upstream-object and snapshot SHA-256 values, HTTP metadata, release/object identifiers, filters, retrieval time, and licenses. `fetch_sources.py --fetch` is an explicit renewal workflow and refuses to overwrite snapshots. Plain `fetch_sources.py` verifies the frozen checksums.",
        "",
        "Source | Release/object | License | Upstream",
        "--- | --- | --- | ---",
    ]
    for row in source_records:
        release = row["release"] + (("; " + row["object_id"]) if row["object_id"] else "")
        lines.append(
            "{} | {} | [{}]({}) | [download]({})".format(
                row["source"], release, row["license"], row["license_url"], row["url"]
            )
        )
    lines += [
        "",
        "The local GFF inputs, source releases, repository object IDs, SHA-256 values, licensing notes, and derivation provenance are frozen in `outputs/INPUT_MANIFEST.tsv`. The universe GFF version is also encoded in the audited filename (`chm13v2.0_RefSeq_Liftoff_v5.2`).",
        "",
        "## Rebuild and validation",
        "",
        "From repository root:",
        "",
        "```bash",
        "python3 paper_prep/_brainstorming/chm13_copy_enrichment/fetch_sources.py",
        "python3 paper_prep/_brainstorming/chm13_copy_enrichment/build_term_maps.py",
        "python3 -m unittest discover -s paper_prep/_brainstorming/chm13_copy_enrichment/tests -v",
        "```",
        "",
        "The build is standard-library-only and deterministic: rerunning against unchanged snapshots reproduces byte-identical gzip outputs (gzip `mtime=0`) and the same `validation_summary.json`.",
        "",
    ]
    (HERE / "TERM_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def build(gff_path, phr_gff_path):
    source_records = read_source_manifest_and_verify()
    loci = read_gff(gff_path)
    phr_loci = read_gff(phr_gff_path)
    universe = {x.copy_id: x for x in loci}
    phr_ids = {x.copy_id for x in phr_loci}
    missing_phr = phr_ids.difference(universe)
    if missing_phr:
        raise ValueError("PHR GFF has copies absent from universe: {}".format(sorted(missing_phr)[:3]))

    hgnc_records = load_hgnc()
    resolver = HGNCResolver(hgnc_records, load_hgnc_withdrawn())
    mappings = {locus.copy_id: map_locus(locus, resolver) for locus in loci}
    rejected_mappings = {
        copy_id: result
        for copy_id, result in mappings.items()
        if result.status.startswith("rejected_")
    }
    if rejected_mappings:
        examples = [
            "{}: {} ({})".format(copy_id, result.status, result.detail)
            for copy_id, result in sorted(rejected_mappings.items())[:5]
        ]
        raise ValueError(
            "ambiguous or conflicting stable mappings rejected; build stopped:\n"
            + "\n".join(examples)
        )
    hgnc_by_id = resolver.records
    go_terms, go_alt_ids, go_edges = parse_obo(SOURCE_DIR / "go-basic_2026-06-15.obo.gz")
    gene2go, go_rejected = load_gene2go(go_alt_ids, go_terms)
    gene2reactome, reactome_terms, reactome_edges = load_reactome()

    accumulator = TermAccumulator()
    group_metadata = {}
    biotype_metadata = {}
    for locus in loci:
        mapping = mappings[locus.copy_id]
        biotype_term = "BIOTYPE:" + locus.gene_biotype
        accumulator.add(locus.copy_id, "biotype", biotype_term, "has_gff_biotype", "GFF3")
        biotype_metadata[biotype_term] = {
            "source": "biotype",
            "term_id": biotype_term,
            "term_name": locus.gene_biotype,
            "namespace": "GFF3_gene_biotype",
            "is_obsolete": "false",
            "description": "Exact gene_biotype value in CHM13v2.0 RefSeq Liftoff v5.2",
        }
        route = "|".join(mapping.routes)
        if mapping.hgnc_id:
            for group_id, group_name in hgnc_groups(hgnc_by_id[mapping.hgnc_id]):
                term_id = "HGNC_GROUP:" + group_id
                accumulator.add(
                    locus.copy_id,
                    "HGNC_group",
                    term_id,
                    "member_of_gene_group",
                    "HGNC",
                    route,
                    mapping.hgnc_id,
                )
                existing = group_metadata.get(term_id)
                if existing and existing["term_name"] != group_name:
                    raise ValueError("HGNC group ID has conflicting names: {}".format(term_id))
                group_metadata[term_id] = {
                    "source": "HGNC_group",
                    "term_id": term_id,
                    "term_name": group_name,
                    "namespace": "gene_group_or_family",
                    "is_obsolete": "false",
                    "description": "HGNC gene group ID {}".format(group_id),
                }
        for gene_id in mapping.gene_ids:
            for term_id, evidence_set in gene2go.get(gene_id, {}).items():
                for evidence in evidence_set:
                    accumulator.add(
                        locus.copy_id,
                        "GO",
                        term_id,
                        "direct_gene_annotation",
                        evidence,
                        route or "direct_geneid",
                        "GeneID:" + gene_id,
                    )
            for term_id, evidence_set in gene2reactome.get(gene_id, {}).items():
                if term_id not in reactome_terms:
                    raise ValueError("Reactome mapping term absent from metadata: {}".format(term_id))
                for evidence in evidence_set:
                    accumulator.add(
                        locus.copy_id,
                        "Reactome",
                        term_id,
                        "participates_in_pathway_or_ancestor",
                        evidence,
                        route or "direct_geneid",
                        "GeneID:" + gene_id,
                    )

    term_rows = accumulator.rows()
    keys = [(x["copy_id"], x["source"], x["term_id"]) for x in term_rows]
    if len(keys) != len(set(keys)):
        raise AssertionError("duplicate copy/source/term keys after collapse")
    if set(x["copy_id"] for x in term_rows).difference(universe):
        raise AssertionError("term row points outside physical universe")

    used_go = {x["term_id"] for x in term_rows if x["source"] == "GO"}
    if used_go.difference(go_terms):
        raise AssertionError("GO term row lacks metadata")
    metadata_map = {}
    metadata_map.update(go_terms)
    metadata_map.update(reactome_terms)
    metadata_map.update(group_metadata)
    metadata_map.update(biotype_metadata)
    metadata = [metadata_map[k] for k in sorted(metadata_map, key=lambda x: (metadata_map[x]["source"], x))]
    edges = [
        {
            "source": source,
            "child_term_id": child,
            "parent_term_id": parent,
            "relation": relation,
        }
        for source, child, parent, relation in sorted(go_edges | reactome_edges)
    ]
    for edge in edges:
        if edge["child_term_id"] not in metadata_map or edge["parent_term_id"] not in metadata_map:
            raise AssertionError("hierarchy edge endpoint lacks term metadata")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    universe_rows = []
    diagnostics = []
    term_counts = collections.Counter((x["copy_id"], x["source"]) for x in term_rows)
    for locus in loci:
        mapping = mappings[locus.copy_id]
        universe_rows.append(
            {
                "copy_id": locus.copy_id,
                "seqid": locus.seqid,
                "start_1based": str(locus.start),
                "end_1based_inclusive": str(locus.end),
                "strand": locus.strand,
                "gff_id": locus.gff_id,
                "liftoff_copy_num_id": locus.liftoff_copy_num_id,
                "gene_name_display_only": locus.gene_name,
                "gene_biotype": locus.gene_biotype,
                "gff_stable_xrefs": stable_xref_display(locus),
                "in_phr_audited_target": "1" if locus.copy_id in phr_ids else "0",
            }
        )
        diagnostics.append(
            {
                "copy_id": locus.copy_id,
                "gene_name_display_only": locus.gene_name,
                "gene_biotype": locus.gene_biotype,
                "stable_mapping_status": mapping.status,
                "stable_mapping_detail": mapping.detail,
                "mapping_routes": "|".join(mapping.routes),
                "resolved_hgnc_id": mapping.hgnc_id,
                "resolved_gene_ids": "|".join(mapping.gene_ids),
                "biotype_term_count": str(term_counts[(locus.copy_id, "biotype")]),
                "go_term_count": str(term_counts[(locus.copy_id, "GO")]),
                "hgnc_group_term_count": str(term_counts[(locus.copy_id, "HGNC_group")]),
                "reactome_term_count": str(term_counts[(locus.copy_id, "Reactome")]),
            }
        )

    unmapped = []
    for locus in loci:
        mapping = mappings[locus.copy_id]
        for source in ("GO", "HGNC_group", "Reactome"):
            if term_counts[(locus.copy_id, source)]:
                continue
            if source == "HGNC_group":
                reason = (
                    "no_unique_hgnc_mapping"
                    if not mapping.hgnc_id
                    else "resolved_hgnc_record_has_no_gene_group"
                )
            elif not mapping.gene_ids:
                reason = "no_usable_stable_geneid"
            else:
                reason = "resolved_geneid_has_no_{}_annotation".format(source.lower())
            unmapped.append(
                {
                    "copy_id": locus.copy_id,
                    "gene_name_display_only": locus.gene_name,
                    "gene_biotype": locus.gene_biotype,
                    "source": source,
                    "reason": reason,
                    "stable_mapping_status": mapping.status,
                }
            )

    coverage = coverage_rows(loci, term_rows, mappings)
    audit = audit_rows(loci, mappings, term_counts)
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "copy_universe.tsv.gz", list(universe_rows[0]), universe_rows
    )
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "copy_to_term.tsv.gz", TERM_COLUMNS, term_rows
    )
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "term_metadata.tsv.gz",
        ["source", "term_id", "term_name", "namespace", "is_obsolete", "description"],
        metadata,
    )
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "term_hierarchy.tsv.gz",
        ["source", "child_term_id", "parent_term_id", "relation"],
        edges,
    )
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "copy_mapping_diagnostics.tsv.gz", list(diagnostics[0]), diagnostics
    )
    deterministic_gzip_dict_writer(
        OUTPUT_DIR / "unmapped_copies.tsv.gz", list(unmapped[0]), unmapped
    )
    write_tsv(
        OUTPUT_DIR / "coverage_by_biotype.tsv",
        ["source", "gene_biotype", "n_loci", "n_mapped_loci", "coverage_percent"],
        coverage,
    )
    write_tsv(
        OUTPUT_DIR / "hand_audit_examples.tsv",
        list(audit[0]),
        audit,
    )
    input_manifest = [
        {
            "input": "CHM13 physical universe GFF3",
            "path": str(gff_path.relative_to(REPO)),
            "source_release": "CHM13v2.0 RefSeq Liftoff v5.2; project import 2026-03-31",
            "object_id": "git-commit:340e61d50dc1e26d527aa2f593c726f63b974530",
            "bytes": str(gff_path.stat().st_size),
            "sha256": sha256(gff_path),
            "license": "project-vendored derived annotation; no additional license declared in file",
            "provenance": "all GFF3 gene features; filename and repository history pin the annotation build",
            "role": "all gene features define the background physical-copy universe",
        },
        {
            "input": "audited PHR target GFF3",
            "path": str(phr_gff_path.relative_to(REPO)),
            "source_release": "derived from pinned CHM13 GFF and PHR intervals; 2026-03-31",
            "object_id": "git-commit:7e2b9cd7dcebcb96f085436cf2c3d0f5f6555334",
            "bytes": str(phr_gff_path.stat().st_size),
            "sha256": sha256(phr_gff_path),
            "license": "project-derived table; inherits source-data terms",
            "provenance": "exact interval-intersection product committed by the PHR analysis",
            "role": "exact subset flag only; never used to assign terms",
        },
    ]
    write_tsv(
        OUTPUT_DIR / "INPUT_MANIFEST.tsv",
        [
            "input",
            "path",
            "source_release",
            "object_id",
            "bytes",
            "sha256",
            "license",
            "provenance",
            "role",
        ],
        input_manifest,
    )
    summary = {
        "schema_version": "1.0",
        "physical_loci": len(loci),
        "unique_copy_ids": len(universe),
        "phr_target_loci": len(phr_ids),
        "biotypes": dict(sorted(collections.Counter(x.gene_biotype for x in loci).items())),
        "mapping_statuses": dict(sorted(collections.Counter(x.status for x in mappings.values()).items())),
        "term_rows": len(term_rows),
        "term_rows_by_source": dict(sorted(collections.Counter(x["source"] for x in term_rows).items())),
        "term_metadata_rows": len(metadata),
        "hierarchy_edges": len(edges),
        "unmapped_diagnostic_rows": len(unmapped),
        "go_rejected_annotations": dict(sorted(go_rejected.items())),
        "validation": {
            "all_term_copy_ids_in_universe": True,
            "unique_copy_source_term_keys": True,
            "phr_target_is_subset": True,
            "all_term_ids_have_metadata": True,
            "all_hierarchy_endpoints_have_metadata": True,
            "symbol_annotation_joins": 0,
            "pseudogene_parent_inheritance_operations": 0,
        },
    }
    (OUTPUT_DIR / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(
        loci,
        len(phr_ids),
        coverage,
        mappings,
        term_rows,
        metadata,
        edges,
        audit,
        source_records,
        go_rejected,
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument("--phr-gff", type=Path, default=DEFAULT_PHR_GFF)
    args = parser.parse_args()
    summary = build(args.gff.resolve(), args.phr_gff.resolve())
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
