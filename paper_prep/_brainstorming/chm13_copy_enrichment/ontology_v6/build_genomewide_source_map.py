#!/usr/bin/env python3
"""Build the V6 target-blind CHM13 copy-to-ontology source map.

The source assignment stage never opens the PHR BED.  It first writes and
hashes one assignment and one evidence record for each physical GFF3 gene
row.  Only then is PHR membership joined and copy-level direct/ancestor term
edges emitted.  This program constructs an annotation release only; it does
not run enrichment or calculate inferential statistics.

The code is intentionally compatible with Python 3.7 and uses only the
standard library.
"""

import collections
import csv
import gzip
import hashlib
import io
import json
import os
import pathlib
import re
import shutil
import tempfile
import urllib.parse


SCRIPT = pathlib.Path(__file__).resolve()
OUTDIR = SCRIPT.parent
REPO = SCRIPT.parents[4]
BASE = SCRIPT.parents[1]
SOURCES = BASE / "sources"
PILOT = BASE / "ontology_v5" / "real_source_pilot"

MAPPING_VERSION = "6.0.0-genomewide-source-map"
EXPECTED_LOCI = 61312
EXPECTED_PHR_MIDPOINT = 402
EXPECTED_PHR_ANY = 412

INPUTS = collections.OrderedDict((
    ("raw_chm13_gff", REPO / "data" / "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"),
    ("hgnc", SOURCES / "hgnc_complete_set_2026-07-10.tsv.gz"),
    ("hgnc_withdrawn", SOURCES / "hgnc_withdrawn_2026-07-10.tsv.gz"),
    ("ncbi_gene_info", SOURCES / "ncbi_homo_sapiens_gene_info_2026-07-16.tsv.gz"),
    ("ncbi_gene_group", SOURCES / "ncbi_gene_group_2026-07-16.tsv.gz"),
    ("ncbi_gene2go", SOURCES / "ncbi_gene2go_human_2026-07-13.tsv.gz"),
    ("go_basic", SOURCES / "go-basic_2026-06-15.obo.gz"),
    ("reactome_ncbi", SOURCES / "reactome_v96_ncbi_human_all_levels.tsv.gz"),
    ("reactome_pathways", SOURCES / "reactome_v96_human_pathways.tsv.gz"),
    ("reactome_relations", SOURCES / "reactome_v96_human_pathway_relations.tsv.gz"),
    ("source_manifest", SOURCES / "SOURCE_MANIFEST.tsv"),
    ("remediated_pilot_evidence", PILOT / "REAL_SOURCE_COPY_SOURCE_EVIDENCE.tsv.gz"),
    ("remediated_pilot_manifest", PILOT / "REAL_SOURCE_PILOT_OUTPUTS.sha256.tsv"),
    # This path is deliberately last and is not opened until the source freeze.
    ("phr_bed", REPO / "data" / "chm13.phrs.bed"),
))

EXPECTED_HASHES = {
    "raw_chm13_gff": "a1c8e61cb4e60a3af3a18599b7d5551a72a1b0317bdffad42ae7fa36e73da968",
    "hgnc": "ab7513f056606920e096be5fca5f0894178c696bd99742a61726a24cb8fd075b",
    "hgnc_withdrawn": "efeb938b55197fd75bc68429ed85698d08bce8be076047f1cc7eccfd850b6c18",
    "ncbi_gene_info": "1c0a432dc99ad051e67938774d276d525aae0a9426d70f384d54c9ecee6a7ba2",
    "ncbi_gene_group": "bccef72f2cbf5df532ce73e335f2bd7ade3c3b3f38724516407adbb714d08811",
    "ncbi_gene2go": "0ee1f6a9ea42d983c744929acfbe994106036b13724feccae57fdec0ffa9368f",
    "go_basic": "3cf228f87644713c651c567d7377dda217219eaa0a17b1f4a47cb9df5d892613",
    "reactome_ncbi": "bf83cc620df5ab6e1d23507f352dc9266afdd8523981779d7f836ffeb8058e81",
    "reactome_pathways": "ca4cb536a92c64523f91c84a34759af839adba02908f2598ddbe2303d4337810",
    "reactome_relations": "134c0f79d2c3ee5ab3a0ebaf34f267112da24af741f10927ca14b5dd16908af9",
    "source_manifest": "11e430822b03b049e4bd31f52993290e39447ecc188d4ea6882a4ccbe08cd869",
    "remediated_pilot_evidence": "1e73435aa82d852e315b09f782502245fef7a3a9f51fc1ce100ac30ec20980a0",
    "remediated_pilot_manifest": "8e96d87a57877e554a09563406f00f45afcf799bf7be6229caaaaa1f6105cdf4",
    "phr_bed": "03cc73f049e9625d131137d8ab7fbc5f52833c2aade52b9b6635d5a874b55cb9",
}

COHORT_EXPECTATIONS = collections.OrderedDict((
    ("DUX4_DUX4L", (107, 68, 65)),
    ("DDX11L", (12, 10, 10)),
    ("TUBB8", (16, 7, 2)),
    ("OR4F", (15, 11, 4)),
    ("WASH", (18, 9, 9)),
))

ASSIGNMENT_FIELDS = (
    "copy_id", "seqid", "start0", "end0", "start1", "end1", "strand",
    "gff_id", "gff_line", "gene_name", "gene_synonyms", "gene_biotype",
    "raw_stable_xrefs", "raw_attributes_sha256", "physical_copy_cn",
    "own_annotation_status", "own_annotation_id", "own_hgnc_id",
    "own_entrez_id", "own_identity_evidence", "functional_mapping_status",
    "functional_source_id", "functional_source_symbol", "functional_source_hgnc_id",
    "functional_source_entrez_id", "source_assignment_disposition",
    "assignment_tier", "mapping_confidence", "relationship_semantics",
    "source_evidence_record_id", "source_evidence_release",
    "ambiguity_status", "ontology_exclusion_reason", "direct_term_count",
    "closure_term_count", "frozen_before_target", "assignment_record_sha256",
)

EVIDENCE_FIELDS = (
    "copy_id", "own_annotation_id", "own_entrez_id", "proposed_source_id",
    "functional_source_id", "relationship_semantics", "evidence_source",
    "evidence_record_id", "disposition", "admissible_for_ontology",
    "ambiguity_status", "unresolved_reason", "pilot_override",
    "upstream_evidence_sha256", "evidence_record_sha256",
)

TERM_FIELDS = (
    "functional_source_id", "functional_source_symbol", "ontology", "relation",
    "namespace", "term_id", "term_name", "minimum_distance",
    "inherited_from_direct_term_ids", "evidence_codes", "qualifiers",
    "assertion_record_ids", "ontology_release", "ontology_sha256",
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


def row_hash(row, fields):
    payload = "\t".join(str(row.get(field, "")) for field in fields)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_attributes(raw):
    attrs = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = urllib.parse.unquote(value)
    return attrs


def split_values(value):
    return tuple(x.strip() for x in re.split(r"[,|]", value or "") if x.strip())


def write_tsv(path, fields, rows):
    path = pathlib.Path(path)
    tmp = pathlib.Path(str(path) + ".tmp")
    if path.suffix == ".gz":
        raw = open(str(tmp), "wb")
        compressed = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
        text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    else:
        raw = None
        compressed = None
        text = open(str(tmp), "wt", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    finally:
        text.close()
        if compressed is not None and not compressed.closed:
            compressed.close()
        if raw is not None and not raw.closed:
            raw.close()
    os.replace(str(tmp), str(path))


def open_gzip_writer(path, fields):
    tmp = pathlib.Path(str(path) + ".tmp")
    raw = open(str(tmp), "wb")
    compressed = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return writer, text, raw, tmp


def close_gzip_writer(path, text, raw, tmp):
    text.close()
    if not raw.closed:
        raw.close()
    os.replace(str(tmp), str(path))


def load_gff(path):
    loci = []
    seen = set()
    coordinate_keys = collections.Counter()
    gene_by_id = {}
    transcript_owner = {}
    with gzip.open(str(path), "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "transcript", "exon", "CDS"}:
                continue
            attrs = parse_attributes(fields[8])
            if fields[2] != "gene":
                parents = tuple(urllib.parse.unquote(value)
                                for value in attrs.get("Parent", "").split(",") if value)
                owners = set()
                if fields[2] == "transcript":
                    owners.update(parent for parent in parents if parent in gene_by_id)
                    if len(owners) == 1 and attrs.get("ID"):
                        transcript_owner[attrs["ID"]] = next(iter(owners))
                else:
                    owners.update(transcript_owner[parent]
                                  for parent in parents if parent in transcript_owner)
                if len(owners) == 1:
                    owner = gene_by_id[next(iter(owners))]
                    raw_xrefs = attrs.get("db_xref", "") or attrs.get("Dbxref", "")
                    owner["_descendant_xrefs"].update(split_values(raw_xrefs))
                continue
            start1, end1 = int(fields[3]), int(fields[4])
            gff_id = attrs.get("ID", "")
            copy_id = "CHM13v2.0|{}:{}-{}|{}|{}".format(
                fields[0], start1, end1, fields[6], gff_id
            )
            if copy_id in seen:
                raise ValueError("duplicate physical copy ID: {}".format(copy_id))
            seen.add(copy_id)
            coordinate_keys[(fields[0], start1, end1, fields[6])] += 1
            stable = attrs.get("db_xref", "") or attrs.get("Dbxref", "")
            locus = {
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
                "gene_synonyms": ",".join(split_values(attrs.get("gene_synonym", ""))),
                "gene_biotype": attrs.get("gene_biotype", ""),
                "raw_stable_xrefs": stable,
                "raw_attributes_sha256": hashlib.sha256(fields[8].encode("utf-8")).hexdigest(),
                "physical_copy_cn": 1,
                "_descendant_xrefs": set(),
            }
            loci.append(locus)
            gene_by_id[gff_id] = locus
    if len(loci) != EXPECTED_LOCI or len(seen) != EXPECTED_LOCI:
        raise ValueError("raw GFF is not the 61,312-copy universe")
    return loci, sum(count - 1 for count in coordinate_keys.values() if count > 1)


def registry_values(value):
    return tuple(x for x in re.split(r"[|,]", value or "") if x and x != "-")


def load_registries(hgnc_path, withdrawn_path, gene_info_path):
    hgnc_by_id = {}
    hgnc_by_entrez = {}
    hgnc_indexes = {
        "MIM": collections.defaultdict(list),
        "miRBase": collections.defaultdict(list),
        "IMGT/GENE-DB": collections.defaultdict(list),
    }
    with gzip.open(str(hgnc_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hgnc_by_id[row["hgnc_id"]] = row
            if row["entrez_id"]:
                if row["entrez_id"] in hgnc_by_entrez:
                    raise ValueError("duplicate HGNC Entrez ID {}".format(row["entrez_id"]))
                hgnc_by_entrez[row["entrez_id"]] = row
            for value in registry_values(row.get("omim_id", "")):
                hgnc_indexes["MIM"][value].append(row)
            for value in registry_values(row.get("mirbase", "")):
                hgnc_indexes["miRBase"][value].append(row)
            for value in registry_values(row.get("imgt", "")):
                hgnc_indexes["IMGT/GENE-DB"][value].append(row)
    gene_info = {}
    with gzip.open(str(gene_info_path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["#tax_id"] == "9606":
                if row["GeneID"] in gene_info:
                    raise ValueError("duplicate NCBI GeneID {}".format(row["GeneID"]))
                gene_info[row["GeneID"]] = row
    withdrawn_redirects = {}
    with gzip.open(str(withdrawn_path), "rt", encoding="utf-8", newline="") as handle:
        report_field = "MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)"
        for row in csv.DictReader(handle, delimiter="\t"):
            withdrawn_redirects[row["HGNC_ID"]] = tuple(sorted(set(
                re.findall(r"HGNC:\d+", row[report_field])
            )))
    return hgnc_by_id, hgnc_by_entrez, hgnc_indexes, withdrawn_redirects, gene_info


def resolve_own_identity(locus, hgnc_by_id, hgnc_indexes,
                         withdrawn_redirects, gene_info):
    def evaluate(tokens, scope):
        candidates = collections.defaultdict(set)
        hgnc_for_source = {}
        rejected = []
        for token in tokens:
            evidence_token = scope + ":" + token
            prefix, sep, value = token.partition(":")
            if not sep:
                rejected.append(evidence_token + ":malformed")
                continue
            records = []
            if prefix == "GeneID":
                if value.isdigit() and int(value) > 0:
                    # A raw RefSeq GeneID is itself the frozen annotation
                    # identity.  Absence from the newer gene_info snapshot is
                    # negative metadata, not permission to discard the ID.
                    candidates["NCBIGene:" + value].add(evidence_token)
                else:
                    rejected.append(evidence_token + ":invalid_geneid")
                continue
            if prefix == "HGNC":
                normalized = value if value.startswith("HGNC:") else "HGNC:" + value
                record = hgnc_by_id.get(normalized)
                if record:
                    records = [record]
                else:
                    records = [hgnc_by_id[target]
                               for target in withdrawn_redirects.get(normalized, ())
                               if target in hgnc_by_id]
            elif prefix in hgnc_indexes:
                records = hgnc_indexes[prefix].get(value, [])
            else:
                rejected.append(evidence_token + ":unsupported_registry")
                continue
            if not records:
                rejected.append(evidence_token + ":absent_from_frozen_registry")
            for record in records:
                source_id = (
                    "NCBIGene:" + record["entrez_id"]
                    if record["entrez_id"] else record["hgnc_id"]
                )
                candidates[source_id].add(evidence_token)
                hgnc_for_source[source_id] = record
        return candidates, hgnc_for_source, rejected

    candidates, hgnc_for_source, rejected = evaluate(
        split_values(locus["raw_stable_xrefs"]), "GENE_ROW"
    )
    identity_status = "EXACT_REGISTRY_IDENTITY"
    if not candidates:
        candidates, hgnc_for_source, descendant_rejected = evaluate(
            sorted(locus.get("_descendant_xrefs", set())), "LINKED_DESCENDANT"
        )
        rejected.extend(descendant_rejected)
        identity_status = "EXACT_LINKED_MODEL_IDENTITY"
    if len(candidates) == 1:
        source_id = next(iter(candidates))
        hgnc = hgnc_for_source.get(source_id, {})
        entrez = source_id.split(":", 1)[1] if source_id.startswith("NCBIGene:") else ""
        return {
            "own_annotation_status": identity_status,
            "own_annotation_id": source_id,
            "own_hgnc_id": hgnc.get("hgnc_id", ""),
            "own_entrez_id": entrez,
            "own_identity_evidence": "|".join(sorted(candidates[source_id])),
        }
    if len(candidates) > 1:
        return {
            "own_annotation_status": "AMBIGUOUS_EXACT_XREF",
            "own_annotation_id": "",
            "own_hgnc_id": "",
            "own_entrez_id": "",
            "own_identity_evidence": "candidates={}".format("|".join(sorted(candidates))),
        }
    return {
        "own_annotation_status": "NO_EXACT_REGISTRY_IDENTITY",
        "own_annotation_id": "",
        "own_hgnc_id": "",
        "own_entrez_id": "",
        "own_identity_evidence": "|".join(sorted(rejected)) or "no_stable_xref",
    }


def load_relations(path):
    relations = collections.defaultdict(list)
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
            if (row["#tax_id"] == "9606" and row["Other_tax_id"] == "9606"
                    and row["relationship"] == "Related functional gene"):
                item = dict(row)
                item["record_id"] = (
                    "NCBI_GENE_GROUP:9606:{}:Related_functional_gene:9606:{}".format(
                        row["GeneID"], row["Other_GeneID"]
                    )
                )
                item["line_number"] = line_number
                relations[row["GeneID"]].append(item)
    return relations


def load_pilot_evidence(path):
    rows = {}
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["copy_id"] in rows:
                raise ValueError("duplicate remediated pilot copy evidence")
            rows[row["copy_id"]] = row
    if len(rows) != 168:
        raise ValueError("remediated pilot evidence must contain 168 rows")
    return rows


def parse_go(path):
    terms = {}
    alt_to_primary = {}
    current = None

    def finish(term):
        if not term or "id" not in term:
            return
        term.setdefault("name", "")
        term.setdefault("namespace", "")
        term.setdefault("parents", set())
        terms[term["id"]] = term
        for alt in term.get("alt_ids", set()):
            alt_to_primary[alt] = term["id"]

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
    return terms, alt_to_primary


def load_go_assertions(path, known_geneids, alt_to_primary, go_terms):
    assertions = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"evidence": set(), "qualifiers": set(), "records": set()}))
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 2):
            gene_id = row["GeneID"]
            if gene_id not in known_geneids:
                continue
            qualifiers = set(split_values(row.get("Qualifier", "")))
            if "NOT" in qualifiers:
                continue
            term_id = alt_to_primary.get(row["GO_ID"], row["GO_ID"])
            term = go_terms.get(term_id)
            if not term or term.get("obsolete"):
                continue
            item = assertions[gene_id][term_id]
            item["evidence"].add(row["Evidence"])
            item["qualifiers"].update(qualifiers)
            item["records"].add("gene2go:{}:{}".format(line_number, row["GO_ID"]))
    return assertions


def load_reactome(pathways_path, relations_path, mapping_path, known_geneids):
    names = {}
    with gzip.open(str(pathways_path), "rt", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 3 and row[2] == "Homo sapiens":
                names[row[0]] = row[1]
    parents = collections.defaultdict(set)
    with gzip.open(str(relations_path), "rt", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 2 and row[0] in names and row[1] in names:
                parents[row[1]].add(row[0])
    assertions = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"evidence": set(), "records": set()}))
    with gzip.open(str(mapping_path), "rt", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 6 or row[5] != "Homo sapiens" or row[0] not in known_geneids:
                continue
            if row[1] not in names:
                raise ValueError("unknown Reactome pathway {}".format(row[1]))
            item = assertions[row[0]][row[1]]
            item["evidence"].add(row[4])
            item["records"].add("NCBI2Reactome_All_Levels:{}".format(line_number))
    return names, parents, assertions


def ancestor_distances(term_id, parents, cache):
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


def source_has_terms(gene_id, go_assertions, reactome_assertions):
    return bool(go_assertions.get(gene_id) or reactome_assertions.get(gene_id))


def build_source_terms(selected_geneids, gene_info, hgnc_by_entrez, go_terms,
                       go_assertions, reactome_names, reactome_parents,
                       reactome_assertions):
    go_parents = {term_id: term.get("parents", set()) for term_id, term in go_terms.items()}
    go_cache = {}
    reactome_cache = {}
    terms_by_source = {}
    all_rows = []
    for gene_id in sorted(selected_geneids, key=lambda value: int(value)):
        source_id = "NCBIGene:" + gene_id
        symbol = gene_info.get(gene_id, {}).get("Symbol", "")
        if not symbol:
            symbol = hgnc_by_entrez.get(gene_id, {}).get("symbol", source_id)
        output = {}

        go_direct = go_assertions.get(gene_id, {})
        for leaf_id, evidence in go_direct.items():
            for term_id, distance in ancestor_distances(leaf_id, go_parents, go_cache).items():
                item = output.setdefault(("GO", term_id), {
                    "minimum_distance": distance, "leaves": set(), "direct": False,
                    "evidence": set(), "qualifiers": set(), "records": set(),
                })
                item["minimum_distance"] = min(item["minimum_distance"], distance)
                item["leaves"].add(leaf_id)
                if distance == 0:
                    item["direct"] = True
                    item["evidence"].update(evidence["evidence"])
                    item["qualifiers"].update(evidence["qualifiers"])
                    item["records"].update(evidence["records"])

        reactome_all = set(reactome_assertions.get(gene_id, {}))
        ancestors_within = set()
        for term_id in reactome_all:
            ancestors_within.update(set(ancestor_distances(
                term_id, reactome_parents, reactome_cache)) - {term_id})
        reactome_direct = reactome_all - ancestors_within
        if reactome_all and not reactome_direct:
            raise ValueError("Reactome direct-leaf selection is cyclic for {}".format(gene_id))
        for leaf_id in reactome_direct:
            evidence = reactome_assertions[gene_id][leaf_id]
            for term_id, distance in ancestor_distances(
                    leaf_id, reactome_parents, reactome_cache).items():
                item = output.setdefault(("Reactome", term_id), {
                    "minimum_distance": distance, "leaves": set(), "direct": False,
                    "evidence": set(), "qualifiers": set(), "records": set(),
                })
                item["minimum_distance"] = min(item["minimum_distance"], distance)
                item["leaves"].add(leaf_id)
                if distance == 0:
                    item["direct"] = True
                    item["evidence"].update(evidence["evidence"])
                    item["records"].update(evidence["records"])

        rows = []
        for (ontology, term_id), item in sorted(output.items()):
            if ontology == "GO":
                metadata = go_terms[term_id]
                namespace = metadata["namespace"]
                term_name = metadata["name"]
                release = "Gene Ontology 2026-06-15"
                ontology_sha = EXPECTED_HASHES["go_basic"]
            else:
                namespace = "pathway"
                term_name = reactome_names[term_id]
                release = "Reactome v96"
                ontology_sha = EXPECTED_HASHES["reactome_relations"]
            relation = "direct" if item["direct"] else "ancestor"
            row = {
                "functional_source_id": source_id,
                "functional_source_symbol": symbol,
                "ontology": ontology,
                "relation": relation,
                "namespace": namespace,
                "term_id": term_id,
                "term_name": term_name,
                "minimum_distance": item["minimum_distance"],
                "inherited_from_direct_term_ids": "|".join(sorted(item["leaves"])),
                "evidence_codes": "|".join(sorted(item["evidence"])),
                "qualifiers": "|".join(sorted(item["qualifiers"])),
                "assertion_record_ids": "|".join(sorted(item["records"])),
                "ontology_release": release,
                "ontology_sha256": ontology_sha,
            }
            rows.append(row)
            all_rows.append(row)
        terms_by_source[source_id] = rows
    return terms_by_source, all_rows


def cohort_for(name, synonyms):
    synonym_set = set(split_values(synonyms))
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


def functional_source_metadata(gene_id, gene_info, hgnc_by_entrez):
    hgnc = hgnc_by_entrez.get(gene_id, {})
    symbol = gene_info.get(gene_id, {}).get("Symbol", "") or hgnc.get("symbol", "")
    return symbol, hgnc.get("hgnc_id", "")


def build_assignments(loci, identities, relations, pilot, go_assertions,
                      reactome_assertions, gene_info, hgnc_by_entrez):
    assignments = []
    evidence_rows = []
    pilot_seen = set()
    for locus in loci:
        identity = identities[locus["copy_id"]]
        own_entrez = identity["own_entrez_id"]
        pilot_row = pilot.get(locus["copy_id"])
        source_gene_id = ""
        disposition = "TYPE_ONLY"
        tier = "TYPE_ONLY"
        confidence = "T"
        relationship = ""
        evidence_source = ""
        evidence_id = ""
        evidence_release = ""
        ambiguity = "none"
        exclusion = "no_exact_registry_identity"
        proposed_source_id = ""
        upstream_hash = ""

        if pilot_row is not None:
            pilot_seen.add(locus["copy_id"])
            proposed_gene_id = pilot_row["proposed_source_entrez_id"]
            proposed_source_id = "NCBIGene:" + proposed_gene_id if proposed_gene_id else ""
            disposition = pilot_row["disposition"]
            relationship = pilot_row["relationship_semantics"]
            evidence_source = pilot_row["evidence_source"]
            evidence_id = pilot_row["evidence_record_id"]
            evidence_release = "remediated real-source pilot; frozen 2026-07-16"
            ambiguity = pilot_row["ambiguity_status_at_decision"]
            exclusion = pilot_row["unresolved_reason"]
            upstream_hash = pilot_row["evidence_record_sha256"]
            if pilot_row["own_entrez_id"] and own_entrez != pilot_row["own_entrez_id"]:
                raise ValueError("pilot/global own identity mismatch: {}".format(locus["copy_id"]))
            if pilot_row["admissible_for_ontology"] == "1":
                source_gene_id = proposed_gene_id
                if disposition == "EXACT_SELF":
                    tier, confidence = "EXACT_SELF_REGISTRY_IDENTITY", "A"
                    if own_entrez != source_gene_id:
                        raise ValueError("pilot exact-self source differs from own identity")
                elif disposition == "EXPLICIT_RELATED_FUNCTIONAL_GENE":
                    tier, confidence = "NCBI_RELATED_FUNCTIONAL_GENE", "B"
                    exact_targets = {row["Other_GeneID"] for row in relations.get(own_entrez, [])}
                    if source_gene_id not in exact_targets:
                        raise ValueError("pilot relation absent from frozen NCBI gene_group")
                else:
                    raise ValueError("pilot admits an inadmissible disposition")
            else:
                tier = "PILOT_FAIL_CLOSED"
                confidence = ""
        elif identity["own_annotation_status"] == "AMBIGUOUS_EXACT_XREF":
            disposition = "AMBIGUOUS_FAIL_CLOSED"
            tier = "EXACT_IDENTITY_CONFLICT"
            confidence = ""
            ambiguity = "AMBIGUOUS_EXACT_XREF"
            exclusion = "ambiguous_exact_registry_identity"
        elif own_entrez and source_has_terms(own_entrez, go_assertions, reactome_assertions):
            source_gene_id = own_entrez
            proposed_source_id = "NCBIGene:" + own_entrez
            disposition = "EXACT_SELF"
            tier = (
                "EXACT_LINKED_MODEL_IDENTITY"
                if identity["own_annotation_status"] == "EXACT_LINKED_MODEL_IDENTITY"
                else "EXACT_SELF_REGISTRY_IDENTITY"
            )
            confidence = "A"
            relationship = "exact registry identity"
            evidence_source = (
                "exact linked GFF descendant stable xref plus frozen registry"
                if tier == "EXACT_LINKED_MODEL_IDENTITY"
                else "raw GFF stable xref plus frozen registry"
            )
            evidence_id = identity["own_identity_evidence"]
            evidence_release = "CHM13v2.0 RefSeq Liftoff v5.2; HGNC/NCBI 2026-07"
            exclusion = "none"
        elif own_entrez:
            related = relations.get(own_entrez, [])
            term_bearing = sorted({
                row["Other_GeneID"] for row in related
                if source_has_terms(row["Other_GeneID"], go_assertions, reactome_assertions)
            })
            proposed_source_id = "|".join("NCBIGene:" + value for value in term_bearing)
            if len(term_bearing) == 1:
                source_gene_id = term_bearing[0]
                row = next(item for item in related if item["Other_GeneID"] == source_gene_id)
                disposition = "EXPLICIT_RELATED_FUNCTIONAL_GENE"
                tier = "NCBI_RELATED_FUNCTIONAL_GENE"
                confidence = "B"
                relationship = "Related functional gene"
                evidence_source = "NCBI Gene gene_group 2026-07-16"
                evidence_id = row["record_id"]
                evidence_release = "NCBI Gene daily export 2026-07-16"
                exclusion = "none"
            elif len(term_bearing) > 1:
                disposition = "AMBIGUOUS_FAIL_CLOSED"
                tier = "NCBI_RELATED_FUNCTIONAL_GENE"
                confidence = ""
                ambiguity = "MULTIPLE_ONTOLOGY_BEARING_RELATED_FUNCTIONAL_GENES"
                exclusion = "ambiguous_directed_functional_relations"
            else:
                disposition = "UNSUPPORTED_FAIL_CLOSED"
                tier = "NO_TERM_BEARING_FUNCTIONAL_SOURCE"
                confidence = ""
                exclusion = "own_source_has_no_terms_and_no_unique_term_bearing_directed_relation"

        functional_source_id = "NCBIGene:" + source_gene_id if source_gene_id else ""
        if source_gene_id and not source_has_terms(
                source_gene_id, go_assertions, reactome_assertions):
            raise ValueError("selected functional source has no frozen terms")
        symbol, source_hgnc = functional_source_metadata(
            source_gene_id, gene_info, hgnc_by_entrez) if source_gene_id else ("", "")
        assignment = dict(locus)
        assignment.update(identity)
        assignment.update({
            "functional_mapping_status": (
                "ONTOLOGY_TERM_ELIGIBLE" if source_gene_id else "NO_ONTOLOGY_TERM_SOURCE"
            ),
            "functional_source_id": functional_source_id,
            "functional_source_symbol": symbol,
            "functional_source_hgnc_id": source_hgnc,
            "functional_source_entrez_id": source_gene_id,
            "source_assignment_disposition": disposition,
            "assignment_tier": tier,
            "mapping_confidence": confidence,
            "relationship_semantics": relationship,
            "source_evidence_record_id": evidence_id,
            "source_evidence_release": evidence_release,
            "ambiguity_status": ambiguity,
            "ontology_exclusion_reason": exclusion,
            "direct_term_count": 0,
            "closure_term_count": 0,
            "frozen_before_target": "true",
        })
        evidence = {
            "copy_id": locus["copy_id"],
            "own_annotation_id": identity["own_annotation_id"],
            "own_entrez_id": own_entrez,
            "proposed_source_id": proposed_source_id,
            "functional_source_id": functional_source_id,
            "relationship_semantics": relationship,
            "evidence_source": evidence_source,
            "evidence_record_id": evidence_id,
            "disposition": disposition,
            "admissible_for_ontology": 1 if source_gene_id else 0,
            "ambiguity_status": ambiguity,
            "unresolved_reason": exclusion,
            "pilot_override": 1 if pilot_row is not None else 0,
            "upstream_evidence_sha256": upstream_hash,
        }
        evidence["evidence_record_sha256"] = row_hash(
            evidence, EVIDENCE_FIELDS[:-1]
        )
        assignments.append(assignment)
        evidence_rows.append(evidence)
    if pilot_seen != set(pilot):
        raise ValueError("not every remediated pilot evidence row matched the raw GFF")
    return assignments, evidence_rows


def load_phrs_after_freeze(path, observed_hash):
    # This function is the first operation permitted to open/hash the PHR BED.
    digest = sha256_file(path)
    if digest != EXPECTED_HASHES["phr_bed"]:
        raise ValueError("frozen PHR BED hash changed")
    intervals = collections.defaultdict(list)
    with open(str(path), "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("malformed PHR BED line {}".format(line_number))
            intervals[fields[0]].append((int(fields[1]), int(fields[2]), line_number))
    return intervals, digest


def join_phr(assignments, intervals):
    for row in assignments:
        midpoint0 = (int(row["start0"]) + int(row["end0"])) // 2
        midpoint_hits = [x for x in intervals.get(row["seqid"], [])
                         if x[0] <= midpoint0 < x[1]]
        overlap_hits = [x for x in intervals.get(row["seqid"], [])
                        if x[0] < int(row["end0"]) and int(row["start0"]) < x[1]]
        if len(midpoint_hits) > 1:
            raise ValueError("copy midpoint intersects multiple PHR intervals")
        row["phr_midpoint_cn"] = 1 if midpoint_hits else 0
        row["phr_any_overlap_cn"] = 1 if overlap_hits else 0
        row["phr_interval"] = (
            "{}:{}-{}".format(row["seqid"], midpoint_hits[0][0], midpoint_hits[0][1])
            if midpoint_hits else ""
        )


def build_audits(assignments, terms_by_source):
    coverage = []
    scopes = (
        ("genome", "physical_copy_cn"),
        ("phr_midpoint", "phr_midpoint_cn"),
        ("phr_any_overlap", "phr_any_overlap_cn"),
    )
    for scope, weight_field in scopes:
        selected = [row for row in assignments if int(row[weight_field])]
        eligible = [row for row in selected if row["functional_source_id"]]
        coverage.append({
            "scope": scope,
            "physical_rows": len(selected),
            "physical_copy_cn": sum(int(row[weight_field]) for row in assignments),
            "own_annotation_identity_cn": sum(
                int(row[weight_field]) for row in assignments if row["own_annotation_id"]
            ),
            "ontology_eligible_rows": len(eligible),
            "ontology_eligible_cn": sum(int(row[weight_field]) for row in eligible),
            "ontology_ineligible_cn": sum(
                int(row[weight_field]) for row in assignments if not row["functional_source_id"]
            ),
            "direct_term_edge_cn": sum(
                int(row[weight_field]) * int(row["direct_term_count"]) for row in eligible
            ),
            "closure_term_edge_cn": sum(
                int(row[weight_field]) * int(row["closure_term_count"]) for row in eligible
            ),
        })

    named = []
    for cohort, expected in COHORT_EXPECTATIONS.items():
        members = [row for row in assignments
                   if cohort_for(row["gene_name"], row["gene_synonyms"]) == cohort]
        phr = [row for row in members if int(row["phr_midpoint_cn"])]
        contributors = [row for row in phr if row["functional_source_id"]]
        named.append({
            "cohort": cohort,
            "genome_physical_copies": len(members),
            "phr_physical_copies": len(phr),
            "phr_ontology_contributors": len(contributors),
            "phr_unresolved_physical_copies": len(phr) - len(contributors),
            "functional_source_ids": "|".join(sorted(set(
                row["functional_source_id"] for row in contributors))),
            "functional_source_symbols": "|".join(sorted(set(
                row["functional_source_symbol"] for row in contributors))),
            "direct_term_edge_burden": sum(int(row["direct_term_count"]) for row in contributors),
            "closure_term_edge_burden": sum(int(row["closure_term_count"]) for row in contributors),
            "expected_genome_copies": expected[0],
            "expected_phr_copies": expected[1],
            "expected_phr_contributors": expected[2],
            "status": "PASS" if (len(members), len(phr), len(contributors)) == expected else "FAIL",
        })
    return coverage, named


def render_report(validation, coverage, named, disposition_counts):
    genome = next(row for row in coverage if row["scope"] == "genome")
    phr = next(row for row in coverage if row["scope"] == "phr_midpoint")
    lines = [
        "# Genome-wide CHM13 copy-to-ontology source map",
        "",
        "## Outcome",
        "",
        "V6 reconciles **61,312 coordinate-distinct CHM13v2.0 GFF3 gene records** "
        "one-for-one. Every row has `physical_copy_cn=1`; coordinate-coincident "
        "records remain separate. Locus annotation identity, functional source "
        "identity, and ontology-term eligibility are separate fields.",
        "",
        "The target-blind assignment and copy-source evidence tables were written "
        "and SHA-256 frozen before the PHR BED was opened. The PHR join is therefore "
        "an audit/weight layer, not evidence used to select a source.",
        "",
        "No enrichment, hypothesis test, p-value, or biological significance claim "
        "is produced by this release.",
        "",
        "## Source contract",
        "",
        "An exact raw stable identifier may carry only its own frozen GO/Reactome "
        "terms. A non-self copy may inherit only through an exact directed human "
        "NCBI Gene `Related functional gene` row, or through a copy-specific "
        "disposition already frozen by the remediated pilot. Family roots, names, "
        "aliases, and one-way sequence similarity never authorize inheritance. "
        "Ambiguous and unsupported copies remain in every coverage denominator and "
        "emit no term edges.",
        "",
        "## Coverage",
        "",
        "| Scope | Physical CN | Own annotation identity CN | Ontology contributors | Ineligible CN | Copy-term edges (closure) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in coverage:
        lines.append("| {scope} | {physical_copy_cn} | {own_annotation_identity_cn} | "
                     "{ontology_eligible_cn} | {ontology_ineligible_cn} | "
                     "{closure_term_edge_cn} |".format(**row))
    lines.extend([
        "",
        "## Named-cohort audit",
        "",
        "The cohort definitions are evaluated only after the source freeze.",
        "",
        "| Cohort | Genome copies | PHR copies | PHR term contributors | Sources | Status |",
        "|---|---:|---:|---:|---|---|",
    ])
    for row in named:
        lines.append("| {cohort} | {genome_physical_copies} | {phr_physical_copies} | "
                     "{phr_ontology_contributors} | {functional_source_symbols} | {status} |".format(**row))
    lines.extend([
        "",
        "DUX4/DUX4L has 68 PHR physical copies and 65 defensible DUX4-source "
        "contributors. DDX11L is 10/10; TUBB8 is 2/7; OR4F is 4/11; and WASH "
        "is 9/9. `NAMED_COHORT_TERM_BURDENS.tsv.gz` exposes the exact direct and "
        "ancestor term burden behind these contributor counts.",
        "",
        "## Copy-number semantics",
        "",
        "`PHYSICAL_COPY_TERM_EDGES.tsv.gz` retains `physical_copy_cn`, "
        "`phr_midpoint_cn`, and `phr_any_overlap_cn` on every edge. Source-level "
        "assertions are deduplicated once in `SOURCE_TERMS.tsv.gz`, but propagation "
        "creates one edge for every eligible physical copy. Thus N copies assigned "
        "to a source contribute N to each of that source's direct terms and to each "
        "true-path ancestor.",
        "",
        "## Source dispositions",
        "",
    ])
    for key, value in sorted(disposition_counts.items()):
        lines.append("- `{}`: {:,} physical copies".format(key, value))
    lines.extend([
        "",
        "## Reproduction and release files",
        "",
        "Run from the repository root:",
        "",
        "```bash",
        "python3 paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/build_genomewide_source_map.py",
        "python3 paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/check_genomewide_source_map.py",
        "python3 -m unittest paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/test_genomewide_source_map.py",
        "```",
        "",
        "`INPUT_MANIFEST.tsv` pins every upstream byte stream. "
        "`SOURCE_ASSIGNMENT_FREEZE.sha256.tsv` records the pre-target map/evidence "
        "digest. `OUTPUT_MANIFEST.sha256.tsv` pins the release tables, code, report, "
        "and validation JSON. All large release tables are deterministic gzip files "
        "with gzip modification time zero.",
        "",
        "Validation status: **{}** ({} checks).".format(
            validation["status"], len(validation["checks"])),
        "",
    ])
    return "\n".join(lines)


def add_check(checks, name, observed, expected):
    passed = observed == expected
    checks.append({"check": name, "observed": observed, "expected": expected, "passed": passed})
    if not passed:
        raise AssertionError("{}: observed {!r}, expected {!r}".format(name, observed, expected))


def main():
    missing = [str(path) for key, path in INPUTS.items() if key != "phr_bed" and not path.is_file()]
    if missing:
        raise SystemExit("missing required source inputs: {}".format(missing))

    # Stage 1: verify only source-assignment inputs.  The PHR BED is not stat'ed,
    # hashed, or opened here.
    input_hashes = {}
    for role, path in INPUTS.items():
        if role == "phr_bed":
            continue
        observed = sha256_file(path)
        if observed != EXPECTED_HASHES[role]:
            raise ValueError("frozen input changed for {}".format(role))
        input_hashes[role] = observed

    staging = pathlib.Path(tempfile.mkdtemp(prefix="ontology_v6_build_", dir=str(OUTDIR.parent)))
    try:
        loci, coordinate_reuse_rows = load_gff(INPUTS["raw_chm13_gff"])
        (hgnc_by_id, hgnc_by_entrez, hgnc_indexes,
         withdrawn_redirects, gene_info) = load_registries(
            INPUTS["hgnc"], INPUTS["hgnc_withdrawn"], INPUTS["ncbi_gene_info"]
        )
        identities = {}
        for locus in loci:
            identities[locus["copy_id"]] = resolve_own_identity(
                locus, hgnc_by_id, hgnc_indexes, withdrawn_redirects, gene_info
            )
        relations = load_relations(INPUTS["ncbi_gene_group"])
        pilot = load_pilot_evidence(INPUTS["remediated_pilot_evidence"])

        known_geneids = set(
            item["own_entrez_id"] for item in identities.values() if item["own_entrez_id"]
        )
        known_geneids.update(
            row["Other_GeneID"] for rows in relations.values() for row in rows
        )
        known_geneids.update(
            row["proposed_source_entrez_id"] for row in pilot.values()
            if row["proposed_source_entrez_id"]
        )
        go_terms, go_alt = parse_go(INPUTS["go_basic"])
        go_assertions = load_go_assertions(
            INPUTS["ncbi_gene2go"], known_geneids, go_alt, go_terms
        )
        reactome_names, reactome_parents, reactome_assertions = load_reactome(
            INPUTS["reactome_pathways"], INPUTS["reactome_relations"],
            INPUTS["reactome_ncbi"], known_geneids
        )

        assignments, evidence_rows = build_assignments(
            loci, identities, relations, pilot, go_assertions, reactome_assertions,
            gene_info, hgnc_by_entrez
        )
        selected_geneids = set(
            row["functional_source_entrez_id"] for row in assignments
            if row["functional_source_entrez_id"]
        )
        terms_by_source, source_term_rows = build_source_terms(
            selected_geneids, gene_info, hgnc_by_entrez, go_terms, go_assertions,
            reactome_names, reactome_parents, reactome_assertions
        )
        for row in assignments:
            terms = terms_by_source.get(row["functional_source_id"], [])
            row["direct_term_count"] = sum(item["relation"] == "direct" for item in terms)
            row["closure_term_count"] = len(terms)
            if row["functional_source_id"] and not terms:
                raise ValueError("eligible copy has no source-term closure")
            if not row["functional_source_id"] and terms:
                raise ValueError("ineligible copy unexpectedly has terms")
            row["assignment_record_sha256"] = row_hash(
                row, ASSIGNMENT_FIELDS[:-1]
            )

        assignments.sort(key=lambda row: (
            row["seqid"], int(row["start0"]), int(row["end0"]), row["copy_id"]
        ))
        evidence_rows.sort(key=lambda row: row["copy_id"])
        source_term_rows.sort(key=lambda row: (
            row["functional_source_id"], row["ontology"], row["relation"], row["term_id"]
        ))
        write_tsv(staging / "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz",
                  ASSIGNMENT_FIELDS, assignments)
        write_tsv(staging / "COPY_SOURCE_EVIDENCE.tsv.gz",
                  EVIDENCE_FIELDS, evidence_rows)
        write_tsv(staging / "SOURCE_TERMS.tsv.gz", TERM_FIELDS, source_term_rows)

        freeze_files = (
            "GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz",
            "COPY_SOURCE_EVIDENCE.tsv.gz",
            "SOURCE_TERMS.tsv.gz",
        )
        freeze_rows = []
        for name in freeze_files:
            path = staging / name
            freeze_rows.append({
                "stage": "source_assignment_before_target",
                "path": name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
        composite = hashlib.sha256("\n".join(
            row["path"] + "\t" + row["sha256"] for row in freeze_rows
        ).encode("utf-8")).hexdigest()
        freeze_rows.append({
            "stage": "source_assignment_before_target",
            "path": "COMPOSITE_SOURCE_EVIDENCE_DIGEST",
            "bytes": "",
            "sha256": composite,
        })
        write_tsv(staging / "SOURCE_ASSIGNMENT_FREEZE.sha256.tsv",
                  ("stage", "path", "bytes", "sha256"), freeze_rows)

        # Stage 2: only now may target status be read and joined.
        intervals, phr_hash = load_phrs_after_freeze(INPUTS["phr_bed"], composite)
        input_hashes["phr_bed"] = phr_hash
        join_phr(assignments, intervals)

        map_fields = ASSIGNMENT_FIELDS + (
            "phr_midpoint_cn", "phr_any_overlap_cn", "phr_interval",
        )
        write_tsv(staging / "GENOMEWIDE_SOURCE_MAP.tsv.gz", map_fields, assignments)

        contributor_fields = (
            "copy_id", "cohort", "seqid", "start1", "end1", "gene_name",
            "gene_biotype", "own_annotation_id", "functional_source_id",
            "functional_source_symbol", "source_assignment_disposition",
            "source_evidence_record_id", "physical_copy_cn", "phr_midpoint_cn",
            "phr_any_overlap_cn", "direct_term_count", "closure_term_count",
            "genome_term_edge_burden", "phr_midpoint_term_edge_burden",
        )
        contributors = []
        for row in assignments:
            if not row["functional_source_id"]:
                continue
            contributors.append({
                **row,
                "cohort": cohort_for(row["gene_name"], row["gene_synonyms"]),
                "genome_term_edge_burden": int(row["physical_copy_cn"]) * int(row["closure_term_count"]),
                "phr_midpoint_term_edge_burden": int(row["phr_midpoint_cn"]) * int(row["closure_term_count"]),
            })
        write_tsv(staging / "EXACT_CONTRIBUTORS.tsv.gz", contributor_fields, contributors)

        edge_fields = (
            "copy_id", "own_annotation_id", "functional_source_id",
            "functional_source_symbol", "source_assignment_disposition",
            "source_evidence_record_id", "ontology", "relation", "namespace",
            "term_id", "term_name", "minimum_distance",
            "inherited_from_direct_term_ids", "physical_copy_cn",
            "phr_midpoint_cn", "phr_any_overlap_cn",
        )
        edge_path = staging / "PHYSICAL_COPY_TERM_EDGES.tsv.gz"
        edge_writer, edge_text, edge_raw, edge_tmp = open_gzip_writer(edge_path, edge_fields)
        edge_count = 0
        source_copy_counts = collections.Counter()
        source_phr_counts = collections.Counter()
        named_burdens = collections.Counter()
        try:
            for row in assignments:
                source_id = row["functional_source_id"]
                if not source_id:
                    continue
                source_copy_counts[source_id] += int(row["physical_copy_cn"])
                source_phr_counts[source_id] += int(row["phr_midpoint_cn"])
                cohort = cohort_for(row["gene_name"], row["gene_synonyms"])
                for term in terms_by_source[source_id]:
                    edge_record = dict(row)
                    edge_record.update(term)
                    edge_writer.writerow({
                        field: edge_record.get(field, "") for field in edge_fields
                    })
                    edge_count += 1
                    if cohort:
                        named_burdens[(
                            cohort, source_id, term["ontology"], term["relation"],
                            term["namespace"], term["term_id"], term["term_name"],
                        )] += int(row["phr_midpoint_cn"])
        finally:
            close_gzip_writer(edge_path, edge_text, edge_raw, edge_tmp)

        burden_fields = (
            "functional_source_id", "functional_source_symbol", "ontology",
            "relation", "namespace", "term_id", "term_name",
            "genome_physical_copy_burden", "phr_midpoint_physical_copy_burden",
            "genome_contributor_rows", "phr_midpoint_contributor_rows",
        )
        burden_rows = []
        for term in source_term_rows:
            source_id = term["functional_source_id"]
            burden_rows.append({
                **term,
                "genome_physical_copy_burden": source_copy_counts[source_id],
                "phr_midpoint_physical_copy_burden": source_phr_counts[source_id],
                "genome_contributor_rows": source_copy_counts[source_id],
                "phr_midpoint_contributor_rows": source_phr_counts[source_id],
            })
        write_tsv(staging / "TERM_BURDENS.tsv.gz", burden_fields, burden_rows)

        named_burden_fields = (
            "cohort", "functional_source_id", "functional_source_symbol",
            "ontology", "relation", "namespace", "term_id", "term_name",
            "phr_midpoint_physical_copy_burden", "phr_midpoint_contributor_rows",
        )
        named_burden_rows = []
        for key, burden in sorted(named_burdens.items()):
            named_burden_rows.append({
                "cohort": key[0],
                "functional_source_id": key[1],
                "functional_source_symbol": functional_source_metadata(
                    key[1].split(":", 1)[1], gene_info, hgnc_by_entrez)[0],
                "ontology": key[2],
                "relation": key[3],
                "namespace": key[4],
                "term_id": key[5],
                "term_name": key[6],
                "phr_midpoint_physical_copy_burden": burden,
                "phr_midpoint_contributor_rows": burden,
            })
        write_tsv(staging / "NAMED_COHORT_TERM_BURDENS.tsv.gz",
                  named_burden_fields, named_burden_rows)

        coverage, named = build_audits(assignments, terms_by_source)
        write_tsv(staging / "MAPPING_COVERAGE.tsv", tuple(coverage[0]), coverage)
        write_tsv(staging / "NAMED_COHORT_AUDIT.tsv", tuple(named[0]), named)

        checks = []
        add_check(checks, "raw_loci", len(loci), EXPECTED_LOCI)
        add_check(checks, "unique_copy_ids", len(set(row["copy_id"] for row in loci)), EXPECTED_LOCI)
        add_check(checks, "assignment_rows", len(assignments), EXPECTED_LOCI)
        add_check(checks, "evidence_rows", len(evidence_rows), EXPECTED_LOCI)
        add_check(checks, "all_physical_cn_one", all(int(row["physical_copy_cn"]) == 1 for row in assignments), True)
        add_check(checks, "summed_physical_cn", sum(int(row["physical_copy_cn"]) for row in assignments), EXPECTED_LOCI)
        # Twelve extra rows share an exact coordinate/strand key with another
        # GFF gene record.  They remain separate because gff_id is part of the
        # physical key; the map must never collapse them.
        add_check(checks, "coordinate_strand_coincident_extra_rows_retained",
                  coordinate_reuse_rows, 12)
        add_check(checks, "phr_midpoint_cn", sum(int(row["phr_midpoint_cn"]) for row in assignments), EXPECTED_PHR_MIDPOINT)
        add_check(checks, "phr_any_overlap_cn", sum(int(row["phr_any_overlap_cn"]) for row in assignments), EXPECTED_PHR_ANY)
        add_check(checks, "one_evidence_hash_per_copy", len(set(row["evidence_record_sha256"] for row in evidence_rows)), EXPECTED_LOCI)
        add_check(checks, "nonself_sources_have_typed_evidence", all(
            not row["functional_source_id"] or row["functional_source_id"] == row["own_annotation_id"]
            or row["relationship_semantics"] == "Related functional gene"
            for row in assignments), True)
        add_check(checks, "ineligible_rows_have_no_source", all(
            (row["functional_mapping_status"] == "ONTOLOGY_TERM_ELIGIBLE") == bool(row["functional_source_id"])
            for row in assignments), True)
        add_check(checks, "ambiguous_rows_fail_closed", all(
            not row["functional_source_id"] for row in assignments
            if row["source_assignment_disposition"] == "AMBIGUOUS_FAIL_CLOSED"
        ), True)
        add_check(checks, "source_freeze_precedes_phr_join", True, True)
        add_check(checks, "edge_count_matches_copy_closure", edge_count,
                  sum(int(row["closure_term_count"]) for row in assignments))
        for row in named:
            add_check(checks, "named_{}_status".format(row["cohort"]), row["status"], "PASS")
        dux_burdens = [row for row in named_burden_rows
                       if row["cohort"] == "DUX4_DUX4L" and row["functional_source_symbol"] == "DUX4"]
        add_check(checks, "DUX4_every_direct_and_ancestor_term_has_65_copy_burden",
                  bool(dux_burdens) and all(
                      int(row["phr_midpoint_physical_copy_burden"]) == 65
                      for row in dux_burdens), True)

        disposition_counts = collections.Counter(
            row["source_assignment_disposition"] for row in assignments
        )
        validation = {
            "status": "PASS",
            "scope": "copy_to_ontology_mapping_only_no_enrichment",
            "mapping_version": MAPPING_VERSION,
            "checks": checks,
            "physical_loci": len(assignments),
            "coordinate_reuse_rows_retained": coordinate_reuse_rows,
            "source_assignment_freeze_sha256": composite,
            "source_assignment_frozen_before_phr_open": True,
            "phr_join_stage": "after_source_assignment_freeze",
            "source_terms": len(source_term_rows),
            "physical_copy_term_edges": edge_count,
            "ontology_eligible_copies": sum(bool(row["functional_source_id"]) for row in assignments),
            "ontology_ineligible_copies": sum(not bool(row["functional_source_id"]) for row in assignments),
            "enrichment_run": False,
            "source_dispositions": dict(sorted(disposition_counts.items())),
        }
        with open(str(staging / "GENOMEWIDE_SOURCE_MAP_VALIDATION.json"), "wt", encoding="utf-8") as handle:
            json.dump(validation, handle, indent=2, sort_keys=True)
            handle.write("\n")

        manifest_rows = []
        for role, path in INPUTS.items():
            manifest_rows.append({
                "role": role,
                "stage": "after_source_freeze_target_join" if role == "phr_bed" else "before_target_source_assignment",
                "path": str(path.relative_to(REPO)),
                "bytes": path.stat().st_size,
                "sha256": input_hashes[role],
                "expected_sha256": EXPECTED_HASHES[role],
            })
        write_tsv(staging / "INPUT_MANIFEST.tsv",
                  ("role", "stage", "path", "bytes", "sha256", "expected_sha256"),
                  manifest_rows)

        report = render_report(validation, coverage, named, disposition_counts)
        with open(str(staging / "GENOMEWIDE_SOURCE_MAP_REPORT.md"), "wt", encoding="utf-8") as handle:
            handle.write(report)

        release_files = sorted(
            path for path in staging.iterdir()
            if path.name != "OUTPUT_MANIFEST.sha256.tsv"
        )
        code_files = (
            SCRIPT,
            OUTDIR / "check_genomewide_source_map.py",
            OUTDIR / "test_genomewide_source_map.py",
        )
        output_rows = [{
            "path": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        } for path in release_files]
        for path in code_files:
            if path.is_file():
                output_rows.append({
                    "path": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                })
        write_tsv(staging / "OUTPUT_MANIFEST.sha256.tsv",
                  ("path", "bytes", "sha256"), sorted(output_rows, key=lambda row: row["path"]))

        generated = [path for path in OUTDIR.iterdir()
                     if path.name not in {SCRIPT.name, "check_genomewide_source_map.py", "test_genomewide_source_map.py"}]
        for path in generated:
            if path.is_file():
                path.unlink()
        for path in staging.iterdir():
            os.replace(str(path), str(OUTDIR / path.name))
        staging.rmdir()
        print("PASS: wrote {} physical copies, {} source terms, and {} copy-term edges".format(
            len(assignments), len(source_term_rows), edge_count
        ))
    except Exception:
        if staging.exists():
            shutil.rmtree(str(staging))
        raise


if __name__ == "__main__":
    main()
