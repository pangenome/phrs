#!/usr/bin/env python3
"""Fetch and freeze the external sources used by build_term_maps.py.

The checked-in snapshots are the reproducible inputs.  This script is an
explicit renewal procedure: it downloads current upstream objects, records
the complete-object SHA-256 before filtering, and writes deterministic gzip
snapshots (mtime=0).  It never silently replaces an existing snapshot.
"""

import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import shutil
import tempfile
import urllib.request


HERE = Path(__file__).resolve().parent
SNAPSHOT_DIR = HERE / "sources"
USER_AGENT = "phrs-chm13-term-map/1.0 (https://github.com/pangenome/phrs)"

SOURCES = {
    "ncbi_gene2go": {
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
        "snapshot": "ncbi_gene2go_human_2026-07-13.tsv.gz",
        "release": "NCBI Gene daily export, 2026-07-13",
        "license": "NCBI molecular data; NCBI places no restrictions on use or distribution",
        "license_url": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
    },
    "hgnc_complete_set": {
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "snapshot": "hgnc_complete_set_2026-07-10.tsv.gz",
        "release": "HGNC complete set, upstream Last-Modified 2026-07-10",
        "license": "HGNC data release policy: no restrictions on access or use",
        "license_url": "https://www.genenames.org/download/",
    },
    "hgnc_withdrawn": {
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt",
        "snapshot": "hgnc_withdrawn_2026-07-10.tsv.gz",
        "release": "HGNC withdrawn reports, upstream Last-Modified 2026-07-10",
        "license": "HGNC data release policy: no restrictions on access or use",
        "license_url": "https://www.genenames.org/download/",
    },
    "go_basic": {
        # The OBO header pins data-version releases/2026-06-15.  The GO dated
        # PURL is published for go.owl, not this go-basic.obo derivative, so
        # the immutable checked-in snapshot and SHA-256 are the object pin.
        "url": "https://current.geneontology.org/ontology/go-basic.obo",
        "snapshot": "go-basic_2026-06-15.obo.gz",
        "release": "Gene Ontology 2026-06-15",
        "object_id": "http://purl.obolibrary.org/obo/go/releases/2026-06-15/go.owl",
        "license": "CC BY 4.0",
        "license_url": "https://geneontology.org/docs/go-citation-policy/",
    },
    "reactome_ncbi": {
        "url": "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
        "snapshot": "reactome_v96_ncbi_human_all_levels.tsv.gz",
        "release": "Reactome version 96 (2026-04-01; files refreshed 2026-06-19)",
        "object_id": "doi:10.5281/zenodo.19581589",
        "license": "CC0 1.0 (Reactome annotation files)",
        "license_url": "https://reactome.org/license",
    },
    "reactome_pathways": {
        "url": "https://reactome.org/download/current/ReactomePathways.txt",
        "snapshot": "reactome_v96_human_pathways.tsv.gz",
        "release": "Reactome version 96 (2026-04-01; files refreshed 2026-06-19)",
        "object_id": "doi:10.5281/zenodo.19581589",
        "license": "CC0 1.0 (Reactome annotation files)",
        "license_url": "https://reactome.org/license",
    },
    "reactome_relations": {
        "url": "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
        "snapshot": "reactome_v96_human_pathway_relations.tsv.gz",
        "release": "Reactome version 96 (2026-04-01; files refreshed 2026-06-19)",
        "object_id": "doi:10.5281/zenodo.19581589",
        "license": "CC0 1.0 (Reactome annotation files)",
        "license_url": "https://reactome.org/license",
    },
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
        return {
            "resolved_url": response.geturl(),
            "last_modified": response.headers.get("Last-Modified", ""),
            "etag": response.headers.get("ETag", ""),
        }


def deterministic_gzip_writer(path):
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, gz


def copy_to_deterministic_gzip(source, destination, source_is_gzip=False):
    count = 0
    raw_out, gz_out = deterministic_gzip_writer(destination)
    try:
        opener = gzip.open if source_is_gzip else open
        with opener(str(source), "rb") as inp:
            for line in inp:
                gz_out.write(line)
                count += 1
    finally:
        gz_out.close()
        raw_out.close()
    return count


def filter_gene2go(source, destination):
    raw_out, gz_out = deterministic_gzip_writer(destination)
    count = 0
    try:
        with gzip.open(str(source), "rt", encoding="utf-8") as inp:
            for line in inp:
                if line.startswith("#tax_id") or line.startswith("9606\t"):
                    gz_out.write(line.encode("utf-8"))
                    if line.startswith("9606\t"):
                        count += 1
    finally:
        gz_out.close()
        raw_out.close()
    return count


def filter_reactome(source, destination, species="Homo sapiens"):
    raw_out, gz_out = deterministic_gzip_writer(destination)
    count = 0
    try:
        with source.open("rt", encoding="utf-8") as inp:
            for line in inp:
                fields = line.rstrip("\n").split("\t")
                if fields and fields[-1] == species:
                    gz_out.write(line.encode("utf-8"))
                    count += 1
    finally:
        gz_out.close()
        raw_out.close()
    return count


def filter_reactome_relations(source, destination, pathway_snapshot):
    pathway_ids = set()
    with gzip.open(str(pathway_snapshot), "rt", encoding="utf-8") as inp:
        for line in inp:
            pathway_ids.add(line.split("\t", 1)[0])
    raw_out, gz_out = deterministic_gzip_writer(destination)
    count = 0
    try:
        with source.open("rt", encoding="utf-8") as inp:
            for line in inp:
                parent, child = line.rstrip("\n").split("\t")[:2]
                if parent in pathway_ids and child in pathway_ids:
                    gz_out.write(line.encode("utf-8"))
                    count += 1
    finally:
        gz_out.close()
        raw_out.close()
    return count


def make_snapshot(key, raw_path, destination):
    if key == "ncbi_gene2go":
        return filter_gene2go(raw_path, destination)
    if key in ("reactome_ncbi", "reactome_pathways"):
        return filter_reactome(raw_path, destination)
    if key == "reactome_relations":
        pathways = SNAPSHOT_DIR / SOURCES["reactome_pathways"]["snapshot"]
        return filter_reactome_relations(raw_path, destination, pathways)
    # Both the HGNC TSV and GO OBO endpoints are plain text.  We freeze their
    # exact bytes inside a deterministic gzip container.
    return copy_to_deterministic_gzip(raw_path, destination, source_is_gzip=False)


def write_manifest(records):
    path = SNAPSHOT_DIR / "SOURCE_MANIFEST.tsv"
    fields = [
        "source",
        "release",
        "object_id",
        "url",
        "resolved_url",
        "upstream_last_modified",
        "upstream_etag",
        "upstream_bytes",
        "upstream_sha256",
        "snapshot",
        "snapshot_filter",
        "snapshot_rows",
        "snapshot_bytes",
        "snapshot_sha256",
        "license",
        "license_url",
        "retrieved_utc",
    ]
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def fetch(retrieved_utc):
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    with tempfile.TemporaryDirectory(prefix="chm13-term-sources-") as temp:
        tempdir = Path(temp)
        order = [
            "ncbi_gene2go",
            "hgnc_complete_set",
            "hgnc_withdrawn",
            "go_basic",
            "reactome_ncbi",
            "reactome_pathways",
            "reactome_relations",
        ]
        for key in order:
            spec = SOURCES[key]
            destination = SNAPSHOT_DIR / spec["snapshot"]
            if destination.exists():
                raise SystemExit("refusing to replace existing snapshot: {}".format(destination))
            raw_path = tempdir / (key + ".raw")
            print("fetching {}".format(key), flush=True)
            headers = download(spec["url"], raw_path)
            rows = make_snapshot(key, raw_path, destination)
            filter_text = {
                "ncbi_gene2go": "tax_id == 9606; header retained",
                "reactome_ncbi": "species == Homo sapiens",
                "reactome_pathways": "species == Homo sapiens",
                "reactome_relations": "both endpoints in frozen human pathway set",
            }.get(key, "none; exact content recompressed with gzip mtime=0")
            records.append(
                {
                    "source": key,
                    "release": spec["release"],
                    "object_id": spec.get("object_id", ""),
                    "url": spec["url"],
                    "resolved_url": headers["resolved_url"],
                    "upstream_last_modified": headers["last_modified"],
                    "upstream_etag": headers["etag"],
                    "upstream_bytes": raw_path.stat().st_size,
                    "upstream_sha256": sha256(raw_path),
                    "snapshot": "sources/" + destination.name,
                    "snapshot_filter": filter_text,
                    "snapshot_rows": rows,
                    "snapshot_bytes": destination.stat().st_size,
                    "snapshot_sha256": sha256(destination),
                    "license": spec["license"],
                    "license_url": spec["license_url"],
                    "retrieved_utc": retrieved_utc,
                }
            )
    write_manifest(records)


def verify():
    manifest = SNAPSHOT_DIR / "SOURCE_MANIFEST.tsv"
    failures = []
    with manifest.open(encoding="utf-8") as inp:
        for row in csv.DictReader(inp, delimiter="\t"):
            path = HERE / row["snapshot"]
            if not path.exists():
                failures.append("missing {}".format(path))
            elif sha256(path) != row["snapshot_sha256"]:
                failures.append("checksum mismatch {}".format(path))
    if failures:
        raise SystemExit("\n".join(failures))
    print("verified all frozen source snapshot checksums")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="renew all snapshots")
    parser.add_argument(
        "--retrieved-utc",
        default="2026-07-13T00:00:00Z",
        help="fixed provenance timestamp written to the manifest",
    )
    args = parser.parse_args()
    if args.fetch:
        fetch(args.retrieved_utc)
    else:
        verify()


if __name__ == "__main__":
    main()
