#!/usr/bin/env python3
"""Create COPY_engine term collections from the frozen coordinate-keyed map.

The annotation build deliberately uses a coordinate-anchored ``copy_id`` while
the independently prepared analysis universe uses the GFF ``locus_id``.  This
bridge performs a strict one-to-one coordinate join; symbols never participate
in it.  Output rows retain both identifiers and the frozen source so every
engine result can be traced back to ``copy_to_term.tsv.gz``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple


COLLECTIONS = ("GO_BP", "GO_CC", "GO_MF", "HGNC_group", "Reactome", "biotype")
GO_NAMESPACE = {
    "biological_process": "GO_BP",
    "cellular_component": "GO_CC",
    "molecular_function": "GO_MF",
}


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    """Prefer a repository-relative record so manifests survive worktree moves."""
    here = Path(__file__).resolve().parent
    try:
        return str(path.resolve().relative_to(here))
    except ValueError:
        return str(path.resolve())


def read_rows(path: Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_gzip_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        text = io.TextIOWrapper(compressed, newline="")
        writer = csv.DictWriter(text, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})
        text.flush()
        text.detach()
        compressed.close()
    os.replace(temporary, path)


def build_crosswalk(loci_path: Path, universe_path: Path) -> Tuple[Dict[str, str], int]:
    """Return copy_id -> locus_id after a bijective physical-coordinate join."""
    locus_by_key: Dict[Tuple[str, str, str, str, str], str] = {}
    locus_ids = set()
    for row in read_rows(loci_path):
        key = (row["chromosome"], row["start1"], row["end1"], row["strand"], row["locus_id"])
        if key in locus_by_key or row["locus_id"] in locus_ids:
            raise ValueError("analysis locus universe is not one-to-one: %s" % (key,))
        locus_by_key[key] = row["locus_id"]
        locus_ids.add(row["locus_id"])

    crosswalk: Dict[str, str] = {}
    matched_loci = set()
    for row in read_rows(universe_path):
        key = (row["seqid"], row["start_1based"], row["end_1based_inclusive"],
               row["strand"], row["gff_id"])
        locus_id = locus_by_key.get(key)
        if locus_id is None:
            raise ValueError("frozen copy has no exact analysis-locus coordinate match: %s" % row["copy_id"])
        if row["copy_id"] in crosswalk or locus_id in matched_loci:
            raise ValueError("copy/locus physical join is not bijective: %s" % row["copy_id"])
        crosswalk[row["copy_id"]] = locus_id
        matched_loci.add(locus_id)
    if matched_loci != locus_ids:
        missing = sorted(locus_ids - matched_loci)[:5]
        raise ValueError("%d analysis loci lack a frozen copy; examples: %s" %
                         (len(locus_ids - matched_loci), missing))
    return crosswalk, len(locus_ids)


def build_collections(loci_path: Path, universe_path: Path, edges_path: Path,
                      metadata_path: Path, output: Path) -> Mapping[str, object]:
    crosswalk, physical_loci = build_crosswalk(loci_path, universe_path)
    metadata = {}
    for row in read_rows(metadata_path):
        key = (row["source"], row["term_id"])
        if key in metadata:
            raise ValueError("duplicate frozen term metadata: %s" % (key,))
        metadata[key] = row

    output.mkdir(parents=True, exist_ok=True)
    by_collection = {name: [] for name in COLLECTIONS}
    seen = set()
    for row in read_rows(edges_path):
        edge = (row["copy_id"], row["source"], row["term_id"])
        if edge in seen:
            raise ValueError("duplicate frozen copy/source/term edge: %s" % (edge,))
        seen.add(edge)
        if row["copy_id"] not in crosswalk:
            raise ValueError("term edge references unknown physical copy: %s" % row["copy_id"])
        meta = metadata.get((row["source"], row["term_id"]))
        if meta is None:
            raise ValueError("term edge lacks frozen metadata: %s" % ((row["source"], row["term_id"]),))
        if row["source"] == "GO":
            try:
                collection = GO_NAMESPACE[meta["namespace"]]
            except KeyError as exc:
                raise ValueError("unknown GO namespace for %s: %s" %
                                 (row["term_id"], meta["namespace"])) from exc
        elif row["source"] in {"HGNC_group", "Reactome", "biotype"}:
            collection = row["source"]
        else:
            raise ValueError("undeclared frozen term source: %s" % row["source"])
        by_collection[collection].append({
            "locus_id": crosswalk[row["copy_id"]],
            "term_id": row["term_id"],
            "term_name": meta["term_name"],
            "frozen_source": row["source"],
            "copy_id": row["copy_id"],
        })

    fields = ("locus_id", "term_id", "term_name", "frozen_source", "copy_id")
    outputs = {}
    for collection in COLLECTIONS:
        rows = sorted(by_collection[collection],
                      key=lambda x: (x["term_id"], x["locus_id"], x["copy_id"]))
        path = output / (collection + ".tsv.gz")
        write_gzip_tsv(path, fields, rows)
        outputs[collection] = {
            "path": portable_path(path), "rows": len(rows),
            "terms": len({row["term_id"] for row in rows}), "sha256": sha256(path),
        }

    manifest = {
        "schema_version": "chm13-engine-terms-v1",
        "join": "exact chromosome,start1,end1,strand,gff_id; symbols excluded",
        "physical_loci": physical_loci,
        "source_edges": len(seen),
        "input_checksums": {
            portable_path(path): sha256(path)
            for path in (loci_path, universe_path, edges_path, metadata_path)
        },
        "collections": outputs,
        "collection_rows_total": sum(value["rows"] for value in outputs.values()),
    }
    if manifest["collection_rows_total"] != len(seen):
        raise AssertionError("collection split lost or duplicated frozen edges")
    manifest_path = output / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def parser() -> argparse.ArgumentParser:
    here = Path(__file__).resolve().parent
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--loci", type=Path, default=here / "analysis_ready/chm13_gene_loci.tsv.gz")
    result.add_argument("--universe", type=Path, default=here / "outputs/copy_universe.tsv.gz")
    result.add_argument("--edges", type=Path, default=here / "outputs/copy_to_term.tsv.gz")
    result.add_argument("--metadata", type=Path, default=here / "outputs/term_metadata.tsv.gz")
    result.add_argument("--output", type=Path, default=here / "engine_terms")
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    manifest = build_collections(args.loci, args.universe, args.edges, args.metadata, args.output)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
