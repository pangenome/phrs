#!/usr/bin/env python3
"""Independent black-box validator for a completed COPY_engine run.

This module intentionally does not import COPY_engine.  It reconstructs every
saved interval from the placement TSVs, performs brute-force coordinate joins,
and compares the resulting burden, denominator, copy-count, and breadth arrays
with the engine caches.  It then independently recomputes empirical p-values,
BH families, collection maxT, and global maxT from those caches.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def rows(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def resolve(recorded, run_dir):
    path = Path(recorded)
    return path if path.is_absolute() else run_dir / path


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    result = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse in range(len(values) - 1, -1, -1):
        index = int(order[reverse])
        running = min(running, values[index] * len(values) / (reverse + 1))
        result[index] = running
    return result


def standardized(observed, null):
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    zobs = np.zeros(len(observed))
    znull = np.zeros_like(null, dtype=float)
    degenerate = np.zeros(len(observed), dtype=bool)
    for column in range(len(observed)):
        values = null[:, column]
        if np.all(values == observed[column]):
            degenerate[column] = True
        elif np.all(values == values[0]) and observed[column] > values[0]:
            zobs[column] = np.inf
        else:
            pooled = np.concatenate(([observed[column]], values))
            sd = float(np.std(pooled, ddof=1))
            if sd == 0:
                degenerate[column] = True
            else:
                mean = float(np.mean(pooled))
                zobs[column] = (observed[column] - mean) / sd
                znull[:, column] = (values - mean) / sd
    return zobs, znull, degenerate


def selected_loci(loci, intervals, assignment):
    selected = set()
    for index, locus in enumerate(loci):
        for arm, start, end in intervals:
            if locus["arm"] != arm:
                continue
            if assignment == "midpoint":
                hit = start <= locus["midpoint"] < end
            else:
                hit = locus["start"] < end and start < locus["end"]
            if hit:
                selected.add(index)
                break
    return sorted(selected)


def recount(loci, selected, term_sets, term_ids):
    counts = np.zeros(len(term_ids), dtype=np.int64)
    arms = [set() for _ in term_ids]
    annotated = 0
    position = {term: index for index, term in enumerate(term_ids)}
    for locus_index in selected:
        locus_id = loci[locus_index]["locus_id"]
        terms = term_sets.get(locus_id, ())
        retained = [term for term in terms if term in position]
        if retained:
            annotated += 1
        for term in retained:
            index = position[term]
            counts[index] += 1
            arms[index].add(loci[locus_index]["arm"])
    return counts, annotated, np.asarray([len(value) for value in arms], dtype=np.int64)


def array_path(run_dir, assignment, collection, statistic, batch):
    name = "%s.%s.%s.%09d-%09d.npy" % (
        assignment, collection, statistic, int(batch["start"]), int(batch["end"]))
    return run_dir / "batches" / name


def validate(run_dir: Path, require_frozen_trace=True):
    run_dir = run_dir.resolve()
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    config = manifest["immutable_configuration"]
    input_paths = {name: resolve(value["path"], run_dir)
                   for name, value in config["inputs"].items()}
    loci = []
    locus_position = {}
    for row in rows(input_paths["loci"]):
        if row["locus_id"] in locus_position:
            raise AssertionError("duplicate physical locus_id")
        locus_position[row["locus_id"]] = len(loci)
        loci.append({"locus_id": row["locus_id"], "arm": row["arm"],
                     "start": int(row["start0"]), "end": int(row["end0"]),
                     "midpoint": int(row["midpoint0"])})

    term_filter = defaultdict(list)
    for row in rows(run_dir / "term_filtering.tsv"):
        if row["status"] == "retained":
            term_filter[row["collection"]].append(row["term_id"])
    collections = {}
    trace_errors = []
    for item in config["collections"]:
        name, path = item["name"], resolve(item["path"], run_dir)
        edges = defaultdict(set)
        seen = set()
        for row in rows(path):
            edge = (row["locus_id"], row["term_id"])
            if edge in seen:
                raise AssertionError("term-map join inflation/duplicate edge: %s %s" % (name, edge))
            seen.add(edge)
            if row["locus_id"] not in locus_position:
                raise AssertionError("term edge outside physical universe: %s" % (edge,))
            if require_frozen_trace and (not row.get("copy_id") or not row.get("frozen_source")):
                trace_errors.append((name, edge))
            edges[row["locus_id"]].add(row["term_id"])
        if trace_errors:
            raise AssertionError("term rows lack copy_id/frozen_source trace fields; first: %s" %
                                 (trace_errors[0],))
        collections[name] = {"terms": sorted(term_filter[name]), "edges": edges, "path": str(path)}

    cached = {}
    placement_rows = 0
    for batch in manifest["batches"]:
        by_replicate = defaultdict(list)
        for row in rows(run_dir / "batches" / batch["placements"]):
            placement_rows += 1
            components = [part.split(":") for part in row["components"].split(";")]
            start = int(row["block_start0"])
            for _identifier, offset, width in components:
                by_replicate[int(row["replicate"])].append(
                    (row["destination_arm"], start + int(offset), start + int(offset) + int(width)))
        expected_replicates = list(range(int(batch["start"]), int(batch["end"]) + 1))
        if sorted(by_replicate) != expected_replicates:
            raise AssertionError("placement batch has missing/duplicate replicate groups")
        for assignment in config["assignments"]:
            burdens = []
            counts_by_collection = {name: [] for name in collections}
            denominators = {name: [] for name in collections}
            breadth = {name: [] for name in collections}
            for replicate in expected_replicates:
                selected = selected_loci(loci, by_replicate[replicate], assignment)
                burdens.append(len(selected))
                for name, collection in collections.items():
                    count, denominator, arm_count = recount(
                        loci, selected, collection["edges"], collection["terms"])
                    counts_by_collection[name].append(count)
                    denominators[name].append(denominator)
                    breadth[name].append(arm_count)
            np.testing.assert_array_equal(
                np.load(array_path(run_dir, assignment, "all", "burden", batch)), burdens)
            cached.setdefault((assignment, "all", "burden"), []).append(np.asarray(burdens))
            for name in collections:
                for statistic, values in (("counts", counts_by_collection[name]),
                                          ("denominators", denominators[name]),
                                          ("breadth", breadth[name])):
                    value = np.asarray(values)
                    np.testing.assert_array_equal(
                        np.load(array_path(run_dir, assignment, name, statistic, batch)), value)
                    cached.setdefault((assignment, name, statistic), []).append(value)

    cached = {key: np.concatenate(value, axis=0) for key, value in cached.items()}
    blocks = []
    for row in rows(run_dir / "placement_blocks.tsv"):
        for _identifier, offset, width in (part.split(":") for part in row["components"].split(";")):
            start = int(row["source_start0"]) + int(offset)
            blocks.append((row["source_arm"], start, start + int(width)))
    observed_selected = {assignment: selected_loci(loci, blocks, assignment)
                         for assignment in config["assignments"]}

    result_rows = list(rows(run_dir / "term_results.tsv"))
    families = defaultdict(list)
    for row in result_rows:
        key = (row["assignment"], row["collection"], row["statistic"])
        families[key].append(row)
        collection = collections[row["collection"]]
        terms = collection["terms"]
        term_index = terms.index(row["term_id"])
        count, denominator, arm_count = recount(
            loci, observed_selected[row["assignment"]], collection["edges"], terms)
        statistic = row["statistic"]
        observed = (count[term_index] if statistic == "copy_burden" else
                    count[term_index] / denominator if statistic == "composition" and denominator else
                    arm_count[term_index])
        if not math.isclose(float(row["observed"]), float(observed), abs_tol=1e-12):
            raise AssertionError("observed recount mismatch: %s/%s" % (row["collection"], row["term_id"]))
        if not (0 <= int(row["observed_copy_count"]) <= int(row["observed_denominator"]) <=
                len(observed_selected[row["assignment"]])):
            raise AssertionError("impossible observed contingency fields")
        null = cached[(row["assignment"], row["collection"],
                       "counts" if statistic != "breadth" else "breadth")][:, term_index].astype(float)
        if statistic == "composition":
            denom = cached[(row["assignment"], row["collection"], "denominators")]
            null = np.divide(null, denom, out=np.zeros_like(null), where=denom > 0)
        exceed = len(null) if observed == 0 else int(np.count_nonzero(null >= observed))
        pvalue = (exceed + 1) / (len(null) + 1)
        if int(row["exceedances"]) != exceed or not math.isclose(float(row["raw_p"]), pvalue):
            raise AssertionError("empirical tail mismatch")

    global_maxima = np.full(int(manifest["completed_permutations"]), -np.inf)
    for key, family in families.items():
        assignment, name, statistic = key
        collection = collections[name]
        terms = collection["terms"]
        family.sort(key=lambda row: terms.index(row["term_id"]))
        reported_p = np.asarray([float(row["raw_p"]) for row in family])
        np.testing.assert_allclose([float(row["bh_q"]) for row in family], bh(reported_p), atol=1e-12)
        count, denominator, arm_count = recount(
            loci, observed_selected[assignment], collection["edges"], terms)
        null = cached[(assignment, name, "counts" if statistic != "breadth" else "breadth")].astype(float)
        observed = count.astype(float) if statistic != "breadth" else arm_count.astype(float)
        if statistic == "composition":
            observed = observed / denominator if denominator else np.zeros_like(observed)
            denom = cached[(assignment, name, "denominators")]
            null = np.divide(null, denom[:, None], out=np.zeros_like(null), where=denom[:, None] > 0)
        zobs, znull, degenerate = standardized(observed, null)
        maxima = np.max(znull, axis=1)
        expected = [(len(maxima) + 1) / (len(maxima) + 1) if degenerate[i] or observed[i] == 0 else
                    (np.count_nonzero(maxima >= zobs[i]) + 1) / (len(maxima) + 1)
                    for i in range(len(observed))]
        np.testing.assert_allclose([float(row["collection_maxT_p"]) for row in family], expected, atol=1e-12)
        if assignment == "midpoint" and manifest["mode"] == "primary":
            global_maxima = np.maximum(global_maxima, maxima)
    if manifest["mode"] == "primary":
        for row in result_rows:
            if row["assignment"] != "midpoint":
                if row["global_maxT_p"]:
                    raise AssertionError("sensitivity assignment leaked into global maxT")
                continue
            expected = 1.0 if row["non_informative"] == "1" or float(row["observed"]) == 0 else (
                np.count_nonzero(global_maxima >= float(row["z_observed"])) + 1) / (len(global_maxima) + 1)
            if not math.isclose(float(row["global_maxT_p"]), float(expected), abs_tol=1e-12):
                raise AssertionError("global maxT mismatch")

    return {
        "schema_version": "chm13-independent-run-validation-v1",
        "run": str(run_dir), "permutations": int(manifest["completed_permutations"]),
        "placement_rows_recounted": placement_rows, "physical_loci": len(loci),
        "collections": {name: {"retained_terms": len(value["terms"]), "source": value["path"]}
                        for name, value in collections.items()},
        "result_rows": len(result_rows), "checks": {
            "coordinate_recount": "pass", "cached_arrays": "pass", "empirical_p": "pass",
            "bh_families": "pass", "collection_maxT": "pass", "global_maxT": "pass",
            "physical_source_trace": "pass" if require_frozen_trace else "not_required",
            "impossible_contingencies": "none",
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--allow-minimal-terms", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = validate(args.run, not args.allow_minimal_terms)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
