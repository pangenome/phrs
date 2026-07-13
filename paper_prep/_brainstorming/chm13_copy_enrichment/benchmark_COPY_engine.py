#!/usr/bin/env python3
"""Reproducible GO-scale throughput and memory benchmark for COPY_engine.py."""

from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np

import COPY_engine as engine


def synthetic_collection(genome: engine.GenomeIndex, terms: int,
                         annotations_per_locus: int) -> engine.Collection:
    """Deterministic, correlated GO-like many-to-many annotation map."""
    locus_terms = []
    counts = np.zeros(terms, dtype=np.int64)
    arm_sets = [set() for _ in range(terms)]
    for locus in genome.loci:
        # Consecutive copies have overlapping memberships, mimicking propagated
        # ancestors and family correlation without constructing a huge TSV.
        base = (locus.index * 104729 + engine.canonical_arm_key(locus.arm)[0] * 1009) % terms
        values = np.unique(np.asarray(
            [(base + offset * 8191) % terms for offset in range(annotations_per_locus)],
            dtype=np.int32))
        locus_terms.append(values)
        counts[values] += 1
        for value in values:
            arm_sets[int(value)].add(locus.arm)
    return engine.Collection(
        "SYNTHETIC_GO", Path("synthetic://go"),
        ["GO:%07d" % index for index in range(terms)],
        ["synthetic term %d" % index for index in range(terms)],
        locus_terms, counts, np.asarray([len(values) for values in arm_sets], dtype=np.int64),
        np.asarray([values.size > 0 for values in locus_terms], dtype=np.bool_), [])


def main() -> int:
    directory = Path(__file__).with_name("analysis_ready")
    parser = argparse.ArgumentParser()
    parser.add_argument("--arms", type=Path, default=directory / "chm13_arm_summary.tsv")
    parser.add_argument("--intervals", type=Path, default=directory / "chm13_phr_intervals.tsv")
    parser.add_argument("--loci", type=Path, default=directory / "chm13_gene_loci.tsv.gz")
    parser.add_argument("--terms", type=int, default=15_000)
    parser.add_argument("--annotations-per-locus", type=int, default=12)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--final-permutations", type=int, default=99_999)
    parser.add_argument("--seed", type=int, default=2026071301)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rss_start = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    setup_start = time.perf_counter()
    arms = engine.load_arms(args.arms)
    intervals = engine.load_intervals(args.intervals, arms)
    loci = engine.load_loci(args.loci, arms)
    genome = engine.GenomeIndex(loci, arms)
    blocks = engine.build_blocks(intervals, arms)
    collection = synthetic_collection(genome, args.terms, args.annotations_per_locus)
    sampler = engine.RegionSampler(blocks, arms, {name: [] for name in arms},
                                   "primary", min_candidates=100)
    setup_seconds = time.perf_counter() - setup_start

    rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(args.seed)))
    count_dtype = np.uint16 if len(loci) <= np.iinfo(np.uint16).max else np.uint32
    breadth_dtype = np.uint8 if len(arms) <= np.iinfo(np.uint8).max else np.uint16
    count_store = {assignment: np.empty((args.permutations, args.terms), dtype=count_dtype)
                   for assignment in ("midpoint", "overlap")}
    breadth_store = {assignment: np.empty((args.permutations, args.terms), dtype=breadth_dtype)
                     for assignment in ("midpoint", "overlap")}
    denominator_store = {assignment: np.empty(args.permutations, dtype=np.int32)
                         for assignment in ("midpoint", "overlap")}
    generation_start = time.perf_counter()
    for replicate in range(args.permutations):
        placements = sampler.sample(rng, replicate + 1)
        placed = engine.placed_intervals(placements, arms)
        for assignment in ("midpoint", "overlap"):
            stats = engine.calculate_stats(genome, [collection], placed, assignment)
            count_store[assignment][replicate] = stats.counts[collection.name]
            breadth_store[assignment][replicate] = stats.breadth[collection.name]
            denominator_store[assignment][replicate] = stats.denominators[collection.name]
    generation_seconds = time.perf_counter() - generation_start

    # Exercise the exact pooled maxT kernel on the complete synthetic family.
    inference_start = time.perf_counter()
    observed = engine.calculate_stats(genome, [collection], intervals, "midpoint")
    family_max = np.full(args.permutations, -np.inf)
    chunk_size = 64
    for start in range(0, args.terms, chunk_size):
        end = min(args.terms, start + chunk_size)
        obs = observed.counts[collection.name][start:end]
        _z_obs, z_null, _degenerate = engine.standardized_values(
            obs, count_store["midpoint"][:, start:end])
        family_max = np.maximum(family_max, np.max(z_null, axis=1))
    inference_seconds = time.perf_counter() - inference_start
    rss_end = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    bytes_per_permutation = sum(array[0].nbytes for array in count_store.values())
    bytes_per_permutation += sum(array[0].nbytes for array in breadth_store.values())
    bytes_per_permutation += sum(array[:1].nbytes for array in denominator_store.values())
    seconds_per_permutation = generation_seconds / args.permutations
    report = {
        "schema_version": engine.SCHEMA_VERSION,
        "benchmark": "GO-scale synthetic term map on frozen CHM13 physical loci",
        "physical_loci": len(loci),
        "target_intervals": len(intervals),
        "placement_blocks": len(blocks),
        "terms": args.terms,
        "annotations_per_locus": args.annotations_per_locus,
        "annotation_edges": sum(values.size for values in collection.locus_terms),
        "benchmark_permutations": args.permutations,
        "assignments": 2,
        "setup_seconds": setup_seconds,
        "generation_seconds": generation_seconds,
        "inference_maxT_kernel_seconds": inference_seconds,
        "seconds_per_joint_permutation_both_assignments": seconds_per_permutation,
        "permutations_per_second": 1.0 / seconds_per_permutation,
        "peak_rss_mib": rss_end / 1024.0,
        "peak_rss_increase_mib": max(0, rss_end - rss_start) / 1024.0,
        "uncompressed_stat_bytes_per_permutation": bytes_per_permutation,
        "final_permutations": args.final_permutations,
        "projected_generation_hours": seconds_per_permutation * args.final_permutations / 3600.0,
        "projected_uncompressed_stat_gib": bytes_per_permutation * args.final_permutations / (1024.0 ** 3),
        "maxT_checksum": float(np.sum(family_max)),
        "numpy_version": np.__version__,
    }
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
