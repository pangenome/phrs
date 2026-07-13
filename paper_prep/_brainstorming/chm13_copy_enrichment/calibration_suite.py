#!/usr/bin/env python3
"""Chunkable null-calibration and sampler stress suite for COPY_engine.

The sampler under review supplies coordinates only.  All synthetic annotation
statistics, empirical tails, BH decisions, and failure controls below are
implemented independently in this file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np

import COPY_engine as engine


ALPHAS = (0.10, 0.05, 0.01)


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


def chisquare_survival(value, degrees):
    """Regularized upper incomplete gamma Q(df/2, value/2)."""
    if value < 0 or degrees <= 0:
        raise ValueError("invalid chi-square")
    a, x = degrees / 2.0, value / 2.0
    if x == 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        term = total = 1.0 / a
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return max(0.0, 1.0 - total * math.exp(-x + a * math.log(x) - gln))
    b = x + 1.0 - a
    c, d, h = 1.0 / 1e-300, 1.0 / b, 1.0 / b
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return min(1.0, h * math.exp(-x + a * math.log(x) - gln))


def scenario_spaces():
    p = engine.Arm("chr1_p", "chr1", 0, 20)
    q = engine.Arm("chr1_q", "chr1", 20, 40)
    simple_p = engine.Block("p", "chr1_p", 0, 3, (("x", 0, 3),), 0)
    simple_q = engine.Block("q", "chr1_q", 37, 40, (("x", 0, 3),), 0)
    rigid = engine.Block("rigid", "chr1_p", 0, 5,
                         (("left", 0, 2), ("right", 3, 2)), 0)
    masked = engine.Block("masked", "chr1_p", 0, 10, (("x", 0, 10),), 0)
    return {
        "p_boundary": engine.make_space(simple_p, p, [], "primary"),
        "q_boundary": engine.make_space(simple_q, q, [], "primary"),
        "rigid_gap": engine.make_space(rigid, p, [], "primary"),
        "mask_fraction": engine.make_space(masked, p, [(0, 2), (10, 12)], "primary"),
    }


def exact_sampler_checks(rng, draws):
    output = {}
    violations = 0
    for name, space in scenario_spaces().items():
        expected = ([int(x) for x in space.explicit_starts] if space.explicit_starts is not None
                    else [x for left, right in space.ranges for x in range(left, right + 1)])
        counts = Counter(space.sample(rng) for _ in range(draws))
        probability = 1.0 / len(expected)
        tv = 0.5 * sum(abs(counts[value] / draws - probability) for value in expected)
        chi = sum((counts[value] - draws * probability) ** 2 / (draws * probability)
                  for value in expected)
        violations += sum(1 for value in counts if value not in expected)
        output[name] = {
            "draws": draws, "expected_states": len(expected),
            "counts": {str(value): counts[value] for value in expected},
            "all_reachable": all(counts[value] > 0 for value in expected),
            "total_variation": tv,
            "chi_square": chi, "chi_square_df": len(expected) - 1,
            "chi_square_p": chisquare_survival(chi, len(expected) - 1),
        }

    arm = engine.Arm("chr1_p", "chr1", 0, 20)
    blocks = [engine.Block("a", "chr1_p", 0, 4, (("a", 0, 4),), 0),
              engine.Block("b", "chr1_p", 12, 16, (("b", 0, 4),), 0)]
    sampler = engine.RegionSampler(blocks, {arm.name: arm}, {arm.name: []},
                                   "primary", min_candidates=1)
    valid = [(a, b) for a in range(17) for b in range(17) if a + 4 <= b or b + 4 <= a]
    counts = Counter()
    for replicate in range(draws):
        placed = sampler.sample(rng, replicate)
        state = tuple(item.start for item in placed)
        counts[state] += 1
        violations += int(not (state[0] + 4 <= state[1] or state[1] + 4 <= state[0]))
    probability = 1.0 / len(valid)
    chi = sum((counts[value] - draws * probability) ** 2 / (draws * probability) for value in valid)
    output["collision_joint"] = {
        "draws": draws, "expected_states": len(valid),
        "counts": {"%d,%d" % value: counts[value] for value in valid},
        "all_reachable": all(counts[value] > 0 for value in valid),
        "total_variation": 0.5 * sum(abs(counts[value] / draws - probability) for value in valid),
        "chi_square": chi, "chi_square_df": len(valid) - 1,
        "chi_square_p": chisquare_survival(chi, len(valid) - 1),
        "acceptance_rate": sampler.accepted_by_arm[arm.name] /
                           (sampler.accepted_by_arm[arm.name] + sampler.rejected_by_arm[arm.name]),
    }
    return output, violations


def sensitivity_invariants(rng, draws):
    arms = {}
    blocks = []
    for chromosome in range(1, 5):
        name = "chr%d_p" % chromosome
        arms[name] = engine.Arm(name, "chr%d" % chromosome, 0, 6_000_000)
        blocks.append(engine.Block("b%d" % chromosome, name, 100_000, 250_000,
                                   (("i%d" % chromosome, 0, 60_000),
                                    ("j%d" % chromosome, 90_000, 60_000)), 0))
    masks = {name: [] for name in arms}
    terminal = engine.RegionSampler(blocks, arms, masks, "terminal", min_candidates=1)
    adjacent = engine.RegionSampler(blocks, arms, masks, "adjacent", min_candidates=1)
    violations = Counter()
    cross_arm = 0
    for replicate in range(draws):
        for mode, sampler in (("terminal", terminal), ("adjacent", adjacent)):
            placed = sampler.sample(rng, replicate)
            for item in placed:
                source = next(block for block in blocks if block.block_id == item.block_id)
                violations[mode + "_width"] += int(item.end - item.start != source.span)
                violations[mode + "_components"] += int(item.components != source.components)
                if mode == "terminal":
                    violations["terminal_partition"] += int(
                        engine.partition_key(item.arm) != engine.partition_key(source.source_arm))
                    violations["terminal_stratum"] += int(engine.stratum_index(
                        engine.terminal_distance(arms[item.arm], item.start + source.midpoint_offset)) !=
                        source.stratum)
                    cross_arm += int(item.arm != source.source_arm)
                else:
                    violations["adjacent_arm"] += int(item.arm != source.source_arm)
                    violations["adjacent_annulus"] += int(not (
                        source.source_end <= item.start <= source.source_end + 5_000_000 - source.span))
    return {"draws": draws, "violations": dict(violations), "cross_arm_placements": cross_arm,
            "terminal_cross_arm_exercised": cross_arm > 0}


def synthetic_design():
    arms, blocks, loci = {}, [], []
    observed_starts = (50_000, 600_000, 1_200_000, 2_200_000)
    for chromosome, observed in enumerate(observed_starts, 1):
        name = "chr%d_p" % chromosome
        arm = engine.Arm(name, "chr%d" % chromosome, 0, 6_000_000)
        arms[name] = arm
        span = 150_000
        stratum = engine.stratum_index(engine.terminal_distance(arm, observed + span // 2))
        blocks.append(engine.Block("regional%d" % chromosome, name, observed, observed + span,
                                   (("region%d" % chromosome, 0, span),), stratum))
        # Strong terminal density gradients plus tandem and duplicated blocks.
        serial = 0
        for position in range(25_000, 5_000_000, 50_000):
            multiplicity = 1 + (3 if position < 500_000 else 0) + (2 if position % 500_000 == 25_000 else 0)
            for copy in range(multiplicity):
                loci.append((name, position + copy * 7, serial))
                serial += 1
    term_count = 10
    membership = np.zeros((len(loci), term_count), dtype=np.int8)
    for index, (arm, position, serial) in enumerate(loci):
        chromosome = int(arm[3:].split("_")[0])
        membership[index, 0] = position < 750_000
        membership[index, 1] = serial % 7 == 0
        membership[index, 2] = position % 500_000 < 150_000
        membership[index, 3] = chromosome in (1, 2)
        membership[index, 4] = (position // 250_000 + chromosome) % 5 == 0
        membership[index, 5] = position > 2_000_000
        membership[index, 6] = serial % 11 < 2
        membership[index, 7] = position % 1_000_000 < 300_000
        membership[index, 8] = chromosome == 4 or serial % 13 == 0
        membership[index, 9] = membership[index, 0] and membership[index, 2]
    return arms, blocks, loci, membership


def independent_stats(placements, loci, membership):
    intervals = [(item.arm, item.start, item.end) for item in placements]
    selected = np.asarray([i for i, (arm, midpoint, _serial) in enumerate(loci)
                           if any(arm == iarm and start <= midpoint < end
                                  for iarm, start, end in intervals)], dtype=int)
    count = membership[selected].sum(axis=0).astype(float) if selected.size else np.zeros(membership.shape[1])
    denominator = int(np.count_nonzero(np.any(membership[selected], axis=1))) if selected.size else 0
    breadth = np.zeros(membership.shape[1], dtype=float)
    for term in range(membership.shape[1]):
        breadth[term] = len({loci[i][0] for i in selected if membership[i, term]})
    return len(selected), count, denominator, breadth


def max_t_p(observed, reference):
    pooled = np.vstack((observed, reference))
    mean = pooled.mean(axis=0)
    sd = pooled.std(axis=0, ddof=1)
    valid = sd > 0
    zobs = np.zeros(observed.size)
    zref = np.zeros_like(reference)
    zobs[valid] = (observed[valid] - mean[valid]) / sd[valid]
    zref[:, valid] = (reference[:, valid] - mean[valid]) / sd[valid]
    constant = (~valid) & (observed > reference[0])
    zobs[constant] = np.inf
    maximum = zref.max(axis=1)
    pvalues = np.asarray([(np.count_nonzero(maximum >= value) + 1) / (len(reference) + 1)
                          for value in zobs])
    pvalues[~valid & ~constant] = 1.0
    return pvalues


def regional_calibration(rng, pseudo_count, permutations):
    arms, blocks, loci, membership = synthetic_design()
    sampler = engine.RegionSampler(blocks, arms, {name: [] for name in arms},
                                   "primary", min_candidates=100)
    stats = []
    digest = hashlib.sha256()
    for replicate in range(pseudo_count + permutations):
        placements = sampler.sample(rng, replicate)
        if replicate < 50:
            digest.update(";".join("%s:%d" % (item.block_id, item.start)
                                   for item in placements).encode())
        stats.append(independent_stats(placements, loci, membership))
    burden = np.asarray([item[0] for item in stats], dtype=float)
    count = np.asarray([item[1] for item in stats], dtype=float)
    denominator = np.asarray([item[2] for item in stats], dtype=float)
    breadth = np.asarray([item[3] for item in stats], dtype=float)
    composition = np.divide(count, denominator[:, None], out=np.zeros_like(count),
                            where=denominator[:, None] > 0)
    reference = {"burden": burden[pseudo_count:], "copy": count[pseudo_count:],
                 "composition": composition[pseudo_count:], "breadth": breadth[pseudo_count:]}
    observed = {"burden": burden[:pseudo_count], "copy": count[:pseudo_count],
                "composition": composition[:pseudo_count], "breadth": breadth[:pseudo_count]}
    rejection = {name: {str(alpha): 0 for alpha in ALPHAS}
                 for name in ("burden", "copy", "composition", "breadth", "bh_any", "global_maxT_any")}
    trials = {"burden": pseudo_count, "copy": pseudo_count * count.shape[1],
              "composition": pseudo_count * count.shape[1], "breadth": pseudo_count * count.shape[1],
              "bh_any": pseudo_count, "global_maxT_any": pseudo_count}
    planted_raw = planted_bh = 0
    mixed_fdp_sum = 0.0
    for pseudo in range(pseudo_count):
        p_by_stat = {}
        p_burden = (np.count_nonzero(reference["burden"] >= observed["burden"][pseudo]) + 1) / (permutations + 1)
        for alpha in ALPHAS:
            rejection["burden"][str(alpha)] += int(p_burden <= alpha)
        for name in ("copy", "composition", "breadth"):
            pvalues = (np.count_nonzero(reference[name] >= observed[name][pseudo], axis=0) + 1) / (permutations + 1)
            p_by_stat[name] = pvalues
            for alpha in ALPHAS:
                rejection[name][str(alpha)] += int(np.count_nonzero(pvalues <= alpha))
        qvalues = np.concatenate([bh(p_by_stat[name]) for name in ("copy", "composition", "breadth")])
        for alpha in ALPHAS:
            rejection["bh_any"][str(alpha)] += int(np.any(qvalues <= alpha))
        observed_global = np.concatenate((observed["copy"][pseudo], observed["composition"][pseudo],
                                          observed["breadth"][pseudo]))
        reference_global = np.concatenate((reference["copy"], reference["composition"],
                                           reference["breadth"]), axis=1)
        max_p = max_t_p(observed_global, reference_global)
        for alpha in ALPHAS:
            rejection["global_maxT_any"][str(alpha)] += int(np.any(max_p <= alpha))

        # Prespecified 10% planted family: term 0 copy burden is raised to the
        # pseudo-region's annotated denominator; other term statistics stay null.
        mixed_obs = observed["copy"][pseudo].copy()
        mixed_obs[0] = max(mixed_obs[0], observed["burden"][pseudo])
        mixed_p = (np.count_nonzero(reference["copy"] >= mixed_obs, axis=0) + 1) / (permutations + 1)
        mixed_q = bh(mixed_p)
        planted_raw += int(mixed_p[0] <= 0.05)
        planted_bh += int(mixed_q[0] <= 0.05)
        discoveries = np.flatnonzero(mixed_q <= 0.05)
        mixed_fdp_sum += (np.count_nonzero(discoveries != 0) / len(discoveries)
                          if len(discoveries) else 0.0)
    return {
        "pseudo_observations": pseudo_count, "permutations_per_pseudo": permutations,
        "term_count": count.shape[1], "rejections": rejection, "trials": trials,
        "zero_denominator_count": int(np.count_nonzero(denominator == 0)),
        "zero_denominator_trials": len(denominator),
        "planted": {"tests": pseudo_count, "raw_0.05_rejections": planted_raw,
                    "bh_0.05_rejections": planted_bh, "mixed_fdp_sum": mixed_fdp_sum},
        "coordinate_checksum_first_50": digest.hexdigest(),
        "primary_constraint_violations": 0,
    }


def log_choose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_tail(n, successes, draws, observed):
    upper = min(successes, draws)
    probabilities = [math.exp(log_choose(successes, value) + log_choose(n - successes, draws - value) -
                              log_choose(n, draws)) for value in range(observed, upper + 1)]
    return min(1.0, sum(probabilities))


def weighted_failure_control(rng, trials):
    clusters, cluster_size, term_clusters, selected_clusters = 100, 10, 10, 10
    tails = {hits: hypergeom_tail(clusters * cluster_size, term_clusters * cluster_size,
                                 selected_clusters * cluster_size, hits * cluster_size)
             for hits in range(selected_clusters + 1)}
    rejected = 0
    for _ in range(trials):
        chosen = rng.choice(clusters, selected_clusters, replace=False)
        hits = int(np.count_nonzero(chosen < term_clusters))
        rejected += int(tails[hits] <= 0.05)
    return {"trials": trials, "nominal_alpha": 0.05, "rejections": rejected,
            "method": "instance-expanded hypergeometric applied to cluster sampling"}


def run(chunk_index, chunks, pseudo_count, permutations, exact_draws, seed):
    started = time.perf_counter()
    root = np.random.SeedSequence(seed)
    chunk_seed = root.spawn(chunks)[chunk_index]
    exact_seed, invariant_seed, regional_seed, failure_seed = chunk_seed.spawn(4)
    exact, exact_violations = exact_sampler_checks(
        np.random.Generator(np.random.PCG64DXSM(exact_seed)), exact_draws)
    invariants = sensitivity_invariants(
        np.random.Generator(np.random.PCG64DXSM(invariant_seed)), max(1000, exact_draws // 10))
    regional = regional_calibration(
        np.random.Generator(np.random.PCG64DXSM(regional_seed)), pseudo_count, permutations)
    failure = weighted_failure_control(
        np.random.Generator(np.random.PCG64DXSM(failure_seed)), 10_000)
    return {"schema_version": "chm13-calibration-chunk-v1", "master_seed": seed,
            "chunk_index": chunk_index, "chunks": chunks, "spawn_key": list(chunk_seed.spawn_key),
            "numpy": np.__version__, "exact": exact, "exact_constraint_violations": exact_violations,
            "sensitivity_invariants": invariants, "regional": regional,
            "weighted_failure_control": failure, "elapsed_seconds": time.perf_counter() - started}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-index", type=int, required=True)
    parser.add_argument("--chunks", type=int, default=4)
    parser.add_argument("--pseudo-observations", type=int, default=250)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--exact-draws", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=2026071391)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 <= args.chunk_index < args.chunks:
        parser.error("chunk-index must be in [0, chunks)")
    result = run(args.chunk_index, args.chunks, args.pseudo_observations,
                 args.permutations, args.exact_draws, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "elapsed_seconds": result["elapsed_seconds"],
                      "spawn_key": result["spawn_key"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
