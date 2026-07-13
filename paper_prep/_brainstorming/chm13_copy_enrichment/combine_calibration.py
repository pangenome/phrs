#!/usr/bin/env python3
"""Combine deterministic Slurm calibration chunks and evaluate hard gates."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from calibration_suite import ALPHAS, chisquare_survival


def betacf(a, b, x):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1e-300 if abs(d) < 1e-300 else d
    d = 1.0 / d
    h = d
    for iteration in range(1, 301):
        twice = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + twice) * (a + twice))
        d = 1.0 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1.0 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / ((a + twice) * (qap + twice))
        d = 1.0 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1.0 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-14:
            return h
    raise ArithmeticError("beta continued fraction did not converge")


def regularized_beta(x, a, b):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    factor = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
                      a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return factor * betacf(a, b, x) / a
    return 1.0 - factor * betacf(b, a, 1 - x) / b


def beta_quantile(probability, a, b):
    lower, upper = 0.0, 1.0
    for _ in range(100):
        middle = (lower + upper) / 2
        if regularized_beta(middle, a, b) < probability:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def interval(successes, trials):
    lower = 0.0 if successes == 0 else beta_quantile(0.025, successes, trials - successes + 1)
    upper = 1.0 if successes == trials else beta_quantile(0.975, successes + 1, trials - successes)
    return lower, upper


def combine(paths):
    chunks = [json.loads(path.read_text()) for path in paths]
    indices = [item["chunk_index"] for item in chunks]
    if len(set(indices)) != len(indices) or sorted(indices) != list(range(chunks[0]["chunks"])):
        raise ValueError("chunks must contain every unique index exactly once")
    if any(item["master_seed"] != chunks[0]["master_seed"] or item["chunks"] != len(chunks)
           for item in chunks):
        raise ValueError("incompatible calibration chunk manifests")

    exact = {}
    exact_pass = True
    for scenario in chunks[0]["exact"]:
        counts = Counter()
        for item in chunks:
            counts.update(item["exact"][scenario]["counts"])
        draws = sum(item["exact"][scenario]["draws"] for item in chunks)
        states = chunks[0]["exact"][scenario]["expected_states"]
        expected = draws / states
        chi = sum((value - expected) ** 2 / expected for value in counts.values())
        tv = 0.5 * sum(abs(value / draws - 1 / states) for value in counts.values())
        record = {"draws": draws, "states": states, "all_reachable": len(counts) == states and min(counts.values()) > 0,
                  "total_variation": tv, "chi_square": chi,
                  "chi_square_p": chisquare_survival(chi, states - 1)}
        record["pass"] = (record["all_reachable"] and tv <= 0.01 and record["chi_square_p"] >= 0.001)
        exact_pass &= record["pass"]
        exact[scenario] = record

    regional_counts = {}
    regional_pass = True
    for family in chunks[0]["regional"]["rejections"]:
        trials = sum(item["regional"]["trials"][family] for item in chunks)
        regional_counts[family] = {}
        for alpha in ALPHAS:
            key = str(alpha)
            rejected = sum(item["regional"]["rejections"][family][key] for item in chunks)
            lower, upper = interval(rejected, trials)
            rate = rejected / trials
            ceiling = max(alpha + 0.02, 1.5 * alpha)
            passed = lower <= alpha and rate <= ceiling
            regional_pass &= passed
            regional_counts[family][key] = {
                "rejections": rejected, "trials": trials, "rate": rate,
                "ci_lower": lower, "ci_upper": upper, "point_tolerance_ceiling": ceiling,
                "pass": passed,
            }

    bh_gate = regional_counts["bh_any"]["0.05"]
    bh_gate["fdr_specific_pass"] = bh_gate["rate"] <= 0.07 and bh_gate["ci_lower"] <= 0.05
    max_gate = regional_counts["global_maxT_any"]["0.05"]
    max_gate["fwer_specific_pass"] = max_gate["rate"] <= 0.07 and max_gate["ci_lower"] <= 0.05
    regional_pass &= bh_gate["fdr_specific_pass"] and max_gate["fwer_specific_pass"]

    zero = sum(item["regional"]["zero_denominator_count"] for item in chunks)
    zero_trials = sum(item["regional"]["zero_denominator_trials"] for item in chunks)
    zero_rate = zero / zero_trials
    regional_pass &= zero_rate <= 0.01
    planted_tests = sum(item["regional"]["planted"]["tests"] for item in chunks)
    planted_raw = sum(item["regional"]["planted"]["raw_0.05_rejections"] for item in chunks)
    planted_bh = sum(item["regional"]["planted"]["bh_0.05_rejections"] for item in chunks)
    mixed_fdr = sum(item["regional"]["planted"]["mixed_fdp_sum"] for item in chunks) / planted_tests
    planted_pass = planted_raw > 0 and planted_bh > 0 and mixed_fdr <= 0.07

    failure_rejected = sum(item["weighted_failure_control"]["rejections"] for item in chunks)
    failure_trials = sum(item["weighted_failure_control"]["trials"] for item in chunks)
    failure_lower, failure_upper = interval(failure_rejected, failure_trials)
    failure_pass = failure_lower > 0.05

    constraint_violations = sum(item["exact_constraint_violations"] +
                                sum(item["sensitivity_invariants"]["violations"].values()) +
                                item["regional"]["primary_constraint_violations"] for item in chunks)
    sensitivity_exercised = all(item["sensitivity_invariants"]["terminal_cross_arm_exercised"] for item in chunks)
    constraints_pass = constraint_violations == 0 and sensitivity_exercised
    report = {
        "schema_version": "chm13-calibration-combined-v1", "master_seed": chunks[0]["master_seed"],
        "chunks": len(chunks), "spawn_keys": [item["spawn_key"] for item in chunks],
        "chunk_files": [str(path) for path in paths], "elapsed_seconds_sum": sum(item["elapsed_seconds"] for item in chunks),
        "tolerance": {
            "exact_sampler": "chi-square p >= 0.001, TV <= 0.01, every state reachable",
            "type_i": "95% lower limit <= alpha and point <= max(alpha+0.02, 1.5*alpha)",
            "BH_FDR_0.05": "point <= 0.07 and 95% lower limit <= 0.05",
            "global_maxT_FWER_0.05": "point <= 0.07 and 95% lower limit <= 0.05",
            "composition_zero_denominator": "rate <= 0.01",
        },
        "exact_sampler": exact,
        "regional_null": regional_counts,
        "zero_denominator": {"count": zero, "trials": zero_trials, "rate": zero_rate,
                             "pass": zero_rate <= 0.01},
        "planted_effect": {"tests": planted_tests, "raw_power_0.05": planted_raw / planted_tests,
                           "bh_power_0.05": planted_bh / planted_tests,
                           "null_term_mean_fdp": mixed_fdr, "pass": planted_pass},
        "weighted_hypergeometric_failure_control": {
            "rejections": failure_rejected, "trials": failure_trials,
            "rate": failure_rejected / failure_trials, "ci_lower": failure_lower,
            "ci_upper": failure_upper, "positive_control_pass": failure_pass},
        "constraint_checks": {"violations": constraint_violations,
                              "terminal_cross_arm_exercised_in_every_chunk": sensitivity_exercised,
                              "pass": constraints_pass},
    }
    report["overall_pass"] = exact_pass and regional_pass and planted_pass and failure_pass and constraints_pass
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chunks", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = combine(args.chunks)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "overall_pass": report["overall_pass"]}))
    return 0 if report["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
