#!/usr/bin/env python3
"""Selectively extend only V6 hypotheses unresolved after the complete screen.

Collection/global maxT remains frozen to the complete 99,999-hypothesis screen;
this program can refine raw/BH uncertainty for the frozen candidate list but a
selected subset is never used to construct a maxT null.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

import v6_engine as V


WORK = V.WORK / "selective_extension"
RELEASE = V.RELEASE
EXTENSION_CHECKPOINTS = (249_999, 599_999, 999_999)


def candidate_indices(release: Path) -> list[int]:
    path = release / "SELECTIVE_EXTENSION_CANDIDATES.tsv"
    if not path.is_file():
        raise FileNotFoundError("finalize must freeze extension candidates first")
    return [int(row["hypothesis_index"]) for row in V.read_rows(path)]


def reduced_arrays(full: V.WeightedTermArrays, active: Sequence[int]) -> V.WeightedTermArrays:
    position = {value: index for index, value in enumerate(active)}
    rows = []
    for terms in full:
        rows.append(np.asarray([position[int(value)] for value in terms if int(value) in position],
                               dtype=np.int32))
    return V.WeightedTermArrays(rows, full.weights.copy())


def initial_values(indices: Sequence[int]) -> dict[int, np.ndarray]:
    manifest = json.loads((V.WORK / "primary/run_manifest.json").read_text())
    return {
        index: V.COUNT_CACHE.load(V.WORK / "primary", manifest, "midpoint", index, index + 1)[:, 0].copy()
        for index in indices
    }


def update_candidate_statistics(
    rows: list[dict[str, object]], values: Mapping[int, np.ndarray], checkpoint: int,
) -> None:
    for index, null in values.items():
        row = rows[index]
        observed = int(row["observed_physical_copy_burden"])
        degenerate = bool(np.all(null == observed))
        if observed == 0 or degenerate:
            exceedances = len(null)
            p_value = 1.0
        else:
            exceedances = int(np.count_nonzero(null >= observed))
            p_value = (exceedances + 1.0) / (len(null) + 1.0)
        mc95 = V.plus_one_interval(exceedances, len(null), 0.05)
        sequential = V.plus_one_interval(exceedances, len(null), 0.0125)
        q025, median, q975 = np.quantile(null, [0.025, 0.5, 0.975])
        mean = float(np.mean(null))
        row.update({
            "null_mean": mean,
            "null_median": float(median),
            "null_q025": float(q025),
            "null_q975": float(q975),
            "null_max": int(np.max(null)),
            "count_difference": float(observed - median),
            "enrichment_ratio": float((observed + 0.5) / (mean + 0.5)),
            "raw_exceedances": exceedances,
            "raw_permutations": len(null),
            "p_empirical": p_value,
            "p_mc95_lower": mc95[0],
            "p_mc95_upper": mc95[1],
            "p_sequential_lower": sequential[0],
            "p_sequential_upper": sequential[1],
            "non_informative_constant": int(degenerate),
            "inference_stage": f"SELECTIVE_RAW_BH_EXTENSION_{checkpoint}",
        })
    # Mixed-depth p-values are valid term estimates; BH is recomputed over the
    # complete frozen family, while maxT remains explicitly at complete B0.
    for collection in V.COLLECTIONS:
        family_indices = [index for index, row in enumerate(rows) if row["collection"] == collection]
        p = [float(rows[index]["p_empirical"]) for index in family_indices]
        lower = [float(rows[index]["p_sequential_lower"]) for index in family_indices]
        upper = [float(rows[index]["p_sequential_upper"]) for index in family_indices]
        bh = V.E.bh_adjust(p)
        by = V.E.by_adjust(p)
        bh_lower = V.E.bh_adjust(lower)
        bh_upper = V.E.bh_adjust(upper)
        for local, index in enumerate(family_indices):
            rows[index]["bh_q"] = float(bh[local])
            rows[index]["by_q"] = float(by[local])
            rows[index]["bh_sequential_lower"] = float(bh_lower[local])
            rows[index]["bh_sequential_upper"] = float(bh_upper[local])
            rows[index]["mc_status"] = V.classify_primary(rows[index])
    for index in values:
        reasons = V.candidate_reasons(rows[index])
        if reasons == ["MAXT_UNRESOLVED_NO_SELECTIVE_EXTENSION"]:
            raise RuntimeError("selective extension cannot turn a row into a maxT candidate")
        rows[index]["stopping_reason"] = (
            "SELECTIVE_EXTENSION_UNRESOLVED:" + ";".join(reasons)
            if reasons and checkpoint < EXTENSION_CHECKPOINTS[-1]
            else "SELECTIVE_EXTENSION_CAP_UNRESOLVED:" + ";".join(reasons)
            if reasons
            else f"SELECTIVE_EXTENSION_RESOLVED_AT_{checkpoint}"
        )


def run(release: Path = RELEASE, batch_size: int = 1_000) -> dict[str, object]:
    initial = json.loads((V.WORK / "primary/run_manifest.json").read_text())
    if int(initial["completed_permutations"]) != V.INITIAL_PERMUTATIONS:
        raise RuntimeError("targeted extension requires the complete 99,999 initial prefix")
    stage = json.loads((release / "STAGE_DECISION.json").read_text())
    selected = candidate_indices(release)
    if len(selected) != int(stage["selective_extension_candidates"]):
        raise RuntimeError("candidate file and stage decision disagree")
    if WORK.exists() and any(WORK.iterdir()):
        raise FileExistsError("refusing to overwrite a selective-extension work directory")
    WORK.mkdir(parents=True, exist_ok=True)
    primary_rows = [
        row for row in V.read_rows(release / "TERM_RESULTS.tsv.gz")
        if row["assignment"] == "midpoint"
    ]
    sensitivity_rows = [
        row for row in V.read_rows(release / "TERM_RESULTS.tsv.gz")
        if row["assignment"] == "overlap"
    ]
    if not selected:
        report = {
            "schema_version": V.SCHEMA_VERSION,
            "status": "NOT_REQUIRED_INITIAL_SCREEN_SUFFICIENT",
            "initial_permutations": V.INITIAL_PERMUTATIONS,
            "selected_hypotheses": 0,
            "additional_placements": 0,
            "maxT_complete_family_permutations": V.INITIAL_PERMUTATIONS,
            "completed_utc": V.utcnow(),
        }
        V.atomic_json(release / "SELECTIVE_EXTENSION_STATUS.json", report)
        return report

    hypotheses = V.load_hypotheses()
    _arms, _intervals, _loci, genome, blocks = V.N.engine_objects()
    arms = _arms
    full_arrays = V.build_locus_term_arrays(genome, hypotheses)
    sampler = V.E.RegionSampler(blocks, arms, V.E.load_masks(None, arms), "primary", min_candidates=100)
    rng = np.random.Generator(np.random.PCG64DXSM())
    rng.bit_generator.state = initial["rng_state"]
    all_values = {index: [values] for index, values in initial_values(selected).items()}
    active = list(selected)
    cursor = V.INITIAL_PERMUTATIONS
    history = []
    manifest = {
        "schema_version": V.SCHEMA_VERSION,
        "created_utc": V.utcnow(),
        "git_commit": V.git_commit(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "seed": V.MASTER_SEED,
        "spawn_key": [0],
        "continued_initial_rng_state": True,
        "initial_manifest_sha256": V.sha256(V.WORK / "primary/run_manifest.json"),
        "candidate_sha256": V.sha256(release / "SELECTIVE_EXTENSION_CANDIDATES.tsv"),
        "maxT_recomputed_on_subset": False,
        "stages": [],
    }
    V.atomic_json(WORK / "run_manifest.json", manifest)

    for checkpoint in EXTENSION_CHECKPOINTS:
        if not active:
            break
        stage_start = cursor + 1
        stage_dir = WORK / f"checkpoint_{checkpoint}"
        stage_dir.mkdir()
        local_arrays = reduced_arrays(full_arrays, active)
        batches = []
        while cursor < checkpoint:
            stop = min(checkpoint, cursor + batch_size)
            started = time.perf_counter()
            placements_batch = []
            counts = np.zeros((stop - cursor, len(active)), dtype=np.uint16)
            for offset, replicate in enumerate(range(cursor + 1, stop + 1)):
                placements = sampler.sample(rng, replicate)
                placements_batch.append(placements)
                selected_loci = V.N.selected_for_placements(placements, arms, genome, "midpoint")
                counts[offset] = V.recount(selected_loci, local_arrays, len(active))
            stem = f"{cursor + 1:09d}-{stop:09d}"
            placement_path = stage_dir / f"placements.{stem}.tsv.gz"
            count_path = stage_dir / f"midpoint.counts.{stem}.npy"
            V.E.save_placement_batch(placement_path, placements_batch, arms)
            V.E.atomic_save_npy(count_path, counts)
            for local, hypothesis_index in enumerate(active):
                all_values[hypothesis_index].append(counts[:, local].copy())
            batches.append({
                "first": cursor + 1, "last": stop,
                "placements": placement_path.name, "placement_sha256": V.sha256(placement_path),
                "counts": count_path.name, "counts_sha256": V.sha256(count_path),
                "elapsed_seconds": time.perf_counter() - started,
            })
            cursor = stop
            manifest["rng_state"] = V.E.jsonable_rng_state(rng.bit_generator.state)
            manifest["completed_permutations"] = cursor
            manifest["updated_utc"] = V.utcnow()
            V.atomic_json(WORK / "run_manifest.json", manifest)
        update_candidate_statistics(primary_rows, {
            index: np.concatenate(all_values[index]) for index in active
        }, checkpoint)
        # BH is a family operation. Re-evaluate every originally selected row
        # after each family update; a previously stopped row may re-enter and
        # simply accumulates another iid subset of the continuing stream.
        resolved = [index for index in selected if not V.candidate_reasons(primary_rows[index])]
        unresolved = [index for index in selected if V.candidate_reasons(primary_rows[index])]
        stage_record = {
            "checkpoint": checkpoint, "first_replicate": stage_start,
            "last_replicate": checkpoint, "active_hypotheses": active,
            "active_hypothesis_sha256": hashlib.sha256(
                ("\n".join(map(str, active)) + "\n").encode()
            ).hexdigest(),
            "resolved_at_checkpoint": resolved, "unresolved_after_checkpoint": unresolved,
            "batches": batches,
        }
        history.append(stage_record)
        manifest["stages"].append(stage_record)
        manifest["rng_state"] = V.E.jsonable_rng_state(rng.bit_generator.state)
        V.atomic_json(WORK / "run_manifest.json", manifest)
        V.write_rows(release / "TERM_RESULTS.tsv.gz", V.RESULT_FIELDS, primary_rows + sensitivity_rows)
        active = unresolved

    total_additional = cursor - V.INITIAL_PERMUTATIONS
    report = {
        "schema_version": V.SCHEMA_VERSION,
        "status": "COMPLETE_RESOLVED" if not active else "COMPLETE_CAP_REACHED_UNRESOLVED",
        "initial_permutations": V.INITIAL_PERMUTATIONS,
        "selected_hypotheses": len(selected),
        "additional_placements": total_additional,
        "final_stream_replicate": cursor,
        "resolved_hypotheses": len(selected) - len(active),
        "unresolved_at_cap": len(active),
        "unresolved_hypothesis_indices": active,
        "checkpoints_completed": [row["checkpoint"] for row in history],
        "maxT_complete_family_permutations": V.INITIAL_PERMUTATIONS,
        "maxT_recomputed_on_selected_subset": False,
        "rng_state_saved": True,
        "work_manifest": str((WORK / "run_manifest.json").relative_to(V.REPO)),
        "completed_utc": V.utcnow(),
    }
    V.atomic_json(release / "SELECTIVE_EXTENSION_STATUS.json", report)
    stage.update({
        "status": report["status"],
        "extension_checkpoints_completed": report["checkpoints_completed"],
        "selective_extension_resolved": report["resolved_hypotheses"],
        "selective_extension_unresolved_at_cap": report["unresolved_at_cap"],
        "additional_joint_placements": total_additional,
        "final_decision_utc": V.utcnow(),
    })
    V.atomic_json(release / "STAGE_DECISION.json", stage)
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, default=RELEASE)
    parser.add_argument("--batch-size", type=int, default=1_000)
    args = parser.parse_args(argv)
    print(json.dumps(run(args.release, args.batch_size), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
