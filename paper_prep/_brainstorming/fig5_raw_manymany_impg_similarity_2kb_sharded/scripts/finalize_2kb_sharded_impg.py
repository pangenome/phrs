#!/usr/bin/env python3
"""Assemble completed IMPG shards and write Fig5 2 kb similarity summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path


OUT = Path(__file__).resolve().parents[1]
CHR_RE = re.compile(r"_chr([0-9XY]+)(?::|$)")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_group(group: str) -> dict[str, str | int]:
    seq, coords = group.rsplit(":", 1)
    start_text, end_text = coords.split("-", 1)
    chrom_match = CHR_RE.search(seq)
    chrom = f"chr{chrom_match.group(1)}" if chrom_match else "unknown"
    return {"seq": seq, "chrom": chrom, "start": int(start_text), "end": int(end_text)}


def read_lengths(target_fasta: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with Path(f"{target_fasta}.fai").open() as handle:
        for line in handle:
            seq, length, *_ = line.rstrip("\n").split("\t")
            lengths[seq] = int(length)
    return lengths


def arm(chrom: str, start: int, end: int, seq_len: int | None) -> str:
    if chrom == "unknown":
        return "unknown"
    if seq_len is None:
        return chrom
    mid = (start + end) // 2
    return f"{chrom}{'p' if mid < seq_len / 2 else 'q'}"


def row_target_side(row: dict[str, str]) -> tuple[dict[str, str | int], dict[str, str | int]] | None:
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    target = {
        "seq": chrom,
        "chrom": parse_group(f"{chrom}:{start}-{end}")["chrom"],
        "start": start,
        "end": end,
    }
    a = parse_group(row["group.a"])
    b = parse_group(row["group.b"])
    if a["seq"] == target["seq"] and a["start"] == target["start"] and a["end"] == target["end"]:
        return target, b
    if b["seq"] == target["seq"] and b["start"] == target["start"] and b["end"] == target["end"]:
        return target, a
    return None


def annotate_signal(row: dict[str, object]) -> str:
    target_arm = str(row["target_arm"])
    other_arm = str(row["other_arm"])
    pair = {target_arm, other_arm}
    acro_p = {"chr13p", "chr14p", "chr15p", "chr21p", "chr22p"}
    if pair == {"chrXq", "chrYq"} or pair == {"chrXp", "chrYp"}:
        return "PAR_control"
    if target_arm in acro_p or other_arm in acro_p:
        return "acrocentric_control"
    if pair == {"chr9q", "chr3q"}:
        return "chr9q_chr3q"
    return "other"


def metadata_for(path: str) -> dict[str, object]:
    p = Path(path)
    if not p.exists():
        return {}
    with p.open() as handle:
        return json.load(handle)


def check_shards(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        path = Path(row["output_tsv_gz"])
        meta = metadata_for(row["metadata_json"])
        status = "OK" if path.exists() and path.stat().st_size > 0 and meta.get("status") == "OK" else "MISSING_OR_INCOMPLETE"
        out.append(
            {
                **row,
                "state": status,
                "output_size_bytes": path.stat().st_size if path.exists() else "",
                "slurm_job_id": meta.get("slurm_job_id", row.get("slurm_array_job_id", "")),
                "slurm_array_task_id": meta.get("slurm_array_task_id", row.get("shard_index", "")),
                "node": meta.get("node", ""),
                "partition": meta.get("partition", row.get("partition", "")),
                "slurm_cpus_per_task": meta.get("slurm_cpus_per_task", row.get("slurm_cpus_per_task", "")),
                "start_utc": meta.get("start_utc", ""),
                "finish_utc": meta.get("finish_utc", ""),
            }
        )
    return out


def assemble_one(method: str, comparison_id: str, rows: list[dict[str, str]], output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    data_rows = 0
    with gzip.open(output, "wt") as dst:
        for row in sorted(rows, key=lambda item: int(item["shard_index"])):
            with gzip.open(row["output_tsv_gz"], "rt") as src:
                header = src.readline()
                if not header:
                    raise RuntimeError(f"Empty shard output: {row['output_tsv_gz']}")
                if not wrote_header:
                    dst.write(header)
                    wrote_header = True
                for line in src:
                    dst.write(line)
                    data_rows += 1
    return data_rows


def summarize_assembled(assembled_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    best_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    for manifest in assembled_rows:
        path = Path(str(manifest["assembled_output_tsv_gz"]))
        lengths = read_lengths(str(manifest["target_fasta"]))
        best_by_window: dict[tuple[str, int, int, str], dict[str, object]] = {}
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sides = row_target_side(row)
                if sides is None:
                    continue
                target, other = sides
                if other["seq"] == target["seq"] and other["start"] == target["start"] and other["end"] == target["end"]:
                    continue
                target_arm = arm(str(target["chrom"]), int(target["start"]), int(target["end"]), lengths.get(str(target["seq"])))
                other_arm = arm(str(other["chrom"]), int(other["start"]), int(other["end"]), lengths.get(str(other["seq"])))
                relation = "same_chromosome" if target["chrom"] == other["chrom"] else "interchromosomal"
                out = {
                    "method": manifest["method"],
                    "comparison_id": manifest["comparison_id"],
                    "target_seq": target["seq"],
                    "target_chrom": target["chrom"],
                    "target_start": target["start"],
                    "target_end": target["end"],
                    "target_arm": target_arm,
                    "other_seq": other["seq"],
                    "other_chrom": other["chrom"],
                    "other_start": other["start"],
                    "other_end": other["end"],
                    "other_arm": other_arm,
                    "relation": relation,
                    "intersection": int(float(row["intersection"])),
                    "jaccard_similarity": float(row["jaccard.similarity"]),
                    "cosine_similarity": float(row["cosine.similarity"]),
                    "dice_similarity": float(row["dice.similarity"]),
                    "estimated_identity": float(row["estimated.identity"]),
                }
                out["signal_class"] = annotate_signal(out)
                edge_rows.append(out)
                key = (str(target["seq"]), int(target["start"]), int(target["end"]), relation)
                previous = best_by_window.get(key)
                if previous is None or float(out["estimated_identity"]) > float(previous["estimated_identity"]):
                    best_by_window[key] = out
        best_rows.extend(best_by_window.values())

    fields = [
        "method",
        "comparison_id",
        "target_seq",
        "target_chrom",
        "target_start",
        "target_end",
        "target_arm",
        "other_seq",
        "other_chrom",
        "other_start",
        "other_end",
        "other_arm",
        "relation",
        "signal_class",
        "intersection",
        "jaccard_similarity",
        "cosine_similarity",
        "dice_similarity",
        "estimated_identity",
    ]
    write_tsv(OUT / "summaries/per_window_target_similarity_support.tsv", best_rows, fields)
    inter = [row for row in edge_rows if row["relation"] == "interchromosomal"]
    write_tsv(OUT / "summaries/all_interchromosomal_targets.tsv", inter, fields)
    top_rows: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for row in sorted(
        inter,
        key=lambda item: (
            str(item["method"]),
            str(item["comparison_id"]),
            str(item["target_seq"]),
            int(item["target_start"]),
            -float(item["estimated_identity"]),
        ),
    ):
        key = (row["method"], row["comparison_id"], row["target_seq"], row["target_start"], row["target_end"])
        if key not in seen:
            seen.add(key)
            top_rows.append(row)
    write_tsv(OUT / "summaries/top_interchromosomal_targets.tsv", top_rows, fields)
    write_tsv(OUT / "summaries/chr9q_chr3q_windows.tsv", [r for r in inter if r["signal_class"] == "chr9q_chr3q"], fields)
    write_tsv(OUT / "summaries/par_controls.tsv", [r for r in inter if r["signal_class"] == "PAR_control"], fields)
    write_tsv(OUT / "summaries/acrocentric_controls.tsv", [r for r in inter if r["signal_class"] == "acrocentric_control"], fields)

    tracks: dict[tuple[object, ...], dict[str, object]] = defaultdict(dict)
    for row in inter:
        key = (row["method"], row["comparison_id"], row["target_seq"], row["target_start"], row["target_end"], row["target_arm"])
        rec = tracks[key]
        rec.update(
            {
                "method": row["method"],
                "comparison_id": row["comparison_id"],
                "target_seq": row["target_seq"],
                "target_start": row["target_start"],
                "target_end": row["target_end"],
                "target_arm": row["target_arm"],
            }
        )
        arms = set(str(rec.get("interchromosomal_arms", "")).split(",")) if rec.get("interchromosomal_arms") else set()
        arms.add(str(row["other_arm"]))
        rec["interchromosomal_arms"] = ",".join(sorted(arms))
        rec["interchromosomal_edge_count"] = int(rec.get("interchromosomal_edge_count", 0)) + 1
        rec["max_estimated_identity"] = max(float(rec.get("max_estimated_identity", 0)), float(row["estimated_identity"]))
    write_tsv(
        OUT / "summaries/full_genome_target_pattern_tracks.tsv",
        list(tracks.values()),
        [
            "method",
            "comparison_id",
            "target_seq",
            "target_start",
            "target_end",
            "target_arm",
            "interchromosomal_arms",
            "interchromosomal_edge_count",
            "max_estimated_identity",
        ],
    )
    return edge_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()

    shard_manifest = OUT / "manifests/shard_manifest.tsv"
    rows = read_tsv(shard_manifest)
    completion = check_shards(rows)
    completion_fields = list(completion[0].keys()) if completion else []
    write_tsv(OUT / "manifests/shard_completion_manifest.tsv", completion, completion_fields)

    incomplete = [row for row in completion if row["state"] != "OK"]
    if incomplete and not args.allow_incomplete:
        raise SystemExit(f"{len(incomplete)} shards are missing or incomplete; see manifests/shard_completion_manifest.tsv")

    complete_rows = [row for row in completion if row["state"] == "OK"]
    assembled_rows: list[dict[str, object]] = []
    by_pair: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in complete_rows:
        by_pair[(str(row["method"]), str(row["comparison_id"]))].append(row)
    expected_pairs = {(str(row["method"]), str(row["comparison_id"])) for row in rows}
    complete_pairs = set(by_pair)
    if expected_pairs != complete_pairs and not args.allow_incomplete:
        raise SystemExit("Not all method/comparison pairs have complete shards")

    for (method, comparison_id), pair_rows in sorted(by_pair.items()):
        expected_count = sum(1 for row in rows if row["method"] == method and row["comparison_id"] == comparison_id)
        if len(pair_rows) != expected_count and not args.allow_incomplete:
            raise SystemExit(f"{method} {comparison_id}: {len(pair_rows)} complete of {expected_count} shards")
        output = OUT / "outputs" / "assembled" / f"{method}.{comparison_id}.full_genome_2kb.impg_similarity.tsv.gz"
        data_rows = assemble_one(method, comparison_id, pair_rows, output)
        first = pair_rows[0]
        assembled_rows.append(
            {
                "method": method,
                "comparison_id": comparison_id,
                "assembled_output_tsv_gz": str(output),
                "completed_shards": len(pair_rows),
                "expected_shards": expected_count,
                "data_rows": data_rows,
                "source_raw_paf": first["source_raw_paf"],
                "impg_alignment_paf": first["impg_alignment_paf"],
                "query_fasta": first["query_fasta"],
                "target_fasta": first["target_fasta"],
                "full_target_bed": first["full_target_bed"],
            }
        )
    write_tsv(
        OUT / "manifests/assembled_outputs.tsv",
        assembled_rows,
        [
            "method",
            "comparison_id",
            "assembled_output_tsv_gz",
            "completed_shards",
            "expected_shards",
            "data_rows",
            "source_raw_paf",
            "impg_alignment_paf",
            "query_fasta",
            "target_fasta",
            "full_target_bed",
        ],
    )
    if len(assembled_rows) == len(expected_pairs):
        edge_rows = summarize_assembled(assembled_rows)
        print(f"Assembled {len(assembled_rows)} outputs and summarized {len(edge_rows)} target-other edges")
    else:
        print(f"Assembled {len(assembled_rows)} incomplete outputs; summaries skipped")


if __name__ == "__main__":
    main()
