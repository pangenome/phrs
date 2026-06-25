#!/usr/bin/env python3
"""Summarize completed IMPG similarity TSV outputs for Fig5 raw many:many scans."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


OUT = Path(__file__).resolve().parent
CHR_RE = re.compile(r"_chr([0-9XY]+)(?::|$)")


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_group(group: str) -> dict[str, str | int]:
    seq, coords = group.rsplit(":", 1)
    start_text, end_text = coords.split("-", 1)
    chrom_match = CHR_RE.search(seq)
    chrom = f"chr{chrom_match.group(1)}" if chrom_match else "unknown"
    return {
        "seq": seq,
        "chrom": chrom,
        "start": int(start_text),
        "end": int(end_text),
    }


def arm(chrom: str, start: int, end: int, seq_len: int | None) -> str:
    if chrom == "unknown":
        return "unknown"
    if seq_len is None:
        return chrom
    mid = (start + end) // 2
    return f"{chrom}{'p' if mid < seq_len / 2 else 'q'}"


def read_lengths(target_fasta: str) -> dict[str, int]:
    fai = Path(f"{target_fasta}.fai")
    lengths: dict[str, int] = {}
    with fai.open() as fh:
        for line in fh:
            seq, length, *_ = line.rstrip("\n").split("\t")
            lengths[seq] = int(length)
    return lengths


def write_table(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def load_manifest() -> list[dict[str, str]]:
    with (OUT / "run_manifest.tsv").open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def row_target_side(row: dict[str, str]) -> tuple[dict[str, str | int], dict[str, str | int]] | None:
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    a = parse_group(row["group.a"])
    b = parse_group(row["group.b"])
    target = {"seq": chrom, "chrom": parse_group(f"{chrom}:{start}-{end}")["chrom"], "start": start, "end": end}
    if a["seq"] == chrom and a["start"] == start and a["end"] == end:
        return target, b
    if b["seq"] == chrom and b["start"] == start and b["end"] == end:
        return target, a
    return None


def summarize_one(manifest_row: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    path = Path(manifest_row["output_tsv_gz"])
    if not path.exists():
        return [], []

    lengths = read_lengths(manifest_row["target_fasta"])
    best_by_window: dict[tuple[str, int, int, str], dict[str, object]] = {}
    all_edges: list[dict[str, object]] = []
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
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
                "method": manifest_row["method"],
                "comparison_id": manifest_row["comparison_id"],
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
            all_edges.append(out)
            key = (str(target["seq"]), int(target["start"]), int(target["end"]), relation)
            previous = best_by_window.get(key)
            if previous is None or float(out["estimated_identity"]) > float(previous["estimated_identity"]):
                best_by_window[key] = out
    return list(best_by_window.values()), all_edges


def annotate_signal(row: dict[str, object]) -> str:
    target_arm = str(row["target_arm"])
    other_arm = str(row["other_arm"])
    pair = {target_arm, other_arm}
    if pair == {"chrXq", "chrYq"} or pair == {"chrXp", "chrYp"}:
        return "PAR_control"
    if target_arm in {"chr13p", "chr14p", "chr15p", "chr21p", "chr22p"} or other_arm in {
        "chr13p",
        "chr14p",
        "chr15p",
        "chr21p",
        "chr22p",
    }:
        return "acrocentric_control"
    if pair == {"chr9q", "chr3q"}:
        return "chr9q_chr3q"
    return "other"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-missing", action="store_true")
    args = parser.parse_args()

    manifest = load_manifest()
    missing = [row["output_tsv_gz"] for row in manifest if not Path(row["output_tsv_gz"]).exists()]
    if missing and not args.allow_missing:
        raise SystemExit("Missing IMPG outputs:\n" + "\n".join(missing))

    best_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    for row in manifest:
        best, edges = summarize_one(row)
        best_rows.extend(best)
        edge_rows.extend(edges)

    for row in best_rows + edge_rows:
        row["signal_class"] = annotate_signal(row)

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
    write_table(OUT / "summary/per_window_target_similarity_support.tsv", best_rows, fields)
    inter = [row for row in edge_rows if row["relation"] == "interchromosomal"]
    write_table(OUT / "summary/all_interchromosomal_targets.tsv", inter, fields)
    top = sorted(
        inter,
        key=lambda r: (str(r["method"]), str(r["comparison_id"]), str(r["target_seq"]), int(r["target_start"]), -float(r["estimated_identity"])),
    )
    seen: set[tuple[object, ...]] = set()
    top_rows: list[dict[str, object]] = []
    for row in top:
        key = (row["method"], row["comparison_id"], row["target_seq"], row["target_start"], row["target_end"])
        if key in seen:
            continue
        seen.add(key)
        top_rows.append(row)
    write_table(OUT / "summary/top_interchromosomal_targets.tsv", top_rows, fields)
    write_table(OUT / "summary/chr9q_chr3q_windows.tsv", [r for r in inter if r["signal_class"] == "chr9q_chr3q"], fields)
    write_table(OUT / "summary/par_controls.tsv", [r for r in inter if r["signal_class"] == "PAR_control"], fields)
    write_table(OUT / "summary/acrocentric_controls.tsv", [r for r in inter if r["signal_class"] == "acrocentric_control"], fields)

    track_counts: dict[tuple[object, ...], dict[str, object]] = defaultdict(dict)
    for row in inter:
        key = (row["method"], row["comparison_id"], row["target_seq"], row["target_start"], row["target_end"], row["target_arm"])
        rec = track_counts[key]
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
    write_table(
        OUT / "summary/full_genome_target_pattern_tracks.tsv",
        list(track_counts.values()),
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
    print(f"Summarized {len(edge_rows)} target-other similarity edges from {len(manifest) - len(missing)} completed outputs")


if __name__ == "__main__":
    main()
