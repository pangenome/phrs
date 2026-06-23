#!/usr/bin/env python3
import csv
import os
import re
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_label(value):
    value = value.replace(":", "to").replace("-", "_")
    value = re.sub(r"[^A-Za-z0-9_.]+", "_", value)
    return value.strip("_")


def main():
    repo = Path(os.environ.get("REPO_ROOT", "/moosefs/erikg/phrs"))
    package = repo / "paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity"
    sources = {row["comparison_id"]: row["source_paf"] for row in read_tsv(package / "config/source_pafs.tsv")}
    windows = read_tsv(package / "config/candidate_windows.tsv")
    matrix = [row for row in read_tsv(package / "config/matrix.tsv") if row.get("include", "yes") == "yes"]

    comparison_order = []
    for row in windows:
        if row["comparison_id"] not in comparison_order:
            comparison_order.append(row["comparison_id"])

    out_path = package / "filter_tasks.tsv"
    fields = [
        "task_index",
        "comparison_id",
        "source_paf",
        "output_paf",
        "cell_id",
        "scaffold_jump",
        "min_aln_length",
        "scoring",
        "num_mappings",
        "scaffold_mass",
        "overlap",
        "command",
    ]
    rows = []
    task_index = 0
    for comparison_id in comparison_order:
        source = sources[comparison_id]
        for cell in matrix:
            task_index += 1
            sj = cell["scaffold_jump"]
            min_len = cell["min_aln_length"]
            scoring = cell["scoring"]
            mappings = cell["num_mappings"]
            mass = cell["scaffold_mass"]
            cell_id = ".".join(
                [
                    f"j{safe_label(sj)}",
                    f"l{safe_label(min_len)}",
                    safe_label(scoring),
                    safe_label(mappings),
                    f"m{safe_label(mass)}",
                ]
            )
            output = package / "filtered_paf" / f"{comparison_id}.{cell_id}.paf.gz"
            command = [
                "/home/erikg/.cargo/bin/sweepga",
                "--num-mappings",
                mappings,
                "--scaffold-jump",
                sj,
                "--scaffold-mass",
                mass,
                "--scoring",
                scoring,
            ]
            if min_len != "default":
                command.extend(["--min-aln-length", min_len])
            if cell["overlap"] != "default":
                command.extend(["--overlap", cell["overlap"]])
            command.extend(["--threads", "${THREADS}", "--output-file", "${WORK}/filtered.paf", "${WORK}/input.paf"])
            rows.append(
                {
                    "task_index": str(task_index),
                    "comparison_id": comparison_id,
                    "source_paf": source,
                    "output_paf": str(output),
                    "cell_id": cell_id,
                    "scaffold_jump": sj,
                    "min_aln_length": min_len,
                    "scoring": scoring,
                    "num_mappings": mappings,
                    "scaffold_mass": mass,
                    "overlap": cell["overlap"],
                    "command": " ".join(command),
                }
            )

    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} tasks to {out_path}")


if __name__ == "__main__":
    main()
