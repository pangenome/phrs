#!/usr/bin/env python3
import csv
import os
import subprocess
import tempfile

from common import PACKAGE_DIR


OUT = os.path.join(PACKAGE_DIR, "summaries", "query_grid_shifted_boundary_audit.tsv")
PAFCHOP = os.environ.get(
    "PAFCHOP_BIN",
    os.path.abspath(os.path.join(PACKAGE_DIR, "..", "..", "..", "target", "release", "pafchop")),
)


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="qgrid_boundary.") as tmp:
        raw = os.path.join(tmp, "shifted.paf")
        summary = os.path.join(tmp, "summary.tsv")
        with open(raw, "w") as fh:
            fh.write("qA\t100\t7\t27\t+\ttA\t1000\t100\t120\t20\t20\t60\tcg:Z:20=\tri:Z:row_a\n")
            fh.write("qA\t100\t12\t32\t+\ttB\t1000\t200\t220\t20\t20\t60\tcg:Z:20=\tri:Z:row_b\n")
        proc = subprocess.run(
            [
                PAFCHOP,
                "--input", raw,
                "--length", "10",
                "--overlap", "0",
                "--chunk-mode", "query-grid",
                "--comparison-id", "shifted_same_query",
                "--summary", summary,
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
    rows = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        tags = {part.split(":", 2)[0]: part for part in f[12:] if ":" in part}
        rows.append({
            "comparison_id": "shifted_same_query_f32",
            "query_name": f[0],
            "target_name": f[5],
            "target_start": f[7],
            "chunk_q_start": f[2],
            "chunk_q_end": f[3],
            "chop_length_bp": "10",
            "chunk_mode": "query-grid",
            "shared_query_grid_boundary": "yes" if f[2] in {"10", "20", "30"} or f[3] in {"10", "20", "30"} else "edge",
            "mode_tag": tags.get("zm", ""),
            "provenance": "synthetic f32 shifted raw mappings row_a 7-27 and row_b 12-32 on qA; both are cut at absolute qA grid boundary 20",
        })
    with open(OUT, "w", newline="") as fh:
        fieldnames = [
            "comparison_id", "query_name", "target_name", "target_start",
            "chunk_q_start", "chunk_q_end", "chop_length_bp", "chunk_mode",
            "shared_query_grid_boundary", "mode_tag", "provenance",
        ]
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)


if __name__ == "__main__":
    main()
