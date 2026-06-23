#!/usr/bin/env python3
import argparse
import csv
import gzip


def read_windows(path, comparison_id):
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["comparison_id"] == comparison_id:
                row["query_start"] = int(row["query_start"])
                row["query_end"] = int(row["query_end"])
                yield row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--raw-paf", required=True)
    args = parser.parse_args()

    windows = list(read_windows(f"{args.panel_dir}/config/panel_windows.tsv", args.comparison_id))
    if not windows:
        raise SystemExit(f"no configured windows for {args.comparison_id}")

    with gzip.open(args.raw_paf, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q_name = fields[0]
            q_start = int(fields[2])
            q_end = int(fields[3])
            for window in windows:
                if q_name != window["query_name"]:
                    continue
                if min(q_end, window["query_end"]) > max(q_start, window["query_start"]):
                    print(line, end="")
                    break


if __name__ == "__main__":
    main()
