#!/usr/bin/env python3
import csv
import gzip
import os


PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COMPARISONS = os.path.join(PACKAGE_DIR, "config", "comparisons.tsv")
PARAMETERS = os.path.join(PACKAGE_DIR, "config", "wfmash_parameters.tsv")
CANDIDATE_WINDOWS = os.path.join(PACKAGE_DIR, "config", "candidate_windows.tsv")
RAW_SUFFIX = ".wfmash.paf.gz"
DEFAULT_RUN_LABEL = os.environ.get("WFMASH_RUN_LABEL", "current_v0.24.2")
SUMMARY_ALL_RUNS = os.environ.get("WFMASH_SUMMARY_ALL_RUNS", "") == "1"


def paf_glob():
    if SUMMARY_ALL_RUNS:
        return os.path.join(PACKAGE_DIR, "raw_paf", "**", "*.paf.gz")
    return os.path.join(PACKAGE_DIR, "raw_paf", DEFAULT_RUN_LABEL, "*.paf.gz")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def paf_records(path):
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            yield fields
