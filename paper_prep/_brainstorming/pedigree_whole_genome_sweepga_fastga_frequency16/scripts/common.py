#!/usr/bin/env python3
import csv
import gzip
import os


PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COMPARISONS = os.path.join(PACKAGE_DIR, "config", "comparisons.tsv")
FILTER_MATRIX = os.path.join(PACKAGE_DIR, "config", "filter_matrix.tsv")
RAW_SUFFIX = ".sweepga_frequency16_many_many_j0.paf.gz"
CHOP_LENGTH = int(os.environ.get("PAF_CHOP_LENGTH", "500000"))
CHOP_OVERLAP = int(os.environ.get("PAF_CHOP_OVERLAP", "0"))


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
