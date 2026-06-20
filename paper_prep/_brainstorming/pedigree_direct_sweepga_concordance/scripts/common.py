#!/usr/bin/env python3
import csv
import gzip
import os
import re


PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FASTA = "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz"
COMPARISONS = os.path.join(PACKAGE_DIR, "config", "comparisons.tsv")
FILTER_MATRIX = os.path.join(PACKAGE_DIR, "config", "filter_matrix.tsv")


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


def sequence_arm(name):
    m = re.search(r"_(chr(?:[0-9]+|X|Y)_[pq]arm)$", name)
    if not m:
        return "unknown"
    return m.group(1).replace("_parm", "p").replace("_qarm", "q")


def sample_hap(name):
    return name.split("#chr", 1)[0]


def native_window(name):
    m = re.search(r"#(chr(?:[0-9]+|X|Y))\.[^:]+:(\d+)-(\d+)_", name)
    if not m:
        return ""
    chrom, start, end_inclusive = m.groups()
    return "%s:%d-%d" % (chrom, int(start), int(end_inclusive) + 1)


def paf_records(path):
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            matches = int(fields[9])
            aln_len = int(fields[10])
            q_len = int(fields[1])
            q_start = int(fields[2])
            q_end = int(fields[3])
            yield {
                "query_name": fields[0],
                "query_length": q_len,
                "query_start": q_start,
                "query_end": q_end,
                "strand": fields[4],
                "target_name": fields[5],
                "target_length": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matches": matches,
                "aln_length": aln_len,
                "mapq": fields[11],
                "identity": (float(matches) / aln_len) if aln_len else 0.0,
                "query_coverage": (float(q_end - q_start) / q_len) if q_len else 0.0,
                "query_arm": sequence_arm(fields[0]),
                "target_arm": sequence_arm(fields[5]),
                "query_hap": sample_hap(fields[0]),
                "target_hap": sample_hap(fields[5]),
                "query_native_window": native_window(fields[0]),
                "target_native_window": native_window(fields[5]),
            }


def overlap_bp(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))
