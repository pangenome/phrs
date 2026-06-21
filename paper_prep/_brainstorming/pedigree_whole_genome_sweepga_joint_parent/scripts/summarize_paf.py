#!/usr/bin/env python3
import glob
import os

from common import PACKAGE_DIR, paf_records, write_tsv


def summarize_file(path):
    records = 0
    query_bp = 0
    aln_bp = 0
    matches = 0
    queries = set()
    targets = set()
    for f in paf_records(path):
        records += 1
        queries.add(f[0])
        targets.add(f[5])
        query_bp += max(0, int(f[3]) - int(f[2]))
        matches += int(f[9])
        aln_bp += int(f[10])
    return {
        "path": path,
        "records": str(records),
        "query_bp": str(query_bp),
        "aln_bp": str(aln_bp),
        "matches": str(matches),
        "mean_identity": ("%.6f" % (float(matches) / aln_bp)) if aln_bp else "0",
        "query_sequences": str(len(queries)),
        "target_sequences": str(len(targets)),
        "gzip_ok": "yes",
    }


def main():
    rows = []
    patterns = [
        os.path.join(PACKAGE_DIR, "raw_paf", "*.paf.gz"),
        os.path.join(PACKAGE_DIR, "chopped_paf", "*.paf.gz"),
        os.path.join(PACKAGE_DIR, "filtered_paf", "*.paf.gz"),
        os.path.join(PACKAGE_DIR, "filtered_paf_unchopped_control", "*.paf.gz"),
    ]
    for pattern in patterns:
        for path in sorted(glob.glob(pattern)):
            rows.append(summarize_file(path))
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"), rows, [
        "path",
        "records",
        "query_bp",
        "aln_bp",
        "matches",
        "mean_identity",
        "query_sequences",
        "target_sequences",
        "gzip_ok",
    ])


if __name__ == "__main__":
    main()
