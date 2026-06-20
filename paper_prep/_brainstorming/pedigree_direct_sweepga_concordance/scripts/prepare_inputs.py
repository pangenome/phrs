#!/usr/bin/env python3
import os
import subprocess
import sys

from common import COMPARISONS, FASTA, PACKAGE_DIR, read_tsv, write_tsv


def read_names(path, required_prefix=None):
    names = []
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if not name:
                continue
            if required_prefix and not name.startswith(required_prefix):
                continue
            names.append(name)
    return names


def write_name_file(path, names):
    with open(path, "w") as fh:
        for name in names:
            fh.write(name + "\n")


def faidx_extract(names_path, fasta_out):
    with open(fasta_out, "w") as out:
        subprocess.check_call(["samtools", "faidx", "-r", names_path, FASTA], stdout=out)


def main():
    rows = []
    for row in read_tsv(COMPARISONS):
        cid = row["comparison_id"]
        query_names = read_names(row["query_list"], row["query_haplotype"])
        target_names = read_names(row["target_list"], row["target_haplotype"])
        if not query_names:
            raise SystemExit("No query names for %s from %s prefix %s" % (cid, row["query_list"], row["query_haplotype"]))
        if not target_names:
            raise SystemExit("No target names for %s from %s prefix %s" % (cid, row["target_list"], row["target_haplotype"]))
        query_names_path = os.path.join(PACKAGE_DIR, "inputs", cid + ".query.names.txt")
        target_names_path = os.path.join(PACKAGE_DIR, "inputs", cid + ".target.names.txt")
        query_fa = os.path.join(PACKAGE_DIR, "inputs", cid + ".query.fa")
        target_fa = os.path.join(PACKAGE_DIR, "inputs", cid + ".target.fa")
        write_name_file(query_names_path, query_names)
        write_name_file(target_names_path, target_names)
        if not os.path.exists(query_fa) or os.path.getsize(query_fa) == 0:
            faidx_extract(query_names_path, query_fa)
        if not os.path.exists(target_fa) or os.path.getsize(target_fa) == 0:
            faidx_extract(target_names_path, target_fa)
        rows.append({
            "comparison_id": cid,
            "query_name_count": str(len(query_names)),
            "target_name_count": str(len(target_names)),
            "query_fasta": query_fa,
            "target_fasta": target_fa,
            "source_fasta": FASTA,
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "input_manifest.tsv"), rows,
              ["comparison_id", "query_name_count", "target_name_count", "query_fasta", "target_fasta", "source_fasta"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
