#!/usr/bin/env python3
import argparse
import gzip
import os
import subprocess
import sys
import tempfile

from common import FILTER_MATRIX, PACKAGE_DIR, read_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--filter-id", required=True)
    args = parser.parse_args()

    filters = {row["filter_id"]: row for row in read_tsv(FILTER_MATRIX)}
    if args.filter_id not in filters:
        raise SystemExit("unknown filter_id %s" % args.filter_id)
    filt = filters[args.filter_id]
    raw = os.path.join(PACKAGE_DIR, "raw_paf", "%s.sweepga_many_many_j0.paf.gz" % args.comparison_id)
    out = os.path.join(PACKAGE_DIR, "filtered_paf", "%s.%s.paf.gz" % (args.comparison_id, args.filter_id))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(raw):
        raise SystemExit("missing raw PAF %s" % raw)

    if args.filter_id.startswith("simple_"):
        min_ident = float(filt["min_identity"])
        min_len = int(filt["min_aln_length"])
        min_qcov = float(filt["min_query_coverage"])
        with gzip.open(raw, "rt") as inp, gzip.open(out, "wt") as oh:
            for line in inp:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                q_len = int(fields[1])
                q_start = int(fields[2])
                q_end = int(fields[3])
                matches = int(fields[9])
                aln_len = int(fields[10])
                identity = float(matches) / aln_len if aln_len else 0.0
                qcov = float(q_end - q_start) / q_len if q_len else 0.0
                if identity >= min_ident and aln_len >= min_len and qcov >= min_qcov:
                    oh.write(line)
        return 0

    sweepga = os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga")
    with tempfile.NamedTemporaryFile(suffix=".paf", delete=False) as tmp:
        tmp_path = tmp.name
    filtered_tmp = tmp_path + ".filtered.paf"
    try:
        with gzip.open(raw, "rb") as inp, open(tmp_path, "wb") as tmp:
            tmp.write(inp.read())
        cmd = [
            sweepga,
            "--num-mappings", filt["num_mappings"],
            "--scaffold-jump", filt["scaffold_jump"],
            "--output-file", filtered_tmp,
            tmp_path,
        ]
        subprocess.check_call(cmd)
        with open(filtered_tmp, "rb") as inp, gzip.open(out, "wb") as oh:
            oh.write(inp.read())
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        try:
            os.unlink(filtered_tmp)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
