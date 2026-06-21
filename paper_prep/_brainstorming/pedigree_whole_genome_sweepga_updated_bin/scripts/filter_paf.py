#!/usr/bin/env python3
import argparse
import gzip
import os
import shutil
import subprocess
import tempfile

from common import FILTER_MATRIX, PACKAGE_DIR, read_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--filter-id", required=True)
    parser.add_argument("--input-dir", default="chopped_paf")
    parser.add_argument("--output-dir", default="filtered_paf")
    parser.add_argument("--chop-length", default=os.environ.get("PAF_CHOP_LENGTH", "500000"))
    parser.add_argument("--overlap", default=os.environ.get("PAF_CHOP_OVERLAP", "0"))
    args = parser.parse_args()

    filters = {row["filter_id"]: row for row in read_tsv(FILTER_MATRIX)}
    if args.filter_id not in filters:
        raise SystemExit("unknown filter_id %s" % args.filter_id)
    filt = filters[args.filter_id]

    src = os.path.join(
        PACKAGE_DIR,
        args.input_dir,
        "%s.chopped_l%s_o%s.paf.gz" % (args.comparison_id, args.chop_length, args.overlap),
    )
    out = os.path.join(PACKAGE_DIR, args.output_dir, "%s.%s.paf.gz" % (args.comparison_id, args.filter_id))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(src):
        raise SystemExit("missing chopped PAF %s" % src)

    if filt["num_mappings"] == "many:many":
        shutil.copyfile(src, out)
        return 0

    sweepga = os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga")
    scratch_base = os.environ.get("SWEEPGA_DEVSHM_BASE") or "/dev/shm"
    if not os.path.isdir(scratch_base) or not os.access(scratch_base, os.W_OK):
        raise SystemExit("required sweepGA filter scratch base is not writable: %s" % scratch_base)
    work = tempfile.mkdtemp(prefix="sweepga_filter.%s." % args.comparison_id, dir=scratch_base)
    try:
        tmp_paf = os.path.join(work, "input.paf")
        filtered_tmp = os.path.join(work, "filtered.paf")
        with gzip.open(src, "rb") as inp, open(tmp_paf, "wb") as tmp:
            shutil.copyfileobj(inp, tmp)
        cmd = [
            sweepga,
            "--num-mappings", filt["num_mappings"],
            "--scaffold-jump", filt["scaffold_jump"],
            "--output-file", filtered_tmp,
            tmp_paf,
        ]
        subprocess.check_call(cmd)
        with open(filtered_tmp, "rb") as inp, gzip.open(out, "wb") as oh:
            shutil.copyfileobj(inp, oh)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
