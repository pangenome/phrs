#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import sys

from common import CHOP_LENGTH, CHOP_OVERLAP, PACKAGE_DIR, RAW_SUFFIX, read_tsv, write_tsv


def open_maybe_gzip(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def split_record(fields, chop_length, overlap):
    q_start = int(fields[2])
    q_end = int(fields[3])
    t_start = int(fields[7])
    t_end = int(fields[8])
    matches = int(fields[9])
    aln_len = int(fields[10])
    span = q_end - q_start
    if span <= 0 or span <= chop_length:
        yield fields
        return

    step = chop_length - overlap
    if step <= 0:
        raise ValueError("chop_length must be greater than overlap")
    target_forward = t_end >= t_start
    emitted = 0
    offset = 0
    while offset < span:
        frag_q_start = q_start + offset
        frag_q_end = min(q_end, frag_q_start + chop_length)
        frag_span = frag_q_end - frag_q_start
        if frag_span <= 0:
            break
        frac_start = float(frag_q_start - q_start) / span
        frac_end = float(frag_q_end - q_start) / span
        if target_forward:
            frag_t_start = int(round(t_start + (t_end - t_start) * frac_start))
            frag_t_end = int(round(t_start + (t_end - t_start) * frac_end))
        else:
            frag_t_start = int(round(t_start - (t_start - t_end) * frac_start))
            frag_t_end = int(round(t_start - (t_start - t_end) * frac_end))
        frag_aln_len = max(1, int(round(aln_len * frag_span / span)))
        frag_matches = min(frag_aln_len, int(round(matches * frag_span / span)))
        out = list(fields)
        out[2] = str(frag_q_start)
        out[3] = str(frag_q_end)
        out[7] = str(frag_t_start)
        out[8] = str(frag_t_end)
        out[9] = str(frag_matches)
        out[10] = str(frag_aln_len)
        out.extend([
            "cg:Z:chopped",
            "zc:i:%d" % emitted,
            "zl:i:%d" % chop_length,
            "zo:i:%d" % overlap,
            "zs:i:%d" % q_start,
            "ze:i:%d" % q_end,
        ])
        yield out
        emitted += 1
        offset += step


def count_records(path):
    n = 0
    bp = 0
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                n += 1
                bp += max(0, int(fields[3]) - int(fields[2]))
    return n, bp


def chop_one(comparison_id, chop_length, overlap):
    raw = os.path.join(PACKAGE_DIR, "raw_paf", comparison_id + RAW_SUFFIX)
    out = os.path.join(PACKAGE_DIR, "chopped_paf", "%s.chopped_l%d_o%d.paf.gz" % (comparison_id, chop_length, overlap))
    if not os.path.exists(raw):
        raise SystemExit("missing raw whole-genome PAF %s" % raw)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    raw_records = 0
    chopped_records = 0
    raw_bp = 0
    chopped_bp = 0
    with gzip.open(raw, "rt") as inp, gzip.open(out, "wt") as oh:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            raw_records += 1
            raw_bp += max(0, int(fields[3]) - int(fields[2]))
            for frag in split_record(fields, chop_length, overlap):
                chopped_records += 1
                chopped_bp += max(0, int(frag[3]) - int(frag[2]))
                oh.write("\t".join(frag) + "\n")
    return {
        "comparison_id": comparison_id,
        "tool": "scripts/chop_paf.py deterministic query-axis splitter",
        "chop_length_bp": str(chop_length),
        "overlap_bp": str(overlap),
        "segmentation_rule": "split each raw PAF row into <= chop_length query-axis fragments; linearly interpolate target coordinates and match/alignment counts",
        "rationale": "bound merged whole-genome PAF intervals before sweepGA joint filtering so long merged paths do not break similarity/path metrics",
        "input_raw_paf": raw,
        "output_chopped_paf": out,
        "raw_records": str(raw_records),
        "chopped_records": str(chopped_records),
        "raw_query_bp": str(raw_bp),
        "chopped_query_bp": str(chopped_bp),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-id")
    parser.add_argument("--chop-length", type=int, default=CHOP_LENGTH)
    parser.add_argument("--overlap", type=int, default=CHOP_OVERLAP)
    args = parser.parse_args()
    comparisons = [r["comparison_id"] for r in read_tsv(os.path.join(PACKAGE_DIR, "config", "comparisons.tsv"))]
    if args.comparison_id:
        comparisons = [args.comparison_id]
    rows = [chop_one(cid, args.chop_length, args.overlap) for cid in comparisons]
    manifest = os.path.join(PACKAGE_DIR, "summaries", "chop_manifest.tsv")
    existing = []
    if os.path.exists(manifest) and args.comparison_id:
        import csv
        with open(manifest, newline="") as fh:
            existing = [r for r in csv.DictReader(fh, delimiter="\t") if r["comparison_id"] != args.comparison_id]
    rows = existing + rows
    write_tsv(manifest, rows, [
        "comparison_id",
        "tool",
        "chop_length_bp",
        "overlap_bp",
        "segmentation_rule",
        "rationale",
        "input_raw_paf",
        "output_chopped_paf",
        "raw_records",
        "chopped_records",
        "raw_query_bp",
        "chopped_query_bp",
    ])
    return 0


if __name__ == "__main__":
    sys.exit(main())
