#!/usr/bin/env python3
import csv
import hashlib
import os
import subprocess
import sys

from common import COMPARISONS, PACKAGE_DIR, write_tsv


LENGTHS = [int(x) for x in os.environ.get("QUERY_GRID_CHOP_LENGTHS", "10000 5000 2000 1000").split()]
FILTER_ID = "one_to_one_ani_o0"
STATUS_DIR = os.path.join(PACKAGE_DIR, "summaries", "query_grid_chop_filter_status")
SLURM_SUMMARY = os.path.join(PACKAGE_DIR, "summaries", "query_grid_chop_filter_slurm.tsv")
MANIFEST = os.path.join(PACKAGE_DIR, "summaries", "query_grid_chop_filter_manifest.tsv")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd):
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def sacct(job_id):
    out = run(["sacct", "-n", "-P", "-j", job_id, "--format=JobID,State,Elapsed,ExitCode,NodeList"])
    for line in out.splitlines():
        parts = line.split("|")
        if parts and parts[0] == job_id:
            return {
                "sacct_state": parts[1] if len(parts) > 1 else "",
                "elapsed": parts[2] if len(parts) > 2 else "",
                "exit_code": parts[3] if len(parts) > 3 else "",
                "node": parts[4] if len(parts) > 4 else "",
            }
    return {"sacct_state": "", "elapsed": "", "exit_code": "", "node": ""}


def main():
    comparisons = [row["comparison_id"] for row in read_tsv(COMPARISONS)]
    status_by_key = {}
    if os.path.isdir(STATUS_DIR):
        for name in sorted(os.listdir(STATUS_DIR)):
            if not name.endswith(".tsv"):
                continue
            rows = read_tsv(os.path.join(STATUS_DIR, name))
            if rows:
                row = rows[0]
                status_by_key[(row["comparison_id"], int(row["chop_length_bp"]))] = row

    slurm_rows = []
    manifest_rows = []
    failures = []
    for cid in comparisons:
        for length in LENGTHS:
            chopped = os.path.join(PACKAGE_DIR, f"chopped_paf_qgrid_l{length}_o0", f"{cid}.chopped_l{length}_o0_query_grid.paf.gz")
            filtered = os.path.join(PACKAGE_DIR, "filtered_paf_chop_sensitivity_query_grid", f"l{length}", f"{cid}.{FILTER_ID}.chopped_l{length}_o0_query_grid.paf.gz")
            chop_summary = os.path.join(PACKAGE_DIR, "summaries", f"pafchop_qgrid_l{length}_o0", f"{cid}.summary.tsv")
            status = status_by_key.get((cid, length), {})
            job_id = status.get("job_id", "")
            sacct_info = sacct(job_id) if job_id and job_id != "manual" else {"sacct_state": "", "elapsed": "", "exit_code": "", "node": ""}

            chopped_ok = os.path.exists(chopped) and run(["pigz", "-t", chopped]) == ""
            filtered_ok = os.path.exists(filtered) and run(["pigz", "-t", filtered]) == ""
            chopped_sha = sha256(chopped) if chopped_ok else ""
            filtered_sha = sha256(filtered) if filtered_ok else ""
            worker_status = status.get("status", "")
            manifest_status = "OK" if chopped_ok and filtered_ok and worker_status in {"OK", "FAILED"} else "MISSING_OR_INVALID"
            if manifest_status != "OK":
                failures.append(f"{cid} length {length}: {manifest_status}")
            elif worker_status == "FAILED":
                status["status"] = "OK_VALIDATED_OUTPUTS"
            slurm_row = dict(status)
            slurm_row.update(sacct_info)
            slurm_rows.append(slurm_row)
            manifest_rows.append({
                "comparison_id": cid,
                "chop_length_bp": str(length),
                "filter_id": FILTER_ID,
                "chunk_mode": "query-grid",
                "output_family": "query_grid_chop_sensitivity",
                "raw_paf": status.get("raw_paf", ""),
                "chopped_paf": chopped,
                "chopped_sha256": chopped_sha,
                "filtered_paf": filtered,
                "filtered_sha256": filtered_sha,
                "chop_summary": chop_summary,
                "num_mappings": "1:1",
                "scaffold_jump": "0",
                "scoring": "ani",
                "filter_overlap": "0",
                "gzip_validation": "pigz -t",
                "status": manifest_status,
            })
            for path, digest in ((chopped, chopped_sha), (filtered, filtered_sha)):
                if digest:
                    with open(path + ".sha256", "w") as fh:
                        fh.write(f"{digest}  {path}\n")

    slurm_fields = [
        "comparison_id", "chop_length_bp", "filter_id", "job_id", "array_task_id", "host",
        "started_utc", "finished_utc", "status", "sacct_state", "elapsed", "exit_code", "node",
        "raw_paf", "chopped_paf", "filtered_paf", "chop_summary", "chunk_mode", "overlap_bp",
        "num_mappings", "scaffold_jump", "scoring", "filter_overlap", "threads", "pigz_compression_level", "scratch_dir",
        "pafchop_bin", "pafchop_sha256", "sweepga_bin", "sweepga_sha256", "pigz_bin", "pigz_sha256",
        "chop_command", "filter_command",
    ]
    manifest_fields = [
        "comparison_id", "chop_length_bp", "filter_id", "chunk_mode", "output_family",
        "raw_paf", "chopped_paf", "chopped_sha256", "filtered_paf", "filtered_sha256",
        "chop_summary", "num_mappings", "scaffold_jump", "scoring", "filter_overlap",
        "gzip_validation", "status",
    ]
    write_tsv(SLURM_SUMMARY, slurm_rows, slurm_fields)
    write_tsv(MANIFEST, manifest_rows, manifest_fields)
    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    print(f"Validated {len(manifest_rows)} query-grid chopped/filter outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
