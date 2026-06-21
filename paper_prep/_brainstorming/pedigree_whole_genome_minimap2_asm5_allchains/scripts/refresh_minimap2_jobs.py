#!/usr/bin/env python3
import os
import subprocess

from common import PACKAGE_DIR, read_tsv, write_tsv


def job_state(job_id):
    try:
        out = subprocess.check_output(
            ["sacct", "-n", "-P", "-j", str(job_id), "--format=JobID,State,ExitCode,Elapsed,NodeList"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return "UNKNOWN"
    states = []
    for line in out.splitlines():
        fields = line.split("|")
        if fields and fields[0] == str(job_id):
            states.append(":".join(fields[1:4]))
    return ";".join(states) if states else "UNKNOWN"


def file_info(path):
    if not os.path.exists(path):
        return "", ""
    size = str(os.path.getsize(path))
    sha_path = path + ".sha256"
    digest = ""
    if os.path.exists(sha_path):
        with open(sha_path) as fh:
            digest = fh.readline().split()[0]
    else:
        digest = subprocess.check_output(["sha256sum", path], text=True).split()[0]
    return size, digest


def main():
    path = os.path.join(PACKAGE_DIR, "summaries", "slurm_jobs.tsv")
    rows = read_tsv(path)
    for row in rows:
        size, digest = file_info(row["output_paf"])
        row["output_size_bytes"] = size
        row["sha256"] = digest
        state = job_state(row["job_id"])
        if size and digest:
            row["status"] = "COMPLETED_OUTPUT_PRESENT"
        elif state != "UNKNOWN":
            row["status"] = state
    write_tsv(path, rows, [
        "run_label",
        "comparison_id",
        "parameter_set",
        "job_id",
        "status",
        "stdout",
        "stderr",
        "output_paf",
        "output_size_bytes",
        "sha256",
        "minimap2_options",
        "scratch_base",
        "minimap2_bin",
    ])


if __name__ == "__main__":
    main()
