#!/usr/bin/env python3
import csv
import os
import subprocess

from common import PACKAGE_DIR, write_tsv


JOBS = os.path.join(PACKAGE_DIR, "summaries", "slurm_jobs.tsv")
FIELDS = [
    "stage",
    "comparison_id",
    "job_id",
    "parameter_set",
    "fastga_frequency",
    "num_mappings",
    "scaffold_jump",
    "sweepga_devshm_base",
    "stdout",
    "command_log",
    "sacct_state",
    "elapsed",
    "exit_code",
    "node",
    "status",
    "note",
]


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
    with open(JOBS, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for row in rows:
        cid = row["comparison_id"]
        job_id = row["job_id"]
        row.update(sacct(job_id))
        command_logs = sorted(
            p for p in os.listdir(os.path.join(PACKAGE_DIR, "logs"))
            if p.startswith(cid + ".sweepga_frequency16_many_many_j0.") and p.endswith(".command.log")
        )
        row["command_log"] = os.path.join(PACKAGE_DIR, "logs", command_logs[-1]) if command_logs else ""
        raw = os.path.join(PACKAGE_DIR, "raw_paf", cid + ".sweepga_frequency16_many_many_j0.paf.gz")
        if os.path.exists(raw) and os.path.getsize(raw) > 0:
            row["status"] = "RAW_PAF_OK"
            row["note"] = "Raw PAF emitted by explicit frequency16 run."
        elif "CANCELLED" in row.get("sacct_state", ""):
            row["status"] = "PATHOLOGICAL_NO_RAW_PAF"
            row["note"] = "Cancelled after active FastGA -f16 runtime with no raw PAF; see pathological_runtime.tsv if present."
        else:
            row["status"] = "NO_RAW_PAF"
            row["note"] = "No raw PAF present at summary time."
        out.append(row)
    write_tsv(JOBS, out, FIELDS)


if __name__ == "__main__":
    main()
