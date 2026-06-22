#!/usr/bin/env python3
import gzip
import hashlib
import os
import shutil
import subprocess

from common import DEFAULT_RUN_LABEL, PACKAGE_DIR, read_tsv, write_tsv


PARAMETER_SET = os.environ.get("MINIMAP2_PARAMETER_SET", "asm5_allchains")
MINIMAP2_SUFFIX = os.environ.get("MINIMAP2_SUFFIX", "minimap2-v2.31-r1302")


def job_state(job_id):
    try:
        out = subprocess.check_output(
            ["sacct", "-n", "-P", "-j", str(job_id), "--format=JobID,State,ExitCode,Elapsed,Timelimit"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return "UNKNOWN"
    for line in out.splitlines():
        fields = line.split("|")
        if fields and fields[0] == str(job_id):
            return ":".join(fields[1:])
    return "UNKNOWN"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return "yes"
    except Exception as exc:
        return "no:%s" % exc


def count_records(path):
    records = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                records += 1
    return records


def destination_path(comparison_id):
    return os.path.join(
        PACKAGE_DIR,
        "raw_paf",
        DEFAULT_RUN_LABEL,
        "%s.%s.%s.paf.gz" % (comparison_id, PARAMETER_SET, MINIMAP2_SUFFIX),
    )


def main():
    slurm_path = os.path.join(PACKAGE_DIR, "summaries", "slurm_jobs.tsv")
    rows = []
    for job in read_tsv(slurm_path):
        source = job["output_paf"]
        dest = destination_path(job["comparison_id"])
        state = job_state(job["job_id"])
        row = {
            "comparison_id": job["comparison_id"],
            "job_id": job["job_id"],
            "sacct_state": state,
            "source_paf": source,
            "dest_paf": dest,
            "source_exists": "yes" if os.path.exists(source) else "no",
            "source_bytes": "",
            "source_gzip_ok": "no",
            "source_records": "0",
            "dest_bytes": "",
            "dest_sha256": "",
            "status": "SOURCE_MISSING",
        }
        if not os.path.exists(source):
            rows.append(row)
            continue
        row["source_bytes"] = str(os.path.getsize(source))
        ok = gzip_ok(source)
        row["source_gzip_ok"] = ok
        if ok != "yes":
            row["status"] = "SOURCE_GZIP_INCOMPLETE"
            rows.append(row)
            continue
        records = count_records(source)
        row["source_records"] = str(records)
        if records == 0:
            row["status"] = "SOURCE_EMPTY"
            rows.append(row)
            continue
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(source, dest)
        digest = sha256(dest)
        with open(dest + ".sha256", "w") as fh:
            fh.write("%s  %s\n" % (digest, dest))
        row["dest_bytes"] = str(os.path.getsize(dest))
        row["dest_sha256"] = digest
        row["status"] = "COPIED"
        rows.append(row)

    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "raw_paf_harvest_manifest.tsv"), rows, [
        "comparison_id",
        "job_id",
        "sacct_state",
        "source_paf",
        "dest_paf",
        "source_exists",
        "source_bytes",
        "source_gzip_ok",
        "source_records",
        "dest_bytes",
        "dest_sha256",
        "status",
    ])


if __name__ == "__main__":
    main()
