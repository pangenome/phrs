#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_sidecar_sha256(path):
    sidecar = path + ".sha256"
    if not os.path.exists(sidecar):
        return ""
    with open(sidecar) as handle:
        return handle.read().strip().split()[0]


def run_filter(raw_paf, output_gz, mode, sweepga, scratch_base, log_path):
    os.makedirs(os.path.dirname(output_gz), exist_ok=True)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    work = tempfile.mkdtemp(prefix="raw_no_chop_filter.", dir=scratch_base)
    started = datetime.now(timezone.utc)
    tmp_paf = os.path.join(work, "filtered.paf")
    cmd = [
        sweepga,
        "--num-mappings", mode["num_mappings"],
        "--scaffold-jump", mode["scaffold_jump"],
        "--scoring", mode["scoring"],
        "--output-file", tmp_paf,
    ]
    with open(log_path, "w") as log:
        log.write("started_utc\t%s\n" % started.isoformat())
        log.write("raw_paf\t%s\n" % raw_paf)
        log.write("output_gz\t%s\n" % output_gz)
        log.write("command\tgzip -dc %s | %s\n" % (raw_paf, " ".join(cmd)))
        log.flush()
        try:
            gzip_proc = subprocess.Popen(["gzip", "-dc", raw_paf], stdout=subprocess.PIPE, stderr=log)
            try:
                proc = subprocess.run(cmd, stdin=gzip_proc.stdout, stdout=log, stderr=log, check=False)
            finally:
                if gzip_proc.stdout is not None:
                    gzip_proc.stdout.close()
            gzip_returncode = gzip_proc.wait()
            log.write("gzip_returncode\t%s\n" % gzip_returncode)
            log.write("sweepga_returncode\t%s\n" % proc.returncode)
            if gzip_returncode != 0:
                raise subprocess.CalledProcessError(gzip_returncode, ["gzip", "-dc", raw_paf])
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, cmd)
            with open(tmp_paf, "rb") as inp, gzip.open(output_gz, "wb", compresslevel=6) as out:
                shutil.copyfileobj(inp, out)
            finished = datetime.now(timezone.utc)
            log.write("finished_utc\t%s\n" % finished.isoformat())
            log.write("elapsed_seconds\t%.3f\n" % (finished - started).total_seconds())
        finally:
            shutil.rmtree(work, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--source-package", required=True)
    parser.add_argument("--sweepga", default=os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga"))
    parser.add_argument("--scratch-base", default=os.environ.get("SWEEPGA_DEVSHM_BASE") or os.environ.get("TMPDIR") or "/tmp")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.sweepga):
        raise SystemExit("missing sweepga binary: %s" % args.sweepga)
    if not os.path.isdir(args.scratch_base) or not os.access(args.scratch_base, os.W_OK):
        raise SystemExit("scratch base is not writable: %s" % args.scratch_base)

    windows = read_tsv(os.path.join(args.panel_dir, "config", "panel_windows.tsv"))
    modes = read_tsv(os.path.join(args.panel_dir, "config", "filter_modes.tsv"))
    comparisons = sorted({row["comparison_id"] for row in windows})
    manifest_rows = []

    for comparison_id in comparisons:
        raw_paf = os.path.join(
            args.source_package,
            "raw_paf",
            "%s.sweepga_frequency16_many_many_j0.paf.gz" % comparison_id,
        )
        if not os.path.exists(raw_paf):
            raise SystemExit("missing raw PAF: %s" % raw_paf)
        for mode in sorted(modes, key=lambda row: int(row["plot_order"])):
            output_gz = os.path.join(args.panel_dir, "filtered_paf", "%s.%s.raw_no_chop.paf.gz" % (comparison_id, mode["filter_mode"]))
            log_path = os.path.join(args.panel_dir, "logs", "%s.%s.raw_no_chop_filter.log" % (comparison_id, mode["filter_mode"]))
            command = "gzip -dc %s | %s --num-mappings %s --scaffold-jump %s --scoring %s --output-file <tmp.paf>" % (
                raw_paf,
                args.sweepga,
                mode["num_mappings"],
                mode["scaffold_jump"],
                mode["scoring"],
            )
            status = "reused"
            if args.force or not os.path.exists(output_gz) or os.path.getsize(output_gz) == 0:
                status = "generated"
                run_filter(raw_paf, output_gz, mode, args.sweepga, args.scratch_base, log_path)
            manifest_rows.append({
                "comparison_id": comparison_id,
                "filter_mode": mode["filter_mode"],
                "filter_label": mode["filter_label"],
                "num_mappings": mode["num_mappings"],
                "scaffold_jump": mode["scaffold_jump"],
                "scoring": mode["scoring"],
                "source_raw_paf": raw_paf,
                "source_raw_paf_sha256": read_sidecar_sha256(raw_paf),
                "filtered_paf": output_gz,
                "filtered_paf_sha256": sha256(output_gz),
                "filter_command": command,
                "filter_log": log_path,
                "status": status,
            })

    write_tsv(
        os.path.join(args.panel_dir, "raw_merge_panel_manifest.tsv"),
        manifest_rows,
        [
            "comparison_id", "filter_mode", "filter_label", "num_mappings", "scaffold_jump",
            "scoring", "source_raw_paf", "source_raw_paf_sha256", "filtered_paf",
            "filtered_paf_sha256", "filter_command", "filter_log", "status",
        ],
    )


if __name__ == "__main__":
    main()
