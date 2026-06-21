#!/usr/bin/env python3
import hashlib
import os
import shutil
import subprocess
from datetime import datetime, timezone

from common import PACKAGE_DIR, write_tsv


def run_text(cmd):
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT).strip()
    except subprocess.CalledProcessError as exc:
        return "ERROR[%s]: %s" % (exc.returncode, exc.output.strip())
    except OSError as exc:
        return "ERROR: %s" % exc


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    explicit = os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga")
    row = {
        "recorded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sweepga_explicit_path": explicit,
        "sweepga_which": shutil.which("sweepga") or "",
        "sweepga_realpath": os.path.realpath(explicit),
        "sweepga_version": run_text([explicit, "--version"]),
        "sweepga_sha256": sha256(explicit),
        "sweepga_help_head": run_text([explicit, "--help"]).replace("\t", " ").replace("\n", "\\n")[:4000],
    }
    os.makedirs(os.path.join(PACKAGE_DIR, "summaries"), exist_ok=True)
    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "sweepga_binary.tsv"),
        [row],
        [
            "recorded_utc",
            "sweepga_explicit_path",
            "sweepga_which",
            "sweepga_realpath",
            "sweepga_version",
            "sweepga_sha256",
            "sweepga_help_head",
        ],
    )


if __name__ == "__main__":
    main()
