#!/usr/bin/env python3
import os
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


def main():
    sweepga = os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga")
    output = run_text([sweepga, "--check-fastga"])
    rows = [{
        "recorded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sweepga_explicit_path": sweepga,
        "command": sweepga + " --check-fastga",
        "check_fastga_output": output.replace("\t", " ").replace("\n", "\\n"),
    }]
    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "fastga_binary.tsv"),
        rows,
        ["recorded_utc", "sweepga_explicit_path", "command", "check_fastga_output"],
    )


if __name__ == "__main__":
    main()
