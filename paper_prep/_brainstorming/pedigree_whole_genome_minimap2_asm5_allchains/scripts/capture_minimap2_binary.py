#!/usr/bin/env python3
import os
import shutil
import subprocess

from common import PACKAGE_DIR, write_tsv


MINIMAP2_BIN = os.environ.get("MINIMAP2_BIN", "/home/erikg/bin/minimap2")
EXPECTED_REALPATH = "/export/local/home/erikg/bin/minimap2-v2.31-r1302"
EXPECTED_VERSION = "2.31-r1302"
EXPECTED_SHA256 = "5a0e9d6b351f1aa5d11a5067bd29a33bc50abe70c51fc9be9e1899ec1643c949"
SOURCE_CHECKOUT = "/home/erikg/minimap2"
SOURCE_TAG = "v2.31"
SOURCE_COMMIT = "3c28777e7e2dcc90f825de1b9f17a89cca7d4452"


def run(args):
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT).strip()


def main():
    os.makedirs(os.path.join(PACKAGE_DIR, "logs"), exist_ok=True)
    explicit_realpath = os.path.realpath(MINIMAP2_BIN)
    version = run([MINIMAP2_BIN, "--version"]).splitlines()[0]
    sha256 = run(["sha256sum", MINIMAP2_BIN]).split()[0]
    help_text = run([MINIMAP2_BIN, "--help"])
    relevant_help = []
    for line in help_text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("-x", "-c", "--eqx", "-P", "--q-occ-frac", "-t")) or "asm5" in stripped:
            relevant_help.append(line)
        if len(relevant_help) == 24:
            break
    help_path = os.path.join(PACKAGE_DIR, "logs", "minimap2_v2.31-r1302.help.txt")
    with open(help_path, "w") as out:
        out.write(help_text)
        if not help_text.endswith("\n"):
            out.write("\n")

    status = "PASS"
    notes = []
    if explicit_realpath != EXPECTED_REALPATH:
        status = "MISMATCH"
        notes.append("realpath expected %s" % EXPECTED_REALPATH)
    if version != EXPECTED_VERSION:
        status = "MISMATCH"
        notes.append("version expected %s" % EXPECTED_VERSION)
    if sha256 != EXPECTED_SHA256:
        status = "MISMATCH"
        notes.append("sha256 expected %s" % EXPECTED_SHA256)

    rows = [{
        "role": "required_source_built_local_binary",
        "which_minimap2": shutil.which("minimap2") or "",
        "explicit_path": MINIMAP2_BIN,
        "realpath": explicit_realpath,
        "version": version,
        "sha256": sha256,
        "expected_realpath": EXPECTED_REALPATH,
        "expected_version": EXPECTED_VERSION,
        "expected_sha256": EXPECTED_SHA256,
        "source_checkout": SOURCE_CHECKOUT,
        "source_tag": SOURCE_TAG,
        "source_commit": SOURCE_COMMIT,
        "help_relevant_lines": " | ".join(relevant_help),
        "help_log": help_path,
        "status": status,
        "note": "; ".join(notes) if notes else "explicit /home/erikg/bin/minimap2 matches task-created provenance",
    }]
    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "minimap2_binary.tsv"),
        rows,
        [
            "role",
            "which_minimap2",
            "explicit_path",
            "realpath",
            "version",
            "sha256",
            "expected_realpath",
            "expected_version",
            "expected_sha256",
            "source_checkout",
            "source_tag",
            "source_commit",
            "help_relevant_lines",
            "help_log",
            "status",
            "note",
        ],
    )
    if status != "PASS":
        raise SystemExit("minimap2 binary provenance mismatch")


if __name__ == "__main__":
    main()
