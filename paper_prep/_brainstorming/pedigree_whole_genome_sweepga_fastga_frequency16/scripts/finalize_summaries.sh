#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

python3 "$SCRIPT_DIR/record_sweepga_binary.py"
python3 "$SCRIPT_DIR/record_fastga_binary.py"
python3 "$SCRIPT_DIR/summarize_paf.py"
python3 "$SCRIPT_DIR/write_raw_chr3_support.py"
python3 "$SCRIPT_DIR/frequency_sensitivity_summary.py"
python3 "$SCRIPT_DIR/write_readme.py"
