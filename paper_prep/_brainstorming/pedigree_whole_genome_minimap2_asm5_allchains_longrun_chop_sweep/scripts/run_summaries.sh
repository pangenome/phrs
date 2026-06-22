#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/refresh_minimap2_jobs.py"
"$SCRIPT_DIR/summarize_minimap2_chop_sweep.py"
