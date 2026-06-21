#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/refresh_wfmash_jobs.py"
"$SCRIPT_DIR/summarize_paf.py"
"$SCRIPT_DIR/summarize_candidate_windows.py"
