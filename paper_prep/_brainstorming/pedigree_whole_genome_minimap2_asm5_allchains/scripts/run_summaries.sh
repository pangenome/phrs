#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/refresh_minimap2_jobs.py"
if compgen -G "$SCRIPT_DIR/../raw_paf/*/*.paf.gz" > /dev/null; then
    "$SCRIPT_DIR/summarize_paf.py"
    "$SCRIPT_DIR/summarize_candidate_windows.py"
    "$SCRIPT_DIR/summarize_chr3_support.py"
else
    "$SCRIPT_DIR/write_pathological_runtime_summaries.py"
fi
