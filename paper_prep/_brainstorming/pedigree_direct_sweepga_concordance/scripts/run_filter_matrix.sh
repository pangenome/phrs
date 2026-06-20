#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

tail -n +2 "$PACKAGE_DIR/config/comparisons.tsv" | cut -f1 | while read -r cid; do
    [[ -n "$cid" ]] || continue
    tail -n +2 "$PACKAGE_DIR/config/filter_matrix.tsv" | cut -f1 | while read -r filter_id; do
        [[ "$filter_id" == "many_many_noscaffold" ]] && continue
        python3 "$SCRIPT_DIR/filter_paf.py" --comparison-id "$cid" --filter-id "$filter_id"
    done
done
