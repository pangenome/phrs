#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

tail -n +2 "$PACKAGE_DIR/config/comparisons.tsv" | cut -f1 | while read -r cid; do
    [[ -n "$cid" ]] || continue
    raw="$PACKAGE_DIR/raw_paf/${cid}.sweepga_many_many_j0.paf.gz"
    many="$PACKAGE_DIR/filtered_paf/${cid}.many_many_noscaffold.paf.gz"
    if [[ -s "$raw" ]]; then
        cp "$raw" "$many"
    fi
    tail -n +2 "$PACKAGE_DIR/config/filter_matrix.tsv" | cut -f1 | while read -r filter_id; do
        [[ "$filter_id" == "many_many_noscaffold" ]] && continue
        python3 "$SCRIPT_DIR/filter_paf.py" --comparison-id "$cid" --filter-id "$filter_id"
    done
done
