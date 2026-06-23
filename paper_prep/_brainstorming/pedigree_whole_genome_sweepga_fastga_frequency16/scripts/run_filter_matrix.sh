#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
COMPARISONS="$PACKAGE_DIR/config/comparisons.tsv"
FILTERS="$PACKAGE_DIR/config/filter_matrix.tsv"
STATUS="$PACKAGE_DIR/summaries/filter_manifest.tsv"
CHOP_MODE="${PAF_CHOP_MODE:-query-grid}"
case "$CHOP_MODE" in
    row-start|query-grid) ;;
    *) echo "invalid PAF_CHOP_MODE: $CHOP_MODE" >&2; exit 2 ;;
esac
MODE_SUFFIX=""
if [[ "$CHOP_MODE" == "query-grid" ]]; then
    MODE_SUFFIX="_query_grid"
fi
OUTPUT_DIR="filtered_paf${MODE_SUFFIX}"
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "required sweepGA filter scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

mkdir -p "$PACKAGE_DIR/$OUTPUT_DIR" "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"
printf "comparison_id\tfilter_id\tinput_paf\toutput_paf\tnum_mappings\tscaffold_jump\tchunk_mode\tsweepga_devshm_base\tstatus\n" > "$STATUS"
input_dir="chopped_paf_l${PAF_CHOP_LENGTH:-10000}_o${PAF_CHOP_OVERLAP:-0}${MODE_SUFFIX}"
if [[ ! -d "$PACKAGE_DIR/$input_dir" ]]; then
    input_dir="chopped_paf"
fi

while IFS=$'\t' read -r cid _rest; do
    [[ "$cid" == "comparison_id" || -z "$cid" ]] && continue
    while IFS=$'\t' read -r filter_id num_mappings scaffold_jump _source_dir _note; do
        [[ "$filter_id" == "filter_id" || -z "$filter_id" ]] && continue
        input="$PACKAGE_DIR/$input_dir/${cid}.chopped_l${PAF_CHOP_LENGTH:-10000}_o${PAF_CHOP_OVERLAP:-0}${MODE_SUFFIX}.paf.gz"
        output="$PACKAGE_DIR/$OUTPUT_DIR/${cid}.${filter_id}.paf.gz"
        echo "filtering comparison=$cid filter=$filter_id input=$input output=$output"
        echo "sweepga_filter_scratch_base=$SCRATCH_BASE"
        SWEEPGA_DEVSHM_BASE="$SCRATCH_BASE" python3 "$SCRIPT_DIR/filter_paf.py" \
            --comparison-id "$cid" \
            --filter-id "$filter_id" \
            --input-dir "$input_dir" \
            --output-dir "$OUTPUT_DIR" \
            --chop-length "${PAF_CHOP_LENGTH:-10000}" \
            --overlap "${PAF_CHOP_OVERLAP:-0}" \
            --chunk-mode "$CHOP_MODE"
        gzip -t "$output"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n" "$cid" "$filter_id" "$input" "$output" "$num_mappings" "$scaffold_jump" "$CHOP_MODE" "$SCRATCH_BASE" >> "$STATUS"
    done < "$FILTERS"
done < "$COMPARISONS"
