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
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "required sweepGA filter scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

mkdir -p "$PACKAGE_DIR/filtered_paf" "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"
printf "comparison_id\tfilter_id\tinput_paf\toutput_paf\tnum_mappings\tscaffold_jump\tsweepga_devshm_base\tstatus\n" > "$STATUS"
input_dir="chopped_paf_l${PAF_CHOP_LENGTH:-10000}_o${PAF_CHOP_OVERLAP:-0}"
if [[ ! -d "$PACKAGE_DIR/$input_dir" ]]; then
    input_dir="chopped_paf"
fi

while IFS=$'\t' read -r cid _rest; do
    [[ "$cid" == "comparison_id" || -z "$cid" ]] && continue
    while IFS=$'\t' read -r filter_id num_mappings scaffold_jump _source_dir _note; do
        [[ "$filter_id" == "filter_id" || -z "$filter_id" ]] && continue
        input="$PACKAGE_DIR/$input_dir/${cid}.chopped_l${PAF_CHOP_LENGTH:-10000}_o${PAF_CHOP_OVERLAP:-0}.paf.gz"
        output="$PACKAGE_DIR/filtered_paf/${cid}.${filter_id}.paf.gz"
        echo "filtering comparison=$cid filter=$filter_id input=$input output=$output"
        echo "sweepga_filter_scratch_base=$SCRATCH_BASE"
        SWEEPGA_DEVSHM_BASE="$SCRATCH_BASE" python3 "$SCRIPT_DIR/filter_paf.py" \
            --comparison-id "$cid" \
            --filter-id "$filter_id" \
            --input-dir "$input_dir" \
            --chop-length "${PAF_CHOP_LENGTH:-10000}" \
            --overlap "${PAF_CHOP_OVERLAP:-0}"
        gzip -t "$output"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n" "$cid" "$filter_id" "$input" "$output" "$num_mappings" "$scaffold_jump" "$SCRATCH_BASE" >> "$STATUS"
    done < "$FILTERS"
done < "$COMPARISONS"
