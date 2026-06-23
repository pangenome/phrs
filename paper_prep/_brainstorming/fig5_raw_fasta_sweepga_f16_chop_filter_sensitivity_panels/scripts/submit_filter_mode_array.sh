#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$PANEL_DIR/../../.." && pwd)}"
SOURCE_PACKAGE="${SOURCE_PACKAGE:-$REPO_ROOT/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16}"
CHOP_LENGTHS="${CHOP_LENGTHS:-10000 5000 2000}"
COMPARISON_IDS="${COMPARISON_IDS:-PAN027pat_vs_PAN011_joint PAN028mat_vs_PAN027_joint}"
TASK_FILE="$PANEL_DIR/filter_tasks.tsv"

mkdir -p "$PANEL_DIR/filtered_paf" "$PANEL_DIR/logs"
printf "comparison_id\tchop_length_bp\tfilter_mode\tfilter_label\tscaffold_jump\tscoring\tnum_mappings\tsource_paf\toutput_paf\n" > "$TASK_FILE"

tail -n +2 "$PANEL_DIR/config/filter_modes.tsv" |
while IFS=$'\t' read -r filter_mode filter_label scaffold_jump scoring num_mappings source; do
    [[ -z "$filter_mode" ]] && continue
    [[ "$source" == "existing" ]] && continue
    for comparison_id in $COMPARISON_IDS; do
        for length in $CHOP_LENGTHS; do
            source_paf="$SOURCE_PACKAGE/chopped_paf_l${length}_o0/${comparison_id}.chopped_l${length}_o0.paf.gz"
            output_paf="$PANEL_DIR/filtered_paf/${comparison_id}.l${length}.${filter_mode}.paf.gz"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$comparison_id" "$length" "$filter_mode" "$filter_label" "$scaffold_jump" "$scoring" "$num_mappings" "$source_paf" "$output_paf" >> "$TASK_FILE"
        done
    done
done

num_tasks="$(($(wc -l < "$TASK_FILE") - 1))"
if [[ "$num_tasks" -le 0 ]]; then
    echo "no filter tasks generated" >&2
    exit 1
fi

sbatch \
    --array="1-${num_tasks}%6" \
    --output="$PANEL_DIR/logs/fig5_filter_modes.%A_%a.out" \
    --error="$PANEL_DIR/logs/fig5_filter_modes.%A_%a.err" \
    --export=ALL,REPO_ROOT="$REPO_ROOT",PANEL_DIR="$PANEL_DIR",SOURCE_PACKAGE="$SOURCE_PACKAGE" \
    "$PANEL_DIR/scripts/run_filter_mode_array.sbatch"
