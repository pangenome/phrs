#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$PANEL_DIR/../../.." && pwd)}"
SOURCE_RAW_DIR="${SOURCE_RAW_DIR:-$REPO_ROOT/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/raw_paf}"
COMPARISON_IDS="${COMPARISON_IDS:-PAN027pat_vs_PAN011_joint PAN028mat_vs_PAN027_joint}"
TASK_FILE="$PANEL_DIR/filter_tasks.tsv"

mkdir -p "$PANEL_DIR/filtered_paf" "$PANEL_DIR/logs"
printf "comparison_id\tfilter_mode\tfilter_label\tscaffold_jump\tscoring\tnum_mappings\tsource_paf\toutput_paf\n" > "$TASK_FILE"

tail -n +2 "$PANEL_DIR/config/filter_modes.tsv" |
while IFS=$'\t' read -r filter_mode filter_label scaffold_jump scoring num_mappings; do
    [[ -z "$filter_mode" ]] && continue
    for comparison_id in $COMPARISON_IDS; do
        source_paf="$SOURCE_RAW_DIR/${comparison_id}.sweepga_frequency16_many_many_j0.paf.gz"
        output_paf="$PANEL_DIR/filtered_paf/${comparison_id}.${filter_mode}.paf.gz"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$comparison_id" "$filter_mode" "$filter_label" "$scaffold_jump" "$scoring" "$num_mappings" "$source_paf" "$output_paf" >> "$TASK_FILE"
    done
done

num_tasks="$(($(wc -l < "$TASK_FILE") - 1))"
if [[ "$num_tasks" -le 0 ]]; then
    echo "no raw merge filter tasks generated" >&2
    exit 1
fi

sbatch \
    --array="1-${num_tasks}%3" \
    --output="$PANEL_DIR/logs/fig5_raw_merge.%A_%a.out" \
    --error="$PANEL_DIR/logs/fig5_raw_merge.%A_%a.err" \
    --export=ALL,REPO_ROOT="$REPO_ROOT",PANEL_DIR="$PANEL_DIR",SOURCE_RAW_DIR="$SOURCE_RAW_DIR" \
    "$PANEL_DIR/scripts/run_raw_merge_array.sbatch"
