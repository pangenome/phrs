#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$PANEL_DIR/../../.." && pwd)"

mkdir -p "$PANEL_DIR/logs"
sbatch \
    --partition=tux \
    --nodelist=tux05 \
    --cpus-per-task=96 \
    --mem=700G \
    --job-name=fig5_f16_panel_96 \
    --output="$PANEL_DIR/logs/fig5_f16_panel_96.%j.out" \
    --error="$PANEL_DIR/logs/fig5_f16_panel_96.%j.err" \
    "$PANEL_DIR/scripts/run_raw_fasta_chopped_panel.sbatch" \
    "$REPO_ROOT" \
    "$PANEL_DIR"
