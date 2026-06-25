#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

python3 scripts/build_length_scaled_tracks.py
Rscript scripts/plot_length_scaled_tracks.R "$here"
bash validate_outputs.sh
