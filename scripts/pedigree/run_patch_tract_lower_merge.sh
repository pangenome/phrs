#!/usr/bin/env bash
set -euo pipefail

# Regenerate the lower-merge WashU untangle intermediates and tract summaries.
# Run from the repository root in the moosefs environment:
#
#   bash scripts/pedigree/run_patch_tract_lower_merge.sh
#
# The BED files are large intermediates and are intentionally not committed.

python3 scripts/pedigree/patch_tract_lengths.py --run-lower-untangle
