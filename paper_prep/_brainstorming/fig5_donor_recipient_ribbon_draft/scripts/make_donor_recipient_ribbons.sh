#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
python3 paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/scripts/plot_donor_recipient_ribbons.py
