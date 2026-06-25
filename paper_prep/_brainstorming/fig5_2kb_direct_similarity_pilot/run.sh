#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
python3 scripts/run_direct_similarity_pilot.py
