#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

mkdir -p "$PACKAGE_DIR/inputs" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/logs"
"$SCRIPT_DIR/prepare_inputs.py" 2>&1 | tee "$PACKAGE_DIR/logs/prepare_inputs.$(date -u +%Y%m%dT%H%M%SZ).log"
