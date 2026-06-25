#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "usage: $0 RAW.paf.gz OUT.paf.gz SUMMARY.tsv COMPARISON_ID [CHOP_LENGTH=10000] [THREADS=8] [CHOP_MODE=row-start|query-grid]" >&2
    exit 2
fi

RAW="$1"
OUT="$2"
SUMMARY="$3"
COMPARISON_ID="$4"
CHOP_LENGTH="${5:-10000}"
THREADS="${6:-${SLURM_CPUS_PER_TASK:-8}}"
CHOP_MODE="${7:-${PAF_CHOP_MODE:-row-start}}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="${PAFCHOP_BIN:-$CRATE_DIR/target/release/pafchop}"

if [[ ! -x "$BIN" ]]; then
    cargo build --release --manifest-path "$CRATE_DIR/Cargo.toml"
fi

mkdir -p "$(dirname "$OUT")" "$(dirname "$SUMMARY")"
TMP_OUT="${OUT}.tmp.$$"
TMP_SUMMARY="${SUMMARY}.tmp.$$"

pigz -dc "$RAW" \
    | "$BIN" --length "$CHOP_LENGTH" --overlap 0 --chunk-mode "$CHOP_MODE" --threads "$THREADS" --comparison-id "$COMPARISON_ID" --summary "$TMP_SUMMARY" \
    | pigz -p "$THREADS" > "$TMP_OUT"

mv "$TMP_OUT" "$OUT"
mv "$TMP_SUMMARY" "$SUMMARY"
gzip -t "$OUT"
sha256sum "$OUT" > "${OUT}.sha256"
