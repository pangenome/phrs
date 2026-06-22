#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
FIXTURE_DIR="$PACKAGE_DIR/synthetic_paf_filter_audit"
OUT_DIR="$FIXTURE_DIR/outputs"
SUMMARY="$PACKAGE_DIR/summaries/sweepga_paf_filter_identity_audit.tsv"
LOG="$PACKAGE_DIR/logs/sweepga_paf_filter_identity_audit.log"

mkdir -p "$OUT_DIR" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/logs"

run_case() {
    local test_id="$1"
    local fixture="$2"
    local expected="$3"
    local interpretation="$4"
    shift 4
    local output="$OUT_DIR/${test_id}.paf"
    local stderr="$OUT_DIR/${test_id}.stderr"
    local cmd=("$SWEEPGA" "$FIXTURE_DIR/$fixture" --output-file "$output" "$@")
    {
        printf '## %s\n' "$test_id"
        printf 'command='
        printf '%q ' "${cmd[@]}"
        printf '\n'
    } >> "$LOG"
    "${cmd[@]}" > "$OUT_DIR/${test_id}.stdout" 2> "$stderr"
    local status=$?
    local observed
    if [[ -s "$output" ]]; then
        observed="$(awk 'BEGIN{ORS=";"} {print $1 ":" $3 "-" $4 ">" $6 ":" $8 "-" $9 ":" $10 "/" $11}' "$output")"
        observed="${observed%;}"
    else
        observed="NO_OUTPUT"
    fi
    local pass="FAIL"
    if [[ "$status" -eq 0 && "$observed" == "$expected" ]]; then
        pass="PASS"
    fi
    local command_text
    command_text="$(printf '%q ' "${cmd[@]}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$test_id" "$fixture" "$command_text" "$expected" "$observed" "$pass" "$interpretation" >> "$SUMMARY"
}

: > "$LOG"
printf 'test_id\tfixture\tcommand\texpected_winner\tobserved_winner\tpass_fail\tinterpretation\n' > "$SUMMARY"

run_case \
    unequal_1to1_default_log_length_ani \
    unequal_length_competition.paf \
    'qA:0-1000>tA:0-1000:800/1000' \
    'Default log-length-ani should prefer the longer lower-identity chunk over the shorter high-identity chunk.' \
    --num-mappings 1:1 --scaffold-jump 0

run_case \
    unequal_1to1_scoring_ani \
    unequal_length_competition.paf \
    'qA:0-100>tA:0-100:99/100' \
    '--scoring ani should rank by per-row identity rather than raw length or length*identity.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    unequal_manymany_default \
    unequal_length_competition.paf \
    'qA:0-1000>tA:0-1000:800/1000;qA:0-100>tA:0-100:99/100' \
    'many:many with scaffolding disabled should preserve both competing rows.' \
    --num-mappings many:many --scaffold-jump 0

run_case \
    unequal_manymany_scoring_ani \
    unequal_length_competition.paf \
    'qA:0-1000>tA:0-1000:800/1000;qA:0-100>tA:0-100:99/100' \
    '--scoring ani does not matter when many:many has no mapping-axis limit and scaffolding is disabled.' \
    --num-mappings many:many --scaffold-jump 0 --scoring ani

run_case \
    equal_1to1_default \
    equal_length_competition.paf \
    'qB:0-100>tB:0-100:90/100' \
    'Equal-length competitors distinguish identity-sensitive scoring without a length confound.' \
    --num-mappings 1:1 --scaffold-jump 0

run_case \
    equal_1to1_scoring_ani \
    equal_length_competition.paf \
    'qB:0-100>tB:0-100:90/100' \
    'ANI scoring should keep the higher identity equal-length row.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    columns_vs_cg_1to1_ani \
    columns_vs_cg_precedence.paf \
    'qC:0-100>tC:0-100:80/100' \
    'cg:Z supplies identity for scoring when present; the output row still preserves the original PAF columns.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    dv_overrides_columns_1to1_ani \
    dv_overrides_columns.paf \
    'qD:0-100>tD:0-100:80/100' \
    'If this passes, dv:f overrides PAF columns 10/11 for scoring identity.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    non_overlapping_1to1_ani \
    non_overlapping_many_many.paf \
    'qE:0-100>tE:0-100:90/100;qE:200-300>tE:200-300:80/100' \
    '1:1 is plane-sweep overlap filtering on query and target axes; non-overlapping rows can both pass.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    non_overlapping_manymany_ani \
    non_overlapping_many_many.paf \
    'qE:0-100>tE:0-100:90/100;qE:200-300>tE:200-300:80/100' \
    'many:many preserves non-overlapping rows when scaffolding is disabled.' \
    --num-mappings many:many --scaffold-jump 0 --scoring ani

run_case \
    scaffold_disabled_row_probe \
    scaffold_chain_probe.paf \
    'qF:0-100>tF:0-100:100/100;qF:120-220>tF:120-220:100/100' \
    'With --scaffold-jump 0, no scaffold-chain merge is applied; row-level plane sweep keeps the two perfect non-overlapping rows under ANI scoring.' \
    --num-mappings 1:1 --scaffold-jump 0 --scoring ani

run_case \
    scaffold_enabled_chain_probe \
    scaffold_chain_probe.paf \
    'qF:0-100>tF:0-100:100/100;qF:120-220>tF:120-220:100/100' \
    'With explicit scaffolding enabled, nearby rows can be handled as a chain and the member rows are emitted.' \
    --num-mappings 1:1 --scaffold-jump 50 --scaffold-mass 100 --scaffold-overlap 0.5 --scoring ani
