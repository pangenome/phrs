#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"

RAW_DIR="${WFMASH_RAW_PAF_DIR:-/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10}"
PAFCHOP_DIR="${PAFCHOP_DIR:-$REPO_ROOT/paper_prep/_brainstorming/pafchop-rs}"
PAFCHOP_BIN="${PAFCHOP_BIN:-$REPO_ROOT/target/release/pafchop}"
SWEEPGA="${SWEEPGA:-$(command -v sweepga)}"
PIGZ="${PIGZ:-$(command -v pigz)}"
LENGTHS="${QUERY_GRID_CHOP_LENGTHS:-10000 5000 2000}"
THREADS="${QUERY_GRID_PIGZ_THREADS:-4}"
SCRATCH_BASE="${QUERY_GRID_SCRATCH_BASE:-/dev/shm}"
FILTER_ID="one_to_one_ani_o0"

if [[ ! -d "$RAW_DIR" ]]; then
    echo "missing raw PAF directory: $RAW_DIR" >&2
    exit 1
fi
if [[ ! -x "$SWEEPGA" ]]; then
    echo "missing executable sweepga: $SWEEPGA" >&2
    exit 1
fi
if [[ ! -x "$PIGZ" ]]; then
    echo "missing executable pigz: $PIGZ" >&2
    exit 1
fi
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi
if [[ ! -x "$PAFCHOP_BIN" ]]; then
    cargo build --release --manifest-path "$PAFCHOP_DIR/Cargo.toml"
fi

OUT_ROOT="$PACKAGE_DIR/query_grid_filter"
STATUS_DIR="$PACKAGE_DIR/summaries/query_grid_filter_status"
LOG_DIR="$PACKAGE_DIR/logs/query_grid_filter"
mkdir -p "$OUT_ROOT" "$STATUS_DIR" "$LOG_DIR"

PAFCHOP_SHA="$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')"
SWEEPGA_SHA="$(sha256sum "$SWEEPGA" | awk '{print $1}')"
PIGZ_SHA="$(sha256sum "$PIGZ" | awk '{print $1}')"

raw_for_comparison() {
    local cid="$1"
    local raw="$RAW_DIR/${cid}.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz"
    if [[ ! -s "$raw" ]]; then
        echo "missing raw PAF for $cid: $raw" >&2
        return 1
    fi
    printf "%s\n" "$raw"
}

verify_all_rows_have_cg() {
    local raw="$1"
    "$PIGZ" -dc "$raw" | awk -F'\t' '
        BEGIN { checked = 0; missing = 0 }
        /^$/ || /^#/ { next }
        {
            checked++
            has_cg = 0
            for (i = 13; i <= NF; i++) {
                if ($i ~ /^cg:Z:/) {
                    has_cg = 1
                    break
                }
            }
            if (!has_cg) {
                missing++
                if (missing <= 5) {
                    printf("missing cg:Z in %s line %d\n", FILENAME, NR) > "/dev/stderr"
                }
            }
        }
        END {
            if (checked == 0) {
                print "no PAF rows checked" > "/dev/stderr"
                exit 3
            }
            if (missing > 0) {
                printf("cg:Z verification failed: %d/%d rows missing cg:Z\n", missing, checked) > "/dev/stderr"
                exit 2
            }
            printf("cg:Z verification passed: %d rows checked\n", checked) > "/dev/stderr"
        }'
}

write_status() {
    local status_file="$1"
    shift
    {
        printf "method\tcomparison_id\tchop_length_bp\tfilter_id\traw_paf\tchopped_paf\tfiltered_paf\tchop_summary\tcg_verification\tpafchop_bin\tpafchop_sha256\tsweepga_bin\tsweepga_sha256\tpigz_bin\tpigz_sha256\tchop_command\tfilter_command\tstarted_utc\tfinished_utc\tstatus\n"
        printf "%s\n" "$*"
    } > "$status_file"
}

while IFS=$'\t' read -r cid _rest; do
    [[ "$cid" == "comparison_id" || -z "$cid" ]] && continue
    raw="$(raw_for_comparison "$cid")"
    verify_log="$LOG_DIR/${cid}.cg_verification.log"
    verify_all_rows_have_cg "$raw" 2> "$verify_log"
    cg_verification="$(tr '\n' ';' < "$verify_log" | sed 's/;$//')"

    for length in $LENGTHS; do
        started_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        chopped_dir="$OUT_ROOT/chopped_paf_qgrid_l${length}_o0"
        filtered_dir="$OUT_ROOT/filtered_paf_qgrid_l${length}_o0"
        summary_dir="$PACKAGE_DIR/summaries/pafchop_wfmash_qgrid_l${length}_o0"
        mkdir -p "$chopped_dir" "$filtered_dir" "$summary_dir"

        chopped="$chopped_dir/${cid}.chopped_l${length}_o0_query_grid.paf.gz"
        filtered="$filtered_dir/${cid}.${FILTER_ID}.chopped_l${length}_o0_query_grid.paf.gz"
        summary="$summary_dir/${cid}.summary.tsv"
        status_file="$STATUS_DIR/${cid}.l${length}.tsv"
        tmp_chopped="${chopped}.tmp.$$"
        tmp_summary="${summary}.tmp.$$"
        scratch="$(mktemp -d "$SCRATCH_BASE/wfmash_qgrid_filter.${cid}.${length}.XXXXXX")"
        tmp_paf="$scratch/input.paf"
        tmp_filtered="$scratch/filtered.paf"
        chop_command="pigz -dc $raw | $PAFCHOP_BIN --length $length --overlap 0 --chunk-mode query-grid --comparison-id $cid --summary $summary | pigz -p $THREADS > $chopped"
        filter_command="$SWEEPGA --num-mappings 1:1 --scaffold-jump 0 --scoring ani --overlap 0 --output-file $tmp_filtered $tmp_paf"
        trap 'rm -rf "$scratch" "$tmp_chopped" "$tmp_summary"' EXIT

        "$PIGZ" -dc "$raw" \
            | "$PAFCHOP_BIN" --length "$length" --overlap 0 --chunk-mode query-grid --comparison-id "$cid" --summary "$tmp_summary" \
            | "$PIGZ" -p "$THREADS" > "$tmp_chopped"
        mv "$tmp_chopped" "$chopped"
        mv "$tmp_summary" "$summary"
        "$PIGZ" -t "$chopped"
        sha256sum "$chopped" > "${chopped}.sha256"

        "$PIGZ" -dc "$chopped" > "$tmp_paf"
        "$SWEEPGA" \
            --num-mappings 1:1 \
            --scaffold-jump 0 \
            --scoring ani \
            --overlap 0 \
            --output-file "$tmp_filtered" \
            "$tmp_paf"
        "$PIGZ" -p "$THREADS" < "$tmp_filtered" > "${filtered}.tmp.$$"
        mv "${filtered}.tmp.$$" "$filtered"
        "$PIGZ" -t "$filtered"
        sha256sum "$filtered" > "${filtered}.sha256"
        rm -rf "$scratch"
        trap - EXIT

        finished_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        write_status "$status_file" "wfmash_p95_updated_bin	$cid	$length	$FILTER_ID	$raw	$chopped	$filtered	$summary	$cg_verification	$PAFCHOP_BIN	$PAFCHOP_SHA	$SWEEPGA	$SWEEPGA_SHA	$PIGZ	$PIGZ_SHA	$chop_command	$filter_command	$started_utc	$finished_utc	OK"
    done
done < "$PACKAGE_DIR/config/comparisons.tsv"

python3 "$SCRIPT_DIR/collect_query_grid_filter_outputs.py"
