#!/usr/bin/env bash
#SBATCH --cpus-per-task=96
#SBATCH --mem=700G
#SBATCH --time=12:00:00

set -euo pipefail

REPO_ROOT="${1:?repo root required}"
PANEL_DIR="${2:?panel dir required}"
SOURCE_PACKAGE="${SOURCE_PACKAGE:-/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16}"
PAFCHOP_BIN="${PAFCHOP_BIN:-/moosefs/erikg/phrs/paper_prep/_brainstorming/pafchop-rs/target/release/pafchop}"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
CHOP_LENGTH="${CHOP_LENGTH:-2000}"
OVERLAP="${OVERLAP:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-96}"
JOB_ID="${SLURM_JOB_ID:-manual}"
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"

if [[ ! -d "$SOURCE_PACKAGE/raw_paf" ]]; then
    echo "missing raw_paf directory in SOURCE_PACKAGE=$SOURCE_PACKAGE" >&2
    exit 1
fi
if [[ ! -x "$PAFCHOP_BIN" ]]; then
    cargo build --release --manifest-path /moosefs/erikg/phrs/paper_prep/_brainstorming/pafchop-rs/Cargo.toml
fi
if [[ ! -x "$PAFCHOP_BIN" ]]; then
    echo "missing pafchop binary: $PAFCHOP_BIN" >&2
    exit 1
fi
if [[ ! -x "$SWEEPGA" ]]; then
    echo "missing sweepga binary: $SWEEPGA" >&2
    exit 1
fi
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

mkdir -p "$PANEL_DIR/work/raw_window_paf" \
         "$PANEL_DIR/work/chopped_paf_l${CHOP_LENGTH}_o${OVERLAP}" \
         "$PANEL_DIR/evidence_paf" \
         "$PANEL_DIR/logs"

exec 9>"$PANEL_DIR/.raw_fasta_chopped_panel.lock"
if ! flock -n 9; then
    echo "another raw FASTA chopped panel job is already running for $PANEL_DIR" >&2
    exit 1
fi

PIPELINE_LOG="$PANEL_DIR/logs/raw_fasta_chopped_panel_pipeline.${JOB_ID}.log"
{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "job_id=$JOB_ID"
    echo "host=$(hostname)"
    echo "repo_root=$REPO_ROOT"
    echo "panel_dir=$PANEL_DIR"
    echo "source_package=$SOURCE_PACKAGE"
    echo "pafchop_bin=$PAFCHOP_BIN"
    echo "pafchop_sha256=$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')"
    echo "sweepga=$SWEEPGA"
    echo "sweepga_sha256=$(sha256sum "$SWEEPGA" | awk '{print $1}')"
    echo "chop_length=$CHOP_LENGTH"
    echo "overlap=$OVERLAP"
    echo "threads=$THREADS"
    echo "scratch_base=$SCRATCH_BASE"
} > "$PIPELINE_LOG"

if [[ -n "${ONLY_COMPARISON:-}" ]]; then
    printf "%s\n" "$ONLY_COMPARISON"
else
    tail -n +2 "$PANEL_DIR/config/panel_windows.tsv" | cut -f2 | sort -u
fi |
while read -r comparison_id; do
    raw="$SOURCE_PACKAGE/raw_paf/${comparison_id}.sweepga_frequency16_many_many_j0.paf.gz"
    raw_window="$PANEL_DIR/work/raw_window_paf/${comparison_id}.raw_window_extract.paf.gz"
    chopped="$PANEL_DIR/work/chopped_paf_l${CHOP_LENGTH}_o${OVERLAP}/${comparison_id}.chopped_l${CHOP_LENGTH}_o${OVERLAP}.paf.gz"
    filtered="$PANEL_DIR/evidence_paf/${comparison_id}.window_extract.one_one_scoring_ani_chopped_l${CHOP_LENGTH}_o${OVERLAP}.paf.gz"
    work="$(mktemp -d "$SCRATCH_BASE/fig5_f16_panel.${comparison_id}.${JOB_ID}.XXXXXX")"
    trap 'rm -rf "$work"' EXIT

    {
        echo "comparison_id=$comparison_id"
        echo "raw_paf=$raw"
        echo "raw_window_extract_paf=$raw_window"
        echo "chopped_paf=$chopped"
        echo "filtered_paf=$filtered"
        date -u '+stage_started_utc=%Y-%m-%dT%H:%M:%SZ'
    } >> "$PIPELINE_LOG"

    pigz -t -p "$THREADS" "$raw"
    if [[ -s "$raw_window" ]] && ! pigz -t -p "$THREADS" "$raw_window"; then
        rm -f "$raw_window"
    fi
    if [[ ! -s "$raw_window" ]]; then
        tmp_raw_window="${raw_window}.tmp.${JOB_ID}"
        rm -f "$tmp_raw_window"
        window_spec="$(
            awk -F '\t' -v cid="$comparison_id" '
                NR > 1 && $2 == cid {
                    printf "%s:%s:%s,", $3, $4, $5
                }
            ' "$PANEL_DIR/config/panel_windows.tsv"
        )"
        pigz -dc -p "$THREADS" "$raw" |
            awk -F '\t' -v spec="$window_spec" '
                BEGIN {
                    n = split(spec, a, ",")
                    for (i = 1; i <= n; i++) {
                        if (a[i] == "") continue
                        split(a[i], b, ":")
                        qs[b[1]] = b[2] + 0
                        qe[b[1]] = b[3] + 0
                    }
                }
                ($1 in qs) && (($4 + 0) > qs[$1]) && (($3 + 0) < qe[$1]) { print }
            ' |
            pigz -p "$THREADS" > "$tmp_raw_window"
        pigz -t -p "$THREADS" "$tmp_raw_window"
        mv "$tmp_raw_window" "$raw_window"
    fi
    pigz -t -p "$THREADS" "$raw_window"

    if [[ -s "$chopped" ]] && ! pigz -t -p "$THREADS" "$chopped"; then
        rm -f "$chopped"
    fi
    if [[ ! -s "$chopped" ]]; then
        tmp_chopped="${chopped}.tmp.${JOB_ID}"
        rm -f "$tmp_chopped"
        pigz -dc -p "$THREADS" "$raw_window" |
            "$PAFCHOP_BIN" \
                --input - \
                --length "$CHOP_LENGTH" \
                --overlap "$OVERLAP" \
                --threads "$THREADS" \
                --comparison-id "$comparison_id" \
                --summary "$PANEL_DIR/work/${comparison_id}.chopped_l${CHOP_LENGTH}_o${OVERLAP}.summary.tsv" |
            pigz -p "$THREADS" > "$tmp_chopped"
        pigz -t -p "$THREADS" "$tmp_chopped"
        mv "$tmp_chopped" "$chopped"
    fi
    pigz -t -p "$THREADS" "$chopped"

    if [[ -s "$filtered" ]] && ! pigz -t -p "$THREADS" "$filtered"; then
        rm -f "$filtered"
    fi
    if [[ ! -s "$filtered" ]]; then
        pigz -dc -p "$THREADS" "$chopped" > "$work/input.chopped.paf"
        "$SWEEPGA" \
            --num-mappings 1:1 \
            --scaffold-jump 0 \
            --scoring ani \
            --output-file "$work/filtered.paf" \
            "$work/input.chopped.paf"
        tmp_filtered="${filtered}.tmp.${JOB_ID}"
        rm -f "$tmp_filtered"
        pigz -p "$THREADS" -c "$work/filtered.paf" > "$tmp_filtered"
        pigz -t -p "$THREADS" "$tmp_filtered"
        mv "$tmp_filtered" "$filtered"
    fi
    pigz -t -p "$THREADS" "$filtered"
    sha256sum "$raw" "$raw_window" "$chopped" "$filtered" >> "$PIPELINE_LOG"
    rm -rf "$work"
    trap - EXIT
done

python3 "$PANEL_DIR/scripts/extract_raw_fasta_chopped_segments.py" \
    --panel-dir "$PANEL_DIR" \
    --source-package "$SOURCE_PACKAGE" \
    --job-id "$JOB_ID" \
    --chop-length "$CHOP_LENGTH" \
    --overlap "$OVERLAP"

Rscript "$PANEL_DIR/scripts/plot_raw_fasta_chopped_panel.R" "$PANEL_DIR"

{
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "segments=$PANEL_DIR/raw_fasta_chopped_panel_segments.tsv"
    echo "provenance=$PANEL_DIR/slurm_jobs.tsv"
    echo "pdf=$PANEL_DIR/fig5_raw_fasta_sweepga_f16_chopped_panels.pdf"
    echo "svg=$PANEL_DIR/fig5_raw_fasta_sweepga_f16_chopped_panels.svg"
} >> "$PIPELINE_LOG"
