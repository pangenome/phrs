#!/usr/bin/env bash
#SBATCH --job-name=fig5_wg_graph_untangle
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH --time=48:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/logs/%x.%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/logs/%x.%j.err

set -euo pipefail

RUN_ID="${1:-wfmash_p95_updated_bin_l2000}"
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
OUT_DIR="${REPO_ROOT}/paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot"
SOURCE_TSV="${OUT_DIR}/config/pilot_sources.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

SEQWISH_BIN="${SEQWISH_BIN:-$(command -v seqwish)}"
ODGI_BIN="${ODGI_BIN:-$(command -v odgi)}"
PIGZ_BIN="${PIGZ_BIN:-$(command -v pigz)}"

mkdir -p "${OUT_DIR}/work/${RUN_ID}" "${OUT_DIR}/summaries" "${OUT_DIR}/logs"
WORK_DIR="${OUT_DIR}/work/${RUN_ID}"
SCRATCH_ROOT="${GRAPH_SCRATCH_ROOT:-${OUT_DIR}/work/tmp}"
mkdir -p "${SCRATCH_ROOT}"
SCRATCH="$(mktemp -d "${SCRATCH_ROOT%/}/fig5_graph_untangle.${RUN_ID}.XXXXXX")"
trap 'rm -rf "${SCRATCH}"' EXIT

row="$(
  awk -F'\t' -v run_id="${RUN_ID}" '
    NR == 1 {
      for (i=1; i<=NF; i++) {
        if ($i == "run_id") {
          run_col=i
        }
      }
      next
    }
    $run_col == run_id {print; found=1}
    END {if (!found) exit 1}
  ' "${SOURCE_TSV}"
)"

header_index() {
  awk -F'\t' -v key="$1" 'NR == 1 {for (i=1; i<=NF; i++) if ($i == key) {print i; exit}}' "${SOURCE_TSV}"
}

get_field() {
  local key="$1"
  local idx
  idx="$(header_index "${key}")"
  printf '%s\n' "${row}" | awk -F'\t' -v idx="${idx}" '{print $idx}'
}

QUERY_FASTA="$(get_field query_fasta)"
TARGET_FASTA="$(get_field target_fasta)"
FILTERED_PAF="$(get_field filtered_paf)"
EXPECTED_PAF_SHA="$(get_field filtered_paf_sha256)"
EXPECTED_ROWS="$(get_field filtered_row_count)"

{
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'seqwish_bin\t%s\n' "${SEQWISH_BIN}"
  "${SEQWISH_BIN}" --version 2>&1 | sed 's/^/seqwish_version\t/'
  printf 'odgi_bin\t%s\n' "${ODGI_BIN}"
  "${ODGI_BIN}" version 2>&1 | sed 's/^/odgi_version\t/' || true
  printf 'pigz_bin\t%s\n' "${PIGZ_BIN}"
  printf 'scratch_root\t%s\n' "${SCRATCH_ROOT}"
  df -h "${SCRATCH_ROOT}" 2>&1 | sed 's/^/scratch_df\t/'
  printf 'query_fasta\t%s\n' "${QUERY_FASTA}"
  printf 'target_fasta\t%s\n' "${TARGET_FASTA}"
  printf 'filtered_paf\t%s\n' "${FILTERED_PAF}"
} > "${OUT_DIR}/summaries/${RUN_ID}.command_manifest.tsv"

cp "${QUERY_FASTA}" "${SCRATCH}/query.fa"
cp "${TARGET_FASTA}" "${SCRATCH}/target.fa"
"${PIGZ_BIN}" -dc "${FILTERED_PAF}" > "${SCRATCH}/filtered.paf"

actual_sha="$(sha256sum "${FILTERED_PAF}" | awk '{print $1}')"
if [[ -n "${EXPECTED_PAF_SHA}" && "${actual_sha}" != "${EXPECTED_PAF_SHA}" ]]; then
  printf 'ERROR: filtered PAF sha256 mismatch: expected %s observed %s\n' "${EXPECTED_PAF_SHA}" "${actual_sha}" >&2
  exit 2
fi

actual_rows="$(wc -l < "${SCRATCH}/filtered.paf" | tr -d ' ')"
if [[ -n "${EXPECTED_ROWS}" && "${actual_rows}" != "${EXPECTED_ROWS}" ]]; then
  printf 'ERROR: filtered PAF row-count mismatch: expected %s observed %s\n' "${EXPECTED_ROWS}" "${actual_rows}" >&2
  exit 2
fi

cat "${SCRATCH}/query.fa" "${SCRATCH}/target.fa" > "${SCRATCH}/seqwish_input.fa"

"${SEQWISH_BIN}" \
  --seqs "${SCRATCH}/seqwish_input.fa" \
  --paf-alns "${SCRATCH}/filtered.paf" \
  --gfa "${SCRATCH}/${RUN_ID}.gfa" \
  --threads "${THREADS}" \
  --temp-dir "${SCRATCH}" \
  --show-progress

"${ODGI_BIN}" build \
  --gfa "${SCRATCH}/${RUN_ID}.gfa" \
  --out "${SCRATCH}/${RUN_ID}.og" \
  --threads "${THREADS}" \
  --progress

"${ODGI_BIN}" paths --idx "${SCRATCH}/${RUN_ID}.og" --list-paths > "${WORK_DIR}/${RUN_ID}.paths.txt"
grep '^PAN027#2#' "${WORK_DIR}/${RUN_ID}.paths.txt" > "${WORK_DIR}/${RUN_ID}.query_paths.txt"
grep '^PAN011#' "${WORK_DIR}/${RUN_ID}.paths.txt" > "${WORK_DIR}/${RUN_ID}.target_paths.txt"

"${ODGI_BIN}" untangle \
  --idx "${SCRATCH}/${RUN_ID}.og" \
  --query-paths "${WORK_DIR}/${RUN_ID}.query_paths.txt" \
  --target-paths "${WORK_DIR}/${RUN_ID}.target_paths.txt" \
  --merge-dist 50000 \
  --n-best 4 \
  --threads "${THREADS}" \
  --progress \
  > "${WORK_DIR}/${RUN_ID}.untangle.bedpe"

"${ODGI_BIN}" untangle \
  --idx "${SCRATCH}/${RUN_ID}.og" \
  --query-paths "${WORK_DIR}/${RUN_ID}.query_paths.txt" \
  --target-paths "${WORK_DIR}/${RUN_ID}.target_paths.txt" \
  --merge-dist 50000 \
  --n-best 4 \
  --paf-output \
  --threads "${THREADS}" \
  --progress \
  > "${WORK_DIR}/${RUN_ID}.untangle.paf"

cp "${SCRATCH}/${RUN_ID}.gfa" "${WORK_DIR}/${RUN_ID}.gfa"
cp "${SCRATCH}/${RUN_ID}.og" "${WORK_DIR}/${RUN_ID}.og"

"${OUT_DIR}/scripts/summarize_untangle_support.py" \
  --run-id "${RUN_ID}" \
  --paths "${WORK_DIR}/${RUN_ID}.paths.txt" \
  --bedpe "${WORK_DIR}/${RUN_ID}.untangle.bedpe" \
  --focus "${OUT_DIR}/config/focus_patterns.tsv" \
  --summary "${OUT_DIR}/summaries/${RUN_ID}.focus_summary.tsv" \
  --segments "${OUT_DIR}/summaries/${RUN_ID}.focus_segments.tsv"

{
  printf 'finished_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'gfa\t%s\n' "${WORK_DIR}/${RUN_ID}.gfa"
  printf 'og\t%s\n' "${WORK_DIR}/${RUN_ID}.og"
  printf 'paths\t%s\n' "${WORK_DIR}/${RUN_ID}.paths.txt"
  printf 'bedpe\t%s\n' "${WORK_DIR}/${RUN_ID}.untangle.bedpe"
  printf 'paf\t%s\n' "${WORK_DIR}/${RUN_ID}.untangle.paf"
  printf 'focus_summary\t%s\n' "${OUT_DIR}/summaries/${RUN_ID}.focus_summary.tsv"
  printf 'focus_segments\t%s\n' "${OUT_DIR}/summaries/${RUN_ID}.focus_segments.tsv"
} >> "${OUT_DIR}/summaries/${RUN_ID}.command_manifest.tsv"
