#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

SCAN_DIR="paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan"
TMP_DIR="$SCAN_DIR/summaries/tmp_worker_bin_support"
OUT_DIR="$SCAN_DIR/summaries"
mkdir -p "$OUT_DIR"

support_manifest="$OUT_DIR/bin_target_support_manifest.tsv"
{
  printf "method_id\tevidence_layer\tcomparison_id\tbin_target_support_tsv\tbytes\n"
  for f in "$TMP_DIR"/*.bin_target_support.tsv; do
    b=$(stat -c '%s' "$f")
    base=$(basename "$f" .bin_target_support.tsv)
    method=${base%%.*}
    rest=${base#*.}
    evidence=${rest%%.*}
    comparison=${rest#*.}
    printf "%s\t%s\t%s\t%s\t%s\n" "$method" "$evidence" "$comparison" "$f" "$b"
  done
} > "$support_manifest"

awk -F'\t' '
BEGIN { OFS=FS }
function side(arm) {
  if (arm ~ /p$/) return "p";
  if (arm ~ /q$/) return "q";
  return "internal";
}
FNR == 1 { next }
{
  method=$1; evidence=$2; comparison=$3; qchrom=$5; qside=side($6);
  tchrom=$11; tarm=$12; class=$13; paf_rows=$14 + 0; aligned=$15 + 0; mean=$17 + 0;
  key=method FS evidence FS comparison FS qchrom FS qside FS tchrom FS tarm FS class;
  bins[key] += 1;
  rows[key] += paf_rows;
  bp[key] += aligned;
  idbp[key] += aligned * mean;

  region="";
  if (qchrom == "chr9" && qside == "q" && tchrom == "chr3" && tarm ~ /q$/) region="chr9q_to_chr3q";
  else if ((qchrom == "chrX" && tchrom == "chrY") || (qchrom == "chrY" && tchrom == "chrX")) region="PAR_XY";
  else if (qchrom ~ /^chr(13|14|15|21|22)$/ && tchrom ~ /^chr(13|14|15|21|22)$/ && qchrom != tchrom && qside == "p" && tarm ~ /p$/) region="acrocentric_p_cross";
  if (region != "") {
    fkey=region FS method FS evidence FS comparison FS qchrom FS tchrom;
    fbins[fkey] += 1;
    fbp[fkey] += aligned;
    fidbp[fkey] += aligned * mean;
  }
}
END {
  print "method_id\tevidence_layer\tcomparison_id\tquery_chrom\tquery_arm_side\ttarget_chrom\ttarget_arm\tsupport_class\tsupport_bins\tpaf_rows\taligned_bp_sum\tmean_identity_weighted\tmean_match_distance" > totals;
  for (key in bp) {
    mean = (bp[key] > 0 ? idbp[key] / bp[key] : 0);
    print key, bins[key], rows[key], bp[key], sprintf("%.6f", mean), sprintf("%.6f", 1 - mean) >> totals;
  }
  print "region\tmethod_id\tevidence_layer\tcomparison_id\tquery_chrom\ttarget_chrom\tsupport_bins\taligned_bp_sum\tmean_identity_weighted\tmean_match_distance" > focal;
  for (key in fbp) {
    mean = (fbp[key] > 0 ? fidbp[key] / fbp[key] : 0);
    print key, fbins[key], fbp[key], sprintf("%.6f", mean), sprintf("%.6f", 1 - mean) >> focal;
  }
}
' totals="$OUT_DIR/target_support_totals.tsv" focal="$OUT_DIR/focal_region_summary.tsv" "$TMP_DIR"/*.bin_target_support.tsv

cat > "$OUT_DIR/resource_usage.tsv" <<EOF
slurm_job_id	hostname	slurm_cpus_per_task	input_paf_count	process_workers	pigz_threads_per_worker	accounted_helper_threads	wall_seconds	python_executable	method_id	evidence_layer	comparison_id	paf_path	input_lines	bin_target_rows	worker_seconds	pigz_threads	reducer	tmp_output
${SLURM_JOB_ID:-not_slurm}	$(hostname)	${SLURM_CPUS_PER_TASK:-48}	12	12	4	48		/usr/bin/python3									compact_recovery	$TMP_DIR
EOF

cat > "$OUT_DIR/parquet_status.tsv" <<EOF
tsv	parquet	status
target_support_totals.tsv	target_support_totals.parquet	SKIP: pyarrow unavailable in default Slurm Python
focal_region_summary.tsv	focal_region_summary.parquet	SKIP: pyarrow unavailable in default Slurm Python
bin_target_support_manifest.tsv	bin_target_support_manifest.parquet	SKIP: pyarrow unavailable in default Slurm Python
EOF
