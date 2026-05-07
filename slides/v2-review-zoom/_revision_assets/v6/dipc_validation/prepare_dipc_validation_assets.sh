#!/usr/bin/env bash
set -euo pipefail

asset_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$asset_dir"

mkdir -p source_pdfs pdf_pngs plots

manifest="source_manifest.tsv"
conversion_log="conversion_log.tsv"

printf "asset\ttype\tsource_path\tstaged_path\tsha256\tbytes\n" > "$manifest"
printf "asset\tsource_pdf\tconverted_png\tppi\tbytes\n" > "$conversion_log"

stage_pdf() {
  local asset="$1"
  local source_pdf="$2"
  local staged_pdf="source_pdfs/${asset}.pdf"
  local png_prefix="pdf_pngs/${asset}"
  local png_path="${png_prefix}.png"

  if [[ ! -s "$source_pdf" ]]; then
    echo "Missing source PDF: $source_pdf" >&2
    exit 1
  fi

  cp -f "$source_pdf" "$staged_pdf"
  guix shell ghostscript poppler -- \
    pdftoppm -r 220 -png -singlefile "$staged_pdf" "$png_prefix"

  local pdf_sha png_sha pdf_bytes png_bytes
  pdf_sha="$(sha256sum "$staged_pdf" | awk '{print $1}')"
  png_sha="$(sha256sum "$png_path" | awk '{print $1}')"
  pdf_bytes="$(stat -c '%s' "$staged_pdf")"
  png_bytes="$(stat -c '%s' "$png_path")"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$asset" "source_pdf" "$source_pdf" "$staged_pdf" "$pdf_sha" "$pdf_bytes" >> "$manifest"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$asset" "converted_png" "$staged_pdf" "$png_path" "$png_sha" "$png_bytes" >> "$manifest"
  printf "%s\t%s\t%s\t%s\t%s\n" "$asset" "$staged_pdf" "$png_path" "220" "$png_bytes" >> "$conversion_log"
}

stage_pdf \
  "gm12878_mantel_scatter" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_scatter.pdf"
stage_pdf \
  "gm12878_radial_community" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf"
stage_pdf \
  "sperm_all20_mantel_scatter" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf"
stage_pdf \
  "sperm_all20_radial_community" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf"
stage_pdf \
  "sperm_all20_by_arm_type_arm" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.arm.pdf"
stage_pdf \
  "sperm_all20_by_arm_type_per_cell" \
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.per-cell.pdf"

Rscript make_dipc_validation_summary_plots.R

find source_pdfs pdf_pngs plots -maxdepth 1 -type f | sort
