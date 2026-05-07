#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

layout_tsv="${LAYOUT_TSV:-/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv}"
component="${COMPONENT:-8}"
out_png="${OUT_PNG:-pggb_graph_2d.png}"
render_log="${RENDER_LOG:-render_log.tsv}"

Rscript "./render_pggb_layout_component8.R" \
  --layout-tsv "${layout_tsv}" \
  --component "${component}" \
  --out "${out_png}" \
  --render-log "${render_log}"
