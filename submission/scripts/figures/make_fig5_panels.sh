#!/usr/bin/env bash
# Generate the three Fig5 ribbon panels (B/C/D) as separate cropped PDF+PNG,
# with NO baked-in panel letter and NO per-panel legend (letters + one shared
# legend are added in paper.tex). Panel A (pedigree) is made by make_fig5_pedigree.R.
# Needs: python3 (stdlib only) + inkscape. Reads vendored data/ inputs.
set -euxo pipefail
cd "$(git rev-parse --show-toplevel)"

PLOTTER=submission/scripts/figures/fig5_ribbon_plotter.py
OUT=submission/fig/MainFigures
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# crop window (SVG user units): keep content y in [CROP_Y0, CROP_Y0+CROP_H]
# and trim empty left margin (content/row-labels start at x~186) to give the
# chromosomes more horizontal space when scaled to \linewidth.
CROP_Y0=75; CROP_H=565; CROP_X0=165; CROP_W=$((3600 - CROP_X0))

gen () {  # id  query_label  h1_label  h2_label  out_base
  local id="$1" q="$2" h1="$3" h2="$4" base="$5"
  echo ">>> panel $base ($id)"
  FIG5_HIDE_LEGEND=1 python3 "$PLOTTER" \
    --comparison-id "$id" \
    --class-winners "data/fig5_${id}.class_winners.impg_similarity.tsv.gz" \
    --query-fai "data/fig5_${id}.query.fa.fai" \
    --target-fai "data/fig5_${id}.target.fa.fai" \
    --output-dir "$TMP" \
    --query-label "$q" --target-h1-label "$h1" --target-h2-label "$h2" \
    --panel-label ""
  local svg="$TMP/${id}.whole_genome_homologous_context_ribbon.svg"
  local cropped="$TMP/${base}.crop.svg"
  # crop the viewBox/height to the content band
  python3 - "$svg" "$cropped" "$CROP_Y0" "$CROP_H" "$CROP_X0" "$CROP_W" <<'PY'
import sys,re
svg,out,y0,h,x0,w=sys.argv[1:7]
s=open(svg).read()
s=re.sub(r'viewBox="0 0 3600 840"', f'viewBox="{x0} {y0} {w} {h}"', s, count=1)
s=re.sub(r'width="3600"',  f'width="{w}"',  s, count=1)
s=re.sub(r'height="840"',  f'height="{h}"', s, count=1)
# the "100%" background rect is anchored at (0,0); after shifting the viewBox it
# no longer covers the whole visible region -> replace with a full-crop white rect
s=s.replace('<rect width="100%" height="100%" fill="white"/>',
            f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="white"/>', 1)
open(out,'w').write(s)
PY
  inkscape --export-type=pdf --export-filename="$OUT/${base}.pdf" "$cropped"
  inkscape --export-type=png -w 3600 --export-filename="$OUT/${base}.png" "$cropped"
}

gen PAN027pat_vs_PAN011_joint "PAN027 paternal" "PAN011 h1" "PAN011 h2" Fig5B_paternal
gen PAN027mat_vs_PAN010_joint "PAN027 maternal" "PAN010 h1" "PAN010 h2" Fig5C_maternal
gen PAN028mat_vs_PAN027_joint "PAN028 maternal" "PAN027 h1" "PAN027 h2" Fig5D_pan028

echo "done: $OUT/Fig5{B_paternal,C_maternal,D_pan028}.{pdf,png}"
