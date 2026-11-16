#!/bin/bash
# Rebuild Fig5 panels A-E in temp tree
# args: run dir with class winners, fai dir
set -euo pipefail
S=$1; FAI=$2; T=$S/tree; R=/home/guarracino/git/phrs
cd $R
rm -rf $T $S/svg
mkdir -p $T/data $T/submission/scripts/figures $T/submission/fig/MainFigures $S/svg
cp submission/scripts/figures/{fig5_ribbon_plotter.py,make_fig5_zoom.R,make_fig5_pedigree.R} $T/submission/scripts/figures/
for c in PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint; do
  cp $S/$c.class_winners.impg_similarity.tsv.gz $T/data/fig5_$c.class_winners.impg_similarity.tsv.gz
  cp $FAI/fig5_$c.query.fa.fai $FAI/fig5_$c.target.fa.fai $T/data/
done
# panels C, D, E as make_fig5_panels.sh
gen () { local id=$1 q=$2 h1=$3 h2=$4 base=$5
  FIG5_HIDE_LEGEND=1 python3 $T/submission/scripts/figures/fig5_ribbon_plotter.py --comparison-id $id \
    --class-winners $T/data/fig5_$id.class_winners.impg_similarity.tsv.gz --query-fai $T/data/fig5_$id.query.fa.fai \
    --target-fai $T/data/fig5_$id.target.fa.fai --output-dir $S/svg --query-label "$q" --target-h1-label "$h1" \
    --target-h2-label "$h2" --panel-label "" 2>&1 | grep -E "segments=" || true
  python3 - $S/svg/$id.whole_genome_homologous_context_ribbon.svg $S/svg/$base.crop.svg 75 565 165 3435 <<'PY'
import sys,re
svg,out,y0,h,x0,w=sys.argv[1:7]
s=open(svg).read()
s=re.sub(r'viewBox="0 0 3600 840"', f'viewBox="{x0} {y0} {w} {h}"', s, count=1)
s=re.sub(r'width="3600"',  f'width="{w}"',  s, count=1)
s=re.sub(r'height="840"',  f'height="{h}"', s, count=1)
s=s.replace('<rect width="100%" height="100%" fill="white"/>', f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="white"/>', 1)
open(out,'w').write(s)
PY
  inkscape --export-type=png -w 3600 --export-filename=$T/submission/fig/MainFigures/$base.png $S/svg/$base.crop.svg >/dev/null 2>&1
}
gen PAN027pat_vs_PAN011_joint "PAN027 paternal" "PAN011 h1" "PAN011 h2" Fig5C_paternal
gen PAN027mat_vs_PAN010_joint "PAN027 maternal" "PAN010 h1" "PAN010 h2" Fig5D_maternal
gen PAN028mat_vs_PAN027_joint "PAN028 maternal" "PAN027 h1" "PAN027 h2" Fig5E_pan028
# panels A and B
Rscript $T/submission/scripts/figures/make_fig5_zoom.R >/dev/null 2>&1; echo "zoom exit=$?"
Rscript $T/submission/scripts/figures/make_fig5_pedigree.R >/dev/null 2>&1; echo "pedigree exit=$?"
echo "=== pixel diff vs paper ==="
for p in Fig5A_pedigree Fig5B_zoom9q3q Fig5C_paternal Fig5D_maternal Fig5E_pan028; do
  printf "%s: " $p; magick compare -metric AE $T/submission/fig/MainFigures/$p.png submission/fig/MainFigures/$p.png null: 2>&1 || true
  printf " of %s\n" "$(magick identify -format '%[fx:w*h]' submission/fig/MainFigures/$p.png)"
done
