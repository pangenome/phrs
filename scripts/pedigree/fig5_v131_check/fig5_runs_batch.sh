#!/bin/bash
# Fig5 drawn runs for all finished f runs, with comparisons
# run on octopus01: bash /tmp/guarracino.fig5_runs_batch.sh ; needs /tmp/guarracino.{fig5_runs.sh,compare_runs.py}
# tags: f32_<v> originals (out_B_v11, out_C_v131, patch_v11/out_P_v11 as f32_v11p), f48_*, f96_*; skips finished runs
set -uo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; R=$D/fig5_f48/runs_check; W=class_winners.impg_similarity.tsv.gz
ids="PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint"
dir_of () { case $1 in f32_v11) echo out_B_v11;; f32_v131) echo out_C_v131;; f32_v11p) echo patch_v11/out_P_v11;;
  f48_v11) echo fig5_f48/out_v11;; f48_v131) echo fig5_f48/out_v131;; f48_v11p) echo fig5_f48/out_v11_patched;;
  f96_v11) echo fig5_f96/out_v11;; f96_v131) echo fig5_f96/out_v131;; f96_v11p) echo fig5_f96/out_v11_patched;;
  lla_*) echo fig5_lla/${1#lla_};; esac; }
fai_of () { case $1 in *_v131) echo inputs_v131;; *_v11p) echo patch_v11/inputs_v11_patched;; *) echo inputs_v11;; esac; }
for t in f32_v11 f32_v131 f32_v11p f48_v11 f48_v131 f48_v11p f96_v11 f96_v131 f96_v11p lla_f32_v11 lla_f32_v131 lla_f32_v11p lla_f48_v11 lla_f48_v131 lla_f48_v11p; do
  for c in $ids; do
    cw=$(dir_of $t)/outputs/$c.$W; r=$R/$t/$c.whole_genome_ribbon_runs.tsv
    [ -s $cw ] || continue; [ -s $r ] && continue
    bash /tmp/guarracino.fig5_runs.sh $cw $c $(fai_of $t) $R/$t 2>&1 | grep -E "segments=|rror|missing"
  done
done
cmp () { # tagA tagB cid
  local a=$R/$1/$3.whole_genome_ribbon_runs.tsv b=$R/$2/$3.whole_genome_ribbon_runs.tsv
  [ -s $a ] && [ -s $b ] || return 0
  echo "== $3 $1 vs $2"; python3 /tmp/guarracino.compare_runs.py $a $b | awk -F'\t' 'NF<9 || $9 !~ /acro_acro/'
  echo "acro_acro only in $2: $(python3 /tmp/guarracino.compare_runs.py $a $b | grep -c '^B_only.*acro_acro'), only in $1: $(python3 /tmp/guarracino.compare_runs.py $a $b | grep -c '^A_only.*acro_acro')"
}
for c in $ids; do for v in v11 v131 v11p; do cmp f32_$v f48_$v $c; cmp f48_$v f96_$v $c; done; done
cmp f48_v11 f48_v11p PAN027pat_vs_PAN011_joint
for c in $ids; do for v in v11 v131 v11p; do cmp f32_$v lla_f32_$v $c; cmp f48_$v lla_f48_$v $c; cmp lla_f32_$v lla_f48_$v $c; done; done
cmp lla_f48_v11 lla_f48_v11p PAN027pat_vs_PAN011_joint
