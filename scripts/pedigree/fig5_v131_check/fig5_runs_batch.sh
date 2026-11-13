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
  lla_*) echo fig5_lla/${1#lla_};; cnd_*) echo fig5_cnd/${1#cnd_};;
  mm2ani_*) echo fig5_mm2/ani_${1#mm2ani_};; mm2lla_*) echo fig5_mm2/lla_${1#mm2lla_};; mm2cnd_*) echo fig5_mm2/cnd_${1#mm2cnd_};;
  mm2cani_*) echo fig5_mm2c/ani_${1#mm2cani_};; mm2ccnd_*) echo fig5_mm2c/cnd_${1#mm2ccnd_};; esac; }
fai_of () { case $1 in *_v131) echo inputs_v131;; *_v11p) echo patch_v11/inputs_v11_patched;; *) echo inputs_v11;; esac; }
for t in f32_v11 f32_v131 f32_v11p f48_v11 f48_v131 f48_v11p f96_v11 f96_v131 f96_v11p lla_f32_v11 lla_f32_v131 lla_f32_v11p lla_f48_v11 lla_f48_v131 lla_f48_v11p cnd_f32_v11 cnd_f32_v131 cnd_f32_v11p cnd_f48_v11 cnd_f48_v131 cnd_f48_v11p \
  mm2ani_v11 mm2ani_v131 mm2ani_v11p mm2lla_v11 mm2lla_v131 mm2lla_v11p mm2cnd_v11 mm2cnd_v131 mm2cnd_v11p mm2cani_v11 mm2cani_v131 mm2cani_v11p mm2ccnd_v11 mm2ccnd_v131 mm2ccnd_v11p; do
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
for c in $ids; do for v in v11 v131 v11p; do cmp lla_f48_$v cnd_f48_$v $c; cmp cnd_f32_$v cnd_f48_$v $c; cmp cnd_f48_$v mm2cnd_$v $c; cmp mm2ani_$v mm2cnd_$v $c; cmp mm2lla_$v mm2cnd_$v $c; done; done
cmp cnd_f48_v11 cnd_f48_v11p PAN027pat_vs_PAN011_joint; cmp mm2cnd_v11 mm2cnd_v11p PAN027pat_vs_PAN011_joint
for c in $ids; do for v in v11 v131 v11p; do cmp cnd_f48_$v mm2ccnd_$v $c; cmp mm2cani_$v mm2ccnd_$v $c; done; done
cmp mm2ccnd_v11 mm2ccnd_v11p PAN027pat_vs_PAN011_joint
