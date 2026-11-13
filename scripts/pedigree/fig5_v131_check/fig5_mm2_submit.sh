#!/bin/bash
# Fig5 with minimap2: align, then steps 5-7 (ani, log-length-ani)
# run on octopus01: bash fig5_mm2_submit.sh ; needs scripts_tmp/{fig5_align_mm2,fig5_similarity_s}.sbatch, minimap2-2.31
# out: $D/fig5_mm2/raw_<v>/raw_paf/, $D/fig5_mm2/<ani|lla>_<v>/ (v = v11, v131, v11p)
set -uo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp; MM=$T/minimap2-2.31/minimap2
P=PAN027pat_vs_PAN011_joint; M=PAN027mat_vs_PAN010_joint; K=PAN028mat_vs_PAN027_joint
sub () { # v cid inputs
  local v=$1 c=$2 in=$D/$3 A S1 S2
  A=$(sbatch --parsable -J mm2-$v-${c%%_vs*} -o logs/mm2.$v.$c.aln.%j.log $T/fig5_align_mm2.sbatch $c $in $D/fig5_mm2/raw_$v $MM)
  S1=$(sbatch --parsable -J mm2ani-$v-${c%%_vs*} -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/mm2.$v.$c.ani.%j.log \
       $T/fig5_similarity_s.sbatch $c $in $D/fig5_mm2/raw_$v/raw_paf/$c.mm2.paf.gz $D/fig5_mm2/ani_$v ani)
  S2=$(sbatch --parsable -J mm2lla-$v-${c%%_vs*} -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/mm2.$v.$c.lla.%j.log \
       $T/fig5_similarity_s.sbatch $c $in $D/fig5_mm2/raw_$v/raw_paf/$c.mm2.paf.gz $D/fig5_mm2/lla_$v log-length-ani)
  echo "$v $c align=$A ani=$S1 lla=$S2"
}
for c in $P $M $K; do sub v11 $c inputs_v11; sub v131 $c inputs_v131; done
sub v11p $P patch_v11/inputs_v11_patched
echo SUBMITDONE
