#!/bin/bash
# Fig5 steps 5-7 with log-length-ani on existing raw PAFs
# run on octopus01: bash fig5_lla_submit.sh ; needs scripts_tmp/fig5_similarity_s.sbatch
# out: $D/fig5_lla/<f>_<v>/ (absolute: the job cds to scratch); skips raw PAFs that are missing or unreadable (out_B_v11 pat, out_C_v131 PAN028mat)
set -uo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D
P=PAN027pat_vs_PAN011_joint; M=PAN027mat_vs_PAN010_joint; K=PAN028mat_vs_PAN027_joint
sub () { # tag cid inputs raw
  [ -s "$4" ] || { echo "skip $1 $2: no $4"; return; }
  echo "$1 $2 $(sbatch --parsable -J lla-$1-${2%%_vs*} -p tux,highmem -c 48 --mem=480G -o logs/lla.$1.$2.%j.log \
    scripts_tmp/fig5_similarity_s.sbatch $2 $D/$3 $D/$4 $D/fig5_lla/$1 log-length-ani)"
}
sub f32_v11 $P inputs_v11 regen_B_v11/raw_paf/$P.frequency32.paf.gz
sub f32_v11 $M inputs_v11 out_B_v11/raw_paf/$M.frequency32.paf.gz
sub f32_v11 $K inputs_v11 out_B_v11/raw_paf/$K.frequency32.paf.gz
sub f32_v131 $P inputs_v131 out_C_v131/raw_paf/$P.frequency32.paf.gz
sub f32_v131 $M inputs_v131 out_C_v131/raw_paf/$M.frequency32.paf.gz
sub f32_v11p $P patch_v11/inputs_v11_patched patch_v11/out_P_v11/raw_paf/$P.frequency32.paf.gz
for v in v11 v131; do for c in $P $M $K; do sub f48_$v $c inputs_$v fig5_f48/out_$v/raw_paf/$c.frequency48.paf.gz; done; done
sub f48_v11p $P patch_v11/inputs_v11_patched fig5_f48/out_v11_patched/raw_paf/$P.frequency48.paf.gz
echo SUBMITDONE
