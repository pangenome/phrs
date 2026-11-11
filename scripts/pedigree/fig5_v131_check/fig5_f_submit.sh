#!/bin/bash
# Submit Fig5 steps 4-7 at several -f, both versions
# run on octopus01: bash fig5_f_submit.sh [FREQS] ; needs scripts_tmp/{fig5_align_f,fig5_similarity}.sbatch
# jobs: align (workers,tux) then similarity (tux,highmem, afterok), 3 comparisons x 2 versions x FREQS + patched v1.1 PAN027pat
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; T=$D/scripts_tmp; cd $D; mkdir -p logs
FREQS=${1:-"48 96"}
CIDS="PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint"
sub () { # tag IN OUT CID F
  local tag=$1 IN=$2 O=$3 cid=$4 F=$5 A S
  A=$(sbatch --parsable -J f$F-$tag-aln -o logs/f$F.$tag.$cid.aln.%j.log $T/fig5_align_f.sbatch $cid $IN $O $F)
  S=$(sbatch --parsable -J f$F-$tag-sim -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/f$F.$tag.$cid.sim.%j.log \
      $T/fig5_similarity.sbatch $cid $IN $O/raw_paf/$cid.frequency$F.paf.gz $O)
  echo "f$F $tag $cid align=$A sim=$S"
}
for F in $FREQS; do
  for v in v11 v131; do for cid in $CIDS; do sub ${v}_${cid%%_vs*} $D/inputs_$v $D/fig5_f$F/out_$v $cid $F; done; done
  sub v11_PAN027pat_patched $D/patch_v11/inputs_v11_patched $D/fig5_f$F/out_v11_patched PAN027pat_vs_PAN011_joint $F
done
echo SUBMITDONE
