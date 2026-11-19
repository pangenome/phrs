#!/bin/bash
# Fig5 KID panel with the real KID.m
# run on octopus01: bash fig5_kidm_submit.sh ; needs scripts_tmp/{fig5_align_f,fig5_align_mm2c,fig5_similarity,fig5_similarity_s}.sbatch
# KID.m per chromosome from kid_phase.sbatch (v1.1 and v1.3.1 agree, = v1.0 trio table in PEDIGREE_NAMING.md):
# hap2 for chr1 3 4 5 7 9 12 14 15 18 21 22 X, hap1 for the rest. The paper query (prepare_inputs.py) is the whole hap1.
# Query contigs keep their names (PAN028#1#chrN or PAN028#2#chrN); target = inputs_<v> target (SELF joint).
# out: $D/fig5_kidm/inputs_<v>/, f48_<v>/ (FastGA -f48, paper ani step), f48cnd_<v>/, mm2c_raw_<v>/, mm2ccnd_<v>/
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp; K=$D/fig5_kidm; mkdir -p $K logs
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; C=PAN028mat_vs_PAN027_joint
SW=$T/sweepga_fix2/sweepga.mode; MM=$T/minimap2-2.31/minimap2; H2=" chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX "
for v in v11 v131; do
  I=$K/inputs_$v; mkdir -p $I
  if [ ! -s $I/$C.query.fa.fai ]; then
    cut -f1 norm_$v/PAN028#1.fa.gz.fai | while read -r c; do n=${c##*#}; case "$H2" in *" $n "*) echo "PAN028#2#$n";; *) echo $c;; esac; done > $I/names.txt
    for h in 1 2; do $S faidx norm_$v/PAN028#$h.fa.gz $(grep "^PAN028#$h#" $I/names.txt); done > $I/$C.query.fa
    $S faidx $I/$C.query.fa; ln -sf $D/inputs_$v/$C.target.fa $I/$C.target.fa; ln -sf $D/inputs_$v/$C.target.fa.fai $I/$C.target.fa.fai
  fi
  echo "$v query: $(wc -l < $I/$C.query.fa.fai) contigs, hap2 $(grep -c '#2#' $I/$C.query.fa.fai), $(awk '{s+=$2} END {print s}' $I/$C.query.fa.fai) bp"
  A=$(sbatch --parsable --exclude=octopus05 -t 3-00:00:00 -J kidm-f48-$v -o logs/kidm.f48.$v.aln.%j.log $T/fig5_align_f.sbatch $C $I $K/f48_$v 48)
  P=$(sbatch --parsable -J kidm-f48sim-$v -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/kidm.f48.$v.sim.%j.log \
      $T/fig5_similarity.sbatch $C $I $K/f48_$v/raw_paf/$C.frequency48.paf.gz $K/f48_$v)
  Q=$(sbatch --parsable -J kidm-f48cnd-$v -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/kidm.f48cnd.$v.%j.log \
      $T/fig5_similarity_s.sbatch $C $I $K/f48_$v/raw_paf/$C.frequency48.paf.gz $K/f48cnd_$v ani $SW candidate)
  B=$(sbatch --parsable --exclude=octopus05 -J kidm-mm2c-$v -o logs/kidm.mm2c.$v.aln.%j.log $T/fig5_align_mm2c.sbatch $C $I $K/mm2c_raw_$v $MM)
  R=$(sbatch --parsable -J kidm-mm2ccnd-$v -p tux,highmem -c 48 --mem=480G -d afterok:$B -o logs/kidm.mm2ccnd.$v.%j.log \
      $T/fig5_similarity_s.sbatch $C $I $K/mm2c_raw_$v/raw_paf/$C.mm2c.paf.gz $K/mm2ccnd_$v ani $SW candidate)
  echo "$v f48 align=$A sim=$P cnd=$Q; mm2c align=$B cnd=$R"
done
echo SUBMITDONE
