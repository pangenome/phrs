#!/bin/bash
# Fig5C on patched DAD v1.3.1
# run on octopus01: bash fig5_v131p_submit.sh ; needs scripts_tmp/{fig5_align_f,fig5_align_mm2c,fig5_similarity,fig5_similarity_s}.sbatch
# Inputs: patch_v131/inputs_v131_patched (patch_dad_chr9.sbatch v131). Same steps as the v1.1 patched runs:
# FastGA -f48 + paper step (fig5_f48/out_v131_patched), + cnd (fig5_cnd/f48_v131p); minimap2 chunks + cnd (fig5_mm2c/cnd_v131p).
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp; C=PAN027pat_vs_PAN011_joint
I=$D/patch_v131/inputs_v131_patched; SW=$T/sweepga_fix2/sweepga.mode; MM=$T/minimap2-2.31/minimap2
echo "query $(wc -l < $I/$C.query.fa.fai) contigs, target $(wc -l < $I/$C.target.fa.fai) contigs"
A=$(sbatch --parsable --exclude=octopus05 -t 3-00:00:00 -J v131p-f48 -o logs/v131p.f48.aln.%j.log $T/fig5_align_f.sbatch $C $I $D/fig5_f48/out_v131_patched 48)
P=$(sbatch --parsable -J v131p-f48sim -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/v131p.f48.sim.%j.log \
    $T/fig5_similarity.sbatch $C $I $D/fig5_f48/out_v131_patched/raw_paf/$C.frequency48.paf.gz $D/fig5_f48/out_v131_patched)
Q=$(sbatch --parsable -J v131p-f48cnd -p tux,highmem -c 48 --mem=480G -d afterok:$A -o logs/v131p.f48cnd.%j.log \
    $T/fig5_similarity_s.sbatch $C $I $D/fig5_f48/out_v131_patched/raw_paf/$C.frequency48.paf.gz $D/fig5_cnd/f48_v131p ani $SW candidate)
B=$(sbatch --parsable --exclude=octopus05 -J v131p-mm2c -o logs/v131p.mm2c.aln.%j.log $T/fig5_align_mm2c.sbatch $C $I $D/fig5_mm2c/raw_v131p $MM)
R=$(sbatch --parsable -J v131p-mm2ccnd -p tux,highmem -c 48 --mem=480G -d afterok:$B -o logs/v131p.mm2ccnd.%j.log \
    $T/fig5_similarity_s.sbatch $C $I $D/fig5_mm2c/raw_v131p/raw_paf/$C.mm2c.paf.gz $D/fig5_mm2c/cnd_v131p ani $SW candidate)
echo "f48 align=$A sim=$P cnd=$Q; mm2c align=$B cnd=$R"
