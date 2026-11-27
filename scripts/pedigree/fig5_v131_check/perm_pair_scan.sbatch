#!/usr/bin/env bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11
# R1-M2 false-positive control: the final Fig5 scan on a permuted parent-child pair, i.e. a child haplotype against the
# parent who did not transmit it (no identical-by-descent homolog in the target).
# usage: sbatch perm_pair_scan.sbatch CID QUERY_CID TARGET_CID
#   query  = fig5_final/inputs_v11f3/QUERY_CID.query.fa  (the final corrected inputs)
#   target = fig5_final/inputs_v11f3/TARGET_CID.target.fa
#   e.g. PAN027pat_vs_PAN010_joint PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint (SELF.p against MOM)
# Inputs are symlinks in fig5_final/inputs_perm (no copies). Same steps and parameters as the real scans: fig5_align_chunk.sbatch
# (--freq 32 --piece 2000000 --batches 8, all groups in this job), fig5_similarity_s.sbatch (ani, N=18). Everything runs in node
# /scratch; only class winners, filtered PAF, beds, summaries and the alignment coverage go to fig5_final/perm/CID/.
# Runs are drawn afterwards with fig5_runs.sh on the copied class winners.
set -euxo pipefail
CID=$1; QC=$2; TC=$3
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; F=$D/fig5_final; SC=$D/scripts_tmp; IN=$F/inputs_perm; O=$F/perm/$CID
log () { echo "[$(date '+%F %T')] $*"; }
log "host=$(hostname) cid=$CID query=$QC target=$TC cpus=$SLURM_CPUS_PER_TASK"
[ ! -e $O ] || { echo "ERROR $O exists"; exit 1; }
mkdir -p $IN $O
for x in query.fa query.fa.fai; do ln -s ../inputs_v11f3/$QC.$x $IN/$CID.$x; done
for x in target.fa target.fa.fai; do ln -s ../inputs_v11f3/$TC.$x $IN/$CID.$x; done
ls -lL $IN/$CID.*
S=/scratch/$USER/tmp/phrs/perm.$SLURM_JOB_ID; mkdir -p $S; trap 'rc=$?; rm -rf $S; echo "exit status $rc"' EXIT
log "step 1 alignment"
bash $SC/fig5_align_chunk.sbatch --cid $CID --in $IN --out $S/out --freq 32 --piece 2000000 --batches 8
RAW=$S/out/raw_paf/$CID.frequency32.piece2000000.batches8.paf.gz
log "step 2 N=18 similarity on $RAW ($(stat -c %s $RAW) bytes)"
bash $SC/fig5_similarity_s.sbatch $CID $IN $RAW $S/n18 ani /home/erikg/.cargo/bin/sweepga "" 18
log "step 3 copy results to $O"
cp $S/out/raw_paf/$CID.frequency32.piece2000000.batches8.coverage.tsv $O/
cp -r $S/n18/outputs $S/n18/summaries $S/n18/beds $S/n18/filtered_paf $O/
ls -la $O/*/
log "DONE $CID"
