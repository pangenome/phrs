#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=300G
#SBATCH -t 1-00:00:00
# FastGA chr5 chain end at one -f
# args: WASHU, OUT_DIR, VERSION (v11|v131), FREQ
# Query PAN027#2#chr5, target joint PAN011 diploid. Same binary and options as fastga_chr5_stop.sbatch.
# out: chain.<v>.f<F>.tsv (longest chain to h2_chr5), fastga.<v>.f<F>.log, chr5.<v>.f<F>.paf.gz
set -euo pipefail
D=$1; OUT=$2; v=$3; F=$4; cd $D; mkdir -p $OUT
B=/home/erikg/.cache/sweepga/0.1.1_fe5037d5e_wd47b7e3e; export PATH=$B:$PATH
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.fgf.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host $(hostname), $v -f$F, scratch $(df -h --output=avail $W | tail -1)"
$S faidx norm_$v/PAN027#2.fa.gz PAN027#2#chr5 > $W/q.fa
cp inputs_$v/PAN027pat_vs_PAN011_joint.target.fa $W/t.fa
t0=$(date +%s)
FastGA -v -k -T$SLURM_CPUS_PER_TASK -P$W -f$F -paf $W/q.fa $W/t.fa > $W/o.paf 2> $W/fg.log || { log "ERROR FastGA"; tail -5 $W/fg.log; exit 1; }
sec=$(( $(date +%s) - t0 ))
cp $W/fg.log $OUT/fastga.$v.f$F.log; gzip -c $W/o.paf > $OUT/chr5.$v.f$F.paf.gz
n=$(wc -l < $W/o.paf)
awk -F'\t' -v v=$v -v f=$F -v s=$sec -v n=$n '$6 ~ /h2_chr5$/ && $11 > b {b=$11; r=sprintf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f", v, f, s, n, $3, $4, $8, $9, $10/$11)} END {print r}' $W/o.paf > $OUT/chain.$v.f$F.tsv
log "done in ${sec}s, $n PAF lines; version freq seconds lines qstart qend tstart tend identity: $(cat $OUT/chain.$v.f$F.tsv)"
