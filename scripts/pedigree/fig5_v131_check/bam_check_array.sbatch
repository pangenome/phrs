#!/bin/bash
#SBATCH -p workers,tux,lowmem
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 12:00:00
# Full read of one BAM per array task
# args: WASHU, BAM_LIST, OUT_DIR ; sbatch --array=1-N%K
# One samtools pass per BAM: records, primary (-F 0x904), primary MAPQ>=20; READ_ERROR on non-zero exit or stderr.
# out: OUT_DIR/<bam>.tsv (bam, status, records, primary, primary_q20, seconds)
set -uo pipefail
D=$1; L=$2; O=$3; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
b=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $L); n=$(basename $b); t0=$(date +%s)
echo "[$(date -Is)] task $SLURM_ARRAY_TASK_ID $b on $(hostname)"
r=$($S view -@ 6 $b 2>$O/$n.err | awk -F'\t' '{n++; f=$2; if (int(f/4)%2==0 && int(f/256)%2==0 && int(f/2048)%2==0) {p++; if ($5>=20) q++}} END {printf "%d\t%d\t%d", n, p, q}'; echo " ${PIPESTATUS[0]}")
e=${r##* }; r=${r% *}
st=OK; { [ $e -ne 0 ] || [ -s $O/$n.err ]; } && st=READ_ERROR
printf "%s\t%s\t%s\t%d\n" $n $st "$r" $(( $(date +%s) - t0 )) > $O/$n.tsv
echo "[$(date -Is)] $n $st $r exit $e"
