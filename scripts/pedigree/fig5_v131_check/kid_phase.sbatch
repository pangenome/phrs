#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 12:00:00
# KID maternal copy per chromosome
# args: WASHU, VERSION (v11|v131), OUT_DIR
# minimap2 2.31 -x asm5: PAN028 hap1 + hap2 (query) on PAN027 hap1 + hap2 (target, SELF).
# Per KID contig: aligned bp and identity (sum matches / sum block length, primary, MAPQ >= 5)
# on the same-name SELF chromosome, per SELF hap. The KID.m copy is the one at transmission grade (>= 0.998).
# out: OUT_DIR/kid_phase.<v>.tsv (kid_contig, self_hap, aligned_bp, identity, kid_len)
set -euo pipefail
D=$1; V=$2; O=$3; cd $D; mkdir -p $O
M=$D/scripts_tmp/minimap2-2.31/minimap2
W=/scratch/$USER.kidphase.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
zcat norm_$V/PAN027#1.fa.gz norm_$V/PAN027#2.fa.gz > $W/self.fa
zcat norm_$V/PAN028#1.fa.gz norm_$V/PAN028#2.fa.gz > $W/kid.fa
log "SELF $(grep -c '>' $W/self.fa) contigs, KID $(grep -c '>' $W/kid.fa) contigs, host $(hostname)"
$M -t 32 -I 8G -x asm5 -c --secondary=no $W/self.fa $W/kid.fa 2> $W/mm2.log > $W/kid_on_self.paf
log "PAF lines $(wc -l < $W/kid_on_self.paf); minimap2 log tail: $(tail -2 $W/mm2.log | tr '\n' ' ')"
cat norm_$V/PAN028#1.fa.gz.fai norm_$V/PAN028#2.fa.gz.fai | cut -f1,2 > $W/kid.len
awk -F'\t' -v OFS='\t' 'NR==FNR {L[$1]=$2; next} $12>=5 {split($1,q,"#"); split($6,t,"#"); if (q[3]!=t[3]) next;
  k=$1"\t"t[2]; m[k]+=$10; b[k]+=$11; a[k]+=$4-$3} END {for (k in m) {split(k,x,"\t"); print k, a[k], sprintf("%.5f", m[k]/b[k]), L[x[1]]}}' \
  $W/kid.len $W/kid_on_self.paf | sort -t'#' -k3,3V -k2,2n > $O/kid_phase.$V.tsv
sed -i '1i kid_contig\tself_hap\taligned_bp\tidentity\tkid_len' $O/kid_phase.$V.tsv
log "rows $(($(wc -l < $O/kid_phase.$V.tsv) - 1))"
cat $O/kid_phase.$V.tsv
