#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 6:00:00
# SELF chr5 end vs the KID chr5q~chr1p tract
# args: WASHU, OUT_DIR ; v1.1
# 1) SELF reads (PAN027.onPAN028) spanning KID tract PAN028#1#chr5:183,044,000-183,074,000 +-5 kb, MAPQ>=20 primary:
#    their alignments in SELF own BAMs (PAN027.v11, ont_full): contig, start, end, MAPQ, flag, clips.
# 2) KID chr5 last 200 kb (tract included) aligned to the SELF diploid (minimap2 asm20, secondaries).
# 3) SELF chr5 ends (h1, h2, last 300 kb) aligned to the KID diploid.
# out: OUT_DIR/self_reads_on_self.tsv, kid_end_on_self.paf, self_ends_on_kid.paf, summary in log
set -euo pipefail
D=$1; O=$2; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
W=/scratch/$USER.selfchr5.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
T0=183044000; T1=183074000
for b in ont_full/SRR*.PAN027.onPAN028.v11.bam; do $S view -q 20 -F 0x904 $b PAN028#1#chr5:$T0-$T1; done \
  | awk -v a=$((T0-5000)) -v z=$((T1+5000)) '{e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; if ($4<=a && e>=z) print $1}' | sort -u > $W/names.txt
log "SELF reads spanning the KID tract +-5 kb: $(wc -l < $W/names.txt)"
printf "read\tbam\tflag\tcontig\tstart\tend\tmapq\tclip_left\tclip_right\n" > $O/self_reads_on_self.tsv
for b in ont_full/SRR*.PAN027.v11.bam; do
  $S view -N $W/names.txt $b | awk -F'\t' -v bam=$(basename $b) -v OFS='\t' '{e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)};
    l=0; r=0; if ($6 ~ /^[0-9]+[SH]/) l=$6+0; if (match($6,/[0-9]+[SH]$/)) r=substr($6,RSTART)+0; print $1, bam, $2, $3, $4, e, $5, l, r}' >> $O/self_reads_on_self.tsv
done
log "their alignments on SELF: $(($(wc -l < $O/self_reads_on_self.tsv) - 1)) records; by contig: $(tail -n +2 $O/self_reads_on_self.tsv | cut -f4 | sort | uniq -c | tr -s ' ' | tr '\n' ';')"
L=$(awk '$1=="PAN028#1#chr5" {print $2}' norm_v11/PAN028#1.fa.gz.fai)
$S faidx norm_v11/PAN028#1.fa.gz PAN028#1#chr5:$((L-200000))-$L > $W/kid_end.fa
zcat norm_v11/PAN027#1.fa.gz norm_v11/PAN027#2.fa.gz > $W/self.fa
$MM -c -x asm20 -N 50 -p 0.3 --secondary=yes -t 16 $W/self.fa $W/kid_end.fa 2> $W/mm.log > $O/kid_end_on_self.paf
log "KID chr5 last 200 kb (start $((L-200000))) on SELF, records >= 5 kb: $(awk -F'\t' '$11>=5000 {printf "q%d-%d %s %s:%d-%d %.4f; ", $3,$4,$5,$6,$8,$9,$10/$11}' $O/kid_end_on_self.paf)"
for h in 1 2; do Lh=$(awk -v c=PAN027#$h#chr5 '$1==c {print $2}' norm_v11/PAN027#$h.fa.gz.fai); $S faidx norm_v11/PAN027#$h.fa.gz PAN027#$h#chr5:$((Lh-300000))-$Lh; done > $W/self_ends.fa
zcat norm_v11/PAN028#1.fa.gz norm_v11/PAN028#2.fa.gz > $W/kid.fa
$MM -c -x asm20 -N 50 -p 0.3 --secondary=yes -t 16 $W/kid.fa $W/self_ends.fa 2>> $W/mm.log > $O/self_ends_on_kid.paf
log "SELF chr5 ends (last 300 kb) on KID, records >= 20 kb: $(awk -F'\t' '$11>=20000 {printf "%s q%d-%d %s %s:%d-%d %.4f; ", $1,$3,$4,$5,$6,$8,$9,$10/$11}' $O/self_ends_on_kid.paf)"
log done
