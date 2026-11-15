#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# Polish a patched segment with HiFi
# args: WASHU, SAMPLE, VERSION (v11), IN_NORM (patched norm dir), SEGMENT (contig:start-end, 1-based), PICK_REGIONS (comma list), OUT_DIR
# 1. HiFi reads of SAMPLE (hifi/bam/*.<sample>.<v>.bam, original assembly) in PICK_REGIONS, as FASTQ.
# 2. Realigned to the patched diploid (IN_NORM) with minimap2 2.31 -ax map-hifi.
# 3. In SEGMENT: bcftools 1.21 call --ploidy 1 (SNV and indel, QUAL >= 20), counts by TYPE only; bcftools consensus writes the
#    polished haplotype FASTA (alleles are never printed). 4. Realign the same reads to the polished diploid and count again.
# out: OUT_DIR/norm/<sample>#<h>.fa.gz (+fai, gzi; other haplotype copied), OUT_DIR/polish.tsv
set -euo pipefail
D=$1; SM=$2; V=$3; IN=$4; SEG=$5; PICK=$6; O=$7; cd $D; mkdir -p $O/norm
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
MM=$D/scripts_tmp/minimap2-2.31/minimap2; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1); T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.hpol.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
C=${SEG%%:*}; H=$(echo $C | cut -d'#' -f2)
log "start $SM $V $SEG, host $(hostname); pick regions: $PICK"
for r in $(echo $PICK | tr ',' ' '); do c=${r%%:*}; e=${r##*-}; L=$(cat $D/norm_$V/${c%#chr*}.fa.gz.fai | awk -v c="$c" '$1==c {print $2}'); [ -n "$L" ] && [ "$e" -le "$L" ] || { log "ERROR region $r: contig length ${L:-missing}"; exit 1; }; done
for b in hifi/bam/*.$SM.$V.bam; do $S view -b -F 0x900 $b $(echo $PICK | tr ',' ' '); done | $S fastq - 2> /dev/null | awk 'NR%4==1 {n=$1} NR%4==1 && seen[n]++ {skip=1; next} NR%4==1 {skip=0} !skip' > $W/reads.fq
log "$SM $V segment $SEG, reads $(awk 'NR%4==2' $W/reads.fq | wc -l) from $PICK"
aln () { # ref_dir tag
  zcat $1/$SM#1.fa.gz $1/$SM#2.fa.gz > $W/ref.$2.fa; $S faidx $W/ref.$2.fa
  $MM -ax map-hifi --eqx -t $T $W/ref.$2.fa $W/reads.fq 2> /dev/null | $S sort -@ 4 -o $W/$2.bam - && $S index $W/$2.bam
  $B mpileup -B -Ou -f $W/ref.$2.fa -r "$SEG" -q 20 --max-depth 2000 -a AD,DP $W/$2.bam 2> /dev/null | $B call -m --ploidy 1 -v -Oz -o $W/$2.vcf.gz 2> /dev/null
  $B index $W/$2.vcf.gz; $B view -i 'QUAL>=20' -Oz -o $W/$2.q20.vcf.gz $W/$2.vcf.gz; $B index $W/$2.q20.vcf.gz
  echo "$2 snv $($B view -H -v snps $W/$2.q20.vcf.gz | wc -l) indel $($B view -H -v indels $W/$2.q20.vcf.gz | wc -l) depth_median $($S depth -a -Q 20 -r "$SEG" $W/$2.bam | cut -f3 | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')"
}
R1=$(aln $IN before); log "$R1"
$S faidx $W/ref.before.fa $(cut -f1 $IN/$SM#$H.fa.gz.fai) | $B consensus -f - $W/before.q20.vcf.gz 2> $W/cons.log > $W/pol.fa
log "consensus: $(tr '\n' ' ' < $W/cons.log | cut -c1-200)"
$BGZ -c $W/pol.fa > $O/norm/$SM#$H.fa.gz; $S faidx $O/norm/$SM#$H.fa.gz
o=$((3 - H)); cp $IN/$SM#$o.fa.gz $IN/$SM#$o.fa.gz.fai $IN/$SM#$o.fa.gz.gzi $O/norm/
R2=$(aln $O/norm after); log "$R2"
printf "sample\tversion\tsegment\tstep\tsnv\tindel\tdepth_median\n" > $O/polish.tsv
for r in "$R1" "$R2"; do set -- $r; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $SM $V "$SEG" $1 $3 $5 $7 >> $O/polish.tsv; done
cat $O/polish.tsv; log "lengths: before $(grep -P "^$C\t" $IN/$SM#$H.fa.gz.fai | cut -f2), after $(grep -P "^$C\t" $O/norm/$SM#$H.fa.gz.fai | cut -f2)"
