#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Polish all applied pieces, check
# usage: sbatch g1_polish_all.sbatch --washu D --sample SM --version V --in-dir APPLY_DIR --out DIR
# APPLY_DIR (g1_apply_all.py): SM#{1,2}.fa, corrections.tsv (applied rows: contig, start/end on norm_V, new_start/new_end).
# 1. reads: own HiFi and ONT primary reads at every applied piece (norm_V coordinates +- 50 kb).
# 2. HiFi on the corrected diploid (minimap2 2.31 map-hifi); bcftools 1.21 call --ploidy 1, QUAL >= 20, only inside the
#    pieces; bcftools consensus on both haplotypes (alleles never printed, counts by TYPE).
# 3. HiFi and ONT (-X ont) on the polished diploid; per piece SNV and indel before and after, on the pre-polish
#    piece windows (polish changes lengths by a few bp; > 1 kb on a contig exits with an error).
# out (new dir only): DIR/norm/SM#{1,2}.fa.gz (+fai, gzi), DIR/polish.tsv (one row per piece)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --in-dir) I=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V I O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O/polish.tsv ] && { log "exists $O/polish.tsv, skip"; exit 0; }
[ -e $O ] && { echo "ERROR $O exists without polish.tsv; use a new dir"; exit 1; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
MM=$D/scripts_tmp/minimap2-2.31/minimap2; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1); T=${SLURM_CPUS_PER_TASK:-4}
W=/scratch/$USER/tmp/phrs/polish.${SLURM_JOB_ID:-$$}; mkdir -p $W $O/norm; trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
awk -F'\t' '$1 == "applied"' $I/corrections.tsv > $W/app.tsv; NP=$(wc -l < $W/app.tsv); [ $NP -gt 0 ] || { echo "ERROR no applied pieces in $I/corrections.tsv"; exit 1; }
awk -F'\t' -v OFS='\t' '{s=$4-50000; print $3, (s < 0 ? 0 : s), $5+50000}' $W/app.tsv | sort -k1,1 -k2,2n > $W/pick.bed
awk -F'\t' -v OFS='\t' '{e=($9 > $8 ? $9 : $8 + 1); print $3, $8, e, $2}' $W/app.tsv | sort -k1,1 -k2,2n > $W/seg.bed
log "$SM $V: $NP applied pieces from $I; $(awk '{s+=$3-$2} END {print s}' $W/seg.bed) bp in pieces"
for s in hifi ont; do
  [ $s = hifi ] && BS=$(ls hifi/bam/*.$SM.$V.bam) || BS=$(ls ont_full/SRR*.$SM.$V.bam)
  # moosefs read errors: retry each BAM
  for b in $BS; do
    for t in 1 2 3; do
      if $S view -F 0x900 -L $W/pick.bed $b > $W/view.sam; then break; fi
      [ $t -lt 3 ] || { echo "ERROR samtools view failed 3 times on $b" >&2; exit 1; }
      log "samtools view failed on $b (try $t of 3), retry in 60 s" >&2; sleep 60
    done
    cat $W/view.sam
  done | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/$s.fq
  log "reads $s: $(awk 'NR%4==2' $W/$s.fq | wc -l)"
done
cat $I/$SM#1.fa $I/$SM#2.fa > $W/ref0.fa; $S faidx $W/ref0.fa
count () { # ref bam set tag -> per piece snv indel depth
  if [ $3 = hifi ]; then X="-B"; else X="-B -X ont"; fi
  $B mpileup $X -Ou -f $1 -R $W/seg.bed -q 20 --max-depth 2000 -a AD,DP $2 2> /dev/null | $B call -m --ploidy 1 -v -Oz -o $W/$4.vcf.gz 2> /dev/null
  $B index -f $W/$4.vcf.gz; $B view -i 'QUAL>=20' -Oz -o $W/$4.q20.vcf.gz $W/$4.vcf.gz; $B index -f $W/$4.q20.vcf.gz
  $B query -f '%CHROM\t%POS\t%TYPE\n' $W/$4.q20.vcf.gz | awk -F'\t' -v OFS='\t' 'NR==FNR {c[NR]=$1; a[NR]=$2; z[NR]=$3; id[NR]=$4; n=NR; next}
    {for (i=1; i<=n; i++) if ($1 == c[i] && $2 > a[i] && $2 <= z[i]) {if ($3 == "SNP") sv[i]++; else ind[i]++}} END {for (i=1; i<=n; i++) print id[i], sv[i]+0, ind[i]+0}' $W/seg.bed - > $W/$4.counts
}
$MM -ax map-hifi -t $T $W/ref0.fa $W/hifi.fq 2> /dev/null | $S sort -@ 4 -o $W/h0.bam - && $S index $W/h0.bam
count $W/ref0.fa $W/h0.bam hifi h0
for h in 1 2; do $S faidx $W/ref0.fa $(grep "^$SM#$h#" $W/ref0.fa.fai | cut -f1) | $B consensus -f - $W/h0.q20.vcf.gz 2> $W/cons.$h.log > $W/pol.$h.fa; done
log "consensus: $(cat $W/cons.*.log | tr '\n' ' ' | cut -c1-300)"
cat $W/pol.1.fa $W/pol.2.fa > $W/ref1.fa; $S faidx $W/ref1.fa
paste <(cut -f1,2 $W/ref0.fa.fai) <(cut -f1,2 $W/ref1.fa.fai) | awk '$1 != $3 || ($2 - $4) ^ 2 > 1e6 {print "ERROR polish changed contig "$1" "$2" -> "$3" "$4 > "/dev/stderr"; exit 1}'
$MM -ax map-hifi -t $T $W/ref1.fa $W/hifi.fq 2> /dev/null | $S sort -@ 4 -o $W/h1.bam - && $S index $W/h1.bam
$MM -ax map-ont -t $T $W/ref1.fa $W/ont.fq 2> /dev/null | $S sort -@ 4 -o $W/o1.bam - && $S index $W/o1.bam
count $W/ref1.fa $W/h1.bam hifi h1; count $W/ref1.fa $W/o1.bam ont o1
for h in 1 2; do $BGZ -@ 4 -c $W/pol.$h.fa > $O/norm/$SM#$h.fa.gz && $S faidx $O/norm/$SM#$h.fa.gz; done
printf "locus\thifi_snv_before\thifi_indel_before\thifi_snv_after\thifi_indel_after\tont_snv_after\n" > $O/polish.tsv
paste $W/h0.counts $W/h1.counts $W/o1.counts | awk -v OFS='\t' '{print $1, $2, $3, $5, $6, $8}' >> $O/polish.tsv
log "totals: $(awk 'NR>1 {a+=$2; b+=$3; c+=$4; d+=$5; e+=$6} END {print "HiFi SNV "a" -> "c", indel "b" -> "d", ONT SNV after "e}' $O/polish.tsv); length change $(awk '{s+=$2} END {print s}' $W/ref1.fa.fai) - $(awk '{s+=$2} END {print s}' $W/ref0.fa.fai) bp"
log done
