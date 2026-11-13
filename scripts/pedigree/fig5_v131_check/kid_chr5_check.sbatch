#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
# KID chr5q~chr1p tract: ONT read support
# args: WASHU, OUT_DIR
# Fig5E call: PAN028#1#chr5 v1.1 183,044,000-183,074,000 (v1.3.1 183,464,000-183,490,000), donor PAN027 h2 chr1.
# Own BAMs (PAN028 reads on PAN028 diploid, ont_full): per 2 kb bin from tract-40 kb to the contig end, depth MAPQ>=20,
# reads spanning tract +-5 kb at MAPQ>=20, reads clipped >= 1 kb with the clip in the bin.
# Cross BAMs (PAN028 reads on PAN027 diploid): for MAPQ>=20 primary reads over the donor PAN027#2#chr1 40,000-80,000,
# chromosomes of their SA alignments. out: kid_chr5.<v>.bins.tsv, kid_chr5.<v>.summary.txt
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p $OUT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
log () { echo "[$(date -Is)] $*"; }
for v in v11 v131; do
  if [ $v = v11 ]; then T0=183044000; T1=183074000; else T0=183464000; T1=183490000; fi
  ls ont_full/SRR*.PAN028.$v.bam > $OUT/own.$v.txt; ls ont_full/SRR*.PAN028.onPAN027.$v.bam > $OUT/cross.$v.txt
  b1=$(head -1 $OUT/own.$v.txt); C=$($S view -H $b1 | awk -F'\t' '$1=="@SQ" {sub("SN:","",$2); print $2}' | grep -E "^PAN028#1#chr5$")
  L=$($S view -H $b1 | awk -F'\t' -v c=$C '$1=="@SQ" && $2=="SN:"c {sub("LN:","",$3); print $3}')
  log "$v: own BAMs $(wc -l < $OUT/own.$v.txt), cross BAMs $(wc -l < $OUT/cross.$v.txt), contig $C length $L, tract $T0-$T1"
  R=$C:$((T0-40000))-$L
  spans=$(for b in $(cat $OUT/own.$v.txt); do $S view -q 20 -F 0x904 $b $C:$T0-$T1 | awk -v a=$((T0-5000)) -v z=$((T1+5000)) '{e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; if ($4<=a && e>=z) n++} END {print n+0}'; done | awk '{s+=$1} END {print s}')
  log "$v: MAPQ>=20 primary reads spanning tract +-5 kb: $spans"
  for b in $(cat $OUT/own.$v.txt); do $S view -q 20 -F 0x904 $b $R; done | awk -v s=$((T0-40000)) -v OFS='\t' '
    {bin=int(($4-s)/2000); if (bin<0) bin=0; c=$6; if (c ~ /^[0-9]+[SH]/ && substr(c,1,match(c,/[SH]/)-1)+0>=1000) cl[bin]++;
     if (c ~ /[0-9]+[SH]$/) {t=c; sub(/.*[MIDNX=]/,"",t); if (t+0>=1000) {e=$4; x=$6; while (match(x,/[0-9]+[MDN=X]/)) {e+=substr(x,RSTART,RLENGTH-1); x=substr(x,RSTART+RLENGTH)}; cl[int((e-s)/2000)]++}}}
    END {for (k in cl) print s+k*2000, cl[k]}' | sort -k1,1n > $OUT/clip.$v.tsv
  $S depth -a -Q 20 -r $R $(cat $OUT/own.$v.txt) | awk -v s=$((T0-40000)) -v OFS='\t' '{d=0; for (i=3;i<=NF;i++) d+=$i; b=s+int(($2-s)/2000)*2000; sum[b]+=d; n[b]++} END {for (b in sum) print b, sum[b]/n[b]}' | sort -k1,1n > $OUT/depth.$v.tsv
  join -a1 -e0 -o 0,1.2,2.2 -t $'\t' <(sort -k1,1 $OUT/depth.$v.tsv) <(sort -k1,1 $OUT/clip.$v.tsv) | sort -k1,1n | awk -v a=$T0 -v z=$T1 -v OFS='\t' 'BEGIN {print "bin_start","mean_depth_q20","clipped_ge_1kb","in_tract"} {print $1, sprintf("%.1f",$2), $3, ($1>=a && $1<z)?"yes":"no"}' > $OUT/kid_chr5.$v.bins.tsv
  log "$v bins: $(wc -l < $OUT/kid_chr5.$v.bins.tsv); tract mean depth $(awk -F'\t' '$4=="yes" {s+=$2; n++} END {printf "%.1f", s/n}' $OUT/kid_chr5.$v.bins.tsv), flanks 40 kb before $(awk -F'\t' 'NR>1 && $4=="no" && $1<'$T0' {s+=$2; n++} END {printf "%.1f", s/n}' $OUT/kid_chr5.$v.bins.tsv), clipped reads in tract bins $(awk -F'\t' '$4=="yes" {s+=$3} END {print s+0}' $OUT/kid_chr5.$v.bins.tsv)"
  for b in $(cat $OUT/cross.$v.txt); do $S view -q 20 -F 0x904 $b PAN027#2#chr1:40000-80000; done | awk -F'\t' '{sa="none"; for (i=12;i<=NF;i++) if ($i ~ /^SA:Z:/) {sa=$i}; if (sa=="none") {print "no_SA"} else {gsub(/^SA:Z:/,"",sa); n=split(sa,a,";"); for (j=1;j<n;j++) {split(a[j],f,","); x=f[1]; sub(/.*#/,"",x); print x}}}' | sort | uniq -c | sort -k1,1nr > $OUT/cross_donor_sa.$v.txt
  log "$v cross: MAPQ>=20 PAN028 reads on PAN027#2#chr1 40-80 kb, SA chromosomes: $(head -8 $OUT/cross_donor_sa.$v.txt | tr -s ' ' | tr '\n' ';')"
done
log done
