#!/bin/bash
#SBATCH -p lowmem,386mem,highmem,1tbmem
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 4:00:00
# DAD chr1p h1 start vs chr9q h2 end
# args: WASHU, OUT_DIR, [SPAN_BP] (default 300000)
# out: chr1_chr9.<v>.paf (ends, asm20, cs, all chains), chr1_chr9.whole.<v>.paf (whole chromosomes), lengths*.tsv
set -euo pipefail
D=$1; OUT=$2; SP=${3:-300000}; cd $D; mkdir -p $OUT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.syn19.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
# whole chromosomes: chr9 h2 vs chr1 h1, no cs
for v in v11 v131; do
  $S faidx norm_$v/PAN011#1.fa.gz "PAN011#1#chr1" > $W/chr1w.fa; $S faidx norm_$v/PAN011#2.fa.gz "PAN011#2#chr9" > $W/chr9w.fa
  $M -cx asm20 --secondary=no -t 16 $W/chr1w.fa $W/chr9w.fa 2> $W/mm2w.log > $OUT/chr1_chr9.whole.$v.paf
  printf "chr9\t%s\nchr1\t%s\n" $(grep -P "^PAN011#2#chr9\t" norm_$v/PAN011#2.fa.gz.fai | cut -f2) $(grep -P "^PAN011#1#chr1\t" norm_$v/PAN011#1.fa.gz.fai | cut -f2) > $OUT/lengths.whole.$v.tsv
  log "$v whole chains $(wc -l < $OUT/chr1_chr9.whole.$v.paf), >= 10 kb $(awk '$11 >= 10000' $OUT/chr1_chr9.whole.$v.paf | wc -l)"
done
for v in v11 v131; do
  L9=$(grep -P "^PAN011#2#chr9\t" norm_$v/PAN011#2.fa.gz.fai | cut -f2)
  $S faidx norm_$v/PAN011#1.fa.gz "PAN011#1#chr1:1-$SP" > $W/chr1.fa
  $S faidx norm_$v/PAN011#2.fa.gz "PAN011#2#chr9:$((L9 - SP + 1))-$L9" > $W/chr9.fa
  $M -cx asm20 --cs -P -t 4 $W/chr1.fa $W/chr9.fa 2> $W/mm2.log > $OUT/chr1_chr9.$v.paf
  printf "chr9_end\t%s\nchr9_window_start\t%s\nspan\t%s\n" $L9 $((L9 - SP + 1)) $SP > $OUT/lengths.$v.tsv
  log "$v chr9 h2 length $L9, chains $(wc -l < $OUT/chr1_chr9.$v.paf): $(awk '{printf "q%s-%s(%s) t%s-%s NM %s; ", $3,$4,$5,$8,$9,$13}' $OUT/chr1_chr9.$v.paf | cut -c1-600)"
done
log done
