#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11
# md5 manifest of the 16 corrected FASTA files (g1_finals.tsv), reads only; small TSV to corrected_release/manifest.tsv
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp
printf "version\tsample\thaplotype\tpath\tbytes\tcontigs\tbp\tmd5\n" > corrected_release/manifest.tsv.part
tail -n +2 $T/g1_finals.tsv | while IFS=$'\t' read -r v sm dir; do for h in 1 2; do f=$D/$dir/$sm#$h.fa.gz
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $v $sm $h $f $(stat -c %s $f) $(wc -l < $f.fai) $(awk '{s += $2} END {printf "%.0f", s}' $f.fai) $(md5sum < $f | cut -c1-32) >> corrected_release/manifest.tsv.part
  echo "[$(date -Is)] $v $sm $h done"; done; done
mv corrected_release/manifest.tsv.part corrected_release/manifest.tsv; cat corrected_release/manifest.tsv
