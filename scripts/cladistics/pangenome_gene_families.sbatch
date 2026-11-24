#!/bin/bash
#SBATCH -J pgfam
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH -p workers
#SBATCH --exclude=octopus05,octopus06
# Gene-window overlaps from the per-haplotype Liftoff GFF3s, one array task per chunk (pangenome_gene_families.py extract)
# Usage: sbatch --array=1-K pangenome_gene_families.sbatch GFF_LIST K OUTDIR
set -euxo pipefail
L=$1; K=$2; O=$3; N=$SLURM_ARRAY_TASK_ID
S=/moosefs/guarracino/HPRCv2/PHR_III/conversion/scripts
mkdir -p $O; hostname; date '+%F %T %Z'
mapfile -t G < <(awk -v n=$N -v k=$K '(NR - n) % k == 0' $L)
echo "${#G[@]} GFF3 files"
python3 $S/pangenome_gene_families.py extract --fai /moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai \
    --gff "${G[@]}" --out $O/chunk_$N.tsv.part
mv $O/chunk_$N.tsv.part $O/chunk_$N.tsv
date '+%F %T %Z'
