#!/bin/bash
#SBATCH -J plant
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -p workers
#SBATCH --exclude=octopus05,octopus06
# Known answer on real data for the conversion scan (R1-M1 test 3b): plant tracts into one real window
# (plant_tracts.py), realign the modified window and, as control, the unmodified one against all 18,827
# windows with the paper's wfmash call (v0.23.0-41, -p 95), then scan both at 1, 2, 5 and 10 kb bins.
# Usage: sbatch plant_run.sbatch RECIPIENT_ARM DONOR_ARM INDEX OUTDIR
set -euxo pipefail
R=$1; D=$2; I=$3; O=$4
B=/moosefs/guarracino/HPRCv2/PHR_III; S=$B/conversion/scripts
FA=$B/hprcv2.1Mb.telo_500kb_trimmed.fa.gz
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash
ST=/gnu/store/02z0lccy5klz70a7z2d5fdaaljp3w3jz-samtools-1.19/bin/samtools
N=${R}_from_${D}_$I
hostname; date '+%F %T %Z'; $WF --version 2>&1 | head -1
T=/scratch/$USER/tmp/phrs/plant.$SLURM_JOB_ID; mkdir -p $T $O; trap 'rm -rf $T' EXIT
python3 $S/plant_tracts.py --list $B/conversion/all/list.txt --fasta $FA --recipient-arm $R --donor-arm $D \
    --index $I --samtools $ST --out $T/$N
NAME=$(head -1 $T/$N.fa | cut -c2-)
$ST faidx $FA "$NAME" > $T/$N.orig.fa
$ST faidx $T/$N.fa; $ST faidx $T/$N.orig.fa  # wfmash needs the target .fai
for k in plant orig; do
    F=$T/$N.fa; [ $k = orig ] && F=$T/$N.orig.fa
    t0=$(date +%s)
    $WF -p 95 -t $SLURM_CPUS_PER_TASK $F $FA --quiet | gzip > $T/$N.$k.paf.gz
    echo "wfmash $k: $(( $(date +%s) - t0 )) s, $(zcat $T/$N.$k.paf.gz | wc -l) PAF lines"
    for b in 1000 2000 5000 10000; do
        python3 $S/conversion_scan.py $T/$N.$k.paf.gz --bin $b --out $T/$N.$k.b$b
    done
done
cp $T/$N.fa $T/$N.plants.tsv $T/$N.*.paf.gz $T/$N.*.b*.tsv $O/
date '+%F %T %Z'
