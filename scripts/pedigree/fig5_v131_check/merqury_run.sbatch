#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05,octopus06
# merqury QV and error k-mers, one diploid
# args: WASHU, SAMPLE, NORM_DIR (e.g. norm_v11), TAG (e.g. v11)
# merqury v1.4.1 (scripts_tmp/merqury-1.4.1), meryl 1.4.1 k=21 from short/meryl/<SAMPLE>.k21.meryl (Illumina, own reads).
# out: WASHU/merqury/<TAG>/<SAMPLE>/ (<SAMPLE>.qv, <SAMPLE>.<hap>.qv, *_only.bed = assembly k-mers absent from reads)
set -euo pipefail
D=$1; SM=$2; N=$3; TAG=$4; cd $D; O=$D/merqury/$TAG/$SM; mkdir -p $O
G=/gnu/store; export MERQURY=$D/scripts_tmp/merqury-1.4.1
export PATH=$G/8zwv41mb9sgb08652jfn9s4h697550v9-meryl-1.4.1/bin:$G/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin:$G/18fiwlwr3d3jglkcypdv5250d098ax12-bedtools-2.31.1/bin:$PATH
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/merqury.XXXX); trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
cp -r short/meryl/$SM.k21.meryl $W/reads.meryl
zcat $N/$SM#1.fa.gz > $W/${SM}_hap1.fasta; zcat $N/$SM#2.fa.gz > $W/${SM}_hap2.fasta
log "$SM $TAG host $(hostname), meryl $(meryl --version 2>&1 | head -1), bedtools $(bedtools --version), Rscript $(which Rscript 2>/dev/null || echo none)"
cd $W; export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
bash $MERQURY/merqury.sh reads.meryl ${SM}_hap1.fasta ${SM}_hap2.fasta $SM > merqury.log 2>&1 || log "WARN merqury exit $? (plots need R; QV files checked below)"
ls $SM*.qv >/dev/null || { tail -30 merqury.log; exit 1; }
cp $SM*.qv $SM*completeness.stats merqury.log $O/ 2>/dev/null || true
for b in *_only.bed; do [ -s $b ] && gzip -c $b > $O/$b.gz; done
for f in $SM*.qv; do log "$f: $(tr '\n' ';' < $f)"; done
log "error k-mer positions: $(for b in *_only.bed; do echo "$b $(wc -l < $b)"; done | tr '\n' ';')"
