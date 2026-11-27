#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 24
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06
# Null Pore-C / CiFi simulation (R1-M4): simulate cis-only concatemers near the HG002 ends (sim_porec_null.py), then
# the paper's pipeline (code_PHR-and-3D.md, HG002 Pore-C and CiFi): pore-c-py digest ENZYME --max_monomers 250 |
# samtools fastq -T '*' | minimap2 -ay -x PRESET --cap-kalloc 100m --cap-sw-mem 50m HG002.fa.gz | pore-c-py annotate
# --monomers --stdout | pairtools parse2 --single-end --readid-transform 'readID.split(":")[0]' --drop-sam --drop-seq
# --expand | pairtools restrict -f FRAGMENTS | pairtools sort | pairtools dedup | cooler cload pairs (5 kb) + zoomify
# --balance. parse2 runs twice on the same annotated alignments: with its default --min-mapq 1 (as in the paper) and
# with --min-mapq 0. Tools from the rebuilt paper env (condatools/porec/1.0-rebuilt). Work in node /scratch.
# usage: sbatch sim_porec_null.sbatch --enzyme NlaIII|DpnII --preset map-ont|map-hifi --site CATG|GATC --domains TSV [--profile PROFILE_TSV_GZ]
#        --fasta EXT_FA_GZ --sim PY --reads-per-kb N --errors "SUB INS DEL" --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --enzyme) EZ=$2;; --preset) PR=$2;; --site) SITE=$2;; --domains) DOM=$2;; --fasta) FA=$2;; --sim) SIM=$2;; --reads-per-kb) RK=$2;; --errors) ER=$2;; --out) O=$2;; --profile) PF=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in EZ PR SITE DOM FA SIM RK ER O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
PF=${PF:-}
log () { echo "[$(date -Is)] $*"; }
export PATH=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin:$PATH
H=/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002; REF=$H/HG002.fa.gz; CS=$H/chrom.sizes; FR=$H/fragments_$EZ.bed
N=${SLURM_CPUS_PER_TASK:-48}; read SUB INS DEL <<< "$ER"
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simporec.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); $N cores; $EZ $PR site $SITE; errors sub $SUB ins $INS del $DEL; $RK reads per kb; $(pore-c-py --version 2>&1 | tail -1); $(minimap2 --version); pairtools $(pairtools --version 2>&1 | awk '{print $NF}'); $(samtools --version | head -1)"
python3 $SIM simulate $DOM $FA $FR $SITE $W/sim --reads-per-kb $RK --sub $SUB --ins $INS --dele $DEL ${PF:+--profile $PF}
pore-c-py digest $W/sim.fastq.gz $EZ --threads $((N / 3)) --max_monomers 250 \
  | samtools fastq -T '*' \
  | minimap2 -ay -t $((N * 2 / 3)) -x $PR --cap-kalloc 100m --cap-sw-mem 50m $REF - \
  | pore-c-py annotate - $W/annotated --monomers --stdout > $W/annotated.stdout
log "annotated alignments: $(stat -c %s $W/annotated.stdout) bytes; first bytes $(head -c 4 $W/annotated.stdout | od -c | head -1)"
for mq in 1 0; do
  pairtools parse2 --single-end --readid-transform 'readID.split(":")[0]' --drop-sam --drop-seq --expand --min-mapq $mq \
      --output-stats $W/parse2_stats.mq$mq.txt -c $CS $W/annotated.stdout \
    | pairtools restrict -f $FR \
    | pairtools sort --nproc $((N / 3)) --tmpdir $W \
    | pairtools dedup --output $W/dedup.mq$mq.pairs
  log "min-mapq $mq: $(grep -c '^[^#]' $W/dedup.mq$mq.pairs) pairs after dedup; $(grep -E '^(total_mapped|pair_types/UU|pair_types/MM)' $W/parse2_stats.mq$mq.txt | tr '\n' ' ')"
  python3 $SIM evaluate $DOM $W/sim.truth.tsv $W/dedup.mq$mq.pairs $W/eval.mq$mq
  cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $CS:5000 $W/dedup.mq$mq.pairs $W/sim.mq$mq.5000.cool
  cooler zoomify $W/sim.mq$mq.5000.cool -r 5000,10000,20000,50000,100000 --balance -o $W/sim.mq$mq.mcool
  cp $W/sim.mq$mq.mcool $W/eval.mq$mq.misassign.tsv $W/parse2_stats.mq$mq.txt $O/
done
cp $W/sim.truth.tsv $O/
log "done: outputs in $O"
