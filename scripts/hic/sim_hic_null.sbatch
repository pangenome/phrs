#!/bin/bash
#SBATCH -p tux
#SBATCH -c 48
#SBATCH -t 1-00:00:00
#SBATCH --exclude=tux06,tux07,tux08,tux09
# Null Hi-C simulation (R1-M4): simulate cis-only read pairs near the ends of SAMPLE (sim_hic_null.py simulate), map
# them with the paper's HiC-Pro 3.0.0 run (same singularity image, config-hicpro.txt with N_CPU set to the allocated
# cores, same bowtie2 index and chrom.sizes), evaluate placement against the truth in the read names, and build the
# mcool exactly as for the real data (cooler cload pairs + zoomify --balance, cooler 0.10.4 env). All work in node
# /scratch; outputs copied to OUT.
# usage: sbatch sim_hic_null.sbatch [--sample HG002|CHM13] --domains TSV --fasta EXT_FA_GZ --sim PY --pairs-per-kb N --out DIR [--profile PROFILE_TSV_GZ]
#        [--read-len L] (default 250: HG002 Hi-C reads are 251 bp; CHM13 Hi-C reads are 150 bp)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
SM=HG002; RL=250
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --domains) DOM=$2;; --fasta) FA=$2;; --sim) SIM=$2;; --pairs-per-kb) PK=$2;; --read-len) RL=$2;; --profile) PF=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in DOM FA SIM PK O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
PF=${PF:-}
log () { echo "[$(date -Is)] $*"; }
H=/moosefs/guarracino/HPRCv2/PHR_III/HiC/$SM; SIF=/moosefs/guarracino/sifs/hicpro-latest.simg; E=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin
N=${SLURM_CPUS_PER_TASK:-48}
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simhic.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "sample $SM; host $(hostname); $N cores; domains $DOM; fasta $FA; $PK pairs per kb; read length $RL"
mkdir -p $W/in/hicpro_input/001 $W/out
$E/python3 $SIM simulate $DOM $FA $W/in/hicpro_input/001/001 --pairs-per-kb $PK --read-len $RL ${PF:+--profile $PF}
cp -r $H/bowtie2_index $W/in/; cp $H/chrom.sizes $W/in/
sed "s|^N_CPU = .*|N_CPU = $N|" $H/config-hicpro.txt > $W/in/config-hicpro.txt
nd=$(awk 'NR==FNR {a[FNR]=$0; next} a[FNR] != $0 {d++; print "config differs: " $0 > "/dev/stderr"} END {print d + 0}' $H/config-hicpro.txt $W/in/config-hicpro.txt)
[ "$nd" -eq 1 ] && grep -q "^N_CPU = $N$" $W/in/config-hicpro.txt || { echo "ERROR config differs beyond N_CPU ($nd lines)"; exit 1; }
log "HiC-Pro start"
singularity exec --bind $W/in:/input:ro --bind $W/out:/output $SIF /bin/bash -c 'export PATH=/usr/local/anaconda/envs/HiC-Pro_v3.0.0/bin:/HiC-Pro_3.0.0/bin:/HiC-Pro_3.0.0/scripts:$PATH; HiC-Pro -i /input/hicpro_input -o /output/hicpro_output -c /input/config-hicpro.txt'
AVP=$W/out/hicpro_output/hic_results/data/001/001.allValidPairs
log "HiC-Pro done: $(wc -l < $AVP) valid pairs"
$E/python3 $SIM evaluate $DOM $AVP $W/eval
awk -F'\t' 'BEGIN{OFS="\t"} {print $2,$3,$5,$6}' $AVP | $E/python3 -m cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 $W/in/chrom.sizes:5000 - $W/sim.5000.cool
$E/python3 -m cooler zoomify $W/sim.5000.cool -r 5000,10000,20000,50000,100000 --balance -o $W/sim.mcool
cp $W/eval.reads.tsv $W/eval.misassign.tsv $W/sim.mcool $O/; cp -r $W/out/hicpro_output/logs $O/hicpro_logs
cp $W/out/hicpro_output/hic_results/data/001/001.allValidPairs $O/sim.allValidPairs
cp -r $W/out/hicpro_output/hic_results/stats $O/hicpro_stats
log "done: outputs in $O"
