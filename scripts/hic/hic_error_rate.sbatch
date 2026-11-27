#!/bin/bash
#SBATCH -c 16
#SBATCH -t 06:00:00
# Is the 0.2% per-base substitution rate of the Hi-C null simulation realistic? Map the first N real R1 reads of SAMPLE
# and N simulated reads (sim_hic_null.py simulate, same read length, --err ERR) with bowtie2 from the paper's HiC-Pro
# image and the paper's HiC-Pro global options (--very-sensitive-local -L 30 --score-min G,20,8 --local --reorder), on
# the paper's bowtie2 index. For primary mapped reads with MAPQ >= 30, and with any MAPQ: mismatches (XM) and gap extensions (XG) per aligned base (CIGAR M),
# soft-clipped bases per read base, and the fraction of reads with 0, 1, 2 and 3 or more mismatches. Writes OUT/rates.tsv.
# usage: sbatch hic_error_rate.sbatch --sample HG002|CHM13 --fastq R1.fastq.gz --read-len L --n N --err ERR --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --fastq) FQ=$2;; --read-len) RL=$2;; --n) NR=$2;; --err) ER=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SM FQ RL NR ER O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; H=/moosefs/guarracino/HPRCv2/PHR_III/HiC/$SM; SIF=/moosefs/guarracino/sifs/hicpro-latest.simg
E=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin; N=${SLURM_CPUS_PER_TASK:?}; s=$(echo $SM | tr A-Z a-z)
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/errrate.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); sample $SM; real $FQ; read length $RL; $NR reads; sim error $ER"
( set +o pipefail; zcat $FQ | head -n $((4 * NR)) ) > $W/real.fq
$E/python3 $Q/sim_hic_null.py simulate $Q/${s}_null.domains.tsv $Q/${s}_null.ext.fa.gz $W/sim --pairs-per-kb 2 --read-len $RL --err $ER
( set +o pipefail; zcat $W/sim_R1.fastq.gz | head -n $((4 * NR)) ) > $W/sim.fq
log "reads: real $(( $(wc -l < $W/real.fq) / 4 )), sim $(( $(wc -l < $W/sim.fq) / 4 ))"
log "real reads: mean error implied by the base qualities $($E/python3 -c "import sys; q=[c - 33 for i, l in enumerate(open(sys.argv[1], 'rb')) if i % 4 == 3 for c in l.rstrip()]; print(round(sum(10 ** (-v / 10) for v in q) / len(q), 6), 'over', len(q), 'bases')" $W/real.fq)"
printf "set\tmin_mapq\treads\taligned_bases\tmismatch_per_base\tgap_ext_per_base\tsoftclip_per_base\tfrac_0mm\tfrac_1mm\tfrac_2mm\tfrac_3plus\n" > $O/rates.tsv
for t in real sim; do
  singularity exec --bind $W:$W --bind /moosefs:/moosefs $SIF /usr/local/anaconda/envs/HiC-Pro_v3.0.0/bin/bowtie2 \
      --very-sensitive-local -L 30 --score-min G,20,8 --local --reorder -p $N -x $H/bowtie2_index/$SM -U $W/$t.fq 2> $W/$t.bt2.log \
    | singularity exec --bind $W:$W $SIF /usr/local/anaconda/envs/HiC-Pro_v3.0.0/bin/samtools view -F 0x904 -o $W/$t.sam -
  for q in 30 0; do awk -v t=$t -v q=$q '$5 >= q' $W/$t.sam \
    | awk -v t=$t -v q=$q 'BEGIN{OFS="\t"} {m=0; s=0; c=$6; while (match(c, /[0-9]+[MIDNSHP=X]/)) {n=substr(c, RSTART, RLENGTH-1); op=substr(c, RSTART+RLENGTH-1, 1); if (op=="M") m+=n; if (op=="S") s+=n; c=substr(c, RSTART+RLENGTH)}
        xm=0; xg=0; for (i=12; i<=NF; i++) {if ($i ~ /^XM:i:/) xm=substr($i, 6); if ($i ~ /^XG:i:/) xg=substr($i, 6)}
        r++; M+=m; S+=s; L+=length($10); X+=xm; G+=xg; k[(xm>=3)?3:xm]++}
        END{print t, q, r, M, X/M, G/M, S/L, k[0]/r, k[1]/r, k[2]/r, k[3]/r}' >> $O/rates.tsv; done
  log "$t: $(tail -1 $O/rates.tsv); bowtie2: $(grep overall $W/$t.bt2.log)"
done
cat $O/rates.tsv
