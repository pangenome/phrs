#!/bin/bash
#SBATCH -c 16
#SBATCH -t 06:00:00
# Are the error rates of the Pore-C / CiFi null simulations realistic (sim_porec_null.py: Pore-C 0.5% substitutions,
# 0.25% insertions, 0.25% deletions; CiFi a tenth of that)? Digest and map the first N real reads and the simulated reads
# with the paper's chain (pore-c-py digest ENZYME --max_monomers 250 | samtools fastq -T '*' | minimap2 -ay -x PRESET
# --cap-kalloc 100m --cap-sw-mem 50m HG002.fa.gz). For primary monomer alignments with MAPQ >= 30, and with any MAPQ:
# gap-compressed divergence (minimap2 de tag) weighted by aligned length, edit distance (NM) per aligned base, and
# soft-clipped bases per monomer base; plus the fraction of monomers mapped. Writes OUT/rates.tsv.
# usage: sbatch porec_error_rate.sbatch --assay porec|cifi --n N --reads-per-kb R --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --assay) AS=$2;; --n) NR=$2;; --reads-per-kb) RK=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in AS NR RK O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
export PATH=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin:$PATH
B=/moosefs/guarracino/HPRCv2; H=$B/PHR_III/HiC/HG002; Q=$B/PHR_III/simhic; N=${SLURM_CPUS_PER_TASK:?}
case $AS in
  porec) EZ=NlaIII; PR=map-ont; SITE=CATG; ER="0.005 0.0025 0.0025"; IN=$B/data/PoreC/HG002/HG002_1_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz;;
  cifi)  EZ=DpnII; PR=map-hifi; SITE=GATC; ER="0.0005 0.00025 0.00025"; IN=$B/data/CiFi/HG002/HG002_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam;;
  *) echo "ERROR assay $AS"; exit 2;;
esac
read SUB INS DEL <<< "$ER"
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/perr.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); assay $AS ($EZ, $PR); first $NR real reads of $IN; sim errors $SUB $INS $DEL at $RK reads per kb"
# the reader stops early, so the upstream gets SIGPIPE: pipefail is off inside these subshells
case $IN in *.bam) ( set +o pipefail; samtools view -h $IN | awk -v n=$NR '/^@/ {print; next} {if (++k > n) exit; print}' ) | samtools view -b -o $W/real.bam -; RIN=$W/real.bam;;
            *) ( set +o pipefail; zcat $IN | head -n $((4 * NR)) ) | gzip -1 > $W/real.fastq.gz; RIN=$W/real.fastq.gz;; esac
python3 $Q/sim_porec_null.py simulate $Q/hg002_null.domains.tsv $Q/hg002_null.ext.fa.gz $H/fragments_$EZ.bed $SITE $W/sim --reads-per-kb $RK --sub $SUB --ins $INS --dele $DEL
printf "set\tmin_mapq\tmonomers\tmapped\tfrac_mapped\taligned_bases\tdivergence\tnm_per_base\tsoftclip_per_base\n" > $O/rates.tsv
for t in real sim; do
  src=$RIN; [ $t = sim ] && src=$W/sim.fastq.gz
  pore-c-py digest $src $EZ --threads 4 --max_monomers 250 | samtools fastq -T '*' > $W/$t.mono.fq
  minimap2 -ay -t $N -x $PR --cap-kalloc 100m --cap-sw-mem 50m $H/HG002.fa.gz $W/$t.mono.fq 2> $W/$t.mm2.log | samtools view -F 0x900 -o $W/$t.sam -
  for q in 30 0; do awk -v t=$t -v q=$q 'BEGIN{OFS="\t"} {n++; if ($3 == "*" || $5 < q) next; r++; m=0; s=0; c=$6
      while (match(c, /[0-9]+[MIDNSHP=X]/)) {x=substr(c, RSTART, RLENGTH-1); op=substr(c, RSTART+RLENGTH-1, 1); if (op=="M") m+=x; if (op=="S") s+=x; c=substr(c, RSTART+RLENGTH)}
      for (i=12; i<=NF; i++) {if ($i ~ /^NM:i:/) nm=substr($i, 6); if ($i ~ /^de:f:/) de=substr($i, 6)}
      M+=m; S+=s; L+=length($10); NMs+=nm; D+=de*m}
      END{print t, q, n, r, r/n, M, D/M, NMs/M, S/L}' $W/$t.sam >> $O/rates.tsv; done
  log "$t: $(grep -c '' $W/$t.mono.fq | awk '{print $1/4}') monomers; $(tail -2 $O/rates.tsv | tr '\n' ' ')"
done
column -t $O/rates.tsv
