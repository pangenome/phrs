#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06
# Reads near contig ends, aligned to a sim assembly
# usage: sbatch g1_sim_reads.sbatch --washu D --sim SIM_DIR --sample SM --version V --window BP --out DIR
# reads: primary HiFi (hifi/bam) and ONT (ont_full) reads of norm_V in the first and last WINDOW bp of every contig, both
#   haplotypes; reads from elsewhere in the genome that would align to an end are not included.
# aligned to the sim diploid (SIM_DIR/mmi/SM.{hifi,ont}.mmi from g1_mmi.sbatch; minimap2 2.31 -a map-hifi / map-ont), sorted, indexed.
# out: DIR/SM.hifi.bam, DIR/SM.ont.bam (+ .bai)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sim) SD=$2;; --sample) SM=$2;; --version) V=$2;; --window) WN=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SD SM V WN O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2; T=${SLURM_CPUS_PER_TASK:-4}
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simreads.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
cat norm_$V/$SM#1.fa.gz.fai norm_$V/$SM#2.fa.gz.fai | awk -F'\t' -v OFS='\t' -v w=$WN '{print $1, 0, ($2 < w ? $2 : w); print $1, ($2 > w ? $2 - w : 0), $2}' | sort -k1,1 -k2,2n > $W/ends.bed
log "$SM $V: $(wc -l < $W/ends.bed) end windows of $WN bp on norm_$V; sim $SD; host $(hostname)"
for s in hifi ont; do
  [ -s $O/$SM.$s.bam.bai ] && { log "exists $O/$SM.$s.bam, skip"; continue; }
  [ -s $SD/mmi/$SM.$s.mmi ] || { echo "ERROR missing $SD/mmi/$SM.$s.mmi"; exit 1; }
  [ $s = hifi ] && BS=$(ls hifi/bam/*.$SM.$V.bam) || BS=$(ls ont_full/SRR*.$SM.$V.bam)
  for b in $BS; do
    for t in 1 2 3; do  # moosefs read errors: retry each BAM
      if $S view -F 0x900 -L $W/ends.bed $b > $W/view.sam; then break; fi
      [ $t -lt 3 ] || { echo "ERROR samtools view failed 3 times on $b" >&2; exit 1; }
      log "samtools view failed on $b (try $t of 3), retry in 60 s" >&2; sleep 60
    done
    cat $W/view.sam
  done | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/$s.fq
  log "reads $s: $(awk 'NR%4==2' $W/$s.fq | wc -l) reads, $(awk 'NR%4==2 {n+=length($0)} END {print n+0}' $W/$s.fq) bp from $(echo $BS | wc -w) BAMs"
  $MM -ax map-$s -t $T $SD/mmi/$SM.$s.mmi $W/$s.fq 2> $W/$s.mm.log | $S sort -@ 4 -m 2G -T $W/$s.sort -o $W/$s.bam -
  tail -3 $W/$s.mm.log; $S index $W/$s.bam; log "aligned $s: $($S view -c -F 0x904 $W/$s.bam) primary mapped, $($S view -c -f 4 $W/$s.bam) unmapped"
  cp $W/$s.bam $O/$SM.$s.bam.part && mv $O/$SM.$s.bam.part $O/$SM.$s.bam && cp $W/$s.bam.bai $O/$SM.$s.bam.bai
  rm -f $W/$s.fq $W/$s.bam $W/$s.bam.bai $W/view.sam
done
log done
