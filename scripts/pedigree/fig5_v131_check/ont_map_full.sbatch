#!/bin/bash
#SBATCH -p workers
#SBATCH -c 48
#SBATCH --mem=240G
#SBATCH -t 10-00:00:00
#SBATCH --exclude=octopus05,octopus06,octopus07,octopus08,octopus09
# One ONT run on whole diploids, no region filter
# args: WASHU, SRR, SAMPLE, MD5, BYTES, URL, SRA_DIR, OUT_DIR, "REF1 REF2 .." (default SAMPLE), ["v11 v131"]
# One BAM per (REF, version): OUT_DIR/<SRR>.<SAMPLE>[.on<REF>].<v>.bam, all reads, primary + supplementary.
# Alignments run one after the other; a BAM that is already in OUT_DIR with its .bai is skipped.
set -euo pipefail
WASHU=$1; SRR=$2; SAMPLE=$3; MD5=$4; BYTES=$5; URL=$6; FQDIR=$7; OUT=$8; REFS=${9:-$3}; VERS=${10:-v11 v131}
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
SRA=/gnu/store/sjzg68akxlirqmm9vmmhrfvxnvdx7fa7-sra-tools-2.10.6/bin
T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.ontfull.$SLURM_JOB_ID; mkdir -p $W $OUT; trap 'rm -rf $W' EXIT
FQ=$W/${SRR}.sra
log () { echo "[$(date -Is)] $*"; }
log "run $SRR ($SAMPLE) on [$REFS] x [$VERS], host $(hostname), scratch free $(df -h --output=avail $W | tail -1)"

# list missing outputs
TODO=()
for REF in $REFS; do for v in $VERS; do
  TAG=$SAMPLE; [ "$REF" = "$SAMPLE" ] || TAG=$SAMPLE.on$REF
  [ -s $OUT/${SRR}.$TAG.$v.bam.bai ] && [ -s $OUT/${SRR}.$TAG.$v.bam ] && { log "skip ${SRR}.$TAG.$v.bam (present)"; continue; }
  TODO+=("$REF:$v:$TAG")
done; done
[ ${#TODO[@]} -gt 0 ] || { log "nothing to do"; exit 0; }

# kept .sra copy, 40 MB/s, md5
[ -s $FQDIR/${SRR}.sra ] && [ "$(stat -c %s $FQDIR/${SRR}.sra)" = "$BYTES" ] || { log "ERROR: no kept $FQDIR/${SRR}.sra of $BYTES bytes ($URL)"; exit 1; }
rsync --bwlimit=40000 $FQDIR/${SRR}.sra $FQ
GOT=$(md5sum $FQ | cut -d' ' -f1); [ "$GOT" = "$MD5" ] || { log "ERROR md5 $GOT != $MD5"; exit 1; }
log "sra copied, md5 ok"

for x in "${TODO[@]}"; do
  REF=${x%%:*}; v=$(echo $x | cut -d: -f2); TAG=${x##*:}; B=$W/${SRR}.$TAG.$v.bam
  for try in 1 2 3 4 5; do
    zcat "$WASHU/norm_$v/$REF#1.fa.gz" "$WASHU/norm_$v/$REF#2.fa.gz" > $W/ref.fa && break
    log "WARN ref $REF $v read failed, try $try"; sleep 120; [ $try = 5 ] && { log "ERROR ref"; exit 1; }
  done
  $MM2 -x map-ont -I 8G -t 24 -d $W/ref.mmi $W/ref.fa 2> $W/idx.log; rm $W/ref.fa
  log "$TAG $v: index $(stat -c %s $W/ref.mmi) bytes"
  $SRA/fastq-dump -Z $FQ 2> $W/dump.log \
  | $MM2 -ax map-ont -Y --secondary=no -K 4G -t $((T - 10)) $W/ref.mmi - 2> $W/map.log \
  | $SAMTOOLS sort -@ 8 -m 3G -T $W/sort -o $B -
  $SAMTOOLS index -@ 8 $B; $SAMTOOLS quickcheck $B
  grep -q "^@SQ" <($SAMTOOLS view -H $B) || { log "ERROR: no @SQ in $B"; exit 1; }
  log "$TAG $v: records $($SAMTOOLS view -@ 8 -c $B), primary $($SAMTOOLS view -@ 8 -c -F 0x904 $B), unmapped $($SAMTOOLS view -@ 8 -c -f 4 $B), $(stat -c %s $B) bytes"
  tail -n 2 $W/map.log; tail -n 3 $W/dump.log
  rsync --bwlimit=40000 $B $OUT/$(basename $B).part.$SLURM_JOB_ID && mv $OUT/$(basename $B).part.$SLURM_JOB_ID $OUT/$(basename $B)
  [ "$(stat -c %s $OUT/$(basename $B))" = "$(stat -c %s $B)" ] || { log "ERROR copy size $B"; exit 1; }
  cp $B.bai $OUT/; rm -f $B $B.bai $W/ref.mmi
  log "$TAG $v: kept $OUT/$(basename $B)"
done
log done
