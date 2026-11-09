#!/bin/bash
#SBATCH -p workers
#SBATCH -c 48
#SBATCH --mem=240G
#SBATCH -t 3-00:00:00
#SBATCH --exclude=octopus05,octopus06,octopus08,octopus09
# Download one ONT .sra, keep it, map
# args: WASHU, SRR, SAMPLE, MD5, BYTES, URL, SRA_DIR, OUT_DIR, [REF]
# REF: assembly to map on (default SAMPLE)
# OUT_DIR/beds/REF.{v11,v131}.bed must exist
set -euo pipefail
WASHU=$1; SRR=$2; SAMPLE=$3; MD5=$4; BYTES=$5; URL=$6; FQDIR=$7; OUT=$8; REF=${9:-$3}
TAG=$SAMPLE; [ "$REF" = "$SAMPLE" ] || TAG=$SAMPLE.on$REF
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
SRA=/gnu/store/sjzg68akxlirqmm9vmmhrfvxnvdx7fa7-sra-tools-2.10.6/bin
W=/scratch/$USER.ont.$SLURM_JOB_ID; mkdir -p $W $FQDIR $OUT; trap 'rm -rf $W' EXIT
FQ=$W/${SRR}.sra
log () { echo "[$(date -Is)] $*"; }
log "reads $SAMPLE on $REF, host $(hostname) scratch free $(df -h --output=avail $W | tail -1) url $URL"

# Diploid refs first, retry moosefs reads
for v in v11 v131; do
  for try in 1 2 3 4 5; do
    zcat "$WASHU/norm_$v/$REF#1.fa.gz" "$WASHU/norm_$v/$REF#2.fa.gz" > $W/ref_$v.fa && cp $OUT/beds/$REF.$v.bed $W/regions_$v.bed && break
    log "WARN ref $v read failed, try $try"; sleep 120
    [ $try = 5 ] && { log "ERROR ref $v"; exit 1; }
  done
  $MM2 -x map-ont -I 8G -t 24 -d $W/ref_$v.mmi $W/ref_$v.fa 2> $W/idx_$v.log
  log "index $v $(stat -c %s $W/ref_$v.mmi) bytes, regions $(wc -l < $W/regions_$v.bed)"
done

# Reuse kept copy, else download
if [ -s $FQDIR/${SRR}.sra ] && [ "$(stat -c %s $FQDIR/${SRR}.sra)" = "$BYTES" ]; then
  rsync --bwlimit=40000 $FQDIR/${SRR}.sra $FQ; log "reused $FQDIR/${SRR}.sra"
fi
python3 - $URL $FQ $BYTES <<'PY'
import os, sys, time, urllib.request
url, path, size = sys.argv[1], sys.argv[2], int(sys.argv[3]); tries = 0
while True:
    done = os.path.getsize(path) if os.path.exists(path) else 0
    if done >= size: break
    try:
        req = urllib.request.Request(url, headers={"Range": f"bytes={done}-"} if done else {})
        with urllib.request.urlopen(req, timeout=120) as r, open(path, 'ab') as o:
            while True:
                b = r.read(1 << 22)
                if not b: break
                o.write(b)
    except Exception as e:
        tries += 1
        sys.stderr.write(f"retry {tries} at {done} bytes: {e!r}\n")
        if tries > 100: raise
        time.sleep(30)
sys.stderr.write(f"downloaded {os.path.getsize(path)} of {size} bytes\n")
PY
GOT=$(md5sum $FQ | cut -d' ' -f1)
[ "$GOT" = "$MD5" ] || { log "ERROR md5 $GOT != $MD5"; exit 1; }
log "md5 ok $GOT"

# Map to both refs, keep call regions
map_one () { local v=$1
  $SRA/fastq-dump -Z $FQ 2> $W/dump_$v.log \
  | $MM2 -ax map-ont -Y --secondary=no -K 4G -t $((SLURM_CPUS_PER_TASK / 2 - 1)) $W/ref_$v.mmi - 2> $W/map_$v.log \
  | $SAMTOOLS view -h -L $W/regions_$v.bed -u - \
  | $SAMTOOLS sort -@ 2 -m 2G -T $W/sort_$v -o $W/${SRR}.$TAG.$v.bam -
}
map_one v11 & P11=$!
map_one v131
wait $P11
for v in v11 v131; do
  B=$W/${SRR}.$TAG.$v.bam; $SAMTOOLS index $B
  log "$v records $($SAMTOOLS view -c $B) primary $($SAMTOOLS view -c -F 0x904 $B)"
  tail -n 2 $W/map_$v.log; tail -n 3 $W/dump_$v.log
  cp $B $B.bai $OUT/
done
# Keep the .sra on moosefs, 40 MB/s
if [ ! -s $FQDIR/${SRR}.sra ]; then
  rsync --bwlimit=40000 $FQ $FQDIR/${SRR}.sra.part.$SLURM_JOB_ID && mv -n $FQDIR/${SRR}.sra.part.$SLURM_JOB_ID $FQDIR/${SRR}.sra; rm -f $FQDIR/${SRR}.sra.part.$SLURM_JOB_ID
fi
[ "$(stat -c %s $FQDIR/${SRR}.sra)" = "$BYTES" ] || { log "ERROR copy size"; exit 1; }
log "kept $FQDIR/${SRR}.sra"
log "done"
