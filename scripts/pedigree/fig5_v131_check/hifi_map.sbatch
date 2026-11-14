#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# HiFi run to own diploid, with sample check
# args: WASHU, RUNS_TSV (SRR, sample, md5, bytes, S3 URL), VERSIONS (e.g. "v11 v131") ; sbatch --array=1-N
# 1. S3 download (resume), md5, .sra kept in WASHU/hifi/sra/.
# 2. Sample check: 20,000 reads (minimap2 -x map-hifi, primary, MAPQ >= 30) on the PAN027 and PAN028 v1.1 diploids
#    and on the expected sample; exact = NM 0. Stop when another sample has more exact reads than the expected one.
# 3. All reads (fastq-dump) to the own diploid of each version: minimap2 2.31 -ax map-hifi --eqx -Y, sorted BAM.
# out: WASHU/hifi/bam/<SRR>.<sample>.<v>.bam, WASHU/hifi/check/<SRR>.who.tsv
set -euo pipefail
D=$1; T=$2; VS=$3; cd $D; mkdir -p hifi/{sra,bam,check}
SRA=/gnu/store/sjzg68akxlirqmm9vmmhrfvxnvdx7fa7-sra-tools-2.10.6/bin; MM=$D/scripts_tmp/minimap2-2.31/minimap2
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; C=$SLURM_CPUS_PER_TASK
read SRR SM MD5 BYTES URL < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $T)
W=/scratch/$USER.hifi.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
log "$SRR expected $SM, host $(hostname), scratch free $(df -h --output=avail /scratch | tail -1)"
K=hifi/sra/$SRR.sra; F=$W/$SRR.sra
if [ -s $K ] && [ "$(stat -c %s $K)" = "$BYTES" ]; then cp $K $F; log "reused $K"; else
python3 - $URL $F $BYTES <<'PY'
import os, sys, time, urllib.request
url, path, size = sys.argv[1], sys.argv[2], int(sys.argv[3]); tries = 0
while True:
    done = os.path.getsize(path) if os.path.exists(path) else 0
    if done >= size: break
    try:
        req = urllib.request.Request(url, headers={"Range": f"bytes={done}-"} if done else {})
        with urllib.request.urlopen(req, timeout=120) as r, open(path, "ab") as o:
            while True:
                b = r.read(1 << 22)
                if not b: break
                o.write(b)
    except Exception as e:
        tries += 1; sys.stderr.write(f"retry {tries} at {done} bytes: {e!r}\n")
        if tries > 100: raise
        time.sleep(30)
PY
fi
GOT=$(md5sum $F | cut -d' ' -f1); [ "$GOT" = "$MD5" ] || { log "ERROR md5 $GOT != $MD5"; exit 1; }
[ -s $K ] || cp $F $K; log "md5 ok, kept $K"
$SRA/fastq-dump -X 20000 -Z $F 2> $W/dump_test.log > $W/test.fq
log "test reads $(awk 'NR%4==2' $W/test.fq | wc -l)"
idx () { [ -s $W/$1.$2.mmi ] || { zcat norm_$2/$1#1.fa.gz norm_$2/$1#2.fa.gz > $W/ref.fa; $MM -x map-hifi -I 8G -t $C -d $W/$1.$2.mmi $W/ref.fa 2> $W/idx.log; rm -f $W/ref.fa; }; }
printf "srr\texpected\tsample\treads_q30\texact\tfrac_exact\n" > hifi/check/$SRR.who.tsv
for c in $(echo PAN027 PAN028 $SM | tr ' ' '\n' | sort -u); do
  idx $c v11
  $MM -ax map-hifi --secondary=no -t $C $W/$c.v11.mmi $W/test.fq 2> /dev/null | $S view -F 0x904 -q 30 - \
    | awk -v s=$SRR -v e=$SM -v c=$c -v OFS='\t' '{n++; for (i=12;i<=NF;i++) if ($i=="NM:i:0") x++} END {print s, e, c, n, x+0, (n ? sprintf("%.4f", x/n) : "NA")}' >> hifi/check/$SRR.who.tsv
  [ $c = $SM ] || rm -f $W/$c.v11.mmi
done
cat hifi/check/$SRR.who.tsv
BEST=$(awk -F'\t' 'NR>1 {print $6"\t"$3}' hifi/check/$SRR.who.tsv | sort -k1,1gr | head -1 | cut -f2)
[ "$BEST" = "$SM" ] || { log "ERROR sample check: best $BEST, expected $SM; not mapped"; exit 1; }
log "sample check ok: $SM"
$SRA/fastq-dump -Z $F 2> $W/dump.log > $W/all.fq; rm -f $F
log "reads $(awk 'NR%4==2' $W/all.fq | wc -l); dump log: $(tail -1 $W/dump.log)"
for v in $VS; do B=hifi/bam/$SRR.$SM.$v.bam
  idx $SM $v
  $MM -ax map-hifi --eqx -Y -t $((C - 4)) $W/$SM.$v.mmi $W/all.fq 2> $W/map.$v.log | $S sort -@ 4 -m 3G -T $W/sort -o $W/out.bam -
  $S index -@ 8 $W/out.bam && $S quickcheck $W/out.bam && cp $W/out.bam $B && cp $W/out.bam.bai $B.bai
  log "$v: $B records $($S view -@ 8 -c $B), primary $($S view -@ 8 -c -F 0x904 $B), unmapped $($S view -@ 8 -c -f 4 $B), $(stat -c %s $B) bytes"
  rm -f $W/$SM.$v.mmi $W/out.bam*
done
log done
