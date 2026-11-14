#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Short reads to meryl k=21, one sample
# args: WASHU, SAMPLE, RUNS_TSV (SRR, sample, HG, md5, bytes, S3 URL)
# Per run: S3 download (resume), md5, .sra kept in WASHU/short/sra/, fastq-dump --split-spot | pigz,
# meryl count k=21. Runs merged with meryl union-sum.
# out: WASHU/short/meryl/<SAMPLE>.k21.meryl, WASHU/short/meryl/<SAMPLE>.k21.stats.txt
set -euo pipefail
D=$1; SM=$2; T=$3; cd $D; mkdir -p short/sra short/meryl
SRA=/gnu/store/sjzg68akxlirqmm9vmmhrfvxnvdx7fa7-sra-tools-2.10.6/bin; MER=/gnu/store/8zwv41mb9sgb08652jfn9s4h697550v9-meryl-1.4.1/bin/meryl
W=/scratch/$USER.shortmer.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
log "$SM host $(hostname), scratch free $(df -h --output=avail /scratch | tail -1), meryl $($MER --version 2>&1 | head -1)"
DBS=()
while read -r SRR S HG MD5 BYTES URL; do
  [ "$S" = "$SM" ] || continue
  K=short/sra/$SRR.sra; F=$W/$SRR.sra
  if [ -s $K ] && [ "$(stat -c %s $K)" = "$BYTES" ]; then cp $K $F; log "$SRR reused kept .sra"; else
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
sys.stderr.write(f"downloaded {os.path.getsize(path)} of {size} bytes\n")
PY
  fi
  GOT=$(md5sum $F | cut -d' ' -f1); [ "$GOT" = "$MD5" ] || { log "ERROR $SRR md5 $GOT != $MD5"; exit 1; }
  [ -s $K ] || cp $F $K; log "$SRR md5 ok, kept $K"
  $SRA/fastq-dump --split-spot -Z $F 2> $W/dump.$SRR.log | pigz -p 8 -1 > $W/$SRR.fq.gz
  log "$SRR fastq.gz $(stat -c %s $W/$SRR.fq.gz) bytes, dump log: $(tail -2 $W/dump.$SRR.log | tr '\n' ' ')"
  rm -f $F
  $MER k=21 threads=$SLURM_CPUS_PER_TASK memory=100g count $W/$SRR.fq.gz output $W/$SRR.meryl > $W/count.$SRR.log 2>&1 || { tail -20 $W/count.$SRR.log; exit 1; }
  rm -f $W/$SRR.fq.gz; DBS+=($W/$SRR.meryl)
  log "$SRR meryl: $($MER statistics $W/$SRR.meryl 2>/dev/null | grep -E '^(Number of|  (unique|distinct|present))' | head -4 | tr -s ' ' | tr '\n' ';')"
done < $T
[ ${#DBS[@]} -ge 1 ] || { log "ERROR no runs for $SM in $T"; exit 1; }
O=short/meryl/$SM.k21.meryl
if [ ${#DBS[@]} -eq 1 ]; then cp -r ${DBS[0]} $O; else $MER threads=$SLURM_CPUS_PER_TASK union-sum ${DBS[@]} output $W/u.meryl > $W/union.log 2>&1 && cp -r $W/u.meryl $O; fi
$MER statistics $O > short/meryl/$SM.k21.stats.txt 2>&1
log "done $O ($(du -sh $O | cut -f1)): $(head -8 short/meryl/$SM.k21.stats.txt | tr -s ' ' | tr '\n' ';')"
