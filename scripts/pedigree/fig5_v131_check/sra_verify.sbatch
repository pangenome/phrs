#!/bin/bash
#SBATCH -p workers
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 2-00:00:00
# md5 check of kept .sra, remove bad ones
# args: WASHU, then run tsv files (srr sample hg md5 bytes ...)
set -uo pipefail
D=$1; shift
log () { echo "[$(date -Is)] $*"; }
cat "$@" | awk -F'\t' '!s[$1]++ {print $1, $3, $4, $5}' > /scratch/$USER.sraverify.$SLURM_JOB_ID.lst
chk () { local srr=$1 hg=$2 md5=$3 bytes=$4 f=$5/ont_sra/$2/$1.sra
  [ -e $f ] || { echo "$srr MISSING"; return; }
  sz=$(stat -c %s $f); [ "$sz" = "$bytes" ] || { echo "$srr BADSIZE $sz"; rm -f $f; return; }
  got=$(md5sum $f 2>/dev/null | cut -d' ' -f1)
  if [ "$got" = "$md5" ]; then echo "$srr OK"; else echo "$srr BAD got=${got:-IOERR}"; rm -f $f; fi; }
export -f chk
log "files $(wc -l < /scratch/$USER.sraverify.$SLURM_JOB_ID.lst), 3 in parallel"
xargs -P3 -L1 bash -c 'chk $0 $1 $2 $3 '"$D" < /scratch/$USER.sraverify.$SLURM_JOB_ID.lst | tee $D/logs/sra_verify.$SLURM_JOB_ID.txt
rm -f /scratch/$USER.sraverify.$SLURM_JOB_ID.lst; log done
