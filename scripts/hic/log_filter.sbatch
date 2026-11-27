#!/bin/bash
#SBATCH -c 2
#SBATCH -t 12:00:00
# Keep the useful lines of logs flooded by pore-c-py digest warnings ("Read <id> has no modified bases.", one per
# monomer): one pass over each LOG writes LOG.filtered (every other line, in order) on node /scratch, then copies it next
# to LOG, and logs the lines kept and dropped. Checks that LOG.filtered ends with the job's last line ("exit status ...")
# and keeps every timestamped step line ("[20..-..-..T"). The original LOG is left for the user to delete.
# usage: sbatch log_filter.sbatch LOG [LOG ...]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/logf.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
for L in "$@"; do
  echo "[$(date -Is)] $L: $(stat -c %s $L) bytes"
  LC_ALL=C awk -v s=$W/counts '/ has no modified bases\.$/ {d++; next} {k++; if (/^\[20[0-9][0-9]-[0-9][0-9]-[0-9][0-9]T/) t++; print} END {print k + 0, d + 0, t + 0 > s}' $L > $W/f.log
  read k d t < $W/counts
  tail -1 $W/f.log | grep -q '^exit status ' || { echo "ERROR $L: the filtered log does not end with the exit status line"; exit 1; }
  cp $W/f.log $L.filtered.part && mv $L.filtered.part $L.filtered
  echo "[$(date -Is)] $L.filtered: $(stat -c %s $L.filtered) bytes; lines kept $k (timestamped $t), dropped $d digest warnings; last line: $(tail -1 $L.filtered)"
done
