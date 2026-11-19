#!/bin/bash
# List my stale scratch dirs per node
# usage: bash scratch_stale.sh --dry-run | --yes NODE:DIR ...   (run on octopus01)
# --dry-run: on each node where I have /scratch/$USER/tmp/phrs, lists locus.<jobid> and mmi.<jobid> dirs whose job is
# not running, and mmi.XXXX dirs (old random names) with size and time; nothing is removed.
# --yes NODE:DIR ...: removes exactly the listed dirs (the user copies them from the dry run).
set -euo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
[ $# -ge 1 ] || { echo "ERROR use --dry-run or --yes NODE:DIR ..."; exit 2; }
if [ "$1" = --dry-run ]; then
  R=" $(squeue -u $USER -h -o '%A' | tr '\n' ' ') "
  for n in $(sinfo -h -N -o '%N' | sort -u); do
    [ "$n" = octopus05 ] && continue
    rc=0; L=$(ssh -o ConnectTimeout=10 -o BatchMode=yes $n "if [ -d /scratch/$USER/tmp/phrs ]; then find /scratch/$USER/tmp/phrs -maxdepth 1 -name 'locus.*' -o -maxdepth 1 -name 'mmi.*'; fi") || rc=$?
    [ $rc -eq 0 ] || { echo "# $n: ssh exit $rc, not checked"; continue; }
    [ -n "$L" ] || continue
    echo "$L" | while read -r d; do
      j=${d##*.}; case "$R" in *" $j "*) continue;; esac; t=""; [[ $j =~ ^[0-9]+$ ]] || t="OWNER-UNKNOWN(random name, may belong to a running job) "
      echo "$t$n:$d $(ssh -n -o ConnectTimeout=10 $n "du -sb $d | cut -f1; stat -c %y $d | cut -c1-16" | tr '\n' ' ')"
    done
  done
  echo "dry-run only, nothing removed"
elif [ "$1" = --yes ]; then
  shift; for x in "$@"; do n=${x%%:*}; d=${x#*:}; case "$d" in /scratch/$USER/tmp/phrs/locus.*|/scratch/$USER/tmp/phrs/mmi.*) ;; *) echo "ERROR refuse $d"; exit 1;; esac
    ssh -o ConnectTimeout=10 $n "rm -rf $d" && echo "removed $n:$d"; done
else echo "ERROR unknown $1"; exit 2; fi
