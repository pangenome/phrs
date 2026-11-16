#!/bin/bash
# Submit whole-genome ONT mapping in lanes
# args: WASHU dir, OUT_DIR
# env LANES="sbatch opts;sbatch opts;.." one job at a time per lane; DEPS="jobid;jobid;.." first dependency per lane (empty = none)
# env SKIP="SRR.w1.v11 .." jobs not to submit (already running)
# Per version (v11 first, then v131): every run on its own diploid, then on its cross references (ont_runs_cross.tsv).
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p logs $OUT
IFS=';' read -ra LANE <<< "${LANES:--p tux -w tux05;-p tux -w tux05;-p workers -w octopus10;-p workers -w octopus11}"
IFS=';' read -ra last <<< "${DEPS:-}"; i=0; NL=${#LANE[@]}
sub () {  # srr sample hg md5 bytes url refs wave version
  case " ${SKIP:-} " in *" $1.w$8.$9 "*) echo "skip $1.w$8.$9"; return;; esac
  k=$((i++ % NL)); dep=${last[$k]:+--dependency=afterany:${last[$k]}}
  id=$(sbatch --parsable ${LANE[$k]} $dep -J of$8$9-$1 -o logs/ontfull.$1.w$8.$9.%j.log scripts_tmp/ont_map_full.sbatch $D $1 $2 $4 $5 $6 $D/ont_sra/$3 $OUT "$7" $9)
  last[$k]=$id; echo "$id $1 $2 $9 wave=$8 refs=[$7] lane=$k ${dep:-}"
}
RUNS=$(cat scripts_tmp/ont_runs.tsv scripts_tmp/ont_runs_grandparents.tsv | sort -u)
for v in v11 v131; do
for s in PAN011 PAN027 PAN028 PAN010; do
  while IFS=$'\t' read -r srr sm hg md5 bytes url; do sub $srr $sm $hg $md5 $bytes $url $sm 1 $v; done < <(echo "$RUNS" | awk -F'\t' -v s=$s '$2 == s')
done
while IFS=$'\t' read -r srr sm hg md5 bytes url; do
  refs=$(awk -F'\t' -v r=$srr '$1 == r {printf "%s ", $7}' scripts_tmp/ont_runs_cross.tsv | sed 's/ $//')
  [ -n "$refs" ] && sub $srr $sm $hg $md5 $bytes $url "$refs" 2 $v
done < <(echo "$RUNS")
done
