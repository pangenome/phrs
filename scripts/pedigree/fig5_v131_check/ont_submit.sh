#!/bin/bash
# Submit one ONT job per run
# args: WASHU dir, runs tsv (default ont_runs.tsv)
# optional 7th column: REF assembly
# env LANES=n lanes, one run per lane
# DEPS="j1 j2 .." first job per lane
# SB_EXTRA: extra sbatch options
set -euo pipefail
D=$1; cd $D; mkdir -p logs ont; i=0; declare -A last lane
while IFS=$'\t' read -r srr s hg md5 bytes url ref; do
  D0=(${DEPS:-}); dep=; [ -n "${lane[$srr]:-}" ] || { lane[$srr]=$((i++ % ${LANES:-1000000})); }; k=${lane[$srr]}; dep=${last[$k]:-${D0[$k]:-}}; dep=${dep:+afterany:$dep}
  id=$(sbatch --parsable ${SB_EXTRA:-} ${dep:+--dependency=$dep} -J ont-$srr${ref:+-$ref} -o logs/ont.$srr.%j.log ont_map.sbatch $D $srr $s $md5 $bytes $url $D/ont_sra/$hg $D/ont $ref)
  last[$k]=$id; echo "$id $srr $ref lane=$k dep=${dep:-none}"
done < ${2:-ont_runs.tsv}
