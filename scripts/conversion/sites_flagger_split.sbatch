#!/usr/bin/env bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05,octopus06
# Flagger overlap of site-level events, significant and non-significant runs apart (sites_cluster_null.py split), same call
# as the earlier sites.flagger.tsv (code.md: conversion_flagger_overlap.py on sites.events.as_tracts.tsv).
# usage (from /moosefs/guarracino/HPRCv2/PHR_III/conversion): sbatch sites_flagger_split.sbatch K15_TSV
set -euxo pipefail
for k in sig nonsig; do
  echo "[$(date '+%F %T')] $k: $(($(wc -l < sites/sites.$k.as_tracts.tsv) - 1)) events"
  python3 scripts/conversion_flagger_overlap.py --flagger-dir flagger/bed --chunks all/out --tracts sites/sites.$k.as_tracts.tsv --communities $1 --out sites/sites.$k.flagger.tsv
  grep -E "^ALL" sites/sites.$k.flagger.tsv
done
echo "[$(date '+%F %T')] DONE"
