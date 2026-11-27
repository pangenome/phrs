#!/bin/bash
#SBATCH -J ks_flanks
#SBATCH -p compute
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o %x.%j.log
# KaryoScope annotation (all six feature sets of HKS_human_CHM13_v2, smoothed BEDs) of telomere-anchored flanks, on gemini.
# Usage: sbatch karyoscope_flanks.sbatch IN.fa OUTDIR
set -euxo pipefail
IN=$1; OUT=$2
E=/tgen_labs/barthel/software/miniforge3_new/envs/karyoscope/bin
DB=/tgen_labs/barthel/software/github/barthel/KaryoScope/resources/registry
export PATH=$E:$PATH
date; hostname; karyoscope --version; md5sum "$IN"; grep -c '>' "$IN"
mkdir -p "$OUT"
karyoscope -v annotate -i "$IN" -o "$OUT" --db-root "$DB" --db HKS_human_CHM13_v2 -t 16 --no-bgzip --no-keep-presmoothed
ls -la "$OUT"; wc -l "$OUT"/*.bed; date
