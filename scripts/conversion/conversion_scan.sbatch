#!/bin/bash
#SBATCH -J convscan
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -p workers
# Gene-conversion scan (R1-M1 test 3b), one array task per chunk of target-centric PAFs.
# Usage: sbatch --array=1-K -o logs/convscan.%A_%a.log conversion_scan.sbatch LIST K OUTDIR
#   LIST   file with one PAF path per line; task N processes lines N, N+K, N+2K, ...
#   K      number of chunks (= array size)
#   OUTDIR where <OUTDIR>/chunk_N.{bins,tracts}.tsv are written
set -euxo pipefail
LIST=$1; K=$2; OUTDIR=$3; N=${SLURM_ARRAY_TASK_ID}
SCRIPT_DIR=/moosefs/guarracino/HPRCv2/PHR_III/conversion/scripts
hostname; date '+%F %T %Z'
TMP=/scratch/$USER/tmp/phrs/$SLURM_JOB_ID.$N
mkdir -p "$TMP" "$OUTDIR"
trap 'rm -rf "$TMP"' EXIT
# stage this chunk's PAFs on node-local disk (moosefs is slow)
awk -v n="$N" -v k="$K" '(NR - n) % k == 0' "$LIST" > "$TMP/paths.txt"
wc -l < "$TMP/paths.txt"
mkdir -p "$TMP/paf"
xargs -a "$TMP/paths.txt" cp -t "$TMP/paf/"
ls "$TMP"/paf/*.paf.gz > "$TMP/local.txt"
python3 "$SCRIPT_DIR/conversion_scan.py" --list "$TMP/local.txt" --out "$TMP/chunk_$N"
mv "$TMP/chunk_$N.bins.tsv" "$TMP/chunk_$N.tracts.tsv" "$OUTDIR/"
date '+%F %T %Z'
