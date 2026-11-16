#!/usr/bin/env bash
#SBATCH --partition=workers
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
# Fig5 step 4: sweepGA/FastGA raw PAF in /scratch
# Args: comparison_id inputs_dir out_dir
set -euo pipefail
CID=$1; IN=$2; OUT=$3
SWEEPGA=/home/erikg/.cargo/bin/sweepga
PIGZ=/usr/bin/pigz
W=/scratch/$USER.fig5align.$SLURM_JOB_ID
mkdir -p "$W" "$OUT/raw_paf"; trap 'rm -rf "$W"' EXIT
echo "[$(date -Is)] host=$(hostname) cid=$CID in=$IN out=$OUT work=$W"
$SWEEPGA --version
# stage inputs to scratch
cp "$IN/$CID.query.fa" "$IN/$CID.target.fa" "$W/"
ls -la "$W"
cd "$W"
echo "[$(date -Is)] sweepga fastga f32 many:many"
$SWEEPGA --fastga --fastga-frequency 32 --num-mappings many:many --scaffold-jump 0 \
  --temp-dir "$W" --output-file "$W/$CID.raw.paf" "$CID.query.fa" "$CID.target.fa"
echo "[$(date -Is)] raw PAF lines=$(wc -l < "$W/$CID.raw.paf")"
$PIGZ -p 32 -c "$W/$CID.raw.paf" > "$W/$CID.frequency32.paf.gz"
cp "$W/$CID.frequency32.paf.gz" "$OUT/raw_paf/"
echo "[$(date -Is)] DONE $OUT/raw_paf/$CID.frequency32.paf.gz"
