#!/usr/bin/env bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
# Fig5 step 4 at one -f, with coverage check
# args: CID, IN_DIR, OUT_DIR, FREQ
# Same as fig5_stageA/fig5_align.sbatch, FREQ is an argument and every query contig >= 10 Mb must be aligned over >= 50% (union).
# Fails (exit 1) when a contig is missing: Run C regeneration 2871643 lost chr1 without error.
# out: OUT_DIR/raw_paf/CID.frequency<FREQ>.paf.gz
set -euo pipefail
CID=$1; IN=$2; OUT=$3; F=$4
SWEEPGA=/home/erikg/.cargo/bin/sweepga
PIGZ=/usr/bin/pigz
W=/scratch/$USER.fig5alnf.$SLURM_JOB_ID
mkdir -p "$W" "$OUT/raw_paf"; trap 'rm -rf "$W"' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host=$(hostname) cid=$CID f=$F in=$IN out=$OUT work=$W"
$SWEEPGA --version
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$W/"
ls -la "$W"; cd "$W"
log "sweepga fastga f$F many:many"
$SWEEPGA --fastga --fastga-frequency $F --num-mappings many:many --scaffold-jump 0 \
  --temp-dir "$W" --output-file "$W/$CID.raw.paf" "$CID.query.fa" "$CID.target.fa"
N=$(wc -l < "$CID.raw.paf")
log "raw PAF lines=$N"
log "coverage check, contig length aligned_union fraction"
cut -f1,3,4 "$CID.raw.paf" | sort -S 8G -k1,1 -k2,2n \
  | awk -F'\t' 'NR==FNR {len[$1]=$2; next} {if ($1!=c) {if (c!="") u[c]=t; c=$1; t=0; e=0} if ($3>e) {s=($2>e)?$2:e; t+=$3-s; e=$3}} END {u[c]=t; bad=0; for (k in len) {f=(k in u)?u[k]/len[k]:0; printf "%s\t%d\t%d\t%.3f\n", k, len[k], u[k]+0, f; if (len[k]>=10000000 && f<0.5) bad=1}; exit bad}' \
    "$CID.query.fa.fai" - | sort -V | tee "$W/cov.tsv" \
  || { log "ERROR contig coverage below 0.5 (see table above)"; exit 1; }
log "coverage ok"
$PIGZ -p 32 -c "$CID.raw.paf" > "$CID.frequency$F.paf.gz"
gzip -t "$CID.frequency$F.paf.gz"
cp "$CID.frequency$F.paf.gz" "$OUT/raw_paf/"; cp cov.tsv "$OUT/raw_paf/$CID.frequency$F.coverage.tsv"
log "DONE $OUT/raw_paf/$CID.frequency$F.paf.gz lines=$N"
