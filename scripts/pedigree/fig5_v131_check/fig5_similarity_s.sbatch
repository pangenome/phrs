#!/usr/bin/env bash
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=48
#SBATCH --mem=480G
#SBATCH --time=1-00:00:00
# Fig5 steps 5-7, scoring as argument
# Args: comparison_id inputs_dir raw_paf out_dir [scoring, default ani]
set -euo pipefail
CID=$1; IN=$2; RAW=$3; OUT=$4; SC=${5:-ani}
SWEEPGA=/home/erikg/.cargo/bin/sweepga
IMPG=/home/erikg/.cargo/bin/impg
BGZIP=$HOME/.local/bin/bgzip
PIGZ=/usr/bin/pigz
FS=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts
CENT=/moosefs/erikg/phrs/data/chm13-annotations.bed
T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.fig5sim.$SLURM_JOB_ID
mkdir -p "$W" "$OUT"/{filtered_paf,beds,outputs,summaries}; trap 'rm -rf "$W"' EXIT
export LC_ALL=C RAYON_NUM_THREADS=$T OMP_NUM_THREADS=$T
echo "[$(date -Is)] host=$(hostname) cid=$CID raw=$RAW out=$OUT work=$W"
$SWEEPGA --version; $IMPG --version
# stage inputs to scratch
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$IN/$CID.target.fa.fai" \
   "$CENT" "$FS/filter_query_windows_by_paf_depth.py" "$FS/filter_impg_similarity_class_winners.py" "$W/"
$PIGZ -dc -p $T "$RAW" > "$W/raw.paf"
echo "[$(date -Is)] raw PAF lines=$(wc -l < "$W/raw.paf")"
cd "$W"
echo "[$(date -Is)] step5 sweepga 10:10 $SC"
$SWEEPGA --num-mappings 10:10 --scaffold-jump 0 --scoring $SC \
  --temp-dir "$W" --output-file filt.paf raw.paf
rm -f raw.paf
echo "[$(date -Is)] filtered PAF lines=$(wc -l < filt.paf)"
$BGZIP -@ $T -c filt.paf > "$CID.10to10.paf.gz"; rm -f filt.paf
echo "[$(date -Is)] step6 2kb depth windows"
python3 filter_query_windows_by_paf_depth.py --paf "$CID.10to10.paf.gz" \
  --query-fai "$CID.query.fa.fai" --centromere-bed chm13-annotations.bed \
  --out-bed "$CID.query_2000bp.bed" --report "$CID.depth_report.tsv" --summary "$CID.depth_summary.tsv" \
  --window-size 2000 --min-depth 1 --max-depth 100 --interchrom-only --pigz-threads $T
echo "[$(date -Is)] windows=$(wc -l < "$CID.query_2000bp.bed")"
echo "[$(date -Is)] step7 impg similarity class winners"
$IMPG similarity --alignment-files "$CID.10to10.paf.gz" --target-bed "$CID.query_2000bp.bed" \
  --sequence-files "$CID.query.fa" "$CID.target.fa" --gfa-engine poa --no-merge \
  --num-mappings many:many --scaffold-jump 0 --threads $T \
  | python3 filter_impg_similarity_class_winners.py --max-candidates 5000 --skip-report "$CID.class_skip.tsv" \
  | $PIGZ -p $T > "$CID.class_winners.impg_similarity.tsv.gz"
gzip -t "$CID.class_winners.impg_similarity.tsv.gz"
echo "[$(date -Is)] class_winners lines=$($PIGZ -dc "$CID.class_winners.impg_similarity.tsv.gz" | wc -l) skipped=$(($(wc -l < "$CID.class_skip.tsv")-1))"
# copy results back to moosefs
cp "$CID.10to10.paf.gz" "$OUT/filtered_paf/"
cp "$CID.query_2000bp.bed" "$OUT/beds/"
cp "$CID.class_winners.impg_similarity.tsv.gz" "$OUT/outputs/"
cp "$CID.depth_report.tsv" "$CID.depth_summary.tsv" "$CID.class_skip.tsv" "$OUT/summaries/"
echo "[$(date -Is)] DONE $CID"
