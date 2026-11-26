#!/bin/bash
# Manual fixes of g1_manual_fixes.tsv for one version and sample, in file order, on top of IN_DIR (g1_rebuild_person.sbatch step 8)
# usage: bash g1_manual_fixes.sh D V SM IN_DIR OUT_PREFIX   (inside a Slurm job: hifi_polish.sbatch uses SLURM_CPUS_PER_TASK)
# end_replace: g1_end_replace.py; then, when g1_end_polish.tsv lists the contig, hifi_polish.sbatch on the inserted piece
#   (junction + 1 to the contig end, 1-based) with that row's HiFi read regions. Without this polish the DAD chr9 h2 tail kept
#   1 indel (3 bp, 2,382 bp from the tip) that the 2026-11-15 manual build (g1/corrected_v11b_pol) had polished.
# tip_trim: tip_trim.py.
# out: OUT_PREFIX<N> per fix (OUT_PREFIX<N>.pol/norm after a polish); the last stdout line is the final dir; log on stderr
set -euxo pipefail
D=$1; V=$2; SM=$3; F=$4; P=$5; T=$D/scripts_tmp; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
log () { echo "[$(date -Is)] $*" >&2; }
i=0
while IFS=$'\t' read -r ver sm step c cf pos src why; do
  [ "$ver" = "$V" ] && [ "$sm" = "$SM" ] || continue
  i=$((i + 1)); log "fix $i: $step $c pos $pos ($why)"
  case $step in
    end_replace)
      IFS=: read -r pf pn pfrom <<< "$src"; python3 $T/g1_end_replace.py $S $BGZ $F $P$i $c $cf $pos $pf $pn $pfrom; F=$P$i
      PK=$(awk -F'\t' -v v=$V -v s=$SM -v c="$c" '$1 == v && $2 == s && $3 == c {print $4}' $T/g1_end_polish.tsv)
      if [ -n "$PK" ]; then
        read -r jc j L < $F/junction.tsv; [ "$jc" = "$c" ] || { echo "ERROR $F/junction.tsv is for $jc, not $c"; exit 1; }
        bash $T/hifi_polish.sbatch $D $SM $V $F "$c:$((j + 1))-$L" "$PK" $P$i.pol >&2; F=$P$i.pol/norm
        log "fix $i polished $c:$((j + 1))-$L: $(tail -n +2 $P$i.pol/polish.tsv | cut -f4-7 | tr '\t\n' ' ;')"
      else log "fix $i: no polish row in g1_end_polish.tsv"; fi;;
    tip_trim) python3 $T/tip_trim.py $S $BGZ $cf $c $pos $(ls $src | paste -sd,) $F $P$i; F=$P$i;;
    *) echo "ERROR unknown step $step"; exit 1;;
  esac
done < <(tail -n +2 $T/g1_manual_fixes.tsv)
log "manual fixes: $i; final $F"; echo $F
