#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --exclude=octopus05
# Does each piece restore the original?
# usage: sbatch g1_sim_check.sbatch --washu D --sim SIM_DIR --sample SM --version V --hap H --out TSV
# Per event and mode with a piece: patched = modified[start-50k, start) + piece + modified[end, end+50k), original =
# norm_V on the matching interval (shift from truth.tsv); minimap2 2.31 asm5 --eqx: edit distance (NM) and length difference.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sim) SD=$2;; --sample) SM=$2;; --version) V=$2;; --hap) H=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SD SM V H OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simchk.XXXX)
printf "event\tmode\tcontig\tpiece_start\tpiece_end\tpatched_bp\toriginal_bp\taligned_bp\tnm\tlength_diff\n" > $OUT
while IFS=$'\t' read -r i n e r; do
  for p in $SD/locus/$i.$n/*.piece.fa; do
    [ -e "$p" ] || continue
    read -r m c a z <<< "$(head -1 $p | tr -d '>' | cut -d' ' -f1-4)"
    read -r o0 o1 m0 m1 <<< "$(awk -F'\t' -v c="$c" -v i=$i 'NR == i + 1 && $3 == c {print $4, $5, $6, $7}' $SD/truth.tsv)"
    [ -n "${m1:-}" ] || { echo "ERROR event $i not found in truth.tsv for $c"; exit 1; }
    { echo ">patched"; { $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((a - 49999))-$a" | tail -n +2; tail -n +2 $p
      $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((z + 1))-$((z + 50000))" | tail -n +2; } | tr -d '\n'; echo; } > $W/p.fa
    $S faidx norm_$V/$SM#$H.fa.gz "$c:$((a - 49999 - (m0 - o0)))-$((z - (m1 - o1) + 50000))" | sed '1s/.*/>original/' > $W/o.fa
    $MM -cx asm5 --eqx $W/o.fa $W/p.fa 2> /dev/null | sort -k10,10nr | head -1 > $W/a.paf
    PL=$(tail -n +2 $W/p.fa | tr -d '\n' | wc -c); OL=$(tail -n +2 $W/o.fa | tr -d '\n' | wc -c)
    AL=$(awk '{print $10 + 0}' $W/a.paf); NM=$(awk '{for (i=13; i<=NF; i++) if ($i ~ /^NM:i:/) print substr($i, 6)}' $W/a.paf)
    printf "$i.$n\t$m\t$c\t$a\t$z\t$PL\t$OL\t${AL:-0}\t${NM:-NA}\t$((PL - OL))\n" >> $OUT
  done
done < $SD/locus_tasks.tsv
cat $OUT
rm -rf "$W"
