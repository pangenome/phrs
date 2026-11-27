#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# PHR FASTA from a calls TSV, with the paper's recipe (code_PHR-and-3D.md, "Find inter-chromosomal regions"): for each
# flank with a region, samtools faidx of flank:(start+1)-end, renamed contig:start-end_chrN_arm (p arm: keep the flank
# start, end = start + region_end; q arm: start = flank start + region_start, keep the flank end); bgzip -l 9, faidx.
# usage: sbatch extract_phr_fasta.sbatch --calls TSV --flanks FA_GZ --out FA (writes FA.gz, FA.gz.fai, FA.gz.gzi)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --calls) REG=$2;; --flanks) FA=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in REG FA OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
[ -e "$OUT" ] || [ -e "$OUT.gz" ] && { echo "ERROR $OUT or $OUT.gz exists"; exit 1; }
set +x
n=0
while read -r seq_name int_start int_end; do
    prefix=$(echo "$seq_name" | cut -d ':' -f 1)
    coords_suffix=$(echo "$seq_name" | cut -d ':' -f 2)
    orig_start=$(echo "$coords_suffix" | cut -d '-' -f 1)
    rest=$(echo "$coords_suffix" | cut -d '-' -f 2)
    orig_end=$(echo "$rest" | cut -d '_' -f 1)
    suffix=$(echo "$rest" | sed 's/^[0-9]*//')
    if [[ "$seq_name" == *"_parm"* ]]; then new_start=$orig_start; new_end=$((orig_start + int_end))
    else new_start=$((orig_start + int_start)); new_end=$orig_end; fi
    new_name="${prefix}:${new_start}-${new_end}${suffix}"
    region="${seq_name}:$((int_start + 1))-${int_end}"
    $S faidx "$FA" "$region" | sed "1s/.*/>$new_name/" >> "$OUT"
    n=$((n + 1))
done < <(sed 1d "$REG" | cut -f 1,4,5 | awk '$2 != "."')
set -x
echo "[$(date -Is)] extracted $n sequences, $(grep -c '^>' $OUT) headers"
$BGZ -@ 4 -l 9 "$OUT"
$S faidx "$OUT.gz"
echo "[$(date -Is)] $OUT.gz: $(wc -l < $OUT.gz.fai) sequences, $(awk '{s+=$2} END {print s}' $OUT.gz.fai) bp"
