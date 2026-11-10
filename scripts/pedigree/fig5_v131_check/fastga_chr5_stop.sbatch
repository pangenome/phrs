#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=300G
#SBATCH -t 1-00:00:00
# Where does FastGA stop on SELF.p chr5 end, by seed frequency
# args: WASHU, OUT_DIR
# Same binary as the Fig5 runs (sweepga 0.1.1 cache). Query = PAN027#2#chr5 only, target = the joint PAN011 diploid of
# each version (inputs_v11, inputs_v131). Runs -f32 (Fig5) and -f255 (maximum), no other option changed (-c 85 -s 1000 -l 100 -i .7).
# out: chr5.<v>.f<F>.paf, summary.tsv (alignments of the last 1.5 Mb of the query: span, target, identity)
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p $OUT
B=/home/erikg/.cache/sweepga/0.1.1_fe5037d5e_wd47b7e3e; export PATH=$B:$PATH
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.fgstop.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host $(hostname), FastGA $($B/FastGA 2>&1 | head -1), scratch $(df -h --output=avail $W | tail -1)"
printf "version\tfreq\tquery_start\tquery_end\tstrand\ttarget\ttarget_start\ttarget_end\tblock\tidentity\n" > $OUT/summary.tsv
for v in v11 v131; do
  $S faidx norm_$v/PAN027#2.fa.gz PAN027#2#chr5 > $W/q.$v.fa
  cp inputs_$v/PAN027pat_vs_PAN011_joint.target.fa $W/t.$v.fa
  for F in 32 255; do
    log "$v -f$F: start"
    FastGA -v -k -T$SLURM_CPUS_PER_TASK -P$W -f$F -paf $W/q.$v.fa $W/t.$v.fa > $W/chr5.$v.f$F.paf 2> $W/fg.$v.f$F.log || { log "ERROR FastGA $v f$F"; tail -5 $W/fg.$v.f$F.log; exit 1; }
    cp $W/chr5.$v.f$F.paf $OUT/; cp $W/fg.$v.f$F.log $OUT/fastga.$v.f$F.log
    L=$(cut -f2 $W/q.$v.fa.fai 2>/dev/null || awk '$1 !~ />/ {n += length($0)} END {print n}' $W/q.$v.fa)
    awk -F'\t' -v v=$v -v f=$F -v L=$L '$4 > L - 1500000 && $11 >= 5000 {t=$6; sub(/PAN011#joint#/, "", t); printf "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%.4f\n", v, f, $3, $4, $5, t, $8, $9, $11, $10/$11}' $W/chr5.$v.f$F.paf | sort -k3,3n >> $OUT/summary.tsv
    log "$v -f$F: $(wc -l < $W/chr5.$v.f$F.paf) PAF lines; query chr5 length $L; longest alignment to h2_chr5: $(awk -F'\t' '$6 ~ /h2_chr5$/ {print $3"-"$4" ("$11" bp)"}' $W/chr5.$v.f$F.paf | sort -t'(' -k2,2nr | head -1)"
  done
done
log done
