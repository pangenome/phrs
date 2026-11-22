#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Largest query piece at which FastGA still aligns the transmitted parent copy to the chromosome end (whole chr8 and chr19
# queries stop 1.9 and 8.5 Mb before the end even vs the transmitted copy alone; 2 Mb end pieces reach the end).
# usage: sbatch fastga_piece_scan.sbatch --washu D --tasks TSV --sizes 5000000,10000000,... --out TSV
# TSV lines: label, query FASTA, target FASTA, query contig, end (p|q), transmitted target contig.
# Per line and size S: the end piece (q: last S bp, p: first S bp) vs the whole transmitted target contig, FastGA -f32 (Fig5
# binary); coverage of the last (q) or first (p) 1 Mb of the piece by all records, and bp uncovered at the end.
# pass = >= 950,000 bp covered and < 10,000 bp uncovered at the end (same as fastga_fscan.sbatch).
# out: TSV label, size, seconds, records, covered_bp_1mb, uncovered_at_end, result
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --tasks) TS=$2;; --sizes) SZ=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D TS SZ O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O ] && { log "exists $O, skip"; exit 0; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
B=/home/erikg/.cache/sweepga/0.1.1_fe5037d5e_wd47b7e3e; export PATH=$B:$PATH; T=${SLURM_CPUS_PER_TASK:-4}
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/pscan.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log "tasks $TS ($(wc -l < $TS) lines), sizes $SZ, host $(hostname)"
printf "label\tsize\tseconds\trecords\tcovered_bp_1mb\tuncovered_at_end\tresult\n" > $W/out.tsv
while IFS=$'\t' read -r LAB QF TF C E TC; do
  [ -n "${TC:-}" ] || { echo "ERROR incomplete task line: $LAB"; exit 1; }
  L=$(awk -F'\t' -v c="$C" '$1 == c {print $2}' $QF.fai); [ -n "$L" ] || { echo "ERROR $C not in $QF.fai"; exit 1; }
  $S faidx $TF "$TC" > $W/t.fa
  for s in $(echo $SZ | tr ',' ' '); do
    [ $s -le $L ] || { log "$LAB size $s > contig length $L, skip"; continue; }
    if [ $E = q ]; then a=$((L - s + 1)); else a=1; fi; b=$((a + s - 1))
    $S faidx $QF "$C:$a-$b" > $W/q.fa
    t0=$(date +%s)
    FastGA -v -T$T -P$W -f32 -paf $W/q.fa $W/t.fa > $W/a.paf 2> $W/fg.log || { tail -5 $W/fg.log; exit 1; }
    sec=$(( $(date +%s) - t0 ))
    # piece coordinates: q end = last 1 Mb [s - 1e6, s), p end = first 1 Mb [0, 1e6)
    read -r cov un < <(awk -F'\t' -v s=$s -v e=$E 'BEGIN {lo = (e == "q") ? s - 1000000 : 0; hi = lo + 1000000} {a = ($3 > lo) ? $3 : lo; b = ($4 < hi) ? $4 : hi; if (b > a) print a "\t" b}' $W/a.paf | sort -k1,1n \
      | awk -F'\t' -v s=$s -v e=$E 'BEGIN {lo = (e == "q") ? s - 1000000 : 0; hi = lo + 1000000; m = lo; x = hi} {if ($2 > m) {c += $2 - (($1 > m) ? $1 : m); m = $2} if ($1 < x) x = $1} END {print c + 0, (e == "q") ? hi - m : x - lo}')
    r=fail; [ $cov -ge 950000 ] && [ $un -lt 10000 ] && r=pass
    printf "$LAB\t$s\t$sec\t$(wc -l < $W/a.paf)\t$cov\t$un\t$r\n" >> $W/out.tsv
    log "$LAB piece $C:$a-$b ($s bp) vs $TC: $sec s, covered $cov bp of the end 1 Mb, uncovered at end $un -> $r"
  done
done < $TS
cp $W/out.tsv $O; cat $O
log done
