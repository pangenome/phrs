#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Smallest FastGA -f at which the transmitted parent copy reaches a chromosome end (Fig5 1p/5q-type stop)
# usage: sbatch --array=1-N%K fastga_fscan.sbatch --washu D --tasks TSV --f-from 32 --f-to 48 --out DIR
# TSV lines: label, query FASTA (Fig5 inputs), target FASTA (joint parent), query contig, end (p|q), transmitted copy (h1|h2)
# Same binary and options as the Fig5 sweepga 0.1.1 runs (fastga_chr5_stop.sbatch): FastGA -k -f F -paf, query = the one
# contig. For F = from, from+1, ...: coverage of the last (q) or first (p) 1 Mb by alignments to the transmitted copy of
# the same chromosome (sweepga_end_reach.py); stops at the first F with >= 950,000 bp covered and < 10,000 bp uncovered at
# the end. The first F (32) shows whether the one-contig run reproduces the whole-genome stop.
# out: DIR/<label>.f<F>.paf.gz, DIR/<label>.tsv (F, seconds, h1/h2 uncovered at end and covered bp, pass)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --tasks) TS=$2;; --f-from) F0=$2;; --f-to) F1=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D TS F0 F1 O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; mkdir -p $O
read -r LAB QF TF C E CP < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $TS)
[ -n "${CP:-}" ] || { echo "ERROR task line $SLURM_ARRAY_TASK_ID incomplete in $TS"; exit 1; }
[ -s $O/$LAB.tsv ] && grep -q -P "\tpass$" $O/$LAB.tsv && { log "exists $O/$LAB.tsv with a pass, skip"; exit 0; }
B=/home/erikg/.cache/sweepga/0.1.1_fe5037d5e_wd47b7e3e; export PATH=$B:$PATH
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/fscan.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
$S faidx $QF "$C" > $W/q.fa; cp $TF $W/t.fa
log "$LAB: query $C ($E end) from $QF, target $TF, transmitted $CP, F $F0..$F1, host $(hostname), FastGA $($B/FastGA 2>&1 | head -1)"
printf "f\tseconds\th1_uncovered_at_end\th1_covered_bp_1mb\th2_uncovered_at_end\th2_covered_bp_1mb\tresult\n" > $O/$LAB.tsv
for F in $(seq $F0 $F1); do
  t0=$(date +%s)
  FastGA -v -k -T$SLURM_CPUS_PER_TASK -P$W -f$F -paf $W/q.fa $W/t.fa > $W/a.paf 2> $W/fg.log || { log "ERROR FastGA f$F"; tail -5 $W/fg.log; exit 1; }
  sec=$(( $(date +%s) - t0 )); gzip -c $W/a.paf > $O/$LAB.f$F.paf.gz
  read -r u1 c1 u2 c2 < <(python3 scripts_tmp/sweepga_end_reach.py $LAB < $W/a.paf 2> /dev/null | awk -F'\t' -v e=$E 'NR>1 && $4==e {print $5, $6, $7, $8}')
  [ -n "${c2:-}" ] || { echo "ERROR no end-reach row for $C $E at f$F"; exit 1; }
  [ $CP = h1 ] && { u=$u1; c=$c1; } || { u=$u2; c=$c2; }
  r=fail; [ $c -ge 950000 ] && [ $u -lt 10000 ] && r=pass
  printf "$F\t$sec\t$u1\t$c1\t$u2\t$c2\t$r\n" >> $O/$LAB.tsv
  log "$LAB f$F ($sec s): transmitted $CP covered $c bp of the end 1 Mb, uncovered at end $u bp -> $r"
  [ $r = pass ] && break
done
cat $O/$LAB.tsv
log done
