#!/usr/bin/env bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11
# Why sweepga's ALNtoPAF step dies (code None = killed by a signal) on one query set: KID chr18 alone (41 pieces of 2 Mb,
# as in fig5_align_chunk.sbatch) died in 2881146; chr6 and chr11 alone passed; RSS 1.2-2.4 GB, so not memory.
# usage: sbatch fig5_alnpaf_diag.sbatch --cid CID --in IN_DIR --piece BP --freq F --contig NAME --out OUT_DIR
# 2881152: sweepga --output-file .1aln exits 1 at once; sweepga PAF with ulimit -s unlimited still dies in ALNtoPAF.
# 1. FastGA -1:x.1aln -T8 -f F, the call sweepga makes (its log)
# 2. ALNtoPAF on x.1aln: -x -T8, -x -T1, no CIGAR -T8, -x -T8 with ulimit -s unlimited; exit (128 + signal), lines, last record
# 3. FastGA -pafx (FastGA writes PAF itself); x.1aln kept in OUT_DIR when 2a fails
# out: OUT_DIR/diag.tsv (test, exit, signal, lines) and OUT_DIR/<test>.paf.gz of the tests that passed; log on stdout
set -euxo pipefail
while [ $# -gt 0 ]; do case $1 in --cid) CID=$2;; --in) IN=$2;; --piece) P=$2;; --freq) F=$2;; --contig) CT=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in CID IN P F CT OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
OUT=$(realpath -m "$OUT"); mkdir -p "$OUT"
SWEEPGA=/home/erikg/.cargo/bin/sweepga; A2P=/home/erikg/.cargo/bin/ALNtoPAF; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/alnpaf.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log "host=$(hostname) cid=$CID contig=$CT piece=$P f=$F out=$OUT work=$W stack limit $(ulimit -s) KB"
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$W/"; cd "$W"
awk -F'\t' -v p=$P -v c="$CT" '$1 == c {n = int(($2 + p - 1) / p); s = int(($2 + n - 1) / n); for (a = 1; a <= $2; a += s) {b = a + s - 1; if (b > $2) b = $2; print $1 ":" a "-" b}}' $CID.query.fa.fai > p.txt
[ -s p.txt ] || { echo "ERROR contig $CT not in $CID.query.fa.fai"; exit 1; }
$S faidx -r p.txt $CID.query.fa > q.fa; log "pieces $(wc -l < p.txt)"
printf "test\texit\tsignal\tpaf_lines\n" > $OUT/diag.tsv
rec () { local s=0; [ $2 -gt 128 ] && s=$(kill -l $(($2 - 128))); printf "%s\t%s\t%s\t%s\n" $1 $2 $s $(wc -l < $1.paf) >> $OUT/diag.tsv
  log "$1: exit $2 signal $s, PAF lines $(wc -l < $1.paf), last record: $(tail -1 $1.paf | cut -f1-12)"
  [ $2 -eq 0 ] && gzip -c $1.paf > $OUT/$1.paf.gz; rm -f $1.paf; }
FG=/home/erikg/.cargo/bin/FastGA; export PATH=/home/erikg/.cargo/bin:$PATH  # FastGA calls FAtoGDB and GIXmake by name (2881165: not found)
log "1 FastGA -1:x.1aln -T8 -f$F (the call sweepga makes, from its log)"; rc=0; $FG -1:x.1aln -T8 -P$W -f$F q.fa $CID.target.fa || rc=$?
log "1 FastGA exit $rc"; printf "fastga_1aln\t%s\tNA\tNA\n" $rc >> $OUT/diag.tsv; [ $rc -eq 0 ] || exit 1; ls -la x.1aln
log "2a ALNtoPAF -x -T8"; rc=0; $A2P -x -T8 x.1aln > a_x_T8.paf || rc=$?; rec a_x_T8 $rc
log "2b ALNtoPAF -x -T1"; rc=0; $A2P -x -T1 x.1aln > b_x_T1.paf || rc=$?; rec b_x_T1 $rc
log "2c ALNtoPAF -T8, no CIGAR"; rc=0; $A2P -T8 x.1aln > c_noCIGAR_T8.paf || rc=$?; rec c_noCIGAR_T8 $rc
log "2d ALNtoPAF -x -T8, stack unlimited"; rc=0; ( ulimit -s unlimited; $A2P -x -T8 x.1aln ) > d_x_T8_stack.paf || rc=$?; rec d_x_T8_stack $rc
log "3 FastGA -pafx -T8 -f$F (PAF streamed by FastGA, no ALNtoPAF)"; rc=0; $FG -pafx -T8 -P$W -f$F q.fa $CID.target.fa > e_fastga_pafx.paf || rc=$?; rec e_fastga_pafx $rc
if grep -q -P "^a_x_T8\t0\t" $OUT/diag.tsv; then log "ALNtoPAF passed: x.1aln not kept"; else cp x.1aln $OUT/x.1aln; log "kept $OUT/x.1aln ($(stat -c %s x.1aln) bytes)"; fi
cat $OUT/diag.tsv
log DONE
