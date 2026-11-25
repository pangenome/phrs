#!/usr/bin/env bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Fig5 step 4 with the query cut into pieces: FastGA stops the transmitted-copy chain before chromosome ends when the query is
# a whole chromosome (chr8q: stops at 144,312,322 even vs DAD h1 alone; a 2 Mb end vs the whole genome reaches the end).
# usage: sbatch fig5_align_chunk.sbatch --cid CID --in IN_DIR --out OUT_DIR --freq F --piece BP [--batches K]
# --batches K: query contigs split into K groups (whole contigs, balanced by bp), one sweepga run per group, PAFs joined
# (SELF.p failed 4 times in sweepga's ALNtoPAF on the whole query: 2879922, 2879985, 2880081; K=1 is the plain run).
# A failing group is run again one contig at a time; finished groups saved in OUT_DIR/raw_paf/<name>.parts/ and reused by a rerun.
# Same sweepga call as fig5_align_f.sbatch. Each query contig is cut into ceil(L/BP) equal pieces (each <= BP, no overlap),
# named CONTIG:A-B (samtools faidx); after alignment the query name, length and coordinates go back to the contig.
# An alignment that crosses a piece boundary becomes two records. Every query contig >= 10 Mb must be aligned over >= 50%.
# out: OUT_DIR/raw_paf/CID.frequency<F>.piece<BP>[.batches<K>].paf.gz and .coverage.tsv
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
K=1
while [ $# -gt 0 ]; do case $1 in --cid) CID=$2;; --in) IN=$2;; --out) OUT=$2;; --freq) F=$2;; --piece) P=$2;; --batches) K=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in CID IN OUT F P; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
OUT=$(realpath -m "$OUT")  # absolute: the script cd's into the work dir before the final cp (2879810 lost 5 h on a relative OUT)
SFX=piece$P; [ $K -gt 1 ] && SFX=$SFX.batches$K
R=$OUT/raw_paf/$CID.frequency$F.$SFX.paf.gz
[ -s $R ] && { log "exists $R, skip"; exit 0; }
SWEEPGA=/home/erikg/.cargo/bin/sweepga; PIGZ=/usr/bin/pigz; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
mkdir -p /scratch/$USER/tmp/phrs "$OUT/raw_paf"; W=$(mktemp -d /scratch/$USER/tmp/phrs/fig5chunk.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log "host=$(hostname) cid=$CID f=$F piece=$P in=$IN out=$OUT work=$W"
$SWEEPGA --version
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$W/"; cd "$W"
awk -F'\t' -v p=$P '{n = int(($2 + p - 1) / p); s = int(($2 + n - 1) / n); for (a = 1; a <= $2; a += s) {b = a + s - 1; if (b > $2) b = $2; print $1 ":" a "-" b}}' $CID.query.fa.fai > pieces.txt
log "query contigs $(wc -l < $CID.query.fa.fai), pieces $(wc -l < pieces.txt), longest piece $(awk -F'[:-]' '{d = $NF - $(NF-1) + 1; if (d > m) m = d} END {print m}' pieces.txt) bp"
$S faidx -r pieces.txt $CID.query.fa > q.fa
[ $(grep -c '>' q.fa) -eq $(wc -l < pieces.txt) ] || { echo "ERROR piece count differs from pieces.txt"; exit 1; }
sort -k2,2nr $CID.query.fa.fai | awk -F'\t' -v k=$K '{g = 0; for (i = 1; i < k; i++) if (b[i] < b[g]) g = i; b[g] += $2; print $1 "\t" g}' > groups.tsv
# finished groups are kept in PD and reused when the same command runs again (KID 2880973, 2881075 lost 5-6 h each)
PD=$OUT/raw_paf/$CID.frequency$F.$SFX.parts; mkdir -p $PD
if [ -s $PD/groups.tsv ]; then cmp groups.tsv $PD/groups.tsv || { echo "ERROR groups differ from $PD/groups.tsv"; exit 1; }; else cp groups.tsv $PD/; fi
run () {  # one sweepga run: TAG QUERY_FA OUT_PAF; logs exit status and peak RSS of FastGA/ALNtoPAF; returns the sweepga status
  ( set +ex; while sleep 15; do ps -C FastGA,ALNtoPAF -o rss=,comm=; done ) > mem.$1.txt 2> /dev/null & local M=$! rc=0  # ps exits 1 when neither runs
  $SWEEPGA --fastga --fastga-frequency $F --num-mappings many:many --scaffold-jump 0 --temp-dir "$W" --output-file $3 $2 "$CID.target.fa" || rc=$?
  kill $M; log "$1 sweepga exit $rc, peak RSS: $(peak mem.$1.txt)"; return $rc
}
peak () { awk '{if ($1 > m[$2]) m[$2] = $1} END {for (p in m) printf "%s %.1f GB ", p, m[p] / 1048576}' $1; }
: > pieces.paf
for g in $(seq 0 $((K - 1))); do
  if [ -s $PD/g$g.paf.gz ]; then gzip -t $PD/g$g.paf.gz; zcat $PD/g$g.paf.gz >> pieces.paf; log "group $g reused from $PD ($(zcat $PD/g$g.paf.gz | wc -l) lines)"; continue; fi
  awk -F'\t' -v g=$g 'NR == FNR {if ($2 == g) c[$1] = 1; next} {i = match($0, /:[0-9]+-[0-9]+$/); if (substr($0, 1, i - 1) in c) print}' groups.tsv pieces.txt > p.$g.txt
  $S faidx -r p.$g.txt q.fa > q.$g.fa
  log "sweepga fastga f$F many:many, group $g of $K: $(cut -f2 groups.tsv | grep -cx $g) contigs ($(awk -F'\t' -v g=$g '$2 == g {printf $1 " "}' groups.tsv)), $(wc -l < p.$g.txt) pieces, $(grep -v '>' q.$g.fa | tr -d '\n' | wc -c) bp"
  # a failing group fails again on every try (KID K=8 group 5, chr6 chr11 chr18: 4 of 4, ALNtoPAF code None, RSS 2.4 GB of 503 GB;
  # the same contigs passed inside other groups with K=4): no retry, run the group's contigs one at a time instead
  if ! run g$g q.$g.fa pieces.$g.paf; then
    log "group $g FAILED; node memory GB: $(free -g | awk '/Mem/ {print "total " $2 ", used " $3 ", available " $7}'); its contigs one at a time"
    : > pieces.$g.paf
    for c in $(awk -F'\t' -v g=$g '$2 == g {print $1}' groups.tsv); do
      awk -v c="$c" '{i = match($0, /:[0-9]+-[0-9]+$/); if (substr($0, 1, i - 1) == c) print}' p.$g.txt > p.$g.c.txt
      $S faidx -r p.$g.c.txt q.fa > q.$g.c.fa
      run "g$g.$c" q.$g.c.fa pieces.$g.c.paf || { echo "ERROR group $g contig $c failed alone ($(wc -l < p.$g.c.txt) pieces)"; exit 1; }
      log "group $g contig $c: $(wc -l < p.$g.c.txt) pieces, PAF lines $(wc -l < pieces.$g.c.paf)"; cat pieces.$g.c.paf >> pieces.$g.paf; rm -f pieces.$g.c.paf q.$g.c.fa
    done
  fi
  $PIGZ -p $SLURM_CPUS_PER_TASK -c pieces.$g.paf > $PD/g$g.paf.gz.part; gzip -t $PD/g$g.paf.gz.part; mv $PD/g$g.paf.gz.part $PD/g$g.paf.gz
  log "group $g PAF lines $(wc -l < pieces.$g.paf), saved $PD/g$g.paf.gz"; cat pieces.$g.paf >> pieces.paf; rm -f pieces.$g.paf q.$g.fa
done
log "piece PAF lines=$(wc -l < pieces.paf)"
awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {L[$1] = $2; next} {i = match($1, /:[0-9]+-[0-9]+$/); if (!i) {print "ERROR query name without piece range: " $1 > "/dev/stderr"; exit 1} c = substr($1, 1, i - 1); split(substr($1, i + 1), r, "-"); if (!(c in L)) {print "ERROR unknown contig " c > "/dev/stderr"; exit 1} o = r[1] - 1; $1 = c; $2 = L[c]; $3 += o; $4 += o; print}' $CID.query.fa.fai pieces.paf > $CID.raw.paf
N=$(wc -l < $CID.raw.paf); [ $N -eq $(wc -l < pieces.paf) ] || { echo "ERROR lift lost lines"; exit 1; }
log "lifted PAF lines=$N; coverage check, contig length aligned_union fraction"
cut -f1,3,4 "$CID.raw.paf" | sort -S 8G -k1,1 -k2,2n \
  | awk -F'\t' 'NR==FNR {len[$1]=$2; next} {if ($1!=c) {if (c!="") u[c]=t; c=$1; t=0; e=0} if ($3>e) {s=($2>e)?$2:e; t+=$3-s; e=$3}} END {u[c]=t; bad=0; for (k in len) {f=(k in u)?u[k]/len[k]:0; printf "%s\t%d\t%d\t%.3f\n", k, len[k], u[k]+0, f; if (len[k]>=10000000 && f<0.5) bad=1}; exit bad}' \
    "$CID.query.fa.fai" - | sort -V | tee cov.tsv \
  || { log "ERROR contig coverage below 0.5 (see table above)"; exit 1; }
$PIGZ -p $SLURM_CPUS_PER_TASK -c "$CID.raw.paf" > out.paf.gz; gzip -t out.paf.gz
cp out.paf.gz $R; cp cov.tsv "$OUT/raw_paf/$CID.frequency$F.$SFX.coverage.tsv"
log "DONE $R lines=$N"
