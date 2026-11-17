#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# One error locus, HiFi and ONT
# usage: sbatch g1_locus.sbatch --washu D --sample SM --version V --error CONTIG:A-Z --flank BP --out DIR [--extra R1,R2]
#   [--asm NORM_DIR (default norm_V)] [--read-error CONTIG:A-Z on norm_V, where the BAM reads sit (default --error)]
# Own reads only, same steps for HiFi (hifi/bam) and ONT (ont_full):
# 1. windows: faulty = ERROR +- FLANK, homolog window on the other haplotype (minimap2 asm5).
# 2. reads: primary reads from both windows (+ extra regions).
# 3. span: reads covering ERROR +- 1 kb in the original BAMs (any MAPQ).
# 4. assemble (hifiasm 0.25.0): mode hifi, mode ont (--ont), mode both (HiFi + --ul ONT).
# 5. candidate per mode: longest contig anchored by exact runs >= 20 kb on both sides of ERROR (g1_locus.py junction),
#    or same_as_assembly when one exact run covers the error.
# 6. compete: reads (+ reads where the new sequence also matches) on the whole diploid + all candidates; per read, best
#    alignment score on assembly interval vs candidate vs elsewhere (g1_locus.py compete).
# accept: more reads vote for the candidate than for the assembly, for HiFi and for ONT. merqury comes after.
# out: DIR/summary.tsv (one row per mode), DIR/<mode>.contigs.fa, DIR/<mode>.piece.fa
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
EX=""; NA=""; RE=""
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --error) ERR=$2;; --flank) F=$2;; --out) O=$2;; --extra) EX=$2;; --asm) NA=$2;; --read-error) RE=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V ERR F O; do [ -n "${!v:-}" ] || { echo "ERROR missing argument for $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O/summary.tsv ] && { log "exists $O/summary.tsv, skip"; exit 0; }
mkdir -p $O; W=/scratch/$USER/tmp/phrs/locus.${SLURM_JOB_ID:-$$}; mkdir -p $W
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
HA=$D/scripts_tmp/hifiasm-0.25.0/hifiasm; PY=$D/scripts_tmp/g1_locus.py; T=${SLURM_CPUS_PER_TASK:-4}
C=${ERR%%:*}; R=${ERR#*:}; A=${R%-*}; Z=${R#*-}; H=$(echo $C | cut -d'#' -f2); OH=$((3 - H)); OC=${C/\#$H\#/\#$OH\#}
NA=${NA:-norm_$V}; RE=${RE:-$ERR}; RR=${RE#*:}; RA=${RR%-*}; RZ=${RR#*-}; [ "${RE%%:*}" = "$C" ] || { echo "ERROR --read-error contig ${RE%%:*} is not $C"; exit 2; }
L=$(awk -v c="$C" '$1==c {print $2}' $NA/$SM#$H.fa.gz.fai); [ -n "$L" ] || { echo "ERROR $C not in $NA"; exit 1; }
RL=$(awk -v c="$C" '$1==c {print $2}' norm_$V/$SM#$H.fa.gz.fai); RWS=$((RA > F ? RA - F : 1)); RWE=$((RZ + F > RL ? RL : RZ + F))
WS=$((A > F ? A - F : 1)); WE=$((Z + F > L ? L : Z + F)); TL=$((WS == 1 ? 1 : 0)); TR=$((WE == L ? 1 : 0))
log "1 windows: $SM $V asm $NA error $ERR (reads at $RE), faulty $C:$WS-$WE (contig $L bp, terminal left $TL right $TR), host $(hostname)"
$S faidx $NA/$SM#$H.fa.gz "$C:$WS-$WE" | sed '1s/.*/>assembly/' > $W/asm.fa; $S faidx $W/asm.fa
$S faidx $NA/$SM#$OH.fa.gz "$OC" > $W/oc.fa
$MM -cx asm5 -t $T $W/oc.fa $W/asm.fa 2> /dev/null | sort -k10,10nr > $W/proj.paf
read OS OE < <(awk -F'\t' 'NR==1 {c=($8+$9)/2} $12 >= 5 && ($8+$9)/2 > c - 1e6 && ($8+$9)/2 < c + 1e6 {s = (s == "" || $8 < s) ? $8 : s; e = ($9 > e) ? $9 : e} END {print s, e}' $W/proj.paf)
[ -n "${OS:-}" ] || { echo "ERROR no homolog window on $OC"; exit 1; }
log "1 windows: homolog $OC:$((OS + 1))-$OE"
for s in hifi ont; do
  [ $s = hifi ] && B=$(ls hifi/bam/*.$SM.$V.bam) || B=$(ls ont_full/SRR*.$SM.$V.bam)
  for b in $B; do $S view -F 0x900 $b "$C:$RWS-$RWE" "$OC:$((OS + 1))-$OE" $(echo $EX | tr ',' ' '); done \
    | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/$s.fq
  N=$(awk 'NR%4==2' $W/$s.fq | wc -l); [ $N -gt 0 ] || { echo "ERROR no $s reads"; exit 1; }
  log "2 reads $s: $N reads, $(awk 'NR%4==2 {s+=length($0)} END {print s}' $W/$s.fq) bp from $(echo $B | wc -w) BAMs"
  SP=$(for b in $B; do $S view -F 0x904 $b "$C:$RA-$RZ"; done | awk -F'\t' -v a=$((RA - 1000)) -v z=$((RZ + 1000 > RL ? RL : RZ + 1000)) \
    '{l=0; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {l+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; if ($4 <= a && $4 + l - 1 >= z) n++} END {print n+0}')
  eval SPAN_$s=$SP; log "3 span $s: $SP reads cover $C:$((RA - 1000))-$((RZ + 1000 > RL ? RL : RZ + 1000)) (read coordinates)"
done
printf "sample\tversion\terror\tmode\treads_hifi\treads_ont\tspan_orig_hifi\tspan_orig_ont\tcontigs\tresult\tcandidate\tasm_start\tasm_end\told_bp\tnew_bp\n" > $W/modes.tsv
: > $W/cand.fa; : > $W/cands.tsv; : > $W/pieces.fa
for m in hifi ont both; do
  mkdir -p $W/$m; cd $W/$m
  case $m in hifi) X="$W/hifi.fq";; ont) X="--ont $W/ont.fq";; both) X="--ul $W/ont.fq $W/hifi.fq";; esac
  $HA -t $T -o loc $X > hifiasm.log 2>&1 || { tail -5 hifiasm.log; echo "ERROR hifiasm $m"; exit 1; }
  for g in loc.bp.hap1.p_ctg loc.bp.hap2.p_ctg loc.bp.p_ctg; do [ -s $g.gfa ] && awk -v p=$m.$g '$1=="S" {print ">"p"."$2"\n"$3}' $g.gfa; done > raw.fa
  cd $D; cp $W/$m/raw.fa $O/$m.contigs.fa; NC=$(awk '/^>/ {n++} END {print n+0}' $W/$m/raw.fa)
  log "4 assemble $m: $NC contigs, $(grep -v '>' $W/$m/raw.fa | tr -d '\n' | wc -c) bp"
  ROW="$SM\t$V\t$ERR\t$m\t$(awk 'NR%4==2' $W/hifi.fq | wc -l)\t$(awk 'NR%4==2' $W/ont.fq | wc -l)\t$SPAN_hifi\t$SPAN_ont\t$NC"; RES=none
  if [ $NC -gt 0 ]; then
    $S faidx $W/$m/raw.fa; $MM -cx asm5 -t $T $W/asm.fa $W/$m/raw.fa 2> /dev/null > $W/$m/o.paf
    awk -F'\t' '{b[$1]+=$10; if ($10 > m[$1]) {m[$1]=$10; s[$1]=$5}} END {for (c in b) print b[c], c, s[c]}' $W/$m/o.paf | sort -k1,1nr > $W/$m/order.txt
    while read -r nb c st; do [ $st = + ] && $S faidx $W/$m/raw.fa "$c" || $S faidx -i $W/$m/raw.fa "$c" | sed "1s/.*/>$c/"; done < $W/$m/order.txt > $W/$m/ctg.fa
    if [ -s $W/$m/ctg.fa ]; then
      $S faidx $W/$m/ctg.fa; $MM -cx asm5 --eqx -t $T $W/asm.fa $W/$m/ctg.fa 2> /dev/null > $W/$m/a.paf
      while read -r nb c st; do
        grep -P "^$c\t" $W/$m/a.paf > $W/$m/c.paf || continue
        QL=$(awk -v c="$c" '$1==c {print $2}' $W/$m/ctg.fa.fai); rc=0
        J=$(python3 $PY junction $W/$m/c.paf $((WE - WS + 1)) $((A - WS)) $((Z - WS + 1)) $TL $TR $QL) || rc=$?
        case $rc in 0) read LT LQ RT RQ <<< "$J"; RES=candidate; CAND=$c; break;; 4) RES=same_as_assembly; CAND=$c; break;; 3) ;; *) echo "ERROR junction rc $rc"; exit 1;; esac
      done < $W/$m/order.txt
    fi
  fi
  if [ $RES = candidate ]; then
    NB=$((RQ - LQ)); log "5 candidate $m: $CAND, assembly $C:$((WS + LT))-$((WS - 1 + RT)) ($((RT - LT)) bp) -> contig $LQ-$RQ ($NB bp)"
    { echo ">cand_$m"; { if [ $LT -gt 0 ]; then $S faidx $W/asm.fa "assembly:1-$LT" | tail -n +2; fi; $S faidx $W/$m/ctg.fa "$CAND:$((LQ + 1))-$RQ" | tail -n +2
      if [ $RT -lt $((WE - WS + 1)) ]; then $S faidx $W/asm.fa "assembly:$((RT + 1))-$((WE - WS + 1))" | tail -n +2; fi; } | tr -d '\n'; echo; } >> $W/cand.fa
    $S faidx $W/$m/ctg.fa "$CAND:$((LQ + 1))-$RQ" | sed "1s/.*/>$m $C $((WS - 1 + LT)) $((WS - 1 + RT)) $CAND/" | tee $O/$m.piece.fa >> $W/pieces.fa
    echo "$m $LT $RT $NB" >> $W/cands.tsv; printf "$ROW\t$RES\t$CAND\t$((WS - 1 + LT))\t$((WS - 1 + RT))\t$((RT - LT))\t$NB\n" >> $W/modes.tsv
  else
    log "5 candidate $m: $RES ${CAND:-}"; printf "$ROW\t$RES\t${CAND:-NA}\tNA\tNA\tNA\tNA\n" >> $W/modes.tsv
  fi
  CAND=""
done
printf "mode" > $W/votes.tsv; for s in hifi ont; do printf "\t${s}_votes_asm\t${s}_votes_cand\t${s}_span_asm\t${s}_span_cand" >> $W/votes.tsv; done; printf "\taccept\n" >> $W/votes.tsv
if [ -s $W/cands.tsv ]; then
  { zcat $NA/$SM#1.fa.gz $NA/$SM#2.fa.gz; cat $W/cand.fa; } > $W/ref.fa
  for s in hifi ont; do $MM -t $T -I 8G -x $([ $s = hifi ] && echo map-hifi || echo map-ont) -d $W/$s.mmi $W/ref.fa 2> $W/$s.idx.log; done
  log "6 compete: indexes of diploid + $(wc -l < $W/cands.tsv) candidates built"
  awk '/^>/ {n=$1; next} {q[n]=q[n] $0} END {for (n in q) for (i=1; i<=length(q[n]); i+=10000) print n"_"i"\n"substr(q[n],i,10000)}' $W/pieces.fa > $W/chunks.fa
  $MM -t $T -N 20 -p 0.5 $W/hifi.mmi $W/chunks.fa 2> /dev/null | awk -F'\t' -v c="$C" -v ws=$WS -v we=$WE '$6 !~ /^cand_/ && !($6 == c && $9 > ws && $8 < we) {s=$8-10000; print $6":"(s < 1 ? 1 : s)"-"$9+10000}' | sort -u > $W/hits.txt
  log "6 compete: new sequence also matches $(wc -l < $W/hits.txt) places elsewhere; their reads join"
  for s in hifi ont; do
    [ $s = hifi ] && B=$(ls hifi/bam/*.$SM.$V.bam) || B=$(ls ont_full/SRR*.$SM.$V.bam)
    if [ -s $W/hits.txt ]; then for b in $B; do $S view -F 0x900 $b $(cat $W/hits.txt); done | awk -F'\t' '{print "@"$1"\n"$10"\n+\n"$11}' >> $W/$s.fq; fi
    awk 'NR%4==1 {k=!s[$1]++} k' $W/$s.fq > $W/$s.u.fq
    $MM -ax $([ $s = hifi ] && echo map-hifi || echo map-ont) -N 20 -p 0.5 -t $T $W/$s.mmi $W/$s.u.fq 2> /dev/null > $W/$s.sam
    python3 $PY compete $W/$s.sam "$C" $((WS - 1)) $W/cands.tsv > $W/$s.votes
    log "6 compete $s: $(awk 'NR%4==1' $W/$s.u.fq | wc -l) reads; per mode votes assembly, candidate, spanning: $(tr '\n' ';' < $W/$s.votes)"
  done
  paste -d' ' $W/hifi.votes $W/ont.votes | awk -v OFS='\t' '$1 != $6 {print "ERROR vote rows differ: "$0 > "/dev/stderr"; exit 1} {print $1, $2, $3, $4, $5, $7, $8, $9, $10, ($3 > $2 && $8 > $7) ? "yes" : "no"}' >> $W/votes.tsv
fi
awk -F'\t' -v OFS='\t' 'NR==FNR {v[$1]=$0; next} {k=(FNR==1 ? "mode" : $4); x=(k in v ? v[k] : "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno"); sub(/^[^\t]*\t/, "", x); print $0, x}' $W/votes.tsv $W/modes.tsv > $O/summary.tsv
cat $O/summary.tsv
rm -rf "$W"
log done
