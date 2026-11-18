#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# One error locus, HiFi and ONT
# usage: sbatch g1_locus.sbatch --washu D --sample SM --version V --error CONTIG:A-Z --flank BP --out DIR [--extra R1,R2]
#   [--asm NORM_DIR (default norm_V)] [--read-error CONTIG:A-Z on norm_V, where the BAM reads sit (default --error)]
#   [--mmi-dir DIR with SM.hifi.mmi and SM.ont.mmi of the --asm diploid] [--rules ex:L,...,id:T,...]
# Own reads only, same steps for HiFi (hifi/bam) and ONT (ont_full):
# 1. windows: faulty = ERROR +- FLANK, homolog window on the other haplotype (minimap2 asm5).
# 2. reads: primary reads from both windows (+ extra regions).
# 3. span: reads covering ERROR +- 1 kb in the original BAMs (any MAPQ).
# 4. assemble (hifiasm 0.25.0): mode hifi, mode ont (--ont), mode both (HiFi + --ul ONT).
# 5. candidates: per mode and anchor rule (--rules, g1_locus.py junction) the first contig, by aligned bases, that
#    anchors both sides; a piece equal to the assembly is "same"; equal pieces are kept once (cand_N).
# 6. votes: reads (+ reads where a new piece also matches elsewhere) on the diploid (--mmi-dir shared index, built
#    here when not given) and on the candidates; per read best AS on assembly interval, candidate, elsewhere.
# accept (strict rule, others offline from votes.<set>.tsv.gz): more votes for the candidate, HiFi and ONT.
# out: DIR/summary.tsv (one row per mode and rule), DIR/cands.tsv, DIR/pieces.fa, DIR/votes.{hifi,ont}.tsv.gz,
#      DIR/<mode>.contigs.fa
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
EX=""; NA=""; RE=""; MD=""; RU=ex:5000,ex:10000,ex:20000,id:0.99,id:0.995,id:0.999
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --error) ERR=$2;; --flank) F=$2;; --out) O=$2;; --extra) EX=$2;; --asm) NA=$2;; --read-error) RE=$2;; --mmi-dir) MD=$2;; --rules) RU=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
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
log "5 rules $RU"
printf "id\tmode\trule\tcontig\tasm_start\tasm_end\told_bp\tnew_bp\tmd5\n" > $W/cands_all.tsv; : > $W/cand.fa; : > $W/cands.tsv; : > $O/pieces.fa
printf "mode\tcontigs\n" > $W/nc.tsv; : > $W/keys.txt
for m in hifi ont both; do
  mkdir -p $W/$m; cd $W/$m
  case $m in hifi) X="$W/hifi.fq";; ont) X="--ont $W/ont.fq";; both) X="--ul $W/ont.fq $W/hifi.fq";; esac
  $HA -t $T -o loc $X > hifiasm.log 2>&1 || { tail -5 hifiasm.log; echo "ERROR hifiasm $m"; exit 1; }
  for g in loc.bp.hap1.p_ctg loc.bp.hap2.p_ctg loc.bp.p_ctg; do [ -s $g.gfa ] && awk -v p=$m.$g '$1=="S" {print ">"p"."$2"\n"$3}' $g.gfa; done > raw.fa
  cd $D; cp $W/$m/raw.fa $O/$m.contigs.fa; NC=$(awk '/^>/ {n++} END {print n+0}' $W/$m/raw.fa); printf "$m\t$NC\n" >> $W/nc.tsv
  log "4 assemble $m: $NC contigs, $(grep -v '>' $W/$m/raw.fa | tr -d '\n' | wc -c) bp"
  [ $NC -gt 0 ] || continue
  $S faidx $W/$m/raw.fa; $MM -cx asm5 -t $T $W/asm.fa $W/$m/raw.fa 2> /dev/null > $W/$m/o.paf
  awk -F'\t' '{b[$1]+=$10; if ($10 > m[$1]) {m[$1]=$10; s[$1]=$5}} END {for (c in b) print b[c], c, s[c]}' $W/$m/o.paf | sort -k1,1nr > $W/$m/order.txt
  [ -s $W/$m/order.txt ] || continue
  while read -r nb c st; do [ $st = + ] && $S faidx $W/$m/raw.fa "$c" || $S faidx -i $W/$m/raw.fa "$c" | sed "1s/.*/>$c/"; done < $W/$m/order.txt > $W/$m/ctg.fa
  $S faidx $W/$m/ctg.fa; $MM -cx asm5 --eqx -t $T $W/asm.fa $W/$m/ctg.fa 2> /dev/null > $W/$m/a.paf
  : > $W/$m/picks.txt
  while read -r nb c st; do
    grep -P "^$c\t" $W/$m/a.paf > $W/$m/c.paf || continue
    python3 $PY junction $W/$m/c.paf $((WE - WS + 1)) $((A - WS)) $((Z - WS + 1)) $TL $TR $(awk -v c="$c" '$1==c {print $2}' $W/$m/ctg.fa.fai) $RU | awk -v c="$c" '$2 != "none" {print c, $0}' >> $W/$m/picks.txt
  done < $W/$m/order.txt
  awk '!s[$2]++' $W/$m/picks.txt > $W/$m/first.txt
  while read -r c rule lt lq rt rq; do
    if [ "$lt" = same ]; then printf "NA\t$m\t$rule\t$c\tNA\tNA\tNA\tNA\tsame\n" >> $W/cands_all.tsv; continue; fi
    $S faidx $W/$m/ctg.fa "$c:$((lq + 1))-$rq" | tail -n +2 | tr -d '\n' > $W/p.seq
    if [ $lt -lt $rt ]; then $S faidx $W/asm.fa "assembly:$((lt + 1))-$rt" | tail -n +2 | tr -d '\n' > $W/a.seq; else : > $W/a.seq; fi
    PM=$(md5sum < $W/p.seq | cut -c1-12); AM=$(md5sum < $W/a.seq | cut -c1-12)
    if [ $PM = $AM ]; then printf "NA\t$m\t$rule\t$c\t$((WS - 1 + lt))\t$((WS - 1 + rt))\t$((rt - lt))\t$((rq - lq))\tsame\n" >> $W/cands_all.tsv; continue; fi
    K="$lt-$rt-$PM"; ID=$(awk -v k=$K '$1 == k {print $2}' $W/keys.txt)
    if [ -z "$ID" ]; then
      ID=cand_$(( $(wc -l < $W/cands.tsv) + 1 )); echo "$K $ID" >> $W/keys.txt; echo "$ID $lt $rt $((rq - lq))" >> $W/cands.tsv
      { echo ">$ID"; { if [ $lt -gt 0 ]; then $S faidx $W/asm.fa "assembly:1-$lt" | tail -n +2; fi; cat $W/p.seq; echo
        if [ $rt -lt $((WE - WS + 1)) ]; then $S faidx $W/asm.fa "assembly:$((rt + 1))-$((WE - WS + 1))" | tail -n +2; fi; } | tr -d '\n'; echo; } >> $W/cand.fa
      { echo ">$ID $C $((WS - 1 + lt)) $((WS - 1 + rt)) $m $rule $c"; cat $W/p.seq; echo; } >> $O/pieces.fa
    fi
    printf "$ID\t$m\t$rule\t$c\t$((WS - 1 + lt))\t$((WS - 1 + rt))\t$((rt - lt))\t$((rq - lq))\t$PM\n" >> $W/cands_all.tsv
  done < $W/$m/first.txt
  log "5 candidates $m: $(awk -v m=$m '$2 == m' $W/cands_all.tsv | cut -f1,3 | tr '\t\n' ': ')"
done
cp $W/cands_all.tsv $O/cands.tsv
log "5 candidates: $(wc -l < $W/cands.tsv) distinct pieces"
printf "id" > $W/votes.tsv; for s in hifi ont; do printf "\t${s}_votes_asm\t${s}_votes_cand\t${s}_span_asm\t${s}_span_cand" >> $W/votes.tsv; done; printf "\n" >> $W/votes.tsv
if [ -s $W/cands.tsv ]; then
  if [ -z "$MD" ]; then
    MD=$W/mmi; mkdir -p $MD; zcat $NA/$SM#1.fa.gz $NA/$SM#2.fa.gz > $W/dip.fa
    for s in hifi ont; do $MM -t $T -I 8G -x $([ $s = hifi ] && echo map-hifi || echo map-ont) -d $MD/$SM.$s.mmi $W/dip.fa 2> $W/$s.idx.log; done
    log "6 votes: no --mmi-dir, diploid indexes built here"
  fi
  for s in hifi ont; do [ -s $MD/$SM.$s.mmi ] || { echo "ERROR missing $MD/$SM.$s.mmi"; exit 1; }; done
  awk '/^>/ {n=substr($1, 2); next} {q[n]=q[n] $0} END {for (n in q) for (i=1; i<=length(q[n]); i+=10000) print ">"n"_"i"\n"substr(q[n],i,10000)}' $O/pieces.fa > $W/chunks.fa
  $MM -t $T -N 20 -p 0.5 $MD/$SM.hifi.mmi $W/chunks.fa 2> /dev/null | awk -F'\t' -v c="$C" -v ws=$WS -v we=$WE '!($6 == c && $9 > ws && $8 < we) {s=$8-10000; print $6":"(s < 1 ? 1 : s)"-"$9+10000}' | sort -u > $W/hits.txt
  log "6 votes: new pieces also match $(wc -l < $W/hits.txt) places elsewhere; their reads join"
  for s in hifi ont; do
    X=$([ $s = hifi ] && echo map-hifi || echo map-ont)
    [ $s = hifi ] && B=$(ls hifi/bam/*.$SM.$V.bam) || B=$(ls ont_full/SRR*.$SM.$V.bam)
    if [ -s $W/hits.txt ]; then for b in $B; do $S view -F 0x900 $b $(cat $W/hits.txt); done | awk -F'\t' '{print "@"$1"\n"$10"\n+\n"$11}' >> $W/$s.fq; fi
    awk 'NR%4==1 {k=!s[$1]++} k' $W/$s.fq > $W/$s.u.fq
    $MM -ax $X -N 20 -p 0.5 -t $T $MD/$SM.$s.mmi $W/$s.u.fq 2> /dev/null > $W/$s.dip.sam
    $MM -ax $X -N 20 -p 0.5 -t $T $W/cand.fa $W/$s.u.fq 2> /dev/null > $W/$s.cand.sam
    python3 $PY votes $W/$s.dip.sam $W/$s.cand.sam "$C" $((WS - 1)) $W/cands.tsv $W/$s.reads.tsv > $W/$s.votes
    gzip -c $W/$s.reads.tsv > $O/votes.$s.tsv.gz
    log "6 votes $s: $(awk 'NR%4==1' $W/$s.u.fq | wc -l) reads; per candidate votes assembly, candidate, spanning: $(tr '\n' ';' < $W/$s.votes)"
  done
  paste -d' ' $W/hifi.votes $W/ont.votes | awk -v OFS='\t' '$1 != $6 {print "ERROR vote rows differ: "$0 > "/dev/stderr"; exit 1} {print $1, $2, $3, $4, $5, $7, $8, $9, $10}' >> $W/votes.tsv
fi
SP="$SM\t$V\t$ERR\t$(awk 'NR%4==2' $W/hifi.fq | wc -l)\t$(awk 'NR%4==2' $W/ont.fq | wc -l)\t$SPAN_hifi\t$SPAN_ont"
awk -F'\t' -v OFS='\t' -v sp="$SP" 'FILENAME ~ /nc.tsv/ {nc[$1]=$2; next} FILENAME ~ /votes.tsv/ {v[$1]=$0; next}
  FNR==1 {print "sample", "version", "error", "reads_hifi", "reads_ont", "span_orig_hifi", "span_orig_ont", "contigs", $0, "hifi_votes_asm", "hifi_votes_cand", "hifi_span_asm", "hifi_span_cand", "ont_votes_asm", "ont_votes_cand", "ont_span_asm", "ont_span_cand", "accept"; next}
  {x=($1 in v) ? v[$1] : "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"; n=split(x, y, "\t"); acc=($1 in v && y[3] > y[2] && y[7] > y[6]) ? "yes" : "no"; sub(/^[^\t]*\t/, "", x); print sp, nc[$2], $0, x, acc}' $W/nc.tsv $W/votes.tsv $W/cands_all.tsv > $O/summary.tsv
cat $O/summary.tsv
rm -rf "$W"
log done
