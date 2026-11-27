#!/bin/bash
#SBATCH -t 5-00:00:00
# Rerun of the paper's HG002 Pore-C or CiFi contact pipeline (code_PHR-and-3D.md, "CiFi pipeline (HG002)" and "Pore-C
# pipeline (HG002)") with two parse2 branches on the same annotated alignments: the paper's call (pairtools parse2
# default --min-mapq 1, unique alignments only) and the same call with --min-mapq 0 (multi-mapping alignments kept), for
# referee R1-M4 / R3-M3. pore-c-py annotate --stdout is split by tee into two FIFOs, each read on stdin by its own
# parse2 | restrict | sort, so both branches see the identical stream and the default branch must reproduce the paper's
# matrix. Pore-C: flow cells 1 and 2 one after the other (digest --threads 32, minimap2 -t 64, sort --nproc 16), then
# pairtools merge --nproc 16 | dedup per branch. CiFi: digest N/3, minimap2 2N/3, sort N/3, then dedup per branch.
# Then cooler cload pairs 5 kb + zoomify --balance per branch. Tools from the rebuilt paper env
# (condatools/porec/1.0-rebuilt, same package list as porec/1.0). Work in node /scratch/$USER/tmp/phrs; outputs to OUT
# (mq1/ and mq0/: mcool and parse2 stats). --head-reads N runs on the first N reads of each input (test).
# pore-c-py digest writes one warning per monomer ("has no modified bases") unless --quiet; with --quiet its records are
# identical (samtools view, 2026-11-27) and only the @PG line differs, which samtools fastq drops. Its stderr goes to
# node /scratch and only counts and the last lines reach the job log (the 2026-11-26 run wrote a 113 GB log on moosefs,
# which kept only the first 28.9 GB).
# The exit code of every stage of the main pipeline is logged.
# usage: sbatch -c 96|48 porec_mapq_rerun.sbatch --assay porec|cifi --out DIR [--head-reads N]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
HR=0
while [ $# -gt 0 ]; do case $1 in --assay) AS=$2;; --out) O=$2;; --head-reads) HR=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in AS O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
export PATH=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin:$PATH
B=/moosefs/guarracino/HPRCv2; H=$B/PHR_III/HiC/HG002; REF=$H/HG002.fa.gz; CS=$H/chrom.sizes; N=${SLURM_CPUS_PER_TASK:?}
case $AS in
  porec) EZ=NlaIII; PR=map-ont; DT=32; MT=64; ST=16
         IN="$B/data/PoreC/HG002/HG002_1_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz $B/data/PoreC/HG002/HG002_2_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz";;
  cifi)  EZ=DpnII; PR=map-hifi; DT=$((N / 3)); MT=$((N * 2 / 3)); ST=$((N / 3)); IN=$B/data/CiFi/HG002/HG002_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam;;
  *) echo "ERROR assay $AS"; exit 2;;
esac
FR=$H/fragments_$EZ.bed
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }
mkdir -p $O/mq1 $O/mq0 /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/rerun_$AS.XXXX)
trap 'rc=$?; [ -z "$(jobs -p)" ] || kill $(jobs -p); rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); $N cores; assay $AS ($EZ, $PR); digest $DT, minimap2 $MT, sort $ST threads; head reads $HR; inputs $IN; $(df -h /scratch | tail -1)"
log "tools: $(pore-c-py --version 2>&1 | tail -1); minimap2 $(minimap2 --version); pairtools $(pairtools --version 2>&1 | awk '{print $NF}'); $(samtools --version | head -1); cooler $(cooler --version 2>&1 | awk '{print $NF}')"
MQ_1=""; MQ_0="--min-mapq 0"   # branch 1 = the paper's call (parse2 default --min-mapq 1)
i=0
for f in $IN; do
  i=$((i + 1)); src=$f
  if [ $HR -gt 0 ]; then
    # test only: the reader stops early, so the upstream gets SIGPIPE; pipefail is off inside these subshells
    case $f in *.bam) ( set +o pipefail; samtools view -h $f | awk -v n=$HR '/^@/ {print; next} {if (++k > n) exit; print}' ) | samtools view -b -o $W/head$i.bam - ; src=$W/head$i.bam;;
                *) ( set +o pipefail; zcat $f | head -n $((4 * HR)) ) | gzip -1 > $W/head$i.fastq.gz; src=$W/head$i.fastq.gz;; esac
    log "input $i: first $HR reads of $f -> $src ($(stat -c %s $src) bytes)"
  fi
  pids=""
  for mq in 1 0; do
    mkfifo $W/ann$i.mq$mq; v=MQ_$mq
    ( pairtools parse2 --single-end --readid-transform 'readID.split(":")[0]' --drop-sam --drop-seq --expand ${!v} \
          --output-stats $W/parse2_stats_$i.mq$mq.txt -c $CS - < $W/ann$i.mq$mq \
        | pairtools restrict -f $FR \
        | pairtools sort --nproc $ST --tmpdir $W -o $W/in$i.mq$mq.sorted.pairs.gz ) &
    pids="$pids $!"
  done
  log "input $i: digest | fastq | minimap2 | annotate | tee into 2 parse2 branches (pids$pids)"
  pore-c-py digest $src $EZ --threads $DT --max_monomers 250 --quiet 2> $W/digest.$i.err \
    | samtools fastq -T '*' \
    | minimap2 -ay -t $MT -x $PR --cap-kalloc 100m --cap-sw-mem 50m $REF - 2> $W/minimap2.$i.err \
    | pore-c-py annotate - $W/annotated_$i --monomers --stdout 2> $W/annotate.$i.err \
    | tee $W/ann$i.mq1 > $W/ann$i.mq0 || { st="${PIPESTATUS[*]}"; log "input $i: main pipeline failed, exit codes digest fastq minimap2 annotate tee: $st"; for e in digest minimap2 annotate; do log "$e stderr, last lines:"; awk '!/has no modified bases/' $W/$e.$i.err | tail -20; done; exit 1; }
  log "input $i: main pipeline done, exit codes ${PIPESTATUS[*]}; digest warnings 'has no modified bases' $(grep -c 'has no modified bases' $W/digest.$i.err), other digest lines $(grep -vc 'has no modified bases' $W/digest.$i.err)"
  awk '!/has no modified bases/' $W/digest.$i.err | tail -20; tail -5 $W/minimap2.$i.err; tail -5 $W/annotate.$i.err
  for p in $pids; do wait $p || { log "input $i: a parse2 branch (pid $p) failed"; exit 1; }; done
  if [ $HR -gt 0 ] && [ $AS = cifi ]; then   # test: the plain paper chain on the same reads must give the same pairs as branch 1
    pore-c-py digest $src $EZ --threads $DT --max_monomers 250 | samtools fastq -T '*' \
      | minimap2 -ay -t $MT -x $PR --cap-kalloc 100m --cap-sw-mem 50m $REF - | pore-c-py annotate - $W/annotated_paper --monomers --stdout \
      | pairtools parse2 --single-end --readid-transform 'readID.split(":")[0]' --drop-sam --drop-seq --expand --output-stats $W/paper_stats.txt -c $CS - \
      | pairtools restrict -f $FR | pairtools sort --nproc $ST --tmpdir $W | pairtools dedup --output $W/paper.dedup.pairs
    pairtools dedup --output $W/test.dedup.pairs $W/in$i.mq1.sorted.pairs.gz
    cmp <(grep -v '^#' $W/paper.dedup.pairs) <(grep -v '^#' $W/test.dedup.pairs)
    log "test: paper chain and branch mq1 give identical pairs ($(grep -c '^[^#]' $W/paper.dedup.pairs))"
  fi
  for mq in 1 0; do log "input $i branch mq$mq: $(zcat $W/in$i.mq$mq.sorted.pairs.gz | grep -c '^[^#]') sorted pairs; $(grep -E '^(total|total_mapped|total_dups)\s' $W/parse2_stats_$i.mq$mq.txt | tr '\n\t' '  ')"; done
  rm -f $W/annotated_$i.* $W/head$i.*
done
for mq in 1 0; do
  if [ $AS = porec ]; then
    pairtools merge $W/in1.mq$mq.sorted.pairs.gz $W/in2.mq$mq.sorted.pairs.gz --nproc 16 --tmpdir $W | pairtools dedup --output $W/dedup.mq$mq.pairs
  else
    pairtools dedup --output $W/dedup.mq$mq.pairs $W/in1.mq$mq.sorted.pairs.gz
  fi
  rm -f $W/in*.mq$mq.sorted.pairs.gz
  log "branch mq$mq: $(grep -c '^[^#]' $W/dedup.mq$mq.pairs) pairs after dedup"
  cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $CS:5000 $W/dedup.mq$mq.pairs $W/$AS.5000.cool
  cooler zoomify $W/$AS.5000.cool -r 5000,10000,20000,50000,100000 --balance -o $W/$AS.mcool
  cp $W/$AS.mcool $W/parse2_stats_*.mq$mq.txt $O/mq$mq/
  rm -f $W/dedup.mq$mq.pairs $W/$AS.5000.cool $W/$AS.mcool
  log "branch mq$mq: $O/mq$mq/$AS.mcool ($(stat -c %s $O/mq$mq/$AS.mcool) bytes)"
done
log "done"
