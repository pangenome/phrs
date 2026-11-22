#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05,octopus06
# Does FastGA miss the transmitted parent copy at a chromosome end also without the other copy (Fig5 stop that -f 32..48
# does not fix: SELF.p chr8q, chr13q, SELF.m chr19q, KID.m chr7p)? Same end region aligned by FastGA (Fig5 binary, -f32)
# and minimap2 asm5 against the transmitted copy alone (t1) and against both copies (t12).
# usage: sbatch fastga_iso.sbatch --washu D --query-fa FA --target-fa FA --q REGION --t-tr REGION --t-other REGION --out TSV
#   [--whole-target 1: t12 = the whole --target-fa (FastGA only; minimap2 not run on it)]
# A REGION can be a whole contig name.
# out: TSV with one line per record >= 20 kb: aligner, target set, query start, end, target, target start, end, matches, block
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
WT=""
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --query-fa) QF=$2;; --target-fa) TF=$2;; --q) Q=$2;; --t-tr) T1=$2;; --t-other) T2=$2;; --out) O=$2;; --whole-target) WT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D QF TF Q T1 T2 O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O ] && { log "exists $O, skip"; exit 0; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
B=/home/erikg/.cache/sweepga/0.1.1_fe5037d5e_wd47b7e3e; export PATH=$B:$PATH; T=${SLURM_CPUS_PER_TASK:-4}  # FastGA calls FAtoGDB from PATH
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/fgiso.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log "query $Q from $QF; transmitted $T1, other $T2 from $TF; host $(hostname)"
$S faidx $QF "$Q" > $W/q.fa; $S faidx $TF "$T1" > $W/t1.fa; $S faidx $TF "$T2" > $W/t2.fa; if [ -n "$WT" ]; then cp $TF $W/t12.fa; else cat $W/t1.fa $W/t2.fa > $W/t12.fa; fi; log "t12: $(grep -c '>' $W/t12.fa) contigs"
printf "aligner\ttargets\tqstart\tqend\ttarget\ttstart\ttend\tmatches\tblock\n" > $W/out.tsv
for t in t1 t12; do
  log "FastGA -f32 vs $t"; $B/FastGA -v -T$T -P$W -f32 -paf $W/q.fa $W/$t.fa > $W/fg.$t.paf 2> $W/fg.$t.log || { tail -5 $W/fg.$t.log; exit 1; }
  if [ $t = t12 ] && [ -n "$WT" ]; then : > $W/mm.$t.paf; log "minimap2 vs whole target not run"; else log "minimap2 asm5 vs $t"; $M -cx asm5 -t$T $W/$t.fa $W/q.fa > $W/mm.$t.paf 2> $W/mm.$t.log || { tail -5 $W/mm.$t.log; exit 1; }; fi
  for a in fg mm; do
    log "$a $t: $(wc -l < $W/$a.$t.paf) records, $(awk -F'\t' '$11 >= 20000' $W/$a.$t.paf | wc -l) >= 20 kb"
    awk -F'\t' -v a=$a -v t=$t 'BEGIN {OFS = "\t"} $11 >= 20000 {print a, t, $3, $4, $6, $8, $9, $10, $11}' $W/$a.$t.paf | sort -k3,3n >> $W/out.tsv
  done
done
cp $W/out.tsv $O; cat $O
log done
