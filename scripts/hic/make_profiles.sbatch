#!/bin/bash
#SBATCH -c 16
#SBATCH -t 12:00:00
# Error profiles of real reads (read_profile.py) for null simulations that look like the real data (R1-M4 / R3-M3).
# Reads are sampled along the whole of each file (every K-th read; the first reads of a file carry about twice the
# mismatches of the following ones) and renamed f<file>_<index> for the position table.
# hic: every K-th read of each R1 file and of its R2 file (FILES lists R1 files; R2 = same name with _R1 -> _R2),
#      mapped single-end with bowtie2 from the paper's HiC-Pro image and the paper's HiC-Pro global options on the
#      paper's index; read_profile.py hic.
# porec: every K-th read of each file, digested and mapped with the paper's chain (pore-c-py digest NlaIII
#      --max_monomers 250 | samtools fastq -T '*' | minimap2 -ay -x map-ont --cap-kalloc 100m --cap-sw-mem 50m
#      HG002.fa.gz); read_profile.py porec.
# Also writes OUT.real_rates.tsv (the rates of the sampled real reads, as hic_error_rate / porec_error_rate compute them)
# and OUT.position.tsv (mismatch or error rate per aligned base by tenth of each file).
# usage: sbatch make_profiles.sbatch --assay hic|porec --sample HG002|CHM13 --files "F1 F2 ..." --every K --out TSV_GZ
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --assay) AS=$2;; --sample) SM=$2;; --files) FS=$2;; --every) K=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in AS SM FS K O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
E=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; H=/moosefs/guarracino/HPRCv2/PHR_III/HiC/$SM
SIF=/moosefs/guarracino/sifs/hicpro-latest.simg; HP=/usr/local/anaconda/envs/HiC-Pro_v3.0.0/bin; N=${SLURM_CPUS_PER_TASK:?}
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $(dirname $O) /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/prof.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); assay $AS; sample $SM; every ${K}th read of each of: $FS"
# the reader stops early, so the upstream gets SIGPIPE: pipefail is off inside these subshells
if [ $AS = hic ]; then
  for m in 1 2; do : > $W/r$m.fq; i=0; for f in $FS; do i=$((i + 1)); g=$f; [ $m = 2 ] && g=${f/_R1/_R2}; [ -e $g ] || { echo "ERROR missing $g"; exit 1; }
    zcat $g | awk -v k=$K -v f=$i 'NR % 4 == 1 {r++; keep = (r % k == 0)} keep {print (NR % 4 == 1 ? "@f" f "_" r : $0)} END {print "file " f ": " r " reads, kept " int(r / k) > "/dev/stderr"}' >> $W/r$m.fq; done
    log "mate $m: $(( $(wc -l < $W/r$m.fq) / 4 )) reads"
    singularity exec --bind $W:$W --bind /moosefs:/moosefs $SIF $HP/bowtie2 --very-sensitive-local -L 30 --score-min G,20,8 --local --reorder \
        -p $N -x $H/bowtie2_index/$SM -U $W/r$m.fq -S $W/r$m.sam 2> $W/r$m.bt2.log
    log "mate $m bowtie2: $(grep overall $W/r$m.bt2.log)"; done
  $E/python3 $Q/read_profile.py hic $W/r1.sam $W/r2.sam $W/p.tsv.gz
  printf "mate\tmin_mapq\treads\taligned_bases\tmismatch_per_base\tgap_ext_per_base\tsoftclip_per_base\tfrac_0mm\tfrac_1mm\tfrac_2mm\tfrac_3plus\n" > $W/rr.tsv
  for m in 1 2; do for q in 30 0; do grep -v '^@' $W/r$m.sam | awk -F'\t' -v t=$m -v q=$q 'BEGIN{OFS="\t"} $3 != "*" && !($2 % 512 >= 256) && $2 < 2048 && $5 >= q {m=0; s=0; c=$6; while (match(c, /[0-9]+[MIDNSHP=X]/)) {n=substr(c, RSTART, RLENGTH-1); op=substr(c, RSTART+RLENGTH-1, 1); if (op=="M") m+=n; if (op=="S") s+=n; c=substr(c, RSTART+RLENGTH)}
      xm=0; xg=0; for (i=12; i<=NF; i++) {if ($i ~ /^XM:i:/) xm=substr($i, 6); if ($i ~ /^XG:i:/) xg=substr($i, 6)}
      r++; M+=m; S+=s; L+=length($10); X+=xm; G+=xg; k[(xm>=3)?3:xm]++} END{print t, q, r, M, X/M, G/M, S/L, k[0]/r, k[1]/r, k[2]/r, k[3]/r}' >> $W/rr.tsv; done; done
  $E/python3 - $W/p.tsv.gz $W/pos.tsv <<'PY'
import collections, gzip, sys
R = collections.defaultdict(list)
with gzip.open(sys.argv[1], "rt") as f:
    f.readline()
    for l in f:
        m, L, q, c5, c3, mm, nm = l.rstrip("\n").split("\t"); fi, idx = nm[1:].split("_")
        R[(m, fi)].append((int(idx), 0 if mm == "." else len(mm.split(",")), int(L) - int(c5) - int(c3), int(c5) + int(c3), int(L)))
with open(sys.argv[2], "w") as o:
    o.write("mate\tfile\ttenth\treads\tmismatch_per_aligned_base\tsoftclip_per_base\n")
    for (m, fi), v in sorted(R.items()):
        v.sort()
        for t in range(10):
            w = v[len(v) * t // 10:len(v) * (t + 1) // 10]
            o.write(f"{m}\t{fi}\t{t + 1}\t{len(w)}\t{sum(x[1] for x in w) / sum(x[2] for x in w):.5f}\t{sum(x[3] for x in w) / sum(x[4] for x in w):.5f}\n")
PY
else
  export PATH=$E:$PATH; : > $W/in.fq; i=0
  for f in $FS; do i=$((i + 1)); zcat $f | awk -v k=$K -v f=$i 'NR % 4 == 1 {r++; keep = (r % k == 0)} keep {print (NR % 4 == 1 ? "@f" f "_" r : $0)} END {print "file " f ": " r " reads, kept " int(r / k) > "/dev/stderr"}' >> $W/in.fq; done
  log "$(( $(wc -l < $W/in.fq) / 4 )) reads"
  pore-c-py digest $W/in.fq NlaIII --threads 4 --max_monomers 250 | samtools fastq -T '*' \
    | minimap2 -ay -t $N -x map-ont --cap-kalloc 100m --cap-sw-mem 50m $H/HG002.fa.gz - 2> $W/mm2.log | samtools view -b -o $W/m.bam -
  python3 $Q/read_profile.py porec $W/m.bam $W/p.tsv.gz
  printf "set\tmin_mapq\tmonomers\tmapped\tfrac_mapped\taligned_bases\tdivergence\tnm_per_base\tsoftclip_per_base\n" > $W/rr.tsv
  for q in 30 0; do samtools view -F 0x900 $W/m.bam | awk -v t=real -v q=$q 'BEGIN{OFS="\t"} {n++; if ($3 == "*" || $5 < q) next; r++; m=0; s=0; c=$6
      while (match(c, /[0-9]+[MIDNSHP=X]/)) {x=substr(c, RSTART, RLENGTH-1); op=substr(c, RSTART+RLENGTH-1, 1); if (op=="M") m+=x; if (op=="S") s+=x; c=substr(c, RSTART+RLENGTH)}
      for (i=12; i<=NF; i++) {if ($i ~ /^NM:i:/) nm=substr($i, 6); if ($i ~ /^de:f:/) de=substr($i, 6)}
      M+=m; S+=s; L+=length($10); NMs+=nm; D+=de*m} END{print t, q, n, r, r/n, M, D/M, NMs/M, S/L}' >> $W/rr.tsv; done
  python3 - $W/p.tsv.gz $W/pos.tsv <<'PY'
import collections, gzip, sys
R = collections.defaultdict(list)
with gzip.open(sys.argv[1], "rt") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t"); fi, idx = x[0][1:].split("_"); a = int(x[1])
        R[fi].append((int(idx), a, a * (float(x[2]) + float(x[3]) + float(x[4]))))
with open(sys.argv[2], "w") as o:
    o.write("file\ttenth\treads\terror_per_aligned_base\n")
    for fi, v in sorted(R.items()):
        v.sort()
        for t in range(10):
            w = v[len(v) * t // 10:len(v) * (t + 1) // 10]
            o.write(f"{fi}\t{t + 1}\t{len(w)}\t{sum(x[2] for x in w) / sum(x[1] for x in w):.5f}\n")
PY
fi
cp $W/rr.tsv ${O%.tsv.gz}.real_rates.tsv; cp $W/pos.tsv ${O%.tsv.gz}.position.tsv
cat ${O%.tsv.gz}.real_rates.tsv ${O%.tsv.gz}.position.tsv
cp $W/p.tsv.gz $O.part && mv $O.part $O
log "done: $O ($(stat -c %s $O) bytes)"
