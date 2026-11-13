#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
# Whole-contig read-vs-assembly scan, one contig per task
# args: WASHU, TASKS (lines: sample version contig), OUT_DIR ; sbatch --array=1-N%K
# Own ONT reads (ont_full/SRR*.<sample>.<version>.bam) on the own assembly, whole contig, 10 kb bins:
# median depth of MAPQ>=20 reads; positions where reads disagree with the assembly = bcftools mpileup -X ont
# (no indel calling: default took 1,679 s per 2 Mb, ont 20 s, same SNVs; wg_speed.sbatch), call --ploidy 1 -v,
# QUAL>=20, counted by TYPE only (no alleles).
# out: OUT_DIR/<sample>.<version>.<contig>.bins.tsv (contig, bin_start, depth_q20_median, snv, indel)
set -euo pipefail
D=$1; T=$2; O=$3; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
read SM V C < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $T); H=${C%#chr*}
OUTF=$O/$SM.$V.$(echo $C | tr '#' '_').bins.tsv
W=/scratch/$USER.wgscan.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
ls ont_full/SRR*.$SM.$V.bam > $W/bams.txt
L=$(awk -v c="$C" '$1==c {print $2}' norm_$V/$H.fa.gz.fai)
echo "[$(date -Is)] task $SLURM_ARRAY_TASK_ID $SM $V $C length $L, BAMs $(wc -l < $W/bams.txt), host $(hostname)"
$S depth -a -Q 20 -r "$C" -f $W/bams.txt | awk -F'\t' '{d=0; for (i=3;i<=NF;i++) d+=$i; print int(($2-1)/10000)"\t"d}' > $W/depth.tsv
$B mpileup -X ont -Ou -f norm_$V/$H.fa.gz -r "$C" -q 20 -B --max-depth 2000 -a AD,DP -b $W/bams.txt 2> $W/mp.log \
  | $B call -m --ploidy 1 -v -Ou 2> $W/call.log | $B view -i 'QUAL>=20' -Ou | $B query -f '%POS\t%TYPE\n' \
  | awk -F'\t' '{b=int(($1-1)/10000); if ($2=="SNP") s[b]++; else i[b]++} END {for (k in s) print k"\tsnv\t"s[k]; for (k in i) print k"\tindel\t"i[k]}' > $W/calls.tsv
python3 - "$C" $W/depth.tsv $W/calls.tsv "$OUTF" <<'PY'
import sys, statistics as st, collections
c, dep, cal, out = sys.argv[1:5]
d = collections.defaultdict(list)
for l in open(dep):
    b, v = l.split(); d[int(b)].append(int(v))
s, i = collections.Counter(), collections.Counter()
for l in open(cal):
    b, t, n = l.split(); (s if t == "snv" else i)[int(b)] = int(n)
with open(out, "w") as o:
    for b in sorted(d):
        o.write(f"{c}\t{b * 10000}\t{st.median(d[b])}\t{s[b]}\t{i[b]}\n")
PY
echo "[$(date -Is)] done: bins $(wc -l < $OUTF), snv $(awk '{s+=$4} END {print s+0}' $OUTF), indel $(awk '{s+=$5} END {print s+0}' $OUTF)"
