#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
# Read support per high-SNV cluster
# args: WASHU, BINS_DIR, OUT_DIR ; sbatch --array=1-8 (PAN010 PAN011 PAN027 PAN028 x v11 v131)
# Clusters: wg_scan bins with >= 10 SNV, merged when <= 20 kb apart. Per cluster: SNV sites (bcftools
# mpileup -X ont, call --ploidy 1, QUAL >= 20) and read fraction against the assembly = AD[1] / (AD[0] + AD[1])
# summed over BAMs (counts only, no alleles). Fraction near 1: assembly error; near 0.5: collapsed copies.
# out: OUT_DIR/<sample>.<v>.clusters.tsv
set -euo pipefail
D=$1; BD=$2; O=$3; cd $D; mkdir -p $O
B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
SM=$(echo PAN010 PAN011 PAN027 PAN028 | cut -d' ' -f$(( (SLURM_ARRAY_TASK_ID - 1) / 2 + 1 ))); V=$( [ $((SLURM_ARRAY_TASK_ID % 2)) -eq 1 ] && echo v11 || echo v131)
W=/scratch/$USER.wgcl.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
ls ont_full/SRR*.$SM.$V.bam > $W/bams.txt; log "$SM $V BAMs $(wc -l < $W/bams.txt), host $(hostname)"
python3 - $BD $SM $V $W/clusters.bed $W/ctx.tsv <<'PY'
import sys, glob, statistics as st
bd, sm, v, bed, ctx = sys.argv[1:6]
with open(bed, "w") as o, open(ctx, "w") as x:
    for f in sorted(glob.glob(f"{bd}/{sm}.{v}.*.bins.tsv")):
        rows = [l.rstrip("\n").split("\t") for l in open(f)]
        L = int(rows[-1][1]) + 10000; med = st.median(float(r[2]) for r in rows); cur = None
        for c, b, d, s, i in rows + [[rows[0][0], str(10**12), "0", "0", "0"]]:
            b, s = int(b), int(s)
            if s >= 10 and cur and b - cur[2] <= 20000: cur[2] = b + 10000; cur[3] += s; cur[4].append(float(d)); continue
            if cur:
                o.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\n")
                x.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\t{cur[3]}\t{st.median(cur[4])}\t{med}\t{min(cur[1], L - cur[2])}\t{L}\n")
            cur = [c, b, b + 10000, s, [float(d)]] if s >= 10 else None
PY
log "clusters $(wc -l < $W/clusters.bed), bp $(awk '{s+=$3-$2} END {print s}' $W/clusters.bed)"
cut -f1 $W/clusters.bed | sort -u | while read -r c; do H=${c%#chr*}
  awk -v c="$c" '$1==c' $W/clusters.bed > $W/r.bed
  $B mpileup -X ont -Ou -f norm_$V/$H.fa.gz -R $W/r.bed -q 20 --max-depth 2000 -a AD,DP -b $W/bams.txt 2>> $W/mp.log \
    | $B call -m --ploidy 1 -v -Ou 2>> $W/call.log | $B view -i 'QUAL>=20' -v snps -Ou | $B query -f '%CHROM\t%POS[\t%AD]\n'
done > $W/sites.tsv
log "SNV sites $(wc -l < $W/sites.tsv)"
python3 - $W/ctx.tsv $W/sites.tsv $O/$SM.$V.clusters.tsv <<'PY'
import sys, bisect, collections as C, statistics as st
ctx, sites, out = sys.argv[1:4]
cl = C.defaultdict(list)
for l in open(ctx): r = l.rstrip("\n").split("\t"); cl[r[0]].append(r)
starts = {c: [int(r[1]) for r in v] for c, v in cl.items()}; fr = C.defaultdict(list)
for l in open(sites):
    r = l.rstrip("\n").split("\t"); c, p = r[0], int(r[1]) - 1; ref = alt = 0
    for a in r[2:]:
        x = a.split(","); ref += int(x[0]) if x[0] != "." else 0; alt += sum(int(y) for y in x[1:] if y != ".")
    k = bisect.bisect_right(starts[c], p) - 1
    if k >= 0 and p < int(cl[c][k][2]) and ref + alt > 0: fr[(c, k)].append(alt / (ref + alt))
with open(out, "w") as o:
    o.write("contig\tstart\tend\tlen\tsnv_bins_sum\tdepth_med\tcontig_depth_med\tdist_to_end\tcontig_len\tsites\tfrac_med\tfrac_ge_0.8\tfrac_0.3_0.7\tclass\n")
    for c, v in cl.items():
        for k, r in enumerate(v):
            f = fr[(c, k)]; n = len(f); L = int(r[2]) - int(r[1])
            m = st.median(f) if n else float("nan"); hi = sum(x >= 0.8 for x in f) / n if n else 0; mid = sum(0.3 <= x <= 0.7 for x in f) / n if n else 0
            cls = "satellite_block" if L > 500000 else ("error_like" if hi >= 0.5 else ("collapse_like" if mid >= 0.5 else "mixed"))
            o.write("\t".join(r[:3] + [str(L)] + r[3:8] + [str(n), f"{m:.3f}", f"{hi:.3f}", f"{mid:.3f}", cls]) + "\n")
PY
log "done: $(awk -F'\t' 'NR>1 {n[$14]++} END {for (k in n) printf "%s %d ", k, n[k]}' $O/$SM.$V.clusters.tsv)"
