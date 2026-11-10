#!/bin/bash
#SBATCH -p workers
#SBATCH --exclude=octopus05,octopus06,octopus08,octopus09
#SBATCH -c 24
#SBATCH --mem=64G
#SBATCH -t 12:00:00
# Read-vs-assembly scan of all chromosome ends
# args: WASHU, OUT_DIR, SAMPLE (PAN011, PAN027, PAN028), [VERSION] (default v11), [END_BP] (default 500000)
# Own ONT reads on the own diploid. Per 10 kb bin of the first and last END_BP of every contig:
# wrong positions (>= 3 MAPQ>=20 reads show a base, >= 80% share one base that is not the assembly base),
# median depth of MAPQ>=20 reads and of all reads.
# out: end_scan.<sample>.<v>.tsv (bins), end_scan.<sample>.<v>.wrong.tsv (positions)
set -euo pipefail
D=$1; OUT=$2; SM=$3; V=${4:-v11}; E=${5:-500000}; cd $D; mkdir -p $OUT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.endscan.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
ls ont/SRR*.$SM.$V.bam > $W/bams.txt
log "host $(hostname), $SM $V, BAMs $(wc -l < $W/bams.txt), end window $E bp"
cat > $W/bin.py <<'PY'
import re, sys, statistics as st, collections
# stdin: mpileup of all BAMs; argv: contig, end label, region start, region end, depth file, out prefix
ctg, lab, s, e, dep, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
q20 = collections.defaultdict(list); wrong = collections.Counter(); wl = []
for l in sys.stdin:
    f = l.rstrip("\n").split("\t"); p = int(f[1]); b = re.sub(r"\^.|\$", "", "".join(f[4::3])); i, x = 0, []
    while i < len(b):
        c = b[i]
        if c in "+-":
            m = re.match(r"[+-](\d+)", b[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        x.append(c); i += 1
    nref = sum(c in ".," for c in x); alts = {a: sum(c.upper() == a for c in x) for a in "ACGT"}
    a = max(alts, key=alts.get); k = (p - 1) // 10000; q20[k].append(len(x))
    if f[2].upper() != "N" and nref + alts[a] >= 3 and alts[a] / (nref + alts[a]) >= 0.8:
        wrong[k] += 1; wl.append(f"{ctg}\t{p}\t{f[2].upper()}\t{a}\t{alts[a]}\t{nref}\n")
allq = collections.defaultdict(list)
for l in open(dep):
    f = l.split("\t"); allq[(int(f[1]) - 1) // 10000].append(sum(map(int, f[2:])))
with open(out + ".bins", "w") as o:
    for k in range((s - 1) // 10000, (e - 1) // 10000 + 1):
        n = min(e, (k + 1) * 10000) - max(s, k * 10000 + 1) + 1; d20 = q20.get(k, []) + [0] * (n - len(q20.get(k, [])))
        o.write(f"{ctg}\t{lab}\t{k * 10000 + 1}\t{wrong[k]}\t{st.median(d20)}\t{st.median(allq.get(k, [0]))}\n")
open(out + ".wrong", "w").writelines(wl)
PY
# one job per contig end
while read -r C L _; do
  H=${C%#chr*}
  echo "$C $H start 1 $((E < L ? E : L))"; echo "$C $H end $((L - E + 1 > 1 ? L - E + 1 : 1)) $L"
done < <(cat norm_$V/$SM#1.fa.gz.fai norm_$V/$SM#2.fa.gz.fai) > $W/jobs.txt
log "contig ends: $(wc -l < $W/jobs.txt)"
export S W V
one () {
  C=$1; H=$2; LAB=$3; A=$4; B=$5; o="$W/$(echo "$C" | tr '#' '_').$LAB"
  $S depth -a -r "$C:$A-$B" -f $W/bams.txt > $o.depth
  $S mpileup -B -Q 0 -q 20 -d 100000 -f norm_$V/$H.fa.gz -b $W/bams.txt -r "$C:$A-$B" 2> $o.log | python3 $W/bin.py "$C" $LAB $A $B $o.depth $o
  rm -f $o.depth
}
export -f one
xargs -a $W/jobs.txt -P $SLURM_CPUS_PER_TASK -L 1 bash -c 'one "$@"' _
[ $(ls $W/*.bins | wc -l) -eq $(wc -l < $W/jobs.txt) ] || { log "ERROR: $(ls $W/*.bins | wc -l) of $(wc -l < $W/jobs.txt) ends done"; exit 1; }
{ printf "contig\tend\tbin_start\twrong_positions\tdepth_mapq20_median\tdepth_all_median\n"; cat $W/*.bins | sort -k1,1V -k3,3n; } > $OUT/end_scan.$SM.$V.tsv
{ printf "contig\tpos\tasm_base\tread_base\tn_read_base\tn_asm_base\n"; cat $W/*.wrong | sort -k1,1V -k2,2n; } > $OUT/end_scan.$SM.$V.wrong.tsv
log "bins $(($(wc -l < $OUT/end_scan.$SM.$V.tsv) - 1)), wrong positions $(($(wc -l < $OUT/end_scan.$SM.$V.wrong.tsv) - 1)), bins with >= 5 wrong: $(awk -F'\t' 'NR > 1 && $4 >= 5' $OUT/end_scan.$SM.$V.tsv | wc -l)"
log "top bins: $(sort -t$'\t' -k4,4nr $OUT/end_scan.$SM.$V.tsv | head -8 | awk -F'\t' '{printf "%s:%s(%s wrong, depth %s/%s); ", $1, $3, $4, $5, $6}')"
log done
