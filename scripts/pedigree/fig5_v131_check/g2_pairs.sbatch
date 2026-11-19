#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05
# Homolog vs non-homolog identity per window
# args: WASHU, TASKS (lines: tag pair chrom), OUT_DIR ; sbatch --array=1-N
# Query chromosome (as g2_homolog.sbatch: pat, mat, kid) vs each parent contig, one minimap2 2.31 run per pair
# (-cx asm5 --eqx), pairs in parallel: no chromosome competes with another; -P on the whole diploid was too slow.
# Per window (10, 20, 50 kb): best identity (=/(=+X+I+D) inside the window, >= 50% of the window aligned)
# to the homolog (same chromosome only) and to any other chromosome (chrY included).
# out: OUT_DIR/win/<tag>.<pair>.<chrom>.w<kb>.tsv (start, end, hom_id, hom_target, other_id, other_target), OUT_DIR/paf/
set -euo pipefail
D=$1; T=$2; O=$3; cd $D; mkdir -p $O/{win,paf}
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
read TAG PAIR C < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $T)
case $TAG in v11|v131) N=norm_$TAG;; v11p) N=patch_v11/norm_v11_patched;; v131p) N=patch_v131/norm_v131_patched;; c11) N=g1/corrected_v11_pol/norm;; c11b) N=g1/corrected_v11b_pol/norm;; c131) N=g1/corrected_v131/norm;; g11) N=g1/apply_v11;; *) echo "unknown tag $TAG"; exit 1;; esac
case $PAIR in pat) Q=PAN027#2; P=PAN011;; mat) Q=PAN027#1; P=PAN010;;
  kid) P=PAN027; case " chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX " in *" $C "*) Q=PAN028#2;; *) Q=PAN028#1;; esac;; esac
NQ=$N; NP=$N
if [ $TAG = g11 ]; then NQ=$N/${Q%%#*}/polished2/norm; NP=$N/$P/polished2/norm; fi
W=/scratch/$USER.g2p.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
$S faidx $NQ/$Q.fa.gz "$Q#$C" > $W/q.fa; $S faidx $W/q.fa
log "$TAG $PAIR $C: $Q#$C vs each $P contig ($NQ, $NP), host $(hostname)"
for h in 1 2; do cut -f1 $NP/$P#$h.fa.gz.fai | sed "s/^/$h\t/"; done > $W/targets.txt
export MM S N NP P W
pair () { t=$2; $S faidx $NP/$P#$1.fa.gz "$t" > $W/t.$t.fa && $MM -cx asm5 --eqx -t 1 $W/t.$t.fa $W/q.fa 2> $W/mm.$t.log > $W/a.$t.paf; rm -f $W/t.$t.fa; }
export -f pair
xargs -P $SLURM_CPUS_PER_TASK -L 1 bash -c 'pair "$0" "$1"' < $W/targets.txt
cat $W/a.*.paf > $W/a.paf
log "pairs $(wc -l < $W/targets.txt), records $(wc -l < $W/a.paf), per target: $(cut -f6 $W/a.paf | sort | uniq -c | sort -rn | head -4 | tr -s ' ' | tr '\n' ';')"
gzip -c $W/a.paf > $O/paf/$TAG.$PAIR.$C.paf.gz
python3 - $W/a.paf $C $O/win/$TAG.$PAIR.$C $(cut -f2 $W/q.fa.fai) <<'PY'
import sys, re, collections
paf, c, out, L = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
hom = {c}  # chrY counts as other, so PAR1 on chrX shows as an X-Y exchange
ws = (10000, 20000, 50000)
best = {w: collections.defaultdict(lambda: [0.0, "NA", 0.0, "NA"]) for w in ws}
n = 0
for l in open(paf):
    f = l.rstrip("\n").split("\t"); qs = int(f[2]); t = f[5]; tc = t.split("#")[-1]
    cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
    ops = [(int(a), b) for a, b in re.findall(r"(\d+)([=XIDM])", cg)]
    if f[4] == "-": ops.reverse()  # cigar follows the target; reversed = query order
    k = 1 if tc in hom else 3
    for w in ws:
        acc = collections.defaultdict(lambda: [0, 0, 0])  # window -> [eq, cols, query bp]
        p = qs
        for m, op in ops:
            if op == "D":
                acc[p // w][1] += m; continue
            e = p + m
            while p < e:
                b = p // w; s = min(e, (b + 1) * w) - p
                a = acc[b]; a[1] += s; a[2] += s
                if op == "=": a[0] += s
                p += s
        for b, (eq, cols, qb) in acc.items():
            if qb < 0.5 * min(w, L - b * w): continue
            idn = eq / cols; r = best[w][b]
            if idn > r[k - 1]: r[k - 1], r[k] = idn, t
    n += 1
for w in ws:
    with open(f"{out}.w{w // 1000}.tsv", "w") as o:
        o.write("start\tend\thom_id\thom_target\tother_id\tother_target\n")
        for b in range(0, (L + w - 1) // w):
            r = best[w][b]
            o.write(f"{b * w}\t{min(L, (b + 1) * w)}\t{r[0]:.5f}\t{r[1]}\t{r[2]:.5f}\t{r[3]}\n")
    x = [r for r in best[w].values()]
    print(f"w{w // 1000} kb: windows {(L + w - 1) // w}, with homolog hit {sum(1 for r in x if r[0] > 0)}, other > homolog {sum(1 for r in x if r[2] > r[0])}", file=sys.stderr)
print(f"alignments {n}", file=sys.stderr)
PY
log done
