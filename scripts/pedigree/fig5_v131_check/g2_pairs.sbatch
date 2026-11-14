#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05
# Homolog vs non-homolog identity per window
# args: WASHU, TASKS (lines: tag pair chrom), OUT_DIR ; sbatch --array=1-N
# Query chromosome (as g2_homolog.sbatch: pat, mat, kid) vs the whole parent diploid with minimap2 2.31
# -cx asm5 --eqx -P: all chains kept, no primary selection, so no chromosome competes with another.
# Per window (10, 20, 50 kb): best identity (=/(=+X+I+D) inside the window, >= 50% of the window aligned)
# to the homolog (same chromosome; chrX also chrY) and to any other chromosome.
# out: OUT_DIR/win/<tag>.<pair>.<chrom>.w<kb>.tsv (start, end, hom_id, hom_target, other_id, other_target), OUT_DIR/paf/
set -euo pipefail
D=$1; T=$2; O=$3; cd $D; mkdir -p $O/{win,paf}
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
read TAG PAIR C < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $T)
case $TAG in v11|v131) N=norm_$TAG;; v11p) N=patch_v11/norm_v11_patched;; v131p) N=patch_v131/norm_v131_patched;; esac
case $PAIR in pat) Q=PAN027#2; P=PAN011;; mat) Q=PAN027#1; P=PAN010;;
  kid) P=PAN027; case " chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX " in *" $C "*) Q=PAN028#2;; *) Q=PAN028#1;; esac;; esac
W=/scratch/$USER.g2p.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
$S faidx $N/$Q.fa.gz "$Q#$C" > $W/q.fa; $S faidx $W/q.fa
zcat $N/$P#1.fa.gz $N/$P#2.fa.gz > $W/t.fa
log "$TAG $PAIR $C: $Q#$C vs $P diploid ($N), host $(hostname)"
$MM -cx asm5 --eqx -P -I 8G -t $SLURM_CPUS_PER_TASK $W/t.fa $W/q.fa 2> $W/mm.log > $W/a.paf
log "records $(wc -l < $W/a.paf), $(grep 'Real time' $W/mm.log)"
gzip -c $W/a.paf > $O/paf/$TAG.$PAIR.$C.paf.gz
python3 - $W/a.paf $C $O/win/$TAG.$PAIR.$C $(cut -f2 $W/q.fa.fai) <<'PY'
import sys, re, collections
paf, c, out, L = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
hom = {c} | ({"chrY"} if c == "chrX" else set())
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
