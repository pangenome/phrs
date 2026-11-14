#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# SELF chromosome vs parent homologs
# args: WASHU, TASKS (lines: tag pair chrom), OUT_DIR ; sbatch --array=1-N
# tag: v11, v131, v11p, v131p (patched DAD chr9); pair: pat (PAN027#2 vs PAN011), mat (PAN027#1 vs PAN010),
# kid (KID.m vs PAN027: hap2 for chr1 3 4 5 7 9 12 14 15 18 21 22 X, hap1 for the rest; kid_phase.sbatch).
# Targets: parent h1 and h2 contigs of the same chromosome (chrX also chrY), one alignment per target.
# Tools: wfmash latest (-n 1) and minimap2 2.31 (-cx asm5 --eqx, primary records).
# Coverage = union of aligned intervals, query and target side; uncovered bases outside N runs are gaps.
# out: OUT_DIR/cov/<tag>.<pair>.<chrom>.tsv, OUT_DIR/gaps/<tag>.<pair>.<chrom>.tsv, OUT_DIR/paf/*.paf.gz
set -euo pipefail
D=$1; T=$2; O=$3; cd $D; mkdir -p $O/{cov,gaps,paf}
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; WF=$D/scripts_tmp/latest/wfmash; export LD_LIBRARY_PATH=$D/scripts_tmp/latest/lib; MM=$D/scripts_tmp/minimap2-2.31/minimap2
read TAG PAIR C < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $T)
case $TAG in v11|v131) N=norm_$TAG;; v11p) N=patch_v11/norm_v11_patched;; v131p) N=patch_v131/norm_v131_patched;; esac
case $PAIR in pat) Q=PAN027#2; P=PAN011;; mat) Q=PAN027#1; P=PAN010;;
  kid) P=PAN027; case " chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX " in *" $C "*) Q=PAN028#2;; *) Q=PAN028#1;; esac;; esac
W=/scratch/$USER.g2h.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
$S faidx $N/$Q.fa.gz "$Q#$C" > $W/q.fa; $S faidx $W/q.fa
TG=$(for h in 1 2; do cut -f1 $N/$P#$h.fa.gz.fai | grep -E "#$C\$$( [ $C = chrX ] && echo '|#chrY$')" | sed "s/^/$h:/"; done)
log "$TAG $PAIR $C query $Q#$C, targets: $(echo $TG | tr '\n' ' '), host $(hostname), wfmash $($WF --version 2>&1 | head -1)"
for x in $TG; do h=${x%%:*}; t=${x#*:}; n=$(echo $t | tr '#' '_')
  $S faidx $N/$P#$h.fa.gz "$t" > $W/t.fa; $S faidx $W/t.fa
  $WF -t $SLURM_CPUS_PER_TASK -n 1 --quiet -B $W $W/t.fa $W/q.fa > $W/wf.$n.paf 2> $W/wf.log || { log "ERROR wfmash $t"; tail -5 $W/wf.log; exit 1; }
  $MM -cx asm5 --eqx -t $SLURM_CPUS_PER_TASK $W/t.fa $W/q.fa 2> $W/mm.log | awk -F'\t' '/tp:A:P/' > $W/mm.$n.paf
  log "$t: wfmash $(wc -l < $W/wf.$n.paf) records, minimap2 primary $(wc -l < $W/mm.$n.paf)"
  cp $W/t.fa $W/t.$n.fa
done
python3 - $W $TAG $PAIR $C $O <<'PY'
import sys, glob, re, os
W, tag, pair, c, O = sys.argv[1:6]
def seq(f): return "".join(l.strip() for l in open(f) if not l.startswith(">"))
def nruns(s): return [(m.start(), m.end()) for m in re.finditer(r"[Nn]+", s)]
def union(iv):
    out = []
    for a, b in sorted(map(tuple, iv)):
        if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out
def gaps(L, cov, ns):
    g, p = [], 0
    for a, b in union(cov + ns) + [[L, L]]:
        if a > p: g.append((p, a))
        p = max(p, b)
    return g
q = seq(f"{W}/q.fa"); qn = nruns(q); rows, gl = [], []
for tf in sorted(glob.glob(f"{W}/t.*.fa")):
    n = tf.split("/t.", 1)[1][:-3]; t = seq(tf); tn = nruns(t)
    for tool in ("wf", "mm"):
        r = [l.split("\t") for l in open(f"{W}/{tool}.{n}.paf")]
        qc = union([(int(x[2]), int(x[3])) for x in r]); tc = union([(int(x[7]), int(x[8])) for x in r])
        for side, L, cov, ns in (("query", len(q), qc, qn), ("target", len(t), tc, tn)):
            g = gaps(L, cov, ns); gb = sum(b - a for a, b in g)
            rows.append([tag, pair, c, tool, n, side, L, sum(b - a for a, b in ns), sum(b - a for a, b in cov), gb,
                         sum(1 for a, b in g if b - a >= 10000), max((b - a for a, b in g), default=0),
                         cov[0][0] if cov else -1, L - cov[-1][1] if cov else -1])
            gl += [[tag, pair, c, tool, n, side, a, b, b - a] for a, b in g if b - a >= 1000]
    os.system(f"gzip -c {W}/wf.{n}.paf > {O}/paf/{tag}.{pair}.{c}.{n}.wfmash.paf.gz; gzip -c {W}/mm.{n}.paf > {O}/paf/{tag}.{pair}.{c}.{n}.mm2.paf.gz")
h = "tag pair chrom tool target side length n_bp covered_bp uncovered_nonN_bp gaps_ge10kb largest_gap start_uncovered end_uncovered".split()
with open(f"{O}/cov/{tag}.{pair}.{c}.tsv", "w") as o:
    o.write("\t".join(h) + "\n"); [o.write("\t".join(map(str, x)) + "\n") for x in rows]
with open(f"{O}/gaps/{tag}.{pair}.{c}.tsv", "w") as o:
    o.write("tag\tpair\tchrom\ttool\ttarget\tside\tstart\tend\tlen\n"); [o.write("\t".join(map(str, x)) + "\n") for x in gl]
for x in rows: print("\t".join(map(str, x)))
PY
log done
