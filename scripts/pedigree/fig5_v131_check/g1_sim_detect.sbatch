#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05,octopus06
# HiFi scan and candidate groups at the event ends of a sim, reads aligned to the sim (g1_sim_reads.sbatch)
# usage: sbatch g1_sim_detect.sbatch --washu D --sim SIM_DIR --sample SM --version V --bam-dir DIR --merge BP[,BP...] --out DIR
# scan (as wg_scan.sbatch with hifi): the 1 Mb at the contig end of each event in SIM_DIR/truth.tsv, 10 kb bins on the contig
#   grid: median MAPQ>=0 depth (samtools depth -a -Q 0), SNV (bcftools mpileup -B -I -q 20, call --ploidy 1, QUAL>=20, by TYPE).
# flag (g1_candidates_hifi.py): depth <= 0.25 x contig median (median of the norm_V whole-contig scan check/wg_scan/bins_hifi_q0)
#   or SNV >= 3 with depth >= 3; groups: flagged bins <= MERGE bp apart; "satellite" when > 50% cenSat (annot_V, mapped back
#   through the event); only "candidate" groups become loci, as in the real run.
# out: DIR/bins/<contig>.<arm>.tsv; DIR/groups.m<MERGE>.tsv; DIR/tasks.m<MERGE>.tsv (event row, <event>.g<k>.m<MERGE>, group
#   interval as error and read interval: reads sit on the sim); DIR/tasks.truth.tsv (true intervals, reads on the sim);
#   DIR/detect.m<MERGE>.tsv (per event: event bp, event bp inside candidate groups, groups)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sim) SD=$2;; --sample) SM=$2;; --version) V=$2;; --bam-dir) BD=$2;; --merge) MG=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SD SM V BD MG O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
[ -s $BD/$SM.hifi.bam.bai ] || { echo "ERROR missing $BD/$SM.hifi.bam.bai"; exit 1; }
mkdir -p $O/bins /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simdet.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
cat $SD/norm/$SM#1.fa.gz.fai $SD/norm/$SM#2.fa.gz.fai > $W/fai
awk -F'\t' -v OFS='\t' 'NR==FNR {n[$1]=$2; next} FNR>1 {a=($6 < 1e6 ? "p" : "q"); s=(a == "p" ? 0 : int((n[$3] > 1e6 ? n[$3] - 1e6 : 0) / 1e4) * 1e4); print $3, a, s, (a == "p" ? (n[$3] < 1e6 ? n[$3] : 1e6) : n[$3])}' $W/fai $SD/truth.tsv | sort -u > $W/ends.tsv
log "$SM $V sim $SD: $(wc -l < $W/ends.tsv) event ends, BAM $BD/$SM.hifi.bam, merge $MG, host $(hostname)"
while IFS=$'\t' read -r c a s e; do
  F=$O/bins/$(echo $c | tr '#' '_').$a.tsv; [ -s $F ] && { log "exists $F, skip"; continue; }
  H=$(echo $c | cut -d'#' -f2)
  $S depth -a -Q 0 -r "$c:$((s + 1))-$e" $BD/$SM.hifi.bam | awk -F'\t' '{print int(($2 - 1) / 10000)"\t"$3}' > $W/depth.tsv
  $B mpileup -B -I -Ou -f $SD/norm/$SM#$H.fa.gz -r "$c:$((s + 1))-$e" -q 20 --max-depth 2000 -a AD,DP $BD/$SM.hifi.bam 2> $W/mp.log \
    | $B call -m --ploidy 1 -v -Ou 2> $W/call.log | $B view -i 'QUAL>=20' -Ou | $B query -f '%POS\t%TYPE\n' \
    | awk -F'\t' '$2 == "SNP" {s[int(($1 - 1) / 10000)]++} END {for (k in s) print k"\t"s[k]}' > $W/snv.tsv
  python3 -c 'import sys, collections, statistics as st
d, s = collections.defaultdict(list), collections.Counter()
for l in open(sys.argv[2]): b, n = l.split(); s[int(b)] = int(n)
for l in open(sys.argv[1]): b, n = l.split(); d[int(b)].append(int(n))
for b in sorted(d): print(f"{sys.argv[3]}\t{b * 10000}\t{st.median(d[b])}\t{s[b]}")' $W/depth.tsv $W/snv.tsv "$c" > $F
  log "$c $a [$s, $e): $(wc -l < $F) bins, SNV $(awk '{n+=$4} END {print n+0}' $F)"
done < $W/ends.tsv
python3 - $D $SD $SM $V $O $MG $W/ends.tsv <<'PY'
import collections, statistics as st, sys
D, SD, SM, V, O, MGS, EN = sys.argv[1:8]
def aname(c):
    s, h, ch = c.split("#"); return f"{s}.{ch}." + ({"1": "maternal", "2": "paternal"}[h] if s == "PAN027" else f"haplotype{h}")
cs = collections.defaultdict(list)
for l in open(f"{D}/annot_{V}/{SM}.polished.diploid.cenSat.bed"):
    t = l.split("\t")
    if len(t) >= 3 and t[1].isdigit(): cs[t[0]].append((int(t[1]), int(t[2])))
tr = [l.rstrip("\n").split("\t") for l in open(f"{SD}/truth.tsv")][1:]
ends = {(c, a) for c, a, s, e in (l.split("\t") for l in open(EN))}
with open(f"{O}/tasks.truth.tsv", "w") as o:
    for i, r in enumerate(tr, 1):
        ms, me = int(r[5]), int(r[6]); b = 10000; a = ms // b * b + 1; z = -(-me // b) * b if me > ms else a + b - 1
        o.write(f"{i}\t{r[0]}_{r[1]}\t{r[2]}:{a}-{z}\t{r[2]}:{a}-{z}\n")
for MG in map(int, MGS.split(",")):
    G, T, R = [], [], []
    for i, r in enumerate(tr, 1):
        c = r[2]; os_, oe, ms, me = map(int, r[3:7]); a = "p" if ms < 1000000 else "q"
        if (c, a) not in ends: sys.exit(f"ERROR event {i} {c} {a} has no scanned end")
        med = st.median(float(x.split("\t")[5]) for x in open(f"{D}/check/wg_scan/bins_hifi_q0/{SM}.{V}.{c.replace('#', '_')}.bins.tsv"))
        b = [(int(t[1]), float(t[2]), int(t[3])) for t in (x.rstrip("\n").split("\t") for x in open(f"{O}/bins/{c.replace('#', '_')}.{a}.tsv"))]
        fl = [x for x in b if x[1] <= 0.25 * med or (x[2] >= 3 and x[1] >= 3)]; g = []
        for x in fl:
            if g and x[0] - g[-1][-1][0] <= MG: g[-1].append(x)
            else: g.append([x])
        back = lambda x: x if x <= ms else (x - (me - ms) + (oe - os_) if x >= me else os_)  # sim -> norm_V coordinate
        cov, ng = 0, 0
        for k, grp in enumerate(g, 1):
            s, e = grp[0][0], grp[-1][0] + 10000; s0, e0 = back(s), back(e)
            sat = sum(max(0, min(e0, y) - max(s0, x)) for x, y in cs[aname(c)]) / max(1, e0 - s0); cl = "satellite" if sat > 0.5 else "candidate"
            ov = max(0, min(e, me) - max(s, ms)); G.append((i, f"{r[0]}_{r[1]}", c, s, e, len(grp), min(x[1] for x in grp), med, sum(x[2] for x in grp), ov, f"{sat:.2f}", cl))
            if cl == "candidate":
                cov += ov; ng += 1; T.append(f"{i}\t{r[0]}_{r[1]}.g{k}.m{MG}\t{c}:{s + 1}-{e}\t{c}:{s + 1}-{e}\n")
        R.append((i, f"{r[0]}_{r[1]}", c, a, me - ms, cov, ng))
    with open(f"{O}/groups.m{MG}.tsv", "w") as o:
        o.write("event_row\tevent\tcontig\tstart\tend\tbins\thifi_depth_min\thifi_depth_contig_median\thifi_snv\tevent_bp_inside\tcensat_frac\tclass\n")
        for x in G: o.write("\t".join(map(str, x)) + "\n")
    open(f"{O}/tasks.m{MG}.tsv", "w").write("".join(T))
    with open(f"{O}/detect.m{MG}.tsv", "w") as o:
        o.write("event_row\tevent\tcontig\tarm\tevent_bp\tevent_bp_in_candidates\tcandidate_groups\n")
        for x in R: o.write("\t".join(map(str, x)) + "\n")
    print(f"merge {MG}: groups {len(G)} (candidate {len(T)}); events with >= 1 bp found {sum(1 for x in R if x[5] > 0)} of {len(R)}; "
          + "; ".join(f"{x[1]} {x[5]}/{x[4]} bp in {x[6]} groups" for x in R), file=sys.stderr)
PY
log done
