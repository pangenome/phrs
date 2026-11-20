#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 0-08:00:00
#SBATCH --exclude=octopus05,octopus06
# cenSat BED from v1.1 to v1.3.1
# Same-md5 contigs: copy rows. Changed contigs: minimap2 -cx asm5 (v1.1 target, v1.3.1 query, same contig name),
# primary records only, each interval projected through the cg CIGAR.
# GAP rows on changed contigs: checked = v1.1 N runs, replaced by v1.3.1 N runs.
# args: --washu DIR --md5 TSV (contig_md5.sbatch) --out-dir DIR ; sbatch --array=1-4 (PAN010 PAN011 PAN027 PAN028)
set -euxo pipefail
trap 'echo "exit $? at line ${LINENO}"' EXIT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
while [ $# -gt 0 ]; do case $1 in --washu) D=$2; shift 2;; --md5) M=$2; shift 2;; --out-dir) O=$2; shift 2;; *) echo "unknown arg $1"; exit 1;; esac; done
[ -n "${D:-}" ] && [ -s "${M:-}" ] && [ -n "${O:-}" ] || { echo "need --washu, --md5 (non-empty), --out-dir"; exit 1; }
SM=$(echo PAN010 PAN011 PAN027 PAN028 | cut -d' ' -f$SLURM_ARRAY_TASK_ID); MM=$D/scripts_tmp/minimap2-2.31/minimap2
log() { echo "[$(date -Is)] $*"; }
OUTF=$O/$SM.polished.diploid.cenSat.bed; [ -s $OUTF ] && { log "skip: $OUTF exists"; exit 0; }
cd $D; mkdir -p $O; W=/scratch/$USER.liftcs.${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W; trap 'echo "exit $? at line ${LINENO}"; rm -rf $W' EXIT TERM
awk -F'\t' -v s=$SM 'index($2,s"#")==1 {m[$1,$2]=$4; c[$2]=1} END {for (k in c) print k"\t"(m["v11",k]==m["v131",k] ? "same" : "changed")}' $M | sort > $W/status.tsv
log "$SM: contigs $(wc -l < $W/status.tsv), changed $(awk '$2=="changed"' $W/status.tsv | wc -l)"
for h in 1 2; do
  awk -v p="$SM#$h#" 'index($1,p)==1 && $2=="changed" {print $1}' $W/status.tsv > $W/ch$h.txt
  [ -s $W/ch$h.txt ] || { log "hap $h: no changed contigs"; : > $W/h$h.paf; : > $W/v11.$h.fa; : > $W/v131.$h.fa; continue; }
  $S faidx -r $W/ch$h.txt norm_v11/$SM#$h.fa.gz > $W/v11.$h.fa; $S faidx -r $W/ch$h.txt norm_v131/$SM#$h.fa.gz > $W/v131.$h.fa
  log "hap $h: minimap2 $(wc -l < $W/ch$h.txt) contigs"
  $MM -t 16 -I 8G -cx asm5 $W/v11.$h.fa $W/v131.$h.fa > $W/h$h.paf
  log "hap $h: PAF lines $(wc -l < $W/h$h.paf)"
done
cat $W/h1.paf $W/h2.paf > $W/all.paf
# N runs of changed contigs (GAP rows)
for V in v11 v131; do cat $W/$V.1.fa $W/$V.2.fa | python3 -c '
import sys, re
def runs(n, s):
    for m in re.finditer(r"[Nn]+", s): print(f"{n}\t{m.start()}\t{m.end()}")
n, b = None, []
for l in sys.stdin:
    if l.startswith(">"):
        if n: runs(n, "".join(b))
        n, b = l[1:].split()[0], []
    else: b.append(l.strip())
if n: runs(n, "".join(b))' > $W/n$V.tsv; log "$V N runs on changed contigs: $(wc -l < $W/n$V.tsv)"; done
python3 - $SM annot_v11/$SM.polished.diploid.cenSat.bed $W/status.tsv $W/all.paf $OUTF.tmp $W/nv11.tsv $W/nv131.tsv <<'PY'
import sys, re, collections
sm, bed, st, paf, out, n11, n131 = sys.argv[1:8]
def aname(c):
    s, h, ch = c.split("#")
    return f"{s}.{ch}." + ({"1": "maternal", "2": "paternal"}[h] if s == "PAN027" else f"haplotype{h}")
status = {aname(l.split()[0]): l.split()[1] for l in open(st)}
alns = collections.defaultdict(list)
for l in open(paf):
    t = l.rstrip("\n").split("\t")
    if t[0] != t[5] or "tp:A:P" not in t[12:]: continue
    cg = next(x[5:] for x in t[12:] if x.startswith("cg:Z:"))
    alns[aname(t[5])].append((int(t[7]), int(t[8]), t[4], int(t[2]), int(t[3]), int(t[1]), re.findall(r"(\d+)([MIDNSHX=])", cg)))
def project(a, s, e):
    ts, te, strand, qs, qe, ql, ops = a
    s, e = max(s, ts), min(e, te)
    if s >= e: return None
    t, q, lo, hi = ts, 0, None, None
    for n, op in ops:
        n = int(n)
        if op in "M=X":
            if lo is None and t + n > s: lo = q + max(0, s - t)
            if t + n >= e: hi = q + (e - t); break
            t += n; q += n
        elif op == "D":
            if lo is None and t + n > s: lo = q
            if t + n >= e: hi = q; break
            t += n
        elif op == "I": q += n
    if lo is None or hi is None or hi <= lo: return None
    return (qs + lo, qs + hi) if strand == "+" else (qe - hi, qe - lo)
def nruns(f): return {(aname(a), int(b), int(e)) for a, b, e in (l.split() for l in open(f))}
nr11, nr131 = nruns(n11), nruns(n131)
c = collections.Counter(); bp = collections.Counter(); rows = []; gap11 = set(); gcol = None
for l in open(bed):
    t = l.rstrip("\n").split("\t")
    if len(t) < 3 or not t[1].isdigit(): continue
    s, e = int(t[1]), int(t[2]); k = status.get(t[0])
    if k is None: print(f"ERROR cenSat contig {t[0]} not in md5 table", file=sys.stderr); sys.exit(1)
    if k == "same": rows.append(t); c["same"] += 1; bp["same_in"] += e - s; bp["same_out"] += e - s; continue
    if t[3] == "GAP": gap11.add((t[0], s, e)); gcol = t[8:]; continue  # replaced by v1.3.1 N runs below
    bp["changed_in"] += e - s; hits = [p for p in (project(a, s, e) for a in alns[t[0]]) if p]
    if not hits: c["lost"] += 1; continue
    c["lifted"] += 1
    for a, b in hits: rows.append([t[0], str(a), str(b)] + t[3:6] + [str(a), str(b)] + t[8:]); bp["changed_out"] += b - a
miss, extra = gap11 - nr11, nr11 - gap11
print(f"{sm}: v1.1 GAP rows on changed contigs {len(gap11)}, N runs {len(nr11)}, GAP not N run {len(miss)}, N run not GAP {len(extra)}")
if miss or extra: print(f"ERROR GAP != N runs on v1.1: {sorted(miss)[:5]} {sorted(extra)[:5]}", file=sys.stderr); sys.exit(1)
for a, s, e in sorted(nr131): rows.append([a, str(s), str(e), "GAP", "0", ".", str(s), str(e)] + (gcol or ["0,0,0"])); c["gap131"] += 1
rows.sort(key=lambda r: (r[0], int(r[1])))
with open(out, "w") as o:
    for r in rows: o.write("\t".join(r) + "\n")
print(f"{sm}: intervals same {c['same']}, lifted {c['lifted']}, lost {c['lost']}, GAP from v1.3.1 N runs {c['gap131']}; bp same {bp['same_in']}, changed {bp['changed_in']} -> {bp['changed_out']}")
PY
mv $OUTF.tmp $OUTF
log "done: $OUTF rows $(wc -l < $OUTF)"
