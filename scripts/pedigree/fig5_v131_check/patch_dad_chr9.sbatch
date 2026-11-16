#!/bin/bash
#SBATCH -p tux
#SBATCH -w tux05
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 12:00:00
# Patch DAD h2 chr9q end from DAD ONT reads
# args: WASHU, OUT_DIR, [VERSION] (default v11), [FILL_TSV]
# FILL_TSV: check/wrong_positions.patched.chr9q_h2_end.<v>.tsv of patch_read_check; fills the N rows only
# Replaces PAN011#2#chr9 from anchor P to the end with the consensus of DAD (PAN011) ONT reads.
# Reads come from the DAD-on-SELF cross BAMs, where SELF.p chr9q (PAN027#2#chr9) is the alignment frame.
# out: patch/ (segment FASTA, patch.tsv), norm_<v>_patched/ (8 haplotypes, PAN011#2 patched),
#      inputs_<v>_patched/ (Fig5 PAN027pat_vs_PAN011_joint), check/ (assembly and read checks)
set -euo pipefail
D=$1; OUT=$2; V=${3:-v11}; FILL=${4:-}; cd $D; mkdir -p $OUT/{patch,norm_${V}_patched,inputs_${V}_patched,check}
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.patch9.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
DC=PAN011#2#chr9; SC=PAN027#2#chr9
MISJOIN=$([ $V = v11 ] && echo 133673490 || echo 133673964)  # nucflag MISJOIN start (v1.1; v1.3.1 shifted +474)
LD=$(grep -P "^$DC\t" norm_$V/PAN011#2.fa.gz.fai | cut -f2); LS=$(grep -P "^$SC\t" norm_$V/PAN027#2.fa.gz.fai | cut -f2)
log "host $(hostname), $V, $DC length $LD, $SC length $LS, anchor must be <= $MISJOIN"

# 1. anchor: exact DAD=SELF run before misjoin
$S faidx norm_$V/PAN011#2.fa.gz "$DC" > $W/dad9.fa; $S faidx $W/dad9.fa
$S faidx norm_$V/PAN027#2.fa.gz "$SC" > $W/self9.fa; $S faidx $W/self9.fa
$S faidx $W/dad9.fa "$DC:$((LD - 400000 + 1))-$LD" > $W/dad_end.fa
$S faidx $W/self9.fa "$SC:$((LS - 400000 + 1))-$LS" > $W/self_end.fa
$M -cx asm5 --cs=short --secondary=no -t 8 $W/dad_end.fa $W/self_end.fa 2> $W/mm_anchor.log > $W/anchor.paf
log "SELF end vs DAD end chains: $(awk '{printf "q%s-%s t%s-%s NM %s; ", $3,$4,$8,$9,$13}' $W/anchor.paf)"
read P Q RUN < <(python3 - $W/anchor.paf $((LD - 400000)) $((LS - 400000)) $MISJOIN <<'PY'
import re, sys
paf, doff, soff, lim = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
best = None  # (dad_end, self_end, run) of the last exact run >= 10 kb, clipped at lim
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    if f[4] != "+": continue
    q, t = int(f[2]), int(f[7]); cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
    for op, s in re.findall(r"([:*+\-])([0-9]+|[a-z]+)", cs):
        if op == ":":
            n = int(s); ts, te = doff + t, doff + t + n  # 0-based DAD run [ts, te)
            e = min(te, lim)
            if e - ts >= 10000 and (best is None or e > best[0]): best = (e, soff + q + (e - ts), e - ts)
            t += n; q += n
        elif op == "*": t += 1; q += 1
        elif op == "-": t += len(s)
        elif op == "+": q += len(s)
if best is None: sys.exit("no exact run >= 10 kb before the misjoin")
print(*best)
PY
)
[ -n "${P:-}" ] && [ -n "${Q:-}" ] || { log "ERROR: no anchor"; exit 1; }
log "anchor: DAD $DC 1-$P kept, SELF frame $SC from $((Q + 1)), exact run ending at anchor $RUN bp"

# 2. DAD reads on SELF.p chr9q end
R="$SC:$((Q - 50000))-$LS"
ls ont/SRR*.PAN011.onPAN027.$V.bam > $W/cross.txt
$S merge -@ 8 -R "$R" -b $W/cross.txt -o $W/cross_all.bam
$S view -@ 8 -b -q 20 -F 0x704 $W/cross_all.bam -o $W/cross.bam; $S index $W/cross.bam
NR=$($S view -c -F 0x800 $W/cross.bam "$SC:$((Q + 1))-$LS")
log "DAD cross BAMs $(wc -l < $W/cross.txt), reads MAPQ>=20 primary over the frame: $NR"
$S depth -a -Q 20 -r "$SC:$((Q + 1))-$LS" $W/cross.bam | awk '{print $3}' | sort -n > $W/depth.txt
DMIN=$(head -1 $W/depth.txt); DMED=$(awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}' $W/depth.txt); DMAX=$(tail -1 $W/depth.txt)
log "frame depth MAPQ>=20: min $DMIN, median $DMED, max $DMAX over $(wc -l < $W/depth.txt) bp"

# 3. consensus over SELF frame
$S consensus -m simple -c 0.6 -d 5 -H 1.0 --min-MQ 20 -a --show-ins yes --show-del no -l 0 \
  -r "$SC:$((Q + 1))-$LS" $W/cross.bam -o $W/cons.fa
python3 - $W/cons.fa $W/seg.fa <<'PY'
import sys
s = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
lead = len(s) - len(s.lstrip("N")); trail = len(s) - len(s.rstrip("N"))
print(f"consensus {len(s)} bp, N total {s.count('N')}, N leading {lead}, N trailing {trail}", file=sys.stderr)
open(sys.argv[2], "w").write(">seg\n" + s.rstrip("N") + "\n")
PY
if [ -n "$FILL" ]; then  # fill N from DAD reads on patched DAD
  python3 - $W/seg.fa "$FILL" $P <<'PY'
import sys
seg, fill, P = sys.argv[1], sys.argv[2], int(sys.argv[3])
s = list(open(seg).read().split("\n")[1]); n = 0
for l in open(fill):
    f = l.split("\t"); i = int(f[0]) - P - 1
    if f[1] != "N": continue  # only rows at N bases
    if s[i] != "N": sys.exit(f"fill position {f[0]} is not N in the segment")
    s[i] = f[2]; n += 1
open(seg, "w").write(">seg\n" + "".join(s) + "\n"); print(f"filled {n} N from {fill}", file=sys.stderr)
PY
fi
NN=$(grep -v '>' $W/seg.fa | tr -cd 'N' | wc -c); NSEG=$(grep -v '>' $W/seg.fa | tr -d '\n' | wc -c)
log "patch segment $NSEG bp after trimming trailing N, internal N $NN"
[ $NN -le 100 ] || { log "ERROR: $NN N in the patch segment"; exit 1; }

# 4. patched chr9 and hap2
{ echo ">$DC"; $S faidx $W/dad9.fa "$DC:1-$P" | grep -v '>' | tr -d '\n'; grep -v '>' $W/seg.fa | tr -d '\n'; echo; } > $W/dad9_patched.fa
$S faidx $W/dad9_patched.fa; LP=$(cut -f2 $W/dad9_patched.fa.fai)
cut -f1 norm_$V/PAN011#2.fa.gz.fai | awk -v c="$DC" '$1 == c {exit} {print}' > $W/before.txt  # keep contig order
cut -f1 norm_$V/PAN011#2.fa.gz.fai | awk -v c="$DC" 'f {print} $1 == c {f = 1}' > $W/after.txt
{ $S faidx norm_$V/PAN011#2.fa.gz -r $W/before.txt; cat $W/dad9_patched.fa; [ -s $W/after.txt ] && $S faidx norm_$V/PAN011#2.fa.gz -r $W/after.txt; } | $BGZ -@ 8 -c > $OUT/norm_${V}_patched/PAN011#2.fa.gz
$S faidx $OUT/norm_${V}_patched/PAN011#2.fa.gz
cmp -s <(cut -f1 norm_$V/PAN011#2.fa.gz.fai) <(cut -f1 $OUT/norm_${V}_patched/PAN011#2.fa.gz.fai) || { log "ERROR: contig order differs"; exit 1; }
cmp -s <(grep -v -P "^$DC\t" norm_$V/PAN011#2.fa.gz.fai | cut -f1,2) <(grep -v -P "^$DC\t" $OUT/norm_${V}_patched/PAN011#2.fa.gz.fai | cut -f1,2) || { log "ERROR: other contigs changed length"; exit 1; }
log "patched PAN011#2: $(wc -l < $OUT/norm_${V}_patched/PAN011#2.fa.gz.fai) contigs (original $(wc -l < norm_$V/PAN011#2.fa.gz.fai)), chr9 $LD -> $LP bp"
for f in norm_$V/*.fa.gz*; do b=$(basename "$f"); [[ $b == PAN011#2.* ]] || ln -sf "$D/$f" "$OUT/norm_${V}_patched/$b"; done
cp $W/seg.fa $OUT/patch/PAN011_h2_chr9_${V}_patch.fa; sed -i "1s|.*|>${DC}_patch_$((P + 1))-$LP|" $OUT/patch/PAN011_h2_chr9_${V}_patch.fa
printf "version\tcontig\tkept_1_to\treplaced_old_start\treplaced_old_end\tnew_length\tsegment_bp\tsegment_N\tframe\tframe_start\tframe_end\treads_mapq20\tdepth_min\tdepth_median\tdepth_max\tsegment_md5\tpatched_chr9_md5\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  $V "$DC" $P $((P + 1)) $LD $LP $NSEG $NN "$SC" $((Q + 1)) $LS $NR $DMIN $DMED $DMAX \
  $(grep -v '>' $W/seg.fa | tr -d '\n' | md5sum | cut -c1-32) $(grep -v '>' $W/dad9_patched.fa | tr -d '\n' | md5sum | cut -c1-32) > $OUT/patch/patch.$V.tsv
cat $OUT/patch/patch.$V.tsv

# 5. Fig5 inputs, pat comparison only
ln -sf $D/inputs_$V/PAN027pat_vs_PAN011_joint.query.fa $OUT/inputs_${V}_patched/; ln -sf $D/inputs_$V/PAN027pat_vs_PAN011_joint.query.fa.fai $OUT/inputs_${V}_patched/
for h in 1 2; do
  f=norm_$V/PAN011#$h.fa.gz; [ $h = 2 ] && f=$OUT/norm_${V}_patched/PAN011#2.fa.gz
  zcat "$f" | sed "s|^>PAN011#${h}#|>PAN011#joint#h${h}_|"
done > $OUT/inputs_${V}_patched/PAN027pat_vs_PAN011_joint.target.fa
$S faidx $OUT/inputs_${V}_patched/PAN027pat_vs_PAN011_joint.target.fa
log "Fig5 target $(wc -l < $OUT/inputs_${V}_patched/PAN027pat_vs_PAN011_joint.target.fa.fai) contigs, $(awk '{s+=$2} END {print s}' $OUT/inputs_${V}_patched/PAN027pat_vs_PAN011_joint.target.fa.fai) bp"

# 6. assembly checks: segment vs frame, old end, chr1p
$S faidx $W/self9.fa "$SC:$((Q + 1))-$LS" > $W/frame.fa
$S faidx $W/dad9.fa "$DC:$((P + 1))-$LD" > $W/oldend.fa
$S faidx norm_$V/PAN011#1.fa.gz "PAN011#1#chr1:1-200000" > $W/chr1p.fa
cs_count () { python3 -c '
import re, sys
for l in sys.stdin:
    f = l.rstrip("\n").split("\t"); cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
    m = sum(int(s) for o, s in re.findall(r"(:)([0-9]+)", cs)); snv = cs.count("*")
    ins = re.findall(r"\+([a-z]+)", cs); dl = re.findall(r"-([a-z]+)", cs)
    print(f"{sys.argv[1]}\t{f[0]}:{int(f[2])+1}-{f[3]}/{f[1]}\t{f[5]}:{int(f[7])+1}-{f[8]}/{f[6]}\t{f[4]}\tmatch={m}\tSNV={snv}\tINS={len(ins)}({sum(map(len, ins))}bp)\tDEL={len(dl)}({sum(map(len, dl))}bp)")
' "$1"; }
for ref in frame oldend chr1p; do
  $M -cx asm5 --cs=short --secondary=no -t 8 $W/$ref.fa $W/seg.fa 2> /dev/null | cs_count "patch_vs_$ref"
done | tee $OUT/check/patch_vs_refs.$V.tsv
# patch vs SELF frame, every difference
$M -cx asm5 --cs=short --secondary=no -t 8 $W/frame.fa $W/seg.fa 2> /dev/null | python3 -c '
import re, sys
l = sys.stdin.readline().split("\t"); t = int(l[7]) + int(sys.argv[1]); cs = next(x for x in l[12:] if x.startswith("cs:Z:"))[5:]
print("frame_pos\ttype\tframe_base\tpatch_base")
for op, s in re.findall(r"([:*+\-])([0-9]+|[a-z]+)", cs):
    if op == ":": t += int(s)
    elif op == "*":
        k = "N" if s[1] == "n" else "SNV"; print(f"{t + 1}\t{k}\t{s[0].upper()}\t{s[1].upper()}"); t += 1
    elif op == "-": print(f"{t + 1}\tDEL\t{s.upper()}\t-"); t += len(s)
    elif op == "+": print(f"{t + 1}\tINS\t-\t{s.upper()}")
' $Q > $OUT/check/patch_vs_frame_diffs.$V.tsv
log "patch vs SELF frame differences: $(tail -n +2 $OUT/check/patch_vs_frame_diffs.$V.tsv | cut -f2 | sort | uniq -c | tr '\n' ' ')"

log done
