#!/bin/bash
#SBATCH -p workers
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# Parent ONT reads at child diagnostic sites
# args: WASHU, OUT_TSV, VERSION, CHILD_REGION, PARENT_REGION, DIAG_TSV, BAM_GLOB
# DIAG_TSV from diag_sites.sbatch (child coords); BAM_GLOB = parent reads on parent assembly
set -euo pipefail
WASHU=$1; OUT=$2; V=$3; CREG=$4; PREG=$5; DIAG=$6; GLOB=$7
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.psites.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
fa () { local r=$1; local h=${r%%#chr*}; $SAMTOOLS faidx "$WASHU/norm_$V/$h.fa.gz" "$r"; }

# Child tract on parent copy, cs tags
fa "$CREG" > $W/c.fa; fa "$PREG" > $W/p.fa
$MM2 -x asm20 -c --cs=short --secondary=no -t 8 $W/p.fa $W/c.fa > $W/cp.paf 2> $W/mm2.log
log "child->parent alignments $(wc -l < $W/cp.paf): $(awk '{printf "%s-%s(%s) t%s-%s id=%.5f; ", $3,$4,$5,$8,$9,$10/$11}' $W/cp.paf)"

# Map child SNV sites to parent coords
python3 - $W/cp.paf "$CREG" "$PREG" "$DIAG" $W/p.fa > $W/sites.tsv <<'PY'
import re, sys
paf, creg, preg, diag, pfa = sys.argv[1:6]
pseq = "".join(l.strip() for l in open(pfa) if not l.startswith(">")).upper()
coff = int(creg.split(":")[1].split("-")[0]) - 1
pctg, prng = preg.rsplit(":", 1); poff = int(prng.split("-")[0]) - 1
want = {}
for l in open(diag):
    f = l.rstrip("\n").split("\t")
    if f[0] == "version": continue
    want[int(f[1])] = f[2]  # child pos -> query_matches (donor/own/neither)
comp = str.maketrans("acgt", "tgca")
n = 0
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    if f[4] != "+": sys.stderr.write("skip minus block\n"); continue
    cs = next(t for t in f[12:] if t.startswith("cs:Z:"))[5:]
    q, t = int(f[2]), int(f[7])
    for op, s in re.findall(r"([:*+\-])([0-9]+|[a-z]+)", cs):
        if op == ":":  # matches: parent base = child base
            for k in range(int(s)):
                cpos = coff + q + k + 1
                if cpos in want:
                    pb = pseq[t + k]; print(f"{pctg}\t{poff + t + k + 1}\t{cpos}\t{want[cpos]}\t{pb}\t{pb}"); n += 1
            q += int(s); t += int(s)
        elif op == "*":
            cpos = coff + q + 1
            if cpos in want:  # parent base s[0], child base s[1]
                print(f"{pctg}\t{poff + t + 1}\t{cpos}\t{want[cpos]}\t{s[0].upper()}\t{s[1].upper()}"); n += 1
            q += 1; t += 1
        elif op == "+": q += len(s)
        elif op == "-": t += len(s)
sys.stderr.write(f"sites mapped {n} of {len(want)}\n")
PY
log "sites $(wc -l < $W/sites.tsv) (child sites on the parent copy; '=' = parent base equals child base), MINQ=${MINQ:-0}"

# Pileup parent reads at parent coords
ls $GLOB > $W/bams.txt; log "BAMs $(wc -l < $W/bams.txt)"
cut -f1,2 $W/sites.tsv > $W/pos.tsv
PCTG=${PREG%:*}
$SAMTOOLS mpileup -B -Q 0 -q ${MINQ:-0} -b $W/bams.txt -l $W/pos.tsv -r "$PCTG" 2> $W/mp.log > $W/pile.tsv || true
python3 - $W/sites.tsv $W/pile.tsv $OUT "$V" <<'PY'
import sys
sites, pile, out, ver = sys.argv[1:5]
pl = {}
for l in open(pile):
    f = l.rstrip("\n").split("\t")
    b = "".join(f[i] for i in range(4, len(f), 3)).upper()  # all BAMs
    pl[int(f[1])] = "".join(ch for ch in b if ch in "ACGT")
rows = []
with open(out, "w") as o:
    o.write("version\tparent_contig\tparent_pos\tchild_pos\tchild_site_class\tparent_asm_base\tchild_base\tdepth\treads_parent_asm_base\treads_child_base\n")
    for l in open(sites):
        c, p, cp, cls, pb, cb = l.rstrip("\n").split("\t"); p = int(p)
        b = pl.get(p, ""); r = (len(b), b.count(pb), b.count(cb))
        rows.append((cls,) + r)
        o.write(f"{ver}\t{c}\t{p}\t{cp}\t{cls}\t{pb}\t{cb}\t{r[0]}\t{r[1]}\t{r[2]}\n")
for cls in sorted({r[0] for r in rows}):
    sel = [r for r in rows if r[0] == cls and r[1] > 0]
    maj_p = sum(1 for r in sel if r[2] > r[3]); maj_c = sum(1 for r in sel if r[3] > r[2])
    dep = sorted(r[1] for r in sel)
    print(f"{ver}\tclass={cls}\tsites={sum(1 for r in rows if r[0] == cls)}\tcovered={len(sel)}\tmedian_depth={dep[len(dep)//2] if dep else 0}"
          f"\treads_majority_parent_asm_base={maj_p}\treads_majority_child_base={maj_c}")
PY
log "done"
