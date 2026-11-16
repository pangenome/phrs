#!/bin/bash
#SBATCH -p workers
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 2:00:00
# Diagnostic sites: own copy vs donor copy
# args: WASHU, OUT tsv, VERSION, QREGION, OWN_REF, DONOR_REF
# refs as HAP#CHR:start-end, e.g. PAN011#2#chr5:180900000-181932123
set -euo pipefail
WASHU=$1; OUT=$2; V=$3; QREG=$4; OWN=$5; DON=$6
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.diag.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
fa () { local r=$1; local h=${r%%#chr*}; $SAMTOOLS faidx "$WASHU/norm_$V/$h.fa.gz" "$r"; }

fa "$QREG" > $W/q.fa; fa "$OWN" > $W/own.fa; fa "$DON" > $W/don.fa
log "V=$V query $QREG $(grep -v '>' $W/q.fa | tr -d '\n' | wc -c) bp; own $OWN; donor $DON"

# Query vs own and donor, cs tags
for t in own don; do
  $MM2 -x asm20 -c --cs=short --secondary=no -t 8 $W/$t.fa $W/q.fa > $W/$t.paf 2> $W/$t.log
  log "$t alignments $(wc -l < $W/$t.paf): $(awk '{printf "%s-%s(%s) id=%.5f; ", $3,$4,$5,$10/$11}' $W/$t.paf)"
done

# Classify sites by query allele
python3 - $W/own.paf $W/don.paf $OUT "$V" "$QREG" <<'PY'
import sys, re
own_paf, don_paf, out, ver, qreg = sys.argv[1:6]
qoff = int(qreg.split(':')[1].split('-')[0]) - 1
COMP = str.maketrans('acgtn', 'tgcan')  # minus-strand cs bases

def diffs(paf):
    """query positions covered and differing"""
    cov, dif = set(), {}
    for l in open(paf):
        f = l.rstrip('\n').split('\t')
        cs = [x for x in f[12:] if x.startswith('cs:Z:')][0][5:]
        rev = f[4] == '-'
        q = int(f[3]) if rev else int(f[2])
        for op, s in re.findall(r'([:*+\-])([0-9]+|[a-z]+)', cs):
            if op == ':':
                n = int(s); r = range(q - n, q) if rev else range(q, q + n)
                cov.update(r); q = q - n if rev else q + n
            elif op == '*':
                p = q - 1 if rev else q
                x, y = (s[0].translate(COMP), s[1].translate(COMP)) if rev else (s[0], s[1])
                cov.add(p); dif[p] = f"SNV {x}>{y}"; q = p if rev else q + 1
            elif op == '+':
                n = len(s); r = range(q - n, q) if rev else range(q, q + n)
                for p in r: cov.add(p); dif[p] = 'INS'
                q = q - n if rev else q + n
            elif op == '-':
                dif[q - 1 if rev else q] = f"DEL{len(s)}"
    return cov, dif

oc, od = diffs(own_paf); dc, dd = diffs(don_paf)
both = oc & dc
own_only = [p for p in od if p in both and p not in dd]   # query = donor, != own
don_only = [p for p in dd if p in both and p not in od]   # query = own, != donor
neither = [p for p in od if p in both and p in dd]
with open(out, 'w') as o:
    o.write("version\tquery_pos\tquery_matches\tvs_own\tvs_donor\n")
    for p in sorted(own_only): o.write(f"{ver}\t{p+qoff+1}\tdonor\t{od[p]}\t.\n")
    for p in sorted(don_only): o.write(f"{ver}\t{p+qoff+1}\town\t.\t{dd[p]}\n")
    for p in sorted(neither):  o.write(f"{ver}\t{p+qoff+1}\tneither\t{od[p]}\t{dd[p]}\n")
print(f"{ver}\tcovered_by_both={len(both)}\tquery_like_own={len(don_only)}\tquery_like_donor={len(own_only)}\tneither={len(neither)}")
PY
log "done"
