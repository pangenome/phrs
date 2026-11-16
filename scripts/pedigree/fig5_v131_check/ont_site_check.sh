#!/bin/bash
# ONT alleles at diagnostic sites
# args: OUT_TSV, BAM, DIAG_TSV, CONTIG, TRACT_START, TRACT_END
# DIAG_TSV from diag_sites.sbatch; SNV sites only
set -euo pipefail
OUT=$1; BAM=$2; DIAG=$3; CTG=$4; TS=$5; TE=$6
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
T=$(mktemp -d); trap 'rm -rf $T' EXIT
echo "[$(date -Is)] $BAM $CTG:$TS-$TE sites from $DIAG" >&2

# Site list: position, assembly base, other copy base
awk -F'\t' -v c="$CTG" 'NR>1 { s = ($4 != ".") ? $4 : $5; if (s ~ /^SNV/) { split(s, a, " "); print c"\t"$2"\t"$3"\t"substr(a[2],1,1)"\t"substr(a[2],3,1)"\t"(($4 != ".") ? "own" : "donor") } }' "$DIAG" > $T/sites.tsv
cut -f1,2 $T/sites.tsv > $T/pos.tsv
echo "sites $(wc -l < $T/sites.tsv)" >&2

# Pileup, no reference, all MAPQ
$SAMTOOLS mpileup -B -Q 0 -q 0 -l $T/pos.tsv -r "$CTG" "$BAM" 2> $T/mp.log | cut -f2,4,5 > $T/pile.tsv

# Spanning reads: 10 kb anchors, primary
SPAN=$($SAMTOOLS view -F 0x904 "$BAM" "$CTG:$TS-$TE" | awk -v s=$TS -v e=$TE '{ l=0; c=$6; while (match(c, /[0-9]+[MDN=X]/)) { l += substr(c, RSTART, RLENGTH-1); c = substr(c, RSTART+RLENGTH) } if ($4 <= s-10000 && $4+l >= e+10000) n++ } END { print n+0 }')

python3 - $T/sites.tsv $T/pile.tsv $OUT "$CTG:$TS-$TE" $SPAN <<'PY'
import sys
sites, pile, out, tract, span = sys.argv[1:6]
pl = {int(l.split('\t')[0]): l.rstrip('\n').split('\t') for l in open(pile)}
rows = []
for l in open(sites):
    c, p, qm, other, asm, kind = l.rstrip('\n').split('\t')
    # cs *xy: x = other copy base, y = assembly base
    p = int(p); depth, n_asm, n_oth = 0, 0, 0
    if p in pl:
        b = ''.join(ch for ch in pl[p][2].upper() if ch in 'ACGT')
        depth = len(b); n_asm = b.count(asm.upper()); n_oth = b.count(other.upper())
    rows.append((p, qm, kind, depth, n_asm, n_oth))
with open(out, 'w') as o:
    o.write("pos\tquery_matches\tvs\tdepth\treads_assembly_base\treads_other_base\n")
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
cov = [r for r in rows if r[3] > 0]
fr = sorted(r[4] / r[3] for r in cov)
maj = sum(1 for r in cov if r[4] > r[5])
med = fr[len(fr)//2] if fr else float('nan')
print(f"{tract}\tsites={len(rows)}\tcovered={len(cov)}\tmedian_depth={sorted(r[3] for r in cov)[len(cov)//2] if cov else 0}\tmedian_frac_assembly_base={med:.3f}\tsites_majority_assembly={maj}\tspanning_reads_10kb={span}")
PY
