#!/bin/bash
#SBATCH -p lowmem,386mem,highmem,1tbmem
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# Read-vs-assembly mismatch scan, chr9q ends
# args: WASHU, OUT_DIR; env PART=ends for the chr1p/chr9q MAPQ 0 and 20 scan
# out: pile.<who>.<v>.tsv (pos depth n_ref alt n_alt), asm.<pair>.<v>.tsv (target pos, query base, target base)
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p $OUT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.asmerr.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }

# Pileup parser: ref matches vs top alt
pile () { python3 -c '
import re, sys
for l in sys.stdin:
    f = l.rstrip("\n").split("\t"); b = re.sub(r"\^.|\$", "", "".join(f[4::3])); i, s = 0, []  # all BAMs
    while i < len(b):
        c = b[i]
        if c in "+-":
            m = re.match(r"[+-](\d+)", b[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        s.append(c); i += 1
    nref = sum(c in ".," for c in s); alts = {x: sum(c.upper() == x for c in s) for x in "ACGT"}
    a = max(alts, key=alts.get); print(f"{f[1]}\t{len(s)}\t{nref}\t{a}\t{alts[a]}")
'; }

# PART=ends: DAD h1 chr1p start and DAD h2 chr9q end, MAPQ>=0 and >=20
if [ "${PART:-main}" = ends ]; then
  for v in v11 v131; do
    for spec in "DAD1p:PAN011#1#chr1:1-330000" "DAD9q:PAN011#2#chr9:END"; do
      n=${spec%%:*}; r=${spec#*:}; C=${r%%:*}
      L=$(grep -P "^$C\t" norm_$v/${C%#chr*}.fa.gz.fai | cut -f2); [ "${r#*:}" = END ] && r="$C:$((L - 330000))-$L"
      $S faidx norm_$v/${C%#chr*}.fa.gz "$C" > $W/ref.fa; $S faidx $W/ref.fa; ls ont/SRR378260[57][0-9].PAN011.$v.bam > $W/bams.txt
      for q in 0 20; do
        $S mpileup -B -Q 0 -q $q -d 100000 -f $W/ref.fa -b $W/bams.txt -r "$r" 2> $W/mp.log | pile > $OUT/pile.$n.q$q.$v.tsv
        log "$n $v $r MAPQ>=$q columns $(wc -l < $OUT/pile.$n.q$q.$v.tsv)"
      done
    done
  done
  log done; exit 0
fi

for v in v11 v131; do
  for who in DAD SELF; do
    if [ $who = DAD ]; then H=PAN011; C="PAN011#2#chr9"; B="ont/SRR378260[57][0-9].PAN011.$v.bam"; else H=PAN027; C="PAN027#2#chr9"; B="ont/SRR3782606[89].PAN027.$v.bam ont/SRR3782607[0-3].PAN027.$v.bam"; fi
    L=$(grep -P "^$C\t" norm_$v/$H#2.fa.gz.fai | cut -f2)
    $S faidx norm_$v/$H#2.fa.gz "$C" > $W/ref.fa; $S faidx $W/ref.fa
    ls $B > $W/bams.txt
    log "$who $v $C length $L, BAMs $(wc -l < $W/bams.txt)"
    $S mpileup -B -Q 0 -q 20 -d 100000 -f $W/ref.fa -b $W/bams.txt -r "$C:$((L - 330000))-$L" 2> $W/mp.log | pile > $OUT/pile.$who.$v.tsv
    log "columns $(wc -l < $OUT/pile.$who.$v.tsv)"
  done
  # Assembly vs assembly on DAD h2 chr9q end
  L9=$(grep -P "^PAN011#2#chr9\t" norm_$v/PAN011#2.fa.gz.fai | cut -f2); LS=$(grep -P "^PAN027#2#chr9\t" norm_$v/PAN027#2.fa.gz.fai | cut -f2)
  $S faidx norm_$v/PAN011#2.fa.gz "PAN011#2#chr9:$((L9 - 330000))-$L9" > $W/t.fa
  $S faidx norm_$v/PAN011#1.fa.gz "PAN011#1#chr1:1-200000" > $W/chr1p.fa
  $S faidx norm_$v/PAN027#2.fa.gz "PAN027#2#chr9:$((LS - 400000))-$LS" > $W/self.fa
  for q in chr1p self; do
    $M -cx asm5 --cs=short --secondary=no -t 8 $W/t.fa $W/$q.fa 2> /dev/null > $W/$q.paf
    log "$q vs DAD h2 chr9q: $(awk '{printf "q%s-%s(%s) t%s-%s NM=%s; ", $3,$4,$5,$8,$9,$13}' $W/$q.paf)"
    python3 - $W/$q.paf $((L9 - 330000)) > $OUT/asm.$q.$v.tsv <<'PY'
import re, sys
paf, off = sys.argv[1], int(sys.argv[2]) - 1
for l in open(paf):
    f = l.rstrip("\n").split("\t"); t = int(f[7]); cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
    print(f"#block\t{off + t + 1}\t{off + int(f[8])}\t{f[4]}")
    for op, s in re.findall(r"([:*+\-])([0-9]+|[a-z]+)", cs):
        if op == ":": t += int(s)
        elif op == "*": print(f"{off + t + 1}\tSNV"); t += 1
        elif op == "-": print(f"{off + t + 1}\tDEL{len(s)}"); t += len(s)
        elif op == "+": print(f"{off + t + 1}\tINS{len(s)}")
PY
  done
done
log done
