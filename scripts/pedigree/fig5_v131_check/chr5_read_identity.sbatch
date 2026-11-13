#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 8:00:00
# Per-read divergence, KID chr5 tract, both assemblies
# args: WASHU, OUT_DIR ; v1.1
# Reads spanning PAN028#1#chr5 183,044,000-183,074,000 +-5 kb at MAPQ>=20 (primary): SELF reads (PAN027.onPAN028)
# and KID reads (PAN028 own). For each read, every alignment record in: KID assembly BAMs and SELF assembly BAMs
# (PAN027 own for SELF reads, PAN028.onPAN027 for KID reads): contig, start, end, MAPQ, flag, de:f, clips.
# out: OUT_DIR/read_identity.tsv (set, read, assembly, contig, start, end, mapq, flag, de, clip_l, clip_r)
set -euo pipefail
D=$1; O=$2; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.chr5id.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
T0=183044000; T1=183074000
span () { for b in "$@"; do $S view -q 20 -F 0x904 $b PAN028#1#chr5:$T0-$T1; done | awk -v a=$((T0-5000)) -v z=$((T1+5000)) '{e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; if ($4<=a && e>=z) print $1}' | sort -u; }
span ont_full/SRR*.PAN027.onPAN028.v11.bam > $W/self.txt; span ont_full/SRR*.PAN028.v11.bam > $W/kid.txt
log "spanning reads: SELF $(wc -l < $W/self.txt), KID $(wc -l < $W/kid.txt)"
rec () { # set assembly names bams...
  local set=$1 asm=$2 names=$3; shift 3
  for b in "$@"; do $S view -N $names $b; done | awk -F'\t' -v s=$set -v a=$asm -v OFS='\t' '{e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)};
    de="NA"; for (i=12;i<=NF;i++) if ($i ~ /^de:f:/) de=substr($i,6); l=0; r=0; if ($6 ~ /^[0-9]+[SH]/) l=$6+0; if (match($6,/[0-9]+[SH]$/)) r=substr($6,RSTART)+0;
    print s, $1, a, $3, $4, e, $5, $2, de, l, r}'
}
printf "set\tread\tassembly\tcontig\tstart\tend\tmapq\tflag\tde\tclip_l\tclip_r\n" > $O/read_identity.tsv
rec SELF KID $W/self.txt ont_full/SRR*.PAN027.onPAN028.v11.bam >> $O/read_identity.tsv
rec SELF SELF $W/self.txt ont_full/SRR*.PAN027.v11.bam >> $O/read_identity.tsv
rec KID KID $W/kid.txt ont_full/SRR*.PAN028.v11.bam >> $O/read_identity.tsv
rec KID SELF $W/kid.txt ont_full/SRR*.PAN028.onPAN027.v11.bam >> $O/read_identity.tsv
log "records: $(($(wc -l < $O/read_identity.tsv) - 1))"
log "primary records (flag without 256, 2048), per set and assembly: contig, n, median de:"
awk -F'\t' 'NR>1 && int($8/256)%2==0 && int($8/2048)%2==0 {k=$1"\t"$3"\t"$4; n[k]++; d[k]=d[k]" "$9} END {for (k in n) print k"\t"n[k]"\t"d[k]}' $O/read_identity.tsv \
  | python3 -c "import sys,statistics as st; [print(l.split('\t')[0], l.split('\t')[1], l.split('\t')[2], l.split('\t')[3], 'median de', st.median(float(x) for x in l.split('\t')[4].split() if x!='NA')) for l in sys.stdin]"
log done
