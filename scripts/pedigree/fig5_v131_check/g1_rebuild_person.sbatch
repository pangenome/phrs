#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 24
#SBATCH --mem=128G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Rebuild one person's final corrected assembly from the original release, in the order used (called by g1_rebuild.sh)
# usage: sbatch g1_rebuild_person.sbatch --washu D --version V --sample SM --runs RUN_DIR --out DIR
# Steps (each script runs inside this job with bash):
#  1 candidates (g1_candidates_hifi.py on the stored HiFi/ONT scan bins and merqury/V) -> task list; must equal g1/apply_V/SM.tasks.tsv
#  1b gap rule: loci whose task window has >= 100 kb of N on norm_V are not applied (large gaps: rDNA, satellites); small
#     gaps may be filled when the reads vote for the new piece (v1.1: 3 filled, 1 removed; v1.3.1: 5 loci skipped) -> gap100k.keep.txt
#  2 G1 round 1: g1_apply.sbatch --keep gap100k.keep.txt (pieces of RUN_DIR = g1_locus_array output, anchor ex:5000, vote more) -> applied
#  3 g1_polish_all.sbatch -> polished;  4 merqury_run.sbatch (tag <name of DIR's parent>.r1);  5 g1_merqury_filter.py (V vs that tag) -> keep list
#  6 G1 round 2: g1_apply.sbatch --keep -> applied2;  7 g1_polish_all.sbatch -> polished2
#  8 g1_manual_fixes.sh: manual fixes of g1_manual_fixes.tsv for V and SM, in file order (end_replace: g1_end_replace.py, then
#    hifi_polish.sbatch on the inserted piece when g1_end_polish.tsv lists it; tip_trim: tip_trim.py)
#  9 per-contig md5 vs the final of g1_finals.tsv (seq_md5.py); a difference exits non-zero
# out: DIR/{cand.tsv,tasks.tsv,applied,polished,merqury_filter.*,applied2,polished2,fix<N>}, DIR/final_dir.txt, DIR/compare.tsv
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --version) V=$2;; --sample) SM=$2;; --runs) R=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D V SM R O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; T=$D/scripts_tmp; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
[ -e $O ] && { echo "ERROR $O exists; use a new dir"; exit 1; }
MT=$(basename $(dirname $O)).r1; [ -e merqury/$MT/$SM ] && { echo "ERROR merqury/$MT/$SM exists"; exit 1; }
mkdir -p $O; log "rebuild $SM $V: runs $R, out $O, host $(hostname), cpus $SLURM_CPUS_PER_TASK"
log "step 1 candidates"
python3 $T/g1_candidates_hifi.py $D $SM $V $O/cand.tsv
awk -F'\t' -v OFS='\t' 'NR > 1 && $11 == "candidate" {n++; r = $1 ":" $2 + 1 "-" $3; print n, "g" n, r, r}' $O/cand.tsv > $O/tasks.tsv
cmp $O/tasks.tsv g1/apply_$V/$SM.tasks.tsv || { echo "ERROR task list differs from g1/apply_$V/$SM.tasks.tsv"; exit 1; }
log "tasks: $(wc -l < $O/tasks.tsv), same as g1/apply_$V/$SM.tasks.tsv; loci with summary in $R: $(ls $R/*/summary.tsv | wc -l)"
log "step 1b gap rule"
while IFS=$'\t' read -r i n r x; do c=$($S faidx norm_$V/${r%%#chr*}.fa.gz "$r" | grep -v '>' | tr -cd 'Nn' | wc -c); if [ $c -lt 100000 ]; then echo "$i.$n"; else log "skip $i.$n $r: $c N" >&2; fi; done < $O/tasks.tsv > $O/gap100k.keep.txt
log "loci kept by the gap rule: $(wc -l < $O/gap100k.keep.txt) of $(wc -l < $O/tasks.tsv)"
log "step 2 G1 round 1 apply"
bash $T/g1_apply.sbatch --washu $D --sample $SM --version $V --runs $R --tasks $O/tasks.tsv --anchor ex:5000 --vote more --out $O/applied --keep $O/gap100k.keep.txt
log "step 3 polish"; bash $T/g1_polish_all.sbatch --washu $D --sample $SM --version $V --in-dir $O/applied --out $O/polished
log "step 4 merqury $MT"; bash $T/merqury_run.sbatch $D $SM $O/polished/norm $MT
log "step 5 merqury filter"; python3 $T/g1_merqury_filter.py $D $SM $V $MT $O/applied/corrections.tsv $O/merqury_filter
log "keep $(wc -l < $O/merqury_filter.keep.txt) of $(awk -F'\t' '$1 == "applied"' $O/applied/corrections.tsv | wc -l) applied pieces"
log "step 6 G1 round 2 apply"
bash $T/g1_apply.sbatch --washu $D --sample $SM --version $V --runs $R --tasks $O/tasks.tsv --anchor ex:5000 --vote more --out $O/applied2 --keep $O/merqury_filter.keep.txt
log "step 7 polish"; bash $T/g1_polish_all.sbatch --washu $D --sample $SM --version $V --in-dir $O/applied2 --out $O/polished2
log "step 8 manual fixes"; F=$(bash $T/g1_manual_fixes.sh $D $V $SM $O/polished2/norm $O/fix | tail -1)
log "final $F"; echo $F > $O/final_dir.txt
log "step 9 compare with the final in g1_finals.tsv"
REF=$(awk -F'\t' -v v=$V -v s=$SM '$1 == v && $2 == s {print $3}' $T/g1_finals.tsv); [ -n "$REF" ] || { echo "ERROR no final for $V $SM"; exit 1; }
python3 $T/seq_md5.py $F $SM | awk -F'\t' -v OFS='\t' '{print $1 "|" $2, $3, $4}' | sort > $O/new.md5
python3 $T/seq_md5.py $REF $SM | awk -F'\t' -v OFS='\t' '{print $1 "|" $2, $3, $4}' | sort > $O/ref.md5
join -t $'\t' -a 1 -a 2 -e NA -o 0,1.2,1.3,2.2,2.3 $O/new.md5 $O/ref.md5 | awk -F'\t' -v OFS='\t' 'BEGIN {print "contig", "new_length", "new_md5", "final_length", "final_md5", "same"} {print $0, ($3 == $5 ? "yes" : "no")}' > $O/compare.tsv
nd=$(awk -F'\t' 'NR > 1 && $6 == "no"' $O/compare.tsv | wc -l)
log "contigs $(($(wc -l < $O/compare.tsv) - 1)), different from $REF: $nd"
[ $nd -eq 0 ] || { awk -F'\t' '$6 == "no"' $O/compare.tsv; echo "ERROR rebuilt assembly differs from $REF"; exit 1; }
log "DONE $SM $V identical to $REF"
