#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11
# Fig. 5d candidate competition on the final corrected v1.1 assemblies (g1_finals.tsv), every candidate end offered at once.
# usage: sbatch fig5d_competition.sbatch OUT_DIR   (run from the washu_v131 dir; OUT_DIR new, on moosefs; work on node scratch)
# Candidates: the terminal 200 kb of each end (q arms as assembled; p arms reverse-complemented so that the telomere is last):
#   SELFp_chr9q (PAN027#2#chr9, the daughter's paternal chr9q), DAD_chr9h2_CORR (final), DAD_chr9h2_ORIG (released norm_v11),
#   DAD_chr9h1, DAD_chr3h2, DAD_chr19h2 (q), DAD_chr1h1, DAD_chr2p_CONTROL (h1, p): the names of the 2026-11-18 run
#   (paper_prep/figures/fig5_concerted/candidate_similarity_9q.tsv), whose candidate FASTA is not on disk.
# Steps: sweepga all-vs-all (many:many, no scaffold jumps, FastGA frequency 5000 so that no candidate is filtered out);
#   impg similarity --gfa-engine poa on 2 kb windows of SELFp_chr9q (as fig5_similarity_s.sbatch); per window and candidate the
#   best dice.similarity and its estimated.identity (block identity in poa mode); means inside and outside the tract.
# Tract: PAN027#2#chr9:136,168,000-136,188,000 of the released v1.1 (response/notes/M2_pedigree_completions.md), taken at the same
#   distance from the chromosome end (the last 200 kb of SELF.p chr9 are identical in the released and corrected v1.1).
# out: OUT_DIR/candidates.fa.fai, cand.paf.gz, similarity.tsv.gz (raw impg), per_window.tsv, candidate_graph_dice_9q.tsv
set -euxo pipefail
OUT=$(realpath -m $1); [ -e $OUT ] && { echo "ERROR $OUT exists"; exit 1; }
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
SWEEPGA=/home/erikg/.cargo/bin/sweepga; IMPG=/home/erikg/.cargo/bin/impg; export PATH=/home/erikg/.cargo/bin:$PATH
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/f5d.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log () { echo "[$(date -Is)] $*"; }
fin () { awk -F'\t' -v s=$1 '$1 == "v11" && $2 == s {print $3}' $T/g1_finals.tsv; }
SELF=$(fin PAN027); DAD=$(fin PAN011); log "host $(hostname); SELF $SELF; DAD $DAD; work $W; $($SWEEPGA --version); $($IMPG --version)"
E=200000
end () {  # FASTA CONTIG NAME p|q -> the 200 kb end as a FASTA record named NAME
  local L; L=$(awk -v c="$2" '$1 == c {print $2}' "$1.fai")
  if [ $4 = q ]; then $S faidx "$1" "$2:$((L - E + 1))-$L"; else $S faidx -i "$1" "$2:1-$E"; fi | sed "1s/.*/>$3/"
  log "$3: $2 ($4 end) of $1, contig length $L" >&2
}
{ end $SELF/PAN027#2.fa.gz 'PAN027#2#chr9' SELFp_chr9q q
  end $DAD/PAN011#2.fa.gz 'PAN011#2#chr9' DAD_chr9h2_CORR q
  end norm_v11/PAN011#2.fa.gz 'PAN011#2#chr9' DAD_chr9h2_ORIG q
  end $DAD/PAN011#1.fa.gz 'PAN011#1#chr9' DAD_chr9h1 q
  end $DAD/PAN011#2.fa.gz 'PAN011#2#chr3' DAD_chr3h2 q
  end $DAD/PAN011#2.fa.gz 'PAN011#2#chr19' DAD_chr19h2 q
  end $DAD/PAN011#1.fa.gz 'PAN011#1#chr1' DAD_chr1h1 p
  end $DAD/PAN011#1.fa.gz 'PAN011#1#chr2' DAD_chr2p_CONTROL p; } > $W/candidates.fa
$S faidx $W/candidates.fa; cat $W/candidates.fa.fai
[ $(wc -l < $W/candidates.fa.fai) -eq 8 ] && [ $(awk -v e=$E '$2 != e' $W/candidates.fa.fai | wc -l) -eq 0 ] || { echo "ERROR not 8 candidates of $E bp"; exit 1; }
LS=$(awk '$1 == "PAN027#2#chr9" {print $2}' norm_v11/PAN027#2.fa.gz.fai); TA=$((136168000 - (LS - E))); TB=$((136188000 - (LS - E)))
[ $TA -ge 0 ] && [ $TB -le $E ] || { echo "ERROR tract $TA-$TB outside the segment"; exit 1; }
log "released SELF.p chr9 length $LS; tract on SELFp_chr9q: $TA-$TB"
cd $W
log "sweepga all-vs-all"; $SWEEPGA --fastga --fastga-frequency 5000 --num-mappings many:many --scaffold-jump 0 --temp-dir $W --output-file cand.paf candidates.fa
log "PAF lines $(wc -l < cand.paf)"; $BGZ -c cand.paf > cand.paf.gz
awk -v e=$E 'BEGIN {OFS = "\t"; for (a = 0; a < e; a += 2000) print "SELFp_chr9q", a, a + 2000}' > win.bed
log "impg similarity poa, windows $(wc -l < win.bed)"
$IMPG similarity --alignment-files cand.paf.gz --target-bed win.bed --sequence-files candidates.fa --gfa-engine poa --no-merge \
  --num-mappings many:many --scaffold-jump 0 --threads $SLURM_CPUS_PER_TASK > sim.tsv
log "similarity lines $(wc -l < sim.tsv)"; head -1 sim.tsv
python3 - sim.tsv $TA $TB <<'PY'
import csv, sys, collections
f, ta, tb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
best = {}  # (window start, candidate) -> (dice, identity): the best SELFp_chr9q-vs-candidate pair in the window
for r in csv.DictReader(open(f), delimiter="\t"):
    a, b = r["group.a"].split(":")[0], r["group.b"].split(":")[0]
    if a != "SELFp_chr9q" or b == "SELFp_chr9q": continue
    k = (int(r["start"]), b); v = (float(r["dice.similarity"]), float(r["estimated.identity"]))
    if k not in best or v > best[k]: best[k] = v
with open("per_window.tsv", "w") as o:
    o.write("window_start\tcandidate\tin_tract\tdice\tblock_identity\n")
    for (w, c), (d, i) in sorted(best.items()): o.write(f"{w}\t{c}\t{int(ta <= w < tb)}\t{d}\t{i}\n")
nwin = collections.Counter((c, ta <= w < tb) for (w, c) in best)
with open("candidate_graph_dice_9q.tsv", "w") as o:
    o.write("candidate\ttract_dice\ttract_block_identity\ttract_windows\toutside_dice\toutside_block_identity\toutside_windows\n")
    for c in sorted({c for (_, c) in best}):
        row = [c]
        for t in (True, False):
            v = [best[(w, cc)] for (w, cc) in best if cc == c and (ta <= w < tb) == t]
            row += [f"{sum(x[0] for x in v) / len(v):.4f}", f"{sum(x[1] for x in v) / len(v):.4f}", str(len(v))] if v else ["NA", "NA", "0"]
        o.write("\t".join(row) + "\n")
print(open("candidate_graph_dice_9q.tsv").read(), file=sys.stderr)
PY
mkdir -p $OUT; cp candidates.fa.fai cand.paf.gz per_window.tsv candidate_graph_dice_9q.tsv $OUT/; gzip -c sim.tsv > $OUT/similarity.tsv.gz
log "DONE $OUT"
