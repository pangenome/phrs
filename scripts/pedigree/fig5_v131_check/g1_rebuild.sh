#!/bin/bash
# Rebuild the 4 final corrected pedigree assemblies (MOM PAN010, DAD PAN011, SELF PAN027, KID PAN028) from the original
# release, in the order used: G1 (2 rounds, merqury filter between them), then the manual fixes of g1_manual_fixes.tsv
# (DAD chr9q end, SELF.m chr22q tail), then a per-contig md5 comparison with the finals of g1_finals.tsv.
# run on octopus01: bash g1_rebuild.sh --version V --out DIR [--loci]
#   default: G1 uses the stored local assemblies g1/apply_V/SM/runs (g1_locus_array output);
#   --loci: reruns them first (g1_locus_array.sbatch, 4 arrays, about 1 day), then uses DIR/SM/runs.
# Needs in D: norm_V, scan bins (check/wg_scan), merqury/V, reads, g1/mmi/V, scripts in scripts_tmp (deployed from the repo).
# One job per person (g1_rebuild_person.sbatch); prints the job IDs.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
LOCI=0
while [ $# -gt 0 ]; do case $1 in --version) V=$2; shift 2;; --out) O=$2; shift 2;; --loci) LOCI=1; shift;; *) echo "ERROR unknown argument $1"; exit 2;; esac; done
for v in V O; do [ -n "${!v:-}" ] || { echo "ERROR missing --${v,,}"; exit 2; }; done
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; T=$D/scripts_tmp; cd $D
[ -e $O ] && { echo "ERROR $O exists; use a new dir"; exit 1; }
mkdir -p $O; O=$(realpath $O)
for SM in PAN010 PAN011 PAN027 PAN028; do
  R=g1/apply_$V/$SM/runs; dep=""
  if [ $LOCI -eq 1 ]; then
    n=$(( ($(wc -l < g1/apply_$V/$SM.tasks.tsv) + 9) / 10 )); R=$O/$SM/runs
    a=$(sbatch --parsable --array=1-$n%3 -J rb_loci_$SM -o $D/logs/rb_loci.$SM.%A_%a.log $T/g1_locus_array.sbatch --washu $D --sample $SM --version $V --tasks g1/apply_$V/$SM.tasks.tsv --flank 400000 --out $R --mmi-dir g1/mmi/$V --chunk 10 --rules ex:5000,ex:20000 --time-limit 7200)
    dep="-d afterok:$a"; echo "$SM loci $a"
  fi
  j=$(sbatch --parsable $dep -J rb_$SM -o $D/logs/rb_person.$SM.$V.%j.log $T/g1_rebuild_person.sbatch --washu $D --version $V --sample $SM --runs $R --out $O/$SM)
  echo "$SM rebuild $j"
done
