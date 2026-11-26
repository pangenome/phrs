#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11
# DAD v1.1 and v1.3.1 on the new finals (chr9 h2 end polished): merqury (work on node scratch, small outputs to merqury/<tag>),
# then g1_changes_report.py to a new prefix (report_dadf); the old reports stay for comparison
set -euxo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; T=$D/scripts_tmp
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
for v in v11 v131; do r=$([ $v = v11 ] && echo v1.1 || echo v1.3.1); F=$D/g1/final_${v}f/PAN011/fix1.pol/norm
  [ -e merqury/g1${v}f/PAN011 ] && { echo "ERROR merqury/g1${v}f/PAN011 exists"; exit 1; }
  echo "[$(date -Is)] merqury $v on $F"; bash $T/merqury_run.sbatch $D PAN011 $F g1${v}f
  mkdir -p corrected_release/$r/report_dadf
  echo "[$(date -Is)] report $v"; python3 $T/g1_changes_report.py $S $D $v PAN011 g1/rebuild_$v/PAN011/applied2 $F $T/g1_manual_fixes.tsv $v g1${v}f corrected_release/$r/report_dadf/$v.PAN011
done
echo "[$(date -Is)] DONE"
