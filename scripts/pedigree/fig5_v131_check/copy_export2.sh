#!/bin/bash
# Save my /export2/scratch data to moosefs
# run on octopus01 (/export2 is local there): nohup bash copy_export2.sh > LOG 2>&1 &
# Copies (rsync -a, no delete, 80 MB/s cap for moosefs) hprc25272.aln.paf.gz, hprcv2-all-vs-all.1Mb.p95, rev2,
# guarracino, chry_rebuild, .chrY143.bps into a new dir, then checks: md5 of the PAF on both sides, rsync --checksum
# dry run per directory (must list no file), file counts and bytes on both sides.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
S=/export2/scratch; O=/moosefs/guarracino/export2_scratch_backup_20261118
log () { echo "[$(date -Is)] $*"; }
[ "$(hostname)" = octopus01 ] || { echo "ERROR run on octopus01, not $(hostname)"; exit 1; }
mkdir -p $O
for x in hprc25272.aln.paf.gz hprcv2-all-vs-all.1Mb.p95 rev2 guarracino chry_rebuild .chrY143.bps; do
  [ -e $S/$x ] || { echo "ERROR missing $S/$x"; exit 1; }
  log "copy $x ($(du -s --apparent-size --block-size=1 $S/$x | cut -f1) bytes)"
  rsync -a --partial --bwlimit=80M $S/$x $O/
  log "copied $x"
done
log "check: md5 of hprc25272.aln.paf.gz"
A=$(md5sum < $S/hprc25272.aln.paf.gz | cut -c1-32); B=$(md5sum < $O/hprc25272.aln.paf.gz | cut -c1-32)
[ "$A" = "$B" ] || { echo "ERROR md5 differs: $A vs $B"; exit 1; }
log "md5 equal $A"
for x in hprcv2-all-vs-all.1Mb.p95 rev2 guarracino chry_rebuild .chrY143.bps; do
  n=$(rsync -a -n --checksum --out-format='%n' $S/$x $O/ | grep -v '/$' | wc -l)
  [ $n -eq 0 ] || { echo "ERROR $x: $n files differ after copy"; exit 1; }
  log "$x: 0 files differ; source $(find $S/$x -type f | wc -l) files $(du -s --apparent-size --block-size=1 $S/$x | cut -f1) bytes, copy $(find $O/$x -type f | wc -l) files $(du -s --apparent-size --block-size=1 $O/$x | cut -f1) bytes"
done
log "all copied and checked into $O"
