#!/bin/bash
# Per finished 5 Mb target: name, CPU (user+sys summed over parts), number of parts, max part RSS (GB).
set -euo pipefail
cd /moosefs/guarracino/HPRCv2/PHR_III/flank5mb/allvsall5mb
for f in *.time; do
  [ -s ${f%.time}.paf.gz ] || continue
  awk -F': ' -v n=${f%.time} '/User time|System time/ {s += $2} /Maximum resident/ {k++; if ($2 > m) m = $2} END {printf "%s\t%.0f\t%d\t%.1f\n", n, s, k, m / 1048576}' $f
done
