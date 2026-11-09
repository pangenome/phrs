#!/bin/bash
# Fig5 class winners: original vs patched DAD
# args: WASHU, [START] (default 136139447, patch frame start on SELF.p chr9, v1.1)
# For each 2 kb window on PAN027#2#chr9 from START to the end, prints the same-chromosome and the best
# interchromosomal identity in Run B (out_B_v11) and in the patched run (patch_v11/out_P_v11).
set -euo pipefail
D=$1; ST=${2:-136139447}; CID=PAN027pat_vs_PAN011_joint
T0=136168001; T1=136188000  # Fig5 tract, SELF.p v1.1
for run in out_B_v11 patch_v11/out_P_v11; do
  f=$D/$run/outputs/$CID.class_winners.impg_similarity.tsv.gz
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  zcat "$f" | awk -F'\t' -v s=$ST -v r=${run%%/*} '$1 == "PAN027#2#chr9" && $2 >= s - 2000 {print r"\t"$2"\t"$13"\t"$4"\t"$12}'
done | python3 -c '
import sys, collections
t0, t1 = int(sys.argv[1]), int(sys.argv[2]); w = collections.defaultdict(dict)
for l in sys.stdin:
    run, start, cls, other, ident = l.rstrip("\n").split("\t")
    w[(run, int(start))][cls] = (other.split(":")[0].replace("PAN011#joint#", ""), float(ident))
print("run\twindow_start\tsame_chrom\tsame_id\tinterchrom\tinter_id\tmargin_same_minus_inter")
summ = collections.Counter()
for (run, s), d in sorted(w.items()):
    sm, it = d.get("same_chrom", ("-", 0.0)), d.get("interchrom", ("-", 0.0))
    m = sm[1] - it[1]; print(f"{run}\t{s}\t{sm[0]}\t{sm[1]:.6f}\t{it[0]}\t{it[1]:.6f}\t{m:+.6f}")
    k = "tract" if t0 - 1 <= s < t1 else "other"; summ[(run, k, "windows")] += 1; summ[(run, k, "inter_wins")] += m < 0
for (run, k, x), n in sorted(summ.items()): print(f"# {run}\t{k}\t{x}\t{n}", file=sys.stderr)
' $T0 $T1
