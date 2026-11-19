#!/bin/bash
# Fig5 9q windows: homolog vs donors
# usage: bash fig5_9q_compare.sh --washu D --start BP --tract A-Z --run LABEL=OUT_DIR [--run LABEL=OUT_DIR ...] --out TSV
# For each 2 kb window on SELF.p chr9 (PAN027#2#chr9) from START to the end, from each run's impg similarity class
# winners: best same-chromosome (DAD chr9) and best interchromosomal (candidate donor) identity, margin = same - inter.
# Tract = the Fig5 9q/3q tract (v1.1 136168001-136188000, v1.3.1 136175732-136195731, minimap2 asm5 20,000 of 20,000).
# Summary on stderr: per run, tract and other windows, windows the donor wins, donors that win.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
R=()
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --start) ST=$2;; --tract) TR=$2;; --run) R+=("$2");; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D ST TR OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
[ ${#R[@]} -gt 0 ] || { echo "ERROR no --run"; exit 2; }
CID=PAN027pat_vs_PAN011_joint
for r in "${R[@]}"; do
  f=$D/${r#*=}/outputs/$CID.class_winners.impg_similarity.tsv.gz
  [ -s "$f" ] || { echo "ERROR missing $f"; exit 1; }
  zcat "$f" | awk -F'\t' -v s=$ST -v r=${r%%=*} '$1 == "PAN027#2#chr9" && $2 >= s - 2000 {print r"\t"$2"\t"$13"\t"$4"\t"$12}'
done | python3 -c '
import sys, collections
t0, t1 = map(int, sys.argv[1].split("-")); w = collections.defaultdict(dict); out = open(sys.argv[2], "w")
for l in sys.stdin:
    run, start, cls, other, ident = l.rstrip("\n").split("\t")
    w[(run, int(start))][cls] = (other.split(":")[0].replace("PAN011#joint#", ""), float(ident))
out.write("run\twindow_start\tin_tract\tsame_chrom\tsame_id\tinterchrom\tinter_id\tmargin_same_minus_inter\n")
summ, don = collections.Counter(), collections.defaultdict(collections.Counter)
for (run, s), d in sorted(w.items()):
    sm, it = d.get("same_chrom", ("-", 0.0)), d.get("interchrom", ("-", 0.0)); m = sm[1] - it[1]
    k = "tract" if t0 - 1 <= s < t1 else "other"; tf = k == "tract"
    out.write(f"{run}\t{s}\t{tf}\t{sm[0]}\t{sm[1]:.6f}\t{it[0]}\t{it[1]:.6f}\t{m:+.6f}\n")
    summ[(run, k, "windows")] += 1; summ[(run, k, "donor_wins")] += m < 0
    if m < 0: don[(run, k)][it[0]] += 1
for (run, k, x), n in sorted(summ.items()): print(f"# {run}\t{k}\t{x}\t{n}", file=sys.stderr)
for (run, k), c in sorted(don.items()): print(f"# {run}\t{k}\twinning donors\t{dict(c.most_common())}", file=sys.stderr)
' "$TR" "$OUT"
echo "[$(date -Is)] wrote $OUT ($(wc -l < $OUT) lines)"
