#!/usr/bin/env python3
"""D-PeerQ3 (orientation half): orientation of within-community arm pairs.

Measures, for the CHM13 telomeric 100 kb arm ends, the strand of every alignment
between arms and reports the fraction of aligned base pairs in the same
TELOMERE-FIRST orientation.

Orientation convention: the flank sequences are stored in forward chromosome
coordinates (q-arm telomere at the high coordinate, p-arm telomere at the low
coordinate). Re-orienting both arms telomere-first therefore maps
  same-type pair (P-P or Q-Q): strand '+' = same telomere-first orientation
  mixed pair  (P-Q):            strand '-' = same telomere-first orientation

Input: arm-arm PAF produced by
  FastGA -T8 -f8192 -s200 -c70 -i0.6 -paf arms.fa arms.fa
where arms.fa holds one 100 kb CHM13 flank per arm (names = chr<NN>_<p|q>).
Usage: peerq3_orientation.py <arms.paf> [--out TSV]
"""
from __future__ import annotations
import argparse, collections, csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ASSIGN = ROOT / "data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("--out", default=str(ROOT / "scripts/cladistics/peerq3_orientation.tsv"))
    a = ap.parse_args()
    comm = {r["ChromArm"]: r["Community"] for r in csv.DictReader(open(ASSIGN), delimiter="\t")}
    typ = {x: ("P" if x.endswith("_p") else "Q") for x in comm}
    bp = collections.defaultdict(lambda: [0, 0])
    for line in open(a.paf):
        f = line.split("\t")
        q, t, strand = f[0], f[5], f[4]
        if q == t or q not in comm or t not in comm:
            continue
        k = (q, t) if q < t else (t, q)
        bp[k][0 if strand == "+" else 1] += int(f[3]) - int(f[2])

    def same_fraction(x, y, plus, minus):
        return plus / (plus + minus) if typ[x] == typ[y] else minus / (plus + minus)

    rows = []
    for i, x in enumerate(sorted(comm)):
        for y in sorted(comm):
            if not x < y or comm[x] != comm[y]:
                continue
            p, m = bp.get((x, y), [0, 0])
            if p + m < 1000:
                continue
            fr = same_fraction(x, y, p, m)
            pairtype = "P-Q" if typ[x] != typ[y] else f"{typ[x]}-{typ[y]}"
            rows.append([comm[x], x, y, pairtype, p + m, f"{fr:.4f}",
                         "same" if fr >= 0.5 else "opposite"])
    with open(a.out, "w") as fh:
        fh.write("community\tarm1\tarm2\ttype\ttotal_bp\tsame_orientation_fraction\tcall\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    mix = [r for r in rows if r[3] == "P-Q"]
    st = [r for r in rows if r[3] != "P-Q"]
    print(f"wrote {a.out}")
    print(f"P-Q within-community pairs alignable (>=1 kb): {len(mix)}/17 "
          f"same {sum(1 for r in mix if r[6]=='same')} "
          f"opposite {sum(1 for r in mix if r[6]=='opposite')}")
    print(f"same-type pairs alignable: {len(st)} "
          f"same {sum(1 for r in st if r[6]=='same')} "
          f"opposite {sum(1 for r in st if r[6]=='opposite')}")


if __name__ == "__main__":
    main()
