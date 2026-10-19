#!/usr/bin/env python3
"""Orientation of shared modules between chromosome ends, relative to the telomere.

Question (reviewer Q3 follow-up): the blocks/modules shared between subtelomeric arm
ends do not match end-to-end and their ORDER along the arm differs. Is their
ORIENTATION relative to the telomere also scrambled, or consistently the same?

Data: HPRCv2 all-vs-all alignments of the telomeric 500 kb flanks
(/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95, one PAF per CHM13 flank). Each
file holds every sampled flank aligned against that flank, so the 46 CHM13 files cover
all arm ends at once.

Orientation rule. Flanks are stored in forward chromosome coordinates, so the telomere
lies at the LOW coordinate for a p-arm and the HIGH coordinate for a q-arm. Re-orienting
both arms telomere-first, an alignment is SAME orientation when
    strand '+' for P-P and Q-Q pairs, strand '-' for P-Q pairs,
and OPPOSITE otherwise. This is per-block: no colinearity or offset matching is required,
because the question is about block orientation, not arm-end colinearity.

Homolog-vs-homolog pairs (same arm, different individual) are excluded - they are
trivially same-orientation and would inflate the statistic.

Usage: python3 scripts/cladistics/arm_module_orientation.py [--min-block 5000]
Writes scripts/cladistics/arm_module_orientation.tsv
"""
from __future__ import annotations
import argparse, collections, glob, gzip, re, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAF_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95")
OUT = ROOT / "scripts/cladistics/arm_module_orientation.tsv"
NAMERE = re.compile(r"_(chr[0-9XY]+_[pq])arm")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-block", type=int, default=5000)
    a = ap.parse_args()
    bp = collections.defaultdict(int)
    files = sorted(glob.glob(str(PAF_DIR / "hprcv2.1Mb.telo_500kb_trimmed.all-vs-CHM13_0_chr*.paf.gz")))
    if not files:
        raise SystemExit(f"no PAFs found under {PAF_DIR}")
    for f in files:
        m = re.search(r"_(chr[0-9XY]+_[pq])arm\.paf\.gz$", f)
        if not m:
            continue
        t = m.group(1)
        with gzip.open(f, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                qm = NAMERE.search(c[0])
                if not qm:
                    continue
                q = qm.group(1)
                if q == t:
                    continue                      # homolog vs homolog: excluded
                L = int(c[3]) - int(c[2])
                if L < a.min_block:
                    continue
                bp[(t, q, c[4])] += L
    typ = lambda x: "P" if x.endswith("_p") else "Q"
    rows = []
    for (t, q) in {(k[0], k[1]) for k in bp}:
        plus, minus = bp[(t, q, "+")], bp[(t, q, "-")]
        tot = plus + minus
        if tot < 20000:
            continue
        same = plus if typ(t) == typ(q) else minus
        rows.append([t, q, f"{typ(t)}-{typ(q)}", same, tot - same, f"{same/tot:.4f}", tot])
    rows.sort(key=lambda r: float(r[5]))
    with open(OUT, "w") as fh:
        fh.write("target_arm\tquery_arm\ttype\tsame_bp\topposite_bp\tsame_fraction\ttotal_bp\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    S = sum(r[3] for r in rows); T = sum(r[6] for r in rows)
    fr = [float(r[5]) for r in rows]
    print(f"arm pairs (>=20 kb aligned, blocks >= {a.min_block} bp): {len(rows)}")
    print(f"same telomere-relative orientation: {100*S/T:.1f}% of {T/1e9:.1f} Gbp")
    print(f"per-pair same-fraction: median {100*st.median(fr):.1f}%, "
          f"min {100*min(fr):.1f}%, pairs <90%: {sum(1 for x in fr if x<0.9)}")
    for tt in ("P-P", "Q-Q", "P-Q"):
        d = [r for r in rows if r[2] in (tt, "Q-P" if tt == "P-Q" else tt)]
        if d:
            s = sum(r[3] for r in d); t_ = sum(r[6] for r in d)
            print(f"  {tt}: {len(d):3d} pairs, same {100*s/t_:.1f}%")
    print("most inverted pairs:")
    for r in rows[:6]:
        print(f"  {r[0]:8s} x {r[1]:8s} ({r[2]})  same {100*float(r[5]):5.1f}%  {r[6]/1e3:8.0f} kb")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
