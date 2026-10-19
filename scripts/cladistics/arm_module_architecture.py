#!/usr/bin/env python3
"""Full architectural pattern of subtelomeric modules.

Three measurements, all from the HPRCv2 all-vs-all flank alignments
(/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95, one PAF per CHM13 500 kb arm end):

1. TILING - every arm end cut into 25 kb bins; each bin labelled with the partner arm that
   contributes the most aligned base pairs to it, aggregated across all sampled individuals.
   This is the "module map" of an arm end.
2. BLOCK SIZE - length distribution of the individual shared blocks, compared for partner
   arms in the same Leiden community versus partners in other communities.
3. INVERSION STRUCTURE - for the arm pairs that are orientation-inverted, whether the
   inversion is one large block (a single structural rearrangement) or many small blocks
   (repeat-family architecture).

Writes scripts/cladistics/arm_module_tiling.tsv, arm_module_blocks.tsv,
arm_module_inversions.tsv.
"""
from __future__ import annotations
import collections, csv, glob, gzip, re, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAF_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95")
ASSIGN = ROOT / "data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
NAMERE = re.compile(r"_(chr[0-9XY]+_[pq])arm")
BIN = 25000


def main():
    comm = {r["ChromArm"]: r["Community"] for r in csv.DictReader(open(ASSIGN), delimiter="\t")}
    bins = collections.defaultdict(lambda: collections.Counter())      # (P, bin) -> partner: bp
    blocks = collections.defaultdict(list)                              # (P, X) -> [lengths]
    inv_detail = collections.defaultdict(list)                          # (P, X) -> [(off, len, strand)]
    files = sorted(glob.glob(str(PAF_DIR / "hprcv2.1Mb.telo_500kb_trimmed.all-vs-CHM13_0_chr*.paf.gz")))
    for f in files:
        m = re.search(r"_(chr[0-9XY]+_[pq])arm\.paf\.gz$", f)
        if not m:
            continue
        P = m.group(1)
        plen = 500000
        with gzip.open(f, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                qm = NAMERE.search(c[0])
                if not qm:
                    continue
                X = qm.group(1)
                ts, te = int(c[7]), int(c[8])
                L = int(c[3]) - int(c[2])
                if L <= 0:
                    continue
                blocks[(P, X)].append(L)
                inv_detail[(P, X)].append((ts / plen, L, c[4]))
                b0, b1 = ts // BIN, (te - 1) // BIN
                share = L / (b1 - b0 + 1)
                for b in range(b0, b1 + 1):
                    bins[(P, b)][X] += share
    # 1. tiling
    with open(ROOT / "scripts/cladistics/arm_module_tiling.tsv", "w") as fh:
        fh.write("arm\tbin\tstart\tend\tpartner_arm\tpartner_community\tbp\n")
        for (P, b), cnt in sorted(bins.items()):
            partner, bp = cnt.most_common(1)[0]
            if partner == P:
                continue
            fh.write(f"{P}\t{b}\t{b*BIN}\t{(b+1)*BIN}\t{partner}\t{comm.get(partner,'NA')}\t{bp:.0f}\n")
    # 2. block sizes: same-community vs different-community partners
    rows = []
    for (P, X), ls in blocks.items():
        if P == X or P not in comm or X not in comm:
            continue
        tot = sum(ls)
        if tot < 20000:
            continue
        rows.append([P, X, comm[P], comm[X], "same" if comm[P] == comm[X] else "diff",
                     len(ls), tot, f"{st.median(ls):.0f}", max(ls)])
    with open(ROOT / "scripts/cladistics/arm_module_blocks.tsv", "w") as fh:
        fh.write("arm\tpartner\tcomm_arm\tcomm_partner\trelation\tn_blocks\ttotal_bp\tmedian_block\tmax_block\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    same = [r for r in rows if r[4] == "same"]
    diff = [r for r in rows if r[4] == "diff"]
    print(f"partner pairs (>=20 kb): {len(rows)}  same-community {len(same)}  different {len(diff)}")
    print(f"  median block length: same {st.median([float(r[7]) for r in same]):.0f} bp, "
          f"different {st.median([float(r[7]) for r in diff]):.0f} bp")
    print(f"  median n_blocks:     same {st.median([int(r[5]) for r in same]):.0f}, "
          f"different {st.median([int(r[5]) for r in diff]):.0f}")
    # 3. inversion structure of the inverted pairs
    out = []
    typ = lambda x: "P" if x.endswith("_p") else "Q"
    for (P, X), det in inv_detail.items():
        if P == X or P not in comm or X not in comm:
            continue
        tot = sum(d[1] for d in det)
        if tot < 20000:
            continue
        plus = sum(d[1] for d in det if d[2] == "+")
        same_bp = plus if typ(P) == typ(X) else tot - plus
        if same_bp / tot > 0.5:
            continue                      # only the inverted pairs
        big = max(det, key=lambda d: d[1])
        out.append([P, X, len(det), tot, f"{same_bp/tot:.3f}", big[1], f"{big[0]:.3f}"])
    out.sort(key=lambda r: -r[3])
    with open(ROOT / "scripts/cladistics/arm_module_inversions.tsv", "w") as fh:
        fh.write("arm\tpartner\tn_blocks\ttotal_bp\tsame_fraction\tlargest_block_bp\tlargest_block_telomere_offset\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"inverted arm pairs (same-fraction <= 0.5, >=20 kb): {len(out)}")
    for r in out[:6]:
        print(f"  {r[0]:8s} x {r[1]:8s} {r[2]:5d} blocks, {r[3]/1e3:8.0f} kb, same {r[4]}, "
              f"largest block {r[5]/1e3:7.0f} kb at telomere offset {r[6]}")


if __name__ == "__main__":
    main()
