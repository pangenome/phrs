#!/usr/bin/env python3
"""Assembled telomere length per chromosome end and its relation to subtelomeric exchange (R1-m12).

Telomere length of a window = the seqtk-telo interval that the 500 kb telomere-anchored window was cut next
to: for a p-arm window starting at contig position a (1-based), the interval [0, a-1) of that contig; for a
q-arm window ending at b, the interval [b, contig length). Windows with no adjacent interval are counted
and reported. Per arm: number of haplotypes and median telomere length. Unit of every test is the arm.
Exchange measures per arm, from the paralogous/allelic divergence table (Extended Data Table 1 input):
  ratio = divergence to the closest arm on another chromosome / allelic divergence of the arm.
We report the Spearman correlation between per-arm median telomere length and that ratio, with and
without the acrocentric, X and Y ends (excluded by Karimian et al. 2024), and the per-arm ranks of the
ends that study found shortest (17p, 20q, 12p) and longest (4q, 12q, 3p).
"""
import argparse, csv, glob, json, os, re, statistics as st, sys
from collections import defaultdict

NAME = re.compile(r"^(.+):(\d+)-(\d+)_(chr[0-9XYM]+)_(parm|qarm)$")


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return num / (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5


def perm_p(x, y, n=10000, seed=7):
    import random
    rng, r0, y2, hit = random.Random(seed), abs(spearman(x, y)), list(y), 0
    for _ in range(n):
        rng.shuffle(y2)
        hit += abs(spearman(x, y2)) >= r0
    return (1 + hit) / (1 + n)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--telo-dir", required=True)
ap.add_argument("--fai", required=True)
ap.add_argument("--het", required=True, help="hprc_het_final.json")
ap.add_argument("--communities", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
telo = defaultdict(list)
files = sorted(glob.glob(os.path.join(a.telo_dir, "*.telo.tsv")))
for f in files:
    for l in open(f):
        c, s, e, L = l.split()[:4]
        telo[c].append((int(s), int(e), int(L)))
sys.stderr.write(f"{len(files)} telo files, {sum(map(len, telo.values()))} intervals on {len(telo)} contigs\n")
per_arm, missing = defaultdict(list), 0
for l in open(a.fai):
    m = NAME.match(l.split("\t")[0])
    if not m:
        raise SystemExit(f"FATAL: cannot parse {l.split()[0]!r}")
    contig, s1, e1, ch, pq = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5)
    arm = ch + ("_p" if pq == "parm" else "_q")
    iv = [x for x in telo.get(contig, []) if (pq == "parm" and x[0] == 0 and x[1] == s1 - 1)
          or (pq == "qarm" and x[1] == x[2] and x[0] == e1)]
    if len(iv) != 1:
        missing += 1
        continue
    per_arm[arm].append(iv[0][1] - iv[0][0])
sys.stderr.write(f"windows with one adjacent telomere interval: {sum(map(len, per_arm.values()))}; without: {missing}\n")
comm = {r["ChromArm"]: r["Community"] for r in csv.DictReader(open(a.communities), delimiter="\t")}
d = json.load(open(a.het))
ale, para = {}, defaultdict(list)
for r in d["pairs"]:
    x, y, v = r["arm_a"], r["arm_b"], float(r["d_median_record"])
    if x == y:
        ale[x] = v
    elif x.rsplit("_", 1)[0] != y.rsplit("_", 1)[0]:
        para[x].append(v); para[y].append(v)
rows = []
for arm, v in per_arm.items():
    ratio = min(para[arm]) / ale[arm] if para.get(arm) and ale.get(arm) else None
    rows.append(dict(arm=arm, community=comm.get(arm, "none"), n=len(v), median_telomere_bp=int(st.median(v)), ratio=ratio))
rows.sort(key=lambda r: r["median_telomere_bp"])
for i, r in enumerate(rows, 1):
    r["rank_short_to_long"] = i
with open(a.out, "w") as fh:
    cols = ["arm", "community", "n", "median_telomere_bp", "rank_short_to_long", "ratio"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join("NA" if r[c] is None else (f"{r[c]:.2f}" if c == "ratio" else str(r[c])) for c in cols) + "\n")
excl = lambda arm: arm in ("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p") or arm.startswith(("chrX", "chrY"))
for lab, sel in (("all arms with a ratio", [r for r in rows if r["ratio"] is not None]),
                 ("without acrocentric p, X and Y", [r for r in rows if r["ratio"] is not None and not excl(r["arm"])])):
    x, y = [r["median_telomere_bp"] for r in sel], [r["ratio"] for r in sel]
    print(f"{lab}: {len(sel)} arms; Spearman(median telomere length, closest-paralog/allelic ratio) = "
          f"{spearman(x, y):.3f}, permutation p = {perm_p(x, y):.3g}")
kept = [r for r in rows if not excl(r["arm"])]
kr = {r["arm"]: i for i, r in enumerate(kept, 1)}
print(f"ranks among {len(kept)} arms without acrocentric p, X and Y (1 = shortest): "
      + ", ".join(f"{e} {kr.get(e, 'NA')}" for e in ("chr17_p", "chr20_q", "chr12_p", "chr4_q", "chr12_q", "chr3_p")))
print("shortest 6: " + ", ".join(f"{r['arm']} ({r['median_telomere_bp']})" for r in kept[:6]))
print("longest 6: " + ", ".join(f"{r['arm']} ({r['median_telomere_bp']})" for r in kept[-6:]))
