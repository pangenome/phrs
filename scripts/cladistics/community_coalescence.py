#!/usr/bin/env python3
"""Best-guess coalescence time for each arm-level subtelomeric community.

METHOD (deliberately crude; stated as a bound, not a measurement)
-----------------------------------------------------------------
For two contemporaneous copies of a subtelomeric module in one genome, the substitution-only
divergence between them is d_sub = X / (= + X), taken over the aligned homologous blocks.
If the two copies last shared an ancestral sequence T generations ago and have accumulated
substitutions independently at rate mu per bp per generation, then

    d_sub ~= 2 * mu * T        =>        T = d_sub / (2 * mu)

mu = 1.2e-8 substitutions/bp/generation (human pedigree average). T is reported in generations
and in Myr at 20 years/generation.

WHY THIS IS AN UPPER BOUND
  Subtelomeric and pericentromeric sequence sits in the male germline, where the mutation rate is
  ~2-4x the genome average, and repeat-rich sequence accumulates substitutions faster still.
  Using a genome-average mu therefore *overestimates* T: with mu = 1.5e-8 or 3e-8 the same
  divergence gives a 20-60% smaller T. The published mutation rate is the conservative choice.

WHY IT IS ALSO CRUDE IN THE OTHER DIRECTION
  Only blocks that align at >= the wfmash identity floor contribute, so highly diverged pairs are
  measured on their most similar blocks only (d_sub biased low -> T biased low). Indels are
  excluded from d_sub by construction (they are not substitutions and mu is a substitution rate),
  which is the right thing for dating but means the number describes only substitutional change.
  Recurrent gene conversion resets segments independently, so T is really "generations since the
  last exchange affecting these blocks" — the distribution across pairs, not any single value,
  is the meaningful object.

INPUTS
  A FASTA of one individual's telomere-anchored 500 kb windows, one per chromosome arm, named
  <arm> (chrN_p / chrN_q), and a wfmash all-vs-all PAF with CIGAR (`cg:Z:`) for that FASTA.
  Reference workflow:

    samtools faidx combined_with_hprcv2.1Mb.telo_500kb_trimmed.fa.gz -r <names for one genome>
    # rename headers to arm names, telomere-first orientation not required (wfmash is strand-aware)
    samtools faidx hg002.win.fa
    wfmash hg002.win.fa hg002.win.fa -p 75 -s 1000 -w 3000 -t 8 -o > hg002.p75.paf

  Community assignments:
    /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv

Usage:
  python3 scripts/cladistics/community_coalescence.py --paf <paf> --fa <fasta> [--min-aln 20000]
"""
import argparse
import collections
import csv
import json
import os
import random
import re
import statistics as st

K15 = ("/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
       "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
MU = 1.2e-8
MU_FAST = 3.0e-8          # plausible subtelomeric/paternal-germline rate
YEARS_PER_GEN = 20
CIGAR_RE = re.compile(r"(\d+)([=XIDM])")


def cigar_counts(cg):
    m = x = i = d = ev = 0
    for n, op in CIGAR_RE.findall(cg):
        n = int(n)
        if op in "=M":
            m += n
        elif op == "X":
            x += n
        elif op == "I":
            i += n
            ev += 1
        elif op == "D":
            d += n
            ev += 1
    return m, x, i, d, ev


def load_pairs(paf, min_aln):
    """Return {(a,b): dict} with accumulated CIGAR counts for a<b, one entry per alignment block
    also kept individually so we can median over blocks rather than over pairs."""
    pairs = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
    blocks = collections.defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, t = f[0], f[5]
            if q == t:
                continue
            cg = [z[5:] for z in f if z.startswith("cg:Z:")]
            if not cg:
                continue
            m, x, i, d, ev = cigar_counts(cg[0])
            if m + x < min_aln:
                continue
            key = (q, t) if q < t else (t, q)
            a = pairs[key]
            for k, v in enumerate((m, x, i, d, ev)):
                a[k] += v
            if x + m > 0:
                blocks[key].append((m, x))
    return pairs, blocks


def boot_median(vals, n=2000, seed=7):
    if len(vals) < 2:
        return (None, None)
    rnd = random.Random(seed)
    meds = []
    for _ in range(n):
        meds.append(st.median([vals[rnd.randrange(len(vals))] for _ in vals]))
    meds.sort()
    return meds[int(0.025 * n)], meds[int(0.975 * n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--fa", required=True)
    ap.add_argument("--min-aln", type=int, default=20000)
    ap.add_argument("--out", default="scripts/cladistics/community_coalescence")
    a = ap.parse_args()

    arms = [l[1:].split()[0] for l in open(a.fa) if l.startswith(">")]
    pairs, blocks = load_pairs(a.paf, a.min_aln)
    comm = {}
    with open(K15) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            comm[r["ChromArm"].strip()] = r["Community"].strip()

    per_pair = []
    for (q, t), (m, x, i, d, ev) in pairs.items():
        alen = m + x
        if alen == 0:
            continue
        per_pair.append({
            "arm_a": q, "arm_b": t,
            "aligned_bp": alen,
            "substitutions": x,
            "d_sub": x / alen,
            "d_lenpen": (x + i + d) / (m + x + i + d) if (m + x + i + d) else None,
            "d_gapcomp": (x + ev) / (m + x + ev) if (m + x + ev) else None,
            "T_gen": (x / alen) / (2 * MU),
            "same_community": comm.get(q) == comm.get(t),
            "community": comm.get(q) if comm.get(q) == comm.get(t) else None,
        })

    rows = []
    for c in sorted({v for v in comm.values()}, key=lambda z: int(z[1:])):
        members = [x for x in arms if comm.get(x) == c]
        ds = [p["d_sub"] for p in per_pair if p["same_community"] and p["community"] == c]
        if not ds:
            rows.append({"community": c, "n_arms": len(members), "arms": ",".join(members),
                         "n_pairs": 0, "status": "no pair aligned >= %d bp" % a.min_aln})
            continue
        md = st.median(ds)
        lo, hi = boot_median(ds)
        rows.append({
            "community": c, "n_arms": len(members), "arms": ",".join(members),
            "n_pairs": len(ds), "status": "ok",
            "median_d_sub": f"{md:.6f}",
            "min_d_sub": f"{min(ds):.6f}", "max_d_sub": f"{max(ds):.6f}",
            "median_T_gen": f"{md / (2 * MU):,.0f}",
            "T_ci95_gen": (f"{lo / (2 * MU):,.0f}-{hi / (2 * MU):,.0f}" if lo is not None else ""),
            "median_T_Myr": f"{md / (2 * MU) * YEARS_PER_GEN / 1e6:.2f}",
            "T_Myr_if_mu3e-8": f"{md / (2 * MU_FAST) * YEARS_PER_GEN / 1e6:.2f}",
        })

    base = a.out
    with open(base + "_pairs.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_pair[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(per_pair, key=lambda p: p["d_sub"]))
    with open(base + "_communities.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"{'comm':5s} {'n_arms':>6s} {'n':>3s} {'med d_sub':>10s} "
          f"{'T (Myr)':>8s} {'if mu=3e-8':>11s}  arms")
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['community']:5s} {r['n_arms']:6d} {0:3d} {'-':>10s} {'-':>8s} {'-':>11s}  "
                  f"{r['arms']}  ({r['status']})")
            continue
        print(f"{r['community']:5s} {r['n_arms']:6d} {r['n_pairs']:3d} {r['median_d_sub']:>10s} "
              f"{r['median_T_Myr']:>8s} {r['T_Myr_if_mu3e-8']:>11s}  {r['arms']}")
    ok = [float(r["median_d_sub"]) for r in rows if r["status"] == "ok"]
    if ok:
        m = st.median(ok)
        print(f"\nacross {len(ok)} communities with data: median d_sub={m:.6f} -> "
              f"T={m / (2 * MU) * YEARS_PER_GEN / 1e6:.2f} Myr "
              f"({m / (2 * MU_FAST) * YEARS_PER_GEN / 1e6:.2f} Myr if mu=3e-8)")
    within = [p["d_sub"] for p in per_pair if p["same_community"]]
    if within:
        print(f"all {len(within)} within-community pairs: median d_sub={st.median(within):.6f} -> "
              f"T={st.median(within) / (2 * MU) * YEARS_PER_GEN / 1e6:.2f} Myr")
    between = [p["d_sub"] for p in per_pair if not p["same_community"]]
    if between:
        print(f"all {len(between)} between-community pairs: median d_sub={st.median(between):.6f} -> "
              f"T={st.median(between) / (2 * MU) * YEARS_PER_GEN / 1e6:.2f} Myr  (a lower bound on the"
              f" age of the duplication that created the families)")
    json.dump({"per_pair": per_pair, "per_community": rows}, open(base + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
