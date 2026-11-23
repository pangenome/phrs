#!/usr/bin/env python3
"""D-PeerQ3: arm-type composition of the 15 arm-level Leiden communities.

For each community of the 41 signal-bearing arms, count within-community arm pairs
by arm type (P-P, Q-Q, P-Q) and compare against a permutation null that keeps the
community sizes and the global count of P and Q arms fixed (arm types are shuffled
across arms, communities fixed). Writes a per-community table and a summary.

Usage:  python3 scripts/cladistics/peerq3_arm_type_asymmetry.py [--nperm 20000]
        [--exclude C7,C14,C15,C6 --out scripts/cladistics/arm_type_asymmetry_r1m7.tsv]
--exclude drops whole communities before counting and before the shuffle (R1-M7).
"""
from __future__ import annotations
import argparse, collections, csv, random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ASSIGN = ROOT / "data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
OUT = ROOT / "scripts/cladistics/peerq3_arm_type_asymmetry.tsv"


def load(exclude=()):
    comm = {}
    with open(ASSIGN) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            comm[r["ChromArm"]] = r["Community"]
    unknown = set(exclude) - set(comm.values())
    if unknown:
        raise SystemExit(f"unknown communities in --exclude: {sorted(unknown)}")
    dropped = sorted(a for a in comm if comm[a] in exclude)
    print(f"arms read: {len(comm)}  dropped (community in --exclude): {len(dropped)} {dropped}")
    comm = {a: c for a, c in comm.items() if c not in exclude}
    arms = sorted(comm)
    typ = {a: ("P" if a.endswith("_p") else "Q") for a in arms}
    return comm, arms, typ


def counts(comm, arms, typ, type_of=None):
    """per-community P-P / Q-Q / P-Q pair counts"""
    by = collections.defaultdict(list)
    for a in arms:
        by[comm[a]].append(a)
    out = {}
    for c, members in by.items():
        t = [type_of[a] if type_of else typ[a] for a in members]
        pp = qq = pq = 0
        for i in range(len(t)):
            for j in range(i + 1, len(t)):
                s = {t[i], t[j]}
                if s == {"P"}:
                    pp += 1
                elif s == {"Q"}:
                    qq += 1
                else:
                    pq += 1
        out[c] = (len(members), pp, qq, pq)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--exclude", default="", help="comma-separated communities to drop")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    comm, arms, typ = load([c for c in args.exclude.split(",") if c])
    obs = counts(comm, arms, typ)
    nP = sum(1 for a in arms if typ[a] == "P")
    nQ = len(arms) - nP
    rng = random.Random(20260)
    perm = collections.defaultdict(lambda: [0, 0, 0])
    gperm = []          # global per-permutation totals
    for _ in range(args.nperm):
        tlist = ["P"] * nP + ["Q"] * nQ
        rng.shuffle(tlist)
        shuffled = dict(zip(arms, tlist))
        c = counts(comm, arms, typ, type_of=shuffled)
        gt = [0, 0, 0]
        for k, (n, pp, qq, pq) in c.items():
            perm[k][0] += pp
            perm[k][1] += qq
            perm[k][2] += pq
            gt[0] += pp; gt[1] += qq; gt[2] += pq
        gperm.append(gt)
    tot = [sum(v[i] for v in obs.values()) for i in range(4)]
    allpairs = len(arms) * (len(arms) - 1) // 2
    same_all = nP * (nP - 1) // 2 + nQ * (nQ - 1) // 2
    print(f"same-type pairs: within communities {tot[1]+tot[2]} of {tot[1]+tot[2]+tot[3]}"
          f" ({(tot[1]+tot[2])/(tot[1]+tot[2]+tot[3]):.3f}); all arm pairs {same_all} of"
          f" {allpairs} ({same_all/allpairs:.3f})")
    print(f"arms: {len(arms)}  P={nP}  Q={nQ}  communities: {len(obs)}")
    print(f"total within-community pairs: {tot[1]+tot[2]+tot[3]}  "
          f"P-P {tot[1]}  Q-Q {tot[2]}  P-Q {tot[3]}")
    exp = [sum(v[i] for v in perm.values()) / args.nperm for i in range(3)]
    print(f"null expectation (arm-type shuffle, {args.nperm} perms): "
          f"P-P {exp[0]:.1f}  Q-Q {exp[1]:.1f}  P-Q {exp[2]:.1f}")
    n = len(gperm)
    pP = sum(1 for g in gperm if g[0] >= tot[1]) / n
    pQ = sum(1 for g in gperm if g[1] >= tot[2]) / n
    pPQ = sum(1 for g in gperm if g[2] <= tot[3]) / n
    print(f"enrichment p (P-P >= obs): {pP:.4f}   (Q-Q >= obs): {pQ:.4f}   "
          f"(P-Q <= obs): {pPQ:.4f}")
    # per-community one-sided p: P-Q pairs <= observed (arm-type clustering)
    perC = collections.defaultdict(list)
    for _ in range(args.nperm):
        tlist = ["P"] * nP + ["Q"] * nQ
        rng.shuffle(tlist)
        c = counts(comm, arms, typ, type_of=dict(zip(arms, tlist)))
        for k, (n, pp, qq, pq) in c.items():
            perC[k].append(pq)
    rows = []
    for c in sorted(obs):
        n, pp, qq, pq = obs[c]
        e = [perm[c][i] / args.nperm for i in range(3)]
        pc = (sum(1 for v in perC[c] if v <= pq) / len(perC[c])) if n > 1 else float("nan")
        rows.append([c, n, pp, qq, pq, f"{e[0]:.2f}", f"{e[1]:.2f}", f"{e[2]:.2f}",
                     "nan" if n <= 1 else f"{pc:.4f}"])
        print(f"  {c:4s} n={n:2d}  P-P {pp:2d} (exp {e[0]:5.2f})  "
              f"Q-Q {qq:2d} (exp {e[1]:5.2f})  P-Q {pq:2d} (exp {e[2]:5.2f})  "
              f"p(P-Q<=obs) {pc:.4f}")
    with open(args.out, "w") as fh:
        fh.write("community\tn_arms\tP-P\tQ-Q\tP-Q\texp_PP\texp_QQ\texp_PQ\tp_PQ_depletion\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
