#!/usr/bin/env python3
"""Population-scale divergence between heterologous subtelomeric windows (the PHR all-vs-all).

DATA — this is the alignment set the PHR calls, arm communities and the similarity matrices
were built from, not a whole-assembly alignment:
  /moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/*.paf.gz
  18,827 target-centric PAFs (one per telomere-anchored 500 kb window instance, 465 assemblies),
  wfmash -p 95, CIGARs included, ~92 GB gzipped. Each sequence in the PAF *is* a window, named
  ..._<chrN>_<parm|qarm>, so arm identity comes straight from the name and no coordinate
  clipping is needed.

METRIC
  d_sub = X / (= + X) from the CIGAR: substitutions per aligned non-indel position. Indels are
  excluded because mu is a substitution rate; they are reported separately via aligned length.
  異  Heterologs = query and target on different chromosomes (non-allelic arm windows).

DEDUPE
  Every unordered window pair appears twice (A->B and B->A). A record is counted only when
  query_name < target_name, which is exact.

Usage:
  python3 scripts/ci/hprc_heterolog_divergence.py --limit 30 --workers 8   # smoke test
  python3 scripts/ci/hprc_heterolog_divergence.py --workers 14 --out PREFIX
"""
import argparse
import glob
import gzip
import json
import os
import re
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

PAFDIR = "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95"
CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
WIN_RE = re.compile(r"(chr[0-9XY]+)_(parm|qarm)$")
CG_TAG = re.compile(r"cg:Z:([^\t\n]+)")
GI_TAG = re.compile(r"gi:f:([0-9.eE+-]+)")
NBINS = 500
MU = 1.2e-8
YEARS = 20


def arm_of(name):
    m = WIN_RE.search(name)
    if not m:
        return None, None
    chrom, kind = m.group(1), m.group(2)
    return chrom, chrom + ("_p" if kind == "parm" else "_q")


def parse_cigar(cg):
    m = x = 0
    for n, op in CIGAR_RE.findall(cg):
        n = int(n)
        if op in "=M":
            m += n
        elif op == "X":
            x += n
    return m, x


def binof(d):
    return NBINS if d >= 0.5 else int(d * 1000)


def _lines(path):
    """Yield lines, tolerating plain-text or truncated gzip files."""
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                yield line
        return
    except (OSError, EOFError, gzip.BadGzipFile):
        pass
    try:
        with open(path, "rt", errors="replace") as fh:
            for line in fh:
                yield line
    except OSError:
        return


def process(path):
    out = defaultdict(lambda: [0, 0, 0, [0] * (NBINS + 1), 0, 1.0, 0.0,
                                [0] * (NBINS + 1), 0])
    for line in _lines(path):
        f = line.split("\t", 13)
        if len(f) < 13:
            continue
        qn, tn = f[0], f[5]
        if qn >= tn:
            continue
        cq, q_arm = arm_of(qn)
        ct, t_arm = arm_of(tn)
        if not q_arm or not t_arm or cq == ct:
            continue                        # heterologs only
        mc = CG_TAG.search(line)
        if mc is None:
            continue
        m, x = parse_cigar(mc.group(1))
        mx = m + x
        if mx < 2000:
            continue
        key = (q_arm, t_arm) if q_arm < t_arm else (t_arm, q_arm)
        a = out[key]
        d = x / mx
        a[0] += 1
        a[1] += mx
        a[2] += x
        a[3][binof(d)] += mx
        if d <= 0.005:
            a[4] += mx
        if d < a[5]:
            a[5] = d
        mg = GI_TAG.search(line)
        if mg is not None:
            dg = 1.0 - float(mg.group(1))
            a[6] += dg * mx
            a[7][binof(dg)] += mx
        else:
            a[8] += 1
    return dict(out)


def merge(parts):
    tot = defaultdict(lambda: [0, 0, 0, [0] * (NBINS + 1), 0, 1.0, 0.0,
                                [0] * (NBINS + 1), 0])
    for p in parts:
        for k, v in p.items():
            t = tot[k]
            t[0] += v[0]
            t[1] += v[1]
            t[2] += v[2]
            for i, c in enumerate(v[3]):
                t[3][i] += c
            t[4] += v[4]
            if v[5] < t[5]:
                t[5] = v[5]
            t[6] += v[6]
            for i, c in enumerate(v[7]):
                t[7][i] += c
            t[8] += v[8]
    return tot


def wquant(hist, q):
    tot = sum(hist)
    if tot == 0:
        return None
    target, run = q * tot, 0
    for i, c in enumerate(hist):
        run += c
        if run >= target:
            return i / 1000.0
    return NBINS / 1000.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="scripts/ci/hprc_heterolog_divergence")
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(PAFDIR, "*.paf.gz")))
    if a.limit:
        import random
        random.seed(7)
        files = random.sample(files, min(a.limit, len(files)))
    t0 = time.time()
    print(f"{len(files)} window PAFs, {a.workers} workers", file=sys.stderr)
    parts = []
    with Pool(a.workers) as pool:
        for i, p in enumerate(pool.imap_unordered(process, files, chunksize=4), 1):
            parts.append(p)
            if i % 200 == 0:
                print(f"  {i}/{len(files)}  {time.time() - t0:.0f}s", file=sys.stderr)
                sys.stderr.flush()
    tot = merge(parts)
    print(f"done {time.time() - t0:.0f}s, {len(tot)} heterologous arm pairs", file=sys.stderr)
    rows = []
    for (x, y), v in tot.items():
        n, mx, xx, hist, bp995, mind, gisum, ghist, n_nogi = v
        if mx < 20000:
            continue
        rows.append({"arm_a": x, "arm_b": y, "n_records": n, "aligned_bp": mx,
                     "d_sub": xx / mx,
                     "d_median_record": wquant(hist, 0.5),
                     "d_q25_record": wquant(hist, 0.25),
                     "d_q75_record": wquant(hist, 0.75),
                     "frac_at_99.5pct": bp995 / mx,
                     "min_record_d": mind,
                     "d_gapcomp": gisum / mx if mx else None,
                     "d_gapcomp_median_record": wquant(ghist, 0.5),
                     "T_Myr_gapcomp": (gisum / mx) / (2 * MU) * YEARS / 1e6 if mx else None,
                     "T_Myr": xx / mx / (2 * MU) * YEARS / 1e6,
                     "T_Myr_at_min": mind / (2 * MU) * YEARS / 1e6})
    rows.sort(key=lambda r: r["d_sub"])
    if not rows:
        print("no pairs aggregated", file=sys.stderr)
        return
    with open(a.out + ".tsv", "w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])
                               for k in r) + "\n")
    allhist = [0] * (NBINS + 1)
    for v in tot.values():
        for i, c in enumerate(v[3]):
            allhist[i] += c
    s = {"n_arm_pairs": len(rows),
         "global_d_sub": sum(v[2] for v in tot.values()) / max(1, sum(v[1] for v in tot.values())),
         "global_d_median_record": wquant(allhist, 0.5),
         "global_d_q05": wquant(allhist, 0.05),
         "global_d_q25": wquant(allhist, 0.25),
         "global_d_q75": wquant(allhist, 0.75),
         "global_d_q95": wquant(allhist, 0.95),
         "global_frac_at_99.5pct": sum(v[4] for v in tot.values()) / max(1, sum(v[1] for v in tot.values())),
         "n_files": len(files), "mu": MU}
    s["global_T_Myr"] = s["global_d_median_record"] / (2 * MU) * YEARS / 1e6
    json.dump({"summary": s, "pairs": rows}, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(s, indent=1))
    print("\nclosest 15 heterologous arm pairs across the cohort:")
    for r in rows[:15]:
        print(f"{r['arm_a'] + '-' + r['arm_b']:20s} {r['aligned_bp']:12d} {r['d_sub']:9.5f} "
              f"{r['d_median_record']:9.5f} {r['frac_at_99.5pct']:8.3f} {r['T_Myr']:7.2f}")


if __name__ == "__main__":
    main()
