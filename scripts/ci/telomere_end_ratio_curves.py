#!/usr/bin/env python3
"""Paralog/allelic divergence ratio vs distance from the telomere, PER ARM and PER END.

Variant of telomere_distance_profile.py that keeps the arm identity (and optionally the individual
end) so that curves can be drawn per chromosome end rather than only as a pooled average.

Two aggregation levels, both in one pass:
  ARM  (arm, bin, class)          -> per-arm median ratio curves (histograms retained)
  END  (sample+haplotype+arm, bin, class) -> per-individual-end weighted ratios (sums only)

Distance bins are 2 kb to 50 kb then coarser; a record is counted only when the query and target
midpoints fall in the same bin; paralogs are inter-chromosomal only.

WARNING baked into the output: within roughly the first 6-10 kb the ALLELIC denominator is invalid
(allelic divergence 0.3-2.5% versus 0.09-0.16% in the interior, impossible as orthologous divergence),
so points inside that zone must not be plotted as if they were measurements. The JSON carries the
per-bin allelic value so a plotting script can mask them automatically.

Original docstring follows.
---

PREDICTION UNDER TEST
  If subtelomeric homogenisation is driven by crossover/gene conversion between ends, it should be
  strongest at the telomere-subtelomere junction and decay inward. So the paralog/allelic ratio
  should approach 1 at distance 0 and rise with distance.

WHY THIS IS THE RIGHT TEST
  A PHR is "any region where any other non-homolog appears in the alignment set", so an arm-level
  average can be carried by a single rare match. Binning by distance tests the mechanism instead:
  every end contributes, at every distance, and the ratio is evaluated position by position.

ORIENTATION
  Each window name carries its original contig interval (`<contig>:<start>-<end>_<chrN>_<p|q>arm`).
  The interval endpoint that coincides with a contig end is the telomere, so:
      p arm -> telomere at sequence position 0      distance = position
      q arm -> telomere at sequence position L      distance = L - position
  (Windows are stored in original contig orientation; the earlier dot-plot work had to
  reverse-complement q arms explicitly, which confirms this.)

BINNING
  A record is assigned by the MIDPOINT of each alignment span, and only when the query and the
  target midpoints fall in the same distance bin - so "paralog divergence at 10-20 kb" means both
  sides are 10-20 kb from their own telomeres. Records are short blocks (a few kb typically); this
  is a stated approximation at block granularity.

OUTPUT
  Per distance bin: allelic (same arm) and paralog (different chromosome) weighted divergence,
  per-record medians, the ratio, and the number of contributing ends.

Usage:
  python3 scripts/ci/telomere_distance_profile.py --limit 3000 --workers 12 --out /tmp/prof
  python3 scripts/ci/telomere_distance_profile.py --workers 14 --out PREFIX
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
SAMPLE_RE = re.compile(r"^([^#]+)#([^#]+)#")
NAME_RE = re.compile(r"^(.*?):(\d+)-(\d+)_(chr[0-9XY]+)_(parm|qarm)$")
CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
CG_TAG = re.compile(r"cg:Z:([^\t\n]+)")

# distance bins in bp (upper bounds), then everything above the last
# fine near the junction (2 kb), coarser inward so every bin stays measurable
EDGES = list(range(2_000, 50_001, 2_000)) + [75_000, 100_000, 150_000, 200_000, 500_000]
NBINS = len(EDGES)
BINS_PER_UNIT = 10000
CHIST = 5000
OUT = "telomere_end_ratio_curves"
MU = 1.2e-8
YPG = 20.0


def parse_name(name):
    """-> (arm, telomere_side) where telomere_side is 'start' for p arms, 'end' for q arms."""
    m = NAME_RE.match(name)
    if not m:
        return None, None
    chrom, kind = m.group(4), m.group(5)
    arm = chrom + ("_p" if kind == "parm" else "_q")
    return arm, ("start" if kind == "parm" else "end")


def dist_mid(pos_mid, seqlen, side):
    d = pos_mid if side == "start" else seqlen - pos_mid
    return max(0.0, d)


def binof_dist(d):
    for i, e in enumerate(EDGES):
        if d < e:
            return i
    return NBINS - 1


def binof_div(d):
    return CHIST if d >= 0.5 else int(d * BINS_PER_UNIT)


def parse_cigar(cg):
    m = x = 0
    for n, op in CIGAR_RE.findall(cg):
        n = int(n)
        if op in "=M":
            m += n
        elif op == "X":
            x += n
    return m, x


def _lines(path):
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
    """{bin_index: {class: [n, mx, x, hist, ends_set_is_not_picklable -> count]}}"""
    out = defaultdict(lambda: [0, 0, 0, [0] * (CHIST + 1)])
    for line in _lines(path):
        f = line.split("\t", 13)
        if len(f) < 13 or "cg:Z:" not in f[12] + (f[13] if len(f) > 13 else ""):
            continue
        qn, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        tn, tlen, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
        q_arm, q_side = parse_name(qn)
        t_arm, t_side = parse_name(tn)
        if not q_arm or not t_arm:
            continue
        dq = dist_mid((qs + qe) / 2.0, qlen, q_side)
        dt = dist_mid((ts + te) / 2.0, tlen, t_side)
        b = binof_dist(dq)
        if b != binof_dist(dt):
            continue                      # both sides must sit in the same distance bin
        mc = CG_TAG.search(line)
        if mc is None:
            continue
        m, x = parse_cigar(mc.group(1))
        mx = m + x
        if mx < 1000:
            continue
        cls = "allelic" if q_arm == t_arm else "paralog"
        if cls == "paralog" and q_arm.rsplit("_", 1)[0] == t_arm.rsplit("_", 1)[0]:
            continue                      # inter-chromosomal paralogs only
        d = x / mx
        a = out[("ARM", t_arm, b, cls)]
        a[0] += 1; a[1] += mx; a[2] += x; a[3][binof_div(d)] += mx
        ms = SAMPLE_RE.match(tn)
        endkey = (ms.group(1) + "#" + ms.group(2)) if ms else "unknown"
        e = out[("END", endkey, t_arm, b, cls)]
        e[0] += 1; e[1] += mx; e[2] += x
        continue
        d = x / mx
        a[0] += 1
        a[1] += mx
        a[2] += x
        a[3][binof_div(d)] += mx
    return dict(out)


def wquant(hist, q):
    tot = sum(hist)
    if tot == 0:
        return None
    tgt, run = q * tot, 0
    for i, c in enumerate(hist):
        run += c
        if run >= tgt:
            return i / BINS_PER_UNIT
    return CHIST / BINS_PER_UNIT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="scripts/ci/telomere_distance_profile")
    global OUT
    a = ap.parse_args()
    OUT = a.out
    files = sorted(glob.glob(os.path.join(PAFDIR, "*.paf.gz")))
    if a.limit:
        import random
        random.seed(11)
        files = random.sample(files, min(a.limit, len(files)))
    t0 = time.time()
    print(f"{len(files)} PAFs, {a.workers} workers", file=sys.stderr)
    tot = defaultdict(lambda: [0, 0, 0, [0] * (CHIST + 1)])
    with Pool(a.workers) as pool:
        for i, part in enumerate(pool.imap_unordered(process, files, chunksize=4), 1):
            for k, v in part.items():
                t = tot[k]
                t[0] += v[0]
                t[1] += v[1]
                t[2] += v[2]
                for j, c in enumerate(v[3]):
                    t[3][j] += c
            if i % 400 == 0:
                print(f"  {i}/{len(files)} {time.time()-t0:.0f}s", file=sys.stderr)
                sys.stderr.flush()
    print(f"done {time.time()-t0:.0f}s", file=sys.stderr)

    def med(v):
        return wquant(v[3], 0.5) if v[1] else None

    arms = sorted({k[1] for k in tot if k[0] == "ARM"})
    out_arms = {}
    for arm in arms:
        rows = []
        for b in range(NBINS):
            lo = 0 if b == 0 else EDGES[b - 1]
            a_ = tot.get(("ARM", arm, b, "allelic"))
            p_ = tot.get(("ARM", arm, b, "paralog"))
            row = {"bin": b, "from_kb": lo / 1000, "to_kb": EDGES[b] / 1000,
                   "allelic_median": med(a_) if a_ else None,
                   "paralog_median": med(p_) if p_ else None,
                   "allelic_bp": a_[1] if a_ else 0, "paralog_bp": p_[1] if p_ else 0}
            if row["allelic_median"] and row["paralog_median"]:
                row["ratio"] = row["paralog_median"] / row["allelic_median"]
            rows.append(row)
        out_arms[arm] = rows
    ends = {}
    for k, v in tot.items():
        if k[0] != "END" or not v[1]:
            continue
        ends.setdefault(k[1] + "|" + k[2], {}).setdefault(k[3], {})[k[4]] = \
            [v[0], v[1], v[2]]
    end_rows = {}
    for eid, bins in ends.items():
        rows = []
        for b in sorted(bins):
            a_ = bins[b].get("allelic"); p_ = bins[b].get("paralog")
            if not a_ or not p_ or not a_[1] or not p_[1]:
                continue
            rows.append({"bin": b, "from_kb": (0 if b == 0 else EDGES[b-1])/1000,
                         "to_kb": EDGES[b]/1000,
                         "allelic": a_[2]/a_[1], "paralog": p_[2]/p_[1],
                         "ratio": (p_[2]/p_[1]) / (a_[2]/a_[1]) if a_[2] else None,
                         "paralog_bp": p_[1]})
        if len(rows) >= 3:
            end_rows[eid] = rows
    import json as _json
    _json.dump({"arms": out_arms, "ends": end_rows, "n_files": len(files),
                "edges_kb": [e/1000 for e in EDGES]}, open(OUT + ".json", "w"))
    print(f"[written] {OUT}.json  arms={len(out_arms)} ends={len(end_rows)}", file=sys.stderr)
    print(f"{'arm':9s} {'usable bins':>11s} {'median ratio >=10kb':>20s}")
    import statistics as _st
    for arm in arms:
        rr = [r["ratio"] for r in out_arms[arm]
              if r.get("ratio") and r["from_kb"] >= 10]
        if rr:
            print(f"{arm:9s} {len(rr):11d} {_st.median(rr):20.2f}")


if __name__ == "__main__":
    main()
