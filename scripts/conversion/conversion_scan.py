#!/usr/bin/env python3
"""Cohort scan for gene-conversion tracts in telomere-anchored subtelomeric windows (R1-M1, test 3b).

INPUT: target-centric PAFs of the PHR all-vs-all (wfmash -p 95, cg:Z CIGAR), one file per haplotype
window, e.g. /moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/*.paf.gz. Every query in a file is
another window aligned onto the target window; names end in _chrN_parm or _chrN_qarm.

PER BIN (default 2 kb of target coordinates) we keep, for every query, the record with the most aligned
target positions (= plus X) in the bin, and its substitution divergence d = X / (= + X). Queries with
less than --min-cov of the bin aligned are ignored in that bin.
  alleles   = queries on the same arm as the target (other haplotypes, the target itself excluded)
  paralogs  = queries on any other arm (other chromosome, or the other arm of the same chromosome)
  dA        = median d over the alleles; dY = median d over the copies of paralogous arm Y, for every
              arm Y with >= --min-arm-copies copies in the bin; Y* = the arm with the smallest dY
A bin is CONVERTED when it has >= --min-alleles alleles and dA - dY* >= --min-delta: the target is
typically closer to the copies of another arm than to the copies of its own arm. Medians on both sides
keep the comparison fair: a minimum over thousands of paralogous copies is ~0 by chance in shared
sequence, a median is not. A conversion carried by fewer than half of the haplotypes still moves dA;
a donor sequence carried by fewer than half of the copies of Y does not move dY and is missed.
TRACT = a maximal run of converted bins. A tract is EMBEDDED when the nearest scored bin on each side
(>= --min-alleles alleles) is not converted, i.e. the target matches its own arm on both sides; tracts at
the window edge or next to unscored bins are reported with embedded = 0.

OUTPUT: <out>.bins.tsv (one row per scored bin) and <out>.tracts.tsv (one row per tract); the log on
stderr gives inputs, parameters, record counts kept and dropped with the reason, and timings.
"""
import argparse, gzip, re, sys, time
from collections import Counter, defaultdict

CIGAR = re.compile(r"(\d+)([=XIDM])")
ARM = re.compile(r"_(chr[0-9XYM]+)_(parm|qarm)$")


def log(msg):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def arm_of(name):
    m = ARM.search(name)
    if not m:
        raise SystemExit(f"FATAL: no _chrN_parm/_qarm suffix in sequence name {name!r}")
    return m.group(1), m.group(1) + "_" + m.group(2)


def bin_counts(tstart, cigar, binsize):
    """{bin: [matches, mismatches]} along target coordinates; I consumes no target, D consumes target unaligned."""
    out = defaultdict(lambda: [0, 0])
    pos = tstart
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op == "M":
            raise SystemExit("FATAL: CIGAR has M; need =/X to count substitutions")
        if op == "I":
            continue
        if op == "D":
            pos += n
            continue
        k = 0 if op == "=" else 1
        while n:
            b = pos // binsize
            take = min(n, (b + 1) * binsize - pos)
            out[b][k] += take
            pos += take
            n -= take
    return out


def med(v):
    v = sorted(v)
    return v[len(v) // 2]


def scan(lines, binsize=2000, min_cov=0.75, min_alleles=20, min_arm_copies=20, min_delta=0.0025):
    stats = Counter()
    target = tlen = None
    best = defaultdict(dict)  # bin -> query -> (aligned, mismatches)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            stats["dropped_short_line"] += 1
            continue
        q, t = f[0], f[5]
        if target is None:
            target, tlen = t, int(f[6])
        elif t != target:
            raise SystemExit(f"FATAL: second target {t!r} in a target-centric PAF (first {target!r})")
        if q == t:
            stats["dropped_self"] += 1
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if cg is None:
            stats["dropped_no_cigar"] += 1
            continue
        stats["records_kept"] += 1
        for b, (m, x) in bin_counts(int(f[7]), cg, binsize).items():
            a = m + x
            if a > best[b].get(q, (0, 0))[0]:
                best[b][q] = (a, x)
    if target is None:
        raise SystemExit("FATAL: no records")
    tchr, tarm = arm_of(target)
    nb = (tlen + binsize - 1) // binsize
    bins = []
    for b in range(nb):
        size = min(binsize, tlen - b * binsize)
        by_arm = defaultdict(list)
        for q, (a, x) in best.get(b, {}).items():
            if a >= min_cov * size:
                by_arm[arm_of(q)[1]].append(x / a)
        al = by_arm.pop(tarm, [])
        row = {"bin": b, "start": b * binsize, "end": b * binsize + size, "n_alleles": len(al),
               "n_paralogs": sum(map(len, by_arm.values())), "scored": int(len(al) >= min_alleles)}
        if al:
            row["d_allele_median"] = med(al)
        if row["n_paralogs"]:
            row["d_paralog_min"] = min(min(v) for v in by_arm.values())
        arms = [(med(v), y, len(v)) for y, v in by_arm.items() if len(v) >= min_arm_copies]
        if arms:
            row["d_arm_median"], row["paralog_arm"], row["n_arm"] = min(arms)
        row["converted"] = int(row["scored"] and bool(arms) and
                               row["d_allele_median"] - row["d_arm_median"] >= min_delta)
        bins.append(row)
    stats["bins_total"], stats["bins_scored"] = nb, sum(r["scored"] for r in bins)
    stats["bins_converted"] = sum(r["converted"] for r in bins)
    tracts, b = [], 0
    while b < nb:
        if not bins[b]["converted"]:
            b += 1
            continue
        e = b
        while e + 1 < nb and bins[e + 1]["converted"]:
            e += 1
        left = next((r for r in reversed(bins[:b]) if r["scored"]), None)
        right = next((r for r in bins[e + 1:] if r["scored"]), None)
        emb = int(left is not None and right is not None and left["bin"] == b - 1 and right["bin"] == e + 1)
        run = bins[b:e + 1]
        tracts.append({"target": target, "arm": tarm, "start": run[0]["start"], "end": run[-1]["end"],
                       "n_bins": len(run), "embedded": emb,
                       "paralog_arm": Counter(r["paralog_arm"] for r in run).most_common(1)[0][0],
                       "d_allele_median": med(r["d_allele_median"] for r in run),
                       "d_arm_median": med(r["d_arm_median"] for r in run),
                       "d_paralog_min": min(r["d_paralog_min"] for r in run)})
        b = e + 1
    stats["tracts"], stats["tracts_embedded"] = len(tracts), sum(t["embedded"] for t in tracts)
    return target, tarm, bins, tracts, stats


BIN_COLS = ["start", "end", "n_alleles", "n_paralogs", "scored", "d_allele_median", "paralog_arm", "n_arm",
            "d_arm_median", "d_paralog_min", "converted"]
TRACT_COLS = ["target", "arm", "start", "end", "n_bins", "embedded", "paralog_arm", "d_allele_median",
              "d_arm_median", "d_paralog_min"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pafs", nargs="*", help="target-centric PAF(.gz) files")
    ap.add_argument("--list", help="file with one PAF path per line; with --chunk N/K process chunk N of K")
    ap.add_argument("--chunk", help="N/K, 1-based")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--bin", type=int, default=2000)
    ap.add_argument("--min-cov", type=float, default=0.75)
    ap.add_argument("--min-alleles", type=int, default=20)
    ap.add_argument("--min-arm-copies", type=int, default=20)
    ap.add_argument("--min-delta", type=float, default=0.0025, help="dA - dY* in substitutions per site")
    a = ap.parse_args()
    paths = list(a.pafs)
    if a.list:
        paths += [l.strip() for l in open(a.list) if l.strip()]
    if a.chunk:
        n, k = map(int, a.chunk.split("/"))
        paths = paths[n - 1::k]
    if not paths:
        raise SystemExit("FATAL: no input PAFs")
    log(f"inputs: {len(paths)} PAFs; bin {a.bin} bp, min_cov {a.min_cov}, min_alleles {a.min_alleles}, "
        f"min_arm_copies {a.min_arm_copies}, min_delta {a.min_delta}; out {a.out}")
    tot = Counter()
    with open(a.out + ".bins.tsv", "w") as fb, open(a.out + ".tracts.tsv", "w") as ft:
        fb.write("target\tarm\t" + "\t".join(BIN_COLS) + "\n")
        ft.write("\t".join(TRACT_COLS) + "\n")
        for i, p in enumerate(paths, 1):
            t0 = time.time()
            op = gzip.open if p.endswith(".gz") else open
            with op(p, "rt") as fh:
                target, tarm, bins, tracts, st = scan(fh, a.bin, a.min_cov, a.min_alleles,
                                                      a.min_arm_copies, a.min_delta)
            for r in bins:
                fb.write(f"{target}\t{tarm}\t" + "\t".join(
                    (f"{r[c]:.5f}" if isinstance(r.get(c), float) else str(r.get(c, "NA"))) for c in BIN_COLS) + "\n")
            for t in tracts:
                ft.write("\t".join((f"{t[c]:.5f}" if isinstance(t[c], float) else str(t[c])) for c in TRACT_COLS) + "\n")
            tot.update(st)
            log(f"{i}/{len(paths)} {target} {tarm}: " + ", ".join(f"{k} {v}" for k, v in sorted(st.items()))
                + f"; {time.time() - t0:.1f} s")
    log("TOTAL " + ", ".join(f"{k} {v}" for k, v in sorted(tot.items())))
    if tot["records_kept"] == 0:
        raise SystemExit("FATAL: 0 records kept from non-empty input")


if __name__ == "__main__":
    main()
