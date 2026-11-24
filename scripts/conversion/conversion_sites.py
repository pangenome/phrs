#!/usr/bin/env python3
"""Site-level scan for conversion tracts between chromosome ends (R1-M1 test 3b, base resolution).

INPUT: target-centric PAFs of the PHR all-vs-all (wfmash -p 95, cg:Z =/X CIGAR), one file per window; query
names end in _chrN_parm or _chrN_qarm. For every query we count each target position once, from its record
with the most aligned target bases (then the next one on positions still free), as covered and as
matching (=) or not (X); positions under D are not covered.
PER POSITION p and arm: n = copies covering p, m = copies matching the target base at p.
  alleles = copies of the target's arm A (other haplotypes); Y = any other arm with >= --min-copies copies.
  A-SITE for Y: >= --min-copies alleles and Y copies, m_A >= --hi * n_A and m_Y <= --lo * n_Y
                (a paralogous sequence variant where the target carries its own arm's base)
  Y-SITE for Y: m_A <= --lo * n_A and m_Y >= --hi * n_Y
                (the target carries the base of arm Y where its own arm carries another one)
TRACT = a run of >= --min-sites Y-sites with no A-site between them; EMBEDDED when A-sites exist on both
sides. The tract spans its first to last Y-site; its bounds lie between the flanking A-sites and the tract.
A conversion shorter than one bin of conversion_scan.py is visible here if it covers >= 2 such sites.
OUTPUT: <out>.tracts.tsv (one row per tract) and <out>.power.tsv (per window and Y: A-sites and Y-sites,
the number of informative positions that make a tract detectable). Log on stderr.
"""
import argparse, gzip, re, sys, time
from collections import defaultdict
import numpy as np

CIGAR = re.compile(r"(\d+)([=XIDM])")
ARM = re.compile(r"_(chr[0-9XYM]+_(?:parm|qarm))$")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def arm_of(n):
    m = ARM.search(n)
    if not m:
        raise SystemExit(f"FATAL: no _chrN_parm/_qarm suffix in {n!r}")
    return m.group(1)


def scan(lines, min_copies=20, hi=0.9, lo=0.1, min_sites=2):
    stats = defaultdict(int)
    target = tlen = None
    recs = defaultdict(list)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            stats["dropped_short_line"] += 1
            continue
        q, t = f[0], f[5]
        if target is None:
            target, tlen = t, int(f[6])
        elif t != target:
            raise SystemExit(f"FATAL: second target {t!r} (first {target!r})")
        if q == t:
            stats["dropped_self"] += 1
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if cg is None:
            stats["dropped_no_cigar"] += 1
            continue
        if "M" in cg:
            raise SystemExit("FATAL: CIGAR has M; need =/X")
        recs[q].append((int(f[8]) - int(f[7]), int(f[7]), cg))
        stats["records_kept"] += 1
    if target is None:
        raise SystemExit("FATAL: no records")
    tarm = arm_of(target)
    cov, mat = {}, {}
    for q, rs in recs.items():
        a = arm_of(q)
        if a not in cov:
            cov[a], mat[a] = np.zeros(tlen, np.uint16), np.zeros(tlen, np.uint16)
        free = np.ones(tlen, bool)
        for _, pos, cg in sorted(rs, reverse=True):
            for n, op in CIGAR.findall(cg):
                n = int(n)
                if op in "=X":
                    f = free[pos:pos + n]
                    cov[a][pos:pos + n] += f
                    if op == "=":
                        mat[a][pos:pos + n] += f
                    free[pos:pos + n] = False
                    pos += n
                elif op == "D":
                    pos += n
        stats["queries"] += 1
    if tarm not in cov:
        return target, tarm, [], [], stats
    nA, mA = cov[tarm].astype(np.int32), mat[tarm].astype(np.int32)
    okA = nA >= min_copies
    a_hi, a_lo = okA & (mA >= hi * nA), okA & (mA <= lo * nA)
    tracts, power = [], []
    for y in sorted(cov):
        if y == tarm:
            continue
        nY, mY = cov[y].astype(np.int32), mat[y].astype(np.int32)
        okY = nY >= min_copies
        if not okY.any():
            continue
        asite = np.flatnonzero(a_hi & okY & (mY <= lo * nY))
        ysite = np.flatnonzero(a_lo & okY & (mY >= hi * nY))
        power.append((y, len(asite), len(ysite)))
        if len(ysite) < min_sites:
            continue
        ev = sorted([(p, 0) for p in asite] + [(p, 1) for p in ysite])
        i = 0
        while i < len(ev):
            if ev[i][1] == 0:
                i += 1
                continue
            j = i
            while j + 1 < len(ev) and ev[j + 1][1] == 1:
                j += 1
            if j - i + 1 >= min_sites:
                left = ev[i - 1][0] if i > 0 else None
                right = ev[j + 1][0] if j + 1 < len(ev) else None
                tracts.append(dict(target=target, arm=tarm, donor_arm=y, first_site=ev[i][0], last_site=ev[j][0],
                                   n_sites=j - i + 1, left_a_site=left, right_a_site=right,
                                   embedded=int(left is not None and right is not None)))
            i = j + 1
    stats["tracts"] = len(tracts)
    stats["tracts_embedded"] = sum(t["embedded"] for t in tracts)
    return target, tarm, tracts, power, stats


TCOLS = ["target", "arm", "donor_arm", "first_site", "last_site", "n_sites", "left_a_site", "right_a_site", "embedded"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pafs", nargs="*")
    ap.add_argument("--list", help="file with one PAF path per line")
    ap.add_argument("--chunk", help="N/K, 1-based")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-copies", type=int, default=20)
    ap.add_argument("--hi", type=float, default=0.9)
    ap.add_argument("--lo", type=float, default=0.1)
    ap.add_argument("--min-sites", type=int, default=2)
    a = ap.parse_args()
    paths = list(a.pafs) + ([l.strip() for l in open(a.list) if l.strip()] if a.list else [])
    if a.chunk:
        n, k = map(int, a.chunk.split("/"))
        paths = paths[n - 1::k]
    if not paths:
        raise SystemExit("FATAL: no input PAFs")
    log(f"inputs: {len(paths)} PAFs; min_copies {a.min_copies}, hi {a.hi}, lo {a.lo}, min_sites {a.min_sites}; out {a.out}")
    tot = defaultdict(int)
    with open(a.out + ".tracts.tsv", "w") as ft, open(a.out + ".power.tsv", "w") as fp:
        ft.write("\t".join(TCOLS) + "\n")
        fp.write("target\tarm\tdonor_arm\ta_sites\ty_sites\n")
        for i, p in enumerate(paths, 1):
            t0 = time.time()
            with (gzip.open if p.endswith(".gz") else open)(p, "rt") as fh:
                target, tarm, tracts, power, st = scan(fh, a.min_copies, a.hi, a.lo, a.min_sites)
            for t in tracts:
                ft.write("\t".join("NA" if t[c] is None else str(t[c]) for c in TCOLS) + "\n")
            for y, na, ny in power:
                fp.write(f"{target}\t{tarm}\t{y}\t{na}\t{ny}\n")
            for k, v in st.items():
                tot[k] += v
            log(f"{i}/{len(paths)} {target}: " + ", ".join(f"{k} {v}" for k, v in sorted(st.items()))
                + f"; {time.time() - t0:.1f} s")
    log("TOTAL " + ", ".join(f"{k} {v}" for k, v in sorted(tot.items())))
    if tot["records_kept"] == 0:
        raise SystemExit("FATAL: 0 records kept from non-empty input")


if __name__ == "__main__":
    main()
