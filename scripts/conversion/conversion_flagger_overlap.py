#!/usr/bin/env python3
"""Overlap of conversion-scan tracts with assembly-QC flags (R1-M1 test 3b control).

Flags: HPRC release-2 HMM-Flagger v1.2.1 HiFi predictions, conservative, canonical, no_Hap (one BED per
sample, PanSN contig names; every interval in it is a non-haploid call: Err, Dup, Col, ...).
Window names are contig:a-b_<arm> with a 1-based, so window offset x is contig position a - 1 + x.
Background = every scored bin of the scan (conversion_scan.py chunk_*.bins.tsv). A tract or bin is
FLAGGED when >= 1 bp of it overlaps a flag. Per community and in total: fraction of scored bins flagged,
fraction of tracts flagged, split into singleton tracts (n_other_individuals = 0) and shared tracts.
"""
import argparse, bisect, csv, glob, os, re, sys, time
from collections import defaultdict

NAME = re.compile(r"^(.+):(\d+)-(\d+)_(chr[0-9XYM]+_(?:parm|qarm))$")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def win(name):
    m = NAME.match(name)
    if not m:
        raise SystemExit(f"FATAL: cannot parse window name {name!r}")
    return m.group(1), int(m.group(2)) - 1, m.group(4)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flagger-dir", required=True, help="<sample>.bed files")
    ap.add_argument("--chunks", required=True, help="directory with chunk_*.bins.tsv")
    ap.add_argument("--tracts", required=True, help="conversion_aggregate.py <out>.tracts.tsv")
    ap.add_argument("--communities", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    comm = {r["ChromArm"].replace("_p", "_parm").replace("_q", "_qarm"): r["Community"]
            for r in csv.DictReader(open(a.communities), delimiter="\t")}
    flags, labels = defaultdict(list), defaultdict(int)
    beds = sorted(glob.glob(os.path.join(a.flagger_dir, "*.bed")))
    for f in beds:
        for line in open(f):
            if line.startswith("track"):
                continue
            c, s, e, lab = line.split("\t")[:4]
            if lab == "Hap":
                raise SystemExit(f"FATAL: Hap interval in a no_Hap BED {f}")
            flags[c].append((int(s), int(e)))
            labels[lab] += 1
    for c in flags:
        flags[c].sort()
    log(f"{len(beds)} Flagger BEDs, {len(flags)} contigs with flags, labels {dict(labels)}")
    samples = {os.path.basename(f)[:-4] for f in beds}

    def flagged(contig, s, e):
        v = flags.get(contig)
        if not v:
            return False
        i = bisect.bisect_left(v, (e, -1))
        return any(x < e and y > s for x, y in v[max(0, i - 50):i])

    bg = defaultdict(lambda: [0, 0])
    nskip = 0
    for f in sorted(glob.glob(os.path.join(a.chunks, "chunk_*.bins.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["scored"] != "1":
                continue
            contig, off, arm = win(r["target"])
            if contig.split("#")[0] not in samples:
                nskip += 1
                continue
            k = comm.get(arm, "none")
            bg[k][0] += 1
            bg[k][1] += flagged(contig, off + int(r["start"]), off + int(r["end"]))
    log(f"background: {sum(v[0] for v in bg.values())} scored bins; {nskip} skipped (sample without Flagger BED)")
    tr = defaultdict(lambda: [0, 0, 0, 0])  # singleton n, flagged, shared n, flagged
    tskip = 0
    for r in csv.DictReader(open(a.tracts), delimiter="\t"):
        contig, off, arm = win(r["target"])
        if contig.split("#")[0] not in samples:
            tskip += 1
            continue
        j = 0 if r["n_other_individuals"] == "0" else 2
        tr[r["community"]][j] += 1
        tr[r["community"]][j + 1] += flagged(contig, off + int(r["start"]), off + int(r["end"]))
    log(f"tracts: {sum(v[0] + v[2] for v in tr.values())}; {tskip} skipped (sample without Flagger BED)")
    with open(a.out, "w") as fh:
        fh.write("community\tbins_scored\tbins_flagged_frac\tsingleton_tracts\tsingleton_flagged_frac\t"
                 "shared_tracts\tshared_flagged_frac\n")
        keys = sorted(set(bg) | set(tr), key=lambda x: (x == "none", x))
        tot = [0, 0, 0, 0, 0, 0]
        for k in keys + ["ALL"]:
            if k == "ALL":
                b, t = tot[:2], tot[2:]
            else:
                b, t = bg[k], tr[k]
                for i, x in enumerate(b + t):
                    tot[i] += x
            fr = lambda n, d: f"{n / d:.4f}" if d else "NA"
            fh.write(f"{k}\t{b[0]}\t{fr(b[1], b[0])}\t{t[0]}\t{fr(t[1], t[0])}\t{t[2]}\t{fr(t[3], t[2])}\n")
    log(f"wrote {a.out}")


if __name__ == "__main__":
    main()
