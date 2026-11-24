#!/usr/bin/env python3
"""Compare conversion-scan runs with different bin sizes on the same windows (R1-M1 test 3b).

Each --run is LABEL=DIR with DIR/chunk_*.{bins,tracts}.tsv from conversion_scan.py. Windows = the targets
present in the FIRST run; other runs are restricted to them (e.g. the full 2 kb run restricted to the 10%
subset). Per run, excluding C7, C14 and C15 (copies as similar as alleles there): windows, scored bp,
converted bp and their fraction, embedded tracts, tract length quartiles. Concordance: the share of the
reference run's converted bp that each run also calls converted, and the reverse.
"""
import argparse, csv, glob, os, re, sys, time
from collections import defaultdict

ARM = re.compile(r"_(chr[0-9XYM]+_(?:parm|qarm))$")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def load(d, keep, skip_arms):
    conv, scored, tracts = defaultdict(list), 0, []
    for f in sorted(glob.glob(os.path.join(d, "chunk_*.bins.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            t = r["target"]
            if (keep is not None and t not in keep) or ARM.search(t).group(1) in skip_arms:
                continue
            if r["scored"] == "1":
                scored += int(r["end"]) - int(r["start"])
                if r["converted"] == "1":
                    conv[t].append((int(r["start"]), int(r["end"])))
    for f in sorted(glob.glob(os.path.join(d, "chunk_*.tracts.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if (keep is None or r["target"] in keep) and r["arm"] not in skip_arms and r["embedded"] == "1":
                tracts.append(int(r["end"]) - int(r["start"]))
    return conv, scored, sorted(tracts)


def overlap(a, b):
    """bp of intervals a covered by intervals b, both lists per target"""
    n = 0
    for t, iv in a.items():
        jv = b.get(t, [])
        for s, e in iv:
            n += sum(max(0, min(e, y) - max(s, x)) for x, y in jv)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="LABEL=DIR; the first run defines the windows")
    ap.add_argument("--reference", required=True, help="label of the run used for concordance")
    ap.add_argument("--communities", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    comm = {r["ChromArm"].replace("_p", "_parm").replace("_q", "_qarm"): r["Community"]
            for r in csv.DictReader(open(a.communities), delimiter="\t")}
    skip = {x for x, c in comm.items() if c in ("C7", "C14", "C15")}
    runs = [x.split("=", 1) for x in a.run]
    keep = {r["target"] for f in glob.glob(os.path.join(runs[0][1], "chunk_*.bins.tsv"))
            for r in csv.DictReader(open(f), delimiter="\t")}
    log(f"{len(keep)} windows from {runs[0][0]}; {len(skip)} arms excluded ({sorted(skip)})")
    res = {}
    for lab, d in runs:
        res[lab] = load(d, keep, skip)
        log(f"{lab}: {res[lab][1]} scored bp, {sum(e - s for v in res[lab][0].values() for s, e in v)} converted bp, "
            f"{len(res[lab][2])} embedded tracts")
    ref = res[a.reference][0]
    rbp = sum(e - s for v in ref.values() for s, e in v)
    q = lambda v, p: v[int(p * (len(v) - 1))] if v else "NA"
    with open(a.out, "w") as fh:
        fh.write("run\twindows_with_calls\tscored_bp\tconverted_bp\tconverted_frac\tembedded_tracts\t"
                 "tract_len_q25\ttract_len_median\ttract_len_q75\tref_converted_bp_also_called\tcalled_bp_also_in_ref\n")
        for lab, _ in runs:
            conv, sc, tr = res[lab]
            cbp = sum(e - s for v in conv.values() for s, e in v)
            fh.write(f"{lab}\t{len(conv)}\t{sc}\t{cbp}\t{cbp / sc:.5f}\t{len(tr)}\t{q(tr, .25)}\t{q(tr, .5)}\t"
                     f"{q(tr, .75)}\t{overlap(ref, conv) / rbp:.4f}\t{(overlap(conv, ref) / cbp) if cbp else 0:.4f}\n")
    log(f"wrote {a.out}")


if __name__ == "__main__":
    main()
