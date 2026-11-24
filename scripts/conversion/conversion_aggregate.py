#!/usr/bin/env python3
"""Aggregate the conversion scan (conversion_scan.py chunk outputs) per community and arm (R1-M1, test 3b).

Per target window: scored and converted 2 kb bins, embedded tracts. Telomere distance of a tract:
p-arm windows start at the telomere (distance = tract start), q-arm windows end at it (distance =
window length - tract end). SHARED tract = an embedded tract whose (arm, donor arm) and telomere-distance
interval overlap (>= --min-overlap of the shorter one) an embedded tract in a window of ANOTHER individual;
a tract seen in one individual only is a SINGLETON (both haplotypes of one person count once), the class
in which assembly errors would fall. C7 (acrocentric short arms) and the pseudoautosomal communities C14
and C15 are reported apart: their copies are as similar as alleles, so a "closer paralog" is expected there.
LOCUS = tracts of the same arm and donor arm merged across individuals where their telomere-distance
intervals overlap or touch.
OUTPUT: <out>.community.tsv, <out>.arm.tsv, <out>.tracts.tsv (embedded tracts with sharing), <out>.loci.tsv.
"""
import argparse, csv, glob, os, re, sys, time
from collections import defaultdict

NAME = re.compile(r"^([^#]+)#([^#]+)#.*:(\d+)-(\d+)_(chr[0-9XYM]+_(?:parm|qarm))$")
SPECIAL = {"C7": "acrocentric", "C14": "PAR", "C15": "PAR"}


def log(msg):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def parse(name):
    m = NAME.match(name)
    if not m:
        raise SystemExit(f"FATAL: cannot parse window name {name!r}")
    s, h, a, b, arm = m.groups()
    return s, h, int(b) - int(a), arm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunks", required=True, help="directory with chunk_*.{bins,tracts}.tsv")
    ap.add_argument("--communities", required=True, help="arm-level Leiden k=15 assignments (ChromArm, Community)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.5)
    a = ap.parse_args()
    comm = {r["ChromArm"].replace("_p", "_parm").replace("_q", "_qarm"): r["Community"]
            for r in csv.DictReader(open(a.communities), delimiter="\t")}
    bfiles = sorted(glob.glob(os.path.join(a.chunks, "chunk_*.bins.tsv")))
    tfiles = sorted(glob.glob(os.path.join(a.chunks, "chunk_*.tracts.tsv")))
    log(f"inputs: {len(bfiles)} bins files, {len(tfiles)} tracts files in {a.chunks}; {len(comm)} arms with community")
    if not bfiles or len(bfiles) != len(tfiles):
        raise SystemExit("FATAL: missing or unpaired chunk files")
    win = {}  # target -> [scored, converted]
    for f in bfiles:
        for r in csv.DictReader(open(f), delimiter="\t"):
            w = win.setdefault(r["target"], [0, 0])
            w[0] += int(r["scored"])
            w[1] += int(r["converted"])
    tracts, n_all = [], 0
    for f in tfiles:
        for r in csv.DictReader(open(f), delimiter="\t"):
            n_all += 1
            if r["embedded"] != "1":
                continue
            s, h, L, arm = parse(r["target"])
            st, en = int(r["start"]), int(r["end"])
            d0, d1 = (st, en) if arm.endswith("parm") else (L - en, L - st)
            tracts.append(dict(r, sample=s, hap=h, d0=d0, d1=d1, length=en - st,
                               community=comm.get(arm, "none")))
    log(f"{len(win)} windows, {n_all} tracts, {len(tracts)} embedded")
    groups = defaultdict(list)
    for i, t in enumerate(tracts):
        groups[(t["arm"], t["paralog_arm"])].append(i)
    for key, idx in groups.items():
        idx.sort(key=lambda i: tracts[i]["d0"])
        for x in idx:
            tx, others = tracts[x], set()
            for y in idx:
                ty = tracts[y]
                if ty["sample"] == tx["sample"]:
                    continue
                ov = min(tx["d1"], ty["d1"]) - max(tx["d0"], ty["d0"])
                if ov >= a.min_overlap * min(tx["length"], ty["length"]):
                    others.add(ty["sample"])
            tx["n_other_individuals"] = len(others)
    percomm, perarm = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))
    for t, (sc, cv) in win.items():
        arm = parse(t)[3]
        for d in (percomm[comm.get(arm, "none")], perarm[arm]):
            d["windows"] += 1
            d["bins_scored"] += sc
            d["bins_converted"] += cv
    for t in tracts:
        for d in (percomm[t["community"]], perarm[t["arm"]]):
            d["tracts_embedded"] += 1
            d["tracts_shared"] += int(t["n_other_individuals"] > 0)
            d["tract_bp"] += t["length"]
    keys = ["windows", "bins_scored", "bins_converted", "tracts_embedded", "tracts_shared", "tract_bp"]
    for path, table, label in ((a.out + ".community.tsv", percomm, "community"), (a.out + ".arm.tsv", perarm, "arm")):
        with open(path, "w") as fh:
            fh.write(f"{label}\tregime\t" + "\t".join(keys) + "\ttracts_per_window\tfrac_bins_converted\n")
            for k in sorted(table):
                d = table[k]
                c = k if label == "community" else comm.get(k, "none")
                fh.write(f"{k}\t{SPECIAL.get(c, 'facultative')}\t" + "\t".join(str(d[x]) for x in keys) +
                         f"\t{d['tracts_embedded'] / max(1, d['windows']):.4f}"
                         f"\t{d['bins_converted'] / max(1, d['bins_scored']):.5f}\n")
    tk = ["target", "sample", "hap", "arm", "community", "paralog_arm", "start", "end", "length", "d0", "d1",
          "d_allele_median", "d_arm_median", "d_paralog_min", "n_other_individuals"]
    with open(a.out + ".tracts.tsv", "w") as fh:
        fh.write("\t".join(tk) + "\n")
        for t in tracts:
            fh.write("\t".join(str(t[k]) for k in tk) + "\n")
    # LOCI: per (arm, donor arm), tracts of all individuals merged where their telomere-distance intervals
    # overlap or touch; n_individuals = distinct samples (both haplotypes of a person count once)
    loci = []
    for (arm, par), idx in sorted(groups.items()):
        cur = None
        for i in sorted(idx, key=lambda i: (tracts[i]["d0"], tracts[i]["d1"])):
            t = tracts[i]
            if cur and t["d0"] <= cur["d1"]:
                cur["d1"] = max(cur["d1"], t["d1"]); cur["samples"].add(t["sample"]); cur["lens"].append(t["length"])
            else:
                cur = dict(community=t["community"], arm=arm, paralog_arm=par, d0=t["d0"], d1=t["d1"],
                           samples={t["sample"]}, lens=[t["length"]])
                loci.append(cur)
    with open(a.out + ".loci.tsv", "w") as fh:
        fh.write("community\tarm\tparalog_arm\td0\td1\tn_tracts\tn_individuals\tmedian_tract_length\n")
        for x in loci:
            fh.write(f"{x['community']}\t{x['arm']}\t{x['paralog_arm']}\t{x['d0']}\t{x['d1']}\t{len(x['lens'])}\t"
                     f"{len(x['samples'])}\t{sorted(x['lens'])[len(x['lens']) // 2]}\n")
    log(f"{len(loci)} loci, {sum(len(x['samples']) == 1 for x in loci)} in one individual")
    log(f"wrote {a.out}.community.tsv, {a.out}.arm.tsv, {a.out}.tracts.tsv, {a.out}.loci.tsv")


if __name__ == "__main__":
    main()
