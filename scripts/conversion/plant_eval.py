#!/usr/bin/env python3
"""Score the planted-tract runs (plant_run.sbatch) for the bin scan and the site scan (R1-M1 test 3b).

For each plant (<name>.plants.tsv, new coordinates) and each scan output of the planted window: FOUND when a
tract (bin scan: <name>.plant.b<size>.tracts.tsv; site scan: <name>.plant.sites.tracts.tsv) overlaps the
planted interval with the plant's donor arm (bin scan: paralog_arm, site scan: donor_arm); we report the
start and end offsets of the best-overlapping tract and whether it is embedded. Calls away from every plant
are counted in the planted and in the control run (the unmodified window), so a plant that creates
spurious calls elsewhere shows up as a difference.
"""
import argparse, csv, glob, os, re, sys


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="plant_run.sbatch OUTDIR")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = []
    for pt in sorted(glob.glob(os.path.join(a.dir, "*.plants.tsv"))):
        name = pt[:-len(".plants.tsv")]
        plants = rows(pt)
        donor = re.search(r"_from_(chr\w+?_[pq]arm)_\d+$", name).group(1)
        runs = [(f"bin {os.path.basename(f).split('.b')[-1].split('.')[0]}", f, "paralog_arm", "start", "end")
                for f in sorted(glob.glob(name + ".plant.b*.tracts.tsv"))]
        if os.path.exists(name + ".plant.sites.tracts.tsv"):
            runs.append(("sites", name + ".plant.sites.tracts.tsv", "donor_arm", "first_site", "last_site"))
        for lab, f, dcol, scol, ecol in runs:
            tr = rows(f)
            ctrl = rows(f.replace(".plant.", ".orig."))
            iv = [(int(p["new_start"]), int(p["new_end"])) for p in plants]
            off = [t for t in tr if not any(int(t[scol]) < e and int(t[ecol]) + (ecol == "last_site") > s for s, e in iv)]
            for p in plants:
                s, e = int(p["new_start"]), int(p["new_end"])
                hit = [t for t in tr if t[dcol] == donor and int(t[scol]) < e and int(t[ecol]) + (ecol == "last_site") > s]
                b = max(hit, key=lambda t: min(e, int(t[ecol])) - max(s, int(t[scol])), default=None)
                out.append(dict(window=os.path.basename(name), donor=donor, scan=lab, length=p["length"],
                                subst=p["subst"], found=int(b is not None),
                                start_offset=int(b[scol]) - s if b else "NA", end_offset=int(b[ecol]) - e if b else "NA",
                                embedded=b["embedded"] if b else "NA", calls_off_plants=len(off),
                                control_calls=len(ctrl) if ctrl is not None else "NA"))
    cols = ["window", "donor", "scan", "length", "subst", "found", "start_offset", "end_offset", "embedded",
            "calls_off_plants", "control_calls"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    sys.stderr.write(f"{len(out)} plant x scan rows; wrote {a.out}\n")


if __name__ == "__main__":
    main()
