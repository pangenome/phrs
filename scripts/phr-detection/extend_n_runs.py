#!/usr/bin/env python3
"""N runs that reach a flank edge (flank_joins_coords.py TSV) continue in the contig, so their far side is not in the flank.
For each such run we read the assembly around it (HPRC release-2 FASTA over HTTPS, paths from the assemblies index CSV,
samtools faidx) and extend it to its full length in the contig. Output: the same TSV with the full N run coordinates for
those rows (all other rows unchanged), and per extended row a log line with the old and new extent and the contig length."""
import argparse, csv, re, subprocess, sys, time

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--coords", required=True)
ap.add_argument("--index", required=True, help="assemblies_release2_v1.0.index.csv")
ap.add_argument("--flanks", help="only these flanks (one per line); default all")
ap.add_argument("--step", type=int, default=1000000)
ap.add_argument("--out", required=True)
a = ap.parse_args()
log = lambda m: (sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"), sys.stderr.flush())
S3 = "https://s3-us-west-2.amazonaws.com/human-pangenomics"
asm = {(r["sample_id"], r["haplotype"]): r["assembly"].replace("s3://human-pangenomics", S3) for r in csv.DictReader(open(a.index))}
keep = set(l.strip() for l in open(a.flanks)) if a.flanks else None
FL = re.compile(r":(\d+)-(\d+)_chr")
fai = {}


def ctg_len(url, ctg):
    if url not in fai:
        r = subprocess.run(["curl", "-sf", url + ".fai"], capture_output=True, text=True, check=True)
        fai[url] = {l.split("\t")[0]: int(l.split("\t")[1]) for l in r.stdout.splitlines()}
    return fai[url][ctg]


def seq(url, ctg, s, e):
    r = subprocess.run(["samtools", "faidx", url, f"{ctg}:{s}-{e}"], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"FATAL: samtools faidx {ctg}:{s}-{e}: {r.stderr[-300:]}")
    return "".join(r.stdout.split("\n")[1:]).upper()


rows = list(csv.DictReader(open(a.coords), delimiter="\t"))
n = 0
for r in rows:
    fa, fb = map(int, FL.search(r["flank"]).groups())
    s, e = int(r["start"]), int(r["end"])
    if r["kind"] != "N" or (s != fa and e != fb) or (keep is not None and r["flank"] not in keep):
        continue
    url = asm[r["sample"], r["hap"]]; L = ctg_len(url, r["contig"]); s0, e0 = s, e
    while s > 1:  # extend left while the base before is N
        lo = max(1, s - a.step); x = seq(url, r["contig"], lo, s - 1); m = re.search("[^N]N*$", x)
        if m:
            s = lo + m.start() + 1; break
        s = lo
    while e < L:
        hi = min(L, e + a.step); x = seq(url, r["contig"], e + 1, hi); m = re.search("[^N]", x)
        if m:
            e = e + 1 + m.start() - 1; break
        e = hi
    r["start"], r["end"] = str(s), str(e); n += 1
    log(f"{r['flank']}: N run {s0}-{e0} ({e0 - s0 + 1} bp in the flank) is {s}-{e} ({e - s + 1} bp) in the contig of {L} bp")
with open(a.out, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
log(f"{n} N runs extended, {len(rows)} rows written to {a.out}")
