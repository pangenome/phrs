#!/usr/bin/env python3
"""Do the p/q labels of the telomere-anchored flanks match where the flanks map on CHM13?

Stage 2 of PIPELINE.md names a flank _parm when the telomere is at the contig start and _qarm when it is at the contig
end. A contig stored q-to-p would get the wrong arm label. Here each 500 kb flank is placed on CHM13 with the
assembly-vs-CHM13 mapping PAFs (wfmash -m, stage 1): mapped bp of the flank region to its own CHM13 chromosome, split
into the p side (target midpoint before the centromere start) and the q side (after the centromere end), and the
strand of those bp. A label AGREES when most own-chromosome bp fall on the labelled side.
Output: TSV per flank and a summary per label class and per contig strand on stderr/stdout."""
import argparse, collections, glob, os, re, sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--fai", required=True, help="fai of the 500 kb flank FASTA")
ap.add_argument("--paf-dir", required=True, help="assembly-vs-chm13 *.map.paf")
ap.add_argument("--annot", required=True, help="CHM13 annotation BED with Centromere rows (data/chm13-annotations.bed)")
ap.add_argument("--out", required=True)
a = ap.parse_args()
cen = {}
for l in open(a.annot):
    c, s, e, k = l.rstrip("\n").split("\t")[:4]
    if k == "Centromere":
        cen[c] = (int(s), int(e))
NAME = re.compile(r"^(.+):(\d+)-(\d+)_(chr[0-9XY]+)_(parm|qarm)$")
fl = collections.defaultdict(list)  # contig -> flanks
for l in open(a.fai):
    n = l.split("\t")[0]
    m = NAME.match(n)
    if not m:
        raise SystemExit(f"FATAL: cannot parse {n}")
    fl[m.group(1)].append((n, int(m.group(2)) - 1, int(m.group(3)), m.group(4), m.group(5)))
sys.stderr.write(f"{sum(map(len, fl.values()))} flanks on {len(fl)} contigs\n")
side = collections.defaultdict(lambda: collections.Counter())
files = sorted(glob.glob(os.path.join(a.paf_dir, "*.map.paf")))
seen = 0
for f in files:
    for l in open(f):
        x = l.split("\t", 9)
        if x[0] not in fl:
            continue
        seen += 1
        qs, qe, strand, t, ts, te = int(x[2]), int(x[3]), x[4], x[5].split("#")[-1], int(x[7]), int(x[8])
        for n, s, e, ch, arm in fl[x[0]]:
            ov = min(qe, e) - max(qs, s)
            if ov <= 0:
                continue
            if t != ch or ch not in cen:
                side[n]["other chromosome"] += ov; continue
            mid = (ts + te) / 2
            k = "p" if mid < cen[ch][0] else ("q" if mid > cen[ch][1] else "centromere")
            side[n][k] += ov; side[n][strand] += ov
sys.stderr.write(f"{len(files)} PAFs, {seen} records on flank contigs\n")
res = collections.Counter()
with open(a.out, "w") as fh:
    fh.write("flank\tchrom\tlabel\tp_bp\tq_bp\tother_chrom_bp\tplus_bp\tminus_bp\tverdict\n")
    for c, v in fl.items():
        for n, s, e, ch, arm in v:
            d = side.get(n, collections.Counter())
            p, q = d["p"], d["q"]
            if p + q == 0:
                ver = "unplaced"
            else:
                ver = "agrees" if (p > q) == (arm == "parm") else "DISAGREES"
            strand = "+" if d["+"] >= d["-"] else "-"
            res[(arm, ver, strand if p + q else ".")] += 1
            fh.write(f"{n}\t{ch}\t{arm}\t{p}\t{q}\t{d['other chromosome']}\t{d['+']}\t{d['-']}\t{ver}\n")
for k in sorted(res):
    print(f"label {k[0]}, {k[1]}, contig strand {k[2]}: {res[k]}")
print(f"total disagreeing: {sum(v for k, v in res.items() if k[1] == 'DISAGREES')} of {sum(res.values())}; unplaced {sum(v for k, v in res.items() if k[1] == 'unplaced')}")
