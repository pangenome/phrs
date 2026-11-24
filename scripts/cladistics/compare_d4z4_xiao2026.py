#!/usr/bin/env python3
"""Compare DUX4/DUX4L copies per 4q/10q window (Liftoff annotation of each HPRC release 2 assembly,
dux4_copies_per_window.tsv) with full D4Z4 repeat units counted on spanning raw reads by Xiao et al. 2026
(Genome Res 36:827, Supplementary Tables S2 and S6, transcribed in d4z4_xiao2026_read_counts.tsv).
Read counts are unphased and labelled by array type (4q or 10q sequence); we compare within the same
label. With as many values on both sides we pair them in sorted order; otherwise each value of the shorter
side is paired with the closest unused value of the other (an allele without spanning reads, or a
haplotype without a 4q or 10q window, has no partner).
"""
import collections, csv, sys
D = "scripts/cladistics/"
ours = collections.defaultdict(list)
for l in open(D + "dux4_copies_per_window.tsv"):
    w, arm, n = l.rstrip("\n").split("\t")
    ours[(w.split("#")[0], arm)].append(int(n))
rows = []
for r in csv.DictReader(open(D + "d4z4_xiao2026_read_counts.tsv"), delimiter="\t"):
    x = sorted(int(v) for v in r["full_units"].split(","))
    o = sorted(ours.get((r["sample"], r["arm_type"]), []))
    if not o:
        sys.stderr.write(f"{r['sample']} {r['arm_type']}: no window in our set\n")
        continue
    if len(x) == len(o):
        pairs = list(zip(x, o))
    else:  # pair each value of the shorter side with the closest unused value of the other side
        few, many = (x, o) if len(x) < len(o) else (o, x)
        left, pairs = list(many), []
        for v in few:
            w = min(left, key=lambda u: abs(u - v))
            left.remove(w)
            pairs.append((v, w) if few is x else (w, v))
    rows += [(r["sample"], r["arm_type"], a, b) for a, b in pairs]
print("sample\tarm\tread_full_units\tannotation_copies\tdifference")
for s, arm, a, b in rows:
    print(f"{s}\t{arm}\t{a}\t{b}\t{b - a:+d}")
d = [b - a for *_, a, b in rows]
print(f"# {len(rows)} arrays in {len({r[0] for r in rows})} individuals; annotation - reads: "
      f"0 in {d.count(0)}, +1 in {d.count(1)}, within 1 unit {sum(abs(v) <= 1 for v in d)}, "
      f"within 2 {sum(abs(v) <= 2 for v in d)}, larger {sum(abs(v) > 2 for v in d)}")
