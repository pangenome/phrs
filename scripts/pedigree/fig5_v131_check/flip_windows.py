#!/usr/bin/env python3
"""Why interchrom calls change between two class-winner tables (A = f32, B = f48)."""
import csv, gzip, statistics as st, sys
from collections import Counter

ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}

def load(p):
    w = {}
    with gzip.open(p, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = (r["chrom"], int(r["start"]))
            w.setdefault(k, {})[r["winner_class"]] = (float(r["estimated.identity"]), int(r["class_candidate_count"]), int(r["raw_candidate_count"]))
    return {k: d for k, d in w.items() if "same_chrom" in d and "interchrom" in d}

def call(d):
    return "inter" if d["interchrom"][0] > d["same_chrom"][0] else "same"

def med(x):
    return f"{st.median(x):.4f}" if x else "NA"

a, b = load(sys.argv[1]), load(sys.argv[2])
com = a.keys() & b.keys()
flip = [k for k in com if call(a[k]) == "same" and call(b[k]) == "inter"]
back = [k for k in com if call(a[k]) == "inter" and call(b[k]) == "same"]
new = [k for k in b.keys() - a.keys() if call(b[k]) == "inter"]
print(f"windows A {len(a)}, B {len(b)}, common {len(com)}; same->inter {len(flip)}, inter->same {len(back)}, inter windows only in B {len(new)}")
chrom = lambda k: k[0].split("#")[-1]
for tag, ks in (("same->inter", flip), ("inter windows only in B", new)):
    c = Counter("acro" if chrom(k) in ACRO else "other" for k in ks)
    t = Counter(chrom(k) for k in ks).most_common(8)
    print(f"{tag}: acro {c['acro']}, other {c['other']}; top chromosomes {t}")
    print(f"  start < 20 Mb: {sum(k[1] < 20_000_000 for k in ks)}, within 1 Mb of query end: n/a")
print("flipped windows, A -> B (median): same_id %s -> %s, inter_id %s -> %s" % (
    med([a[k]["same_chrom"][0] for k in flip]), med([b[k]["same_chrom"][0] for k in flip]),
    med([a[k]["interchrom"][0] for k in flip]), med([b[k]["interchrom"][0] for k in flip])))
print("flipped windows, candidates per class A -> B (median): same %s -> %s, inter %s -> %s" % (
    st.median([a[k]["same_chrom"][1] for k in flip]) if flip else "NA", st.median([b[k]["same_chrom"][1] for k in flip]) if flip else "NA",
    st.median([a[k]["interchrom"][1] for k in flip]) if flip else "NA", st.median([b[k]["interchrom"][1] for k in flip]) if flip else "NA"))
rest = [k for k in com if k not in set(flip)]
print("all other common windows, candidates per class A -> B (median): same %s -> %s, inter %s -> %s" % (
    st.median([a[k]["same_chrom"][1] for k in rest]), st.median([b[k]["same_chrom"][1] for k in rest]),
    st.median([a[k]["interchrom"][1] for k in rest]), st.median([b[k]["interchrom"][1] for k in rest])))
up = sum(b[k]["interchrom"][0] > a[k]["interchrom"][0] for k in flip)
dn = sum(b[k]["same_chrom"][0] < a[k]["same_chrom"][0] for k in flip)
print(f"flipped windows where inter_id rose {up}, same_id fell {dn} (of {len(flip)})")
