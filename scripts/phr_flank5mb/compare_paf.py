#!/usr/bin/env python3
"""Compare two PAF files record by record (all fields except the chain id ch:Z, which numbers chains within one run).
Records whose query (column 1) or target (column 6) is in EXCLUDE are dropped from both before the comparison, so a
run on the 18,827 flanks can be compared with a run on the 18,799 flanks without the 28 removed ends. Prints one line
per pair: records in each file, dropped by EXCLUDE, identical, only in REF, only in NEW, identical multiset.
usage: compare_paf.py EXCLUDE.txt REF.paf.gz NEW.paf.gz [REF2 NEW2 ...]"""
import gzip, sys

excl = set(l.strip() for l in open(sys.argv[1]) if l.strip())
pairs = sys.argv[2:]
if not pairs or len(pairs) % 2:
    sys.exit("ERROR give REF NEW pairs")


def load(p):
    keep, drop = [], 0
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as f:
        for l in f:
            x = l.rstrip("\n").split("\t")
            if x[0] in excl or x[5] in excl:
                drop += 1
                continue
            keep.append("\t".join(v for v in x if not v.startswith("ch:Z:")))
    return sorted(keep), drop


bad = 0
for ref, new in zip(pairs[::2], pairs[1::2]):
    a, da = load(ref)
    b, db = load(new)
    sa, sb = set(a), set(b)
    same = a == b
    bad += not same
    print(f"{new}\treference {len(a)} kept ({da} dropped), new {len(b)} kept ({db} dropped); identical {len(sa & sb)}, "
          f"only reference {len(sa - sb)}, only new {len(sb - sa)}; identical multiset {same}", flush=True)
print(f"{len(pairs) // 2 - bad} of {len(pairs) // 2} pairs identical")
