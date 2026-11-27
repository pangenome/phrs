#!/usr/bin/env python3
"""Recheck the files compare_envs.py called different: replace each env path by a placeholder and collapse NUL runs
(conda pads a replaced binary prefix with NULs, and the two paths differ in length); report what still differs.
usage: recheck_diff.py ORIGINAL REBUILT COMPARE_TSV"""
import collections, os, re, sys
a, b, t = sys.argv[1:4]
v = [l.split("\t")[0] for l in open(t) if l.rstrip("\n").endswith("\tdifferent")]
norm = lambda d, p: re.sub(rb"\0+", b"\0", d.replace(p.encode(), b"<P>"))
still = []
for p in v:
    da, db = open(os.path.join(a, p), "rb").read(), open(os.path.join(b, p), "rb").read()
    if norm(da, a) != norm(db, b): still.append((p, len(da), len(db)))
print(f"{len(v)} different; {len(v) - len(still)} identical after prefix and NUL padding; still different {len(still)}")
c = collections.Counter("/".join(p.split("/")[:2]) for p, _, _ in still); print(c.most_common(10))
for p, la, lb in still[:12]: print(f"  {p}: {la} vs {lb} bytes")
