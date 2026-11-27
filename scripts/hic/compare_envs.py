#!/usr/bin/env python3
"""Compare two conda env trees file by file (a damaged original and its rebuild from the exact package list).
Skips conda-meta/, __pycache__/ and *.pyc. A file is: identical (same bytes); identical after prefix (same bytes once
the original env path is replaced by the rebuilt one, as conda writes the env path into text files and pads it with
NUL in binaries); different; missing in the rebuild; or only in the rebuild (expected for the files lost from the
original). Writes OUT.tsv (path, status) and prints counts per status, plus the first different files.
usage: compare_envs.py ORIGINAL REBUILT OUT"""
import hashlib, os, sys
a, b, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3]
def files(root):
    r = set()
    for d, dirs, fs in os.walk(root):
        dirs[:] = [x for x in dirs if x not in ("conda-meta", "__pycache__")]
        for f in fs:
            p = os.path.join(d, f)
            if not f.endswith(".pyc") and (os.path.isfile(p) or os.path.islink(p)): r.add(os.path.relpath(p, root))
    return r
def norm(data, prefix):
    if prefix.encode() not in data: return data
    # binaries: conda pads the replaced prefix with NUL bytes to the placeholder length, so compare up to the NUL run
    return b"\0".join(s.replace(prefix.encode(), b"<PREFIX>").rstrip(b"\0") for s in data.split(b"\0"))
fa, fb = files(a), files(b)
st, rows = {}, []
for p in sorted(fa | fb):
    if p not in fb: s = "missing_in_rebuild"
    elif p not in fa: s = "only_in_rebuild"
    else:
        pa, pb = os.path.join(a, p), os.path.join(b, p)
        if os.path.islink(pa) or os.path.islink(pb):
            s = "identical" if os.path.islink(pa) and os.path.islink(pb) and os.readlink(pa).replace(a, b) == os.readlink(pb) else "different"
        else:
            da, db = open(pa, "rb").read(), open(pb, "rb").read()
            s = "identical" if da == db else ("identical_after_prefix" if norm(da, a) == norm(db, b) else "different")
    st[s] = st.get(s, 0) + 1; rows.append((p, s))
with open(out, "w") as o:
    o.write("path\tstatus\n")
    for p, s in rows: o.write(f"{p}\t{s}\n")
print(f"{a} vs {b}: " + "; ".join(f"{k} {v}" for k, v in sorted(st.items())))
for k in ("different", "missing_in_rebuild", "only_in_rebuild"):
    v = [p for p, s in rows if s == k]
    print(f"{k} ({len(v)}): {v[:15]}")
