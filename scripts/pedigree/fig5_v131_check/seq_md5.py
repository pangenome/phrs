#!/usr/bin/env python3
"""md5 of every contig sequence (line breaks ignored, case kept) of one person's diploid assembly.

usage: seq_md5.py DIR SAMPLE > TSV   (DIR/SAMPLE#{1,2}.fa.gz)
out: file, contig, length, md5
"""
import gzip, hashlib, sys
d, sm = sys.argv[1:3]
for h in (1, 2):
    fn = f"{sm}#{h}.fa.gz"; c = None
    def out():
        if c is not None: print(f"{fn}\t{c}\t{n}\t{m.hexdigest()}")
    for l in gzip.open(f"{d}/{fn}", "rt"):
        if l[0] == ">": out(); c = l[1:].split()[0]; n = 0; m = hashlib.md5(); continue
        s = l.rstrip("\n"); n += len(s); m.update(s.encode())
    out()
