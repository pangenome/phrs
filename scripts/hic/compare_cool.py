#!/usr/bin/env python3
"""Compare two cooler matrices at one resolution pixel by pixel (bins, then bin1_id, bin2_id, count in chunks).
Prints the number of pixels in each, the pixels that differ and the total counts; exits 1 if they differ.
usage: compare_cool.py MCOOL_A MCOOL_B [--resolution 5000] [--chunk 20000000]"""
import argparse, sys
import cooler
ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--resolution", type=int, default=5000); ap.add_argument("--chunk", type=int, default=20000000); a = ap.parse_args()
A, B = (cooler.Cooler(f"{p}::resolutions/{a.resolution}") for p in (a.a, a.b))
if not A.bins()[["chrom", "start", "end"]][:].equals(B.bins()[["chrom", "start", "end"]][:]): sys.exit(f"FATAL different bins: {a.a} {a.b}")
na, nb = A.info["nnz"], B.info["nnz"]; diff = 0; ta = tb = 0
if na == nb:
    for lo in range(0, na, a.chunk):
        x, y = (c.pixels()[lo:lo + a.chunk][["bin1_id", "bin2_id", "count"]].values for c in (A, B))
        diff += int((x != y).any(axis=1).sum()); ta += int(x[:, 2].sum()); tb += int(y[:, 2].sum())
print(f"{a.a} vs {a.b} at {a.resolution} bp: pixels {na} vs {nb}; differing pixels {diff if na == nb else 'NA (different nnz)'}; total counts {ta} vs {tb}")
sys.exit(0 if na == nb and diff == 0 else 1)
