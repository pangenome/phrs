#!/usr/bin/env python3
"""Every sequence of FASTA A must be in FASTA B with the same bases (both faidx-indexed, bgzip or plain).
Prints counts of identical, different and missing sequences; exits 1 if any differs or is missing."""
import hashlib, subprocess, sys

S, A, B = sys.argv[1:4]
names = [l.split("\t")[0] for l in open(A + ".fai")]
inb = {l.split("\t")[0] for l in open(B + ".fai")}


def md5(fa, n):
    r = subprocess.run([S, "faidx", fa, n], capture_output=True, text=True, check=True)
    return hashlib.md5("".join(r.stdout.split("\n")[1:]).upper().encode()).hexdigest()


same = diff = miss = 0
for n in names:
    if n not in inb:
        miss += 1; print(f"missing in B: {n}"); continue
    if md5(A, n) == md5(B, n):
        same += 1
    else:
        diff += 1; print(f"differs: {n}")
print(f"{len(names)} sequences in A: identical {same}, different {diff}, missing in B {miss}")
sys.exit(1 if diff or miss else 0)
