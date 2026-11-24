#!/usr/bin/env python3
"""Trim extra sequence beyond the telomere at one q end (terminal rule; the read vote cannot reject extra tip sequence).

usage: tip_trim.py SAMTOOLS BGZIP NORM_FA CONTIG KEEP_TO BAM[,BAM...] IN_DIR OUT_DIR
Rule (2026-11-23, user yes), checked on NORM_FA (the assembly the BAMs are aligned to), keep CONTIG 1..KEEP_TO:
  1. the 2,000 bp before KEEP_TO are >= 90% TTAGGG (canonical telomere);
  2. the removed tail is < 50% TTAGGG (not canonical telomere);
  3. >= 20 primary ONT reads that cross a unique anchor 30 kb before the contig end end within KEEP_TO +- 1,500 bp.
Then IN_DIR/<hap fasta> is trimmed at the same place: the 2,000 bp before KEEP_TO must occur once in the last 100 kb of
CONTIG in IN_DIR, and the IN_DIR tail after it must equal the NORM_FA tail. Both haplotype files go to OUT_DIR (new).
A failed rule check exits non-zero before writing; a failed tail match exits non-zero (OUT_DIR then incomplete). Log on stderr.
"""
import os, re, subprocess, sys
S, BG, nf, C, K, bams, ind, outd = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6].split(","), sys.argv[7], sys.argv[8]
def log(m): print(m, file=sys.stderr, flush=True)
def fai(fa): return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fa + ".fai")}
def get(fa, r): return "".join(subprocess.run([S, "faidx", fa, r], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
L = fai(nf)[C]; tail = get(nf, f"{C}:{K + 1}-{L}"); keep2k = get(nf, f"{C}:{K - 1999}-{K}")
p1 = keep2k.count("TTAGGG") * 6 / 2000; p2 = tail.count("TTAGGG") * 6 / len(tail)
log(f"{C} in {nf}: length {L}, keep 1..{K}, remove {L - K} bp; TTAGGG before {p1:.3f} (>= 0.9), tail {p2:.3f} (< 0.5)")
A = L - 30000; ends = []
for b in bams:
    for l in subprocess.run([S, "view", "-F", "0x904", b, f"{C}:{A}-{A}"], capture_output=True, text=True, check=True).stdout.splitlines():
        f = l.split("\t"); r, q = int(f[3]), 0
        for k, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            k = int(k)
            if op in "M=X":
                if r + k > A: q += max(0, A - r); break
                r += k; q += k
            elif op in "DN": r += k
            elif op in "IS": q += k
        else: continue
        ends.append(A + len(f[9]) - q)
near = sum(abs(e - K) <= 1500 for e in ends)
log(f"reads crossing anchor {A}: {len(ends)}; ending within {K} +- 1500: {near} (>= 20); ending past {K + 1500}: {sum(e > K + 1500 for e in ends)}")
if not (p1 >= 0.9 and p2 < 0.5 and near >= 20): sys.exit("ERROR rule not met, nothing written")
os.makedirs(outd, exist_ok=False)
for fn in sorted(os.listdir(ind)):
    if not fn.endswith(".fa.gz"): continue
    src = f"{ind}/{fn}"; dst = f"{outd}/{fn}"; n = fai(src)
    with subprocess.Popen([BG, "-c"], stdin=subprocess.PIPE, stdout=open(dst, "wb"), text=True) as z:
        for c in n:
            s = get(src, c)
            if c == C:
                i = s.rfind(keep2k, max(0, len(s) - 100000))
                if i < 0 or s.count(keep2k, max(0, len(s) - 100000)) != 1: sys.exit(f"ERROR {fn}: anchor 2 kb found {s.count(keep2k, max(0, len(s) - 100000))} times in the last 100 kb")
                if s[i + 2000:] != tail: sys.exit(f"ERROR {fn}: tail differs from {nf} ({len(s) - i - 2000} vs {len(tail)} bp)")
                log(f"{fn}: {C} {len(s)} -> {i + 2000} bp (removed {len(s) - i - 2000})"); s = s[:i + 2000]
            z.stdin.write(f">{c}\n" + "\n".join(s[j:j + 80] for j in range(0, len(s), 80)) + "\n")
    subprocess.run([S, "faidx", dst], check=True); log(f"wrote {dst}: {len(n)} contigs")
