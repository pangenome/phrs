#!/usr/bin/env python3
"""Telomere layout at every contig end: is the pure telomere at the tip, or is there a non-telomeric tail beyond it?

usage: end_telo_profile.py SAMTOOLS NORM_DIR SAMPLE[,SAMPLE...] > TSV
Per contig end, the last (q) or first (p) 50 kb in 500 bp windows, from the tip inward; motif TTAGGG at q, CCCTAA at p;
pure = window with motif x 6 >= 90% of 500 bp. tail = bp from the tip to the first pure window (0 = telomere at the tip);
pure_bp = pure windows in the run that starts there (gaps of one window allowed); tail_motif_pct = motif % in the tail;
tail_AT_pct = A+T % in the tail.
out: sample, contig, end, tail_bp, pure_bp, tail_motif_pct, tail_AT_pct (NA when no pure window in 50 kb)
"""
import subprocess, sys
S, nd, sms = sys.argv[1], sys.argv[2], sys.argv[3].split(","); W, E = 500, 50000
print("sample\tcontig\tend\ttail_bp\tpure_bp\ttail_motif_pct\ttail_AT_pct")
for sm in sms:
    for h in (1, 2):
        fa = f"{nd}/{sm}#{h}.fa.gz"
        for l in open(fa + ".fai"):
            c, L = l.split("\t")[0], int(l.split("\t")[1]); n = min(E, L)
            for end, a, m in (("p", 1, "CCCTAA"), ("q", L - n + 1, "TTAGGG")):
                s = "".join(subprocess.run([S, "faidx", fa, f"{c}:{a}-{a + n - 1}"], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
                if end == "q": s = s[::-1]; m = m[::-1]  # walk from the tip inward on the reversed string
                w = [s[k:k + W] for k in range(0, len(s) - W + 1, W)]
                pure = [x.count(m) * 6 >= 0.9 * W for x in w]
                if True not in pure: print(f"{sm}\t{c}\t{end}\tNA\t0\tNA\tNA"); continue
                i = pure.index(True); j = i; gap = 0
                while j < len(w) and (pure[j] or (gap == 0 and j + 1 < len(w) and pure[j + 1])):
                    gap += not pure[j]; j += 1
                t = s[:i * W]
                tm = f"{100 * t.count(m) * 6 / len(t):.0f}" if t else "NA"
                at = f"{100 * (t.count('A') + t.count('T')) / len(t):.0f}" if t else "NA"
                print(f"{sm}\t{c}\t{end}\t{i * W}\t{sum(pure[i:j]) * W}\t{tm}\t{at}")
