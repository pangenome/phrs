#!/usr/bin/env python3
"""Helpers for g1_locus.sbatch.

junction PAF WLEN EA EZ TERM_L TERM_R QLEN: contig (+ strand) vs faulty window with cg tags; prints lt lq rt rq.
  Left: end of the last exact run >= 20 kb ending <= EA (window start if TERM_L). Right: start of the first exact run
  >= 20 kb starting >= EZ (window end and contig end if TERM_R). Exit 4 when one exact run covers the whole error
  (contig = assembly), exit 3 when a side has no anchor.
compete SAM CONTIG OFFSET CANDS: reads on the whole diploid + candidates (minimap2 -N 20 -p 0.5). Per mode in CANDS
  (mode lt rt nb): best alignment score (AS) on the assembly interval, on the new part of cand_<mode>, and elsewhere.
  A read votes for the candidate when its AS there beats both others, for the assembly likewise; span = a voting read
  covering its interval +- 1 kb. Reads with the same best AS in two places (paralogs, flanks) do not vote.
  Prints per mode: mode votes_asm votes_cand span_asm span_cand.
"""
import re, sys

def runs(paf):
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        if f[4] != "+": continue
        q, t = int(f[2]), int(f[7])
        for n, op in re.findall(r"(\d+)([=XID])", next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))):
            n = int(n)
            if op == "=" and n >= 20000: yield t, t + n, q, q + n
            if op in "=XD": t += n
            if op in "=XI": q += n

cmd = sys.argv[1]
if cmd == "junction":
    paf, wl, ea, ez, tl, tr, ql = sys.argv[2], *map(int, sys.argv[3:9])
    r = list(runs(paf))
    if any(t0 <= ea and t1 >= (wl if tr else ez) for t0, t1, q0, q1 in r): sys.exit(4)
    L = (0, 0) if tl else max(((t1, q1) for t0, t1, q0, q1 in r if t1 <= ea), default=None)
    R = (wl, ql) if tr else min(((t0, q0) for t0, t1, q0, q1 in r if t0 >= ez), default=None)
    if not L or not R or L[1] >= R[1]: sys.exit(3)
    print(L[0], L[1], R[0], R[1])
elif cmd == "compete":
    sam, ctg, off, cf = sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    cands = [(m, int(lt), int(rt), int(nb)) for m, lt, rt, nb in (l.split() for l in open(cf))]
    reads, ln = {}, {}
    for l in open(sam):
        f = l.split("\t")
        if l.startswith("@SQ"): ln[f[1][3:]] = int(f[2][3:])
        if l[0] == "@" or int(f[1]) & 0x4: continue
        p0 = int(f[3]) - 1; p1 = p0 + sum(int(x) for x, op in re.findall(r"(\d+)([MDN=X])", f[5]))
        a = next(int(x[5:]) for x in f[11:] if x.startswith("AS:i:"))
        reads.setdefault(f[0], []).append((f[2], p0, p1, a))
    for m, lt, rt, nb in cands:
        iv = {"asm": (ctg, off + lt, off + rt), "cand": (f"cand_{m}", lt, lt + nb)}
        v = {"asm": [0, 0], "cand": [0, 0]}
        for al in reads.values():
            best = {"asm": -1, "cand": -1, "else": -1}; cov = {"asm": False, "cand": False}
            for t, p0, p1, a in al:
                k = next((k for k, (c, s0, s1) in iv.items() if t == c and p0 < s1 and p1 > s0), None)
                if k is None and t.startswith("cand_"): continue
                k = k or "else"; best[k] = max(best[k], a)
                if k != "else":
                    c, s0, s1 = iv[k]; cov[k] |= p0 <= max(s0 - 1000, 0) and p1 >= min(s1 + 1000, ln[c])
            for k, o in (("asm", "cand"), ("cand", "asm")):
                if best[k] > max(best[o], best["else"]): v[k][0] += 1; v[k][1] += cov[k]
        print(m, v["asm"][0], v["cand"][0], v["asm"][1], v["cand"][1])
else:
    sys.exit(f"ERROR unknown command {cmd}")
