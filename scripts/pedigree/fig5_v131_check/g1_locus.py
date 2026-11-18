#!/usr/bin/env python3
"""Helpers for g1_locus.sbatch.

junction PAF WLEN EA EZ TERM_L TERM_R QLEN RULES: contig (+ strand) vs faulty window (asm5 --eqx, cg tags).
  Prints one line per rule: rule lt lq rt rq, or rule same (one anchor covers the whole error), or rule none.
  ex:L  exact runs >= L bp. Left = end of the run that starts before EA and ends last before EZ; right = start of the
        first run after the left junction that ends after EZ. (v2 asked the left run to end before EA: the 10 kb
        error bins are wider than the error, so a run reaching into the bin was dropped; G1r 2026-11-17.)
  id:T  alignment records with identity >= T (= / (= X I D)). Left = a record with >= 20 kb before EA, junction at
        min(EA, record end) projected to the contig; right = a record with >= 20 kb after EZ, junction at
        max(EZ, record start). One record over the whole error: the piece is the contig copy of [EA, EZ).
  TERM_L / TERM_R: window starts at the contig start / ends at the contig end; that side goes to the contig end.
votes DIP_SAM CAND_SAM CONTIG OFFSET CANDS OUT_TSV: reads on the whole diploid (shared index) and on the candidates
  (small index), same minimap2 preset, -N and -p. Per candidate (id lt rt nb): best AS on the assembly interval
  [OFFSET+lt, OFFSET+rt) of CONTIG, on the new part [lt, lt+nb) of cand id, and elsewhere in the diploid.
  OUT_TSV rows: cand, read, as_asm, as_cand, as_else, cov_asm, cov_cand (cov = covers its interval +- 1 kb), only
  reads that touch either interval. Prints per candidate: id votes_asm votes_cand span_asm span_cand, where a read
  votes for one side when its AS there beats both others (ties do not vote).
"""
import re, sys

def recs(paf):
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        if f[4] != "+": continue
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([=XID])", next(x[5:] for x in f[12:] if x.startswith("cg:Z:")))]
        yield int(f[7]), int(f[8]), int(f[2]), ops

def exact(rs, L):
    for t, t1, q, ops in rs:
        for n, op in ops:
            if op == "=" and n >= L: yield t, t + n, q, q + n
            if op in "=XD": t += n
            if op in "=XI": q += n

def proj(t, q, ops, x):
    for n, op in ops:
        if op in "=X" and t + n > x: return q + (x - t)
        if op == "D" and t + n > x: return q
        if op in "=XD": t += n
        if op in "=XI": q += n
    return q

def ident(ops):
    m = sum(n for n, op in ops if op == "="); return m / max(sum(n for n, op in ops), 1)

cmd = sys.argv[1]
if cmd == "junction":
    paf, wl, ea, ez, tl, tr, ql, rules = sys.argv[2], *map(int, sys.argv[3:9]), sys.argv[9].split(",")
    rs = list(recs(paf)); end = wl if tr else ez
    for rule in rules:
        kind, v = rule.split(":"); L = R = None
        if kind == "ex":
            r = list(exact(rs, int(v)))
            if any(t0 <= ea and t1 >= end for t0, t1, q0, q1 in r): print(rule, "same"); continue
            L = (0, 0) if tl else max(((t1, q1) for t0, t1, q0, q1 in r if t0 < ea and t1 < ez), default=None)
            if L: R = (wl, ql) if tr else min(((t0, q0) for t0, t1, q0, q1 in r if t0 >= L[0] and t1 > ez), default=None)
        elif kind == "id":
            g = [x for x in rs if ident(x[3]) >= float(v)]
            if tl: L = (0, 0)
            else:
                c = [x for x in g if x[0] <= ea - 20000]
                if c: x = max(c, key=lambda x: x[1]); j = min(ea, x[1]); L = (j, proj(x[0], x[2], x[3], j))
            if tr: R = (wl, ql)
            else:
                c = [x for x in g if x[1] >= ez + 20000]
                if c: x = min(c, key=lambda x: x[0]); j = max(ez, x[0]); R = (j, proj(x[0], x[2], x[3], j))
        else: sys.exit(f"ERROR unknown rule {rule}")
        if not L or not R or L[1] >= R[1] or L[0] >= R[0]: print(rule, "none")
        else: print(rule, L[0], L[1], R[0], R[1])
elif cmd == "votes":
    dsam, csam, ctg, off, cf, out = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6], sys.argv[7]
    cands = [(c, int(lt), int(rt), int(nb)) for c, lt, rt, nb in (l.split() for l in open(cf))]
    reads, ln = {}, {}
    for f in (dsam, csam):
        for l in open(f):
            x = l.split("\t")
            if l.startswith("@SQ"): ln[x[1][3:]] = int(x[2][3:])
            if l[0] == "@" or int(x[1]) & 0x4: continue
            p0 = int(x[3]) - 1; p1 = p0 + sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", x[5]) )
            reads.setdefault(x[0], []).append((x[2], p0, p1, next(int(y[5:]) for y in x[11:] if y.startswith("AS:i:"))))
    with open(out, "w") as o:
        for c, lt, rt, nb in cands:
            iv = {"asm": (ctg, off + lt, off + rt), "cand": (c, lt, lt + nb)}; v = [0, 0, 0, 0]
            for rn, al in reads.items():
                b = {"asm": -1, "cand": -1, "else": -1}; cov = {"asm": 0, "cand": 0}
                for t, p0, p1, a in al:
                    k = next((k for k, (cc, s0, s1) in iv.items() if t == cc and p0 < s1 and p1 > s0), None)
                    if k is None and t.startswith("cand_"): continue
                    k = k or "else"; b[k] = max(b[k], a)
                    if k != "else": cc, s0, s1 = iv[k]; cov[k] |= p0 <= max(s0 - 1000, 0) and p1 >= min(s1 + 1000, ln[cc])
                if b["asm"] < 0 and b["cand"] < 0: continue
                o.write(f"{c}\t{rn}\t{b['asm']}\t{b['cand']}\t{b['else']}\t{cov['asm']}\t{cov['cand']}\n")
                if b["asm"] > max(b["cand"], b["else"]): v[0] += 1; v[2] += cov["asm"]
                if b["cand"] > max(b["asm"], b["else"]): v[1] += 1; v[3] += cov["cand"]
            print(c, *v)
else:
    sys.exit(f"ERROR unknown command {cmd}")
