#!/usr/bin/env python3
"""Crossovers from homolog alignments, vs published blocks.

usage: g2_crossover.py HOMOLOG_DIR ANNOT_DIR TAG PAIR WINDOW_BP DELTA OUT_PREFIX [MIN_WINDOWS]
Per query chromosome (pat: PAN027#2 vs PAN011, mat: PAN027#1 vs PAN010, kid: KID.m vs PAN027) and window: identity
(= / (= + X + I + D), >= 50% of the window aligned) to parent h1 and h2 from g2_homolog.sbatch minimap2 PAFs (primary).
A window is h1 or h2 when the identity gap is >= DELTA, else tie. Runs of h1/h2 windows (ties skipped) give segments;
a change of haplotype between segments is a crossover. Published blocks (annot_v11/{maternal,paternal}.blocks...) give the
same per window on the query, from the parent haplotype of each block (they cover only part of each chromosome); windows
covered by blocks of both parent haplotypes are dropped (2026-11-21); kid: only blocks on KID's maternal copy.
Match = the two crossover intervals overlap within 1 Mb; an unmatched crossover of ours is a conflict when published
blocks lie within 5 Mb on both sides, else untestable.
out: OUT_PREFIX.segments.tsv, OUT_PREFIX.crossovers.tsv; summary on stdout.
"""
import collections, glob, gzip, re, sys

hd, ad, tag, pair, W, dl, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
MINW = int(sys.argv[8]) if len(sys.argv) > 8 else 2  # segments shorter than MINW windows are dropped (ours only)
KIDM = set("chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX".split())
CH = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
def win_id(paf):
    acc = collections.defaultdict(lambda: [0, 0, 0]); L = 0
    for l in gzip.open(paf, "rt"):
        f = l.rstrip("\n").split("\t"); L = int(f[1]); p = int(f[2])
        ops = [(int(a), b) for a, b in re.findall(r"(\d+)([=XID])", next(x[5:] for x in f[12:] if x.startswith("cg:Z:")))]
        if f[4] == "-": ops.reverse()
        for m, op in ops:
            if op == "D": acc[p // W][1] += m; continue
            e = p + m
            while p < e:
                b = p // W; s = min(e, (b + 1) * W) - p; a = acc[b]; a[1] += s; a[2] += s
                if op == "=": a[0] += s
                p += s
    return {b: eq / cols for b, (eq, cols, qb) in acc.items() if cols and qb >= 0.5 * W}, L
def segs(calls):
    s = []
    for b, h in calls:
        if s and s[-1][2] == h: s[-1][1] = b
        else: s.append([b, b, h])
    return s
# published per-window haplotype on the query
qpat = {"pat": ("paternal", "PAN027", "paternal"), "mat": ("maternal", "PAN027", "maternal")}
pub = collections.defaultdict(lambda: collections.defaultdict(set))  # window -> parent haplotypes of the published blocks covering it
if pair in qpat:
    fn, qs, ql = qpat[pair]
    for l in open(f"{ad}/{fn}.blocks.formatted.filtered.tsv"):
        t = l.rstrip("\n").split("\t"); g, q = t[0], t[1]
        h = re.search(r"haplotype(\d)", g).group(1); c, rg = q.split(":"); ch = c.split(".")[1]; a, b = map(int, rg.split("-"))
        for k in range(a // W, b // W + 1): pub[ch][k].add("h" + h)
else:
    skip = collections.Counter()
    for fn, h in (("maternal", "1"), ("paternal", "2")):  # KID.m blocks from SELF.m (PAN027#1) or SELF.p (PAN027#2)
        for l in open(f"{ad}/{fn}.blocks.formatted.filtered.tsv"):
            t = l.rstrip("\n").split("\t")
            if len(t) < 3 or not t[2].startswith("PAN028"): continue
            c, rg = t[2].split(":"); ch = c.split(".")[1]; a, b = map(int, rg.split("-"))
            # only blocks on KID's maternal copy (the query); blocks on her other copy (from PARTNER) are not KID.m
            if re.search(r"haplotype(\d)", c).group(1) != ("2" if ch in KIDM else "1"): skip[fn] += 1; continue
            for k in range(a // W, b // W + 1): pub[ch][k].add("h" + h)
    print(f"published blocks on KID's other copy, skipped: {dict(skip)}", file=sys.stderr)
amb = {ch: sum(1 for v in d.values() if len(v) > 1) for ch, d in pub.items()}  # windows covered by blocks of both haplotypes: no call
print(f"published windows covered by both haplotypes (dropped): {sum(amb.values())} ({ {c: n for c, n in amb.items() if n} })", file=sys.stderr)
pub = {ch: {k: next(iter(v)) for k, v in d.items() if len(v) == 1} for ch, d in pub.items()}
P = {"pat": "PAN011", "mat": "PAN010", "kid": "PAN027"}[pair]
so, xo = open(f"{out}.segments.tsv", "w"), open(f"{out}.crossovers.tsv", "w")
so.write("chrom\tsource\tstart\tend\thap\twindows\n"); xo.write("chrom\tsource\tleft_end\tright_start\tfrom\tto\tstatus\n")
tot = collections.Counter()
for ch in CH:
    f = {h: glob.glob(f"{hd}/paf/{tag}.{pair}.{ch}.{P}_{h}_{ch}.mm2.paf.gz") for h in "12"}
    if not f["1"] or not f["2"]: continue
    i1, L = win_id(f["1"][0]); i2, _ = win_id(f["2"][0])
    calls = []
    for b in sorted(set(i1) & set(i2)):
        d = i1[b] - i2[b]
        if abs(d) >= dl: calls.append((b, "h1" if d > 0 else "h2"))
    ours = segs(calls); pubs = segs(sorted(pub[ch].items())) if pub.get(ch) else []
    xs = {}
    for name, sg in (("ours", ours), ("published", pubs)):
        sg = [x for x in sg if x[1] - x[0] + 1 >= (MINW if name == "ours" else 2)]
        merged = []
        for x in sg:
            if merged and merged[-1][2] == x[2]: merged[-1][1] = x[1]
            else: merged.append(list(x))
        for x in merged: so.write(f"{ch}\t{name}\t{x[0] * W}\t{(x[1] + 1) * W}\t{x[2]}\t{x[1] - x[0] + 1}\n")
        xs[name] = [((merged[i][1] + 1) * W, merged[i + 1][0] * W, merged[i][2], merged[i + 1][2]) for i in range(len(merged) - 1)]
    cov = [(x[0] * W, (x[1] + 1) * W) for x in pubs]
    for name, other in (("ours", "published"), ("published", "ours")):
        for a, b, h1, h2 in xs[name]:
            hit = any(a - 1e6 <= d and c <= b + 1e6 for c, d, _, _ in xs[other])
            if name == "ours" and not hit:  # testable only if published blocks lie on both sides within 5 Mb
                l = any(e <= a and a - e <= 5e6 for st, e in cov); r = any(st >= b and st - b <= 5e6 for st, e in cov)
                st_ = "conflict" if (l and r) else "untestable"
            else: st_ = "match" if hit else "missed"
            tot[(name, st_)] += 1
            xo.write(f"{ch}\t{name}\t{a}\t{b}\t{h1}\t{h2}\t{st_}\n")
    tot["chrom"] += 1
print(f"{tag} {pair} window {W} delta {dl} min {MINW}: chromosomes {tot['chrom']}; ours: match {tot[('ours', 'match')]}, "
      f"conflict {tot[('ours', 'conflict')]}, untestable {tot[('ours', 'untestable')]}; published: match {tot[('published', 'match')]}, missed {tot[('published', 'missed')]}")
