#!/usr/bin/env python3
"""Candidate assembly errors from the HiFi scan.

usage: g1_candidates_hifi.py WASHU SAMPLE VERSION OUT_TSV
Rules (fixed 2026-11-16 before the HiFi scan finished):
  a 10 kb bin is flagged when HiFi depth <= 0.25 x the contig median, or HiFi SNV >= 3 with depth >= 3;
  depth = MAPQ>=0 depth (bins_hifi_q0). Change 2026-11-16 after the first scan: MAPQ>=20 depth was 0 in
  hap1 = hap2 stretches (diploid mapping, MAPQ 1-5), so 22% of bins were flagged. Thresholds unchanged.
  flagged bins <= 20 kb apart are merged; a group is "satellite" when > 50% of it is in cenSat (annot_v11).
Per group: bins, min depth, SNV sum, ONT SNV sum (check/wg_scan/bins_ont), merqury error k-mers (merqury/<v>/<s>).
"""
import bisect, collections, gzip, os, statistics as st, sys

D, SM, V, out = sys.argv[1:5]
def aname(c):
    s, h, ch = c.split("#")
    return f"{s}.{ch}." + ({"1": "maternal", "2": "paternal"}[h] if s == "PAN027" else f"haplotype{h}")
cs = collections.defaultdict(list)
for l in open(f"{D}/annot_v11/{SM}.polished.diploid.cenSat.bed"):
    t = l.split("\t")
    if len(t) >= 3 and t[1].isdigit(): cs[t[0]].append((int(t[1]), int(t[2])))
mq = collections.defaultdict(list)
for h in (1, 2):
    f = f"{D}/merqury/{V}/{SM}/{SM}_hap{h}_only.bed.gz"
    if os.path.exists(f):
        for l in gzip.open(f, "rt"): t = l.split("\t"); mq[t[0]].append(int(t[1]))
rows, nb, nf = [], 0, 0
for h in (1, 2):
    for l in open(f"{D}/norm_{V}/{SM}#{h}.fa.gz.fai"):
        c = l.split("\t")[0]; tag = c.replace("#", "_")
        hb = f"{D}/check/wg_scan/bins_hifi_q0/{SM}.{V}.{tag}.bins.tsv"; ob = f"{D}/check/wg_scan/bins_ont/{SM}.{V}.{tag}.bins.tsv"
        if not os.path.exists(hb): sys.exit(f"ERROR missing {hb}")
        b = [(int(r[1]), float(r[5]), int(r[3])) for r in (x.rstrip("\n").split("\t") for x in open(hb))]
        ont = {int(r[1]): int(r[3]) for r in (x.rstrip("\n").split("\t") for x in open(ob))} if os.path.exists(ob) else {}
        med = st.median(d for _, d, _ in b); nb += len(b)
        fl = [x for x in b if x[1] <= 0.25 * med or (x[2] >= 3 and x[1] >= 3)]; nf += len(fl)
        g = []
        for x in fl:
            if g and x[0] - g[-1][-1][0] <= 20000: g[-1].append(x)
            else: g.append([x])
        an = aname(c); p = sorted(mq.get(c, []))
        for grp in g:
            a, z = grp[0][0], grp[-1][0] + 10000
            sat = sum(max(0, min(z, e) - max(a, s)) for s, e in cs[an]) / (z - a)
            k = bisect.bisect_left(p, z) - bisect.bisect_left(p, a)
            rows.append((c, a, z, len(grp), min(x[1] for x in grp), med, sum(x[2] for x in grp),
                         sum(ont.get(x[0], 0) for x in grp), k, f"{sat:.2f}", "satellite" if sat > 0.5 else "candidate"))
with open(out, "w") as o:
    o.write("contig\tstart\tend\tbins\thifi_depth_min\thifi_depth_contig_median\thifi_snv\tont_snv\tmerqury_kmers\tcensat_frac\tclass\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
c = collections.Counter(r[10] for r in rows)
print(f"{SM} {V}: bins {nb}, flagged {nf}, groups {len(rows)} (candidate {c['candidate']}, satellite {c['satellite']}), "
      f"candidate bp {sum(r[2] - r[1] for r in rows if r[10] == 'candidate')}")
