#!/usr/bin/env python3
"""Depth calibration of the mapping artefact (referee R1-M4 / R3-M3): how many of the contacts observed between the PHR
windows of two chromosome ends can mapping errors alone explain?

In a null simulation (sim_hic_null.py, sim_porec_null.py) every read comes from one end, so every raw contact between
the PHR windows of two ends on different chromosomes (X_sim) is a placement error. The simulation and the real data
went through the same pipeline and are binned on the same chromosomes, so the artefact expected in the real data is
X_sim scaled by the real-to-simulated depth of the ends: s_A = C_real(A) / C_sim(A), where C is the raw contact count
inside the PHR window of end A with itself, excluding bin pairs closer than MIN_DIAG bins (the simulation has no
contacts below 1 kb). An artefact between A and B comes from reads of A or of B, so the expected count lies between
E_lo = X_sim * min(s_A, s_B) and E_hi = X_sim * max(s_A, s_B). The observed count O is the raw real contact count
between the two PHR windows. Raw counts at RESOLUTION (5 kb), no balancing, both matrices. An end with no simulated
cis contact in its PHR window (possible with MAPQ filters) has no scale: a pair with one such end uses the scale of the
other end for both bounds, a pair with two such ends gets E = 0 when X_sim = 0 and otherwise the smallest and the
largest scale over all ends as its bounds; the ends and pairs concerned are logged.

Outputs: OUT.pairs.tsv (per pair of ends on different chromosomes, homologs of the two haplotypes count as the same
chromosome: O, X_sim, s_A, s_B, E_lo, E_hi, bin pairs); OUT.arms.tsv (per pair of arms: sums over the ends, densities
per bin pair, arm mean Jaccard from the paper's PAIR_CORR table, 0 when the arm pair is absent); OUT.summary.tsv
(totals, pairs with E_hi >= O, Spearman of arm Jaccard vs observed density and vs the density left after subtracting
E_hi, floored at 0).
usage: calibrate_artefact.py DOMAINS_TSV SIM_MCOOL REAL_MCOOL PAIR_CORR_TSV OUT [--resolution 5000] [--min-diag 2]"""
import argparse, collections, re, sys
import cooler, numpy as np
from scipy.stats import spearmanr
ap = argparse.ArgumentParser()
for k in ("domains", "sim", "real", "paircorr", "out"): ap.add_argument(k)
ap.add_argument("--resolution", type=int, default=5000); ap.add_argument("--min-diag", type=int, default=2)
a = ap.parse_args()
log = lambda *m: print(*m, file=sys.stderr, flush=True)
D = [dict(zip(h, l.rstrip("\n").split("\t"))) for h in [open(a.domains).readline().rstrip("\n").split("\t")] for l in open(a.domains) if not l.startswith("end\t")]
hap = lambda x: re.sub(r"_(MATERNAL|PATERNAL)", "", x)
C = {k: cooler.Cooler(f"{p}::resolutions/{a.resolution}") for k, p in (("sim", a.sim), ("real", a.real))}
if dict(C["sim"].chromsizes) != dict(C["real"].chromsizes): sys.exit("FATAL simulated and real matrices have different chromosomes")
log(f"domains {a.domains}: {len(D)} ends; sim {a.sim}; real {a.real}; resolution {a.resolution}; min diag {a.min_diag}")
reg = {d["end"]: f"{d['chrom']}:{int(d['phr_start'])}-{int(d['phr_end'])}" for d in D}
def cis(clr, r):
    m = clr.matrix(balance=False, sparse=True).fetch(r).tocoo(); k = m.col - m.row >= a.min_diag
    return float(m.data[k].sum()), int(m.shape[0])
s, nb = {}, {}
for e, r in reg.items():
    (cs, nb[e]), (cr, _) = cis(C["sim"], r), cis(C["real"], r)
    s[e] = cr / cs if cs > 0 else float("nan")
nos = sorted(e for e in s if s[e] != s[e])
log(f"ends without simulated cis contacts in the PHR window (no scale): {len(nos)} {nos}")
v = [x for x in s.values() if x == x]
log("scale s = real / simulated cis contacts in the PHR window: median %.4g, range %.4g-%.4g" % (np.median(v), min(v), max(v)))
def bounds(xs, sa, sb):
    w = [x for x in (sa, sb) if x == x]
    if w: return xs * min(w), xs * max(w)
    if xs > 0: both[0] += 1; both[1] += xs; return xs * min(v), xs * max(v)
    return 0.0, 0.0
both = [0, 0.0]
half = collections.Counter()
tot = lambda clr, x, y: float(clr.matrix(balance=False, sparse=True).fetch(reg[x], reg[y]).sum())
ends = sorted(reg); rows = []
for i, x in enumerate(ends):
    for y in ends[i + 1:]:
        cx, cy = hap(x).rsplit("_", 1)[0], hap(y).rsplit("_", 1)[0]
        if cx == cy: continue
        o, xs = tot(C["real"], x, y), tot(C["sim"], x, y)
        half[sum(v != v for v in (s[x], s[y]))] += 1
        rows.append((x, y, o, xs, s[x], s[y], *bounds(xs, s[x], s[y]), nb[x] * nb[y]))
log(f"end pairs by number of ends without a scale: {dict(half)}; pairs with two such ends and simulated artefact: {both[0]}"
    f" ({both[1]:.0f} simulated contacts, bounded by the smallest and largest scale)")
with open(a.out + ".pairs.tsv", "w") as f:
    f.write("end_a\tend_b\tobserved\tsim_artefact\tscale_a\tscale_b\texpected_lo\texpected_hi\tbin_pairs\n")
    for r in rows: f.write("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r) + "\n")
jac = {}
with open(a.paircorr) as f:
    h = f.readline().rstrip("\n").split("\t")
    for l in f:
        r = dict(zip(h, l.rstrip("\n").split("\t"))); jac[frozenset((r["arm_a"], r["arm_b"]))] = float(r["mean_jaccard"])
arm = collections.defaultdict(lambda: np.zeros(4))
for x, y, o, xs, _, _, el, eh, bp in rows:
    k = tuple(sorted((hap(x), hap(y))))
    if k[0] != k[1]: arm[k] += (o, el, eh, bp)
with open(a.out + ".arms.tsv", "w") as f:
    f.write("arm_a\tarm_b\tmean_jaccard\tobserved\texpected_lo\texpected_hi\tbin_pairs\tobserved_density\texcess_density_hi\n")
    A = []
    for (p, q), (o, el, eh, bp) in sorted(arm.items()):
        j = jac.get(frozenset((p, q)), 0.0); A.append((j, o / bp, max(o - eh, 0) / bp))
        f.write(f"{p}\t{q}\t{j:.6g}\t{o:.6g}\t{el:.6g}\t{eh:.6g}\t{bp:.0f}\t{o / bp:.6g}\t{max(o - eh, 0) / bp:.6g}\n")
O, EL, EH = (sum(r[k] for r in rows) for k in (2, 6, 7))
A = np.array(A); top = A[:, 0] > 0
with open(a.out + ".summary.tsv", "w") as f:
    f.write("measure\tvalue\n")
    for k, v in (("end_pairs", len(rows)), ("observed_total", O), ("expected_lo_total", EL), ("expected_hi_total", EH),
                 ("expected_lo_over_observed", EL / O), ("expected_hi_over_observed", EH / O),
                 ("end_pairs_expected_hi_ge_observed", sum(r[7] >= r[2] for r in rows)),
                 ("arm_pairs", len(A)), ("arm_pairs_jaccard_gt0", int(top.sum())),
                 ("rho_jaccard_observed", spearmanr(A[:, 0], A[:, 1])[0]), ("rho_jaccard_excess_hi", spearmanr(A[:, 0], A[:, 2])[0]),
                 ("rho_jaccard_observed_jaccard_gt0", spearmanr(A[top, 0], A[top, 1])[0]), ("rho_jaccard_excess_hi_jaccard_gt0", spearmanr(A[top, 0], A[top, 2])[0])):
        f.write(f"{k}\t{v:.6g}\n" if isinstance(v, float) else f"{k}\t{v}\n")
log(open(a.out + ".summary.tsv").read())
