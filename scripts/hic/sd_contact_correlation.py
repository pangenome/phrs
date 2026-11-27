#!/usr/bin/env python3
"""Similarity vs contact for PHR loci and matched non-subtelomeric SD loci (referee R1-M5), with the paper's per-pair
code: extract_hic_contact and load_similarity_matrix imported from sequence_hic_correlation.py (balanced contact summed
over the region pair, divided by the number of valid bin pairs; odgi Jaccard, 0 when the pair is absent; loci shorter
than 5 kb skipped, pairs with no valid bin skipped, as --min-region-bp 5000 and n_bins == 0 in the paper script).
For each set (PHR loci, control loci from matched_sd_select.py) all pairs of loci on different chromosomes (homologs of
the two haplotypes count as the same chromosome) are scored, as in the paper's sequence-level analysis; Spearman rho
of Jaccard vs contact per set, also restricted to pairs with Jaccard > 0. Uncertainty: 1,000 bootstrap resamples of
loci (all pairs among the resampled loci), giving a 95% interval for each rho and for rho(PHR) - rho(control).
usage: sd_contact_correlation.py PAPER_SCRIPT_DIR PAIRS_TSV NAMES_TSV SIMILARITY_TSV_GZ MCOOL RESOLUTION OUT [--seed 20261126]"""
import argparse, collections, logging, random, sys
import numpy as np
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
for k in ("scripts", "pairs", "names", "sim", "mcool", "resolution", "out"): ap.add_argument(k)
ap.add_argument("--seed", type=int, default=20261126); ap.add_argument("--boot", type=int, default=1000)
a = ap.parse_args()
sys.path.insert(0, a.scripts)
from sequence_hic_correlation import extract_hic_contact, load_similarity_matrix  # the paper's functions
import cooler
logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s"); log = logging.getLogger()
base = lambda c: c.rsplit("_", 1)[0] if c.endswith(("_MATERNAL", "_PATERNAL")) else c
pansn = {l.split("\t")[1].strip(): l.split("\t")[0] for l in open(a.names)}  # "chr:start1-end" -> L<n>#0#...
loci = collections.defaultdict(dict)  # set -> {pansn: (chrom, start0, end)}
with open(a.pairs) as f:
    h = f.readline().rstrip("\n").split("\t")
    for l in f:
        r = dict(zip(h, l.rstrip("\n").split("\t"))); st = "PHR" if r["set"] == "PHR" else "control"
        for s in "AB":
            c, s0, e = r["chr" + s], int(r["start" + s]), int(r["end" + s])
            loci[st][pansn[f"{c}:{s0 + 1}-{e}"]] = (c, s0, e)
clr = cooler.Cooler(f"{a.mcool}::resolutions/{a.resolution}")
log.info(f"{a.mcool} at {a.resolution} bp; loci PHR {len(loci['PHR'])}, control {len(loci['control'])}")
keep = {st: {n: v for n, v in d.items() if v[2] - v[1] >= 5000 and v[0] in clr.chromnames} for st, d in loci.items()}
log.info(f"loci >= 5 kb on cooler chromosomes: PHR {len(keep['PHR'])}, control {len(keep['control'])}")
sim = load_similarity_matrix(a.sim, seq_set=set(keep["PHR"]) | set(keep["control"]))
rows = []
for st, d in keep.items():
    ns = sorted(d); nz = 0
    for i, x in enumerate(ns):
        for y in ns[i + 1:]:
            if base(d[x][0]) == base(d[y][0]): continue
            c, nb = extract_hic_contact(clr, f"{d[x][0]}:{d[x][1]}-{d[x][2]}", f"{d[y][0]}:{d[y][1]}-{d[y][2]}")
            if nb == 0: nz += 1; continue
            rows.append((st, x, y, sim.get((x, y), 0.0), c / nb, nb))
    log.info(f"{st}: {sum(r[0] == st for r in rows)} scored pairs; {nz} pairs with no valid bin skipped")
with open(a.out + ".pairs.tsv", "w") as o:
    o.write("set\tlocus_a\tlocus_b\tjaccard\tcontact_norm\tvalid_bins\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
def rho(rs, nonzero=False):
    v = [(j, c) for _, _, _, j, c, _ in rs if j > 0 or not nonzero]
    return spearmanr([j for j, _ in v], [c for _, c in v])[0] if len(v) > 2 else float("nan"), len(v)
by = {st: [r for r in rows if r[0] == st] for st in keep}
rng = random.Random(a.seed); boot = collections.defaultdict(list)
for _ in range(a.boot):
    res = {}
    for st, d in keep.items():
        pick = collections.Counter(rng.choice(sorted(d)) for _ in d)
        rs = [r for r in by[st] for _ in range(pick[r[1]] * pick[r[2]])]
        res[st] = rho(rs)[0]
    boot["PHR"].append(res["PHR"]); boot["control"].append(res["control"]); boot["diff"].append(res["PHR"] - res["control"])
ci = lambda v: tuple(np.nanpercentile(v, [2.5, 97.5]))
with open(a.out + ".summary.tsv", "w") as o:
    o.write("set\tpairs\trho\tci_low\tci_high\tpairs_jaccard_gt0\trho_jaccard_gt0\n")
    for st in ("PHR", "control"):
        r0, n0 = rho(by[st]); r1, n1 = rho(by[st], True); lo, hi = ci(boot[st])
        o.write(f"{st}\t{n0}\t{r0:.4f}\t{lo:.4f}\t{hi:.4f}\t{n1}\t{r1:.4f}\n")
    lo, hi = ci(boot["diff"]); d0 = rho(by["PHR"])[0] - rho(by["control"])[0]
    o.write(f"PHR_minus_control\tNA\t{d0:.4f}\t{lo:.4f}\t{hi:.4f}\tNA\tNA\n")
log.info("summary:\n" + open(a.out + ".summary.tsv").read())
