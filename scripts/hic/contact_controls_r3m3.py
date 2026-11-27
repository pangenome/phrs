#!/usr/bin/env python3
"""R3-M3: does subtelomeric sequence similarity still predict contact after controlling for arm type, telomere
position, chromosome size and the general contact propensity of each end, with the dependence among pairs that
share an end kept in the null?

Inputs: the arm-level distance matrix (1 - mean inter-arm Jaccard; data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv) and
an arm-level contact matrix (community_based/<res>/<sample>_contact_matrix.tsv). Haplotype-resolved matrices
(chrN_MATERNAL_p ...) are collapsed to arms by the mean over the haplotype combinations of the two arms.
Pairs: the upper triangle over the arms present in both matrices (as the paper's Mantel test); --inter drops pairs
of the same chromosome.
Pair covariates:
  arm type       both p or both q (1) vs one of each (0); acrocentric p arms in the pair (0-2); X/Y arms (0-2)
  telomere pos.  log10 median PHR length of each arm across haplotypes (paper calls), summed over the pair
  chrom. size    log10 CHM13 length of each chromosome (cytoband ends): sum, and absolute difference
  propensity     log10 of each arm's total contact with all other arms (row sum without the diagonal), summed
Statistics (Spearman throughout; every p is from --perm permutations of the arm labels of the similarity matrix,
rows and columns together, which keeps the dependence among pairs that share an arm):
  mantel     rho(similarity, contact)
  partial    partial Mantel: ranks of similarity and contact each regressed on all covariates (OLS on ranks),
             correlation of the residuals (Smouse et al. 1986); in each permutation the permuted similarity is
             residualized again
  oe         rho(similarity, contact / expected), expected_ij = r_i r_j / sum(r) (degree-corrected propensity model)
  strata     mantel within same-type pairs, within p-q pairs, and without acrocentric p and X/Y arms
"""
import argparse, csv, re, statistics as st, sys
import numpy as np
from scipy.stats import rankdata

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dist", required=True)
ap.add_argument("--contact", required=True, nargs="+", help="LABEL=PATH")
ap.add_argument("--cytobands", required=True, help="CHM13 cytoband BED; chromosome length = last band end (chm13.chrom.sizes lacks chrY)")
ap.add_argument("--calls", required=True, help="paper PHR calls (all-vs-all.1Mb.p95.id95.len.tsv)")
ap.add_argument("--perm", type=int, default=10000)
ap.add_argument("--seed", type=int, default=20261125)
ap.add_argument("--inter", action="store_true", help="drop pairs of the same chromosome")
ap.add_argument("--out", required=True)
a = ap.parse_args()


def load(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    lab = rows[0][1:]
    M = np.array([[float(v) for v in r[1:]] for r in rows[1:]])
    if [r[0] for r in rows[1:]] != lab:
        raise SystemExit(f"FATAL: row and column labels differ in {path}")
    return lab, (M + M.T) / 2


def collapse(lab, M):
    """chrN_MATERNAL_p -> chrN_p, mean over haplotype combinations; plain labels are kept"""
    arm = [re.sub(r"_(MATERNAL|PATERNAL)_", "_", x) for x in lab]
    if arm == lab:
        return lab, M
    u = sorted(set(arm), key=arm.index)
    C = np.zeros((len(u), len(u)))
    for i, x in enumerate(u):
        for j, y in enumerate(u):
            ii = [k for k, z in enumerate(arm) if z == x]; jj = [k for k, z in enumerate(arm) if z == y]
            C[i, j] = np.mean([M[p, q] for p in ii for q in jj if p != q]) if i != j else 0.0
    return u, C


size = {}
for l in open(a.cytobands):
    f = l.split("\t"); size[f[0]] = max(size.get(f[0], 0), int(f[2]))
L = {}
for r in csv.DictReader(open(a.calls), delimiter="\t"):
    if r["region_start"] != ".":
        m = re.search(r"_(chr[0-9XY]+)_(p|q)arm$", r["seq"])
        L.setdefault(f"{m.group(1)}_{m.group(2)}", []).append(int(r["region_end"]) - int(r["region_start"]))
medL = {k: st.median(v) for k, v in L.items()}
dl, D = load(a.dist)
ACRO = {"chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p"}
rng = np.random.default_rng(a.seed)
out = open(a.out, "w")
out.write("contact\tstatistic\tsubset\tn_arms\tn_pairs\trho\tp_perm\n")


def spearman(x, y):
    return np.corrcoef(rankdata(x), rankdata(y))[0, 1]


def resid(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


for spec in a.contact:
    name, path = spec.split("=", 1)
    cl, C = collapse(*load(path))
    arms = [x for x in cl if x in dl and x in medL]
    drop = sorted(set(cl) - set(arms))
    S = 1 - D[np.ix_([dl.index(x) for x in arms], [dl.index(x) for x in arms])]
    C = C[np.ix_([cl.index(x) for x in arms], [cl.index(x) for x in arms])]
    n = len(arms)
    chrom = [x.rsplit("_", 1)[0] for x in arms]
    I, J = np.triu_indices(n, k=1)
    keep = np.array([not (a.inter and chrom[i] == chrom[j]) for i, j in zip(I, J)])
    I, J = I[keep], J[keep]
    r = C.sum(1) - np.diag(C)
    z = [arms[k] for k in np.where(r <= 0)[0]]
    base = spearman(S[I, J], C[I, J])  # the paper's Mantel, every arm kept
    nb = np.array([spearman(S[np.ix_(p, p)][I, J], C[I, J]) for p in (rng.permutation(n) for _ in range(a.perm))])
    pb = (1 + (nb >= base).sum()) / (1 + a.perm)
    out.write(f"{name}\tmantel\tall arms\t{n}\t{len(I)}\t{base:.3f}\t{pb:.2g}\n")
    print(f"{name}\tmantel\tall arms\tarms {n}\tpairs {len(I)}\trho {base:.3f}\tp {pb:.2g}")
    if z:  # the propensity covariate needs contact: those arms leave the controlled statistics, and we say so
        sys.stderr.write(f"{name}: arms with no contact to any other arm, dropped from the controlled statistics: {z}\n")
        k = [i for i in range(n) if arms[i] not in z]
        arms = [arms[i] for i in k]; S = S[np.ix_(k, k)]; C = C[np.ix_(k, k)]; n = len(arms)
        chrom = [x.rsplit("_", 1)[0] for x in arms]
        I, J = np.triu_indices(n, k=1)
        keep = np.array([not (a.inter and chrom[i] == chrom[j]) for i, j in zip(I, J)])
        I, J = I[keep], J[keep]
        r = C.sum(1) - np.diag(C)
    typ = np.array([x[-1] for x in arms]); acro = np.array([x in ACRO for x in arms]); sex = np.array([x[:4] in ("chrX", "chrY") for x in arms])
    ls = np.log10([size[c] for c in chrom]); lp = np.log10([medL[x] for x in arms]); lr = np.log10(r)
    X = np.column_stack([np.ones(len(I)), typ[I] == typ[J], acro[I].astype(int) + acro[J], sex[I].astype(int) + sex[J],
                         lp[I] + lp[J], ls[I] + ls[J], np.abs(ls[I] - ls[J]), lr[I] + lr[J]]).astype(float)
    E = np.outer(r, r) / r.sum()
    c, s, oe = C[I, J], S[I, J], (C / E)[I, J]
    rc = rankdata(c); res_c = resid(rc, X)
    def stats(Sm):
        sv = Sm[I, J]
        rs = rankdata(sv)
        return spearman(sv, c), np.corrcoef(resid(rs, X), res_c)[0, 1], spearman(sv, oe)
    obs = stats(S)
    null = np.array([stats(S[np.ix_(p, p)]) for p in (rng.permutation(n) for _ in range(a.perm))])
    sys.stderr.write(f"{name}: {n} arms (not in the similarity matrix or without PHR calls: {drop}), {len(I)} pairs\n")
    for k, lab in enumerate(("mantel", "partial", "oe")):
        p = (1 + (null[:, k] >= obs[k]).sum()) / (1 + a.perm)
        out.write(f"{name}\t{lab}\tall\t{n}\t{len(I)}\t{obs[k]:.3f}\t{p:.2g}\n")
        print(f"{name}\t{lab}\tall\tarms {n}\tpairs {len(I)}\trho {obs[k]:.3f}\tp {p:.2g}")
    for j, cov in enumerate(("same arm type", "acrocentric p arms", "X/Y arms", "log median PHR length", "log chromosome sizes",
                             "|log size difference|", "log contact propensity"), 1):
        print(f"   covariate {cov}: rho with contact {spearman(X[:, j], c):.3f}, with similarity {spearman(X[:, j], s):.3f}")
    for lab, sel in (("same-type pairs", typ[I] == typ[J]), ("p-q pairs", typ[I] != typ[J]),
                     ("without acrocentric p and X/Y", ~(acro[I] | acro[J] | sex[I] | sex[J]))):
        if sel.sum() < 10:
            continue
        o = spearman(s[sel], c[sel])
        nl = []
        for _ in range(a.perm):
            p_ = rng.permutation(n); sp = S[np.ix_(p_, p_)][I, J]
            nl.append(spearman(sp[sel], c[sel]))
        p = (1 + (np.array(nl) >= o).sum()) / (1 + a.perm)
        out.write(f"{name}\tmantel\t{lab}\t{n}\t{int(sel.sum())}\t{o:.3f}\t{p:.2g}\n")
        print(f"{name}\tmantel\t{lab}\tpairs {int(sel.sum())}\trho {o:.3f}\tp {p:.2g}")
out.close()
