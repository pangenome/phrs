#!/usr/bin/env python3
"""Clustering null for the site-level conversion scan (conversion_sites.py, R1-M1 test 3b).

For one target window and one donor arm Y the scan orders the informative positions along the target: A-sites (the target
carries its own arm's base where Y differs) and Y-sites (the target carries Y's base where its own arm differs). A tract is
a maximal run of >= 2 consecutive Y-sites. If Y-sites were independent of each other (sequencing or assembly errors,
recurrent mutation), their order among the A-sites would be random. Under that null, for nA A-sites and nY Y-sites
(n = nA + nY), the expected numbers are exact:
  runs of exactly L Y-sites  = (n-L-1) nA (nA-1) (nY)_L / (n)_{L+2}  +  2 nA (nY)_L / (n)_{L+1}   [L < n; (x)_k falling]
  embedded runs (A on both sides) = the first term only;  Y-Y adjacencies = nY (nY-1) / n.
We compare them with the observed runs (tracts.tsv) per (target, donor arm), summed over all pairs, by community and by nY.
usage: sites_cluster_null.py OUT_DIR K15_TSV OUT_PREFIX [--events EVENTS_TSV] [--selftest]
  OUT_DIR: chunk_*.power.tsv and chunk_*.tracts.tsv of conversion_sites.py. --selftest: check the formulas by shuffling.
Per run: significant when the expected number of runs at least as long, in its pair, is <= 0.05; the sum of those expected
  numbers over all pairs bounds the false significant runs. EVENTS_TSV (conversion_sites_aggregate.py events): an event is
  kept when one of its donor runs is significant.
out: OUT_PREFIX.summary.tsv (observed vs expected), OUT_PREFIX.runs.tsv, OUT_PREFIX.events_significant.tsv; log on stderr.
"""
import collections, csv, glob, math, random, sys, time
def log(m): print(f"[{time.strftime('%F %T')}] {m}", file=sys.stderr, flush=True)
LMAX = 12  # exact terms for L = 2..LMAX-1; ">= LMAX" = total minus those

def lff(x, k):  # log falling factorial (x)_k, -inf if k > x
    return -math.inf if k > x else math.lgamma(x + 1) - math.lgamma(x - k + 1)

def expected(nA, nY):
    n = nA + nY; e = [0.0] * (LMAX + 1); emb = 0.0
    if nY < 2: return e, emb, 0.0
    if nA == 0: e[min(nY, LMAX)] = 1.0; return e, emb, nY - 1.0
    tot = nY * (nY - 1) * (nA + 1) / (n * (n - 1))  # all runs >= 2
    for L in range(2, min(nY, LMAX - 1) + 1):
        a = (n - L - 1) * math.exp(math.log(nA) + (math.log(nA - 1) if nA > 1 else -math.inf) + lff(nY, L) - lff(n, L + 2)) if nA > 1 and n - L - 1 > 0 else 0.0
        b = 2 * math.exp(math.log(nA) + lff(nY, L) - lff(n, L + 1))
        e[L] = a + b; emb += a
    e[LMAX] = max(tot - sum(e[2:LMAX]), 0.0)
    if nY >= LMAX:  # embedded runs >= LMAX
        for L in range(LMAX, nY + 1):
            if n - L - 1 <= 0 or nA < 2: break
            t = (n - L - 1) * math.exp(math.log(nA) + math.log(nA - 1) + lff(nY, L) - lff(n, L + 2))
            emb += t
            if t < 1e-12 * max(emb, 1e-300): break
    return e, emb, nY * (nY - 1) / n

def observed_runs(seq):
    runs, i = [], 0
    while i < len(seq):
        if seq[i]:
            j = i
            while j + 1 < len(seq) and seq[j + 1]: j += 1
            if j - i + 1 >= 2: runs.append((j - i + 1, i > 0 and j + 1 < len(seq)))
            i = j + 1
        else: i += 1
    return runs

if "--selftest" in sys.argv:
    random.seed(1)
    for nA, nY, it in [(300, 4, 200000), (50, 10, 100000), (5, 8, 100000), (100, 30, 50000), (2, 20, 50000)]:
        s = [1] * nY + [0] * nA; c = collections.Counter(); emb = adj = 0
        for _ in range(it):
            random.shuffle(s); r = observed_runs(s)
            for L, e in r: c[min(L, LMAX)] += 1; emb += e
            adj += sum(L - 1 for L, _ in r)
        e, ee, ea = expected(nA, nY)
        log(f"selftest nA={nA} nY={nY}: runs>=2 sim {sum(c.values()) / it:.5f} exact {sum(e[2:]):.5f}; L=2 sim {c[2] / it:.5f} exact {e[2]:.5f}; "
            f"embedded sim {emb / it:.5f} exact {ee:.5f}; adjacencies sim {adj / it:.5f} exact {ea:.5f}")
    sys.exit(0)

D, K15, OUT = sys.argv[1:4]
com = {r["ChromArm"].replace("_p", "_parm").replace("_q", "_qarm"): r["Community"] for r in csv.DictReader(open(K15), delimiter="\t")}
pw = {}
for f in sorted(glob.glob(f"{D}/chunk_*.power.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        k = (r["target"], r["donor_arm"])
        if k in pw: sys.exit(f"ERROR duplicate power row {k}")
        pw[k] = (r["arm"], int(r["a_sites"]), int(r["y_sites"]))
obs = collections.defaultdict(list)
for f in sorted(glob.glob(f"{D}/chunk_*.tracts.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        obs[(r["target"], r["donor_arm"])].append((int(r["n_sites"]), r["embedded"] == "1"))
log(f"power rows {len(pw)} from {len(glob.glob(f'{D}/chunk_*.power.tsv'))} chunks; tract rows {sum(map(len, obs.values()))} in {len(obs)} pairs")
miss = [k for k in obs if k not in pw]
if miss: sys.exit(f"ERROR {len(miss)} tract pairs without a power row, e.g. {miss[:2]}")
for k, v in obs.items():
    if sum(L for L, _ in v) > pw[k][2]: sys.exit(f"ERROR {k}: tract sites {sum(L for L, _ in v)} > y_sites {pw[k][2]}")
def nyclass(ny): return "nY 2-4" if ny <= 4 else "nY 5-9" if ny <= 9 else "nY 10-49" if ny <= 49 else "nY >= 50"
keys = ["all"]; agg = collections.defaultdict(lambda: collections.defaultdict(float))
for k, (arm, nA, nY) in pw.items():
    e, emb, adj = expected(nA, nY); o = obs.get(k, [])
    for g in ("all", f"community {com.get(arm, 'none')}", nyclass(nY) if nY >= 2 else "nY 0-1"):
        a = agg[g]; a["pairs"] += 1; a["y_sites"] += nY; a["a_sites"] += nA
        a["exp_runs2"] += sum(e[2:]); a["exp_runs3"] += sum(e[3:]); a["exp_runs5"] += sum(e[5:]); a["exp_emb2"] += emb; a["exp_adj"] += adj
        a["obs_runs2"] += len(o); a["obs_runs3"] += sum(L >= 3 for L, _ in o); a["obs_runs5"] += sum(L >= 5 for L, _ in o)
        a["obs_emb2"] += sum(x for _, x in o); a["obs_adj"] += sum(L - 1 for L, _ in o); a["obs_sites_in_runs"] += sum(L for L, _ in o)
cols = ["group", "pairs", "a_sites", "y_sites", "obs_runs2", "exp_runs2", "obs_emb2", "exp_emb2", "obs_runs3", "exp_runs3",
        "obs_runs5", "exp_runs5", "obs_adj", "exp_adj", "obs_sites_in_runs"]
order = ["all"] + sorted(g for g in agg if g.startswith("nY")) + sorted((g for g in agg if g.startswith("community")), key=lambda g: (len(g), g))
with open(f"{OUT}.summary.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for g in order:
        a = agg[g]; fo.write("\t".join([g] + [f"{a[c]:.0f}" if c.startswith(("obs", "pairs", "a_", "y_")) else f"{a[c]:.2f}" for c in cols[1:]]) + "\n")
        log(f"{g}: pairs {a['pairs']:.0f}, Y-sites {a['y_sites']:.0f}; runs>=2 observed {a['obs_runs2']:.0f} expected {a['exp_runs2']:.1f}; "
            f"embedded {a['obs_emb2']:.0f} vs {a['exp_emb2']:.1f}; runs>=3 {a['obs_runs3']:.0f} vs {a['exp_runs3']:.1f}; runs>=5 {a['obs_runs5']:.0f} vs {a['exp_runs5']:.1f}; "
            f"Y-Y adjacencies {a['obs_adj']:.0f} vs {a['exp_adj']:.1f}")
log(f"wrote {OUT}.summary.tsv")
# Per run: E_ge = expected number of runs of >= L Y-sites in its pair under the null; a run is significant when E_ge <= ALPHA.
# Expected false significant runs = sum over pairs of E_ge(L_min), L_min = the shortest significant length of that pair.
ALPHA = 0.05
def e_exact(nA, nY, L):
    n = nA + nY
    if L > nY: return 0.0
    if L == n: return 1.0
    a = (n - L - 1) * math.exp(math.log(nA) + math.log(nA - 1) + lff(nY, L) - lff(n, L + 2)) if nA > 1 and n - L - 1 > 0 else 0.0
    return a + 2 * math.exp(math.log(nA) + lff(nY, L) - lff(n, L + 1))
def lmin_tail(nA, nY):  # (L_min, E_ge(L_min)); tail summed directly from L upward to avoid cancellation
    if nY < 2: return None, 0.0
    if nA == 0: return None, 0.0  # one run of all sites: no A-site, never called embedded, not tested
    tail = [0.0] * (nY + 2)
    for L in range(nY, 1, -1): tail[L] = tail[L + 1] + e_exact(nA, nY, L)
    for L in range(2, nY + 1):
        if tail[L] <= ALPHA: return L, tail[L]
    return None, 0.0
lm, efalse, sig = {}, 0.0, set()
for k, (arm, nA, nY) in pw.items():
    if nY >= 2:
        L, t = lmin_tail(nA, nY); lm[k] = L; efalse += t
nrun = nsig = 0
with open(f"{OUT}.runs.tsv", "w") as fo:
    fo.write("target\tdonor_arm\ta_sites\ty_sites\tn_sites\tembedded\tL_min\tsignificant\n")
    for f in sorted(glob.glob(f"{D}/chunk_*.tracts.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            k = (r["target"], r["donor_arm"]); _, nA, nY = pw[k]; L = lm.get(k); s = int(L is not None and int(r["n_sites"]) >= L)
            nrun += 1; nsig += s
            if s: sig.add((r["target"], r["donor_arm"], int(r["first_site"]), int(r["last_site"])))
            fo.write(f"{r['target']}\t{r['donor_arm']}\t{nA}\t{nY}\t{r['n_sites']}\t{r['embedded']}\t{L if L else 'NA'}\t{s}\n")
log(f"runs {nrun}, significant (E_ge <= {ALPHA} in their pair) {nsig}; expected false significant runs <= {efalse:.1f} (FDR <= {efalse / max(nsig, 1):.3f})")
EV = sys.argv[sys.argv.index("--events") + 1] if "--events" in sys.argv else None
if EV:
    byt = collections.defaultdict(list)
    for t, y, a, b in sig: byt[t].append((y, a, b))
    rows = list(csv.DictReader(open(EV), delimiter="\t")); c = collections.Counter(); keep = []
    for r in rows:
        ds = {x.split(":")[0] for x in r["donors"].split(",")}; a, b = int(r["start"]), int(r["end"])
        s = any(y in ds and a <= p and q < b for y, p, q in byt.get(r["target"], []))
        c[(r["community"], s)] += 1; c[("all", s)] += 1; c[(f"cls {r['cls']}", s)] += 1
        c[("one individual" if r["n_individuals_locus"] == "1" else "shared", s)] += 1
        c[(f"bin overlap {'yes' if r['bin_tract_overlap'] != '0' else 'no'}", s)] += 1
        if s: keep.append(r)
    for g in ["all", "one individual", "shared", "bin overlap yes", "bin overlap no", "cls embedded", "cls terminal", "cls open"] + \
             sorted({x for x, _ in c if x.startswith("C")}, key=lambda g: (len(g), g)):
        log(f"events {g}: {c[(g, True)] + c[(g, False)]}, with a significant run {c[(g, True)]}")
    with open(f"{OUT}.events_significant.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(keep)
    log(f"wrote {OUT}.events_significant.tsv ({len(keep)} events)")
