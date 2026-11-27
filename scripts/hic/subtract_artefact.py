#!/usr/bin/env python3
"""Real contact matrix minus the mapping artefact expected at its depth (referee R1-M4 / R3-M3), at the base resolution,
so that the paper's statistics can run on it after cooler zoomify --balance.

The null simulation (sim_hic_null.py, sim_porec_null.py) draws every read from the extended domain of one chromosome
end, so every simulated contact that is not inside a single end's domain is a placement error. Each such simulated
pixel is scaled to the real depth with s = real / simulated raw cis contacts inside the PHR window of an end (bin pairs
>= MIN_DIAG apart, as calibrate_artefact.py) and subtracted from the real pixel, floored at 0. Scale per pixel (the
conservative upper bound, the largest plausible artefact): both bins in domains of two different ends, max(s_A, s_B);
one bin in the domain of end A, s_A; no bin in a domain, the largest s over the ends. An end without simulated cis
contacts has no scale and is skipped in max(); a pixel whose ends all lack a scale gets the largest s over the ends
(logged). Writes OUT_COOL (float
counts, same bins) and logs the simulated and expected artefact per class, the real total, the amount removed and the
amount lost to the floor, and the real contacts left between different ends' domains.
--factor F multiplies the expected artefact (0: nothing removed, a control of the pipeline itself). --scale-resolution
sets the resolution of the scale (default 5000, as calibrate_artefact.py) when the subtraction runs at a coarser one.
usage: subtract_artefact.py DOMAINS_TSV SIM_MCOOL REAL_MCOOL OUT_COOL [--resolution 5000] [--min-diag 2] [--chunk 20000000] [--factor 1]
       [--scale-resolution 5000]"""
import argparse, sys
import cooler, numpy as np, pandas as pd
ap = argparse.ArgumentParser()
for k in ("domains", "sim", "real", "out"): ap.add_argument(k)
for k, v in (("--resolution", 5000), ("--min-diag", 2), ("--chunk", 20000000)): ap.add_argument(k, type=int, default=v)
ap.add_argument("--factor", type=float, default=1.0); ap.add_argument("--scale-resolution", type=int, default=5000); a = ap.parse_args()
log = lambda *m: print(*m, file=sys.stderr, flush=True)
D = pd.read_csv(a.domains, sep="\t")
S, R = (cooler.Cooler(f"{p}::resolutions/{a.resolution}") for p in (a.sim, a.real))
S5, R5 = (cooler.Cooler(f"{p}::resolutions/{a.scale_resolution}") for p in (a.sim, a.real))
if dict(S.chromsizes) != dict(R.chromsizes) or S.binsize != R.binsize: sys.exit("FATAL simulated and real matrices have different bins")
nb = R.info["nbins"]; own = np.full(nb, -1)
for k, d in D.iterrows():
    lo, hi = R.extent(f"{d['chrom']}:{d['ext_start']}-{d['ext_end']}")
    if (own[lo:hi] >= 0).any(): sys.exit(f"FATAL domain of {d['end']} overlaps another end's domain")
    own[lo:hi] = k
def cis(c, d):
    m = c.matrix(balance=False, sparse=True).fetch(f"{d['chrom']}:{d['phr_start']}-{d['phr_end']}").tocoo()
    return float(m.data[m.col - m.row >= a.min_diag].sum())
s = np.array([(lambda x, y: y / x if x > 0 else np.nan)(cis(S5, d), cis(R5, d)) for _, d in D.iterrows()])
log(f"{a.real} minus artefact from {a.sim} at {a.resolution} bp (scale at {a.scale_resolution} bp): {len(D)} ends, {int((own >= 0).sum())} of {nb} bins in domains;"
    f" scale median {np.nanmedian(s):.4g}, range {np.nanmin(s):.4g}-{np.nanmax(s):.4g}; ends without scale {[D['end'][k] for k in np.flatnonzero(np.isnan(s))]}")
p = S.pixels()[:]; ea, eb = own[p["bin1_id"].values], own[p["bin2_id"].values]
cls = np.where((ea >= 0) & (ea == eb), 0, np.where((ea >= 0) & (eb >= 0), 1, np.where((ea >= 0) | (eb >= 0), 2, 3)))
sa = np.where(ea >= 0, s[np.maximum(ea, 0)], np.nan); sb = np.where(eb >= 0, s[np.maximum(eb, 0)], np.nan)
with np.errstate(all="ignore"): sc = np.where(cls == 3, np.nanmax(s), np.fmax(sa, sb))
art = cls > 0; ns = art & np.isnan(sc)
log(f"artefact pixels whose ends lack a scale (given the largest scale): {int(ns.sum())} pixels, {int(p['count'].values[ns].sum())} simulated contacts")
sc = np.where(ns, np.nanmax(s), sc)
E = p["count"].values[art] * sc[art] * a.factor
log(f"factor applied to the expected artefact: {a.factor}")
names = {1: "two different ends", 2: "one end and outside the domains", 3: "outside the domains"}
for c in (1, 2, 3):
    k = cls[art] == c; log(f"class {names[c]}: {int(k.sum())} simulated pixels, {int(p['count'].values[art][k].sum())} simulated contacts, expected {E[k].sum():.1f}")
log(f"simulated contacts inside one end's domain (kept, not artefact): {int(p['count'].values[cls == 0].sum())}")
key = p["bin1_id"].values[art].astype(np.int64) * nb + p["bin2_id"].values[art]; o = np.argsort(key); key, E = key[o], E[o]
st = dict(real=0.0, removed=0.0, floor=0.0, left_between=0.0, matched=0)
def chunks():
    for lo in range(0, R.info["nnz"], a.chunk):
        c = R.pixels()[lo:lo + a.chunk]; x = c["count"].values.astype(float); st["real"] += x.sum()
        q = c["bin1_id"].values.astype(np.int64) * nb + c["bin2_id"].values
        i = np.minimum(np.searchsorted(key, q), len(key) - 1); m = key[i] == q
        y = x.copy(); y[m] = np.maximum(x[m] - E[i[m]], 0)
        st["removed"] += (x - y).sum(); st["floor"] += np.maximum(E[i[m]] - x[m], 0).sum(); st["matched"] += int(m.sum())
        b1, b2 = own[c["bin1_id"].values], own[c["bin2_id"].values]; st["left_between"] += y[(b1 >= 0) & (b2 >= 0) & (b1 != b2)].sum()
        log(f"pixels {lo}-{lo + len(c)} of {R.info['nnz']}")
        yield pd.DataFrame({"bin1_id": c["bin1_id"].values, "bin2_id": c["bin2_id"].values, "count": y})
cooler.create_cooler(a.out, R.bins()[["chrom", "start", "end"]][:], chunks(), dtypes={"count": np.float64}, ordered=True,
                     assembly=R.info.get("genome-assembly"), metadata={"source": a.real, "artefact_from": a.sim})
log(f"done {a.out}: real contacts {st['real']:.0f}; expected artefact {E.sum():.1f} on {len(E)} pixels, {st['matched']} present in the real"
    f" matrix; removed {st['removed']:.1f} ({st['removed'] / st['real']:.6f}); not removable (floor at 0) {st['floor']:.1f};"
    f" real contacts left between different ends' domains {st['left_between']:.1f}")
