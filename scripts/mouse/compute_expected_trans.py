#!/usr/bin/env python3
"""
Step 1 of the proper bouquet test (run on a host that mounts /moosefs and has the
per-stage mouse Hi-C coolers + cooltools).

Computes the genuine trans expected with cooltools.expected_trans: for each
chromosome pair, the mean balanced contact across the whole inter-chromosomal
block.  That is the "expected" a real observed/expected needs -- the column in the
seqlevel TSVs (hic_contact_norm = raw / hic_bins) is only a contact DENSITY, not
an O/E.  With this table, step 2 forms the true per-PHR-pair O/E exactly:

    O/E_mean(pair) = hic_contact_norm(pair) / expected_trans[chr_a, chr_b]

(exact, because the trans expected is constant per chromosome pair).

Output (one file per stage x resolution), written to OUT_DIR:
    expected_trans_<stage>_<res>bp.tsv   columns: region1  region2  expected

Run:  python3 scripts/mouse/compute_expected_trans.py
Needs: cooler, cooltools, pandas.  Fill COOLERS below with the real paths.
"""

import os
import sys
import pandas as pd
import cooler
import cooltools

# ----------------------------- CONFIG (EDIT) -----------------------------
# Per-stage cooler. Point at .mcool (multi-res) or single-res .cool.
# These are the Zuo-2021 mouse meiotic Hi-C maps used upstream; put the real
# /moosefs paths here (the same maps that produced the seqlevel TSVs).
COOLERS = {
    "leptotene": "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/hic/leptotene.mcool",
    "zygotene":  "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/hic/zygotene.mcool",
    "pachytene": "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/hic/pachytene.mcool",
    "diplotene": "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/hic/diplotene.mcool",
}
RESOLUTIONS = [5000, 10000, 20000, 50000]      # fine scales where the bouquet lives + 50 kb for cross-check
OUT_DIR = os.environ.get("OUT_DIR", "data/mouse_meiosis_sweep/expected_trans")
WEIGHT_NAME = os.environ.get("WEIGHT_NAME", "weight")   # balancing column in the coolers
NPROC = int(os.environ.get("NPROC", "4"))
# -------------------------------------------------------------------------


def log(msg):
    print(f"[expected_trans] {msg}", file=sys.stderr, flush=True)


def open_cooler(path, res):
    """Open a single-resolution cooler from an .mcool (preferred) or .cool."""
    if "::" in path:
        return cooler.Cooler(path)
    if path.endswith(".mcool"):
        uri = f"{path}::/resolutions/{res}"
        log(f"open {uri}")
        return cooler.Cooler(uri)
    clr = cooler.Cooler(path)
    if clr.binsize != res:
        raise ValueError(f"{path} binsize={clr.binsize} != requested {res}; give an .mcool")
    log(f"open {path} (binsize {clr.binsize})")
    return clr


def chrom_view(clr):
    """Whole-chromosome viewframe."""
    try:
        from cooltools.lib.common import make_cooler_view
        return make_cooler_view(clr)
    except Exception:                       # fallback: build it by hand
        return pd.DataFrame({
            "chrom": clr.chromnames,
            "start": 0,
            "end": [clr.chromsizes[c] for c in clr.chromnames],
            "name": clr.chromnames,
        })


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for stage, path in COOLERS.items():
        if not (os.path.exists(path) or "::" in path):
            log(f"WARNING: missing cooler for {stage}: {path} -- skipping")
            continue
        for res in RESOLUTIONS:
            log(f"stage={stage} res={res}: computing expected_trans (nproc={NPROC})")
            clr = open_cooler(path, res)
            view = chrom_view(clr)
            exp = cooltools.expected_trans(
                clr, view_df=view, clr_weight_name=WEIGHT_NAME, nproc=NPROC)
            # cooltools returns region1, region2, n_valid, count.sum, balanced.sum, balanced.avg
            avg_col = "balanced.avg" if "balanced.avg" in exp.columns else \
                      [c for c in exp.columns if c.endswith(".avg")][0]
            out = exp[["region1", "region2", avg_col]].rename(columns={avg_col: "expected"})
            out = out[out["region1"] != out["region2"]]            # trans only
            fp = os.path.join(OUT_DIR, f"expected_trans_{stage}_{res}bp.tsv")
            out.to_csv(fp, sep="\t", index=False)
            log(f"  wrote {fp}  ({len(out)} chrom pairs; "
                f"expected range {out.expected.min():.3g}..{out.expected.max():.3g})")
    log("done")


if __name__ == "__main__":
    main()
