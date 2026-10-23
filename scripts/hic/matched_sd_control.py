#!/usr/bin/env python3
"""M5 -- specificity control for the sequence-to-contact association:
is the arm-level similarity/contact relationship a property of the TELOMERIC end of an arm, or
does between-arm shared sequence predict contact equally well from INTERSTITIAL windows?

Referee #1, major point 5 (verbatim):
  "A control is needed to make the relationship specific to chromosome ends. As presented the
   result is compatible with a generic property of high-identity segmental duplications. The
   decisive control is non-subtelomeric duplications, pericentromeric and interstitial, matched
   on identity, length and copy number and run through the identical pipeline. ... If
   subtelomeric pairs show a steeper relationship the compartment claim is supported, but if not,
   the result concerns duplication and mappability rather than chromosome ends."

WHAT THIS SCRIPT IS, AND IS NOT
-------------------------------
This is a **window-class control**, not the referee's literal "matched segmental-duplication
census". The same 41 arms and the same arm-level Hi-C contact matrix are used throughout; only
the sequence windows that define the similarity axis change. See the notes file
`response/notes/M5_matched_sd_control.md` for the explicit statement of what the literal control
would still require.

TWO SUBSTRATES (both already on disk, one streaming pass each, no new alignments):

 AXIS 1 -- arm level, within the telomere-anchored flank (primary; arm-labelled partners)
   Substrate: all-vs-all.1Mb.p95.id95.len.per_window.tsv (the PHR call's own 5 kb window table;
   one row per window per flank sequence, with the arms it aligns to at >=95% identity).
   The table covers the inner 500 kb of the 1 Mb telomere-anchored flank (dist_from_telomere
   0..~495 kb), so the furthest available band is 400-500 kb from the telomere.
   Axis T = outermost 100 kb band, Axis I = innermost 100 kb band, plus the full 5-band gradient.
   Between-arm similarity for a band:
       share(A,B) = 0.5 * [ #{A-windows listing B}/#{A-windows} + #{B-windows listing A}/#{B-windows} ]
   plus a stricter variant `both_directions` = min of the two directional fractions.

 AXIS 2 -- chromosome level, genome-wide (secondary; truly interstitial windows)
   Substrate: hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz (31,202 CHM13 100 kb windows,
   cross-chromosome haplotype counts at >=98% identity).
   Axis T = windows <=100 kb from a telomere, Axis I = windows >=1 Mb from BOTH telomeres
   (the literal reviewer geometry, which axis 1 cannot reach).
   Similarity is defined by the same share() formula on partner chromosomes; chromosome-level
   rather than arm-level because that table resolves partner chromosomes but not partner arms.

CONTACT AXIS (identical for every window class)
   The published arm-level (axis 1) / chromosome-level (axis 2) mean balanced contact from the
   v5 pair tables: data/human_{CHM13_hic,HG002_hic,HG002_porec,HG002_cifi}_50000bp_seqlevel.tsv.
   Three contact regimes spanning two mapper families: HiC-Pro+bowtie2 (CHM13, HG002) and
   minimap2+pairtools (HG002 Pore-C, HG002 CiFi).

STATISTICS
   Mantel correlation (Spearman) between arm-pair similarity and arm-pair contact, with
   - ties handled by AVERAGE RANKS (scipy rankdata, method="average"); no argsort-of-argsort,
   - 10,000 arm-label row/column permutations for p,
   - 2,000 arm-block bootstrap replicates (resampling arms, which keeps the dependence induced by
     one arm appearing in many pairs) for the 95% CI,
   - a paired bootstrap CI for rho(T) - rho(I) on the same resampled arm sets.

OUTPUTS
   scripts/hic/matched_sd_control.tsv   tidy results + coverage table
   scripts/hic/matched_sd_control.json  full record (matrices, coverage, statistics, provenance)

USAGE
   python3 scripts/hic/matched_sd_control.py            # everything, writes both outputs
   python3 scripts/hic/matched_sd_control.py --perms 2000 --boot 500
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
MOOSEFS = Path("/moosefs/guarracino/HPRCv2/PHR_III")

PER_WINDOW = DATA / "all-vs-all.1Mb.p95.id95.len.per_window.tsv"
if not PER_WINDOW.exists():                       # not vendored -> read-only from moosefs
    PER_WINDOW = MOOSEFS / "all-vs-all.1Mb.p95.id95.len.per_window.tsv"
WINDOW_TABLE = DATA / "hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz"
CHROM_SIZES_FROM_FAI = MOOSEFS / "HiC/CHM13/CHM13.fa.gz.fai"
CONTACT_TABLES = {
    # label: (path, mapper / pipeline, note)
    "CHM13_HiC": (DATA / "human_CHM13_hic_50000bp_seqlevel.tsv",
                  "HiC-Pro + bowtie2 (MIN_MAPQ=0, RM_MULTI=0)", "primary"),
    "HG002_HiC": (DATA / "human_HG002_hic_50000bp_seqlevel.tsv",
                  "HiC-Pro + bowtie2 (MIN_MAPQ=0, RM_MULTI=0)", ""),
    "HG002_PoreC": (DATA / "human_HG002_porec_50000bp_seqlevel.tsv",
                    "minimap2 + pairtools (no MAPQ filter)", ""),
    "HG002_CiFi": (DATA / "human_HG002_cifi_50000bp_seqlevel.tsv",
                   "minimap2 + pairtools (no MAPQ filter)", ""),
}
OUT_TSV = ROOT / "scripts/hic/matched_sd_control.tsv"
OUT_JSON = ROOT / "scripts/hic/matched_sd_control.json"

MIN_WINDOWS = 20          # windows an arm needs inside a band to enter that band's matrix
MIN_HAP = 2               # haplotype support for a partner chromosome to count (axis 2)
# (name, lo, hi, dup_only): hi=None means "distance >= lo"; dup_only=True keeps only windows
# that actually carry a cross-chromosome duplicate, i.e. the closest available analogue of the
# referee's "interstitial segmental duplications" set.
GW_BANDS_AXIS2 = [
    ("T_telomeric_100kb", 0, 100_000, False),
    ("T_telomeric_500kb", 0, 500_000, False),
    ("T_telomeric_100kb_duponly", 0, 100_000, True),
    ("I_interstitial_1Mb", 1_000_000, None, False),
    ("I_interstitial_2Mb", 2_000_000, None, False),
    ("I_interstitial_1Mb_duponly", 1_000_000, None, True),
]
SEED = 20260518
N_PERM = 10_000
N_BOOT = 2_000

# axis 1 bands (dist_from_telomere, bp)
# NOTE (see the containment check in build_axis1): the per-window table contains windows inside
# the *called PHR intervals*, not the whole flank, so these bands are PHR-internal depth bands,
# not "distance from the telomere into the interstitial genome".
BANDS = [
    ("PHR_0_100kb_from_telomere", 0, 100_000),
    ("PHR_100_200kb_from_telomere", 100_000, 200_000),
    ("PHR_200_300kb_from_telomere", 200_000, 300_000),
    ("PHR_300_400kb_from_telomere", 300_000, 400_000),
    ("PHR_400_500kb_from_telomere", 400_000, 500_000),
]


def die(msg: str) -> "NoReturn":
    sys.stderr.write(f"\nFATAL: {msg}\n")
    raise SystemExit(1)


def log(msg: str) -> None:
    print(msg, flush=True)


ARM_RE = re.compile(r"(chr(?:[0-9]+|[XYM]))[_]?([pq])")


def norm_arm(tok: str) -> str:
    """'chr10q' / 'chr10_q' -> 'chr10_q'."""
    m = ARM_RE.match(tok.strip())
    if not m:
        die(f"cannot parse an arm label from {tok!r}")
    return f"{m.group(1)}_{m.group(2)}"


def arm_from_seq(seq: str) -> str:
    """'CHM13#0#chr10:134254995-134754994_chr10_qarm' -> 'chr10_q'."""
    m = re.search(r"_(chr(?:[0-9]+|[XYM]))_([pq])arm$", seq.strip())
    if not m:
        die(f"cannot parse an arm label from sequence name {seq!r}")
    return f"{m.group(1)}_{m.group(2)}"


# ============================================================ axis 1: within-flank windows
def phr_regions(path: Path | None = None) -> dict:
    """seq -> (region_start, region_end) relative to the flank start, from the PHR-level table."""
    path = path or (DATA / "all-vs-all.1Mb.p95.id95.len.tsv")
    if not path.exists():
        path = MOOSEFS / "all-vs-all.1Mb.p95.id95.len.tsv"
    if not path.exists():
        die(f"missing PHR-level table {path}")
    out = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            rs, re_ = c[idx["region_start"]].strip(), c[idx["region_end"]].strip()
            if rs in (".", "") or re_ in (".", ""):
                continue
            out[c[idx["seq"]]] = (int(rs), int(re_))
    return out

def build_axis1(path: Path = PER_WINDOW, min_windows: int = MIN_WINDOWS) -> dict:
    """One streaming pass over the per-window table -> per-band between-arm sharing."""
    if not path.exists():
        die(f"missing per-window table {path}")
    band_idx = {}
    for k, (name, lo, hi) in enumerate(BANDS):
        band_idx[(lo, hi)] = k
    phr = phr_regions()
    n_rows = 0
    n_inside_phr = 0
    n_unknown_seq = 0
    max_dist = defaultdict(int)                             # arm -> max dist_from_telomere seen
    cov = defaultdict(int)                                  # (arm, band) -> n windows
    share = defaultdict(int)                                # (band, arm, partner) -> n windows
    n_partners = defaultdict(list)                          # (arm, band) -> partner counts
    passes = defaultdict(int)
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        need = ["seq", "arm", "window_start", "window_end", "dist_from_telomere",
                "arms_involved", "passes_threshold"]
        if not set(need).issubset(hdr):
            die(f"unexpected per-window header: {hdr}")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            n_rows += 1
            c = line.rstrip("\n").split("\t")
            d = int(c[idx["dist_from_telomere"]])
            b = band_idx.get((d // 100_000 * 100_000, d // 100_000 * 100_000 + 100_000))
            if b is None:
                continue
            arm = arm_from_seq(c[idx["seq"]])
            seq = c[idx["seq"]]
            ext = phr.get(seq)
            ws, we = int(c[idx["window_start"]]), int(c[idx["window_end"]])
            if ext is None:
                n_unknown_seq += 1
            elif ext[0] <= ws and we <= ext[1]:
                n_inside_phr += 1
            if d > max_dist[arm]:
                max_dist[arm] = d
            cov[(arm, b)] += 1
            if c[idx["passes_threshold"]].strip() == "True":
                passes[(arm, b)] += 1
            raw = c[idx["arms_involved"]].strip()
            partners = [] if raw in ("", ".") else [p for p in raw.split(",") if p and p != "."]
            n_partners[(arm, b)].append(len(partners))
            for p in partners:
                pa = norm_arm(p)
                if pa != arm:
                    share[(b, arm, pa)] += 1
    log(f"axis1: {n_rows} window rows from {path.name}, "
        f"{len({a for a, _ in cov})} arms, bands {[b[0] for b in BANDS]}")
    log(f"axis1 containment: {n_inside_phr}/{n_rows - n_unknown_seq} windows lie inside the "
        f"sequence's own called PHR interval ({100.0 * n_inside_phr / max(n_rows - n_unknown_seq, 1):.2f}%); "
        f"{n_unknown_seq} rows belong to sequences with no called PHR")

    # per-band similarity matrices
    out = {}
    for b, (name, lo, hi) in enumerate(BANDS):
        arms = sorted({a for (a, bb) in cov if bb == b and cov[(a, bb)] >= min_windows})
        dropped = sorted({a for (a, bb) in cov if bb == b} - set(arms))
        sim = {}
        sim_bd = {}
        for i, a in enumerate(arms):
            for c in arms[i + 1:]:
                fa = share.get((b, a, c), 0) / cov[(a, b)]
                fc = share.get((b, c, a), 0) / cov[(c, b)]
                sim[(a, c)] = 0.5 * (fa + fc)
                sim_bd[(a, c)] = min(fa, fc)
        med_partners = {a: float(np.median(n_partners[(a, b)])) for a in arms}
        out[name] = dict(band=name, lo=lo, hi=hi, arms=arms, dropped_arms=dropped,
                         cov={a: cov[(a, b)] for a in arms},
                         cov_dropped={a: cov[(a, b)] for a in dropped},
                         passes={a: passes[(a, b)] for a in arms},
                         median_partner_arms=med_partners,
                         sim=sim, sim_both=sim_bd)
        log(f"  {name:22s} arms={len(arms):3d} (dropped {len(dropped)}: {','.join(dropped) or '-'}) "
            f"median windows/arm={np.median(list(out[name]['cov'].values())):.0f} "
            f"median partner arms/window={np.median(list(med_partners.values())):.0f} "
            f"nonzero pairs={sum(1 for v in sim.values() if v > 0)}/{len(sim)}")
    out["_containment"] = dict(n_rows=n_rows, n_inside_phr=n_inside_phr,
                              n_unknown_seq=n_unknown_seq,
                              max_dist_from_telomere_per_arm=dict(max_dist))
    return out


# ============================================================ axis 2: genome-wide windows
def build_axis2(path: Path = WINDOW_TABLE) -> dict:
    """Genome-wide 100 kb CHM13 windows: telomeric bands vs >=1/2 Mb-interstitial bands,
    chromosome-level between-chromosome sharing."""
    if not path.exists():
        die(f"missing genome-wide window table {path}")
    sizes = {}
    for line in open(CHROM_SIZES_FROM_FAI):
        c = line.split("\t")
        sizes[c[0]] = int(c[1])
    labels = [b[0] for b in GW_BANDS_AXIS2]
    hits = {t: defaultdict(int) for t in labels}
    cov = {t: Counter() for t in labels}
    n_rows = 0
    n_dropped_chr = 0
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            n_rows += 1
            c = line.rstrip("\n").split("\t")
            chrom = c[idx["chrom"]].split("#")[-1]
            if chrom not in sizes:
                n_dropped_chr += 1
                continue
            s, e = int(c[idx["start"]]), int(c[idx["end"]])
            d_telo = min(s, sizes[chrom] - e)
            for name, lo, hi, dup_only in GW_BANDS_AXIS2:
                if hi is not None:
                    if not (d_telo <= hi and lo == 0):
                        continue
                else:
                    if not (d_telo >= lo):
                        continue
                field = c[idx["chroms-num_haplotypes"]]
                partners = []
                if field:
                    for kv in field.split(","):
                        pc, _, v = kv.rpartition("-")
                        if pc and pc != chrom and int(v) >= MIN_HAP:
                            partners.append(pc)
                if dup_only and not partners:
                    continue                 # keep only duplicated windows
                cov[name][chrom] += 1
                for pc in partners:
                    hits[name][(chrom, pc)] += 1
    log(f"axis2: {n_rows} windows from {path.name} ({n_dropped_chr} non-CHM13-chrom rows "
        f"dropped, i.e. chrY which the masked-Y fasta lacks)")
    out = {}
    for name, lo, hi, dup_only in GW_BANDS_AXIS2:
        cls = name
        chroms = sorted(cov[cls])
        sim = {}
        for i, a in enumerate(chroms):
            for b in chroms[i + 1:]:
                fa = hits[cls].get((a, b), 0) / cov[cls][a]
                fb = hits[cls].get((b, a), 0) / cov[cls][b]
                sim[(a, b)] = 0.5 * (fa + fb)
        vals = np.array(list(sim.values()))
        out[cls] = dict(band=cls, lo=lo, hi=hi, chroms=chroms, cov=dict(cov[cls]), sim=sim,
                        n_nonzero=int((vals > 0).sum()),
                        median_sim=float(np.median(vals)) if len(vals) else float("nan"),
                        max_sim=float(vals.max()) if len(vals) else float("nan"))
        log(f"  {cls:22s} chromosomes={len(chroms):2d} windows={sum(cov[cls].values()):6d} "
            f"nonzero pairs={out[cls]['n_nonzero']:3d}/{len(sim):3d} "
            f"sim[max]={out[cls]['max_sim']:.3f} sim[median]={out[cls]['median_sim']:.4f}")
    return out


# ============================================================ contact axis
def read_contact(label: str) -> dict:
    path, mapper, note = CONTACT_TABLES[label]
    if not path.exists():
        die(f"missing contact table {path}")
    acc = defaultdict(list)
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            a, b = c[idx["arm_a"]], c[idx["arm_b"]]
            v = float(c[idx["hic_contact_norm"]])
            acc[(a, b) if a < b else (b, a)].append(v)
    return {k: float(np.mean(v)) for k, v in acc.items()}


def contact_to_chrom(arm_contact: dict) -> dict:
    acc = defaultdict(list)
    for (a, b), v in arm_contact.items():
        ca, cb = a.split("_")[0], b.split("_")[0]
        if ca == cb:
            continue
        acc[(ca, cb) if ca < cb else (cb, ca)].append(v)
    return {k: float(np.mean(v)) for k, v in acc.items()}


# ============================================================ statistics
def _avg_ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks (ties get the mean rank). Explicitly NOT argsort-of-argsort."""
    from scipy.stats import rankdata
    return rankdata(x, method="average")


def mantel(sim: dict, con: dict, n_perm: int = N_PERM, n_boot: int = N_BOOT,
           seed: int = SEED) -> dict:
    """Mantel test between a similarity and a contact table, keyed by unordered node pairs.

    rho is Spearman on arm-pair values with average ranks; the null is generated by permuting
    node labels (row and column permutation jointly), which preserves the dependency structure of
    the matrix; the CI comes from resampling nodes with replacement.
    """
    pairs = sorted(set(sim) & set(con))
    if len(pairs) < 10:
        return dict(rho=float("nan"), p_perm=float("nan"), ci_lo=float("nan"),
                    ci_hi=float("nan"), n_pairs=len(pairs), note="too few pairs")
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    S = np.full((n, n), np.nan)
    C = np.full((n, n), np.nan)
    for (a, b) in pairs:
        i, j = nidx[a], nidx[b]
        S[i, j] = S[j, i] = sim[(a, b)]
        C[i, j] = C[j, i] = con[(a, b)]
    iu = np.triu_indices(n, 1)
    mask = np.isfinite(S[iu]) & np.isfinite(C[iu])
    xs, ys = S[iu][mask], C[iu][mask]
    rx, ry = _avg_ranks(xs), _avg_ranks(ys)

    def rho_of(a, b):
        a = a - a.mean()
        b = b - b.mean()
        da, db = a.std(), b.std()
        if da == 0 or db == 0:
            return float("nan")
        return float(np.mean((a / da) * (b / db)))

    obs = rho_of(rx, ry)
    rng = np.random.default_rng(seed)
    # label permutation null
    null = np.empty(n_perm)
    for k in range(n_perm):
        p = rng.permutation(n)
        Sp = S[np.ix_(p, p)]
        x = Sp[iu][mask]
        null[k] = rho_of(_avg_ranks(x), ry)
    cnt = int(np.sum(null >= obs))
    p_perm = (cnt + 1) / (n_perm + 1)
    # node bootstrap
    boot = np.empty(n_boot)
    for k in range(n_boot):
        draw = rng.integers(0, n, n)
        Sp = S[np.ix_(draw, draw)]
        Cp = C[np.ix_(draw, draw)]
        m = np.isfinite(Sp[iu]) & np.isfinite(Cp[iu])
        if m.sum() < 10:
            boot[k] = np.nan
            continue
        boot[k] = rho_of(_avg_ranks(Sp[iu][m]), _avg_ranks(Cp[iu][m]))
    ok = boot[np.isfinite(boot)]
    lo, hi = (float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))) \
        if len(ok) > 10 else (float("nan"), float("nan"))
    uniq_s = len(set(xs.tolist()))
    uniq_c = len(set(ys.tolist()))
    return dict(rho=obs, p_perm=float(p_perm), ci_lo=lo, ci_hi=hi,
                n_pairs=int(m.sum()), n_nodes=n,
                n_distinct_sim=uniq_s, n_distinct_con=uniq_c,
                n_tied_sim=int(m.sum()) - uniq_s, n_tied_con=int(m.sum()) - uniq_c,
                n_zero_sim=int(np.sum(xs == 0)), n_zero_con=int(np.sum(ys == 0)),
                null_mean=float(null.mean()), null_sd=float(null.std()),
                n_perm=n_perm, n_boot=int(len(ok)), pct_perm_ge=float(cnt / n_perm))


def paired_delta(sim_t: dict, sim_i: dict, con: dict, n_boot: int = N_BOOT,
                 seed: int = SEED) -> dict:
    """Bootstrapped CI of rho(T) - rho(I) on the SAME resampled node sets."""
    pairs = sorted(set(sim_t) & set(sim_i) & set(con))
    if len(pairs) < 10:
        return dict(delta=float("nan"), lo=float("nan"), hi=float("nan"), p=float("nan"))
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    def mat(d):
        M = np.full((n, n), np.nan)
        for (a, b) in pairs:
            i, j = nidx[a], nidx[b]
            M[i, j] = M[j, i] = d[(a, b)]
        return M

    ST, SI, C = mat(sim_t), mat(sim_i), mat(con)
    iu = np.triu_indices(n, 1)
    mask = np.isfinite(ST[iu]) & np.isfinite(SI[iu]) & np.isfinite(C[iu])

    def rho(a, b):
        a = _avg_ranks(a)
        b = _avg_ranks(b)
        a, b = a - a.mean(), b - b.mean()
        if a.std() == 0 or b.std() == 0:
            return float("nan")
        return float(np.mean((a / a.std()) * (b / b.std())))

    obs = rho(ST[iu][mask], C[iu][mask]) - rho(SI[iu][mask], C[iu][mask])

    def rho_full(M):
        x, y = M[iu][mask], C[iu][mask]
        return rho(x, y)

    rng = np.random.default_rng(seed)
    ds = []
    for _ in range(n_boot):
        draw = rng.integers(0, n, n)
        A = ST[np.ix_(draw, draw)]
        B = SI[np.ix_(draw, draw)]
        Cc = C[np.ix_(draw, draw)]
        m = np.isfinite(A[iu]) & np.isfinite(B[iu]) & np.isfinite(Cc[iu])
        if m.sum() < 10:
            continue
        r1 = rho(A[iu][m], Cc[iu][m])
        r2 = rho(B[iu][m], Cc[iu][m])
        if np.isfinite(r1) and np.isfinite(r2):
            ds.append(r1 - r2)
    ds = np.array(ds)
    if len(ds) < 10:
        return dict(delta=float(obs), lo=float("nan"), hi=float("nan"), p=float("nan"))
    return dict(delta=float(obs), lo=float(np.percentile(ds, 2.5)),
                hi=float(np.percentile(ds, 97.5)),
                p=float(2 * min((ds <= 0).mean(), (ds >= 0).mean())), n_boot=int(len(ds)))


def restrict_nonzero(sim: dict, con: dict) -> tuple:
    """Keep only node pairs with a nonzero similarity (sensitivity for cells whose x axis is
    dominated by ties at zero)."""
    keep = {k: v for k, v in sim.items() if v > 0 and k in con}
    return keep, {k: con[k] for k in keep}


def gate_published_pair_rho() -> dict:
    """Reproduction gate: the published CHM13 pair-level Spearman is 0.716."""
    from scipy.stats import spearmanr
    rows = []
    with open(CONTACT_TABLES["CHM13_HiC"][0]) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            rows.append((float(c[idx["jaccard"]]), float(c[idx["hic_contact_norm"]])))
    a = np.array(rows)
    r = spearmanr(a[:, 0], a[:, 1])
    out = dict(gate="G1_published_pair_rho", n_pairs=int(len(a)), rho=float(r.statistic),
               p=float(r.pvalue), published=0.716, passed=bool(abs(r.statistic - 0.716) < 0.002))
    log(f"GATE G1: published CHM13 pair-level rho reproduced = {r.statistic:.5f} "
        f"(published 0.716) on {len(a)} pairs -> {'PASS' if out['passed'] else 'FAIL'}")
    return out


def gate_axis1_vs_published(pubdist: Path, band: dict) -> dict:
    """Sanity gate: the telomeric-band arm similarity must agree in rank with the published
    arm-level Jaccard distance matrix (1 - dist), otherwise the window aggregation is wrong."""
    from scipy.stats import spearmanr
    if not pubdist.exists():
        return dict(gate="G_axis1_vs_published", skipped=True, reason=f"missing {pubdist}")
    D = {}
    with open(pubdist) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        names = hdr[1:]
        for line in f:
            c = line.rstrip("\n").split("\t")
            for n2, v in zip(names, c[1:]):
                if v in ("", "NA", "nan"):
                    continue
                k = (c[0], n2) if c[0] < n2 else (n2, c[0])
                D[k] = float(v)
    pairs = sorted(set(D) & set(band["sim"]))
    x = [band["sim"][p] for p in pairs]
    y = [1.0 - D[p] for p in pairs]
    r = spearmanr(x, y)
    out = dict(gate="G_axis1_vs_published_arm_jaccard", n_pairs=len(pairs),
               rho=float(r.statistic), p=float(r.pvalue),
               passed=bool(r.statistic > 0.5),
               note="telomeric-band share vs 1 - published arm-level Jaccard distance")
    log(f"GATE G: telomeric-band similarity vs published arm matrix: rho={r.statistic:+.4f} "
        f"on {len(pairs)} arm pairs -> {'PASS' if out['passed'] else 'FAIL'}")
    return out


def gate_axis2_vs_axis1(axis1_band: dict, axis2_band: dict) -> dict:
    """Sanity gate: the genome-wide telomeric, chromosome-level sharing must agree in rank with
    the arm-level telomeric sharing aggregated to chromosome pairs. If it does not, the two
    substrates are not measuring the same thing and the axis-2 contrast is uninterpretable."""
    from scipy.stats import spearmanr
    acc = defaultdict(list)
    for (a, b), v in axis1_band["sim"].items():
        ca, cb = a.split("_")[0], b.split("_")[0]
        if ca == cb:
            continue
        acc[(ca, cb) if ca < cb else (cb, ca)].append(v)
    agg = {k: float(np.mean(v)) for k, v in acc.items()}
    pairs = sorted(set(agg) & set(axis2_band["sim"]))
    x = [agg[p] for p in pairs]
    y = [axis2_band["sim"][p] for p in pairs]
    r = spearmanr(x, y)
    out = dict(gate="G_axis2_vs_axis1_chromosome_sharing", n_pairs=len(pairs),
               rho=float(r.statistic), p=float(r.pvalue), passed=bool(r.statistic > 0.5))
    log(f"GATE G: axis2 telomeric sharing vs axis1 telomeric sharing aggregated to chromosomes: "
        f"rho={r.statistic:+.4f} on {len(pairs)} chromosome pairs -> "
        f"{'PASS' if out['passed'] else 'FAIL'}")
    return out


# ============================================================ driver
ROW_COLS = ["analysis", "substrate", "axis", "band", "dist_from_telomere_bp", "contact_regime",
            "mapper_pipeline", "similarity_def", "n_nodes", "n_pairs", "n_distinct_sim",
            "n_distinct_con", "n_tied_sim", "n_tied_con", "n_zero_sim", "n_zero_con",
            "rho", "p_perm", "ci_lo", "ci_hi", "null_mean", "n_perm", "n_boot",
            "delta_T_minus_I", "delta_ci_lo", "delta_ci_hi", "delta_p", "note"]


def _fmt(v):
    if isinstance(v, float):
        return "nan" if not np.isfinite(v) else f"{v:.6g}"
    return "" if v is None else str(v)


def _row(**kw) -> dict:
    r = {c: "" for c in ROW_COLS}
    r.update({k: v for k, v in kw.items() if k in r})
    return r


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perms", type=int, default=N_PERM)
    ap.add_argument("--boot", type=int, default=N_BOOT)
    ap.add_argument("--min-windows", type=int, default=MIN_WINDOWS)
    a = ap.parse_args()

    rows = []
    record = {}

    gates = [gate_published_pair_rho()]
    axis1 = build_axis1(min_windows=a.min_windows)
    gates.append(gate_axis1_vs_published(DATA / "hprcv2.1Mb.subtelo.arm_dist_matrix.tsv",
                                         axis1["PHR_0_100kb_from_telomere"]))
    record["gates"] = gates

    # ---------------- coverage table (axis 1) ----------------
    for name, d in axis1.items():
        if name.startswith("_"):
            continue
        if name != "PHR_0_100kb_from_telomere":
            continue
        for arm in sorted(set(d["cov"]) | set(d["cov_dropped"])):
            n = d["cov"].get(arm, d["cov_dropped"].get(arm, 0))
            rows.append(_row(analysis="coverage_axis1", substrate="per_window_5kb",
                             axis="arm", band=name, dist_from_telomere_bp=f"{d['lo']}-{d['hi']}",
                             n_pairs=n, note=f"arm={arm}; passes_threshold="
                                            f"{d['passes'].get(arm, '')};"
                                            f"median_partner_arms="
                                            f"{d['median_partner_arms'].get(arm, 'NA')}"))
    # coverage for every band, one row per band summarising
    for name, d in axis1.items():
        if name.startswith("_"):
            continue
        rows.append(_row(analysis="band_summary_axis1", substrate="per_window_5kb", axis="arm",
                         band=name, dist_from_telomere_bp=f"{d['lo']}-{d['hi']}",
                         n_nodes=len(d["arms"]), n_pairs=len(d["sim"]),
                         n_distinct_sim=len(set(d["sim"].values())),
                         note=f"arms={len(d['arms'])};dropped={','.join(d['dropped_arms']) or '-'};"
                              f"median_windows_per_arm="
                              f"{int(np.median(list(d['cov'].values())))};"
                              f"median_partner_arms_per_window="
                              f"{np.median(list(d['median_partner_arms'].values())):.0f}"))
    cov_axis1 = {n: dict(n_arms=len(d["arms"]), dropped=d["dropped_arms"],
                         windows_per_arm=d["cov"], windows_per_arm_dropped=d["cov_dropped"],
                         median_partner_arms=d["median_partner_arms"])
                 for n, d in axis1.items() if not n.startswith("_")}
    record["axis1_coverage"] = cov_axis1
    record["axis1_containment"] = axis1["_containment"]

    # ---------------- axis 1 associations ----------------
    for label, (path, mapper, note) in CONTACT_TABLES.items():
        con = read_contact(label)
        for name, d in axis1.items():
            if name.startswith("_"):
                continue
            for sdef, key in (("share_symmetric", "sim"), ("share_both_directions", "sim_both")):
                m = mantel(d[key], con, n_perm=a.perms, n_boot=a.boot)
                rows.append(_row(analysis="mantel", substrate="per_window_5kb", axis="arm",
                                 band=name, dist_from_telomere_bp=f"{d['lo']}-{d['hi']}",
                                 contact_regime=label, mapper_pipeline=mapper,
                                 similarity_def=sdef, **m))
                log(f"  axis1 {name:22s} {label:11s} {sdef:22s} n={m['n_pairs']:4d} "
                    f"rho={m['rho']:+.4f} p={m['p_perm']:.4g} "
                    f"CI[{m['ci_lo']:+.3f},{m['ci_hi']:+.3f}]")
                if sdef == "share_symmetric":
                    sn, cn_ = restrict_nonzero(d[key], con)
                    if len(sn) >= 10:
                        mn = mantel(sn, cn_, n_perm=min(a.perms, 2000), n_boot=a.boot)
                        rows.append(_row(analysis="mantel_nonzero_sim",
                                         substrate="per_window_5kb", axis="arm", band=name,
                                         dist_from_telomere_bp=f"{d['lo']}-{d['hi']}",
                                         contact_regime=label, mapper_pipeline=mapper,
                                         similarity_def=sdef + "_nonzero_only", **mn))
                        log(f"    nonzero-only n={mn['n_pairs']:4d} rho={mn['rho']:+.4f} "
                            f"p={mn['p_perm']:.4g} CI[{mn['ci_lo']:+.3f},{mn['ci_hi']:+.3f}]")
        # paired delta T vs I, same contact regime
        for sdef, key in (("share_symmetric", "sim"), ("share_both_directions", "sim_both")):
            dl = paired_delta(axis1["PHR_0_100kb_from_telomere"][key],
                              axis1["PHR_400_500kb_from_telomere"][key], con, n_boot=a.boot)
            rows.append(_row(analysis="delta_T_vs_I", substrate="per_window_5kb", axis="arm",
                             band="PHRdepth(0-100kb) - PHRdepth(400-500kb)",
                             contact_regime=label, mapper_pipeline=mapper, similarity_def=sdef,
                             delta_T_minus_I=dl["delta"], delta_ci_lo=dl["lo"],
                             delta_ci_hi=dl["hi"], delta_p=dl["p"], n_boot=dl.get("n_boot"),
                             note="paired arm bootstrap on the arms present in both PHR-depth bands; positive = outer association steeper"))
            log(f"  DELTA {label:11s} {sdef:22s} rho_T - rho_I = {dl['delta']:+.4f} "
                f"CI[{dl['lo']:+.3f},{dl['hi']:+.3f}] p={dl['p']:.4g}")

    # ---------------- axis 2 (genome-wide, chromosome level) ----------------
    axis2 = build_axis2()
    gates.append(gate_axis2_vs_axis1(axis1["PHR_0_100kb_from_telomere"],
                                     axis2["T_telomeric_500kb"]))
    record["axis2_coverage"] = {k: dict(n_chroms=len(d["chroms"]), windows=d["cov"])
                                for k, d in axis2.items()}
    for name, d in axis2.items():
        rows.append(_row(analysis="band_summary_axis2", substrate="CHM13_100kb_genomewide",
                         axis="chromosome", band=name,
                         dist_from_telomere_bp=(f"0-{d['hi']}" if d["hi"]
                                                else f">={d['lo']}"),
                         n_nodes=len(d["chroms"]), n_pairs=len(d["sim"]),
                         n_distinct_sim=len(set(d["sim"].values())),
                         note=f"windows={sum(d['cov'].values())};nonzero_pairs="
                              f"{d['n_nonzero']}/{len(d['sim'])};"
                              f"median_sim={np.median(list(d['sim'].values())):.4f}"))
    for label, (path, mapper, note) in CONTACT_TABLES.items():
        con_arm = read_contact(label)
        con_chr = contact_to_chrom(con_arm)
        for name, d in axis2.items():
            dib = (f"0-{d['hi']}" if d["hi"] else f">={d['lo']}")
            m = mantel(d["sim"], con_chr, n_perm=a.perms, n_boot=a.boot)
            rows.append(_row(analysis="mantel", substrate="CHM13_100kb_genomewide",
                             axis="chromosome", band=name, dist_from_telomere_bp=dib,
                             contact_regime=label, mapper_pipeline=mapper,
                             similarity_def="share_symmetric", **m))
            log(f"  axis2 {name:22s} {label:11s} n={m['n_pairs']:4d} rho={m['rho']:+.4f} "
                f"p={m['p_perm']:.4g} CI[{m['ci_lo']:+.3f},{m['ci_hi']:+.3f}]")
            sn, cn_ = restrict_nonzero(d["sim"], con_chr)
            if len(sn) >= 10:
                mn = mantel(sn, cn_, n_perm=min(a.perms, 2000), n_boot=a.boot)
                rows.append(_row(analysis="mantel_nonzero_sim",
                                 substrate="CHM13_100kb_genomewide", axis="chromosome",
                                 band=name, dist_from_telomere_bp=dib, contact_regime=label,
                                 mapper_pipeline=mapper,
                                 similarity_def="share_symmetric_nonzero_only", **mn))
                log(f"    nonzero-only n={mn['n_pairs']:4d} rho={mn['rho']:+.4f} "
                    f"p={mn['p_perm']:.4g} CI[{mn['ci_lo']:+.3f},{mn['ci_hi']:+.3f}]")
        # sensitivity: exclude acrocentric + sex chromosomes (nucleolar association, arm type)
        drop = {"chr13", "chr14", "chr15", "chr21", "chr22", "chrX", "chrY"}
        for name in ("T_telomeric_500kb", "I_interstitial_1Mb",
                     "T_telomeric_100kb_duponly", "I_interstitial_1Mb_duponly"):
            d = axis2[name]
            keep = set(d["chroms"]) - drop
            sub = {k: v for k, v in d["sim"].items() if k[0] in keep and k[1] in keep}
            m = mantel(sub, con_chr, n_perm=a.perms, n_boot=a.boot)
            rows.append(_row(analysis="mantel_no_acro_sex", substrate="CHM13_100kb_genomewide",
                             axis="chromosome", band=name,
                             dist_from_telomere_bp=(f"0-{d['hi']}" if d["hi"]
                                                    else f">={d['lo']}"),
                             contact_regime=label, mapper_pipeline=mapper,
                             similarity_def="share_symmetric",
                             note="chr13,14,15,21,22,X,Y excluded", **m))
            log(f"  axis2-noacro {name:22s} {label:11s} n={m['n_pairs']:3d} "
                f"rho={m['rho']:+.4f} p={m['p_perm']:.4g} "
                f"CI[{m['ci_lo']:+.3f},{m['ci_hi']:+.3f}]")
        dl2 = paired_delta(axis2["T_telomeric_100kb_duponly"]["sim"],
                           axis2["I_interstitial_1Mb_duponly"]["sim"], con_chr, n_boot=a.boot)
        rows.append(_row(analysis="delta_T_vs_I_duponly", substrate="CHM13_100kb_genomewide",
                         axis="chromosome", band="T(<=100kb,dup) - I(>=1Mb,dup)",
                         contact_regime=label, mapper_pipeline=mapper,
                         similarity_def="share_symmetric", delta_T_minus_I=dl2["delta"],
                         delta_ci_lo=dl2["lo"], delta_ci_hi=dl2["hi"], delta_p=dl2["p"],
                         n_boot=dl2.get("n_boot"),
                         note="both bands restricted to windows carrying a cross-chromosome "
                              "duplicate; paired chromosome bootstrap"))
        log(f"  DELTA dup-only {label:11s} rho_T - rho_I = {dl2['delta']:+.4f} "
            f"CI[{dl2['lo']:+.3f},{dl2['hi']:+.3f}] p={dl2['p']:.4g}")
        dl = paired_delta(axis2["T_telomeric_100kb"]["sim"],
                          axis2["I_interstitial_1Mb"]["sim"], con_chr, n_boot=a.boot)
        rows.append(_row(analysis="delta_T_vs_I", substrate="CHM13_100kb_genomewide",
                         axis="chromosome", band="T(<=100kb) - I(>=1Mb)",
                         contact_regime=label, mapper_pipeline=mapper,
                         similarity_def="share_symmetric", delta_T_minus_I=dl["delta"],
                         delta_ci_lo=dl["lo"], delta_ci_hi=dl["hi"], delta_p=dl["p"],
                         n_boot=dl.get("n_boot"),
                         note="paired chromosome bootstrap; positive = telomeric steeper"))
        log(f"  DELTA axis2 {label:11s} rho_T - rho_I = {dl['delta']:+.4f} "
            f"CI[{dl['lo']:+.3f},{dl['hi']:+.3f}] p={dl['p']:.4g}")

    # ---------------- write ----------------
    with open(OUT_TSV, "w") as f:
        f.write("\t".join(ROW_COLS) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r[c]) for c in ROW_COLS) + "\n")
    log(f"wrote {len(rows)} rows -> {OUT_TSV}")

    old = {}
    if OUT_JSON.exists():
        try:
            old = json.loads(OUT_JSON.read_text())
        except json.JSONDecodeError:
            old = {}
    old.update(record)
    for stale in ("sequence_axis", "window_classes", "extraction", "alignment", "contact",
                  "oracle", "liftover", "statistics"):
        old.pop(stale, None)          # keys produced by the earlier, abandoned wfmash design
    old["provenance"] = {
        "per_window_table": str(PER_WINDOW),
        "genomewide_window_table": str(WINDOW_TABLE),
        "contact_tables": {k: str(v[0]) for k, v in CONTACT_TABLES.items()},
        "mappers": {k: v[1] for k, v in CONTACT_TABLES.items()},
        "min_windows_per_arm_per_band": a.min_windows,
        "min_haplotypes_axis2": MIN_HAP,
        "n_perm": a.perms, "n_boot": a.boot, "seed": SEED,
        "python": sys.version.split()[0], "numpy": np.__version__,
        "statistics": "Spearman on arm-pair values, average ranks (ties averaged); "
                      "p from arm-label row/column permutation; CI from arm bootstrap",
    }
    old["axis1_matrices"] = {n: {f"{k[0]}|{k[1]}": v for k, v in d["sim"].items()}
                             for n, d in axis1.items() if not n.startswith("_")}
    old["axis2_matrices"] = {n: {f"{k[0]}|{k[1]}": v for k, v in d["sim"].items()}
                             for n, d in axis2.items()}
    OUT_JSON.write_text(json.dumps(old, indent=2, default=str) + "\n")
    log(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
