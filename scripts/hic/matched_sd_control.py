#!/usr/bin/env python3
"""M5: the decisive specificity control -- are subtelomeres special, or is the
sequence-to-contact relationship a generic property of high-identity segmental duplications?

Referee #1, major point 5 (verbatim):
  "A control is needed to make the relationship specific to chromosome ends. As presented the
   result is compatible with a generic property of high-identity segmental duplications. The
   decisive control is non-subtelomeric duplications, pericentromeric and interstitial, matched
   on identity, length and copy number and run through the identical pipeline. ... If
   subtelomeric pairs show a steeper relationship the compartment claim is supported, but if
   not, the result concerns duplication and mappability rather than chromosome ends."

Referee #3, major point 3 asks for the same control from the mappability side (restrict to
uniquely mappable/anchored reads, aggregate-level analysis, control for arm type, telomere
position, chromosome size and general chromosome-end contact propensity).

DESIGN
------
Units are CHM13 v2.0 100 kb windows (the window grid of the published Fig. 1A analysis), assigned
to four mutually exclusive classes by position and by whether they carry cross-chromosome
high-identity (>=98%) duplicated content in a genome-wide 465-haplotype alignment of CHM13:

  T     subtelomeric PHR windows          : inside a CHM13 PHR interval, cross-chr duplicated
  P     pericentromeric SD windows        : <=1 Mb of the centromere, outside PHRs/inner genome
  I     interstitial SD windows           : >1 Mb from both telomere and centromere, duplicated
  Utelo subtelomeric UNIQUE windows       : in the 1 Mb telomere flank but not in a PHR, no
                                            cross-chromosome duplicate (the paper's flank control)
  Ugen  genome-wide UNIQUE windows        : unique, sampled away from telomere and centromere

Because every unit is exactly 100 kb, the estimator is length-matched by construction. Windows are
additionally matched (caliper, 1:1, without replacement) on
  * duplicated-sequence identity (bp-weighted mean identity of cross-chromosome alignments),
  * copy number (number of distinct partner chromosomes carrying >=MIN_HAP copies),
  * shared-sequence support (number of merged alignments in the window),
which are the three axes the referee names.

Sequence axis (x): an alignment-block Jaccard computed by ONE code path for every class --
  windows are extracted from CHM13, aligned all-vs-all with wfmash at the same identity
  threshold used for the PHR call (p95), and for each window pair
      J(a,b) = shared_bp / (len_a + len_b - shared_bp)
  with shared_bp = mean(union of a's positions aligned to b, union of b's positions aligned to a).
  This is the continuous counterpart of the odgi graph-node Jaccard the paper uses; on PHR pairs
  it reproduces the published Jaccard ranking (gate G4).

Contact axis (y): the paper's exact contract (gate G2, bit-exact) --
      sum of cooler balanced pixels over the window-pair block / number of non-NaN bin pairs,
  from the HiC-Pro mcool of the sample at 50 kb (primary) and 20 kb, plus HG002 Pore-C and CiFi
  (minimap2/pairtools) as the second and third aligner regimes. An observed/expected variant
  divides by the chromosome-pair baseline computed from all same-chromosome-pair windows, which
  removes the chromosome-territory confound the referee asks about (arm type, chromosome size,
  general chromosome-end contact propensity).

STATISTICS
----------
Per class and per regime: Spearman rho with 10,000-permutation p and a block bootstrap 95% CI
(resampling windows, i.e. whole rows/columns of the pair design, which respects the fact that one
window appears in many pairs). Class contrast: bootstrap CI on delta-rho (T minus matched SD) plus
a Fisher-z test. Confidence stratification: identity tier and repeat content
(RepeatMasker fraction of the window), reported separately rather than pooled.

OUTPUTS
-------
  scripts/hic/matched_sd_control.tsv   tidy results table (one row per class x regime x stratum)
  scripts/hic/matched_sd_control.json  full record: gates, class counts, matching balance, all
                                       statistics, per-class pair counts, software/dataset provenance

Intermediate files (window FASTA, wfmash PAF, pair tables) are written to --workdir (default
/tmp/m5_sd_work) and are reproducible by re-running this script; nothing is written into the repo
except the two result files above.

USAGE
-----
  python3 scripts/hic/matched_sd_control.py gates          # G1-G3 reproduction gates
  python3 scripts/hic/matched_sd_control.py windows        # build + match window classes
  python3 scripts/hic/matched_sd_control.py extract        # samtools faidx window sequences
  python3 scripts/hic/matched_sd_control.py align          # wfmash all-vs-all of the window set
  python3 scripts/hic/matched_sd_control.py jaccard        # PAF -> per-pair Jaccard
  python3 scripts/hic/matched_sd_control.py contact        # mcool -> per-pair contact (all regimes)
  python3 scripts/hic/matched_sd_control.py stats          # matching-aware associations + contrasts
  python3 scripts/hic/matched_sd_control.py all            # everything in order
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
WORK = Path("/tmp/m5_sd_work")

# ---------------------------------------------------------------- data provenance (named files)
WINDOW_TABLE = ROOT / "data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz"
PHR_BED = ROOT / "data/chm13.phrs.bed"
ANNOT_BED = ROOT / "data/chm13-annotations.bed"
CHM13_FAI = Path("/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/CHM13.fa.gz.fai")
CHM13_FA = Path("/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/CHM13.fa.gz")
RM_BED = Path("/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/"
              "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz")
HIC_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/HiC")
MCools = {
    "CHM13_HiC": HIC_DIR / "CHM13/hicpro_output/cool/CHM13.mcool",
    "HG002_HiC": HIC_DIR / "HG002/hicpro_output/cool/HG002.mcool",
    "HG002_PoreC": HIC_DIR / "HG002/porec_output/porec.mcool",
    "HG002_CiFi": HIC_DIR / "HG002/cifi_output/cifi.mcool",
}
PUB_CHM13_PAIRS = ROOT / "data/human_CHM13_hic_50000bp_seqlevel.tsv"
ARM_DIST = ROOT / "data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"

OUT_TSV = ROOT / "scripts/hic/matched_sd_control.tsv"
OUT_JSON = ROOT / "scripts/hic/matched_sd_control.json"

WIN = 100_000
MIN_HAP = 2               # haplotypes needed for a partner chromosome to count as a copy
IDENT_MIN = 0.98          # window-table identity floor (id098)
TELO_BAND = 1_000_000
CENT_BAND = 1_000_000
WFMASH_P = 95             # identity threshold, same as the PHR call (p95)
WFMASH_S = 20000          # segment length
SEED = 20260518


def die(msg: str) -> "NoReturn":
    sys.stderr.write(f"\nFATAL: {msg}\n")
    raise SystemExit(1)


def log(msg: str) -> None:
    print(msg, flush=True)


# =============================================================== gates
def gate_contact_contract(n_check: int = 60, resolution: int = 50000) -> dict:
    """G2: recompute the published hic_contact_norm values bit-exactly from the mcool.

    This pins the contact extraction code path to the published pipeline:
      contact_norm = nansum(balanced submatrix) / n_nonnan_bin_pairs
    """
    import cooler
    import pandas as pd

    df = pd.read_csv(PUB_CHM13_PAIRS, sep="\t")
    clr = cooler.Cooler(f"{MCools['CHM13_HiC']}::resolutions/{resolution}")

    def coords(name: str):
        rest = name.split("#", 2)[2]
        m = re.match(r"[^:]+:(\d+)-(\d+)", rest)
        if not m:
            die(f"cannot parse coordinates from {name}")
        return int(m.group(1)), int(m.group(2))

    exact = mismatch = skipped = 0
    maxdiff = 0.0
    for k in range(min(n_check, len(df))):
        r = df.iloc[k]
        sa, sb = coords(r.seq_a), coords(r.seq_b)
        A = f"{r.chr_a}:{sa[0]}-{sa[1]}"
        B = f"{r.chr_b}:{sb[0]}-{sb[1]}"
        mat = clr.matrix(balance=True).fetch(A, B)
        valid = ~np.isnan(mat)
        if not valid.any():
            skipped += 1
            continue
        rec = float(np.nansum(mat)) / int(valid.sum())
        d = abs(rec - r.hic_contact_norm)
        maxdiff = max(maxdiff, d)
        if d == 0.0:
            exact += 1
        else:
            mismatch += 1
            log(f"  G2 MISMATCH {A} {B}: recomputed {rec!r} published {r.hic_contact_norm!r}")
    out = dict(gate="G2_contact_contract", file=str(PUB_CHM13_PAIRS),
               checked=n_check, exact=exact, mismatch=mismatch, skipped=skipped,
               max_abs_diff=maxdiff, resolution=resolution,
               mcool=str(MCools["CHM13_HiC"]))
    log(f"G2 contact contract: {exact}/{exact + mismatch} bit-exact, max|diff|={maxdiff} "
        f"(skipped {skipped})")
    if mismatch:
        die("G2 failed: contact extraction does not reproduce the published values")
    return out


def gate_published_rho() -> dict:
    """G1: reproduce the published pair-level Spearman (CHM13 Hi-C, 50 kb) = 0.716."""
    import pandas as pd
    from scipy.stats import spearmanr

    df = pd.read_csv(PUB_CHM13_PAIRS, sep="\t")
    rho = spearmanr(df["jaccard"], df["hic_contact_norm"])
    out = dict(gate="G1_published_pair_rho", n_pairs=int(len(df)),
               rho=float(rho.statistic), p=float(rho.pvalue), published_rho=0.716)
    log(f"G1 published pair-level rho: {rho.statistic:.5f} (published 0.716) on {len(df)} pairs")
    return out


def gate_arm_mantel(n_perm: int = 2000) -> dict:
    """G3: arm-level Mantel from the published tables. Documented as a partial reproduction."""
    import itertools
    import pandas as pd
    from scipy.stats import spearmanr

    df = pd.read_csv(PUB_CHM13_PAIRS, sep="\t")
    D = pd.read_csv(ARM_DIST, sep="\t", index_col=0)
    arms = list(D.columns)
    n = len(arms)
    Dm = np.array(D.values, dtype=float)
    contact = df.groupby([df.arm_a, df.arm_b])["hic_contact_norm"].mean()
    cm = {tuple(sorted(k)): float(v) for k, v in contact.items()}
    xs, ys = [], []
    for i, j in itertools.combinations(range(n), 2):
        k = tuple(sorted((arms[i], arms[j])))
        if k in cm and np.isfinite(Dm[i, j]):
            xs.append(1.0 - Dm[i, j])
            ys.append(cm[k])
    obs = float(spearmanr(xs, ys).statistic)
    M = np.where(np.isnan(Dm), np.nan, 1.0 - Dm)
    C = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            k = tuple(sorted((arms[i], arms[j])))
            if k in cm:
                C[i, j] = cm[k]
    iup = np.triu_indices(n, 1)
    mask = np.isfinite(M) & np.isfinite(C)
    mask = mask[iup]
    rng = np.random.default_rng(SEED)
    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        Mp = M[np.ix_(p, p)]
        r = spearmanr(Mp[iup][mask], C[iup][mask]).statistic
        cnt += (r >= obs)
    pval = (cnt + 1) / (n_perm + 1)
    out = dict(gate="G3_arm_mantel_documented", n_arms=n, n_arm_pairs=len(xs),
               rho=obs, p_perm=pval, n_perm=n_perm,
               published_rho=0.656, published_n_arms=38, reproduced=False,
               note=("The report's arm-level Mantel (0.656, n=38 arms) is not reproducible from "
                     "the exported pair table; the paper builds its arm contact matrix inside "
                     "analyze_hic_communities.py from matrix-pixel sums. Variant aggregations "
                     "give 0.725-0.750. M5 therefore uses the pair-level Spearman (exactly "
                     "reproducible, 0.716) as the primary comparand and reports an arm-level "
                     "Mantel computed identically for all classes as a secondary statistic."))
    log(f"G3 arm-level Mantel (documented non-reproduction): rho={obs:.4f} p={pval:.4f} "
        f"n_arms={n} n_pairs={len(xs)} vs published 0.656 (n=38)")
    return out


# =============================================================== windows
def load_chrom_sizes() -> dict:
    sizes = {}
    for line in open(CHM13_FAI):
        c = line.split("\t")
        sizes[c[0]] = int(c[1])
    return sizes


def load_phr() -> list:
    out = []
    for line in open(PHR_BED):
        c = line.rstrip("\n").split("\t")
        out.append((c[0], int(c[1]), int(c[2])))
    return out


def load_centromeres() -> dict:
    cent = defaultdict(list)
    for line in open(ANNOT_BED):
        c = line.rstrip("\n").split("\t")
        if c[3].strip() == "Centromere":
            cent[c[0]].append((int(c[1]), int(c[2])))
    return dict(cent)


def dist_to_interval(s: int, e: int, ivs: list) -> float:
    """Distance in bp from [s,e) to the nearest of ivs (0 if overlapping)."""
    if not ivs:
        return -1
    best = math.inf
    for (a, b) in ivs:
        if s < b and e > a:
            return 0.0
        best = min(best, abs(a - e), abs(s - b))
    return best


def parse_haplotype_list(field: str) -> dict:
    """'chr1-389,chr2-465' -> {'chr1': 389, 'chr2': 465}"""
    d = {}
    if not field:
        return d
    for kv in field.split(","):
        k, _, v = kv.rpartition("-")
        if not k:
            die(f"cannot parse chromosome-haplotype field entry {kv!r}")
        d[k] = int(v)
    return d


def parse_identity_list(field: str) -> dict:
    d = {}
    if not field:
        return d
    for kv in field.split(","):
        k, _, v = kv.rpartition("-")
        if not k:
            die(f"cannot parse chromosome-identity field entry {kv!r}")
        d[k] = float(v)
    return d


def build_windows() -> dict:
    """Classify every CHM13 100 kb window and write the class table."""
    WORK.mkdir(parents=True, exist_ok=True)
    sizes = load_chrom_sizes()
    phr = load_phr()
    cent = load_centromeres()

    rows = []
    n_table = 0
    n_drop_chry = 0
    with gzip.open(WINDOW_TABLE, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        want = ["chrom", "start", "end", "num_alignments", "num_alignments_merged",
                "num_haplotypes", "num_samples", "chroms-num_haplotypes", "chroms-avg_identity"]
        if hdr[:len(want)] != want:
            die(f"window-table header changed: {hdr[:9]}")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            n_table += 1
            chrom = c[0]
            chr_ = chrom.split("#")[-1]
            if chr_ not in sizes:
                n_drop_chry += 1
                continue
            s, e = int(c[1]), int(c[2])
            if e - s != WIN:
                # terminal window of the table grid: anchor it at the chromosome end so that
                # every window in the set is exactly WIN bp (length-matched by construction)
                # while keeping the telomere-proximal content.
                if sizes[chr_] < WIN:
                    continue
                e = e if e - s == WIN else min(sizes[chr_], e)
                s = max(0, e - WIN)
                e = s + WIN
            haps = parse_haplotype_list(c[7])
            idents = parse_identity_list(c[8])
            cross = {k: v for k, v in haps.items() if k != chr_ and v >= MIN_HAP}
            cross_any = {k: v for k, v in haps.items() if k != chr_ and v >= 1}
            shared = len(cross) >= 1
            if cross:
                wsum = sum(cross.values())
                ident_mean = sum(idents[k] * cross[k] for k in cross if k in idents) / max(wsum, 1)
                ident_max = max(idents[k] for k in cross if k in idents)
            else:
                ident_mean = float("nan")
                ident_max = float("nan")
            d_telo = min(s, sizes[chr_] - e)
            d_cen = dist_to_interval(s, e, cent.get(chr_, []))
            in_phr = dist_to_interval(s, e, [(iv[1], iv[2]) for iv in phr if iv[0] == chr_]) == 0.0
            if in_phr and shared:
                cls = "T"
            elif in_phr and not shared:
                cls = "T_nosharing"
            elif d_telo <= TELO_BAND and not in_phr and shared:
                cls = "telo_band_shared_not_PHR"
            elif d_telo <= TELO_BAND and not in_phr and not shared:
                cls = "Utelo"
            elif d_telo > TELO_BAND and d_cen <= CENT_BAND and shared:
                cls = "P"
            elif d_telo > TELO_BAND and d_cen <= CENT_BAND and not shared:
                cls = "cen_band_unique"
            elif shared:
                cls = "I"
            else:
                cls = "Ugen"
            rows.append(dict(
                chrom=chr_, start=s, end=e, cls=cls,
                n_haps=int(c[5]), n_merged=int(c[4]),
                cn=len(cross), cn_any=len(cross_any),
                hap_strong=sum(cross.values()),
                hap_any=sum(cross_any.values()),
                ident_mean=ident_mean, ident_max=ident_max,
                dist_telo=d_telo, dist_cen=d_cen, in_phr=in_phr))

    if not rows:
        die("no usable windows in the window table")
    counts = Counter(r["cls"] for r in rows)
    log(f"window table: {n_table} rows, {n_drop_chry} dropped (not in the masked-Y CHM13 "
        f"fasta, i.e. chrY), {len(rows)} classified")
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        log(f"  {k:26s} {v:6d}")

    with open(WORK / "windows_all.tsv", "w") as f:
        f.write("chrom\tstart\tend\tcls\tn_haps\tn_merged\tcn\tcn_any\thap_strong\thap_any\t"
                "ident_mean\tident_max\tdist_telo\tdist_cen\tin_phr\n")
        for r in rows:
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cls']}\t{r['n_haps']}\t"
                    f"{r['n_merged']}\t{r['cn']}\t{r['cn_any']}\t{r['hap_strong']}\t"
                    f"{r['hap_any']}\t{r['ident_mean']:.6f}\t{r['ident_max']:.6f}\t"
                    f"{int(r['dist_telo'])}\t{int(r['dist_cen'])}\t{int(r['in_phr'])}\n")
    return dict(counts=dict(counts), n_rows=len(rows), n_table=n_table,
                n_drop_not_in_fasta=n_drop_chry, window_bp=WIN,
                min_haplotypes_per_partner=MIN_HAP,
                source=str(WINDOW_TABLE))


# =============================================================== window set + extraction
CLASSES = ["T", "P", "I", "Utelo", "Ucen", "Ugen"]
CLASS_CAP = 470
CLASS_ALIAS = {"cen_band_unique": "Ucen"}
CLASS_LABEL = {
    "T": "subtelomeric PHR windows",
    "P": "pericentromeric SD windows",
    "I": "interstitial SD windows",
    "Utelo": "subtelomeric unique windows (flank control)",
    "Ucen": "pericentromeric unique windows (position-matched control)",
    "Ugen": "genome-wide unique windows",
}


def windows_from_table() -> list:
    rows = []
    with open(WORK / "windows_all.tsv") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = dict(zip(hdr, line.rstrip("\n").split("\t")))
            for k in ("start", "end", "n_haps", "n_merged", "cn", "cn_any", "hap_strong",
                      "hap_any", "dist_telo", "dist_cen"):
                c[k] = int(c[k])
            for k in ("ident_mean", "ident_max"):
                c[k] = float("nan") if c[k] == "nan" else float(c[k])
            rows.append(c)
    return rows


def select_windows(rows: list, cap: int = CLASS_CAP) -> list:
    """Deterministic, seeded subsample per class (all windows of a class if under the cap).

    Sampling strata: chromosome, so that a subsample spans the genome rather than a few arms.
    """
    rng = np.random.default_rng(SEED)
    by_cls = defaultdict(list)

    for r in rows:
        if r["chrom"] == "chrM":
            continue
        cls = CLASS_ALIAS.get(r["cls"], r["cls"])
        if cls in CLASSES:
            by_cls[cls].append(r)
    out = []
    for cls in CLASSES:
        cand = sorted(by_cls[cls], key=lambda r: (r["chrom"], r["start"]))
        if len(cand) <= cap:
            sel = cand
        else:
            # round-robin over chromosomes for even genomic coverage
            per_chr = defaultdict(list)
            for r in cand:
                per_chr[r["chrom"]].append(r)
            for ch in per_chr:
                rng.shuffle(per_chr[ch])
            order = sorted(per_chr)
            rng.shuffle(order)
            sel = []
            while len(sel) < cap:
                progressed = False
                for ch in order:
                    if per_chr[ch] and len(sel) < cap:
                        sel.append(per_chr[ch].pop())
                        progressed = True
                if not progressed:
                    break
            sel.sort(key=lambda r: (r["chrom"], r["start"]))
        for r in sel:
            rec = dict(r)
            rec["cls"] = cls
            out.append(rec)
    return out


def write_window_set(windows: list) -> dict:
    """Name the windows, write FASTA regions + the window->interval map."""
    WORK.mkdir(parents=True, exist_ok=True)
    with open(WORK / "window_set.tsv", "w") as f, open(WORK / "regions.txt", "w") as g:
        f.write("win_id\tchrom\tstart\tend\tcls\tcn\thap_strong\tident_mean\tn_merged\t"
                "dist_telo\tdist_cen\n")
        for i, r in enumerate(windows):
            r["win_id"] = f"CHM13#0#{r['chrom']}_{r['start'] + 1}_{r['end']}"
            f.write(f"{r['win_id']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cls']}\t"
                    f"{r['cn']}\t{r['hap_strong']}\t{r['ident_mean']:.6f}\t{r['n_merged']}\t"
                    f"{r['dist_telo']}\t{r['dist_cen']}\n")
            g.write(f"{r['chrom']}:{r['start'] + 1}-{r['end']}\n")
    counts = Counter(r["cls"] for r in windows)
    log(f"window set: {len(windows)} windows -> {dict(counts)}")
    return dict(n=len(windows), counts=dict(counts), cap=CLASS_CAP, seed=SEED)


def stage_windows() -> dict:
    meta = build_windows()
    rows = windows_from_table()
    sel = select_windows(rows)
    meta["window_set"] = write_window_set(sel)
    return meta


def stage_extract() -> dict:
    """Extract the window sequences from the CHM13 fasta (masked-Y, the Hi-C coordinate space)."""
    if not (WORK / "regions.txt").exists():
        die("run `windows` first")
    fa = WORK / "windows.fa"
    cmd = ["samtools", "faidx", str(CHM13_FA), "-r", str(WORK / "regions.txt"), "-o", str(fa)]
    log(" ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, check=True)
    # verify header order matches the region file (samtools >= 1.16 preserves it)
    regions = [l.strip() for l in open(WORK / "regions.txt")]
    hdrs = []
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                hdrs.append(line[1:].strip())
    if hdrs != regions:
        die(f"sequence order mismatch: {len(hdrs)} headers vs {len(regions)} regions")
    lens = {}
    name = None
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
                lens[name] = 0
            else:
                lens[name] += len(line.strip())
    bad = {k: v for k, v in lens.items() if v != WIN}
    if bad:
        die(f"{len(bad)} windows are not {WIN} bp, e.g. {list(bad.items())[:3]}")
    out = dict(fasta=str(fa), n_sequences=len(hdrs), all_lengths_equal=WIN)
    log(f"extracted {len(hdrs)} sequences, all {WIN} bp: {fa}")
    return out


def stage_align(threads: int = 96) -> dict:
    """wfmash all-vs-all of the pooled window set under two alignment regimes."""
    fa = WORK / "windows.fa"
    if not fa.exists():
        die("run `extract` first")
    # `-f -M` (no filtering, no merging of fragment mappings) is deliberate: the shared-sequence
    # estimator needs *all* paralogous mappings of every segment, not one best mapping per
    # segment. The paper's PHR all-vs-all is an all-vs-ONE run per target sequence, in which a
    # query segment has only one possible target, so no cross-target competition occurs; `-f -M`
    # in all-vs-all reproduces that completeness. Measured effect on a 30-window pilot: default
    # filtering 36 records / 142 kb, `-M` 128 / 318 kb, `-f` 229 / 505 kb, `-f -M` 351 / 732 kb.
    regimes = {
        # regime A: the identity regime of the PHR call's all-vs-all (all-vs-all.p95)
        "p95": ["-p", "95", "-L", "-f", "-M"],
        # regime B: the exact wfmash invocation pggb uses to build the graph that the paper's
        # odgi Jaccard is computed on (-p 90 -s 1000 -l 3000 -n 10 -k 15), i.e. the *filtered*
        # mapping regime with at most 10 mappings per segment. This is the alignment layer the
        # published Jaccard actually rests on, so it is the second alignment regime for the
        # sequence axis.
        "p90_pggb": ["-p", "90", "-s", "1000", "-l", "3000", "-n", "10", "-k", "15", "-L"],
    }
    out = {}
    for name, opts in regimes.items():
        paf = WORK / f"windows.{name}.paf"
        part = WORK / f"windows.{name}.paf.part"
        cmd = ["wfmash", str(fa), *opts, "-t", str(threads), "-B", str(WORK), "--quiet"]
        log(" ".join(shlex.quote(c) for c in cmd))
        with open(part, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
        part.replace(paf)          # atomic: a partial PAF never looks complete
        n = sum(1 for _ in open(paf))
        if n == 0:
            die(f"wfmash regime {name} produced 0 mappings")
        out[name] = dict(paf=str(paf), n_mappings=n, wfmash_opts=opts)
        log(f"  regime {name}: {n} mappings -> {paf}")
    return out


# =============================================================== sequence axis (Jaccard)
CPP_JACCARD = r"""
// Aggregate a PAF (all-vs-all of fixed-size windows) into per-window-pair shared sequence.
//
// For window pair (a,b) the shared sequence is the mean of the base pairs of a covered by
// alignments to b and of b covered by alignments to a; the Jaccard is
//     shared / (len_a + len_b - shared)
// which is the continuous counterpart of an odgi graph-node Jaccard between two sequences.
// Unions are exact interval unions (merged blocks), not bin approximations.
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

struct Rec {
    uint32_t a, b, s, e;
    uint32_t bl, mt;   // alignment block length / matches (identity weight)
    uint8_t side;      // 0 = interval on a, 1 = interval on b
    uint8_t hi;        // 1 if this block identity >= 0.99
    bool operator<(const Rec& o) const {
        if (a != o.a) return a < o.a;
        if (b != o.b) return b < o.b;
        if (side != o.side) return side < o.side;
        if (s != o.s) return s < o.s;
        return e < o.e;
    }
};

struct Union { uint64_t all = 0, hi = 0; };

static Union merge_intervals(const std::vector<Rec>& v) {
    Union u;
    uint32_t cur_s = 0, cur_e = 0; bool open = false;
    uint32_t hi_s = 0, hi_e = 0; bool hi_open = false;
    for (const auto& r : v) {
        if (!open) { cur_s = r.s; cur_e = r.e; open = true; }
        else if (r.s <= cur_e) { if (r.e > cur_e) cur_e = r.e; }
        else { u.all += cur_e - cur_s; cur_s = r.s; cur_e = r.e; }
        if (r.hi) {
            if (!hi_open) { hi_s = r.s; hi_e = r.e; hi_open = true; }
            else if (r.s <= hi_e) { if (r.e > hi_e) hi_e = r.e; }
            else { u.hi += hi_e - hi_s; hi_s = r.s; hi_e = r.e; }
        }
    }
    if (open) u.all += cur_e - cur_s;
    if (hi_open) u.hi += hi_e - hi_s;
    return u;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s <paf> <fai> <out.tsv>\n", argv[0]);
        return 2;
    }
    // names + lengths
    std::unordered_map<std::string, uint32_t> id;
    std::vector<std::string> names;
    std::vector<uint32_t> lens;
    {
        FILE* f = fopen(argv[2], "r");
        if (!f) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
        char line[8192];
        while (fgets(line, sizeof line, f)) {
            char nm[4096]; unsigned long long len;
            if (sscanf(line, "%4095s\t%llu", nm, &len) != 2) continue;
            id.emplace(std::string(nm), (uint32_t)names.size());
            names.push_back(nm);
            lens.push_back((uint32_t)len);
        }
        fclose(f);
    }
    std::vector<Rec> recs;
    {
        FILE* f = fopen(argv[1], "r");
        if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
        char *line = nullptr; size_t cap = 0; ssize_t n;
        char qn[4096], tn[4096];
        unsigned long long qlen, qs, qe, tlen, ts, te, matches, blocklen;
        while ((n = getline(&line, &cap, f)) > 0) {
            int got = sscanf(line, "%4095s\t%llu\t%llu\t%llu\t%*c\t%4095s\t%llu\t%llu\t%llu\t%llu\t%llu",
                             qn, &qlen, &qs, &qe, tn, &tlen, &ts, &te, &matches, &blocklen);
            if (got != 10) { fprintf(stderr, "malformed PAF line: %.120s\n", line); return 1; }
            auto ia = id.find(std::string(qn));
            auto ib = id.find(std::string(tn));
            if (ia == id.end() || ib == id.end()) {
                fprintf(stderr, "PAF name not in fai: %s / %s\n", qn, tn);
                return 1;
            }
            uint32_t a = ia->second, b = ib->second;
            if (a == b) continue;                       // skip self mappings
            if (a > b) { std::swap(a, b); std::swap(qs, ts); std::swap(qe, te); }
            uint8_t hi = (blocklen > 0 && (double)matches / (double)blocklen >= 0.99) ? 1 : 0;
            recs.push_back({a, b, (uint32_t)qs, (uint32_t)qe, (uint32_t)blocklen,
                            (uint32_t)matches, 0, hi});
            recs.push_back({a, b, (uint32_t)ts, (uint32_t)te, (uint32_t)blocklen,
                            (uint32_t)matches, 1, hi});
        }
        free(line);
        fclose(f);
    }
    std::sort(recs.begin(), recs.end());
    FILE* out = fopen(argv[3], "w");
    if (!out) { fprintf(stderr, "cannot write %s\n", argv[3]); return 1; }
    fprintf(out, "a_id\tb_id\tshared_a\tshared_b\tshared_a_hi\tshared_b_hi\tn_blocks\t"
                 "sum_blocklen\tsum_matches\n");
    size_t i = 0;
    while (i < recs.size()) {
        uint32_t a = recs[i].a, b = recs[i].b;
        std::vector<Rec> sideA, sideB;
        size_t j = i;
        while (j < recs.size() && recs[j].a == a && recs[j].b == b) {
            (recs[j].side == 0 ? sideA : sideB).push_back(recs[j]);
            ++j;
        }
        Union ua = merge_intervals(sideA), ub = merge_intervals(sideB);
        // n_blocks / identity: count query-side records only
        uint64_t nb = sideA.size(), sb = 0, sm = 0;
        for (const auto& r : sideA) { sb += r.bl; sm += r.mt; }
        fprintf(out, "%u\t%u\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                a, b, (unsigned long long)ua.all, (unsigned long long)ub.all,
                (unsigned long long)ua.hi, (unsigned long long)ub.hi,
                (unsigned long long)nb, (unsigned long long)sb, (unsigned long long)sm);
        i = j;
    }
    fclose(out);
    return 0;
}
"""


def compile_jaccard(work: Path = WORK) -> Path:
    src = work / "paf_jaccard.cpp"
    exe = work / "paf_jaccard"
    src.write_text(CPP_JACCARD)
    cmd = ["g++", "-O3", "-std=c++17", "-o", str(exe), str(src)]
    log(" ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, check=True)
    return exe


def stage_jaccard(regime: str = "p95") -> dict:
    """PAF -> per-window-pair shared sequence (identical code path for every class)."""
    paf = WORK / f"windows.{regime}.paf"
    fai = WORK / "windows.fa.fai"
    if not paf.exists():
        die(f"missing {paf}; run `align` first")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(WORK / "windows.fa")], check=True)
    exe = compile_jaccard()
    out = WORK / f"pairs.{regime}.tsv"
    subprocess.run([str(exe), str(paf), str(fai), str(out)], check=True)

    # The FASTA sequence names are '#'-free samtools region names (wfmash groups mappings by the
    # last '#' prefix and then skips self-group mappings, so PanSN-style names would suppress the
    # paralogous mappings we need). Map them back to win_id.
    region2win = {}
    with open(WORK / "window_set.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            c = dict(zip(h, line.rstrip("\n").split("\t")))
            region2win[f"{c['chrom']}:{int(c['start']) + 1}-{c['end']}"] = c["win_id"]
    fai_names = [l.split("\t")[0] for l in open(fai)]
    missing = [n for n in fai_names if n not in region2win]
    if missing:
        die(f"{len(missing)} FASTA names have no win_id mapping, e.g. {missing[:3]}")
    fai_names = [region2win[n] for n in fai_names]
    fai_lens = [int(l.split("\t")[1]) for l in open(fai)]
    fin = WORK / f"pairs.{regime}.labelled.tsv"
    n_lab = 0
    with open(out) as f, open(fin, "w") as g:
        g.write("win_a\twin_b\tlen_a\tlen_b\tshared\tjaccard\tidentity\tn_blocks\t"
                "jaccard_hi\n")
        hdr = f.readline()
        if not hdr.startswith("a_id"):
            die("unexpected jaccard output header")
        for line in f:
            c = line.rstrip("\n").split("\t")
            a, b = int(c[0]), int(c[1])
            ua, ub = int(c[2]), int(c[3])
            ua_hi, ub_hi = int(c[4]), int(c[5])
            nb, sb, sm = int(c[6]), int(c[7]), int(c[8])
            la, lb = fai_lens[a], fai_lens[b]
            shared = 0.5 * (ua + ub)
            jac = shared / (la + lb - shared) if (la + lb - shared) > 0 else 0.0
            shared_hi = 0.5 * (ua_hi + ub_hi)
            jac_hi = shared_hi / (la + lb - shared_hi) if (la + lb - shared_hi) > 0 else 0.0
            ident = sm / sb if sb else float("nan")
            g.write(f"{fai_names[a]}\t{fai_names[b]}\t{la}\t{lb}\t{shared:.1f}\t{jac:.9f}\t"
                    f"{ident:.6f}\t{nb}\t{jac_hi:.9f}\n")
            n_lab += 1
    log(f"jaccard[{regime}]: {n_lab} window pairs -> {fin}")
    return dict(paf=str(paf), labelled=str(fin), n_pairs=n_lab)


# =============================================================== contact axis
REGIMES = [
    # (label, mcool, resolution, coordinate space, mapper/fragment pipeline)
    ("CHM13_HiC_50kb", "CHM13_HiC", 50000, "CHM13", "HiC-Pro + bowtie2 (MIN_MAPQ=0, RM_MULTI=0)"),
    ("CHM13_HiC_20kb", "CHM13_HiC", 20000, "CHM13", "HiC-Pro + bowtie2 (MIN_MAPQ=0, RM_MULTI=0)"),
    ("HG002_HiC_50kb", "HG002_HiC", 50000, "HG002", "HiC-Pro + bowtie2 (MIN_MAPQ=0, RM_MULTI=0)"),
    ("HG002_PoreC_50kb", "HG002_PoreC", 50000, "HG002", "minimap2 + pairtools (no MAPQ filter)"),
    ("HG002_CiFi_50kb", "HG002_CiFi", 50000, "HG002", "minimap2 + pairtools (no MAPQ filter)"),
]


def load_window_set() -> dict:
    rows = {}
    with open(WORK / "window_set.tsv") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = dict(zip(hdr, line.rstrip("\n").split("\t")))
            for k in ("start", "end", "cn", "hap_strong", "n_merged", "dist_telo", "dist_cen"):
                c[k] = int(c[k])
            c["ident_mean"] = float(c["ident_mean"])
            rows[c["win_id"]] = c
    return rows


def enumerate_pairs(windows: dict) -> list:
    """All inter-chromosomal window pairs, per class (the published pair-set construction)."""
    by_cls = defaultdict(list)
    for w in windows.values():
        by_cls[w["cls"]].append(w["win_id"])
    pairs = []
    for cls, ids in by_cls.items():
        ids = sorted(ids)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = ids[i], ids[j]
                if windows[a]["chrom"] == windows[b]["chrom"]:
                    continue
                pairs.append((cls, a, b))
    log("pairs per class: " +
        str({c: sum(1 for p in pairs if p[0] == c) for c in sorted(by_cls)}))
    return pairs


# ---------------------------------------------------------------- HG002 liftover
def stage_liftover(threads: int = 96) -> dict:
    """Map the CHM13 windows onto HG002's haplotype-resolved chromosomes (for the second and
    third contact-aligner regimes, whose mcools live in HG002 coordinate space).

    One wfmash run, query = our CHM13 windows, target = the HG002 hap fasta that the HG002
    Hi-C/Pore-C/CiFi mcools were built from (chrN_PATERNAL / chrN_MATERNAL naming).
    """
    fa = WORK / "windows.fa"
    if not fa.exists():
        die("run `extract` first")
    hg = HIC_DIR / "HG002/HG002.fa.gz"
    if not hg.exists():
        die(f"missing HG002 fasta {hg}")
    paf = WORK / "windows_vs_HG002.paf"
    cmd = ["wfmash", str(hg), str(fa), "-p", "95", "-t", str(threads), "-B", str(WORK),
           "--quiet"]
    log(" ".join(shlex.quote(c) for c in cmd))
    with open(paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)
    n = sum(1 for _ in open(paf))
    if n == 0:
        die("HG002 liftover produced 0 mappings")
    # best interval per (window, haplotype chromosome): span of alignments, requiring that the
    # covered query fraction be >= 20% and that the target's base chromosome match the window's
    fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(WORK / "windows.fa.fai")}
    wins = load_window_set()
    best = {}
    cover = defaultdict(int)
    with open(paf) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            qn, qs, qe = c[0], int(c[2]), int(c[3])
            tn, ts, te = c[5], int(c[7]), int(c[8])
            if qn not in wins:
                die(f"unexpected query {qn} in liftover PAF")
            w = wins[qn]
            base = tn.split("_")[0]
            if base != w["chrom"]:
                continue
            key = (qn, tn)
            if key not in best:
                best[key] = [ts, te, qs, qe]
            else:
                best[key][0] = min(best[key][0], ts)
                best[key][1] = max(best[key][1], te)
                best[key][2] += (qe - qs)
            cover[key] += (qe - qs)
    out = WORK / "hg002_liftover.tsv"
    n_out = 0
    with open(out, "w") as g:
        g.write("win_id\thap_chrom\tstart\tend\tcovered_bp\n")
        for (qn, tn), v in sorted(best.items()):
            w = wins[qn]
            cov = min(cover[(qn, tn)], fai[qn])
            if cov < 0.2 * fai[qn]:
                continue
            g.write(f"{qn}\t{tn}\t{v[0]}\t{v[1]}\t{cov}\n")
            n_out += 1
    log(f"liftover: {n} mappings -> {n_out} window/haplotype-chromosome intervals ({out})")
    return dict(paf=str(paf), table=str(out), n_mappings=n, n_intervals=n_out,
                hg002_fasta=str(hg))


# ---------------------------------------------------------------- contact extraction
_worker_state = {}


def _init_worker(mcool_key: str, resolution: int):
    import cooler
    _worker_state["clr"] = cooler.Cooler(f"{MCools[mcool_key]}::resolutions/{resolution}")


def _fetch(task):
    import numpy as np
    clr = _worker_state["clr"]
    cls, a, b, ra, rb = task
    out = []
    for reg_a, reg_b in zip(ra, rb):
        try:
            mat = clr.matrix(balance=True).fetch(reg_a, reg_b)
        except Exception as exc:  # out-of-range or unknown chromosome
            die(f"cooler fetch failed for {reg_a} x {reg_b}: {exc}")
        valid = ~np.isnan(mat)
        if valid.any():
            s = float(np.nansum(mat))
            out.append((s, int(valid.sum())))
        else:
            out.append((0.0, 0))
    return (cls, a, b, out)


def stage_contact(regimes=None, procs: int = 48, liftover: bool = True) -> dict:
    import multiprocessing as mp

    if regimes is None:
        regimes = REGIMES
    windows = load_window_set()
    pairs = enumerate_pairs(windows)

    # CHM13 region strings
    def chm13_region(w):
        return f"{w['chrom']}:{w['start']}-{w['end']}"

    # HG002 liftover table
    hg_int = defaultdict(list)          # win_id -> [(hap_chrom, start, end)]
    hg_path = WORK / "hg002_liftover.tsv"
    if liftover:
        if not hg_path.exists():
            die("missing hg002_liftover.tsv; run `liftover` first")
        for line in open(hg_path):
            c = line.rstrip("\n").split("\t")
            hg_int[c[0]].append((c[1], int(c[2]), int(c[3])))

    meta = {}
    for label, mcool_key, res, space, mapper in regimes:
        out = WORK / f"contact.{label}.tsv"
        if space == "CHM13":
            tasks = [(cls, a, b, [chm13_region(windows[a])], [chm13_region(windows[b])])
                     for cls, a, b in pairs]
        else:
            tasks = []
            for cls, a, b in pairs:
                if a not in hg_int or b not in hg_int:
                    continue
                ia, ib = hg_int[a], hg_int[b]
                ra = [f"{h}:{s}-{e}" for h, s, e in ia]
                rb = [f"{h}:{s}-{e}" for h, s, e in ib]
                tasks.append((cls, a, b, ra, rb))
        if not tasks:
            die(f"no pairs for regime {label}")
        log(f"contact[{label}]: {len(tasks)} pairs x up to {max(len(t[3]) for t in tasks)}"
            f" haplotype intervals, {procs} workers, resolution {res}")
        rows = []
        with mp.Pool(processes=procs, initializer=_init_worker,
                     initargs=(mcool_key, res)) as pool:
            for k, (cls, a, b, res_list) in enumerate(
                    pool.imap_unordered(_fetch, tasks, chunksize=64)):
                tot = sum(s for s, _ in res_list)
                nb = sum(n for _, n in res_list)
                if nb == 0:
                    continue                      # all bins ICE-masked: missing, not zero
                rows.append((cls, a, b, tot, nb, tot / nb, len(res_list)))
                if (k + 1) % 50000 == 0:
                    log(f"  {k + 1}/{len(tasks)}")
        with open(out, "w") as f:
            f.write("cls\twin_a\twin_b\tcontact_sum\tn_bins\tcontact_norm\tn_regions\n")
            for r in rows:
                f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.10g}\t{r[4]}\t{r[5]:.10g}\t{r[6]}\n")
        meta[label] = dict(contact=str(out), n_pairs=len(rows), n_tasks=len(tasks),
                           resolution=res, space=space, mapper=mapper,
                           mcool=str(MCools[mcool_key]))
        log(f"  wrote {len(rows)} pairs -> {out}")
    return meta


# =============================================================== statistics
def load_rm_fraction() -> dict:
    """Per-window RepeatMasker footprint (all classes, and Satellite class separately)."""
    if not RM_BED.exists():
        die(f"missing RepeatMasker bed {RM_BED}")
    windows = load_window_set()
    idx = defaultdict(list)                     # chrom -> [(start, end, is_sat)]
    cov = {w: [0, 0] for w in windows}
    # windows are 100 kb; bucket them by chromosome for a sweep
    by_chrom = defaultdict(list)
    for w, r in windows.items():
        by_chrom[r["chrom"]].append((r["start"], r["end"], w))
    for ch in by_chrom:
        by_chrom[ch].sort()
    with gzip.open(RM_BED, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            chrom = c[0].split("#")[-1]
            if chrom not in by_chrom:
                continue
            s, e = int(c[1]), int(c[2])
            sat = 1 if (c[6] == "Satellite" or c[6] == "Simple_repeat") else 0
            for (ws, we, w) in by_chrom[chrom]:
                if we <= s:
                    continue
                if ws >= e:
                    break
                ov = max(0, min(e, we) - max(s, ws))
                if ov:
                    idx[(w, sat)].append(ov)
    out = {}
    for w, r in windows.items():
        tot = sum(idx[(w, 0)]) + sum(idx[(w, 1)])
        sat = sum(idx[(w, 1)])
        L = r["end"] - r["start"]
        out[w] = dict(repeat_frac=min(1.0, tot / L), satellite_frac=min(1.0, sat / L))
    return out


def _ranks(x: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    r = rankdata(x, method="average")
    r = r - r.mean()
    sd = r.std()
    return r / sd if sd > 0 else np.zeros_like(r)


def _rho(a_std: np.ndarray, b_std: np.ndarray) -> float:
    if len(a_std) < 3:
        return float("nan")
    return float(np.dot(a_std, b_std) / len(a_std))


def assoc(jac, con, chr_a, chr_b, n_perm=5000, n_boot=1000, seed=SEED) -> dict:
    """Spearman rho with permutation p and a chromosome-block bootstrap 95% CI.

    The block bootstrap resamples chromosomes (not pairs), so it keeps the dependence induced by
    the fact that one window/arm occurs in many pairs -- the objection both referees raise.
    """
    jac = np.asarray(jac, dtype=float)
    con = np.asarray(con, dtype=float)
    chr_a = np.asarray(chr_a)
    chr_b = np.asarray(chr_b)
    if len(jac) < 10:
        return dict(rho=float("nan"), p_perm=float("nan"), ci_lo=float("nan"),
                    ci_hi=float("nan"), n=len(jac))
    rj = _ranks(jac)
    rc = _ranks(con)
    obs = _rho(rj, rc)
    rng = np.random.default_rng(seed)
    # permutation null on the contact ranks
    cnt = 0
    null = np.empty(n_perm)
    # centre for a faster dot-product null: rj.sum()==0 and rc.sum()==0
    for k in range(n_perm):
        p = rng.permutation(len(rc))
        r = float(np.dot(rj, rc[p]) / len(rc))
        null[k] = r
        cnt += (r >= obs)
    p_perm = (cnt + 1) / (n_perm + 1)
    # chromosome-block bootstrap
    chars = sorted(set(chr_a.tolist()) | set(chr_b.tolist()))
    order = {c: i for i, c in enumerate(chars)}
    ia = np.array([order[c] for c in chr_a])
    ib = np.array([order[c] for c in chr_b])
    lo = np.minimum(ia, ib)
    hi = np.maximum(ia, ib)
    blocks = defaultdict(list)
    for k in range(len(jac)):
        blocks[(lo[k], hi[k])].append(k)
    blocks = {k: np.array(v) for k, v in blocks.items()}
    boot = np.empty(n_boot)
    m = len(chars)
    for b in range(n_boot):
        draw = rng.integers(0, m, m)
        idx_parts = []
        for x in draw:
            for y in draw:
                key = (min(x, y), max(x, y))
                if key in blocks:
                    idx_parts.append(blocks[key])
        if not idx_parts:
            boot[b] = np.nan
            continue
        idx = np.concatenate(idx_parts)
        boot[b] = _rho(rj[idx], rc[idx])
    ok = boot[np.isfinite(boot)]
    ci_lo, ci_hi = (float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))) \
        if len(ok) > 10 else (float("nan"), float("nan"))
    return dict(rho=obs, p_perm=float(p_perm), ci_lo=ci_lo, ci_hi=ci_hi, n=int(len(jac)),
                null_mean=float(null.mean()), null_sd=float(null.std()))


def caliper_match(source: list, candidates: list, wins: dict, variables: list,
                  caliper: float = 0.5) -> tuple:
    """1:1 greedy nearest-neighbour matching without replacement on standardized covariates."""
    def vec(w):
        v = []
        for var in variables:
            x = wins[w][var]
            if var == "cn":
                v.append(math.log1p(x))
            elif var == "hap_strong":
                v.append(math.log1p(x))
            else:
                v.append(x)
        return v

    pool = source + candidates
    X = np.array([vec(w) for w in pool], dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    s_idx = list(range(len(source)))
    c_idx = list(range(len(source), len(pool)))
    # most-constrained first: fewest candidates within caliper
    def n_within(i):
        d = np.linalg.norm(Xs[c_idx] - Xs[i], axis=1)
        return int((d <= caliper).sum())
    s_order = sorted(s_idx, key=n_within)
    used = set()
    pairs = []
    for i in s_order:
        best, bd = None, math.inf
        for j in c_idx:
            if j in used:
                continue
            d = float(np.linalg.norm(Xs[i] - Xs[j]))
            if d < bd:
                best, bd = j, d
        if best is None or bd > caliper:
            continue
        used.add(best)
        pairs.append((pool[i], pool[best], bd))
    matched_src = [a for a, _, _ in pairs]
    matched_cand = [b for _, b, _ in pairs]
    return matched_src, matched_cand, pairs


# ---------------------------------------------------------------- joining
def _key(a: str, b: str) -> tuple:
    return (a, b) if a <= b else (b, a)


def load_jaccard(seq_regime: str) -> dict:
    path = WORK / f"pairs.{seq_regime}.labelled.tsv"
    if not path.exists():
        die(f"missing {path}; run `jaccard --regime {seq_regime}` first")
    out = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            out[_key(c[idx["win_a"]], c[idx["win_b"]])] = dict(
                jaccard=float(c[idx["jaccard"]]),
                jaccard_hi=float(c[idx["jaccard_hi"]]),
                identity=float(c[idx["identity"]]) if c[idx["identity"]] != "nan"
                else float("nan"),
                n_blocks=int(c[idx["n_blocks"]]))
    return out


def load_contact(label: str) -> dict:
    path = WORK / f"contact.{label}.tsv"
    if not path.exists():
        die(f"missing {path}; run `contact` first")
    out = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            out[_key(c[idx["win_a"]], c[idx["win_b"]])] = float(c[idx["contact_norm"]])
    return out


def build_design(windows: dict, jac: dict, con: dict) -> dict:
    """Assemble per-class arrays; pairs missing from the contact table are dropped (ICE-masked),
    pairs missing from the alignment table enter with Jaccard 0 (as in the published analysis)."""
    pairs = enumerate_pairs(windows)
    design = defaultdict(lambda: defaultdict(list))
    n_zero = Counter()
    n_dropped = Counter()
    for cls, a, b in pairs:
        k = _key(a, b)
        if k not in con:
            n_dropped[cls] += 1
            continue
        j = jac.get(k)
        jv = 0.0
        jh = 0.0
        ident = float("nan")
        nblk = 0
        if j is None:
            n_zero[cls] += 1
        else:
            jv = j["jaccard"]
            jh = j["jaccard_hi"]
            ident = j["identity"]
            nblk = j["n_blocks"]
            if jv == 0.0:
                n_zero[cls] += 1
        d = design[cls]
        d["jac"].append(jv)
        d["jac_hi"].append(jh)
        d["con"].append(con[k])
        d["ident"].append(ident)
        d["nblk"].append(nblk)
        d["win_a"].append(a)
        d["win_b"].append(b)
        d["chr_a"].append(windows[a]["chrom"])
        d["chr_b"].append(windows[b]["chrom"])
        d["rm_a"].append(None)      # filled by caller
    return dict(design=design, n_zero=dict(n_zero), n_dropped=dict(n_dropped),
                n_pairs_total=len(pairs))


def add_oe(design: dict) -> None:
    """Observed/expected contact per class: divide by the chromosome-pair baseline computed over
    ALL pairs in the design (the referee's 'general chromosome-end contact propensity' and
    'chromosome size' control)."""
    pooled = defaultdict(list)
    for cls, d in design.items():
        for k in range(len(d["con"])):
            pooled[(min(d["chr_a"][k], d["chr_b"][k]), max(d["chr_a"][k], d["chr_b"][k]))] \
                .append(d["con"][k])
    base = {k: float(np.mean(v)) for k, v in pooled.items()}
    for cls, d in design.items():
        oe = []
        for k in range(len(d["con"])):
            key = (min(d["chr_a"][k], d["chr_b"][k]), max(d["chr_a"][k], d["chr_b"][k]))
            b = base.get(key, 0.0)
            oe.append(d["con"][k] / b if b > 0 else 0.0)
        d["con_oe"] = oe
        d["baseline_n"] = len(pooled)


def class_mantel(d: dict) -> dict:
    """Chromosome-level Mantel analogue of the paper's arm-level Mantel (secondary statistic)."""
    pa = defaultdict(list)
    ca = defaultdict(list)
    for k in range(len(d["con"])):
        key = (min(d["chr_a"][k], d["chr_b"][k]), max(d["chr_a"][k], d["chr_b"][k]))
        pa[key].append(d["jac"][k])
        ca[key].append(d["con"][k])
    keys = sorted(pa)
    if len(keys) < 10:
        return dict(rho=float("nan"), p_perm=float("nan"), n_cell=len(keys))
    x = np.array([float(np.mean(pa[k])) for k in keys])
    y = np.array([float(np.mean(ca[k])) for k in keys])
    return assoc(x, y, [k[0] for k in keys], [k[1] for k in keys], n_perm=2000, n_boot=500)


def fisher_z_diff(rho1: float, n1: int, rho2: float, n2: int) -> dict:
    if not (np.isfinite(rho1) and np.isfinite(rho2)) or n1 < 4 or n2 < 4:
        return dict(z=float("nan"), p=float("nan"))
    z1, z2 = math.atanh(max(min(rho1, 0.999999), -0.999999)), \
        math.atanh(max(min(rho2, 0.999999), -0.999999))
    se = math.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z = (z1 - z2) / se
    from scipy.stats import norm
    return dict(z=float(z), p=float(2 * norm.sf(abs(z))))


# ---------------------------------------------------------------- stage: statistics
ROW_COLS = ["analysis", "class", "class_label", "seq_regime", "contact_regime", "metric",
            "stratum", "n_pairs", "n_windows", "n_jaccard_zero", "rho", "p_perm",
            "ci_lo", "ci_hi", "null_mean", "n_perm", "n_boot", "delta_vs",
            "delta_ci_lo", "delta_ci_hi", "delta_p_boot", "fisher_z_p",
            "matched_n_windows", "note"]


def _row(**kw) -> dict:
    out = {c: "" for c in ROW_COLS}
    out.update({k: v for k, v in kw.items() if k in out})
    return out


def _fmt(v):
    if isinstance(v, float):
        if not np.isfinite(v):
            return "nan"
        return f"{v:.6g}"
    return str(v)


def _subset(d: dict, keep: set) -> dict:
    out = {}
    for k, v in d.items():
        out[k] = [x for x, w in zip(v, d["win_a"]) if w in keep] if k != "win_a" else \
            [w for w in v if w in keep]
    if "win_a" in d:
        out["win_b"] = [w for w in d["win_b"] if w in keep]
    return out


def _windows_in(d: dict) -> list:
    return sorted(set(d["win_a"]) | set(d["win_b"]))


def _blocks_of(chr_a, chr_b, chrom_index):
    lo = np.array([min(chrom_index[a], chrom_index[b]) for a, b in zip(chr_a, chr_b)])
    hi = np.array([max(chrom_index[a], chrom_index[b]) for a, b in zip(chr_a, chr_b)])
    blocks = defaultdict(list)
    for k in range(len(lo)):
        blocks[(lo[k], hi[k])].append(k)
    return {k: np.array(v) for k, v in blocks.items()}


def delta_rho_bootstrap(d1: dict, d2: dict, metric: str, n_boot: int = 1000,
                        seed: int = SEED) -> dict:
    """Paired chromosome-block bootstrap of rho(d1) - rho(d2)."""
    chroms = sorted(set(d1["chr_a"]) | set(d1["chr_b"]) | set(d2["chr_a"]) | set(d2["chr_b"]))
    ci = {c: i for i, c in enumerate(chroms)}
    b1 = _blocks_of(d1["chr_a"], d1["chr_b"], ci)
    b2 = _blocks_of(d2["chr_a"], d2["chr_b"], ci)
    rj1, rc1 = _ranks(np.array(d1["jac"])), _ranks(np.array(d1[metric]))
    rj2, rc2 = _ranks(np.array(d2["jac"])), _ranks(np.array(d2[metric]))
    rng = np.random.default_rng(seed)
    m = len(chroms)
    deltas = []
    for _ in range(n_boot):
        draw = rng.integers(0, m, m)
        i1, i2 = [], []
        for x in draw:
            for y in draw:
                key = (min(x, y), max(x, y))
                if key in b1:
                    i1.append(b1[key])
                if key in b2:
                    i2.append(b2[key])
        if not i1 or not i2:
            continue
        a1 = np.concatenate(i1)
        a2 = np.concatenate(i2)
        r1 = _rho(rj1[a1], rc1[a1])
        r2 = _rho(rj2[a2], rc2[a2])
        if np.isfinite(r1) and np.isfinite(r2):
            deltas.append(r1 - r2)
    deltas = np.array(deltas)
    if len(deltas) < 10:
        return dict(delta=float("nan"), lo=float("nan"), hi=float("nan"),
                    p=float("nan"), n_boot=len(deltas))
    obs = _rho(rj1, rc1) - _rho(rj2, rc2)
    return dict(delta=float(obs), lo=float(np.percentile(deltas, 2.5)),
                hi=float(np.percentile(deltas, 97.5)),
                p=float(2 * min((deltas <= 0).mean(), (deltas >= 0).mean())),
                n_boot=int(len(deltas)))


def stage_stats(seq_regimes=("p95", "p90_pggb"),
                contact_labels=("CHM13_HiC_50kb", "CHM13_HiC_20kb", "HG002_HiC_50kb",
                                "HG002_PoreC_50kb", "HG002_CiFi_50kb"),
                match_vars=("ident_mean", "cn"), n_perm: int = 5000,
                n_boot: int = 1000, primary: str = "CHM13_HiC_50kb") -> dict:
    windows = load_window_set()
    rm = load_rm_fraction()
    jac_cache = {sr: load_jaccard(sr) for sr in seq_regimes if (WORK / f"pairs.{sr}.labelled.tsv").exists()}
    if not jac_cache:
        die("no alignment pair tables found; run `jaccard`")
    rows = []
    detail = {"matching": {}, "classes": {}, "contact_regimes": list(jac_cache and contact_labels)}

    # ---------- matching (identical for every regime, so computed once) ----------
    by_cls = defaultdict(list)
    for w, r in windows.items():
        by_cls[r["cls"]].append(w)
    for w in windows:
        windows[w]["repeat_frac"] = rm[w]["repeat_frac"]
        windows[w]["satellite_frac"] = rm[w]["satellite_frac"]
    t_wins = sorted(by_cls["T"])
    match_sets = {}
    for target_cls, cands in (("P", sorted(by_cls["P"])), ("I", sorted(by_cls["I"])),
                              ("SD", sorted(by_cls["P"]) + sorted(by_cls["I"]))):
        ms, mc, pairs = caliper_match(t_wins, cands, windows, list(match_vars))
        match_sets[target_cls] = (set(ms), set(mc))
        bal = {}
        for var in ("ident_mean", "cn", "hap_strong", "repeat_frac"):
            def stats_for(ids):
                v = np.array([windows[w][var] if var != "ident_mean" else
                              (0.0 if not np.isfinite(windows[w][var]) else windows[w][var])
                              for w in ids], dtype=float)
                return float(v.mean()) if len(v) else float("nan"), \
                    float(v.std()) if len(v) else float("nan")
            t_mu, t_sd = stats_for(sorted(ms))
            c_mu, c_sd = stats_for(sorted(mc))
            t0, _ = stats_for(t_wins)
            c0, _ = stats_for(cands)
            pooled_sd = math.sqrt((t_sd ** 2 + c_sd ** 2) / 2) if np.isfinite(t_sd) else float("nan")
            bal[var] = dict(T_mean_matched=t_mu, ctrl_mean_matched=c_mu,
                            T_mean_all=t0, ctrl_mean_all=c0,
                            smd_matched=(t_mu - c_mu) / pooled_sd if pooled_sd else float("nan"),
                            smd_unmatched=(t0 - c0) / pooled_sd if pooled_sd else float("nan"))
        detail["matching"][target_cls] = dict(n_matched_pairs=len(pairs), balance=bal,
                                              n_T_matched=len(ms), n_ctrl_matched=len(mc),
                                              n_T_total=len(t_wins), n_ctrl_total=len(cands),
                                              caliper_sd=0.5, variables=list(match_vars))
        log(f"matching T~{target_cls}: {len(ms)} T windows <-> {len(mc)} control windows "
            f"(of {len(t_wins)} / {len(cands)})")
        for var, b in bal.items():
            log(f"   {var:13s} T {b['T_mean_all']:.4g}->{b['T_mean_matched']:.4g}  "
                f"ctrl {b['ctrl_mean_all']:.4g}->{b['ctrl_mean_matched']:.4g}  "
                f"SMD {b['smd_unmatched']:+.3f} -> {b['smd_matched']:+.3f}")

    # ---------- associations per contact regime x sequence regime x class ----------
    for label in contact_labels:
        if not (WORK / f"contact.{label}.tsv").exists():
            log(f"  (skipping contact regime {label}: no table)")
            continue
        con = load_contact(label)
        for sr, jac in jac_cache.items():
            dd = build_design(windows, jac, con)
            design = dd["design"]
            add_oe(design)
            for cls in sorted(design):
                d = design[cls]
                for metric, mlabel in (("con", "balanced_contact"), ("con_oe", "contact_O_over_E")):
                    a = assoc(d["jac"], d[metric], d["chr_a"], d["chr_b"],
                              n_perm=n_perm, n_boot=n_boot)
                    # repeat-methylation-style confidence stratum: RepeatMasker footprint
                    rmvec = np.array([
                        (windows[d["win_a"][k]]["repeat_frac"] +
                         windows[d["win_b"][k]]["repeat_frac"]) / 2
                        for k in range(len(d["jac"]))])
                    rows.append(_row(
                        analysis="class_association", **{"class": cls},
                        class_label=CLASS_LABEL.get(cls, cls), seq_regime=sr,
                        contact_regime=label, metric=mlabel, stratum="all",
                        n_pairs=a["n"], n_windows=len(_windows_in(d)),
                        n_jaccard_zero=dd["n_zero"].get(cls, 0),
                        rho=a["rho"], p_perm=a["p_perm"], ci_lo=a["ci_lo"], ci_hi=a["ci_hi"],
                        null_mean=a.get("null_mean", float("nan")), n_perm=n_perm,
                        n_boot=n_boot,
                        note=f"n_zero_contact_dropped={dd['n_dropped'].get(cls, 0)};"
                             f"median_repeat_frac={float(np.median(rmvec)):.3g}"))
                    # strata (primary contact regime only, to keep the table readable)
                    if label == primary and metric == "con":
                        ident = np.array(d["ident"], dtype=float)
                        for tier, m in (
                                ("identity>=0.99", ident >= 0.99),
                                ("identity 0.98-0.99", (ident >= 0.98) & (ident < 0.99)),
                                ("no_alignment", ~np.isfinite(ident))):
                            if m.sum() < 30:
                                continue
                            sub = {k: [x for x, keep in zip(v, m) if keep] if k != "win_a"
                                   else [] for k, v in d.items()}
                            sub["win_a"] = [x for x, keep in zip(d["win_a"], m) if keep]
                            sub["win_b"] = [x for x, keep in zip(d["win_b"], m) if keep]
                            sub["chr_a"] = [x for x, keep in zip(d["chr_a"], m) if keep]
                            sub["chr_b"] = [x for x, keep in zip(d["chr_b"], m) if keep]
                            aa = assoc(sub["jac"], sub["con"], sub["chr_a"], sub["chr_b"],
                                       n_perm=2000, n_boot=500)
                            rows.append(_row(
                                analysis="class_association", **{"class": cls},
                                class_label=CLASS_LABEL.get(cls, cls), seq_regime=sr,
                                contact_regime=label, metric=mlabel, stratum=tier,
                                n_pairs=aa["n"], n_windows=len(_windows_in(sub)),
                                rho=aa["rho"], p_perm=aa["p_perm"], ci_lo=aa["ci_lo"],
                                ci_hi=aa["ci_hi"], n_perm=2000, n_boot=500))
                        med = float(np.median(rmvec))
                        for tier, m in (("repeat_frac<median", rmvec < med),
                                        ("repeat_frac>=median", rmvec >= med)):
                            if m.sum() < 30:
                                continue
                            sub = {k: [x for x, keep in zip(v, m) if keep] if k != "win_a"
                                   else [] for k, v in d.items()}
                            sub["win_a"] = [x for x, keep in zip(d["win_a"], m) if keep]
                            sub["win_b"] = [x for x, keep in zip(d["win_b"], m) if keep]
                            sub["chr_a"] = [x for x, keep in zip(d["chr_a"], m) if keep]
                            sub["chr_b"] = [x for x, keep in zip(d["chr_b"], m) if keep]
                            aa = assoc(sub["jac"], sub["con"], sub["chr_a"], sub["chr_b"],
                                       n_perm=2000, n_boot=500)
                            rows.append(_row(
                                analysis="class_association", **{"class": cls},
                                class_label=CLASS_LABEL.get(cls, cls), seq_regime=sr,
                                contact_regime=label, metric=mlabel, stratum=tier,
                                n_pairs=aa["n"], n_windows=len(_windows_in(sub)),
                                rho=aa["rho"], p_perm=aa["p_perm"], ci_lo=aa["ci_lo"],
                                ci_hi=aa["ci_hi"], n_perm=2000, n_boot=500))
                # chromosome-level Mantel (secondary, = the paper's arm-level statistic shape)
                mant = class_mantel(d)
                rows.append(_row(
                    analysis="chromosome_mantel", **{"class": cls},
                    class_label=CLASS_LABEL.get(cls, cls), seq_regime=sr,
                    contact_regime=label, metric="balanced_contact", stratum="all",
                    n_pairs=mant.get("n_cell", 0), rho=mant["rho"], p_perm=mant["p_perm"],
                    n_perm=2000, n_boot=500, note="aggregated to chromosome pairs"))

            # ---------- matched contrast: same pipeline, matched windows only ----------
            for target_cls, (mT, mC) in match_sets.items():
                dT = _subset(design["T"], mT)
                dC = _subset(design[target_cls], mC)
                if len(dT["jac"]) < 30 or len(dC["jac"]) < 30:
                    continue
                for metric, mlabel in (("con", "balanced_contact"),
                                       ("con_oe", "contact_O_over_E")):
                    aT = assoc(dT["jac"], dT[metric], dT["chr_a"], dT["chr_b"],
                               n_perm=n_perm, n_boot=n_boot)
                    aC = assoc(dC["jac"], dC[metric], dC["chr_a"], dC["chr_b"],
                               n_perm=n_perm, n_boot=n_boot)
                    dl = delta_rho_bootstrap(dT, dC, metric, n_boot=n_boot)
                    fz = fisher_z_diff(aT["rho"], aT["n"], aC["rho"], aC["n"])
                    for cls, aa, nn in ((f"T~{target_cls}:T", aT, len(mT)),
                                        (f"T~{target_cls}:ctrl", aC, len(mC))):
                        rows.append(_row(
                            analysis="matched_association", **{"class": cls},
                            class_label=("subtelomeric PHR windows (matched)" if ":T" in cls
                                         else f"{CLASS_LABEL.get(target_cls, target_cls)} "
                                              f"(matched)"),
                            seq_regime=sr, contact_regime=label, metric=mlabel, stratum="all",
                            n_pairs=aa["n"], n_windows=nn, rho=aa["rho"], p_perm=aa["p_perm"],
                            ci_lo=aa["ci_lo"], ci_hi=aa["ci_hi"], n_perm=n_perm, n_boot=n_boot,
                            delta_vs=target_cls, delta_ci_lo=dl["lo"], delta_ci_hi=dl["hi"],
                            delta_p_boot=dl["p"], fisher_z_p=fz["p"],
                            matched_n_windows=nn,
                            note=f"delta_rho={dl['delta']:.4g} (T minus {target_cls})"))
    per_contact = None
    tsv = OUT_TSV
    with open(tsv, "w") as f:
        f.write("\t".join(ROW_COLS) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r[c]) for c in ROW_COLS) + "\n")
    log(f"wrote {len(rows)} result rows -> {tsv}")
    detail["n_rows"] = len(rows)
    return detail


# =============================================================== oracle (tiny-scale check)
def stage_oracle(workdir: Path = None) -> dict:
    """Brute-force verification of the Jaccard estimator on hand-built sequences with known
    overlap, plus an independent per-base recomputation from the same PAF."""
    d = (workdir or (WORK / "oracle"))
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1234)
    L = 100_000
    base = "".join(rng.choice(list("ACGT"), L))
    mut = list(base)
    for pos in rng.choice(L, 200, replace=False):        # ~2% substitutions
        mut[pos] = rng.choice([c for c in "ACGT" if c != mut[pos]])
    w1 = base
    w2 = "".join(mut)
    w3 = "".join(rng.choice(list("ACGT"), L))
    w4 = base[: L // 2] + w3[: L // 2]
    w5 = base                                            # exact duplicate of w1
    # NOTE: names deliberately carry no '#'. wfmash groups mappings by the last '#'-delimited
    # prefix (-Y '#') and then skips self-group mappings, so PanSN-style names would make every
    # synthetic sequence its own target group's sibling and suppress the very homology we are
    # testing. The production window names are likewise '#'-free (samtools region names).
    names = ["w1", "w2", "w3", "w4", "w5"]
    seqs = [w1, w2, w3, w4, w5]
    fa = d / "oracle.fa"
    with open(fa, "w") as f:
        for n, s in zip(names, seqs):
            f.write(f">{n}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fa)], check=True)
    paf = d / "oracle.paf"
    with open(paf, "w") as fh:
        subprocess.run(["wfmash", str(fa), "-p", "90", "-L", "-f", "-M", "-t", "4",
                        "-B", str(d), "--quiet"], stdout=fh, check=True)
    exe = compile_jaccard(d)
    out = d / "oracle_pairs.tsv"
    subprocess.run([str(exe), str(paf), str(fa) + ".fai", str(out)], check=True)
    fai = [l.split("\t")[0] for l in open(str(fa) + ".fai")]
    lens = [int(l.split("\t")[1]) for l in open(str(fa) + ".fai")]
    got = {}
    with open(out) as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            a, b = fai[int(c[0])], fai[int(c[1])]
            shared = 0.5 * (int(c[2]) + int(c[3]))
            got[_key(a, b)] = dict(shared=shared, jac=shared / (lens[int(c[0])] + lens[int(c[1])]
                                                              - shared))
    # independent per-base recomputation from the same PAF, per pair
    L2 = {n: len(sq) for n, sq in zip(names, seqs)}
    pairpos = defaultdict(lambda: {n: None for n in names})
    with open(paf) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            q, t = c[0], c[5]
            qs, qe, ts, te = int(c[2]), int(c[3]), int(c[7]), int(c[8])
            if q > t:
                q, t, qs, qe, ts, te = t, q, ts, te, qs, qe
            key = (q, t)
            slot = pairpos[key]
            if slot[q] is None:
                slot[q] = np.zeros(L2[q], dtype=bool)
                slot[t] = np.zeros(L2[t], dtype=bool)
            slot[q][qs:qe] = True
            slot[t][ts:te] = True
    checks = []
    expected = {
        ("w1", "w5"): (1.0, 0.02, "identical sequences -> Jaccard 1"),
        ("w1", "w3"): (0.0, 0.02, "unrelated sequences -> Jaccard 0"),
        ("w1", "w4"): (1.0 / 3.0, 0.05, "half-shared -> Jaccard 1/3"),
        ("w1", "w2"): (1.0, 0.02, "~2% divergent -> Jaccard ~1"),
    }
    ok = True
    for (a, b), (exp, tol, why) in expected.items():
        v = got.get(_key(a, b), dict(jac=0.0))["jac"]
        passed = abs(v - exp) <= tol
        ok &= passed
        checks.append(dict(pair=f"{a}|{b}", expected=exp, got=v, tol=tol, passed=bool(passed),
                           why=why))
        log(f"  oracle {a} vs {b}: expected {exp:.4f} got {v:.4f} "
            f"({'PASS' if passed else 'FAIL'}) -- {why}")
    # integer equality against the independent per-base union
    for (a, b), v in got.items():
        slot = pairpos[(a, b)]
        brute = 0.5 * (slot[a].sum() + slot[b].sum())
        if abs(brute - v["shared"]) > 1e-9:
            ok = False
            checks.append(dict(pair=f"{a}|{b}", brute_force=float(brute),
                               estimator=float(v["shared"]), passed=False,
                               why="per-base union disagrees with interval-union estimator"))
            log(f"  oracle brute-force mismatch {a} vs {b}: per-base {brute} vs interval "
                f"{v['shared']}")
    if not ok:
        die("oracle FAILED: the Jaccard estimator does not behave as specified")
    log(f"oracle: {len(checks)} checks, all passed ({paf})")
    return dict(oracle_fasta=str(fa), paf=str(paf), n_pairs=len(got), checks=checks)


# =============================================================== main
def merge_json(record: dict) -> None:
    old = {}
    if OUT_JSON.exists():
        try:
            old = json.loads(OUT_JSON.read_text())
        except json.JSONDecodeError:
            old = {}
    old.update(record)
    old.setdefault("provenance", {})
    old["provenance"].update({
        "window_table": str(WINDOW_TABLE),
        "phr_bed": str(PHR_BED),
        "annotations_bed": str(ANNOT_BED),
        "chm13_fasta": str(CHM13_FA),
        "repeatmasker_bed": str(RM_BED),
        "mcools": {k: str(v) for k, v in MCools.items()},
        "published_pair_table": str(PUB_CHM13_PAIRS),
        "arm_dist_matrix": str(ARM_DIST),
        "workdir": str(WORK),
        "wfmash": (subprocess.run(["wfmash", "--version"], capture_output=True,
                                 text=True).stdout.strip() or "unknown")
        if shutil.which("wfmash") else "missing",
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    })
    OUT_JSON.write_text(json.dumps(old, indent=2, default=str) + "\n")
    log(f"wrote {OUT_JSON}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["gates", "oracle", "windows", "extract", "align", "jaccard",
                                      "liftover", "contact", "stats", "all"])
    ap.add_argument("--threads", type=int, default=96)
    ap.add_argument("--procs", type=int, default=48)
    ap.add_argument("--regime", default=None, help="sequence regime for `jaccard`")
    ap.add_argument("--perms", type=int, default=5000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--skip-p90", action="store_true",
                    help="skip the second (p90/pggb) alignment regime")
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)

    if a.stage == "gates":
        merge_json(dict(gates=[gate_published_rho(), gate_contact_contract(),
                               gate_arm_mantel()]))
        return
    if a.stage == "oracle":
        merge_json(dict(oracle=stage_oracle()))
        return
    if a.stage == "windows":
        merge_json(dict(window_classes=stage_windows()))
        return
    if a.stage == "extract":
        merge_json(dict(extraction=stage_extract()))
        return
    if a.stage == "align":
        merge_json(dict(alignment=stage_align(threads=a.threads)))
        return
    if a.stage == "jaccard":
        regs = [a.regime] if a.regime else ["p95", "p90_pggb"]
        merge_json(dict(sequence_axis={r: stage_jaccard(r) for r in regs}))
        return
    if a.stage == "liftover":
        merge_json(dict(liftover=stage_liftover(threads=a.threads)))
        return
    if a.stage == "contact":
        merge_json(dict(contact=stage_contact(procs=a.procs)))
        return
    if a.stage == "stats":
        merge_json(dict(statistics=stage_stats(n_perm=a.perms, n_boot=a.boot)))
        return
    if a.stage == "all":
        rec = dict(gates=[gate_published_rho(), gate_contact_contract(), gate_arm_mantel()],
                   oracle=stage_oracle(),
                   window_classes=stage_windows(),
                   extraction=stage_extract(),
                   liftover=stage_liftover(threads=a.threads),
                   alignment=stage_align(threads=a.threads),
                   contact=stage_contact(procs=a.procs))
        regs = ["p95"] if a.skip_p90 else ["p95", "p90_pggb"]
        rec["sequence_axis"] = {r: stage_jaccard(r) for r in regs}
        rec["statistics"] = stage_stats(seq_regimes=tuple(regs), n_perm=a.perms, n_boot=a.boot)
        merge_json(rec)
        return
    die(f"unknown stage {a.stage}")


if __name__ == "__main__":
    main()
