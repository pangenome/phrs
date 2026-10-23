#!/usr/bin/env python3
"""M6: mouse p-arm satellite confound -- quantification and sensitivity of the arm-level
Mantel association to satellite content.

Referee #1, major point 6:
  "Mouse chromosomes are acro/telocentric with essentially no short arm, so a 1 Mb
   telomere-anchored 'p-arm' flank is largely composed of minor and major satellite, which
   are near-identical between chromosomes, placing these windows at once at the top of the
   Jaccard axis and at the maximum of the multi-mapping problem."

What this script does
---------------------
(a) Annotates satellite content of every 1 Mb telomere-anchored mouse flank used by the PHR
    call, at 5 kb resolution, from two independent sources:

    SOURCE 1 (ANNOTATION, not a proxy): blastn of the flank sequences against the curated
    Dfam 3.8 mouse satellite consensus families shipped inside RepeatMasker 4.1.7p1
    (RepeatMasker.lib).  Dfam 3.8 in this local library has a family explicitly named for
    mouse major satellite (GSAT_MM) but *no* family explicitly named "minor satellite" for
    Mus; the minor-satellite figure reported here therefore comes from MurSAT1 / MMSAT4 /
    CENSAT_MC and is labelled PROXY everywhere it is used.  Per-family coverage is kept
    separate so major and minor can be read independently.

    SOURCE 2 (PROXY, composition only, no annotation at all): tandem-periodicity detector.
    A 5 kb bin is PROXY-satellite-like when the fraction of positions i with s[i]==s[i+lag]
    exceeds 0.5 for some lag in 100..260 bp (major satellite is a 234 bp tandem repeat,
    minor satellite a 120 bp one), while the 6 bp match fraction (telomeric CCCTAA repeat)
    is below 0.8.  Computed in C (embedded below, compiled at run time) -- no Python over
    the sequence.

(b) Re-runs the arm-level Mantel test of scripts/hic/mantel_mouse_d_m5.py
    (imported, NOT modified) with the arm-pair set stratified by the satellite fraction of
    the two flanking windows, at 20 kb and 50 kb, for all four meiotic stages.
    Also re-derives the sequence-similarity axis from the raw wfmash PAFs so it can be
    recomputed with satellite bins excluded (the published arm matrix is an odgi graph-node
    Jaccard and cannot be recomputed without rebuilding the pggb graph), and runs the
    Mantel on that satellite-masked sequence axis.

    NOT done here (see notes): the contact axis is not re-masked from raw Hi-C.  That needs
    the .mcool files read per-bin; h5py/cooler are absent from the system python and the
    per-stage balanced contact matrices are 0.3-0.5 Gpixels each.  The contact values are
    therefore re-used exactly as in mantel_mouse_d_m5.py (mean normalised contact per arm
    pair aggregated from data/mouse_<stage>_phr_<res>bp_seqlevel.tsv).

(c) Reports whether the stage peak (pachytene at 20 kb in D-M5) and the p-arm vs q-arm
    contrast survive satellite masking, with permutation tests.

Verification pyramid
--------------------
Tiny-scale brute-force oracle first (--selftest, run automatically):
  * interval-union coverage vs a base-by-base bitmask oracle;
  * the Mantel statistic from mantel_mouse_d_m5.mantel vs scipy spearmanr;
  * a synthetic perfect-correlation case with a known permutation p.

Usage:
  python3 scripts/hic/mouse_satellite_sensitivity.py [--perms 10000] [--skip-blast]
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
M5_PATH = Path(__file__).resolve().parent / "mantel_mouse_d_m5.py"
OUT = Path(__file__).resolve().parent / "mouse_satellite_sensitivity.json"

MOUSE = Path("/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T")
FLANK_FA = MOUSE / "subtelo_1Mb/mouse_subtelo_1Mb_trimmed.fa.gz"
PAF_DIR = MOUSE / "subtelo_1Mb/all-vs-all.p95"
PHR_TSV = MOUSE / "subtelo_1Mb/mouse.all-vs-all.p95.id95.len.tsv"
DIST = MOUSE / "similarity_1Mb/mouse.dist_matrix.tsv"
RMSK_LIB = Path("/moosefs/guarracino/condatools/repeatmasker/4.1.7p1/share/RepeatMasker"
                "/Libraries/RepeatMasker.lib")
BLASTN = Path("/moosefs/guarracino/condatools/repeatmasker/4.1.7p1/bin/blastn")

STAGES = ["leptotene", "zygotene", "pachytene", "diplotene"]
RESOLUTIONS = (20000, 50000)
BIN = 5000

# Dfam families used.  Major satellite: GSAT_MM (the Dfam mouse major-satellite family).
# "Minor-satellite PROXY": MurSAT1 (120 bp, @Muridae -- the 120 bp mouse minor-satellite
# period), MMSAT4, CENSAT_MC.  CYRA11_Mm is Y-chromosome satellite and is reported apart.
MAJOR = ["GSAT_MM"]
MINOR_PROXY = ["MurSAT1", "MMSAT4", "CENSAT_MC"]
OTHER_SAT = ["SYNREP_MM", "IMPB_01", "ZP3AR"]
Y_SAT = ["CYRA11_Mm"]
ALL_FAMILIES = MAJOR + MINOR_PROXY + OTHER_SAT + Y_SAT

C_SRC = r"""
/* Per-5kb-bin tandem-periodicity / base-composition scanner for M6.
 * Usage: binscan <fasta> <bin>  ->  TSV on stdout:
 *   seq bin_start bin_end len N_frac AT_frac lag6 best_lag best_match m120 m234
 * best_match = max over lag in [100,260] of #{i: s[i]==s[i+lag]} / (len-lag)
 * Independent of any repeat annotation: pure sequence composition.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLEN (4u*1024u*1024u)
#define LAGLO 100
#define LAGHI 260

static char *buf;
static size_t bufn;

static void die(const char *m) { fprintf(stderr, "binscan: %s\n", m); exit(1); }

static void emit(const char *name, size_t off, size_t len, size_t bin)
{
    size_t n = len < bin ? len : bin;
    size_t i, lag, a = 0, c = 0, g = 0, t = 0, nn = 0;
    double best = -1.0; int bestlag = -1;
    double m120 = -1.0, m234 = -1.0, lag6 = -1.0;
    if (n < LAGHI + 2) return;
    for (i = 0; i < n; i++) {
        switch (buf[off + i]) {
            case 'A': a++; break; case 'C': c++; break;
            case 'G': g++; break; case 'T': t++; break;
            default: nn++; break;
        }
    }
    {   /* telomeric-repeat indicator: 6 bp periodicity */
        size_t m = 0, tot = 0;
        for (i = 0; i + 6 < n; i++) { tot++; if (buf[off+i] == buf[off+i+6]) m++; }
        lag6 = tot ? (double)m / (double)tot : -1.0;
    }
    for (lag = LAGLO; lag <= LAGHI; lag++) {
        size_t m = 0, tot = 0;
        double f;
        if (n <= lag) continue;
        for (i = 0; i + lag < n; i++) { tot++; if (buf[off+i] == buf[off+i+lag]) m++; }
        f = (double)m / (double)tot;
        if (lag == 120) m120 = f;
        if (lag == 234) m234 = f;
        if (f > best) { best = f; bestlag = (int)lag; }
    }
    printf("%s\t%zu\t%zu\t%zu\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n",
           name, off, off + n, n,
           n ? (double)nn / (double)n : 0.0,
           n ? (double)(a + t) / (double)n : 0.0,
           lag6, bestlag, best, m120, m234);
}

int main(int argc, char **argv)
{
    FILE *f;
    size_t bin;
    char line[70000];
    char *name = NULL;
    size_t off = 0;

    if (argc != 3) die("usage: binscan <fasta> <bin>");
    bin = (size_t)strtoull(argv[2], NULL, 10);
    if (bin < LAGHI + 2) die("bin too small");
    f = fopen(argv[1], "r");
    if (!f) die("cannot open fasta");
    buf = malloc(MAXLEN);
    if (!buf) die("oom");
    bufn = 0;
    name = malloc(4096); if (!name) die("oom");
    name[0] = '\0';

    while (fgets(line, sizeof line, f)) {
        if (line[0] == '>') {
            size_t bo;
            if (name[0]) for (bo = 0; bo < bufn; bo += bin) emit(name, bo, bufn - bo, bin);
            bufn = 0;
            line[strcspn(line, "\r\n")] = '\0';
            snprintf(name, 4096, "%s", line + 1);
            off = 0;
        } else {
            size_t k, L = strlen(line);
            for (k = 0; k < L; k++) {
                char ch = line[k];
                if (ch == '\n' || ch == '\r') continue;
                if (bufn >= MAXLEN) die("sequence longer than MAXLEN");
                buf[bufn++] = (ch >= 'a' && ch <= 'z') ? (char)(ch - 32) : ch;
            }
        }
    }
    if (name[0]) for (off = 0; off < bufn; off += bin) emit(name, off, bufn - off, bin);
    fclose(f);
    free(buf); free(name);
    return 0;
}
"""


# --------------------------------------------------------------------------- helpers

def run(cmd, **kw):
    """Run a command, fail loud."""
    print("$", " ".join(str(c) for c in cmd), flush=True)
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:])
        sys.stderr.write(r.stderr[-4000:])
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(str(c) for c in cmd)}")
    return r.stdout


def load_m5():
    """Import scripts/hic/mantel_mouse_d_m5.py without modifying it."""
    spec = importlib.util.spec_from_file_location("mantel_mouse_d_m5", M5_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def merge_len(iv):
    """Total length of the union of half-open integer intervals."""
    if not iv:
        return 0
    iv = sorted(iv)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        elif e > ce:
            ce = e
    return tot + (ce - cs)


def merge_minus_len(iv, mask):
    """Total length of union(iv) after subtracting union(mask) (both half-open)."""
    if not iv:
        return 0
    return merge_len(iv) - merge_len(intersect_intervals(iv, mask))


def intersect_intervals(a, b):
    """Intersection of two interval lists (each may be internally overlapping)."""
    if not a or not b:
        return []
    a = sorted(a)
    b = sorted(b)
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _ranks(v):
    return np.argsort(np.argsort(np.asarray(v)))


def boot_rho_ci(xs, ys, rng, B=1000):
    """Percentile bootstrap CI of the rank correlation, resampling arm pairs.

    The CI in mantel_mouse_d_m5.mantel is the spread of the PERMUTATION NULL, not a CI of
    rho; this is the latter and is computed here without touching that file.
    """
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    n = len(x)
    if n < 4:
        return None
    vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        xx, yy = x[idx], y[idx]
        if np.ptp(xx) == 0 or np.ptp(yy) == 0:
            continue
        vals.append(float(np.corrcoef(_ranks(xx), _ranks(yy))[0, 1]))
    if len(vals) < 50:
        return None
    return [round(float(np.percentile(vals, 2.5)), 4),
            round(float(np.percentile(vals, 97.5)), 4)]


def paf_parse_columns(c):
    """Return (qname, qstart, qend, tname, matches, blocklen, gi) from a split PAF line.

    PAF column order (1-based): 1 qname, 2 qlen, 3 qstart, 4 qend, 5 strand, 6 tname,
    7 tlen, 8 tstart, 9 tend, 10 matches, 11 blocklen, 12 mapq.  Off-by-one here silently
    produced empty interval lists (caught by the selftest below, which pins the fields to a
    real record from mouse.all-vs-B6_1_chr10_parm_8494_1000000.paf.gz).
    """
    gi = None
    for tag in c[12:]:
        if tag.startswith("gi:f:"):
            gi = float(tag[5:])
            break
    return (c[0], int(c[2]), int(c[3]), c[5], int(c[9]), int(c[10]), gi)


def paf_identity_ok(matches, blocklen, gi=None, min_identity=0.95):
    """True when an alignment's identity reaches the threshold.

    wfmash (run with -p 95) reports the filtering identity as the gap-compressed gi:f tag;
    matches/blocklen is lower because blocklen counts indels.  Prefer gi:f when present --
    using matches/blocklen alone silently discarded ~20% of the qualifying alignments
    (caught by comparing per-pair record counts against the gi filter).
    """
    if gi is not None:
        return gi >= min_identity
    if blocklen <= 0:
        return False
    return (matches / blocklen) >= min_identity


def read_fasta_gz(path):
    """Stream a bgzipped FASTA (plain gzip compatible) -> dict name->sequence."""
    seqs = {}
    name = None
    chunks = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def flank_of(seq_name):
    """'B6#1#chr10_parm:8494-1000000:1-990000' -> 'B6#1#chr10_parm:8494-1000000'."""
    return seq_name.rsplit(":", 1)[0]


def end_of(flank):
    """'B6#1#chr10_parm:8494-1000000' -> ('B6#1#chr10', 'parm')."""
    stem, _ = flank.split(":")
    chrom, kind = stem.rsplit("_", 1)
    return chrom, kind


# --------------------------------------------------------------------------- (a) annotation

def extract_satellite_library(out_fa):
    """Pull the curated mouse satellite consensus sequences out of RepeatMasker.lib."""
    fam = {}
    name = None
    hdr = None
    chunks = []
    with open(RMSK_LIB) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    fam[name] = ("".join(chunks), hdr)
                hdr = line[1:].strip()
                name = hdr.split("#")[0].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        fam[name] = ("".join(chunks), hdr)
    missing = [k for k in ALL_FAMILIES if k not in fam]
    if missing:
        raise RuntimeError(f"satellite families missing from {RMSK_LIB}: {missing}")
    with open(out_fa, "w") as f:
        for k in ALL_FAMILIES:
            f.write(f">{fam[k][1]}\n{fam[k][0]}\n")
    return {k: len(fam[k][0]) for k in ALL_FAMILIES}, {k: fam[k][1] for k in ALL_FAMILIES}


def blast_satellites(flank_fa, lib_fa, out_tsv, tmpdir):
    """blastn each flank against the tiny curated satellite library; per-family hits."""
    seqs = {}
    name = None
    for line in open(flank_fa):
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    parts = []
    for i, (nm, chunks) in enumerate(seqs.items()):
        p = tmpdir / f"shard_{i:03d}.fa"
        with open(p, "w") as f:
            f.write(f">{nm}\n" + "".join(chunks) + "\n")
        parts.append(p)
    outs = []
    for i, p in enumerate(parts):
        o = tmpdir / f"shard_{i:03d}.hits"
        run([BLASTN, "-query", p, "-subject", lib_fa, "-task", "blastn",
             "-word_size", "11", "-dust", "no", "-evalue", "1e-10",
             "-perc_identity", "70", "-max_target_seqs", "5000",
             "-outfmt", "6 qseqid qstart qend sseqid pident length", "-out", o])
        outs.append(o)
    with open(out_tsv, "w") as fo:
        for o in outs:
            fo.write(o.read_text())
    return out_tsv


def hits_to_intervals(tsv, min_len=100, min_pid=75.0):
    """flank -> family -> list of (start,end) 0-based half-open intervals."""
    out = {}
    with open(tsv) as f:
        for line in f:
            q, qs, qe, s, pid, ln = line.rstrip("\n").split("\t")
            if float(ln) < min_len or float(pid) < min_pid:
                continue
            fam = s.split("#")[0]
            out.setdefault(q, {}).setdefault(fam, []).append((int(qs) - 1, int(qe)))
    return out


def bin_intervals(intervals, bin=BIN):
    """Set of bin indices touched by an interval list."""
    b = set()
    for s, e in intervals:
        for i in range(s // bin, (e - 1) // bin + 1):
            b.add(i)
    return b


def binwise_coverage(intervals, binlen, bin=BIN):
    """Per-bin covered fraction (0..1) from possibly overlapping intervals.

    Union semantics: each bin's intervals are merged before measuring coverage.  Verified
    bin-by-bin against a base-by-base bitmask oracle in selftest().
    """
    n = (binlen + bin - 1) // bin
    if n == 0:
        return np.zeros(0)
    sizes = np.full(n, bin, dtype=np.int64)
    if binlen % bin:
        sizes[-1] = binlen % bin
    out = np.zeros(n)
    buckets = {}
    for s, e in intervals:
        while s < e:
            bi = s // bin
            stop = min(e, (bi + 1) * bin)
            buckets.setdefault(bi, []).append((s, stop))
            s = stop
    for bi, iv in buckets.items():
        out[bi] = merge_len(iv) / sizes[bi]
    return out


# --------------------------------------------------------------------------- selftest

def selftest(m5):
    """Brute-force oracles at tiny scale.  Aborts the run on any mismatch."""
    rng = np.random.default_rng(7)
    # 1. interval union vs base-by-base bitmask
    for _ in range(200):
        L = int(rng.integers(20, 200))
        n = int(rng.integers(1, 8))
        iv = []
        for _ in range(n):
            s = int(rng.integers(0, L))
            e = int(rng.integers(s + 1, L + 1))
            iv.append((s, e))
        mask = np.zeros(L, dtype=bool)
        for s, e in iv:
            mask[s:e] = True
        got = merge_len(iv)
        want = int(mask.sum())
        if got != want:
            raise AssertionError(f"interval-union oracle: got {got} want {want} for {iv}")
        cv = binwise_coverage(iv, L, bin=7)
        fine = np.zeros((L + 6) // 7)
        for i in range(L):
            if mask[i]:
                fine[i // 7] += 1
        sizes = np.array([min(7, L - 7 * i) for i in range(len(fine))])
        if not np.allclose(cv, fine / sizes):
            raise AssertionError("binwise_coverage oracle mismatch")
    # 2. Mantel statistic vs scipy spearmanr, and known p for perfect correlation
    from scipy.stats import spearmanr
    for _ in range(20):
        x = rng.normal(size=40)
        y = 0.5 * x + rng.normal(size=40)
        obs, p, lo, hi, mu = m5.mantel(x, y, 500, np.random.default_rng(11))
        ref = spearmanr(x, y).statistic
        if abs(obs - ref) > 1e-9:
            raise AssertionError(f"mantel rho {obs} != spearman {ref}")
    obs, p, _, _, _ = m5.mantel(np.arange(30.0), np.arange(30.0), 200,
                                np.random.default_rng(3))
    if abs(obs - 1.0) > 1e-12 or p > 0.02:
        raise AssertionError(f"perfect-correlation oracle: rho={obs} p={p}")
    # 3. PAF identity filter orientation (matches/blocklen, NOT divergence)
    assert paf_identity_ok(6807, 6981)              # 97.5% identity -> keep
    assert paf_identity_ok(9500, 10000)             # exactly 95% -> keep
    assert not paf_identity_ok(9400, 10000)         # 94% -> drop
    assert not paf_identity_ok(0, 0)                # degenerate -> drop
    assert paf_identity_ok(9817, 10484, 0.99282)    # gi:f decides when present
    assert not paf_identity_ok(9817, 10484, 0.93)
    # 4. PAF column parsing against a verbatim record
    rec = ("B6#1#chr6_parm:21150-1000000\t978851\t0\t6981\t+\t"
           "B6#1#chr10_parm:8494-1000000\t991507\t0\t6835\t6807\t6981\t16\t"
           "gi:f:0.995758\tbi:f:0.974796").split("\t")
    q, qs, qe, t, mm, blk, gi = paf_parse_columns(rec)
    assert q == "B6#1#chr6_parm:21150-1000000" and t == "B6#1#chr10_parm:8494-1000000"
    assert (qs, qe) == (0, 6981), (qs, qe)
    assert (mm, blk) == (6807, 6981) and abs(gi - 0.995758) < 1e-9
    assert qs < qe and paf_identity_ok(mm, blk, gi)
    # 5. end/flank naming parser
    assert flank_of("B6#1#chr10_parm:8494-1000000:1-990000") == "B6#1#chr10_parm:8494-1000000"
    assert end_of("B6#1#chr10_parm:8494-1000000") == ("B6#1#chr10", "parm")
    return True


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--skip-blast", action="store_true")
    a = ap.parse_args()

    m5 = load_m5()
    if not selftest(m5):
        raise SystemExit("selftest failed")
    print("selftest: interval-union + binwise-coverage + Mantel-statistic oracles OK")

    tmp = Path(a.tmpdir or tempfile.mkdtemp(prefix="m6_mouse_sat_"))
    tmp.mkdir(parents=True, exist_ok=True)
    print("scratch dir:", tmp)

    for p in (FLANK_FA, PHR_TSV, DIST, RMSK_LIB):
        if not p.exists():
            raise SystemExit(f"missing required input: {p}")

    # ---- flanks
    seqs = read_fasta_gz(FLANK_FA)
    print(f"flanks: {len(seqs)} sequences, {sum(len(s) for s in seqs.values())} bp")

    # ---- (a1) Dfam-derived satellite annotation
    lib_fa = tmp / "mouse_sat.lib.fa"
    famlen, famhdr = extract_satellite_library(lib_fa)
    print("satellite library families:", famlen)
    hits_tsv = tmp / "sat_hits.tsv"
    if a.skip_blast and hits_tsv.exists():
        print("re-using", hits_tsv)
    else:
        if not BLASTN.exists():
            raise SystemExit(f"blastn not found: {BLASTN}")
        flank_fa = tmp / "flanks.fa"
        with open(flank_fa, "w") as f:
            for nm, sq in seqs.items():
                f.write(f">{nm}\n{sq}\n")
        blast_satellites(flank_fa, lib_fa, hits_tsv, tmp)
    hits = hits_to_intervals(hits_tsv)
    print(f"flanks with satellite hits: {len(hits)} / {len(seqs)}")

    # ---- (a2) composition PROXY in C
    c_src = tmp / "binscan.c"
    c_bin = tmp / "binscan"
    c_src.write_text(C_SRC)
    run(["cc", "-O2", "-o", c_bin, c_src])
    flanks_fa = tmp / "flanks.fa"
    if not flanks_fa.exists():
        with open(flanks_fa, "w") as f:
            for nm, sq in seqs.items():
                f.write(f">{nm}\n{sq}\n")
    comp_out = run([c_bin, flanks_fa, str(BIN)])
    print("composition bins produced:", len(comp_out.splitlines()))

    proxy = {}
    for line in comp_out.splitlines():
        nm, s, e, ln, nf, at, lag6, bl, bm, m120, m234 = line.split("\t")
        proxy.setdefault(nm, {})[int(s) // BIN] = (
            dict(start=int(s), end=int(e), len=int(ln), n_frac=float(nf),
                 at_frac=float(at), lag6=float(lag6), best_lag=int(bl),
                 best_match=float(bm), m120=float(m120),
                 m234=float(m234),
                 proxy_sat=(float(bm) >= 0.5 and float(lag6) < 0.8)))

    # ---- per-flank satellite summary at 5 kb
    fam_cov = {}
    per_flank = {}
    for nm, sq in sorted(seqs.items()):
        L = len(sq)
        fam_iv = hits.get(nm, {})
        cov = {}
        for fam in ALL_FAMILIES:
            dv = fam_iv.get(fam, [])
            cov[fam] = binwise_coverage(dv, L) if dv else np.zeros((L + BIN - 1) // BIN)
        union_iv = []
        for fam in ALL_FAMILIES:
            union_iv.extend(fam_iv.get(fam, []))
        union_cov = binwise_coverage(union_iv, L)
        major_iv = []
        for fam in MAJOR:
            major_iv.extend(fam_iv.get(fam, []))
        minor_iv = []
        for fam in MINOR_PROXY:
            minor_iv.extend(fam_iv.get(fam, []))
        major_cov = binwise_coverage(major_iv, L)
        minor_cov = binwise_coverage(minor_iv, L)
        p_ = proxy.get(nm, {})
        p_list = [p_[i] for i in sorted(p_)]
        proxy_sat_frac = (sum(b["len"] for b in p_list if b["proxy_sat"]) /
                          max(1, sum(b["len"] for b in p_list))) if p_list \
            else float("nan")
        per_flank[nm] = dict(
            length=L,
            n_bins=len(union_cov),
            n_proxy_bins=len(p_list),
            sat_bp_union=int(round(float(union_cov.sum()) * BIN)),
            sat_frac_union=float(union_cov.mean()),
            sat_frac_major_GSAT_MM=float(major_cov.mean()),
            sat_frac_minor_PROXY=float(minor_cov.mean()),
            frac_bins_sat50=float((union_cov >= 0.5).mean()),
            frac_bins_sat25=float((union_cov >= 0.25).mean()),
            proxy_periodic_bin_frac=proxy_sat_frac,
            proxy_best_lag_median=float(np.median([b["best_lag"] for b in p_list]))
            if p_list else None,
            at_frac=float(np.mean([b["at_frac"] * b["len"] for b in p_list]) /
                          max(1, sum(b["len"] for b in p_list))) if p_list else None,
            # agreement between the two independent sources, bin level
            agree_bins=float(np.mean([(b["proxy_sat"] == (union_cov[i] >= 0.5))
                                      for i in sorted(p_) for b in [p_[i]]]))
            if p_list else None,
            union_cov=union_cov, proxy=p_,
        )
        fam_cov[nm] = cov

    # arm-level satellite fraction: mean over the haplotype flanks of that chromosome end
    by_end = {}
    for nm, d in per_flank.items():
        by_end.setdefault(end_of(nm), []).append(d)
    end_sat = {}
    for (chrom, kind), ds in by_end.items():
        end_sat[(chrom, kind)] = dict(
            sat_frac_union=float(np.mean([d["sat_frac_union"] for d in ds])),
            sat_frac_major=float(np.mean([d["sat_frac_major_GSAT_MM"] for d in ds])),
            sat_frac_minor_proxy=float(np.mean([d["sat_frac_minor_PROXY"] for d in ds])),
            frac_bins_sat50=float(np.mean([d["frac_bins_sat50"] for d in ds])),
            proxy_periodic=float(np.mean([d["proxy_periodic_bin_frac"] for d in ds])),
            n_flanks=len(ds))
    print(f"chromosome ends annotated: {len(end_sat)}")

    # ---- flank -> arm name, from the seqlevel tables (authoritative for the Mantel)
    flank2arm = {}
    for stage in STAGES:
        for res in RESOLUTIONS:
            f = ROOT / f"data/mouse_{stage}_phr_{res}bp_seqlevel.tsv"
            if not f.exists():
                continue
            with open(f) as fh:
                next(fh)
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    flank2arm[flank_of(c[0])] = c[3]
                    flank2arm[flank_of(c[1])] = c[5]
    if not flank2arm:
        raise SystemExit("no seqlevel tables found under data/")
    print(f"flank->arm entries: {len(flank2arm)}")

    arm_sat = {}
    for fl, arm in flank2arm.items():
        chrom, kind = end_of(fl)
        if (chrom, kind) not in end_sat:
            raise SystemExit(f"flank {fl} has no satellite annotation")
        arm_sat.setdefault(arm, []).append(end_sat[(chrom, kind)])
    arm_sat_frac = {k: float(np.mean([d["sat_frac_union"] for d in v]))
                    for k, v in arm_sat.items()}
    arm_sat_major = {k: float(np.mean([d["sat_frac_major"] for d in v]))
                     for k, v in arm_sat.items()}
    arm_sat_minor = {k: float(np.mean([d["sat_frac_minor_proxy"] for d in v]))
                     for k, v in arm_sat.items()}

    # ---- sequence axis re-derived from the PAFs (maskable, unlike the odgi graph Jaccard)
    print("re-deriving sequence-sharing similarity from", PAF_DIR)
    phr = {}
    with open(PHR_TSV) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        i_seq, i_st, i_en = hdr.index("seq"), hdr.index("region_start"), hdr.index("region_end")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[i_st] == "." or c[i_en] == ".":
                phr[c[i_seq]] = []
            else:
                phr[c[i_seq]] = [(int(c[i_st]), int(c[i_en]))]
    paf_shared = {}   # (q,t) -> list of query-coordinate intervals (identity>=0.95)
    nfiles = 0
    for pf in sorted(PAF_DIR.glob("*.paf.gz")):
        nfiles += 1
        with gzip.open(pf, "rt") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                q, qs, qe, strand, t = c[0], int(c[1]), int(c[2]), c[4], c[5]
                gi = None
                for tag in c[12:]:
                    if tag.startswith("gi:f:"):
                        gi = float(tag[5:])
                        break
                if not paf_identity_ok(int(c[9]), int(c[10]), gi):
                    continue
                paf_shared.setdefault((q, t), []).append((qs, qe))
    print(f"paf files read: {nfiles}, ordered flank pairs with >=95% alignments: {len(paf_shared)}")

    def sat_intervals(fl):
        """5 kb bins of this flank that are >=50% satellite, as flank-local intervals."""
        d = per_flank[fl]
        cov = d["union_cov"]
        iv = []
        for i, v in enumerate(cov):
            if v >= 0.5:
                iv.append((i * BIN, min((i + 1) * BIN, d["length"])))
        return iv

    sat_iv = {fl: sat_intervals(fl) for fl in per_flank}

    def pair_shared(fla, flb, masked):
        """Symmetric shared bases between two flanks' PHRs; optionally satellite-masked."""
        ab = paf_shared.get((fla, flb), [])
        ba = paf_shared.get((flb, fla), [])
        pa = phr.get(fla, [])
        pb = phr.get(flb, [])
        if not ab and not ba:
            return 0.0, 0, 0
        # shared bases restricted to each flank's called PHR region, symmetrised as the
        # mean of the two alignment directions
        ivq = intersect_intervals(ab, pa)
        ivr = intersect_intervals(ba, pb)
        la = merge_len(pa)
        lb = merge_len(pb)
        if masked:
            # a shared base is dropped when either side of the alignment falls in a
            # satellite bin, and the denominators lose the satellite bins too
            ivq = merge_len(ivq) - merge_len(intersect_intervals(ivq, sat_iv[fla]))
            ivr = merge_len(ivr) - merge_len(intersect_intervals(ivr, sat_iv[flb]))
            la = max(0, la - merge_len(intersect_intervals(pa, sat_iv[fla])))
            lb = max(0, lb - merge_len(intersect_intervals(pb, sat_iv[flb])))
        sh = max(0.0, 0.5 * ((ivq if masked else merge_len(ivq)) +
                             (ivr if masked else merge_len(ivr))))
        if la + lb - sh <= 0:
            return 0.0, la, lb
        return sh / (la + lb - sh), la, lb

    arms, dist = m5.read_dist()
    all_fl = sorted(per_flank)
    for tag, masked in (("unmasked", False), ("masked", True)):
        acc = {}
        for i in range(len(all_fl)):
            for j in range(i, len(all_fl)):
                fla, flb = all_fl[i], all_fl[j]
                arma, armb = flank2arm.get(fla), flank2arm.get(flb)
                if arma is None or armb is None or arma == armb:
                    continue
                jac, la, lb = pair_shared(fla, flb, masked)
                key = (arma, armb) if arma < armb else (armb, arma)
                acc.setdefault(key, []).append(jac)
        seqjac = {k: float(np.mean(v)) for k, v in acc.items()}
        if tag == "unmasked":
            seqjac_unmasked = seqjac
        else:
            seqjac_masked = seqjac

    # validate the derived axis against the published odgi graph Jaccard
    both = [k for k in seqjac_unmasked if k in dist]
    from scipy.stats import spearmanr, pearsonr
    x = [seqjac_unmasked[k] for k in both]
    y = [1.0 - dist[k] for k in both]
    val = dict(n=len(both),
               spearman=float(spearmanr(x, y).statistic),
               pearson=float(pearsonr(x, y)[0]))
    print(f"derived Jaccard vs published odgi Jaccard: n={val['n']} "
          f"spearman={val['spearman']:.3f} pearson={val['pearson']:.3f}")
    if val["spearman"] < 0.5:
        print("WARNING: derived sequence axis is a weak stand-in for the published graph "
              "Jaccard; masked-axis results are reported but must be read as an "
              "approximation (see notes).")

    # ---- (b) Mantel: baseline, stratified, masked axis
    rng = np.random.default_rng(20260518)
    base = {}
    strat = {}
    pair_counts = {}
    for stage in STAGES:
        for res in RESOLUTIONS:
            con = m5.read_contact(stage, res)
            if not con:
                print(f"no contact table for {stage} {res}")
                continue
            keys = [k for k in con if k in dist and k[0] != k[1]]
            if len(keys) < 20:
                print(f"too few pairs for {stage} {res}")
                continue
            xs = [1.0 - dist[k] for k in keys]
            ys = [con[k] for k in keys]
            obs, p, lo, hi, mu = m5.mantel(xs, ys, a.perms, rng)
            base[f"{stage}_{res//1000}kb"] = dict(
                rho=round(obs, 4), p_perm=p, null_mean=round(mu, 4),
                ci95=[[round(lo, 4), round(hi, 4)]][0],
                boot_ci95_rho=boot_rho_ci(xs, ys, np.random.default_rng(99)),
                n_arm_pairs=len(keys))

            # stratification variable: satellite fraction of the pair =
            #   max over the two arms (a pair is "high satellite" if either flank is)
            pair_sat = {k: max(arm_sat_frac.get(k[0], np.nan),
                               arm_sat_frac.get(k[1], np.nan)) for k in keys}
            pair_sat_major = {k: max(arm_sat_major.get(k[0], np.nan),
                                     arm_sat_major.get(k[1], np.nan)) for k in keys}
            for thr in (0.10, 0.25, 0.50):
                lo_k = [k for k in keys if pair_sat[k] < thr]
                hi_k = [k for k in keys if pair_sat[k] >= thr]
                for lab, kk in (("low", lo_k), ("high", hi_k)):
                    if len(kk) < 15:
                        strat[f"{stage}_{res//1000}kb_sat{'<=' if lab=='low' else '>'}"
                              f"{thr}"] = dict(n_arm_pairs=len(kk), skipped="too few pairs")
                        continue
                    xa = [1.0 - dist[k] for k in kk]
                    ya = [con[k] for k in kk]
                    o2, p2, l2, h2, m2 = m5.mantel(xa, ya, a.perms, rng)
                    strat[f"{stage}_{res//1000}kb_sat_{lab}_{int(thr*100)}"] = dict(
                        rho=round(o2, 4), p_perm=p2, null_mean=round(m2, 4),
                        ci95=[round(l2, 4), round(h2, 4)], n_arm_pairs=len(kk),
                        threshold=thr,
                        boot_ci95_rho=boot_rho_ci(xa, ya, np.random.default_rng(98)),
                        mean_sat=round(float(np.mean([pair_sat[k] for k in kk])), 4))
            # major-satellite-only stratification (GSAT_MM only, no minor proxy)
            lo_k = [k for k in keys if pair_sat_major[k] < 0.05]
            hi_k = [k for k in keys if pair_sat_major[k] >= 0.05]
            for lab, kk in (("low", lo_k), ("high", hi_k)):
                if len(kk) < 15:
                    continue
                xa = [1.0 - dist[k] for k in kk]
                ya = [con[k] for k in kk]
                o2, p2, l2, h2, m2 = m5.mantel(xa, ya, a.perms, rng)
                strat[f"{stage}_{res//1000}kb_majorSat_{lab}_5"] = dict(
                    rho=round(o2, 4), p_perm=p2, n_arm_pairs=len(kk),
                    ci95=[round(l2, 4), round(h2, 4)],
                    boot_ci95_rho=boot_rho_ci(xa, ya, np.random.default_rng(97)))
            # masked sequence axis (satellite bins removed from both PHR sets)
            ke = [k for k in keys if k in seqjac_masked]
            if len(ke) >= 20:
                xs2 = [max(0.0, seqjac_masked[k]) for k in ke]
                ys2 = [con[k] for k in ke]
                o3, p3, l3, h3, m3 = m5.mantel(xs2, ys2, a.perms, rng)
                xs1 = [max(0.0, seqjac_unmasked[k]) for k in ke]
                o1, p1, l1, h1, m1 = m5.mantel(xs1, ys2, a.perms, rng)
                strat[f"{stage}_{res//1000}kb_axis_maskedSat"] = dict(
                    rho=round(o3, 4), p_perm=p3, n_arm_pairs=len(ke),
                    ci95=[round(l3, 4), round(h3, 4)],
                    boot_ci95_rho=boot_rho_ci(xs2, ys2, np.random.default_rng(96)),
                    rho_derived_axis_unmasked=round(o1, 4),
                    p_derived_axis_unmasked=p1,
                    boot_ci95_rho_derived_unmasked=boot_rho_ci(
                        xs1, ys2, np.random.default_rng(95)))
            pair_counts[f"{stage}_{res//1000}kb"] = len(keys)

    # ---- (c) p-arm vs q-arm contrast, overall and after satellite masking
    def arm_type(arm):
        return "p" if arm.endswith("_p") else "q"

    pq = {}
    for stage in STAGES:
        for res in RESOLUTIONS:
            con = m5.read_contact(stage, res)
            if not con:
                continue
            keys = [k for k in con if k in dist and k[0] != k[1]]
            keys = [k for k in keys if arm_type(k[0]) and arm_type(k[1])]
            cls = {}
            for k in keys:
                t = "".join(sorted([arm_type(k[0]), arm_type(k[1])]))
                cls.setdefault(t, []).append(con[k])
            # stage-resolved rho within arm-type classes, satellite-masked (>=50% sat
            # bins removed from the sequence axis)
            row = {}
            for t in ("pp", "pq", "qq"):
                kk = [k for k in keys if "".join(sorted([arm_type(k[0]),
                                                         arm_type(k[1])])) == t]
                row[t] = dict(n=len(kk),
                              mean_contact=round(float(np.mean([con[k] for k in kk])), 6)
                              if kk else None)
                if len(kk) >= 15:
                    xa = [max(0.0, seqjac_masked.get(k, 0.0)) for k in kk]
                    ya = [con[k] for k in kk]
                    o, p, *_ = m5.mantel(xa, ya, a.perms, rng)
                    row[t]["rho_masked_axis"] = round(o, 4)
                    row[t]["p_masked_axis"] = p
            # permutation test of the pp-vs-qq mean-contact contrast, permuting arm-type
            # labels among arms (preserves the arm-level dependence structure)
            all_arms = sorted({k[0] for k in keys} | {k[1] for k in keys})
            obs_pp = float(np.mean(cls.get("pp", [np.nan])))
            obs_qq = float(np.mean(cls.get("qq", [np.nan])))
            obs_d = obs_pp - obs_qq
            n_p = sum(1 for ar in all_arms if arm_type(ar) == "p")
            rng2 = np.random.default_rng(4242)
            cnt = 0
            null = []
            for _ in range(a.perms):
                lab = np.array(["q"] * len(all_arms))
                lab[rng2.choice(len(all_arms), n_p, replace=False)] = "p"
                lt = dict(zip(all_arms, lab))
                vv = {}
                for k in keys:
                    t = "".join(sorted([lt[k[0]], lt[k[1]]]))
                    vv.setdefault(t, []).append(con[k])
                if "pp" not in vv or "qq" not in vv:
                    continue
                d = float(np.mean(vv["pp"]) - np.mean(vv["qq"]))
                null.append(d)
                if abs(d) >= abs(obs_d):
                    cnt += 1
            row["pp_minus_qq"] = dict(observed=round(obs_d, 6),
                                      p_perm_armlabel=(cnt + 1) / (len(null) + 1),
                                      n_perms=len(null))
            pq[f"{stage}_{res//1000}kb"] = row

    # independent satellite annotation of the acrocentric/telocentric ends:
    # which end of each chromosome carries the satellite
    ends_tab = {}
    for (chrom, kind), d in sorted(end_sat.items()):
        ends_tab.setdefault(chrom, {})[kind] = round(d["sat_frac_union"], 4)

    # ---- verdicts (strictly derived from the numbers above)
    def peak_stage(tab, suffix, field="rho"):
        vals = {s: tab.get(f"{s}_{suffix}", {}).get(field) for s in STAGES}
        vals = {k: v for k, v in vals.items() if isinstance(v, (int, float))}
        if not vals:
            return None
        return max(vals, key=vals.get)

    verdict = {
        "peak_stage_20kb_baseline": peak_stage(base, "20kb"),
        "peak_stage_50kb_baseline": peak_stage(base, "50kb"),
        "peak_stage_20kb_low_sat(both<0.25)": peak_stage(
            {k.replace("_sat_low_25", ""): v for k, v in strat.items()
             if k.endswith("_sat_low_25")}, "20kb"),
        "peak_stage_20kb_high_sat(both>=0.25)": peak_stage(
            {k.replace("_sat_high_25", ""): v for k, v in strat.items()
             if k.endswith("_sat_high_25")}, "20kb"),
        "peak_stage_20kb_masked_axis": peak_stage(
            {k.replace("_axis_maskedSat", ""): v for k, v in strat.items()
             if k.endswith("_axis_maskedSat")}, "20kb"),
    }

    out = dict(
        item="R1-M6 mouse p-arm satellite confound",
        referee_verbatim=(
            "Mouse chromosomes are acro/telocentric with essentially no short arm, so a 1 Mb "
            "telomere-anchored 'p-arm' flank is largely composed of minor and major satellite, "
            "which are near-identical between chromosomes, placing these windows at once at "
            "the top of the Jaccard axis and at the maximum of the multi-mapping problem."),
        inputs=dict(flanks=str(FLANK_FA), phr=str(PHR_TSV), paf=str(PAF_DIR),
                    dist_matrix=str(DIST), repeat_library=str(RMSK_LIB),
                    blastn=str(BLASTN), bin=BIN,
                    mantel_functions=str(M5_PATH),
                    contact_tables="data/mouse_<stage>_phr_<res>bp_seqlevel.tsv"),
        satellite_annotation=dict(
            label=("major satellite GSAT_MM = Dfam-curated ANNOTATION; "
                   "minor satellite = PROXY (no Dfam 3.8 family named 'minor satellite' "
                   "for Mus in this library; MurSAT1+MMSAT4+CENSAT_MC used instead); "
                   "periodicity detector = PROXY"),
            families=famlen, family_headers=famhdr,
            hit_filters="blastn length>=100 bp, identity>=75%; bin covered if >=50%",
            per_flank={k: {kk: vv for kk, vv in v.items()
                           if kk not in ("union_cov", "proxy")}
                       for k, v in per_flank.items()},
            per_arm_sat_frac={k: round(v, 4) for k, v in sorted(arm_sat_frac.items())},
            per_arm_major={k: round(v, 4) for k, v in sorted(arm_sat_major.items())},
            per_arm_minor_proxy={k: round(v, 4) for k, v in sorted(arm_sat_minor.items())},
            per_chromosome_end_sat_frac=ends_tab,
        ),
        sequence_axis=dict(
            description=("arm-level sequence-sharing Jaccard re-derived from the wfmash "
                         "all-vs-all PAFs (>=95% identity blocks intersected with the called "
                         "PHR regions), unmasked and with >=50%-satellite 5 kb bins removed"),
            validation_vs_published_odgi_jaccard=val,
        ),
        mantel=dict(perms=a.perms, baseline_reproduced=base, stratified=strat,
                    n_pairs=pair_counts),
        p_vs_q=pq,
        verdicts=verdict,
        not_established=[
            "contact axis is NOT satellite-masked here: re-masking contacts needs the "
            "zuo2021_<stage>.mcool per-bin matrices (h5py/cooler absent from system python; "
            "0.3-0.5 Gpixel per stage) plus the paternal/maternal haplotype mapping",
            "no full RepeatMasker run was possible (RepeatMasker 4.1.7p1 is missing "
            "PubRef.pm/FastaDB.pm on this box); only the curated Dfam satellite consensus "
            "sequences were used, so non-satellite repeat content is unmeasured",
            "the published arm-level Jaccard is an odgi pangenome-graph node Jaccard and "
            "cannot be recomputed with satellite bins excluded without rebuilding the pggb "
            "graph; the masked-axis result uses the re-derived PAF-based Jaccard instead",
        ],
    )
    OUT.write_text(json.dumps(out, indent=2, sort_keys=False) + "\n")
    print("wrote", OUT)

    # ---- console summary
    print("\n== per-arm satellite fraction (5 kb bins, >=50% covered, Dfam-derived) ==")
    for k in sorted(arm_sat_frac):
        print(f"  {k:10s} union={arm_sat_frac[k]:.3f} major={arm_sat_major[k]:.3f} "
              f"minor(PROXY)={arm_sat_minor[k]:.3f}")
    print("\n== Mantel baseline (should reproduce scripts/hic/mantel_mouse_d_m5.json) ==")
    for k, v in base.items():
        print(f"  {k:16s} rho={v['rho']:+.3f} p={v['p_perm']:.4g} n={v['n_arm_pairs']}")
    print("\n== stratified (satellite fraction of flank pair, max of the two arms) ==")
    for k in sorted(strat):
        v = strat[k]
        if "rho" in v:
            print(f"  {k:34s} rho={v['rho']:+.3f} p={v['p_perm']:.4g} n={v['n_arm_pairs']}")
        else:
            print(f"  {k:34s} n={v.get('n_arm_pairs')} ({v.get('skipped','')})")
    print("\n== p vs q ==")
    for k, v in pq.items():
        pp, qq = v["pp"], v["qq"]
        print(f"  {k:14s} pp(n={pp['n']},mean={pp['mean_contact']},"
              f"rhoMasked={pp.get('rho_masked_axis')}) "
              f"qq(n={qq['n']},mean={qq['mean_contact']},"
              f"rhoMasked={qq.get('rho_masked_axis')}) "
              f"pp-qq={v['pp_minus_qq']['observed']:+.2e} "
              f"p={v['pp_minus_qq']['p_perm_armlabel']:.4g}")
    print("\n== verdicts ==")
    for k, v in verdict.items():
        print(f"  {k}: {v}")
    print("\nscratch kept at:", tmp)


if __name__ == "__main__":
    main()
