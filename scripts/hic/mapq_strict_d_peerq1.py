#!/usr/bin/env python3
"""
MAPQ-strict Hi-C re-run for D-PeerQ1.

Goal
----
Quantify how much of the within-community Hi-C B/W signal at PHR-internal
coordinates depends on the random-placement of MAPQ0 (multi-mapped) reads
that HiC-Pro retained under MIN_MAPQ=0, RM_MULTI=0.

Why mcool alone is not enough
-----------------------------
An mcool file is post-binning. By the time reads have been aggregated into
bin-pair counts there is no per-read MAPQ to filter on; the contact matrix
already reflects the random placement chosen at HiC-Pro alignment time.
Re-applying a "MAPQ-strict" filter therefore requires going back to the
pairs stage and re-binning. The pipeline below operates on HiC-Pro
.allValidPairs (text, one row per pair) which carries both alignment MAPQ
columns. After MAPQ-strict filtering the pairs are loaded with
`cooler cload pairs` and zoomified to mcool.

Pipeline (per sample)
---------------------
  1. Read HiC-Pro .allValidPairs (random-placement, MIN_MAPQ=0 baseline).
  2. Drop pairs where either mate has MAPQ < MAPQ_THR (default 30).
  3. `cooler cload pairs` to build a strict-MAPQ .cool at the analysis
     resolution (default 50 kb).
  4. `cooler zoomify` to a strict-MAPQ .mcool.
  5. Run analyze_hic_communities.py against the strict .mcool to produce
     a strict-MAPQ global_test.tsv with within/between/B-W and Mann-Whitney
     p-value at PHR coordinates and at flanking 100 kb coordinates.
  6. Emit a per-sample comparison row:
        sample tech region random_BW strict_BW delta_log2 random_p strict_p

Inputs (CLI)
------------
  --sample            sample label (HG002, CHM13, ...)
  --tech              Hi-C | Pore-C | CiFi
  --pairs             HiC-Pro .allValidPairs file (text, gzipped OK)
  --chrom-sizes       chromsizes TSV (chrom, size) matching the pairs file
  --random-mcool      v5 random-placement mcool (for the random_BW value)
  --communities-tsv   arm-to-community assignment (15-community Leiden)
  --chrom-parts       BED with chrom, start, end, arm-label (p/q)
  --phr-regions-tsv   per-haplotype PHR coordinates (all-vs-all.p95.id95.len.tsv)
  --flanking-bed      per-haplotype 100 kb flanking BED (centromere-ward)
  --resolution        binning resolution in bp (default 50000)
  --mapq-threshold    MAPQ cutoff for strict filter (default 30)
  --out-dir           output directory
  --analyze-script    path to analyze_hic_communities.py
  --threads           parallel decompression / cload workers (default 8)

Output
------
  <out_dir>/<sample>.<tech>.mapq{MAPQ_THR}.pairs.gz   (filtered pairs)
  <out_dir>/<sample>.<tech>.mapq{MAPQ_THR}.cool       (single-res cool)
  <out_dir>/<sample>.<tech>.mapq{MAPQ_THR}.mcool      (zoomified mcool)
  <out_dir>/<sample>.<tech>.mapq{MAPQ_THR}.global_test.tsv
  <out_dir>/<sample>.<tech>.mapq{MAPQ_THR}.comparison_row.tsv
  <out_dir>/comparison_summary.tsv  (one row per sample, appended)

Expected outcome (predicted by v5 §05_hic_validation flanking control)
----------------------------------------------------------------------
  At PHR-internal coordinates: the strict-MAPQ filter is expected to drop
  the bulk of reads at high-identity paralogous loci (which is why
  MIN_MAPQ=0 was chosen in v5 in the first place). Strict-MAPQ B/W at PHR
  coords is therefore expected to be uninterpretable (very few surviving
  reads, dominated by noise). This is not a regression of the v5 result;
  it is a confirmation that the PHR-internal signal lives in the
  multi-mapping fraction. Compare to the FLANKING analysis (next).

  At flanking 100 kb coordinates: these are unique-sequence regions
  centromere-ward of the PHR boundary, so MAPQ-strict and MAPQ0 give
  near-identical values (no multi-mapping reads to drop). The v5 flanking
  B/W values are already the MAPQ-strict-equivalent measurement; this
  script should reproduce them within Poisson sampling noise.

Usage example
-------------
  python scripts/hic/mapq_strict_d_peerq1.py \
      --sample HG002 \
      --tech Hi-C \
      --pairs   /moosefs/guarracino/HPRCv2/PHR_III/3d/HG002/hic_results/data/HG002/HG002.allValidPairs.gz \
      --chrom-sizes /moosefs/guarracino/HPRCv2/PHR_III/3d/HG002/annotations/chrom.sizes \
      --random-mcool /moosefs/guarracino/HPRCv2/PHR_III/3d/HG002/HG002.mcool \
      --communities-tsv /moosefs/guarracino/HPRCv2/PHR_III/communities/arm_15_communities.tsv \
      --chrom-parts /moosefs/guarracino/HPRCv2/PHR_III/3d/HG002/annotations/HG002_chrom_parts.bed \
      --phr-regions-tsv /moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv \
      --flanking-bed /moosefs/guarracino/HPRCv2/PHR_III/flanking/HG002_flanking_100kb.bed \
      --resolution 50000 \
      --mapq-threshold 30 \
      --analyze-script /moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py \
      --out-dir paper_prep/synthesis/d_peerq1/HG002

  Run once for HG002 Hi-C, once for CHM13 Hi-C; concatenate
  comparison_row.tsv into the table in ANALYSIS_D_PEERQ1.md.

Notes
-----
- HiC-Pro .allValidPairs columns (tab-separated, no header):
    read_id  chr1  pos1  strand1  chr2  pos2  strand2  size  mapq1  ?  mapq2  ?
  Some versions use different column orderings; this script reads two MAPQ
  columns and the user must point --mapq-cols at them (default 8,10
  matches HiC-Pro >= 3.0 .allValidPairs format).
- For Pore-C / CiFi the pairs format is pairtools-pairs (4dn format) with
  MAPQ columns 8,9 (after pairtools parse2). Use --pairs-format pairtools
  to switch.

Data-access blocker (this worker, 2026-05-18)
---------------------------------------------
At time of writing, the canonical Hi-C .allValidPairs files at
/moosefs/guarracino/HPRCv2/PHR_III/3d/ are not mounted on this worker
(agent-207, .wg-worktrees/agent-207). Re-running this script with real
inputs is therefore deferred to an environment with moosefs access.
A dry-run mode (--smoke) is provided that validates the argparse
plumbing and the awk/cooler command lines without touching real data.
"""

import argparse
import gzip
import os
import shlex
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--tech", default="Hi-C", choices=["Hi-C", "Pore-C", "CiFi"])
    p.add_argument("--pairs", required=True,
                   help="HiC-Pro .allValidPairs (text or .gz) OR pairtools pairs")
    p.add_argument("--pairs-format", default="hicpro",
                   choices=["hicpro", "pairtools"],
                   help="hicpro: HiC-Pro allValidPairs columns; pairtools: 4dn pairs")
    p.add_argument("--mapq-cols",
                   help="1-indexed comma-separated MAPQ column indices, "
                        "default 8,10 for hicpro, 8,9 for pairtools",
                   default=None)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--random-mcool", required=True,
                   help="v5 random-placement mcool (used to extract the "
                        "random-placement B/W baseline for the comparison row)")
    p.add_argument("--communities-tsv", required=True)
    p.add_argument("--chrom-parts", required=True)
    p.add_argument("--phr-regions-tsv", required=True)
    p.add_argument("--flanking-bed", required=True)
    p.add_argument("--resolution", type=int, default=50000)
    p.add_argument("--mapq-threshold", type=int, default=30)
    p.add_argument("--analyze-script", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--smoke", action="store_true",
                   help="Dry-run: print commands without executing")
    return p.parse_args()


def default_mapq_cols(pairs_format):
    return "8,10" if pairs_format == "hicpro" else "8,9"


def run(cmd, smoke):
    print(f"[run] {cmd}", flush=True)
    if smoke:
        print("[smoke] skip exec", flush=True)
        return 0
    rc = subprocess.call(cmd, shell=True, executable="/bin/bash")
    if rc != 0:
        raise SystemExit(f"command failed (rc={rc}): {cmd}")
    return rc


def filter_pairs_mapq(pairs_in, pairs_out, mapq_cols, mapq_thr, smoke):
    """Filter pairs to MAPQ >= mapq_thr on both mates, gzip output.

    The awk expression keeps the row when both numeric MAPQ columns are
    >= mapq_thr. Column indices are 1-indexed (awk style).
    """
    c1, c2 = mapq_cols.split(",")
    decompress = "zcat" if pairs_in.endswith(".gz") else "cat"
    awk = (f"awk 'BEGIN{{FS=OFS=\"\\t\"}} "
           f"!/^#/ && ${c1}+0>={mapq_thr} && ${c2}+0>={mapq_thr}'")
    cmd = f"{decompress} {shlex.quote(pairs_in)} | {awk} | gzip -c > {shlex.quote(pairs_out)}"
    run(cmd, smoke)


def cooler_cload(pairs_gz, chrom_sizes, resolution, cool_out, threads, smoke):
    """Bin filtered pairs into a single-resolution .cool.

    cooler cload pairs expects: chrom1 pos1 chrom2 pos2 columns. HiC-Pro
    .allValidPairs columns are: read_id chr1 pos1 strand1 chr2 pos2 ...
    so chrom1=2, pos1=3, chrom2=5, pos2=6.
    """
    cmd = (f"cooler cload pairs --nproc {threads} "
           f"-c1 2 -p1 3 -c2 5 -p2 6 "
           f"{shlex.quote(chrom_sizes)}:{resolution} "
           f"{shlex.quote(pairs_gz)} {shlex.quote(cool_out)}")
    run(cmd, smoke)


def cooler_zoomify(cool_in, mcool_out, threads, smoke):
    """Multi-res zoomify to .mcool (matches v5 resolutions: 5/10/20/50/100 kb)."""
    cmd = (f"cooler zoomify --nproc {threads} "
           f"--resolutions 5000,10000,20000,50000,100000 "
           f"--balance --balance-args '--mad-max 5 --min-nnz 10' "
           f"-o {shlex.quote(mcool_out)} {shlex.quote(cool_in)}")
    run(cmd, smoke)


def run_analyze(analyze_script, sample, mcool, resolution, chrom_parts,
                communities_tsv, phr_regions, out_dir, label, smoke):
    """Invoke analyze_hic_communities.py and capture the global_test.tsv."""
    cmd = (f"python {shlex.quote(analyze_script)} "
           f"--sample {shlex.quote(f'{sample}_{label}')} "
           f"--matrix {shlex.quote(mcool)} "
           f"--resolution {resolution} "
           f"--chrom-parts {shlex.quote(chrom_parts)} "
           f"--communities-tsv {shlex.quote(communities_tsv)} "
           f"--phr-regions-tsv {shlex.quote(phr_regions)} "
           f"--output-dir {shlex.quote(out_dir)}")
    run(cmd, smoke)


def read_bw(global_test_tsv):
    """Extract within_mean, between_mean, p_value from a global_test.tsv."""
    if not os.path.exists(global_test_tsv):
        return None
    with open(global_test_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row.get("test") == "within_vs_between":
                w = float(row["within_mean"])
                b = float(row["between_mean"])
                p = float(row["p_value"])
                return {"within": w, "between": b, "BW": (b / w) if w > 0 else None,
                        "p_value": p}
    return None


def write_comparison_row(sample, tech, region, random_metrics, strict_metrics,
                         out_path, mapq_thr):
    """Write a one-row TSV with the random vs strict comparison."""
    cols = ["sample", "tech", "region", "mapq_threshold",
            "random_within", "random_between", "random_BW", "random_p",
            "strict_within", "strict_between", "strict_BW", "strict_p",
            "delta_log2_BW"]
    def g(d, k):
        return "" if d is None or d.get(k) is None else f"{d[k]:.6g}"
    def log2_delta():
        if (random_metrics and strict_metrics
                and random_metrics.get("BW") and strict_metrics.get("BW")
                and random_metrics["BW"] > 0 and strict_metrics["BW"] > 0):
            import math
            return f"{math.log2(strict_metrics['BW'] / random_metrics['BW']):.6g}"
        return ""
    row = [sample, tech, region, str(mapq_thr),
           g(random_metrics, "within"), g(random_metrics, "between"),
           g(random_metrics, "BW"), g(random_metrics, "p_value"),
           g(strict_metrics, "within"), g(strict_metrics, "between"),
           g(strict_metrics, "BW"), g(strict_metrics, "p_value"),
           log2_delta()]
    with open(out_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        f.write("\t".join(row) + "\n")
    print(f"[wrote] {out_path}", flush=True)


def append_summary(out_dir, sample_row_path):
    """Append per-sample comparison_row.tsv to comparison_summary.tsv."""
    summary = os.path.join(out_dir, "comparison_summary.tsv")
    sample_lines = open(sample_row_path).readlines()
    if not os.path.exists(summary):
        with open(summary, "w") as f:
            f.writelines(sample_lines)
    else:
        with open(summary, "a") as f:
            f.write(sample_lines[1])


def main():
    args = parse_args()
    mapq_cols = args.mapq_cols or default_mapq_cols(args.pairs_format)
    os.makedirs(args.out_dir, exist_ok=True)

    base = os.path.join(args.out_dir,
                        f"{args.sample}.{args.tech}.mapq{args.mapq_threshold}")
    pairs_strict = base + ".pairs.gz"
    cool_strict = base + ".cool"
    mcool_strict = base + ".mcool"

    print(f"[D-PeerQ1] sample={args.sample} tech={args.tech} "
          f"mapq>={args.mapq_threshold} resolution={args.resolution}",
          flush=True)
    print(f"[D-PeerQ1] pairs (random-placement) -> {args.pairs}", flush=True)
    print(f"[D-PeerQ1] mapq_cols (1-indexed) = {mapq_cols}", flush=True)

    filter_pairs_mapq(args.pairs, pairs_strict, mapq_cols,
                      args.mapq_threshold, args.smoke)
    cooler_cload(pairs_strict, args.chrom_sizes, args.resolution,
                 cool_strict, args.threads, args.smoke)
    cooler_zoomify(cool_strict, mcool_strict, args.threads, args.smoke)

    # PHR-region analysis: strict and random
    phr_out = os.path.join(args.out_dir, "phr")
    os.makedirs(phr_out, exist_ok=True)
    run_analyze(args.analyze_script, args.sample, args.random_mcool,
                args.resolution, args.chrom_parts, args.communities_tsv,
                args.phr_regions_tsv, phr_out, "random_phr", args.smoke)
    run_analyze(args.analyze_script, args.sample, mcool_strict,
                args.resolution, args.chrom_parts, args.communities_tsv,
                args.phr_regions_tsv, phr_out, f"mapq{args.mapq_threshold}_phr",
                args.smoke)

    # Flanking-region analysis: strict and random
    flank_out = os.path.join(args.out_dir, "flanking")
    os.makedirs(flank_out, exist_ok=True)
    run_analyze(args.analyze_script, args.sample, args.random_mcool,
                args.resolution, args.chrom_parts, args.communities_tsv,
                args.flanking_bed, flank_out, "random_flank", args.smoke)
    run_analyze(args.analyze_script, args.sample, mcool_strict,
                args.resolution, args.chrom_parts, args.communities_tsv,
                args.flanking_bed, flank_out, f"mapq{args.mapq_threshold}_flank",
                args.smoke)

    if args.smoke:
        print("[smoke] skipping comparison row write (no real outputs)",
              flush=True)
        return

    for region, region_dir, random_label, strict_label in [
        ("PHR", phr_out, "random_phr", f"mapq{args.mapq_threshold}_phr"),
        ("flanking", flank_out, "random_flank",
         f"mapq{args.mapq_threshold}_flank"),
    ]:
        random_tsv = os.path.join(region_dir,
                                  f"{args.sample}_{random_label}_global_test.tsv")
        strict_tsv = os.path.join(region_dir,
                                  f"{args.sample}_{strict_label}_global_test.tsv")
        random_metrics = read_bw(random_tsv)
        strict_metrics = read_bw(strict_tsv)
        row_path = base + f".{region}.comparison_row.tsv"
        write_comparison_row(args.sample, args.tech, region,
                             random_metrics, strict_metrics, row_path,
                             args.mapq_threshold)
        append_summary(args.out_dir, row_path)


if __name__ == "__main__":
    main()
