#!/usr/bin/env python3
"""Fig5 drawn runs under alternative run-calling rules, v1.1 vs v1.3.1.

Rules (all keep the paper filter: run bp >= 10 kb, mean inter identity >= 0.95):
  paper      same donor haplotype+chrom, query and donor touching (gap 0)
  gap1k      same donor haplotype+chrom, query and donor gap <= 1 kb
  gap2k      same donor haplotype+chrom, query and donor gap <= 2 kb
  gap3k      same donor haplotype+chrom, query and donor gap <= 3 kb
  gap4k      same donor haplotype+chrom, query and donor gap <= 4 kb
  anyhap     same donor chrom (h1 or h2), query gap 0, donor ignored
  anyhap4k   same donor chrom (h1 or h2), query gap <= 4 kb, donor ignored

usage: run_sensitivity.py PLOTTER.py LABEL=CLASS_WINNERS.tsv.gz [...]
"""
import importlib.util, sys
from collections import defaultdict

spec = importlib.util.spec_from_file_location("plotter", sys.argv[1])
P = importlib.util.module_from_spec(spec); sys.modules["plotter"] = P; spec.loader.exec_module(P)
print("label\trule\tcategory\tquery\tbp\twindows\tdonor\tmean_same\tmean_inter")
FOCUS = ["chr9_chr3_candidate", "chr5_chr1_candidate", "chr21_chr4_candidate", "PAR_XY"]


def query_only_runs(segs, gap):
    """Runs on query contiguity, any donor haplotype"""
    by = defaultdict(list)
    for s in segs:
        by[(s.query_seq, s.target_chrom)].append(s)
    runs = []
    for key, ss in by.items():
        cur = None
        for s in sorted(ss, key=lambda s: s.query_start):
            if cur and s.query_start - cur.query_end <= gap:
                P.extend_run(cur, s); cur.haps.add(s.donor_haplotype)
            else:
                cur = P.start_run(s, "run", defaultdict(int)); cur.haps = {s.donor_haplotype}; runs.append(cur)
    return runs


def drawn(runs):
    return [r for r in runs if r.bp >= 10_000 and r.mean_inter_identity >= 0.95]


for arg in sys.argv[2:]:
    label, path = arg.split("=", 1)
    segs = P.read_segments(path)
    rules = {}
    P.CONTIGUOUS_MERGE_GAP_BP = 0; rules["paper"] = P.group_runs(segs)
    for g in (1, 2, 3, 4):
        P.CONTIGUOUS_MERGE_GAP_BP = g * 1000; rules[f"gap{g}k"] = P.group_runs(segs)
    P.CONTIGUOUS_MERGE_GAP_BP = 0
    rules["anyhap"] = query_only_runs(segs, 0)
    rules["anyhap4k"] = query_only_runs(segs, 4000)
    sys.stderr.write(f"{label}: segments={len(segs)} " + " ".join(f"{n}={len(drawn(r))}" for n, r in rules.items()) + "\n")
    for name, runs in rules.items():
        for r in sorted((r for r in drawn(runs) if r.category in FOCUS), key=lambda r: (r.category, r.query_start)):
            haps = ",".join(sorted(getattr(r, "haps", {r.donor_haplotype})))
            print(f"{label}\t{name}\t{r.category}\t{r.query_chrom}:{r.query_start}-{r.query_end}\t{r.bp}\t{r.windows}\t{r.target_chrom}:{haps}\t{r.mean_same_identity:.6f}\t{r.mean_inter_identity:.6f}")
