#!/usr/bin/env python3
"""Known-answer test for junction_kmers.py: a 60 bp segment present only in the query is found as
specific k-mers; reads carrying it forward and reverse-complemented are both counted; reads without it
give zero. Run: python3 scripts/conversion/test_junction_kmers.py WORKDIR"""
import os, random, subprocess, sys

wd = sys.argv[1]
os.makedirs(wd, exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(3)
rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
left, right = rnd(500), rnd(500)
ins = rnd(60)
while ins[0] == right[0] or ins[-1] == left[-1]:  # otherwise the edge k-mers also exist in left+right
    ins = rnd(60)
open(f"{wd}/q.fa", "w").write(f">q\n{left}{ins}{right}\n")
open(f"{wd}/ex.fa", "w").write(f">ex\n{left}{right}\n")
reads = [left[-80:] + ins + right[:80], rc(left[-70:] + ins + right[:70]), left[:200], right[100:300]]
open(f"{wd}/reads.fa", "w").write("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
open(f"{wd}/neg.fa", "w").write("".join(f">n{i}\n{r}\n" for i, r in enumerate(reads[2:])))
py = [sys.executable, os.path.join(here, "junction_kmers.py")]
subprocess.run(py + ["specific", "--query-fasta", f"{wd}/q.fa", "--exclude-fasta", f"{wd}/ex.fa", "--k", "31",
                     "--out", f"{wd}/spec.tsv"], check=True)
nspec = len(open(f"{wd}/spec.tsv").read().splitlines()) - 1
subprocess.run(py + ["count", "--kmers", f"{wd}/spec.tsv", "--reads", f"{wd}/reads.fa", f"{wd}/neg.fa",
                     "--out", f"{wd}/count.tsv"], check=True)
rows = [l.split("\t") for l in open(f"{wd}/count.tsv").read().splitlines()[1:]]
ok = nspec == 60 + 31 - 1 and int(rows[0][4]) == 2 * nspec and int(rows[1][4]) == 0
print(f"{'PASS' if ok else 'FAIL'} specific k-mers {nspec} (expect 90); occurrences in 2 carrier reads "
      f"{rows[0][4]} (expect {2 * nspec}); in non-carriers {rows[1][4]} (expect 0)")
sys.exit(0 if ok else 1)
