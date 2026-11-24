#!/usr/bin/env python3
"""Known-answer test for junction_kmers.py: a 60 bp segment present only in the query is found as
specific k-mers; reads carrying it forward and reverse-complemented are both counted; reads without it
give zero. Run: python3 scripts/conversion/test_junction_kmers.py WORKDIR"""
import os, random, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

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
# cooccur: A = k-mers of the insert's left half, B = of its right half; fragment r0 (one read) and
# pair p (left half in mate 1, right half in mate 2, as BAM) carry both, n0 carries only A.
import pysam
open(f"{wd}/a.tsv", "w").write("kmer\tx\n" + "".join(f"{x}\t0\n" for x in set(
    __import__("junction_kmers").canon_kmers(ins[:31], 31))))
open(f"{wd}/b.tsv", "w").write("kmer\tx\n" + "".join(f"{x}\t0\n" for x in set(
    __import__("junction_kmers").canon_kmers(ins[29:], 31))))
hdr = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "q", "LN": 2000}]})
with pysam.AlignmentFile(f"{wd}/frag.bam", "wb", header=hdr) as bam:
    for name, seq in (("r0", left[-40:] + ins + right[:40]), ("p", left[-60:] + ins[:31]), ("p", rc(ins[29:] + right[:60])),
                      ("n0", left[-60:] + ins[:31])):
        r = pysam.AlignedSegment(hdr)
        r.query_name, r.query_sequence, r.flag, r.reference_id, r.reference_start = name, seq, 4, -1, -1
        bam.write(r)
subprocess.run(py + ["cooccur", "--kmers-a", f"{wd}/a.tsv", "--kmers-b", f"{wd}/b.tsv", "--reads", f"{wd}/frag.bam",
                     "--out", f"{wd}/co.tsv"], check=True)
row = open(f"{wd}/co.tsv").read().splitlines()[1].split("\t")
ok2 = row[2:] == ["3", "2", "2"]
print(f"{'PASS' if ok2 else 'FAIL'} cooccur fragments with A, B, A and B: {row[2:]} (expect 3, 2, 2)")
sys.exit(0 if ok and ok2 else 1)
