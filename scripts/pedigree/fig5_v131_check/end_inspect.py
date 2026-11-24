#!/usr/bin/env python3
"""What do the reads and the Illumina k-mers say about one contig end? (ends with many merqury error k-mers)

usage: end_inspect.py SAMTOOLS NORM_FA FINAL_FA MERQ_ONLY_BED_GZ CONTIG END(p|q) ONT_BAM[,...] HIFI_BAM[,...] > TSV
Distances are bp from the tip. FINAL_FA (the corrected assembly, scored by merqury) gives error k-mers and the telomere
motif (CCCTAA at p, TTAGGG at q); NORM_FA (the assembly the BAMs are aligned to) gives the read data. The last 25 kb must
be identical in both, otherwise exit non-zero. Per 1 kb bin from the tip (0-20 kb), per technology (primary alignments):
  depth = aligned bases / 1,000; mapq0 = reads with MAPQ 0 overlapping the bin;
  stop = reads whose aligned part stops in the bin on the tip side; clip = of those, soft clip >= 500 bp on the tip side
  (the read goes on but does not match the assembly); clip_tel = of those clips, >= 50% motif (the read goes on in telomere);
  tel_end = of all stopping reads, the last 1 kb of the read on the tip side is >= 50% motif (the read ends in telomere);
  past = of the clips, the soft clip is longer than the distance to the tip (the read reaches past the contig end).
out: contig, end, tech, bin_kb, err_kmers, motif_pct, depth, mapq0, stop, clip, clip_tel, tel_end, past
"""
import collections, gzip, re, subprocess, sys
S, nf, ff, bed, C, E, ont, hifi = sys.argv[1:9]; W, N = 1000, 20
m = "CCCTAA" if E == "p" else "TTAGGG"
def log(x): print(x, file=sys.stderr, flush=True)
def L(fa): return next(int(l.split("\t")[1]) for l in open(fa + ".fai") if l.split("\t")[0] == C)
def get(fa, a, b): return "".join(subprocess.run([S, "faidx", fa, f"{C}:{a}-{b}"], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
Ln, Lf = L(nf), L(ff); K = 25000
en = get(nf, 1, K) if E == "p" else get(nf, Ln - K + 1, Ln); ef = get(ff, 1, K) if E == "p" else get(ff, Lf - K + 1, Lf)
log(f"{C} {E}: norm length {Ln}, final length {Lf}, last {K} bp identical: {en == ef}")
if en != ef: sys.exit(f"ERROR {C} {E}: the {K} bp end differs between {nf} and {ff}")
def dist(x, L0): return x - 1 if E == "p" else L0 - x  # 1-based position -> bp from the tip
seq = ef if E == "p" else ef[::-1]; mm = m if E == "p" else m[::-1]
mot = [round(100 * seq[k * W:(k + 1) * W].count(mm) * 6 / W) for k in range(N)]
err = collections.Counter()
for l in gzip.open(bed, "rt"):
    f = l.split("\t")
    if f[0] == C:
        d = dist(int(f[1]) + 1, Lf)
        if 0 <= d < N * W: err[d // W] += 1
log(f"error k-mers in the first {N} kb: {sum(err.values())}")
reg = f"{C}:1-{N * W + 5000}" if E == "p" else f"{C}:{Ln - N * W - 5000}-{Ln}"
for tech, bams in (("ont", ont), ("hifi", hifi)):
    dep, mq0, stop, clip, ctel, tel, past = (collections.Counter() for _ in range(7)); n = 0
    for b in bams.split(","):
        for l in subprocess.run([S, "view", "-F", "0x904", b, reg], capture_output=True, text=True, check=True).stdout.splitlines():
            f = l.split("\t"); r0 = r = int(f[3]); cig = re.findall(r"(\d+)([MIDNSHP=X])", f[5]); n += 1
            for k, op in cig:
                if op in "M=XDN": r += int(k)
            a, z = (dist(r0, Ln), dist(r - 1, Ln)) if E == "p" else (dist(r - 1, Ln), dist(r0, Ln))  # a: tip side
            for k in range(max(0, a // W), min(N, z // W + 1)):
                dep[k] += max(0, min(z + 1, (k + 1) * W) - max(a, k * W)) / W
                if f[4] == "0": mq0[k] += 1
            sc = int(cig[0][0]) if E == "p" and cig[0][1] == "S" else int(cig[-1][0]) if E == "q" and cig[-1][1] == "S" else 0
            s = f[9]; tip = s[:sc] if E == "p" else s[len(s) - sc:]; t1 = s[:1000] if E == "p" else s[-1000:]
            k = a // W
            if k >= N: continue
            stop[k] += 1; tel[k] += t1.count(m) * 6 >= 0.5 * len(t1)
            if sc >= 500: clip[k] += 1; ctel[k] += tip.count(m) * 6 >= 0.5 * len(tip); past[k] += sc > a
    log(f"{tech}: {n} primary reads in {reg}")
    for k in range(N):
        print(f"{C}\t{E}\t{tech}\t{k}\t{err[k]}\t{mot[k]}\t{dep[k]:.1f}\t{mq0[k]}\t{stop[k]}\t{clip[k]}\t{ctel[k]}\t{tel[k]}\t{past[k]}")
