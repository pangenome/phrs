#!/usr/bin/env python3
"""Data for the images of the two chromosome-end fixes (v1.3.1): DAD h2 chr9q and SELF.m chr22q. Aggregates only.

usage: fixes_figdata.py SAMTOOLS MINIMAP2 D OUT_DIR   (run in D = washu_v131; OUT_DIR new)
chr9q: SELF.p chr9 (PAN027#2#chr9 of g1/final_v131d, the copy SELF got from DAD) last 400 kb in 10 kb windows, aligned
  (minimap2 -cx asm5) to the last 700 kb of DAD h2 chr9 in norm_v131 (original) and g1/final_v131c (corrected); per window
  the best record: identity = matching bases / alignment block length, query bases aligned.
chr22q: SELF.m chr22 (PAN027#1#chr22), last 30 kb of norm_v131 in 500 bp bins (bin 0 = the original contig end):
  TTAGGG % of the bin; ONT read ends (primary reads that cross 30 kb before the end, end projected with soft clips as in
  tip_trim.py); merqury error k-mer starts in the original (merqury/v131) and in the corrected assembly (g1v131d, same
  coordinates up to the cut at KEEP_TO, checked).
out: OUT_DIR/chr9q_windows.tsv, OUT_DIR/chr22q_bins.tsv; log on stderr.
"""
import glob, gzip, os, re, subprocess, sys
S, MM, D, O = sys.argv[1:5]; os.chdir(D); os.makedirs(O, exist_ok=False)
def log(m): print(m, file=sys.stderr, flush=True)
def fai(fa, c): return next(int(l.split("\t")[1]) for l in open(fa + ".fai") if l.split("\t")[0] == c)
def get(fa, c, a, b): return "".join(subprocess.run([S, "faidx", fa, f"{c}:{a + 1}-{b}"], capture_output=True, text=True, check=True).stdout.split("\n")[1:])
# chr9q
q, Lq = "g1/final_v131d/PAN027/norm/PAN027#2.fa.gz", fai("g1/final_v131d/PAN027/norm/PAN027#2.fa.gz", "PAN027#2#chr9")
with open(f"{O}/q.fa", "w") as f:
    for a in range(Lq - 400000, Lq, 10000): f.write(f">w{Lq - a}\n{get(q, 'PAN027#2#chr9', a, a + 10000)}\n")
best = {}
for lab, t in (("original", "norm_v131/PAN011#2.fa.gz"), ("corrected", "g1/final_v131c/PAN011/norm/PAN011#2.fa.gz")):
    Lt = fai(t, "PAN011#2#chr9")
    with open(f"{O}/t.fa", "w") as f: f.write(f">t\n{get(t, 'PAN011#2#chr9', Lt - 700000, Lt)}\n")
    out = subprocess.run([MM, "-cx", "asm5", "-t", "4", f"{O}/t.fa", f"{O}/q.fa"], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        x = l.split("\t"); w = int(x[0][1:]); k = (lab, w)
        if k not in best or int(x[9]) > best[k][0]: best[k] = (int(x[9]), int(x[10]), int(x[3]) - int(x[2]))
    log(f"chr9q {lab}: DAD h2 chr9 length {Lt}, windows aligned {sum(1 for k in best if k[0] == lab)} of 40")
with open(f"{O}/chr9q_windows.tsv", "w") as f:
    f.write("window_end_from_tip_bp\tidentity_original\tquery_aligned_original\tidentity_corrected\tquery_aligned_corrected\n")
    for w in range(10000, 400001, 10000):
        r = [best.get((lab, w)) for lab in ("original", "corrected")]
        f.write(f"{w}\t" + "\t".join(f"{m / b:.5f}\t{qa}" if m else "NA\t0" for m, b, qa in (x or (0, 1, 0) for x in r)) + "\n")
os.remove(f"{O}/q.fa"); os.remove(f"{O}/t.fa")
# chr22q
C, K, nf, ff = "PAN027#1#chr22", 48463777, "norm_v131/PAN027#1.fa.gz", "g1/final_v131d/PAN027/norm/PAN027#1.fa.gz"
L, Lf = fai(nf, C), fai(ff, C)
if Lf != K or get(ff, C, K - 2000, K) != get(nf, C, K - 2000, K): sys.exit(f"ERROR corrected {C} is not the original cut at {K}")
W, N = 500, 60; s = get(nf, C, L - W * N, L).upper()
tel = [s[len(s) - (i + 1) * W:len(s) - i * W].count("TTAGGG") * 6 * 100 / W for i in range(N)]
A, ends = L - 30000, [0] * N; past = 0
for b in sorted(glob.glob("ont_full/*.PAN027.v131.bam")):
    for l in subprocess.run([S, "view", "-F", "0x904", b, f"{C}:{A}-{A}"], capture_output=True, text=True, check=True).stdout.splitlines():
        x = l.split("\t"); r, qp = int(x[3]), 0
        for k, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]):
            k = int(k)
            if op in "M=X":
                if r + k > A: qp += max(0, A - r); break
                r += k; qp += k
            elif op in "DN": r += k
            elif op in "IS": qp += k
        else: continue
        e = A + len(x[9]) - qp; i = (L - e) // W
        if e > L: past += 1
        elif 0 <= i < N: ends[i] += 1
def km(tag, Lc):
    v = [0] * N
    for l in gzip.open(f"merqury/{tag}/PAN027/PAN027_hap1_only.bed.gz", "rt"):
        x = l.split("\t")
        if x[0] == C:
            i = (L - 1 - int(x[1])) // W
            if 0 <= i < N: v[i] += 1
    return v
ko, kf = km("v131", L), km("g1v131d", Lf)
with open(f"{O}/chr22q_bins.tsv", "w") as f:
    f.write("bin_from_original_end\tbin_start_bp_from_end\tTTAGGG_pct\tont_read_ends\tmerqury_kmers_original\tmerqury_kmers_corrected\tremoved_by_trim\n")
    for i in range(N): f.write(f"{i}\t{i * W}\t{tel[i]:.1f}\t{ends[i]}\t{ko[i]}\t{kf[i]}\t{int((i + 1) * W <= L - K)}\n")
log(f"chr22q: original length {L}, cut at {K} ({L - K} bp removed); ONT reads crossing {A}: {sum(ends) + past} ({past} end past the original end, {sum(ends)} in the last {W * N} bp); merqury k-mers in the last {W * N} bp: {sum(ko)} -> {sum(kf)}")
