#!/usr/bin/env python3
"""Incremental PHR call for one telomere-anchored flank (prototype, to compare with the all-vs-all calls).

The all-vs-all caller (scripts/phr-detection/find-multichr-regions-incremental.py) counts, per 5 kb window, the impg
records (>= --min-identity, >= --min-output-length) per partner chromosome on a bidirectional index of all all-vs-one
PAFs, and a window passes when >= --min-chrs chromosomes have >= --min-count records (the window's own record counts
once for its own chromosome). Counts only grow when partners are added, so a window that passes stays passing.
Here we add partners one haplotype at a time (seeded random order): each round aligns the windows still undecided
(merged, padded by --pad on both sides) against the flanks of that haplotype in both directions, as the all-vs-one runs
do (partner flank as target with our pieces as queries, one wfmash run per partner; our pieces as target with all its
partner flanks as queries), lifts the coordinates back onto the flank, indexes the round with impg, adds the records of
each undecided window, and drops the windows that pass. After the last round the caller's rule gives the region: from
the telomere to the last passing window before the first run of --stop failing windows.
Output: one TSV row per window (pass, round that decided it) and a summary line; verbose log on stderr."""
import argparse, collections, os, random, re, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--end", required=True, help="flank name, e.g. HG00097#1#CM094072.1:4314-504313_chr21_parm")
ap.add_argument("--flank-dir", required=True, help="one FASTA per flank (<sanitized name>.fa, with .fai) and haplotypes.tsv")
ap.add_argument("--wfmash", required=True)
ap.add_argument("--impg", required=True)
ap.add_argument("--samtools", required=True)
ap.add_argument("--seed", type=int, default=20261125)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--pad", type=int, default=20000)
ap.add_argument("--window", type=int, default=5000)
ap.add_argument("--max", type=int, default=500000)
ap.add_argument("--min-identity", default="0.95")
ap.add_argument("--min-output-length", default="3000")
ap.add_argument("--min-count", type=int, default=5)
ap.add_argument("--min-chrs", type=int, default=2)
ap.add_argument("--stop", type=int, default=4)
ap.add_argument("--mode", choices=["frontier", "all"], default="frontier", help="frontier: only the windows the caller would scan; all: every window")
ap.add_argument("--max-rounds", type=int, default=0, help="stop after this many haplotypes (pilot); 0 = all")
ap.add_argument("--out", required=True)
a = ap.parse_args()
T = a.end
san = lambda s: s.translate(str.maketrans(":@#-", "____"))
chrom = lambda s: s.rsplit("_", 2)[-2]
t0 = time.time()
log = lambda m: (sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {T}: {m}\n"), sys.stderr.flush())  # python 3.7 buffers stderr in a file


def run(cmd, out=None):
    with open(out, "w") if out else open(os.devnull, "w") as fh:
        r = subprocess.run(cmd, stdout=fh if out else subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise SystemExit(f"FATAL: {' '.join(cmd)} exit {r.returncode}: {r.stderr[-1500:]}")
    return r.stdout


haps = collections.defaultdict(list)
for l in open(os.path.join(a.flank_dir, "haplotypes.tsv")):
    h, s = l.rstrip("\n").split("\t")
    haps[h].append(s)
tfa = os.path.join(a.flank_dir, san(T) + ".fa")
L = int(open(tfa + ".fai").readline().split("\t")[1])
own = chrom(T)
order = sorted(haps)
random.Random(a.seed).shuffle(order)
if a.max_rounds:
    order = order[:a.max_rounds]
mx = min(a.max, L)
wins = [(i * a.window, i * a.window + a.window) for i in range((mx - a.window) // a.window + 1)]
cnt = {w: collections.Counter({own: 1}) for w in wins}  # the window's own record
passed, hidx = {}, {w: 0 for w in wins}  # haplotypes already checked per window
tel = wins if T.endswith("_parm") else wins[::-1]  # telomere order
work_bp, rounds, H = 0, 0, len(order)
log(f"length {L}, {len(wins)} windows, {H} haplotypes, seed {a.seed}, pad {a.pad}, mode {a.mode}, own chromosome {own}")
W = tempfile.mkdtemp(prefix="phrinc.", dir=os.environ.get("TMPDIR", "/tmp"))
ok = lambda c: sum(v >= a.min_count for v in c.values()) >= a.min_chrs


def active():
    """frontier: the windows the caller would scan if every undecided window failed (telomere to the first run of
    --stop windows that do not pass); all: every window"""
    if a.mode == "all":
        return tel
    A, fail = [], 0
    for w in tel:
        A.append(w)
        fail = 0 if w in passed else fail + 1
        if fail >= a.stop:
            break
    return A


while True:
    und = [w for w in active() if w not in passed and hidx[w] < H]
    if not und:
        break
    i = min(hidx[w] for w in und)
    und = sorted(w for w in und if hidx[w] == i)  # windows that still have to see haplotype order[i]
    h = order[i]
    rounds += 1
    iv = []
    for s, e in und:
        s, e = max(0, s - a.pad), min(L, e + a.pad)
        if iv and s <= iv[-1][1]:
            iv[-1][1] = max(iv[-1][1], e)
        else:
            iv.append([s, e])
    pieces = os.path.join(W, "pieces.fa")
    with open(pieces, "w") as fh:
        for s, e in iv:
            seq = run([a.samtools, "faidx", tfa, f"{T}:{s + 1}-{e}"]).split("\n", 1)[1]
            fh.write(f">P{s}\n{seq}")
    run([a.samtools, "faidx", pieces])
    partners = [u for u in haps[h] if u != T]
    pfa = os.path.join(W, "partners.fa")
    with open(pfa, "w") as fh:
        for u in partners:
            shutil.copyfileobj(open(os.path.join(a.flank_dir, san(u) + ".fa")), fh)
    run([a.samtools, "faidx", pfa])
    pbp = sum(e - s for s, e in iv)
    work_bp += pbp * len(partners)
    # direction 1: partner flank as target, our pieces as queries (one run per partner, as the all-vs-one runs)
    def d1(u):
        o = os.path.join(W, f"d1.{san(u)}.paf")
        run([a.wfmash, "-p", "95", "-t", "1", os.path.join(a.flank_dir, san(u) + ".fa"), pieces], o)
        return o
    with ThreadPoolExecutor(a.threads) as ex:
        d1files = list(ex.map(d1, partners))
    # direction 2: our pieces as target, the partner flanks as queries
    d2file = os.path.join(W, "d2.paf")
    run([a.wfmash, "-p", "95", "-t", str(a.threads), pieces, pfa], d2file)
    rp = os.path.join(W, "round.paf")
    n = 0
    with open(rp, "w") as fh:  # lift P<offset> back to flank coordinates
        for f, side in [(x, "q") for x in d1files] + [(d2file, "t")]:
            for l in open(f):
                x = l.rstrip("\n").split("\t")
                j = 0 if side == "q" else 5
                off = int(x[j][1:])
                x[j], x[j + 1], x[j + 2], x[j + 3] = T, str(L), str(int(x[j + 2]) + off), str(int(x[j + 3]) + off)
                fh.write("\t".join(x) + "\n"); n += 1
    newly = 0
    if n:
        lst = os.path.join(W, "list"); open(lst, "w").write(rp + "\n")
        idx = os.path.join(W, "round.impg")
        if os.path.exists(idx):
            os.remove(idx)
        run([a.impg, "index", "--alignment-list", lst, "-i", idx, "-t", "1"])
        for w in und:
            res = run([a.impg, "query", "--alignment-list", lst, "-i", idx, "-r", f"{T}:{w[0]}-{w[1]}",
                       "--min-identity", a.min_identity, "--min-output-length", a.min_output_length, "--threads", "1", "-v", "0"])
            for line in res.splitlines():
                nm = line.split("\t")[0]
                if nm != T:  # the own record is already counted
                    cnt[w][chrom(nm)] += 1
    for w in und:
        hidx[w] += 1
        if ok(cnt[w]):
            passed[w] = hidx[w]; newly += 1
    log(f"round {rounds}: haplotype {i + 1}/{H} {h}, {len(und)} windows in {len(iv)} pieces ({pbp} bp), {len(partners)} partner flanks, "
        f"{n} records, {newly} windows pass, elapsed {time.time() - t0:.0f} s")
shutil.rmtree(W)
last, fail = -1, 0
for i, w in enumerate(tel):
    if w in passed:
        last, fail = i, 0
    else:
        fail += 1
        if fail >= a.stop:
            break
region = (last + 1) * a.window
full_bp = mx * sum(len(v) for v in haps.values())
with open(a.out, "w") as fh:
    fh.write("seq\twindow_start\tdist_from_telomere\tpasses\thaplotypes_at_pass\thaplotypes_checked\tchrs_ge_min\n")
    for w in tel:
        d = w[0] if T.endswith("_parm") else mx - w[1]
        fh.write(f"{T}\t{w[0]}\t{d}\t{w in passed}\t{passed.get(w, 'NA')}\t{hidx[w]}\t"
                 f"{','.join(sorted(c for c, v in cnt[w].items() if v >= a.min_count)) or '.'}\n")
    fh.write(f"#SUMMARY\t{T}\tregion_length={region}\twindows_pass={len(passed)}/{len(wins)}\trounds={rounds}\t"
             f"work_bp={work_bp}\tfull_bp={full_bp}\tseconds={time.time() - t0:.0f}\n")
log(f"done: region {region} bp, {len(passed)} of {len(wins)} windows pass, work {work_bp / full_bp:.3f} of the full scan, "
    f"{time.time() - t0:.0f} s")
