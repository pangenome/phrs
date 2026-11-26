#!/usr/bin/env python3
"""Replace a chromosome end with a piece of an own contig (manual fix, g1_manual_fixes.tsv step end_replace).

usage: g1_end_replace.py SAMTOOLS BGZIP IN_DIR OUT_DIR CONTIG COORD_FA POS PIECE_FA PIECE_NAME PIECE_FROM
POS is on COORD_FA (the assembly where the fix was found): CONTIG keeps its first POS bases (0-based junction).
The 20 kb of COORD_FA before POS must occur once in CONTIG of IN_DIR; there the end is cut and replaced by
PIECE_NAME[PIECE_FROM:] of PIECE_FA (same result as g1_dad_chr9_final.sbatch when IN_DIR does not change those 20 kb).
Other contigs and the other haplotype are copied unchanged, contig order kept. OUT_DIR must not exist. Log on stderr.
Also writes OUT_DIR/junction.tsv: contig, 0-based junction in OUT_DIR, new contig length.
"""
import os, shutil, subprocess, sys
S, BG, ind, outd, C, cf, pos, pf, pn, pfrom = sys.argv[1:11]; pos, pfrom = int(pos), int(pfrom)
def log(m): print(m, file=sys.stderr, flush=True)
def get(fa, r): return "".join(subprocess.run([S, "faidx", fa, r], capture_output=True, text=True, check=True).stdout.split("\n")[1:])
sm, h = C.split("#")[0], C.split("#")[1]; fn = f"{sm}#{h}.fa.gz"
r, g, p = get(cf, C), get(f"{ind}/{fn}", C), get(pf, pn)
a = r[pos - 20000:pos]; n = g.count(a)
if n != 1: sys.exit(f"ERROR the 20 kb before {pos} on {cf} occurs {n} times in {C} of {ind}")
j = g.find(a) + 20000; new = g[:j] + p[pfrom:]
log(f"{C}: junction {pos} on {cf}, {j} in {ind}; {ind} length {len(g)} -> {len(new)} ({len(g) - j} bp end replaced by {len(p) - pfrom} bp of {pn}[{pfrom}:])")
bad = len(new) - sum(new.upper().count(x) for x in "ACGTN")
if bad: sys.exit(f"ERROR {bad} non-ACGTN letters in the new {C}")
os.makedirs(outd, exist_ok=False)
names = [l.split("\t")[0] for l in open(f"{ind}/{fn}.fai")]
with subprocess.Popen([BG, "-c"], stdin=subprocess.PIPE, stdout=open(f"{outd}/{fn}", "wb"), text=True) as z:
    for c in names:
        s = new if c == C else get(f"{ind}/{fn}", c)
        z.stdin.write(f">{c}\n" + "\n".join(s[k:k + 80] for k in range(0, len(s), 80)) + "\n")
subprocess.run([S, "faidx", f"{outd}/{fn}"], check=True)
open(f"{outd}/junction.tsv", "w").write(f"{C}\t{j}\t{len(new)}\n")  # 0-based junction, new length: read by g1_manual_fixes.sh
for f in os.listdir(ind):
    if f.startswith(f"{sm}#") and not f.startswith(fn): shutil.copy(f"{ind}/{f}", outd)
log(f"wrote {outd}: {fn} ({len(names)} contigs) + {sorted(f for f in os.listdir(outd) if not f.startswith(fn))}")
