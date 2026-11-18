#!/usr/bin/env python3
"""Apply accepted G1 pieces to one person's assembly.

usage: g1_apply_all.py NORM_DIR SAMPLE RUN_DIR ANCHOR VOTE OUT_DIR
Per locus RUN_DIR/<id>.<name>/ (g1_locus.sbatch): candidates of rule ANCHOR (any mode); a candidate passes when the
VOTE rule (more, ratio2, binom as in g1_grid_score.py) holds for HiFi and ONT; the passing one with most HiFi + ONT
candidate votes is applied. Overlapping accepted pieces on the same contig: the one with more votes is kept, the other
is listed as dropped. Pieces are applied right to left, so coordinates in corrections.tsv are on the input assembly.
out: OUT_DIR/SAMPLE#{1,2}.fa (all contigs), OUT_DIR/corrections.tsv (applied and dropped, with votes and md5).
"""
import collections, glob, gzip, hashlib, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g1_grid_score import passes, votes  # same vote counting as the grid

na, sm, rd, ar, vr, out = sys.argv[1:7]
acc = []
for sf in sorted(glob.glob(f"{rd}/*/summary.tsv")):
    ev = os.path.dirname(sf); rows = [l.rstrip("\n").split("\t") for l in open(sf)]
    rows = [dict(zip(rows[0], r)) for r in rows[1:]]; v = votes(ev); best = None
    for r in rows:
        if r["rule"] != ar or r["id"] == "NA": continue
        h, o = v["hifi"][r["id"]], v["ont"][r["id"]]
        if passes(vr, *h) and passes(vr, *o) and (best is None or h[1] + o[1] > best[0]): best = (h[1] + o[1], r, h, o)
    if best is None: continue
    n, r, h, o = best
    seq = None
    for l in open(f"{ev}/pieces.fa"):
        if l.startswith(">"): cur = l[1:].split()[0]
        elif cur == r["id"]: seq = l.strip()
    if seq is None: sys.exit(f"ERROR {r['id']} not in {ev}/pieces.fa")
    c = r["error"].split(":")[0]
    acc.append(dict(locus=os.path.basename(ev), contig=c, start=int(r["asm_start"]), end=int(r["asm_end"]), seq=seq,
                    mode=r["mode"], votes=n, hifi=f"{h[0]}/{h[1]}", ont=f"{o[0]}/{o[1]}"))
print(f"{sm}: {len(glob.glob(rd + '/*/summary.tsv'))} loci, {len(acc)} accepted ({ar}, {vr})", file=sys.stderr)
by = collections.defaultdict(list)
for a in acc: by[a["contig"]].append(a)
keep, drop = [], []
for c, l in by.items():
    for a in sorted(l, key=lambda a: -a["votes"]):
        (drop if any(a["start"] < b["end"] and b["start"] < a["end"] for b in keep if b["contig"] == c) else keep).append(a)
print(f"{sm}: kept {len(keep)}, dropped {len(drop)} (overlap with a piece with more votes)", file=sys.stderr)
os.makedirs(out, exist_ok=True)
for h in (1, 2):
    with gzip.open(f"{na}/{sm}#{h}.fa.gz", "rt") as f, open(f"{out}/{sm}#{h}.fa", "w") as o:
        name, parts = None, []
        def flush():
            if name is None: return
            s = "".join(parts)
            for a in sorted((a for a in keep if a["contig"] == name), key=lambda a: -a["start"]):
                s = s[:a["start"]] + a["seq"] + s[a["end"]:]
            o.write(f">{name}\n{s}\n")
        for l in f:
            if l.startswith(">"): flush(); name, parts = l[1:].split()[0], []
            else: parts.append(l.strip())
        flush()
with open(f"{out}/corrections.tsv", "w") as o:
    o.write("status\tlocus\tcontig\tstart\tend\told_bp\tnew_bp\tmode\tanchor\tvote\thifi_votes_asm_cand\tont_votes_asm_cand\tpiece_md5\n")
    for st, l in (("applied", keep), ("dropped", drop)):
        for a in sorted(l, key=lambda a: (a["contig"], a["start"])):
            o.write(f"{st}\t{a['locus']}\t{a['contig']}\t{a['start']}\t{a['end']}\t{a['end'] - a['start']}\t{len(a['seq'])}\t{a['mode']}\t{ar}\t{vr}\t{a['hifi']}\t{a['ont']}\t{hashlib.md5(a['seq'].encode()).hexdigest()[:12]}\n")
print(f"{sm}: wrote {out}", file=sys.stderr)
