#!/usr/bin/env python3
"""Telomeric blocks that KaryoScope finds in the flanks (subtelomeric set: canonical_telomere or noncanonical_telomere,
merged within --gap bp, >= --min-len bp) and what their sequence is. Per block: bp from the flank's telomere end (near
edge), length, bp labelled canonical and noncanonical, the strand (G-rich TTAGGG or C-rich CCCTAA), and on the G strand the
share of bases in canonical TTAGGG, in telomere variant repeats (TVR: TCAGGG, TGAGGG, TTGGGG, TAAGGG, CTAGGG, TTCGGG,
GTAGGG, TTAGGGG counted as TTAGGG plus one G) and in neither; the bp of TAR1 (KaryoScope) within --tar bp of the block.
Output TSV per block plus a summary by flank verdict (--table) and by position (at the end: < --end bp from the telomere)."""
import argparse, collections, csv, re

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--bed", required=True, help="KaryoScope subtelomeric smoothed BED")
ap.add_argument("--fasta", required=True)
ap.add_argument("--table", required=True)
ap.add_argument("--gap", type=int, default=500)
ap.add_argument("--min-len", type=int, default=1000)
ap.add_argument("--tar", type=int, default=5000)
ap.add_argument("--end", type=int, default=2000)
ap.add_argument("--out", required=True)
a = ap.parse_args()
TVR = ["TCAGGG", "TGAGGG", "TTGGGG", "TAAGGG", "CTAGGG", "TTCGGG", "GTAGGG"]
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
fa, name, seq = {}, None, []
for l in open(a.fasta):
    if l.startswith(">"):
        if name:
            fa[name] = "".join(seq)
        name, seq = l[1:].split()[0], []
    else:
        seq.append(l.strip().upper())
fa[name] = "".join(seq)
verdict = {r["flank"]: r["verdict"] for r in csv.DictReader(open(a.table), delimiter="\t")}
iv, tar = collections.defaultdict(list), collections.defaultdict(list)
for l in open(a.bed):
    c, s, e, f = l.rstrip("\n").split("\t")[:4]
    if f in ("canonical_telomere", "noncanonical_telomere"):
        iv[c].append([int(s), int(e), {f: int(e) - int(s)}])
    elif f == "TAR1":
        tar[c].append((int(s), int(e)))


def count(s, m):  # overlapping occurrences
    return len(re.findall(f"(?={m})", s))


rows = []
for c, v in iv.items():
    v.sort(); bl = []
    for s, e, d in v:
        if bl and s - bl[-1][1] <= a.gap:
            bl[-1][1] = max(bl[-1][1], e)
            for k, x in d.items():
                bl[-1][2][k] = bl[-1][2].get(k, 0) + x
        else:
            bl.append([s, e, dict(d)])
    L = len(fa[c]); q = c.endswith("_qarm") or c not in verdict
    for s, e, d in bl:
        if e - s < a.min_len:
            continue
        x = fa[c][s:e]; g, cc = count(x, "TTAGGG"), count(x, "CCCTAA")
        gs = x if g >= cc else rc(x); n = len(x)
        can = 6 * count(gs, "TTAGGG") / n; tv = 6 * sum(count(gs, m) for m in TVR) / n
        d0 = L - e if q else s
        tb = sum(max(0, min(te, e + a.tar) - max(ts, s - a.tar)) for ts, te in tar[c])
        rows.append([c, verdict.get(c, "deep"), d0, n, d.get("canonical_telomere", 0), d.get("noncanonical_telomere", 0),
                     "G" if g >= cc else "C", round(can, 3), round(tv, 3), round(max(0, 1 - can - tv), 3), tb])
with open(a.out, "w") as fh:
    fh.write("flank\tverdict\tdist_from_telomere\tlength\tbp_canonical_label\tbp_noncanonical_label\tstrand\tshare_TTAGGG\tshare_TVR\tshare_other\tTAR1_bp_within\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
print(f"{len(rows)} blocks >= {a.min_len} bp in {len(set(r[0] for r in rows))} flanks")
for pos, sel in (("at the end", lambda r: r[2] < a.end), ("inside", lambda r: r[2] >= a.end)):
    for vd in ("wrong", "wrong_tiny", "undecided", "correct", "deep"):
        v = [r for r in rows if r[1] == vd and sel(r)]
        if v:
            med = lambda i: sorted(r[i] for r in v)[len(v) // 2]
            print(f"{pos}\t{vd}\tblocks {len(v)}\tflanks {len(set(r[0] for r in v))}\tmedian length {med(3)}\tmedian TTAGGG {med(7)}\t"
                  f"median TVR {med(8)}\tmedian other {med(9)}\twith TAR1 within {a.tar} bp {sum(r[10] > 0 for r in v)}\t"
                  f"labelled mostly noncanonical {sum(r[5] > r[4] for r in v)}")
