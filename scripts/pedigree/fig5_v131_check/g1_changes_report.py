#!/usr/bin/env python3
"""What we changed in one person's assembly, and where: one row per change, one summary row per haplotype.

usage: g1_changes_report.py SAMTOOLS D VERSION SAMPLE APPLIED2_DIR FINAL_DIR FIXES_TSV MQ_ORIG_TAG MQ_FINAL_TAG OUT_PREFIX
Inputs in D: norm_VERSION/SAMPLE#{1,2}.fa.gz (original release), APPLIED2_DIR (G1 round 2: corrections.tsv and the pre-polish
SAMPLE#{1,2}.fa), FINAL_DIR/SAMPLE#{1,2}.fa.gz (corrected), FIXES_TSV (g1_manual_fixes.tsv), merqury/<tag>/SAMPLE/
SAMPLE_hap{1,2}_only.bed.gz (error k-mers = assembly 21-mers absent from the Illumina reads) and SAMPLE.qv.
Coordinates are 0-based, end excluded. G1 piece on the final assembly: polish edits only inside pieces (a few bp), so the
sequence before and after the piece in the pre-polish assembly is searched in the final contig within +- 2 kb of the
expected position (pre-polish position + shift of the previous piece on the contig), with anchors of 1, 5 and 20 kb in turn;
the first anchor found once gives the coordinate, otherwise NA and a log line. end_replace: [pos, end) -> [junction, end); tip_trim: [keep, end) -> [end, end).
merqury k-mers per change: error k-mer starts in the interval +- 1 kb (same window as g1_merqury_filter.py).
out: OUT_PREFIX.changes.tsv, OUT_PREFIX.summary.tsv; log on stderr.
"""
import bisect, collections, gzip, re, subprocess, sys
S, D, V, SM, A2, FD, FX, MQO, MQF, OUT = sys.argv[1:11]; ANCH, TOL = (1000, 5000, 20000), 2000
def log(m): print(m, file=sys.stderr, flush=True)
def get(fa, c, a=None, b=None):
    if a is not None and b <= a: return ""  # empty interval (G1 insertion)
    r = c if a is None else f"{c}:{a + 1}-{b}"
    return "".join(subprocess.run([S, "faidx", fa, r], capture_output=True, text=True, check=True).stdout.split("\n")[1:])
def fai(fa): return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fa + ".fai")}
def nN(s): return s.count("N") + s.count("n")
def beds(tag):
    p = collections.defaultdict(list)
    for h in (1, 2):
        for l in gzip.open(f"{D}/merqury/{tag}/{SM}/{SM}_hap{h}_only.bed.gz", "rt"):
            f = l.split("\t"); p[f[0]].append(int(f[1]))
    for v in p.values(): v.sort()
    return p
def kmers(p, c, a, b): v = p.get(c, []); return bisect.bisect_left(v, b + 1000) - bisect.bisect_left(v, max(0, a - 1000))
def qv(tag):  # SAMPLE.qv lines: assembly, asm-only k-mers, total k-mers, QV, error rate
    return {l.split("\t")[0]: (int(l.split("\t")[1]), float(l.split("\t")[3])) for l in open(f"{D}/merqury/{tag}/{SM}/{SM}.qv")}
orig = {h: f"{D}/norm_{V}/{SM}#{h}.fa.gz" for h in (1, 2)}; pre = {h: f"{A2}/{SM}#{h}.fa" for h in (1, 2)}
fin = {h: f"{FD}/{SM}#{h}.fa.gz" for h in (1, 2)}
lo, lf = {h: fai(orig[h]) for h in (1, 2)}, {h: fai(fin[h]) for h in (1, 2)}
po, pf = beds(MQO), beds(MQF)
log(f"{SM} {V}: original {orig[1]}, pre-polish {A2}, final {FD}, merqury {MQO} -> {MQF}")
cols = ["version", "sample", "haplotype", "contig", "step", "change_id", "orig_start", "orig_end", "orig_bp", "orig_N",
        "final_start", "final_end", "final_bp", "final_N", "mode", "hifi_votes_old_new", "ont_votes_old_new",
        "merqury_kmers_orig", "merqury_kmers_final", "reason"]
rows, na = [], 0
last = [None, 0]
def lift(s, fa, c, p, exp, left):  # final position of pre-polish position p (expected exp); left: anchor ends at p
    n = []
    for L in ANCH:
        an = get(fa, c, max(0, p - L), p) if left else get(fa, c, p, p + L)
        if len(an) < min(L, 1000): break
        lo_, hi_ = exp - TOL - (len(an) if left else 0), exp + TOL + (0 if left else len(an))
        w = s[max(0, lo_):hi_]; k = [max(0, lo_) + m.start() for m in re.finditer(f"(?={re.escape(an)})", w)]
        n.append(f"{len(an)} bp: {len(k)}")
        if len(k) == 1: return (k[0] + len(an) if left else k[0]), ", ".join(n)
    return None, ", ".join(n)
def find(s, a, lo_, hi_):  # unique position of a in s[lo_:hi_]
    w = s[max(0, lo_):hi_]; n = w.count(a)
    return (max(0, lo_) + w.find(a), n) if n == 1 else (None, n)
fs = {}
def final_seq(h, c):  # one contig cached
    if (h, c) not in fs: fs.clear(); fs[(h, c)] = get(fin[h], c)
    return fs[(h, c)]
fixed = {(int(l.split("\t")[3].split("#")[1]), l.split("\t")[3]) for l in list(open(FX))[1:] if l.split("\t")[0] == V and l.split("\t")[1] == SM}
for l in open(f"{A2}/corrections.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "applied": continue
    loc, c, a, b, ns, ne = f[1], f[2], int(f[3]), int(f[4]), int(f[7]), int(f[8]); h = int(c.split("#")[1])
    s = final_seq(h, c); n1 = n2 = "contig end"
    if c != last[0]: last[:] = [c, 0]  # shift final - pre-polish, from the previous resolved piece on the contig
    fa_ = 0 if a == 0 else None  # piece starts at the contig start
    if fa_ is None: fa_, n1 = lift(s, pre[h], c, ns, ns + last[1], True)
    if b == lo[h][c]:  # piece ends at the contig end
        fe_ = len(s)
        if (h, c) in fixed: log(f"NOTE {loc} {c} ends at the contig end, which a manual fix changes later")
    else: fe_, n2 = lift(s, pre[h], c, ne, ne + (fa_ - ns if fa_ is not None else last[1]), False)
    if fa_ is None or fe_ is None:
        na += 1; log(f"NA {loc} {c}: left anchor {n1}, right anchor {n2} (occurrences within +- {TOL} bp of the expected position)"); fa_ = fe_ = None
    else: last[1] = fe_ - ne
    on = nN(get(orig[h], c, a, b))
    fN = nN(s[fa_:fe_]) if fa_ is not None else "NA"
    rows.append([V, SM, h, c, "G1", loc, a, b, b - a, on, fa_ if fa_ is not None else "NA", fe_ if fe_ is not None else "NA",
                 fe_ - fa_ if fa_ is not None else "NA", fN, f[9], f[12], f[13], kmers(po, c, a, b),
                 kmers(pf, c, fa_, fe_) if fa_ is not None else "NA", "HiFi and ONT reads vote for the reassembled piece (hifiasm, anchors >= 5 kb exact)"])
for l in list(open(FX))[1:]:
    ver, sm, step, c, cf, pos, src, why = l.rstrip("\n").split("\t")
    if ver != V or sm != SM: continue
    h = int(c.split("#")[1]); pos = int(pos); o = get(orig[h], c).upper(); s = final_seq(h, c).upper(); Lo, Lf = len(o), len(s)
    if step == "end_replace":
        an = get(f"{D}/{cf}", c, pos - 20000, pos).upper(); x, n1 = find(o, an, 0, Lo); y, n2 = find(s, an, 0, Lf)
        if x is None or y is None: sys.exit(f"ERROR {step} {c}: anchor found {n1} times in the original, {n2} in the final")
        a, b, fa_, fe_ = x + 20000, Lo, y + 20000, Lf
    elif step == "tip_trim":
        an = get(f"{D}/{cf}", c, pos - 2000, pos).upper(); x, n1 = find(o, an, Lo - 100000, Lo)
        if x is None or s[-2000:] != an: sys.exit(f"ERROR {step} {c}: anchor found {n1} times in the original end; final ends with it: {s[-2000:] == an}")
        a, b, fa_, fe_ = x + 2000, Lo, Lf, Lf
    else: sys.exit(f"ERROR unknown step {step}")
    rows.append([V, SM, h, c, step, f"fix:{step}", a, b, b - a, nN(o[a:b]), fa_, fe_, fe_ - fa_, nN(s[fa_:fe_]), "NA", "NA", "NA",
                 kmers(po, c, a, b), kmers(pf, c, fa_, fe_), why])
rows.sort(key=lambda r: (r[2], [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", r[3])], r[6]))
with open(f"{OUT}.changes.tsv", "w") as o:
    o.write("\t".join(cols) + "\n"); o.writelines("\t".join(map(str, r)) + "\n" for r in rows)
qo, qf = qv(MQO), qv(MQF)
with open(f"{OUT}.summary.tsv", "w") as o:
    o.write("version\tsample\thaplotype\tcontigs\torig_bp\tfinal_bp\torig_N\tfinal_N\torig_N_pct\tfinal_N_pct\torig_N_runs\tfinal_N_runs\t"
            "final_non_ACGTN\tmerqury_kmers_orig\tmerqury_kmers_final\tQV_orig\tQV_final\tchanges_G1\tchanges_manual\tNA_coordinates\n")
    for h in (1, 2):
        st = {}
        for k, fa in (("o", orig[h]), ("f", fin[h])):
            n = runs = bad = L = 0
            for l in gzip.open(fa, "rt"):
                if l[0] == ">": prev = ""; continue
                s = l.rstrip("\n"); L += len(s); n += nN(s); runs += len(re.findall(r"[Nn]+", s)) - (prev in ("N", "n") and s[:1] in ("N", "n"))
                bad += len(s) - sum(s.upper().count(x) for x in "ACGTN"); prev = s[-1:]
            st[k] = (L, n, runs, bad)
        g = sum(r[2] == h and r[4] == "G1" for r in rows); m = sum(r[2] == h and r[4] != "G1" for r in rows)
        key = f"{SM}_hap{h}"
        o.write("\t".join(map(str, [V, SM, h, len(lf[h]), st["o"][0], st["f"][0], st["o"][1], st["f"][1],
                f"{100 * st['o'][1] / st['o'][0]:.4f}", f"{100 * st['f'][1] / st['f'][0]:.4f}", st["o"][2], st["f"][2], st["f"][3],
                qo[key][0], qf[key][0], qo[key][1], qf[key][1], g, m, sum(r[2] == h and r[10] == "NA" for r in rows)])) + "\n")
        log(f"hap{h}: N {st['o'][1]} -> {st['f'][1]}, N runs {st['o'][2]} -> {st['f'][2]}, non-ACGTN final {st['f'][3]}, G1 {g}, manual {m}")
log(f"wrote {OUT}.changes.tsv ({len(rows)} changes, {na} with NA final coordinates) and {OUT}.summary.tsv")
