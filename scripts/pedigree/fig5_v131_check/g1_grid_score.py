#!/usr/bin/env python3
"""Score G1 settings offline.

usage: g1_grid_score.py OUT_TSV SET... ; SET = label:kind:grid_dir (window runs grid_dir/f<bp>/)
Margins: a read votes only when its best AS beats the next by >= margin x best AS (MARGINS; 0 = original).
kind = sim (run_dir events <i>.<type>_<len>, check_tsv from g1_sim_check.sbatch) or holdout (DAD: chr9, chr7, chr3,
chr10 = real errors, chr21ctl = clean). Sim needs grid_dir/f<bp>/check.tsv (g1_sim_check.sbatch). Per window run and setting (anchor rule x vote rule x NM tolerance), per event:
candidates of that anchor rule (any assembly mode); a candidate passes when the vote rule holds for HiFi and for ONT
(votes from votes.<set>.tsv.gz: a read votes for a side when its AS beats the two others). Among passing candidates the
one with the most HiFi + ONT candidate votes is applied.
Vote rules: more (cand > asm), ratio2 (cand >= 2 x asm and cand >= 3), binom (one-sided binomial P(X >= cand) < 0.01,
n = cand + asm, p = 0.5); ontlead: ONT more, and HiFi not against (HiFi cand >= asm). Sim: right = applied on an error and the patched window equals the original within the
tolerance (NM <= tol, |length difference| <= tol, matches >= original length - tol); wrong = applied on a control, or
applied on an error without restoring it; score = right - wrong. Holdout: fixed real errors and changed controls.
"""
import collections, glob, gzip, math, os, sys

RULES = ["ex:5000", "ex:10000", "ex:20000", "id:0.99", "id:0.995", "id:0.999"]
VOTES = ["more", "ratio2", "binom", "ontlead"]; TOLS = [10, 50, 200]; MARGINS = [0.0, 0.001, 0.002, 0.005, 0.01]
def binom_p(c, n):
    lc = lambda k: math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) - n * math.log(2)  # python 3.6: no math.comb
    return sum(math.exp(lc(k)) for k in range(c, n + 1)) if n else 1.0
def passes(rule, a, c):
    if rule == "more": return c > a
    if rule == "ratio2": return c >= 2 * a and c >= 3
    return binom_p(c, a + c) < 0.01
def votes(ev, margin=0.0):
    v = {}
    for s in ("hifi", "ont"):
        f = f"{ev}/votes.{s}.tsv.gz"
        cnt = collections.defaultdict(lambda: [0, 0])
        if os.path.exists(f):
            for l in gzip.open(f, "rt"):
                cid, rn, aa, ac, ae = l.split("\t")[:5]; aa, ac, ae = int(aa), int(ac), int(ae)
                if aa > max(ac, ae) and aa - max(ac, ae) >= margin * aa: cnt[cid][0] += 1
                if ac > max(aa, ae) and ac - max(aa, ae) >= margin * ac: cnt[cid][1] += 1
        v[s] = cnt
    return v
if __name__ == "__main__":  # g1_apply_all.py imports passes and votes
    out = open(sys.argv[1], "w"); out.write("set\twindow\tanchor\tvote\tmargin\ttol\tevents\tright\twrong\tmissed\tscore\tdetail\n")
    for spec in sys.argv[2:]:
        p = spec.split(":"); label, kind, rd = p[0], p[1], p[2]
        for wd in sorted(glob.glob(f"{rd}/f*")):
            chk = {}
            if kind == "sim":
                for l in list(open(f"{wd}/check.tsv"))[1:]:
                    t = l.rstrip("\n").split("\t"); chk[(t[0], t[1])] = (int(t[6]), int(t[7]), t[8], int(t[9]))
            evs = sorted(glob.glob(f"{wd}/*/summary.tsv"))
            data = {}
            for sf in evs:
                ev = os.path.dirname(sf); rows = [l.rstrip("\n").split("\t") for l in open(sf)]
                h = rows[0]; rows = [dict(zip(h, r)) for r in rows[1:]]
                data[os.path.basename(ev)] = (rows, {m: votes(ev, m) for m in MARGINS})
            for ar, vr, mg in ((a, b, c) for a in RULES for b in VOTES for c in MARGINS):
                    for tol in (TOLS if kind == "sim" else [0]):
                        right = wrong = miss = 0; det = []
                        for name, (rows, vm) in data.items():
                            best = None; v = vm[mg]
                            for r in rows:
                                if r["rule"] != ar or r["id"] == "NA": continue
                                h, o = v["hifi"][r["id"]], v["ont"][r["id"]]
                                ok = passes("more", *o) and h[1] >= h[0] if vr == "ontlead" else passes(vr, *h) and passes(vr, *o)
                                if ok and (best is None or h[1] + o[1] > best[1]): best = (r["id"], h[1] + o[1])
                            if kind == "sim":
                                err = "none_" not in name
                                if best is None: miss += err; continue
                                ol, mt, nm, ld = chk[(name, best[0])]
                                ok = err and nm != "NA" and int(nm) <= tol and abs(ld) <= tol and mt >= ol - tol
                                right += ok; wrong += not ok; det.append(f"{name}={'R' if ok else 'W'}")
                            else:
                                err = name.split(".")[1] != "chr21ctl"
                                if best is None: miss += err; continue
                                right += err; wrong += not err; det.append(f"{name}={'fix' if err else 'W'}")
                        out.write(f"{label}\t{os.path.basename(wd)}\t{ar}\t{vr}\t{mg}\t{tol}\t{len(data)}\t{right}\t{wrong}\t{miss}\t{right - wrong}\t{','.join(det)}\n")
    out.close()
    print(open(sys.argv[1]).read())
