#!/usr/bin/env python3
"""Alignment check of site-level conversion events (conversion_sites.py): are the event's donor-type sites still there when
the copies are aligned again with a different aligner (minimap2) instead of wfmash?

usage:
  sites_realign_check.py sample EVENTS_TSV SIGNIFICANT_TSV N_SIG_ONE N_SIG_SHARED N_NONSIG SEED OUT_TSV
      random events: significant ones found in one individual, significant shared ones, non-significant ones (all strata
      limited to spans <= 20 kb); OUT_TSV adds the column 'stratum'.
  sites_realign_check.py regions SAMPLE_TSV PAF_DIR PAD OUT_PREFIX
      per event, the region R = [start - PAD, end + PAD) of the target window, and for every copy of the window's own arm and
      of the primary donor arm the query interval that its wfmash record aligns to R (record chosen as conversion_sites.py
      does: most aligned target bases first). out: OUT_PREFIX.regions.txt (samtools faidx regions), OUT_PREFIX.copies.tsv
  sites_realign_check.py check SAMPLE_TSV PAF_DIR COPIES_TSV SEQ_FA MINIMAP2 PAD OUT_PREFIX
      wfmash classes recomputed exactly as conversion_sites.py (whole window); minimap2 -cx asm20 --eqx of each copy
      sequence against the target region, best-scoring primary record per copy; sites with the same thresholds
      (>= 20 copies per arm, hi 0.9, lo 0.1). An event is REPRODUCED when >= 2 of its wfmash donor-type sites are donor-type
      sites again and no own-arm site falls between its first and last site. out: OUT_PREFIX.events.tsv, log on stderr.
PAF_DIR holds hprcv2.1Mb.telo_500kb_trimmed.all-vs-<target with #, : and - as _>.paf.gz.
"""
import collections, csv, gzip, os, random, re, subprocess, sys, tempfile, time
import numpy as np
CIG = re.compile(r"(\d+)([=XIDM])"); ARM = re.compile(r"_(chr[0-9XYM]+_(?:parm|qarm))$")
MINC, HI, LO = 20, 0.9, 0.1
def log(m): print(f"[{time.strftime('%F %T')}] {m}", file=sys.stderr, flush=True)
def arm(n): return ARM.search(n).group(1)
def paf_path(d, t): return f"{d}/hprcv2.1Mb.telo_500kb_trimmed.all-vs-{re.sub(r'[#:-]', '_', t)}.paf.gz"

def read_paf(path, arms):
    recs = collections.defaultdict(list); target = tlen = None
    for l in gzip.open(path, "rt"):
        f = l.rstrip("\n").split("\t"); q, t = f[0], f[5]
        if target is None: target, tlen = t, int(f[6])
        if q == t or arm(q) not in arms: continue
        cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        recs[q].append((int(f[8]) - int(f[7]), int(f[7]), cg, int(f[1]), int(f[2]), int(f[3]), f[4]))
    return target, tlen, recs

def classes(cov, mat, own, y):  # positions as arrays: A-sites and Y-sites for donor y
    nA, mA, nY, mY = (x.astype(np.int64) for x in (cov[own], mat[own], cov[y], mat[y]))
    ok = (nA >= MINC) & (nY >= MINC)
    return np.flatnonzero(ok & (mA >= HI * nA) & (mY <= LO * nY)), np.flatnonzero(ok & (mA <= LO * nA) & (mY >= HI * nY))

def wfmash_arrays(tlen, recs, keep):  # exactly conversion_sites.scan: per query, records by aligned length, first come first served
    cov, mat, src = {}, {}, {}
    for q, rs in recs.items():
        a = arm(q)
        if a not in cov: cov[a], mat[a] = np.zeros(tlen, np.uint16), np.zeros(tlen, np.uint16)
        free = np.ones(tlen, bool)
        for r in sorted(rs, key=lambda x: (x[0], x[1], x[2]), reverse=True):
            pos = r[1]; ca, cb = keep
            for n, op in CIG.findall(r[2]):
                n = int(n)
                if op in "=X":
                    f = free[pos:pos + n]; cov[a][pos:pos + n] += f
                    if op == "=": mat[a][pos:pos + n] += f
                    if pos < cb and pos + n > ca and f[max(ca - pos, 0):min(cb - pos, n)].any(): src.setdefault(q, r)
                    free[pos:pos + n] = False; pos += n
                elif op == "D": pos += n
    return cov, mat, src

def query_interval(r, a, b):  # query forward-strand interval aligned to target [a, b) by record r
    tl, ts, cg, ql, qs, qe, strand = r; t, q = ts, 0; lo = hi = None
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in "=X":
            s, e = max(t, a), min(t + n, b)
            if s < e:
                q0, q1 = q + (s - t), q + (e - t); lo = q0 if lo is None else lo; hi = q1
            t += n; q += n
        elif op == "D": t += n
        elif op == "I":
            if a <= t < b and lo is not None: hi = q + n
            q += n
    if lo is None: return None
    return (qs + lo, qs + hi) if strand == "+" else (qe - hi, qe - lo)

mode = sys.argv[1]
if mode == "sample":
    ev, sg, n1, n2, n0, seed, out = sys.argv[2], sys.argv[3], *map(int, sys.argv[4:7]), int(sys.argv[7]), sys.argv[8]
    rows = list(csv.DictReader(open(ev), delimiter="\t")); sig = {(r["target"], r["start"], r["end"]) for r in csv.DictReader(open(sg), delimiter="\t")}
    st = collections.defaultdict(list)
    for r in rows:
        if int(r["span"]) > 20000: continue
        s = (r["target"], r["start"], r["end"]) in sig
        st["sig_one" if s and r["n_individuals_locus"] == "1" else "sig_shared" if s else "nonsig"].append(r)
    random.seed(seed); pick = []
    for k, n in (("sig_one", n1), ("sig_shared", n2), ("nonsig", n0)):
        log(f"stratum {k}: {len(st[k])} events with span <= 20 kb, sampling {n}")
        for r in random.sample(st[k], n): r["stratum"] = k; pick.append(r)
    with open(out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(pick[0]), delimiter="\t"); w.writeheader(); w.writerows(pick)
    log(f"wrote {out}: {len(pick)} events, seed {seed}")
elif mode == "regions":
    smp, pd, pad, out = sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    rows = list(csv.DictReader(open(smp), delimiter="\t")); reg = []; cp = []
    for i, r in enumerate(rows):
        t = r["target"]; own, y = r["arm"], r["primary_donor"]
        target, tlen, recs = read_paf(paf_path(pd, t), {own, y})
        if target != t: sys.exit(f"ERROR {t}: PAF target {target}")
        a, b = max(int(r["start"]) - pad, 0), min(int(r["end"]) + pad, tlen)
        _, _, src = wfmash_arrays(tlen, recs, (a, b))
        reg.append(f"{t}:{a + 1}-{b}"); n = collections.Counter()
        for q, rec in src.items():
            iv = query_interval(rec, a, b)
            if iv is None or iv[1] - iv[0] < 50: n["too_short"] += 1; continue
            reg.append(f"{q}:{iv[0] + 1}-{iv[1]}"); cp.append((i, t, q, arm(q), rec[6], f"{q}:{iv[0] + 1}-{iv[1]}")); n[arm(q)] += 1
        log(f"{i + 1}/{len(rows)} {t} {r['start']}-{r['end']} ({r['stratum']}): region {a}-{b}, copies {dict(n)}")
    open(f"{out}.regions.txt", "w").write("\n".join(dict.fromkeys(reg)) + "\n")
    with open(f"{out}.copies.tsv", "w") as fo:
        fo.write("event\ttarget\tquery\tarm\tstrand\tregion\n"); fo.writelines("\t".join(map(str, c)) + "\n" for c in cp)
    log(f"wrote {out}.regions.txt ({len(dict.fromkeys(reg))} regions) and {out}.copies.tsv ({len(cp)} copies)")
elif mode == "check":
    smp, pd, cpt, fa, mm2, pad, out = sys.argv[2:9]; pad = int(pad)
    seqs = {}; name = None
    for l in open(fa):
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(l.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    cps = collections.defaultdict(list)
    for r in csv.DictReader(open(cpt), delimiter="\t"): cps[int(r["event"])].append(r)
    rows = list(csv.DictReader(open(smp), delimiter="\t")); res = []
    tmp = tempfile.mkdtemp(prefix="realign_", dir=os.path.dirname(os.path.abspath(out)))
    for i, r in enumerate(rows):
        t = r["target"]; own, y = r["arm"], r["primary_donor"]; s0, e0 = int(r["start"]), int(r["end"])
        target, tlen, recs = read_paf(paf_path(pd, t), {own, y})
        a, b = max(s0 - pad, 0), min(e0 + pad, tlen)
        cov, mat, _ = wfmash_arrays(tlen, recs, (a, b))
        A0, Y0 = classes(cov, mat, own, y)
        y0 = Y0[(Y0 >= s0) & (Y0 < e0)]; a0 = A0[(A0 >= s0) & (A0 < e0)]
        tr = f"{t}:{a + 1}-{b}"; open(f"{tmp}/t.fa", "w").write(f">{tr}\n{seqs[tr]}\n")
        with open(f"{tmp}/q.fa", "w") as fo:
            for c in cps[i]: fo.write(f">{c['query']}\n{seqs[c['region']]}\n")
        p = subprocess.run([mm2, "-cx", "asm20", "--eqx", "--secondary=no", "-t", "4", f"{tmp}/t.fa", f"{tmp}/q.fa"], capture_output=True, text=True, check=True)
        best = {}
        for l in p.stdout.splitlines():
            f = l.split("\t"); sc = int(next(x[5:] for x in f[12:] if x.startswith("AS:i:")))
            if f[0] not in best or sc > best[f[0]][0]: best[f[0]] = (sc, int(f[7]), next(x[5:] for x in f[12:] if x.startswith("cg:Z:")))
        L = b - a; cov2 = {own: np.zeros(L, np.uint16), y: np.zeros(L, np.uint16)}; mat2 = {own: np.zeros(L, np.uint16), y: np.zeros(L, np.uint16)}
        for c in cps[i]:
            if c["query"] not in best: continue
            _, pos, cg = best[c["query"]]
            for n, op in CIG.findall(cg):
                n = int(n)
                if op in "=X":
                    cov2[c["arm"]][pos:pos + n] += 1
                    if op == "=": mat2[c["arm"]][pos:pos + n] += 1
                    pos += n
                elif op == "D": pos += n
        A1, Y1 = classes(cov2, mat2, own, y); A1 += a; Y1 += a
        f0, l0 = (int(y0.min()), int(y0.max())) if len(y0) else (s0, e0 - 1)
        kept = int(np.isin(y0, Y1).sum()); a_in = int(((A1 > f0) & (A1 < l0)).sum())
        rep = int(kept >= 2 and a_in == 0)
        na = sum(1 for c in cps[i] if c["arm"] == own); ny = sum(1 for c in cps[i] if c["arm"] == y)
        aa = sum(1 for c in cps[i] if c["arm"] == own and c["query"] in best); ay = sum(1 for c in cps[i] if c["arm"] == y and c["query"] in best)
        res.append(dict(event=i, stratum=r["stratum"], target=t, start=s0, end=e0, community=r["community"], primary_donor=y,
                        n_sites_event=r["n_sites"], wfmash_y_sites=len(y0), wfmash_a_sites_inside=len(a0), copies_own=na, copies_donor=ny,
                        minimap2_aligned_own=aa, minimap2_aligned_donor=ay, minimap2_y_sites_kept=kept,
                        minimap2_y_sites_in_span=int(((Y1 >= s0) & (Y1 < e0)).sum()), minimap2_a_sites_between=a_in, reproduced=rep))
        log(f"{i + 1}/{len(rows)} {r['stratum']} {t} {s0}-{e0}: wfmash Y-sites {len(y0)} (event n_sites {r['n_sites']}), A inside {len(a0)}; "
            f"copies own {aa}/{na}, donor {ay}/{ny} aligned by minimap2; Y-sites kept {kept}, A-sites between {a_in} -> reproduced {rep}")
    with open(f"{out}.events.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(res[0]), delimiter="\t"); w.writeheader(); w.writerows(res)
    for k in ("sig_one", "sig_shared", "nonsig"):
        s = [x for x in res if x["stratum"] == k]
        if s: log(f"{k}: {sum(x['reproduced'] for x in s)} of {len(s)} reproduced; wfmash Y-sites {sum(x['wfmash_y_sites'] for x in s)}, kept {sum(x['minimap2_y_sites_kept'] for x in s)}")
    log(f"wrote {out}.events.tsv")
