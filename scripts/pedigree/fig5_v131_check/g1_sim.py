#!/usr/bin/env python3
"""Put known errors into one haplotype.

usage: g1_sim.py WASHU SAMPLE VERSION HAP SEED OUT_DIR
Clean 10 kb bins: HiFi and ONT Q0 scans unflagged (depth > 0.25 x contig median, SNV < 3) for +- 300 kb, no cenSat.
Types x sizes (2, 10, 50 kb), one each, >= 2 Mb apart, >= 1 Mb from ends: del (remove), dup (tandem copy), foreign
(replace with a clean segment of another contig), snv (1 change per 1,000 bp); term_del removes the last L bp of a
contig 20 kb before its end (the last 20 kb, telomere, kept): read flags in the last 30 kb and the subtelomeric
satellites are allowed there, since every contig end has them, but no other flag in the last 1 Mb; none = 10 clean windows left as they are (controls).
out: OUT_DIR/norm/SAMPLE#{1,2}.fa (HAP modified, other copied), OUT_DIR/truth.tsv (original and modified coordinates).
"""
import collections, gzip, os, random, statistics as st, sys

D, SM, V, H, seed, O = sys.argv[1:7]; H = int(H); rng = random.Random(int(seed))
def read(f):
    name, seq = None, []
    for l in gzip.open(f, "rt"):
        if l[0] == ">":
            if name: yield name, "".join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    yield name, "".join(seq)
def aname(c):
    s, h, ch = c.split("#")
    return f"{s}.{ch}." + ({"1": "maternal", "2": "paternal"}[h] if s == "PAN027" else f"haplotype{h}")
cs = collections.defaultdict(list)
for l in open(f"{D}/annot_v11/{SM}.polished.diploid.cenSat.bed"):
    t = l.split("\t")
    if len(t) >= 3 and t[1].isdigit(): cs[t[0]].append((int(t[1]), int(t[2])))
seqs = dict(read(f"{D}/norm_{V}/{SM}#{H}.fa.gz"))
clean, flag = {}, {}
for c in seqs:
    bad = set()
    for sc in ("hifi", "ont"):
        b = [(int(r[1]), float(r[5]), int(r[3])) for r in (x.split("\t") for x in open(f"{D}/check/wg_scan/bins_{sc}_q0/{SM}.{V}.{c.replace('#', '_')}.bins.tsv"))]
        med = st.median(d for _, d, _ in b)
        bad |= {p for p, d, s in b if d <= 0.25 * med or (s >= 3 and d >= 3)}
    flag[c] = set(bad)
    bad |= {p for s, e in cs[aname(c)] for p in range(s // 10000 * 10000, e, 10000)}
    clean[c] = {p for p in range(0, len(seqs[c]), 10000) if not any(q in bad for q in range(p - 300000, p + 310000, 10000))}
    print(f"{c}: length {len(seqs[c])}, clean bins {len(clean[c])}", file=sys.stderr)
used = collections.defaultdict(list)
def pick(L, term=False):
    for _ in range(100000):
        c = rng.choice(sorted(seqs)); n = len(seqs[c])
        p = (n - 20000 - L) if term else rng.randrange(1000000, n - 1000000 - L)
        p0 = p // 10000 * 10000
        if term and any(q in flag[c] for q in range(max(n - 1000000, 0) // 10000 * 10000, n - 30000, 10000)): continue
        if not term and not all(q in clean[c] for q in range(p0, p + L, 10000)): continue
        if any(abs(p - u) < 2000000 for u in used[c]): continue
        used[c].append(p); return c, p
    sys.exit(f"ERROR no clean place for L={L} term={term}")
ev = []
for L in (2000, 10000, 50000):
    for t in ("del", "dup", "foreign", "snv", "term_del"):
        c, p = pick(L, t == "term_del"); ev.append((t, L, c, p))
for _ in range(10): c, p = pick(50000); ev.append(("none", 50000, c, p))
new, truth = dict(seqs), []
for c in seqs:
    s, shift = seqs[c], 0; parts, last = [], 0
    for t, L, cc, p in sorted((e for e in ev if e[2] == c), key=lambda e: e[3]):
        parts.append(s[last:p]); seg = s[p:p + L]; ms = p + shift
        if t in ("del", "term_del"): rep = ""
        elif t == "dup": rep = seg + seg
        elif t == "foreign":
            fc, fp = pick(L); rep = seqs[fc][fp:fp + L]
        elif t == "snv":
            rep = list(seg)
            for i in range(rng.randrange(1000), L, 1000): rep[i] = rng.choice([b for b in "ACGT" if b != rep[i].upper()])
            rep = "".join(rep)
        else: rep = seg
        parts.append(rep); last = p + L; shift += len(rep) - L
        truth.append((t, L, c, p, p + L, ms, ms + len(rep)))
    parts.append(s[last:]); new[c] = "".join(parts)
os.makedirs(f"{O}/norm", exist_ok=True)
with open(f"{O}/norm/{SM}#{H}.fa", "w") as o:
    for c in seqs: o.write(f">{c}\n{new[c]}\n")
with open(f"{O}/truth.tsv", "w") as o:
    o.write("type\tlength\tcontig\torig_start\torig_end\tmod_start\tmod_end\n")
    for r in sorted(truth, key=lambda r: (r[2], r[3])): o.write("\t".join(map(str, r)) + "\n")
print(f"{SM}#{H} {V}: {len(truth)} events ({collections.Counter(r[0] for r in truth)}), seed {seed}", file=sys.stderr)
