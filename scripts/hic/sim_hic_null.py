#!/usr/bin/env python3
"""Read-level Hi-C simulation for referee R1-M4 / R3-M3: how often do reads from one chromosome end land on another
end, and how much end-to-end contact does that alone create?

Null model: only cis contacts near each end, no true contact between ends. For each end e (PHR window from
hg002_subtelomeric_regions.bed plus its 100 kb flank window), read 1 starts uniformly in the end domain (PHR + flank);
read 2 lies at distance s from read 1, with P(s) proportional to 1/s on [MIN_S, MAX_S], towards the telomere or the
centromere at random, and must stay inside the extended domain (end domain + MAX_S inward) and the chromosome. Strands
are random; reads have uniform substitution errors at rate ERR and quality F, or, with --profile (read_profile.py hic),
each read takes the profile of a random real read of the same mate: its length, base qualities and mismatch positions,
random bases over its 5' soft clip, and over its 3' soft clip the sequence of the partner locus (the other read's
position, opposite strand), as when a read runs through a ligation junction. Read names carry the truth:
s<n>:<chr1>:<pos1>:<strand1>:<chr2>:<pos2>:<strand2> (1-based leftmost reference position of each read).

subcommands
  regions  PHR_BED FLANK_BED CHROM_SIZES OUT_PREFIX [--max-s 500000]
           writes OUT.regions.txt (samtools regions of the extended domains) and OUT.domains.tsv
  simulate DOMAINS_TSV EXTENDED_FA OUT_PREFIX [--pairs-per-kb 8] [--read-len 250] [--err 0.002] [--min-s 1000]
           [--max-s 500000] [--seed 20261126] [--profile PROFILES_TSV_GZ]   writes OUT_R1.fastq.gz, OUT_R2.fastq.gz
  evaluate DOMAINS_TSV ALL_VALID_PAIRS OUT_PREFIX [--tol 600]  (HiC-Pro reports a 5' end, up to read length plus
           soft clips away from the leftmost position)
           writes OUT.reads.tsv (per read: true end and window, assigned end and window, correct within TOL bp) and
           OUT.misassign.tsv (per ordered pair of ends: reads from the PHR window of end A assigned to the PHR window
           of end B, and the rate over reads from A that HiC-Pro kept)"""
import argparse, collections, gzip, math, sys
import numpy as np

log = lambda *m: print(*m, file=sys.stderr, flush=True)
ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
r = sp.add_parser("regions"); r.add_argument("phr"); r.add_argument("flank"); r.add_argument("sizes"); r.add_argument("out")
r.add_argument("--max-s", type=int, default=500000)
s = sp.add_parser("simulate"); s.add_argument("domains"); s.add_argument("fa"); s.add_argument("out")
s.add_argument("--pairs-per-kb", type=float, default=8); s.add_argument("--read-len", type=int, default=250)
s.add_argument("--err", type=float, default=0.002); s.add_argument("--min-s", type=int, default=1000)
s.add_argument("--max-s", type=int, default=500000); s.add_argument("--seed", type=int, default=20261126)
s.add_argument("--profile", help="real-read profiles (read_profile.py hic); replaces --read-len, --err and the F qualities")
e = sp.add_parser("evaluate"); e.add_argument("domains"); e.add_argument("avp"); e.add_argument("out")
e.add_argument("--tol", type=int, default=600)
a = ap.parse_args()
COLS = ["end", "chrom", "arm", "phr_start", "phr_end", "flank_start", "flank_end", "ext_start", "ext_end"]


def read_domains(p):
    d = [dict(zip(COLS, l.rstrip("\n").split("\t"))) for l in open(p) if not l.startswith("end\t")]
    for x in d:
        for k in COLS[3:]: x[k] = int(x[k])
    return d


if a.cmd == "regions":
    size = {l.split()[0]: int(l.split()[1]) for l in open(a.sizes)}
    fl = {l.split("\t")[3].strip(): l.split("\t")[:3] for l in open(a.flank)}
    rows = []
    for l in open(a.phr):
        c, ps, pe, name = l.rstrip("\n").split("\t")[:4]; ps, pe = int(ps), int(pe)
        fc, fs, fe = fl[name]; fs, fe = int(fs), int(fe)
        if fc != c: sys.exit(f"FATAL flank of {name} on {fc}")
        arm = name.rsplit("_", 1)[1]
        es, ee = (min(ps, fs), min(size[c], max(pe, fe) + a.max_s)) if arm == "p" else (max(0, min(ps, fs) - a.max_s), max(pe, fe))
        rows.append([name, c, arm, ps, pe, fs, fe, es, ee])
    with open(a.out + ".domains.tsv", "w") as o:
        o.write("\t".join(COLS) + "\n")
        for x in rows: o.write("\t".join(map(str, x)) + "\n")
    with open(a.out + ".regions.txt", "w") as o:
        for x in rows: o.write(f"{x[1]}:{x[7] + 1}-{x[8]}\n")
    log(f"{len(rows)} ends; extended domains {sum(x[8] - x[7] for x in rows)} bp; PHR {sum(x[4] - x[3] for x in rows)} bp;"
        f" flank {sum(x[6] - x[5] for x in rows)} bp; max s {a.max_s}")

elif a.cmd == "simulate":
    D = read_domains(a.domains); rng = np.random.default_rng(a.seed)
    seqs, name, buf = {}, None, []
    with (gzip.open(a.fa, "rt") if a.fa.endswith(".gz") else open(a.fa)) as f:
        for l in f:
            if l.startswith(">"):
                if name: seqs[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    seqs[name] = "".join(buf).upper()
    comp = str.maketrans("ACGTN", "TGCAN"); B = np.frombuffer(b"ACGT", dtype=np.uint8); L = a.read_len; P = None
    if a.profile:
        P = {1: [], 2: []}
        with gzip.open(a.profile, "rt") as f:
            if not f.readline().startswith("mate\t"): sys.exit("FATAL profile header")
            for l in f:
                m, n_, q, c5, c3, mm = l.rstrip("\n").split("\t")[:6]
                P[int(m)].append((int(n_), q, int(c5), int(c3), np.array([int(v) for v in mm.split(",")] if mm != "." else [], dtype=int)))
        L = max(x[0] for v in P.values() for x in v)
        log(f"profiles {a.profile}: mate 1 {len(P[1])}, mate 2 {len(P[2])}; longest read {L}; mean mismatches per aligned base"
            f" {sum(len(x[4]) for v in P.values() for x in v) / sum(x[0] - x[2] - x[3] for v in P.values() for x in v):.5f}")
    log(f"simulate: {len(D)} ends, {len(seqs)} sequences, {a.pairs_per_kb} pairs per kb of end domain, read length {L},"
        f" error {'from profiles' if P else a.err}, s in [{a.min_s}, {a.max_s}] with P(s) ~ 1/s, seed {a.seed}")
    n, drop = 0, collections.Counter()
    with gzip.open(a.out + "_R1.fastq.gz", "wt", compresslevel=1) as o1, gzip.open(a.out + "_R2.fastq.gz", "wt", compresslevel=1) as o2:
        for d in D:
            sq = seqs[f"{d['chrom']}:{d['ext_start'] + 1}-{d['ext_end']}"]; off = d["ext_start"]
            if len(sq) != d["ext_end"] - d["ext_start"]: sys.exit(f"FATAL length of {d['end']}")
            lo, hi = min(d["phr_start"], d["flank_start"]), max(d["phr_end"], d["flank_end"])
            want = int(round((hi - lo) / 1000 * a.pairs_per_kb)); made = 0
            while made < want:
                x = int(rng.integers(lo, hi - L)); sd = int(round(math.exp(rng.uniform(math.log(a.min_s), math.log(a.max_s)))))
                y = x + sd if rng.random() < 0.5 else x - sd
                if y < d["ext_start"] or y + L > d["ext_end"]: drop["partner outside extended domain"] += 1; continue
                pr = [P[k][int(rng.integers(len(P[k])))] for k in (1, 2)] if P else [(L, "F" * L, 0, 0, None)] * 2
                rr = []
                for p, q in zip((x, y), pr):
                    t = sq[p - off:p - off + q[0]]; st = "+" if rng.random() < 0.5 else "-"
                    if st == "-": t = t.translate(comp)[::-1]
                    rr.append((p, st, t))
                if any("N" in t for _, _, t in rr): drop["N in read"] += 1; continue
                out = []
                for i, ((p, st, t), (Li, q, c5, c3, mm)) in enumerate(zip(rr, pr)):
                    b = np.frombuffer(t.encode(), dtype=np.uint8).copy()
                    if P:
                        pj, sj = rr[1 - i][0], rr[1 - i][1]; seg = sq[pj - off:pj - off + c3]
                        if c3: b[Li - c3:] = np.frombuffer((seg.translate(comp)[::-1] if sj == "+" else seg).encode(), dtype=np.uint8)
                        if c5: b[:c5] = B[rng.integers(0, 4, c5)]
                        m = np.zeros(Li, dtype=bool); m[mm] = True
                    else: m = rng.random(L) < a.err
                    if m.any(): b[m] = B[(np.searchsorted(B, b[m]) + rng.integers(1, 4, m.sum())) % 4]  # a different base
                    out.append((b.tobytes().decode(), q))
                if any("N" in r_ for r_, _ in out): drop["N in the partner sequence of a 3' clip"] += 1; continue
                rn = f"s{n}:{d['chrom']}:{rr[0][0] + 1}:{rr[0][1]}:{d['chrom']}:{rr[1][0] + 1}:{rr[1][1]}"
                for (r_, q), o in zip(out, (o1, o2)): o.write(f"@{rn}\n{r_}\n+\n{q}\n")
                n += 1; made += 1
    log(f"simulate: wrote {n} read pairs to {a.out}_R1/_R2.fastq.gz; dropped draws {dict(drop)}")

else:
    D = read_domains(a.domains)
    win = collections.defaultdict(list)
    for d in D:
        win[d["chrom"]] += [(d["phr_start"], d["phr_end"], d["end"], "phr"), (d["flank_start"], d["flank_end"], d["end"], "flank")]
    def where(c, p):
        for s0, e0, nm, w in win.get(c, []):
            if s0 <= p < e0: return nm, w
        return ".", "other"
    kept = collections.Counter(); mis = collections.Counter(); n = 0; correct = 0
    with open(a.avp) as f, open(a.out + ".reads.tsv", "w") as o:
        o.write("read\tside\ttrue_end\ttrue_window\tassigned_chrom\tassigned_pos\tassigned_end\tassigned_window\tcorrect\n")
        for l in f:
            x = l.rstrip("\n").split("\t"); t = x[0].split(":")
            truth = [(t[1], int(t[2])), (t[4], int(t[5]))]
            # allValidPairs does not keep R1/R2 order: take the side order with the smaller summed distance to the truth
            asg = [(x[1], int(x[2])), (x[4], int(x[5]))]
            dist = lambda u, v: abs(u[1] - v[1]) if u[0] == v[0] else 10 ** 12
            if dist(asg[1], truth[0]) + dist(asg[0], truth[1]) < dist(asg[0], truth[0]) + dist(asg[1], truth[1]): asg = asg[::-1]
            n += 1
            for side in (0, 1):
                te, tw = where(*truth[side]); ae, aw = where(*asg[side])
                ok = dist(asg[side], truth[side]) <= a.tol
                correct += ok
                o.write(f"{t[0]}\t{side + 1}\t{te}\t{tw}\t{asg[side][0]}\t{asg[side][1]}\t{ae}\t{aw}\t{int(ok)}\n")
                if tw == "phr":
                    kept[te] += 1
                    if aw == "phr" and ae != te: mis[(te, ae)] += 1
    with open(a.out + ".misassign.tsv", "w") as o:
        o.write("end_a\tend_b\treads_a_to_b\treads_a_kept\trate\n")
        for (ea, eb), v in sorted(mis.items()): o.write(f"{ea}\t{eb}\t{v}\t{kept[ea]}\t{v / kept[ea]:.6g}\n")
    log(f"evaluate: {n} valid pairs, {2 * n} reads, {correct} placed within tolerance ({correct / max(1, 2 * n):.4f});"
        f" PHR-window reads moved to another end's PHR window {sum(mis.values())} of {sum(kept.values())}")
