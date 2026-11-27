#!/usr/bin/env python3
"""Read-level Pore-C / CiFi simulation for referee R1-M4 / R3-M3: concatemers whose monomers all come from one
chromosome end, so any contact between two ends after mapping is a placement error.

Null model: for each end (PHR window plus 100 kb flank, as sim_hic_null.py domains), a read has K >= 2 monomers
(K - 2 ~ Poisson(MEAN - 2), at most 20). The first monomer is the restriction fragment (FRAGMENTS_BED, from cooler
digest) containing a uniform position of the end domain; each other one is the fragment containing a position at
distance s from it, P(s) proportional to 1/s on [MIN_S, MAX_S], inside the extended domain. Monomers are stripped of
the enzyme site at their ends, oriented at random and joined by one site, as in a ligation product, so that pore-c-py
digest splits them again. Errors per base: substitution, insertion, deletion at SUB, INS, DEL, or, with --profile
(read_profile.py porec), each read takes the substitution, insertion and deletion rates of a random real read and each
monomer gets random bases over the fractions of its length that a random real monomer had soft-clipped at its 5' and 3'
ends. Read names s<n>; the truth (end and monomer coordinates) goes to OUT.truth.tsv.

subcommands
  simulate DOMAINS_TSV EXTENDED_FA FRAGMENTS_BED SITE OUT [--reads-per-kb 20] [--mean-monomers 4] [--sub 0.005]
           [--ins 0.0025] [--dele 0.0025] [--min-s 1000] [--max-s 500000] [--seed 20261126] [--profile TSV_GZ]
           writes OUT.fastq.gz and OUT.truth.tsv
  evaluate DOMAINS_TSV TRUTH_TSV PAIRS OUT [--tol 1000]
           per pair side of the pairs file (readID chrom1 pos1 chrom2 pos2 ...): the read's true end, the end window
           (PHR, flank, other) where it was placed, and whether it lies within TOL bp of a true monomer of the read.
           Writes OUT.misassign.tsv (true end A, placed in the PHR window of end B != A: sides, all sides of A's
           reads, rate) and logs the totals."""
import argparse, bisect, collections, gzip, math, re, sys
import numpy as np

log = lambda *m: print(*m, file=sys.stderr, flush=True)
ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
s = sp.add_parser("simulate")
for k in ("domains", "fa", "frags", "site", "out"): s.add_argument(k)
for k, v in (("--reads-per-kb", 20.0), ("--mean-monomers", 4.0), ("--sub", 0.005), ("--ins", 0.0025), ("--dele", 0.0025)): s.add_argument(k, type=float, default=v)
for k, v in (("--min-s", 1000), ("--max-s", 500000), ("--seed", 20261126)): s.add_argument(k, type=int, default=v)
s.add_argument("--profile", help="real-read profiles (read_profile.py porec); replaces --sub, --ins and --dele")
e = sp.add_parser("evaluate")
for k in ("domains", "truth", "pairs", "out"): e.add_argument(k)
e.add_argument("--tol", type=int, default=1000)
a = ap.parse_args()
COLS = ["end", "chrom", "arm", "phr_start", "phr_end", "flank_start", "flank_end", "ext_start", "ext_end"]
D = [dict(zip(COLS, l.rstrip("\n").split("\t"))) for l in open(a.domains) if not l.startswith("end\t")]
for d in D:
    for k in COLS[3:]: d[k] = int(d[k])

if a.cmd == "simulate":
    rng = np.random.default_rng(a.seed); comp = str.maketrans("ACGTN", "TGCAN"); site = a.site.upper()
    seqs, name, buf = {}, None, []
    with gzip.open(a.fa, "rt") if a.fa.endswith(".gz") else open(a.fa) as f:
        for l in f:
            if l.startswith(">"):
                if name: seqs[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    seqs[name] = "".join(buf).upper()
    want = collections.defaultdict(list)
    for d in D: want[d["chrom"]].append((d["ext_start"], d["ext_end"]))
    fr = collections.defaultdict(list)
    with open(a.frags) as f:
        for l in f:
            c, s0, e0 = l.split("\t")[:3]
            if c in want:
                s0, e0 = int(s0), int(e0)
                if any(s0 >= x and e0 <= y for x, y in want[c]): fr[c].append((s0, e0))
    for c in fr: fr[c].sort()
    fs = {c: [x[0] for x in v] for c, v in fr.items()}
    log(f"simulate: {len(D)} ends, {sum(map(len, fr.values()))} {site} fragments inside the extended domains; {a.reads_per_kb} reads"
        f" per kb, mean monomers {a.mean_monomers}, errors sub {a.sub} ins {a.ins} del {a.dele}, s in [{a.min_s}, {a.max_s}], seed {a.seed}")
    B = np.frombuffer(b"ACGT", dtype=np.uint8)
    PR, PC = None, None
    if a.profile:
        PR, PC = [], []
        with gzip.open(a.profile, "rt") as f:
            if not f.readline().startswith("read\t"): sys.exit("FATAL profile header")
            for l in f:
                x = l.rstrip("\n").split("\t"); PR.append((float(x[2]), float(x[3]), float(x[4])))
                PC += [tuple(map(float, c.split("/"))) for c in x[5].split(",")]
        log(f"profiles {a.profile}: {len(PR)} reads, {len(PC)} monomers; mean rates sub {np.mean([r[0] for r in PR]):.5f},"
            f" ins {np.mean([r[1] for r in PR]):.5f}, del {np.mean([r[2] for r in PR]):.5f}; mean clipped fraction 5' {np.mean([c[0] for c in PC]):.4f}, 3' {np.mean([c[1] for c in PC]):.4f}")
    def mutate(t, rs, ri, rd):
        b = np.frombuffer(t.encode(), dtype=np.uint8).copy(); u = rng.random(len(b))
        sub = u < rs; b[sub] = B[(np.searchsorted(B, b[sub]) + rng.integers(1, 4, sub.sum())) % 4]
        keep = ~((u >= rs) & (u < rs + rd)); ins = (u >= rs + rd) & (u < rs + rd + ri)
        out = np.insert(b, np.nonzero(ins)[0] + 1, B[rng.integers(0, 4, ins.sum())])
        mask = np.insert(keep, np.nonzero(ins)[0] + 1, True)
        return out[mask].tobytes().decode()
    n, drop, nmono = 0, collections.Counter(), 0
    with gzip.open(a.out + ".fastq.gz", "wt", compresslevel=1) as o, open(a.out + ".truth.tsv", "w") as t:
        t.write("read\tend\tmonomers\n")
        for d in D:
            c, off = d["chrom"], d["ext_start"]; sq = seqs[f"{c}:{d['ext_start'] + 1}-{d['ext_end']}"]
            if len(sq) != d["ext_end"] - d["ext_start"]: sys.exit(f"FATAL length of {d['end']}")
            lo, hi = min(d["phr_start"], d["flank_start"]), max(d["phr_end"], d["flank_end"])
            def frag(p):
                j = bisect.bisect_right(fs[c], p) - 1
                return fr[c][j] if j >= 0 and fr[c][j][0] <= p < fr[c][j][1] else None
            for _ in range(int(round((hi - lo) / 1000 * a.reads_per_kb))):
                k = min(20, 2 + int(rng.poisson(max(0.0, a.mean_monomers - 2))))
                f0 = frag(int(rng.integers(lo, hi)))
                if f0 is None: drop["no fragment at the first position"] += 1; continue
                mono = [f0]
                while len(mono) < k:
                    sd = int(round(math.exp(rng.uniform(math.log(a.min_s), math.log(a.max_s)))))
                    p = f0[0] + (sd if rng.random() < 0.5 else -sd)
                    if p < d["ext_start"] or p >= d["ext_end"]: drop["partner outside extended domain"] += 1; continue
                    f1 = frag(p)
                    if f1 is None or f1 in mono: drop["no new fragment at the partner position"] += 1; continue
                    mono.append(f1)
                parts, desc = [], []
                for s0, e0 in mono:
                    x = sq[s0 - off:e0 - off]
                    x = x[len(site):] if x.startswith(site) else x
                    x = x[:-len(site)] if x.endswith(site) else x
                    st = "+" if rng.random() < 0.5 else "-"
                    x = x if st == "+" else x.translate(comp)[::-1]
                    if PC:
                        f5, f3 = PC[int(rng.integers(len(PC)))]; k5, k3 = int(round(f5 * len(x))), int(round(f3 * len(x)))
                        if k5 + k3 < len(x): x = "".join("ACGT"[v] for v in rng.integers(0, 4, k5)) + x[k5:len(x) - k3] + "".join("ACGT"[v] for v in rng.integers(0, 4, k3))
                        else: x = "".join("ACGT"[v] for v in rng.integers(0, 4, len(x)))
                    parts.append(x); desc.append(f"{c}:{s0}-{e0}:{st}")
                r = mutate(site.join(parts), *(PR[int(rng.integers(len(PR)))] if PR else (a.sub, a.ins, a.dele)))
                if "N" in r: drop["N in read"] += 1; continue
                o.write(f"@s{n}\n{r}\n+\n{'5' * len(r)}\n"); t.write(f"s{n}\t{d['end']}\t{','.join(desc)}\n")
                n += 1; nmono += len(mono)
    log(f"simulate: wrote {n} reads with {nmono} monomers to {a.out}.fastq.gz; dropped draws {dict(drop)}")

else:
    win = collections.defaultdict(list)
    for d in D: win[d["chrom"]] += [(d["phr_start"], d["phr_end"], d["end"], "phr"), (d["flank_start"], d["flank_end"], d["end"], "flank")]
    def where(c, p):
        for s0, e0, nm, w in win.get(c, []):
            if s0 <= p < e0: return nm, w
        return ".", "other"
    arm = lambda x: re.sub(r"_(MATERNAL|PATERNAL)", "", x)
    tr = {}
    for l in open(a.truth):
        if l.startswith("read\t"): continue
        r, en, m = l.rstrip("\n").split("\t")
        tr[r] = (en, [(x.split(":")[0], int(x.split(":")[1].split("-")[0]), int(x.split(":")[1].split("-")[1])) for x in m.split(",")])
    sides, mis, homo, ok, tot = collections.Counter(), collections.Counter(), 0, 0, 0
    opener = gzip.open if a.pairs.endswith(".gz") else open
    with opener(a.pairs, "rt") as f:
        for l in f:
            if l.startswith("#"): continue
            x = l.rstrip("\n").split("\t"); en, mono = tr[x[0]]
            for c, p in ((x[1], int(x[2])), (x[3], int(x[4]))):
                if c == "!": continue
                tot += 1; sides[en] += 1
                ok += any(c == mc and ms - a.tol <= p <= me + a.tol for mc, ms, me in mono)
                ae, aw = where(c, p)
                if aw == "phr" and ae != en:
                    if arm(ae) == arm(en): homo += 1
                    else: mis[(en, ae)] += 1
    with open(a.out + ".misassign.tsv", "w") as o:
        o.write("end_a\tend_b\tsides_a_to_b\tsides_a\trate\n")
        for (ea, eb), v in sorted(mis.items()): o.write(f"{ea}\t{eb}\t{v}\t{sides[ea]}\t{v / sides[ea]:.6g}\n")
    log(f"evaluate: {tot} pair sides; {ok} within {a.tol} bp of a true monomer ({ok / max(1, tot):.4f}); placed in the PHR window of"
        f" the homologous end {homo}; of another arm {sum(mis.values())} ({sum(mis.values()) / max(1, tot):.4f})")
