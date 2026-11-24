#!/usr/bin/env python3
"""Plant known conversion tracts into real subtelomeric windows (R1-M1 test 3b, known answer on real data).

Windows of --recipient-arm (CHM13 excluded) are visited in a seeded random order. For each window we read
its target-centric PAF from the PHR all-vs-all and try to place one tract per --lengths value. A site for
length L is accepted when
  (i)   one record of a donor-arm window of ANOTHER individual covers L plus --flank bp on each side,
        with no insertion or deletion >= --max-indel bp inside that span;
  (ii)  >= --min-informative of the 1 kb bins of the site have >= 20 alleles and >= 20 donor-arm copies
        (best record per query, >= 75% of the bin aligned, as in conversion_scan.py) and median d(donor
        copies) - median d(alleles) >= --min-gap, so a converted copy would differ from its own arm;
  (iii) the site is >= --spacing bp from any other planted site.
Each length is tried on its own (longest first); the --index-th window where >= --min-placed lengths get a
site is used. Inside each site the recipient sequence is
replaced by the aligned donor sequence. OUTPUT: <out>.fa (the modified window, same name) and
<out>.plants.tsv (planted intervals in the new and old coordinates, donor, substitutions introduced).
"""
import argparse, gzip, random, re, subprocess, sys, time
from collections import defaultdict
sys.path.insert(0, __import__("os").path.dirname(__file__))
from conversion_scan import CIGAR, arm_of, bin_counts, med

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def fetch(samtools, fasta, name):
    out = subprocess.run([samtools, "faidx", fasta, name], check=True, capture_output=True, text=True).stdout
    seq = "".join(out.split("\n")[1:])
    if not seq:
        raise SystemExit(f"FATAL: empty sequence for {name}")
    return seq


def walk(rec):
    """[(op, n, t, k)]: op runs with target position t and offset k into the aligned query (as aligned)."""
    t, k, out = rec["tstart"], 0, []
    for n, op in CIGAR.findall(rec["cg"]):
        n = int(n)
        out.append((op, n, t, k))
        t += n if op in "=XD" else 0
        k += n if op in "=XI" else 0
    return out


def t2k(ops, t):
    for op, n, t0, k0 in ops:
        if op in "=X" and t0 <= t <= t0 + n:
            return k0 + (t - t0)
    raise SystemExit(f"FATAL: target position {t} is not in a =/X run")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", required=True, help="target-centric PAF paths, one per line")
    ap.add_argument("--fasta", required=True, help="bgzipped windows FASTA with .fai and .gzi")
    ap.add_argument("--recipient-arm", required=True)
    ap.add_argument("--donor-arm", required=True)
    ap.add_argument("--lengths", default="500,1000,2000,5000,10000")
    ap.add_argument("--index", type=int, default=1, help="use the index-th window where all lengths fit")
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--spacing", type=int, default=20000)
    ap.add_argument("--max-indel", type=int, default=50)
    ap.add_argument("--min-gap", type=float, default=0.003)
    ap.add_argument("--min-informative", type=float, default=0.75)
    ap.add_argument("--min-placed", type=int, default=3)
    ap.add_argument("--min-copies", type=int, default=20)
    ap.add_argument("--max-windows", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    L = sorted((int(x) for x in a.lengths.split(",")), reverse=True)
    log(f"recipient {a.recipient_arm}, donor {a.donor_arm}, lengths {L}, index {a.index}, flank {a.flank}, "
        f"spacing {a.spacing}, max_indel {a.max_indel}, min_gap {a.min_gap}, min_informative {a.min_informative}, "
        f"min_placed {a.min_placed}, seed {a.seed}")
    paths = [p.strip() for p in open(a.list) if p.strip().endswith(f"_{a.recipient_arm}.paf.gz") and "CHM13" not in p]
    random.Random(a.seed).shuffle(paths)
    log(f"{len(paths)} candidate {a.recipient_arm} windows")
    hits = 0
    for p in paths[:a.max_windows]:
        recs, best = [], defaultdict(dict)
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                q, t = f[0], f[5]
                cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
                if q == t or cg is None:
                    continue
                for b, (m, x) in bin_counts(int(f[7]), cg, 1000).items():
                    if m + x > best[b].get(q, (0, 0))[0]:
                        best[b][q] = (m + x, x)
                if arm_of(q)[1] == a.donor_arm and q.split("#")[0] != t.split("#")[0]:
                    recs.append(dict(q=q, qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=t,
                                     tstart=int(f[7]), tend=int(f[8]), cg=cg))
        target = t
        ok_bin = {}
        for b, d in best.items():
            al = [x / n for q, (n, x) in d.items() if n >= 750 and arm_of(q)[1] == a.recipient_arm]
            do = [x / n for q, (n, x) in d.items() if n >= 750 and arm_of(q)[1] == a.donor_arm]
            if len(al) >= a.min_copies and len(do) >= a.min_copies and med(do) - med(al) >= a.min_gap:
                ok_bin[b] = (med(al), med(do))
        recs.sort(key=lambda r: r["tend"] - r["tstart"], reverse=True)
        placed = []
        for ln in L:
            site = None
            for r in recs[:50]:
                ops = walk(r)
                big = [(t0, t0 + (n if op == "D" else 0)) for op, n, t0, k0 in ops if op in "ID" and n >= a.max_indel]
                for s in range(r["tstart"] + a.flank, r["tend"] - a.flank - ln, 1000):
                    e = s + ln
                    if any(abs(s - ps) < a.spacing + ln + pl for ps, pl, _ in placed):
                        continue
                    if any(s - a.flank <= y and x <= e + a.flank for x, y in big):
                        continue
                    bs = range(s // 1000, (e - 1) // 1000 + 1)
                    if sum(b in ok_bin for b in bs) < a.min_informative * len(bs):
                        continue
                    site = (s, ln, r)
                    break
                if site:
                    break
            if site:
                placed.append(site)
        log(f"{p.split('/')[-1]}: {len(recs)} donor records, {len(ok_bin)} informative 1 kb bins, "
            f"{len(placed)}/{len(L)} lengths placed")
        if len(placed) >= a.min_placed:
            hits += 1
            if hits == a.index:
                break
    else:
        raise SystemExit(f"FATAL: only {hits} windows with >= {a.min_placed} lengths placed, index {a.index} requested")
    rseq = fetch(a.samtools, a.fasta, target)
    rows = []
    for s, ln, r in placed:
        ops = walk(r)
        k0, k1 = t2k(ops, s), t2k(ops, s + ln)
        dseq = fetch(a.samtools, a.fasta, r["q"])
        if len(dseq) != r["qlen"]:
            raise SystemExit(f"FATAL: {r['q']} length {len(dseq)} != PAF {r['qlen']}")
        if r["strand"] == "+":
            seg = dseq[r["qs"] + k0:r["qs"] + k1]
        else:
            seg = dseq[r["qe"] - k1:r["qe"] - k0].translate(COMP)[::-1]
        sub = sum(min(t0 + n, s + ln) - max(t0, s) for op, n, t0, _ in ops if op == "X" and t0 < s + ln and t0 + n > s)
        ind = sum(1 for op, n, t0, _ in ops if op in "ID" and s <= t0 < s + ln)
        bins = [ok_bin[b] for b in range(s // 1000, (s + ln - 1) // 1000 + 1) if b in ok_bin]
        rows.append(dict(old_start=s, old_end=s + ln, seg=seg, donor=r["q"], strand=r["strand"], subst=sub, indels=ind,
                         d_allele=sum(x for x, _ in bins) / len(bins), d_donor=sum(y for _, y in bins) / len(bins)))
    rows.sort(key=lambda x: x["old_start"])
    new, prev, shift = [], 0, 0
    for x in rows:
        new.append(rseq[prev:x["old_start"]])
        x["new_start"] = x["old_start"] + shift
        new.append(x["seg"])
        shift += len(x["seg"]) - (x["old_end"] - x["old_start"])
        x["new_end"] = x["old_end"] + shift
        prev = x["old_end"]
    new.append(rseq[prev:])
    mseq = "".join(new)
    for x in rows:
        if mseq[x["new_start"]:x["new_end"]] != x["seg"]:
            raise SystemExit("FATAL: planted segment not found at its new coordinates")
    with open(a.out + ".fa", "w") as fh:
        fh.write(f">{target}\n")
        for i in range(0, len(mseq), 80):
            fh.write(mseq[i:i + 80] + "\n")
    cols = ["new_start", "new_end", "old_start", "old_end", "donor", "strand", "subst", "indels", "d_allele", "d_donor"]
    with open(a.out + ".plants.tsv", "w") as fh:
        fh.write("target\tlength\t" + "\t".join(cols) + "\n")
        for x in rows:
            fh.write(f"{target}\t{x['old_end'] - x['old_start']}\t" + "\t".join(
                f"{x[c]:.5f}" if isinstance(x[c], float) else str(x[c]) for c in cols) + "\n")
            log(f"planted {x['old_end'] - x['old_start']} bp at {x['old_start']}-{x['old_end']} -> "
                f"{x['new_start']}-{x['new_end']} from {x['donor']} ({x['strand']}): {x['subst']} substitutions, "
                f"{x['indels']} indels; bins d_allele {x['d_allele']:.4f}, d_donor {x['d_donor']:.4f}")
    log(f"wrote {a.out}.fa ({len(rseq)} -> {len(mseq)} bp) and {a.out}.plants.tsv")


if __name__ == "__main__":
    main()
