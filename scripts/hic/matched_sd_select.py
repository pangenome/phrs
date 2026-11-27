#!/usr/bin/env python3
"""Matched non-subtelomeric segmental-duplication (SD) control for referee R1-M5 / R1-M9.

From the CHM13 SEDEF table (one row per pair and orientation; kept: original == True, inter-chromosomal, chr1-22/X/Y)
we merge the copies of each chromosome into loci and aggregate the alignments into locus pairs:
  identity  = aligned-length-weighted fracMatch of the alignments linking the two loci,
  length    = summed aligned length (aln_len),
  copies    = mean number of distinct other chromosomes that the two loci share SDs with (as column 4 of the PHR BED).
Classes: a locus is PHR if it overlaps a PHR interval (PHR_BED); non-subtelomeric if it lies at least END_DIST from
both chromosome ends; among those, acrocentric if it lies on the p arm (ends before the centromere) of a chromosome in
NOR_CHROMS, pericentromeric if it overlaps the centromere (CEN_BED) extended by PERI_DIST, else interstitial. The SEDEF
telo/peri/acro flags are not used (the HG002 v1.1 file has 0 acro and 53 peri loci, against 44 and 693 in CHM13). PHR pairs = both loci PHR; control candidates = both loci non-subtelomeric.
Each PHR pair gets K controls by greedy nearest-neighbour matching without replacement on (identity, log10 length,
log2 copies), standardized by the PHR-pair SD, within a caliper; PHR pairs are visited in a seeded random order.
Writes OUT.pairs.tsv (PHR pairs and their controls) and OUT.balance.tsv (means and standardized mean differences
before and after matching), and logs counts kept and dropped with the reason.
usage: matched_sd_select.py SEDEF_FULL_BED CHROM_SIZES PHR_BED CEN_BED NOR_CHROMS OUT [--end-dist 1000000]
       [--peri-dist 5000000] [--k 3] [--caliper 0.25] [--seed 20261126] [--no-acro]"""
import argparse, collections, math, random, sys

ap = argparse.ArgumentParser()
ap.add_argument("sedef"); ap.add_argument("sizes"); ap.add_argument("phr"); ap.add_argument("cen"); ap.add_argument("nor"); ap.add_argument("out")
ap.add_argument("--peri-dist", type=int, default=5000000)
ap.add_argument("--end-dist", type=int, default=1000000); ap.add_argument("--k", type=int, default=3)
ap.add_argument("--caliper", type=float, default=0.25); ap.add_argument("--seed", type=int, default=20261126)
ap.add_argument("--no-acro", action="store_true", help="drop control loci with the SEDEF acro flag")
a = ap.parse_args()
log = lambda *m: print(*m, file=sys.stderr, flush=True)
SEDEF_COLS = ("chr1 start1 end1 name fakeScore strand1 start1 end1 color chr2 start2 end2 score strand2 max_len aln_len indel_a "
              "indel_b alnB matchB mismatchB transitionsB transversions fracMatch fracMatchIndel jck k2K aln_gaps uppercaseA "
              "uppercaseB uppercaseMatches aln_matches aln_mismatches aln_gap_bases filter_score sat_bases unique_id original "
              "telo peri acro telo2 peri2 acro2").split()  # the CHM13 header; the HG002 v1.1 file has the same 44 columns, no header
base = lambda c: c.rsplit("_", 1)[0] if c.endswith(("_MATERNAL", "_PATERNAL")) else c  # chr1_MATERNAL and chr1_PATERNAL are homologs
CHR = {f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]}
size = {l.split()[0]: int(l.split()[1]) for l in open(a.sizes) if base(l.split()[0]) in CHR}
CHR = set(size)  # contigs of the contact reference (CHM13 Hi-C: chrY absent); SDs elsewhere are dropped
phr = collections.defaultdict(list)
for l in open(a.phr):
    c, s, e = l.split("\t")[:3]; phr[c].append((int(s), int(e)))
cen = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in open(a.cen)}
nor = {l.strip() for l in open(a.nor) if l.strip()}
miss = [c for c in size if c not in cen]
if miss: sys.exit(f"FATAL no centromere for {miss}")
log(f"inputs: {a.sedef}; {len(size)} chromosomes; {sum(map(len, phr.values()))} PHR intervals; {len(cen)} centromeres; NOR chromosomes {sorted(nor)}; end distance {a.end_dist}; peri distance {a.peri_dist};"
    f" k {a.k}; caliper {a.caliper} SD; seed {a.seed}; no-acro {a.no_acro}")

rows, drop = [], collections.Counter()
with open(a.sedef) as f:
    first = f.readline()
    h = first.lstrip("#").rstrip("\n").split("\t") if first.startswith("#") else SEDEF_COLS
    if len(h) != 44: sys.exit(f"FATAL expected 44 SEDEF columns, got {len(h)}")
    ix = {k: i for i, k in reversed(list(enumerate(h)))}
    for l in ([] if first.startswith("#") else [first]) + list(f):
        x = l.rstrip("\n").split("\t")
        if x[ix["original"]] != "True": drop["mirror row"] += 1; continue
        c1, c2 = x[0], x[ix["chr2"]]
        if c1 not in CHR or c2 not in CHR: drop["chromosome not in CHROM_SIZES (chr1-22/X/Y only)"] += 1; continue
        if base(c1) == base(c2): drop["same chromosome (incl. the other haplotype's homolog)"] += 1; continue
        rows.append((c1, int(x[1]), int(x[2]), c2, int(x[ix["start2"]]), int(x[ix["end2"]]), int(x[ix["aln_len"]]),
                     float(x[ix["fracMatch"]]), x[39] == "1", x[40] == "1", x[42] == "1", x[43] == "1"))
log(f"SEDEF pairs kept {len(rows)}; dropped {dict(drop)}")

# merge copies into loci per chromosome
iv = collections.defaultdict(list)
for r in rows:
    iv[r[0]].append((r[1], r[2])); iv[r[3]].append((r[4], r[5]))
loci, lid = {}, collections.defaultdict(list)  # loci[id] = (chr, start, end); lid[chr] = sorted (start, end, id)
for c, v in iv.items():
    v.sort(); cs, ce = v[0]
    for s, e in v[1:] + [(math.inf, math.inf)]:
        if s > ce:
            i = len(loci); loci[i] = (c, cs, ce); lid[c].append((cs, ce, i)); cs, ce = s, e
        else: ce = max(ce, e)
import bisect
starts = {c: [t[0] for t in v] for c, v in lid.items()}
def locus(c, s):
    j = bisect.bisect_right(starts[c], s) - 1; cs, ce, i = lid[c][j]; assert cs <= s < ce or cs == s; return i
peri, acro = collections.defaultdict(bool), collections.defaultdict(bool)
pairs = collections.defaultdict(lambda: [0, 0.0, 0])  # (la, lb) -> [aln_len, sum aln_len * fracMatch, n]
partners = collections.defaultdict(set)
for c1, s1, e1, c2, s2, e2, al, fm, p1, a1, p2, a2 in rows:
    la, lb = locus(c1, s1), locus(c2, s2)
    peri[la] |= p1; acro[la] |= a1; peri[lb] |= p2; acro[lb] |= a2
    k = (min(la, lb), max(la, lb)); pairs[k][0] += al; pairs[k][1] += al * fm; pairs[k][2] += 1
    partners[la].add(lb); partners[lb].add(la)
log(f"loci {len(loci)}; locus pairs {len(pairs)}")

def cls(i):
    c, s, e = loci[i]
    if any(s < pe and ps < e for ps, pe in phr.get(c, [])): return "PHR"
    if s >= a.end_dist and e <= size[c] - a.end_dist:
        cs, ce = cen[c]
        if c in nor and e <= cs: return "acro"
        return "peri" if s < ce + a.peri_dist and cs - a.peri_dist < e else "interstitial"
    return "end_other"
lc = {i: cls(i) for i in loci}
log(f"locus classes: {dict(collections.Counter(lc.values()))}")
pchr = {i: {loci[j][0] for j in v} for i, v in partners.items()}
feat = lambda k: (pairs[k][1] / pairs[k][0], math.log10(pairs[k][0]), math.log2((len(pchr[k[0]]) + len(pchr[k[1]])) / 2))
phrp = [k for k in pairs if lc[k[0]] == "PHR" and lc[k[1]] == "PHR"]
okc = {"peri", "interstitial"} | (set() if a.no_acro else {"acro"})
cand = [k for k in pairs if lc[k[0]] in okc and lc[k[1]] in okc]
sp = set(phrp); mixed = sum(1 for k in pairs if k not in sp and (lc[k[0]] == "PHR" or lc[k[1]] == "PHR"))
log(f"PHR-PHR locus pairs {len(phrp)}; non-subtelomeric candidate pairs {len(cand)}; PHR pairs with a non-PHR partner {mixed} (not used)")
if not phrp or not cand: sys.exit("FATAL: no PHR pairs or no candidates")
F = {k: feat(k) for k in phrp + cand}
mu = [sum(F[k][j] for k in phrp) / len(phrp) for j in range(3)]
sd = [math.sqrt(sum((F[k][j] - mu[j]) ** 2 for k in phrp) / max(1, len(phrp) - 1)) or 1.0 for j in range(3)]
rng = random.Random(a.seed); order = phrp[:]; rng.shuffle(order); free = set(cand); match = {}
for p in order:
    d = sorted((max(abs(F[p][j] - F[q][j]) / sd[j] for j in range(3)), q) for q in free)
    match[p] = [q for dist, q in d[:a.k] if dist <= a.caliper]; free -= set(match[p])
nm = collections.Counter(len(v) for v in match.values())
log(f"matched: PHR pairs by number of controls {dict(sorted(nm.items()))}; controls used {sum(nm[k] * k for k in nm)}")
ctrl = [q for v in match.values() for q in v]; mphr = [p for p, v in match.items() if v]; unm = [p for p, v in match.items() if not v]
with open(a.out + ".balance.tsv", "w") as o:
    o.write("set\tn\tidentity_mean\tlog10_len_mean\tlog2_copies_mean\tsmd_identity\tsmd_log10_len\tsmd_log2_copies\n")
    for name, S in (("PHR_all", phrp), ("PHR_matched", mphr), ("PHR_unmatched", unm), ("candidates", cand), ("matched_controls", ctrl)):
        m = [sum(F[k][j] for k in S) / len(S) for j in range(3)] if S else [math.nan] * 3
        o.write(f"{name}\t{len(S)}\t" + "\t".join(f"{v:.4f}" for v in m) + "\t" + "\t".join(f"{(m[j] - mu[j]) / sd[j]:.3f}" for j in range(3)) + "\n")
fmt = lambda i: f"{loci[i][0]}\t{loci[i][1]}\t{loci[i][2]}\t{lc[i]}\t{len(pchr[i])}"
with open(a.out + ".pairs.tsv", "w") as o:
    o.write("set\tphr_pair\tchrA\tstartA\tendA\tclassA\tpartner_chrsA\tchrB\tstartB\tendB\tclassB\tpartner_chrsB\taln_len\tidentity\tn_alignments\n")
    for n, p in enumerate(sorted(phrp, key=lambda k: loci[k[0]])):
        for setn, k in [("PHR", p)] + [("control", q) for q in match[p]]:
            o.write(f"{setn}\t{n}\t{fmt(k[0])}\t{fmt(k[1])}\t{pairs[k][0]}\t{pairs[k][1] / pairs[k][0]:.6f}\t{pairs[k][2]}\n")
log(f"wrote {a.out}.pairs.tsv and {a.out}.balance.tsv")
