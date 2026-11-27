#!/usr/bin/env python3
"""Pedigree statistics on the drawn non-homologous runs of the final corrected assemblies (replaces the statistics of the
withdrawn odgi-untangle patch set: within-community fraction, assembly QC of called intervals, PAR1 and acrocentric support).

usage: pedigree_run_stats.py D K15_TSV SUMMARY_DIR OUT_PREFIX RUNS_TSV:LABEL:QUERY_SAMPLE:QUERY_HAP:DONOR_SAMPLE ...
  D: washu_v131 dir. K15_TSV: hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv (arm -> Leiden community).
  RUNS_TSV: fig5_runs.sh output (*.whole_genome_ribbon_runs.tsv); QUERY_HAP: 1, 2 or kid (KID.m query: hap2 for chr1 3 4 5 7 9
  12 14 15 18 21 22 X, hap1 otherwise, as fig5_inputs_final.sbatch builds it).
Arms: p if the run ends before the first active-HOR block of its own contig (annot_v11/SAMPLE.polished.cenSat.active_hor.bed;
  haplotype1/maternal = #1, haplotype2/paternal = #2), q if it starts after the last one, otherwise 'cen' (between blocks).
Within-community: query arm and donor arm in the same Leiden community (arms missing from K15_TSV: no community).
Null: 10,000 permutations of the donor arms over all runs (keeps the per-arm counts on both sides), seed 1.
Merqury QC: error k-mer starts (assembly 21-mers absent from the person's Illumina reads) inside the query and the donor
  interval, on the final assemblies: merqury tags g1v11b (PAN010, PAN028), g1v11f (PAN011), g1v11d (PAN027); expected =
  the haplotype's genome-wide rate x interval length (rates from SUMMARY_DIR/v11.SAMPLE.summary.tsv, g1_changes_report.py).
out: OUT_PREFIX.runs.tsv (one row per run), OUT_PREFIX.summary.tsv; log on stderr.
"""
import bisect, collections, csv, glob, gzip, math, os, random, sys
D, K15, SD, OUT = sys.argv[1:5]; specs = sys.argv[5:]
def log(m): print(m, file=sys.stderr, flush=True)
TAG = {"PAN010": "g1v11b", "PAN028": "g1v11b", "PAN011": "g1v11f", "PAN027": "g1v11d"}
KID2 = {f"chr{c}" for c in (1, 3, 4, 5, 7, 9, 12, 14, 15, 18, 21, 22, "X")}
com = {r["ChromArm"]: r["Community"] for r in csv.DictReader(open(K15), delimiter="\t")}
cen = {}; HN = {"haplotype1": "1", "maternal": "1", "haplotype2": "2", "paternal": "2"}
for sm in TAG:
    for l in open(f"{D}/annot_v11/{sm}.polished.cenSat.active_hor.bed"):
        f = l.split("\t"); _, ch, hn = f[0].split("."); k = f"{sm}#{HN[hn]}#{ch}"; a, b = int(f[1]), int(f[2])
        cen[k] = (min(a, cen.get(k, (a, b))[0]), max(b, cen.get(k, (a, b))[1]))
mq = {}
def kmers(sm, contig):
    h = contig.split("#")[1]
    if (sm, h) not in mq:
        p = collections.defaultdict(list)
        for l in gzip.open(f"{D}/merqury/{TAG[sm]}/{sm}/{sm}_hap{h}_only.bed.gz", "rt"):
            f = l.split("\t"); p[f[0]].append(int(f[1]))
        for v in p.values(): v.sort()
        mq[(sm, h)] = p; log(f"merqury {TAG[sm]} {sm} hap{h}: {sum(len(v) for v in p.values())} error k-mers")
    v = mq[(sm, h)].get(contig, []); return v
rate = {}
for f in glob.glob(f"{SD}/v11.*.summary.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"): rate[(r["sample"], r["haplotype"])] = int(r["merqury_kmers_final"]) / int(r["final_bp"])
if len(rate) != 8: sys.exit(f"ERROR {len(rate)} haplotype rates, expected 8")
def arm(contig, a, b):
    s, e = cen[contig]; chrom = contig.split("#")[2]
    return f"{chrom}_p" if b <= s else f"{chrom}_q" if a >= e else f"{chrom}_cen"
rows = []
for sp in specs:
    path, lab, qs, qh, ds = sp.split(":")
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = r["query_chrom"]; h = ("2" if c in KID2 else "1") if qh == "kid" else qh
        qc = f"{qs}#{h}#{c}"; qa, qb = int(r["query_start"]), int(r["query_end"])
        dseq = r["donor_seq"]; dh = dseq.split("#joint#h")[1].split("_")[0]; dchr = dseq.split("_")[-1]
        dc = f"{ds}#{dh}#{dchr}"; da, db = int(r["donor_start"]), int(r["donor_end"])
        qarm, darm = arm(qc, qa, qb), arm(dc, da, db)
        kq = kmers(qs, qc); kd = kmers(ds, dc)
        oq = bisect.bisect_left(kq, qb) - bisect.bisect_left(kq, qa); od = bisect.bisect_left(kd, db) - bisect.bisect_left(kd, da)
        eq = rate[(qs, h)] * (qb - qa); ed = rate[(ds, dh)] * (db - da)
        rows.append(dict(transmission=lab, run_id=r["run_id"], category=r["category"], query=f"{qc}:{qa}-{qb}", donor=f"{dc}:{da}-{db}",
                         bp=int(r["bp"]), query_arm=qarm, donor_arm=darm, query_community=com.get(qarm, "none"),
                         donor_community=com.get(darm, "none"), within=int(com.get(qarm, "a") == com.get(darm, "b")),
                         mean_same_identity=r["mean_same_identity"], mean_inter_identity=r["mean_inter_identity"],
                         kmers_query=oq, kmers_query_expected=f"{eq:.2f}", kmers_donor=od, kmers_donor_expected=f"{ed:.2f}"))
    log(f"{lab}: {sum(1 for x in rows if x['transmission'] == lab)} runs from {path}")
with open(f"{OUT}.runs.tsv", "w") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
n = len(rows); k = sum(x["within"] for x in rows)
def wilson(k, n, z=1.959964):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h
random.seed(1); qa_ = [x["query_arm"] for x in rows]; da_ = [x["donor_arm"] for x in rows]; null = []
for _ in range(10000):
    random.shuffle(da_); null.append(sum(com.get(a, "a") == com.get(b, "b") for a, b in zip(qa_, da_)) / n)
null.sort(); ge = sum(v >= k / n for v in null)
lo, hi = wilson(k, n)
out = [("runs", n), ("runs_bp", sum(x["bp"] for x in rows)), ("within_community", k), ("within_fraction", f"{k / n:.4f}"),
       ("within_wilson95", f"{lo:.4f}-{hi:.4f}"), ("null_mean", f"{sum(null) / len(null):.4f}"),
       ("null_95", f"{null[249]:.4f}-{null[9749]:.4f}"), ("null_p_ge_observed", f"{(ge + 1) / 10001:.6f}")]
for lab in sorted({x["transmission"] for x in rows}):
    s = [x for x in rows if x["transmission"] == lab]
    out += [(f"{lab}_runs", len(s)), (f"{lab}_within", sum(x["within"] for x in s)),
            (f"{lab}_categories", ";".join(f"{c}={n_}" for c, n_ in sorted(collections.Counter(x["category"] for x in s).items())))]
for side in ("query", "donor"):
    o_ = sum(x[f"kmers_{side}"] for x in rows); e_ = sum(float(x[f"kmers_{side}_expected"]) for x in rows)
    out += [(f"merqury_kmers_{side}", o_), (f"merqury_kmers_{side}_expected", f"{e_:.1f}"),
            (f"runs_with_kmers_{side}", sum(x[f"kmers_{side}"] > 0 for x in rows))]
acro = [x for x in rows if x["category"] == "acro_acro"]; fl = collections.defaultdict(list)
for x in acro: fl[(x["transmission"], x["query_arm"])].append(x["donor_arm"])
out += [("acro_flanks", len(fl)), ("acro_flank_runs_min_max", f"{min(len(v) for v in fl.values())}-{max(len(v) for v in fl.values())}" if fl else "NA"),
        ("acro_flank_donor_arms_min_max", f"{min(len(set(v)) for v in fl.values())}-{max(len(set(v)) for v in fl.values())}" if fl else "NA")]
with open(f"{OUT}.summary.tsv", "w") as o:
    o.write("statistic\tvalue\n"); o.writelines(f"{a}\t{b}\n" for a, b in out)
for a, b in out: log(f"{a}\t{b}")
for x in rows:
    if not x["within"] or x["query_arm"].endswith("_cen") or x["donor_arm"].endswith("_cen"):
        log(f"NOT WITHIN or CEN: {x['transmission']} {x['run_id']} {x['category']} {x['query']} ({x['query_arm']}, {x['query_community']}) -> {x['donor']} ({x['donor_arm']}, {x['donor_community']})")
