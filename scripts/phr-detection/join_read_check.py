#!/usr/bin/env python3
"""Long-read test of candidate joins (flank_joins_coords.py TSV) with the HPRC release-2 read alignments to the same
assemblies, read over HTTPS from S3 (Flagger supplementary: read_alignments/ont/*.bam, else hifi/*.bam).

Per feature [start, end] (contig, 1-based) we fetch primary alignments in [start - --window, end + --window] and count:
  span      MAPQ >= --mapq, start < start - --anchor and end > end + --anchor (the two sides are adjacent in the reads)
  stop      end from --near bp before the feature to its end with a right soft clip >= --min-clip, or start from its start to
            --near bp after it with a left clip >= --min-clip; for each, the share of the clipped bases in telomeric 6-mers
  cover     reads with >= 1 kb on the 20 kb before the feature (left) and after it (right), MAPQ >= --mapq and any MAPQ
Verdict, first that applies: SPANNED (>= 2 spanning at MAPQ >= --mapq), SPANNED_LOW_MAPQ (>= 2 spanning at any MAPQ),
TELOMERE_STOP (no spanning read, >= 2 stopping reads whose clip is >= 50% telomeric: the chromosome ends at the feature),
ONE_SIDE_NO_READS (<= 2 reads of any MAPQ on one side, >= 10 on the other: no read supports that side),
NOT_SPANNED (no spanning read, >= 5 reads on both sides), LOW_COVERAGE.
With --merge, features of a flank closer than that are checked as one join (from the first start to the last end).
Samples run in parallel, the features of one sample in order (the BAM index is fetched once per sample)."""
import argparse, collections, csv, os, re, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--features", required=True)
ap.add_argument("--hifi-index", required=True, help="flagger_hifi_processing_metadata_v0.1.csv")
ap.add_argument("--workdir", required=True)
ap.add_argument("--window", type=int, default=30000)
ap.add_argument("--anchor", type=int, default=5000)
ap.add_argument("--near", type=int, default=3000)
ap.add_argument("--min-clip", type=int, default=2000)
ap.add_argument("--mapq", type=int, default=20)
ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--merge", type=int, default=0, help="merge features of a flank closer than this (bp) into one join")
ap.add_argument("--out", required=True)
a = ap.parse_args()
S3 = "https://s3-us-west-2.amazonaws.com/human-pangenomics"
QC = "submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC"
OP = re.compile(r"(\d+)([MIDNSHP=X])")
tel = lambda s: (len(re.findall("TTAGGG", s)) + len(re.findall("CCCTAA", s))) * 6 / max(len(s), 1)
log = lambda m: (sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"), sys.stderr.flush())
hifi = {r["sample_id"]: r["bam"].replace("s3://human-pangenomics", S3) for r in csv.DictReader(open(a.hifi_index))}


def ont_bam(sample):
    x = urllib.request.urlopen(f"{S3}?list-type=2&prefix={QC}/{sample}/hprc_r2/assembly_qc/read_alignments/ont/", timeout=60).read().decode()
    keys = [k for k in re.findall(r"<Key>([^<]+)</Key>", x) if k.endswith(".bam")]
    return [f"{S3}/{k}" for k in keys]


def check(bam, contig, s, e):
    r = subprocess.run(["samtools", "view", bam, f"{contig}:{max(1, s - a.window)}-{e + a.window}"],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"samtools view exit {r.returncode}: {r.stderr[-300:]}")
    span = span_any = stop = stop_tel = coverL = coverR = allL = allR = 0
    sa_hits = collections.Counter()
    for l in r.stdout.splitlines():
        f = l.split("\t")
        flag, pos, mq, cig, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
        if flag & 0x904:
            continue
        ops = [(int(k), o) for k, o in OP.findall(cig)]
        end = pos + sum(k for k, o in ops if o in "MDN=X") - 1
        lc = ops[0][0] if ops[0][1] == "S" else 0; rc = ops[-1][0] if ops[-1][1] == "S" else 0
        ovL = min(end, s - 1) - max(pos, s - 20000) + 1  # bp of the read on the 20 kb before the feature
        ovR = min(end, e + 20000) - max(pos, e + 1) + 1
        anc = min(a.anchor, s - 1)
        if pos < s - anc and end > e + a.anchor:
            span_any += 1; span += mq >= a.mapq
        if ovL >= 1000:
            allL += 1; coverL += mq >= a.mapq
        if ovR >= 1000:
            allR += 1; coverR += mq >= a.mapq
        clip = None
        if s - a.near <= end <= e and rc >= a.min_clip:  # ends at or inside the join (e.g. in its telomeric array)
            clip = seq[len(seq) - rc:]
        elif s <= pos <= e + a.near and lc >= a.min_clip:
            clip = seq[:lc]
        if clip:
            stop += 1; stop_tel += tel(clip) >= 0.5
            sa = next((t[5:] for t in f[11:] if t.startswith("SA:Z:")), "")
            for x in sa.split(";"):
                if x and int(x.split(",")[4]) >= a.mapq:
                    sa_hits[x.split(",")[0].split("#")[-1]] += 1
    if span >= 2:
        v = "SPANNED"
    elif span_any >= 2:
        v = "SPANNED_LOW_MAPQ"
    elif span_any == 0 and stop_tel >= 2:
        v = "TELOMERE_STOP"
    elif (allL <= 2 and allR >= 10) or (allR <= 2 and allL >= 10):
        v = "ONE_SIDE_NO_READS"
    elif span_any == 0 and allL >= 5 and allR >= 5:
        v = "NOT_SPANNED"
    else:
        v = "LOW_COVERAGE"
    return span, span_any, stop, stop_tel, coverL, coverR, allL, allR, ";".join(f"{k}:{n}" for k, n in sa_hits.most_common(3)) or ".", v


rows = list(csv.DictReader(open(a.features), delimiter="\t"))
if a.merge:  # one join per run of features closer than --merge: its side windows then hold no other feature
    g, n0 = collections.defaultdict(list), len(rows)
    for r in rows:
        g[r["flank"]].append(r)
    rows = []
    for fl, v in g.items():
        cur = None
        for r in sorted(v, key=lambda r: int(r["start"])):
            if cur and int(r["start"]) - int(cur["end"]) <= a.merge:
                cur["end"] = str(max(int(cur["end"]), int(r["end"]))); cur["kind"] += "+" + r["kind"]
                cur["strand_vs_end"] += "+" + r["strand_vs_end"]
                cur["dist_from_telomere"] = str(min(int(cur["dist_from_telomere"]), int(r["dist_from_telomere"])))
            else:
                rows.append(cur) if cur else None; cur = dict(r)
        rows.append(cur)
    log(f"{n0} features merged within {a.merge} bp into {len(rows)} joins")
by = collections.defaultdict(list)
for r in rows:
    by[r["sample"]].append(r)
os.makedirs(a.workdir, exist_ok=True)


def run_sample(sm):
    d = os.path.join(a.workdir, sm); os.makedirs(d, exist_ok=True); os.chdir(d) if False else None
    res = []
    for k in range(4):
        try:
            onts = ont_bam(sm); break
        except OSError as ex:
            if k == 3:
                raise SystemExit(f"FATAL: S3 listing for {sm}: {ex}")
            time.sleep(30 * (k + 1))
    bams = [("ONT", b) for b in onts] + ([("HiFi", hifi[sm])] if sm in hifi else [])
    if not bams:
        raise SystemExit(f"FATAL: no ONT or HiFi BAM for {sm}")
    for r in by[sm]:
        for tech, bam in bams[:1] + ([b for b in bams if b[0] == "HiFi"] if bams[0][0] == "ONT" else []):
            for k in range(4):
                try:
                    m = check(bam, r["contig"], int(r["start"]), int(r["end"])); break
                except (RuntimeError, subprocess.TimeoutExpired) as ex:
                    if k == 3:
                        raise SystemExit(f"FATAL: {sm} {r['contig']}:{r['start']}: {ex}")
                    log(f"retry {k + 1}: {sm} {r['contig']}:{r['start']}: {ex}"); time.sleep(30 * (k + 1))
            res.append((r, tech, os.path.basename(bam), m))
            log(f"{sm} {r['kind']} {r['contig']}:{r['start']}-{r['end']} {tech}: span {m[0]} (any MAPQ {m[1]}), stop {m[2]} (telomeric {m[3]}), cover MAPQ>={a.mapq} {m[4]}/{m[5]}, any {m[6]}/{m[7]}, clipped parts also on {m[8]} -> {m[9]}")
    return res


cwd = os.getcwd(); os.chdir(a.workdir)  # htslib saves remote BAM indexes in the working directory
with ThreadPoolExecutor(a.jobs) as ex:
    allres = [x for v in ex.map(run_sample, sorted(by)) for x in v]
os.chdir(cwd)
with open(a.out, "w") as fh:
    fh.write("flank\tsample\tcontig\tstart\tend\tkind\tstrand_vs_end\tdist_from_telomere\ttech\tbam\tspan\tspan_any_mapq\tstop\tstop_telomeric\tcover_left\tcover_right\tany_left\tany_right\tclip_sa\tverdict\n")
    for r, tech, bam, m in allres:
        fh.write("\t".join([r["flank"], r["sample"], r["contig"], r["start"], r["end"], r["kind"], r["strand_vs_end"],
                            r["dist_from_telomere"], tech, bam] + [str(x) for x in m]) + "\n")
c = collections.Counter((r["kind"], tech, m[9]) for r, tech, bam, m in allres)
for k in sorted(c):
    print(f"{k[0]}\t{k[1]}\t{k[2]}\t{c[k]}")
