#!/usr/bin/env python3
"""Whole-flank long-read check: for every base x of each flank (names contig:a-b_chrN_parm|qarm, 1-based), the number
of primary alignments that span it with >= --anchor bp on both sides (any MAPQ), from the HPRC release-2 ONT read
alignments to the same assembly (Flagger supplementary, over HTTPS). A BREAK is a run of bases, >= --anchor bp from the
flank edges, that no read spans: the two sides are never adjacent in a read, the signature of a false join. For each
break we give its distance from the telomere end, its length, and the reads covering the base --anchor bp before and
after it (telomere side / far side). Output per flank: reads, median coverage, minimum span, breaks, seconds, bytes."""
import argparse, collections, os, re, subprocess, sys, time, urllib.request
import numpy as np
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--flanks", required=True, help="one flank name per line")
ap.add_argument("--workdir", required=True)
ap.add_argument("--anchor", type=int, default=5000)
ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--out", required=True)
a = ap.parse_args()
S3 = "https://s3-us-west-2.amazonaws.com/human-pangenomics"
QC = "submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC"
NAME = re.compile(r"^(([^#]+)#[^#]+#[^:]+):(\d+)-(\d+)_(chr[0-9XY]+)_(parm|qarm)$")
OP = re.compile(r"(\d+)([MIDNSHP=X])")
log = lambda m: (sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"), sys.stderr.flush())


def ont_bam(sample):
    x = urllib.request.urlopen(f"{S3}?list-type=2&prefix={QC}/{sample}/hprc_r2/assembly_qc/read_alignments/ont/", timeout=60).read().decode()
    keys = [k for k in re.findall(r"<Key>([^<]+)</Key>", x) if k.endswith(".bam")]
    if len(keys) != 1:
        raise SystemExit(f"FATAL: {sample}: {len(keys)} ONT BAMs")
    return f"{S3}/{keys[0]}"


def scan(bam, nm):
    ctg, smp, fa, fb, arm = NAME.match(nm).group(1, 2, 3, 4, 6)
    fa, fb = int(fa), int(fb); L = fb - fa + 1; t0 = time.time()
    r = subprocess.run(["samtools", "view", bam, f"{ctg}:{max(1, fa - a.anchor)}-{fb + a.anchor}"], capture_output=True, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"samtools view exit {r.returncode}: {r.stderr[-300:].decode()}")
    span, cov, n = np.zeros(L + 1, np.int32), np.zeros(L + 1, np.int32), 0
    for l in r.stdout.decode().splitlines():
        f = l.split("\t", 6)
        if int(f[1]) & 0x904:
            continue
        n += 1; pos = int(f[3])
        end = pos + sum(int(k) for k, o in OP.findall(f[5]) if o in "MDN=X") - 1
        for arr, s, e in ((cov, pos, end), (span, pos + a.anchor, end - a.anchor)):
            s, e = max(s, fa) - fa, min(e, fb) - fa  # 0-based flank indexes, inclusive
            if s <= e:
                arr[s] += 1; arr[e + 1] -= 1
    span, cov = np.cumsum(span)[:L], np.cumsum(cov)[:L]
    brk = np.zeros(L, bool); brk[a.anchor:L - a.anchor] = span[a.anchor:L - a.anchor] == 0
    runs, x = [], 0
    while x < L:
        if brk[x]:
            y = x
            while y + 1 < L and brk[y + 1]:
                y += 1
            c0, c1 = cov[max(0, x - a.anchor)], cov[min(L - 1, y + a.anchor)]
            d, ct, cf = (x, c0, c1) if arm == "parm" else (L - 1 - y, c1, c0)  # bp from the telomere end to the near edge
            runs.append(f"{d}:{y - x + 1}:{ct}/{cf}"); x = y + 1
        else:
            x += 1
    inner = span[a.anchor:L - a.anchor]
    return [nm, smp, str(n), str(int(np.median(cov))), str(int(inner.min())), str(len(runs)), ",".join(runs) or ".",
            f"{time.time() - t0:.1f}", str(len(r.stdout))]


names = [l.strip() for l in open(a.flanks) if l.strip()]
by = collections.defaultdict(list)
for nm in names:
    if not NAME.match(nm):
        raise SystemExit(f"FATAL: cannot parse {nm}")
    by[NAME.match(nm).group(2)].append(nm)
log(f"{len(names)} flanks from {len(by)} samples, anchor {a.anchor}, jobs {a.jobs}")
a.out = os.path.abspath(a.out); os.makedirs(a.workdir, exist_ok=True); os.chdir(a.workdir)  # htslib saves remote BAM indexes here


def run_sample(sm):
    bam = ont_bam(sm); res = []
    for nm in by[sm]:
        for k in range(4):
            try:
                res.append(scan(bam, nm)); break
            except (RuntimeError, subprocess.TimeoutExpired) as ex:
                if k == 3:
                    raise SystemExit(f"FATAL: {nm}: {ex}")
                log(f"retry {k + 1}: {nm}: {ex}"); time.sleep(30 * (k + 1))
        log("\t".join(res[-1]))
    return res


with ThreadPoolExecutor(a.jobs) as ex:
    rows = [x for v in ex.map(run_sample, sorted(by)) for x in v]
with open(a.out, "w") as fh:
    fh.write("flank\tsample\treads\tmedian_cov\tmin_span\tbreaks\tbreaks(dist_from_telomere:len:cov_telomere_side/cov_far_side)\tseconds\tsam_bytes\n")
    fh.writelines("\t".join(r) + "\n" for r in rows)
log(f"done: {sum(r[5] != '0' for r in rows)} of {len(rows)} flanks with a break")
