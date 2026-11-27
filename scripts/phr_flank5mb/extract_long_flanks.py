#!/usr/bin/env python3
"""Long telomere-anchored flanks (e.g. 5 Mb) with the same contigs and telomere-trimmed starts as the 500 kb set,
extended inward to --flank bp (capped at the contig end): p arm a..a+F-1, q arm b-F+1..b, named
contig:start-end_chrN_arm (1-based). Reads from moosefs can fail now and then (EIO): every read of an assembly index
and every samtools faidx is retried (--retries, waiting --wait s x attempt), and each extracted sequence must have the
planned length, so a failed or short read never reaches the output. With --urls the assemblies are read over HTTPS
(e.g. the HPRC S3 bucket) instead of --fasta-dir: each assembly's .fai and .gzi are downloaded once to --idx-dir and
samtools faidx reads the remote bgzipped FASTA with them (set CURL_CA_BUNDLE if libcurl lacks CA certificates). With
--check-fasta every long flank must contain its 500 kb flank at the telomere end (p: prefix, q: suffix). Output: --out
(plain FASTA, plan order) and --names (long name, 500 kb name)."""
import argparse, glob, os, re, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--flank-fai", required=True, help="fai of the 500 kb flank FASTA")
src = ap.add_mutually_exclusive_group(required=True)
src.add_argument("--fasta-dir")
src.add_argument("--urls", help="file with one bgzipped FASTA URL per line")
ap.add_argument("--idx-dir", help="local copies of the remote .fai and .gzi (with --urls)")
ap.add_argument("--check-fasta", help="500 kb flank FASTA (bgzipped, indexed): each long flank must contain its flank")
ap.add_argument("--flank", type=int, required=True)
ap.add_argument("--arms", required=True, help="space-separated, e.g. 'chr13_parm chrX_qarm'")
ap.add_argument("--samtools", required=True)
ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--retries", type=int, default=6)
ap.add_argument("--wait", type=int, default=30)
ap.add_argument("--out", required=True)
ap.add_argument("--names", required=True)
a = ap.parse_args()
log = lambda m: (sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n"), sys.stderr.flush())  # python 3.7 buffers a file
retried = 0


def retry(what, fn):
    global retried
    for k in range(1, a.retries + 1):
        try:
            return fn()
        except (OSError, RuntimeError) as e:
            if k == a.retries:
                raise SystemExit(f"FATAL: {what} failed {k} times: {e}")
            retried += 1
            log(f"retry {k}/{a.retries - 1} in {a.wait * k} s: {what}: {e}")
            time.sleep(a.wait * k)


def read_lines(p):
    with open(p) as fh:
        return fh.read().splitlines()


contig, idx = {}, {}
if a.urls:
    if not a.idx_dir:
        raise SystemExit("FATAL: --urls needs --idx-dir")
    os.makedirs(a.idx_dir, exist_ok=True)
    fastas = [u.strip() for u in open(a.urls) if u.strip()]

    def fetch(u, dst):
        def once():
            urllib.request.urlretrieve(u, dst + ".part")
            if os.path.getsize(dst + ".part") == 0:
                raise RuntimeError("empty file")
            os.replace(dst + ".part", dst)
        if not os.path.exists(dst):
            retry(f"download {u}", once)
    for f in fastas:
        idx[f] = tuple(os.path.join(a.idx_dir, os.path.basename(f) + e) for e in (".fai", ".gzi"))
    with ThreadPoolExecutor(a.jobs) as ex:  # list() re-raises a failed download here
        list(ex.map(lambda x: fetch(*x), [(f + e, idx[f][k]) for f in fastas for k, e in enumerate((".fai", ".gzi"))]))
    log(f"indexes of {len(fastas)} remote FASTAs in {a.idx_dir}")
else:
    fastas = retry(f"list {a.fasta_dir}", lambda: sorted(glob.glob(os.path.join(a.fasta_dir, "*.fa.gz"))))
for f in fastas:
    fai = idx[f][0] if a.urls else f + ".fai"
    if not a.urls and not os.path.exists(fai):
        continue
    for l in retry(f"read {fai}", lambda: read_lines(fai)):
        n, L = l.split("\t")[:2]
        contig[n] = (f, int(L))
log(f"contigs indexed: {len(contig)} from {len(fastas)} FASTA files")
arms = set(a.arms.split())
plan, capped = [], 0
for l in retry(f"read {a.flank_fai}", lambda: read_lines(a.flank_fai)):
    name = l.split("\t")[0]
    m = re.match(r"^(.+):(\d+)-(\d+)_(chr\w+?)_(parm|qarm)$", name)
    if not m:
        raise SystemExit(f"FATAL: cannot parse {name}")
    ctg, s, e, ch, arm = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5)
    if f"{ch}_{arm}" not in arms:
        continue
    if ctg not in contig:
        raise SystemExit(f"FATAL: contig {ctg} in no FASTA index")
    f, L = contig[ctg]
    na, nb = (s, min(s + a.flank - 1, L)) if arm == "parm" else (max(e - a.flank + 1, 1), e)
    capped += nb - na + 1 < a.flank
    plan.append((f, f"{ctg}:{na}-{nb}", f"{ctg}:{na}-{nb}_{ch}_{arm}", name, nb - na + 1))
log(f"plan: {len(plan)} flanks, {capped} shorter than {a.flank} bp (contig end)")


def extract(row):
    f, reg, new, old, want = row
    def once():
        cmd = [a.samtools, "faidx"] + (["--fai-idx", idx[f][0], "--gzi-idx", idx[f][1]] if a.urls else []) + [f, reg]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"exit {r.returncode}: {r.stderr.strip()[-300:]}")
        seq = "".join(r.stdout.split("\n")[1:])
        if len(seq) != want:
            raise RuntimeError(f"got {len(seq)} bp, planned {want}")
        return seq
    seq = retry(f"samtools faidx {f} {reg}", once)
    if a.check_fasta:
        r = subprocess.run([a.samtools, "faidx", a.check_fasta, old], capture_output=True, text=True)
        short = "".join(r.stdout.split("\n")[1:])
        if r.returncode != 0 or not short:
            raise SystemExit(f"FATAL: cannot read {old} from {a.check_fasta}: {r.stderr.strip()[-300:]}")
        if (seq[:len(short)] if old.endswith("_parm") else seq[-len(short):]) != short:
            raise SystemExit(f"FATAL: {new} does not contain {old} at its telomere end")
    return new, seq


bp = 0
with open(a.out, "w") as out, open(a.names, "w") as nm, ThreadPoolExecutor(a.jobs) as ex:
    for i in range(0, len(plan), 4 * a.jobs):  # chunks keep memory bounded and the output in plan order
        chunk = plan[i:i + 4 * a.jobs]
        for (new, seq), row in zip(ex.map(extract, chunk), chunk):
            out.write(f">{new}\n{seq}\n"); nm.write(f"{new}\t{row[3]}\n"); bp += len(seq)
        if (i // (4 * a.jobs)) % 20 == 0:
            log(f"extracted {min(i + 4 * a.jobs, len(plan))} of {len(plan)} flanks, {bp} bp, retries so far {retried}")
log(f"done: {len(plan)} flanks, {bp} bp, retries {retried}" + (f"; all contain their 500 kb flank ({a.check_fasta})" if a.check_fasta else ""))
