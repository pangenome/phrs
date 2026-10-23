#!/usr/bin/env python3
"""Per-family test of the homogenisation signature requested by Referee #1 (major point 1).

Referee #1: "Concerted evolution is never tested. It requires that paralogous copies at
different ends be more similar than their divergence times predict, i.e. actively
homogenized; high identity alone does not show this, since recent duplication gives the
same signature. The data allow a direct test: ... a per-family comparison of
within-individual divergence between non-homologous copies against between-individual
divergence at the same end ..."

This script implements the *per-family version* of the allele-vs-paralog test that the
manuscript already reports in aggregate (end-to-end-report/report/04_heterogeneity.md,
section "Allele vs paralog distance"):

  allelic distance  = distance between the two haplotypes of ONE individual at ONE arm
                      (maternal vs paternal copy of the same arm)
  paralog distance  = distance from that sequence to the CLOSEST sequence on a DIFFERENT
                      arm of the SAME community, within the same individual

  homogenisation signature: paralog distance < allelic distance.

Families (names are the ones Referee #1 lists): D4Z4/DUX4L, TUBB8/TUBB8B, OR4F, SEPTIN14.

Pipeline
--------
1. Locate every annotated copy of each family in CHM13 from the RefSeq Liftoff GFF3
   (data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz), assign each copy to a chromosome arm
   using the centromere intervals in data/chm13-annotations.bed, and keep the
   subtelomeric copies (midpoint within 1 Mb of that arm's distal end).
2. Arms -> communities from data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv.
3. Two sequence substrates, one window per haplotype per arm:
     P = the telomere-anchored PHR windows
         (/moosefs/.../PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz) -- the exact
         substrate of the aggregate analysis in report section 04, i.e. the windows that
         carry the duplicons.
     F = the 100 kb unique-sequence FLANKS that abut those PHR windows on the centromere
         side (/moosefs/.../PHR_III/flanking/hprcv2.1Mb.flanking_100kb.fa.gz) -- the
         flanking control substrate; it is adjacent to, not inside, the PHR.
4. Distance = Mash distance (k=21, sketch 100000). Every parameter above is echoed into
   the output so the number is traceable.
5. Per (family, substrate, arm pair in the same community) report n, medians, the
   fraction of matched units where the paralog is closer than the alleles, an exact
   binomial sign test and a one-sided Wilcoxon signed-rank test.

Everything reported is computed and printed; nothing is estimated or assumed.

Usage
-----
  concerted_evolution_families.py                 # full run, writes the TSV to stdout
  concerted_evolution_families.py --self-test     # oracle + adversarial unit tests only
  concerted_evolution_families.py --workscratch /scratch/erikg/celfam
"""

import argparse
import gzip
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

GFF3 = os.path.join(REPO, "data", "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
CHROM_SIZES = os.path.join(REPO, "data", "chm13.chrom.sizes")
ANNOT_BED = os.path.join(REPO, "data", "chm13-annotations.bed")
COMMUNITIES = os.path.join(
    REPO, "data", "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")

PHR_FA = ("/moosefs/guarracino/HPRCv2/PHR_III/"
          "hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz")
FLANK_FA = ("/moosefs/guarracino/HPRCv2/PHR_III/flanking/"
            "hprcv2.1Mb.flanking_100kb.fa.gz")

SUBSTRATES = {"P": PHR_FA, "F": FLANK_FA}

# A copy counts as subtelomeric when its midpoint is within this distance of the arm tip.
SUBTELOMERIC_MAX_DIST = 1000000

MASH_K = 21
MASH_SKETCH = 20000

# Family definitions: (family label, gene-name regex, human-readable scope note).
FAMILIES = [
    ("D4Z4_DUX4L",
     re.compile(r"^(DUX4|DUX4L[0-9]+)$"),
     "D4Z4 macrosatellite unit genes (DUX4 + DUX4L copies)"),
    ("TUBB8_TUBB8B",
     re.compile(r"^(TUBB8|TUBB8B)$"),
     "TUBB8 and TUBB8B"),
    ("OR4F",
     re.compile(r"^OR4F[0-9]*[A-Z]*$"),
     "OR4F gene family incl. OR4FnP pseudogene copies"),
    ("SEPTIN14",
     re.compile(r"^(SEPTIN14|SEPTIN14P[0-9]+)$"),
     "SEPTIN14 and SEPTIN14P copies"),
]

SAMPLE_NAME_RE = re.compile(r"^([^#]+)#([0-9]+)#(.*)$")


def die(msg):
    sys.stderr.write("FATAL: %s\n" % msg)
    sys.exit(1)


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------

def run(cmd, timeout=3600, stdout=subprocess.PIPE):
    """Run a command, fail loud on non-zero status or timeout."""
    try:
        res = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE,
                             timeout=timeout)
    except subprocess.TimeoutExpired:
        die("timeout (>%ss): %s" % (timeout, " ".join(cmd)))
    if res.returncode != 0:
        die("command failed (%d): %s\n%s" % (res.returncode, " ".join(cmd),
                                             res.stderr.decode("utf8", "replace")[-2000:]))
    return res


def fai_records(fa):
    """Read a .fai next to a (bgzip) FASTA. Returns list of (name, length)."""
    fai = fa + ".fai"
    if not os.path.exists(fai):
        die("missing .fai for %s" % fa)
    out = []
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1])))
    return out


def arm_of_window_name(name):
    """Extract chrN_parm / chrN_qarm from either substrate's window name."""
    m = re.search(r"(chr[0-9XY]+_[pq]arm)", name)
    if not m:
        die("cannot parse arm from window name: %r" % name)
    return m.group(1)


def sample_hap_of_window_name(name):
    m = SAMPLE_NAME_RE.match(name)
    if not m:
        die("cannot parse sample/hap from window name: %r" % name)
    return m.group(1), int(m.group(2))


def arm_to_community_label(arm):
    """'chr4_qarm' -> 'chr4_q' (the label used by the community assignment table)."""
    return arm.replace("_parm", "_p").replace("_qarm", "_q")


def window_key(name):
    """(sample, hap, arm) for a window sequence name."""
    s, h = sample_hap_of_window_name(name)
    return (s, h, arm_of_window_name(name))


# --------------------------------------------------------------------------------------
# Step 1: family copies in CHM13
# --------------------------------------------------------------------------------------

def load_sizes():
    sizes = {}
    with open(CHROM_SIZES) as fh:
        for line in fh:
            c, s = line.split()
            sizes[c] = int(s)
    return sizes


def load_centromeres():
    cent = {}
    with open(ANNOT_BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == "Centromere":
                cent[f[0]] = (int(f[1]), int(f[2]))
    return cent


def load_communities():
    """arm label ('chr4_q') -> community id ('C1')."""
    out = {}
    with open(COMMUNITIES) as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            out[f[0]] = f[1]
    return out


def build_family_copies():
    """Return list of dicts, one per annotated copy of any target family."""
    sizes = load_sizes()
    cent = load_centromeres()
    copies = []
    with gzip.open(GFF3, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            gene = attrs.get("gene", "")
            fam = None
            for label, rx, _note in FAMILIES:
                if rx.match(gene):
                    fam = label
                    break
            if fam is None:
                continue
            chrom = f[0]
            if chrom not in cent or chrom not in sizes:
                continue
            start, end = int(f[3]), int(f[4])
            cs, ce = cent[chrom]
            L = sizes[chrom]
            if end <= cs:
                arm = chrom + "_parm"
                tip_dist = start
            elif start >= ce:
                arm = chrom + "_qarm"
                tip_dist = L - end
            else:
                arm = chrom + "_c"
                tip_dist = min(start, L - end)
            mid = (start + end) // 2
            if arm.endswith("_parm"):
                sub = mid < SUBTELOMERIC_MAX_DIST
            elif arm.endswith("_qarm"):
                sub = (L - mid) < SUBTELOMERIC_MAX_DIST
            else:
                sub = False
            # distance from the copy midpoint to the distal (telomeric) end of its arm
            tip_dist = mid if arm.endswith("_parm") else (L - mid)
            copies.append(dict(family=fam, chrom=chrom, arm=arm, start=start, end=end,
                               mid=mid, strand=f[6], gene=gene,
                               tip_dist=tip_dist, subtelomeric=sub))
    copies.sort(key=lambda r: (r["family"], r["arm"], r["start"]))
    return copies


def family_arms(copies):
    """family -> sorted list of arms carrying >=1 subtelomeric copy."""
    out = {}
    for c in copies:
        if c["subtelomeric"]:
            out.setdefault(c["family"], set()).add(c["arm"])
    return {k: sorted(v) for k, v in out.items()}


# --------------------------------------------------------------------------------------
# Step 2/3: windows from a substrate
# --------------------------------------------------------------------------------------

def extract_windows(fa, wanted_arms, workdir, tag, timeout=600):
    """Write a FASTA of all windows on `wanted_arms` (uppercased, headers unchanged).

    Returns (path, {header: length}). Fails loud if the substrate is unreadable.
    """
    recs = fai_records(fa)
    keep = [(n, ln) for (n, ln) in recs if arm_of_window_name(n) in wanted_arms]
    if not keep:
        die("no windows found on %s in %s" % (sorted(wanted_arms), fa))

    regions = os.path.join(workdir, "regions_%s.txt" % tag)
    with open(regions, "w") as fh:
        for n, _ln in keep:
            fh.write(n + "\n")

    raw = os.path.join(workdir, "raw_%s.fa" % tag)
    res = run(["samtools", "faidx", "-r", regions, fa], timeout=timeout)
    text = res.stdout.decode("ascii", "replace")
    if "returning empty sequence" in res.stderr.decode("utf8", "replace"):
        die("substrate %s returned empty sequence(s); treating as unreadable" % fa)

    lengths = {}
    with open(raw, "w") as out:
        hdr = None
        buf = []
        for line in text.splitlines():
            if line.startswith(">"):
                if hdr is not None:
                    seq = "".join(buf)
                    out.write(">%s\n%s\n" % (hdr, seq))
                    lengths[hdr] = len(seq)
                    WINDOW_NAMES[window_key(hdr)] = hdr
                hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line.upper())
        if hdr is not None:
            seq = "".join(buf)
            out.write(">%s\n%s\n" % (hdr, seq))
            lengths[hdr] = len(seq)
            WINDOW_NAMES[window_key(hdr)] = hdr
    if not lengths:
        die("substrate %s produced no sequences" % fa)
    return raw, lengths


def window_span_of_chm13(fa, arm):
    """CHM13 window name and (start,end) on the chromosome for one arm, or None."""
    for n, _ln in fai_records(fa):
        if n.startswith("CHM13#0#") and arm_of_window_name(n) == arm:
            m = re.search(r":([0-9]+)-([0-9]+)_chr", n)
            if not m:
                return None
            return n, int(m.group(1)), int(m.group(2))
    return None


# --------------------------------------------------------------------------------------
# Step 4: mash distances
# --------------------------------------------------------------------------------------

def read_fasta(path):
    """Minimal FASTA reader -> list of (header, sequence)."""
    out = []
    hdr = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(buf)))
                hdr = line[1:]
                buf = []
            else:
                buf.append(line)
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    return out


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for h, sq in recs:
            fh.write(">%s\n%s\n" % (h, sq))


def mash_sketch(fa, out_prefix):
    run(["mash", "sketch", "-i", "-k", str(MASH_K), "-s", str(MASH_SKETCH),
         "-o", out_prefix, fa], timeout=3600)


def mash_dist_pairwise(msh_path):
    """All pairwise distances inside one sketch group -> dict[(name1,name2)] = rec."""
    res = run(["mash", "dist", msh_path, msh_path], timeout=3600)
    out = {}
    for line in res.stdout.decode("ascii", "replace").splitlines():
        f = line.split("\t")
        if len(f) != 5:
            die("unexpected mash dist line (%d fields): %r" % (len(f), line))
        if f[0] == f[1]:
            continue
        sh, tot = f[4].split("/")
        out[(f[0], f[1])] = (float(f[2]), int(sh), int(tot))
    return out


def same_sample_distances(raw_fa, workdir, tag):
    """Distances between the windows of each individual, computed one sample at a time.

    Mash's all-vs-all is quadratic in the number of windows, so sketching and comparing
    per individual keeps the work proportional to the number of pairs actually needed.
    """
    recs = read_fasta(raw_fa)
    by_sample = {}
    for h, sq in recs:
        by_sample.setdefault(sample_hap_of_window_name(h)[0], []).append((h, sq))
    same = {}
    for smp in sorted(by_sample):
        fa = os.path.join(workdir, "smp_%s_%s.fa" % (tag, smp))
        write_fasta(fa, by_sample[smp])
        sk = os.path.join(workdir, "sk_%s_%s" % (tag, smp))
        mash_sketch(fa, sk)
        for (a, b), rec in mash_dist_pairwise(sk + ".msh").items():
            k1, k2 = window_key(a), window_key(b)
            key = (k1, k2) if k1 <= k2 else (k2, k1)
            same[key] = rec
        os.unlink(fa)
    sys.stderr.write("  within-individual windows: %d sequences, %d pairs\n"
                     % (len(recs), len(same)))
    return same


def between_individual_distances(raw_fa, workdir, tag, arms, max_haps=120, seed=20260501):
    """Reference scale: distances between DIFFERENT individuals at the same arm."""
    import random as _random
    rng = _random.Random(seed)
    recs = read_fasta(raw_fa)
    per_arm = {}
    for h, sq in recs:
        per_arm.setdefault(arm_of_window_name(h), []).append((h, sq))
    out = {}
    for arm in arms:
        lst = per_arm.get(arm, [])
        if len(lst) < 3:
            continue
        if len(lst) > max_haps:
            lst = rng.sample(lst, max_haps)
        fa = os.path.join(workdir, "ba_%s_%s.fa" % (tag, arm))
        write_fasta(fa, lst)
        sk = os.path.join(workdir, "bask_%s_%s" % (tag, arm))
        mash_sketch(fa, sk)
        vals = []
        for (a, b), rec in mash_dist_pairwise(sk + ".msh").items():
            if sample_hap_of_window_name(a)[0] == sample_hap_of_window_name(b)[0]:
                continue
            vals.append(rec)
        if vals:
            out[arm] = vals
        os.unlink(fa)
    sys.stderr.write("  between-individual reference arms: %d\n" % len(out))
    return out


def jaccard_from_shared(shared, total):
    """Jaccard estimate from a Mash dist line.

    Mash prints "shared/total" where total is the size of the union of the two
    sketches (verified against an exact canonical-k-mer oracle in --self-test), so
    the Jaccard estimate is simply shared/total.
    """
    if total <= 0:
        return float("nan")
    if shared > total:
        die("mash reported more shared hashes than union size")
    return shared / float(total)


# --------------------------------------------------------------------------------------
# Step 5: allelic / paralog classification
# --------------------------------------------------------------------------------------

WINDOW_NAMES = {}


def window_name(sample, hap, arm):
    """FASTA header of the window for (sample, hap, arm) in the current substrate."""
    try:
        return WINDOW_NAMES[(sample, hap, arm)]
    except KeyError:
        die("unknown window (sample=%s hap=%s arm=%s)" % (sample, hap, arm))


def window_distance(same_sample, s, h1, a1, h2, a2):
    k1, k2 = (s, h1, a1), (s, h2, a2)
    key = (k1, k2) if k1 <= k2 else (k2, k1)
    return same_sample.get(key)


def classify(same_sample, arms, communities, haps=(1, 2)):
    """Build allelic and paralog records from the same-individual distance table.

    allelic[(sample, arm)]       = (distance, shared, total)  hap1 vs hap2 at that arm
    paralog[(sample, armA, armB)] = (distance, shared, total)  closest cross-arm copy,
                                   minimum over the 2x2 haplotype combinations, and only
                                   for arm pairs that share a community.
    """
    idx = {}
    for (k1, k2) in same_sample:
        idx[k1[0] + "|" + k1[2]] = True
        idx[k2[0] + "|" + k2[2]] = True

    def get(s, h1, a1, h2, a2):
        k1, k2 = (s, h1, a1), (s, h2, a2)
        key = (k1, k2) if k1 <= k2 else (k2, k1)
        return same_sample.get(key)

    samples = sorted({s for s in (k.split("|")[0] for k in idx)})
    allelic = {}
    paralog = {}
    for s in samples:
        for a in arms:
            if not any(("%s|%s" % (s, a)) == k for k in idx):
                continue
            d = get(s, 1, a, 2, a)
            if d is not None:
                allelic[(s, a)] = d
        for i in range(len(arms)):
            for j in range(i + 1, len(arms)):
                a, b = arms[i], arms[j]
                ca = communities.get(arm_to_community_label(a))
                cb = communities.get(arm_to_community_label(b))
                if ca is None or cb is None or ca != cb:
                    continue
                best = None
                for ha in haps:
                    for hb in haps:
                        d = get(s, ha, a, hb, b)
                        if d is not None and (best is None or d[0] < best[0]):
                            best = d
                if best is not None:
                    paralog[(s, a, b)] = best
    return allelic, paralog


# --------------------------------------------------------------------------------------
# Step 6: statistics
# --------------------------------------------------------------------------------------

def sign_test(n_better, n_worse):
    """Exact two-sided binomial sign test on discordant pairs."""
    from scipy.stats import binomtest
    n = n_better + n_worse
    if n == 0:
        return float("nan")
    return binomtest(n_better, n, 0.5, alternative="two-sided").pvalue


def wilcoxon_less(x, y):
    """One-sided Wilcoxon signed-rank: H1 = x tends to be smaller than y."""
    from scipy.stats import wilcoxon
    if len(x) < 10:
        return float("nan")
    try:
        return wilcoxon(x, y, alternative="less", zero_method="wilcox").pvalue
    except ValueError:
        return float("nan")


def run_self_test():
    """Brute-force oracle at tiny scale plus adversarial tests of the pairing logic."""
    print("== self-test 1: mash distance vs exact Jaccard on tiny synthetic sequences")
    rng = random.Random(7)
    alphabet = "ACGT"
    base = "".join(rng.choice(alphabet) for _ in range(5000))
    # derived sequences with known edit content
    mut = list(base)
    for pos in range(0, 5000, 200):        # 25 substitutions = 0.5% divergence
        mut[pos] = rng.choice([c for c in alphabet if c != mut[pos]])
    seqs = {"s_ref": base,
            "s_ident": base,
            "s_mut5": "".join(mut),
            "s_unrel": "".join(rng.choice(alphabet) for _ in range(5000))}
    tmp = tempfile.mkdtemp(prefix="celfam_selftest_")
    try:
        fa = os.path.join(tmp, "t.fa")
        with open(fa, "w") as fh:
            for k, v in seqs.items():
                fh.write(">%s\n%s\n" % (k, v))
        msh = os.path.join(tmp, "sst")
        run(["mash", "sketch", "-i", "-k", str(MASH_K), "-s", str(MASH_SKETCH),
             "-o", msh, fa])
        res = run(["mash", "dist", msh + ".msh", msh + ".msh"])
        d = {}
        for line in res.stdout.decode("ascii", "replace").splitlines():
            f = line.split("\t")
            assert len(f) == 5, line
            if f[0] == f[1]:
                continue
            sh, tot = f[4].split("/")
            d[(f[0], f[1])] = (float(f[2]), int(sh), int(tot))

        comp = str.maketrans("ACGT", "TGCA")

        def canon(kmer):
            rc = kmer.translate(comp)[::-1]
            return kmer if kmer <= rc else rc

        def exact_jaccard(x, y, k=21):
            # mash canonicalises k-mers by default, so the oracle must too
            kx = {canon(x[i:i + k]) for i in range(len(x) - k + 1)}
            ky = {canon(y[i:i + k]) for i in range(len(y) - k + 1)}
            return len(kx & ky) / len(kx | ky)

        for (a, b) in (("s_ref", "s_ident"), ("s_ref", "s_mut5"), ("s_ref", "s_unrel")):
            dd = d[(a, b)]
            ej = exact_jaccard(seqs[a], seqs[b])
            mj = jaccard_from_shared(dd[1], dd[2])
            print("   %-8s vs %-8s mash_d=%.6f  exact_J=%.5f  mash_J=%.5f"
                  % (a, b, dd[0], ej, mj))
            assert abs(ej - mj) < 0.06, ("mash Jaccard disagrees with exact Jaccard", a, b, ej, mj)
        assert d[("s_ref", "s_ident")][0] == 0.0
        assert d[("s_ref", "s_ident")][1] == d[("s_ref", "s_ident")][2]

        print("== self-test 2: allelic/paralog classification on a synthetic panel")
        # Two arms in community C1, one arm in C9 (must not form paralog pairs).
        comm = {"chr4_q": "C1", "chr10_q": "C1", "chr1_p": "C9"}
        same = {}

        def put(s, h1, a1, h2, a2, val):
            k1, k2 = (s, h1, a1), (s, h2, a2)
            key = (k1, k2) if k1 <= k2 else (k2, k1)
            same[key] = (val, 10, 100)

        # S1: paralog much closer than alleles -> counts as paralog closer
        put("S1", 1, "chr4_qarm", 2, "chr4_qarm", 0.10)     # allelic 4q
        put("S1", 1, "chr10_qarm", 2, "chr10_qarm", 0.10)   # allelic 10q
        put("S1", 1, "chr4_qarm", 2, "chr10_qarm", 0.01)    # closest cross-arm
        # S2: paralog farther than alleles -> not counted
        put("S2", 1, "chr4_qarm", 2, "chr4_qarm", 0.01)
        put("S2", 1, "chr10_qarm", 2, "chr10_qarm", 0.01)
        put("S2", 1, "chr4_qarm", 1, "chr10_qarm", 0.50)
        # S3: only one haplotype at 10q -> no allelic for 10q, but paralog still defined
        put("S3", 1, "chr4_qarm", 2, "chr4_qarm", 0.05)
        put("S3", 1, "chr4_qarm", 1, "chr10_qarm", 0.20)
        # chr1_p (C9) share must never enter the paralog dict
        put("S1", 1, "chr4_qarm", 1, "chr1_parm", 0.0001)

        allelic, paralog = classify(same, ["chr4_qarm", "chr10_qarm", "chr1_parm"], comm)
        assert ("S3", "chr10_qarm") not in allelic, "single-haplotype arm leaked an allelic pair"
        assert ("S1", "chr4_qarm") in allelic and allelic[("S1", "chr4_qarm")][0] == 0.10
        assert ("S1", "chr4_qarm", "chr10_qarm") in paralog
        assert paralog[("S1", "chr4_qarm", "chr10_qarm")][0] == 0.01, \
            "closest-cross-arm selection failed (must take min over 2x2 haplotypes)"
        assert ("S1", "chr1_parm", "chr4_qarm") not in paralog and \
               ("S1", "chr4_qarm", "chr1_parm") not in paralog, \
            "cross-community pair leaked into paralog set"
        print("   allelic=%d paralog=%d (expected 5 and 3)" % (len(allelic), len(paralog)))
        assert len(allelic) == 5, sorted(allelic)
        assert len(paralog) == 3, sorted(paralog)

        print("== self-test 3: sign test and Wilcoxon on a known-direction toy vector")
        p_sign = sign_test(10, 0)
        assert p_sign < 0.01, p_sign
        p_w = wilcoxon_less([0.1] * 12, [0.5] * 12)
        assert p_w < 1e-3, p_w
        print("   sign p=%.3g  wilcoxon p=%.3g" % (p_sign, p_w))
        print("ALL SELF-TESTS PASSED")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------------------
# Alignment cross-check: does the Mash (k-mer) distance track base-level alignment identity?
# --------------------------------------------------------------------------------------

def alignment_identity(fa, ida, idb, workdir, tag):
    """Base-level identity of the best alignment between two windows, via wfmash.

    Returns (identity, aligned_bp) where identity = matches / alignment length over all
    alignment blocks of the selected mapping, computed from the =/X CIGAR.
    """
    tmp = os.path.join(workdir, "ac_%s.fa" % tag)
    run(["samtools", "faidx", fa, ida, idb], stdout=open(tmp, "wb"))
    # wfmash needs a .fai next to the query file
    run(["samtools", "faidx", tmp])
    map_paf = os.path.join(workdir, "ac_%s.map.paf" % tag)
    with open(map_paf, "wb") as fh:
        run(["wfmash", "-m", "-s", "1000", "-w", "1000", "-p", "60", "-n", "20",
             tmp, tmp], stdout=fh, timeout=600)
    aln_paf = os.path.join(workdir, "ac_%s.aln.paf" % tag)
    with open(aln_paf, "wb") as fh:
        run(["wfmash", "-s", "1000", "-w", "10000", "-p", "60",
             "-i", map_paf, tmp, tmp], stdout=fh, timeout=600)
    best = None
    with open(aln_paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == f[5]:
                continue
            cg = None
            for tagf in f[12:]:
                if tagf.startswith("cg:Z:"):
                    cg = tagf[5:]
            if cg is None:
                continue
            m = sum(int(n) for n, op in re.findall(r"(\d+)([=XID])", cg) if op == "=")
            tot = sum(int(n) for n, op in re.findall(r"(\d+)([=XID])", cg))
            if tot == 0:
                continue
            if best is None or tot > best[1]:
                best = (m / tot, tot)
    if best is None:
        return None
    return best


def alignment_crosscheck(pool, n, workdir, seed=20260501):
    from scipy.stats import spearmanr
    rng = random.Random(seed)
    sample = rng.sample(pool, min(n, len(pool)))
    rows = []
    for i, (fa, ida, idb, mdist, mjac) in enumerate(sample):
        res = alignment_identity(fa, ida, idb, workdir, "%d" % i)
        if res is None:
            sys.stderr.write("[align-check] no alignment for %s vs %s\n" % (ida, idb))
            continue
        rows.append((mdist, mjac, res[0], res[1]))
    if len(rows) < 5:
        sys.stderr.write("[align-check] only %d usable pairs; skipping correlation\n"
                         % len(rows))
        return
    md = [r[0] for r in rows]
    ai = [r[2] for r in rows]
    rho = spearmanr(md, ai)
    out = os.path.join(workdir, "align_check.tsv")
    with open(out, "w") as fh:
        fh.write("mash_dist\tmash_jaccard\taln_identity\taln_bp\n")
        for r in rows:
            fh.write("%.6f\t%.6f\t%.6f\t%d\n" % r)
    sys.stderr.write("[align-check] n=%d  spearman(mash_dist, aln_identity)=%.4f p=%.3g\n"
                     "  median aln_identity=%.4f  median mash_jaccard=%.4f  -> %s\n"
                     % (len(rows), rho.statistic, rho.pvalue,
                        sorted(ai)[len(ai) // 2], sorted(r[1] for r in rows)[len(rows) // 2],
                        out))
    assert rho.statistic < 0, ("Mash distance should anti-correlate with alignment "
                               "identity", rho.statistic)


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--workscratch", default=None,
                    help="scratch directory for extracted windows (default: mktemp)")
    ap.add_argument("--substrates", default="P,F",
                    help="comma separated subset of substrates to run (P,F)")
    ap.add_argument("--families", default=",".join(f[0] for f in FAMILIES))
    ap.add_argument("--align-check", type=int, default=0,
                    help="N random real pairs to cross-check against wfmash alignment identity")
    args = ap.parse_args()

    if args.self_test:
        run_self_test()
        return

    fams = [f for f in FAMILIES if f[0] in args.families.split(",")]
    subs = [s for s in args.substrates.split(",") if s]

    workdir = args.workscratch or tempfile.mkdtemp(prefix="celfam_")
    os.makedirs(workdir, exist_ok=True)

    copies = build_family_copies()
    arms_by_family = family_arms(copies)
    communities = load_communities()

    meta = []
    for label, _rx, note in fams:
        arms = arms_by_family.get(label, [])
        meta.append((label, note, arms))
        sys.stderr.write("[family] %-14s arms=%s\n" % (label, ",".join(arms) or "-"))

    # --- header of the deliverable table ----------------------------------------------
    cols = [
        "family", "substrate", "row_type", "community", "arm_a", "arm_b",
        "n_allelic_pairs", "median_allelic_dist", "median_allelic_jaccard",
        "n_paralog_pairs", "median_paralog_dist", "median_paralog_jaccard",
        "n_between_pairs", "median_between_individual_dist",
        "n_matched_units", "n_paralog_closer", "n_allelic_closer",
        "frac_paralog_closer", "sign_test_p", "wilcoxon_p",
        "chm13_copies_in_arm", "chm13_copies_in_chm13_window",
        "chm13_window_bp", "mash_k", "mash_sketch",
    ]
    print("\t".join(cols))

    import statistics

    def emit(**kw):
        print("\t".join(str(kw.get(c, "")) for c in cols))

    align_pool = []   # (dist, jaccard, aln_identity) collected for the cross-check

    for label, note, arms in meta:
        if not arms:
            emit(family=label, substrate="-", row_type="NOT_RUN",
                 chm13_copies_in_arm=0)
            continue
        for sub in subs:
            fa = SUBSTRATES[sub]
            if not os.path.exists(fa):
                sys.stderr.write("[WARN] missing substrate %s for %s\n" % (sub, label))
                emit(family=label, substrate=sub, row_type="NOT_RUN_MISSING_SUBSTRATE")
                continue
            tag = "%s_%s" % (label, sub)
            raw, lengths = extract_windows(fa, set(arms), workdir, tag)
            sys.stderr.write("[%s %s] windows=%d bp=%d\n"
                             % (label, sub, len(lengths), sum(lengths.values())))
            same_sample = same_sample_distances(raw, workdir, tag)
            between = between_individual_distances(raw, workdir, tag, arms)
            allelic, paralog = classify(same_sample, arms, communities)

            # per-arm CHM13 copy accounting (how much of the family the window holds)
            copies_in_arm = {}
            copies_in_win = {}
            win_bp = {}
            for a in arms:
                copies_in_arm[a] = sum(1 for c in copies
                                       if c["family"] == label and c["arm"] == a)
                sp = window_span_of_chm13(fa, a)
                if sp is None:
                    copies_in_win[a] = 0
                    win_bp[a] = 0
                    continue
                _n, ws, we = sp
                win_bp[a] = we - ws
                copies_in_win[a] = sum(
                    1 for c in copies
                    if c["family"] == label and c["arm"] == a
                    and c["start"] >= ws and c["end"] <= we)

            # matched units: (sample, armA, armB) with both allelic values available
            per_arm_pair = {}
            for (s, a, b), d in paralog.items():
                ka, kb = (s, a), (s, b)
                al = []
                if ka in allelic:
                    al.append(allelic[ka][0])
                if kb in allelic:
                    al.append(allelic[kb][0])
                if len(al) == 2:
                    per_arm_pair.setdefault((a, b), []).append((d[0], sum(al) / 2.0))
                    for h1 in (1, 2):
                        for h2 in (1, 2):
                            dd = window_distance(same_sample, s, h1, a, h2, b)
                            if dd is not None:
                                align_pool.append(
                                    (raw, window_name(s, h1, a), window_name(s, h2, b),
                                     dd[0], jaccard_from_shared(dd[1], dd[2])))

            # family aggregate
            al_d = [v[0] for v in allelic.values()]
            pa_d = [v[0] for v in paralog.values()]
            al_j = [jaccard_from_shared(v[1], v[2]) for v in allelic.values()]
            pa_j = [jaccard_from_shared(v[1], v[2]) for v in paralog.values()]
            be_d = [x[0] for v in between.values() for x in v]
            matched = [(p, a) for lst in per_arm_pair.values() for (p, a) in lst]
            n_pc = sum(1 for p, a in matched if p < a)
            n_ac = sum(1 for p, a in matched if a < p)
            emit(family=label, substrate=sub, row_type="family_aggregate",
                 community="all",
                 arm_a=";".join(arms), arm_b="",
                 n_allelic_pairs=len(al_d),
                 median_allelic_dist="%.6g" % statistics.median(al_d) if al_d else "",
                 median_allelic_jaccard="%.6g" % statistics.median(al_j) if al_j else "",
                 n_paralog_pairs=len(pa_d),
                 median_paralog_dist="%.6g" % statistics.median(pa_d) if pa_d else "",
                 median_paralog_jaccard="%.6g" % statistics.median(pa_j) if pa_j else "",
                 n_between_pairs=len(be_d),
                 median_between_individual_dist=("%.6g" % statistics.median(be_d)) if be_d else "",
                 n_matched_units=len(matched), n_paralog_closer=n_pc, n_allelic_closer=n_ac,
                 frac_paralog_closer=("%.4f" % (n_pc / len(matched))) if matched else "",
                 sign_test_p="%.6g" % sign_test(n_pc, n_ac) if matched else "",
                 wilcoxon_p=("%.6g" % wilcoxon_less([p for p, _ in matched],
                                                    [a for _, a in matched])) if matched else "",
                 chm13_copies_in_arm=sum(copies_in_arm.values()),
                 chm13_copies_in_chm13_window=sum(copies_in_win.values()),
                 chm13_window_bp=sum(win_bp.values()),
                 mash_k=MASH_K, mash_sketch=MASH_SKETCH)

            # per arm-pair rows
            for (a, b), lst in sorted(per_arm_pair.items()):
                al = [allelic[(s, a)][0] for (s, aa, bb) in paralog if (aa, bb) == (a, b)
                      and (s, a) in allelic]
                al += [allelic[(s, b)][0] for (s, aa, bb) in paralog if (aa, bb) == (a, b)
                       and (s, b) in allelic]
                alj = [jaccard_from_shared(*allelic[(s, a)][1:])
                       for (s, aa, bb) in paralog if (aa, bb) == (a, b) and (s, a) in allelic]
                alj += [jaccard_from_shared(*allelic[(s, b)][1:])
                        for (s, aa, bb) in paralog if (aa, bb) == (a, b) and (s, b) in allelic]
                pa = [d[0] for (s, aa, bb), d in paralog.items() if (aa, bb) == (a, b)]
                paj = [jaccard_from_shared(d[1], d[2])
                       for (s, aa, bb), d in paralog.items() if (aa, bb) == (a, b)]
                pc = sum(1 for p, q in lst if p < q)
                ac = sum(1 for p, q in lst if q < p)
                emit(family=label, substrate=sub, row_type="arm_pair",
                     community=communities.get(arm_to_community_label(a), ""),
                     arm_a=a, arm_b=b,
                     n_allelic_pairs=len(al),
                     median_allelic_dist="%.6g" % statistics.median(al) if al else "",
                     median_allelic_jaccard="%.6g" % statistics.median(alj) if alj else "",
                     n_paralog_pairs=len(pa),
                     median_paralog_dist="%.6g" % statistics.median(pa) if pa else "",
                     median_paralog_jaccard="%.6g" % statistics.median(paj) if paj else "",
                     n_matched_units=len(lst), n_paralog_closer=pc, n_allelic_closer=ac,
                     frac_paralog_closer=("%.4f" % (pc / len(lst))) if lst else "",
                     sign_test_p="%.6g" % sign_test(pc, ac) if lst else "",
                     wilcoxon_p=("%.6g" % wilcoxon_less([p for p, _ in lst],
                                                        [q for _, q in lst])) if lst else "",
                     chm13_copies_in_arm=copies_in_arm[a] + copies_in_arm[b],
                     chm13_copies_in_chm13_window=copies_in_win[a] + copies_in_win[b],
                     chm13_window_bp=win_bp[a] + win_bp[b],
                     mash_k=MASH_K, mash_sketch=MASH_SKETCH)

    # ---- optional alignment cross-check (mash distance vs wfmash identity) -----------
    if args.align_check:
        alignment_crosscheck(align_pool, args.align_check, workdir)
    sys.stderr.write("[workdir] %s\n" % workdir)


if __name__ == "__main__":
    main()
