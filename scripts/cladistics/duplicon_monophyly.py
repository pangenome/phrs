#!/usr/bin/env python3
"""Duplicon monophyly test -- Referee #1, major point 1 (Nature review, PEER_REVIEW.md).

Referee #1: "Concerted evolution is never tested. It requires that paralogous copies at
different ends be more similar than their divergence times predict, i.e. actively
homogenized; high identity alone does not show this, since recent duplication gives the
same signature. The data allow a direct test: duplicon phylogenies (D4Z4/DUX4L,
TUBB8/TUBB8B, OR4F, SEPTIN14) asking whether copies cluster by end or interdigitate
across ends within a haplotype ..."

This script builds copy-level phylogenies and asks whether copies from the same
chromosome arm form a monophyletic clade.

WHAT THE TOPOLOGY CAN AND CANNOT ESTABLISH
------------------------------------------
* Reciprocal monophyly of the arm classes over the whole haplotype panel  =>  the copies
  on arm A and arm B last shared an ancestor before the sampled haplotypes' common
  ancestor.  No between-arm homogenisation within the sampled genealogy; ancient
  duplication + drift explains the data.
* Arm classes para-/polyphyletic while haplotypes are monophyletic  =>  copies track the
  haplotypes, i.e. a haplotype's copies on different arms are closer to each other than
  to their orthologues on other haplotypes.  That is the referee's "interdigitate across
  ends within a haplotype", the topology produced by between-arm homogenisation.  It is
  ALSO what a single recent cross-arm duplication followed by drift produces: the
  topology cannot date an event, count events, or infer direction.  Both hypotheses are
  excluded by reciprocal monophyly; paraphyly excludes neither.
* A homogenised state shared by ALL sampled haplotypes bounds the last event as older
  than the sample's common ancestor.  A state that tracks haplotypes does not.

METHOD (all parameters are echoed into the JSON output)
-------------------------------------------------------
1. Window substrates (read-only, per haplotype, one window per chromosome arm):
     COPIES : /moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz
              telomere-anchored, PHR-agnostic, max 500 kb, one record per
              (sample, haplotype, arm).  PHR-agnostic on purpose: no ascertainment on
              PHR calling and the identical window definition for every haplotype.
     FLANKS : /moosefs/guarracino/HPRCv2/PHR_III/flanking/hprcv2.1Mb.flanking_100kb.fa.gz
              the 100 kb centromere-proximal neighbours of the PHR windows -- the
              control substrate.
   The flank file cannot carry the copies (its CHM13 4q window ends at 193392741, the
   D4Z4 array is chr4:193389387-193553139); it is used for the (c) machinery control.
2. Copy discovery.  A family reference sequence is cut from CHM13 (coordinates fixed
   below, taken from data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz).  FastGA aligns the
   reference against every selected window; alignments of >= 500 bp and >= 0.80 identity
   are de-duplicated to <-20 % target overlap, giving discrete copy occurrences per
   window, capped at 10 per window (evenly spaced along the array) to keep trees
   tractable.  Two aligner regimes are always run: FastGA -f512 (sensitive) and
   FastGA -f8192 (strict seed cutoff).
3. Alignment frame.  Every copy is projected onto the reference coordinates from the
   long-form cs tag of its own FastGA alignment, giving a reference-anchored alignment
   matrix (single-reference star alignment).  Reference columns absent from a copy are
   MISSING; reference columns deleted in a copy are GAP (a difference); unaligned target
   insertions are dropped.  Nothing is imputed.
4. Distances.  p-distance over pairwise-complete reference columns (gap-gap equal,
   gap-base different); pairs with < 200 comparable columns take a documented
   saturation cap of 1.0 and are counted.
5. Trees.  Neighbour-joining (Saitou & Nei 1987) implemented here in numpy, on the point
   distance matrix and on >= 100 alignment-BLOCK bootstrap replicates (50-column blocks
   resampled with replacement).  Every run also recomputes the tree on high-confidence
   blocks only (blocks with >= 80 % of units represented) and on high-confidence units
   only (reference alignment identity >= 0.90) -- the required confidence stratification.
6. Monophyly statistics per label (chromosome arm and, separately, haplotype):
      mono        : does some tree edge separate exactly this label's copies? (1/0)
      extra_clades: minimum number of all-same-label clades needed to cover the label
                    (= 1 iff monophyletic; greedy set cover, exact for the ==1 test)
      frac_in_largest_same_label_clade : |largest all-label clade| / |label|
   reported for the point tree and as a fraction over the bootstrap replicates, plus a
   label-permutation null on the point tree (same label sizes, N_PERM draws).
7. Distance-class table: mean/median p-distance for same-arm/same-haplotype,
   same-arm/different-haplotype, different-arm/same-haplotype and
   different-arm/different-haplotype copy pairs.  This is the referee's second suggested
   test (within-individual non-homologous divergence vs between-individual divergence at
   the same end) at copy resolution, on the same alignments as the tree.
8. Control (c), flanks: for every arm of every family the identical discovery+projection+
   block-bootstrap+NJ machinery is run on the 100 kb flank windows.  The flank is unique
   single-copy sequence, so the testable truth there is ALLELIC monophyly (each
   individual's two haplotypes must be sisters).  Cross-arm flank alignments are counted
   and are expected to be zero, which is why a cross-chromosome flank tree cannot be
   built at all; that is also why the flank substrate cannot test homogenisation (cf.
   response/notes/M1_concerted_evolution.md).

VERIFICATION
------------
`--self-test` runs brute-force oracles and adversarial cases at tiny scale: cs/CIGAR
projection against hand-computed alignments, NJ against Biopython's implementation and
against analytically known trees, monophyly/extra_clades against hand-built trees
(monophyletic, 2-clade paraphyletic, polyphyletic, all-tips, singleton, tie), block
bootstrap recovery of a known split and rejection of a shuffled label, dedupe logic, and
the saturated-pair cap.

Usage
-----
  duplicon_monophyly.py                     # full run, writes the TSV and JSON
  duplicon_monophyly.py --self-test         # oracle + adversarial unit tests only
  duplicon_monophyly.py --workscratch DIR   # where intermediate FASTAs live
"""

import argparse
import gzip
import json
import math
import os
import random
import re
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter, OrderedDict, defaultdict

import numpy as np

# --------------------------------------------------------------------------------------
# Configuration.  Everything here is echoed into the JSON output.
# --------------------------------------------------------------------------------------

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FASTGA = os.path.expanduser("~/.cargo/bin/FastGA")

COPIES_FA = ("/moosefs/guarracino/HPRCv2/PHR_III/"
             "hprcv2.1Mb.telo_500kb_trimmed.fa.gz")
FLANK_FA = ("/moosefs/guarracino/HPRCv2/PHR_III/flanking/"
            "hprcv2.1Mb.flanking_100kb.fa.gz")
SUPERPOP_TSV = os.path.join(
    REPO, "slides", "v2-review-zoom", "_revision_assets", "v9",
    "labels_superpop_stats_polish", "nearest_same_superpop_mds_distances.tsv")
GFF3 = os.path.join(REPO, "data", "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
COMMUNITIES = os.path.join(REPO, "data", "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
ANNOT_BED = os.path.join(REPO, "data", "chm13-annotations.bed")

# --- aligner / discovery parameters ---
REGIMES = [512, 8192]          # FastGA -f (adaptive seed count cutoff)
MIN_ALN = 500                  # FastGA -l500 and minimum unit length
FASTGA_MIN_ID = 0.70           # FastGA -i0.7
MIN_UNIT_ID = 0.80             # identity filter applied to discovered units
DEDUPE_MAX_OVERLAP = 0.20      # target-overlap fraction that makes two hits the same copy
MAX_COPIES_PER_WINDOW = 10     # evenly spaced subsample of an array
FASTGA_THREADS = 32
FASTGA_TIMEOUT = 7200

# --- tree / statistics parameters ---
BLOCK = 50                     # alignment block size for the bootstrap (columns)
MIN_BLOCK_COMPLETENESS = 0.80  # a "high-confidence" block has >=80% of units represented
HIGHCONF_UNIT_ID = 0.90        # a "high-confidence" unit aligned at >=90% identity
MIN_COMPARABLE_COLS = 200      # below this a pair takes the saturation cap
P_DIST_CAP = 1.0               # documented saturation cap
GAP_AS_MISSING = True          # complete deletion: gaps are not scored as differences
N_BOOT = 100
N_PERM = 100
SEED = 20261023

# --- sampling ---
MAX_HAPS = 20                  # haplotypes per family when a common set exists
MIN_COMMON_HAPS = 8            # below this, fall back to per-arm sampling
MAX_HAPS_PER_ARM = 8           # per-arm mode (acrocentric control needs this)
REF_MAX = 4000                 # reference frame is clipped to this many bp

# Families.  `ref` = (chrom, start, end) of the CHM13 reference unit; `arms` = the
# chromosome arms whose copies are compared; `query` = the CHM13 window used to cut the
# reference; `comm` = the PHR community (data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv).
#
# D4Z4_DUX4L  : reference = the D4Z4 repeat unit around DUX4L6 (chr4:193495415-193496699
#               is the annotated gene; 3.3 kb is the tandem period measured from FastGA).
# TUBB8_TUBB8B: reference = the first 4 kb of TUBB8 (chr10:39416-69557 is the gene span;
#               the 4 kb from the gene start was verified to hit the TUBB8B block on 18p).
# OR4F        : reference = the first 4 kb of OR4F3 (chr5:181905487-181933378 gene span).
# SEPTIN14    : reference = SEPTIN14P5 (chr5:181897110-181899694, 2585 bp).
# ACRO_P_REPEAT: reference chosen at run time, see pick_acro_reference().
FAMILIES = [
    dict(label="D4Z4_DUX4L", comm="C1",
         ref=("chr4", 193494500, 193497800),
         ref_note="D4Z4 tandem unit centred on DUX4L6 (GFF3 chr4:193495415-193496699)",
         arms=["chr4_qarm", "chr10_qarm"]),
    dict(label="TUBB8_TUBB8B", comm="C2",
         ref=("chr10", 39416, 43416),
         ref_note="first 4 kb of TUBB8 (GFF3 gene span chr10:39416-69557)",
         arms=["chr10_parm", "chr18_parm"]),
    dict(label="OR4F", comm="C11",
         ref=("chr5", 181905487, 181909487),
         ref_note="first 4 kb of OR4F3 (GFF3 gene span chr5:181905487-181933378)",
         arms=["chr5_qarm", "chr6_qarm"]),
    dict(label="SEPTIN14", comm="C3/C11",
         ref=("chr5", 181897110, 181899694),
         ref_note="SEPTIN14P5 (GFF3 chr5:181897110-181899694)",
         arms=["chr1_parm", "chr2_qarm", "chr3_qarm", "chr5_qarm", "chr7_parm", "chr9_qarm"]),
    dict(label="ACRO_P_REPEAT", comm="C7", acro=True,
         ref=None,
         ref_note="data-driven: most self-similar 2 kb tile of the CHM13 chr13_parm window",
         arms=["chr13_parm", "chr14_parm", "chr15_parm", "chr21_parm", "chr22_parm"]),
]

ACRO_PICK_ARM = "chr13_parm"
ACRO_TILE = 2000
ACRO_KMER = 40
ACRO_MIN_SHARED = 0.50

ARM_RX = re.compile(r"(chr[0-9XY]+_[pq]arm)")


# --------------------------------------------------------------------------------------
# Infrastructure
# --------------------------------------------------------------------------------------

def die(msg):
    sys.stderr.write("FATAL: %s\n" % msg)
    sys.stderr.flush()
    sys.exit(1)


def log(msg):
    sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    sys.stderr.flush()


def run(cmd, timeout=FASTGA_TIMEOUT, stdin_text=None):
    t0 = time.time()
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             timeout=timeout, input=stdin_text)
    except subprocess.TimeoutExpired:
        die("timeout (>%ss): %s" % (timeout, " ".join(cmd)))
    if res.returncode != 0:
        die("command failed (%d): %s\n%s" % (res.returncode, " ".join(cmd),
                                             res.stderr.decode("utf8", "replace")[-4000:]))
    log("ran (%0.1fs): %s" % (time.time() - t0, " ".join(cmd)))
    return res


def tool_version(cmd):
    try:
        res = subprocess.run(cmd, capture_output=True, timeout=120)
    except Exception as exc:                                  # noqa: BLE001
        return "unavailable (%s)" % exc
    txt = (res.stdout + res.stderr).decode("utf8", "replace").strip().splitlines()
    return txt[0] if txt else "no output"


def read_fai(fa):
    fai = fa + ".fai"
    if not os.path.exists(fai):
        die("missing .fai for %s" % fa)
    out = []
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1])))
    return out


def parse_window_name(name):
    """(sample, hap, arm, dict) for a substrate record name. Dies if unparseable."""
    parts = name.split("#")
    if len(parts) < 3:
        die("cannot parse window name: %r" % name)
    m = ARM_RX.search(parts[2])
    if not m:
        die("cannot parse arm from window name: %r" % name)
    arm = m.group(1)
    return parts[0], parts[1], arm, {"contig": parts[2].split(":")[0],
                                     "chrom": arm.split("_")[0],
                                     "raw": name}


def read_fasta(path):
    out = OrderedDict()
    hdr = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    out[hdr] = "".join(buf)
                hdr = line[1:]
                buf = []
            else:
                buf.append(line)
    if hdr is not None:
        out[hdr] = "".join(buf)
    return out


def write_fasta(path, items):
    with open(path, "w") as fh:
        for name, seq in items:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def load_superpops():
    """(sample, hap) -> superpopulation, from the repo-local label table."""
    sup = {}
    with open(SUPERPOP_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_s = header.index("Sample")
            i_h = header.index("Haplotype")
            i_p = header.index("Superpopulation")
        except ValueError:
            die("unexpected header in %s: %r" % (SUPERPOP_TSV, header))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sup[(f[i_s], f[i_h])] = f[i_p]
    if not sup:
        die("empty superpopulation table")
    return sup


def load_communities():
    out = {}
    with open(COMMUNITIES) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = f[1]
    return out


def stratified_pick(cands, n, sup):
    """Deterministically pick <=n (sample, hap) keys, round-robin across superpopulations."""
    by_sup = defaultdict(list)
    for k in cands:
        by_sup[sup.get(k, "NA")].append(k)
    for v in by_sup.values():
        v.sort()
    order = sorted(by_sup)
    picked = []
    idx = 0
    while len(picked) < n:
        added = False
        for s in order:
            if idx < len(by_sup[s]) and len(picked) < n:
                picked.append(by_sup[s][idx])
                added = True
        if not added:
            break
        idx += 1
    picked.sort()
    return picked


# --------------------------------------------------------------------------------------
# FastGA wrapper
# --------------------------------------------------------------------------------------

CS_OPS = set("=*+-")


def fastga(query_fa, target_fa, seed_cutoff, out_paf, extra=None):
    cmd = [FASTGA, "-pafS", "-l%d" % MIN_ALN, "-i%g" % FASTGA_MIN_ID,
           "-f%d" % seed_cutoff, "-T%d" % FASTGA_THREADS]
    if extra:
        cmd += extra
    cmd += [query_fa, target_fa]
    res = run(cmd)
    with open(out_paf, "wb") as fh:
        fh.write(res.stdout)
    return out_paf


def parse_cs(cs, convention="query_plus"):
    """Parse a long-form cs tag into a list of (op, seq) with op in {=,*,+,-}.

    convention="query_plus": '+' consumes query bases not in target, '-' consumes target
    bases not in query.  convention="query_minus" is the swap (used if the assertion in
    alignment_columns() fails under the default, which cannot be assumed blindly).
    """
    ops = []
    i = 0
    n = len(cs)
    while i < n:
        ch = cs[i]
        if ch not in CS_OPS:
            die("unexpected character %r in cs tag at %d: %r" % (ch, i, cs[:80]))
        if ch == "*":
            if i + 2 >= n:
                die("truncated mismatch in cs tag: %r" % cs[max(0, i - 20):i + 20])
            ops.append(("*", cs[i + 1:i + 3]))
            i += 3
            continue
        j = i + 1
        while j < n and cs[j] not in CS_OPS:
            j += 1
        ops.append((ch, cs[i + 1:j]))
        i = j
    return ops


def cs_consumption_convention(ops, convention):
    q = t = 0
    for op, s in ops:
        if op == "=":
            q += len(s)
            t += len(s)
        elif op == "*":
            q += 1
            t += 1
        elif op == "+":
            if convention == "query_plus":
                q += len(s)
            else:
                t += len(s)
        else:                                                  # '-'
            if convention == "query_plus":
                t += len(s)
            else:
                q += len(s)
    return q, t


CS_CONVENTIONS = ["query_plus", "query_minus"]


def alignment_columns(rec, target_seq, convention, l_ref):
    """Project one alignment onto reference (query) coordinates.

    Returns dict with:
      refcols    : list of reference columns that are ALIGNED (base or gap)
      codes      : list of base codes (0..3 = ACGT, 4 = gap, 5 = missing) parallel to refcols
      identity, alnlen, nmatch, tstart, tend, strand, tname, qlen...
    Raises ValueError if the cs tag does not reconcile with qstart/qend/tstart/tend.
    """
    ops = parse_cs(rec["cs"])
    q_cons, t_cons = cs_consumption_convention(ops, convention)
    qspan = rec["qend"] - rec["qstart"]
    tspan = rec["tend"] - rec["tstart"]
    if q_cons != qspan or t_cons != tspan:
        raise ValueError("cs consumption (%d,%d) != spans (%d,%d)"
                         % (q_cons, t_cons, qspan, tspan))

    qstart, qend = rec["qstart"], rec["qend"]
    strand = rec["strand"]
    nb = 0                       # query bases consumed so far
    tp = rec["tstart"]           # target position, ascending
    refcols = []
    codes = []

    def refcol_of(nb):
        return (qend - 1 - nb) if strand == "-" else (qstart + nb)

    base_code = {"A": 0, "C": 1, "G": 2, "T": 3}
    comp_base = {"A": "T", "C": "G", "G": "C", "T": "A"}

    for op, s in ops:
        if op == "=" or op == "*":
            # '=' is an identical run; '*' is exactly one query/target base pair.  Both
            # advance query and target together.  The BASE is taken from the target
            # FASTA at the walked position, never from the cs string: that makes the
            # projection independent of the aligner's base ordering inside a mismatch,
            # which is an aligner convention, not a fact about the data.
            n = len(s) if op == "=" else 1
            for _k in range(n):
                rc = refcol_of(nb)
                if rc < 0 or rc >= l_ref:
                    raise ValueError("reference column %d outside frame" % rc)
                if tp < 0 or tp >= len(target_seq):
                    raise ValueError("target position %d outside the target" % tp)
                # An inverted copy is homologous to the reverse complement of the
                # reference, so its bases are COMPLEMENTED to be comparable on the
                # reference frame.  Without this, an inverted copy lands on the frame as
                # its own reverse complement and every column "differs" (measured: a
                # SEPTIN14 chr1p copy vs a chr3q copy gave p = 1.0000 over 2583 columns
                # before the fix).
                tb = target_seq[tp].upper()
                if strand == "-":
                    tb = comp_base.get(tb, "N")
                refcols.append(rc)
                codes.append(base_code.get(tb, 5))
                nb += 1
                tp += 1
        elif op == "+":                                    # query bases not in target
            n = len(s)
            for k in range(n):
                rc = refcol_of(nb)
                if rc < 0 or rc >= l_ref:
                    raise ValueError("reference column %d outside frame" % rc)
                refcols.append(rc)
                codes.append(4)                            # deleted in this copy
                nb += 1
        else:                                              # target bases not in query
            tp += len(s)
    if tp != rec["tend"]:
        raise ValueError("target walk ended at %d != tend %d" % (tp, rec["tend"]))
    return dict(refcols=refcols, codes=codes)


def parse_paf(path):
    """Parse a FastGA -pafS file into records (with cs)."""
    recs = []
    n_bad = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                n_bad += 1
                continue
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            if cs is None:
                n_bad += 1
                continue
            recs.append(dict(qname=f[0], qlen=int(f[1]), qstart=int(f[2]),
                             qend=int(f[3]), strand=f[4], tname=f[5], tlen=int(f[6]),
                             tstart=int(f[7]), tend=int(f[8]), nmatch=int(f[9]),
                             alnlen=int(f[10]), cs=cs))
    if n_bad:
        log("  %d alignment line(s) dropped (no cs / too few fields)" % n_bad)
    return recs


# --------------------------------------------------------------------------------------
# Substrates
# --------------------------------------------------------------------------------------

class Substrate:
    def __init__(self, path, tag):
        self.path = path
        self.tag = tag
        self.records = read_fai(path)              # [(name, length)]
        self.by_arm = defaultdict(list)            # arm -> [(sample, hap, name, length)]
        self.seq_cache = {}
        for name, length in self.records:
            smp, hap, arm, _meta = parse_window_name(name)
            self.by_arm[arm].append((smp, hap, name, length))
        for arm in self.by_arm:
            self.by_arm[arm].sort(key=lambda r: (r[0], r[1]))
        self.name_of = {}
        for arm in self.by_arm:
            for smp, hap, name, length in self.by_arm[arm]:
                self.name_of[(smp, hap, arm)] = name

    def haps(self, arm):
        return [(s, h) for s, h, _n, _l in self.by_arm.get(arm, [])]

    def extract(self, names, out_fa):
        """Extract records by exact name via samtools faidx -r; uppercased."""
        regions = out_fa + ".regions"
        with open(regions, "w") as fh:
            for n in names:
                fh.write(n + "\n")
        res = run(["samtools", "faidx", "-r", regions, self.path])
        recs = OrderedDict()
        hdr = None
        buf = []
        for line in res.stdout.decode("ascii", "replace").splitlines():
            if line.startswith(">"):
                if hdr is not None:
                    recs[hdr] = "".join(buf)
                hdr = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
        if hdr is not None:
            recs[hdr] = "".join(buf)
        got = set(recs)
        missing = [n for n in names if n not in got]
        if missing:
            die("samtools faidx returned no sequence for %d record(s), e.g. %r "
                "(substrate %s)" % (len(missing), missing[:3], self.path))
        write_fasta(out_fa, list(recs.items()))
        for k, v in recs.items():
            self.seq_cache[k] = v
        return out_fa


def chm13_window(sub, arm):
    for smp, hap, name, length in sub.by_arm.get(arm, []):
        if smp == "CHM13":
            return name, length
    die("no CHM13 window for arm %s in %s" % (arm, sub.path))


def window_coords(name):
    """(contig, start_1based, end_inclusive) parsed from a window record name."""
    m = re.search(r":([0-9]+)-([0-9]+)_", name)
    if not m:
        die("cannot parse coordinates from window name %r" % name)
    return int(m.group(1)), int(m.group(2))


def cut_from_window(seq, win_name, start1, end1):
    """Cut 1-based inclusive CHM13 coordinates out of a window record."""
    ws, we = window_coords(win_name)
    if start1 < ws or end1 > we:
        die("interval %d-%d not inside window %s (%d-%d)"
            % (start1, end1, win_name, ws, we))
    return seq[start1 - ws:end1 - ws + 1]


# --------------------------------------------------------------------------------------
# Copy discovery
# --------------------------------------------------------------------------------------

def discover_copies(ref_fa, target_fa, paf_path, max_copies, target_seqs, l_ref):
    """Turn FastGA alignments into copy units.

    Returns (units, stats) where each unit is a dict with refcols/codes and metadata.
    """
    recs = parse_paf(paf_path)
    conv = None
    by_target = defaultdict(list)
    n_proj_fail = 0
    for rec in recs:
        if rec["alnlen"] < MIN_ALN:
            continue
        ident = rec["nmatch"] / float(rec["alnlen"])
        if ident < MIN_UNIT_ID:
            continue
        tseq = target_seqs.get(rec["tname"])
        if tseq is None:
            die("alignment against unknown target %r" % rec["tname"])
        # determine the cs convention once, from a real alignment, then assert it always
        if conv is None:
            for c in CS_CONVENTIONS:
                try:
                    proj = alignment_columns(rec, tseq, c, l_ref)
                    conv = c
                    break
                except ValueError:
                    proj = None
            if conv is None:
                die("cannot reconcile the cs tag with the PAF spans under either "
                    "convention for %r; refusing to guess" % rec["tname"])
        else:
            try:
                proj = alignment_columns(rec, tseq, conv, l_ref)
            except ValueError:
                n_proj_fail += 1
                continue
        by_target[rec["tname"]].append(dict(
            tstart=rec["tstart"], tend=rec["tend"], ident=ident, alnlen=rec["alnlen"],
            nmatch=rec["nmatch"], refcols=proj["refcols"], codes=proj["codes"],
            strand=rec["strand"]))
    if n_proj_fail:
        log("  %d alignment(s) failed projection under the established convention"
            % n_proj_fail)

    units = []
    n_dedup = 0
    for tname in sorted(by_target):
        hits = sorted(by_target[tname], key=lambda h: (-h["nmatch"], h["tstart"]))
        kept = []
        for h in hits:
            if h["tend"] - h["tstart"] < MIN_ALN:
                continue
            bad = False
            for k in kept:
                ov = min(h["tend"], k["tend"]) - max(h["tstart"], k["tstart"])
                if ov > DEDUPE_MAX_OVERLAP * min(h["tend"] - h["tstart"],
                                                 k["tend"] - k["tstart"]):
                    bad = True
                    break
            if bad:
                n_dedup += 1
                continue
            kept.append(h)
        kept.sort(key=lambda h: h["tstart"])
        if len(kept) > max_copies:
            idx = np.linspace(0, len(kept) - 1, max_copies).round().astype(int)
            idx = sorted(set(int(i) for i in idx))
            kept = [kept[i] for i in idx]
        smp, hap, arm, _meta = parse_window_name(tname)
        for h in kept:
            units.append(dict(target=tname, sample=smp, hap=hap, arm=arm,
                              chrom=arm.split("_")[0],
                              tstart=h["tstart"], tend=h["tend"],
                              ident=h["ident"], alnlen=h["alnlen"],
                              strand=h["strand"], refcols=h["refcols"],
                              codes=h["codes"]))
    per_target = Counter(r["tname"] for r in recs
                         if r["alnlen"] >= MIN_ALN
                         and r["nmatch"] / float(r["alnlen"]) >= MIN_UNIT_ID)
    stats = dict(n_alignments=len(recs), n_deduped=n_dedup,
                 cs_convention=conv, n_projection_failures=n_proj_fail,
                 n_targets=len(per_target),
                 median_qualifying_alignments_per_target=(
                     float(np.median(list(per_target.values())))
                     if per_target else 0.0),
                 max_qualifying_alignments_per_target=(
                     int(max(per_target.values())) if per_target else 0))
    return units, stats


def build_matrix(units, l_ref):
    """Reference-anchored alignment matrix: int8 (U x L), codes 0..3 ACGT, 4 gap, 5 missing."""
    M = np.full((len(units), l_ref), 5, dtype=np.int8)
    for i, u in enumerate(units):
        for rc, code in zip(u["refcols"], u["codes"]):
            M[i, rc] = code
    return M


# --------------------------------------------------------------------------------------
# Distances
# --------------------------------------------------------------------------------------

def pdists_from_columns(M, cols, cap=P_DIST_CAP, min_cols=MIN_COMPARABLE_COLS,
                        gap_as_missing=GAP_AS_MISSING):
    """Pairwise p-distance over the given reference columns (pairwise complete).

    With gap_as_missing (the default and the primary model) a column is comparable only
    when BOTH copies carry a sequenced base, so the distance is substitutions per aligned
    base -- complete deletion.  This matters here: tandem-array copies differ by indels of
    10-100 bp, and scoring an indel column as one difference per deleted base would let a
    single indel outvote tens of substitutions.  With gap_as_missing=False, gap-gap counts
    as equal and gap-base as a difference (documented alternative model).
    Pairs with < min_cols comparable columns take the documented cap.
    """
    sub = M[:, cols]
    if gap_as_missing:
        present = ((sub >= 0) & (sub <= 3)).astype(np.float32)
        n_cmp = present @ present.T
        n_eq = np.zeros(n_cmp.shape, dtype=np.float32)
        for code in (0, 1, 2, 3):
            onehot = (sub == code).astype(np.float32)
            n_eq += onehot @ onehot.T
    else:
        present = (sub != 5).astype(np.float32)
        n_cmp = present @ present.T
        n_eq = np.zeros(n_cmp.shape, dtype=np.float32)
        for code in (0, 1, 2, 3, 4):
            onehot = (sub == code).astype(np.float32)
            n_eq += onehot @ onehot.T
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (n_cmp - n_eq) / n_cmp
    undefined = n_cmp < min_cols
    d[undefined] = cap
    d[np.isnan(d)] = cap
    np.fill_diagonal(d, 0.0)
    return d.astype(np.float64), n_cmp


# --------------------------------------------------------------------------------------
# Neighbour joining (Saitou & Nei 1987), numpy
# --------------------------------------------------------------------------------------

def nj_tree(D):
    """Neighbour-joining on a full distance matrix.

    Returns adjacency dict {node_id: set(node_id)} where leaves are 0..n-1.
    The final three nodes are joined in a trifurcation (unrooted tree).
    """
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n < 2:
        die("NJ needs >=2 taxa, got %d" % n)
    if n == 2:
        return {0: {1}, 1: {0}}
    size = 2 * n - 2
    DW = np.zeros((size, size), dtype=float)
    DW[:n, :n] = D
    adj = defaultdict(set)
    active = list(range(n))
    next_id = n
    while len(active) > 3:
        m = len(active)
        sub = DW[np.ix_(active, active)]
        r = sub.sum(axis=1)
        Q = (m - 2.0) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        flat = int(np.argmin(Q))
        i, j = divmod(flat, m)
        if i > j:
            i, j = j, i
        a, b = active[i], active[j]
        dab = DW[a, b]
        new = next_id
        next_id += 1
        adj[a].add(new)
        adj[new].add(a)
        adj[b].add(new)
        adj[new].add(b)
        for k in active:
            if k == a or k == b:
                continue
            dk = 0.5 * (DW[a, k] + DW[b, k] - dab)
            DW[new, k] = dk
            DW[k, new] = dk
        active = [x for x in active if x not in (a, b)] + [new]
        DW[a, :] = 0.0
        DW[:, a] = 0.0
        DW[b, :] = 0.0
        DW[:, b] = 0.0
    a, b, c = active
    centre = next_id
    for x in (a, b, c):
        adj[x].add(centre)
        adj[centre].add(x)
    return dict(adj)


def edge_clades(adj, n_tips):
    """All bipartition sides induced by the tree edges, as integer bitmasks.

    One rooted traversal computes the tip mask below every node, so each edge yields
    (mask_below_child, complement).  Rooted at an internal node, so every edge is one
    parent-child pair.  Bitmasks keep the per-bootstrap cost linear (a Python-set BFS per
    edge made 400-taxon bootstraps 10x slower for identical results).
    """
    if not adj:
        return set()
    internal = [n for n in adj if n >= n_tips]
    root = internal[0] if internal else next(iter(adj))
    parent = {root: None}
    order = []
    stack = [root]
    while stack:
        x = stack.pop()
        order.append(x)
        for y in adj[x]:
            if y not in parent:
                parent[y] = x
                stack.append(y)
    if len(order) != len(adj):
        die("tree traversal did not reach every node (%d of %d); the topology has a "
            "cycle" % (len(order), len(adj)))
    mask = {}
    for x in reversed(order):
        if x < n_tips:
            mask[x] = 1 << x
        else:
            m = 0
            for y in adj[x]:
                if parent.get(y) == x:
                    m |= mask[y]
            mask[x] = m
    full = (1 << n_tips) - 1
    clades = set()
    for x, par in parent.items():
        if par is None:
            continue
        m = mask[x]
        if 0 < m < full:
            clades.add(m)
            clades.add(full & ~m)
    return clades


def label_stats(clades, label_tips, n_tips):
    """Monophyly statistics for one label (integer bitmask interface).

    Returns (mono, extra_clades, frac_in_largest_same_label_clade).
    `extra_clades` is a greedy set-cover count of all-same-label clades needed to cover the
    label; the ==1 test is exact (a cover of size 1 exists iff such a clade exists).
    Greedy is an upper bound in general, which the notes state.
    """
    L = label_tips
    size = bin(L).count("1")
    if size == 0:
        return (None, None, None)
    if size == n_tips:
        return (None, None, 1.0)               # trivially every tip: uninformative
    covers = []
    best = 0
    for cl in clades:
        if cl & ~L:
            continue                            # not fully inside the label
        k = bin(cl).count("1")
        if k > best:
            best = k
        covers.append(cl)
    frac = best / float(size)
    if best == size:
        return (1, 1, 1.0)
    uncovered = L
    count = 0
    while uncovered:
        pick = None
        bestk = 0
        for cl in covers:
            k = bin(cl & uncovered).count("1")
            if k > bestk:
                bestk = k
                pick = cl
        if bestk == 0:
            count += bin(uncovered).count("1")  # lone tips
            break
        uncovered &= ~pick
        count += 1
    return (0, count, frac)


# --------------------------------------------------------------------------------------
# Bootstrap + statistics driver
# --------------------------------------------------------------------------------------

def block_index(l_ref, block=BLOCK):
    blocks = []
    for s in range(0, l_ref, block):
        blocks.append(list(range(s, min(l_ref, s + block))))
    return blocks


def usable_blocks(M, block=BLOCK, min_completeness=0.0):
    """Blocks (as column lists) that carry data for >=2 units.

    With min_completeness>0 a block is kept only when >= min_completeness of the units have
    a SEQUENCED BASE (not missing, not a gap) in it -- the high-confidence-block
    stratification.  Gap-aware on purpose: without it, indel-heavy tandem-array blocks pass
    any completeness threshold and the stratum is degenerate (measured, not assumed).
    Returns (blocks, diagnostics).
    """
    present = (M != 5)
    seqd = (M >= 0) & (M <= 3)
    n_units = M.shape[0]
    out = []
    fracs = []
    for cols in block_index(M.shape[1], block):
        n_data = present[:, cols].sum(axis=0)
        if (n_data >= 2).sum() < 2:
            continue
        frac = seqd[:, cols].mean()
        fracs.append(frac)
        if frac < min_completeness:
            continue
        out.append(cols)
    return out, dict(n_blocks_total=len(fracs),
                     completeness_mean=(float(np.mean(fracs)) if fracs else None),
                     completeness_min=(float(np.min(fracs)) if fracs else None),
                     completeness_median=(float(np.median(fracs)) if fracs else None))


def highconf_block_set(M, min_completeness=MIN_BLOCK_COMPLETENESS):
    """Blocks for the high-confidence stratification, with a data-driven threshold.

    A fixed threshold can be degenerate (measured on the D4Z4 array: every block had
    completeness >= 0.887, so a 0.80 cut kept all of them and the stratum was identical to
    the unfiltered one).  The threshold used is therefore max(min_completeness, median
    block completeness), i.e. the better half of the blocks, and the rule and the realised
    threshold are recorded in the output.
    """
    ab, abinfo = usable_blocks(M)
    med = abinfo.get("completeness_median")
    thr = min_completeness
    rule = "fixed"
    if med is not None and med > thr:
        thr = med
        rule = "max(fixed, median block completeness)"
    hb, hbinfo = usable_blocks(M, min_completeness=thr)
    hbinfo = dict(hbinfo, threshold=thr, rule=rule)
    return ab, abinfo, hb, hbinfo


def analyse_tree_set(M, units, level_key, n_boot=N_BOOT, seed=SEED,
                     blocks=None, tag=""):
    """Bootstrap + point-tree monophyly analysis at one label level.

    level_key: function unit -> label string.
    Returns dict with per-label statistics and the point tree.
    """
    n = M.shape[0]
    labels = [level_key(u) for u in units]
    lab_to_tips = defaultdict(int)
    for i, l in enumerate(labels):
        lab_to_tips[l] |= (1 << i)
    if blocks is None:
        blocks, _bd = usable_blocks(M)
    if len(blocks) < 5:
        return dict(error="too few usable blocks (%d)" % len(blocks), n_units=n)
    all_cols = [c for b in blocks for c in b]

    Dpoint, ncmp = pdists_from_columns(M, all_cols)
    adj = nj_tree(Dpoint)
    clades = edge_clades(adj, n)
    point = {}
    for lab, tips in sorted(lab_to_tips.items()):
        mono, extra, frac = label_stats(clades, tips, n)
        point[lab] = dict(n_copies=bin(tips).count("1"), mono=mono, extra_clades=extra,
                          frac_largest=frac)

    rng = random.Random(seed)
    boot_mono = defaultdict(list)
    boot_extra = defaultdict(list)
    for _ in range(n_boot):
        cols = []
        for _k in range(len(blocks)):
            cols.extend(blocks[rng.randrange(len(blocks))])
        Db, _ = pdists_from_columns(M, cols)
        adjb = nj_tree(Db)
        clb = edge_clades(adjb, n)
        for lab, tips in lab_to_tips.items():
            mono, extra, _frac = label_stats(clb, tips, n)
            if mono is None:
                continue
            boot_mono[lab].append(mono)
            boot_extra[lab].append(extra)

    perm_labels = sorted(lab_to_tips)
    perm_sizes = [bin(lab_to_tips[l]).count("1") for l in perm_labels]
    prng = random.Random(seed + 1)
    perm_mono = defaultdict(list)
    for _ in range(N_PERM):
        idx = list(range(n))
        prng.shuffle(idx)
        pos = 0
        for lab, size in zip(perm_labels, perm_sizes):
            m = 0
            for i in idx[pos:pos + size]:
                m |= (1 << i)
            pos += size
            mono, _e, _f = label_stats(clades, m, n)
            if mono is not None:
                perm_mono[lab].append(mono)

    out = {"n_units": n, "n_blocks": len(blocks), "n_columns": len(all_cols),
           "n_taxa": n, "labels": {}, "tag": tag,
           "mean_pairwise_pdist": float(np.mean(Dpoint[np.triu_indices(n, 1)])),
           "frac_pairs_capped": float(np.mean(ncmp[np.triu_indices(n, 1)] <
                                              MIN_COMPARABLE_COLS)),
           "point_newick": to_newick(adj, n, units)}
    for lab in sorted(lab_to_tips):
        bm = boot_mono.get(lab, [])
        be = boot_extra.get(lab, [])
        pm = perm_mono.get(lab, [])
        p = point[lab]
        out["labels"][lab] = dict(
            n_copies=p["n_copies"],
            n_haplotypes=len({(u["sample"], u["hap"]) for i, u in enumerate(units)
                              if level_key(u) == lab}),
            mono_point=p["mono"],
            extra_clades_point=p["extra_clades"],
            frac_largest_same_label_clade=p["frac_largest"],
            boot_monophyly_frac=(float(np.mean(bm)) if bm else None),
            boot_extra_clades_median=(float(np.median(be)) if be else None),
            boot_extra_clades_max=(int(max(be)) if be else None),
            perm_monophyly_rate=(float(np.mean(pm)) if pm else None),
        )
    return out


def consensus_matrix(M, units):
    """Collapse each window's copies to one majority-rule representative per column.

    The "array state" of a window: majority base per reference column over that window's
    copies, ignoring missing entries (ties -> the smaller code, deterministically);
    columns with no data at all stay missing.  This gives the array-level tree, which is
    the classical duplicon-phylogeny question (do the chrA arrays form a clade?) without
    the unit-level turnover noise that dominates a copy-level topology.
    """
    groups = OrderedDict()
    for i, u in enumerate(units):
        groups.setdefault((u["sample"], u["hap"], u["arm"]), []).append(i)
    ncol = M.shape[1]
    reps = []
    C = np.full((len(groups), ncol), 5, dtype=np.int8)
    for gi, (key, rows) in enumerate(groups.items()):
        sub = M[rows, :]
        counts = np.zeros((5, ncol), dtype=np.int16)
        for code in range(5):
            counts[code] = (sub == code).sum(axis=0)
        data = counts.sum(axis=0) > 0
        arg = counts.argmax(axis=0).astype(np.int8)
        C[gi, data] = arg[data]
        reps.append(dict(sample=key[0], hap=key[1], arm=key[2], n_copies=len(rows)))
    return C, reps


def paralog_vs_ortholog(M, units, cols=None):
    """Referee #1's second suggested test, at copy resolution and on the same alignments.

    For every copy i in window (haplotype h, arm a):
      paralog_d(i) = distance to the CLOSEST copy of h on a DIFFERENT arm
      ortholog_d(i) = distance to the same arm a in OTHER haplotypes (median)
    Homogenisation between arms means paralog_d < ortholog_d.
    """
    if cols is None:
        cols = list(range(M.shape[1]))
    D, ncmp = pdists_from_columns(M, cols)
    ok = ncmp >= MIN_COMPARABLE_COLS
    by_win = defaultdict(list)
    for i, u in enumerate(units):
        by_win[(u["sample"], u["hap"], u["arm"])].append(i)
    haps = sorted({(u["sample"], u["hap"]) for u in units})
    arms = sorted({u["arm"] for u in units})
    pairs = []
    for (smp, hap, arm) in by_win:
        others = [by_win[(smp, hap, a2)] for a2 in arms
                  if a2 != arm and (smp, hap, a2) in by_win]
        if not others:
            continue
        other_rows = [j for lst in others for j in lst]
        for i in by_win[(smp, hap, arm)]:
            pd = [D[i, j] for j in other_rows if ok[i, j]]
            od = [D[i, j] for (s2, h2) in haps if (s2, h2) != (smp, hap)
                  for j in by_win.get((s2, h2, arm), []) if ok[i, j]]
            if pd and len(od) >= 3:
                pairs.append((min(pd), float(np.median(od))))
    if len(pairs) < 10:
        return dict(n=len(pairs), note="too few matched copies")
    p = np.array([x for x, _y in pairs])
    o = np.array([y for _x, y in pairs])
    from scipy.stats import wilcoxon
    try:
        pval = float(wilcoxon(p, o, alternative="less", zero_method="wilcox").pvalue)
    except ValueError:
        pval = float("nan")
    return dict(n=int(len(pairs)),
                median_paralog=float(np.median(p)),
                median_ortholog=float(np.median(o)),
                ratio_paralog_over_ortholog=float(np.median(p) / np.median(o)),
                frac_paralog_farther=float(np.mean(p > o)),
                wilcoxon_p_paralog_less=pval)


def arm_pair_distances(M, units, cols=None):
    """Median p-distance for every arm pair and within each arm.

    This is the direct "which ends are close" readout behind the topology: the diagonal is
    within-arm (orthologous/paralogous copies of one end), an off-diagonal cell is between
    two ends.  Capped pairs (too few comparable columns) are excluded and counted.
    """
    if cols is None:
        cols = list(range(M.shape[1]))
    D, ncmp = pdists_from_columns(M, cols)
    ok = ncmp >= MIN_COMPARABLE_COLS
    acc = defaultdict(list)
    n_cap = Counter()
    arms = sorted({u["arm"] for u in units})
    for i in range(len(units)):
        for j in range(i + 1, len(units)):
            a, b = sorted((units[i]["arm"], units[j]["arm"]))
            key = "%s|%s" % (a, b)
            if not ok[i, j]:
                n_cap[key] += 1
                continue
            acc[key].append(float(D[i, j]))
    out = {}
    for key in sorted(acc) + [k for k in sorted(n_cap) if k not in acc]:
        v = np.array(acc.get(key, []), dtype=float)
        out[key] = dict(n=int(v.size),
                        median=(float(np.median(v)) if v.size else None),
                        mean=(float(v.mean()) if v.size else None),
                        n_capped=int(n_cap.get(key, 0)))
    return out


def nearest_neighbour_summary(D, units, ok, ncmp):
    """Nearest-neighbour composition: a robust companion to the strict clade test.

    For every unit, its nearest neighbour in the distance matrix.  This is interpretable
    where a strict clade test is binary: if a family's copies are separated by arm, then
    nearly every copy's nearest neighbour is on the same arm; the haplotype fractions say
    whether that same-arm neighbour is the individual's own array (inter-arm
    homogenisation) or an unrelated haplotype's array (arm-specific divergence).
    """
    n = len(units)
    same_arm = 0
    same_hap = 0
    same_arm_other_hap = 0
    counted = 0
    detail = []
    for i in range(n):
        best = None
        for j in range(n):
            if i == j or not ok[i, j]:
                continue
            if best is None or D[i, j] < D[i, best]:
                best = j
        if best is None:
            continue
        counted += 1
        ui, uj = units[i], units[best]
        sa = ui["arm"] == uj["arm"]
        sh = (ui["sample"], ui["hap"]) == (uj["sample"], uj["hap"])
        same_arm += int(sa)
        same_hap += int(sh)
        same_arm_other_hap += int(sa and not sh)
        detail.append(dict(unit="%s#%s|%s" % (ui["sample"], ui["hap"], ui["arm"]),
                           nearest="%s#%s|%s" % (uj["sample"], uj["hap"], uj["arm"]),
                           distance=float(D[i, best])))
    if counted == 0:
        return dict(n=0)
    return dict(n=counted,
                frac_nearest_same_arm=same_arm / float(counted),
                frac_nearest_same_haplotype=same_hap / float(counted),
                frac_nearest_same_arm_other_haplotype=same_arm_other_hap / float(counted),
                detail=detail)


def to_newick(adj, n_tips, units=None):
    """Rooted-at-an-arbitrary-edge Newick string of the unrooted tree (for the JSON)."""
    if not adj:
        return ""
    # pick an edge to root: the edge whose removal balances the two sides best
    clades = edge_clades(adj, n_tips)
    if not clades:
        return ""
    root_edge = None
    best = -1
    for (u, v) in sorted({(min(a, b), max(a, b)) for a, s in adj.items() for b in s}):
        stack = [v]
        seen = {u}
        tips = 0
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            if x < n_tips:
                tips += 1
            stack.extend(y for y in adj[x] if y not in seen)
        score = min(tips, n_tips - tips)
        if score > best:
            best = score
            root_edge = (u, v, tips)

    def name_of(i):
        if units is None:
            return "t%d" % i
        u = units[i]
        return "%s|%s_%s|c%d" % (u["sample"], u["hap"], u["arm"], i)

    def rec(node, parent):
        if node < n_tips:
            return name_of(node)
        parts = [rec(c, node) for c in sorted(adj[node]) if c != parent]
        return "(" + ",".join(parts) + ")"

    u, v, _t = root_edge
    return "(%s,%s);" % (rec(u, None), rec(v, None))


def distance_classes(M, units, cols=None):
    """p-distance table by (same arm?, same haplotype?) -- the referee's second test."""
    if cols is None:
        cols = list(range(M.shape[1]))
    D, ncmp = pdists_from_columns(M, cols)
    n = M.shape[0]
    acc = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            ui, uj = units[i], units[j]
            same_arm = ui["arm"] == uj["arm"]
            same_hap = (ui["sample"], ui["hap"]) == (uj["sample"], uj["hap"])
            if same_hap and ui["hap"] != uj["hap"]:
                same_hap = False
            key = ("same_arm_same_hap" if same_arm and same_hap else
                   "same_arm_diff_hap" if same_arm else
                   "diff_arm_same_hap" if same_hap else "diff_arm_diff_hap")
            if ncmp[i, j] >= MIN_COMPARABLE_COLS:
                acc[key].append(float(D[i, j]))
    out = {}
    for key in ("same_arm_same_hap", "same_arm_diff_hap",
                "diff_arm_same_hap", "diff_arm_diff_hap"):
        v = np.array(acc.get(key, []), dtype=float)
        out[key] = dict(n=int(v.size),
                        mean=(float(v.mean()) if v.size else None),
                        median=(float(np.median(v)) if v.size else None))
    return out


# --------------------------------------------------------------------------------------
# Pairwise validation of the reference-frame p-distance
# --------------------------------------------------------------------------------------

def pairwise_crosscheck(units, M, L_ref, n_pairs=100, seed=SEED):
    """Compare reference-frame p-distances with direct pairwise FastGA p-distances."""
    rng = random.Random(seed)
    n = len(units)
    if n < 4:
        return dict(n=0)
    pairs = set()
    while len(pairs) < min(n_pairs, n * (n - 1) // 2):
        i = rng.randrange(n)
        j = rng.randrange(n)
        if i != j:
            pairs.add((min(i, j), max(i, j)))
    pairs = sorted(pairs)
    tmpd = os.path.join(WORKDIR, "xcheck")
    os.makedirs(tmpd, exist_ok=True)
    items = []
    for i, u in enumerate(units):
        items.append(("u%d" % i, u["seq"]))
    fa = os.path.join(tmpd, "units.fa")
    write_fasta(fa, items)
    paf = os.path.join(tmpd, "units.paf")
    fastga(fa, fa, REGIMES[0], paf)
    recs = parse_paf(paf)
    obs = {}
    for rec in recs:
        if not (rec["qname"].startswith("u") and rec["tname"].startswith("u")):
            continue
        try:
            qi = int(rec["qname"][1:])
            ti = int(rec["tname"][1:])
        except ValueError:
            continue
        key = (min(qi, ti), max(qi, ti))
        if key not in pairs:
            continue
        if rec["alnlen"] < MIN_ALN:
            continue
        # the SAME distance model as the reference-frame distance: substitutions per
        # column where both copies carry a base (complete deletion, indels not scored)
        ops = parse_cs(rec["cs"])
        n_eq = sum(len(sq) for op, sq in ops if op == "=")
        n_mm = sum(1 for op, _sq in ops if op == "*")
        if n_eq + n_mm < MIN_COMPARABLE_COLS:
            continue
        val = n_mm / float(n_eq + n_mm)
        if key not in obs or (n_eq + n_mm) > obs[key][1]:
            obs[key] = (val, n_eq + n_mm)
    D, _ = pdists_from_columns(M, list(range(L_ref)))
    xs, ys = [], []
    for k in pairs:
        if k in obs:
            xs.append(obs[k][0])
            ys.append(float(D[k[0], k[1]]))
    if len(xs) < 10:
        return dict(n=len(xs), note="too few pairwise alignments")
    xs = np.array(xs)
    ys = np.array(ys)
    from scipy.stats import spearmanr, pearsonr
    return dict(n=int(len(xs)),
                spearman=float(spearmanr(xs, ys).correlation),
                pearson=float(pearsonr(xs, ys)[0]),
                max_abs_dev=float(np.max(np.abs(xs - ys))),
                median_abs_dev=float(np.median(np.abs(xs - ys))))


# --------------------------------------------------------------------------------------
# CHM13 annotation cross-check
# --------------------------------------------------------------------------------------

def chm13_annotated_copies(family, ref_arm):
    """Annotated CHM13 gene copies of a family, with arm assignment and tip distance.

    Released as a real-data oracle: the number of discovered CHM13 copies in the CHM13
    window of `ref_arm` must be of the same order as the annotated count.
    """
    grows = []
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
            grows.append((f[0], int(f[3]), int(f[4]), attrs.get("gene", "")))
    return grows


# --------------------------------------------------------------------------------------
# Reference picking for the acrocentric control
# --------------------------------------------------------------------------------------

def pick_acro_reference(sub, arms, workdir, tile=ACRO_TILE):
    """Choose the reference unit for the acrocentric positive control, by alignment.

    The rDNA arrays are NOT reachable from a telomere-anchored window (CHM13 RNA18SN1 sits at
    chr13:11731724 while the acrocentric PHR blocks in data/chm13-annotations.bed run to
    chr13:13562256, and these windows stop at 500 kb), so the C7 control has to use the
    acrocentric p-arm repeat family that IS in the windows.

    The unit is chosen from the data, not from an annotation: every window of every
    acrocentric arm is cut into `tile`-bp tiles, all tiles are aligned against all those
    windows, and the tile whose homologs are found on the most DISTINCT acrocentric arms wins
    (ties: more homologs, then the earlier arm, then the lower coordinate).

    Two earlier scorings were measured and rejected: (i) scoring tiles on a 40-mer containment
    rule, (ii) scoring tiles of the chr13p window only.  Both can pick a unit that is private
    to one arm -- measured on this substrate, the CHM13 chr13p 500 kb window is a chr13-private
    divergent 36-mer satellite and shares almost nothing with the other acrocentric p arms
    (a FastGA all-vs-all of the five CHM13 acrocentric 500 kb windows at -l500 -i0.8 found
    2,537-3,398 cross-window alignments for every pair among chr14p/chr15p/chr21p/chr22p but
    only 5-18 for any pair involving chr13p).
    """
    panel_names = [chm13_window(sub, a)[0] for a in arms]
    panel_fa = os.path.join(workdir, "acro_pick_panel.fa")
    sub.extract(panel_names, panel_fa)
    panel = read_fasta(panel_fa)
    tiles = []
    for wi, wname in enumerate(panel_names):
        seq = panel[wname]
        for st in range(0, max(1, len(seq) - tile + 1), tile):
            tiles.append((wi, st, seq[st:st + tile]))
    if len(tiles) < 20:
        die("acrocentric pick panel is shorter than 20 tiles")
    qfa = os.path.join(workdir, "acro_pick_tiles.fa")
    write_fasta(qfa, [("W%dT%d" % (wi, st), sq) for wi, st, sq in tiles])
    paf = os.path.join(workdir, "acro_pick_tiles.paf")
    fastga(qfa, panel_fa, REGIMES[0], paf)
    hits = Counter()
    arms_hit = defaultdict(set)
    for rec in parse_paf(paf):
        if rec["alnlen"] < MIN_ALN:
            continue
        if rec["nmatch"] / float(rec["alnlen"]) < MIN_UNIT_ID:
            continue
        wi = panel_names.index(rec["tname"]) if rec["tname"] in panel_names else None
        if wi is None:
            die("acrocentric panel alignment against an unknown target %r"
                % rec["tname"])
        q_wi = int(rec["qname"][1:].split("T")[0])
        q_st = int(rec["qname"].split("T")[1])
        if wi == q_wi and abs(rec["tstart"] - q_st) < tile and \
                rec["nmatch"] == tile and rec["alnlen"] == tile:
            continue                      # the tile's own perfect self-hit
        hits[rec["qname"]] += 1
        arms_hit[rec["qname"]].add(parse_window_name(panel_names[wi])[2])
    best_key = None
    best = None
    for (wi, st, sq) in tiles:
        nm = "W%dT%d" % (wi, st)
        key = (len(arms_hit.get(nm, ())), hits.get(nm, 0), -wi, -st)
        if best_key is None or key > best_key:
            best_key = key
            best = (wi, st)
    if best is None:
        die("acrocentric reference pick failed (%s)" % ",".join(arms))
    wi, st = best
    nm = "W%dT%d" % (wi, st)
    return dict(window_idx=wi, window=panel_names[wi], arm=arms[wi], start=st, tile=tile,
                n_homologs=hits.get(nm, 0), n_arms_hit=len(arms_hit.get(nm, ())),
                arms_hit=sorted(arms_hit.get(nm, ())),
                median_homologs=float(np.median([hits.get("W%dT%d" % (w, t), 0)
                                                 for (w, t, _s) in tiles])),
                n_tiles=len(tiles), panel_windows=panel_names)


# --------------------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------------------

WORKDIR = None


def family_run(fam, copies_sub, flank_sub, sup, out):
    label = fam["label"]
    log("=== family %s (community %s), arms %s" % (label, fam["comm"], ",".join(fam["arms"])))

    # ---- reference sequence -----------------------------------------------------------
    if fam.get("acro"):
        pickinfo = pick_acro_reference(copies_sub, fam["arms"], WORKDIR)
        name = pickinfo["window"]
        start, tilelen = pickinfo["start"], pickinfo["tile"]
        win_fa = os.path.join(WORKDIR, "ref_%s.fa" % label)
        copies_sub.extract([name], win_fa)
        seq = read_fasta(win_fa)[name]
        ref_seq = seq[start:start + tilelen]
        ref_name = "REF_%s_%s_%d-%d" % (label, pickinfo["arm"], start,
                                        start + tilelen)
        ref_meta = dict(source="CHM13 %s window offset %d..%d" % (pickinfo["arm"], start,
                                                                 start + tilelen),
                        chosen_from_arm=pickinfo["arm"], tile=tilelen, pick=pickinfo)
    else:
        chrom, s1, e1 = fam["ref"]
        arm_for_ref = None
        for arm in fam["arms"]:
            if arm.split("_")[0] == chrom:
                arm_for_ref = arm
                break
        if arm_for_ref is None:
            arm_for_ref = fam["arms"][0]
        name, _len = chm13_window(copies_sub, arm_for_ref)
        win_fa = os.path.join(WORKDIR, "ref_%s.fa" % label)
        copies_sub.extract([name], win_fa)
        seq = read_fasta(win_fa)[name]
        e1 = min(e1, s1 + REF_MAX - 1)
        ref_seq = cut_from_window(seq, name, s1, e1)
        ref_name = "REF_%s_%s_%d-%d" % (label, chrom, s1, e1)
        ref_meta = dict(source=("CHM13 %s %s:%d-%d" % (name.split("#")[-1].split("_")[0],
                                                      chrom, s1, e1)),
                        note=fam["ref_note"], window=name)
    L_ref = len(ref_seq)
    if L_ref < MIN_ALN:
        die("reference frame for %s is only %d bp (< %d)" % (label, L_ref, MIN_ALN))
    ref_fa = os.path.join(WORKDIR, "query_%s.fa" % label)
    write_fasta(ref_fa, [(ref_name, ref_seq)])
    out["reference"] = dict(name=ref_name, length=L_ref, sequence_md5=md5of(ref_seq),
                            params=ref_meta)
    log("  reference %s, %d bp" % (ref_name, L_ref))

    # ---- haplotype selection ---------------------------------------------------------
    cand = {arm: set(copies_sub.haps(arm)) for arm in fam["arms"]}
    for arm in fam["arms"]:
        if not cand[arm]:
            die("no windows for arm %s in %s" % (arm, copies_sub.path))
    common = set.intersection(*[cand[a] for a in fam["arms"]])
    common = {k for k in common if k in sup}          # only labelled haplotypes
    if len(common) >= MIN_COMMON_HAPS:
        mode = "common_haplotypes"
        picked = stratified_pick(sorted(common), MAX_HAPS, sup)
        haps = {arm: picked for arm in fam["arms"]}
    else:
        mode = "per_arm"
        haps = {arm: stratified_pick(sorted(cand[arm] & set(sup)), MAX_HAPS_PER_ARM, sup)
                for arm in fam["arms"]}
    out["sampling"] = dict(mode=mode, n_common_haplotypes=len(common),
                           haplotypes={a: ["%s#%s" % (s, h) for s, h in v]
                                       for a, v in haps.items()},
                           superpopulation_of_haplotype={
                               "%s#%s" % (s, h): sup.get((s, h), "NA")
                               for a in fam["arms"] for s, h in haps[a]})
    for arm in fam["arms"]:
        log("  %s: %d haplotype(s)" % (arm, len(haps[arm])))

    # ---- extract windows --------------------------------------------------------------
    wanted = []
    for arm in fam["arms"]:
        for s, h in haps[arm]:
            wanted.append(copies_sub.name_of[(s, h, arm)])
    target_fa = os.path.join(WORKDIR, "target_%s.fa" % label)
    copies_sub.extract(sorted(wanted), target_fa)
    target_seqs = read_fasta(target_fa)
    out["target_windows"] = dict(n=len(target_seqs),
                                 bp=int(sum(len(v) for v in target_seqs.values())))

    for regime in REGIMES:
        rtag = "f%d" % regime
        paf = os.path.join(WORKDIR, "paf_%s_%s.paf" % (label, rtag))
        t0 = time.time()
        fastga(ref_fa, target_fa, regime, paf)
        units, dstat = discover_copies(ref_fa, target_fa, paf, MAX_COPIES_PER_WINDOW,
                                      target_seqs, L_ref)
        # attach sequences (needed for the pairwise cross-check)
        for u in units:
            u["seq"] = target_seqs[u["target"]][u["tstart"]:u["tend"]]
        wall = time.time() - t0
        log("  regime %s: %d copies from %d windows (%.1f s)"
            % (rtag, len(units), len(target_seqs), wall))
        rkey = "regime_%s" % rtag
        out.setdefault("regimes", {})[rkey] = dict(
            fastga_filter=regime, discovery=dstat, wall_seconds=round(wall, 2),
            copies_per_arm=dict(Counter(u["arm"] for u in units)),
            copies_per_haplotype=dict(Counter("%s#%s" % (u["sample"], u["hap"])
                                              for u in units)),
            identity_quantiles=quantiles([u["ident"] for u in units]),
            length_quantiles=quantiles([u["tend"] - u["tstart"] for u in units]),
        )
        if len(units) < 4:
            out["regimes"][rkey]["error"] = "fewer than 4 copies discovered"
            continue
        M = build_matrix(units, L_ref)
        n_kept_cols = int((M != 5).any(axis=0).sum())
        out["regimes"][rkey]["reference_columns_with_data"] = n_kept_cols

        # strata
        ab, abinfo, hb, hbinfo = highconf_block_set(M)
        strata = [("all_blocks", units, M, ab), ("highconf_blocks", units, M, hb)]
        out["regimes"][rkey]["block_strata"] = dict(all_blocks=abinfo,
                                                    highconf_blocks=hbinfo)
        keep = [i for i, u in enumerate(units) if u["ident"] >= HIGHCONF_UNIT_ID]
        if len(keep) >= 6:
            hu, huinfo = usable_blocks(M[keep, :])
            strata.append(("highconf_units", [units[i] for i in keep], M[keep, :], hu))
            out["regimes"][rkey]["block_strata"]["highconf_units"] = huinfo
        out["regimes"][rkey]["strata"] = {}
        for sname, usub, Msub, blocks in strata:
            if len(usub) < 4 or len(blocks) < 5:
                out["regimes"][rkey]["strata"][sname] = dict(
                    error="insufficient data (units=%d blocks=%d)"
                          % (len(usub), len(blocks)))
                continue
            arm_res = analyse_tree_set(Msub, usub, lambda u: u["arm"],
                                       tag="%s|%s|%s|arm" % (label, rtag, sname),
                                       blocks=blocks)
            hap_res = analyse_tree_set(Msub, usub,
                                       lambda u: "%s#%s" % (u["sample"], u["hap"]),
                                       tag="%s|%s|%s|hap" % (label, rtag, sname),
                                       blocks=blocks)
            # array level: one majority-rule representative per window
            C, reps = consensus_matrix(Msub, usub)
            cblocks, _cbinfo = usable_blocks(C)
            array_arm = array_hap = None
            array_nn = None
            if len(reps) >= 6 and len(cblocks) >= 5:
                Darr, ncarr = pdists_from_columns(C, [c for b in cblocks for c in b])
                array_nn = nearest_neighbour_summary(Darr, reps, ncarr >=
                                                     MIN_COMPARABLE_COLS, ncarr)
                array_arm = analyse_tree_set(C, reps, lambda u: u["arm"],
                                             tag="%s|%s|%s|array_arm"
                                                 % (label, rtag, sname),
                                             blocks=cblocks)
                array_hap = analyse_tree_set(C, reps,
                                             lambda u: "%s#%s" % (u["sample"],
                                                                  u["hap"]),
                                             tag="%s|%s|%s|array_hap"
                                                 % (label, rtag, sname),
                                             blocks=cblocks)
            dcls = distance_classes(Msub, usub)
            pvo = paralog_vs_ortholog(Msub, usub)
            apd = arm_pair_distances(Msub, usub, [c for b in blocks for c in b])
            Dnn, ncmpnn = pdists_from_columns(Msub, [c for b in blocks for c in b])
            oknn = ncmpnn >= MIN_COMPARABLE_COLS
            nn_unit = nearest_neighbour_summary(Dnn, usub, oknn, ncmpnn)
            out["regimes"][rkey]["strata"][sname] = dict(
                n_units=len(usub), n_blocks=len(blocks), n_windows=len(reps),
                arm=arm_res, haplotype=hap_res,
                array_arm=array_arm, array_haplotype=array_hap,
                distance_classes=dcls, paralog_vs_ortholog=pvo,
                nearest_neighbour_unit=nn_unit,
                nearest_neighbour_array=array_nn,
                arm_pair_distances=apd)
            log("    %-15s units=%-4d blocks=%-3d arm-mono=%s hap-mono=%s "
                "array-arm-mono=%s  paralog/ortholog=%.2f"
                % (sname, len(usub), len(blocks),
                   fmt_mono(arm_res), fmt_mono(hap_res), fmt_mono(array_arm),
                   (pvo.get("ratio_paralog_over_ortholog") or float("nan"))))
            if sname == "all_blocks":
                out["regimes"][rkey]["pairwise_crosscheck"] = pairwise_crosscheck(
                    usub, Msub, L_ref)
        # persist the point trees compactly (newick is inside the per-level result)

    # ---- flank control ----------------------------------------------------------------
    out["flank_control"] = flank_control(fam, flank_sub, sup, haps, ref_arm=(
        ACRO_PICK_ARM if fam.get("acro") else
        [a for a in fam["arms"] if a.split("_")[0] == fam["ref"][0]][0]))
    return out


def fmt_mono(res):
    if not res or "labels" not in res:
        return "n/a"
    n = len(res["labels"])
    mono = sum(1 for v in res["labels"].values() if v["mono_point"] == 1)
    return "%d/%d" % (mono, n)


def flank_control(fam, flank_sub, sup, haps, ref_arm):
    """Control (c): the identical machinery on the 100 kb flank windows.

    The flank is unique single-copy sequence.  Two facts are measured:
      * cross-arm alignability -- the reference of one arm against the other arms' flank
        windows.  Expected 0: cross-arm flank sequence is not homologous, so a
        cross-chromosome flank tree cannot be built at all (this is also why the flank
        substrate cannot test homogenisation).
      * allelic monophyly -- within one arm the flank windows are orthologous, so each
        individual's two haplotypes must be sisters.  This is the true-answer case that
        validates the alignment/projection/bootstrap/monophyly code on real data.
    """
    res = dict(note="flanks are unique sequence; see docstring", regimes={})
    arms = fam["arms"]
    for regime in REGIMES:
        rtag = "f%d" % regime
        rres = dict(arms={}, cross_arm_alignments={})
        for arm in arms:
            cand = sorted(set(flank_sub.haps(arm)) & set(sup))
            if len(cand) < 6:
                rres["arms"][arm] = dict(error="only %d labelled flanks" % len(cand))
                continue
            picked = stratified_pick(cand, MAX_HAPS, sup)
            wanted = [flank_sub.name_of[(s, h, arm)] for s, h in picked]
            # add other arms' flanks to the same target so cross-arm hits are visible
            other = []
            for a2 in arms:
                if a2 == arm:
                    continue
                c2 = sorted(set(flank_sub.haps(a2)) & set(sup))
                for s, h in stratified_pick(c2, 4, sup):
                    other.append(flank_sub.name_of[(s, h, a2)])
            target_fa = os.path.join(WORKDIR, "flank_target_%s_%s_%s.fa"
                                     % (fam["label"], arm, rtag))
            flank_sub.extract(sorted(set(wanted + other)), target_fa)
            tseq = read_fasta(target_fa)
            chm_name, _l = chm13_window(flank_sub, arm)
            fname = os.path.join(WORKDIR, "flank_ref_%s_%s.fa" % (fam["label"], arm))
            flank_sub.extract([chm_name], fname)
            fseq = read_fasta(fname)[chm_name].upper()
            mid = len(fseq) // 2
            ref = fseq[max(0, mid - REF_MAX // 2):mid + REF_MAX // 2]
            L = len(ref)
            qfa = os.path.join(WORKDIR, "flank_q_%s_%s.fa" % (fam["label"], arm))
            write_fasta(qfa, [("REF_%s_flank" % arm, ref)])
            paf = os.path.join(WORKDIR, "flank_%s_%s_%s.paf" % (fam["label"], arm, rtag))
            fastga(qfa, target_fa, regime, paf)
            units, dstat = discover_copies(qfa, target_fa, paf, 1, tseq, L)
            for u in units:
                u["seq"] = tseq[u["target"]][u["tstart"]:u["tend"]]
            n_cross = sum(1 for u in units if u["arm"] != arm)
            rres["cross_arm_alignments"][arm] = n_cross
            units = [u for u in units if u["arm"] == arm]
            if len(units) < 6:
                rres["arms"][arm] = dict(error="%d flank copies discovered" % len(units),
                                         discovery=dstat)
                continue
            M = build_matrix(units, L)
            blocks, _binfo = usable_blocks(M)
            if len(blocks) < 5:
                rres["arms"][arm] = dict(error="too few flank blocks")
                continue
            ar = analyse_tree_set(M, units, lambda u: "%s#%s" % (u["sample"], u["hap"]),
                                  tag="%s|%s|flank_%s|hap" % (fam["label"], rtag, arm),
                                  blocks=blocks)
            n_lab = len(ar.get("labels", {}))
            n_mono = sum(1 for v in ar.get("labels", {}).values() if v["mono_point"] == 1)
            rres["arms"][arm] = dict(
                discovery=dstat, n_copies=len(units), n_blocks=len(blocks),
                n_individuals_labelled=n_lab, n_allelic_pairs_monophyletic=n_mono,
                allelic_monophyly_frac=(n_mono / float(n_lab) if n_lab else None),
                boot_monophyly_median=float(np.median(
                    [v["boot_monophyly_frac"] for v in ar["labels"].values()
                     if v["boot_monophyly_frac"] is not None])) if n_lab else None,
                labels=ar.get("labels"))
            log("    flank %s %s: %d copies, allelic monophyly %d/%d"
                % (arm, rtag, len(units), n_mono, n_lab))
        res["regimes"][rtag] = rres
    return res


def quantiles(vals):
    if not vals:
        return None
    a = np.array(vals, dtype=float)
    return dict(n=int(a.size), min=float(a.min()), q25=float(np.percentile(a, 25)),
                median=float(np.median(a)), q75=float(np.percentile(a, 75)),
                max=float(a.max()), mean=float(a.mean()))


def md5of(text):
    import hashlib
    return hashlib.md5(text.encode("ascii", "replace")).hexdigest()


def md5file(path):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------

TSV_COLUMNS = ["record_type", "family", "community", "regime", "stratum", "level",
               "label", "n_copies", "n_haplotypes", "mono_point", "extra_clades_point",
               "frac_largest_same_label_clade", "n_boot", "boot_monophyly_frac",
               "boot_extra_clades_median", "perm_monophyly_rate", "n_pairs",
               "pdist_mean", "pdist_median", "d_paralog", "d_ortholog", "ratio",
               "pvalue", "note"]

LEVELS = ["arm", "haplotype", "array_arm", "array_haplotype"]


def emit_tsv(report, path):
    rows = []

    def row(**kw):
        rows.append({c: kw.get(c, "") for c in TSV_COLUMNS})

    def fmt(v, spec="%.6f"):
        return (spec % v) if isinstance(v, (int, float)) else ("" if v is None else v)

    for fam_label, fam_res in report["families"].items():
        comm = fam_res["community"]
        if "regimes" not in fam_res:
            row(record_type="family_error", family=fam_label, community=comm,
                note=str(fam_res.get("error", "")))
            continue
        row(record_type="reference", family=fam_label, community=comm,
            n_copies=fam_res["reference"]["length"],
            note=fam_res["reference"]["name"] + " : " +
                 str(fam_res["reference"]["params"].get("source")
                     or fam_res["reference"]["params"].get("note")))
        for rkey, rres in sorted(fam_res["regimes"].items()):
            regime = rres.get("fastga_filter", "")
            d = rres.get("discovery", {})
            row(record_type="discovery", family=fam_label, community=comm, regime=regime,
                n_copies=sum(rres.get("copies_per_arm", {}).values()),
                note=("alignments=%s deduped=%s projection_failures=%s cs=%s "
                      "windows=%s"
                      % (d.get("n_alignments"), d.get("n_deduped"),
                         d.get("n_projection_failures"), d.get("cs_convention"),
                         fam_res.get("target_windows", {}).get("n"))))
            for arm, n in sorted(rres.get("copies_per_arm", {}).items()):
                row(record_type="copy_count", family=fam_label, community=comm,
                    regime=regime, label=arm, n_copies=n)
            if rres.get("pairwise_crosscheck", {}).get("n"):
                xc = rres["pairwise_crosscheck"]
                row(record_type="validation", family=fam_label, community=comm,
                    regime=regime, n_pairs=xc["n"],
                    ratio=fmt(xc["spearman"], "%.4f"),
                    pvalue=fmt(xc["pearson"], "%.4f"),
                    note=("reference-frame p-distance vs direct FastGA p-distance: "
                          "spearman=ratio col, pearson=pvalue col, "
                          "median_abs_dev=%.4f max_abs_dev=%.4f"
                          % (xc["median_abs_dev"], xc["max_abs_dev"])))
            for sname, sres in sorted(rres.get("strata", {}).items()):
                if "arm" not in sres:
                    row(record_type="stratum_error", family=fam_label, community=comm,
                        regime=regime, stratum=sname, note=sres.get("error", ""))
                    continue
                for level in LEVELS:
                    lvl = sres.get(level)
                    if not lvl or "labels" not in lvl:
                        continue
                    for lab, st in sorted(lvl["labels"].items()):
                        row(record_type="label_stat", family=fam_label, community=comm,
                            regime=regime, stratum=sname, level=level, label=lab,
                            n_copies=st["n_copies"],
                            n_haplotypes=st["n_haplotypes"],
                            mono_point=st["mono_point"],
                            extra_clades_point=st["extra_clades_point"],
                            frac_largest_same_label_clade=fmt(
                                st["frac_largest_same_label_clade"], "%.4f"),
                            n_boot=N_BOOT,
                            boot_monophyly_frac=fmt(st["boot_monophyly_frac"], "%.4f"),
                            boot_extra_clades_median=fmt(
                                st["boot_extra_clades_median"], "%.2f"),
                            perm_monophyly_rate=fmt(st["perm_monophyly_rate"], "%.4f"))
                    labs = lvl["labels"]
                    n_lab = len(labs)
                    n_mono = sum(1 for v in labs.values() if v["mono_point"] == 1)
                    n_bmono = sum(1 for v in labs.values()
                                  if (v["boot_monophyly_frac"] or 0.0) >= 0.95)
                    row(record_type="verdict", family=fam_label, community=comm,
                        regime=regime, stratum=sname, level=level,
                        n_copies=sres.get("n_units"),
                        n_haplotypes=n_lab,
                        mono_point="%d/%d" % (n_mono, n_lab),
                        boot_monophyly_frac="%d/%d" % (n_bmono, n_lab),
                        note="labels monophyletic: point tree / bootstrap>=0.95")
                for pair, dv in (sres.get("arm_pair_distances") or {}).items():
                    row(record_type="arm_pair_distance", family=fam_label,
                        community=comm, regime=regime, stratum=sname,
                        level="unit_pair", label=pair, n_pairs=dv["n"],
                        pdist_median=fmt(dv["median"]),
                        pdist_mean=fmt(dv["mean"]),
                        note="n_capped=%d" % dv["n_capped"])
                for dc, dv in sres["distance_classes"].items():
                    row(record_type="distance_class", family=fam_label, community=comm,
                        regime=regime, stratum=sname, level="unit_pair", label=dc,
                        n_pairs=dv["n"], pdist_mean=fmt(dv["mean"]),
                        pdist_median=fmt(dv["median"]))
                nnu = sres.get("nearest_neighbour_unit") or {}
                if nnu.get("n"):
                    row(record_type="nearest_neighbour", family=fam_label,
                        community=comm, regime=regime, stratum=sname, level="unit",
                        n_pairs=nnu["n"],
                        frac_largest_same_label_clade=fmt(
                            nnu["frac_nearest_same_arm"], "%.3f"),
                        ratio=fmt(nnu["frac_nearest_same_arm_other_haplotype"], "%.3f"),
                        pvalue=fmt(nnu["frac_nearest_same_haplotype"], "%.3f"),
                        note=("frac_largest_same_label_clade = frac of copies whose "
                              "nearest neighbour is on the same arm; ratio = ... same arm "
                              "other haplotype; pvalue = ... same haplotype"))
                nna = sres.get("nearest_neighbour_array") or {}
                if nna.get("n"):
                    row(record_type="nearest_neighbour", family=fam_label,
                        community=comm, regime=regime, stratum=sname, level="array",
                        n_pairs=nna["n"],
                        frac_largest_same_label_clade=fmt(
                            nna["frac_nearest_same_arm"], "%.3f"),
                        ratio=fmt(nna["frac_nearest_same_arm_other_haplotype"], "%.3f"),
                        pvalue=fmt(nna["frac_nearest_same_haplotype"], "%.3f"),
                        note=("array level; frac_largest_same_label_clade = frac of "
                              "arrays whose nearest array is on the same arm; ratio = "
                              "...same arm other haplotype; pvalue = ...same haplotype"))
                pvo = sres.get("paralog_vs_ortholog") or {}
                if pvo.get("n"):
                    row(record_type="paralog_vs_ortholog", family=fam_label,
                        community=comm, regime=regime, stratum=sname,
                        level="copy", n_pairs=pvo["n"],
                        d_paralog=fmt(pvo["median_paralog"]),
                        d_ortholog=fmt(pvo["median_ortholog"]),
                        ratio=fmt(pvo["ratio_paralog_over_ortholog"], "%.3f"),
                        pvalue=fmt(pvo["wilcoxon_p_paralog_less"], "%.3g"),
                        note=("d_paralog = closest cross-arm copy in the SAME haplotype; "
                              "d_ortholog = median distance to the same arm in OTHER "
                              "haplotypes; ratio = d_paralog/d_ortholog; pvalue = "
                              "one-sided Wilcoxon on %d matched copies (H1: paralog "
                              "closer); frac with paralog farther = %.3f"
                              % (pvo["n"], pvo["frac_paralog_farther"])))
        fl = fam_res.get("flank_control", {})
        for rtag, rres in sorted(fl.get("regimes", {}).items()):
            regime = rtag[1:]
            for arm, ares in sorted(rres.get("arms", {}).items()):
                row(record_type="flank_control", family=fam_label, community=comm,
                    regime=regime, level="array_haplotype", label=arm,
                    n_copies=ares.get("n_copies", ""),
                    mono_point=("%s/%s" % (ares.get("n_allelic_pairs_monophyletic", ""),
                                           ares.get("n_individuals_labelled", ""))),
                    boot_monophyly_frac=fmt(ares.get("boot_monophyly_median"), "%.3f"),
                    note=("allelic (haplotype-pair) monophyly in 100 kb flank windows; "
                          "cross-arm alignments seen for this arm=%s; %s"
                          % (rres.get("cross_arm_alignments", {}).get(arm, ""),
                             ares.get("error", ""))))
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in TSV_COLUMNS) + "\n")
    log("wrote %s (%d rows)" % (path, len(rows)))
    return len(rows)


# --------------------------------------------------------------------------------------
# Self-tests
# --------------------------------------------------------------------------------------

def kmer_jaccard(x, y, k=21):
    comp = str.maketrans("ACGT", "TGCA")

    def canon(s):
        rc = s.translate(comp)[::-1]
        return s if s <= rc else rc

    kx = {canon(x[i:i + k]) for i in range(len(x) - k + 1)}
    ky = {canon(y[i:i + k]) for i in range(len(y) - k + 1)}
    return len(kx & ky) / float(len(kx | ky))


def run_self_tests():
    ok = 0
    tmp = os.path.join(WORKDIR or "/tmp", "selftest")
    os.makedirs(tmp, exist_ok=True)

    print("== self-test 1: NJ against Biopython and against an analytically known tree")
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    rng = random.Random(11)
    for trial in range(20):
        n = rng.randint(4, 9)
        pts = [(rng.random(), rng.random()) for _ in range(n)]
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                D[i, j] = math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])
        adj = nj_tree(D)
        mine = edge_clades(adj, n)
        dm = DistanceMatrix(names=["t%d" % i for i in range(n)],
                            matrix=[[float(D[i, j]) for j in range(i + 1)] for i in range(n)])
        bp = DistanceTreeConstructor().nj(dm)
        theirs = set()
        full = (1 << n) - 1
        for clade in bp.get_nonterminals():
            tips = [int(t.name[1:]) for t in clade.get_terminals()]
            m = sum(1 << i for i in tips)
            if 0 < m < full:
                theirs.add(m)
                theirs.add(full & ~m)
        # pendant-edge bipartitions are not emitted by Biopython's nonterminal clades;
        # compare the informative (internal) bipartitions only
        mine = {c for c in mine if 2 <= bin(c).count("1") <= n - 2}
        theirs = {c for c in theirs if 2 <= bin(c).count("1") <= n - 2}
        if mine != theirs:
            print("   trial %d: MISMATCH mine-theirs=%r theirs-mine=%r"
                  % (trial, mine - theirs, theirs - mine))
            raise AssertionError("NJ disagrees with Biopython")
        ok += 1
    print("   20/20 random Euclidean matrices: my NJ clade set == Biopython NJ clade set")

    print("== self-test 2: monophyly statistics on hand-built trees (adversarial)")
    # Unrooted tree (((0,1),(2,3)),(4,5)) with 6 tips:
    #   clades = {0,1}, {2,3}, {0,1,2,3}, {4,5} and their complements
    adj = {0: {6}, 1: {6}, 6: {0, 1, 8}, 8: {6, 7, 9}, 7: {2, 3, 8},
           2: {7}, 3: {7}, 4: {9}, 5: {9}, 9: {4, 5, 8}}
    def msk(*ids):
        return sum(1 << i for i in ids)
    cl = edge_clades(adj, 6)
    for want in ((0, 1), (2, 3), (0, 1, 2, 3), (4, 5)):
        assert msk(*want) in cl, want
    assert label_stats(cl, msk(0, 1), 6) == (1, 1, 1.0), label_stats(cl, msk(0, 1), 6)
    assert label_stats(cl, msk(0, 1, 2, 3), 6) == (1, 1, 1.0)
    # paraphyletic: {0,1,2} needs {0,1} + the lone tip 2 -> 2 extra clades
    m, e, f = label_stats(cl, msk(0, 1, 2), 6)
    assert m == 0 and e == 2, (m, e, f)
    assert abs(f - 2 / 3.0) < 1e-9, f
    # polyphyletic: {0,2} picks tips from two clades, no shared clade
    m, e, f = label_stats(cl, msk(0, 2), 6)
    assert m == 0 and e == 2, (m, e, f)
    # a caterpillar where an adjacent pair is a clade and a split pair is not
    adj2 = {0: {6}, 1: {6}, 6: {0, 1, 7}, 7: {6, 2, 8}, 8: {7, 3, 9}, 9: {8, 4, 5},
            2: {7}, 3: {8}, 4: {9}, 5: {9}}
    cl2 = edge_clades(adj2, 6)
    assert label_stats(cl2, msk(4, 5), 6)[0] == 1
    assert label_stats(cl2, msk(2, 3), 6)[0] == 0
    # label == all tips: uninformative
    assert label_stats(cl2, msk(*range(6)), 6)[0] is None
    # singleton label is trivially monophyletic
    m, e, f = label_stats(cl2, msk(0), 6)
    assert m == 1 and e == 1, (m, e, f)
    print("   monophyly / extra_clades correct for clade, paraphyly, polyphyly, caterpillar,"
          " all-tips, singleton")
    ok += 1

    print("== self-test 3: cs projection against a hand-computed alignment")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    rng = random.Random(3)
    ref = "".join(rng.choice("ACGT") for _ in range(3000))
    target = list(ref)
    # 3 substitutions, one 10 bp deletion in the target, one 5 bp insertion
    sub_positions = (100, 500, 2500)
    for pos in sub_positions:
        target[pos] = {"A": "C", "C": "G", "G": "T", "T": "A"}[target[pos]]
    del target[2000:2010]                     # target lacks 10 ref bases
    target[1200:1200] = list("ACGTA")         # target has 5 extra bases
    target = "".join(target)
    write_fasta(os.path.join(tmp, "q.fa"), [("REFX", ref)])
    write_fasta(os.path.join(tmp, "t.fa"), [("TGT", target)])
    paf = os.path.join(tmp, "x.paf")
    fastga(os.path.join(tmp, "q.fa"), os.path.join(tmp, "t.fa"), 512, paf)
    recs = parse_paf(paf)
    assert len(recs) == 1, [r["tname"] for r in recs]
    rec = recs[0]
    conv = None
    for c in CS_CONVENTIONS:
        try:
            proj = alignment_columns(rec, target, c, len(ref))
            conv = c
            break
        except ValueError:
            pass
    assert conv is not None, "cs could not be reconciled on a clean synthetic alignment"
    m = np.full(len(ref), 5, dtype=np.int8)
    for rc, code in zip(proj["refcols"], proj["codes"]):
        m[rc] = code
    n_gap = int((m == 4).sum())
    base_code = {"A": 0, "C": 1, "G": 2, "T": 3}
    comp_base = {"A": "T", "C": "G", "G": "C", "T": "A"}
    # the target's 10-base deletion must appear as gap columns for exactly those columns
    got_gap_cols = set(np.where(m == 4)[0].tolist())
    assert got_gap_cols == set(range(2000, 2010)), sorted(got_gap_cols)
    # every substituted reference column must carry the TARGET base, and only those
    diff = [i for i in range(len(ref))
            if m[i] in (0, 1, 2, 3) and m[i] != base_code[ref[i]]]
    assert diff == list(sub_positions), diff
    # the 5-base target insertion has no reference column and must not create one
    assert m.shape[0] == len(ref)
    print("   projection exact: gap columns == %r, differing columns == %r, convention=%s"
          % (sorted(got_gap_cols), diff, conv))
    ok += 1

    print("== self-test 4: p-distance, cap logic and bootstrap recovery")
    # reference of 1000 bp; 8 copies from 2 labels with a KNOWN split
    ref = "".join(rng.choice("ACGT") for _ in range(1000))

    def mutate(seq, positions):
        s = list(seq)
        for p in positions:
            s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
        return "".join(s)

    rng_p = [sorted(rng.sample(range(1000), 20)) for _ in range(8)]
    seqs = [mutate(ref, rng_p[0]), mutate(ref, rng_p[0]),   # label X, close pair
            mutate(ref, rng_p[1]), mutate(ref, rng_p[1]),   # label X
            mutate(ref, rng_p[2]), mutate(ref, rng_p[2]),   # label Y
            mutate(ref, rng_p[3]), mutate(ref, rng_p[3])]   # label Y
    # make X far from Y: give X an extra 60 substitutions
    seqs = [mutate(s, range(0, 120)) if i < 4 else s for i, s in enumerate(seqs)]
    units = [dict(arm=("X" if i < 4 else "Y"), sample="s%d" % i, hap="1",
                  seq=s, ident=1.0, tstart=0, tend=len(s)) for i, s in enumerate(seqs)]
    M = np.full((8, 1000), 5, dtype=np.int8)
    for i, u in enumerate(units):
        M[i, :] = ["ACGT".index(c) for c in u["seq"]]
    D, ncmp = pdists_from_columns(M, list(range(1000)))
    within = [D[0, 1], D[2, 3], D[4, 5], D[6, 7]]
    across = [D[i, j] for i in range(4) for j in range(4, 8)]
    assert max(within) < min(across), (within, across)
    res = analyse_tree_set(M, units, lambda u: u["arm"], n_boot=50, seed=7)
    assert res["labels"]["X"]["mono_point"] == 1 and res["labels"]["Y"]["mono_point"] == 1
    assert res["labels"]["X"]["boot_monophyly_frac"] >= 0.95, res["labels"]["X"]
    print("   known split recovered in %.2f of 50 replicates; within=%.4f across=%.4f"
          % (res["labels"]["X"]["boot_monophyly_frac"], max(within), min(across)))
    # adversarial: shuffled labels must NOT look monophyletic
    units2 = [dict(u, arm=("X", "Y", "X", "Y", "X", "Y", "X", "Y")[i])
              for i, u in enumerate(units)]
    res2 = analyse_tree_set(M, units2, lambda u: u["arm"], n_boot=20, seed=8)
    assert res2["labels"]["X"]["mono_point"] == 0, res2["labels"]["X"]
    print("   shuffled (adversarial) labels are not monophyletic; extra_clades=%d"
          % res2["labels"]["X"]["extra_clades_point"])
    # cap logic
    Msparse = np.full((4, 1000), 5, dtype=np.int8)
    Msparse[0, :300] = 0
    Msparse[1, :300] = 0
    Msparse[2, 500:800] = 1
    Msparse[3, 500:800] = 1
    Ds, ncmps = pdists_from_columns(Msparse, list(range(1000)))
    # no shared columns -> below MIN_COMPARABLE_COLS -> documented cap
    assert Ds[0, 2] == P_DIST_CAP and ncmps[0, 2] == 0, (Ds[0, 2], ncmps[0, 2])
    assert Ds[2, 3] == 0.0 and ncmps[2, 3] == 300, (Ds[2, 3], ncmps[2, 3])
    assert Ds[0, 1] == 0.0 and ncmps[0, 1] == 300
    # a single planted difference in a 300-column overlap gives p = 1/300
    Mdiff = Msparse.copy()
    Mdiff[1, 0] = 3
    Dd, _n = pdists_from_columns(Mdiff, list(range(1000)))
    # float32 matmul accumulation: agreement to ~1e-7 is expected and irrelevant
    assert abs(Dd[0, 1] - 1.0 / 300.0) < 1e-6, Dd[0, 1]
    print("   saturated pairs take the documented cap (%.1f) and are flagged" % P_DIST_CAP)
    ok += 1

    print("== self-test 5: dedupe logic (adversarial overlaps)")
    hits = [dict(tstart=0, tend=1000, nmatch=1000), dict(tstart=500, tend=1500, nmatch=900),
            dict(tstart=2000, tend=3000, nmatch=1000)]
    kept = []
    for h in sorted(hits, key=lambda x: (-x["nmatch"], x["tstart"])):
        if any(min(h["tend"], k["tend"]) - max(h["tstart"], k["tstart"])
               > DEDUPE_MAX_OVERLAP * min(h["tend"] - h["tstart"], k["tend"] - k["tstart"])
               for k in kept):
            continue
        kept.append(h)
    assert len(kept) == 2 and kept[1]["tstart"] == 2000, kept
    print("   overlapping hit removed, non-overlapping hit kept")
    ok += 1

    print("== self-test 6: acrocentric reference picker on a synthetic tandem array")
    unit = "".join(rng.choice("ACGT") for _ in range(ACRO_TILE))
    fl1, fl2 = 6000, 48000
    arr = ("".join(rng.choice("ACGT") for _ in range(fl1)) + (unit * 8) +
           "".join(rng.choice("ACGT") for _ in range(fl2)))
    d = os.path.join(tmp, "acropick_fa")
    os.makedirs(d, exist_ok=True)
    hp = os.path.join(tmp, "acro_win.fa")
    write_fasta(hp, [("CHM13#0#chr13:2544-502543_chr13_parm", arr)])
    hp2 = os.path.join(tmp, "acro_win_copy.fa")
    write_fasta(hp2, [("CHM13#0#chr13:2544-502543_chr13_parm", arr)])

    class _FakeSub:
        path = "synthetic"
        by_arm = {"chr13_parm": [("CHM13", "0",
                                  "CHM13#0#chr13:2544-502543_chr13_parm", len(arr))]}

        def extract(self, names, out_fa):
            shutil.copy(hp, out_fa)
            return out_fa
    info = pick_acro_reference(_FakeSub(), ["chr13_parm"], tmp)
    st, tl, nh, nt = info["start"], info["tile"], info["n_homologs"], info["n_tiles"]
    assert nh >= 5, (nh, info)
    assert info["median_homologs"] < nh, info
    assert fl1 - ACRO_TILE <= st <= fl1 + 8 * ACRO_TILE, (st, nh)
    print("   tandem array located: tile at %d with %d homologs (median %s over %d tiles)"
          % (st, nh, info["median_homologs"], nt))
    ok += 1

    print("== self-test 7: inverted copies must land in the reference orientation")
    comp = str.maketrans("ACGT", "TGCA")
    ref = "".join(rng.choice("ACGT") for _ in range(3000))   # own reference
    fwd = ref
    inv = ref.translate(comp)[::-1]            # reverse complement: a '-' strand copy
    write_fasta(os.path.join(tmp, "q7.fa"), [("REFX", ref)])
    # FastGA emits nothing when the target record is byte-identical to the query, so both
    # copies carry a few substitutions (the point of the test is the orientation, not the
    # identity)
    def _sub(seq_, positions):
        out = list(seq_)
        for pp in positions:
            out[pp] = {"A": "C", "C": "G", "G": "T", "T": "A"}[out[pp]]
        return "".join(out)

    sub_at = (250, 900, 1700, 2400)
    fwd_m = _sub(fwd, sub_at)
    inv_m = _sub(inv, sub_at)[::-1].translate(comp) if False else None
    # mutate the same reference positions, then take the reverse complement
    inv_m = _sub(ref, sub_at).translate(comp)[::-1]
    write_fasta(os.path.join(tmp, "inv_t.fa"), [("FWD", fwd_m), ("INV", inv_m)])
    paf_inv = os.path.join(tmp, "inv.paf")
    fastga(os.path.join(tmp, "q7.fa"), os.path.join(tmp, "inv_t.fa"), 512, paf_inv)
    recs_inv = parse_paf(paf_inv)
    strands = {}
    rows = {}
    for rec in recs_inv:
        for c in CS_CONVENTIONS:
            try:
                proj = alignment_columns(rec,
                                         inv_m if rec["tname"] == "INV" else fwd_m, c,
                                         len(ref))
                break
            except ValueError:
                proj = None
        assert proj is not None
        strands[rec["tname"]] = rec["strand"]
        m = np.full(len(ref), 5, dtype=np.int8)
        for rc, code in zip(proj["refcols"], proj["codes"]):
            m[rc] = code
        rows[rec["tname"]] = m
    assert strands.get("INV") == "-", strands
    assert strands.get("FWD") == "+", strands
    both = (rows["FWD"] >= 0) & (rows["FWD"] <= 3) & (rows["INV"] >= 0) & (rows["INV"] <= 3)
    assert both.sum() > 0.9 * len(ref), both.sum()
    assert np.all(rows["FWD"][both] == rows["INV"][both]), \
        "an inverted copy was not projected into the reference orientation"
    print("   reverse-complement copy (strand %s) projects identically to the forward copy "
          "over %d columns" % (strands["INV"], int(both.sum())))
    ok += 1

    print("\nALL SELF-TESTS PASSED (%d groups)" % ok)


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------

def main():
    global WORKDIR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--workscratch", default="/scratch/erikg/m1c")
    ap.add_argument("--out-tsv", default=os.path.join(REPO, "scripts", "cladistics",
                                                     "duplicon_monophyly.tsv"))
    ap.add_argument("--out-json", default=os.path.join(REPO, "scripts", "cladistics",
                                                      "duplicon_monophyly.json"))
    ap.add_argument("--only-family", default=None)
    args = ap.parse_args()

    WORKDIR = args.workscratch
    os.makedirs(WORKDIR, exist_ok=True)

    if args.self_test:
        run_self_tests()
        return

    for path in (COPIES_FA, COPIES_FA + ".fai", FLANK_FA, FLANK_FA + ".fai",
                 SUPERPOP_TSV, COMMUNITIES, GFF3):
        if not os.path.exists(path):
            die("missing required input: %s" % path)

    t_start = time.time()
    report = dict(
        title="Duplicon monophyly test (Referee #1 major point 1)",
        generated=time.strftime("%Y-%m-%dT%H:%M:%S"),
        host=os.uname().nodename,
        cwd=os.getcwd(),
        inputs={},
        parameters=dict(
            REGIMES=REGIMES, MIN_ALN=MIN_ALN, FASTGA_MIN_ID=FASTGA_MIN_ID,
            MIN_UNIT_ID=MIN_UNIT_ID, DEDUPE_MAX_OVERLAP=DEDUPE_MAX_OVERLAP,
            MAX_COPIES_PER_WINDOW=MAX_COPIES_PER_WINDOW, BLOCK=BLOCK,
            MIN_BLOCK_COMPLETENESS=MIN_BLOCK_COMPLETENESS,
            HIGHCONF_UNIT_ID=HIGHCONF_UNIT_ID,
            MIN_COMPARABLE_COLS=MIN_COMPARABLE_COLS, P_DIST_CAP=P_DIST_CAP,
            N_BOOT=N_BOOT, N_PERM=N_PERM, SEED=SEED, MAX_HAPS=MAX_HAPS,
            MIN_COMMON_HAPS=MIN_COMMON_HAPS, MAX_HAPS_PER_ARM=MAX_HAPS_PER_ARM,
            REF_MAX=REF_MAX, FASTGA_THREADS=FASTGA_THREADS),
        tools=dict(fastga=tool_version([FASTGA]),
                   samtools=tool_version(["samtools", "--version"]),
                   python=sys.version.split()[0], numpy=np.__version__),
        families={})

    for key, path in (("copies_fasta", COPIES_FA), ("flank_fasta", FLANK_FA),
                      ("superpop_tsv", SUPERPOP_TSV), ("gff3", GFF3),
                      ("communities", COMMUNITIES), ("annotations_bed", ANNOT_BED)):
        report["inputs"][key] = dict(path=path, exists=os.path.exists(path))
    report["inputs"]["copies_fasta"]["n_windows"] = len(read_fai(COPIES_FA))
    report["inputs"]["flank_fasta"]["n_windows"] = len(read_fai(FLANK_FA))

    sup = load_superpops()
    comms = load_communities()
    copies_sub = Substrate(COPIES_FA, "copies")
    flank_sub = Substrate(FLANK_FA, "flanks")
    report["inputs"]["copies_fasta"]["arms"] = {
        k: len(v) for k, v in sorted(copies_sub.by_arm.items())}
    log("loaded substrates: %d copy windows, %d flank windows, %d labelled haplotypes"
        % (len(copies_sub.records), len(flank_sub.records), len(sup)))

    for fam in FAMILIES:
        if args.only_family and fam["label"] != args.only_family:
            continue
        fout = dict(label=fam["label"], community=fam["comm"], arms=fam["arms"],
                    community_of_arms={a.replace("_parm", "_p").replace("_qarm", "_q"):
                                       comms.get(a.replace("_parm", "_p")
                                                 .replace("_qarm", "_q"))
                                       for a in fam["arms"]})
        try:
            family_run(fam, copies_sub, flank_sub, sup, fout)
        except SystemExit:
            raise
        except Exception as exc:                              # noqa: BLE001
            import traceback
            fout["error"] = "%s: %s" % (type(exc).__name__, exc)
            fout["traceback"] = traceback.format_exc()
            log("FAMILY %s FAILED: %s" % (fam["label"], exc))
        report["families"][fam["label"]] = fout

    report["wall_seconds_total"] = round(time.time() - t_start, 1)
    report["max_rss_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1)
    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=1, sort_keys=True, default=str)
    log("wrote %s" % args.out_json)
    emit_tsv(report, args.out_tsv)
    log("total wall %.1f s" % report["wall_seconds_total"])


if __name__ == "__main__":
    main()
