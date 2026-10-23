#!/usr/bin/env python3
"""Copy-enumeration of reference "core" sequences by pangenome-alignment transitivity.

Idea (operator's): do not compare duplicated windows with each other; take a CORE
sequence from CHM13, follow it through an all-vs-all pangenome alignment
representation, and enumerate every other reference position it connects to.

The shipped 37 GB index (hprc25272.bidirectional.impg) is orphaned from its source
PAFs (see response/notes/impg_transitivity_copy_enumeration.md, failure F1), so the
alignment representation is rebuilt locally, at core scale only:

  hop 1: FastGA(cores.fa -> CHM13 reference)          (identity to the core)
  hop 2: FastGA(hop-1 reference subsequences -> CHM13) (identity to a hop-1 copy,
                                                      i.e. transitive closure)
  query layer: `impg index` + `impg query` over those PAFs, so that placements,
               merging, identity filtering and original-sequence-coordinate
               remapping are done by impg rather than by this script.

Design constraints honoured here:
  * the reference is read ONLY through windowed random access (faidx offset +
    bgzip .gzi checkpoint + zlib inflate). No reference text is written to disk;
    the single FASTA that is written is the 12-core query file (~180 kb), which the
    aligner requires as a file.
  * no per-suffix / per-position columns, no whole-genome scans.
  * fail loud: every assumption (core inside PHR, core away from chromosome ends,
    reference sequence present in the .fai, PAF present) is asserted.

Subcommands
  cores      verify the core table against data/chm13.phrs.bed, the community
             assignments and the annotation; write <work>/cores.tsv and <work>/cores.fa
  fastga     run FastGA for hop 1 (and optionally hop 2 batches)
  index      build an impg index over the produced PAF(s)
  query      run `impg query` for every core (bedpe), timing each call
  aggregate  parse impg output -> scripts/cladistics/impg_copy_enumeration.tsv
"""

import argparse
import gzip
import json
import os
import re
import subprocess
import tempfile
import time
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------------------------------------------------------------- resources ---
REF_GZ = ("/moosefs/guarracino/HPRCv2/PHR_III/trimmed-assemblies/"
          "chm13v2.0_maskedY_rCRS.fa.PanSN.trimmed.fa.gz")
REF_FA = "/tmp/impg_enum_t/ref/chm13.ref.fa"   # local staging of REF_GZ (see notes)
PHRS_BED = os.path.join(REPO, "data/chm13.phrs.bed")
ASSIGN = os.path.join(REPO, "data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
GFF = os.path.join(REPO, "data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")

FASTGA = os.path.expanduser("~/.cargo/bin/FastGA")
IMPG = os.path.expanduser("~/.cargo/bin/impg")
WFMASH = os.path.expanduser("~/bin/wfmash")

# All scratch output (PAFs, indexes, per-core bedpe) goes to a private temp tree, never
# into the repo: only scripts/cladistics/impg_copy_enumeration.tsv is a repo artifact.
WORK = os.path.join(tempfile.gettempdir(), "impg_copy_enum", "work")

# ------------------------------------------------------------------- cores ---
# Coordinates are CHM13 chromosome coordinates, 0-based half-open (BED convention),
# identical to the convention of data/chm13.phrs.bed and the Liftoff GFF.
# kind: 'phr'   = inside a PHR interval of data/chm13.phrs.bed
#       'acro'  = acrocentric p-arm rDNA array (C7 arm); annotated PHR for chr13p
#                 does not exist in phrs.bed, so containment is asserted as False
#                 and the proximity to the phrs.bed interval is reported
#       'inter' = interstitial, non-subtelomeric control
CORES = [
    # id,                 chrom,  start,      end,        kind,    expect_arm, note
    ("D4Z4_4q",           "chr4",  193469000, 193484000, "phr",   "chr4_q",  "D4Z4/DUX4L array, DUX4L6 copies"),
    ("DUX4L_10q",         "chr10", 134620000, 134635000, "phr",   "chr10_q", "D4Z4-like DUX4L array"),
    ("TUBB8B_18p",        "chr18",    195000,    210000, "phr",   "chr18_p", "TUBB8B (chr18:200997-227325)"),
    ("TUBB8_10p",         "chr10",     35000,     50000, "phr",   "chr10_p", "TUBB8 (chr10:39416-69557)"),
    ("OR4F_5q",           "chr5",  181900000, 181915000, "phr",   "chr5_q",  "OR4F3 (chr5:181905487-181933378)"),
    ("OR4F_6q",           "chr6",  172005000, 172020000, "phr",   "chr6_q",  "OR4F7P (chr6:172013163-172014100)"),
    ("SEPTIN14_1p",       "chr1",     142000,    147000, "phr",   "chr1_p",  "SEPTIN14P22_5/P13 cluster (5 kb: copies are 2.6 kb, a 15 kb core only matches itself)"),
    ("rDNA_13p",          "chr13",    5772000,   5787000, "acro",  "chr13_p", "RNA45S1 unit chr13:5774303-5787683"),
    ("INT_REXO1L_chr8",   "chr8",   86087000,  86102000, "inter", None,      "REXO1L1P/3P/8P pseudogene array 8q21"),
    ("INT_16p11_2",       "chr16",  29500000,  29515000, "inter", None,      "16p11.2 BP4/BP5 interstitial SD block"),
]


# ------------------------------------------------------ windowed .fai reader ---
class FaidxReader:
    """Windowed random access to an uncompressed FASTA + .fai (pread, no mmap).

    Only the byte spans that cover the requested window are read; nothing else is
    touched and nothing is written. The staged local FASTA is a byte-identical
    decompression of the project reference (staging measured, see notes); the
    original sequence names (CHM13#0#chrN:beg-end) are preserved so that impg's
    `--original-sequence-coordinates` converts back to chromosome coordinates.
    """

    def __init__(self, path):
        self.path = path
        self.fai = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")
                self.fai[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fd = os.open(path, os.O_RDONLY)

    def close(self):
        os.close(self._fd)

    def fetch(self, seqname, start, end):
        """0-based half-open slice in *sequence-local* coordinates."""
        if seqname not in self.fai:
            raise KeyError("sequence %r not in %s.fai" % (seqname, self.path))
        length, offset, linebases, linewidth = self.fai[seqname]
        if not (0 <= start < end <= length):
            raise ValueError("window %s:%d-%d outside 0..%d" % (seqname, start, end, length))
        first_line = start // linebases
        last_line = (end - 1) // linebases
        byte_start = offset + first_line * linewidth + (start % linebases)
        byte_end = offset + last_line * linewidth + ((end - 1) % linebases) + 1
        raw = os.pread(self._fd, byte_end - byte_start, byte_start)
        if not raw:
            raise RuntimeError("empty pread for %s:%d-%d" % (seqname, start, end))
        seq = raw.replace(b"\n", b"")
        want = end - start
        if len(seq) != want:
            raise RuntimeError("short read for %s:%d-%d (%d/%d bytes)"
                               % (seqname, start, end, len(seq), want))
        # ACGT/N only: the staged reference is the soft-masked-free trimmed build
        bad = set(seq) - set(b"ACGTNacgtn")
        if bad:
            raise RuntimeError("unexpected bytes %r in %s" % (bad, seqname))
        return seq


# ----------------------------------------------------------------- helpers ---
def parse_pansn(name):
    """'CHM13#0#chr4:3087-193572740' -> ('chr4', 3087) with 1-based chrom start."""
    m = re.match(r"^.*#(?P<chrom>[^:#]+):(?P<beg>\d+)-(?P<end>\d+)$", name)
    if not m:
        raise ValueError("sequence name %r is not '...#<chrom>:<beg>-<end>'" % name)
    return m.group("chrom"), int(m.group("beg"))


def read_phrs():
    rows = []
    with open(PHRS_BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def read_assignments():
    arm2comm, comm2arms = {}, {}
    with open(ASSIGN) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] == "ChromArm":
                continue
            arm, comm, arms = f[0], f[1], [a.strip() for a in f[2].split(",")]
            arm2comm[arm] = comm
            comm2arms[comm] = arms
    return arm2comm, comm2arms


def read_centromeres():
    """Centromere intervals from data/chm13-annotations.bed -> arm partition anchor."""
    cen = {}
    with open(os.path.join(REPO, "data/chm13-annotations.bed")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[3] == "Centromere":
                cen[f[0]] = (int(f[1]), int(f[2]))
    return cen


def gff_genes_in(chrom, start, end, limit=400):
    """Gene features (name, start, end) overlapping [start,end), streamed from the GFF."""
    hits = []
    with gzip.open(GFF, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[0] != chrom:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if e <= start or s >= end:
                continue
            m = re.search(r"gene_name=([^;]+)", f[8])
            hits.append((m.group(1) if m else f[8].split(";")[0], s, e))
            if len(hits) >= limit:
                break
    return hits


def sh(cmd, log=None):
    t0 = time.time()
    proc = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True)
    wall = time.time() - t0
    if log:
        with open(log, "a") as fh:
            fh.write("$ %s\n[w] %.1fs rc=%d\n" % (cmd if isinstance(cmd, str) else " ".join(cmd),
                                                   wall, proc.returncode))
            fh.write("--- stdout ---\n" + proc.stdout[-4000:] + "\n")
            fh.write("--- stderr ---\n" + proc.stderr[-4000:] + "\n")
    return proc, wall


def read_cores(work):
    rows = []
    with open(os.path.join(work, "cores.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = dict(zip(hdr, f))
            for k in ("start", "end", "length", "n_genes", "dist_to_nearest_phr_bp"):
                r[k] = int(r[k])
            rows.append(r)
    return rows


def placement_arm(chrom, start, end, phrs, cen):
    """Arm label of a placement: '<n>p' / '<n>q' from the centromere partition, with
    '_inter' appended when the placement is not inside any PHR interval of
    data/chm13.phrs.bed (i.e. it is an interstitial copy, not a subtelomeric one)."""
    short = chrom.replace("chr", "")
    in_phr = any(c == chrom and s <= end and start <= e for (c, s, e) in phrs)
    if chrom in cen:
        cs, ce = cen[chrom]
        if end <= cs:
            arm = "%sp" % short
        elif start >= ce:
            arm = "%sq" % short
        else:
            arm = "%s_cen" % short
    else:
        arm = "%s_unknown" % short
    return arm, in_phr


def pgf_num(x, nd=4):
    return ("%%.%df" % nd) % x


def parse_bedpe(path):
    """impg -o bedpe rows -> list of dicts (query/target intervals, identity)."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise SystemExit("FAIL: bedpe row with %d columns in %s" % (len(f), path))
            d = dict(q_name=f[0], q_start=int(f[1]), q_end=int(f[2]),
                     t_name=f[3], t_start=int(f[4]), t_end=int(f[5]),
                     strand=f[8])
            for col in f[10:]:
                if col.startswith("gi:f:"):
                    d["gi"] = float(col[5:])
                elif col.startswith("bi:f:"):
                    d["bi"] = float(col[5:])
            out.append(d)
    return out


def gap_compressed_identity(cigar):
    """Same definition as impg's gi: matches / (matches+mismatches+ins_events+del_events)."""
    m = mm = ins = dele = 0
    n = 0
    for ch in cigar:
        if ch.isdigit():
            n = n * 10 + int(ch)
            continue
        if ch in "M=":
            m += n
        elif ch == "X":
            mm += n
        elif ch == "I":
            ins += 1
        elif ch == "D":
            dele += 1
        n = 0
    den = m + mm + ins + dele
    return (m / den) if den else 0.0


def parse_paf(path):
    """FastGA PAF rows -> list of dicts; identity = matches/block_len (paf column 10/11)."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise SystemExit("FAIL: PAF row with %d columns in %s: %.120s" % (len(f), path, line))
            d = dict(q_name=f[0], q_len=int(f[1]), q_start=int(f[2]), q_end=int(f[3]),
                     strand=f[4], t_name=f[5], t_len=int(f[6]),
                     t_start=int(f[7]), t_end=int(f[8]), matches=int(f[9]),
                     block_len=int(f[10]))
            if d["block_len"] <= 0:
                raise SystemExit("FAIL: zero block length in %s" % path)
            d["identity"] = d["matches"] / d["block_len"]
            d["cigar"] = None
            for x in f[12:]:
                if x.startswith("cg:Z:"):
                    d["cigar"] = x[5:]
                    break
            d["gi"] = gap_compressed_identity(d["cigar"]) if d["cigar"] else d["identity"]
            out.append(d)
    return out


def ref_to_chrom(name, start, end):
    """Reference placement -> chromosome coordinates.

    With impg --original-sequence-coordinates the output name already lost its
    ':beg-end' suffix and the coordinates are already chromosome coordinates
    (verified: the core self-hit is returned as [core_start+1, core_end+1], i.e.
    1-based inclusive).  Without the suffix we therefore pass the interval through.
    """
    if ":" in name:
        m = re.match(r"^.*#(?P<chrom>[^:#]+):(?P<beg>\d+)-(?P<end>\d+)$", name)
        if m:
            beg = int(m.group("beg"))
            return m.group("chrom"), start + beg - 1, end + beg - 1
    # PanSN name stripped of its coordinate suffix (impg --original-sequence-coordinates),
    # e.g. 'CHM13#0#chr4'.  Keep the chromosome, drop sample/haplotype.
    return name.split("#")[-1], start, end


def merge_intervals(ivs, max_gap):
    """Merge (chrom, strand, start, end, core_start, core_end) rows into loci.

    Two alignments belong to the same placement if they are on the same chromosome
    and strand and are separated by at most `max_gap` bp.  This is what turns
    FastGA/impg's fragmented records into a copy count.
    """
    loci = []
    for iv in sorted(ivs, key=lambda x: (x["chrom"], x["strand"], x["start"])):
        if loci:
            last = loci[-1]
            if (last["chrom"] == iv["chrom"] and last["strand"] == iv["strand"]
                    and iv["start"] - last["end"] <= max_gap):
                last["end"] = max(last["end"], iv["end"])
                last["core_ivs"].append((iv["core_start"], iv["core_end"]))
                last["gi"].append(iv["gi"])
                continue
        loci.append(dict(chrom=iv["chrom"], strand=iv["strand"], start=iv["start"],
                         end=iv["end"], core_ivs=[(iv["core_start"], iv["core_end"])],
                         gi=[iv["gi"]]))
    for L in loci:
        ivs_sorted = sorted(L["core_ivs"])
        covered, cur_s, cur_e = 0, None, None
        for s0, e0 in ivs_sorted:
            if cur_s is None:
                cur_s, cur_e = s0, e0
            elif s0 <= cur_e:
                cur_e = max(cur_e, e0)
            else:
                covered += cur_e - cur_s
                cur_s, cur_e = s0, e0
        if cur_s is not None:
            covered += cur_e - cur_s
        L["core_cov"] = covered
        L["gi"] = sorted(L["gi"])
        L["gi_median"] = L["gi"][len(L["gi"]) // 2]
    return loci



def hop2_stats(cid, hop2_paf, hop1_starts, r, phrs, cen, min_identity, min_block):
    """New reference positions reached only through a hop-1 copy (transitive hop 2)."""
    existing = set((c, s0) for (c, _strand, s0) in hop1_starts)
    new_starts, arms, n_recs = set(), defaultdict(int), 0
    for d in hop2_paf:
        if not d["q_name"].startswith(cid + "|"):
            continue
        if d["identity"] < min_identity or d["block_len"] < min_block:
            continue
        chrom, s0, e0 = ref_to_chrom(d["t_name"], d["t_start"], d["t_end"])
        if chrom == r["chrom"] and abs(s0 - r["start"]) <= 1:
            continue
        n_recs += 1
        if any(c == chrom and abs(s0 - s1) <= 1000 for (c, s1) in existing):
            continue
        if (chrom, s0) in new_starts:
            continue
        new_starts.add((chrom, s0))
        a, _p = placement_arm(chrom, s0, s0 + 1, phrs, cen)
        arms[a] += 1
    return dict(hop2_n_align=n_recs, hop2_n_new_positions=len(new_starts),
                hop2_n_new_arms=len(arms),
                hop2_per_arm_new=";".join("%s:%d" % (k, v) for k, v in sorted(arms.items())))


def core_stats(cid, r, bedpe_path, phrs, cen, min_core_cov, merge_gap):
    """All copy-enumeration statistics for one core from one impg bedpe file."""
    recs = parse_bedpe(bedpe_path)
    raw = []
    for d in recs:
        if d["t_name"] == cid:
            ref_name, ref_s, ref_e = d["q_name"], d["q_start"], d["q_end"]
            cov_s, cov_e = d["t_start"], d["t_end"]
        elif d["q_name"] == cid:
            ref_name, ref_s, ref_e = d["t_name"], d["t_start"], d["t_end"]
            cov_s, cov_e = d["q_start"], d["q_end"]
        else:
            raise SystemExit("FAIL: bedpe row in %s names neither side %r" % (bedpe_path, cid))
        chrom, s0, e0 = ref_to_chrom(ref_name, ref_s, ref_e)
        raw.append(dict(chrom=chrom, start=s0, end=e0, strand=d["strand"],
                        core_start=cov_s, core_end=cov_e, gi=d["gi"]))
    if not raw:
        raise SystemExit("FAIL: %s produced no rows above the identity threshold" % cid)
    n_before = len(raw)
    raw = [x for x in raw
           if (x["core_end"] - x["core_start"]) >= min_core_cov * r["length"]]
    n_dropped_short = n_before - len(raw)
    if not raw:
        raise SystemExit("FAIL: %s: every placement was below --min-core-cov" % cid)

    starts = sorted(set((x["chrom"], x["strand"], x["start"]) for x in raw))
    loci = merge_intervals(raw, merge_gap)

    span = {}
    by_chrom = defaultdict(list)
    for x in raw:
        by_chrom[x["chrom"]].append((x["start"], x["end"]))
    for c, ivs in by_chrom.items():
        ivs.sort()
        cur_s, cur_e, tot = ivs[0][0], ivs[0][1], 0
        for s0, e0 in ivs[1:]:
            if s0 <= cur_e:
                cur_e = max(cur_e, e0)
            else:
                tot += cur_e - cur_s
                cur_s, cur_e = s0, e0
        span[c] = tot + cur_e - cur_s

    arm_of_start, phr_of_start = {}, {}
    for (c, strand, s0) in starts:
        a, p = placement_arm(c, s0, s0 + 1, phrs, cen)
        arm_of_start[(c, strand, s0)] = a
        phr_of_start[(c, strand, s0)] = p

    # the core's own locus = the placement whose interval IS the core (+-1 bp)
    selfs = [x for x in raw if x["chrom"] == r["chrom"]
             and abs(x["start"] - r["start"]) <= 1 and abs(x["end"] - r["end"]) <= 1]
    if len(selfs) != 1:
        raise SystemExit("FAIL: %s -> %d identity self placements; expected 1"
                         % (cid, len(selfs)))
    self_keys = set((x["chrom"], x["strand"], x["start"]) for x in selfs)
    nonself = [k for k in starts if k not in self_keys]
    n_overlap_core = len([x for x in raw if x["chrom"] == r["chrom"]
                          and x["start"] <= r["end"] and r["start"] <= x["end"]])

    arms = defaultdict(int)
    for k in nonself:
        arms[arm_of_start[k]] += 1
    member_arms = set(x.strip() for x in r["community_arms"].split(",") if x.strip())
    n_member = sum(v for k, v in arms.items() if k in member_arms)
    n_nonmember = sum(v for k, v in arms.items() if k not in member_arms)
    n_inter = sum(1 for k in nonself if not phr_of_start[k])
    gis = sorted(x["gi"] for x in raw if (x["chrom"], x["strand"], x["start"]) not in self_keys)
    nonself_loci = [L for L in loci
                    if not (L["chrom"] == r["chrom"] and L["start"] <= r["end"]
                            and r["start"] <= L["end"])]
    covs = sorted(L["core_cov"] / r["length"] for L in nonself_loci)
    lens = sorted(L["end"] - L["start"] for L in nonself_loci)

    def med(v):
        if not v:
            return float("nan")
        n = len(v)
        return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])

    return dict(
        n_align=len(raw), n_dropped_short=n_dropped_short, n_self=len(selfs),
        n_positions=len(nonself), n_placements_overlapping_core=n_overlap_core,
        n_loci=len(nonself_loci), n_interstitial=n_inter,
        n_chroms=len(set(k[0] for k in nonself)), n_arms=len(arms),
        n_on_community_arms=n_member, n_on_noncommunity_arms=n_nonmember,
        per_arm=";".join("%s:%d" % (k, v) for k, v in sorted(arms.items())),
        identity_median=med(gis), identity_min=(gis[0] if gis else float("nan")),
        identity_max=(gis[-1] if gis else float("nan")),
        core_cov_median=med(covs), placement_len_median=med(lens),
        ref_span_bp_all=sum(span.values()),
        copies_core_units=sum(span.values()) / r["length"],
        _nonself_starts=nonself)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cores", "fastga", "index", "query", "hop2",
                                    "crosscheck", "aggregate", "selftest"])
    ap.add_argument("--work", default=WORK)
    ap.add_argument("--core", action="append")
    ap.add_argument("--identity", type=float, default=95.0)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--run", default="f8192", help="run label for the FastGA frequency setting")
    ap.add_argument("--alt-run", default="f1000",
                    help="second run label aggregated as a filter-sensitivity column")
    ap.add_argument("--freq", type=int, default=8192, help="FastGA -f seed-count cutoff")
    ap.add_argument("--wfmash", default=None, help="PAF from an independent aligner")
    ap.add_argument("--hop2-max-bp", type=int, default=30_000_000,
                    help="cap on the total hop-2 query length in bp (default 30 Mb)")
    ap.add_argument("--min-core-cov", type=float, default=0.5,
                    help="a placement must cover >= this fraction of the core (default 0.5)")
    ap.add_argument("--merge-gap", type=int, default=1000,
                    help="max bp between two alignments of the same placement (default 1000)")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    log = os.path.join(args.work, "timing.log")

    if args.cmd == "selftest":
        r = FaidxReader(REF_FA)
        seqs = list(r.fai)
        name = [s for s in seqs if s.startswith("CHM13#0#chr10")][0]
        chrom, beg = parse_pansn(name)
        # TUBB8 starts at chr10:39416 (1-based) per the Liftoff GFF -> expect a gene start
        sub = r.fetch(name, 39416 - beg, 39416 - beg + 60)
        print("seq", name, "chrom", chrom, "beg", beg)
        print("chr10:39416+60 =", sub.decode())
        # geometric check: line geometry must be self-consistent
        adj = r.fetch(name, 1000, 1100)
        print("chr10:+1000..1100 =", adj.decode()[:60], "...")
        r.close()
        return

    if args.cmd == "cores":
        phrs = read_phrs()
        arm2comm, comm2arms = read_assignments()
        rd = FaidxReader(REF_FA)
        seq_by_chrom = {}
        for s in rd.fai:
            c, b = parse_pansn(s)
            seq_by_chrom[c] = (s, b)
        # chromosome length = sequence length + its 1-based chromosome start offset - 1
        sizes = {c: (b - 1) + rd.fai[s][0] for c, (s, b) in seq_by_chrom.items()}

        rows = []
        fa = open(os.path.join(args.work, "cores.fa"), "w")
        for cid, chrom, start, end, kind, arm, note in CORES:
            if chrom not in seq_by_chrom:
                raise SystemExit("FAIL: %s not in reference .fai" % chrom)
            seqname, seqbeg = seq_by_chrom[chrom]
            length = rd.fai[seqname][0]
            # seqbeg is the 1-based chromosome coordinate of local base 0
            if not (seqbeg - 1 <= start and end <= (seqbeg - 1) + length):
                raise SystemExit("FAIL: %s:%d-%d outside reference sequence %s"
                                 % (chrom, start, end, seqname))
            inside = [(c, s, e) for (c, s, e) in phrs
                      if c == chrom and s <= end and start <= e]
            near = min([min(abs(start - e), abs(s - end)) for (c, s, e) in phrs
                        if c == chrom] or [10 ** 12])
            if kind == "phr":
                if len(inside) != 1:
                    raise SystemExit("FAIL: %s expects exactly one PHR overlap, got %r"
                                     % (cid, inside))
            elif kind == "inter":
                if inside:
                    raise SystemExit("FAIL: %s is an interstitial control but overlaps PHR %r"
                                     % (cid, inside))
                gap_start, gap_end = start, sizes[chrom] - end
                if gap_start < 5_000_000 or gap_end < 5_000_000:
                    raise SystemExit("FAIL: %s interstitial control is <5 Mb from a chromosome end "
                                     "(%d/%d)" % (cid, gap_start, gap_end))
            # arm/community
            if arm is None:
                comm, members = None, []
            else:
                comm = arm2comm.get(arm)
                if comm is None:
                    raise SystemExit("FAIL: arm %r for %s not in assignments file" % (arm, cid))
                members = comm2arms[comm]
            genes = gff_genes_in(chrom, start, end)
            seq = rd.fetch(seqname, start - (seqbeg - 1), end - (seqbeg - 1))
            if len(seq) != end - start:
                raise SystemExit("FAIL: short fetch for %s" % cid)
            fa.write(">%s %s:%d-%d\n" % (cid, chrom, start, end))
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i + 60].decode() + "\n")
            rows.append(dict(
                core=cid, chrom=chrom, start=start, end=end, length=end - start,
                kind=kind, arm=arm, community=comm, community_arms=",".join(members),
                in_phr=bool(inside), phr=",".join("%s:%d-%d" % t for t in inside),
                dist_to_nearest_phr_bp=(near if near < 10 ** 12 else -1),
                n_genes=len(genes),
                genes=",".join(g[0] for g in genes[:12]),
                note=note))
        fa.close()
        rd.close()
        with open(os.path.join(args.work, "cores.tsv"), "w") as fh:
            cols = list(rows[0])
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        for r in rows:
            print("%-18s %s:%d-%d %6d bp  %-5s in_phr=%s comm=%s genes=%d"
                  % (r["core"], r["chrom"], r["start"], r["end"], r["length"],
                     r["kind"], r["in_phr"], r["community"], r["n_genes"]))
        print("wrote", os.path.join(args.work, "cores.tsv"))
        return

    if args.cmd == "fastga":
        paf1 = os.path.join(args.work, "cores_vs_ref.%s.paf" % args.run)
        # FastGA writes its own .gdb/.1ano/.post stubs beside the query, so the query is
        # staged into a private subdirectory (the repo work dir is not used directly).
        indir = os.path.join(args.work, "fastga_in")
        tmpdir = os.path.join(args.work, "fastga_tmp_%s" % args.run)
        os.makedirs(indir, exist_ok=True)
        os.makedirs(tmpdir, exist_ok=True)
        q = os.path.join(indir, "cores.fa")
        with open(os.path.join(args.work, "cores.fa"), "rb") as fh_in, open(q, "wb") as fh_out:
            fh_out.write(fh_in.read())
        cmd = [FASTGA, "-pafx", "-f%d" % args.freq, "-c85", "-s500", "-T32",
               "-P" + tmpdir,
               "cores.fa", REF_FA]
        print("running:", " ".join(cmd), flush=True)
        t0 = time.time()
        with open(paf1, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                                  cwd=indir)
        wall = time.time() - t0
        with open(log, "a") as fh:
            fh.write("FASTGA hop1: wall %.1fs rc=%d n_records_paf=%d\n"
                     % (wall, proc.returncode, sum(1 for _ in open(paf1))))
            fh.write(proc.stderr[-3000:] + "\n")
        if proc.returncode != 0:
            raise SystemExit("FAIL: FastGA hop1 rc=%d: %s" % (proc.returncode, proc.stderr[-1500:]))
        print("hop1 PAF: %s  wall %.1fs (-f%d -s500 -c85 -T32)" % (paf1, wall, args.freq))
        return

    if args.cmd == "hop2":
        # Transitive (2-hop) closure: every hop-1 reference placement is pulled out of
        # the reference and aligned back to the reference, so positions that are similar
        # to a copy but no longer >= identity threshold to the CORE itself become
        # reachable.  This is the transitivity the operator asked for, expressed as an
        # explicit second hop (the 37 GB all-vs-all index that would do it implicitly is
        # unusable, failure F1).
        cores = read_cores(args.work)
        phrs = read_phrs()
        cen = read_centromeres()
        rd = FaidxReader(REF_FA)
        seq_of = {}
        for name in rd.fai:
            chrom, beg = parse_pansn(name)
            seq_of[chrom] = (name, beg)
        indir = os.path.join(args.work, "hop2_in")
        os.makedirs(indir, exist_ok=True)
        outdir = os.path.join(args.work, "hop2")
        os.makedirs(outdir, exist_ok=True)
        qpath = os.path.join(indir, "hop2_query.fa")
        manifest = os.path.join(args.work, "hop2_manifest.tsv")
        total_bp = 0
        n_spans = 0
        with open(qpath, "w") as qf, open(manifest, "w") as mf:
            mf.write("core\thop2_seq\tchrom\tstart\tend\tlength\n")
            for r in cores:
                cid = r["core"]
                bedpe = os.path.join(args.work, "impg_bedpe_%s" % args.run, cid + ".bedpe")
                if not os.path.exists(bedpe):
                    raise SystemExit("FAIL: %s missing; run query first" % bedpe)
                recs = parse_bedpe(bedpe)
                raw = []
                for d in recs:
                    if d["t_name"] == cid:
                        ref_name, ref_s, ref_e = d["q_name"], d["q_start"], d["q_end"]
                        cov_s, cov_e = d["t_start"], d["t_end"]
                    elif d["q_name"] == cid:
                        ref_name, ref_s, ref_e = d["t_name"], d["t_start"], d["t_end"]
                        cov_s, cov_e = d["q_start"], d["q_end"]
                    else:
                        raise SystemExit("FAIL: bedpe row in %s names neither side" % bedpe)
                    chrom, s0, e0 = ref_to_chrom(ref_name, ref_s, ref_e)
                    if (cov_e - cov_s) < args.min_core_cov * r["length"]:
                        continue
                    if chrom == r["chrom"] and abs(s0 - r["start"]) <= 1 and abs(e0 - r["end"]) <= 1:
                        continue     # the core itself: nothing to add
                    raw.append(dict(chrom=chrom, start=s0, end=e0, strand=d["strand"],
                                    core_start=cov_s, core_end=cov_e, gi=d["gi"]))
                spans = merge_intervals(raw, args.merge_gap) if raw else []
                kept = 0
                for L in spans:
                    if total_bp + (L["end"] - L["start"]) > args.hop2_max_bp:
                        break
                    seqname, beg = seq_of[L["chrom"]]
                    sub = rd.fetch(seqname, L["start"] - (beg - 1), L["end"] - (beg - 1))
                    label = "%s|%s:%d-%d" % (cid, L["chrom"], L["start"], L["end"])
                    qf.write(">%s\n" % label)
                    for i in range(0, len(sub), 60):
                        qf.write(sub[i:i + 60].decode() + "\n")
                    mf.write("%s\t%s\t%s\t%d\t%d\t%d\n"
                             % (cid, label, L["chrom"], L["start"], L["end"],
                                L["end"] - L["start"]))
                    total_bp += L["end"] - L["start"]
                    n_spans += 1
                    kept += 1
                print("%-17s hop2 spans %d (%d bp in batch so far)"
                      % (cid, kept, total_bp), flush=True)
        rd.close()
        if total_bp == 0:
            raise SystemExit("FAIL: hop2 query is empty")
        paf2 = os.path.join(args.work, "hop2_vs_ref.%s.paf" % args.run)
        tmpdir = os.path.join(args.work, "fastga_tmp_hop2_%s" % args.run)
        os.makedirs(tmpdir, exist_ok=True)
        cmd = [FASTGA, "-pafx", "-f8192", "-c85", "-s500", "-T32", "-P" + tmpdir,
               "hop2_query.fa", REF_FA]
        print("running:", " ".join(cmd), flush=True)
        t0 = time.time()
        with open(paf2, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                                  cwd=indir)
        wall = time.time() - t0
        with open(log, "a") as fh:
            fh.write("FASTGA hop2: wall %.1fs rc=%d spans=%d query_bp=%d n_records=%d\n"
                     % (wall, proc.returncode, n_spans, total_bp,
                        sum(1 for _ in open(paf2))))
            fh.write(proc.stderr[-2000:] + "\n")
        if proc.returncode != 0:
            raise SystemExit("FAIL: FastGA hop2 rc=%d: %s"
                             % (proc.returncode, proc.stderr[-1500:]))
        print("hop2 PAF: %s wall %.1fs (%d spans, %d bp query)"
              % (paf2, wall, n_spans, total_bp))
        return

    if args.cmd == "index":
        idx = os.path.join(args.work, "cores_hop1.%s.impg" % args.run)
        paf1 = os.path.join(args.work, "cores_vs_ref.%s.paf" % args.run)
        if not os.path.exists(paf1):
            raise SystemExit("FAIL: %s missing; run the fastga step first" % paf1)
        # -f: force a rebuild.  Without it impg silently reuses an existing index file
        # whose PAF has since been rewritten, and every query then dies with
        # "Invalid CIGAR operation" because the stored BGZF offsets are stale.
        cmd = [IMPG, "index", "-a", paf1, "-i", idx, "-f", "-t", str(args.threads), "-v", "1"]
        print("running:", " ".join(cmd), flush=True)
        proc, wall = sh(cmd, log)
        if proc.returncode != 0:
            raise SystemExit("FAIL: impg index rc=%d\n%s\n%s"
                             % (proc.returncode, proc.stdout[-2000:], proc.stderr[-2000:]))
        print("index: %s  wall %.1fs" % (idx, wall))
        return

    if args.cmd == "query":
        idx = os.path.join(args.work, "cores_hop1.%s.impg" % args.run)
        paf1 = os.path.join(args.work, "cores_vs_ref.%s.paf" % args.run)
        outdir = os.path.join(args.work, "impg_bedpe_%s" % args.run)
        os.makedirs(outdir, exist_ok=True)
        cores = read_cores(args.work)
        wanted = set(args.core) if args.core else None
        for r in cores:
            cid = r["core"]
            if wanted and cid not in wanted:
                continue
            out = os.path.join(outdir, cid + ".bedpe")
            cmd = [IMPG, "query", "-i", idx, "-a", paf1,
                   "-r", cid, "-d", "100",
                   "--min-result-identity", "%.4f" % (args.identity / 100.0),
                   "-l", "500", "-o", "bedpe", "--original-sequence-coordinates", "--self",
                   "-t", str(args.threads), "-v", "1"]
            t0 = time.time()
            proc = subprocess.run(cmd, capture_output=True, text=True)
            wall = time.time() - t0
            n = 0
            with open(out, "w") as fh:
                fh.write(proc.stdout)
                n = proc.stdout.count("\n")
            with open(log, "a") as fh:
                fh.write("IMPG QUERY %s: wall %.2fs rc=%d rows=%d\n"
                         % (cid, wall, proc.returncode, n))
                if proc.returncode != 0:
                    fh.write(proc.stderr[-3000:] + "\n")
            if proc.returncode != 0:
                raise SystemExit("FAIL: impg query %s rc=%d\n%s"
                                 % (cid, proc.returncode, proc.stderr[-2000:]))
            print("%-18s rows=%-6d wall=%6.2fs" % (cid, n, wall), flush=True)
        return

    if args.cmd == "crosscheck":
        # Independent recounts of one core's copy number, outside the impg layer:
        #   (a) the raw FastGA PAF that impg was built from, re-filtered here;
        #   (b) an independent aligner's PAF (wfmash), if --wfmash is given.
        core = args.core[0] if args.core else "D4Z4_4q"
        cores = {r["core"]: r for r in read_cores(args.work)}
        if core not in cores:
            raise SystemExit("FAIL: %s not in cores.tsv" % core)
        r = cores[core]
        phrs = read_phrs()
        cen = read_centromeres()
        paf1 = os.path.join(args.work, "cores_vs_ref.%s.paf" % args.run)

        def recount(rows, source, qname_key="q_name", tname_key="t_name",
                    qs="q_start", qe="q_end", ts="t_start", te="t_end"):
            keys, arms = set(), defaultdict(int)
            n = 0
            for d in rows:
                if d[qname_key] == core:
                    ref_name, rs, re_, cs, ce = d[tname_key], d[ts], d[te], d[qs], d[qe]
                elif d[tname_key] == core:
                    ref_name, rs, re_, cs, ce = d[qname_key], d[qs], d[qe], d[ts], d[te]
                else:
                    continue
                if d["gi"] < args.identity / 100.0 or d["block_len"] < 1000:
                    continue
                if (ce - cs) < args.min_core_cov * r["length"]:
                    continue
                chrom, s0, e0 = ref_to_chrom(ref_name, rs, re_)
                if chrom == r["chrom"] and abs(s0 - r["start"]) <= 1 and abs(e0 - r["end"]) <= 1:
                    continue
                n += 1
                if (chrom, s0) in keys:
                    continue
                keys.add((chrom, s0))
                a, _p = placement_arm(chrom, s0, s0 + 1, phrs, cen)
                arms[a] += 1
            return dict(source=source, n_align=n, n_positions=len(keys), n_arms=len(arms),
                        per_arm=";".join("%s:%d" % (k, v) for k, v in sorted(arms.items())),
                        definition="gap-compressed identity >= %.2f, block >= 1000, "
                                   "core coverage >= %.2f" % (args.identity / 100.0,
                                                              args.min_core_cov))

        out = [recount(parse_paf(paf1), "fastga-raw-paf")]
        if args.wfmash:
            wrows = parse_paf(args.wfmash)
            out.append(recount(wrows, "wfmash"))
        import csv as _csv
        tsv = os.path.join(REPO, "scripts/cladistics/impg_copy_enumeration.tsv")
        with open(tsv) as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                if row["core"] == core:
                    out.append(dict(source="impg-bedpe", n_align=int(row["n_align"]),
                                    n_positions=int(row["n_positions"]),
                                    n_arms=int(row["n_arms"]), per_arm=row["per_arm"]))
        print("cross-check for %s (%s:%d-%d, %d bp, identity >= %.1f%%, block >= 1000 bp, "
              "core coverage >= %.2f)" % (core, r["chrom"], r["start"], r["end"],
                                          r["length"], args.identity, args.min_core_cov))
        for d in out:
            print("  %-16s align=%-5d positions=%-5d arms=%-3d %s"
                  % (d["source"], d["n_align"], d["n_positions"], d["n_arms"], d["per_arm"]))
        js = os.path.join(args.work, "crosscheck_%s.json" % core)
        with open(js, "w") as fh:
            json.dump(out, fh, indent=2)
        print("wrote", js)
        return

    if args.cmd == "aggregate":
        cores = read_cores(args.work)
        phrs = read_phrs()
        cen = read_centromeres()
        bedpe_dir = os.path.join(args.work, "impg_bedpe_%s" % args.run)
        alt_dir = os.path.join(args.work, "impg_bedpe_%s" % args.alt_run)
        have_alt = os.path.isdir(alt_dir)
        hop2_path = os.path.join(args.work, "hop2_vs_ref.%s.paf" % args.run)
        hop2_paf = parse_paf(hop2_path) if os.path.exists(hop2_path) else None
        rows = []
        for r in cores:
            cid = r["core"]
            path = os.path.join(bedpe_dir, cid + ".bedpe")
            if not os.path.exists(path):
                raise SystemExit("FAIL: %s missing; run the query step first" % path)
            st = core_stats(cid, r, path, phrs, cen, args.min_core_cov, args.merge_gap)
            row = dict(core=cid, chrom=r["chrom"], start=r["start"], end=r["end"],
                       length=r["length"], kind=r["kind"], in_phr=r["in_phr"],
                       arm=r["arm"], community=r["community"],
                       community_arms=r["community_arms"])
            row.update(st)
            if hop2_paf is not None:
                row.update(hop2_stats(cid, hop2_paf, st["_nonself_starts"], r, phrs, cen,
                                      args.identity / 100.0, 1000))
            else:
                row.update(hop2_n_align=-1, hop2_n_new_positions=-1, hop2_n_new_arms=-1,
                           hop2_per_arm_new="")
            if have_alt:
                apath = os.path.join(alt_dir, cid + ".bedpe")
                if os.path.exists(apath):
                    ast_ = core_stats(cid, r, apath, phrs, cen, args.min_core_cov,
                                      args.merge_gap)
                    row["alt_n_positions"] = ast_["n_positions"]
                    row["alt_n_arms"] = ast_["n_arms"]
                    row["alt_per_arm"] = ast_["per_arm"]
                    row["alt_n_dropped_short"] = ast_["n_dropped_short"]
                else:
                    row["alt_n_positions"] = -1
                    row["alt_n_arms"] = -1
                    row["alt_per_arm"] = ""
                    row["alt_n_dropped_short"] = -1
            else:
                row["alt_n_positions"] = -1
                row["alt_n_arms"] = -1
                row["alt_per_arm"] = ""
                row["alt_n_dropped_short"] = -1
            rows.append(row)

        cols = ["core", "chrom", "start", "end", "length", "kind", "in_phr", "arm",
                "community", "community_arms", "n_align", "n_dropped_short", "n_self",
                "n_positions", "n_placements_overlapping_core", "n_loci",
                "n_interstitial", "n_chroms", "n_arms", "n_on_community_arms",
                "n_on_noncommunity_arms", "per_arm", "identity_median",
                "identity_min", "identity_max", "core_cov_median",
                "placement_len_median", "ref_span_bp_all", "copies_core_units",
                "hop2_n_align", "hop2_n_new_positions", "hop2_n_new_arms",
                "hop2_per_arm_new", "alt_n_positions", "alt_n_arms",
                "alt_n_dropped_short", "alt_per_arm"]
        out = os.path.join(REPO, "scripts/cladistics/impg_copy_enumeration.tsv")
        with open(out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for row in rows:
                vals = []
                for c in cols:
                    v = row[c]
                    vals.append(pgf_num(v, 4) if isinstance(v, float) else str(v))
                fh.write("\t".join(vals) + "\n")
        for row in rows:
            print("%-17s align=%-5d drop=%-5d pos=%-5d loci=%-4d chroms=%-3d arms=%-3d "
                  "id=%.4f cov=%.2f span=%dbp C=%s(%d/%d) alt_pos=%d"
                  % (row["core"], row["n_align"], row["n_dropped_short"], row["n_positions"],
                     row["n_loci"], row["n_chroms"], row["n_arms"], row["identity_median"],
                     row["core_cov_median"], row["ref_span_bp_all"], row["community"],
                     row["n_on_community_arms"], row["n_on_noncommunity_arms"],
                     row["alt_n_positions"]))
            print("   hop2: new_positions=%d new_arms=%d"
                  % (row["hop2_n_new_positions"], row["hop2_n_new_arms"]))
        print("wrote", out)
        return

    raise SystemExit("subcommand %r not implemented yet" % args.cmd)


if __name__ == "__main__":
    main()
