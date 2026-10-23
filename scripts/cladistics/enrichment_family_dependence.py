#!/usr/bin/env python3
"""Family dependence of the CHM13 PHR copy-ontology enrichment (Referee #1, major point 9).

Referee #1, major point 9 (quoted verbatim from PEER_REVIEW.md):

  "Of 402 gene copies with midpoints in CHM13 PHRs, 215 are unannotated and excluded,
   leaving K = 187 across 31,235 term combinations, so over half the PHR gene content is
   invisible. The rest are near-identical members of a few amplified families, so one
   family drives a term. Thus, the 25-fold 'spermatid development' enrichment rests on 19
   copies of a single SEPTIN14-associated duplicon. The whole genome is also the wrong
   comparator ... Please repeat against identity- and length-matched non-subtelomeric
   duplications and with copies collapsed to families, reporting how many of the 209 terms
   survive."

Scope of THIS script (the arithmetic + the family-dependence half of point 9; the
matched-duplication control is a separate deliverable):

  (1) Verify or correct the referee's arithmetic, from raw repository artefacts:
        402 PHR-midpoint gene copies          -> independent recount from the GFF3 + BED
        215 excluded / 187 retained (K)       -> pinned ledger MAPPING_COVERAGE.tsv
        31,235 term combinations              -> row count of TERM_RESULTS.tsv.gz
        209 supported terms                   -> recomputed from the frozen p-values
      Also test the referee's phrasing "215 are unannotated": all 402 PHR copies have
      their own annotation identity; the 215 are copies with no *defensible ontology
      source*, which is a different (and weaker) statement than "unannotated".

  (2) Leave-largest-family-out for every significantly enriched term: how many of a term's
      PHR copies come from its single largest contributing family, and does the enrichment
      survive when that family's copies are removed? The referee's specific claim about
      "spermatid development" (GO:0007286, 19 copies, ~25-fold, SEPTIN14) is tested
      exactly.

  (3) The invisible-content fraction and what it does to the interpretation.

Two assignment regimes and two confidence stratifications (the enrichment analogue of the
"two aligner regimes" rule -- there is no alignment in this analysis):
  regime "midpoint"    : primary frozen partition (K = 187 PHR copies)
  regime "overlap"     : paired any-overlap sensitivity partition (K = 193 PHR copies)
  family definition "source" : functional_source_id (the ontology-inheritance family root)
  family definition "locus"  : gene_name (RefSeq Liftoff locus symbol)
  confidence columns: how many of a term's PHR copies rest on inherited
  (NCBI Gene "Related functional gene") annotation rather than exact self annotation.

Nothing here is estimated. Every number printed is computed from a file whose path,
byte size and SHA-256 are echoed. Any failed expectation aborts the run.

Usage
-----
  enrichment_family_dependence.py                 # full run, writes the TSV
  enrichment_family_dependence.py --output -      # full run, TSV to stdout
  enrichment_family_dependence.py --self-test     # exact-tail oracle + adversarial cases
"""

import argparse
import gzip
import hashlib
import itertools
import math
import os
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ONT_TABLE = os.path.join(REPO, "submission/supporting_material/ED_Fig2_supported_ontology_terms.tsv")
ONTV7 = os.path.join(REPO, "paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7")
TERM_RESULTS = os.path.join(ONTV7, "TERM_RESULTS.tsv.gz")
CONTRIBUTORS = os.path.join(ONTV7, "EXACT_TERM_CONTRIBUTORS.tsv.gz")
COVERAGE = os.path.join(ONTV7, "MAPPING_COVERAGE.tsv")
COHORT_AUDIT = os.path.join(ONTV7, "NAMED_COHORT_AUDIT.tsv")
GFF3 = os.path.join(REPO, "data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
PHR_BED = os.path.join(REPO, "data/chm13.phrs.bed")

DEFAULT_OUT = os.path.join(REPO, "scripts/cladistics/enrichment_family_dependence.tsv")

# Referee #1's stated quantities; each is checked against a computed value below.
REFEREE = {
    "phr_midpoint_copies": 402,
    "phr_excluded": 215,
    "K": 187,
    "term_combinations": 31235,
    "supported_terms": 209,
    "spermatid_copies": 19,
    "spermatid_fold": 25.0,
}

SELF_TEST_LOG = []


def fail(msg):
    raise SystemExit("FATAL: %s" % msg)


def note(msg):
    SELF_TEST_LOG.append(msg)
    print(msg, file=sys.stderr)


# --------------------------------------------------------------------------------------
# exact hypergeometric upper tail
# --------------------------------------------------------------------------------------
def hypergeom_upper(a, K, M, N):
    """Exact P(X >= a) for X ~ Hypergeometric(N, M, K), as an exact Fraction.

    a observed term-bearing copies among K PHR copies in a fixed population of N eligible
    copies containing M term-bearing copies (the frozen V7 finite-population test).
    """
    if K < 0 or M < 0 or N < 0 or K > N or M > N:
        fail("invalid hypergeometric parameters N=%r K=%r M=%r" % (N, K, M))
    # X is the number of term-bearing copies drawn into the K-copy PHR set: 0 <= X <= min(K, M).
    if a <= 0:
        return Fraction(1)
    if K == 0:
        return Fraction(0)  # no PHR copy can be term-bearing
    hi = min(K, M)
    if a > hi:
        return Fraction(0)
    total = math.comb(N, K)
    num = 0
    for x in range(a, hi + 1):
        if K - x > N - M:
            continue
        num += math.comb(M, x) * math.comb(N - M, K - x)
    return Fraction(num, total)


def brute_force_upper(a, K, M, N):
    """Oracle: enumerate every K-subset of N items of which M are marked."""
    hits = 0
    total = 0
    for subset in itertools.combinations(range(N), K):
        total += 1
        if sum(1 for i in subset if i < M) >= a:
            hits += 1
    if total != math.comb(N, K):
        fail("oracle enumeration count %d != C(%d,%d)" % (total, N, K))
    return Fraction(hits, total)


def run_self_test():
    print("== exact upper-tail oracle (tiny scale) ==", file=sys.stderr)
    cases = [
        (0, 0, 0, 1),   # K = 0 -> p = 1
        (1, 1, 1, 1),   # degenerate certainty
        (0, 3, 4, 10),
        (1, 3, 4, 10),
        (2, 3, 4, 10),
        (3, 3, 4, 10),
        (3, 3, 4, 10),
        (2, 5, 5, 12),
        (1, 1, 10, 10),  # M = N -> any PHR copy is term-bearing
        (4, 4, 10, 30),
        (5, 2, 3, 9),   # a > min(K, M) -> 0
        (0, 4, 0, 10),  # no term-bearing copy anywhere -> 1
    ]
    for a, K, M, N in cases:
        exact = hypergeom_upper(a, K, M, N)
        oracle = brute_force_upper(a, K, M, N)
        if exact != oracle:
            fail("exact tail %s != brute force %s for a=%d K=%d M=%d N=%d"
                 % (exact, oracle, a, K, M, N))
        note("  ok  N=%-3d K=%-2d M=%-2d a=%-2d  exact=%-22s oracle=%s"
             % (N, K, M, a, exact, oracle))

    # adversarial: exhaustive over all tiny parameter sets
    bad = 0
    for N in range(1, 13):
        for K in range(0, min(N, 4) + 1):
            for M in range(0, N + 1):
                for a in range(0, K + 2):
                    exact = hypergeom_upper(a, K, M, N)
                    oracle = brute_force_upper(a, K, M, N)
                    if exact != oracle:
                        bad += 1
                        note("  MISMATCH N=%d K=%d M=%d a=%d: %s vs %s"
                             % (N, K, M, a, exact, oracle))
    if bad:
        fail("%d exhaustive tiny-scale oracle mismatches" % bad)
    note("  ok  exhaustive tiny-scale sweep: every (N<=12, K<=4, M, a) matches the oracle")

    # monotonicity / limits on a realistic scale
    N, K, M = 31966, 187, 146
    prev = None
    for a in range(0, K + 1):
        p = float(hypergeom_upper(a, K, M, N))
        if prev is not None and p > prev + 0.0:
            fail("upper tail not monotone non-increasing at a=%d" % a)
        prev = p
    note("  ok  upper tail monotone non-increasing over a=0..K at N=%d K=%d M=%d" % (N, K, M))
    note("  187/31966 population: p(a=19, M=146) = %.6e" % float(hypergeom_upper(19, 187, 146, 31966)))
    note("  187/31966 population: p(a=1, M=1)     = %.6e" % float(hypergeom_upper(1, 187, 1, 31966)))
    print("== self-test PASS ==", file=sys.stderr)


# --------------------------------------------------------------------------------------
# multiplicity: BH within collection, Holm global (definitions are recovered from the
# frozen table rather than assumed; see verify_multiplicity_definition)
# --------------------------------------------------------------------------------------
def bh_qvalues(pvals):
    """Benjamini-Hochberg step-up q values."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = pvals[i] * n / rank
        running = min(running, val)
        q[i] = min(1.0, running)
    return q


def holm_pvalues(pvals):
    """Holm step-down adjusted p values (monotone enforced)."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    running = 0.0
    for rank in range(1, n + 1):
        i = order[rank - 1]
        val = (n - rank + 1) * pvals[i]
        running = max(running, val)
        adj[i] = min(1.0, running)
    return adj


# --------------------------------------------------------------------------------------
# absolute-path / provenance helpers
# --------------------------------------------------------------------------------------
def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(paths):
    out = []
    for p in paths:
        if not os.path.exists(p):
            fail("required input is missing: %s" % p)
        out.append((os.path.relpath(p, REPO), os.path.getsize(p), sha256_of(p)))
    return out


def read_tsv_gz(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            yield dict(zip(header, line.split("\t")))


# --------------------------------------------------------------------------------------
# step 1: independent recount of PHR gene copies from the raw GFF3 + BED
# --------------------------------------------------------------------------------------
def load_phr_bed(path):
    intervals = defaultdict(list)
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            intervals[f[0]].append((int(f[1]), int(f[2])))
            n += 1
    for c in intervals:
        intervals[c].sort()
    if n == 0:
        fail("empty PHR BED: %s" % path)
    return intervals, n


def load_gff3_genes(path):
    genes = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            genes.append({
                "seqid": f[0],
                "start0": int(f[3]) - 1,
                "end0": int(f[4]),
                "biotype": attrs.get("gene_biotype", ""),
                "gene_name": attrs.get("gene_name", attrs.get("gene", "")),
                "gene_id": attrs.get("ID", ""),
            })
    if not genes:
        fail("no gene records parsed from %s" % path)
    return genes


def midpoint_hits(genes, intervals):
    """Return (n_midpoint, n_any_overlap) for gene records against the PHR intervals."""
    n_mid = 0
    n_ovl = 0
    mid_hit = []
    for g in genes:
        mid = (g["start0"] + g["end0"]) // 2
        hit_mid = False
        hit_ovl = False
        for s, e in intervals.get(g["seqid"], ()):
            if s <= mid < e:
                hit_mid = True
            if g["start0"] < e and g["end0"] > s:
                hit_ovl = True
            if hit_mid and hit_ovl:
                break
        if hit_mid:
            n_mid += 1
            mid_hit.append(g)
        if hit_ovl:
            n_ovl += 1
    return n_mid, n_ovl, mid_hit


# --------------------------------------------------------------------------------------
# step 6 machinery: leave-largest-family-out
# --------------------------------------------------------------------------------------
def tops(counter, n):
    """Deterministic 'most common': count desc, key asc (avoids hash-order ties)."""
    return sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[:n]


def largest_family(counts):
    """Deterministic: max copy count, ties broken by family id; returns (fam, n, n_tied)."""
    if not counts:
        return None, 0, 0
    best = max(counts.values())
    tied = sorted([f for f, v in counts.items() if v == best])
    return tied[0], best, len(tied)


def main():
    t0 = time.time()
    print("== inputs (relative path, bytes, sha256) ==", file=sys.stderr)
    prov = provenance([ONT_TABLE, TERM_RESULTS, CONTRIBUTORS, COVERAGE, COHORT_AUDIT,
                       GFF3, PHR_BED])
    for rel, size, digest in prov:
        note("  %-95s %10d  %s" % (rel, size, digest))

    # ---------------- (1) arithmetic ----------------
    note("")
    note("== (1) referee arithmetic ==")
    intervals, n_phr_intervals = load_phr_bed(PHR_BED)
    note("  PHR BED intervals: %d" % n_phr_intervals)
    genes = load_gff3_genes(GFF3)
    note("  CHM13v2.0 GFF3 gene records: %d" % len(genes))
    n_mid, n_ovl, mid_genes = midpoint_hits(genes, intervals)
    note("  independent recount, gene midpoints inside PHR BED: %d (referee: %d)"
         % (n_mid, REFEREE["phr_midpoint_copies"]))
    note("  independent recount, any overlap with PHR BED:     %d" % n_ovl)
    if n_mid != REFEREE["phr_midpoint_copies"]:
        fail("PHR-midpoint recount %d != referee %d" % (n_mid, REFEREE["phr_midpoint_copies"]))

    cov = list(read_tsv_gz(COVERAGE))
    cov_mid = [r for r in cov if r["assignment"] == "midpoint"]
    cov_ovl = [r for r in cov if r["assignment"] == "overlap"]
    if len(cov_mid) != 1 or len(cov_ovl) != 1:
        fail("MAPPING_COVERAGE.tsv should have one midpoint and one overlap row")
    cm, co = cov_mid[0], cov_ovl[0]
    led = {
        "phr_raw_mid": int(cm["phr_raw_physical_copy_cn"]),
        "phr_elig_mid": int(cm["phr_ontology_eligible_cn"]),
        "phr_inelig_mid": int(cm["phr_ontology_ineligible_cn"]),
        "nonphr_elig_mid": int(cm["nonphr_ontology_eligible_cn"]),
        "univ_mid": int(cm["ontology_eligible_cn"]),
        "raw_total": int(cm["raw_physical_copy_cn"]),
        "inelig_total": int(cm["ontology_ineligible_cn"]),
        "phr_raw_ovl": int(co["phr_raw_physical_copy_cn"]),
        "phr_elig_ovl": int(co["phr_ontology_eligible_cn"]),
        "phr_inelig_ovl": int(co["phr_ontology_ineligible_cn"]),
        "nonphr_elig_ovl": int(co["nonphr_ontology_eligible_cn"]),
    }
    note("  MAPPING_COVERAGE.tsv midpoint ledger: raw PHR %d, eligible PHR %d, ineligible PHR %d"
         % (led["phr_raw_mid"], led["phr_elig_mid"], led["phr_inelig_mid"]))
    note("  MAPPING_COVERAGE.tsv overlap  ledger: raw PHR %d, eligible PHR %d, ineligible PHR %d"
         % (led["phr_raw_ovl"], led["phr_elig_ovl"], led["phr_inelig_ovl"]))
    if led["phr_raw_mid"] != n_mid:
        fail("ledger raw PHR %d != my recount %d" % (led["phr_raw_mid"], n_mid))
    if led["phr_raw_ovl"] != n_ovl:
        fail("ledger overlap PHR %d != my recount %d" % (led["phr_raw_ovl"], n_ovl))
    if led["phr_elig_mid"] != REFEREE["K"]:
        fail("ledger K %d != referee %d" % (led["phr_elig_mid"], REFEREE["K"]))
    if led["phr_raw_mid"] - led["phr_elig_mid"] != REFEREE["phr_excluded"]:
        fail("ledger excluded count %d != referee %d"
             % (led["phr_raw_mid"] - led["phr_elig_mid"], REFEREE["phr_excluded"]))
    note("  VERIFIED: %d raw PHR copies - %d excluded = %d retained (K); excluded = %d (referee: %d)"
         % (led["phr_raw_mid"], led["phr_inelig_mid"], led["phr_elig_mid"],
            led["phr_inelig_mid"], REFEREE["phr_excluded"]))

    # term results
    rows = {"midpoint": [], "overlap": []}
    for r in read_tsv_gz(TERM_RESULTS):
        a = r["assignment"]
        if a not in rows:
            continue
        rows[a].append(r)
    for regime in ("midpoint", "overlap"):
        note("  TERM_RESULTS.tsv.gz %s rows: %d (referee hypothesis count: %d)"
             % (regime, len(rows[regime]), REFEREE["term_combinations"]))
    if len(rows["midpoint"]) != REFEREE["term_combinations"]:
        fail("midpoint hypothesis count %d != referee %d"
             % (len(rows["midpoint"]), REFEREE["term_combinations"]))
    if len(rows["overlap"]) != REFEREE["term_combinations"]:
        fail("overlap hypothesis count %d != referee %d"
             % (len(rows["overlap"]), REFEREE["term_combinations"]))

    for regime, exp_K, exp_non, exp_univ in (("midpoint", led["phr_elig_mid"], led["nonphr_elig_mid"], led["univ_mid"]),
                                             ("overlap", led["phr_elig_ovl"], led["nonphr_elig_ovl"], led["univ_mid"])):
        Ks = {int(r["K_phr_total"]) for r in rows[regime]}
        ns = {int(r["nonphr_total"]) for r in rows[regime]}
        us = {int(r["N_eligible_total"]) for r in rows[regime]}
        if Ks != {exp_K} or ns != {exp_non} or us != {exp_univ}:
            fail("%s regime universe inconsistent: K=%s non-PHR=%s N=%s (expected %d/%d/%d)"
                 % (regime, Ks, ns, us, exp_K, exp_non, exp_univ))
        note("  %s universe constant: K=%d non-PHR=%d N=%d" % (regime, exp_K, exp_non, exp_univ))

    # ---------------- multiplicity definition recovery + 209 gate ----------------
    note("")
    note("== (1b) reproduce the frozen 209 supported terms from the frozen p-values ==")
    for regime in ("midpoint", "overlap"):
        rs = rows[regime]
        pv = [float(r["p_exact_upper"]) for r in rs]
        by_col = defaultdict(list)
        for i, r in enumerate(rs):
            by_col[r["collection"]].append(i)
        q = [0.0] * len(rs)
        for col, idx in by_col.items():
            sub = bh_qvalues([pv[i] for i in idx])
            for j, i in enumerate(idx):
                q[i] = sub[j]
        holm = holm_pvalues(pv)
        max_q_diff = max(abs(q[i] - float(rs[i]["bh_q_within_collection"])) for i in range(len(rs)))
        max_h_diff = max(abs(holm[i] - float(rs[i]["holm_p_global"])) for i in range(len(rs)))
        note("  %s: max |recomputed BH q - frozen bh_q| = %.3e ; max |recomputed Holm - frozen holm| = %.3e"
             % (regime, max_q_diff, max_h_diff))
        if max_q_diff > 1e-12 or max_h_diff > 1e-12:
            fail("%s multiplicity definitions do not reproduce the frozen columns "
                 "(q diff %.3e, holm diff %.3e)" % (regime, max_q_diff, max_h_diff))
        n_sup = sum(1 for i, r in enumerate(rs)
                    if q[i] <= 0.05 and holm[i] <= 0.05)
        support_col = "primary_support" if regime == "midpoint" else "sensitivity_support"
        n_sup_frozen = sum(1 for r in rs if r[support_col] == "1")
        note("  %s: recomputed support (BH<=0.05 and Holm<=0.05) = %d ; frozen %s=1: %d"
             % (regime, n_sup, support_col, n_sup_frozen))
        if n_sup != n_sup_frozen:
            fail("%s recomputed support count %d != frozen %d" % (regime, n_sup, n_sup_frozen))
        if regime == "midpoint" and n_sup != REFEREE["supported_terms"]:
            fail("midpoint support count %d != referee 209" % n_sup)
        if regime == "midpoint":
            note("  VERIFIED: 31,235 hypotheses, %d supported terms (referee: 209)" % n_sup)

    # published supporting-material table must be the midpoint supported subset
    ont = list(read_tsv_gz(ONT_TABLE))
    note("  ED_Fig2_supported_ontology_terms.tsv rows: %d" % len(ont))
    if len(ont) != REFEREE["supported_terms"]:
        fail("supporting table has %d rows, expected %d" % (len(ont), REFEREE["supported_terms"]))
    key = lambda r: (r["collection"], r["relation"], r["term_id"])
    frozen_mid_by_key = {key(r): r for r in rows["midpoint"]}
    for r in ont:
        fr = frozen_mid_by_key.get(key(r))
        if fr is None:
            fail("supported term %s absent from TERM_RESULTS midpoint rows" % (key(r),))
        for col in ("a_T", "b_T", "c_T", "d_T", "K_phr_total", "nonphr_total", "N_eligible_total"):
            if int(r[col]) != int(fr[col]):
                fail("2x2 mismatch for %s column %s: table %s vs frozen %s"
                     % (key(r), col, r[col], fr[col]))
        if float(r["p_exact_upper"]) != float(fr["p_exact_upper"]):
            fail("p mismatch for %s" % (key(r),))
    note("  VERIFIED: all 209 supporting-table rows agree with the frozen midpoint rows "
         "(a,b,c,d,K,non-PHR,N,p)")

    # ---------------- per-copy contributor ledger ----------------
    note("")
    note("== contributor ledger ==")
    # cell -> (collection, relation, term_id) -> Counter(copy_id)
    terms_mid = defaultdict(Counter)   # term key -> copy_id -> n (n always 1 per row)
    terms_ovl = defaultdict(Counter)
    copy_meta = {}                     # copy_id -> dict of family/confidence fields
    copy_part_mid = set()
    copy_part_ovl = set()
    n_ledger = 0
    dup_guard = set()
    for r in read_tsv_gz(CONTRIBUTORS):
        n_ledger += 1
        k = (r["collection"], r["relation"], r["term_id"])
        cid = r["copy_id"]
        if r["v7_midpoint_count_cell"] == "a_T":
            if (k, cid) in dup_guard:
                fail("duplicate ledger row for %s %s" % (k, cid))
            dup_guard.add((k, cid))
            terms_mid[k][cid] += 1
            copy_part_mid.add(cid)
        if r["v7_overlap_count_cell"] == "a_T":
            terms_ovl[k][cid] += 1
            copy_part_ovl.add(cid)
        meta = copy_meta.get(cid)
        new = {
            "functional_source_id": r["functional_source_id"],
            "functional_source_symbol": r["functional_source_symbol"],
            "gene_name": r["gene_name"],
            "source_assignment_disposition": r["source_assignment_disposition"],
            "mapping_confidence": r["mapping_confidence"],
            "seqid": r["seqid"],
            "start0": r["start0"],
            "end0": r["end0"],
            "gene_biotype": r["gene_biotype"],
        }
        if meta is None:
            copy_meta[cid] = new
        elif meta != new:
            fail("copy %s has inconsistent metadata across ledger rows:\n  %r\n  %r"
                 % (cid, meta, new))
    note("  ledger rows: %d" % n_ledger)
    note("  distinct PHR copies, midpoint partition: %d (ledger K=%d)"
         % (len(copy_part_mid), led["phr_elig_mid"]))
    note("  distinct PHR copies, overlap  partition: %d (ledger K=%d)"
         % (len(copy_part_ovl), led["phr_elig_ovl"]))
    if len(copy_part_mid) != led["phr_elig_mid"]:
        fail("midpoint partition copy count %d != %d" % (len(copy_part_mid), led["phr_elig_mid"]))
    if len(copy_part_ovl) != led["phr_elig_ovl"]:
        fail("overlap partition copy count %d != %d" % (len(copy_part_ovl), led["phr_elig_ovl"]))

    # ledger completeness: every frozen a_T must be reproduced by contributor rows
    for regime, table in (("midpoint", terms_mid), ("overlap", terms_ovl)):
        bad = 0
        for r in rows[regime]:
            k = key(r)
            if int(r["a_T"]) != len(table.get(k, {})):
                bad += 1
                if bad < 5:
                    note("    mismatch %s %s: frozen a_T=%s ledger=%d"
                         % (regime, k, r["a_T"], len(table.get(k, {}))))
        if bad:
            fail("%s: %d terms where ledger contributor count != frozen a_T" % (regime, bad))
        tot = sum(len(t) for t in table.values())
        note("  %s: ledger reproduces a_T for all %d hypotheses; total PHR term edges %d"
             % (regime, len(rows[regime]), tot))

    # ---------------- p-value implementation gate at fixture scale ----------------
    note("")
    note("== exact-tail implementation gate against the frozen p_exact_upper ==")
    for regime in ("midpoint", "overlap"):
        worst = 0.0
        worst_key = None
        n_checked = 0
        for r in rows[regime]:
            a = int(r["a_T"])
            if a <= 0:
                continue
            K = int(r["K_phr_total"])
            M = int(r["M_T_genome_burden"])
            N = int(r["N_eligible_total"])
            mine = float(hypergeom_upper(a, K, M, N))
            theirs = float(r["p_exact_upper"])
            n_checked += 1
            if theirs == 0.0:
                d = abs(mine)
            else:
                d = abs(mine - theirs) / theirs
            if d > worst:
                worst = d
                worst_key = (key(r), a, K, M, N, mine, theirs)
        note("  %s: %d rows with a_T>0, max relative deviation %.3e" % (regime, n_checked, worst))
        if worst > 1e-9:
            fail("%s exact-tail implementation deviates from the frozen p (%.3e) at %r"
                 % (regime, worst, worst_key))

    # ---------------- (2) leave-largest-family-out ----------------
    note("")
    note("== (2) family dependence of the 209 supported terms ==")
    for regime in ("midpoint", "overlap"):
        part = copy_part_mid if regime == "midpoint" else copy_part_ovl
        note("  %s PHR universe: %d copies" % (regime, len(part)))

    # PHR universe composition by family, for reporting
    for defn, field in (("source", "functional_source_id"), ("locus", "gene_name")):
        cnt = Counter(copy_meta[c][field] for c in copy_part_mid)
        note("  midpoint PHR universe by %s family (%d distinct):" % (defn, len(cnt)))
        for fam, n in tops(cnt, 12):
            note("    %-20s %3d  (%.1f%% of 187)" % (fam, n, 100.0 * n / len(copy_part_mid)))
        conf = Counter((copy_meta[c]["source_assignment_disposition"],
                        copy_meta[c]["mapping_confidence"]) for c in copy_part_mid)
        note("    annotation confidence (disposition, tier letter) in the midpoint PHR universe: %s"
             % dict(conf))
        syms = Counter(copy_meta[c]["functional_source_symbol"] for c in copy_part_mid)
        if defn == "source":
            note("    symbols: " + ", ".join("%s=%d" % (s, n) for s, n in tops(syms, 12)))

    out_rows = []
    summary = {}
    regime_params = {
        "midpoint": dict(K=led["phr_elig_mid"], nonphr=led["nonphr_elig_mid"], N=led["univ_mid"],
                         part=copy_part_mid, table=terms_mid, n_supported=REFEREE["supported_terms"]),
        "overlap": dict(K=led["phr_elig_ovl"], nonphr=led["nonphr_elig_ovl"], N=led["univ_mid"],
                        part=copy_part_ovl, table=terms_ovl, n_supported=None),
    }

    for regime in ("midpoint", "overlap"):
        P = regime_params[regime]
        K0, N, part, table = P["K"], P["N"], P["part"], P["table"]
        rs = rows[regime]
        support_col = "primary_support" if regime == "midpoint" else "sensitivity_support"
        nsup = sum(1 for r in rs if r[support_col] == "1")
        if P["n_supported"] is not None and nsup != P["n_supported"]:
            fail("%s supported count %d != %d" % (regime, nsup, P["n_supported"]))

        for defn, field in (("source", "functional_source_id"), ("locus", "gene_name")):
            fam_of_copy = {c: copy_meta[c][field] for c in part}
            K_by_fam = Counter(fam_of_copy.values())

            # p_after for every hypothesis (a_T == 0 -> exactly 1)
            p_after = [1.0] * len(rs)
            lifted = {}
            for i, r in enumerate(rs):
                a = int(r["a_T"])
                if a <= 0:
                    continue
                M = int(r["M_T_genome_burden"])
                counts = Counter()
                for cid in table.get(key(r), {}):
                    counts[fam_of_copy[cid]] += 1
                fam, aF, ntied = largest_family(counts)
                if fam is None:
                    fail("%s %s %s: a_T=%d but no contributor families" % (regime, defn, key(r), a))
                KF = K_by_fam[fam]
                lifted[key(r)] = (fam, aF, ntied, KF)
                p_after[i] = float(hypergeom_upper(a - aF, K0 - KF, M, N))

            # multiplicity on the leave-out p values, identical definitions
            by_col = defaultdict(list)
            for i, r in enumerate(rs):
                by_col[r["collection"]].append(i)
            q_after = [0.0] * len(rs)
            for col, idx in by_col.items():
                sub = bh_qvalues([p_after[i] for i in idx])
                for j, i in enumerate(idx):
                    q_after[i] = sub[j]
            holm_after = holm_pvalues(p_after)

            n_surv = 0
            n_surv_nom = 0
            n_surv_bh = 0
            removed_fams = Counter()
            for i, r in enumerate(rs):
                if r[support_col] != "1":
                    continue
                k = key(r)
                a = int(r["a_T"])
                b = int(r["b_T"])
                c = int(r["c_T"])
                M = int(r["M_T_genome_burden"])
                fam, aF, ntied, KF = lifted[k]
                after = (a - aF, K0 - KF)
                surv_nom = p_after[i] <= 0.05
                surv_bh = q_after[i] <= 0.05
                surv = surv_bh and holm_after[i] <= 0.05
                n_surv_nom += surv_nom
                n_surv_bh += surv_bh
                n_surv += surv
                removed_fams[copy_meta_family_symbol(fam, field, part, copy_meta)] += 1

                # annotation-confidence stratification for this term
                a_self = sum(1 for cid in table.get(k, {})
                             if copy_meta[cid]["source_assignment_disposition"] == "EXACT_SELF")
                a_rel = a - a_self

                fold_before = (a / K0) / (c / P["nonphr"]) if c > 0 else float("inf")
                Ka = after[1]
                fold_after = ((a - aF) / Ka) / (c / P["nonphr"]) if (c > 0 and Ka > 0) else float("inf")

                out_rows.append({
                    "regime": regime,
                    "family_definition": defn,
                    "collection": r["collection"],
                    "relation": r["relation"],
                    "term_id": r["term_id"],
                    "term_name": r["term_name"],
                    "copies_in_phr_a": a,
                    "K_phr_total": K0,
                    "copies_in_nonphr_c": c,
                    "term_bearing_genome_M": M,
                    "universe_N": N,
                    "fold_before": fold_before,
                    "p_before": float(r["p_exact_upper"]),
                    "bh_q_before": float(r["bh_q_within_collection"]),
                    "holm_before": float(r["holm_p_global"]),
                    "largest_family": fam,
                    "largest_family_symbol": (copy_meta_family_symbol(fam, field, part, copy_meta)
                                              if field == "functional_source_id" else fam),
                    "largest_family_copies_in_term": aF,
                    "largest_family_share_of_term_copies": aF / a,
                    "largest_family_share_of_phr_universe": KF / K0,
                    "largest_family_tied_at_max": ntied,
                    "phr_copies_of_family_in_universe": KF,
                    "a_after_leaveout": a - aF,
                    "K_after_leaveout": Ka,
                    "fold_after_leaveout": fold_after,
                    "p_after_leaveout": p_after[i],
                    "bh_q_after_leaveout": q_after[i],
                    "holm_after_leaveout": holm_after[i],
                    "survives_nominal_0.05": int(surv_nom),
                    "survives_bh_within_collection_0.05": int(surv_bh),
                    "survives_primary_rule_after": int(surv),
                    "term_copies_exact_self_annotation": a_self,
                    "term_copies_inherited_annotation": a_rel,
                })
            summary[(regime, defn)] = dict(n_supported=nsup, survive=n_surv,
                                           survive_nom=n_surv_nom, survive_bh=n_surv_bh,
                                           removed=removed_fams)
            note("  %s / %-6s family: of %d supported terms -> survive primary rule after "
                 "leave-largest-family-out: %d (nominal 0.05: %d, BH-within-collection 0.05: %d)"
                 % (regime, defn, nsup, n_surv, n_surv_nom, n_surv_bh))
            note("    largest families removed (by number of supported terms they dominate): "
                 + ", ".join("%s=%d" % (f, n) for f, n in tops(removed_fams, 8)))

    # ---------------- (2b) coherent re-analyses: remove named family everywhere --------
    note("")
    note("== (2b) coherent re-analyses: remove one family, then the k largest families, "
         "from EVERY term ==")
    removal_summary = {}
    for defn, field in (("source", "functional_source_id"), ("locus", "gene_name")):
        for regime in ("midpoint", "overlap"):
            P = regime_params[regime]
            K0, N, part, table = P["K"], P["N"], P["part"], P["table"]
            rs = rows[regime]
            support_col = "primary_support" if regime == "midpoint" else "sensitivity_support"
            fam_of_copy = {c: copy_meta[c][field] for c in part}
            K_by_fam = Counter(fam_of_copy.values())
            aF_by_term = {}
            for r in rs:
                if int(r["a_T"]) > 0:
                    aF_by_term[key(r)] = Counter(fam_of_copy[c] for c in table.get(key(r), {}))
            nsup = sum(1 for r in rs if r[support_col] == "1")
            supported_idx = [i for i, r in enumerate(rs) if r[support_col] == "1"]

            def reanalyse(removed):
                pv = [1.0] * len(rs)
                KF = sum(K_by_fam[f] for f in removed)
                if K0 - KF <= 0:
                    return None
                for i, r in enumerate(rs):
                    a = int(r["a_T"])
                    if a <= 0:
                        continue
                    M = int(r["M_T_genome_burden"])
                    aF = sum(n for f, n in aF_by_term[key(r)].items() if f in removed)
                    pv[i] = float(hypergeom_upper(a - aF, K0 - KF, M, N))
                by_col = defaultdict(list)
                for i, r in enumerate(rs):
                    by_col[r["collection"]].append(i)
                q = [0.0] * len(rs)
                for col, idx in by_col.items():
                    sub = bh_qvalues([pv[i] for i in idx])
                    for j, i in enumerate(idx):
                        q[i] = sub[j]
                return pv, q, holm_pvalues(pv)

            # deterministic family ranking: PHR copy count desc, then family id
            ranked = sorted([f for f in K_by_fam if K_by_fam[f] >= 2],
                            key=lambda f: (-K_by_fam[f], f))
            plans = [("single:" + f, {f}) for f in ranked[:8]]
            for k in range(1, 9):
                if k <= len(ranked):
                    plans.append(("cumulative_top%d:%s" % (k, "+".join(ranked[:k])), set(ranked[:k])))
            for label, removed in plans:
                res = reanalyse(removed)
                if res is None:
                    note("  %-6s %-6s %-40s : PHR copies of removed families exhausted the PHR set"
                         % (regime, defn, label))
                    continue
                pv, q, holm = res
                surv = sum(1 for i in supported_idx if q[i] <= 0.05 and holm[i] <= 0.05)
                best_i = min(supported_idx, key=lambda i: pv[i])
                removal_summary[(regime, defn, label)] = (surv, nsup, pv[best_i], key(rs[best_i]),
                                                          q[best_i], holm[best_i])
                note("  %-6s %-6s %-40s : %3d/%3d supported terms survive | best p=%9.3e q=%9.3e holm=%9.3e (%s)"
                     % (regime, defn, label, surv, nsup, pv[best_i], q[best_i], holm[best_i],
                        rs[best_i]["term_id"]))

    # ---------------- referee's specific spermatid claim ----------------
    note("")
    note("== (2b) the referee's 'spermatid development' claim ==")
    sp_key = None
    for r in rows["midpoint"]:
        if r["term_id"] == "GO:0007286" and r["relation"] == "direct":
            sp_key = key(r)
            break
    if sp_key is None:
        fail("GO:0007286 (spermatid development, direct) not found in TERM_RESULTS midpoint rows")
    sp = frozen_mid_by_key[sp_key]
    note("  GO:0007286 %s: a_T=%s / K=%s ; c_T=%s / non-PHR=%s ; fold=%.3f ; p=%.6e ; holm=%.3e"
         % (sp["term_name"], sp["a_T"], sp["K_phr_total"], sp["c_T"], sp["nonphr_total"],
            float(sp["fold_enrichment"]), float(sp["p_exact_upper"]), float(sp["holm_p_global"])))
    if int(sp["a_T"]) != REFEREE["spermatid_copies"]:
        fail("spermatid a_T %s != referee %d" % (sp["a_T"], REFEREE["spermatid_copies"]))
    if abs(float(sp["fold_enrichment"]) - REFEREE["spermatid_fold"]) > 1.0:
        fail("spermatid fold %.3f not near referee's ~25" % float(sp["fold_enrichment"]))
    sp_copies = sorted(table_copies(terms_mid, sp_key))
    fams = Counter(copy_meta[c]["functional_source_id"] for c in sp_copies)
    syms = Counter(copy_meta[c]["functional_source_symbol"] for c in sp_copies)
    locs = Counter(copy_meta[c]["gene_name"] for c in sp_copies)
    chrs = Counter(copy_meta[c]["seqid"] for c in sp_copies)
    note("  its %d PHR copies: functional sources %s" % (len(sp_copies), dict(fams)))
    note("  symbols %s" % dict(syms))
    note("  loci (gene_name) %s" % dict(locs))
    note("  chromosomes/seqids %s" % dict(chrs))
    note("  distinct genomic loci: %d  -> the referee's 'single duplicon' is %s"
         % (len(locs), "exact" if len(locs) == 1 else "not literally correct (%d loci, one family)" % len(locs)))
    if len(fams) != 1:
        fail("spermatid copies come from %d functional sources, expected 1 (referee claim)" % len(fams))

    # ---------------- (3) invisible content ----------------
    note("")
    note("== (3) invisible PHR gene content ==")
    inelig = led["phr_inelig_mid"]
    note("  PHR gene copies (midpoints): %d ; with a defensible ontology source: %d ; without: %d (%.2f%%)"
         % (led["phr_raw_mid"], led["phr_elig_mid"], inelig, 100.0 * inelig / led["phr_raw_mid"]))
    note("  genome-wide gene copies: %d ; eligible %d ; ineligible %d (%.2f%%)"
         % (led["raw_total"], led["univ_mid"], led["inelig_total"],
            100.0 * led["inelig_total"] / led["raw_total"]))
    elig_copies_coords = {(copy_meta[c]["seqid"], int(copy_meta[c]["start0"]), int(copy_meta[c]["end0"]))
                          for c in copy_part_mid}
    gff_by_coord = {(g["seqid"], g["start0"], g["end0"]): g for g in mid_genes}
    missing = [c for c in elig_copies_coords if c not in gff_by_coord]
    if missing:
        fail("%d eligible PHR copy coordinates from the ledger do not match a GFF3 gene record: %r"
             % (len(missing), missing[:3]))
    bio_402 = Counter(g["biotype"] for g in mid_genes)
    bio_187 = Counter(gff_by_coord[c]["biotype"] for c in elig_copies_coords)
    inelig_coords = [(g["seqid"], g["start0"], g["end0"]) for g in mid_genes
                     if (g["seqid"], g["start0"], g["end0"]) not in elig_copies_coords]
    if len(inelig_coords) != inelig:
        fail("ineligible copy recount %d != ledger %d" % (len(inelig_coords), inelig))
    bio_215 = Counter(gff_by_coord[c]["biotype"] for c in inelig_coords)
    for lab, cc, tot in (("all 402 PHR copies", bio_402, led["phr_raw_mid"]),
                         ("187 ontology-eligible", bio_187, led["phr_elig_mid"]),
                         ("215 ontology-ineligible", bio_215, inelig)):
        note("  biotypes, %s: %s" % (lab, ", ".join("%s=%d (%.1f%%)"
             % (bt, n, 100.0 * n / tot) for bt, n in tops(cc, len(cc)))))
    def merged_bp(coords):
        tot = 0
        for seq in {c[0] for c in coords}:
            iv = sorted((c[1], c[2]) for c in coords if c[0] == seq)
            cur_s = cur_e = None
            for s0, e0 in iv:
                if cur_s is None:
                    cur_s, cur_e = s0, e0
                elif s0 <= cur_e:
                    cur_e = max(cur_e, e0)
                else:
                    tot += cur_e - cur_s
                    cur_s, cur_e = s0, e0
            if cur_s is not None:
                tot += cur_e - cur_s
        return tot

    phr_bp = sum(e - s for seq, ivs in intervals.items() for s, e in ivs)
    elig_bp = merged_bp(list(elig_copies_coords))
    inel_bp = merged_bp(inelig_coords)
    note("  PHR bases: %d ; gene-body bases covered by the 187 eligible copies: %d (%.1f%%) ; "
         "by the 215 ineligible copies: %d (%.1f%%)"
         % (phr_bp, elig_bp, 100.0 * elig_bp / phr_bp, inel_bp, 100.0 * inel_bp / phr_bp))
    for sym in ("DUX4", "SEPTIN14", "WBP1L", "RPL23A", "DDX11", "WASHC1", "CIC", "GTF2I"):
        cc = Counter(copy_meta[c]["gene_biotype"] for c in copy_part_mid
                     if copy_meta[c]["functional_source_symbol"] == sym)
        n = sum(cc.values())
        note("  biotypes of the %d PHR copies of the %s source family: %s"
             % (n, sym, ", ".join("%s=%d" % kv for kv in tops(cc, len(cc)))))
    loc_like = sum(1 for c in inelig_coords if gff_by_coord[c]["gene_name"].startswith("LOC"))
    note("  of the 215 ineligible copies, %d (%.1f%%) carry a LOC* placeholder symbol"
         % (loc_like, 100.0 * loc_like / inelig))
    elig_loc = sum(1 for c in elig_copies_coords if gff_by_coord[c]["gene_name"].startswith("LOC"))
    note("  of the 187 eligible copies, %d (%.1f%%) carry a LOC* placeholder symbol"
         % (elig_loc, 100.0 * elig_loc / led["phr_elig_mid"]))

    # ---------------- write TSV ----------------
    out = opt.output if opt.output else DEFAULT_OUT
    cols = ["regime", "family_definition", "collection", "relation", "term_id", "term_name",
            "copies_in_phr_a", "K_phr_total", "copies_in_nonphr_c", "term_bearing_genome_M",
            "universe_N", "fold_before", "p_before", "bh_q_before", "holm_before",
            "largest_family", "largest_family_symbol", "largest_family_copies_in_term",
            "largest_family_share_of_term_copies", "largest_family_share_of_phr_universe",
            "largest_family_tied_at_max", "phr_copies_of_family_in_universe",
            "a_after_leaveout", "K_after_leaveout", "fold_after_leaveout", "p_after_leaveout",
            "bh_q_after_leaveout", "holm_after_leaveout", "survives_nominal_0.05",
            "survives_bh_within_collection_0.05", "survives_primary_rule_after",
            "term_copies_exact_self_annotation", "term_copies_inherited_annotation"]
    out_rows.sort(key=lambda r: (r["regime"], r["family_definition"], r["holm_before"],
                                 r["collection"], r["term_id"]))
    lines = ["\t".join(cols)]
    for r in out_rows:
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                vals.append("%.6e" % v if (c.startswith("p_") or c.startswith("bh_") or c.startswith("holm_"))
                            else "%.6g" % v)
            else:
                vals.append(str(v))
        lines.append("\t".join(vals))
    blob = "\n".join(lines) + "\n"
    if out == "-":
        sys.stdout.write(blob)
    else:
        with open(out, "w") as fh:
            fh.write(blob)
        note("")
        note("wrote %s (%d data rows)" % (os.path.relpath(out, REPO), len(out_rows)))

    note("")
    note("== wall time: %.1f s ==" % (time.time() - t0))


def copy_meta_family_symbol(fam, field, part, copy_meta):
    if field == "functional_source_id":
        for c in part:
            if copy_meta[c]["functional_source_id"] == fam:
                return copy_meta[c]["functional_source_symbol"]
        return ""
    return fam


def table_copies(table, k):
    return set(table.get(k, {}))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="", help="TSV path ('-' for stdout)")
    ap.add_argument("--self-test", action="store_true")
    opt = ap.parse_args()
    if opt.self_test:
        run_self_test()
        sys.exit(0)
    main()
