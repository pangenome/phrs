#!/usr/bin/env python3
"""Audit the numbers quoted in the reviewer-response material against their source artifacts.

Motivation: three defects in this project were all the same class - a number quoted in prose
diverging from the file that produced it (mixed-resolution B/W 0.027 vs 0.0031; hap1 vs
per-chromosome maternal copies; zygotene vs pachytene peak). This script recomputes the
canonical value from the source artifact for every headline statistic and checks that the
response documents contain it.

Usage: python3 scripts/ci/audit_response_numbers.py [--docs RESPONSE_LETTER.md ...]
Exit code 1 if any documented value is missing or disagrees with its source.
"""
from __future__ import annotations
import argparse, csv, json, math, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def wilson(k, n, z=1.959963985):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h


def collect():
    """canonical statistic -> (formatted string that must appear, source file)"""
    out = {}
    # pedigree fraction + Wilson CI (recomputed here, so the doc cannot drift from the counts)
    for k, n, tag in [(494, 538, "published"), (293, 310, "corrected")]:
        p, lo, hi = wilson(k, n)
        out[f"pedigree_wilson_{tag}"] = (f"{p:.3f}", f"[{lo:.3f}, {hi:.3f}]", "counts in RESPONSE_LETTER.md")
    # D-M4 null
    j = json.loads((ROOT / "paper_prep/synthesis/ANALYSIS_D_M4.json").read_text())
    txt = json.dumps(j, indent=1)
    for needle in ("0.6772", "0.8396"):
        out[f"dm4_null_{needle}"] = (needle, "", "ANALYSIS_D_M4.json")
    # D-M5 mouse arm-level Mantel
    m = json.loads((ROOT / "scripts/hic/mantel_mouse_d_m5.json").read_text())
    for stage in ("zygotene", "pachytene"):
        v = m.get(f"{stage}_20kb", {}).get("rho")
        if v is not None:
            out[f"dm5_{stage}_20kb"] = (f"{abs(v):.2f}", "", "mantel_mouse_d_m5.json")
    # D-M9 character bootstrap
    d9 = {}
    for line in (ROOT / "scripts/cladistics/d_m9/d_m9_support_summary.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        d9[f[0]] = (f[4], f[5])
    for clade in ("PAR1", "PAR2", "P10_18", "ACRO_p", "TIGHT_q", "DUX4"):
        if clade in d9:
            out[f"dm9_{clade}_upgma"] = (f"{float(d9[clade][1]):.1f}", "", "d_m9_support_summary.tsv")
    # Q3 arm-type asymmetry
    pq = list(csv.DictReader(open(ROOT / "scripts/cladistics/peerq3_arm_type_asymmetry.tsv"), delimiter="\t"))
    c6 = next(r for r in pq if r["community"] == "C6")
    c7 = next(r for r in pq if r["community"] == "C7")
    out["q3_C6_QQ"] = (c6["Q-Q"], "", "peerq3_arm_type_asymmetry.tsv")
    out["q3_C7_PP"] = (c7["P-P"], "", "peerq3_arm_type_asymmetry.tsv")
    # module architecture
    ori = list(csv.DictReader(open(ROOT / "scripts/cladistics/arm_module_orientation.tsv"), delimiter="\t"))
    ordr = list(csv.DictReader(open(ROOT / "scripts/cladistics/arm_module_order.tsv"), delimiter="\t"))
    tot = sum(int(r["same_bp"]) + int(r["opposite_bp"]) for r in ori)
    same = sum(int(r["same_bp"]) for r in ori) * 100.0 / tot if tot else 0
    out["arch_orientation_pct"] = (f"{same:.1f}", "", "arm_module_orientation.tsv")
    js = [float(r["repertoire_jaccard"]) for r in ordr if r["community_relation"] == "same"]
    jd = [float(r["repertoire_jaccard"]) for r in ordr if r["community_relation"] == "diff"]
    med = lambda v: sorted(v)[len(v) // 2]
    out["arch_repertoire_same"] = (f"{med(js):.2f}", "", "arm_module_order.tsv")
    out["arch_repertoire_diff"] = (f"{med(jd):.2f}", "", "arm_module_order.tsv")
    # F_ST matched baseline
    mt = (ROOT / "scripts/popgen/matched_fst_d_m6.tsv")
    if mt.exists():
        rows = list(csv.DictReader(open(mt), delimiter="\t"))
        for r in rows:
            if r.get("superpop_pair", "").startswith("AFR-"):
                out[f"fst_{r['superpop_pair']}"] = (f"{float(r['ratio_subtelo_over_matched']):.2f}", "", mt.name)
    # threshold sensitivity: chr9q hap1 exclusion coverage at every f
    th = list(csv.DictReader(open(ROOT / "scripts/pedigree/threshold_sensitivity_9q3q.tsv"), delimiter="\t"))
    cov = {r["chr9q_hap1_cov"] for r in th}
    out["threshold_hap1_excluded_everywhere"] = ("yes" if all(float(c) < 10 for c in cov) else "NO", "", "threshold_sensitivity_9q3q.tsv")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", nargs="*", default=[
        "RESPONSE_LETTER.md", "RESPONSE_STATS_AND_BIB.md", "RESPONSE_Q3_DRAFT.md",
        "RESPONSE_Q1_Q2_DRAFT.md"])
    a = ap.parse_args()
    canon = collect()
    texts = {d: (ROOT / d).read_text() for d in a.docs if (ROOT / d).exists()}
    missing = []
    print(f"{'statistic':34s} {'value':>10s} {'where quoted':>28s}  source")
    for name, (val, extra, src) in canon.items():
        def found(text, v):
            if v in text:
                return True
            # numeric tolerance: the prose may round (0.677 vs 0.6772)
            try:
                target = float(v)
            except ValueError:
                return False
            for m in re.finditer(r"[-+]?\d+\.\d+", text):
                if abs(float(m.group(0)) - target) <= 0.005 * max(1.0, abs(target)):
                    return True
            return False

        hits = [d for d, t in texts.items() if found(t, val) and (not extra or extra in t)]
        status = "ok" if hits else "MISSING"
        if not hits:
            missing.append((name, val, extra, src))
        print(f"{name:34s} {val:>10s} {(hits[0] if hits else '-'):>28s}  {src}"
              + (f"  [{extra}]" if extra else ""))
    print()
    if missing:
        print("VALUES NOT FOUND IN THE DOCUMENTS (fix the prose or the source):")
        for n, v, e, s in missing:
            print(f"  {n}: expected '{v}' (from {s})")
        return 1
    print("all canonical values appear in the response documents")
    return 0


if __name__ == "__main__":
    sys.exit(main())
