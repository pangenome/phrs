#!/usr/bin/env python3
"""GENECONV on duplicon copy alignments (R1-M1, test 3a), with cross-arm classification of fragments.

INPUT: FASTA alignments written by duplicon_monophyly.py --msa-dir (<FAMILY>_<regime>_{units,arrays}.fa),
names <arm>_<sample>_<hap>[_u<i>], e.g. 4q_HG00099_1_u3; codes ACGT, '-' gap, 'N' missing.
STEPS per alignment: keep sequences with >= --min-complete non-N columns and within --max-cons-dist of their own
arm consensus (a divergent outlier copy otherwise yields fragments wherever it is not divergent), drop columns that are N or '-'
in every kept sequence, run GENECONV (patched 1.81a, scripts/conversion/geneconv_1.81a_va_copy.patch)
in its own directory, parse the .frags table, and classify every globally significant inner fragment
(GI, simulation P corrected over all pairs < 0.05) by pair: different arms (an exchange between chromosome
ends), same arm in different haplotypes, or same arm in the same haplotype (within one array).
OUTPUT: <out>/summary.tsv (one row per alignment) and <out>/fragments.tsv (one row per GI fragment).
"""
import argparse, os, subprocess, sys, time
from collections import Counter

FRAG_COLS = ["type", "names", "sim_p", "bc_ka_p", "begin", "end", "len", "npoly", "ndif", "totdifs", "mismat_pen"]


def log(msg):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seqs:
                raise SystemExit(f"FATAL: duplicate name {name} in {path}")
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    seqs = {k: "".join(v) for k, v in seqs.items()}
    if len({len(v) for v in seqs.values()}) != 1:
        raise SystemExit(f"FATAL: sequences of different length in {path}")
    return seqs


def pair_class(a, b):
    pa, pb = a.split("_"), b.split("_")
    if pa[0] != pb[0]:
        return "cross_arm"
    return "same_arm_same_hap" if pa[1:3] == pb[1:3] else "same_arm_diff_hap"


def run_one(msa, geneconv, outdir, min_complete, gopts, reuse=False, max_cons=0.0):
    tag = os.path.basename(msa).rsplit(".", 1)[0] + ("_" + gopts.replace("/", "").replace(" ", "")
                                                     if gopts else "")
    wd = os.path.join(outdir, tag)
    os.makedirs(wd, exist_ok=True)
    seqs = read_fasta(msa)
    L = len(next(iter(seqs.values())))
    keep = {k: v for k, v in seqs.items() if sum(c != "N" for c in v) >= min_complete * L}
    log(f"{tag}: {len(seqs)} sequences x {L} columns; kept {len(keep)} with >= {min_complete:.2f} non-N, "
        f"dropped {len(seqs) - len(keep)} (missing data)")
    row = {"alignment": tag, "n_in": len(seqs), "n_kept": len(keep), "columns": L,
           "arms_kept": ";".join(f"{a}:{n}" for a, n in sorted(Counter(k.split("_")[0] for k in keep).items()))}
    if max_cons > 0 and keep:
        cons = {}
        for arm in {k.split("_")[0] for k in keep}:
            vs = [v for k, v in keep.items() if k.split("_")[0] == arm]
            cons[arm] = []
            for j in range(L):
                c = Counter(v[j] for v in vs if v[j] in "ACGT")
                cons[arm].append(c.most_common(1)[0][0] if c else "N")
        far = {}
        for k, v in keep.items():
            cs = cons[k.split("_")[0]]
            n = sum(1 for j in range(L) if v[j] in "ACGT" and cs[j] != "N")
            d = sum(1 for j in range(L) if v[j] in "ACGT" and cs[j] != "N" and v[j] != cs[j]) / max(1, n)
            if d > max_cons:
                far[k] = d
        for k in far:
            del keep[k]
        log(f"{tag}: dropped {len(far)} copies farther than {max_cons} from their own arm consensus: "
            + ", ".join(f"{k} {d:.4f}" for k, d in sorted(far.items())))
        row["n_far_dropped"] = len(far)
        row["n_kept"] = len(keep)
    if len({k.split("_")[0] for k in keep}) < 2:
        log(f"{tag}: fewer than 2 arms after filtering, not run")
        return dict(row, status="not_run_lt2_arms"), []
    cols = [j for j in range(L) if any(v[j] not in "N-" for v in keep.values())]
    fa = os.path.join(wd, "aln.fasta")
    with open(fa, "w") as fh:
        for k, v in keep.items():
            fh.write(f">{k}\n{''.join(v[j] for j in cols)}\n")
    cmd = [geneconv, "aln.fasta", "/w123", "/lp", "/es"] + gopts.split()
    t0 = time.time()
    if reuse and os.path.exists(os.path.join(wd, "aln.frags")):
        log(f"{tag}: --reuse, parsing the existing {wd}/aln.frags")
    else:
        res = subprocess.run(cmd, cwd=wd, capture_output=True, text=True)
        open(os.path.join(wd, "geneconv.log"), "w").write(res.stdout + res.stderr)
        if res.returncode != 0:
            raise SystemExit(f"FATAL: {' '.join(cmd)} exit {res.returncode} in {wd}")
    frags_path = os.path.join(wd, "aln.frags")
    npoly = None
    frags = []
    for line in open(frags_path):
        if line.startswith("#"):
            if "polymorphisms to be permuted" in line:
                npoly = int(line.split()[1].replace(",", ""))
            continue
        f = line.split()
        if len(f) < len(FRAG_COLS):
            continue
        d = dict(zip(FRAG_COLS, f))
        nm = d["names"].split(";")
        d["pair_class"] = pair_class(*nm) if len(nm) == 2 else "outer"
        d["alignment"] = tag
        frags.append(d)
    gi = [d for d in frags if d["type"] == "GI"]
    c = Counter(d["pair_class"] for d in gi)
    row.update(status="ok", columns_used=len(cols), polymorphisms=npoly, n_GI=len(gi),
               GI_cross_arm=c["cross_arm"], GI_same_arm_diff_hap=c["same_arm_diff_hap"],
               GI_same_arm_same_hap=c["same_arm_same_hap"], seconds=round(time.time() - t0, 1))
    log(f"{tag}: {len(cols)} columns, {npoly} polymorphisms, GI {len(gi)} (cross-arm {c['cross_arm']}, "
        f"same arm other haplotype {c['same_arm_diff_hap']}, same array {c['same_arm_same_hap']}); "
        f"{row['seconds']} s")
    return row, gi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("msas", nargs="+")
    ap.add_argument("--geneconv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-complete", type=float, default=0.9)
    ap.add_argument("--gopts", action="append", default=None,
                    help="extra GENECONV options per run, e.g. '' and '/g1' (repeatable)")
    ap.add_argument("--reuse", action="store_true", help="parse an existing aln.frags instead of rerunning")
    ap.add_argument("--max-cons-dist", type=float, default=0.02,
                    help="drop copies whose p-distance to their own arm consensus exceeds this (0 = keep all)")
    a = ap.parse_args()
    gopts = a.gopts if a.gopts is not None else [""]
    os.makedirs(a.out, exist_ok=True)
    log(f"inputs: {len(a.msas)} alignments; geneconv {a.geneconv}; min_complete {a.min_complete}; "
        f"runs per alignment: {gopts}")
    rows, allf = [], []
    for m in a.msas:
        for g in gopts:
            r, f = run_one(m, a.geneconv, a.out, a.min_complete, g, a.reuse, a.max_cons_dist)
            rows.append(r)
            allf += f
    keys = ["alignment", "status", "n_in", "n_far_dropped", "n_kept", "arms_kept", "columns", "columns_used", "polymorphisms",
            "n_GI", "GI_cross_arm", "GI_same_arm_diff_hap", "GI_same_arm_same_hap", "seconds"]
    with open(os.path.join(a.out, "summary.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "NA")) for k in keys) + "\n")
    fk = ["alignment", "pair_class"] + FRAG_COLS
    with open(os.path.join(a.out, "fragments.tsv"), "w") as fh:
        fh.write("\t".join(fk) + "\n")
        for d in allf:
            fh.write("\t".join(d[k] for k in fk) + "\n")
    log(f"wrote {a.out}/summary.tsv ({len(rows)} rows) and fragments.tsv ({len(allf)} GI fragments)")


if __name__ == "__main__":
    main()
