#!/usr/bin/env python3
"""IGV-style screenshots of the 9q/3q tract: DAD and SELF reads on DAD and SELF assemblies, v1.1 and v1.3.1.

usage: screenshot.py DATA_DIR OUT_DIR [FLANK_BP]
  DATA_DIR: copy of check/igv_slices80 (igv_slices.sbatch, flank 80000); FLANK_BP default 60000 each side of the Fig5 tract
Writes 9q3q_<READS>reads_on_<ASM>.<v>.png (8 files), with flagger/nucflag tracks for v1.1 (flagger_v11/). A base is coloured only where >= 3 MAPQ >= 20 reads
show a base and >= 80% of them share one base that differs from the assembly.
"""
import sys
import numpy as np, pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

d, out = sys.argv[1:3]; FL = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
BASEC = {"A": "#00A000", "C": "#0000FF", "G": "#D17105", "T": "#FF0000"}
NAME = {"DAD": "DAD (PAN011)", "SELF": "SELF (PAN027)"}
CTG = {"DAD": "PAN011#2#chr9", "SELF": "PAN027#2#chr9"}
QCN = {"DAD": ("PAN011", "PAN011.chr9.haplotype2"), "SELF": ("PAN027", "PAN027.chr9.paternal")}
QCC = {"Dup": "#E41A1C", "Col": "#984EA3", "Err": "#000000", "Hap": "#4DAF4A", "MISJOIN": "#A65628", "COLLAPSE": "#F781BF", "COLLAPSE_OTHER": "#FF7F00", "COLLAPSE_VAR": "#FFD92F"}
QCDIR = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "flagger_v11")
TRACT = {("DAD", "v11"): (133702042, 133722041), ("DAD", "v131"): (133702516, 133722515),
         ("SELF", "v11"): (136168001, 136188000), ("SELF", "v131"): (136175732, 136195731)}

def disagree(bam, ctg, s, e, seq):
    """columns where MAPQ>=20 reads reject the assembly base"""
    cov = np.array(bam.count_coverage(ctg, s, e, quality_threshold=0, read_callback=lambda r: r.mapping_quality >= 20 and not r.is_secondary))
    flag = {}
    for i, rb in enumerate(seq):
        if rb not in "ACGT": continue
        n, b = max((cov[k, i], x) for k, x in enumerate("ACGT") if x != rb); nr = cov["ACGT".index(rb), i]
        if nr + n >= 3 and n / (nr + n) >= 0.8: flag[s + i] = b
    return flag, cov.sum(0)

def bases(r, cols):
    """read base at sorted reference columns"""
    o, ref, q, k = {}, r.reference_start, 0, 0
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            while k < len(cols) and cols[k] < ref: k += 1
            while k < len(cols) and cols[k] < ref + n: o[cols[k]] = r.query_sequence[q + cols[k] - ref]; k += 1
            ref += n; q += n
        elif op in (2, 3): ref += n
        elif op in (1, 4): q += n
    return o

def shot(reads_of, asm, v):
    ctg = CTG[asm]; t0, t1 = TRACT[(asm, v)]
    fa = pysam.FastaFile(f"{d}/{ctg.rsplit('#', 1)[0]}.chr9.{v}.fa.gz"); L = fa.get_reference_length(ctg)
    s, e = max(0, t0 - 1 - FL), min(L, t1 + FL); seq = fa.fetch(ctg, s, e).upper()
    bam = pysam.AlignmentFile(f"{d}/{reads_of}reads_on_{asm}.{v}.bam")
    flag, tot = disagree(bam, ctg, s, e, seq); cols = sorted(flag)
    reads = [r for r in bam.fetch(ctg, s, e) if not r.is_secondary and r.query_sequence]
    rows, ends = [], []
    for grp in (sorted([r for r in reads if r.mapping_quality >= 20], key=lambda r: r.reference_start),
                sorted([r for r in reads if r.mapping_quality < 20], key=lambda r: r.reference_start)):
        first = len(ends)
        for r in grp:  # greedy packing, MAPQ>=20 rows on top
            for i in range(first, len(ends)):
                if ends[i] + 200 < r.reference_start: ends[i] = r.reference_end; rows[i].append(r); break
            else: ends.append(r.reference_end); rows.append([r])
    nhi = sum(1 for row in rows if row[0].mapping_quality >= 20); n = len(rows)
    hr = max(1.2, 0.028 * n)  # tracks panel fixed in inches, reads panel by rows
    fig, (at, ax) = plt.subplots(2, 1, figsize=(14, 2.3 + hr), sharex=True, gridspec_kw={"height_ratios": [1.75, hr], "hspace": 0.03})
    lab_x = s - (e - s) * 0.005
    sm, qn = QCN[asm]  # QC tracks, v1.1 only
    for j, (lab, fn) in enumerate((("flagger HiFi", f"flagger.{sm}.hifi.bed"), ("flagger ONT", f"flagger.{sm}.ONT.bed"), ("nucflag HiFi", f"nucflag.hifi.{sm}.bed"))):
        yy = 8.4 - j * 1.2
        at.text(lab_x, yy + 0.45, lab, fontsize=7.5, va="center", ha="right", color="#555555")
        if v != "v11": at.text(s, yy + 0.45, " not available for v1.3.1", fontsize=7.5, va="center", color="#555555"); continue
        for l in open(f"{QCDIR}/{fn}"):
            f = l.rstrip("\n").split("\t")
            if f[0] == qn and int(f[2]) > s and int(f[1]) < e:
                a, b = max(int(f[1]), s), min(int(f[2]), e)
                at.add_patch(Rectangle((a, yy), b - a, 0.9, fc=QCC.get(f[3], "#555555"), ec="none"))
                at.text(a, yy + 0.45, f" {f[3]}", fontsize=7.5, va="center", color="white", clip_on=True)
    for l in open(f"{d}/regions.{asm}.{v}.bed"):  # regions of interest
        c, a, b, nm = l.rstrip("\n").split("\t"); fig5 = "Fig5" in nm; yy = 4.7 if fig5 else 3.5
        at.add_patch(Rectangle((int(a), yy), int(b) - int(a), 0.9, fc="#E07B00" if fig5 else "#7B3294", ec="none"))
        at.text(lab_x, yy + 0.45, "9q/3q call (Fig5)" if fig5 else "wrong in DAD asm", fontsize=7.5, va="center", ha="right", color="#555555")
        at.text(max(int(a), s), yy + 0.45, " " + nm, fontsize=7.5, va="center", color="white", clip_on=True)
    xs = np.arange(s, e); top = max(tot.max(), 1)  # coverage track
    at.fill_between(xs, 0, 2.8 * tot / top, color="#BDBDBD", lw=0, step="mid")
    for p, b in flag.items(): at.plot([p, p], [0, 2.8 * tot[p - s] / top], color=BASEC[b], lw=0.7)
    at.text(lab_x, 1.4, f"coverage (max {top})", fontsize=7.5, va="center", ha="right", color="#555555")
    at.set_ylim(0, 9.5); at.set_yticks([]); at.spines[["left", "right", "top", "bottom"]].set_visible(False); at.tick_params(bottom=False)
    for y, row in enumerate(rows):
        yy = n - 1 - y
        for r in row:
            hi = r.mapping_quality >= 20
            ax.plot([max(r.reference_start, s), min(r.reference_end, e)], [yy, yy], color="#8C8C8C" if hi else "#DADADA", lw=1.2, solid_capstyle="butt")
            for p, b in bases(r, cols).items():
                if b != seq[p - s]: ax.plot([p, p], [yy - 0.42, yy + 0.42], color=BASEC.get(b, "#000000"), lw=0.7, alpha=1 if hi else 0.4)
    if nhi < n: ax.axhline(n - nhi - 0.5, color="#555555", lw=0.6, ls="--")
    ax.text(lab_x, n - nhi / 2, "MAPQ >= 20", fontsize=7.5, va="center", ha="right", color="#555555")
    if nhi < n: ax.text(lab_x, (n - nhi) / 2, "MAPQ < 20", fontsize=7.5, va="center", ha="right", color="#555555")
    ax.set_xlim(s, e); ax.set_ylim(-1, n); ax.set_yticks([]); ax.spines[["left", "right", "top"]].set_visible(False)
    tk = np.linspace(s, e, 7).astype(int); ax.set_xticks(tk); ax.set_xticklabels([f"{x:,}" for x in tk], fontsize=8)
    ax.set_xlabel(f"{ctg} (the chr9q end is at the right)", fontsize=9)
    ver = "v1.1" if v == "v11" else "v1.3.1"; kind = "same person" if reads_of == asm else "cross"
    at.set_title(f"{NAME[reads_of]} ONT reads on the {NAME[asm]} assembly, {ver} ({kind})\n"
                 f"{len(reads)} reads ({sum(r.mapping_quality >= 20 for r in reads)} with MAPQ >= 20); "
                 f"{len(flag)} positions where the reads disagree with the assembly; window: {(t0 - 1 - s) // 1000} kb before the Fig5 tract to " + ("the chr9q end" if e == L else f"{(e - t1) // 1000} kb after it"),
                 fontsize=10.5, loc="left")
    hs = [Patch(fc="#E07B00", label="9q/3q exchange called in Fig5")] + ([Patch(fc="#7B3294", label="wrong in the DAD assembly: copy of chr1p h1")] if asm == "DAD" else [])
    qc = [Patch(fc=QCC[k], label=k) for k in ("Dup", "Col", "Err", "MISJOIN")] if v == "v11" else []
    ax.legend(handles=hs + [Patch(fc=BASEC[x], label=x) for x in "ACGT"] + qc, loc="upper left", bbox_to_anchor=(0, -0.45 / hr), ncol=10, fontsize=8.5, frameon=False)
    fn = f"{out}/9q3q_{reads_of}reads_on_{asm}.{v}.png"; fig.savefig(fn, dpi=150, facecolor="white", bbox_inches="tight"); plt.close(fig)
    sys.stderr.write(f"{fn}\treads={len(reads)}\tdisagree={len(flag)}\twindow={ctg}:{s + 1}-{e}\n")

for v in ("v11", "v131"):
    for asm in ("DAD", "SELF"):
        for reads_of in ("DAD", "SELF"): shot(reads_of, asm, v)
