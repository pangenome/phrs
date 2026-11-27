#!/usr/bin/env python3
"""KaryoScope view of the flanks with candidate joins (flank_joins_read_check.tsv) plus longer ends (records of --fasta not
in the table: contig:start-end, q ends). One track per flank, telomere on the left: the KaryoScope HKS_human_CHM13_v2
chromosome of origin (smoothed, KaryoScope colors; labels shared by several chromosomes in light blue), telomeric repeat
(subtelomeric set: canonical or noncanonical telomere) in black and TAR1 in red with a black edge (drawn from the BED
intervals merged within 500 bp, blocks >= 500 bp, drawn at least 1/1000 of the axis wide), N runs hatched (KaryoScope fills them with the labels around them). The square left of the track has the color of the flank's own chromosome. Above: candidate joins (triangle
= N run, diamond = telomeric array) colored by the long-read verdict. Below: the paper PHR. Each bin gets the label of
its bp majority. Pages: wrong join, undecided, correct, deep ends."""
import argparse, collections, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--bed-prefix", required=True, help="<prefix>.<feature_set>.smoothed.bed")
ap.add_argument("--colors", required=True)
ap.add_argument("--table", required=True)
ap.add_argument("--fasta", required=True, help="FASTA of the table flanks and the deep ends")
ap.add_argument("--out", required=True, help="output prefix (.pdf and .pageN.png)")
a = ap.parse_args()
SETS = ["chromosome", "subtelomeric"]
TEL = {"canonical_telomere", "noncanonical_telomere"}
VCOL = {"spanned": "#1a9850", "end_at_join": "#d73027", "no_reads_tel": "#d73027", "no_reads_far": "#fdae61",
        "not_spanned": "#fdae61"}
col = {}
for l in open(a.colors):
    fs, ft, c = l.rstrip("\n").split("\t")
    if fs != "feature_set":
        col[fs, ft] = c
for fs in SETS:
    col[fs, "novel"] = "#FFFFFF"; col[fs, "."] = "#FFFFFF"

lines = open(a.table).read().splitlines(); hdr = lines[0].split("\t")
tab = {x[0]: dict(zip(hdr, x)) for x in (l.split("\t") for l in lines[1:])}
fa, name, seq = {}, None, []
for l in open(a.fasta):
    if l.startswith(">"):
        if name:
            fa[name] = "".join(seq)
        name, seq = l[1:].split()[0], []
    else:
        seq.append(l.strip().upper())
fa[name] = "".join(seq)
missing = set(tab) - set(fa)
if missing:
    raise SystemExit(f"FATAL: {len(missing)} table flanks not in --fasta")
deep = {k: v for k, v in fa.items() if k not in tab}
nruns = {k: [(m.start(), m.end()) for m in re.finditer("N+", v)] for k, v in fa.items()}
NAME = re.compile(r"^([^#]+#[^#]+)#[^:]+:(\d+)-(\d+)(?:_(chr[0-9XY]+)_(parm|qarm))?$")
DEEP_LAB = {"HG03521#2": "chr3q", "HG02922#1": "chrXq"}
length = {nm: len(s) for nm, s in fa.items()}
print(f"{len(tab)} flanks in the table, {len(deep)} deep ends")

labs = {fs: {} for fs in SETS}
iv = {fs: collections.defaultdict(list) for fs in SETS}
for fs in SETS:
    n = 0
    for l in open(f"{a.bed_prefix}.{fs}.smoothed.bed"):
        c, s, e, f = l.split("\t")[:4]
        if c in length:
            iv[fs][c].append((int(s), int(e), labs[fs].setdefault(f.strip(), len(labs[fs])))); n += 1
    miss = [k for k in labs[fs] if (fs, k) not in col]
    print(f"{fs}: {n} intervals, {len(labs[fs])} labels, without a KaryoScope color: {miss or 'none'}")
    if miss:
        raise SystemExit(f"FATAL: labels without a color in {fs}: {miss}")


def binned(nm, fs, B):
    L = length[nm]; nl = len(labs[fs]) + 1; arr = np.full(L, nl - 1, np.int32)  # last id = not annotated
    for s, e, k in iv[fs][nm]:
        arr[s:e] = k
    nb = L // B; arr = arr[:nb * B]
    cnt = np.bincount(np.arange(nb * B) // B * nl + arr, minlength=nb * nl).reshape(nb, nl)
    b = cnt.argmax(1)
    q = nm.endswith("_qarm") or nm in deep
    return b[::-1] if q else b  # telomere on the left


inv = {fs: {v: k for k, v in labs[fs].items()} for fs in SETS}


def segments(b, B, keep):
    runs = np.flatnonzero(np.diff(np.r_[-1, b, -1]) != 0); segs = collections.defaultdict(list)
    for i0, i1 in zip(runs[:-1], runs[1:]):
        if keep(b[i0]):
            segs[b[i0]].append((i0 * B, (i1 - i0) * B))
    return segs


def draw(ax, nm, y, B, label, own, phr_bp, marks):
    used = set(); h = 0.5
    for k, v in segments(binned(nm, "chromosome", B), B, lambda k: True).items():
        f = inv["chromosome"].get(k, "."); ax.broken_barh(v, (y - h / 2, h), facecolors=col["chromosome", f], edgecolor="none")
        used.add(f)
    L = length[nm]; q = nm.endswith("_qarm") or nm in deep; wmin = (8e6 if nm in deep else 5e5) / 1000
    for lab, fc, ec in (("TAR1", "#FF0000", "black"), ("TEL", "black", "none")):
        k = labs["subtelomeric"]; ids = {k[x] for x in (TEL if lab == "TEL" else {lab}) if x in k}
        bl = []  # intervals merged within 500 bp, blocks >= 500 bp (most k-mer hits are < 100 bp fragments)
        for s, e, i in sorted(x for x in iv["subtelomeric"][nm] if x[2] in ids):
            if bl and s - bl[-1][1] <= 500:
                bl[-1][1] = max(bl[-1][1], e)
            else:
                bl.append([s, e])
        v = [((L - e if q else s), max(e - s, wmin)) for s, e in bl if e - s >= 500]
        ax.broken_barh(v, (y - h / 2, h), facecolors=fc, edgecolor=ec, linewidth=0.4)
    for s0, e0 in nruns[nm]:
        ax.add_patch(plt.Rectangle((L - e0 if q else s0, y - h / 2), e0 - s0, h, fc="white", ec="black", lw=0.3, hatch="/////"))
    if phr_bp:
        ax.plot([0, phr_bp], [y - 0.38, y - 0.38], color="black", lw=1.0, solid_capstyle="butt")
    for kind, x, v in marks:
        ax.plot(x, y + 0.36, marker="v" if "N" in kind else "D", ms=5 if "N" in kind else 4, color=VCOL[v], mec="black", mew=0.4)
    ax.add_patch(plt.Rectangle((-0.012, y - h / 2), 0.008, h, transform=ax.get_yaxis_transform(), clip_on=False,
                               fc=col["chromosome", own], ec="black", lw=0.3))
    ax.text(-0.016, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.5)
    return used


def short(nm):
    m = NAME.match(nm)
    return f"{m.group(1)} {m.group(4)[3:]}{m.group(5)[0]}" if m.group(4) else f"{m.group(1)} {DEEP_LAB[m.group(1)][3:]}"


chrkey = lambda nm: (lambda c: (int(c) if c.isdigit() else 23 + "XY".index(c)))(NAME.match(nm).group(4)[3:])
groups = [("wrong", "Flanks with a wrong join (reads stop in telomere at it, or no read on its telomere side)"),
          ("undecided", "Flanks with an unresolved join (no read spans it, nothing shows it is wrong)"),
          ("correct", "Flanks whose joins long reads span")]
pdf = PdfPages(f"{a.out}.pdf"); page = 0
for g, title in groups + [("deep", "Two chimeras past 500 kb (8 Mb from the telomere)")]:
    if g == "deep":
        rows = sorted(deep); B = 4000; xmax = 8e6
    else:
        rows = sorted((nm for nm in tab if tab[nm]["verdict"] in (g, g + "_tiny")), key=lambda nm: (chrkey(nm), nm)); B = 250; xmax = 5e5
    W, rows_h, top, leg_h = 12.5, 0.27 * len(rows) + 0.5, 0.75, 1.0  # inches
    H = rows_h + top + leg_h + 0.45
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([3.0 / W, (leg_h + 0.45) / H, (W - 3.2) / W, rows_h / H])
    used = set()
    for i, nm in enumerate(rows):
        y = -i * 1.0
        if g == "deep":
            sq = deep[nm]; L = len(sq)
            marks = [("N", L - mm.end(), "end_at_join") for mm in re.finditer("N+", sq)]
            lab = f"{short(nm)}  HPRC T2T: {'scaffold' if nm.startswith('HG03521') else 'no'}"
            used |= draw(ax, nm, y, B, lab, DEEP_LAB[NAME.match(nm).group(1)][:-1], 0, marks)
        else:
            r = tab[nm]
            marks = [(j.split("@")[0], float(j.split("@")[1].split(":")[0]) * 1000, j.split(":")[-1]) for j in r["joins(kind@kb:len:class)"].split(";")]
            lv = {"contig": "T2T contig", "scaffold": "T2T scaffold"}.get(r["hprc_t2t_level"], "not T2T")
            tiny = "  (piece < 10 kb)" if r["verdict"] == "wrong_tiny" else ""
            used |= draw(ax, nm, y, B, f"{short(nm)}  {lv}  PHR {int(r['paper_phr_bp']) // 1000} kb{tiny}", NAME.match(nm).group(4),
                         int(r["paper_phr_bp"]), marks)
    ax.set_xlim(0, xmax); ax.set_ylim(-(len(rows) - 1) - 0.6, 0.6); ax.set_yticks([])
    ticks = np.arange(0, xmax + 1, xmax / 10)
    ax.set_xticks(ticks); ax.set_xticklabels([f"{t / 1e6:g} Mb" if g == "deep" else f"{t / 1e3:g}" for t in ticks], fontsize=7)
    ax.tick_params(axis="x", top=True, labeltop=True)
    ax.set_xlabel("distance from the telomere (kb)" if g != "deep" else "distance from the telomere", fontsize=8)
    for s in ("left", "right"):
        ax.spines[s].set_visible(False)
    page += 1
    ax.set_title(f"{title}: {len(rows)}\nKaryoScope chromosome of origin, telomeric repeat black, TAR1 red, N runs hatched; square = the flank's own chromosome", fontsize=9, loc="left")
    ck = lambda f: (int(f[3:]) if f[3:].isdigit() else 23 + "XY".index(f[3:])) if f.startswith("chr") else 99
    spec = sorted((f for f in used if f.startswith("chr")), key=ck)
    hs = [Patch(fc=col["chromosome", f], ec="black", lw=0.3, label=f) for f in spec]
    hs += [Patch(fc="#B0C4DE", ec="black", lw=0.3, label="shared by several chromosomes"), Patch(fc="white", ec="black", lw=0.3, label="novel (not in CHM13)"),
           Patch(fc="black", ec="black", lw=0.3, label="telomeric repeat"), Patch(fc="#FF0000", ec="black", lw=0.4, label="TAR1")]
    hs += [Line2D([], [], marker="v", ls="", color=c, mec="black", label=t) for t, c in
           (("join spanned by reads", "#1a9850"), ("join wrong", "#d73027"), ("join unresolved", "#fdae61"))]
    hs += [Line2D([], [], marker="D", ls="", color="white", mec="black", label="telomeric array (triangle: N run)"),
           Patch(fc="white", ec="black", lw=0.3, hatch="/////", label="N run"),
           Line2D([], [], color="black", lw=1.2, label="paper PHR")]
    fig.legend(handles=hs, loc="lower left", ncol=9, fontsize=6.5, frameon=False, bbox_to_anchor=(0.01, 0.005), handlelength=1.2, columnspacing=0.8)
    pdf.savefig(fig); fig.savefig(f"{a.out}.page{page}.png", dpi=110); plt.close(fig)
    print(f"page {page}: {g}, {len(rows)} rows, bin {B} bp")
pdf.close()
