#!/usr/bin/env python3
"""Extended Data Fig. 4: pairwise dot plots of the telomere-first 60 kb ends at the chr9q candidate locus, on the corrected
assemblies plus the released DAD chr9q haplotype 2.

usage: ed4_pairwise_dotplot.py FASTA PAF TRACT_START TRACT_END OUT_PREFIX
  FASTA: the six 60 kb windows (q ends reverse-complemented so that the telomere is at position 0; ed4_windows.sh,
  data/ed4_pairwise_windows.fa): SELF.p.chr9q, DAD.h2.chr9q.corrected, DAD.h2.chr9q.released, DAD.h2.chr3q, DAD.h1.chr3q, DAD.h1.chr9q.
  PAF: FastGA -T8 -f8192 -s200 -c70 -i0.6 -pafx FASTA FASTA (all against all, the settings of the earlier panel).
  TRACT_START/END: the called tract on SELF.p.chr9q in these telomere-first coordinates (shaded in its row and column).
Each cell: x = column sequence, y = row sequence, one line per alignment record, coloured by identity (PAF matches / block
length) on a 98-100% scale; records below 98% in grey. out: OUT_PREFIX.pdf and .png.
"""
import sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
fa, paf, ta, tb, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
L = {}
for l in open(fa):
    if l.startswith(">"): n = l[1:].split()[0]; L[n] = 0
    else: L[n] += len(l.strip())
order = ["SELF.p.chr9q", "DAD.h2.chr9q.corrected", "DAD.h2.chr9q.released", "DAD.h2.chr3q", "DAD.h1.chr3q", "DAD.h1.chr9q"]
lab = {"SELF.p.chr9q": "SELF.p chr9q", "DAD.h2.chr9q.corrected": "DAD.h2 chr9q\ncorrected", "DAD.h2.chr9q.released": "DAD.h2 chr9q\nreleased",
       "DAD.h2.chr3q": "DAD.h2 chr3q", "DAD.h1.chr3q": "DAD.h1 chr3q", "DAD.h1.chr9q": "DAD.h1 chr9q"}
assert sorted(order) == sorted(L), (sorted(L), order)
rec = {}
for l in open(paf):
    f = l.split("\t"); q, t = f[0], f[5]
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8]); idn = int(f[9]) / int(f[10])
    seg = ((ts, qs), (te, qe)) if f[4] == "+" else ((ts, qe), (te, qs))
    rec.setdefault((q, t), []).append((seg, idn))
    rec.setdefault((t, q, "swap"), []).append((((seg[0][1], seg[0][0]), (seg[1][1], seg[1][0])), idn))
n = len(order); cm = plt.get_cmap("turbo"); norm = Normalize(0.98, 1.0)
fig, axs = plt.subplots(n, n, figsize=(7.2, 7.6), gridspec_kw=dict(wspace=0.06, hspace=0.06, bottom=0.12, top=0.9, left=0.16, right=0.98))
for i, r in enumerate(order):
    for j, c in enumerate(order):
        ax = axs[i][j]; ax.set_xlim(0, L[c]); ax.set_ylim(0, L[r]); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_color("#c9c8c4"); sp.set_linewidth(0.6)
        if r == "SELF.p.chr9q": ax.axhspan(ta, tb, color="#eb6834", alpha=0.12, lw=0)
        if c == "SELF.p.chr9q": ax.axvspan(ta, tb, color="#eb6834", alpha=0.12, lw=0)
        rs = rec.get((r, c)) or rec.get((r, c, "swap")) or []
        if rs:
            cols = [cm(norm(v)) if v >= 0.98 else (0.6, 0.6, 0.6, 1) for _, v in rs]
            ax.add_collection(LineCollection([s for s, _ in rs], colors=cols, linewidths=1.1))
        if j == 0: ax.set_ylabel(lab[r], rotation=0, ha="right", va="center", fontsize=8.5, labelpad=4)
        if i == 0: ax.set_title(lab[c], fontsize=8.5, rotation=0, va="bottom")
fig.text(0.57, 0.955, "telomere-first 60 kb ends (x: column, y: row); shaded: former 9q/3q tract on SELF.p chr9q",
         ha="center", fontsize=8, color="#52514e")
cax = fig.add_axes([0.35, 0.06, 0.45, 0.018])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal")
cb.set_label("alignment identity (grey: below 98%)", fontsize=8); cb.ax.tick_params(labelsize=7)
fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=250)
print(f"wrote {out}.pdf/.png; records {sum(len(v) for k, v in rec.items() if len(k) == 2)}", file=sys.stderr)
