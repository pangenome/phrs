#!/usr/bin/env python3
"""16:9 schematic of the chr10q D4Z4 array length in DAD h2 and in the copy SELF got from him (v1.3.1).

usage: plot_d4z4_chr10q.py SPAN_READS_TSV OUT_PNG
SPAN_READS_TSV (check/d4z4_chr10/span_reads.tsv): one row per primary ONT read that spans +- 20 kb around the array
(DAD h2 at PAN011#2#chr10:135,057,730; SELF paternal at PAN027#2#chr10:134,417,447), net inserted minus deleted bases
against its own assembly. SELF's assembled array is 39,698 bp longer than DAD h2's, so SELF reads are shifted by 39,698 bp
to share DAD h2's assembly as the reference. Unit counts: D4Z4 unit 3.3 kb; DAD h2 assembled array 8 units (25 kb),
DAD h1 41 units (135 kb), from the copy-count scan of 2026-11-22.
"""
import csv, random, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
T, OUT = sys.argv[1:3]
INK, INK2, GRID, SURF, BLUE, GRAY, ORANGE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb", "#2a78d6", "#a3a29e", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF})
rows = list(csv.DictReader(open(T), delimiter="\t"))
dad = [int(r["net_indel_bp"]) for r in rows if r["sample"] == "DAD_h2"]
slf = [int(r["net_indel_bp"]) + 39698 for r in rows if r["sample"] == "SELF_paternal"]
nl = sum(d > 20000 for d in dad); long_kb = sorted(d for d in dad if d > 20000)
fig = plt.figure(figsize=(13.333, 7.5))
fig.text(0.04, 0.93, "chr10q D4Z4: the daughter's longer array is already in her father", fontsize=24, weight="bold")
fig.text(0.04, 0.885, "D4Z4 = 3.3 kb repeat unit next to the chr10q telomere (same family as the 4q FSHD array); ONT reads, v1.3.1", fontsize=14, color=INK2)
ax = fig.add_axes([0.03, 0.19, 0.52, 0.64]); ax.set_xlim(0, 52); ax.set_ylim(0, 10); ax.axis("off")
def chrom(y, units, label, sub, col):
    ax.add_patch(Rectangle((0.5, y - 0.18), 8, 0.36, fc="#d8d6d0", ec="none")); ax.text(4.5, y - 0.55, "unique sequence", ha="center", fontsize=9.5, color=INK2)
    for i in range(units): ax.add_patch(Rectangle((8.8 + i * 0.95, y - 0.32), 0.8, 0.64, fc=col, ec="white", lw=0.8))
    x = 8.8 + units * 0.95; ax.add_patch(Rectangle((x + 0.1, y - 0.18), 1.6, 0.36, fc=INK, ec="none")); ax.text(x + 0.9, y - 0.55, "telomere", ha="center", fontsize=9.5, color=INK2)
    ax.text(0.5, y + 0.55, label, fontsize=13, weight="bold"); ax.text(x + 3.6, y, sub, va="center", fontsize=11, color=INK2)
chrom(8.8, 8, "DAD h2 (the copy he passed on): assembly", f"8 units, 25 kb\n{len(dad) - nl} of {len(dad)} reads", GRAY)
chrom(6.7, 20, "DAD h2: the other half of his reads", f"about 20 units\n{nl} of {len(dad)} reads,\n+{long_kb[0] / 1000:.1f} to +{long_kb[-1] / 1000:.1f} kb", BLUE)
chrom(4.4, 20, "SELF paternal chr10q (from DAD h2)", f"about 20 units\n{len(slf)} of {len(slf)} reads", BLUE)
ax.add_patch(FancyArrowPatch((14, 6.2), (14, 5.15), arrowstyle="-|>", mutation_scale=18, color=INK2, lw=1.6)); ax.text(14.6, 5.7, "transmitted", fontsize=11, color=INK2)
ax.text(0.5, 2.3, "DAD h1 (not passed on) has 41 units (135 kb).", fontsize=11.5, color=INK2)
ax.text(0.5, 1.3, "Both lengths come from the same sequencing runs and all 16 reads match h2 better than h1:", fontsize=11.5, color=INK2)
ax.text(0.5, 0.6, "DAD carries two lengths of this array on h2 (mosaic); the assembly kept the short one.", fontsize=11.5, color=INK2)
b = fig.add_axes([0.63, 0.25, 0.33, 0.55]); random.seed(1)
for k, (v, lab_) in enumerate(((dad, f"DAD h2 reads (n={len(dad)})"), (slf, f"SELF paternal reads (n={len(slf)})"))):
    y = 1 - k
    b.scatter([x / 1000 for x in v], [y + random.uniform(-0.18, 0.18) for _ in v], s=55, color=[BLUE if x > 20000 else GRAY for x in v], ec="white", lw=0.8, zorder=3)
b.set_yticks([1, 0]); b.set_yticklabels([f"DAD h2\n(n={len(dad)})", f"SELF paternal\n(n={len(slf)})"]); b.set_ylim(-0.6, 1.6)
b.set_xlabel("array length vs DAD h2 assembly (kb)\nper ONT read spanning the array"); b.set_xlim(-8, 46)
b.axvline(0, color=GRAY, lw=1.2, ls=":"); b.axvline(39.698, color=BLUE, lw=1.2, ls=":")
b.text(0, 1.52, "DAD assembly", ha="center", fontsize=10.5, color=INK2); b.text(39.698, 1.52, "SELF assembly", ha="center", fontsize=10.5, color=INK2)
b.grid(axis="x", color=GRID, lw=0.8); b.set_axisbelow(True)
for sp in ("top", "right", "left"): b.spines[sp].set_visible(False)
b.tick_params(axis="y", length=0)
fig.text(0.04, 0.035, "So SELF did not gain the 12 extra units at transmission: she inherited a length her father already carries. In the parent-child scan\n"
         "of transmitted ends this is the only difference (SELF paternal chr10q 134.41-134.46 Mb) besides the acrocentric short arms, PAR1,\n"
         "the fixed chr9q error and 2-window runs at telomeres.", fontsize=12, color=INK, linespacing=1.35)
fig.savefig(OUT, dpi=200, facecolor=SURF); print(f"wrote {OUT}: DAD h2 {len(dad)} reads ({nl} long), SELF {len(slf)} reads", file=sys.stderr)
