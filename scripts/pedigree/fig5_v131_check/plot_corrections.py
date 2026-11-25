#!/usr/bin/env python3
"""16:9 images that explain the corrections of the WashU pedigree assemblies (v1.1 and v1.3.1).

usage: plot_corrections.py RESULTS_DIR OUT_DIR
RESULTS_DIR = results/corrected_release: v1.1/ and v1.3.1/ (g1_changes_report.py tables), original_contig_lengths.tsv,
figdata/ (fixes_figdata.py; the end-fix image is skipped when it is missing).
out: OUT_DIR/1_method.png, 2_quality.png, 3_gaps.png, 4_where_v1.1.png, 4_where_v1.3.1.png, 5_end_fixes.png (13.33 x 7.5 in, 200 dpi)
"""
import collections, csv, os, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
R, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#a3a29e"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF})
W, H, DPI = 13.333, 7.5, 200
PEOPLE = [("PAN010", "MOM"), ("PAN011", "DAD"), ("PAN027", "SELF"), ("PAN028", "KID")]
VERS = [("v11", "v1.1"), ("v131", "v1.3.1")]
def hapname(sm, h):
    n = dict(PEOPLE)[sm]
    return f"{n} {'maternal' if h == 1 else 'paternal'}" if sm == "PAN027" else f"{n} h{h}"
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))
summ = {v: [r for sm, _ in PEOPLE for r in tsv(f"{R}/{d}/{v}.{sm}.summary.tsv")] for v, d in VERS}
chg = {v: [r for sm, _ in PEOPLE for r in tsv(f"{R}/{d}/{v}.{sm}.changes.tsv")] for v, d in VERS}
lens = collections.defaultdict(dict)
for r in tsv(f"{R}/original_contig_lengths.tsv"): lens[(r["version"], r["sample"], int(r["haplotype"]))][r["contig"].split("#")[-1]] = int(r["length"])
def save(fig, name):
    fig.savefig(f"{O}/{name}", dpi=DPI, facecolor=SURF); plt.close(fig); print(f"wrote {O}/{name}", file=sys.stderr)
def title(fig, t, s):
    fig.text(0.04, 0.93, t, fontsize=24, weight="bold", color=INK); fig.text(0.04, 0.885, s, fontsize=14, color=INK2)
# 1 method
fig = plt.figure(figsize=(W, H)); title(fig, "How we corrected the four assemblies", "Same steps and the same script for v1.1 and v1.3.1, each person checked with its own reads")
n_reg = {v: sum(1 for r in chg[v]) for v, _ in VERS}
tasks = {"v11": 481 + 385 + 218 + 484, "v131": 486 + 385 + 214 + 486}
g1 = {v: sum(r["step"] == "G1" for r in chg[v]) for v, _ in VERS}
steps = [("1  Find suspect regions", "HiFi and ONT depth and mismatches\nin 10 kb bins, plus merqury\nerror k-mers (Illumina)", f"{tasks['v11']:,} regions (v1.1)\n{tasks['v131']:,} regions (v1.3.1)"),
         ("2  Reassemble locally", "own HiFi and ONT reads within\n400 kb of each region, hifiasm\nwith HiFi, ONT and both", ""),
         ("3  Let the reads vote", "replace only when both HiFi and\nONT reads fit the new piece better;\ngaps >= 100 kb (rDNA) left alone", ""),
         ("4  Polish and check", "HiFi polishing; drop any piece\nthat adds merqury error k-mers;\napply the kept pieces again", f"{g1['v11']} pieces kept (v1.1)\n{g1['v131']} pieces kept (v1.3.1)"),
         ("5  Fix two chromosome ends", "DAD h2 chr9q: end from a DAD-only\nhifiasm contig; SELF maternal\nchr22q: 10.4 kb tail trimmed", "the same 2 fixes\nin both versions"),
         ("6  Corrected assemblies", "one script rebuilds all 8 from the\nrelease, identical base by base\nto the files we share", f"{n_reg['v11']} changes (v1.1)\n{n_reg['v131']} changes (v1.3.1)")]
ax = fig.add_axes([0.03, 0.03, 0.94, 0.80]); ax.set_xlim(0, 3); ax.set_ylim(0, 2); ax.axis("off")
for i, (t, d, n) in enumerate(steps):
    r, cidx = divmod(i, 3); x, y0 = cidx + 0.08, 1.06 - r * 1.02; c = ORANGE if i == 4 else BLUE
    ax.add_patch(FancyBboxPatch((x, y0), 0.84, 0.86, boxstyle="round,pad=0.0,rounding_size=0.03", fc="white", ec=c, lw=2))
    ax.add_patch(Rectangle((x, y0 + 0.72), 0.84, 0.14, fc=c, ec=c))
    ax.text(x + 0.42, y0 + 0.79, t, ha="center", va="center", fontsize=15, color="white", weight="bold")
    ax.text(x + 0.42, y0 + 0.50, d, ha="center", va="center", fontsize=13, color=INK, linespacing=1.4)
    if n: ax.text(x + 0.42, y0 + 0.15, n, ha="center", va="center", fontsize=13, color=INK, weight="bold", linespacing=1.4)
    if cidx < 2: ax.annotate("", xy=(x + 0.91, y0 + 0.43), xytext=(x + 0.85, y0 + 0.43), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.8))
ax.annotate("", xy=(0.5, 0.915), xytext=(2.5, 1.05), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.8))
save(fig, "1_method.png")
# 2 quality
import matplotlib.transforms as mtr
fig = plt.figure(figsize=(W, H)); title(fig, "Base accuracy goes up in all 16 haplotypes", "merqury QV from Illumina 21-mers (k-mers in the assembly but not in the reads); +1 QV = 21% fewer errors")
for j, (v, lab) in enumerate(VERS):
    ax = fig.add_axes([0.13 + j * 0.44, 0.14, 0.235, 0.66]); rows = summ[v]; tr = mtr.blended_transform_factory(ax.transAxes, ax.transData)
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i; a_, b_ = float(r["QV_orig"]), float(r["QV_final"])
        ax.plot([a_, b_], [y, y], color=GRID, lw=3, solid_capstyle="round", zorder=1)
        ax.scatter([a_], [y], s=70, fc="white", ec=GRAY, lw=2, zorder=2); ax.scatter([b_], [y], s=80, color=BLUE, ec="white", lw=1.5, zorder=3)
        ko, kf = int(r["merqury_kmers_orig"]), int(r["merqury_kmers_final"])
        ax.text(1.03, y, f"{a_:.2f} → {b_:.2f}   {100 * (kf - ko) / ko:+.1f}%", transform=tr, va="center", ha="left", fontsize=11, color=INK2)
    ax.text(1.03, len(rows) - 0.4, "QV                     error k-mers", transform=tr, ha="left", fontsize=10.5, color=INK2)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([hapname(r["sample"], int(r["haplotype"])) for r in rows][::-1] if j == 0 else [])
    ax.tick_params(axis="y", length=0); ax.set_xlim(56.5, 61.2); ax.set_xlabel("QV"); ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_title(lab, loc="left", fontsize=16, weight="bold", color=INK)
fig.legend(handles=[plt.scatter([], [], s=70, fc="white", ec=GRAY, lw=2, label="original release"), plt.scatter([], [], s=80, color=BLUE, label="corrected")],
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=12)
save(fig, "2_quality.png")
# 3 gaps
fig = plt.figure(figsize=(W, H)); title(fig, "Gaps: small ones closed when the reads span them", "Bases called N and number of gaps (N runs) per haplotype, original release -> corrected; gaps >= 100 kb (rDNA) left as they are")
ax = fig.add_axes([0.04, 0.17, 0.92, 0.66]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 9.4)
cols = [(0.0, "", "left"), (0.30, "v1.1  N bases (%)", "center"), (0.47, "v1.1  gaps", "center"), (0.70, "v1.3.1  N bases (%)", "center"), (0.87, "v1.3.1  gaps", "center")]
for x, t, ha in cols: ax.text(x, 9.0, t, ha=ha, va="center", fontsize=13, weight="bold", color=INK)
ax.plot([0, 1], [8.6, 8.6], color=GRID, lw=1)
for i in range(8):
    y = 7.9 - i; ax.text(0.0, y, hapname(summ["v11"][i]["sample"], int(summ["v11"][i]["haplotype"])), ha="left", va="center", fontsize=13, color=INK)
    for k, v in enumerate(("v11", "v131")):
        r = summ[v][i]; ch = r["orig_N"] != r["final_N"]
        for x, t in ((0.30 + 0.40 * k, f"{r['orig_N_pct']} -> {r['final_N_pct']}"), (0.47 + 0.40 * k, f"{r['orig_N_runs']} -> {r['final_N_runs']}")):
            ax.text(x, y, t, ha="center", va="center", fontsize=13, color=BLUE if ch else INK2, weight="bold" if ch else "normal")
    if i < 7: ax.plot([0, 1], [y - 0.5, y - 0.5], color="#efeeea", lw=0.8)
filled = [r for r in chg["v11"] if r["step"] == "G1" and int(r["orig_N"]) > int(r["final_N"])]
fig.text(0.04, 0.10, "Closed in v1.1 (blue): " + "; ".join(f"{hapname(r['sample'], int(r['haplotype']))} {r['contig'].split('#')[-1]}: {int(r['orig_N']):,} N in {int(r['orig_bp']):,} bp -> {int(r['final_bp']):,} bp with {r['final_N']} N" for r in filled),
         fontsize=11, color=INK2, wrap=True)
fig.text(0.04, 0.045, "The new sequence matches the independent v1.3.1 release at the same place: DAD chr9 70,702 of 70,711 bases, KID chr19 52,576 of 52,580, KID chr14 98.8%.", fontsize=11, color=INK2)
save(fig, "3_gaps.png")
# 4 where
CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
for v, lab in VERS:
    fig = plt.figure(figsize=(W, H)); n = collections.Counter(r["step"] for r in chg[v])
    title(fig, f"Where we changed the sequence ({lab})", f"{n['G1']} G1 pieces (each tick = one piece, drawn at least 1.5 Mb wide) and 2 chromosome-end fixes, on the original coordinates")
    for k, (sm, nm) in enumerate(PEOPLE):
        ax = fig.add_axes([0.04 + k * 0.24, 0.07, 0.21, 0.70]); ax.set_title(nm, fontsize=16, weight="bold", loc="left")
        chrs = [c for c in CHR if c in lens[(v, sm, 1)] or c in lens[(v, sm, 2)]]
        for i, c in enumerate(chrs):
            for h in (1, 2):
                L = lens[(v, sm, h)].get(c)
                if L is None: continue
                y = len(chrs) - 1 - i + (0.17 if h == 1 else -0.17)
                ax.add_patch(Rectangle((0, y - 0.12), L / 1e6, 0.24, fc="#eceae6", ec="none"))
        for r in chg[v]:
            if r["sample"] != sm: continue
            c = r["contig"].split("#")[-1]; i = chrs.index(c); h = int(r["haplotype"]); y = len(chrs) - 1 - i + (0.17 if h == 1 else -0.17)
            a, b = int(r["orig_start"]) / 1e6, int(r["orig_end"]) / 1e6; col = ORANGE if r["step"] != "G1" else AQUA if int(r["orig_N"]) > int(r["final_N"]) else BLUE
            ax.add_patch(Rectangle((a, y - 0.17), max(b - a, 1.5), 0.34, fc=col, ec="none", zorder=3 + (col != BLUE)))
        ax.set_ylim(-0.7, len(chrs) - 0.3); ax.set_xlim(0, 260); ax.set_yticks(range(len(chrs))); ax.set_yticklabels([c[3:] for c in chrs][::-1], fontsize=9.5)
        ax.set_xlabel("position (Mb)", fontsize=11); ax.tick_params(axis="x", labelsize=10)
        for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
        ax.tick_params(axis="y", length=0)
    hs = [Rectangle((0, 0), 1, 1, color=BLUE, label="G1 piece"), Rectangle((0, 0), 1, 1, color=ORANGE, label="end fix (DAD chr9q, SELF chr22q)"),
          Rectangle((0, 0), 1, 1, color="#eceae6", label="chromosome (upper bar h1 or maternal, lower h2 or paternal)")]
    if v == "v11": hs.insert(1, Rectangle((0, 0), 1, 1, color=AQUA, label="G1 piece that closed a gap"))
    fig.legend(handles=hs, loc="upper left", ncol=len(hs), frameon=False, fontsize=11, bbox_to_anchor=(0.035, 0.87))
    save(fig, f"4_where_{lab}.png")
# 5 end fixes
F = f"{R}/figdata"
if os.path.exists(f"{F}/chr9q_windows.tsv"):
    fx = {(r["version"], r["sample"], r["step"]): r for r in tsv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "g1_manual_fixes.tsv"))}
    cut = (lens[("v131", "PAN027", 1)]["chr22"] - int(fx[("v131", "PAN027", "tip_trim")]["pos"])) / 1000
    fig = plt.figure(figsize=(W, H)); title(fig, "The two chromosome-end fixes (same in v1.1 and v1.3.1)", "Both found at chromosome ends by comparing parent and child; v1.3.1 shown")
    a9 = fig.add_axes([0.08, 0.12, 0.37, 0.56]); w = tsv(f"{F}/chr9q_windows.tsv")
    x = [-(int(r["window_end_from_tip_bp"]) - 5000) / 1000 for r in w]
    for col, key, lab_ in ((GRAY, "original", "DAD h2 as released"), (BLUE, "corrected", "DAD h2 corrected")):
        pts = [(xx, 100 * float(r[f"identity_{key}"])) for xx, r in zip(x, w) if r[f"identity_{key}"] != "NA"]
        a9.plot([p_[0] for p_ in pts], [p_[1] for p_ in pts], color=col, lw=2, marker="o", ms=6, label=lab_, mec="white", mew=1)
    a9.set_xlabel("position to the chr9q tip (kb), 10 kb windows"); a9.set_ylabel("identity to DAD h2 (%)")
    a9.set_ylim(97.8, 100.15); a9.grid(color=GRID, lw=0.8); a9.set_axisbelow(True); a9.legend(frameon=False, fontsize=11.5, loc="lower left")
    for sp in ("top", "right"): a9.spines[sp].set_visible(False)
    lo9 = min(100 * float(r["identity_original"]) for r in w if r["identity_original"] != "NA")
    fig.text(0.04, 0.81, "DAD h2 chr9q end", fontsize=15, weight="bold")
    fig.text(0.04, 0.715, f"SELF's paternal chr9 (10 kb windows) vs DAD h2: the released end\ndid not match (down to {lo9:.1f}% identity in the last 70 kb). We replaced\nthe last 150,287 bp with 147,586 bp of a hifiasm contig from DAD reads only", fontsize=10.5, color=INK2, linespacing=1.35)
    b = [r for r in tsv(f"{F}/chr22q_bins.tsv") if int(r["bin_start_bp_from_end"]) < 22000]; xb = [-(int(r["bin_start_bp_from_end"]) + 250) / 1000 for r in b]
    ax3 = [fig.add_axes([0.57, 0.53 - k * 0.195, 0.39, 0.165]) for k in range(3)]
    ax3[0].bar(xb, [float(r["TTAGGG_pct"]) for r in b], width=0.42, color=BLUE); ax3[0].set_ylabel("TTAGGG %", fontsize=11)
    ax3[1].bar(xb, [int(r["ont_read_ends"]) for r in b], width=0.42, color=BLUE); ax3[1].set_ylabel("ONT reads\nending here", fontsize=11)
    ax3[2].bar([xx - 0.11 for xx in xb], [int(r["merqury_kmers_original"]) for r in b], width=0.22, color=GRAY, label="released")
    ax3[2].bar([xx + 0.11 for xx in xb], [int(r["merqury_kmers_corrected"]) for r in b], width=0.22, color=BLUE, label="corrected")
    ax3[2].set_ylabel("merqury\nerror k-mers", fontsize=11); ax3[2].legend(frameon=False, fontsize=10.5, loc="upper left")
    for k, a_ in enumerate(ax3):
        a_.axvline(-cut, color=ORANGE, lw=2); a_.set_xlim(-22, 0.3); a_.grid(axis="y", color=GRID, lw=0.8); a_.set_axisbelow(True); a_.tick_params(labelsize=10)
        for sp in ("top", "right"): a_.spines[sp].set_visible(False)
        if k < 2: a_.tick_params(labelbottom=False)
    ax3[2].set_xlabel("position to the released chr22q end (kb), 500 bp bins")
    ax3[0].set_ylim(0, 125); ax3[0].text(-cut + 0.2, 122, f"cut: last {cut * 1000:,.0f} bp removed", color=ORANGE, fontsize=11, va="top")
    fig.text(0.54, 0.81, "SELF maternal chr22q end", fontsize=15, weight="bold")
    ko = sum(int(r["merqury_kmers_original"]) for r in b if r["removed_by_trim"] == "1")
    fig.text(0.54, 0.715, f"The canonical telomere ends {cut:.1f} kb before the released end, and 43 of the 76\nONT reads that cover the last 30 kb end within 1.5 kb of it. The tail after it is\nnot telomere and holds {ko:,} merqury error k-mers; we trimmed it", fontsize=10.5, color=INK2, linespacing=1.35)
    save(fig, "5_end_fixes.png")
