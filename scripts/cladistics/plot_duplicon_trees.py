#!/usr/bin/env python3
"""Draw all duplicon copy trees of duplicon_monophyly.py (R1-M1 test 1) in one compact figure.

Input: one or more duplicon_monophyly.json; we use the point neighbour-joining topology of each family
(regime, stratum and label set as arguments). No branch lengths are stored, so each tree is a cladogram,
drawn sideways with its copies along the x axis; the strip under it gives each copy's chromosome arm.
The root is an arbitrary edge, not a biological root. The JSON Newick written before the to_newick fix is
"(T_u,T_v);" where T_u is the whole tree hung from u (every tip once) and T_v repeats it from v: for that
format we keep T_u. Either way the tip count must equal n_taxa. Families are sorted by the share of their
arms whose copies form one clade (point tree), most separated first. One block of colour per arm in the
strip = copies cluster by end; interleaved colours = copies of different ends in the same clade.
Output: <out>.pdf and <out>.png.
"""
import argparse, json, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--json", nargs="+", default=["scripts/cladistics/duplicon_monophyly.json"],
                help="PATH or PATH:FAM1,FAM2 to draw only those families of that file")
ap.add_argument("--cols", type=int, default=2)
ap.add_argument("--min-arms", type=int, default=2, help="skip families whose copies were found on fewer arms")
ap.add_argument("--controls", default="ACRO_P_REPEAT,AMD1,GTPB,PLCXD1,PPP2R3,SHOX,SPRY3,VAM,WASIR1",
                help="families marked as positive controls (acrocentric and pseudoautosomal)")
ap.add_argument("--regime", default="regime_f8192")
ap.add_argument("--stratum", default="all_blocks")
ap.add_argument("--labels", default="arm")
ap.add_argument("--out", required=True)
a = ap.parse_args()
TOK = re.compile(r"\(|\)|,|:[-0-9.eE]+|;|[^(),:;]+")
COLS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
        "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7",
        "#dbdb8d", "#9edae5"]  # tab20: distinct colours for up to 20 arms


def parse(s):
    """Newick (topology, optional lengths ignored) -> nested lists; leaf = name, node = [children]."""
    toks, i = [t for t in TOK.findall(s.strip()) if not t.startswith(":")], 0

    def node():
        nonlocal i
        if toks[i] == "(":
            i += 1
            kids = [node()]
            while toks[i] == ",":
                i += 1
                kids.append(node())
            if toks[i] != ")":
                raise SystemExit(f"FATAL: newick parse error at token {i}: {toks[i]!r}")
            i += 1
            if i < len(toks) and toks[i] not in "(),;":
                i += 1  # internal node label
            return kids
        i += 1
        return toks[i - 1]
    return node()


def ntip(t):
    return 1 if isinstance(t, str) else sum(ntip(k) for k in t)


def draw(t, tips, segs):
    """cladogram sideways: tips at y = 0 in order along x, node y = its height; returns (x, y)"""
    if isinstance(t, str):
        tips.append(t)
        return len(tips) - 1, 0
    kids = [draw(k, tips, segs) for k in t]
    y = 1 + max(yk for _, yk in kids)
    for xk, yk in kids:
        segs.append(((xk, xk), (yk, y)))
    segs.append(((min(x for x, _ in kids), max(x for x, _ in kids)), (y, y)))
    return sum(x for x, _ in kids) / len(kids), y


fams = []
ctrl = set(a.controls.split(","))
for arg in a.json:
    path, _, keep = arg.partition(".json:")
    path = path if path.endswith(".json") else path + ".json"
    d = json.load(open(path))
    for fam, v in d["families"].items():
        if keep and fam not in keep.split(","):
            continue
        if "error" in v:
            sys.stderr.write(f"{path}: {fam} skipped, run error: {v['error'][:80]}\n")
            continue
        rec = v["regimes"][a.regime]["strata"][a.stratum][a.labels]
        tree = parse(rec["point_newick"])
        if ntip(tree) == 2 * rec["n_taxa"] and ntip(tree[0]) == rec["n_taxa"]:
            tree = tree[0]  # Newick written before the to_newick fix: the tree twice
        lab = rec["labels"]
        if len(lab) < a.min_arms:
            sys.stderr.write(f"{path}: {fam} skipped, copies on {len(lab)} arm(s)\n")
            continue
        mono = sum(x["mono_point"] == 1 for x in lab.values())
        fams.append((mono / len(lab), fam, tree, rec["n_taxa"], mono, len(lab)))
        sys.stderr.write(f"{path}: {fam} {rec['n_taxa']} copies, {mono}/{len(lab)} arms in one clade\n")
fams.sort(key=lambda f: (f[1] in ctrl, -f[0], f[1]))
plt.rcParams.update({"font.size": 6})
nrow = -(-len(fams) // a.cols)
fig = plt.figure(figsize=(7.2, 0.72 * nrow + 0.1))
outer = fig.add_gridspec(1, a.cols, left=0.01, right=0.99, top=0.985, bottom=0.0, wspace=0.06)
cells = [outer[c].subgridspec(3 * nrow, 1, height_ratios=[1, 0.16, 0.62] * nrow, hspace=0) for c in range(a.cols)]
for k, (_, fam, tree, n, mono, narm) in enumerate(fams):
    g = cells[k // nrow]
    r = k % nrow
    tips, segs = [], []
    draw(tree, tips, segs)
    if len(tips) != n or len(set(tips)) != n:
        raise SystemExit(f"FATAL: {fam}: {len(tips)} tips ({len(set(tips))} distinct) vs n_taxa {n}")
    arm = [t.split("|")[1].split("_", 1)[1] for t in tips]
    arms = sorted(set(arm), key=lambda x: (int(re.sub(r"\D", "", x.split("_")[0]) or 99), x))
    if len(arms) > len(COLS):
        raise SystemExit(f"FATAL: {fam} has {len(arms)} arms, {len(COLS)} colours")
    cm = {x: COLS[i] for i, x in enumerate(arms)}
    ax, axs = fig.add_subplot(g[3 * r]), fig.add_subplot(g[3 * r + 1])
    for xs, ys in segs:
        ax.plot(xs, ys, color="#444", lw=0.25)
    ax.set_xlim(-0.5, n - 0.5); axs.set_xlim(-0.5, n - 0.5)
    axs.bar(range(n), [1] * n, width=1.0, color=[cm[x] for x in arm], linewidth=0)
    for z in (ax, axs):
        z.set_xticks([]); z.set_yticks([])
        for sp in z.spines.values():
            sp.set_visible(False)
    ax.set_title(f"{fam.replace('_', ' ')}{' (control)' if fam in ctrl else ''}: {n} copies, "
                 f"{mono} of {narm} arms in one clade", loc="left", fontsize=6.5, pad=1.5)
    short = lambda x: x.replace("chr", "").replace("_parm", "p").replace("_qarm", "q")
    for i, x in enumerate(arms):
        axs.text(0.066 * i, -0.6, short(x), color=cm[x], transform=axs.transAxes,
                 ha="left", va="top", fontweight="bold", fontsize=5.5)
fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=200)
sys.stderr.write(f"wrote {a.out}.pdf/.png\n")
