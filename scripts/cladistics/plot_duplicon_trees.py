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
ap.add_argument("--json", nargs="+", default=["scripts/cladistics/duplicon_monophyly.json"])
ap.add_argument("--regime", default="regime_f8192")
ap.add_argument("--stratum", default="all_blocks")
ap.add_argument("--labels", default="arm")
ap.add_argument("--out", required=True)
a = ap.parse_args()
TOK = re.compile(r"\(|\)|,|:[-0-9.eE]+|;|[^(),:;]+")
COLS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
        "#17becf"]  # tab10: distinct colours for up to 10 arms


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
for path in a.json:
    d = json.load(open(path))
    for fam, v in d["families"].items():
        rec = v["regimes"][a.regime]["strata"][a.stratum][a.labels]
        tree = parse(rec["point_newick"])
        if ntip(tree) == 2 * rec["n_taxa"] and ntip(tree[0]) == rec["n_taxa"]:
            tree = tree[0]  # Newick written before the to_newick fix: the tree twice
        lab = rec["labels"]
        mono = sum(x["mono_point"] == 1 for x in lab.values())
        fams.append((mono / len(lab), fam, tree, rec["n_taxa"], mono, len(lab)))
        sys.stderr.write(f"{path}: {fam} {rec['n_taxa']} copies, {mono}/{len(lab)} arms in one clade\n")
fams.sort(key=lambda f: (-f[0], f[1]))
plt.rcParams.update({"font.size": 7})
fig = plt.figure(figsize=(7.2, 0.78 * len(fams) + 0.4))
gs = fig.add_gridspec(2 * len(fams), 1, height_ratios=[1, 0.16] * len(fams), hspace=0.05,
                      left=0.16, right=0.76, top=0.99, bottom=0.01)
for k, (_, fam, tree, n, mono, narm) in enumerate(fams):
    tips, segs = [], []
    draw(tree, tips, segs)
    if len(tips) != n or len(set(tips)) != n:
        raise SystemExit(f"FATAL: {fam}: {len(tips)} tips ({len(set(tips))} distinct) vs n_taxa {n}")
    arm = [t.split("|")[1].split("_", 1)[1] for t in tips]
    arms = sorted(set(arm), key=lambda x: (int(re.sub(r"\D", "", x.split("_")[0]) or 99), x))
    if len(arms) > len(COLS):
        raise SystemExit(f"FATAL: {fam} has {len(arms)} arms, {len(COLS)} colours")
    cm = {x: COLS[i] for i, x in enumerate(arms)}
    ax, axs = fig.add_subplot(gs[2 * k]), fig.add_subplot(gs[2 * k + 1])
    for xs, ys in segs:
        ax.plot(xs, ys, color="#444", lw=0.3)
    ax.set_xlim(-0.5, n - 0.5); axs.set_xlim(-0.5, n - 0.5)
    axs.bar(range(n), [1] * n, width=1.0, color=[cm[x] for x in arm], linewidth=0)
    for z in (ax, axs):
        z.set_xticks([]); z.set_yticks([])
        for sp in z.spines.values():
            sp.set_visible(False)
    ax.text(-0.02, 0.35, f"{fam.replace('_', ' ')}\n{n} copies\n{mono}/{narm} arms one clade",
            transform=ax.transAxes, ha="right", va="center")
    short = lambda x: x.replace("chr", "").replace("_parm", "p").replace("_qarm", "q")
    for i, x in enumerate(arms):
        axs.text(1.01 + 0.06 * (i % 4), 1.4 - 1.7 * (i // 4), short(x), color=cm[x], transform=axs.transAxes,
                 ha="left", va="center", fontweight="bold")
fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=200)
sys.stderr.write(f"wrote {a.out}.pdf/.png\n")
