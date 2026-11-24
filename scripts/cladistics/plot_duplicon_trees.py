#!/usr/bin/env python3
"""Draw the duplicon copy trees of duplicon_monophyly.py (R1-M1 test 1), one panel per family.

Input: the point neighbour-joining topologies stored in duplicon_monophyly.json (regime, stratum and label
set as arguments; no branch lengths are stored, so we draw cladograms with aligned tips). Tips are coloured
by chromosome arm; the root is an arbitrary edge, not a biological root. The JSON Newick written before
the to_newick fix is "(T_u,T_v);" where T_u is the whole tree hung from u (every tip once) and T_v repeats
it from v: for that format we keep T_u. Either way the tip count must equal n_taxa.
A family whose copies cluster by end shows one colour per clade; interleaved colours are copies from
different ends in the same clade. Output: <out>.pdf and <out>.png.
"""
import argparse, json, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--json", default="scripts/cladistics/duplicon_monophyly.json")
ap.add_argument("--regime", default="regime_f8192")
ap.add_argument("--stratum", default="all_blocks")
ap.add_argument("--labels", default="arm")
ap.add_argument("--out", required=True)
a = ap.parse_args()
TOK = re.compile(r"\(|\)|,|:[-0-9.eE]+|;|[^(),:;]+")


def parse(s):
    """Newick -> nested lists; leaf = (name, length), node = (children, length)."""
    toks, i = TOK.findall(s.strip()), 0

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
            name = None
            if i < len(toks) and toks[i] not in "(),;" and not toks[i].startswith(":"):
                i += 1
            ln = float(toks[i][1:]) if i < len(toks) and toks[i].startswith(":") else 0.0
            i += 1 if i < len(toks) and toks[i].startswith(":") else 0
            return (kids, ln)
        name = toks[i]
        i += 1
        ln = float(toks[i][1:]) if toks[i].startswith(":") else 0.0
        i += 1 if toks[i].startswith(":") else 0
        return (name, ln)
    return node()


def height(t):
    return 0 if isinstance(t[0], str) else 1 + max(height(k) for k in t[0])


def layout(t, ys, out):
    """cladogram: x = -height, tips at 0; returns (x, y) of the node"""
    lab, _ = t
    if isinstance(lab, str):
        y = len(ys)
        ys.append(lab)
        out["tips"].append((0, y, lab))
        return 0, y
    x = -height(t)
    kids = [layout(k, ys, out) for k in lab]
    for xk, yk in kids:
        out["h"].append((x, xk, yk))
    ky = [yk for _, yk in kids]
    out["v"].append((x, min(ky), max(ky)))
    return x, (min(ky) + max(ky)) / 2


d = json.load(open(a.json))
fams = list(d["families"])
cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
        "#17becf"]  # tab10: distinct colours for up to 10 arms
fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 9))
for ax, fam in zip(axes, fams):
    rec = d["families"][fam]["regimes"][a.regime]["strata"][a.stratum][a.labels]
    tree = parse(rec["point_newick"])
    ntip = lambda t: 1 if isinstance(t[0], str) else sum(ntip(k) for k in t[0])
    if ntip(tree) == 2 * rec["n_taxa"] and ntip(tree[0][0]) == rec["n_taxa"]:
        tree = tree[0][0]  # Newick written before the to_newick fix: the tree twice
    out = {"tips": [], "h": [], "v": []}
    layout(tree, [], out)
    names = [t[2] for t in out["tips"]]
    if len(names) != rec["n_taxa"] or len(set(names)) != len(names):
        raise SystemExit(f"FATAL: {fam}: {len(names)} tips ({len(set(names))} distinct) vs n_taxa {rec['n_taxa']}")
    arms = sorted({t[2].split("|")[1].split("_", 1)[1] for t in out["tips"]})
    if len(arms) > len(cols):
        raise SystemExit(f"FATAL: {fam} has {len(arms)} arms, {len(cols)} colours")
    cm = {x: cols[i] for i, x in enumerate(arms)}
    for x0, x1, y in out["h"]:
        ax.plot([x0, x1], [y, y], color="#555", lw=0.6)
    for x, y0, y1 in out["v"]:
        ax.plot([x, x], [y0, y1], color="#555", lw=0.6)
    for x, y, lab in out["tips"]:
        ax.scatter([x], [y], s=14, color=cm[lab.split("|")[1].split("_", 1)[1]], zorder=3)
    for x in arms:
        ax.scatter([], [], s=30, color=cm[x], label=x.replace("_parm", "p").replace("_qarm", "q").replace("chr", ""))
    ax.legend(loc="lower left", fontsize=10, frameon=False, ncol=1 if len(arms) < 4 else 2)
    ax.set_title(f"{fam.replace('_', ' ')}\n{len(out['tips'])} copies", fontsize=12, loc="left")
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(False)
    sys.stderr.write(f"{fam}: {len(out['tips'])} tips, arms {arms}\n")
fig.suptitle(f"Duplicon copy trees (neighbour joining; {a.regime.replace('regime_', 'FastGA -')}, "
             f"{a.stratum}); colour = chromosome arm; topology only", fontsize=13)
fig.tight_layout()
fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=150)
sys.stderr.write(f"wrote {a.out}.pdf/.png\n")
