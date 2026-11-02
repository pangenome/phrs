#!/usr/bin/env python3
"""Plot paralog/allelic divergence ratio against distance from the telomere-subtelomere junction."""
import argparse, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("--json", required=True)
ap.add_argument("--out", default="response/notes/telomere_profile"); a = ap.parse_args()
d = json.load(open(a.json))["rows"]
mid = [(r["from_kb"] + r["to_kb"]) / 2 for r in d]
al = [r["allelic_median"] for r in d]; pa = [r["paralog_median"] for r in d]
ra = [r.get("ratio") for r in d]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
a0 = ax[0]
a0.plot(mid, [x * 100 for x in al], "o-", ms=3, color="#4c72b0", label="allelic (same arm)")
a0.plot(mid, [x * 100 for x in pa], "o-", ms=3, color="#757575", label="paralogous (other chromosome)")
a0.axvspan(0, 10, color="#fde0c5", zorder=0, label="first 10 kb")
a0.set_xscale("log"); a0.set_xlabel("distance from telomere (kb)"); a0.set_ylabel("median substitution divergence (%)")
a0.set_title("A  both classes diverge more at the tip", fontsize=9, loc="left"); a0.legend(fontsize=6.5)
a1 = ax[1]
a1.plot(mid, ra, "o-", ms=3, color="#D95F02")
a1.axhline(1, color="k", lw=1, ls=":")
a1.axvspan(0, 10, color="#fde0c5", zorder=0)
a1.set_xscale("log"); a1.set_xlabel("distance from telomere (kb)")
a1.set_ylabel("paralog / allelic divergence ratio")
a1.set_title("B  ratio approaches 1 at the junction", fontsize=9, loc="left")
a1.text(0.55, 1.15, "allelic-like", fontsize=6.5, color="#666666")
for x, y in ((1, ra[0]), (500, ra[-1])):
    a1.annotate(f"{y:.2f}x", (x, y), textcoords="offset points", xytext=(4, 5), fontsize=6.5)
plt.tight_layout()
for ext in ("pdf", "png"): fig.savefig(a.out + "." + ext, dpi=190, bbox_inches="tight")
print("wrote", a.out)
print("first 10 kb ratios:", [f"{r:.2f}" for r in ra[:5]])
print("interior 50-200 kb ratios:", [f"{r:.2f}" for r in ra[25:29] if r])
