#!/usr/bin/env python3
"""Fig5 plotter with every inter-chromosomal run outside acro_acro annotated.
usage: fig5_annotate.py PLOTTER.py [--min-kb N] -- <plotter arguments>
Runs the unchanged plotter (same runs, same ribbons); only draw_callouts picks all non-acro_acro runs, not four fixed events."""
import importlib.util, sys

plotter, rest = sys.argv[1], sys.argv[2:]
min_bp = 0
if rest[:1] == ["--min-kb"]:
    min_bp, rest = int(float(rest[1]) * 1000), rest[2:]
rest = rest[1:] if rest[:1] == ["--"] else rest
spec = importlib.util.spec_from_file_location("fig5_plotter", plotter)
P = importlib.util.module_from_spec(spec)
sys.modules["fig5_plotter"] = P
spec.loader.exec_module(P)

def merge(rs, gap=300_000):
    out = []
    for r in sorted(rs, key=lambda r: (r.query_chrom, r.target_chrom, r.donor_haplotype, r.query_start)):
        g = out[-1] if out else None
        if g and (g["q"], g["t"], g["h"]) == (r.query_chrom, r.target_chrom, r.donor_haplotype) and r.query_start - g["e"] <= gap:
            g["e"], g["bp"], g["n"] = max(g["e"], r.query_end), g["bp"] + r.bp, g["n"] + 1
        else:
            out.append({"q": r.query_chrom, "t": r.target_chrom, "h": r.donor_haplotype, "s": r.query_start, "e": r.query_end,
                        "bp": r.bp, "n": 1, "r": r})
    return out

def callouts_all(svg, runs, query_layout):
    sel = merge([r for r in runs if r.category != "acro_acro" and r.bp >= min_bp])
    sel.sort(key=lambda g: P.x_for(query_layout, g["q"], (g["s"] + g["e"]) // 2))
    print(f"annotating {len(sel)} groups of {sum(g['n'] for g in sel)} runs, all except acro_acro, of {len(runs)} drawn runs", file=sys.stderr)
    lanes = []
    for g in sel:
        run = g["r"]
        x0, x1 = P.interval_x(query_layout, g["q"], g["s"], g["e"])
        c = (x0 + x1) / 2
        w = max(170, len(P.label_text(run)) * 15)
        for i, last in enumerate(lanes):
            if c - w / 2 > last + 24:
                lane, lanes[i] = i, c + w / 2
                break
        else:
            lane = len(lanes)
            lanes.append(c + w / 2)
        y = P.Y_QUERY - 140 - lane * 80
        col = P.COLORS.get(run.category, P.COLORS["other_nonacro"])
        svg.line(c, P.Y_QUERY - 9, c, y + 8, col, 1.4, 0.9)
        svg.text(c, y, P.label_text(run), 27, "700", col, "middle")
        svg.text(c, y + int(29 * P.FONT_SCALE), f"{g['bp'] / 1000:.0f} kb, {g['h']}" + (f", {g['n']} runs" if g["n"] > 1 else ""), 24, "600", P.TEXT, "middle")

P.draw_callouts = callouts_all
sys.argv = [plotter] + rest
P.main()
