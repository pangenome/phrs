#!/usr/bin/env python3
"""Fig. 3 panels: the May 2026 UCSC panels with the cut side label "EDEF Segmental Dups" completed to "SEDEF Segmental Dups".

usage (repo root): python3 submission/scripts/figures/fix_fig3_sedef_label.py
UCSC drew the 20-character label in a 20-character column, so the first "S" falls left of the image (only its right 4 columns
remain). A new UCSC render draws the text heavier (2026-11-27), so we keep the May pixels and change only that letter:
in: slides/chm13-phr-ucsc-browser/_assets/ucsc/panels/<window>.png (May render, 1 px #cfcfcf border; the paper used them
  without the border); we drop the border, add M = 6 white columns on the left, and copy the "S" of the side label "Scale"
  (same image, same font and size: columns 107-115, rows 5-14) to the "S" of "SEDEF" (rows 85-94, 111 columns to the left).
Check: the 4 visible columns of the cut "S" equal the matching columns of the "Scale" S, pixel for pixel, in every panel.
out: submission/fig/MainFigures/<panel>.png (no border, as before); log on stderr.
"""
import sys
from PIL import Image
M, DX, DY = 6, 111, 80
P = [("Fig3a_C1_chr4q", "chr4q_chr4_193304946_193574945"), ("Fig3a_C1_chr10q", "chr10q_chr10_134488135_134758134"),
     ("Fig3b_C2_chr10p", "chr10p_chr10_1_105002"), ("Fig3b_C2_chr18p", "chr18p_chr18_1_397502"),
     ("Fig3c_C11_chr5q", "chr5q_chr5_181752940_182045439"), ("Fig3c_C11_chr6q", "chr6q_chr6_171909129_172126628")]
for name, win in P:
    a = Image.open(f"slides/chm13-phr-ucsc-browser/_assets/ucsc/panels/{win}.png").convert("RGB")
    im = a.crop((1, 1, a.width - 1, a.height - 1))
    bad = [(x, y) for x in range(112, 116) for y in range(5, 15) if im.getpixel((x, y)) != im.getpixel((x - DX, y + DY))]
    if bad: sys.exit(f"ERROR {name}: visible part of the cut S differs from the Scale S at {bad[:5]}")
    out = Image.new("RGB", (im.width + M, im.height), (255, 255, 255)); out.paste(im, (M, 0))
    for x in range(107, 116):
        for y in range(5, 15): out.putpixel((x - DX + M, y + DY), im.getpixel((x, y)))
    out.save(f"submission/fig/MainFigures/{name}.png")
    print(f"{name}: May {a.width}x{a.height} -> without border {im.width}x{im.height} -> {out.width}x{out.height}, S completed", file=sys.stderr)
