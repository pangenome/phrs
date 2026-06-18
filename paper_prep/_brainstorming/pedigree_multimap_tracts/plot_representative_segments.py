#!/usr/bin/env python3
import csv
import html
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "representative_segments.tsv").open(), delimiter="\t"))
if not rows:
    raise SystemExit("representative_segments.tsv is empty")

colors = {
    "unique_donor_haplotype": "#1b9e77",
    "unique_donor_arm_multi_haplotype": "#7570b3",
    "same_community_ambiguous": "#e6ab02",
    "cross_community_ambiguous": "#d95f02",
    "unresolved_no_call": "#666666",
}
keys = []
for row in rows:
    key = (row["setting"], row["pair"], row["query_name"], row["tract_start"], row["tract_end"])
    if key not in keys:
        keys.append(key)

width = 1200
left = 90
right = 500
row_h = 28
top = 52
height = top + row_h * len(keys) + 70
min_x = min(int(row["tract_start"]) for row in rows)
max_x = max(int(row["tract_end"]) for row in rows)
span = max(1, max_x - min_x)

def sx(value):
    return left + (int(value) - min_x) / span * (width - left - right)

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
    "<style>text{font-family:Arial,sans-serif;font-size:11px}.title{font-size:16px;font-weight:bold}.axis{fill:#333}.label{font-size:9px}</style>",
    f'<text class="title" x="{left}" y="24">Representative multimap-aware WashU untangle tract calls</text>',
    f'<line x1="{left}" x2="{width-right}" y1="{height-42}" y2="{height-42}" stroke="#333" stroke-width="1"/>',
    f'<text class="axis" x="{left}" y="{height-20}">child/query coordinate within 500 kb flank (bp)</text>',
]
for y, key in enumerate(keys):
    group = [row for row in rows if (row["setting"], row["pair"], row["query_name"], row["tract_start"], row["tract_end"]) == key]
    tract_start = int(key[3])
    tract_end = int(key[4])
    cy = top + y * row_h
    parts.append(f'<line x1="{sx(tract_start):.1f}" x2="{sx(tract_end):.1f}" y1="{cy}" y2="{cy}" stroke="#222" stroke-width="1" opacity="0.5"/>')
    for row in group:
        start = int(row["segment_start"])
        end = int(row["segment_end"])
        cls = row["segment_resolvability"]
        parts.append(
            f'<rect x="{sx(start):.1f}" y="{cy-6}" width="{max(1, sx(end)-sx(start)):.1f}" height="12" '
            f'fill="{colors.get(cls, "#999999")}"><title>{html.escape(cls)}; {start}-{end}</title></rect>'
        )
    label = f"{key[0]} {key[1]}:{key[2].split('#')[-1]} {tract_end - tract_start} bp"
    parts.append(f'<text class="label" x="{width-right+12}" y="{cy+4}">{html.escape(label)}</text>')

legend_y = height - 60
legend_x = left
for cls, color in colors.items():
    parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="12" height="12" fill="{color}"/>')
    parts.append(f'<text class="label" x="{legend_x+17}" y="{legend_y+10}">{html.escape(cls)}</text>')
    legend_x += 225
    if legend_x > width - 260:
        legend_x = left
        legend_y += 18
parts.append("</svg>")
(HERE / "representative_segments.svg").write_text("\n".join(parts))
