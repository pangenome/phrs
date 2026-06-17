#!/usr/bin/env python3
"""Build a community-sorted CHM13/hs1 UCSC browser review deck.

This reuses the cached UCSC browser panels from render_ucsc_browser_suite.py and
emits a compact Typst deck with full-width panels stacked when they fit.
"""

import csv
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
SUITE_DIR = REPO_ROOT / "slides" / "chm13-phr-ucsc-browser"
MANIFEST = SUITE_DIR / "manifest.tsv"
ASSIGNMENTS = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/similarity/"
    "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
)
COMMUNITY_ORDER = (
    REPO_ROOT
    / "slides"
    / "v2-review-zoom"
    / "_revision_assets"
    / "v5"
    / "07a_tree_then_community_heatmap"
    / "arm_order_community.tsv"
)

OUT_TYP = SUITE_DIR / "chm13_phr_ucsc_browser_suite_by_community.typ"
OUT_LAYOUT = SUITE_DIR / "manifest_by_community.tsv"

CONTENT_WIDTH_IN = 13.33 - (2 * 0.18)
CONTENT_HEIGHT_IN = 7.5 - (2 * 0.14)
PANEL_TEXT_IN = 0.65
PANEL_GAP_IN = 0.06


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def community_number(name):
    return int(name[1:])


def arm_sort_key(chrom_arm):
    chrom, arm = chrom_arm.split("_")
    tag = chrom.replace("chr", "")
    if tag == "X":
        chrom_num = 23
    elif tag == "Y":
        chrom_num = 24
    else:
        chrom_num = int(tag)
    return chrom_num, arm


def escape_typst_string(value):
    return value.replace("\\", "\\\\").replace('"', '\\"')


def relative_image_path(path_text):
    prefix = "slides/chm13-phr-ucsc-browser/"
    if path_text.startswith(prefix):
        return path_text[len(prefix) :]
    return path_text


def fmt_int(value):
    return "{:,}".format(int(value))


def fmt_kbp(value):
    bp = int(value)
    kbp = bp / 1000.0
    if kbp >= 100:
        return "{:.0f} kbp".format(kbp)
    if kbp >= 10:
        return "{:.1f} kbp".format(kbp).replace(".0 ", " ")
    return "{:.2f} kbp".format(kbp).replace(".00 ", " ")


def panel_height_in(row):
    width = float(row["image_width_px"])
    height = float(row["image_height_px"])
    return (height / width * CONTENT_WIDTH_IN) + PANEL_TEXT_IN


def paginate(rows):
    pages = []
    current = []
    used = 0.0
    for row in rows:
        block_h = panel_height_in(row)
        gap = PANEL_GAP_IN if current else 0.0
        if current and used + gap + block_h > CONTENT_HEIGHT_IN:
            pages.append(current)
            current = []
            used = 0.0
            gap = 0.0
        current.append(row)
        used += gap + block_h
    if current:
        pages.append(current)
    return pages


def make_typst_panel(row):
    label = row["label"]
    chrom_arm = row["chrom_arm"]
    community = row["Community"]
    members = row["Arms"]
    region = row["ucsc_position_display"]
    details = "PHR {}; window {}; terminal gap {}; C arms {}".format(
        fmt_kbp(row["phr_inclusive_bp"]),
        fmt_kbp(row["browser_window_bp"]),
        fmt_kbp(row["terminal_gap_bp"]),
        members,
    )
    image_path = relative_image_path(row["image_path"])
    url = row["ucsc_url"]
    return (
        '#browser-panel("{}", "{}", "{}", "{}", "{}", "{}", "{}")\n'.format(
            escape_typst_string(community),
            escape_typst_string(label),
            escape_typst_string(region),
            escape_typst_string(details),
            escape_typst_string(members),
            escape_typst_string(image_path),
            escape_typst_string(url),
        )
    )


def title_page(rows, assignments):
    by_community = {}
    rendered_by_community = {}
    for row in assignments:
        by_community.setdefault(row["Community"], []).append(row["ChromArm"])
    for row in rows:
        rendered_by_community.setdefault(row["Community"], []).append(row["chrom_arm"])

    missing_lines = []
    for community in sorted(by_community, key=community_number):
        missing = sorted(
            set(by_community[community]) - set(rendered_by_community.get(community, [])),
            key=arm_sort_key,
        )
        if missing:
            missing_lines.append(
                "{} missing from chm13.phrs.bed: {}".format(
                    community, ", ".join(m.replace("_", "") for m in missing)
                )
            )
    missing_text = "; ".join(missing_lines) if missing_lines else "none"
    return """
#align(center + horizon)[
  #block(width: 11.4in)[
    #align(center)[#text(size: 25pt, weight: "bold", fill: col-title)[CHM13/hs1 PHR Browser Panels by Community]]
    #v(0.22in)
    #text(size: 13pt)[37 cached UCSC Genome Browser panels sorted by the HPRCv2 arm-level Leiden C1-C15 community assignments. Panels keep the same full-width image scale as the original browser suite; multiple panels are stacked only when they fit at that scale.]
    #v(0.16in)
    #text(size: 10.5pt, fill: col-muted)[Every panel includes a visible, clickable UCSC URL. Browser windows are terminal-anchored 1.5x zoom-outs around the CHM13 PHR interval.]
    #v(0.16in)
    #text(size: 10pt, fill: col-muted)[Rendered arms absent from the CHM13 BED remain absent here: #missing.]
    #v(0.16in)
    #text(size: 9.2pt, fill: col-muted)[Sources: `manifest.tsv`; `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`; `arm_order_community.tsv`. Layout audit: `manifest_by_community.tsv`.]
  ]
]
""".replace("#missing", escape_typst_string(missing_text))


def write_typst(rows, pages):
    lines = [
        "// Community-sorted CHM13/hs1 UCSC browser slide suite.",
        "// Generated by _scripts/render_ucsc_browser_suite_by_community.py.",
        "",
        "#set page(",
        "  width: 13.33in,",
        "  height: 7.5in,",
        "  margin: (x: 0.18in, y: 0.14in),",
        ")",
        '#set text(font: "DejaVu Sans", size: 9pt, lang: "en")',
        "#set par(justify: false, leading: 0.62em)",
        "",
        '#let col-title = rgb("#1a3a6b")',
        '#let col-muted = rgb("#555555")',
        '#let col-link = rgb("#135f9f")',
        '#let col-rule = rgb("#ccd6e0")',
        '#let col-chip = rgb("#e7eef7")',
        "",
        "#let browser-panel(community, label, region, details, members, img, url) = {",
        "  box(width: 100%)[",
        "    #grid(",
        "      columns: (0.48in, 0.78in, 1fr, auto),",
        "      column-gutter: 0.08in,",
        "      align: (left, left, left, right),",
        "      [#box(inset: (x: 0.05in, y: 0.012in), radius: 1.5pt, fill: col-chip)[#text(size: 8pt, weight: \"bold\", fill: col-title)[#community]]],",
        "      [#text(size: 10pt, weight: \"semibold\", fill: col-title)[#label]],",
        "      [#text(size: 7.1pt, fill: col-muted)[#details]],",
        "      [#text(size: 7.3pt, fill: col-muted)[#region]],",
        "    )",
        "    #v(0.006in)",
        "    #text(size: 5.25pt, fill: col-link)[#link(url)[#url]]",
        "    #v(0.01in)",
        "    #image(img, width: 100%)",
        "  ]",
        "}",
        "",
    ]

    assignments = read_tsv(ASSIGNMENTS)
    lines.append(title_page(rows, assignments))
    lines.append("#pagebreak()\n")
    for page_index, page_rows in enumerate(pages, start=1):
        for panel_index, row in enumerate(page_rows):
            if panel_index:
                lines.append("#v(0.07in)\n")
            lines.append(make_typst_panel(row))
        if page_index != len(pages):
            lines.append("#pagebreak()\n")
    OUT_TYP.write_text("\n".join(lines))


def write_layout(rows):
    fields = [
        "layout_page",
        "layout_panel_on_page",
        "layout_height_in",
        "Community",
        "label",
        "chrom_arm",
        "ucsc_position_display",
        "phr_inclusive_bp",
        "browser_window_bp",
        "terminal_gap_bp",
        "image_path",
        "ucsc_url",
    ]
    with OUT_LAYOUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {field: row.get(field, "") for field in fields}
            out["layout_height_in"] = "{:.3f}".format(panel_height_in(row))
            writer.writerow(out)


def main():
    manifest_rows = read_tsv(MANIFEST)
    assignments = {row["ChromArm"]: row for row in read_tsv(ASSIGNMENTS)}
    order_rows = read_tsv(COMMUNITY_ORDER)
    order = {row["ChromArm"]: int(row["position_left_to_right"]) for row in order_rows}

    enriched = []
    for row in manifest_rows:
        chrom_arm = row["chrom_arm"]
        assignment = assignments[chrom_arm]
        merged = dict(row)
        merged.update(assignment)
        merged["community_num"] = community_number(assignment["Community"])
        merged["community_order"] = order.get(chrom_arm, 10_000)
        enriched.append(merged)

    enriched.sort(
        key=lambda row: (
            row["community_num"],
            row["community_order"],
            arm_sort_key(row["chrom_arm"]),
        )
    )

    pages = paginate(enriched)
    flat = []
    for page_num, page_rows in enumerate(pages, start=1):
        for panel_num, row in enumerate(page_rows, start=1):
            row["layout_page"] = str(page_num)
            row["layout_panel_on_page"] = str(panel_num)
            flat.append(row)

    write_typst(flat, pages)
    write_layout(flat)
    print("wrote {}".format(OUT_TYP))
    print("wrote {}".format(OUT_LAYOUT))
    print("content pages: {}".format(len(pages)))


if __name__ == "__main__":
    main()
