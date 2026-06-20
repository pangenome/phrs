#!/usr/bin/env python3
import os
import sys

from common import PACKAGE_DIR, read_tsv


COLORS = {
    "agree": "#2b8cbe",
    "discordant_target": "#d95f0e",
    "inconclusive_partial_overlap": "#756bb1",
    "no_direct_overlap": "#969696",
}


def esc(text):
    return (text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_simple_pdf(path, title, lines):
    def pdf_escape(text):
        return str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    content = ["BT", "/F1 16 Tf", "50 760 Td", "(%s) Tj" % pdf_escape(title), "/F1 9 Tf"]
    y_step = 14
    for line in lines[:48]:
        content.append("0 -%d Td" % y_step)
        content.append("(%s) Tj" % pdf_escape(line[:140]))
    content.append("ET")
    stream = "\n".join(content).encode("ascii", "replace")
    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objects.append(b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 842 595] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream")
    with open(path, "wb") as fh:
        fh.write(b"%PDF-1.4\n")
        offsets = [0]
        for i, obj in enumerate(objects, 1):
            offsets.append(fh.tell())
            fh.write(("%d 0 obj\n" % i).encode("ascii"))
            fh.write(obj)
            fh.write(b"\nendobj\n")
        xref = fh.tell()
        fh.write(("xref\n0 %d\n" % (len(objects) + 1)).encode("ascii"))
        fh.write(b"0000000000 65535 f \n")
        for off in offsets[1:]:
            fh.write(("%010d 00000 n \n" % off).encode("ascii"))
        fh.write(("trailer << /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objects) + 1, xref)).encode("ascii"))


def main():
    path = os.path.join(PACKAGE_DIR, "summaries", "direct_vs_graph_concordance.tsv")
    rows = read_tsv(path)
    width = 1500
    row_h = 24
    left = 360
    height = 90 + row_h * len(rows)
    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (width, height, width, height))
    svg.append('<style>text{font-family:Arial,sans-serif;font-size:12px}.h{font-weight:bold;font-size:14px}.small{font-size:10px}</style>')
    svg.append('<text class="h" x="20" y="28">Direct sweepGA parental haplotype concordance</text>')
    svg.append('<text x="20" y="50">Rows are graph-derived selected segments; bar length is direct PAF overlap with the graph interval.</text>')
    max_len = 1
    for r in rows:
        a, b = [int(x) for x in r["query_local_interval"].split("-")]
        max_len = max(max_len, b - a)
    y = 80
    for i, r in enumerate(rows):
        status = r["direct_support"]
        color = COLORS.get(status, "#969696")
        ov = int(r["overlap_bp"] or 0)
        bar_w = int((width - left - 300) * float(ov) / max_len)
        svg.append('<text x="20" y="%d">%s</text>' % (y + 14, esc(r["event_id"][:42])))
        svg.append('<text class="small" x="20" y="%d">%s | %s</text>' % (y + 28, esc(r["event_role"]), esc(r["expected_target_haplotype"] + " " + r["expected_target_arm"])))
        svg.append('<rect x="%d" y="%d" width="%d" height="14" fill="%s"/>' % (left, y + 4, max(1, bar_w), color))
        svg.append('<text x="%d" y="%d">%s; direct %s %s; ov=%s bp; %s</text>' %
                   (left + bar_w + 8, y + 16, esc(status), esc(r["direct_target_haplotype"]), esc(r["direct_target_arm"]), esc(r["overlap_bp"]), esc(r.get("evidence_source_recommendation", ""))))
        y += row_h
    svg.append('</svg>')
    out_svg = os.path.join(PACKAGE_DIR, "plots", "direct_sweepga_concordance_focused.svg")
    with open(out_svg, "w") as fh:
        fh.write("\n".join(svg) + "\n")
    legacy_svg = os.path.join(PACKAGE_DIR, "plots", "direct_sweepga_concordance_review.svg")
    with open(legacy_svg, "w") as fh:
        fh.write("\n".join(svg) + "\n")

    focused_lines = []
    for r in rows:
        focused_lines.append("%s | %s | %s | expected %s %s | direct %s %s | ov=%s | %s" %
                             (r["event_id"], r["event_role"], r["direct_support"],
                              r["expected_target_haplotype"], r["expected_target_arm"],
                              r["direct_target_haplotype"], r["direct_target_arm"],
                              r["overlap_bp"], r.get("evidence_source_recommendation", "")))
    write_simple_pdf(os.path.join(PACKAGE_DIR, "plots", "direct_sweepga_concordance_focused.pdf"),
                     "Focused direct sweepGA concordance", focused_lines)

    summary = read_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"))
    raw = [r for r in summary if "/raw_paf/" in r["file"]]
    fg_width = 1200
    fg_height = 120 + 38 * len(raw)
    fg = []
    fg.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (fg_width, fg_height, fg_width, fg_height))
    fg.append('<style>text{font-family:Arial,sans-serif;font-size:12px}.h{font-weight:bold;font-size:14px}</style>')
    fg.append('<text class="h" x="20" y="28">Direct sweepGA full-window overview</text>')
    fg.append('<text x="20" y="50">Raw many:many/no-scaffold PAF summaries by child-vs-parent-haplotype comparison.</text>')
    max_bp = max([int(r["query_aligned_bp_sum"]) for r in raw] or [1])
    y = 82
    overview_lines = []
    for r in raw:
        name = os.path.basename(r["file"]).split(".sweepga_", 1)[0]
        total = int(r["query_aligned_bp_sum"])
        inter = int(r["inter_arm_query_bp_sum"])
        bar = int(600 * float(total) / max_bp)
        fg.append('<text x="20" y="%d">%s</text>' % (y + 12, esc(name)))
        fg.append('<rect x="330" y="%d" width="%d" height="14" fill="#3182bd"/>' % (y, max(1, bar)))
        fg.append('<text x="%d" y="%d">%s records; %.1f Mb aligned; %.1f Mb inter-arm</text>' %
                  (340 + bar, y + 12, esc(r["records"]), total / 1000000.0, inter / 1000000.0))
        overview_lines.append("%s | records=%s | query_aligned_bp_sum=%s | inter_arm_query_bp_sum=%s" %
                              (name, r["records"], r["query_aligned_bp_sum"], r["inter_arm_query_bp_sum"]))
        y += 38
    fg.append('</svg>')
    full_svg = os.path.join(PACKAGE_DIR, "plots", "direct_sweepga_full_genome_overview.svg")
    with open(full_svg, "w") as fh:
        fh.write("\n".join(fg) + "\n")
    write_simple_pdf(os.path.join(PACKAGE_DIR, "plots", "direct_sweepga_full_genome_overview.pdf"),
                     "Direct sweepGA full-window overview", overview_lines)
    return 0


if __name__ == "__main__":
    sys.exit(main())
