#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
segments_path <- file.path(panel_dir, "raw_fasta_chopped_panel_segments.tsv")
summary_path <- file.path(panel_dir, "raw_fasta_chopped_panel_summary.tsv")

segments <- read.delim(segments_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(segments) == 0) {
  stop("no segments to plot")
}

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)
summary <- summary[match(event_order, summary$event_id), ]
summary <- summary[!is.na(summary$event_id), ]

target_palette <- c(
  chr3 = "#D95F02",
  chr9 = "#1B9E77",
  chrX = "#7570B3",
  chrY = "#E7298A",
  other = "#8A8A8A"
)
target_col <- function(x) {
  ifelse(x %in% names(target_palette), target_palette[x], target_palette[["other"]])
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.0, 8.2, 2.4, 1.0), oma = c(0.8, 0, 1.5, 0), xaxs = "i", yaxs = "i")
  plot(
    NA,
    xlim = c(0, 500000),
    ylim = c(0.35, nrow(summary) + 0.75),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(1, at = seq(0, 500000, by = 100000), labels = seq(0, 500, by = 100), cex.axis = 0.75)
  mtext("Query-window coordinate (kb)", side = 1, line = 2.0, cex = 0.8)
  title("Fig5 raw FASTA SweepGA/FastGA f16 evidence: 2 kb chopped, 1:1 ANI-filtered PAF", cex.main = 0.9)

  for (i in seq_len(nrow(summary))) {
    event_id <- summary$event_id[i]
    y <- nrow(summary) - i + 1
    rows <- segments[segments$event_id == event_id, ]
    rect(0, y - 0.12, 500000, y + 0.12, col = "#F3F3F3", border = "#CCCCCC")
    if (nrow(rows) > 0) {
      rows <- rows[order(rows$query_rel_start, rows$query_rel_end), ]
      for (j in seq_len(nrow(rows))) {
        col <- target_col(rows$target_chrom[j])
        border <- if (rows$is_expected_target[j] == "yes") "#111111" else NA
        rect(
          rows$query_rel_start[j],
          y - 0.18,
          rows$query_rel_end[j],
          y + 0.18,
          col = col,
          border = border,
          lwd = if (rows$is_expected_target[j] == "yes") 0.55 else 0
        )
      }
    }
    label <- sub(" -> ", " ->\n", summary$event_id[i], fixed = TRUE)
    label <- switch(
      event_id,
      PAR1_XY_positive_control = "PAR1 X/Y\npositive control",
      PAN027_chr9q_chr3q_PHR_candidate = "PAN027\nchr9q -> chr3q",
      PAN028_chr9q_chr3q_PHR_candidate = "PAN028\nchr9q -> chr3q",
      label
    )
    text(-25000, y, label, adj = 1, cex = 0.78, xpd = NA)
    text(
      505000,
      y,
      paste0("expected rows: ", summary$expected_target_rows[i]),
      adj = 0,
      cex = 0.62,
      xpd = NA,
      col = "#333333"
    )
  }

  legend(
    "bottom",
    inset = -0.34,
    horiz = TRUE,
    xpd = NA,
    bty = "n",
    cex = 0.72,
    legend = c("chr3 donor", "chr9 context", "chrX", "chrY", "other"),
    fill = unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")]),
    border = NA
  )
}

pdf(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_chopped_panels.pdf"), width = 8.2, height = 3.5, useDingbats = FALSE)
draw_panel()
dev.off()

svg_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_svg_panel <- function(path) {
  width <- 820
  height <- 350
  left <- 185
  right <- 75
  top <- 58
  track_w <- width - left - right
  row_gap <- 62
  y0 <- top + 24
  xscale <- function(x) left + (x / 500000) * track_w
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', width, height, width, height),
    '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
    '<text x="410" y="28" text-anchor="middle" font-family="Arial, sans-serif" font-size="13">Fig5 raw FASTA SweepGA/FastGA f16 evidence: 2 kb chopped, 1:1 ANI-filtered PAF</text>'
  )
  for (tick in seq(0, 500000, by = 100000)) {
    x <- xscale(tick)
    lines <- c(lines,
      sprintf('<line x1="%.2f" y1="286" x2="%.2f" y2="292" stroke="#333" stroke-width="1"/>', x, x),
      sprintf('<text x="%.2f" y="307" text-anchor="middle" font-family="Arial, sans-serif" font-size="10">%d</text>', x, tick / 1000)
    )
  }
  lines <- c(lines, '<text x="465" y="330" text-anchor="middle" font-family="Arial, sans-serif" font-size="11">Query-window coordinate (kb)</text>')
  for (i in seq_len(nrow(summary))) {
    event_id <- summary$event_id[i]
    y <- y0 + (i - 1) * row_gap
    label <- switch(
      event_id,
      PAR1_XY_positive_control = "PAR1 X/Y positive control",
      PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
      PAN028_chr9q_chr3q_PHR_candidate = "PAN028 chr9q -> chr3q",
      event_id
    )
    lines <- c(lines,
      sprintf('<text x="172" y="%.2f" text-anchor="end" dominant-baseline="middle" font-family="Arial, sans-serif" font-size="12">%s</text>', y, svg_escape(label)),
      sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="20" fill="#F3F3F3" stroke="#CCCCCC" stroke-width="1"/>', left, y - 10, track_w)
    )
    rows <- segments[segments$event_id == event_id, ]
    if (nrow(rows) > 0) {
      rows <- rows[order(rows$query_rel_start, rows$query_rel_end), ]
      for (j in seq_len(nrow(rows))) {
        x1 <- xscale(rows$query_rel_start[j])
        x2 <- xscale(rows$query_rel_end[j])
        fill <- target_col(rows$target_chrom[j])
        stroke <- if (rows$is_expected_target[j] == "yes") "#111111" else fill
        sw <- if (rows$is_expected_target[j] == "yes") "0.7" else "0"
        lines <- c(lines, sprintf(
          '<rect x="%.2f" y="%.2f" width="%.2f" height="28" fill="%s" stroke="%s" stroke-width="%s"/>',
          x1, y - 14, max(0.6, x2 - x1), fill, stroke, sw
        ))
      }
    }
    lines <- c(lines, sprintf(
      '<text x="755" y="%.2f" dominant-baseline="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">expected rows: %s</text>',
      y, summary$expected_target_rows[i]
    ))
  }
  legend <- data.frame(
    label = c("chr3 donor", "chr9 context", "chrX", "chrY", "other"),
    col = unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")]),
    stringsAsFactors = FALSE
  )
  lx <- 210
  for (i in seq_len(nrow(legend))) {
    x <- lx + (i - 1) * 105
    lines <- c(lines,
      sprintf('<rect x="%d" y="318" width="14" height="10" fill="%s"/>', x, legend$col[i]),
      sprintf('<text x="%d" y="327" font-family="Arial, sans-serif" font-size="10">%s</text>', x + 20, svg_escape(legend$label[i]))
    )
  }
  lines <- c(lines, '</svg>')
  writeLines(lines, path)
}

write_svg_panel(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_chopped_panels.svg"))
