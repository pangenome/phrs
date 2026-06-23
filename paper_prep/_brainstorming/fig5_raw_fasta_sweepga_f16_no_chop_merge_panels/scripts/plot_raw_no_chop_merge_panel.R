#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
segments_path <- file.path(panel_dir, "raw_merge_panel_segments.tsv")
summary_path <- file.path(panel_dir, "raw_merge_panel_summary.tsv")

segments <- read.delim(segments_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)
mode_order <- c("no_merge_ani", "merge50k_ani", "merge50k_log_length_ani")
summary <- summary[match(
  paste(rep(event_order, each = length(mode_order)), rep(mode_order, times = length(event_order))),
  paste(summary$event_id, summary$filter_mode)
), ]
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
event_label <- function(event_id) {
  switch(
    event_id,
    PAR1_XY_positive_control = "PAR1 X/Y\npositive control",
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027\nchr9q -> chr3q",
    PAN028_chr9q_chr3q_PHR_candidate = "PAN028\nchr9q -> chr3q",
    event_id
  )
}
mode_label <- function(mode) {
  switch(
    mode,
    no_merge_ani = "1:1 no merge\nANI",
    merge50k_ani = "1:1 50 kb merge\nANI",
    merge50k_log_length_ani = "1:1 50 kb merge\nlog-length ANI",
    mode
  )
}
fmt_mb <- function(x) {
  ifelse(x == 0, "0", sprintf("%.1f", x / 1e6))
}

draw_base_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(matrix(seq_len(length(event_order) * length(mode_order)), nrow = length(event_order), byrow = TRUE))
  par(oma = c(3.0, 7.2, 2.5, 0.8), mar = c(2.3, 0.8, 2.0, 0.6), xaxs = "i", yaxs = "i")

  for (event_id in event_order) {
    for (mode in mode_order) {
      srow <- summary[summary$event_id == event_id & summary$filter_mode == mode, ]
      if (nrow(srow) != 1) {
        plot.new()
        next
      }
      x0 <- as.numeric(srow$window_start)
      x1 <- as.numeric(srow$window_end)
      rows <- segments[segments$event_id == event_id & segments$filter_mode == mode, ]
      plot(
        NA,
        xlim = c(x0, x1),
        ylim = c(0, 1),
        axes = FALSE,
        xlab = "",
        ylab = ""
      )
      rect(x0, 0.38, x1, 0.62, col = "#F3F3F3", border = "#CCCCCC")
      if (nrow(rows) > 0) {
        rows <- rows[order(rows$query_clip_start, rows$query_clip_end), ]
        for (j in seq_len(nrow(rows))) {
          fill <- target_col(rows$target_chrom[j])
          border <- if (rows$is_expected_target[j] == "yes") "#111111" else NA
          rect(
            as.numeric(rows$query_clip_start[j]), 0.30,
            as.numeric(rows$query_clip_end[j]), 0.70,
            col = fill,
            border = border,
            lwd = if (rows$is_expected_target[j] == "yes") 0.55 else 0
          )
        }
      }
      axis(1, at = pretty(c(x0, x1), n = 3), labels = fmt_mb(pretty(c(x0, x1), n = 3)), cex.axis = 0.62)
      box(col = "#DDDDDD")
      if (event_id == event_order[1]) {
        title(mode_label(mode), cex.main = 0.8, line = 0.4)
      }
      if (mode == mode_order[1]) {
        mtext(event_label(event_id), side = 2, line = 4.2, cex = 0.78, las = 1)
      }
      text(
        x0 + 0.02 * (x1 - x0), 0.87,
        sprintf("expected rows %s; union %s bp", srow$expected_target_rows, srow$union_expected_overlap_bp),
        adj = c(0, 0.5),
        cex = 0.54,
        col = "#333333"
      )
      if (!is.na(srow$chr3_survival_status) && srow$chr3_survival_status != "NA") {
        text(
          x1 - 0.02 * (x1 - x0), 0.14,
          srow$chr3_survival_status,
          adj = c(1, 0.5),
          cex = 0.54,
          col = if (srow$chr3_survival_status == "CHR3_SURVIVES") "#7A3B00" else "#666666"
        )
      }
    }
  }
  mtext("Absolute query chromosome coordinate (Mb)", side = 1, outer = TRUE, line = 1.4, cex = 0.82)
  mtext("Fig5 raw FASTA SweepGA/FastGA f16 no-chop 1:1 filtering and scaffold-merge comparison", side = 3, outer = TRUE, line = 0.8, cex = 0.95)
  legend(
    "bottom",
    inset = -0.18,
    horiz = TRUE,
    xpd = NA,
    bty = "n",
    cex = 0.66,
    legend = c("chr3 donor", "chr9 context", "chrX", "chrY", "other"),
    fill = unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")]),
    border = NA
  )
}

pdf(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.pdf"), width = 10.0, height = 5.6, useDingbats = FALSE)
draw_base_panel()
dev.off()

svg_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_svg_panel <- function(path) {
  width <- 1000
  height <- 560
  left <- 135
  top <- 70
  cell_w <- 260
  cell_h <- 112
  row_gap <- 28
  col_gap <- 25
  track_h <- 24
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', width, height, width, height),
    '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
    '<text x="500" y="28" text-anchor="middle" font-family="Arial, sans-serif" font-size="14">Fig5 raw FASTA SweepGA/FastGA f16 no-chop 1:1 filtering and scaffold-merge comparison</text>'
  )
  for (ci in seq_along(mode_order)) {
    x <- left + (ci - 1) * (cell_w + col_gap) + cell_w / 2
    label <- gsub("\n", " ", mode_label(mode_order[ci]), fixed = TRUE)
    lines <- c(lines, sprintf('<text x="%.2f" y="55" text-anchor="middle" font-family="Arial, sans-serif" font-size="12">%s</text>', x, svg_escape(label)))
  }
  for (ri in seq_along(event_order)) {
    event_id <- event_order[ri]
    y <- top + (ri - 1) * (cell_h + row_gap)
    lines <- c(lines, sprintf(
      '<text x="122" y="%.2f" text-anchor="end" dominant-baseline="middle" font-family="Arial, sans-serif" font-size="12">%s</text>',
      y + 42, svg_escape(gsub("\n", " ", event_label(event_id), fixed = TRUE))
    ))
    for (ci in seq_along(mode_order)) {
      mode <- mode_order[ci]
      x <- left + (ci - 1) * (cell_w + col_gap)
      srow <- summary[summary$event_id == event_id & summary$filter_mode == mode, ]
      x0 <- as.numeric(srow$window_start)
      x1 <- as.numeric(srow$window_end)
      xscale <- function(v) x + ((v - x0) / (x1 - x0)) * cell_w
      track_y <- y + 32
      lines <- c(lines,
        sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%d" fill="#F3F3F3" stroke="#CCCCCC" stroke-width="1"/>', x, track_y, cell_w, track_h),
        sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%d" fill="none" stroke="#DDDDDD" stroke-width="1"/>', x, y, cell_w, cell_h)
      )
      rows <- segments[segments$event_id == event_id & segments$filter_mode == mode, ]
      if (nrow(rows) > 0) {
        rows <- rows[order(rows$query_clip_start, rows$query_clip_end), ]
        for (j in seq_len(nrow(rows))) {
          x_start <- xscale(as.numeric(rows$query_clip_start[j]))
          x_end <- xscale(as.numeric(rows$query_clip_end[j]))
          fill <- target_col(rows$target_chrom[j])
          stroke <- if (rows$is_expected_target[j] == "yes") "#111111" else fill
          sw <- if (rows$is_expected_target[j] == "yes") "0.7" else "0"
          lines <- c(lines, sprintf(
            '<rect x="%.2f" y="%.2f" width="%.2f" height="32" fill="%s" stroke="%s" stroke-width="%s"/>',
            x_start, track_y - 4, max(0.6, x_end - x_start), fill, stroke, sw
          ))
        }
      }
      tick_vals <- pretty(c(x0, x1), n = 3)
      tick_vals <- tick_vals[tick_vals >= x0 & tick_vals <= x1]
      for (tick in tick_vals) {
        tx <- xscale(tick)
        lines <- c(lines,
          sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#333" stroke-width="1"/>', tx, track_y + track_h + 4, tx, track_y + track_h + 10),
          sprintf('<text x="%.2f" y="%.2f" text-anchor="middle" font-family="Arial, sans-serif" font-size="9">%s</text>', tx, track_y + track_h + 23, svg_escape(fmt_mb(tick)))
        )
      }
      lines <- c(lines, sprintf(
        '<text x="%.2f" y="%.2f" font-family="Arial, sans-serif" font-size="9" fill="#333">expected rows %s; union %s bp</text>',
        x + 7, y + 16, srow$expected_target_rows, srow$union_expected_overlap_bp
      ))
      if (!is.na(srow$chr3_survival_status) && srow$chr3_survival_status != "NA") {
        fill <- if (srow$chr3_survival_status == "CHR3_SURVIVES") "#7A3B00" else "#666666"
        lines <- c(lines, sprintf(
          '<text x="%.2f" y="%.2f" text-anchor="end" font-family="Arial, sans-serif" font-size="9" fill="%s">%s</text>',
          x + cell_w - 7, y + cell_h - 10, fill, srow$chr3_survival_status
        ))
      }
    }
  }
  lines <- c(lines, '<text x="525" y="492" text-anchor="middle" font-family="Arial, sans-serif" font-size="11">Absolute query chromosome coordinate (Mb)</text>')
  legend <- data.frame(
    label = c("chr3 donor", "chr9 context", "chrX", "chrY", "other"),
    col = unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")]),
    stringsAsFactors = FALSE
  )
  lx <- 250
  for (i in seq_len(nrow(legend))) {
    x <- lx + (i - 1) * 110
    lines <- c(lines,
      sprintf('<rect x="%d" y="515" width="14" height="10" fill="%s"/>', x, legend$col[i]),
      sprintf('<text x="%d" y="524" font-family="Arial, sans-serif" font-size="10">%s</text>', x + 20, svg_escape(legend$label[i]))
    )
  }
  lines <- c(lines, '</svg>')
  writeLines(lines, path)
}

write_svg_panel(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.svg"))
