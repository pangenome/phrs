#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
segments <- read.delim(file.path(panel_dir, "raw_merge_panel_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(panel_dir, "raw_merge_panel_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(summary) == 0) {
  stop("no summary rows to plot")
}

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)
mode_order <- c("raw_no_merge_ani", "raw_merge50k_ani", "raw_merge50k_log_length_ani")

summary$event_id <- factor(summary$event_id, levels = event_order)
summary$filter_mode <- factor(summary$filter_mode, levels = mode_order)
summary <- summary[order(summary$filter_mode, summary$event_id), ]
if (nrow(segments) > 0) {
  segments$event_id <- factor(segments$event_id, levels = event_order)
  segments$filter_mode <- factor(segments$filter_mode, levels = mode_order)
}

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
    as.character(event_id),
    PAR1_XY_positive_control = "PAR1 X/Y positive control",
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
    PAN028_chr9q_chr3q_PHR_candidate = "PAN028 chr9q -> chr3q",
    as.character(event_id)
  )
}
mode_label_short <- function(mode_id) {
  switch(
    as.character(mode_id),
    raw_no_merge_ani = "no chop\nno merge; ANI",
    raw_merge50k_ani = "no chop\n50 kb merge; ANI",
    raw_merge50k_log_length_ani = "no chop\n50 kb merge; log-length ANI",
    as.character(mode_id)
  )
}
fmt_coord <- function(x) {
  if (x >= 1e6) sprintf("%.2f Mb", x / 1e6) else sprintf("%d kb", round(x / 1000))
}
draw_legend <- function(x0, y0) {
  labels <- c("chr3 donor", "chr9 context", "chrX", "chrY", "other")
  cols <- unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")])
  x_pos <- c(x0, x0 + 0.22, x0 + 0.41, x0 + 0.58, x0 + 0.75)
  for (i in seq_along(labels)) {
    x <- x_pos[i]
    rect(x, y0 - 0.08, x + 0.018, y0 + 0.08, col = cols[i], border = NA, xpd = NA)
    text(x + 0.024, y0, labels[i], adj = 0, cex = 0.72, xpd = NA)
  }
}

draw_page <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(2.0, 8.0, 3.0, 2.2), oma = c(2.0, 0.2, 1.0, 0.2), xaxs = "i", yaxs = "i")
  n_tracks <- length(mode_order) * length(event_order)
  plot(NA, xlim = c(-0.38, 1.30), ylim = c(0.3, n_tracks + 1.2), axes = FALSE, xlab = "", ylab = "")
  title("Raw FASTA f16 whole-genome PAF -> 1:1 filtered without chopping", cex.main = 0.92)
  draw_legend(0.08, n_tracks + 0.82)

  row_i <- 0
  for (mode_id in mode_order) {
    mode_summary <- summary[summary$filter_mode == mode_id, ]
    mode_label <- mode_summary$filter_label[1]
    for (event_id in event_order) {
      row_i <- row_i + 1
      y <- n_tracks - row_i + 1
      row <- mode_summary[as.character(mode_summary$event_id) == event_id, ]
      if (nrow(row) == 0) next
      w_start <- as.numeric(row$window_start[1])
      w_end <- as.numeric(row$window_end[1])
      xscale <- function(x) (x - w_start) / (w_end - w_start)
      rect(0, y - 0.16, 1, y + 0.16, col = "#F4F4F4", border = "#CCCCCC")
      rows <- segments[as.character(segments$filter_mode) == mode_id & as.character(segments$event_id) == event_id, ]
      if (nrow(rows) > 0) {
        rows <- rows[order(rows$query_clip_start, rows$query_clip_end), ]
        for (j in seq_len(nrow(rows))) {
          x1 <- max(0, min(1, xscale(as.numeric(rows$query_clip_start[j]))))
          x2 <- max(0, min(1, xscale(as.numeric(rows$query_clip_end[j]))))
          col <- target_col(rows$target_chrom[j])
          border <- if (rows$is_expected_target[j] == "yes") "#111111" else NA
          rect(x1, y - 0.22, x2, y + 0.22, col = col, border = border, lwd = if (rows$is_expected_target[j] == "yes") 0.55 else 0)
        }
      }
      label <- sprintf("%s\n%s", event_label(event_id), mode_label_short(mode_id))
      text(-0.035, y, label, adj = 1, cex = 0.62)
      text(1.015, y, sprintf("exp %s; %.1f kb", row$expected_target_rows[1], as.numeric(row$union_expected_overlap_bp[1]) / 1000), adj = 0, cex = 0.58, col = "#333333")
      ticks <- pretty(c(w_start, w_end), n = 4)
      ticks <- ticks[ticks >= w_start & ticks <= w_end]
      for (tick in ticks) {
        x <- xscale(tick)
        segments(x, y - 0.36, x, y - 0.29, col = "#444444", lwd = 0.6)
        text(x, y - 0.45, fmt_coord(tick), cex = 0.48, adj = c(0.5, 1))
      }
    }
  }
  mtext("Absolute query chromosome coordinate per 500 kb display window", side = 1, outer = TRUE, line = 0.7, cex = 0.8)
}

pdf(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.pdf"), width = 12.5, height = 7.6, useDingbats = FALSE)
draw_page()
dev.off()

svg(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.svg"), width = 12.5, height = 7.6, onefile = TRUE)
draw_page()
dev.off()

preview_dir <- file.path(panel_dir, "preview_png")
dir.create(preview_dir, showWarnings = FALSE)
png(file.path(preview_dir, "fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.png"), width = 2500, height = 1520, res = 200, type = "cairo")
draw_page()
dev.off()
