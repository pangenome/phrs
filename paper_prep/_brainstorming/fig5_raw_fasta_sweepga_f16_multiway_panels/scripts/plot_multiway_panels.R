#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
panel_dir <- if (length(args) >= 1) args[[1]] else dirname(dirname(normalizePath(sys.frame(1)$ofile %||% ".")))

support_path <- file.path(panel_dir, "multiway_candidate_support.tsv")
summary_path <- file.path(panel_dir, "multiway_candidate_summary.tsv")
support <- read.delim(support_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)
layer_order <- c(
  "raw_many_many_whole_genome",
  "chopped_many_many",
  "chopped_four_many",
  "chopped_one_one"
)

support <- support[support$event_id %in% event_order & support$source_layer %in% layer_order, ]
summary <- summary[summary$event_id %in% event_order & summary$source_layer %in% layer_order, ]
support$event_id <- factor(support$event_id, levels = event_order)
support$source_layer <- factor(support$source_layer, levels = layer_order)
summary$event_id <- factor(summary$event_id, levels = event_order)
summary$source_layer <- factor(summary$source_layer, levels = layer_order)

target_palette <- c(
  chr3 = "#D55E00",
  chr9 = "#009E73",
  chrX = "#0072B2",
  chrY = "#CC79A7",
  other = "#7F7F7F"
)

target_col <- function(bucket) {
  out <- target_palette["other"]
  if (bucket %in% names(target_palette)) {
    out <- target_palette[bucket]
  }
  unname(out)
}

event_label <- function(event_id) {
  switch(
    as.character(event_id),
    PAR1_XY_positive_control = "PAR1 X/Y control",
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
    PAN028_chr9q_chr3q_PHR_candidate = "PAN028 chr9q -> chr3q",
    as.character(event_id)
  )
}

layer_label <- function(layer) {
  switch(
    as.character(layer),
    raw_many_many_whole_genome = "raw many:many\nwhole-genome",
    chopped_many_many = "chopped\nmany:many",
    chopped_four_many = "chopped\nfour:many",
    chopped_one_one = "chopped\none:one",
    as.character(layer)
  )
}

fmt_coord <- function(x) {
  sprintf("%.3f Mb", as.numeric(x) / 1e6)
}

draw_legend <- function(x0, y0) {
  labels <- c("chr3", "chr9", "chrX", "chrY", "other")
  for (i in seq_along(labels)) {
    x <- x0 + (i - 1) * 0.16
    rect(x, y0 - 0.09, x + 0.022, y0 + 0.09, col = target_col(labels[i]), border = NA, xpd = NA)
    text(x + 0.028, y0, labels[i], adj = 0, cex = 0.68, xpd = NA)
  }
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(2.4, 7.2, 3.0, 4.4), oma = c(1.4, 0.2, 1.0, 0.2), xaxs = "i", yaxs = "i")
  n_tracks <- length(event_order) * length(layer_order)
  plot(
    NA,
    xlim = c(-0.40, 1.42),
    ylim = c(0.2, n_tracks + 1.45),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  title("Fig5 raw FASTA SweepGA f16 multiway support before and after chopping/filtering", cex.main = 0.92)
  text(-0.40, n_tracks + 1.18, "Absolute query chromosome coordinates", adj = 0, cex = 0.72)
  draw_legend(0.22, n_tracks + 1.18)

  row_i <- 0
  for (event_id in event_order) {
    event_summary <- summary[as.character(summary$event_id) == event_id, ]
    if (nrow(event_summary) == 0) {
      next
    }
    w_start <- as.numeric(event_summary$window_start[1])
    w_end <- as.numeric(event_summary$window_end[1])
    xscale <- function(x) (as.numeric(x) - w_start) / (w_end - w_start)
    for (layer in layer_order) {
      row_i <- row_i + 1
      y <- n_tracks - row_i + 1
      row_summary <- event_summary[as.character(event_summary$source_layer) == layer, ]
      rect(0, y - 0.18, 1, y + 0.18, col = "#F4F4F4", border = "#C8C8C8", lwd = 0.6)
      rows <- support[
        as.character(support$event_id) == event_id &
          as.character(support$source_layer) == layer,
      ]
      if (nrow(rows) > 0) {
        rows <- rows[order(as.numeric(rows$query_clip_start), rows$target_bucket, decreasing = FALSE), ]
        lane_count <- min(8, max(1, ceiling(log10(nrow(rows) + 1))))
        for (j in seq_len(nrow(rows))) {
          x1 <- max(0, min(1, xscale(rows$query_clip_start[j])))
          x2 <- max(0, min(1, xscale(rows$query_clip_end[j])))
          lane <- ((j - 1) %% lane_count) + 1
          y1 <- y - 0.24 + (lane - 1) * (0.48 / lane_count)
          y2 <- y - 0.24 + lane * (0.48 / lane_count) - 0.01
          border <- if (rows$is_expected_target[j] == "yes") "#111111" else NA
          rect(x1, y1, x2, y2, col = target_col(rows$target_bucket[j]), border = border, lwd = 0.45)
        }
      }
      if (layer == layer_order[1]) {
        text(-0.055, y + 0.18, event_label(event_id), adj = 1, cex = 0.72, font = 2)
      }
      text(-0.055, y - 0.08, layer_label(layer), adj = 1, cex = 0.56)
      if (nrow(row_summary) > 0) {
        right_label <- sprintf(
          "rows %s | union %.1f kb | mult %s",
          row_summary$row_count[1],
          as.numeric(row_summary$query_union_coverage_bp[1]) / 1000,
          row_summary$row_multiplicity_per_query_bp[1]
        )
        text(1.02, y, right_label, adj = 0, cex = 0.50)
      }
    }
    ticks <- pretty(c(w_start, w_end), n = 4)
    ticks <- ticks[ticks >= w_start & ticks <= w_end]
    axis_y <- n_tracks - row_i + 0.44
    for (tick in ticks) {
      x <- xscale(tick)
      segments(x, axis_y, x, axis_y + 0.08, col = "#444444", lwd = 0.5)
      text(x, axis_y - 0.06, fmt_coord(tick), cex = 0.45, adj = c(0.5, 1))
    }
  }
  mtext("Raw many:many rows are the source-of-truth multiway support; chopped many:many/four:many/one:one rows are comparison layers.", side = 1, outer = TRUE, line = 0.4, cex = 0.72)
}

pdf(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_multiway_panels.pdf"), width = 13.2, height = 8.2, useDingbats = FALSE)
draw_panel()
dev.off()

svg(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_multiway_panels.svg"), width = 13.2, height = 8.2, onefile = TRUE)
draw_panel()
dev.off()

png(file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_multiway_panels.png"), width = 2640, height = 1640, res = 200, type = "cairo")
draw_panel()
dev.off()

png(file.path(panel_dir, "preview_png", "fig5_raw_fasta_sweepga_f16_multiway_panels.png"), width = 2640, height = 1640, res = 200, type = "cairo")
draw_panel()
dev.off()
