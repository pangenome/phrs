#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
prefix <- "fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels"
segments_path <- file.path(panel_dir, "query_grid_panel_segments.tsv")
summary_path <- file.path(panel_dir, "query_grid_panel_summary.tsv")

segments <- read.delim(segments_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(summary) == 0) {
  stop("no summary rows to plot")
}

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)
length_order <- c(10000, 5000, 2000)
summary$event_id <- factor(summary$event_id, levels = event_order)
summary$chop_length_bp <- as.integer(summary$chop_length_bp)
summary <- summary[order(summary$event_id, summary$chop_length_bp), ]

if (nrow(segments) > 0) {
  segments$event_id <- factor(segments$event_id, levels = event_order)
  segments$chop_length_bp <- as.integer(segments$chop_length_bp)
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
    PAR1_XY_positive_control = "PAR1 X/Y control",
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
    PAN028_chr9q_chr3q_PHR_candidate = "PAN028 chr9q -> chr3q",
    as.character(event_id)
  )
}

fmt_coord <- function(x) {
  ifelse(x >= 1e6, sprintf("%.2f Mb", x / 1e6), sprintf("%d kb", round(x / 1000)))
}

fmt_len <- function(x) {
  ifelse(x >= 1000, paste0(x / 1000, " kb"), paste0(x, " bp"))
}

draw_legend <- function() {
  par(xpd = NA)
  legend_labels <- c("chr3 donor", "chr9 context", "chrX", "chrY", "other")
  legend_cols <- unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")])
  x0 <- 0.13
  y0 <- 10.28
  gap <- 0.175
  for (i in seq_along(legend_labels)) {
    x <- x0 + (i - 1) * gap
    rect(x, y0 - 0.070, x + 0.020, y0 + 0.070, col = legend_cols[i], border = NA)
    text(x + 0.026, y0, legend_labels[i], adj = 0, cex = 0.62)
  }
  par(xpd = FALSE)
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.1, 8.9, 3.9, 4.3), xaxs = "i", yaxs = "i")
  n_tracks <- length(event_order) * length(length_order)
  plot(
    NA,
    xlim = c(-0.52, 1.30),
    ylim = c(0.20, n_tracks + 1.55),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  title("Fig5 query-grid chopped raw FASTA SweepGA/FastGA f16 panels", cex.main = 1.02, line = 2.25)
  mtext(
    "1:1 SweepGA filter: --overlap 0, --scoring ani, --scaffold-jump 0; query-grid chunks only (10 kb, 5 kb, 2 kb)",
    side = 3,
    line = 0.55,
    cex = 0.73
  )
  draw_legend()

  row_i <- 0
  for (event_id in event_order) {
    for (length_bp in length_order) {
      row_i <- row_i + 1
      y <- n_tracks - row_i + 1
      row <- summary[as.character(summary$event_id) == event_id & summary$chop_length_bp == length_bp, ]
      if (nrow(row) == 0) {
        next
      }
      w_start <- as.numeric(row$window_start[1])
      w_end <- as.numeric(row$window_end[1])
      xscale <- function(x) (x - w_start) / (w_end - w_start)
      if (length_bp == 10000) {
        segments(-0.49, y + 0.54, 1.25, y + 0.54, col = "#D0D0D0", lwd = 0.8)
      }
      rect(0, y - 0.18, 1, y + 0.18, col = "#F3F3F3", border = "#C7C7C7", lwd = 0.7)
      rows <- segments[
        as.character(segments$event_id) == event_id &
          segments$chop_length_bp == length_bp,
      ]
      if (nrow(rows) > 0) {
        rows <- rows[order(rows$query_clip_start, rows$query_clip_end, rows$target_chrom), ]
        for (j in seq_len(nrow(rows))) {
          x1 <- max(0, min(1, xscale(as.numeric(rows$query_clip_start[j]))))
          x2 <- max(0, min(1, xscale(as.numeric(rows$query_clip_end[j]))))
          col <- target_col(rows$target_chrom[j])
          border <- if (rows$is_expected_target[j] == "yes") "#111111" else NA
          rect(
            x1,
            y - 0.24,
            x2,
            y + 0.24,
            col = col,
            border = border,
            lwd = if (rows$is_expected_target[j] == "yes") 0.45 else 0
          )
        }
      }
      label <- sprintf("%s\n%s query grid", event_label(event_id), fmt_len(length_bp))
      text(-0.040, y, label, adj = 1, cex = 0.69)
      text(
        1.025,
        y,
        sprintf(
          "expected rows %s; union %.1f kb%s",
          row$expected_target_rows[1],
          as.numeric(row$union_expected_overlap_bp[1]) / 1000,
          ifelse(
            nzchar(row$audit_chr3_query_redundant_bp[1]),
            sprintf("\nq-grid audit: %s chr3 redundant bp", row$audit_chr3_query_redundant_bp[1]),
            ""
          )
        ),
        adj = 0,
        cex = 0.50,
        col = "#333333"
      )
      ticks <- pretty(c(w_start, w_end), n = 4)
      ticks <- ticks[ticks >= w_start & ticks <= w_end]
      for (tick in ticks) {
        x <- xscale(tick)
        segments(x, y - 0.36, x, y - 0.29, col = "#444444", lwd = 0.55)
        text(x, y - 0.45, fmt_coord(tick), cex = 0.45, adj = c(0.5, 1))
      }
    }
  }
  mtext("Genomic query coordinates within each 500 kb display window", side = 1, line = 2.15, cex = 0.78)
}

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 13.5, height = 8.1, useDingbats = FALSE)
draw_panel()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 13.5, height = 8.1, onefile = TRUE)
draw_panel()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 2700, height = 1620, res = 200, type = "cairo")
draw_panel()
dev.off()
