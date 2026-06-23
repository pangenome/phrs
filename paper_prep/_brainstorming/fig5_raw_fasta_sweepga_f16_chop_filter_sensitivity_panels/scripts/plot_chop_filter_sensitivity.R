#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
segments <- read.delim(file.path(panel_dir, "chop_filter_panel_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(panel_dir, "chop_filter_panel_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(segments) == 0) {
  stop("no segments to plot")
}

event_levels <- c("PAN027_chr9q_chr3q_PHR_candidate", "PAN028_chr9q_chr3q_PHR_candidate")
mode_levels <- c("no_merge_ani", "no_merge_log_length_ani", "scaffold50k_ani", "scaffold50k_log_length_ani")
chop_levels <- c(10000, 5000, 2000)
summary <- summary[summary$event_id %in% event_levels, ]
summary$event_id <- factor(summary$event_id, levels = event_levels)
summary$filter_mode <- factor(summary$filter_mode, levels = mode_levels)
summary <- summary[order(summary$event_id, summary$chop_length_bp, summary$filter_mode), ]

target_palette <- c(
  chr3 = "#C43B2B",
  chr9 = "#8FA2B2",
  chr16 = "#5B8E5A",
  chr19 = "#B5833A",
  chr1 = "#8E6AAE",
  other = "#B8B8B8"
)

target_col <- function(chrom) {
  ifelse(chrom %in% names(target_palette), target_palette[chrom], target_palette[["other"]])
}

event_label <- function(event_id) {
  switch(
    as.character(event_id),
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027 paternal chr9q candidate",
    PAN028_chr9q_chr3q_PHR_candidate = "PAN028 maternal chr9q candidate",
    as.character(event_id)
  )
}

fmt_bp <- function(x) {
  x <- as.integer(x)
  ifelse(x >= 1000, paste0(round(x / 1000, 1), " kb"), paste0(x, " bp"))
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.8, 7.2, 3.0, 1.8), oma = c(0.5, 0.5, 2.0, 0.5), xaxs = "i", yaxs = "i")
  layout(matrix(seq_along(event_levels), ncol = 1), heights = c(1, 1))
  for (event_id in event_levels) {
    rows_event <- summary[as.character(summary$event_id) == event_id, ]
    x0 <- min(rows_event$window_start_abs)
    x1 <- max(rows_event$window_end_abs)
    y_max <- length(mode_levels) * length(chop_levels)
    plot(
      NA,
      xlim = c(x0, x1),
      ylim = c(0.25, y_max + 0.85),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    axis_at <- pretty(c(x0, x1), n = 5)
    axis(1, at = axis_at, labels = sprintf("%.2f", axis_at / 1e6), cex.axis = 0.72)
    mtext("Absolute query chromosome coordinate (Mb)", side = 1, line = 2.3, cex = 0.78)
    title(event_label(event_id), cex.main = 0.9, font.main = 2)
    rect(x0, 0.25, x1, y_max + 0.85, col = "#FFFFFF", border = NA)
    abline(v = axis_at, col = "#EFEFEF", lwd = 0.6)
    y <- y_max
    y_labels <- character()
    y_pos <- numeric()
    for (chop in chop_levels) {
      for (mode in mode_levels) {
        s <- rows_event[rows_event$chop_length_bp == chop & as.character(rows_event$filter_mode) == mode, ]
        if (nrow(s) != 1) {
          next
        }
        rect(x0, y - 0.22, x1, y + 0.22, col = "#F6F6F6", border = "#D1D1D1", lwd = 0.5)
        seg <- segments[
          segments$event_id == event_id &
            segments$chop_length_bp == chop &
            segments$filter_mode == mode,
        ]
        if (nrow(seg) > 0) {
          seg <- seg[order(seg$query_start_abs, seg$query_end_abs), ]
          for (i in seq_len(nrow(seg))) {
            col <- target_col(seg$target_chrom[i])
            border <- if (seg$target_chrom[i] == "chr3") "#111111" else NA
            rect(
              seg$query_start_abs[i],
              y - 0.28,
              seg$query_end_abs[i],
              y + 0.28,
              col = col,
              border = border,
              lwd = if (seg$target_chrom[i] == "chr3") 0.45 else 0
            )
          }
        }
        text(
          x1 + (x1 - x0) * 0.012,
          y,
          paste0("chr3 ", fmt_bp(s$chr3_query_union_bp), " (", s$chr3_rows, " rows)"),
          adj = 0,
          cex = 0.52,
          xpd = NA,
          col = target_palette[["chr3"]]
        )
        y_labels <- c(y_labels, paste0(chop / 1000, " kb  ", s$filter_label))
        y_pos <- c(y_pos, y)
        y <- y - 1
      }
    }
    axis(2, at = y_pos, labels = y_labels, las = 1, tick = FALSE, cex.axis = 0.52)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  title("Fig5 raw FASTA SweepGA/FastGA f16 chop/filter sensitivity", outer = TRUE, cex.main = 0.95, font.main = 2)
  legend(
    "bottom",
    horiz = TRUE,
    bty = "n",
    inset = 0.015,
    cex = 0.65,
    legend = c("chr3 donor", "chr9 context", "chr16", "chr19", "chr1", "other"),
    fill = unname(target_palette[c("chr3", "chr9", "chr16", "chr19", "chr1", "other")]),
    border = NA
  )
}

pdf(
  file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.pdf"),
  width = 11.2,
  height = 7.2,
  useDingbats = FALSE
)
draw_panel()
dev.off()

svg(
  file.path(panel_dir, "fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.svg"),
  width = 11.2,
  height = 7.2,
  pointsize = 12
)
draw_panel()
dev.off()
