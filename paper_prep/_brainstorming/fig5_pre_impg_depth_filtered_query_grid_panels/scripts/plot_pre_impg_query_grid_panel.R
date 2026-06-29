#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels"
prefix <- "fig5_pre_impg_depth_filtered_query_grid_panels"
segments <- read.delim(file.path(panel_dir, "pre_impg_query_grid_panel_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(panel_dir, "pre_impg_query_grid_panel_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

event_order <- c("PAR1_XY_positive_control", "PAN027_chr9q_chr3q_PHR_candidate")
basis_order <- c("1:1", "4:4", "10:10")

event_label <- function(event_id) {
  switch(
    as.character(event_id),
    PAR1_XY_positive_control = "PAR1 X/Y control",
    PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
    as.character(event_id)
  )
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

fmt_coord <- function(x) {
  ifelse(x >= 1e6, sprintf("%.3f Mb", x / 1e6), sprintf("%d kb", round(x / 1000)))
}

fmt_bp <- function(x) {
  x <- as.numeric(x)
  ifelse(x >= 1000, sprintf("%.1f kb", x / 1000), sprintf("%d bp", round(x)))
}

summarize_target_union <- function(value) {
  if (is.na(value) || !nzchar(value)) return("")
  parts <- strsplit(value, ";", fixed = TRUE)[[1]]
  parsed <- do.call(rbind, lapply(parts, function(part) {
    kv <- strsplit(part, ":", fixed = TRUE)[[1]]
    if (length(kv) != 2) return(NULL)
    data.frame(chrom = kv[1], bp = as.numeric(kv[2]), stringsAsFactors = FALSE)
  }))
  if (is.null(parsed) || nrow(parsed) == 0) return("")
  parsed <- parsed[order(-parsed$bp, parsed$chrom), ]
  parsed <- parsed[seq_len(min(2, nrow(parsed))), ]
  paste(paste0(parsed$chrom, " ", fmt_bp(parsed$bp)), collapse = ", ")
}

draw_legend <- function() {
  par(xpd = NA)
  labels <- c("chr3 donor", "chr9 context", "chrX", "chrY", "other")
  cols <- unname(target_palette[c("chr3", "chr9", "chrX", "chrY", "other")])
  x0 <- 0.10
  y0 <- length(event_order) * length(basis_order) + 1.10
  gap <- 0.175
  for (i in seq_along(labels)) {
    x <- x0 + (i - 1) * gap
    rect(x, y0 - 0.070, x + 0.020, y0 + 0.070, col = cols[i], border = NA)
    text(x + 0.026, y0, labels[i], adj = 0, cex = 0.64)
  }
  par(xpd = FALSE)
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.1, 8.9, 3.9, 4.7), xaxs = "i", yaxs = "i")
  n_tracks <- length(event_order) * length(basis_order)
  plot(NA, xlim = c(-0.52, 1.55), ylim = c(0.20, n_tracks + 1.45), axes = FALSE, xlab = "", ylab = "")
  title("Fig5 pre-IMPG depth-filtered SweepGA f32 query-grid panels", cex.main = 1.02, line = 2.25)
  mtext(
    "IMPG best-all class winners over query-space 2 kb windows; outlined blocks are windows where interchromosomal similarity beats the best homolog",
    side = 3,
    line = 0.55,
    cex = 0.70
  )
  draw_legend()

  row_i <- 0
  for (event_id in event_order) {
    for (basis in basis_order) {
      row_i <- row_i + 1
      y <- n_tracks - row_i + 1
      row <- summary[summary$event_id == event_id & summary$basis == basis, ]
      if (nrow(row) == 0) next
      w_start <- as.numeric(row$window_start[1])
      w_end <- as.numeric(row$window_end[1])
      xscale <- function(x) (x - w_start) / (w_end - w_start)
      if (basis == basis_order[1]) {
        segments(-0.49, y + 0.54, 1.27, y + 0.54, col = "#D0D0D0", lwd = 0.8)
      }
      rect(0, y - 0.18, 1, y + 0.18, col = "#F3F3F3", border = "#C7C7C7", lwd = 0.7)
      rows <- segments[segments$event_id == event_id & segments$basis == basis, ]
      if (nrow(rows) > 0) {
        rows <- rows[order(rows$query_clip_start, rows$query_clip_end, rows$target_bucket), ]
        for (j in seq_len(nrow(rows))) {
          x1 <- max(0, min(1, xscale(as.numeric(rows$query_clip_start[j]))))
          x2 <- max(0, min(1, xscale(as.numeric(rows$query_clip_end[j]))))
          col <- target_col(rows$target_bucket[j])
          border <- if (rows$interchrom_beats_homolog[j] == "yes") "#111111" else NA
          rect(x1, y - 0.24, x2, y + 0.24, col = adjustcolor(col, alpha.f = if (rows$interchrom_beats_homolog[j] == "yes") 0.95 else 0.28), border = border, lwd = if (rows$interchrom_beats_homolog[j] == "yes") 0.55 else 0)
        }
      }
      label <- sprintf("%s\n%s", event_label(event_id), basis)
      text(-0.040, y, label, adj = 1, cex = 0.69)
      text(
        1.025,
        y,
        sprintf(
          "interchrom wins %s; union %.1f kb\ntop: %s",
          row$interchrom_win_rows[1],
          as.numeric(row$union_interchrom_win_bp[1]) / 1000,
          summarize_target_union(row$target_union_overlap_bp[1])
        ),
        adj = 0,
        cex = 0.47,
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

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 13.5, height = 6.6, useDingbats = FALSE)
draw_panel()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 13.5, height = 6.6, onefile = TRUE)
draw_panel()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 2700, height = 1320, res = 200, type = "cairo")
draw_panel()
dev.off()
