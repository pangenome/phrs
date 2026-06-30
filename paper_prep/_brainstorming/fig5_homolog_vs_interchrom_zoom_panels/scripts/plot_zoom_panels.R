#!/usr/bin/env Rscript

out_dir <- "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels"
segments <- read.delim(file.path(out_dir, "zoom_window_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(out_dir, "zoom_panel_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

if (!"target_haplotype" %in% names(segments)) {
  segments$target_haplotype <- "NA"
}
summary <- summary[order(summary$panel_order), ]
segments$relative_start_kb <- segments$relative_start / 1000
segments$relative_end_kb <- segments$relative_end / 1000
segments$relative_mid_kb <- (segments$relative_start + segments$relative_end) / 2000

target_cols <- c(
  chrY = "#E7298A",
  chr1 = "#4E79A7",
  chr3 = "#D95F02",
  other = "#B8B8B8"
)

target_col <- function(x) {
  ifelse(x %in% names(target_cols), target_cols[x], target_cols[["other"]])
}

panel_coord_label <- function(row) {
  len <- as.numeric(row$query_length)
  zoom <- as.numeric(row$zoom_bp)
  if (row$arm == "p") {
    sprintf("%s p: %.3f-%.3f Mb", row$query_chrom, 0, zoom / 1e6)
  } else {
    sprintf("%s q: %.3f-%.3f Mb", row$query_chrom, (len - zoom) / 1e6, len / 1e6)
  }
}

row_coord_endpoints <- function(row) {
  len <- as.numeric(row$query_length)
  zoom <- as.numeric(row$zoom_bp)
  if (row$arm == "p") {
    return(c(sprintf("%.3f Mb", 0), sprintf("%.3f Mb", zoom / 1e6)))
  }
  c(sprintf("%.3f Mb", (len - zoom) / 1e6), sprintf("%.3f Mb", len / 1e6))
}

draw_break_glyph <- function(x, y, direction = c("right", "left")) {
  direction <- match.arg(direction)
  slant <- if (direction == "right") 1 else -1
  xpd_old <- par("xpd")
  par(xpd = NA)
  for (offset in c(-2.4, 2.4)) {
    segments(
      x + offset - slant * 2.2, y - 0.16,
      x + offset + slant * 2.2, y + 0.16,
      col = "#555555",
      lwd = 1.3,
      lend = "butt"
    )
  }
  par(xpd = xpd_old)
}

draw_coord_endpoints <- function(x0, x1, y, labels, cex = 0.42) {
  text(x0, y, labels[1], cex = cex, col = "#555555", adj = c(0, 1), xpd = NA)
  text(x1, y, labels[2], cex = cex, col = "#555555", adj = c(1, 1), xpd = NA)
}

draw_scale_bar <- function(x0, y, width_kb = 100) {
  x1 <- x0 + width_kb
  segments(x0, y, x1, y, col = "#444444", lwd = 1.1, lend = "butt")
  segments(c(x0, x1), y - 0.055, c(x0, x1), y + 0.055, col = "#444444", lwd = 1.1)
  text((x0 + x1) / 2, y - 0.10, sprintf("%d kb", width_kb), cex = 0.58, col = "#333333", adj = c(0.5, 1), xpd = NA)
}

draw_phr_span <- function(x0, x1, y) {
  y0 <- y + 0.145
  segments(x0, y0, x1, y0, col = "#555555", lwd = 0.85, lend = "butt")
  segments(c(x0, x1), y0 - 0.045, c(x0, x1), y0 + 0.045, col = "#555555", lwd = 0.85, lend = "butt")
}

phr_interval_kb <- function(row) {
  if (!all(c("phr_full_start", "phr_full_end") %in% names(row))) {
    return(NULL)
  }
  full_start <- suppressWarnings(as.numeric(row$phr_full_start))
  full_end <- suppressWarnings(as.numeric(row$phr_full_end))
  if (is.na(full_start) || is.na(full_end) || full_end <= full_start) {
    return(NULL)
  }
  len <- as.numeric(row$query_length)
  zoom <- as.numeric(row$zoom_bp)
  if (row$arm == "p") {
    x0 <- full_start
    x1 <- full_end
  } else {
    window_start <- len - zoom
    x0 <- full_start - window_start
    x1 <- full_end - window_start
  }
  x0 <- max(0, x0)
  x1 <- min(zoom, x1)
  if (x1 <= x0) {
    return(NULL)
  }
  c(x0 / 1000, x1 / 1000)
}

draw_zoom <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(2.4, 6.6, 3.3, 4.2), xaxs = "i", yaxs = "i")
  track_x0 <- 0
  track_x1 <- 500
  n <- nrow(summary)
  plot_rows <- summary
  row_step <- 0.66
  plot_rows$plot_y <- 0.74 + rev(seq_len(n) - 1) * row_step
  header_y <- max(plot_rows$plot_y) + 0.42
  legend_y <- header_y + 0.28
  x_center <- (track_x0 + track_x1) / 2
  x_half_span <- 380
  x_center_offset <- -28

  plot(
    NA,
    xlim = x_center + x_center_offset + c(-x_half_span, x_half_span),
    ylim = c(0.03, legend_y + 0.23),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  title("PAN027 paternal hap2 (PAN027#2) vs father PAN011 joint haplotypes (PAN011#1 + PAN011#2)", line = 2.15, cex.main = 0.98)
  mtext(
    "Filled 2 kb windows show where the best interchromosomal IMPG match beats the same-chromosome/homolog match; thin brackets mark population-derived PHR intervals",
    side = 3,
    line = 0.55,
    cex = 0.60
  )
  text(
    (track_x0 + track_x1) / 2,
    header_y,
    "500 kb subtelomeric windows (chromosome coordinates; p telomere left, q telomere right)",
    cex = 0.62,
    font = 2
  )

  p_ticks <- seq(0, 500, by = 100)
  q_ticks <- seq(0, 500, by = 100)
  label_x <- -38

  for (i in seq_len(nrow(plot_rows))) {
    row <- plot_rows[i, ]
    y <- row$plot_y
    is_p <- row$arm == "p"
    x0 <- track_x0
    x1 <- track_x1

    rect(x0, y - 0.075, x1, y + 0.075, col = "#F5F5F5", border = "#C7C7C7", lwd = 0.75)
    if (is_p) {
      segments(track_x0 + p_ticks, y - 0.075, track_x0 + p_ticks, y + 0.075, col = "#EFEFEF", lwd = 0.5)
      draw_break_glyph(track_x1 + 15, y, "right")
      text(label_x, y + 0.065, row$panel_label, adj = 1, cex = 0.66, font = 2, xpd = NA)
      text(label_x, y - 0.115, panel_coord_label(row), adj = 1, cex = 0.49, col = "#444444", xpd = NA)
    } else {
      segments(track_x0 + q_ticks, y - 0.075, track_x0 + q_ticks, y + 0.075, col = "#EFEFEF", lwd = 0.5)
      draw_break_glyph(track_x0 - 15, y, "left")
      text(label_x, y + 0.065, row$panel_label, adj = 1, cex = 0.66, font = 2, xpd = NA)
      text(label_x, y - 0.115, panel_coord_label(row), adj = 1, cex = 0.49, col = "#444444", xpd = NA)
    }
    draw_coord_endpoints(x0, x1, y - 0.19, row_coord_endpoints(row), cex = 0.38)

    rows <- segments[segments$panel_id == row$panel_id, ]
    if (nrow(rows) > 0) {
      if (is_p) {
        rows$plot_start_kb <- x0 + rows$relative_start_kb
        rows$plot_end_kb <- x0 + rows$relative_end_kb
      } else {
        rows$plot_start_kb <- x0 + as.numeric(row$zoom_bp) / 1000 - rows$relative_end_kb
        rows$plot_end_kb <- x0 + as.numeric(row$zoom_bp) / 1000 - rows$relative_start_kb
      }
      rows <- rows[order(rows$plot_start_kb, rows$plot_end_kb, rows$target_chrom, rows$target_haplotype), ]

      for (j in seq_len(nrow(rows))) {
        rect(
          rows$plot_start_kb[j],
          y - 0.105,
          rows$plot_end_kb[j],
          y + 0.105,
          col = target_col(rows$target_bucket[j]),
          border = NA
        )
      }
    }
    phr <- phr_interval_kb(row)
    if (!is.null(phr)) {
      draw_phr_span(phr[1], phr[2], y)
    }
  }

  legend_targets <- unique(segments$target_bucket)
  legend_targets <- legend_targets[order(match(legend_targets, c("chrY", "chr1", "chr3", "other")))]
  legend_targets <- legend_targets[!is.na(legend_targets)]
  legend_labels <- c(chrY = "target chrY", chr1 = "target chr1", chr3 = "target chr3", other = "other target")
  if (length(legend_targets) > 0) {
    legend_step <- 82
    x0 <- ((track_x0 + track_x1) - (length(legend_targets) - 1) * legend_step) / 2
    y0 <- legend_y
    for (k in seq_along(legend_targets)) {
      x <- x0 + (k - 1) * legend_step
      rect(x, y0 - 0.075, x + 5, y0 + 0.075, col = target_col(legend_targets[k]), border = NA)
      text(x + 6.5, y0, legend_labels[legend_targets[k]], adj = 0, cex = 0.55)
    }
  }
  draw_scale_bar((track_x0 + track_x1 - 100) / 2, 0.19, 100)
}

pdf(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.pdf"), width = 13.8, height = 4.6, useDingbats = FALSE)
draw_zoom()
dev.off()

png(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.png"), width = 2760, height = 920, res = 200, type = "cairo")
draw_zoom()
dev.off()

svg(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.svg"), width = 13.8, height = 4.6, onefile = TRUE)
draw_zoom()
dev.off()
