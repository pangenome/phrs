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
  chr1 = "#4E79A7", chr2 = "#A0CBE8", chr3 = "#D95F02", chr4 = "#FFBE7D",
  chr5 = "#59A14F", chr6 = "#8CD17D", chr7 = "#B6992D", chr8 = "#F1CE63",
  chr9 = "#1B9E77", chr10 = "#86BCB6", chr11 = "#E15759", chr12 = "#FF9D9A",
  chr13 = "#79706E", chr14 = "#8A8A8A", chr15 = "#9C755F", chr16 = "#D7B5A6",
  chr17 = "#B07AA1", chr18 = "#D4A6C8", chr19 = "#2F4B7C", chr20 = "#665191",
  chr21 = "#A05195", chr22 = "#D45087", chrX = "#7570B3", chrY = "#E7298A",
  other = "#2C7FB8"
)

target_col <- function(x) {
  ifelse(x %in% names(target_cols), target_cols[x], target_cols[["other"]])
}

fmt_targets <- function(value) {
  if (is.na(value) || !nzchar(value)) return("")
  parts <- strsplit(value, ";", fixed = TRUE)[[1]]
  pieces <- vapply(parts, function(part) {
    kv <- strsplit(part, ":", fixed = TRUE)[[1]]
    if (length(kv) != 2) return("")
    sprintf("%s %.0f kb", kv[1], as.numeric(kv[2]) / 1000)
  }, character(1))
  paste(pieces[seq_len(min(3, length(pieces)))], collapse = "\n")
}

target_run_label <- function(target_chrom, target_haplotype) {
  if (is.na(target_haplotype) || target_haplotype == "NA" || !nzchar(target_haplotype)) {
    return(target_chrom)
  }
  sprintf("%s %s", target_chrom, target_haplotype)
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
  for (offset in c(-3.2, 3.2)) {
    segments(
      x + offset - slant * 2.8, y - 0.32,
      x + offset + slant * 2.8, y + 0.32,
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

draw_zoom <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.0, 6.6, 3.7, 5.4), xaxs = "i", yaxs = "i")
  n <- nrow(summary)
  p_x0 <- 0
  p_x1 <- 500
  q_x0 <- 650
  q_x1 <- 1150

  plot(NA, xlim = c(-88, 1285), ylim = c(0.04, n + 1.12), axes = FALSE, xlab = "", ylab = "")
  title("PAN027 paternal hap2 subtelomeric homolog-vs-interchrom zooms", line = 2.15, cex.main = 1.02)
  mtext(
    "Query PAN027#2 paternal haplotype vs PAN011 father joint haplotypes; filled 2 kb windows: best interchromosomal IMPG similarity beats best same-chromosome/homolog match",
    side = 3,
    line = 0.55,
    cex = 0.62
  )
  text((p_x0 + p_x1) / 2, n + 0.60, "p-arm telomeric windows (chromosome coordinates)", cex = 0.66, font = 2)
  text((q_x0 + q_x1) / 2, n + 0.60, "q-arm telomeric windows (chromosome coordinates)", cex = 0.66, font = 2)

  p_ticks <- seq(0, 500, by = 100)
  q_ticks <- seq(0, 500, by = 100)

  for (i in seq_len(n)) {
    row <- summary[i, ]
    y <- n - i + 1
    is_p <- row$arm == "p"
    x0 <- if (is_p) p_x0 else q_x0
    x1 <- if (is_p) p_x1 else q_x1

    rect(x0, y - 0.18, x1, y + 0.18, col = "#F5F5F5", border = "#C7C7C7", lwd = 0.8)
    if (is_p) {
      segments(p_ticks, y - 0.18, p_ticks, y + 0.18, col = "#EFEFEF", lwd = 0.55)
      draw_break_glyph(p_x1 + 15, y, "right")
      text(-8, y + 0.11, row$panel_label, adj = 1, cex = 0.69, font = 2, xpd = NA)
      text(-8, y - 0.18, panel_coord_label(row), adj = 1, cex = 0.52, col = "#444444", xpd = NA)
    } else {
      segments(q_x0 + q_ticks, y - 0.18, q_x0 + q_ticks, y + 0.18, col = "#EFEFEF", lwd = 0.55)
      draw_break_glyph(q_x0 - 15, y, "left")
      text(q_x0 - 26, y + 0.11, row$panel_label, adj = 1, cex = 0.69, font = 2, xpd = NA)
      text(q_x0 - 26, y - 0.18, panel_coord_label(row), adj = 1, cex = 0.52, col = "#444444", xpd = NA)
    }
    draw_coord_endpoints(x0, x1, y - 0.43, row_coord_endpoints(row))

    rows <- segments[segments$panel_id == row$panel_id, ]
    if (nrow(rows) > 0) {
      if (is_p) {
        rows$plot_start_kb <- x0 + rows$relative_start_kb
        rows$plot_end_kb <- x0 + rows$relative_end_kb
      } else {
        rows$plot_start_kb <- q_x0 + as.numeric(row$zoom_bp) / 1000 - rows$relative_end_kb
        rows$plot_end_kb <- q_x0 + as.numeric(row$zoom_bp) / 1000 - rows$relative_start_kb
      }
      rows <- rows[order(rows$plot_start_kb, rows$plot_end_kb, rows$target_chrom, rows$target_haplotype), ]

      for (j in seq_len(nrow(rows))) {
        rect(
          rows$plot_start_kb[j],
          y - 0.30,
          rows$plot_end_kb[j],
          y + 0.30,
          col = target_col(rows$target_chrom[j]),
          border = NA
        )
      }

      # Label contiguous target runs of at least 6 kb, plus every chr3/chrY run.
      run_start <- rows$plot_start_kb[1]
      run_end <- rows$plot_end_kb[1]
      run_target <- rows$target_chrom[1]
      run_haplotype <- rows$target_haplotype[1]
      emit_run <- function(start_kb, end_kb, target, haplotype, y) {
        width <- end_kb - start_kb
        if (width >= 12 || (target %in% c("chr3", "chrY") && width >= 10)) {
          label <- target_run_label(target, haplotype)
          text((start_kb + end_kb) / 2, y + 0.43, sprintf("%s %.0f kb", label, width), cex = 0.40, col = "#222222")
        }
      }
      if (nrow(rows) > 1) {
        for (j in 2:nrow(rows)) {
          contiguous <- rows$target_chrom[j] == run_target &&
            rows$target_haplotype[j] == run_haplotype &&
            rows$plot_start_kb[j] <= run_end + 0.001
          if (contiguous) {
            run_end <- max(run_end, rows$plot_end_kb[j])
          } else {
            emit_run(run_start, run_end, run_target, run_haplotype, y)
            run_start <- rows$plot_start_kb[j]
            run_end <- rows$plot_end_kb[j]
            run_target <- rows$target_chrom[j]
            run_haplotype <- rows$target_haplotype[j]
          }
        }
      }
      emit_run(run_start, run_end, run_target, run_haplotype, y)
    }

    text(
      if (is_p) p_x1 + 32 else q_x1 + 8,
      y,
      sprintf("inter>same %.0f kb\n%s", as.numeric(row$winning_bp) / 1000, fmt_targets(row$top_targets)),
      adj = 0,
      cex = 0.46,
      col = "#333333"
    )
  }

  legend_targets <- unique(segments$target_chrom)
  legend_targets <- legend_targets[order(match(legend_targets, c("chrY", "chr3", "chr1", "chr13", "chr14", "chr15", "chr21", "chr22")), legend_targets)]
  legend_targets <- legend_targets[!is.na(legend_targets)]
  legend_targets <- legend_targets[seq_len(min(length(legend_targets), 10))]
  if (length(legend_targets) > 0) {
    x0 <- 235
    y0 <- n + 0.88
    for (k in seq_along(legend_targets)) {
      x <- x0 + (k - 1) * 45
      rect(x, y0 - 0.075, x + 5, y0 + 0.075, col = target_col(legend_targets[k]), border = NA)
      text(x + 6.5, y0, legend_targets[k], adj = 0, cex = 0.55)
    }
  }
  draw_scale_bar((p_x1 + q_x0 - 100) / 2, 0.23, 100)
}

pdf(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.pdf"), width = 13.8, height = 6.8, useDingbats = FALSE)
draw_zoom()
dev.off()

png(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.png"), width = 2760, height = 1360, res = 200, type = "cairo")
draw_zoom()
dev.off()

svg(file.path(out_dir, "fig5_homolog_vs_interchrom_zoom_panels.svg"), width = 13.8, height = 6.8, onefile = TRUE)
draw_zoom()
dev.off()
