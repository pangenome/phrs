#!/usr/bin/env Rscript

out_dir <- "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels"
segments <- read.delim(file.path(out_dir, "zoom_window_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(out_dir, "zoom_panel_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

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

panel_coord_label <- function(row) {
  len <- as.numeric(row$query_length)
  zoom <- as.numeric(row$zoom_bp)
  if (row$arm == "p") {
    sprintf("%s p: 0-%.0f kb", row$query_chrom, zoom / 1000)
  } else {
    sprintf("%s q: %.3f-%.3f Mb", row$query_chrom, (len - zoom) / 1e6, len / 1e6)
  }
}

draw_zoom <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(4.0, 7.2, 3.7, 5.1), xaxs = "i", yaxs = "i")
  n <- nrow(summary)
  plot(NA, xlim = c(-92, 660), ylim = c(0.35, n + 1.0), axes = FALSE, xlab = "", ylab = "")
  title("Subtelomeric homolog-vs-interchrom zooms", line = 2.15, cex.main = 1.08)
  mtext("Filled 2 kb windows: best interchromosomal IMPG similarity beats best same-chromosome/homolog match; x-axis is distance inward from the telomere", side = 3, line = 0.55, cex = 0.66)
  axis(1, at = seq(0, 500, by = 100), labels = paste0(seq(0, 500, by = 100), " kb"), cex.axis = 0.67)
  abline(v = seq(0, 500, by = 100), col = "#EFEFEF", lwd = 0.55)

  for (i in seq_len(n)) {
    row <- summary[i, ]
    y <- n - i + 1
    rect(0, y - 0.18, 500, y + 0.18, col = "#F5F5F5", border = "#C7C7C7", lwd = 0.8)
    rows <- segments[segments$panel_id == row$panel_id, ]
    if (nrow(rows) > 0) {
      rows <- rows[order(rows$relative_start, rows$relative_end, rows$target_chrom), ]
      for (j in seq_len(nrow(rows))) {
        rect(
          rows$relative_start_kb[j],
          y - 0.30,
          rows$relative_end_kb[j],
          y + 0.30,
          col = target_col(rows$target_chrom[j]),
          border = NA
        )
      }

      # Label contiguous target runs of at least 6 kb, plus every chr3/chrY run.
      run_start <- rows$relative_start_kb[1]
      run_end <- rows$relative_end_kb[1]
      run_target <- rows$target_chrom[1]
      emit_run <- function(start_kb, end_kb, target, y) {
        width <- end_kb - start_kb
                  if (width >= 12 || (target %in% c("chr3", "chrY") && width >= 10)) {
          text((start_kb + end_kb) / 2, y + 0.43, sprintf("%s %.0f kb", target, width), cex = 0.42, col = "#222222")
        }
      }
      if (nrow(rows) > 1) {
        for (j in 2:nrow(rows)) {
          contiguous <- rows$target_chrom[j] == run_target && rows$relative_start_kb[j] <= run_end + 0.001
          if (contiguous) {
            run_end <- max(run_end, rows$relative_end_kb[j])
          } else {
            emit_run(run_start, run_end, run_target, y)
            run_start <- rows$relative_start_kb[j]
            run_end <- rows$relative_end_kb[j]
            run_target <- rows$target_chrom[j]
          }
        }
      }
      emit_run(run_start, run_end, run_target, y)
    }
    text(-8, y + 0.11, row$panel_label, adj = 1, cex = 0.69, font = 2, xpd = NA)
    text(-8, y - 0.18, panel_coord_label(row), adj = 1, cex = 0.52, col = "#444444", xpd = NA)
    text(
      506,
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
    x0 <- 0
    y0 <- n + 0.58
    for (k in seq_along(legend_targets)) {
      x <- x0 + (k - 1) * 45
      rect(x, y0 - 0.075, x + 5, y0 + 0.075, col = target_col(legend_targets[k]), border = NA)
      text(x + 6.5, y0, legend_targets[k], adj = 0, cex = 0.55)
    }
  }
  mtext("Distance from telomere into chromosome arm", side = 1, line = 2.7, cex = 0.78)
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
