#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels"
prefix <- "fig5_homolog_vs_interchrom_zoom_panels"

segments <- read.delim(file.path(panel_dir, "telomeric_interchrom_winner_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(file.path(panel_dir, "telomeric_interchrom_winner_arm_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
selected <- read.delim(file.path(panel_dir, "selected_telomeric_arms.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

target_cols <- c(
  chr1 = "#4E79A7", chr2 = "#A0CBE8", chr3 = "#D95F02", chr4 = "#FFBE7D",
  chr5 = "#59A14F", chr6 = "#8CD17D", chr7 = "#B6992D", chr8 = "#F1CE63",
  chr9 = "#1B9E77", chr10 = "#86BCB6", chr11 = "#E15759", chr12 = "#FF9D9A",
  chr13 = "#79706E", chr14 = "#BAB0AC", chr15 = "#9C755F", chr16 = "#D7B5A6",
  chr17 = "#B07AA1", chr18 = "#D4A6C8", chr19 = "#2F4B7C", chr20 = "#665191",
  chr21 = "#A05195", chr22 = "#D45087", chrX = "#7570B3", chrY = "#E7298A",
  other = "#8A8A8A"
)

target_col <- function(chrom) {
  chrom <- as.character(chrom)
  out <- unname(target_cols[chrom])
  out[is.na(out)] <- unname(target_cols["other"])
  out
}

arm_label <- function(chrom, arm) {
  if (chrom == "chrX" && arm == "p") return("PAR1 Xp")
  if (chrom == "chr9" && arm == "q") return("chr9q / chr3q")
  paste0(chrom, arm)
}

draw_panels <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(4.0, 9.2, 3.6, 6.4), xaxs = "i", yaxs = "i")

  n <- nrow(selected)
  plot(
    NA,
    xlim = c(-110, 680),
    ylim = c(0.12, n + 0.86),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  title("PAN027 paternal telomeric homolog-vs-interchrom winners", cex.main = 1.05, line = 2.15)
  mtext(
    "Each colored block is one 2 kb query window where best interchromosomal IMPG similarity beats the best same-chromosome/homolog match.",
    side = 3,
    line = 0.55,
    cex = 0.70
  )

  axis(
    1,
    at = seq(0, 500, by = 100),
    labels = seq(0, 500, by = 100),
    tick = FALSE,
    cex.axis = 0.72,
    line = 0.3
  )
  mtext("Distance across distal 500 kb arm window, kb", side = 1, line = 2.35, cex = 0.78)

  for (i in seq_len(n)) {
    chrom <- selected$query_chrom[i]
    arm <- selected$arm[i]
    y <- n - i + 1
    row_summary <- summary[summary$query_chrom == chrom & summary$arm == arm, ][1, ]
    row_segments <- segments[segments$query_chrom == chrom & segments$arm == arm, ]
    display_start <- as.numeric(row_segments$display_start[1])

    rect(0, y - 0.16, 500, y + 0.16, col = "#F2F2F2", border = "#BDBDBD", lwd = 0.75)
    if (nrow(row_segments) > 0) {
      row_segments <- row_segments[order(row_segments$query_start, row_segments$query_end, row_segments$target_chrom), ]
      for (j in seq_len(nrow(row_segments))) {
        x1 <- (as.numeric(row_segments$query_start[j]) - display_start) / 1000
        x2 <- (as.numeric(row_segments$query_end[j]) - display_start) / 1000
        border <- if (chrom == "chr9" && arm == "q" && row_segments$target_chrom[j] == "chr3") "#111111" else NA
        rect(
          max(0, x1),
          y - 0.24,
          min(500, x2),
          y + 0.24,
          col = target_col(row_segments$target_chrom[j]),
          border = border,
          lwd = 0.35
        )
      }
    }

    text(-6, y, arm_label(chrom, arm), adj = 1, cex = 0.86, font = 2)
    text(
      514,
      y,
      sprintf(
        "%.0f kb; n=%d; delta=%.3f",
        as.numeric(row_summary$winner_bp) / 1000,
        as.integer(row_summary$winner_windows),
        as.numeric(row_summary$mean_inter_minus_same_identity)
      ),
      adj = 0,
      cex = 0.61,
      col = "#333333"
    )
  }

  legend_targets <- unique(segments$target_chrom)
  preferred <- c("chr3", "chrY", "chr1", "chr13", "chr14", "chr15", "chr21", "chr22")
  legend_targets <- c(intersect(preferred, legend_targets), sort(setdiff(legend_targets, preferred)))
  legend(
    0,
    0.58,
    legend = legend_targets,
    fill = target_col(legend_targets),
    border = NA,
    title = "Winning target",
    bty = "n",
    horiz = TRUE,
    cex = 0.62,
    title.cex = 0.66
  )
}

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 12.6, height = 5.9, useDingbats = FALSE)
draw_panels()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 2520, height = 1180, res = 200, type = "cairo")
draw_panels()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 12.6, height = 5.9, onefile = TRUE)
draw_panels()
dev.off()
