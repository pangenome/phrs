#!/usr/bin/env Rscript
# Extended Data Figure 3 — Annotation: TAR1, internal (TTAGGG)n islands, telomere length
# Inputs: see paper_prep/figures/ed3/sources.tsv
# Output: figure_ed3.pdf, figure_ed3.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III"
OUT  <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "paper_prep/figures/ed3"

# ---------- Load data ----------
tar1_arm <- read_tsv(file.path(ROOT, "enrichment", "community_tar1_by_arm.tsv"),
                     show_col_types = FALSE)
length_dist <- read_tsv(file.path(ROOT, "ttaggg_analysis", "length_distribution.tsv"),
                        show_col_types = FALSE)
motif_comp <- read_tsv(file.path(ROOT, "ttaggg_analysis", "motif_composition.tsv"),
                       show_col_types = FALSE)
telo_comm <- read_tsv(file.path(ROOT, "telomere_length_by_community.tsv"),
                      show_col_types = FALSE, comment = "#")
tar1_pos <- read_tsv(file.path(ROOT, "enrichment", "tar1_positional",
                               "tar1_positional_per_arm.tsv"),
                     show_col_types = FALSE)

# ---------- Helpers ----------
classify_arm <- function(arm) {
  is_par <- arm %in% c("chrX_parm", "chrY_parm")
  is_acro_p <- arm %in% c("chr13_parm", "chr14_parm", "chr15_parm",
                          "chr21_parm", "chr22_parm")
  ifelse(is_par, "PAR1",
         ifelse(is_acro_p, "Acrocentric p", "Autosomal/sex"))
}

palette_arm <- c("PAR1" = "#d62728",
                 "Acrocentric p" = "#ff7f0e",
                 "Autosomal/sex" = "#1f77b4")

short_arm <- function(a) sub("_parm", "p", sub("_qarm", "q", a))

# ---------- Panel A: TAR1 prevalence per arm ----------
panel_a <- function() {
  d <- tar1_arm %>%
    mutate(class = classify_arm(chr_arm),
           short = short_arm(chr_arm)) %>%
    arrange(pct_with_tar1)

  par(mar = c(4.6, 2.7, 2.4, 0.6))
  cols <- palette_arm[d$class]
  bp <- barplot(d$pct_with_tar1, horiz = FALSE, las = 2,
                names.arg = d$short, col = cols, border = NA,
                ylab = "TAR1 prevalence (%)",
                ylim = c(0, 105),
                cex.names = 0.45, cex.axis = 0.7, cex.lab = 0.85,
                main = "")
  mtext("a   TAR1 prevalence per arm (n=41)", side = 3, line = 0.6,
        adj = 0, cex = 0.85, font = 2)
  abline(h = c(50, 100), col = "grey80", lty = 3)
  legend("bottomright", legend = names(palette_arm), fill = palette_arm,
         border = NA, cex = 0.6, bty = "n", inset = c(0.02, 0.02))
}

# ---------- Panel B: TTAGGG island length + motif composition ----------
panel_b <- function() {
  par(mar = c(4.0, 4.2, 2.4, 4.4))
  d <- length_dist
  d$bin <- factor(d$bin, levels = c("50-74","75-99","100-149","150-199",
                                    "200-299","300-499","500-1000"))
  bp <- barplot(d$pct, names.arg = as.character(d$bin), las = 2,
                col = "#4c72b0", border = NA,
                ylab = "% of islands",
                ylim = c(0, 60),
                cex.names = 0.7, cex.axis = 0.75, cex.lab = 0.85,
                main = "")
  mtext("b   (TTAGGG)n island length + motif composition",
        side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)
  text(bp, d$pct + 1.5, sprintf("%.1f%%", d$pct), cex = 0.6)
  text(par("usr")[2] * 0.98, 56, "n = 18,352 islands", adj = c(1, 1),
       cex = 0.65, col = "grey25")
  mtext("Length (bp)", side = 1, line = 2.7, cex = 0.75)

  # Motif composition stacked bar drawn in plot coords (no fig override)
  m <- motif_comp[order(-motif_comp$hexamer_fraction), ]
  m$pct <- m$hexamer_fraction * 100
  cum <- c(0, cumsum(m$pct))
  cols_m <- c("TTAGGG" = "#2ca02c", "TGAGGG" = "#9467bd",
              "TTGGGG" = "#bcbd22", "TCAGGG" = "#17becf")
  usr <- par("usr")
  xmin <- usr[1] + 0.45 * (usr[2] - usr[1])
  xmax <- usr[2] - 0.02 * (usr[2] - usr[1])
  ybot <- 38; ytop <- 47
  text(xmin, ytop + 1.6, "Hexamer composition (296,406 hexamers)",
       adj = c(0, 0), cex = 0.62, font = 2)
  for (i in seq_len(nrow(m))) {
    x0 <- xmin + (cum[i] / 100) * (xmax - xmin)
    x1 <- xmin + (cum[i + 1] / 100) * (xmax - xmin)
    rect(x0, ybot, x1, ytop, col = cols_m[m$motif[i]], border = "white", lwd = 0.6)
    if (m$pct[i] > 8) {
      text((x0 + x1) / 2, (ybot + ytop) / 2,
           sprintf("%s\n%.1f%%", m$motif[i], m$pct[i]),
           cex = 0.5, col = "white")
    }
  }
  text(xmin, ybot - 1, "0%", adj = c(0.5, 1), cex = 0.55)
  text(xmax, ybot - 1, "100%", adj = c(0.5, 1), cex = 0.55)
}

# ---------- Panel C: Terminal telomere length by community ----------
panel_c <- function() {
  d <- telo_comm %>%
    arrange(median_telo_bp) %>%
    mutate(community = factor(community, levels = community))

  par(mar = c(4.2, 4.4, 2.4, 0.8))
  ymin <- min(d$mean_telo_bp - d$sd_telo_bp) - 200
  ymax <- max(d$mean_telo_bp + d$sd_telo_bp) + 400
  plot(seq_len(nrow(d)), d$median_telo_bp, type = "n",
       xlim = c(0.5, nrow(d) + 0.5),
       ylim = c(ymin, ymax),
       xaxt = "n", xlab = "", ylab = "Telomere length (bp)",
       cex.axis = 0.75, cex.lab = 0.85)
  axis(1, at = seq_len(nrow(d)), labels = as.character(d$community),
       las = 2, cex.axis = 0.7)
  arrows(seq_len(nrow(d)), d$mean_telo_bp - d$sd_telo_bp,
         seq_len(nrow(d)), d$mean_telo_bp + d$sd_telo_bp,
         angle = 90, code = 3, length = 0.02, col = "grey55")
  points(seq_len(nrow(d)), d$mean_telo_bp, pch = 23, bg = "grey90",
         col = "grey45", cex = 0.65)
  points(seq_len(nrow(d)), d$median_telo_bp, pch = 21, bg = "#d62728",
         col = "white", cex = 1.05)

  mtext("c   Terminal telomere length by community", side = 3, line = 0.6,
        adj = 0, cex = 0.85, font = 2)
  mtext("Community (ordered by median)", side = 1, line = 3.0, cex = 0.75)
  legend("topleft",
         legend = c("median", "mean ± SD"),
         pch = c(21, 23), pt.bg = c("#d62728", "grey90"),
         col = c("white", "grey45"), pt.cex = c(1.0, 0.7),
         cex = 0.65, bty = "n", inset = c(0.02, 0.05))
  txt <- sprintf("Kruskal-Wallis H = 100.89, p = 3.2e-15\nN = %d sequences across %d communities",
                 sum(d$n_sequences), nrow(d))
  text(par("usr")[2] - 0.5, ymax - 200, txt, adj = c(1, 1),
       cex = 0.65, col = "grey25")
}

# ---------- Panel D: TAR1 positional per arm ----------
panel_d <- function() {
  d <- tar1_pos %>%
    mutate(class = classify_arm(chr_arm),
           short = short_arm(chr_arm),
           dist_kb = pmax(median_dist_from_telo_kb, 0.05)) %>%
    arrange(dist_kb)

  par(mar = c(4.6, 4.4, 2.4, 0.6))
  cols <- palette_arm[d$class]
  bp <- barplot(log10(d$dist_kb), names.arg = d$short, las = 2, horiz = FALSE,
                col = cols, border = NA,
                ylab = expression("Median TAR1 distance from telomere (log"[10]*" kb)"),
                cex.names = 0.45, cex.axis = 0.7, cex.lab = 0.85,
                ylim = c(-1.5, 2.6),
                main = "")
  mtext("d   Per-arm TAR1 distance-from-telomere",
        side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)
  axis(2, at = log10(c(0.1, 1, 10, 100, 250)),
       labels = c("0.1", "1", "10", "100", "250"),
       cex.axis = 0.7)
  abline(h = log10(c(10, 25)), col = "grey80", lty = 3)
  legend("topleft", legend = c("≤10 kb of telomere (66.9% of TAR1)", "≤25 kb (70.3%)"),
         col = "grey55", lty = 3, cex = 0.6, bty = "n", inset = c(0.02, 0.02))
}

# ---------- Render ----------
render_to <- function(dev_open, dev_close) {
  dev_open()
  layout(matrix(1:4, nrow = 2, byrow = TRUE),
         widths = c(1.0, 1.0), heights = c(1.0, 1.0))
  panel_a()
  panel_b()
  panel_c()
  panel_d()
  dev_close()
}

pdf_path <- file.path(OUT, "figure_ed3.pdf")
png_path <- file.path(OUT, "figure_ed3.png")

render_to(function() pdf(pdf_path, width = 11, height = 8.5),
          dev.off)
render_to(function() png(png_path, width = 11 * 200, height = 8.5 * 200,
                         res = 200, type = "cairo"),
          dev.off)

cat("Wrote", pdf_path, "and", png_path, "\n")
