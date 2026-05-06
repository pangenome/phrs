#!/usr/bin/env Rscript
# BoG v2 slide 14 — gene biology aside (DUX4, OR4F, TAR1)
# Three-panel summary, 1 row x 3 cols, designed to render at slide size (16:9 ~ wide).
#
# Panel A (DUX4): per-arm DUX4L copy distribution. C1 (chr4_q/chr10_q) carries the
#   D4Z4 macrosatellite (median 22 copies); other arms hit only 0-2 DUX4 copies.
# Panel B (OR4F): per-arm pseudogene fraction sweep, 11.1% (chr7_p) -> 99.8% (chr15_q).
# Panel C (TAR1): per-arm TAR1 prevalence, highlighting PAR1 (chrX_p / chrY_p ~0.5%).
#
# Inputs:
#   /moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv
#   /moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv
#   /moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv
# Output: slides/v2/slide_14_gene_biology.pdf and .png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PLOTS_ROOT  <- "/moosefs/guarracino/HPRCv2/PHR_III/plots"
ENRICH_ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III/enrichment"
OUT <- dirname(sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "slides/v2"

dux4 <- read_tsv(file.path(PLOTS_ROOT, "d4z4_dux4l_by_community.tsv"),
                 show_col_types = FALSE)
or4f <- read_csv(file.path(PLOTS_ROOT, "or4f_pseudogene_fraction.csv"),
                 show_col_types = FALSE)
tar1 <- read_tsv(file.path(ENRICH_ROOT, "community_tar1_by_arm.tsv"),
                 show_col_types = FALSE)

short_arm <- function(a) sub("_parm", "p", sub("_qarm", "q", a))

# ---------- Panel A: DUX4 copy distribution by arm ----------
panel_a <- function() {
  d <- dux4 %>%
    mutate(arm_short = sub("_q$", "q", sub("_p$", "p", ChromArm)),
           is_d4z4 = arm_short %in% c("chr4q", "chr10q"))
  arms_order <- d %>%
    group_by(arm_short) %>%
    summarise(med = median(n_dux4l), .groups = "drop") %>%
    arrange(med) %>%
    pull(arm_short)
  d$arm_short <- factor(d$arm_short, levels = arms_order)

  par(mar = c(4.6, 4.4, 3.0, 1.0))
  cols <- ifelse(levels(d$arm_short) %in% c("chr4q", "chr10q"),
                 "#d62728", "#7f7f7f")
  boxplot(n_dux4l ~ arm_short, data = d,
          horizontal = FALSE, las = 2,
          col = cols, border = "grey30", outcol = "grey50",
          ylab = "DUX4L copies per haplotype (C1 only)",
          xlab = "",
          ylim = c(0, max(d$n_dux4l) * 1.18),
          cex.axis = 0.85, cex.lab = 0.85, outpch = 16, outcex = 0.5,
          main = "")
  mtext("a   DUX4: D4Z4 array lives only at chr4q / chr10q",
        side = 3, line = 1.2, adj = 0, cex = 0.95, font = 2)
  abline(h = 22, col = "#d62728", lty = 2, lwd = 1.2)
  text(0.6, 24, "C1 median = 22", adj = c(0, 0),
       cex = 0.75, col = "#d62728")
  text(0.6, par("usr")[4] * 0.97,
       "DUX4 annotated on 18 q-arms,\nbut only chr4q/chr10q (C1) carry\nthe full D4Z4 macrosatellite array.\nOther arms: 0-2 copies only.",
       adj = c(0, 1), cex = 0.78, col = "grey25")
}

# ---------- Panel B: OR4F pseudogenisation gradient ----------
panel_b <- function() {
  d <- or4f %>%
    arrange(pseudo_frac) %>%
    mutate(short = short_arm(chr_arm),
           pct = pseudo_frac * 100,
           is_extreme = chr_arm %in% c("chr7_parm", "chr15_qarm"))

  par(mar = c(4.6, 4.4, 3.0, 1.0))
  cols <- ifelse(d$is_extreme, "#d62728", "#7f7f7f")
  bp <- barplot(d$pct, names.arg = d$short, horiz = FALSE, las = 2,
                col = cols, border = NA,
                ylab = "OR4F pseudogene fraction (%)",
                ylim = c(0, 110),
                cex.names = 0.7, cex.axis = 0.75, cex.lab = 0.85)
  mtext("b   OR4F: gradient of decay across arms",
        side = 3, line = 1.2, adj = 0, cex = 0.95, font = 2)
  abline(h = 62.1, col = "#d62728", lty = 2, lwd = 1.2)
  text(par("usr")[2], 65, " mean = 62.1%", adj = c(1, 0),
       cex = 0.7, col = "#d62728")
  for (i in seq_along(d$is_extreme)) {
    if (d$is_extreme[i]) {
      text(bp[i], d$pct[i] + 5, sprintf("%.1f%%", d$pct[i]),
           cex = 0.75, col = "#d62728", font = 2)
    }
  }
  text(par("usr")[2] * 0.98, 100,
       sprintf("n = %d arms\nN = %d OR4F entries", nrow(d), sum(d$total)),
       adj = c(1, 1), cex = 0.65, col = "grey25")
}

# ---------- Panel C: TAR1 prevalence per arm ----------
panel_c <- function() {
  d <- tar1 %>%
    arrange(pct_with_tar1) %>%
    mutate(short = short_arm(chr_arm),
           is_par1 = chr_arm %in% c("chrX_parm", "chrY_parm"))

  par(mar = c(4.6, 4.4, 3.0, 1.0))
  cols <- ifelse(d$is_par1, "#1f77b4", "#7f7f7f")
  bp <- barplot(d$pct_with_tar1, names.arg = d$short, horiz = FALSE, las = 2,
                col = cols, border = NA,
                ylab = "Sequences carrying TAR1 (%)",
                ylim = c(0, 110),
                cex.names = 0.55, cex.axis = 0.75, cex.lab = 0.85)
  mtext("c   TAR1: near-universal except PAR1",
        side = 3, line = 1.2, adj = 0, cex = 0.95, font = 2)
  abline(h = 94.6, col = "grey30", lty = 2, lwd = 1.0)
  text(par("usr")[2], 97, " all-PHR mean = 94.6%", adj = c(1, 0),
       cex = 0.65, col = "grey30")
  par1_idx <- which(d$is_par1)
  if (length(par1_idx) > 0) {
    label_text <- paste(sprintf("%s %.1f%%", d$short[par1_idx],
                                d$pct_with_tar1[par1_idx]), collapse = ", ")
    text(bp[par1_idx[1]], 18,
         paste0("PAR1: ", label_text),
         adj = c(0, 0.5), cex = 0.75, col = "#1f77b4", font = 2)
  }
  text(par("usr")[2] * 0.98, 50,
       "PAR1 (chrXp/chrYp) has obligate\nmeiotic crossover -- no satellite-\nmediated exchange anchor needed.",
       adj = c(1, 1), cex = 0.7, col = "#1f77b4")
}

render_to <- function(dev_open, dev_close) {
  dev_open()
  layout(matrix(1:3, nrow = 1, byrow = TRUE))
  panel_a()
  panel_b()
  panel_c()
  dev_close()
}

pdf_path <- file.path(OUT, "slide_14_gene_biology.pdf")
png_path <- file.path(OUT, "slide_14_gene_biology.png")

render_to(function() pdf(pdf_path, width = 16, height = 5.2), dev.off)
render_to(function() png(png_path, width = 16 * 200, height = 5.2 * 200,
                         res = 200, type = "cairo"), dev.off)

cat("Wrote", pdf_path, "and", png_path, "\n")
