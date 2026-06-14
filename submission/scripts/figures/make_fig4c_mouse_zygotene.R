#!/usr/bin/env Rscript
# Main Figure 4c -- mouse meiotic Hi-C.
# Left: all-points scatter, one dot per inter-chromosomal mouse PHR sequence pair
#       at zygotene (length-normalised Hi-C contact vs PHR Jaccard; no averaging,
#       matching Fig 4a). Right: Mantel rho across meiotic prophase stages
#       (arm-level matrix test; per-stage values from Methods; zygotene = bouquet).
# Base R only. Paths resolve from the script's location, so run from anywhere:
#   Rscript submission/scripts/figures/make_fig4c_mouse_zygotene.R
# Input  (override dir with DATA_DIR=...):
#   data/mouse_zygotene_phr_50000bp_seqlevel.tsv
#   columns: ... chr_a arm_a chr_b arm_b ... jaccard hic_contact_raw hic_contact_norm hic_bins
# Output (override dir with OUT_DIR=...):
#   submission/fig/MainFigures/Fig4c_mouse_zygotene.{png,pdf}

# Resolve the repo root from this script's own location
# (submission/scripts/figures/), so it runs from any working directory.
.cmd_args  <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
script_dir <- if (length(.this_file)) normalizePath(dirname(.this_file)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir <- Sys.getenv("DATA_DIR", file.path(repo_root, "data"))
out_dir  <- Sys.getenv("OUT_DIR",  file.path(repo_root, "submission/fig/MainFigures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_rho <- function(x) formatC(x, format = "f", digits = 3)
fmt_p   <- function(p) if (is.na(p) || p == 0) "<1e-300" else
                       formatC(p, format = "e", digits = 1)

d <- read.delim(file.path(data_dir, "mouse_zygotene_phr_50000bp_seqlevel.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)
d <- d[d$chr_a != d$chr_b &
       !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]   # inter-chromosomal PHR pairs
ct  <- suppressWarnings(cor.test(d$jaccard, d$hic_contact_norm,
                                 method = "spearman", exact = FALSE))
rho <- unname(ct$estimate); n <- nrow(d)

# Per-stage Mantel rho (manuscript values; zygotene is the bouquet stage).
stages <- data.frame(
  stage = factor(c("lepto", "zygo", "pachy", "diplo"),
                 levels = c("lepto", "zygo", "pachy", "diplo")),
  rho   = c(0.687, 0.718, 0.683, 0.577))

draw <- function() {
  par(mfrow = c(1, 2), mar = c(5.1, 5.1, 4.2, 1.2), family = "sans")

  # left: all-points zygotene scatter (one dot per inter-chromosomal PHR pair)
  x <- d$jaccard; y <- d$hic_contact_norm
  pos <- y[y > 0]; floor_y <- if (length(pos)) min(pos) / 2 else 1e-9
  y_plot <- pmax(y, floor_y)
  plot(x, y_plot, log = "y", pch = 21, xlim = c(0, 1),
       bg = adjustcolor("#1f77b4", alpha.f = 0.22),
       col = adjustcolor("#1f1f1f", alpha.f = 0.12), lwd = 0.25, cex = 0.6,
       xlab = "PHR-pair Jaccard similarity",
       ylab = "Zygotene Hi-C contact (length-normalised, log)",
       main = "Mouse: sequence-similar PHR pairs contact more",
       cex.main = 1.10, cex.lab = 1.18, cex.axis = 1.05)
  grid(col = "#e6e6e6", lwd = 0.7)
  if (length(unique(x)) > 2) {
    fit <- lm(log10(y_plot) ~ x)
    xs <- seq(min(x), max(x), length.out = 100)
    lines(xs, 10 ^ predict(fit, data.frame(x = xs)), col = "#111111", lwd = 1.35)
  }
  legend("topleft", bty = "n", cex = 0.95, text.col = "#222222",
         legend = c(sprintf("n = %s PHR pairs", format(n, big.mark = ",")),
                    sprintf("pointwise Spearman rho = %s", fmt_rho(rho)),
                    sprintf("p = %s", fmt_p(ct$p.value))))

  # right: stage trajectory
  plot(seq_len(nrow(stages)), stages$rho, type = "b", pch = 21,
       bg = ifelse(stages$stage == "zygo", "#d62728", "#1f77b4"),
       col = "#111111", lwd = 1.2, cex = 1.6, xaxt = "n", ylim = c(0.55, 0.78),
       xlab = "Meiotic prophase stage", ylab = "Mantel rho (sequence vs Hi-C)",
       main = "Stage trajectory (bouquet peak)", cex.main = 1.10, cex.lab = 1.18,
       cex.axis = 1.05)
  axis(1, at = seq_len(nrow(stages)), labels = levels(stages$stage), cex.axis = 1.05)
  text(seq_len(nrow(stages)), stages$rho + 0.017,
       sprintf("%.3f", stages$rho), cex = 0.90)
  text(2, 0.760, "bouquet", col = "#d62728", font = 2, cex = 0.95)
}

png(file.path(out_dir, "Fig4c_mouse_zygotene.png"),
    width = 1950, height = 1050, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "Fig4c_mouse_zygotene.pdf"),
    width = 10.83, height = 5.83); draw(); dev.off()

cat(sprintf("mouse zygotene  n=%d  rho=%s  p=%s\n", n, fmt_rho(rho), fmt_p(ct$p.value)))
cat("wrote ", file.path(out_dir, "Fig4c_mouse_zygotene.{png,pdf}"), "\n", sep = "")
