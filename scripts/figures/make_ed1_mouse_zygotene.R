#!/usr/bin/env Rscript
# Extended Data Fig. 1 -- mouse meiotic Hi-C.
# Left: per arm-pair zygotene Hi-C contact (log) vs mean PHR Jaccard similarity
#       (pointwise Spearman). Right: Mantel rho across meiotic prophase stages
#       (per-stage values from the manuscript/Methods; zygotene = bouquet).
# Base R only. Run from the repo root:
#   Rscript scripts/figures/make_ed1_mouse_zygotene.R
# Input  (override dir with DATA_DIR=...):
#   data/zuo2021_zygotene_phr_pair_correlation.tsv
#   columns: arm_a  arm_b  mean_jaccard  hic_contact
# Output (override dir with OUT_DIR=...):
#   paper_prep/submission/fig/ExtendedDataFigures/ED_Fig1_mouse_zygotene.{png,pdf}

data_dir <- Sys.getenv("DATA_DIR", "data")
out_dir  <- Sys.getenv("OUT_DIR",  "paper_prep/submission/fig/ExtendedDataFigures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_rho <- function(x) formatC(x, format = "f", digits = 3)
fmt_p   <- function(p) if (is.na(p) || p == 0) "<1e-300" else
                       formatC(p, format = "e", digits = 1)

d <- read.delim(file.path(data_dir, "zuo2021_zygotene_phr_pair_correlation.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)
chrom <- function(x) sub("_[pq]$", "", x)
d <- d[chrom(d$arm_a) != chrom(d$arm_b) &
       !is.na(d$mean_jaccard) & !is.na(d$hic_contact) & d$hic_contact > 0, ]
ct  <- suppressWarnings(cor.test(d$mean_jaccard, d$hic_contact,
                                 method = "spearman", exact = FALSE))
rho <- unname(ct$estimate); n <- nrow(d)

# Per-stage Mantel rho (manuscript values; zygotene is the bouquet stage).
stages <- data.frame(
  stage = factor(c("lepto", "zygo", "pachy", "diplo"),
                 levels = c("lepto", "zygo", "pachy", "diplo")),
  rho   = c(0.687, 0.718, 0.683, 0.577))

draw <- function() {
  par(mfrow = c(1, 2), mar = c(5.1, 5.1, 4.2, 1.2), family = "sans")

  # left: zygotene scatter
  y <- d$hic_contact; floor_y <- min(y[y > 0]) / 2
  plot(d$mean_jaccard, pmax(y, floor_y), log = "y", pch = 21,
       bg = adjustcolor("#1f77b4", alpha.f = 0.45),
       col = adjustcolor("#1f1f1f", alpha.f = 0.30), lwd = 0.35, cex = 0.9,
       xlab = "Mean PHR Jaccard similarity per arm pair",
       ylab = "Zygotene Hi-C contact (log scale)",
       main = "Mouse: sequence-similar subtelomeres contact more",
       cex.main = 0.98, cex.lab = 0.92, cex.axis = 0.80)
  grid(col = "#e6e6e6", lwd = 0.7)
  fit <- lm(log10(pmax(y, floor_y)) ~ d$mean_jaccard)
  xs <- seq(min(d$mean_jaccard), max(d$mean_jaccard), length.out = 100)
  lines(xs, 10 ^ (coef(fit)[1] + coef(fit)[2] * xs), col = "#111111", lwd = 1.35)
  legend("topleft", bty = "n", cex = 0.74, text.col = "#222222",
         legend = c(sprintf("n = %s pairs", format(n, big.mark = ",")),
                    sprintf("Spearman rho = %s", fmt_rho(rho)),
                    sprintf("p = %s", fmt_p(ct$p.value))))

  # right: stage trajectory
  plot(seq_len(nrow(stages)), stages$rho, type = "b", pch = 21,
       bg = ifelse(stages$stage == "zygo", "#d62728", "#1f77b4"),
       col = "#111111", lwd = 1.2, cex = 1.6, xaxt = "n", ylim = c(0.55, 0.75),
       xlab = "Meiotic prophase stage", ylab = "Mantel rho (sequence vs Hi-C)",
       main = "Stage trajectory (bouquet peak)", cex.main = 0.98, cex.lab = 0.92)
  axis(1, at = seq_len(nrow(stages)), labels = levels(stages$stage))
  text(seq_len(nrow(stages)), stages$rho + 0.013,
       sprintf("%.3f", stages$rho), cex = 0.72)
  text(2, 0.732, "bouquet", col = "#d62728", font = 2, cex = 0.8)
}

png(file.path(out_dir, "ED_Fig1_mouse_zygotene.png"),
    width = 1800, height = 1050, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "ED_Fig1_mouse_zygotene.pdf"),
    width = 10, height = 5.83); draw(); dev.off()

cat(sprintf("mouse zygotene  n=%d  rho=%s  p=%s\n", n, fmt_rho(rho), fmt_p(ct$p.value)))
cat("wrote ", file.path(out_dir, "ED_Fig1_mouse_zygotene.{png,pdf}"), "\n", sep = "")
