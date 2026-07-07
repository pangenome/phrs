#!/usr/bin/env Rscript
# Main Figure 4c -- mouse meiotic Hi-C, all-points (per-PHR-pair), 20 kbp, 1 Mb flank.
# Left: zygotene scatter, one dot per inter-chromosomal mouse PHR sequence pair
#       (length-normalised Hi-C contact vs PHR Jaccard; no averaging, as Fig 4a).
# Right: per-stage pointwise Spearman trajectory, the SAME per-PHR-pair
#       statistic computed here from the four stage files. At 20 kbp it is
#       strongest at zygotene (rho 0.614) and lowest at diplotene (rho 0.245);
#       reported as descriptive point estimates, no formal stage-contrast test.
# Base R only. Paths resolve from the script's location, so run from anywhere:
#   Rscript submission/scripts/figures/make_fig4c_mouse_zygotene.R
# Input  (override dir with DATA_DIR=...):
#   data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_20000bp_seqlevel.tsv
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

RES <- "20000bp"   # per-PHR-pair contacts at 20 kbp (1 Mb flank)
stage_keys <- c(lepto = "leptotene", zygo = "zygotene",
                pachy = "pachytene", diplo = "diplotene")
load_stage <- function(s) {
  d <- read.delim(file.path(data_dir,
                  sprintf("mouse_%s_phr_%s_seqlevel.tsv", s, RES)),
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
}
dl  <- lapply(stage_keys, load_stage)
cts <- lapply(dl, function(d) suppressWarnings(
         cor.test(d$jaccard, d$hic_contact_norm, method = "spearman", exact = FALSE)))

# per-stage trajectory = the same per-PHR-pair Spearman, one value per stage
stages <- data.frame(
  stage = factor(names(stage_keys), levels = names(stage_keys)),
  rho   = vapply(cts, function(c) unname(c$estimate), numeric(1)))

# zygotene is the left-hand scatter
d   <- dl[["zygo"]]
ct  <- cts[["zygo"]]
rho <- unname(ct$estimate); n <- nrow(d)

draw <- function() {
  par(mfrow = c(1, 2), mar = c(5.6, 5.8, 1.4, 1.0), mgp = c(3.6, 1.0, 0),
      family = "sans")

  # left: all-points zygotene scatter (no title; inset stats bottom-right as in 4a)
  x <- d$jaccard; y <- d$hic_contact_norm
  pos <- y[y > 0]; floor_y <- if (length(pos)) min(pos) / 2 else 1e-9
  y_plot <- pmax(y, floor_y)
  plot(x, y_plot, log = "y", pch = 21, xlim = c(0, 1),
       bg = adjustcolor("#1f77b4", alpha.f = 0.22),
       col = adjustcolor("#1f1f1f", alpha.f = 0.12), lwd = 0.25, cex = 0.6,
       xlab = "PHR-pair Jaccard similarity",
       ylab = "Zygotene 3D contact frequency",
       main = "", cex.lab = 1.55, cex.axis = 1.4)
  grid(col = "#e6e6e6", lwd = 0.7)
  if (length(unique(x)) > 2) {
    fit <- lm(log10(y_plot) ~ x)
    xs <- seq(min(x), max(x), length.out = 100)
    lines(xs, 10 ^ predict(fit, data.frame(x = xs)), col = "#111111", lwd = 1.35)
  }
  legend("bottomright", inset = c(0.02, 0.04), bty = "n", cex = 1.15,
         text.col = "#222222",
         legend = c(sprintf("n = %s PHR pairs", format(n, big.mark = ",")),
                    sprintf("Spearman rho = %s",
                            fmt_rho(rho))))
  legend("topleft", legend = "y axis: log scale; 0 shown at floor", bty = "n",
         cex = 1.1, text.col = "black", text.font = 3, inset = c(-0.04, -0.01))

  # right: per-stage trajectory (same per-PHR-pair Spearman statistic; no title)
  yl <- range(stages$rho)
  plot(seq_len(nrow(stages)), stages$rho, type = "b", pch = 21,
       bg = "#1f77b4",
       col = "#111111", lwd = 1.4, cex = 2.1, xaxt = "n",
       xlim = c(0.55, nrow(stages) + 0.45),
       ylim = c(yl[1] - 0.06, yl[2] + 0.07),
       xlab = "Meiotic prophase stage", ylab = "Per-pair Spearman rho",
       main = "", cex.lab = 1.55, cex.axis = 1.4)
  axis(1, at = seq_len(nrow(stages)), labels = FALSE)
  mtext(unname(stage_keys[levels(stages$stage)]), side = 1,
        at = seq_len(nrow(stages)), line = 0.9, cex = 1.05)
  text(seq_len(nrow(stages)), stages$rho + 0.030,
       sprintf("%.3f", stages$rho), cex = 1.2, xpd = NA)
}

png(file.path(out_dir, "Fig4c_mouse_zygotene.png"),
    width = 1950, height = 1050, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "Fig4c_mouse_zygotene.pdf"),
    width = 10.83, height = 5.83); draw(); dev.off()

cat(sprintf("mouse zygotene  n=%d  descriptive pointwise rho=%s\n", n, fmt_rho(rho)))
cat("wrote ", file.path(out_dir, "Fig4c_mouse_zygotene.{png,pdf}"), "\n", sep = "")
