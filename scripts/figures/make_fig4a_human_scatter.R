#!/usr/bin/env Rscript
# Main Figure 4a -- Human sequence similarity vs 3D contact at 50 kb.
# Two panels (HG002 Pore-C, CHM13 Hi-C): one dot per inter-chromosomal arm pair,
# 3D contact (log) vs mean PHR Jaccard similarity, pointwise Spearman.
# Base R only. Run from the repo root:
#   Rscript scripts/figures/make_fig4a_human_scatter.R
# Inputs  (override dir with DATA_DIR=...):
#   data/hg002_porec_phr_pair_correlation.tsv
#   data/chm13_phr_pair_correlation.tsv
#   columns: arm_a  arm_b  mean_jaccard  n_pairs  hic_contact
# Outputs (override dir with OUT_DIR=...):
#   paper_prep/submission/fig/MainFigures/Fig4a_human_scatter.png  and  .pdf

data_dir <- Sys.getenv("DATA_DIR", "data")
out_dir  <- Sys.getenv("OUT_DIR",  "paper_prep/submission/fig/MainFigures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chrom_from_arm <- function(x) sub("_[pq]$", "", x)
fmt_rho <- function(x) formatC(x, format = "f", digits = 3)
fmt_p   <- function(p) if (is.na(p) || p == 0) "<1e-300" else
                       formatC(p, format = "e", digits = 1)

sources <- list(
  list(dataset = "HG002 Pore-C 50 kbp",
       path = file.path(data_dir, "hg002_porec_phr_pair_correlation.tsv"),
       color = "#007a78"),
  list(dataset = "CHM13 Hi-C 50 kbp",
       path = file.path(data_dir, "chm13_phr_pair_correlation.tsv"),
       color = "#b64b2a")
)

load_arm_pair <- function(path) {
  d <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  d <- d[chrom_from_arm(d$arm_a) != chrom_from_arm(d$arm_b), ]   # inter-chromosomal only
  d[!is.na(d$mean_jaccard) & !is.na(d$hic_contact), ]
}

items <- lapply(sources, function(src) {
  d  <- load_arm_pair(src$path)
  ct <- suppressWarnings(cor.test(d$mean_jaccard, d$hic_contact,
                                  method = "spearman", exact = FALSE))
  list(src = src, d = d, n = nrow(d),
       rho = unname(ct$estimate), p = ct$p.value)
})

plot_one <- function(it) {
  x <- it$d$mean_jaccard; y <- it$d$hic_contact
  pos <- y[y > 0]; floor_y <- if (length(pos)) min(pos) / 2 else 1e-12
  y_plot <- pmax(y, floor_y)
  plot(x, y_plot, log = "y", pch = 21,
       bg = adjustcolor(it$src$color, alpha.f = 0.42),
       col = adjustcolor("#1f1f1f", alpha.f = 0.32),
       lwd = 0.35, cex = 0.88, xaxs = "i", xlim = c(0, 1), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       main = it$src$dataset,
       cex.main = 1.32, cex.lab = 1.36, cex.axis = 1.12)
  grid(col = "#e6e6e6", lwd = 0.7)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0), cex.axis = 1.12)
  yd <- seq(ceiling(par("usr")[3]), floor(par("usr")[4]))
  axis(2, at = 10^yd, las = 1, cex.axis = 1.12,
       labels = vapply(yd, function(d)
                       if (d >= 0) formatC(10^d, format = "d")
                       else if (d >= -3) formatC(10^d, format = "f", digits = -d)
                       else sprintf("1e%03d", d), ""))
  if (length(unique(x)) > 2) {
    fit <- lm(log10(y_plot) ~ x)
    xs <- seq(min(x), max(x), length.out = 100)
    lines(xs, 10 ^ predict(fit, data.frame(x = xs)), col = "#111111", lwd = 1.35)
  }
  legend("bottomright", inset = c(0.01, 0.09), bty = "n", cex = 1.02,
         text.col = "#222222",
         legend = c(sprintf("n = %s arm pairs", format(it$n, big.mark = ",")),
                    sprintf("pointwise Spearman rho = %s", fmt_rho(it$rho)),
                    sprintf("p = %s", fmt_p(it$p))))
  legend("topleft", legend = "y axis: log scale; 0 shown at floor", bty = "n",
         cex = 0.74, text.col = "black", text.font = 3, inset = c(-0.02, 0.0))
}

draw <- function() {
  par(mfrow = c(1, 2), oma = c(2.2, 2.6, 0.2, 0.2),
      mar = c(2.3, 3.6, 1.8, 0.9), mgp = c(2.4, 0.7, 0), family = "sans")
  invisible(lapply(items, plot_one))
  mtext("Mean PHR Jaccard similarity per arm pair", side = 1, outer = TRUE,
        line = 0.8, cex = 1.45)
  mtext("3D contact frequency", side = 2, outer = TRUE, line = 0.9, cex = 1.45)
}

png(file.path(out_dir, "Fig4a_human_scatter.png"),
    width = 1760, height = 780, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "Fig4a_human_scatter.pdf"),
    width = 9.78, height = 4.33); draw(); dev.off()

for (it in items)
  cat(sprintf("%-20s n=%d  rho=%s  p=%s\n",
              it$src$dataset, it$n, fmt_rho(it$rho), fmt_p(it$p)))
cat("wrote ", file.path(out_dir, "Fig4a_human_scatter.{png,pdf}"), "\n", sep = "")
