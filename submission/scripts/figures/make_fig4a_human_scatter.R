#!/usr/bin/env Rscript
# Main Figure 4a -- one dot per inter-chromosomal PHR sequence pair
# (single-sample: each sample's own PHRs vs that sample's contact), 50 kbp.
# x = PHR-pair Jaccard similarity; y = length-normalized 3D contact (log).
# Base R only. Paths resolve from the script's location, so run from anywhere:
#   Rscript submission/scripts/figures/make_fig4a_human_scatter.R
# Inputs  (override dir with DATA_DIR=...):
#   data/human_HG002_porec_50000bp_seqlevel.tsv
#   data/human_CHM13_hic_50000bp_seqlevel.tsv
#   columns: ... chr_a arm_a chr_b arm_b ... jaccard hic_contact_raw hic_contact_norm hic_bins
# Outputs (override dir with OUT_DIR=...):
#   submission/fig/MainFigures/Fig4a_human_scatter.{png,pdf}

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

sources <- list(
  list(dataset = "HG002 Pore-C 50 kbp",
       path = file.path(data_dir, "human_HG002_porec_50000bp_seqlevel.tsv"),
       color = "#007a78"),
  list(dataset = "CHM13 Hi-C 50 kbp",
       path = file.path(data_dir, "human_CHM13_hic_50000bp_seqlevel.tsv"),
       color = "#b64b2a")
)

load_seq_pair <- function(path) {
  d <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  d <- d[d$chr_a != d$chr_b, ]                                   # inter-chromosomal only
  d[!is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
}

items <- lapply(sources, function(src) {
  d  <- load_seq_pair(src$path)
  ct <- suppressWarnings(cor.test(d$jaccard, d$hic_contact_norm,
                                  method = "spearman", exact = FALSE))
  list(src = src, d = d, n = nrow(d),
       rho = unname(ct$estimate), p = ct$p.value)
})

plot_one <- function(it) {
  x <- it$d$jaccard; y <- it$d$hic_contact_norm
  pos <- y[y > 0]; floor_y <- if (length(pos)) min(pos) / 2 else 1e-12
  y_plot <- pmax(y, floor_y)
  plot(x, y_plot, log = "y", pch = 21,
       bg = adjustcolor(it$src$color, alpha.f = 0.22),
       col = adjustcolor("#1f1f1f", alpha.f = 0.12),
       lwd = 0.25, cex = 0.55, xaxs = "i", xlim = c(0, 1), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       main = it$src$dataset,
       cex.main = 1.56, cex.lab = 1.56, cex.axis = 1.36)
  grid(col = "#e6e6e6", lwd = 0.7)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0), cex.axis = 1.36)
  yd <- seq(ceiling(par("usr")[3]), floor(par("usr")[4]))
  axis(2, at = 10^yd, las = 1, cex.axis = 1.36,
       labels = vapply(yd, function(d)
                       if (d >= 0) formatC(10^d, format = "d")
                       else if (d >= -3) formatC(10^d, format = "f", digits = -d)
                       else sprintf("1e%03d", d), ""))
  if (length(unique(x)) > 2) {
    fit <- lm(log10(y_plot) ~ x)
    xs <- seq(min(x), max(x), length.out = 100)
    lines(xs, 10 ^ predict(fit, data.frame(x = xs)), col = "#111111", lwd = 1.35)
  }
  legend("bottomright", inset = c(0.01, 0.09), bty = "n", cex = 1.22,
         text.col = "#222222",
         legend = c(sprintf("n = %s PHR pairs", format(it$n, big.mark = ",")),
                    sprintf("pointwise Spearman rho = %s", fmt_rho(it$rho)),
                    sprintf("p = %s", fmt_p(it$p))))
  legend("topleft", legend = "y axis: log scale; 0 shown at floor", bty = "n",
         cex = 0.84, text.col = "black", text.font = 3, inset = c(-0.05, -0.015))
}

draw <- function() {
  par(mfrow = c(1, 2), oma = c(2.6, 2.0, 0.2, 0.2),
      mar = c(2.7, 4.0, 2.1, 0.7), mgp = c(2.4, 0.7, 0), family = "sans")
  invisible(lapply(items, plot_one))
  mtext("PHR sequence-pair Jaccard similarity", side = 1, outer = TRUE,
        line = 0.8, cex = 1.68)
  mtext("3D contact frequency", side = 2, outer = TRUE, line = 0.6, cex = 1.68)
}

png(file.path(out_dir, "Fig4a_human_scatter.png"),
    width = 1840, height = 880, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "Fig4a_human_scatter.pdf"),
    width = 10.22, height = 4.89); draw(); dev.off()

for (it in items)
  cat(sprintf("%-20s n=%d  rho=%s  p=%s\n",
              it$src$dataset, it$n, fmt_rho(it$rho), fmt_p(it$p)))
cat("wrote ", file.path(out_dir, "Fig4a_human_scatter.{png,pdf}"), "\n", sep = "")
