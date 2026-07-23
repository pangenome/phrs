#!/usr/bin/env Rscript
# Extended Data Fig. 1 -- replication of the Fig 4a sequence<->3D coupling across
# additional human contact assays and a second genome, all at 20 kbp, all-points
# (one dot per inter-chromosomal PHR sequence pair; single-sample, no averaging):
#   panel 1  CHM13 Hi-C        (second genome, Hi-C)
#   panel 2  HG002 Hi-C        (same genome as Fig 4a Pore-C, different assay)
#   panel 3  HG002 CiFi        (third assay; sparse, hence weaker)
# x = PHR-pair Jaccard similarity; y = mean contact per bin-pair (density, log).
# Base R only. Paths resolve from the script's location, so run from anywhere:
#   Rscript submission/scripts/figures/make_ed1_human_contacts.R
# Input  (override dir with DATA_DIR=...):
#   data/human_CHM13_hic_20000bp_seqlevel.tsv
#   data/human_HG002_hic_20000bp_seqlevel.tsv
#   data/human_HG002_cifi_20000bp_seqlevel.tsv
# Outputs (override dir with OUT_DIR=...):
#   submission/fig/ExtendedDataFigures/ED_Fig1_human_contacts.{png,pdf}

.cmd_args  <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
script_dir <- if (length(.this_file)) normalizePath(dirname(.this_file)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir <- Sys.getenv("DATA_DIR", file.path(repo_root, "data"))
out_dir  <- Sys.getenv("OUT_DIR",  file.path(repo_root, "submission/fig/ExtendedDataFigures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_rho <- function(x) formatC(x, format = "f", digits = 3)

sources <- list(
  list(dataset = "CHM13 Hi-C", path = file.path(data_dir, "human_CHM13_hic_20000bp_seqlevel.tsv"),
       color = "#b64b2a"),
  list(dataset = "HG002 Hi-C", path = file.path(data_dir, "human_HG002_hic_20000bp_seqlevel.tsv"),
       color = "#2a6fb6"),
  list(dataset = "HG002 CiFi", path = file.path(data_dir, "human_HG002_cifi_20000bp_seqlevel.tsv"),
       color = "#7a4fb6")
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
  list(src = src, d = d, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value)
})

plot_one <- function(it) {
  x <- it$d$jaccard; y <- it$d$hic_contact_norm
  pos <- y[y > 0]; floor_y <- if (length(pos)) min(pos) / 2 else 1e-12
  y_plot <- pmax(y, floor_y)
  plot(x, y_plot, log = "y", pch = 21,
       bg = adjustcolor(it$src$color, alpha.f = 0.22),
       col = adjustcolor("#1f1f1f", alpha.f = 0.12),
       lwd = 0.25, cex = 0.55, xaxs = "i", xlim = c(0, 1), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = it$src$dataset, cex.main = 2.2,
       cex.lab = 1.56, cex.axis = 1.75)
  grid(col = "#e6e6e6", lwd = 0.7)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = 1.75)
  yd <- seq(ceiling(par("usr")[3]), floor(par("usr")[4]))
  axis(2, at = 10^yd, las = 1, cex.axis = 1.75,
       labels = vapply(yd, function(d)
                       if (d >= 0) formatC(10^d, format = "d")
                       else if (d >= -3) formatC(10^d, format = "f", digits = -d)
                       else sprintf("1e%03d", d), ""))
  if (length(unique(x)) > 2) {
    fit <- lm(log10(y_plot) ~ x)
    xs <- seq(min(x), max(x), length.out = 100)
    lines(xs, 10 ^ predict(fit, data.frame(x = xs)), col = "#111111", lwd = 1.35)
  }
  legend("bottomright", inset = c(0.01, 0.05), bty = "n", cex = 1.55,
         text.col = "#222222",
         legend = c(sprintf("n = %s PHR pairs", format(it$n, big.mark = ",")),
                    sprintf("Spearman rho = %s",
                            fmt_rho(it$rho))))
}

draw <- function() {
  par(mfrow = c(1, 3), oma = c(3.0, 3.0, 0.4, 0.4),
      mar = c(2.7, 4.7, 2.2, 1.4), mgp = c(2.4, 0.7, 0), family = "sans")
  invisible(lapply(items, plot_one))
  mtext("PHR sequence-pair Jaccard similarity", side = 1, outer = TRUE,
        line = 1.0, cex = 1.5)
  mtext("3D contact frequency", side = 2, outer = TRUE, line = 1.0, cex = 1.5)
}

png(file.path(out_dir, "ED_Fig1_human_contacts.png"),
    width = 2700, height = 900, res = 200, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "ED_Fig1_human_contacts.pdf"),
    width = 13.5, height = 4.5); draw(); dev.off()

for (it in items)
  cat(sprintf("%-12s n=%d  descriptive pointwise rho=%s\n",
              it$src$dataset, it$n, fmt_rho(it$rho)))
cat("wrote ", file.path(out_dir, "ED_Fig1_human_contacts.{png,pdf}"), "\n", sep = "")
