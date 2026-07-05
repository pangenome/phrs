#!/usr/bin/env Rscript
# DEBUG: panel C (mouse) reproduced at every Hi-C resolution.
# One row per resolution (5/10/20/50/100 kb): LEFT = zygotene scatter
# (PHR-pair Jaccard vs density contact, hic_contact_norm, log y); RIGHT = the
# per-stage per-PHR-pair Spearman trajectory (lepto/zygo/pachy/diplo).
# 1 Mb flank.  Shows where the zygotene bouquet peak is present (it is a
# fine-resolution feature).
# Input: data/mouse_meiosis_sweep/seqlevel/1Mb/mouse_<stage>_phr_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): fig4c_resolution_debug.{png,pdf}
# Run: Rscript scripts/mouse/fig4c_resolution_debug.R

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
dir_in <- file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel/1Mb")
out_dir <- Sys.getenv("OUT_DIR", "/tmp")
res_list <- c(5000, 10000, 20000, 50000, 100000)
stages <- c(lepto="leptotene", zygo="zygotene", pachy="pachytene", diplo="diplotene")

load_stage <- function(stage, res) {
  f <- list.files(dir_in, pattern = sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", stage, res),
                  full.names = TRUE)[1]
  d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
}
sp <- function(a, b) suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE))
fmt <- function(x) formatC(x, format = "f", digits = 3)

draw <- function() {
  par(mfcol = c(2, length(res_list)), mar = c(4.0, 4.3, 2.2, 0.8),
      mgp = c(2.5, 0.8, 0), family = "sans")
  for (res in res_list) {
    dl  <- lapply(stages, load_stage, res = res)
    rho <- vapply(dl, function(d) unname(sp(d$jaccard, d$hic_contact_norm)$estimate), numeric(1))
    # --- top: zygotene scatter ---
    d <- dl[["zygo"]]; x <- d$jaccard; y <- d$hic_contact_norm
    pos <- y[y > 0]; fl <- if (length(pos)) min(pos)/2 else 1e-9; yp <- pmax(y, fl)
    plot(x, yp, log = "y", pch = 21, xlim = c(0, 1),
         bg = adjustcolor("#1f77b4", 0.22), col = adjustcolor("#1f1f1f", 0.12),
         lwd = 0.2, cex = 0.5, xlab = "PHR Jaccard", ylab = "zygotene contact",
         main = sprintf("%d kb", res/1000), cex.main = 1.5, cex.lab = 1.1)
    grid(col = "#ededed")
    if (length(unique(x)) > 2) { fit <- lm(log10(yp) ~ x); xs <- seq(0, 1, length.out = 80)
      lines(xs, 10^predict(fit, data.frame(x = xs)), col = "#111", lwd = 1.2) }
    legend("bottomright", bty = "n", cex = 1.05, text.col = "#222",
           legend = c(sprintf("rho = %s", fmt(rho["zygo"])), sprintf("n = %d", nrow(d))))
    # --- bottom: per-stage trajectory ---
    yl <- range(rho)
    plot(seq_along(rho), rho, type = "b", pch = 21,
         bg = c("#1f77b4","#d62728","#1f77b4","#1f77b4"), col = "#111", lwd = 1.3, cex = 1.9,
         xaxt = "n", ylim = c(yl[1]-0.06, yl[2]+0.09), xlab = "stage",
         ylab = "per-pair rho", main = "", cex.lab = 1.1)
    axis(1, at = seq_along(rho), labels = names(stages), cex.axis = 0.95)
    text(seq_along(rho), rho + 0.035, fmt(rho), cex = 0.95, xpd = NA)
    pk <- which.max(rho)
    text(pk, rho[pk] + 0.075, "peak", col = "#d62728", font = 2, cex = 1.0, xpd = NA)
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(out_dir, "fig4c_resolution_debug.png"), width = 2600, height = 1100,
    res = 150, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "fig4c_resolution_debug.pdf"), width = 17.3, height = 7.3); draw(); dev.off()
cat("wrote ", file.path(out_dir, "fig4c_resolution_debug.{png,pdf}"), "\n", sep = "")
