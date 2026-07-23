#!/usr/bin/env Rscript
# Fig 5 panel B -- identity-crossover zoom(s) for detected exchange events.
# Default run writes the chosen 9q/3q panel as the paper asset
# fig/MainFigures/Fig5B_zoom9q3q.{pdf,png}. Passing --out=<dir> instead writes a
# review PNG per candidate (Fig5B_zoom_<id>.png) for choosing among events.
# Along a query subtelomere, plot best-match sequence identity to its own
# chromosome (same_chrom winner) versus to the donor / best non-homologous
# chromosome (interchrom winner). The dip in own-chromosome identity across the
# called tract, where the donor stays near-identical, is the exchange signal.
# Reads the vendored 10:10 IMPG class-winner tables.

suppressPackageStartupMessages({library(ggplot2); library(data.table)})

.args <- commandArgs(trailingOnly = FALSE)
.f    <- sub("^--file=", "", .args[grep("^--file=", .args)])
root  <- if (length(.f)) normalizePath(file.path(dirname(.f), "..", "..", "..")) else getwd()
mfig  <- file.path(root, "submission/fig/MainFigures")
rarg  <- .args[grep("^--out=", .args)]
REVIEW <- length(rarg) > 0
outdir <- if (REVIEW) sub("^--out=", "", rarg) else mfig
pat <- file.path(root, "data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz")
p28 <- file.path(root, "data/fig5_PAN028mat_vs_PAN027_joint.class_winners.impg_similarity.tsv.gz")
TEXT <- "#202124"; GREY <- "#6f6f6f"

cfgs <- list(
  list(id="9q3q", file=pat, qchr="chr9",  lo=136130000, hi=136190000, tlo=136164000, thi=136188000,
       donor="PAN011 h2 chr3q", own="PAN011 h2 chr9q",  col="#D95F02", ylo=0.97, xlab="PAN027 paternal chr9q position (Mb)",  tract="9q/3q tract, 20 kb"),
  list(id="5q1p", file=pat, qchr="chr5",  lo=182020000, hi=182100000, tlo=182052000, thi=182080000,
       donor="chr1p", own="chr5",  col="#1F77B4", ylo=0.95, xlab="PAN027 paternal chr5q position (Mb)",  tract="5q/1p tract, 28 kb"),
  list(id="PAR1", file=pat, qchr="chrX",  lo=0,        hi=200000,    tlo=14000,     thi=156000,
       donor="chrY",  own="chrX",  col="#E7298A", ylo=0.20, xlab="PAN027 paternal chrX PAR1 position (Mb)", tract="PAR1, 138 kb"),
  list(id="21p4p", file=p28, qchr="chr21", lo=6240000, hi=6340000,   tlo=6276000,   thi=6316000,
       donor="chr4p", own="chr21", col="#2CA02C", ylo=0.40, xlab="PAN028 maternal chr21p position (Mb)", tract="21p/4p tract, 40 kb")
)

build <- function(cfg) {
  dt <- fread(cmd = sprintf("zcat %s", shQuote(cfg$file)))
  d  <- dt[query_chrom == cfg$qchr & start >= cfg$lo & end <= cfg$hi,
           .(mid = (start + 1000) / 1e6, winner_class, id = estimated.identity)]
  d[, series := ifelse(winner_class == "same_chrom", cfg$own, cfg$donor)]
  d$series <- factor(d$series, levels = c(cfg$donor, cfg$own))
  COL <- setNames(c(cfg$col, GREY), c(cfg$donor, cfg$own))
  yann <- cfg$ylo + 0.05 * (1 - cfg$ylo)
  ggplot(d, aes(mid, id, color = series)) +
    annotate("rect", xmin = cfg$tlo/1e6, xmax = cfg$thi/1e6, ymin = -Inf, ymax = Inf,
             fill = cfg$col, alpha = 0.08) +
    annotate("text", x = (cfg$tlo + cfg$thi)/2/1e6, y = yann,
             label = cfg$tract, size = 5.2, color = cfg$col, fontface = "bold") +
    geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
    scale_color_manual(values = COL, name = NULL) +
    scale_y_continuous("best-match\nidentity", limits = c(cfg$ylo, 1.002)) +
    scale_x_continuous(cfg$xlab) +
    theme_classic(base_size = 16) +
    theme(legend.position = c(0.02, 0.34), legend.justification = c(0, 1),
          legend.key.height = unit(15, "pt"), legend.background = element_blank(),
          legend.text = element_text(size = 16, color = TEXT),
          axis.title = element_text(size = 15, color = TEXT),
          axis.text  = element_text(size = 13, color = TEXT),
          plot.margin = margin(8, 12, 6, 6))
}

if (REVIEW) {
  for (cfg in cfgs) {
    ggsave(file.path(outdir, sprintf("Fig5B_zoom_%s.png", cfg$id)), build(cfg),
           width = 6.2, height = 3.1, dpi = 300, bg = "white")
    cat("wrote review", cfg$id, "\n")
  }
} else {
  p <- build(cfgs[[1]])   # 9q/3q is the chosen paper panel
  ggsave(file.path(mfig, "Fig5B_zoom9q3q.pdf"), p, width = 6.0, height = 2.4, bg = "white")
  ggsave(file.path(mfig, "Fig5B_zoom9q3q.png"), p, width = 6.0, height = 2.4, dpi = 300, bg = "white")
  cat("wrote canonical", file.path(mfig, "Fig5B_zoom9q3q.{pdf,png}"), "\n")
}
