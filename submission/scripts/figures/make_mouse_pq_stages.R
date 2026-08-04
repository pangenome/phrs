#!/usr/bin/env Rscript
# NOT A MANUSCRIPT FIGURE. Exploratory companion to Fig 4c: the per-stage
# per-PHR-pair Spearman trajectory split into p-arm-only and q-arm-only pairs,
# stacked vertically, styled to match the Fig 4c trajectory panel.
# Restricts to inter-chromosomal pairs where BOTH PHRs are on p-arms (top) or
# BOTH on q-arms (bottom). Same input files and rho computation as Fig 4c.
#   Rscript submission/scripts/figures/make_mouse_pq_stages.R
# Input  (override dir with DATA_DIR=...):
#   data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_20000bp_seqlevel.tsv
# Output (override dir with OUT_DIR=...): mouse_pq_stages.{png,pdf} at repo root.

.cmd_args  <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
script_dir <- if (length(.this_file)) normalizePath(dirname(.this_file)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir <- Sys.getenv("DATA_DIR", file.path(repo_root, "data"))
out_dir  <- Sys.getenv("OUT_DIR",  repo_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

RES <- "20000bp"
stage_keys <- c(lepto = "leptotene", zygo = "zygotene",
                pachy = "pachytene", diplo = "diplotene")

load_stage <- function(s) {
  d <- read.delim(file.path(data_dir,
                  sprintf("mouse_%s_phr_%s_seqlevel.tsv", s, RES)),
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
}
dl <- lapply(stage_keys, load_stage)

# arm class from the trailing character of arm_a / arm_b (e.g. "chr10_p" -> "p")
arm_of <- function(v) substr(v, nchar(v), nchar(v))
rho_arm <- function(d, a) {
  s <- d[arm_of(d$arm_a) == a & arm_of(d$arm_b) == a, ]      # both PHRs same arm
  r <- suppressWarnings(cor.test(s$jaccard, s$hic_contact_norm,
                                 method = "spearman", exact = FALSE))
  c(rho = unname(r$estimate), n = nrow(s))
}
P <- sapply(dl, rho_arm, a = "p")
Q <- sapply(dl, rho_arm, a = "q")

BLUE <- "#1f77b4"
YLIM <- c(-0.02, 0.60)                    # shared across panels for comparison

panel <- function(rhos, nlab, show_x) {
  plot(seq_len(4), rhos, type = "b", pch = 21, bg = BLUE,
       col = "#111111", lwd = 1.4, cex = 2.1, xaxt = "n", yaxt = "n",
       xlim = c(0.55, 4.45), ylim = YLIM,
       xlab = if (show_x) "Meiotic prophase stage" else "",
       ylab = "Spearman rho", main = "", cex.lab = 1.92)
  axis(2, at = seq(0, 0.6, 0.1), las = 1, cex.axis = 1.68)
  axis(1, at = seq_len(4), labels = FALSE)
  if (show_x)
    mtext(unname(stage_keys), side = 1, at = seq_len(4), line = 0.9, cex = 1.3)
  text(seq_len(4) + c(0, 0, 0, 0.16), rhos + 0.05,
       sprintf("%.3f", rhos), cex = 1.5, xpd = NA)
  mtext(nlab, side = 3, line = 0.25, cex = 1.45, font = 2)
}

draw <- function() {
  par(mfrow = c(2, 1), family = "sans",
      mgp = c(3.6, 1.0, 0),
      mar = c(2.4, 5.8, 2.6, 1.0))
  panel(P["rho", ], sprintf("p-arm PHR pairs  (n = %d)", P["n", 1]), FALSE)
  par(mar = c(4.8, 5.8, 2.6, 1.0))
  panel(Q["rho", ], sprintf("q-arm PHR pairs  (n = %d)", Q["n", 1]), TRUE)
}

png(file.path(out_dir, "mouse_pq_stages.png"),
    width = 1100, height = 1250, res = 180, type = "cairo"); draw(); dev.off()
pdf(file.path(out_dir, "mouse_pq_stages.pdf"),
    width = 6.1, height = 6.9); draw(); dev.off()

cat("p-arm rho:", sprintf("%.3f", P["rho", ]), " n=", P["n", 1], "\n")
cat("q-arm rho:", sprintf("%.3f", Q["rho", ]), " n=", Q["n", 1], "\n")
cat("wrote ", file.path(out_dir, "mouse_pq_stages.{png,pdf}"), "\n", sep = "")
