#!/usr/bin/env Rscript
# Fig 5 panel A -- pedigree of the WashU trio-of-trios comparisons.
# Standalone piece (no panel letter; the "A" is added in paper.tex via overpic).
# Output: submission/fig/MainFigures/Fig5A_pedigree.{pdf,png}, 3600 x 300 units
# to stack cleanly above the B/C/D ribbon panels (same 3600-wide coordinate).

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.f    <- sub("^--file=", "", .args[grep("^--file=", .args)])
root  <- if (length(.f)) normalizePath(file.path(dirname(.f), "..", "..", "..")) else getwd()
out   <- file.path(root, "submission/fig/MainFigures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

TEXT <- "#202124"
GG_TS <- function(svg) svg * 1.55 / 11.856

# tight framing window around the pedigree only (data units; y grows downward)
XLO <- 55; XHI <- 515; YLO <- -14; YHI <- 235
G1Y <- 60; G2Y <- 130; G3Y <- 200   # three generation rows (compact vertical spacing)

# Faithful WashU trio-of-trios topology (Cechova et al.):
#   gen-1  PAN011 (grandfather, square)  x  PAN010 (grandmother, circle)
#   gen-2  their daughter PAN027 (mother, circle)  x  unsampled father (open square)
#   gen-3  granddaughter PAN028 (circle), child of PAN027 x the unsampled father
# labels are sample IDs only (roles are given in the caption) to keep the panel compact
PAN011_X <- 205; PAN010_X <- 430; UNS_X <- 100
G1_MX <- (PAN011_X + PAN010_X) / 2   # gen-1 couple midpoint (descent drops to PAN027)
G2_MX <- (UNS_X + G1_MX) / 2         # gen-2 couple midpoint (descent drops to PAN028)

ped <- data.frame(
  x     = c(PAN011_X, PAN010_X, G1_MX, G2_MX),
  y     = c(G1Y,      G1Y,      G2Y,   G3Y),
  shape = c(22,       21,       21,    21),
  lab   = c("Father\n(PAN011)", "Mother\n(PAN010)", "Daughter\n(PAN027)", "Granddaughter\n(PAN028)"),
  labx  = c(PAN011_X, PAN010_X, G1_MX + 34, G2_MX + 34),
  laby  = c(G1Y - 48, G1Y - 48, G2Y, G3Y),
  hj    = c(0.5, 0.5, 0, 0))
# unsampled father of the granddaughter (open square, not sequenced)
partner <- data.frame(x = UNS_X, y = G2Y)
# gen-1 marriage bar, descent to PAN027, gen-2 marriage bar (PAN027 x partner), descent to PAN028
lines <- data.frame(
  x    = c(PAN011_X, G1_MX, UNS_X, G2_MX),
  xend = c(PAN010_X, G1_MX, G1_MX, G2_MX),
  y    = c(G1Y,      G1Y,   G2Y,   G2Y),
  yend = c(G1Y,      G2Y,   G2Y,   G3Y))

p <- ggplot() +
  geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend),
               color = "#333333", linewidth = 0.8, lineend = "round") +
  geom_point(data = partner, aes(x = x, y = y), shape = 22, fill = "grey88",
             color = "grey55", size = 9, stroke = 1.1) +
  geom_point(data = ped, aes(x = x, y = y, shape = shape), fill = "white",
             color = "black", size = 9, stroke = 1.4) +
  geom_text(data = ped, aes(x = labx, y = laby, label = lab, hjust = hj),
            size = GG_TS(24), fontface = "bold", color = TEXT, lineheight = 0.85) +
  geom_text(data = partner, aes(x = x, y = y + 42), label = "unsampled",
            size = GG_TS(22), color = "grey45", fontface = "italic") +
  scale_shape_identity() + scale_fill_identity() +
  scale_y_reverse(limits = c(YHI, YLO), expand = c(0, 0)) +
  scale_x_continuous(limits = c(XLO, XHI), expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_void() + theme(legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out, "Fig5A_pedigree.pdf"), p,
       width = 2.9, height = 2.9 * (YHI - YLO) / (XHI - XLO), bg = "white")
ggsave(file.path(out, "Fig5A_pedigree.png"), p,
       width = 2.9, height = 2.9 * (YHI - YLO) / (XHI - XLO), dpi = 300, bg = "white")
cat("wrote", file.path(out, "Fig5A_pedigree.{pdf,png}"), "\n")
