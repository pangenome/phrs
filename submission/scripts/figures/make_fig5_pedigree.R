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
COL_B <- "#D95F02"; COL_C <- "#6f6f6f"; COL_D <- "#1F77B4"   # panel-identifier colours
GG_TS <- function(svg) svg * 1.55 / 11.856
CW <- 3600; CH <- 350
PAN028_X <- 1150             # PAN028 pulled in close to PAN027 (short descent branch)
PED_CX <- (560 + PAN028_X) / 2   # horizontal centre of the pedigree drawing

# pedigree symbols (SVG-like coords, y grows downward)
ped <- data.frame(x = c(560, 1010, 785, PAN028_X), y = c(155, 155, 240, 240),
                  shape = c(22, 21, 21, 21),
                  lab = c("PAN011 father", "PAN010 mother", "PAN027", "PAN028"))
# first generation (top row) labelled ABOVE the shapes, children BELOW
ped$laby <- ifelse(ped$y < 200, ped$y - 52, ped$y + 60)
lines <- data.frame(x = c(560, 785, 785), xend = c(1010, 785, PAN028_X),
                    y = c(155, 155, 240), yend = c(155, 240, 240))
bcd <- data.frame(y = c(140, 195, 250), let = c("B", "C", "D"),
  txt = c("PAN027 paternal haplotype vs PAN011 (father)",
          "PAN027 maternal haplotype vs PAN010 (mother)",
          "PAN028 maternal haplotype vs PAN027 (mother)"))

p <- ggplot() +
  geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend),
               color = "#333333", linewidth = 0.7, lineend = "round") +
  geom_point(data = ped, aes(x = x, y = y, shape = shape), fill = "white",
             color = "black", size = 6, stroke = 1.2) +
  geom_text(aes(x = PED_CX, y = 40, label = "Pedigree comparisons"),
            size = GG_TS(40), fontface = "bold", color = TEXT, hjust = 0.5) +
  geom_text(data = ped, aes(x = x, y = laby, label = lab),
            size = GG_TS(23), fontface = "bold", color = TEXT) +
  geom_text(data = bcd, aes(x = 1440, y = y, label = let),
            size = GG_TS(27), fontface = "bold", color = TEXT, hjust = 0) +
  geom_text(data = bcd, aes(x = 1495, y = y, label = txt),
            size = GG_TS(25), color = TEXT, hjust = 0) +
  scale_shape_identity() + scale_fill_identity() +
  scale_y_reverse(limits = c(CH, 0), expand = c(0, 0)) +
  # x-window shifted so the drawn content (centre ~1404) sits at the canvas centre
  scale_x_continuous(limits = c(1404 - CW / 2, 1404 + CW / 2), expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_void() + theme(legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out, "Fig5A_pedigree.pdf"), p, width = 12, height = CW / CH * 0 + CH / 300, bg = "white")
ggsave(file.path(out, "Fig5A_pedigree.png"), p, width = 12, height = CH / 300, dpi = 300, bg = "white")
cat("wrote", file.path(out, "Fig5A_pedigree.{pdf,png}"), "\n")
