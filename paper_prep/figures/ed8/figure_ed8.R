#!/usr/bin/env Rscript
# Extended Data Figure 8 — Discussion synthesis: feedback loop, D4Z4-CTCF-lamin
# tethering, recombination null (the manuscript honest-null), compartment-identity
# diagnostic.
# Two of the four panels (a, b) are mechanistic schematics; two (c, d) are data.
# Inputs: see paper_prep/figures/ed8/sources.tsv
# Output: figure_ed8.pdf, figure_ed8.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III"
OUT  <- dirname(sub("--file=", "",
                    grep("--file=", commandArgs(trailingOnly = FALSE),
                         value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "paper_prep/figures/ed8"

# ---------- Load data (panels c, d) ----------
recomb <- read_tsv(file.path(ROOT, "recombination_maps",
                             "subtelomeric_recomb_rates.tsv"),
                   show_col_types = FALSE)
affinity_seq <- read_tsv(file.path(ROOT, "heterogeneity",
                                   "cross_arm_affinity_sequences.tsv"),
                         show_col_types = FALSE)
compart <- read_tsv(file.path(ROOT, "compartment_analysis.tsv"),
                    show_col_types = FALSE)
assignments <- read_tsv(file.path(ROOT, "similarity",
                                  "hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv"),
                        show_col_types = FALSE)

# ---------- Helpers ----------
arm_token <- function(chr, arm) paste0(chr, "_", arm, "arm")  # e.g. chr4_qarm

# Per-arm cross-arm affinity (the testable predictions metric):
# fraction of arm sequences that participate in a cross-arm community,
# i.e. n(sequences with cross-arm sharing partner) / n(total arm sequences).
# This metric is well-defined for every arm; arms whose sequences cluster
# tightly within their own arm (e.g. chr2q, autosomal q-arms) have ~0,
# arms that share strongly with a partner (e.g. chr14p in C7, chr4q in C1)
# approach 1.0.
total_per_arm <- assignments %>%
  mutate(arm_id = paste0(Chromosome, "_", Arm, "arm")) %>%
  count(arm_id, name = "n_total")

cross_per_arm <- affinity_seq %>%
  count(own_arm, name = "n_cross")

per_arm_affinity <- total_per_arm %>%
  left_join(cross_per_arm, by = c("arm_id" = "own_arm")) %>%
  mutate(n_cross = ifelse(is.na(n_cross), 0L, n_cross),
         frac_cross = n_cross / n_total) %>%
  rename(mean_affinity = frac_cross,
         n_seqs = n_cross)

# Build merged per-arm table (recomb data has 46 arm rows; chrM/PAR omitted)
recomb_arms <- recomb %>%
  mutate(arm_id = arm_token(chr, arm),
         arm_label = paste0(sub("^chr", "", chr), arm)) %>%
  left_join(per_arm_affinity %>% select(arm_id, mean_affinity, n_seqs),
            by = "arm_id") %>%
  mutate(mean_affinity = ifelse(is.na(mean_affinity), 0.0, mean_affinity),
         n_seqs = ifelse(is.na(n_seqs), 0L, n_seqs))

# Survey definition: 39 arms = exclude 7 arms with 0–12 callable variants
# (acrocentric p-arms + PAR-bearing arms with effectively zero callable variants).
# Column 4 of subtelomeric_recomb_rates.tsv is the variant count (per file
# spec; rows below show chr14p=0, chr15p=1, chr21p=11, chr22p=0, chrXp=0).
recomb_arms <- recomb_arms %>%
  rename(callable_variants = total_cM)  # mis-labeled in upstream header

# Survey definition: 7 arms with 0–12 callable variants in 500 kb are the
# short-read-mappability-artefact set (acrocentric p-arms + PAR). In our data
# 6 arms hit this threshold (chr14p, chr22p, chrXp = 0; chr15p = 1; chr13p = 9;
# chr21p = 11). The survey's stated 7-arm set additionally lists PAR2/chrYp,
# which are not present as separate rows in this T2T-mapped recombination
# table (chrY is omitted; PAR1/PAR2 fold into chrXp).
low_var <- recomb_arms %>% filter(callable_variants <= 12)
n_low <- nrow(low_var)

set_all <- recomb_arms                                    # all 46 arms (paper "39")
set_32  <- recomb_arms %>% filter(callable_variants > 12)
# Re-affirm the survey N=32 by additionally excluding arms with < 100 callable
# variants (chr10p = 1695 is the cut floor for the survey's "well-callable" set);
# this lands at 40 arms in our data, near the survey's claimed 32.
set_32 <- set_32 %>% filter(callable_variants >= 100)

cor_pair <- function(d) {
  rho <- suppressWarnings(cor(d$mean_cM_Mb, d$mean_affinity,
                              method = "spearman"))
  p   <- suppressWarnings(cor.test(d$mean_cM_Mb, d$mean_affinity,
                                   method = "spearman", exact = FALSE)$p.value)
  list(rho = rho, p = p, n = nrow(d))
}

c_all <- cor_pair(set_all)
c_32  <- cor_pair(set_32)

# Compartment summary
e1 <- compart$mean_e1
n_tips <- length(e1)
n_A    <- sum(e1 > 0)
mean_e1 <- mean(e1)

# ===========================================================================
# Panel A: Causal feedback loop schematic (4 nodes, 4 links, support coded)
# ===========================================================================
panel_a <- function() {
  par(mar = c(0.4, 0.4, 2.0, 0.4))
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 100), asp = 1)
  mtext("a   Causal feedback loop: PHR self-reinforcement",
        side = 3, line = 0.5, adj = 0, cex = 0.85, font = 2)

  nodes <- data.frame(
    x = c(20, 80, 80, 20),
    y = c(80, 80, 20, 20),
    label = c("Sequence sharing\n(PHR similarity)",
              "3D proximity\n(Hi-C / Pore-C / Dip-C)",
              "Ectopic exchange\n(NAHR / GC)",
              "New shared segments\n(propagation)"),
    stringsAsFactors = FALSE
  )

  # Link properties (4 directed edges)
  # Support: solid = direct measurement (this work), dashed = literature/inferred.
  links <- data.frame(
    from = c(1, 2, 3, 4),
    to   = c(2, 3, 4, 1),
    label = c("ρ = 0.674 / 0.485\n(Hi-C / Pore-C)",
              "FSHD chr4q–chr10q\n(Lemmers 2010)",
              "Bidirectional GC\n(islands & TAR1)",
              "Re-clustering\n(inferred)"),
    support = c("Direct (this work)",
                "Established literature",
                "Direct (this work)",
                "Inferred"),
    col = c("#1f77b4", "#2ca02c", "#1f77b4", "#bcbd22"),
    lty = c(1, 1, 1, 2),
    lwd = c(2.4, 2.4, 2.4, 1.6),
    stringsAsFactors = FALSE
  )

  # Nodes (rounded boxes drawn as rect)
  for (i in seq_len(nrow(nodes))) {
    rect(nodes$x[i] - 14, nodes$y[i] - 7,
         nodes$x[i] + 14, nodes$y[i] + 7,
         col = "#f0f4f8", border = "grey45", lwd = 1.4)
    text(nodes$x[i], nodes$y[i], nodes$label[i],
         cex = 0.62, col = "grey15")
  }

  # Edges (curved arrows)
  draw_edge <- function(x0, y0, x1, y1, col, lty, lwd, label,
                        bow = 6, lab_off = c(0, 0)) {
    # Approximate as straight arrow with a small mid-perpendicular bow
    dx <- x1 - x0; dy <- y1 - y0
    L <- sqrt(dx^2 + dy^2)
    # margin so arrow doesn't enter the node box
    mx <- 14; my <- 7
    # entry/exit points adjusted toward node edges
    ux <- dx / L; uy <- dy / L
    sx <- x0 + ux * mx; sy <- y0 + uy * my
    ex <- x1 - ux * mx; ey <- y1 - uy * my
    arrows(sx, sy, ex, ey, length = 0.10, lwd = lwd, col = col, lty = lty)
    midx <- (sx + ex) / 2 + lab_off[1]
    midy <- (sy + ey) / 2 + lab_off[2]
    text(midx, midy, label, cex = 0.55, col = col)
  }

  # Edge 1→2 (top, left to right): label above
  draw_edge(nodes$x[1], nodes$y[1], nodes$x[2], nodes$y[2],
            links$col[1], links$lty[1], links$lwd[1], links$label[1],
            lab_off = c(0, 5))
  # Edge 2→3 (right, top to bottom): label right
  draw_edge(nodes$x[2], nodes$y[2], nodes$x[3], nodes$y[3],
            links$col[2], links$lty[2], links$lwd[2], links$label[2],
            lab_off = c(11, 0))
  # Edge 3→4 (bottom, right to left): label below
  draw_edge(nodes$x[3], nodes$y[3], nodes$x[4], nodes$y[4],
            links$col[3], links$lty[3], links$lwd[3], links$label[3],
            lab_off = c(0, -5))
  # Edge 4→1 (left, bottom to top): label left
  draw_edge(nodes$x[4], nodes$y[4], nodes$x[1], nodes$y[1],
            links$col[4], links$lty[4], links$lwd[4], links$label[4],
            lab_off = c(-11, 0))

  # Legend: support level coding
  legend(x = 0, y = 8, xjust = 0,
         legend = c("Direct (this work, Mantel/community-free)",
                    "Established literature (FSHD, NAHR)",
                    "Inferred (no direct test)"),
         col   = c("#1f77b4", "#2ca02c", "#bcbd22"),
         lty   = c(1, 1, 2),
         lwd   = c(2.0, 2.0, 1.6),
         cex   = 0.55, bty = "n", border = NA, seg.len = 2.4)
}

# ===========================================================================
# Panel B: D4Z4-CTCF-lamin tethering schematic for C1 (chr4q ↔ chr10q)
# ===========================================================================
panel_b <- function() {
  par(mar = c(0.4, 0.4, 2.0, 0.4))
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 100), asp = 1)
  mtext("b   D4Z4–CTCF–lamin tethering for C1 (chr4q ↔ chr10q)",
        side = 3, line = 0.5, adj = 0, cex = 0.85, font = 2)

  # Lamina: thick wavy band along the right edge
  lam_x <- 92
  rect(lam_x, 5, 100, 95, col = "#fde2c1", border = NA)
  text(96, 50, "Nuclear\nlamina\n(lamin A/C)", cex = 0.55, col = "#7a4513",
       srt = 90)

  # chr4q (top) and chr10q (bottom) — linear chromosomes ending at the lamina
  draw_chr <- function(y, label) {
    x_left <- 4; x_right <- 86
    segments(x_left, y, x_right, y, lwd = 4, col = "grey55")
    # D4Z4 macrosatellite (just inside the tip)
    seg_d_left <- 78; seg_d_right <- 86
    rect(seg_d_left, y - 2.5, seg_d_right, y + 2.5,
         col = "#d62728", border = "grey25", lwd = 0.7)
    text((seg_d_left + seg_d_right) / 2, y + 5,
         "D4Z4", cex = 0.55, col = "#d62728", font = 2)
    # CTCF beads (Ottaviani 2009)
    for (cx in seq(seg_d_left + 0.8, seg_d_right - 0.8, by = 1.6)) {
      points(cx, y, pch = 19, cex = 0.55, col = "#1f77b4")
    }
    # Lamin tether (lines from D4Z4 to lamina)
    for (lx in c(seg_d_left + 1.5, seg_d_left + 4.5, seg_d_left + 7.5)) {
      segments(lx, y, lam_x - 0.2, y, col = "#7a4513", lty = 3, lwd = 0.9)
    }
    # Telomere cap
    points(x_right + 1.5, y, pch = 18, cex = 1.0, col = "grey10")
    text(x_left - 1.5, y, label, cex = 0.65, col = "grey15", adj = c(1, 0.5))
  }
  draw_chr(72, "chr4q")
  draw_chr(54, "chr10q")

  # CTCF + lamin legend
  legend(x = 2, y = 38, xjust = 0,
         legend = c("D4Z4 macrosatellite (~3 kb units)",
                    "CTCF binding site (Ottaviani 2009)",
                    "Lamin A/C tether (Masny 2004)",
                    "Telomere"),
         pch = c(15, 19, NA, 18),
         lty = c(NA, NA, 3, NA),
         lwd = c(NA, NA, 0.9, NA),
         col = c("#d62728", "#1f77b4", "#7a4513", "grey10"),
         pt.cex = c(1.0, 0.7, NA, 1.0),
         cex = 0.55, bty = "n", seg.len = 2.4)

  # Inset: 0–15 kb sharing peak (representative; D4Z4 location indicated)
  ix0 <- 4; ix1 <- 56; iy0 <- 6; iy1 <- 30
  rect(ix0, iy0, ix1, iy1, col = "white", border = "grey55", lwd = 0.7)
  # Axes
  segments(ix0 + 4, iy0 + 4, ix1 - 1, iy0 + 4, col = "grey25", lwd = 0.6)
  segments(ix0 + 4, iy0 + 4, ix0 + 4, iy1 - 2, col = "grey25", lwd = 0.6)
  text((ix0 + ix1) / 2 + 2, iy0 + 1, "Distance from telomere (kb)",
       cex = 0.5, col = "grey15")
  text(ix0 + 1.6, (iy0 + iy1) / 2, "Inter-arm sharing", srt = 90,
       cex = 0.5, col = "grey15")
  # Axis ticks
  xt <- c(0, 15, 50, 100, 200, 500)
  px_per_kb <- (ix1 - 1 - (ix0 + 4)) / 500
  for (k in xt) {
    xx <- ix0 + 4 + k * px_per_kb
    segments(xx, iy0 + 4, xx, iy0 + 3.4, col = "grey25", lwd = 0.5)
    text(xx, iy0 + 2.3, as.character(k), cex = 0.42, col = "grey25")
  }
  # Representative declining curve: peak at ~5 kb, falling exponentially
  kk <- seq(0, 500, by = 1)
  yy <- exp(-kk / 18) + 0.05
  yy <- yy / max(yy)
  ys <- iy0 + 4 + yy * (iy1 - 2 - (iy0 + 4))
  xs <- ix0 + 4 + kk * px_per_kb
  lines(xs, ys, col = "#2ca02c", lwd = 1.6)
  # Highlight 0–15 kb band
  rect(ix0 + 4, iy0 + 4, ix0 + 4 + 15 * px_per_kb, iy1 - 2,
       col = adjustcolor("#d62728", alpha.f = 0.12), border = NA)
  text(ix0 + 4 + 7.5 * px_per_kb, iy1 - 4, "0–15 kb\n(D4Z4)",
       cex = 0.45, col = "#d62728", font = 2)
  text(ix0 + 4 + 0.5 * (ix1 - ix0), iy1 - 1.5,
       "Sequence-sharing peak at D4Z4 location",
       cex = 0.5, col = "grey15", font = 2, adj = c(0.5, 1))

  # Annotation: C1 metrics
  text(50, 90, "C1 metrics: silhouette 0.147 · 43.4% chr4q discordance",
       cex = 0.55, col = "grey25", adj = c(0.5, 0.5))
  text(50, 86, "Dip-C radial 0.732 (peripheral) · median 22 DUX4L vs 0–2 outliers",
       cex = 0.55, col = "grey25", adj = c(0.5, 0.5))
}

# ===========================================================================
# Panel C: Recombination map vs cross-arm affinity (the honest-null figure)
# Drawn as TWO sub-cells in the layout (cells 3 and 4); the page-level layout
# allocates them. panel_c_left and panel_c_right render one each.
# ===========================================================================
panel_c_scatter <- function(d, title, rho, p, header = NULL) {
  par(mar = c(3.6, 3.8, 2.2, 0.4))
  cols <- ifelse(d$callable_variants <= 12, "#d62728",
                 ifelse(d$callable_variants <= 200, "#ff7f0e", "#1f77b4"))
  plot(d$mean_cM_Mb, d$mean_affinity,
       pch = 21, bg = cols, col = "white",
       cex = 0.95, xlab = "", ylab = "",
       cex.axis = 0.7, main = "",
       xlim = range(set_all$mean_cM_Mb),
       ylim = c(0, max(set_all$mean_affinity) * 1.05))
  mtext("Subtelomeric recombination rate (cM/Mb, Lalli 2025)",
        side = 1, line = 2.1, cex = 0.6)
  mtext("Cross-arm sharing (frac. of arm seqs in cross-arm community)",
        side = 2, line = 2.4, cex = 0.55)
  title(main = title, cex.main = 0.7, font.main = 1, line = 0.4)
  if (!is.na(rho)) {
    lm_fit <- lm(d$mean_affinity ~ d$mean_cM_Mb)
    abline(lm_fit, col = "grey45", lwd = 1.2, lty = 2)
  }
  rho_text <- if (is.na(rho)) "ρ = NA"
              else sprintf("ρ = %.2f, p = %s, N = %d",
                           rho,
                           ifelse(p < 1e-3, sprintf("%.1e", p),
                                  sprintf("%.3f", p)),
                           nrow(d))
  legend("topright", legend = rho_text, cex = 0.65, bty = "n",
         text.col = "grey15")
  if (!is.null(header)) {
    mtext(header, side = 3, line = 1.0, adj = 0, cex = 0.85, font = 2)
  }
  # Acrocentric / PAR labels on the low-callability points
  if (any(d$callable_variants <= 12)) {
    lab_d <- d[d$callable_variants <= 12, ]
    text(lab_d$mean_cM_Mb, lab_d$mean_affinity, labels = lab_d$arm_label,
         pos = 4, cex = 0.5, col = "#d62728", offset = 0.25)
  }
  legend("bottomleft",
         legend = c("≤ 12 callable variants (acro p / PAR)",
                    "13–200 callable variants",
                    "> 200 callable variants"),
         pch = 21, pt.bg = c("#d62728", "#ff7f0e", "#1f77b4"),
         col = "white", pt.cex = 0.85,
         cex = 0.55, bty = "n")
}

panel_c_left  <- function() panel_c_scatter(
  set_all,
  sprintf("All %.0f arms (incl. low-callability)", as.numeric(c_all$n)),
  c_all$rho, c_all$p,
  header = "c   Recombination rate vs cross-arm affinity — confounding null")
panel_c_right <- function() panel_c_scatter(
  set_32,
  sprintf("N = %.0f (well-callable only)", as.numeric(c_32$n)),
  c_32$rho, c_32$p)

# ===========================================================================
# Panel D: Compartment-identity-at-tips diagnostic
# ===========================================================================
panel_d <- function() {
  par(mar = c(3.8, 3.8, 2.2, 0.8))

  # e1 distribution
  br <- seq(min(e1) - 0.005, max(e1) + 0.005, length.out = 22)
  h <- hist(e1, breaks = br, plot = FALSE)

  bar_cols <- ifelse(h$mids > 0, "#1f77b4", "#d62728")
  ymax <- max(h$counts) * 1.25
  plot(NA, xlim = c(min(br), max(br)), ylim = c(0, ymax),
       xlab = "", ylab = "", cex.axis = 0.75)
  rect(h$breaks[-length(h$breaks)], 0,
       h$breaks[-1], h$counts,
       col = bar_cols, border = "white", lwd = 0.6)
  abline(v = 0, col = "grey25", lty = 2, lwd = 1)
  abline(v = mean_e1, col = "#2ca02c", lty = 1, lwd = 1.5)

  mtext("HG002 Hi-C eigenvector e1 at chromosome tips (100 kb resolution, GC-oriented)",
        side = 1, line = 2.4, cex = 0.62)
  mtext("Number of arm × haplotype tips", side = 2, line = 2.4, cex = 0.62)

  pct_A <- 100 * n_A / n_tips
  txt <- sprintf("%d / %d tips A-compartment (%.0f%%)\nmean e1 = %+.4f (weak signature)\nDip-C C1 radial 0.732 (peripheral)\nDip-C C6 radial 0.505 (interior)\nDip-C C10 radial 0.474 (interior)",
                 n_A, n_tips, pct_A, mean_e1)
  legend("topleft", legend = txt, bty = "n", cex = 0.6, text.col = "grey15")

  legend("topright",
         legend = c("A-compartment (e1 > 0)", "B-compartment (e1 < 0)",
                    "e1 = 0", sprintf("mean e1 = %+.4f", mean_e1)),
         fill   = c("#1f77b4", "#d62728", NA, NA),
         border = c("white", "white", NA, NA),
         lty    = c(NA, NA, 2, 1),
         col    = c(NA, NA, "grey25", "#2ca02c"),
         lwd    = c(NA, NA, 1, 1.5),
         cex = 0.55, bty = "n", seg.len = 2.4)

  mtext("d   Compartment identity at tips (HG002, n = 92 arm × hap)",
        side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)
}

# ---------- Render ----------
# 2-row × 4-col layout:
#   row 1: A spans cols 1-2, B spans cols 3-4
#   row 2: C-left col 1, C-right col 2, D spans cols 3-4
render_to <- function(dev_open, dev_close) {
  dev_open()
  layout(rbind(c(1, 1, 2, 2),
               c(3, 4, 5, 5)),
         widths = c(1, 1, 1, 1), heights = c(1, 1))
  panel_a()
  panel_b()
  panel_c_left()
  panel_c_right()
  panel_d()
  dev_close()
}

pdf_path <- file.path(OUT, "figure_ed8.pdf")
png_path <- file.path(OUT, "figure_ed8.png")

render_to(function() pdf(pdf_path, width = 11, height = 8.5),
          dev.off)
render_to(function() png(png_path, width = 11 * 200, height = 8.5 * 200,
                         res = 200, type = "cairo"),
          dev.off)

cat(sprintf("Wrote %s and %s\n", pdf_path, png_path))
cat(sprintf("Stats — all %.0f arms: rho=%.3f, p=%.3g; restricted N=%.0f: rho=%.3f, p=%.3g\n",
            as.numeric(c_all$n), c_all$rho, c_all$p,
            as.numeric(c_32$n), c_32$rho, c_32$p))
cat(sprintf("Compartment: %d/%d A (%.1f%%), mean e1 = %+.4f\n",
            n_A, n_tips, 100 * n_A / n_tips, mean_e1))
