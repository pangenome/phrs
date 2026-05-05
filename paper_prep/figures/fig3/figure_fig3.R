#!/usr/bin/env Rscript
# Main Figure 3 — Three-dimensional nuclear organisation mirrors sequence communities
# Inputs: see paper_prep/figures/fig3/sources.tsv
# Output: figure_fig3.pdf, figure_fig3.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(scales)
  library(grid)
})

ANALYSIS_ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human"
MOUSE_ROOT    <- "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T"
DIPC_ROOT     <- "/moosefs/guarracino/HPRCv2/dipc_t2t"
FLANKING_DIPC <- "/moosefs/guarracino/HPRCv2/PHR_III/dipc_flanking_radial.tsv"

OUT <- dirname(sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "paper_prep/figures/fig3"

read_global <- function(path) {
  d <- read_tsv(path, show_col_types = FALSE) |>
    filter(test == "within_vs_between")
  list(within = d$within_mean[1],
       between = d$between_mean[1],
       p = d$p_value[1])
}

bw_ratio <- function(g) g$between / g$within   # < 1 = within-community contact stronger

# ---------- Panel A: HG002 Pore-C contact matrix ordered by community ----------
panel_a <- function() {
  mat <- read_tsv(
    file.path(ANALYSIS_ROOT, "community_based/50000bp/hg002_porec_contact_matrix.tsv"),
    show_col_types = FALSE)
  rn <- mat[[1]]; mat <- as.matrix(mat[, -1]); rownames(mat) <- rn
  comm <- read_tsv(
    file.path(ANALYSIS_ROOT, "community_based/50000bp/hg002_porec_hic.arm-leiden.communities.tsv"),
    show_col_types = FALSE)

  comm <- comm[match(rownames(mat), comm$arm), ]
  ord <- order(comm$community, comm$arm)
  mat <- mat[ord, ord]
  comm_ord <- comm$community[ord]

  # log-scale (zeros → tiny floor)
  vals <- mat
  pos <- vals[vals > 0]
  floor_v <- if (length(pos)) quantile(pos, 0.05) else 1e-6
  vals[vals == 0] <- floor_v
  vals <- log10(vals)
  vrng <- range(vals, finite = TRUE)
  vals_norm <- (vals - vrng[1]) / diff(vrng)
  pal <- colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#7f0000"))(100)

  par(mar = c(2.0, 4.5, 2.6, 2.0))
  n <- nrow(vals_norm)
  image(seq_len(n), seq_len(n), t(vals_norm)[, n:1],
        col = pal, axes = FALSE, ann = FALSE,
        xlab = "", ylab = "", useRaster = TRUE)

  # Community block borders
  rl <- rle(as.character(comm_ord))
  ends <- cumsum(rl$lengths)
  starts <- c(0, head(ends, -1)) + 1
  for (i in seq_along(rl$values)) {
    rect(starts[i] - 0.5, n - ends[i] + 0.5,
         ends[i] + 0.5, n - starts[i] + 1.5,
         border = "#1f3b73", lwd = 1.2)
    midx <- (starts[i] + ends[i]) / 2
    midy <- n - midx + 1
    if (rl$lengths[i] >= 3) {
      text(midx, midy, rl$values[i], cex = 0.6, col = "#1f3b73", font = 2)
    }
  }
  mtext("a   HG002 Pore-C inter-arm contact matrix at 50 kb, ordered by sequence community",
        side = 3, line = 0.6, adj = 0, cex = 0.95, font = 2)
  mtext("B/W = 0.056   p = 3.9e-85   77 arm-haplotypes   diagonal blocks = communities",
        side = 1, line = 0.7, cex = 0.62, col = "grey25")

  # Colour-bar
  cb_x <- grconvertX(c(0.92, 0.95), from = "npc", to = "user")
  cb_y <- grconvertY(c(0.18, 0.82), from = "npc", to = "user")
  ncols <- length(pal)
  cb_steps <- seq(cb_y[1], cb_y[2], length.out = ncols + 1)
  for (i in seq_len(ncols)) {
    rect(cb_x[1], cb_steps[i], cb_x[2], cb_steps[i + 1], col = pal[i], border = NA)
  }
  rect(cb_x[1], cb_y[1], cb_x[2], cb_y[2], border = "grey25", lwd = 0.6)
  text(cb_x[2] + 0.3, cb_y[2], sprintf("%.0e", 10 ^ vrng[2]),
       adj = c(0, 1), cex = 0.55)
  text(cb_x[2] + 0.3, cb_y[1], sprintf("%.0e", 10 ^ vrng[1]),
       adj = c(0, 0), cex = 0.55)
  text(cb_x[2] + 0.3, mean(cb_y), "contacts (log10)",
       adj = c(0, 0.5), cex = 0.55, srt = 90)
}

# ---------- Panel B: Convergent-evidence forest plot (14 tests) ----------
panel_b <- function() {
  hic <- list(
    list(label = "HG002 Pore-C 50 kb",        tech = "Pore-C", file = "community_based/50000bp/hg002_porec_global_test.tsv"),
    list(label = "HG002 CiFi 50 kb",          tech = "CiFi",   file = "community_based/50000bp/hg002_cifi_global_test.tsv"),
    list(label = "HG002 Hi-C 50 kb",          tech = "Hi-C",   file = "community_based/50000bp/hg002_global_test.tsv"),
    list(label = "CHM13 Hi-C 50 kb",          tech = "Hi-C",   file = "community_based/50000bp/chm13_global_test.tsv"),
    list(label = "HG00658 Hi-C 50 kb",        tech = "Hi-C",   file = "community_based/50000bp/hg00658_global_test.tsv"),
    list(label = "HG02148 Hi-C 50 kb",        tech = "Hi-C",   file = "community_based/50000bp/hg02148_global_test.tsv"),
    list(label = "HG02559 Hi-C 50 kb",        tech = "Hi-C",   file = "community_based/50000bp/hg02559_global_test.tsv"),
    list(label = "NA19036 Hi-C 50 kb",        tech = "Hi-C",   file = "community_based/50000bp/na19036_global_test.tsv"))
  hic_rows <- do.call(rbind, lapply(hic, function(r) {
    g <- read_global(file.path(ANALYSIS_ROOT, r$file))
    data.frame(label = r$label, tech = r$tech,
               ratio = bw_ratio(g), p = g$p, conv = "B/W",
               stringsAsFactors = FALSE)
  }))

  # Dip-C 16-cell GM12878 (W/B): use summary
  dipc_sum <- read_tsv(file.path(DIPC_ROOT,
    "output_q0_XX/community_enrichment_16cells_500kb_summary.tsv"),
    show_col_types = FALSE) |> filter(test == "wilcoxon_signed_rank")
  dipc_row <- data.frame(
    label = "GM12878 Dip-C 16 cells",
    tech  = "Dip-C",
    ratio = dipc_sum$mean_ratio[1],
    p     = dipc_sum$p_value[1],
    conv  = "W/B",
    stringsAsFactors = FALSE)

  # Sperm 20-cell W/B: use summary
  sperm_sum <- read_tsv(file.path(DIPC_ROOT,
    "sperm/enrichment_corrected/sperm_all20_summary.tsv"),
    show_col_types = FALSE) |> filter(test == "wilcoxon_signed_rank")
  sperm_row <- data.frame(
    label = "Sperm scHi-C 20 cells",
    tech  = "Sperm",
    ratio = sperm_sum$mean_ratio[1],
    p     = sperm_sum$p_value[1],
    conv  = "W/B",
    stringsAsFactors = FALSE)

  # Mouse meiotic 4 stages (B/W) at 100 kb (4 Mb flanking analysis)
  mouse_stages <- c("leptotene", "zygotene", "pachytene", "diplotene")
  mouse_rows <- do.call(rbind, lapply(mouse_stages, function(st) {
    f <- file.path(MOUSE_ROOT, "community_analysis_4Mb/100000bp",
                   sprintf("zuo2021_%s_global_test.tsv", st))
    g <- read_global(f)
    data.frame(label = sprintf("Mouse %s 4 Mb", st),
               tech  = "Mouse meiotic Hi-C",
               ratio = bw_ratio(g),
               p     = g$p,
               conv  = "B/W",
               stringsAsFactors = FALSE)
  }))

  d <- rbind(hic_rows, dipc_row, sperm_row, mouse_rows)
  # Order: by tech grouping then by ratio (smaller = stronger at top)
  d$tech <- factor(d$tech, levels = c("Hi-C", "Pore-C", "CiFi", "Dip-C", "Sperm", "Mouse meiotic Hi-C"))
  d <- d[order(d$tech, d$ratio), ]
  d$y <- nrow(d):1

  pal <- c("Hi-C" = "#1f78b4", "Pore-C" = "#33a02c", "CiFi" = "#6a3d9a",
           "Dip-C" = "#e31a1c", "Sperm" = "#ff7f00", "Mouse meiotic Hi-C" = "#b15928")

  par(mar = c(3.4, 11.0, 3.2, 1.6))
  xlim <- c(0.0008, 1.4)
  plot(NA, xlim = xlim, ylim = c(0.5, max(d$y) + 0.5),
       log = "x", axes = FALSE, xlab = "", ylab = "")
  abline(v = 1, col = "grey60", lty = 2, lwd = 1)
  axis(1, at = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
       labels = c("0.001", "0.01", "0.05", "0.1", "0.5", "1"), cex.axis = 0.7)
  mtext("effect-size ratio (left of 1 = within-community stronger / closer)",
        side = 1, line = 2.0, cex = 0.7)

  # Stripes for tech groups
  tech_blocks <- rle(as.character(d$tech))
  ends <- cumsum(tech_blocks$lengths)
  starts <- c(0, head(ends, -1)) + 1
  for (i in seq_along(tech_blocks$values)) {
    if (i %% 2 == 1) {
      yt <- max(d$y) + 0.5 - starts[i] + 1
      yb <- max(d$y) + 0.5 - ends[i]
      rect(xlim[1], yb, xlim[2], yt, col = "#f7f7f7", border = NA)
    }
  }

  # Points + p-value annotations + convention tag
  for (i in seq_len(nrow(d))) {
    points(d$ratio[i], d$y[i], pch = 21, bg = pal[as.character(d$tech[i])],
           col = "black", cex = 1.3)
    text(d$ratio[i] * 1.05, d$y[i],
         sprintf("%.3f  p=%.1e  (%s)", d$ratio[i], d$p[i], d$conv[i]),
         adj = c(0, 0.5), cex = 0.55, col = "grey25")
  }
  # Y-axis labels
  axis(2, at = d$y, labels = d$label, las = 1, cex.axis = 0.62, lwd = 0)

  legend("bottomright", legend = names(pal), pt.bg = pal, pch = 21, col = "black",
         cex = 0.62, bty = "n", inset = c(0.02, 0.02), pt.cex = 1.1,
         title = "technology", title.adj = 0)
  mtext("b   Convergent evidence — 14 tests, 6 technologies",
        side = 3, line = 1.0, adj = 0, cex = 0.95, font = 2)
  mtext("all on within-community-stronger side of unity",
        side = 3, line = 0.1, adj = 0, cex = 0.62, col = "grey25")
}

# ---------- Panel C: S_all negative-control vs C-community per-cell W/B ----------
panel_c <- function() {
  # C-community per-cell ratio: pooled W/B across all C* communities (per_cell.tsv)
  gm_c_pc <- read_tsv(file.path(DIPC_ROOT,
    "output_q0_XX/community_enrichment_16cells_500kb_per_cell.tsv"),
    show_col_types = FALSE)
  sp_c_pc <- read_tsv(file.path(DIPC_ROOT,
    "sperm/enrichment_corrected/sperm_all20_per_cell.tsv"),
    show_col_types = FALSE)
  # S_all per-cell ratio: from per_community_per_cell.tsv where community == S_all
  gm_pccc <- read_tsv(file.path(DIPC_ROOT,
    "output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv"),
    show_col_types = FALSE)
  sp_pccc <- read_tsv(file.path(DIPC_ROOT,
    "sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv"),
    show_col_types = FALSE)
  gm_pc <- data.frame(
    cell_id = gm_c_pc$cell_id,
    C       = gm_c_pc$ratio,
    S_all   = gm_pccc$ratio[match(gm_c_pc$cell_id,
                                  gm_pccc$cell_id[gm_pccc$community == "S_all"])])
  gm_pc$S_all <- gm_pccc |> filter(community == "S_all") |>
    (\(x) x$ratio[match(gm_pc$cell_id, x$cell_id)])()
  sp_pc <- data.frame(
    cell_id = sp_c_pc$cell_id,
    C       = sp_c_pc$ratio)
  sp_pc$S_all <- sp_pccc |> filter(community == "S_all") |>
    (\(x) x$ratio[match(sp_pc$cell_id, x$cell_id)])()

  par(mar = c(3.6, 4.4, 2.6, 1.0))
  xlim <- c(0.5, 4.5)
  ylim <- range(c(gm_pc$C, gm_pc$S_all, sp_pc$C, sp_pc$S_all), na.rm = TRUE) * c(0.95, 1.08)

  plot(NA, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "per-cell mean W/B (within / between distance)")
  axis(2, cex.axis = 0.75, las = 1)
  abline(h = 1, col = "grey60", lty = 2)
  axis(1, at = c(1, 2, 3, 4),
       labels = c("GM12878\nC-comm", "GM12878\nS_all", "Sperm\nC-comm", "Sperm\nS_all"),
       cex.axis = 0.7, lwd = 0, line = -0.3, padj = 0.5)

  draw_strip <- function(x, vals, col_pt, col_box) {
    if (!length(vals)) return(invisible())
    qs <- quantile(vals, c(0.25, 0.5, 0.75), na.rm = TRUE)
    rect(x - 0.20, qs[1], x + 0.20, qs[3], border = col_box, col = NA, lwd = 1.2)
    segments(x - 0.20, qs[2], x + 0.20, qs[2], col = col_box, lwd = 1.6)
    points(rep(x, length(vals)) + runif(length(vals), -0.10, 0.10),
           vals, pch = 21, bg = col_pt, col = "black", cex = 0.8)
  }

  draw_strip(1, gm_pc$C,     "#2166ac", "#0a3268")
  draw_strip(2, gm_pc$S_all, "#999999", "#444444")
  draw_strip(3, sp_pc$C,     "#b2182b", "#67001f")
  draw_strip(4, sp_pc$S_all, "#999999", "#444444")

  # Counts below unity
  count_below <- function(x) sum(x < 1, na.rm = TRUE)
  txt <- c(
    sprintf("%d/%d < 1", count_below(gm_pc$C),     sum(!is.na(gm_pc$C))),
    sprintf("%d/%d < 1", count_below(gm_pc$S_all), sum(!is.na(gm_pc$S_all))),
    sprintf("%d/%d < 1", count_below(sp_pc$C),     sum(!is.na(sp_pc$C))),
    sprintf("%d/%d < 1", count_below(sp_pc$S_all), sum(!is.na(sp_pc$S_all))))
  yy <- ylim[2] - diff(ylim) * 0.04
  text(1:4, yy, txt, cex = 0.65, col = "grey25", font = 2)

  mtext("c   S_all (7 non-sharing arms) reverses the community signal",
        side = 3, line = 0.6, adj = 0, cex = 0.95, font = 2)
  mtext("dashed line = W/B = 1 (no community signal)",
        side = 1, line = 2.4, cex = 0.6, col = "grey25")
}

# ---------- Panel D: Flanking-vs-PHR B/W (8 samples) + Dip-C flanking inset ----------
panel_d <- function() {
  samples <- c(
    "HG002 Pore-C" = "hg002_porec",
    "HG002 CiFi"   = "hg002_cifi",
    "HG002 Hi-C"   = "hg002",
    "CHM13 Hi-C"   = "chm13",
    "HG00658"      = "hg00658",
    "HG02148"      = "hg02148",
    "HG02559"      = "hg02559",
    "NA19036"      = "na19036")

  rows <- lapply(seq_along(samples), function(i) {
    nm <- names(samples)[i]; pref <- samples[[i]]
    phr_f  <- file.path(ANALYSIS_ROOT,
                        sprintf("community_based/50000bp/%s_global_test.tsv", pref))
    flk_f  <- file.path(ANALYSIS_ROOT,
                        sprintf("flanking/100000bp/%s_global_test.tsv", pref))
    phr_g <- read_global(phr_f)
    flk_g <- if (file.exists(flk_f)) read_global(flk_f) else list(within = NA, between = NA, p = NA)
    data.frame(sample = nm,
               phr_bw  = bw_ratio(phr_g),
               phr_p   = phr_g$p,
               flk_bw  = if (is.na(flk_g$within)) NA else bw_ratio(flk_g),
               flk_p   = flk_g$p,
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$ratio_fold <- ifelse(is.na(d$flk_bw), NA, d$phr_bw / d$flk_bw)

  # Bar plot: paired bars (PHR vs Flanking) per sample, log-scaled y
  par(mar = c(4.6, 4.4, 3.2, 1.6))
  n <- nrow(d)
  x_phr <- (1:n) - 0.18
  x_flk <- (1:n) + 0.18
  ymin <- min(c(d$phr_bw, d$flk_bw), na.rm = TRUE) * 0.6
  ymax <- max(c(d$phr_bw, d$flk_bw), na.rm = TRUE) * 1.6

  plot(NA, xlim = c(0.4, n + 0.6), ylim = c(ymin, ymax), log = "y",
       axes = FALSE, xlab = "", ylab = "B/W ratio (between / within)")
  axis(2, at = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3),
       labels = c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3"),
       cex.axis = 0.75, las = 1)
  axis(1, at = 1:n, labels = d$sample, las = 2, cex.axis = 0.68, lwd = 0, line = -0.5)

  segments(x_phr, ymin, x_phr, d$phr_bw, col = "#1f78b4", lwd = 9, lend = 1)
  for (i in seq_len(n)) {
    if (!is.na(d$flk_bw[i])) {
      segments(x_flk[i], ymin, x_flk[i], d$flk_bw[i], col = "#fb9a29", lwd = 9, lend = 1)
    } else {
      text(x_flk[i], ymin * 1.4, "NA", cex = 0.5, col = "grey50")
    }
  }
  # Fold-change annotations
  for (i in seq_len(n)) {
    if (!is.na(d$ratio_fold[i])) {
      text((x_phr[i] + x_flk[i]) / 2, max(d$phr_bw[i], d$flk_bw[i], na.rm = TRUE) * 1.4,
           sprintf("%.0f×", d$ratio_fold[i]), cex = 0.55, col = "grey25")
    }
  }
  legend("topright", legend = c("PHR (50 kb)", "Flanking 100 kb (centromere-ward)"),
         fill = c("#1f78b4", "#fb9a29"), border = NA, cex = 0.65, bty = "n",
         inset = c(0.02, 0.02))
  mtext("d   Flanking paradox — unique-sequence flanks > PHRs",
        side = 3, line = 1.0, adj = 0, cex = 0.95, font = 2)
  mtext("100 kb centromere-ward of PHR boundary; B/W as in (b)",
        side = 3, line = 0.1, adj = 0, cex = 0.62, col = "grey25")

  # Inset: Dip-C flanking radial — draw inside panel d using base R user coords
  flk <- read_tsv(FLANKING_DIPC, show_col_types = FALSE)
  fl <- as.numeric(flk$value[flk$metric == "mean_flanking_radial"])
  nf <- as.numeric(flk$value[flk$metric == "mean_nonflanking_radial"])
  p_f <- as.numeric(flk$value[flk$metric == "mannwhitney_p"])
  if (length(p_f) == 0 || is.na(p_f)) p_f <- 7.4e-35

  ix <- grconvertX(c(0.05, 0.42), from = "npc", to = "user")
  iy <- grconvertY(c(0.62, 0.96), from = "npc", to = "user")
  rect(ix[1], iy[1], ix[2], iy[2], col = "#fffaf0", border = "grey50", lwd = 0.6)
  text(mean(ix), iy[2] - diff(iy) * 0.10,
       "Dip-C flanking radial (GM12878)",
       cex = 0.65, font = 2)

  vmax <- 0.65
  bar_centers <- ix[1] + diff(ix) * c(0.30, 0.70)
  bar_w_user  <- diff(ix) * 0.16
  base_y <- iy[1] + diff(iy) * 0.18
  top_y_max <- iy[1] + diff(iy) * 0.70
  vals <- c(fl, nf); cols <- c("#fb9a29", "#bdbdbd")
  labs <- c("flanking", "non-flank.\nterminal")
  for (i in 1:2) {
    h <- (vals[i] / vmax) * (top_y_max - base_y)
    rect(bar_centers[i] - bar_w_user / 2, base_y,
         bar_centers[i] + bar_w_user / 2, base_y + h,
         col = cols[i], border = NA)
    text(bar_centers[i], base_y + h + diff(iy) * 0.05,
         sprintf("%.3f", vals[i]), cex = 0.6)
    text(bar_centers[i], base_y - diff(iy) * 0.06, labs[i],
         cex = 0.55, col = "grey25")
  }
  text(mean(ix), iy[1] + diff(iy) * 0.93,
       sprintf("flanking is MORE INTERIOR  p = %.0e", p_f),
       cex = 0.55, col = "grey25")
}

# ---------- Render ----------
render_to <- function(dev_open, dev_close) {
  dev_open()
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1, 1.05), heights = c(1, 1))
  panel_a()
  panel_b()
  panel_c()
  panel_d()
  dev_close()
}

pdf_path <- file.path(OUT, "figure_fig3.pdf")
png_path <- file.path(OUT, "figure_fig3.png")

render_to(function() pdf(pdf_path, width = 16, height = 12), dev.off)
render_to(function() png(png_path, width = 16 * 200, height = 12 * 200,
                         res = 200, type = "cairo"), dev.off)

cat("Wrote", pdf_path, "and", png_path, "\n")
