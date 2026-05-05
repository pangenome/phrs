#!/usr/bin/env Rscript
# Extended Data Figure 5 — Multi-resolution + confound robustness for Hi-C
# Inputs: see paper_prep/figures/ed5/sources.tsv
# Output: figure_ed5.pdf, figure_ed5.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(grid)
})

ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human"
OUT  <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "paper_prep/figures/ed5"

# ---------- Configuration ----------
RES_LIST <- c("5000bp", "10000bp", "20000bp", "50000bp", "100000bp")
RES_KB   <- c(5, 10, 20, 50, 100)

HPRC_SAMPLES <- c("chm13", "hg002", "hg002_cifi", "hg002_porec",
                  "hg00658", "hg02148", "hg02559", "na19036")
HPRC_LABEL <- c(chm13 = "CHM13",
                hg002 = "HG002 Hi-C",
                hg002_cifi = "HG002 CiFi",
                hg002_porec = "HG002 Pore-C",
                hg00658 = "HG00658",
                hg02148 = "HG02148",
                hg02559 = "HG02559",
                na19036 = "NA19036")

NO_ACRO_SAMPLES <- c("chm13", "hg002", "hg002_porec",
                     "hg00658", "hg02148", "hg02559", "na19036")  # no cifi

RPE_SAMPLES <- c("rpe1_async_cifi", "rpe1_async_porec", "rpe1_mitotic_cifi")
RPE_LABEL <- c(rpe1_async_cifi  = "RPE-1 Async CiFi",
               rpe1_async_porec = "RPE-1 Async Pore-C",
               rpe1_mitotic_cifi = "RPE-1 Mitotic CiFi")

# Sample colors (Tableau-ish, distinct)
sample_colors <- c(
  "CHM13"          = "#4e79a7",
  "HG002 Hi-C"     = "#f28e2b",
  "HG002 CiFi"     = "#e15759",
  "HG002 Pore-C"   = "#76b7b2",
  "HG00658"        = "#59a14f",
  "HG02148"        = "#edc948",
  "HG02559"        = "#b07aa1",
  "NA19036"        = "#9c755f"
)

# ---------- Loaders ----------
read_global <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d
}

extract_global <- function(d, what = c("bw", "mantel")) {
  what <- match.arg(what)
  if (is.null(d)) return(NA_real_)
  if (what == "bw") {
    r <- d[d$test == "within_vs_between", ]
    if (nrow(r) == 0) return(NA_real_)
    wm <- as.numeric(r$within_mean[1]); bm <- as.numeric(r$between_mean[1])
    if (is.na(wm) || is.na(bm) || bm <= 0) return(NA_real_)
    wm / bm
  } else {
    r <- d[d$test == "mantel", ]
    if (nrow(r) == 0) return(NA_real_)
    as.numeric(r$U_statistic[1])  # ρ stored in U_statistic column
  }
}

# ---------- Build ED5a/ED5b matrices ----------
build_global_grid <- function(samples, base_dir) {
  out <- expand.grid(sample = samples, res = RES_LIST,
                     stringsAsFactors = FALSE)
  out$bw     <- NA_real_
  out$mantel <- NA_real_
  for (i in seq_len(nrow(out))) {
    p <- file.path(base_dir, out$res[i],
                   paste0(out$sample[i], "_global_test.tsv"))
    g <- read_global(p)
    out$bw[i]     <- extract_global(g, "bw")
    out$mantel[i] <- extract_global(g, "mantel")
  }
  out$res_kb <- as.integer(sub("bp", "", out$res)) / 1000L
  out
}

cat("Loading community_based grids...\n")
cb <- build_global_grid(HPRC_SAMPLES, file.path(ROOT, "community_based"))
cat("Loading no_acrocentric grids...\n")
na <- build_global_grid(NO_ACRO_SAMPLES, file.path(ROOT, "no_acrocentric"))

# ---------- ED5c: O/E within vs between ----------
parse_arms <- function(s) strsplit(s, ";")[[1]]

oe_within_between <- function(sample, base_dir = file.path(ROOT, "community_based"),
                              res = "50000bp") {
  oe_path  <- file.path(base_dir, res, paste0(sample, "_oe_matrix.tsv"))
  bs_path  <- file.path(base_dir, res, paste0(sample, "_community_bootstrap_tests.tsv"))
  if (!file.exists(oe_path) || !file.exists(bs_path)) return(NULL)
  oe <- as.matrix(read.table(oe_path, sep = "\t", header = TRUE,
                             row.names = 1, check.names = FALSE))
  bs <- read_tsv(bs_path, show_col_types = FALSE)
  bs <- bs[grepl("^C", bs$community) & bs$type == "multi-arm", ]
  arm_to_comm <- list()
  for (i in seq_len(nrow(bs))) {
    for (a in parse_arms(bs$arms[i])) arm_to_comm[[a]] <- bs$community[i]
  }
  arms <- intersect(rownames(oe), names(arm_to_comm))
  if (length(arms) < 4) return(NULL)
  M <- oe[arms, arms]
  comm <- unlist(arm_to_comm[arms])
  # Same-chromosome (cis) excluded: identify chromosome from arm name
  chrom <- sub("_.*", "", arms)
  within  <- c(); between <- c()
  n <- length(arms)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (chrom[i] == chrom[j]) next  # exclude cis
      v <- M[i, j]
      if (is.na(v)) next
      if (comm[i] == comm[j]) within <- c(within, v) else between <- c(between, v)
    }
  }
  list(within = within, between = between)
}

cat("Computing O/E within vs between for ED5c...\n")
oe_data <- list()
for (s in HPRC_SAMPLES) {
  res <- oe_within_between(s)
  if (!is.null(res)) oe_data[[s]] <- res
}

# ---------- ED5d: per-community reproducibility ----------
build_comm_grid <- function() {
  comms <- paste0("C", 1:15)
  datasets <- c(HPRC_SAMPLES, RPE_SAMPLES)
  enrichment <- matrix(NA_real_, nrow = length(comms), ncol = length(datasets),
                       dimnames = list(comms, datasets))
  qval <- matrix(NA_real_, nrow = length(comms), ncol = length(datasets),
                 dimnames = list(comms, datasets))
  for (s in HPRC_SAMPLES) {
    p <- file.path(ROOT, "community_based", "50000bp",
                   paste0(s, "_community_bootstrap_tests.tsv"))
    if (!file.exists(p)) next
    d <- read_tsv(p, show_col_types = FALSE)
    d <- d[d$type == "multi-arm" & d$community %in% comms, ]
    for (i in seq_len(nrow(d))) {
      ratio <- d$observed_mean_contact[i] / d$random_mean[i]
      enrichment[d$community[i], s] <- ratio
      qval[d$community[i], s] <- d$p_adjusted_bh[i]
    }
  }
  for (s in RPE_SAMPLES) {
    p <- file.path(ROOT, "community_based", "RPE1", "50000bp",
                   paste0(s, "_community_bootstrap_tests.tsv"))
    if (!file.exists(p)) next
    d <- read_tsv(p, show_col_types = FALSE)
    d <- d[d$type == "multi-arm" & d$community %in% comms, ]
    for (i in seq_len(nrow(d))) {
      ratio <- d$observed_mean_contact[i] / d$random_mean[i]
      enrichment[d$community[i], s] <- ratio
      qval[d$community[i], s] <- d$p_adjusted_bh[i]
    }
  }
  list(enrichment = enrichment, qval = qval, datasets = datasets)
}

cat("Building ED5d community grid...\n")
comm_grid <- build_comm_grid()

# ---------- Plotting helpers ----------
panel_label <- function(lbl, x = 0.01, y = 0.97, cex = 0.8) {
  mtext(lbl, side = 3, line = 0.6, adj = 0, cex = cex, font = 2)
}

# ---------- Panel A: B/W ratio across resolutions ----------
panel_a <- function() {
  par(mar = c(3.6, 3.4, 2.2, 0.6))
  ymin <- 1; ymax <- max(cb$bw, na.rm = TRUE) * 1.2
  plot(NA, xlim = range(RES_KB), ylim = c(ymin, ymax),
       log = "xy", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "l")
  axis(1, at = RES_KB, labels = RES_KB, cex.axis = 0.7)
  ats <- c(1, 3, 10, 30, 100, 300)
  axis(2, at = ats, labels = ats, las = 1, cex.axis = 0.7)
  mtext("Hi-C resolution (kb)", side = 1, line = 2.1, cex = 0.75)
  mtext("Within / Between contact", side = 2, line = 2.2, cex = 0.75)
  abline(h = 1, col = "grey60", lty = 3)
  for (s in HPRC_SAMPLES) {
    sub <- cb[cb$sample == s & !is.na(cb$bw), ]
    sub <- sub[order(sub$res_kb), ]
    if (nrow(sub) == 0) next
    col <- sample_colors[HPRC_LABEL[s]]
    lines(sub$res_kb, sub$bw, col = col, lwd = 1.6)
    points(sub$res_kb, sub$bw, col = col, bg = col, pch = 21, cex = 0.7)
  }
  legend("topright", legend = HPRC_LABEL[HPRC_SAMPLES],
         col = sample_colors[HPRC_LABEL[HPRC_SAMPLES]],
         lwd = 1.6, pch = 21, pt.bg = sample_colors[HPRC_LABEL[HPRC_SAMPLES]],
         pt.cex = 0.7, cex = 0.55, bty = "n", ncol = 2)
  panel_label("a  Within/Between contact across 5 mcool resolutions")
}

# ---------- Panel B: Mantel ρ before/after acro+sex exclusion ----------
panel_b <- function() {
  par(mar = c(3.6, 3.4, 2.2, 0.6))
  ord <- NO_ACRO_SAMPLES
  full   <- cb[cb$sample %in% ord & cb$res == "50000bp", c("sample", "mantel")]
  excl   <- na[na$sample %in% ord & na$res == "50000bp", c("sample", "mantel")]
  m <- merge(full, excl, by = "sample", suffixes = c("_full", "_noacro"))
  m <- m[match(ord, m$sample), ]
  ymin <- 0; ymax <- 1.0
  plot(NA, xlim = c(ymin, ymax), ylim = c(ymin, ymax),
       xlab = "", ylab = "", bty = "l", xaxs = "i", yaxs = "i",
       cex.axis = 0.7, las = 1)
  abline(0, 1, col = "grey60", lty = 3)
  for (i in seq_len(nrow(m))) {
    s <- m$sample[i]
    col <- sample_colors[HPRC_LABEL[s]]
    points(m$mantel_full[i], m$mantel_noacro[i], pch = 21,
           col = "grey20", bg = col, cex = 1.5, lwd = 0.7)
  }
  mtext("Mantel ρ — full arm set", side = 1, line = 2.1, cex = 0.75)
  mtext("Mantel ρ — no acrocentric+sex", side = 2, line = 2.2, cex = 0.75)
  # Annotate two notable cases (CHM13, HG02148)
  for (i in seq_len(nrow(m))) {
    s <- m$sample[i]
    if (s %in% c("hg02148", "chm13")) {
      lab <- HPRC_LABEL[s]
      text(m$mantel_full[i] + 0.02, m$mantel_noacro[i] - 0.04,
           lab, cex = 0.55, adj = c(0, 0.5), col = "grey20")
    }
  }
  panel_label("b  Mantel ρ — full vs no acro/sex (50 kb)")
}

# ---------- Panel C: O/E within vs between, 8 datasets ----------
panel_c <- function() {
  par(mar = c(4.8, 3.4, 2.2, 0.6))
  samples <- names(oe_data)
  # Mean O/E within vs between (both include zero pairs — sparse-aware)
  mean_w <- sapply(samples, function(s) mean(oe_data[[s]]$within))
  mean_b <- sapply(samples, function(s) mean(oe_data[[s]]$between))
  ratios <- mean_w / mean_b
  ng <- length(samples)
  xpos_w <- seq_len(ng) - 0.20
  xpos_b <- seq_len(ng) + 0.20
  ymax <- max(c(mean_w, mean_b), na.rm = TRUE) * 1.55
  ymin <- min(c(mean_w, mean_b)[c(mean_w, mean_b) > 0], na.rm = TRUE) * 0.7
  plot(NA, xlim = c(0.4, ng + 0.6), ylim = c(ymin, ymax), log = "y",
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "l")
  axis(1, at = seq_len(ng), labels = HPRC_LABEL[samples], las = 2, cex.axis = 0.6)
  yat <- c(0.001, 0.003, 0.01, 0.03, 0.1)
  axis(2, at = yat, labels = c("0.001", "0.003", "0.01", "0.03", "0.1"),
       las = 1, cex.axis = 0.7)
  mtext("Mean O/E contact (log)", side = 2, line = 2.2, cex = 0.75)
  bw <- 0.34
  for (i in seq_along(samples)) {
    rect(xpos_w[i] - bw/2, ymin, xpos_w[i] + bw/2, mean_w[i],
         col = adjustcolor("#1f77b4", 0.85), border = "#0d3a66", lwd = 0.6)
    rect(xpos_b[i] - bw/2, ymin, xpos_b[i] + bw/2, mean_b[i],
         col = adjustcolor("#999999", 0.85), border = "#555555", lwd = 0.6)
    text(i, ymax * 0.85, sprintf("%.1f×", ratios[i]),
         cex = 0.6, col = "grey20")
  }
  legend("bottomleft", legend = c("Within community", "Between"),
         fill = c(adjustcolor("#1f77b4", 0.85), adjustcolor("#999999", 0.85)),
         border = c("#0d3a66", "#555555"),
         bty = "n", cex = 0.6)
  panel_label("c  O/E contact (mean): within vs between sequence community")
}

# ---------- Panel D: Per-community reproducibility heatmap ----------
panel_d <- function() {
  par(mar = c(5.2, 3.6, 2.2, 3.2))
  E <- comm_grid$enrichment
  Q <- comm_grid$qval
  # log2 enrichment, capped
  L <- log2(pmax(E, 1e-3))
  L[!is.finite(L)] <- NA
  cap <- 6
  L[L > cap]  <- cap
  L[L < -cap] <- -cap
  ds <- comm_grid$datasets
  ds_lab <- c(HPRC_LABEL[HPRC_SAMPLES], RPE_LABEL[RPE_SAMPLES])
  comms <- rownames(E)
  nr <- nrow(E); nc <- ncol(E)
  # Image coordinate convention: rows = y, cols = x
  plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       xaxs = "i", yaxs = "i")
  pal <- colorRampPalette(c("#053061", "#4393c3", "#f7f7f7",
                            "#d6604d", "#67001f"))(101)
  zlim <- c(-cap, cap)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      v <- L[i, j]
      if (is.na(v)) {
        col <- "grey85"
      } else {
        idx <- round((v - zlim[1]) / (zlim[2] - zlim[1]) * 100) + 1
        idx <- max(1, min(101, idx))
        col <- pal[idx]
      }
      rect(j - 0.5, nr - i + 0.5, j + 0.5, nr - i + 1.5,
           col = col, border = "white", lwd = 0.4)
      q <- Q[i, j]
      if (!is.na(q)) {
        if (q < 0.001) sym <- "**"
        else if (q < 0.05) sym <- "*"
        else sym <- ""
        if (sym != "") {
          text(j, nr - i + 1, sym, cex = 0.55,
               col = ifelse(abs(v) > 3, "white", "black"))
        }
      }
    }
  }
  axis(1, at = seq_len(nc), labels = ds_lab, las = 2, cex.axis = 0.55, tick = FALSE)
  axis(2, at = seq(nr, 1), labels = comms, las = 1, cex.axis = 0.6, tick = FALSE)
  # Color bar
  usr <- par("usr")
  bar_x <- usr[2] + 0.4
  bar_w <- 0.25
  bar_y0 <- 1; bar_y1 <- nr
  ny <- 60
  ysteps <- seq(bar_y0, bar_y1, length.out = ny + 1)
  vsteps <- seq(zlim[1], zlim[2], length.out = ny + 1)
  for (k in seq_len(ny)) {
    idx <- round((vsteps[k] - zlim[1]) / (zlim[2] - zlim[1]) * 100) + 1
    rect(bar_x, ysteps[k], bar_x + bar_w, ysteps[k + 1],
         col = pal[idx], border = NA, xpd = TRUE)
  }
  rect(bar_x, bar_y0, bar_x + bar_w, bar_y1, border = "grey40", lwd = 0.4, xpd = TRUE)
  bar_at_v <- seq(-cap, cap, by = 3)
  bar_at_y <- bar_y0 + (bar_at_v - zlim[1]) / (zlim[2] - zlim[1]) * (bar_y1 - bar_y0)
  for (k in seq_along(bar_at_v)) {
    text(bar_x + bar_w + 0.1, bar_at_y[k], bar_at_v[k],
         cex = 0.55, adj = 0, xpd = TRUE)
  }
  text(bar_x + bar_w/2, bar_y1 + 0.7, expression(log[2]~"obs/null"),
       cex = 0.6, xpd = TRUE)
  panel_label("d  Per-community enrichment (15 × 11; * q<0.05, ** q<0.001)")
}

# ---------- Render ----------
render <- function(device, path, w, h) {
  if (device == "pdf") {
    cairo_pdf(path, width = w, height = h)
  } else {
    png(path, width = w * 150, height = h * 150, res = 150, type = "cairo")
  }
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
         widths = c(1, 1), heights = c(1, 1.05))
  par(family = "sans", cex = 0.85)
  panel_a()
  panel_b()
  panel_c()
  panel_d()
  invisible(dev.off())
}

W <- 10.5; H <- 9.0
render("pdf", file.path(OUT, "figure_ed5.pdf"), W, H)
render("png", file.path(OUT, "figure_ed5.png"), W, H)

# ---------- Print summary stats for caption ----------
cat("\n=== Summary statistics ===\n")
cat("ED5a — W/B range across all (sample × resolution) cells:\n")
bw_all <- cb$bw[!is.na(cb$bw)]
cat(sprintf("  min=%.2f, median=%.2f, max=%.2f, n=%d\n",
            min(bw_all), median(bw_all), max(bw_all), length(bw_all)))

cat("ED5b — Mantel ρ (50 kb) full vs no-acro:\n")
m50 <- merge(cb[cb$res == "50000bp", c("sample", "mantel")],
             na[na$res == "50000bp", c("sample", "mantel")],
             by = "sample", suffixes = c("_full", "_noacro"))
print(m50)

cat("ED5c — within/between O/E ratios (means, zeros included):\n")
for (s in names(oe_data)) {
  w <- oe_data[[s]]$within; b <- oe_data[[s]]$between
  cat(sprintf("  %s: within mean=%.4f (n=%d), between mean=%.4f (n=%d), ratio=%.1f×\n",
              s, mean(w), length(w), mean(b), length(b),
              mean(w)/mean(b)))
}

cat("ED5d — significant communities (q<0.05) per dataset:\n")
sig_counts <- colSums(comm_grid$qval < 0.05, na.rm = TRUE)
nz_counts  <- colSums(!is.na(comm_grid$qval))
for (k in seq_along(sig_counts)) {
  cat(sprintf("  %s: %d/%d communities significant\n",
              names(sig_counts)[k], sig_counts[k], nz_counts[k]))
}

cat("Done.\n")
