#!/usr/bin/env Rscript

out_dir <- "slides/v2-review-zoom/_revision_assets/hic_visual_redesign"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

human_root <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human"
dipc_root <- "/moosefs/guarracino/HPRCv2/dipc_t2t"
mouse_root <- "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T"

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, ...)
}

extract_mantel <- function(path) {
  d <- read_tsv(path)
  as.numeric(d$U_statistic[d$test == "mantel"][1])
}

draw_note <- function(text, x0, y0, x1, y1, fill = "#fff7e6",
                      border = "#9a6b00", cex = 1.0,
                      font = 2, col = "#222222") {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = 1.1)
  text((x0 + x1) / 2, (y0 + y1) / 2, text, cex = cex, font = font,
       col = col)
}

render_10a <- function() {
  mat_path <- file.path(human_root, "community_based/50000bp",
                        "hg002_porec_contact_matrix.tsv")
  comm_path <- file.path(human_root, "community_based/50000bp",
                         "hg002_porec_hic.arm-leiden.communities.tsv")
  mat_df <- read_tsv(mat_path)
  rn <- mat_df[[1]]
  vals <- as.matrix(mat_df[, -1])
  rownames(vals) <- rn
  comm <- read_tsv(comm_path)
  comm <- comm[match(rownames(vals), comm$arm), ]
  ord <- order(comm$community, comm$arm)
  vals <- vals[ord, ord]
  comm_ord <- comm$community[ord]

  pos <- vals[vals > 0]
  floor_v <- if (length(pos)) stats::quantile(pos, 0.05) else 1e-6
  vals[vals == 0] <- floor_v
  log_vals <- log10(vals)
  zlim <- range(log_vals, finite = TRUE)
  norm <- (log_vals - zlim[1]) / diff(zlim)
  pal <- grDevices::colorRampPalette(
    c("white", "#fee0d2", "#fc9272", "#de2d26", "#7f0000"))(100)
  cols <- pal[pmax(1, pmin(100, floor(norm * 99) + 1))]
  dim(cols) <- dim(norm)

  png(file.path(out_dir, "slide_10a_square_matrix_candidate.png"),
      width = 1800, height = 1800, res = 180, type = "cairo")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(0.78, 0.22))
  par(mar = c(5.0, 5.0, 5.0, 0.8), pty = "s", family = "sans")
  n <- nrow(vals)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
       xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  rasterImage(as.raster(t(cols)[, n:1]), 0.5, 0.5, n + 0.5, n + 0.5,
              interpolate = FALSE)
  box(col = "#444444", lwd = 1.0)

  rl <- rle(as.character(comm_ord))
  ends <- cumsum(rl$lengths)
  starts <- c(0, head(ends, -1)) + 1
  for (i in seq_along(rl$values)) {
    rect(starts[i] - 0.5, n - ends[i] + 0.5,
         ends[i] + 0.5, n - starts[i] + 1.5,
         border = "#0b3c78", lwd = 1.25)
    if (rl$lengths[i] >= 3) {
      midx <- (starts[i] + ends[i]) / 2
      midy <- n - midx + 1
      text(midx, midy, rl$values[i], cex = 0.78,
           col = "#0b3c78", font = 2)
    }
  }

  axis(1, at = c(1, ceiling(n / 2), n),
       labels = c("first", "community ordered", "last"),
       cex.axis = 0.85, tick = FALSE)
  axis(2, at = c(1, ceiling(n / 2), n),
       labels = c("last", "", "first"),
       cex.axis = 0.85, tick = FALSE, las = 1)
  mtext("community ordered", side = 2, line = 3.0, cex = 0.80)
  title("HG002 Pore-C contact matrix - square 1:1 layout",
        cex.main = 1.15, font.main = 2)
  mtext("77 arm-haplotypes ordered by sequence community",
        side = 3, line = 1.0, cex = 0.80, col = "#333333")
  mtext("B/W = 0.056    p = 3.9e-85    diagonal blocks = communities",
        side = 1, line = 3.0, cex = 0.84, col = "#333333")

  par(mar = c(5.0, 1.0, 5.0, 4.0), pty = "m")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  y <- seq(0.10, 0.90, length.out = length(pal) + 1)
  for (i in seq_along(pal)) {
    rect(0.20, y[i], 0.34, y[i + 1], col = pal[i], border = NA)
  }
  rect(0.20, 0.10, 0.34, 0.90, border = "#555555", lwd = 0.8)
  text(0.38, 0.90, sprintf("%.0e", 10 ^ zlim[2]), adj = c(0, 0.5),
       cex = 0.70)
  text(0.38, 0.10, sprintf("%.0e", 10 ^ zlim[1]), adj = c(0, 0.5),
       cex = 0.70)
  text(0.73, 0.50, "contacts (log10)", srt = 90, cex = 0.72)
}

render_10b <- function() {
  samples <- c("chm13", "hg002", "hg002_porec",
               "hg00658", "hg02148", "hg02559", "na19036")
  labels <- c(chm13 = "CHM13",
              hg002 = "HG002 Hi-C",
              hg002_porec = "HG002 Pore-C",
              hg00658 = "HG00658",
              hg02148 = "HG02148",
              hg02559 = "HG02559",
              na19036 = "NA19036")
  cols <- c("#4e79a7", "#f28e2b", "#76b7b2", "#59a14f",
            "#edc948", "#b07aa1", "#9c755f")
  full <- vapply(samples, function(s) {
    extract_mantel(file.path(human_root, "community_based/50000bp",
                             paste0(s, "_global_test.tsv")))
  }, numeric(1))
  excl <- vapply(samples, function(s) {
    extract_mantel(file.path(human_root, "no_acrocentric/50000bp",
                             paste0(s, "_global_test.tsv")))
  }, numeric(1))

  png(file.path(out_dir, "slide_10b_mantel_exclusions_clarity.png"),
      width = 1800, height = 1200, res = 180, type = "cairo")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(5.6, 6.0, 5.0, 2.2), family = "sans")
  plot(NA, xlim = c(0, 0.92), ylim = c(0, 0.92), xaxs = "i", yaxs = "i",
       xlab = "Mantel rho, full arm set",
       ylab = "Mantel rho, excluding acrocentric + sex arms",
       main = "Mantel exclusions: points above x=y mean the signal gets cleaner",
       cex.main = 1.0, cex.lab = 0.95, cex.axis = 0.84, las = 1)
  grid(col = "#e6e6e6", lwd = 0.7)
  abline(0, 1, col = "#555555", lwd = 1.4, lty = 2)
  polygon(c(0, 0.92, 0.92), c(0, 0.92, 0.92),
          col = grDevices::adjustcolor("#dff0d8", 0.28), border = NA)
  abline(0, 1, col = "#555555", lwd = 1.4, lty = 2)
  arrows(full, full, full, excl, angle = 18, length = 0.06,
         col = "#666666", lwd = 1.0)
  points(full, excl, pch = 21, bg = cols, col = "#222222", cex = 1.45,
         lwd = 0.8)
  for (i in seq_along(samples)) {
    dx <- ifelse(samples[i] %in% c("hg02148", "na19036"), 0.02, 0.018)
    dy <- ifelse(samples[i] %in% c("chm13", "hg002"), -0.025, 0.018)
    text(full[i] + dx, excl[i] + dy, labels[samples[i]], cex = 0.68,
         adj = c(0, 0.5), col = "#222222")
  }
  text(0.12, 0.83, "above x=y:\nexclusion increases rho", adj = c(0, 1),
       cex = 0.88, font = 2, col = "#2f5d31")
  text(0.88, 0.08, "below x=y would mean\nexclusion weakens rho",
       adj = c(1, 0), cex = 0.72, col = "#666666")
  mtext("n = 7; HG002 CiFi was not run in the no-acrocentric control",
        side = 1, line = 4.2, cex = 0.70, col = "#555555")
}

render_11 <- function() {
  gm_c <- read_tsv(file.path(dipc_root,
    "output_q0_XX/community_enrichment_16cells_500kb_per_cell.tsv"))
  gm_p <- read_tsv(file.path(dipc_root,
    "output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv"))
  sp_c <- read_tsv(file.path(dipc_root,
    "sperm/enrichment_corrected/sperm_all20_per_cell.tsv"))
  sp_p <- read_tsv(file.path(dipc_root,
    "sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv"))

  s_all <- function(per_comm, cells) {
    x <- per_comm[per_comm$community == "S_all", ]
    x$ratio[match(cells, x$cell_id)]
  }
  d <- data.frame(
    group = factor(rep(c("GM12878\nC communities", "GM12878\nS_all control",
                         "Sperm\nC communities", "Sperm\nS_all control"),
                       c(nrow(gm_c), nrow(gm_c), nrow(sp_c), nrow(sp_c))),
                   levels = c("GM12878\nC communities", "GM12878\nS_all control",
                              "Sperm\nC communities", "Sperm\nS_all control")),
    ratio = c(gm_c$ratio, s_all(gm_p, gm_c$cell_id),
              sp_c$ratio, s_all(sp_p, sp_c$cell_id)),
    class = rep(c("community", "control", "community", "control"),
                c(nrow(gm_c), nrow(gm_c), nrow(sp_c), nrow(sp_c)))
  )
  set.seed(12)
  xpos <- as.numeric(d$group) + stats::runif(nrow(d), -0.10, 0.10)
  cols <- ifelse(d$class == "community", "#2b6cb0", "#999999")

  png(file.path(out_dir, "slide_11_single_cell_purpose_candidate.png"),
      width = 1800, height = 1200, res = 180, type = "cairo")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(0.72, 0.28))
  par(mar = c(5.8, 5.0, 5.0, 1.0), family = "sans")
  plot(NA, xlim = c(0.45, 4.55), ylim = range(d$ratio, na.rm = TRUE) * c(0.92, 1.08),
       xaxt = "n", xlab = "", ylab = "per-cell W/B distance ratio",
       main = "Single-cell 3D: the purpose is to rule out bulk averaging",
       cex.main = 0.98, cex.lab = 0.92, cex.axis = 0.82, las = 1)
  abline(h = 1, lty = 2, col = "#777777")
  axis(1, at = 1:4, labels = levels(d$group), tick = FALSE, cex.axis = 0.72)
  for (g in 1:4) {
    vals <- d$ratio[as.numeric(d$group) == g]
    qs <- stats::quantile(vals, c(0.25, 0.50, 0.75), na.rm = TRUE)
    rect(g - 0.20, qs[1], g + 0.20, qs[3], border = "#333333",
         col = NA, lwd = 1.1)
    segments(g - 0.20, qs[2], g + 0.20, qs[2], lwd = 1.6)
  }
  points(xpos, d$ratio, pch = 21, bg = cols, col = "#222222", cex = 0.85,
         lwd = 0.6)
  title(sub = "Below 1 = same-community arms are closer; S_all has zero-sharing arms and should not cluster",
        cex.sub = 0.72, line = 4.4)
  counts <- tapply(d$ratio < 1, d$group, function(x) sprintf("%d/%d < 1", sum(x), length(x)))
  text(1:4, par("usr")[4] - diff(par("usr")[3:4]) * 0.06, counts,
       cex = 0.72, font = 2, col = "#333333")

  par(mar = c(5.8, 0.6, 5.0, 0.6))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  draw_note("Recommendation:\nkeep Fig. 3c,\nbut retitle it as a test", 0.06, 0.70, 0.94, 0.95,
            cex = 0.80)
  text(0.08, 0.60, "Purpose", adj = c(0, 0.5), font = 2, cex = 0.82)
  text(0.08, 0.52,
       "Bulk Hi-C could be an average.\nPer-cell Dip-C/sperm asks\nwhether individual nuclei show\nsame-community proximity.",
       adj = c(0, 0.5), cex = 0.68)
  text(0.08, 0.32, "Why S_all matters", adj = c(0, 0.5),
       font = 2, cex = 0.82)
  text(0.08, 0.23,
       "Zero-sharing arms move the\nopposite way, so this is not\njust chromosome territory crowding.",
       adj = c(0, 0.5), cex = 0.68)
}

render_12 <- function() {
  d <- read_tsv(file.path(mouse_root, "community_analysis_1Mb/50000bp",
                          "zuo2021_zygotene_phr_pair_correlation.tsv"))
  d$chr_a <- sub("_[pq]$", "", d$arm_a)
  d$chr_b <- sub("_[pq]$", "", d$arm_b)
  d <- d[d$chr_a != d$chr_b & !is.na(d$mean_jaccard) &
           !is.na(d$hic_contact) & d$hic_contact > 0, ]
  rho <- stats::cor(d$mean_jaccard, d$hic_contact, method = "spearman")
  pval <- suppressWarnings(stats::cor.test(d$mean_jaccard, d$hic_contact,
                                           method = "spearman"))$p.value
  stages <- data.frame(stage = factor(c("leptotene", "zygotene",
                                        "pachytene", "diplotene"),
                                      levels = c("leptotene", "zygotene",
                                                 "pachytene", "diplotene")),
                       rho = c(0.687, 0.718, 0.683, 0.577))

  png(file.path(out_dir, "slide_12_mouse_zygotene_trajectory_pairing.png"),
      width = 1800, height = 1200, res = 180, type = "cairo")
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2, 1, 3), 2, 2, byrow = TRUE),
         widths = c(0.64, 0.36), heights = c(0.40, 0.60))
  par(mar = c(5.2, 5.2, 4.6, 1.2), family = "sans")
  plot(d$mean_jaccard, d$hic_contact, log = "y", pch = 21,
       bg = grDevices::adjustcolor("#4575b4", 0.55), col = "#27517a",
       xlab = "Mean PHR Jaccard similarity",
       ylab = "Hi-C contact (zygotene, log scale)",
       main = "Mouse zygotene: similar subtelomeres contact more",
       cex.main = 0.98, cex.lab = 0.90, cex.axis = 0.80)
  fit <- stats::lm(log10(hic_contact) ~ mean_jaccard, data = d)
  xs <- seq(min(d$mean_jaccard), max(d$mean_jaccard), length.out = 100)
  lines(xs, 10 ^ stats::predict(fit, newdata = data.frame(mean_jaccard = xs)),
        col = "#111111", lwd = 1.5)
  legend("topleft",
         legend = sprintf("Spearman rho = %.3f\np = %.1e\nn = %d pairs",
                          rho, pval, nrow(d)),
         bty = "n", cex = 0.78)

  par(mar = c(2.0, 1.2, 4.6, 1.2))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  title("Readable framing", cex.main = 0.92, font.main = 2)
  text(0.02, 0.80,
       "Use panel d for the\nzygotene scatter, but pair it\nwith the four-stage trajectory.\nThe trajectory explains why\nzygotene is being shown.",
       adj = c(0, 1), cex = 0.70)
  draw_note("One sentence on slide:\nZygotene is the bouquet stage,\nwhen telomeres cluster at the\nnuclear envelope.", 0.04, 0.08, 0.96, 0.42,
            cex = 0.63)

  par(mar = c(5.2, 4.4, 2.8, 1.0))
  plot(seq_along(stages$stage), stages$rho, type = "b", pch = 21,
       bg = ifelse(stages$stage == "zygotene", "#d62728", "#1f77b4"),
       col = "#333333", lwd = 1.2, cex = 1.4,
       xaxt = "n", xlim = c(0.85, 4.25), ylim = c(0.54, 0.75),
       xlab = "meiotic prophase stage",
       ylab = "Mantel rho",
       main = "Compact stage trajectory", cex.main = 0.88,
       cex.lab = 0.78, cex.axis = 0.72)
  axis(1, at = seq_along(stages$stage),
       labels = c("lepto", "zygo", "pachy", "diplo"), cex.axis = 0.72)
  text(seq_along(stages$stage), stages$rho + 0.012,
       sprintf("%.3f", stages$rho), cex = 0.62)
  text(2, 0.742, "bouquet", cex = 0.68, font = 2, col = "#d62728")
}

render_10a()
render_10b()
render_11()
render_12()

message("Wrote Hi-C redesign candidates to ", out_dir)
