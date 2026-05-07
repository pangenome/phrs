#!/usr/bin/env Rscript

out_dir <- "slides/v2-review-zoom/_revision_assets/v3/11_wb_labels"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dipc_root <- "/moosefs/guarracino/HPRCv2/dipc_t2t"

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE)
}

extract_s_all <- function(per_comm, cells) {
  s_all <- per_comm[per_comm$community == "S_all", ]
  s_all$ratio[match(cells, s_all$cell_id)]
}

gm_c <- read_tsv(file.path(
  dipc_root,
  "output_q0_XX/community_enrichment_16cells_500kb_per_cell.tsv"
))
gm_p <- read_tsv(file.path(
  dipc_root,
  "output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv"
))
sp_c <- read_tsv(file.path(
  dipc_root,
  "sperm/enrichment_corrected/sperm_all20_per_cell.tsv"
))
sp_p <- read_tsv(file.path(
  dipc_root,
  "sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv"
))

groups <- c(
  "GM12878 Dip-C\nC communities",
  "GM12878 Dip-C\nS_all zero-sharing\ncontrol",
  "Sperm scHi-C\nC communities",
  "Sperm scHi-C\nS_all zero-sharing\ncontrol"
)
summary_groups <- c(
  "GM12878 Dip-C C communities",
  "GM12878 Dip-C S_all zero-sharing control",
  "Sperm scHi-C C communities",
  "Sperm scHi-C S_all zero-sharing control"
)

d <- data.frame(
  group = factor(rep(groups, c(nrow(gm_c), nrow(gm_c), nrow(sp_c), nrow(sp_c))),
                 levels = groups),
  ratio = c(gm_c$ratio, extract_s_all(gm_p, gm_c$cell_id),
            sp_c$ratio, extract_s_all(sp_p, sp_c$cell_id)),
  class = rep(c("community", "control", "community", "control"),
              c(nrow(gm_c), nrow(gm_c), nrow(sp_c), nrow(sp_c)))
)

count_below_one <- function(vals) {
  sprintf("%d/%d lower than 1", sum(vals < 1, na.rm = TRUE),
          sum(!is.na(vals)))
}

group_counts <- tapply(d$ratio, d$group, count_below_one)
group_means <- tapply(d$ratio, d$group, function(x) mean(x, na.rm = TRUE))

summary_table <- data.frame(
  group = summary_groups,
  plotted_mean_within_over_between_3d_distance = as.numeric(group_means[groups]),
  cells_lower_than_one = as.character(group_counts[groups]),
  stringsAsFactors = FALSE
)
write.table(
  summary_table,
  file = file.path(out_dir, "slide11_explicit_distance_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

render_candidate <- function(device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(file.path(out_dir, "slide11_explicit_distance_labels_candidate.png"),
        width = 1800, height = 1200, res = 180, type = "cairo")
  } else {
    pdf(file.path(out_dir, "slide11_explicit_distance_labels_candidate.pdf"),
        width = 10, height = 6.6667, useDingbats = FALSE)
  }

  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(0.68, 0.32))

  par(mar = c(6.4, 5.9, 5.2, 1.0), family = "sans")
  ylim <- range(d$ratio, na.rm = TRUE) * c(0.90, 1.10)
  plot(NA, xlim = c(0.45, 4.55), ylim = ylim,
       xaxt = "n",
       xlab = "",
       ylab = "Within-community 3D distance / between-community 3D distance",
       main = "Single-cell 3D distance: spell out the ratio",
       cex.main = 1.02,
       cex.lab = 0.90,
       cex.axis = 0.80,
       las = 1)
  grid(nx = NA, ny = NULL, col = "#e6e6e6", lwd = 0.7)
  abline(h = 1, lty = 2, col = "#5f5f5f", lwd = 1.4)
  axis(1, at = 1:4, labels = groups, tick = FALSE,
       cex.axis = 0.60, line = 0.1)

  set.seed(1107)
  xpos <- as.numeric(d$group) + stats::runif(nrow(d), -0.10, 0.10)
  point_cols <- ifelse(d$class == "community", "#2f6fb0", "#9b9b9b")
  box_cols <- ifelse(d$class == "community", "#174b7a", "#555555")

  for (g in seq_along(groups)) {
    vals <- d$ratio[as.numeric(d$group) == g]
    qs <- stats::quantile(vals, c(0.25, 0.50, 0.75), na.rm = TRUE)
    rect(g - 0.22, qs[1], g + 0.22, qs[3],
         border = box_cols[match(g, as.numeric(d$group))][1],
         col = NA,
         lwd = 1.1)
    segments(g - 0.22, qs[2], g + 0.22, qs[2],
             col = box_cols[match(g, as.numeric(d$group))][1],
             lwd = 1.6)
  }

  points(xpos, d$ratio, pch = 21, bg = point_cols, col = "#1f1f1f",
         cex = 0.84, lwd = 0.55)
  text(1:4,
       par("usr")[4] - diff(par("usr")[3:4]) * 0.062,
       group_counts,
       cex = 0.61,
       font = 2,
       col = "#262626")

  text(0.50, 1.02, "1 = no distance difference",
       adj = c(0, 0.5), cex = 0.67, col = "#4c4c4c")
  text(0.50, 0.92, "lower than 1 = same-community arms are closer",
       adj = c(0, 0.5), cex = 0.67, font = 2, col = "#174b7a")

  title(sub = "Distance metric, not contact; lower than 1 means closer; S_all is the zero-sharing negative control.",
        cex.sub = 0.62, line = 4.8)

  par(mar = c(6.4, 0.8, 5.2, 1.0), family = "sans")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)

  rect(0.04, 0.68, 0.96, 0.96, col = "#f5f7fb", border = "#6d7f99",
       lwd = 1.0)
  text(0.08, 0.91, "Use explicit labels", adj = c(0, 0.5),
       font = 2, cex = 0.83, col = "#1f2d3d")
  text(0.08, 0.845,
       "mean within-community\n3D distance",
       adj = c(0, 0.5), cex = 0.60, col = "#1f2d3d")
  text(0.08, 0.765, "divided by", adj = c(0, 0.5),
       cex = 0.62, col = "#1f2d3d")
  text(0.08, 0.715,
       "mean between-community\n3D distance",
       adj = c(0, 0.5), cex = 0.60, col = "#1f2d3d")

  rect(0.04, 0.40, 0.96, 0.66, col = "#fff8eb", border = "#9a6b00",
       lwd = 1.0)
  text(0.08, 0.60, "Direction", adj = c(0, 0.5),
       font = 2, cex = 0.83, col = "#3c2b00")
  text(0.08, 0.50,
       "Lower than 1 means the\nsame-community arms are\ncloser in the 3D model.",
       adj = c(0, 0.5), cex = 0.67, col = "#3c2b00")

  rect(0.04, 0.11, 0.96, 0.35, col = "#f4f4f4", border = "#777777",
       lwd = 1.0)
  text(0.08, 0.30, "S_all control", adj = c(0, 0.5),
       font = 2, cex = 0.83, col = "#333333")
  text(0.08, 0.21,
       "Seven arms with zero\nsubtelomeric sequence sharing.\nThey move the opposite way.",
       adj = c(0, 0.5), cex = 0.67, col = "#333333")
}

render_candidate("png")
render_candidate("pdf")
