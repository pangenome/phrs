#!/usr/bin/env Rscript

# Candidate slide-09 community MDS asset for review-zoom-09-mds-community-leading.
# It reads Andrea's cached classical MDS / PCoA coordinates and the arm-level
# Leiden k=15 assignments, then produces a talk-facing plot that labels the six
# abstract-named communities directly on the MDS.

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
mds_rds <- file.path(root, "hprcv2.1Mb.subtelo.full_mds.rds")
assign_tsv <- file.path(root, "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")

out_dir <- "slides/v2-review-zoom/_revision_assets/09_mds_community"
out_base <- file.path(out_dir, "candidate_labeled_mds_community")

fit <- readRDS(mds_rds)
var_explained <- fit$eig / sum(abs(fit$eig)) * 100

points <- as.data.frame(fit$points[, 1:2])
names(points) <- c("x", "y")
points$Name <- rownames(fit$points)
points$ChromArm <- sub(".*_(chr[^_]+)_(p|q)arm$", "\\1_\\2", points$Name)

arm_assignments <- read.delim(assign_tsv, stringsAsFactors = FALSE)
points <- merge(points, arm_assignments[, c("ChromArm", "Community", "Arms")],
                by = "ChromArm", all.x = TRUE, sort = FALSE)

named_order <- c("C15", "C14", "C7", "C2", "C6", "C1")
legend_labels <- c(
  "C15 PAR1 (Xp/Yp)",
  "C14 PAR2 (Xq/Yq)",
  "C7 acrocentric p-arms",
  "C2 10p-18p",
  "C6 tight q-arm clade",
  "C1 DUX4 / D4Z4",
  "Other Leiden communities"
)
names(legend_labels) <- c(named_order, "Other")

# Match the slide-07 NJ-tree palette for the six named abstract clades.
point_palette <- c(
  C15 = "#E41A1C",  # PAR1
  C14 = "#377EB8",  # PAR2
  C7 = "#4DAF4A",   # acrocentric p-arms
  C2 = "#FF7F00",   # 10p-18p
  C6 = "#984EA3",   # tight q-arm clade
  C1 = "#A65628",   # DUX4
  Other = "#BDBDBD"
)

points$Legend <- ifelse(points$Community %in% named_order, points$Community, "Other")
points$Legend <- factor(points$Legend, levels = c(named_order, "Other"))

label_df <- data.frame(
  Community = named_order,
  label = c(
    "C15 PAR1\nXp / Yp",
    "C14 PAR2\nXq / Yq",
    "C7 acrocentric p\n13p/14p/15p/21p/22p",
    "C2 10p-18p",
    "C6 tight q-arm\n1q/13q/17q/19q/21q/22q",
    "C1 DUX4/D4Z4\n4q / 10q"
  ),
  label_x = c(-0.040, 0.115, 0.185, -0.085, 0.325, 0.285),
  label_y = c(-0.145, -0.265, -0.335, -0.250, 0.310, -0.020),
  stringsAsFactors = FALSE
)

centroids <- aggregate(cbind(x, y) ~ Community, points[points$Community %in% named_order, ],
                       FUN = median)
label_df <- merge(label_df, centroids, by = "Community", sort = FALSE)
label_df$Legend <- factor(label_df$Community, levels = c(named_order, "Other"))

label_fill <- vapply(point_palette, function(color) {
  adjustcolor(color, alpha.f = 0.18)
}, character(1))

p <- ggplot(points, aes(x = x, y = y)) +
  geom_point(data = points[points$Legend == "Other", ],
             aes(color = Legend), size = 0.8, alpha = 0.28) +
  geom_point(data = points[points$Legend != "Other", ],
             aes(color = Legend), size = 1.15, alpha = 0.72) +
  geom_segment(data = label_df,
               aes(x = x, y = y, xend = label_x, yend = label_y, color = Legend),
               linewidth = 0.55, alpha = 0.85) +
  geom_label(data = label_df,
             aes(x = label_x, y = label_y, label = label, fill = Legend),
             color = "black", fontface = "bold", size = 3.0,
             lineheight = 0.92, label.size = 0.35,
             label.r = unit(0.12, "lines")) +
  scale_color_manual(
    values = point_palette,
    breaks = c(named_order, "Other"),
    labels = legend_labels,
    name = "Community story"
  ) +
  scale_fill_manual(values = label_fill, guide = "none") +
  coord_cartesian(xlim = c(-0.43, 0.58), ylim = c(-0.48, 0.42), expand = FALSE) +
  labs(
    title = "MDS / PCoA - Leiden communities with named clades",
    subtitle = "Classical MDS on 1 - Jaccard; direct labels use the slide-07 named-clade palette",
    x = sprintf("MDS dimension 1 (%.2f%%)", var_explained[1]),
    y = sprintf("MDS dimension 2 (%.2f%%)", var_explained[2])
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 12, 10, 10)
  )

ggsave(paste0(out_base, ".png"), p, width = 14, height = 8, dpi = 300)
ggsave(paste0(out_base, ".pdf"), p, width = 14, height = 8)

message("Wrote ", paste0(out_base, ".png"))
message("Wrote ", paste0(out_base, ".pdf"))
