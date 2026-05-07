#!/usr/bin/env Rscript

# Slide 09 v3 asset: square equal-scale MDS / PCoA with all Leiden C1-C15
# communities labeled. The coordinate source is Andrea's cached cmdscale output;
# the community source is the canonical arm-level Leiden k=15 assignment TSV.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

similarity_dir <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
mds_rds <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.full_mds.rds")
assign_tsv <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")

out_dir <- "slides/v2-review-zoom/_revision_assets/v3/09_all_communities_1to1"
out_base <- file.path(out_dir, "mds_pcoa_all_communities_1to1")
label_tsv <- file.path(out_dir, "label_positions.tsv")
validation_tsv <- file.path(out_dir, "validation_summary.tsv")

fit <- readRDS(mds_rds)
if (!is.list(fit) || is.null(fit$points) || is.null(fit$eig)) {
  stop("Expected full_mds.rds to contain a cmdscale-style list with points and eig")
}
if (ncol(fit$points) < 2) {
  stop("MDS / PCoA point matrix must contain at least two dimensions")
}

var_explained <- fit$eig / sum(abs(fit$eig)) * 100

points <- as.data.frame(fit$points[, 1:2])
names(points) <- c("x", "y")
points$Name <- rownames(fit$points)
points$ChromArm <- sub(".*_(chr[^_]+)_(p|q)arm$", "\\1_\\2", points$Name)
points$Arm <- sub("^.*_([pq])$", "\\1", points$ChromArm)

arm_assignments <- read.delim(assign_tsv, stringsAsFactors = FALSE)
required_assignment_cols <- c("ChromArm", "Community", "Arms")
missing_cols <- setdiff(required_assignment_cols, names(arm_assignments))
if (length(missing_cols) > 0) {
  stop("Assignment TSV is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

points <- merge(
  points,
  arm_assignments[, required_assignment_cols],
  by = "ChromArm",
  all.x = TRUE,
  sort = FALSE
)

if (any(is.na(points$Community))) {
  missing_arms <- sort(unique(points$ChromArm[is.na(points$Community)]))
  stop("Missing community assignments for arm(s): ", paste(missing_arms, collapse = ", "))
}

community_levels <- paste0("C", 1:15)
observed_communities <- unique(points$Community)
observed_communities <- observed_communities[order(as.integer(sub("^C", "", observed_communities)))]
if (!identical(observed_communities, community_levels)) {
  stop(
    "Expected communities C1-C15; observed: ",
    paste(observed_communities, collapse = ", ")
  )
}
points$Community <- factor(points$Community, levels = community_levels)

# Preserve the v2 / slide-07 named-clade palette for the abstract-anchored
# communities. Additional communities use muted, distinct colors so the all-C1-C15
# callouts remain legible without changing those six established meanings.
community_colors <- c(
  C1 = "#A65628",  # DUX4 / D4Z4, v2 named-clade brown
  C2 = "#FF7F00",  # 10p-18p, v2 named-clade orange
  C3 = "#CC79A7",
  C4 = "#009E73",
  C5 = "#56B4E9",
  C6 = "#984EA3",  # tight q-arm clade, v2 named-clade purple
  C7 = "#4DAF4A",  # acrocentric p-arms, v2 named-clade green
  C8 = "#666666",
  C9 = "#F0E442",
  C10 = "#0072B2",
  C11 = "#8DA0CB",
  C12 = "#E78AC3",
  C13 = "#BDBDBD",
  C14 = "#377EB8", # PAR2, v2 named-clade blue
  C15 = "#E41A1C"  # PAR1, v2 named-clade red
)

community_arms <- unique(arm_assignments[, c("Community", "Arms")])
community_arms$num <- as.integer(sub("^C", "", community_arms$Community))
community_arms <- community_arms[order(community_arms$num), ]

centroids <- aggregate(
  cbind(x, y) ~ Community,
  data = points,
  FUN = median
)
centroids$Community <- as.character(centroids$Community)
names(centroids)[names(centroids) == "x"] <- "anchor_x"
names(centroids)[names(centroids) == "y"] <- "anchor_y"

counts <- aggregate(Name ~ Community, data = points, FUN = length)
names(counts)[names(counts) == "Name"] <- "n_points"
counts$Community <- as.character(counts$Community)

arm_counts <- aggregate(ChromArm ~ Community, data = arm_assignments, FUN = length)
names(arm_counts)[names(arm_counts) == "ChromArm"] <- "n_arms"

# Fixed callout centers. Labels sit on the plot perimeter and the dense
# lower-central C2/C4/C7/C9/C10/C13/C14/C15 region is deliberately fanned out.
label_positions <- data.frame(
  Community = community_levels,
  label_x = c(
    0.52, -0.18, -0.50, 0.05, -0.42,
    0.46, 0.30, -0.50, 0.52, 0.52,
    -0.50, -0.50, 0.52, 0.52, 0.52
  ),
  label_y = c(
    0.18, -0.55, 0.27, -0.55, -0.55,
    0.50, -0.55, -0.08, -0.30, 0.06,
    0.47, 0.09, -0.42, -0.06, -0.18
  ),
  stringsAsFactors = FALSE
)

label_df <- Reduce(
  function(left, right) merge(left, right, by = "Community", sort = FALSE),
  list(label_positions, community_arms[, c("Community", "Arms")], centroids, counts, arm_counts)
)
label_df$num <- as.integer(sub("^C", "", label_df$Community))
label_df <- label_df[order(label_df$num), ]
label_df$Community <- factor(label_df$Community, levels = community_levels)

wrap_arms <- function(arms_string) {
  arms <- strsplit(arms_string, ", ", fixed = TRUE)[[1]]
  if (length(arms) <= 3) {
    return(paste(arms, collapse = " "))
  }
  split_at <- ceiling(length(arms) / 2)
  paste(
    paste(arms[seq_len(split_at)], collapse = " "),
    paste(arms[(split_at + 1):length(arms)], collapse = " "),
    sep = "\n"
  )
}

label_df$label <- paste0(
  label_df$Community,
  "\n",
  vapply(label_df$Arms, wrap_arms, character(1))
)
label_df$fill <- adjustcolor(community_colors[as.character(label_df$Community)], alpha.f = 0.18)

label_table <- label_df
label_table$label <- gsub("\n", " | ", label_table$label, fixed = TRUE)
write.table(
  label_table[, c(
    "Community", "Arms", "n_arms", "n_points",
    "anchor_x", "anchor_y", "label_x", "label_y", "label"
  )],
  file = label_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

axis_limits <- c(-0.68, 0.68)

p <- ggplot(points, aes(x = x, y = y, color = Community)) +
  geom_point(size = 0.72, alpha = 0.38) +
  geom_segment(
    data = label_df,
    aes(
      x = anchor_x,
      y = anchor_y,
      xend = label_x,
      yend = label_y,
      color = Community
    ),
    inherit.aes = FALSE,
    linewidth = 0.42,
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = label_df,
    aes(x = anchor_x, y = anchor_y, fill = Community),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.35,
    color = "black",
    alpha = 0.92
  ) +
  geom_label(
    data = label_df,
    aes(x = label_x, y = label_y, label = label, fill = Community),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 3.15,
    lineheight = 0.88,
    label.size = 0.32,
    label.r = unit(0.12, "lines")
  ) +
  scale_color_manual(values = community_colors, guide = "none") +
  scale_fill_manual(
    values = setNames(
      adjustcolor(community_colors, alpha.f = 0.18),
      names(community_colors)
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = axis_limits,
    breaks = seq(-0.5, 0.5, by = 0.25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = axis_limits,
    breaks = seq(-0.5, 0.5, by = 0.25),
    expand = c(0, 0)
  ) +
  coord_fixed(ratio = 1, xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "MDS / PCoA: all Leiden communities labeled",
    subtitle = "Classical MDS on 1 - Jaccard; labels are C1-C15 arm-level Leiden communities; axes use equal scale",
    x = sprintf("MDS dimension 1 (%.2f%%)", var_explained[1]),
    y = sprintf("MDS dimension 2 (%.2f%%)", var_explained[2])
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 10.8, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.28),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave(paste0(out_base, ".png"), p, width = 12, height = 12, dpi = 300)
ggsave(paste0(out_base, ".pdf"), p, width = 12, height = 12)

validation <- data.frame(
  check = c(
    "labels_present",
    "coordinate_source",
    "assignment_source",
    "terminology",
    "coord_fixed_ratio",
    "axis_limit_x",
    "axis_limit_y",
    "output_inches"
  ),
  value = c(
    paste(as.character(label_df$Community), collapse = ","),
    mds_rds,
    assign_tsv,
    "MDS / PCoA, not PCA",
    "1",
    paste(axis_limits, collapse = ","),
    paste(axis_limits, collapse = ","),
    "12 x 12"
  )
)

write.table(
  validation,
  file = validation_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote ", paste0(out_base, ".png"))
message("Wrote ", paste0(out_base, ".pdf"))
message("Wrote ", label_tsv)
message("Wrote ", validation_tsv)
