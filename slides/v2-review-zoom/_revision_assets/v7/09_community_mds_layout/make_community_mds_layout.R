#!/usr/bin/env Rscript

# Slide 09 v7 asset: MDS / PCoA scatter using the same 12 x 10 inch
# plot grammar as slide 08a, but colored and labeled by arm-level C1-C15
# Leiden community. The coordinates are the cached cmdscale output from
# the HPRCv2 graph-path Jaccard distance workflow.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

repo_root <- tryCatch(
  normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)[1]),
  error = function(e) normalizePath(getwd())
)

path_from_root <- function(...) file.path(repo_root, ...)

similarity_dir <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
mds_rds <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.full_mds.rds")
assign_tsv <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
v3_label_tsv <- path_from_root(
  "slides/v2-review-zoom/_revision_assets/v3/09_all_communities_1to1/label_positions.tsv"
)
s08a_asset <- path_from_root("slides/v2-review-zoom/_typst/assets/s08a_mds_chrom.png")

out_dir <- path_from_root("slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(out_dir, "mds_community_layout")
label_tsv <- file.path(out_dir, "label_positions.tsv")
validation_tsv <- file.path(out_dir, "validation_summary.tsv")

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

assert_file(mds_rds, "MDS coordinate RDS")
assert_file(assign_tsv, "Arm-level Leiden assignment TSV")
assert_file(v3_label_tsv, "v3 all-community label-position TSV")
assert_file(s08a_asset, "slide 08a layout reference PNG")

community_number <- function(x) {
  suppressWarnings(as.integer(sub("^C", "", as.character(x))))
}

short_arms <- function(arms_string) {
  gsub("chr", "", gsub("_", "", arms_string), fixed = TRUE)
}

wrap_arms <- function(arms_string) {
  arms <- strsplit(arms_string, ", ", fixed = TRUE)[[1]]
  arms <- short_arms(arms)
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

fit <- readRDS(mds_rds)
if (!is.list(fit) || is.null(fit$points) || is.null(fit$eig)) {
  stop("Expected full_mds.rds to contain a cmdscale-style list with points and eig", call. = FALSE)
}
if (ncol(fit$points) < 2) {
  stop("MDS point matrix must contain at least two dimensions", call. = FALSE)
}

var_explained <- fit$eig / sum(abs(fit$eig)) * 100

points <- as.data.frame(fit$points[, 1:2])
names(points) <- c("x", "y")
points$Name <- rownames(fit$points)
points$ChromArm <- sub(".*_(chr[^_]+)_(p|q)arm$", "\\1_\\2", points$Name)
points$Arm <- sub("^.*_([pq])$", "\\1", points$ChromArm)

assignments <- read.delim(assign_tsv, stringsAsFactors = FALSE)
required_assignment_cols <- c("ChromArm", "Community", "Arms")
missing_cols <- setdiff(required_assignment_cols, names(assignments))
if (length(missing_cols) > 0) {
  stop("Assignment TSV is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

points <- merge(
  points,
  assignments[, required_assignment_cols],
  by = "ChromArm",
  all.x = TRUE,
  sort = FALSE
)

if (any(is.na(points$Community))) {
  missing_arms <- sort(unique(points$ChromArm[is.na(points$Community)]))
  stop("Missing community assignments for arm(s): ", paste(missing_arms, collapse = ", "), call. = FALSE)
}

community_levels <- paste0("C", 1:15)
observed_communities <- sort(unique(as.character(points$Community)), method = "radix")
observed_communities <- observed_communities[order(community_number(observed_communities))]
if (!identical(observed_communities, community_levels)) {
  stop("Expected communities C1-C15; observed: ", paste(observed_communities, collapse = ", "), call. = FALSE)
}

points$Community <- factor(points$Community, levels = community_levels)
points$Arm <- factor(points$Arm, levels = c("p", "q"))

# C1/C2/C6/C7/C14/C15 retain the named-clade colors used by earlier review
# slides. Other C labels use muted distinct colors from the existing v3
# all-community slide so C1-C15 remain stable across revisions.
community_colors <- c(
  C1 = "#A65628",
  C2 = "#FF7F00",
  C3 = "#CC79A7",
  C4 = "#009E73",
  C5 = "#56B4E9",
  C6 = "#984EA3",
  C7 = "#4DAF4A",
  C8 = "#666666",
  C9 = "#F0E442",
  C10 = "#0072B2",
  C11 = "#8DA0CB",
  C12 = "#E78AC3",
  C13 = "#BDBDBD",
  C14 = "#377EB8",
  C15 = "#E41A1C"
)

label_positions <- read.delim(v3_label_tsv, stringsAsFactors = FALSE)
required_label_cols <- c("Community", "label_x", "label_y")
missing_label_cols <- setdiff(required_label_cols, names(label_positions))
if (length(missing_label_cols) > 0) {
  stop("v3 label TSV is missing required column(s): ", paste(missing_label_cols, collapse = ", "), call. = FALSE)
}
label_positions <- label_positions[, required_label_cols]
label_positions$Community <- as.character(label_positions$Community)
label_positions <- label_positions[order(community_number(label_positions$Community)), ]
if (!identical(label_positions$Community, community_levels)) {
  stop("v3 label TSV must contain exactly C1-C15 in numeric order", call. = FALSE)
}

community_arms <- unique(assignments[, c("Community", "Arms")])
community_arms$Community <- as.character(community_arms$Community)
community_arms$num <- community_number(community_arms$Community)
community_arms <- community_arms[order(community_arms$num), ]

centroids <- aggregate(cbind(x, y) ~ Community, data = points, FUN = median)
centroids$Community <- as.character(centroids$Community)
names(centroids)[names(centroids) == "x"] <- "anchor_x"
names(centroids)[names(centroids) == "y"] <- "anchor_y"

point_counts <- aggregate(Name ~ Community, data = points, FUN = length)
point_counts$Community <- as.character(point_counts$Community)
names(point_counts)[names(point_counts) == "Name"] <- "n_points"

arm_counts <- aggregate(ChromArm ~ Community, data = assignments, FUN = length)
arm_counts$Community <- as.character(arm_counts$Community)
names(arm_counts)[names(arm_counts) == "ChromArm"] <- "n_arms"

label_df <- Reduce(
  function(left, right) merge(left, right, by = "Community", sort = FALSE),
  list(label_positions, community_arms[, c("Community", "Arms")], centroids, point_counts, arm_counts)
)
label_df$num <- community_number(label_df$Community)
label_df <- label_df[order(label_df$num), ]
if (!identical(label_df$Community, community_levels)) {
  stop("Label data must contain exactly C1-C15 in numeric order", call. = FALSE)
}

label_df$label <- paste0(
  label_df$Community,
  "\n",
  vapply(label_df$Arms, wrap_arms, character(1))
)
label_df$Community <- factor(label_df$Community, levels = community_levels)

legend_labels <- setNames(
  paste0(
    community_levels,
    " (",
    arm_counts$n_arms[match(community_levels, arm_counts$Community)],
    " arms)"
  ),
  community_levels
)

axis_limits <- c(-0.68, 0.68)
axis_breaks <- seq(-0.5, 0.5, by = 0.25)

p <- ggplot(points, aes(x = x, y = y, color = Community, shape = Arm)) +
  geom_point(size = 1.75, alpha = 0.58, stroke = 0) +
  geom_segment(
    data = label_df,
    aes(x = anchor_x, y = anchor_y, xend = label_x, yend = label_y, color = Community),
    inherit.aes = FALSE,
    linewidth = 0.55,
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = label_df,
    aes(x = anchor_x, y = anchor_y, color = Community),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2.7,
    stroke = 0.55,
    alpha = 0.96
  ) +
  geom_label(
    data = label_df,
    aes(x = label_x, y = label_y, label = label, fill = Community),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 3.6,
    lineheight = 0.88,
    label.size = 0.42,
    label.r = unit(0.12, "lines")
  ) +
  scale_color_manual(
    values = community_colors,
    breaks = community_levels,
    labels = legend_labels,
    name = "Community"
  ) +
  scale_fill_manual(
    values = setNames(adjustcolor(community_colors, alpha.f = 0.18), names(community_colors)),
    guide = "none"
  ) +
  scale_shape_manual(values = c(p = 16, q = 17), name = "Arm") +
  scale_x_continuous(limits = axis_limits, breaks = axis_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = axis_limits, breaks = axis_breaks, expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "hprcv2.1Mb.subtelo - Full MDS colored by Leiden community",
    subtitle = paste(
      "Classical MDS on 1 - Jaccard; C1-C15 are arm-level Leiden communities",
      "assigned from graph-path Jaccard distances, not 3D/gene labels",
      sep = "\n"
    ),
    x = sprintf("Dimension 1 (%.2f%%)", var_explained[1]),
    y = sprintf("Dimension 2 (%.2f%%)", var_explained[2])
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    plot.subtitle = element_text(hjust = 0.5, size = 12.6, lineheight = 0.95),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12.4),
    legend.key.height = unit(0.33, "in"),
    legend.key.width = unit(0.32, "in"),
    legend.box = "vertical",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.34),
    panel.border = element_rect(color = "#303030", fill = NA, linewidth = 0.55),
    plot.margin = margin(9, 13, 9, 12)
  ) +
  guides(
    color = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1, shape = 16)),
    shape = guide_legend(override.aes = list(size = 5, alpha = 1, color = "#4D4D4D"))
  )

ggsave(paste0(out_base, ".png"), p, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(paste0(out_base, ".pdf"), p, width = 12, height = 10, bg = "white")

label_table <- label_df
label_table$Community <- as.character(label_table$Community)
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

validation <- data.frame(
  check = c(
    "coordinate_method",
    "coordinate_source",
    "assignment_source",
    "layout_reference",
    "output_plot_inches",
    "output_png_pixels",
    "coord_fixed_ratio",
    "axis_limit_x",
    "axis_limit_y",
    "point_count",
    "labels_present",
    "support_text",
    "pca_used"
  ),
  value = c(
    "Classical MDS / cmdscale on 1 - graph-path Jaccard distances",
    mds_rds,
    assign_tsv,
    "slides/v2-review-zoom/_typst/assets/s08a_mds_chrom.png; 12 x 10 inch frame, right legend, theme_bw grammar",
    "12 x 10",
    "3600 x 3000",
    "1",
    paste(axis_limits, collapse = ","),
    paste(axis_limits, collapse = ","),
    as.character(nrow(points)),
    paste(as.character(label_df$Community), collapse = ","),
    "C1-C15 assigned from arm-level graph-path Jaccard distances, not 3D/gene labels",
    "no"
  ),
  stringsAsFactors = FALSE
)
write.table(validation, file = validation_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote ", paste0(out_base, ".png"))
message("Wrote ", paste0(out_base, ".pdf"))
message("Wrote ", label_tsv)
message("Wrote ", validation_tsv)
