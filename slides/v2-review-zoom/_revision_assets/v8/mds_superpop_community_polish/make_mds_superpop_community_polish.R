#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(grid)
})

repo_root <- tryCatch(
  normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)[1]),
  error = function(e) normalizePath(getwd())
)

path_from_root <- function(...) file.path(repo_root, ...)

out_dir <- path_from_root("slides/v2-review-zoom/_revision_assets/v8/mds_superpop_community_polish")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

similarity_dir <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
mds_rds <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.full_mds.rds")
seq_assignment_tsv <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv")
arm_assignment_tsv <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
upstream_superpop_png <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png")
upstream_chromosome_png <- file.path(similarity_dir, "hprcv2.1Mb.subtelo.mds.color-by-chromosome.png")
slide08a_asset <- path_from_root("slides/v2-review-zoom/_typst/assets/s08a_mds_chrom.png")
slide08b_asset <- path_from_root("slides/v2-review-zoom/_typst/assets/s08b_mds_superpop.png")
source_script <- "/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R"

superpop_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")
superpop_colors <- c(
  "AFR" = "#FFCD33",
  "AMR" = "#ED1E24",
  "EAS" = "#108C44",
  "EUR" = "#6AA5CD",
  "SAS" = "#9B59B6"
)

arm_shapes <- c("p" = 16, "q" = 17)

community_levels <- paste0("C", 1:15)
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

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
format_sci <- function(x, digits = 2) formatC(x, format = "e", digits = digits)
format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 1e-3, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
  )
}
format_signed_sci <- function(x, digits = 2) {
  paste0(ifelse(x > 0, "+", ""), format_sci(x, digits))
}
format_signed_num <- function(x, digits = 2) {
  paste0(ifelse(x > 0, "+", ""), format_num(x, digits))
}

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  line <- system2("sha256sum", path, stdout = TRUE)
  strsplit(line, "[[:space:]]+")[[1]][1]
}

community_number <- function(x) suppressWarnings(as.integer(sub("^C", "", as.character(x))))

short_arms <- function(arms_string) {
  gsub("chr", "", gsub("_", "", arms_string, fixed = TRUE), fixed = TRUE)
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

pick_existing_dims <- function(dat, dims) dims[dims %in% colnames(dat)]

nearest_by_superpop <- function(points, dims) {
  if (length(dims) < 2) {
    stop("Nearest-neighbor distance requires at least two MDS dimensions.", call. = FALSE)
  }

  out <- lapply(superpop_order, function(superpop) {
    dat <- points %>%
      filter(Superpopulation == superpop) %>%
      arrange(Name)

    if (nrow(dat) < 2) {
      stop(sprintf("Superpopulation %s has fewer than two points.", superpop), call. = FALSE)
    }

    mat <- as.matrix(dat[, dims, drop = FALSE])
    dm <- as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
    diag(dm) <- Inf

    nearest_idx <- max.col(-dm, ties.method = "first")
    nearest_dist <- dm[cbind(seq_len(nrow(dm)), nearest_idx)]
    nearest <- dat[nearest_idx, , drop = FALSE]

    tibble(
      Name = dat$Name,
      nearest_same_superpop_name = nearest$Name,
      nearest_same_superpop_distance_d1_d2 = nearest_dist,
      nearest_same_superpop_sample = nearest$Sample,
      nearest_same_superpop_chromarm = nearest$ChromArm
    )
  })

  bind_rows(out)
}

summarize_nearest <- function(dat) {
  dat %>%
    group_by(Superpopulation) %>%
    summarise(
      n_points = n(),
      n_samples = n_distinct(Sample),
      min_nearest_distance = min(nearest_same_superpop_distance_d1_d2),
      q05_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.05)),
      q25_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.25)),
      median_nearest_distance = median(nearest_same_superpop_distance_d1_d2),
      mean_nearest_distance = mean(nearest_same_superpop_distance_d1_d2),
      q75_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.75)),
      q90_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.90)),
      q95_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.95)),
      q99_nearest_distance = unname(quantile(nearest_same_superpop_distance_d1_d2, 0.99)),
      max_nearest_distance = max(nearest_same_superpop_distance_d1_d2),
      sd_nearest_distance = sd(nearest_same_superpop_distance_d1_d2),
      n_exact_zero = sum(nearest_same_superpop_distance_d1_d2 == 0),
      n_le_1e_minus_7 = sum(nearest_same_superpop_distance_d1_d2 <= 1e-7),
      .groups = "drop"
    ) %>%
    mutate(Superpopulation = as.character(Superpopulation)) %>%
    arrange(match(Superpopulation, superpop_order))
}

pairwise_stats <- function(dat) {
  combos <- combn(superpop_order, 2, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    g1 <- pair[[1]]
    g2 <- pair[[2]]
    x <- dat$nearest_same_superpop_distance_d1_d2[dat$Superpopulation == g1]
    y <- dat$nearest_same_superpop_distance_d1_d2[dat$Superpopulation == g2]
    nx <- length(x)
    ny <- length(y)

    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = FALSE))
    ranks <- rank(c(x, y), ties.method = "average")
    u <- sum(ranks[seq_len(nx)]) - nx * (nx + 1) / 2
    cliff_delta <- 2 * u / (nx * ny) - 1

    tibble(
      group1 = g1,
      group2 = g2,
      n_group1 = nx,
      n_group2 = ny,
      median_group1 = median(x),
      median_group2 = median(y),
      median_difference_group1_minus_group2 = median(x) - median(y),
      mean_group1 = mean(x),
      mean_group2 = mean(y),
      mean_difference_group1_minus_group2 = mean(x) - mean(y),
      wilcoxon_w = unname(wt$statistic),
      p_value = wt$p.value,
      cliff_delta_group1_vs_group2 = cliff_delta
    )
  })

  bind_rows(rows) %>%
    mutate(
      p_adj_method = "BH",
      p_adj_bh = p.adjust(p_value, method = "BH"),
      cliff_delta_magnitude = case_when(
        abs(cliff_delta_group1_vs_group2) < 0.147 ~ "negligible",
        abs(cliff_delta_group1_vs_group2) < 0.33 ~ "small",
        abs(cliff_delta_group1_vs_group2) < 0.474 ~ "medium",
        TRUE ~ "large"
      )
    ) %>%
    arrange(p_adj_bh, desc(abs(cliff_delta_group1_vs_group2)), group1, group2)
}

assert_file(mds_rds, "MDS coordinate RDS")
assert_file(seq_assignment_tsv, "Sequence assignment TSV")
assert_file(arm_assignment_tsv, "Arm assignment TSV")
assert_file(slide08a_asset, "Slide 08a reference PNG")

fit <- readRDS(mds_rds)
if (!is.list(fit) || is.null(fit$points) || is.null(fit$eig)) {
  stop("Expected full_mds.rds to contain a cmdscale-style list with points and eig", call. = FALSE)
}
if (ncol(fit$points) < 2) {
  stop("MDS point matrix must contain at least two dimensions", call. = FALSE)
}

coords <- as.data.frame(fit$points)
colnames(coords) <- paste0("D", seq_len(ncol(coords)))
coords$Name <- rownames(fit$points)
coordinate_cols <- pick_existing_dims(coords, paste0("D", 1:5))
var_explained <- fit$eig / sum(abs(fit$eig)) * 100

labels <- read_tsv(seq_assignment_tsv, col_types = cols(.default = col_character())) %>%
  select(Name, Sample, Haplotype, Chromosome, Arm, ChromArm, Superpopulation)

all_points <- coords %>%
  inner_join(labels, by = "Name") %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    Arm = factor(Arm, levels = c("p", "q")),
    ChromArm = factor(ChromArm)
  )

stopifnot("Every MDS row must have a matching label row" = nrow(all_points) == nrow(coords))
stopifnot("Every plotted point must have a continental superpopulation label" = !any(is.na(all_points$Superpopulation)))
stopifnot("Displayed MDS dimensions D1 and D2 must exist" = all(c("D1", "D2") %in% colnames(all_points)))

mds_x_limits <- range(all_points$D1, na.rm = TRUE)
mds_y_limits <- range(all_points$D2, na.rm = TRUE)

superpop_plot <- ggplot(
  all_points,
  aes(x = D1, y = D2, color = Superpopulation, shape = Arm)
) +
  geom_point(size = 2, alpha = 0.7) +
  xlab(paste0("Dimension 1 (", round(var_explained[1], 2), "%)")) +
  ylab(paste0("Dimension 2 (", round(var_explained[2], 2), "%)")) +
  ggtitle("hprcv2.1Mb.subtelo - Full MDS (distances preserved)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  scale_color_manual(values = superpop_colors, na.value = "grey50", name = "Superpopulation") +
  scale_shape_manual(values = arm_shapes, na.value = 1, name = "Arm") +
  guides(
    color = guide_legend(ncol = 2, override.aes = list(size = 5, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 5, alpha = 1, color = "grey30"))
  )

ggsave(
  file.path(out_dir, "superpopulation_mds_08a_matched.png"),
  superpop_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "superpopulation_mds_08a_matched.pdf"),
  superpop_plot,
  width = 12,
  height = 10,
  bg = "white"
)

nearest_2d <- nearest_by_superpop(all_points, c("D1", "D2"))
nearest_table <- all_points %>%
  select(Name, Sample, Haplotype, Chromosome, Arm, ChromArm, Superpopulation, all_of(coordinate_cols)) %>%
  mutate(
    Superpopulation = as.character(Superpopulation),
    Arm = as.character(Arm),
    ChromArm = as.character(ChromArm)
  ) %>%
  left_join(nearest_2d, by = "Name") %>%
  arrange(match(Superpopulation, superpop_order), Name)

stopifnot("Nearest-neighbor table must contain one row per MDS point" = nrow(nearest_table) == nrow(all_points))
stopifnot("The nearest-neighbor operation must exclude self" =
            !any(nearest_table$Name == nearest_table$nearest_same_superpop_name))
stopifnot("Nearest-neighbor operation must produce finite D1-D2 distances" =
            all(is.finite(nearest_table$nearest_same_superpop_distance_d1_d2)))

summary_table <- summarize_nearest(nearest_table)
pairwise_table <- pairwise_stats(nearest_table)
significance_table <- pairwise_table %>%
  select(group1, group2, n_group1, n_group2, wilcoxon_w, p_value, p_adj_method, p_adj_bh) %>%
  arrange(p_adj_bh, group1, group2)
effect_table <- pairwise_table %>%
  select(
    group1, group2, n_group1, n_group2,
    median_group1, median_group2, median_difference_group1_minus_group2,
    mean_group1, mean_group2, mean_difference_group1_minus_group2,
    cliff_delta_group1_vs_group2, cliff_delta_magnitude,
    p_adj_bh
  ) %>%
  arrange(desc(abs(cliff_delta_group1_vs_group2)), p_adj_bh, group1, group2)

kw <- kruskal.test(nearest_same_superpop_distance_d1_d2 ~ Superpopulation, data = nearest_table)
global_tests <- tibble(
  test = "Kruskal-Wallis rank-sum",
  metric = "nearest same-superpopulation Euclidean distance in displayed D1-D2 MDS space",
  statistic = unname(kw$statistic),
  parameter = unname(kw$parameter),
  p_value = kw$p.value,
  p_adj_method = "not applicable; single global test"
)

write_tsv(nearest_table, file.path(out_dir, "nearest_same_superpop_mds_distances.tsv"))
write_tsv(summary_table, file.path(out_dir, "nearest_same_superpop_mds_summary.tsv"))
write_tsv(significance_table, file.path(out_dir, "nearest_same_superpop_pairwise_wilcoxon.tsv"))
write_tsv(effect_table, file.path(out_dir, "nearest_same_superpop_effect_sizes.tsv"))
write_tsv(global_tests, file.path(out_dir, "nearest_same_superpop_global_tests.tsv"))

axis_label_vec <- summary_table %>%
  mutate(label = paste0(Superpopulation, "\nn=", comma(n_points))) %>%
  select(Superpopulation, label)
axis_label_vec <- setNames(axis_label_vec$label, axis_label_vec$Superpopulation)

plot_distances <- nearest_table %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    superpop_index = as.numeric(Superpopulation)
  )

summary_2d_plot <- summary_table %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    superpop_index = as.numeric(Superpopulation),
    median_label = paste0("median ", format_sci(median_nearest_distance, 2))
  )

top_contrasts <- pairwise_table %>%
  filter(p_adj_bh < 0.05) %>%
  arrange(desc(abs(cliff_delta_group1_vs_group2)), p_adj_bh, desc(abs(median_difference_group1_minus_group2))) %>%
  slice_head(n = 3)

contrast_lines <- if (nrow(top_contrasts) == 0) {
  "No BH-significant pairwise Wilcoxon contrasts"
} else {
  paste0(
    top_contrasts$group1,
    "-",
    top_contrasts$group2,
    ": d_med=",
    format_signed_sci(top_contrasts$median_difference_group1_minus_group2, 2),
    ", q=",
    format_p(top_contrasts$p_adj_bh),
    ", Cliff=",
    format_signed_num(top_contrasts$cliff_delta_group1_vs_group2, 2)
  )
}
annotation_label <- paste(
  paste0("Global KW p=", format_p(kw$p.value)),
  "Largest BH-significant effects:",
  paste(contrast_lines, collapse = "\n"),
  sep = "\n"
)

x_upper <- max(plot_distances$nearest_same_superpop_distance_d1_d2) * 1.35
x_upper <- max(x_upper, 0.12)

set.seed(8081)
distance_plot <- ggplot(
  plot_distances,
  aes(x = nearest_same_superpop_distance_d1_d2, y = superpop_index)
) +
  geom_boxplot(
    aes(group = Superpopulation, fill = Superpopulation),
    orientation = "y",
    width = 0.34,
    outlier.shape = NA,
    alpha = 0.62,
    linewidth = 0.36,
    color = "grey20"
  ) +
  geom_jitter(
    aes(color = Superpopulation),
    width = 0,
    height = 0.10,
    alpha = 0.16,
    size = 0.44,
    stroke = 0
  ) +
  geom_point(
    data = summary_2d_plot,
    aes(x = median_nearest_distance, y = superpop_index, color = Superpopulation),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    size = 2.65,
    stroke = 0.70
  ) +
  geom_label(
    data = summary_2d_plot,
    aes(x = x_upper, y = superpop_index, label = median_label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3.0,
    color = "grey15",
    fill = "white",
    alpha = 0.94,
    label.size = 0.12,
    label.padding = grid::unit(0.11, "lines")
  ) +
  annotate(
    "label",
    x = x_upper,
    y = 5.92,
    label = annotation_label,
    hjust = 1,
    vjust = 1,
    size = 3.15,
    lineheight = 0.94,
    color = "grey12",
    fill = "white",
    alpha = 0.96,
    label.size = 0.18,
    label.padding = grid::unit(0.15, "lines")
  ) +
  scale_fill_manual(values = superpop_colors, guide = "none") +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(superpop_order),
    labels = axis_label_vec[superpop_order],
    limits = c(0.55, 6.05),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-4),
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1"),
    expand = expansion(mult = c(0.01, 0.035))
  ) +
  coord_cartesian(xlim = c(0, x_upper), clip = "off") +
  labs(
    title = "Nearest same-superpopulation neighbor in displayed MDS space",
    subtitle = "Boxplots plus jitter replace the long-tailed violins. Each point contributes one nearest other same-superpopulation D1-D2 distance.",
    x = "Nearest same-superpopulation MDS distance (D1-D2 Euclidean)",
    y = NULL,
    caption = paste(
      "Metric: nearest other point from the same continental superpopulation in displayed D1-D2 MDS space; self excluded.",
      "No centroid, no all-pairwise summary. P-values: Kruskal-Wallis plus pairwise Wilcoxon, BH corrected; effect size: Cliff's delta.",
      sep = "\n"
    )
  ) +
  theme_classic(base_size = 12.5) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.4, color = "grey25"),
    plot.caption = element_text(size = 8.25, color = "grey35", hjust = 0),
    axis.title = element_text(face = "bold", size = 11.2),
    axis.text = element_text(size = 10.2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 18, 10, 12)
  )

ggsave(
  file.path(out_dir, "nearest_same_superpop_distance_boxplot.png"),
  distance_plot,
  width = 12.2,
  height = 6.7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "nearest_same_superpop_distance_boxplot.pdf"),
  distance_plot,
  width = 12.2,
  height = 6.7,
  bg = "white"
)

arm_assignments <- read_tsv(arm_assignment_tsv, col_types = cols(.default = col_character()))
required_assignment_cols <- c("ChromArm", "Community", "Arms")
missing_cols <- setdiff(required_assignment_cols, names(arm_assignments))
if (length(missing_cols) > 0) {
  stop("Arm assignment TSV is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

community_points <- all_points %>%
  transmute(
    Name,
    x = D1,
    y = D2,
    ChromArm = as.character(ChromArm),
    Arm = factor(as.character(Arm), levels = c("p", "q"))
  ) %>%
  left_join(arm_assignments[, required_assignment_cols], by = "ChromArm")

if (any(is.na(community_points$Community))) {
  missing_arms <- sort(unique(community_points$ChromArm[is.na(community_points$Community)]))
  stop("Missing community assignments for arm(s): ", paste(missing_arms, collapse = ", "), call. = FALSE)
}

observed_communities <- sort(unique(as.character(community_points$Community)), method = "radix")
observed_communities <- observed_communities[order(community_number(observed_communities))]
if (!identical(observed_communities, community_levels)) {
  stop("Expected communities C1-C15; observed: ", paste(observed_communities, collapse = ", "), call. = FALSE)
}

community_points$Community <- factor(community_points$Community, levels = community_levels)

centroids <- aggregate(cbind(x, y) ~ Community, data = community_points, FUN = median)
centroids$Community <- as.character(centroids$Community)
names(centroids)[names(centroids) == "x"] <- "anchor_x"
names(centroids)[names(centroids) == "y"] <- "anchor_y"

point_counts <- aggregate(Name ~ Community, data = community_points, FUN = length)
point_counts$Community <- as.character(point_counts$Community)
names(point_counts)[names(point_counts) == "Name"] <- "n_points"

arm_counts <- aggregate(ChromArm ~ Community, data = arm_assignments, FUN = length)
arm_counts$Community <- as.character(arm_counts$Community)
names(arm_counts)[names(arm_counts) == "ChromArm"] <- "n_arms"

community_arms <- unique(arm_assignments[, c("Community", "Arms")])
community_arms$Community <- as.character(community_arms$Community)
community_arms$num <- community_number(community_arms$Community)
community_arms <- community_arms[order(community_arms$num), ]

manual_label_positions <- tibble::tribble(
  ~Community, ~label_x, ~label_y,
  "C1", 0.52, 0.17,
  "C2", -0.18, -0.55,
  "C3", -0.50, 0.27,
  "C4", 0.05, -0.55,
  "C5", -0.42, -0.55,
  "C6", 0.46, 0.50,
  "C7", 0.29, -0.55,
  "C8", -0.50, -0.09,
  "C9", 0.52, -0.31,
  "C10", 0.52, 0.06,
  "C11", -0.50, 0.47,
  "C12", -0.50, 0.09,
  "C13", 0.52, -0.44,
  "C14", 0.52, -0.07,
  "C15", 0.52, -0.19
)

label_df <- Reduce(
  function(left, right) merge(left, right, by = "Community", sort = FALSE),
  list(manual_label_positions, community_arms[, c("Community", "Arms")], centroids, point_counts, arm_counts)
)
label_df$num <- community_number(label_df$Community)
label_df <- label_df[order(label_df$num), ]
if (!identical(label_df$Community, community_levels)) {
  stop("Label data must contain exactly C1-C15 in numeric order", call. = FALSE)
}

label_df$label <- paste0(label_df$Community, "\n", vapply(label_df$Arms, wrap_arms, character(1)))
label_df$label_position_strategy <- "manual edge labels with leader lines; direct labels replace dense community legend"
label_df$Community <- factor(label_df$Community, levels = community_levels)

community_axis_limits <- c(-0.68, 0.68)
community_axis_breaks <- seq(-0.5, 0.5, by = 0.25)

community_plot <- ggplot(community_points, aes(x = x, y = y, color = Community, shape = Arm)) +
  geom_point(size = 1.75, alpha = 0.58, stroke = 0) +
  geom_segment(
    data = label_df,
    aes(x = anchor_x, y = anchor_y, xend = label_x, yend = label_y, color = Community),
    inherit.aes = FALSE,
    linewidth = 0.52,
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = label_df,
    aes(x = anchor_x, y = anchor_y, color = Community),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2.55,
    stroke = 0.52,
    alpha = 0.96
  ) +
  geom_label(
    data = label_df,
    aes(x = label_x, y = label_y, label = label, fill = Community),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 3.25,
    lineheight = 0.88,
    label.size = 0.38,
    label.r = grid::unit(0.12, "lines")
  ) +
  scale_color_manual(values = community_colors, guide = "none") +
  scale_fill_manual(
    values = setNames(adjustcolor(community_colors, alpha.f = 0.18), names(community_colors)),
    guide = "none"
  ) +
  scale_shape_manual(values = arm_shapes, name = "Arm") +
  scale_x_continuous(limits = community_axis_limits, breaks = community_axis_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = community_axis_limits, breaks = community_axis_breaks, expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = community_axis_limits, ylim = community_axis_limits, clip = "off") +
  labs(
    title = "hprcv2.1Mb.subtelo - Full MDS colored by Leiden community",
    subtitle = paste(
      "Classical MDS on 1 - Jaccard; C1-C15 are arm-level Leiden communities",
      "Direct labels use non-overlapping manual positions with leader lines; not PCA",
      sep = "\n"
    ),
    x = sprintf("Dimension 1 (%.2f%%)", var_explained[1]),
    y = sprintf("Dimension 2 (%.2f%%)", var_explained[2])
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 23),
    plot.subtitle = element_text(hjust = 0.5, size = 12.1, lineheight = 0.95),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12.4),
    legend.key.height = grid::unit(0.33, "in"),
    legend.key.width = grid::unit(0.32, "in"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.34),
    panel.border = element_rect(color = "#303030", fill = NA, linewidth = 0.55),
    plot.margin = margin(9, 13, 9, 12)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 5, alpha = 1, color = "#4D4D4D")))

ggsave(
  file.path(out_dir, "community_mds_labeled.png"),
  community_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "community_mds_labeled.pdf"),
  community_plot,
  width = 12,
  height = 10,
  bg = "white"
)

label_table <- label_df
label_table$Community <- as.character(label_table$Community)
label_table$label <- gsub("\n", " | ", label_table$label, fixed = TRUE)
write_tsv(
  label_table[, c(
    "Community", "Arms", "n_arms", "n_points",
    "anchor_x", "anchor_y", "label_x", "label_y",
    "label", "label_position_strategy"
  )],
  file.path(out_dir, "community_mds_label_positions.tsv")
)

source_manifest <- tibble(
  source = c(
    "mds_rds",
    "seq_assignment_tsv",
    "arm_assignment_tsv",
    "upstream_superpop_png",
    "upstream_chromosome_png",
    "slide08a_asset",
    "slide08b_asset",
    "source_script"
  ),
  path = c(
    mds_rds,
    seq_assignment_tsv,
    arm_assignment_tsv,
    upstream_superpop_png,
    upstream_chromosome_png,
    slide08a_asset,
    slide08b_asset,
    source_script
  )
) %>%
  mutate(exists = file.exists(path), sha256 = vapply(path, sha256, character(1)))
write_tsv(source_manifest, file.path(out_dir, "source_manifest.tsv"))

validation_summary <- tibble(
  check = c(
    "coordinate_method",
    "coordinate_source",
    "point_count",
    "available_mds_dimensions",
    "displayed_dimensions",
    "d1_range",
    "d2_range",
    "d1_percent",
    "d2_percent",
    "slide08b_superpop_plot_inches",
    "slide08b_superpop_png_pixels",
    "slide08b_superpop_style_reference",
    "nearest_distance_metric",
    "nearest_distance_self_excluded",
    "nearest_distance_global_test",
    "nearest_distance_pairwise_test",
    "nearest_distance_effect_size",
    "slide09_community_plot_inches",
    "slide09_community_png_pixels",
    "slide09_coord_fixed_ratio",
    "slide09_axis_limit_x",
    "slide09_axis_limit_y",
    "slide09_labels_present",
    "slide09_pca_used"
  ),
  value = c(
    "Classical MDS / cmdscale on 1 - graph-path Jaccard distances",
    mds_rds,
    as.character(nrow(all_points)),
    as.character(length(coordinate_cols)),
    "D1,D2",
    paste(format_num(mds_x_limits, 6), collapse = ","),
    paste(format_num(mds_y_limits, 6), collapse = ","),
    format_num(var_explained[1], 4),
    format_num(var_explained[2], 4),
    "12 x 10",
    "3600 x 3000",
    "slides/v2-review-zoom/_typst/assets/s08a_mds_chrom.png; upstream 12x10 theme_bw MDS grammar, size 2 points, alpha 0.7, right legend",
    "For each point, Euclidean distance in displayed D1-D2 MDS space to nearest other point from same continental superpopulation",
    as.character(!any(nearest_table$Name == nearest_table$nearest_same_superpop_name)),
    paste0("Kruskal-Wallis p=", format_p(kw$p.value)),
    "Pairwise Wilcoxon rank-sum, exact=FALSE, continuity correction disabled, BH adjusted",
    "Median difference and Cliff's delta / rank-biserial from Mann-Whitney U",
    "12 x 10",
    "3600 x 3000",
    "1",
    paste(community_axis_limits, collapse = ","),
    paste(community_axis_limits, collapse = ","),
    paste(as.character(label_df$Community), collapse = ","),
    "no"
  )
)
write_tsv(validation_summary, file.path(out_dir, "validation_summary.tsv"))

summary_rows <- summary_table %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", comma(n_points),
      " | ", comma(n_samples),
      " | ", format_sci(median_nearest_distance, 2),
      " | ", format_sci(mean_nearest_distance, 2),
      " | ", format_sci(q25_nearest_distance, 2), "-", format_sci(q75_nearest_distance, 2),
      " | ", format_sci(q95_nearest_distance, 2),
      " | ", format_sci(max_nearest_distance, 2),
      " |"
    )
  ) %>%
  pull(row)

top_effect_rows <- effect_table %>%
  filter(p_adj_bh < 0.05) %>%
  arrange(desc(abs(cliff_delta_group1_vs_group2)), p_adj_bh) %>%
  slice_head(n = 5) %>%
  transmute(
    row = paste0(
      "| ", group1, " vs ", group2,
      " | ", format_signed_sci(median_difference_group1_minus_group2, 2),
      " | ", format_signed_num(cliff_delta_group1_vs_group2, 3),
      " | ", format_p(p_adj_bh),
      " |"
    )
  ) %>%
  pull(row)

if (length(top_effect_rows) == 0) {
  top_effect_rows <- "| none | NA | NA | NA |"
}

readme <- c(
  "# V8 MDS Superpopulation and Community Polish",
  "",
  "Task: `review-zoom-v8-mds-superpop-community-polish`.",
  "",
  "## Deliverables",
  "",
  "- `superpopulation_mds_08a_matched.png` / `.pdf`: slide 08b replacement, colored by continental superpopulation but rendered in the 12 x 10 slide 08a MDS grammar.",
  "- `nearest_same_superpop_distance_boxplot.png` / `.pdf`: slide 08b.1 replacement using boxplots plus jitter and concise global/pairwise statistical annotations.",
  "- `nearest_same_superpop_mds_distances.tsv`: raw per-point nearest same-superpopulation D1-D2 MDS distance table.",
  "- `nearest_same_superpop_mds_summary.tsv`: robust per-superpopulation summaries.",
  "- `nearest_same_superpop_pairwise_wilcoxon.tsv`: pairwise Wilcoxon rank-sum tests with BH correction.",
  "- `nearest_same_superpop_effect_sizes.tsv`: pairwise median differences and Cliff's delta effect sizes.",
  "- `nearest_same_superpop_global_tests.tsv`: Kruskal-Wallis global test.",
  "- `community_mds_labeled.png` / `.pdf`: slide 09 replacement, colored by C1-C15 arm-level Leiden community with non-overlapping direct labels and leader lines.",
  "- `community_mds_label_positions.tsv`: explicit community label anchors and label positions.",
  "- `source_manifest.tsv`, `validation_summary.tsv`, `VALIDATION.md`, and `SLIDE_PATCH.md`: provenance, checks, and fan-in instructions.",
  "",
  "## Metric",
  "",
  "The slide 08b.1 metric is exactly the nearest same-superpopulation neighbor distance in the displayed MDS panel. For each sequence-level subtelomeric point, the script computes Euclidean distances in dimensions D1 and D2 to all other points from the same continental superpopulation, excludes self by setting the diagonal to `Inf`, and retains the minimum.",
  "",
  "This is not a centroid distance, not a within-superpopulation all-pairwise distribution, and not an average against all in-group points.",
  "",
  "## Main 2D Summary",
  "",
  "| Superpop | MDS points | Samples | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Max |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  summary_rows,
  "",
  paste0("Global Kruskal-Wallis p-value: `", format_p(kw$p.value), "`."),
  "",
  "Largest BH-significant pairwise effects:",
  "",
  "| Contrast | Median difference | Cliff's delta | BH q |",
  "|---|---:|---:|---:|",
  top_effect_rows,
  "",
  "## Reproduction",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript slides/v2-review-zoom/_revision_assets/v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R",
  "```"
)
write_lines(readme, file.path(out_dir, "README.md"))

validation <- c(
  "# Validation",
  "",
  paste0("Generated by `make_mds_superpop_community_polish.R` on ", format(Sys.Date(), "%Y-%m-%d"), " UTC."),
  "",
  "## Source Checks",
  "",
  paste0("- MDS coordinate rows: `", nrow(all_points), "`."),
  paste0("- Coordinate dimensions available in cached RDS and exported here: `", length(coordinate_cols), "`."),
  paste0("- D1 range: `", paste(format_num(mds_x_limits, 6), collapse = ", "), "`."),
  paste0("- D2 range: `", paste(format_num(mds_y_limits, 6), collapse = ", "), "`."),
  paste0("- D1/D2 percentages: `", format_num(var_explained[1], 2), "%`, `", format_num(var_explained[2], 2), "%`."),
  "- Superpopulation labels: no missing labels after joining the sequence-level assignment TSV.",
  "- Slide 08b replacement is saved at 12 x 10 inches / 3600 x 3000 px to match slide 08a's upstream MDS render format.",
  "- Slide 09 replacement is saved at 12 x 10 inches / 3600 x 3000 px, with `coord_fixed(ratio = 1)` and identical x/y limits.",
  "",
  "## Nearest-Distance Metric Checks",
  "",
  "- Main metric uses Euclidean distance in displayed MDS dimensions `D1` and `D2`.",
  "- For each superpopulation, the script builds the within-superpopulation distance matrix, sets the diagonal to `Inf`, and keeps one minimum distance per point.",
  "- Self matches are explicitly rejected in the generated table.",
  "- There is exactly one nearest-neighbor row per MDS point.",
  "- No centroid distance, all-pairwise summary, or in-group average is used.",
  "",
  "## Statistical Checks",
  "",
  paste0("- Global Kruskal-Wallis p-value: `", format_p(kw$p.value), "`."),
  "- Pairwise tests are Wilcoxon rank-sum tests with BH correction.",
  "- Effect sizes are median differences plus Cliff's delta / rank-biserial effect sizes from the Mann-Whitney U statistic.",
  "",
  "## Community Label Checks",
  "",
  "- Community coordinates are the same cached MDS D1-D2 coordinates, not PCA.",
  "- C1-C15 are direct labels with leader lines and manual edge positions recorded in `community_mds_label_positions.tsv`.",
  "- The plot subtitle and recommended deck footer both state MDS and not PCA."
)
write_lines(validation, file.path(out_dir, "VALIDATION.md"))

slide_patch <- c(
  "# Slide Patch Recommendation",
  "",
  "Task: `review-zoom-v8-mds-superpop-community-polish`.",
  "",
  "Do not edit the Typst deck in this task. These are the recommended replacements for `review-zoom-v8-fanin-render`.",
  "",
  "## Slide 08b",
  "",
  "Replace the current slide 08b figure asset with:",
  "",
  "```typst",
  "#figure-slide(",
  "  \"08b\",",
  "  \"MDS colored by continental superpopulation\",",
  "  \"../_revision_assets/v8/mds_superpop_community_polish/superpopulation_mds_08a_matched.png\",",
  "  source: \"v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS D1-D2 from hprcv2.1Mb.subtelo.full_mds.rds; 08a-matched 12x10 MDS render style; color is continental superpopulation\",",
  ")",
  "```",
  "",
  "Exact caption/source line to use if the slide is captioned outside `figure-slide`:",
  "",
  "`MDS D1-D2 from hprcv2.1Mb.subtelo.full_mds.rds; rendered in the same 12 x 10 theme_bw MDS format as slide 08a; color is continental superpopulation and shape is p/q arm.`",
  "",
  "## Slide 08b.1",
  "",
  "Replace the current slide 08b.1 figure asset with:",
  "",
  "```typst",
  "#figure-slide(",
  "  \"08b.1\",",
  "  \"Nearest same-superpopulation neighbor in MDS space\",",
  "  \"../_revision_assets/v8/mds_superpop_community_polish/nearest_same_superpop_distance_boxplot.png\",",
  "  source: \"v8/mds_superpop_community_polish/nearest_same_superpop_mds_distances.tsv; metric is nearest other same-superpopulation Euclidean distance in displayed D1-D2 MDS, self excluded; Wilcoxon BH and Cliff delta TSVs\",",
  ")",
  "```",
  "",
  "Exact caption/source line:",
  "",
  "`For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Self is excluded; this is not a centroid metric and not an all-pairwise within-population distance. Pairwise Wilcoxon tests use BH correction; effect sizes are median differences and Cliff's delta.`",
  "",
  "Traceability files for slide 08b.1:",
  "",
  "- `nearest_same_superpop_mds_distances.tsv`",
  "- `nearest_same_superpop_mds_summary.tsv`",
  "- `nearest_same_superpop_pairwise_wilcoxon.tsv`",
  "- `nearest_same_superpop_effect_sizes.tsv`",
  "- `nearest_same_superpop_global_tests.tsv`",
  "",
  "## Slide 09",
  "",
  "Replace the current slide 09 figure asset with:",
  "",
  "```typst",
  "#figure-slide(",
  "  \"09\",",
  "  \"MDS: all Leiden communities C1-C15 labeled\",",
  "  \"../_revision_assets/v8/mds_superpop_community_polish/community_mds_labeled.png\",",
  "  source: \"v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS from hprcv2.1Mb.subtelo.full_mds.rds; C1-C15 from arm-level graph-path Jaccard Leiden assignments; direct labels from community_mds_label_positions.tsv; 1:1 axes; not PCA\",",
  ")",
  "```",
  "",
  "Exact caption/source line:",
  "",
  "`Classical MDS on 1 - graph-path Jaccard distances; C1-C15 are arm-level Leiden communities. Direct community labels and leader lines use positions recorded in community_mds_label_positions.tsv. This is MDS, not PCA.`",
  "",
  "## Required Caveats",
  "",
  "- Slide 08b.1 must retain the nearest-neighbor wording: nearest other same-superpopulation point, self excluded.",
  "- Do not describe slide 08b.1 as centroid distance, RMS radius, all-pairwise distance, or average in-group distance.",
  "- Slide 09 must be labeled as MDS in title/footer/caption; do not call it PCA."
)
write_lines(slide_patch, file.path(out_dir, "SLIDE_PATCH.md"))

message("Wrote v8 MDS polish assets to ", out_dir)
