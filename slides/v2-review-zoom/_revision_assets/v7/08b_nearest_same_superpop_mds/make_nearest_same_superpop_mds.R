#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

out_dir <- "slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mds_rds <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds"
label_export_tsv <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv"
sample_metadata_tsv <- "/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv"
source_png <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png"
current_slide_png <- "slides/v2-review-zoom/_typst/assets/s08b_mds_superpop.png"
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

format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
format_sci <- function(x, digits = 2) formatC(x, format = "e", digits = digits)
format_ci <- function(lo, hi, digits = 2) paste0(format_sci(lo, digits), "-", format_sci(hi, digits))

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  line <- system2("sha256sum", path, stdout = TRUE)
  strsplit(line, "[[:space:]]+")[[1]][1]
}

pick_existing_dims <- function(dat, dims) {
  dims[dims %in% colnames(dat)]
}

nearest_by_superpop <- function(points, dims, suffix) {
  if (length(dims) < 2) {
    stop("Nearest-neighbor distance requires at least two MDS dimensions.")
  }

  out <- lapply(superpop_order, function(superpop) {
    dat <- points %>%
      filter(Superpopulation == superpop) %>%
      arrange(Name)

    if (nrow(dat) < 2) {
      stop(sprintf("Superpopulation %s has fewer than two points.", superpop))
    }

    mat <- as.matrix(dat[, dims, drop = FALSE])
    dm <- as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
    diag(dm) <- Inf

    nearest_idx <- max.col(-dm, ties.method = "first")
    nearest_dist <- dm[cbind(seq_len(nrow(dm)), nearest_idx)]
    nearest <- dat[nearest_idx, , drop = FALSE]

    tibble(
      Name = dat$Name,
      !!paste0("nearest_same_superpop_name_", suffix) := nearest$Name,
      !!paste0("nearest_same_superpop_distance_", suffix) := nearest_dist,
      !!paste0("nearest_same_superpop_sample_", suffix) := nearest$Sample,
      !!paste0("nearest_same_superpop_chromarm_", suffix) := nearest$ChromArm
    )
  })

  bind_rows(out)
}

summarize_nearest <- function(dat, distance_col, dimensionality) {
  dat %>%
    group_by(Superpopulation) %>%
    summarise(
      dimensionality = dimensionality,
      n_points = n(),
      n_samples = n_distinct(Sample),
      min_nearest_distance = min(.data[[distance_col]]),
      q05_nearest_distance = unname(quantile(.data[[distance_col]], 0.05)),
      q25_nearest_distance = unname(quantile(.data[[distance_col]], 0.25)),
      median_nearest_distance = median(.data[[distance_col]]),
      mean_nearest_distance = mean(.data[[distance_col]]),
      q75_nearest_distance = unname(quantile(.data[[distance_col]], 0.75)),
      q95_nearest_distance = unname(quantile(.data[[distance_col]], 0.95)),
      q99_nearest_distance = unname(quantile(.data[[distance_col]], 0.99)),
      max_nearest_distance = max(.data[[distance_col]]),
      sd_nearest_distance = sd(.data[[distance_col]]),
      n_exact_zero = sum(.data[[distance_col]] == 0),
      n_le_1e_minus_7 = sum(.data[[distance_col]] <= 1e-7),
      .groups = "drop"
    ) %>%
    mutate(
      Superpopulation = as.character(Superpopulation),
      dimensionality = as.character(dimensionality)
    ) %>%
    arrange(match(Superpopulation, superpop_order))
}

mds <- readRDS(mds_rds)
coords <- as.data.frame(mds$points)
colnames(coords) <- paste0("D", seq_len(ncol(coords)))
coords$Name <- rownames(mds$points)

labels <- read_tsv(label_export_tsv, col_types = cols(.default = col_character())) %>%
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

var_explained <- mds$eig / sum(abs(mds$eig)) * 100
source_png_sha <- sha256(source_png)
current_slide_png_sha <- sha256(current_slide_png)
if (!is.na(source_png_sha) && !is.na(current_slide_png_sha)) {
  stopifnot("Current slide 08b PNG should match upstream source PNG" = identical(source_png_sha, current_slide_png_sha))
}

nearest_2d <- nearest_by_superpop(all_points, c("D1", "D2"), "2d")
five_d_dims <- pick_existing_dims(all_points, paste0("D", 1:5))
nearest_5d <- nearest_by_superpop(all_points, five_d_dims, "5d")

coordinate_cols <- pick_existing_dims(all_points, paste0("D", 1:5))
nearest_table <- all_points %>%
  select(Name, Sample, Haplotype, Chromosome, Arm, ChromArm, Superpopulation, all_of(coordinate_cols)) %>%
  mutate(
    Superpopulation = as.character(Superpopulation),
    Arm = as.character(Arm),
    ChromArm = as.character(ChromArm)
  ) %>%
  left_join(nearest_2d, by = "Name") %>%
  left_join(nearest_5d, by = "Name") %>%
  arrange(match(Superpopulation, superpop_order), Name)

stopifnot("Nearest-neighbor table must contain one row per MDS point" = nrow(nearest_table) == nrow(all_points))
stopifnot("The nearest-neighbor operation must exclude self" =
            !any(nearest_table$Name == nearest_table$nearest_same_superpop_name_2d))
stopifnot("Nearest-neighbor operation must produce finite D1-D2 distances" =
            all(is.finite(nearest_table$nearest_same_superpop_distance_2d)))

summary_2d <- summarize_nearest(
  nearest_table,
  "nearest_same_superpop_distance_2d",
  "2D displayed MDS"
)
summary_5d <- summarize_nearest(
  nearest_table,
  "nearest_same_superpop_distance_5d",
  "5D cached MDS sensitivity"
)
summary_table <- bind_rows(summary_2d, summary_5d) %>%
  arrange(match(dimensionality, c("2D displayed MDS", "5D cached MDS sensitivity")),
          match(Superpopulation, superpop_order))

write_tsv(nearest_table, file.path(out_dir, "nearest_same_superpop_mds_distances.tsv"))
write_tsv(summary_table, file.path(out_dir, "nearest_same_superpop_mds_summary.tsv"))

source_manifest <- tibble(
  source = c("mds_rds", "label_export_tsv", "sample_metadata_tsv", "source_png", "current_slide_png", "source_script"),
  path = c(mds_rds, label_export_tsv, sample_metadata_tsv, source_png, current_slide_png, source_script),
  exists = file.exists(path),
  sha256 = vapply(path, sha256, character(1))
)
write_tsv(source_manifest, file.path(out_dir, "source_manifest.tsv"))

point_counts <- all_points %>%
  count(Superpopulation, name = "n_points") %>%
  mutate(Superpopulation = as.character(Superpopulation))
sample_counts <- all_points %>%
  distinct(Sample, Superpopulation) %>%
  count(Superpopulation, name = "n_samples") %>%
  mutate(Superpopulation = as.character(Superpopulation))

legend_labels <- point_counts %>%
  left_join(sample_counts, by = "Superpopulation") %>%
  mutate(label = paste0(Superpopulation, " (", comma(n_points), " points)")) %>%
  select(Superpopulation, label)
legend_label_vec <- setNames(legend_labels$label, legend_labels$Superpopulation)

axis_label_vec <- point_counts %>%
  mutate(label = paste0(Superpopulation, "\nn=", comma(n_points))) %>%
  select(Superpopulation, label)
axis_label_vec <- setNames(axis_label_vec$label, axis_label_vec$Superpopulation)

scatter_data <- all_points %>%
  mutate(Superpopulation = factor(Superpopulation, levels = superpop_order))

scatter_plot <- ggplot(scatter_data, aes(x = D1, y = D2, color = Superpopulation, shape = Arm)) +
  geom_point(size = 1.65, alpha = 0.72, stroke = 0.12) +
  scale_color_manual(
    values = superpop_colors,
    breaks = superpop_order,
    labels = legend_label_vec[superpop_order],
    name = "Superpopulation"
  ) +
  scale_shape_manual(values = arm_shapes, name = "Arm", na.value = 1) +
  coord_fixed(ratio = 1) +
  labs(
    title = "hprcv2.1Mb.subtelo - Full MDS (distances preserved)",
    subtitle = "Color = continental superpopulation; shape = chromosome arm. Displayed axes are D1-D2 from classical MDS / PCoA.",
    x = paste0("Dimension 1 (", format_num(var_explained[1], 2), "%)"),
    y = paste0("Dimension 2 (", format_num(var_explained[2], 2), "%)"),
    caption = "Same source coordinates and 1000 Genomes-style colors as slide 08b; panel is rendered with 1:1 MDS-axis scaling."
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 1, override.aes = list(size = 4.2, alpha = 1, shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(size = 4.2, alpha = 1, color = "grey30"))
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 11.5, color = "grey25"),
    plot.caption = element_text(size = 8.7, color = "grey35", hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 13.5),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14.5),
    axis.text = element_text(size = 12.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    plot.margin = margin(9, 14, 9, 12)
  )

ggsave(
  file.path(out_dir, "superpopulation_mds_original_style.png"),
  scatter_plot,
  width = 11.6,
  height = 7.2,
  dpi = 300
)
ggsave(
  file.path(out_dir, "superpopulation_mds_original_style.pdf"),
  scatter_plot,
  width = 11.6,
  height = 7.2
)

plot_distances <- nearest_table %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    superpop_axis = factor(Superpopulation, levels = superpop_order)
  )

summary_2d_plot <- summary_2d %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    median_label = paste0("median ", format_sci(median_nearest_distance, 2))
  )
x_upper <- max(plot_distances$nearest_same_superpop_distance_2d) * 1.18
x_upper <- max(x_upper, 0.12)

set.seed(807)
distance_plot <- ggplot(
  plot_distances,
  aes(x = nearest_same_superpop_distance_2d, y = Superpopulation)
) +
  geom_violin(
    aes(fill = Superpopulation),
    orientation = "y",
    alpha = 0.20,
    width = 0.82,
    trim = FALSE,
    linewidth = 0.28,
    color = "grey35"
  ) +
  geom_boxplot(
    aes(fill = Superpopulation),
    orientation = "y",
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.62,
    linewidth = 0.34,
    color = "grey20"
  ) +
  geom_jitter(
    aes(color = Superpopulation),
    width = 0,
    height = 0.10,
    alpha = 0.20,
    size = 0.48,
    stroke = 0
  ) +
  geom_point(
    data = summary_2d_plot,
    aes(x = median_nearest_distance, y = Superpopulation, color = Superpopulation),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    size = 2.5,
    stroke = 0.65
  ) +
  geom_label(
    data = summary_2d_plot,
    aes(
      x = x_upper,
      y = Superpopulation,
      label = median_label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3.0,
    color = "grey15",
    fill = "white",
    alpha = 0.92,
    label.size = 0.12,
    label.padding = grid::unit(0.11, "lines")
  ) +
  scale_fill_manual(values = superpop_colors, guide = "none") +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_y_discrete(labels = axis_label_vec[superpop_order]) +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-4),
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1"),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  coord_cartesian(xlim = c(0, x_upper)) +
  labs(
    title = "Nearest same-superpopulation neighbor in displayed MDS space",
    subtitle = "Each sequence-level subtelomeric point contributes one D1-D2 Euclidean distance to its nearest other point from the same continental superpopulation.",
    x = "Nearest same-superpopulation MDS distance (D1-D2)",
    y = NULL,
    caption = "Lower values mean each subtelomere has a closer same-superpopulation neighbor in sequence-similarity MDS space.\nSelf is excluded; no centroid/all-pairwise average. Pseudo-log x-axis shows near-zero values and sparse outliers."
  ) +
  theme_classic(base_size = 12.5) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.5, color = "grey25"),
    plot.caption = element_text(size = 8.3, color = "grey35", hjust = 0),
    axis.title = element_text(face = "bold", size = 11.2),
    axis.text = element_text(size = 10.2),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 18, 10, 12)
  )

ggsave(
  file.path(out_dir, "nearest_same_superpop_distance_distribution.png"),
  distance_plot,
  width = 12.2,
  height = 6.7,
  dpi = 300
)
ggsave(
  file.path(out_dir, "nearest_same_superpop_distance_distribution.pdf"),
  distance_plot,
  width = 12.2,
  height = 6.7
)

ordered_2d <- summary_2d %>%
  arrange(match(Superpopulation, superpop_order))
ordered_5d <- summary_5d %>%
  arrange(match(Superpopulation, superpop_order))

results_rows <- ordered_2d %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", comma(n_points),
      " | ", comma(n_samples),
      " | ", format_sci(median_nearest_distance, 2),
      " | ", format_sci(mean_nearest_distance, 2),
      " | ", format_ci(q25_nearest_distance, q75_nearest_distance, 2),
      " | ", format_sci(q95_nearest_distance, 2),
      " | ", format_sci(q99_nearest_distance, 2),
      " | ", format_sci(max_nearest_distance, 2),
      " | ", comma(n_exact_zero),
      " |"
    )
  ) %>%
  pull(row)

sensitivity_rows <- ordered_5d %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", format_sci(median_nearest_distance, 2),
      " | ", format_sci(mean_nearest_distance, 2),
      " | ", format_ci(q25_nearest_distance, q75_nearest_distance, 2),
      " | ", format_sci(q95_nearest_distance, 2),
      " | ", format_sci(max_nearest_distance, 2),
      " |"
    )
  ) %>%
  pull(row)

sample_count_text <- sample_counts %>%
  arrange(match(Superpopulation, superpop_order)) %>%
  mutate(label = paste0(Superpopulation, "=", n_samples)) %>%
  pull(label) %>%
  paste(collapse = ", ")

point_count_text <- point_counts %>%
  arrange(match(Superpopulation, superpop_order)) %>%
  mutate(label = paste0(Superpopulation, "=", n_points)) %>%
  pull(label) %>%
  paste(collapse = ", ")

readme <- c(
  "# Slide 08b Nearest Same-Superpopulation MDS Distance",
  "",
  paste0("Task: `review-zoom-v7-slide08b-nearest-superpop-mds`."),
  "",
  "## Deliverables",
  "",
  "- `superpopulation_mds_original_style.png` / `.pdf`: replacement source-style MDS scatter plot using the original D1-D2 coordinates, original continental superpopulation colors, p/q arm shapes, and 1:1 MDS-axis scaling.",
  "- `nearest_same_superpop_distance_distribution.png` / `.pdf`: violin/box/jitter plot of each point's nearest same-superpopulation neighbor distance in displayed D1-D2 MDS space.",
  "- `nearest_same_superpop_mds_distances.tsv`: row-level nearest-neighbor table, including each point's 2D nearest same-superpopulation neighbor and a secondary 5D sensitivity nearest neighbor.",
  "- `nearest_same_superpop_mds_summary.tsv`: per-superpopulation summaries for the 2D slide metric plus the optional 5D sensitivity check.",
  "- `source_manifest.tsv`: input paths, existence checks, and SHA-256 hashes where available.",
  "- `VALIDATION.md`: explicit validation note for metric choice and source checks.",
  "- `SLIDE_PATCH.md`: recommended deck insertion guidance for the v7 fan-in renderer.",
  "- `make_nearest_same_superpop_mds.R`: reproducible generator for all files in this directory.",
  "",
  "## Metric",
  "",
  "The slide metric is nearest same-superpopulation neighbor distance in the displayed MDS panel.",
  "",
  "For each sequence-level subtelomeric MDS point with a continental superpopulation label, the script computes Euclidean distance in dimensions D1 and D2 to every other point from the same superpopulation, excludes the point itself by setting the diagonal distance to `Inf`, and keeps the minimum distance and nearest-neighbor identity.",
  "",
  "This metric is deliberately not a centroid distance, not a within-population all-pairwise distribution, and not an average against all in-group points. Pairwise distances are only used internally to identify each point's nearest eligible neighbor.",
  "",
  "Lower values mean each subtelomere has a closer same-superpopulation neighbor in sequence-similarity MDS space.",
  "",
  "## Source Confirmation",
  "",
  paste0("- Coordinate source: `", mds_rds, "` (`", nrow(all_points), "` rows x `", length(coordinate_cols), "` exported dimensions used here; cached R object stores `", ncol(coords) - 1, "` MDS dimensions)."),
  paste0("- Direct row-label source: `", label_export_tsv, "`."),
  paste0("- Original sample metadata source used by the upstream plot script: `", sample_metadata_tsv, "`."),
  paste0("- Original plotting script: `", source_script, "`; the MDS call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`."),
  paste0("- Current slide 08b asset: `", current_slide_png, "`."),
  paste0("- Exact upstream PNG: `", source_png, "`."),
  paste0("- Upstream PNG SHA-256: `", source_png_sha, "`."),
  paste0("- Current slide PNG SHA-256: `", current_slide_png_sha, "`; the script asserts identity when both files are present."),
  "",
  paste0("Terminology: this is classical MDS / PCoA on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses D1 and D2. Axis percentages are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = ", format_num(var_explained[1], 2), "% and D2 = ", format_num(var_explained[2], 2), "%."),
  "",
  "## Main 2D Results",
  "",
  "| Superpop | MDS points | Samples | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Q99 | Max | Exact zero distances |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  results_rows,
  "",
  "The nearest-neighbor values are small because many subtelomeric points have an extremely close same-superpopulation point in the displayed D1-D2 MDS panel. The violin/box/jitter figure uses a pseudo-log distance axis to show both near-zero distances and sparse outliers without changing the raw TSV values.",
  "",
  "## 5D Sensitivity",
  "",
  "The slide figure should use the displayed 2D MDS metric. A secondary nearest-neighbor check in all five cached MDS dimensions is included in the TSV for auditability only.",
  "",
  "| Superpop | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Max |",
  "|---|---:|---:|---:|---:|---:|",
  sensitivity_rows,
  "",
  "## Sample Sizes",
  "",
  paste0("- Sequence-level MDS point counts: `", point_count_text, "`."),
  paste0("- Distinct sample counts: `", sample_count_text, "`."),
  "",
  "## Reproduction",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/make_nearest_same_superpop_mds.R",
  "```",
  "",
  "The script regenerates the TSVs, PNG/PDF figures, source manifest, README, validation note, and slide patch note in this directory.",
  "",
  "## Validation Note",
  "",
  "No centroid, RMS-radius, covariance-area, hull-area, all-pairwise summary, or in-group average metric is used as the slide metric. The code path for the slide distribution is the pointwise nearest-neighbor minimum in displayed D1-D2 MDS space."
)
write_lines(readme, file.path(out_dir, "README.md"))

validation <- c(
  "# Validation",
  "",
  paste0("Generated by `make_nearest_same_superpop_mds.R` on ", format(Sys.Date(), "%Y-%m-%d"), " UTC."),
  "",
  "## Source Checks",
  "",
  paste0("- MDS coordinate rows: `", nrow(all_points), "`."),
  paste0("- Coordinate dimensions available in the cached RDS: `", ncol(coords) - 1, "`."),
  paste0("- Superpopulation labels: no missing labels after joining `", label_export_tsv, "`."),
  paste0("- Current slide 08b PNG SHA-256: `", current_slide_png_sha, "`."),
  paste0("- Upstream source PNG SHA-256: `", source_png_sha, "`."),
  "- The script asserts the current slide PNG and upstream source PNG are byte-identical when both are present.",
  "",
  "## Metric Checks",
  "",
  "- Main slide metric uses Euclidean distance in displayed MDS dimensions `D1` and `D2`.",
  "- For each superpopulation, the script builds the within-superpopulation distance matrix, sets the diagonal to `Inf`, and retains one minimum distance per point.",
  "- Self matches are explicitly rejected in the generated table.",
  "- There is exactly one nearest-neighbor row per MDS point.",
  "- The optional 5D columns are labeled as sensitivity output and are not the recommended slide metric.",
  "",
  "## Explicit Non-Use",
  "",
  "- No centroid distance is computed for the slide metric.",
  "- No all-pairwise within-population distribution is reported as the slide metric.",
  "- No point is averaged against all same-superpopulation points.",
  "- Pairwise distances are only an internal exact search mechanism for the nearest-neighbor minimum."
)
write_lines(validation, file.path(out_dir, "VALIDATION.md"))

slide_patch <- c(
  "# Slide Patch Recommendation",
  "",
  "Task: `review-zoom-v7-slide08b-nearest-superpop-mds`.",
  "",
  "Do not integrate this directly in this task. This note is for `review-zoom-v7-fanin-render`.",
  "",
  "## Recommended Deck Insertion",
  "",
  "Split the current slide 08b content into two slides.",
  "",
  "### Slide 08b.1 - Source MDS View",
  "",
  "- Use asset: `slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/superpopulation_mds_original_style.png`.",
  "- Suggested title: `MDS colored by continental superpopulation`.",
  "- Suggested note: `Same D1-D2 MDS coordinates as the original slide 08b view; color is continental superpopulation and shape is p/q arm. Axes are rendered 1:1.`",
  "- Keep this as the familiar orientation/context slide before the metric slide.",
  "",
  "### Slide 08b.2 - Nearest Same-Superpopulation Distance",
  "",
  "- Use asset: `slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/nearest_same_superpop_distance_distribution.png`.",
  "- Suggested title: `Nearest same-superpopulation neighbor in MDS space`.",
  "- Suggested note: `For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Lower values mean each subtelomere has a closer same-superpopulation neighbor in sequence-similarity MDS space.`",
  "- Include the sample-size labels already shown on the y-axis.",
  "- If space is tight, keep only the title, the figure, and the lower-values interpretation note.",
  "",
  "## Required Caveat",
  "",
  "State or preserve this caveat in speaker notes: `Nearest same-superpopulation neighbor only; self is excluded. This is not a centroid metric, not an all-pairwise within-population distance plot, and not an average against all in-group points.`",
  "",
  "## Assets",
  "",
  "- `superpopulation_mds_original_style.png` / `.pdf`",
  "- `nearest_same_superpop_distance_distribution.png` / `.pdf`",
  "- `nearest_same_superpop_mds_distances.tsv`",
  "- `nearest_same_superpop_mds_summary.tsv`",
  "- `README.md`",
  "- `VALIDATION.md`"
)
write_lines(slide_patch, file.path(out_dir, "SLIDE_PATCH.md"))

message("Wrote nearest same-superpopulation MDS outputs to ", out_dir)
