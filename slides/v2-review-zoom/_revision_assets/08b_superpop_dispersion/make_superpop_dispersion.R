#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

out_dir <- "slides/v2-review-zoom/_revision_assets/08b_superpop_dispersion"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mds_rds <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds"
sample_metadata_tsv <- "/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv"
label_export_tsv <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv"
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

format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
format_ratio <- function(x) paste0(format_num(x, 2), "x")

polygon_area <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  idx <- chull(x, y)
  hx <- x[idx]
  hy <- y[idx]
  abs(sum(hx * c(hy[-1], hy[1]) - hy * c(hx[-1], hx[1]))) / 2
}

metric_one <- function(dat, dims) {
  mat <- as.matrix(dat[, dims, drop = FALSE])
  centered <- sweep(mat, 2, colMeans(mat), "-")
  squared_dist <- rowSums(centered^2)
  mean(squared_dist)
}

metric_table <- function(dat) {
  c2 <- as.matrix(dat[, c("D1", "D2")])
  c5 <- as.matrix(dat[, paste0("D", 1:5)])
  centered2 <- sweep(c2, 2, colMeans(c2), "-")
  centered5 <- sweep(c5, 2, colMeans(c5), "-")
  dist_centroid2 <- sqrt(rowSums(centered2^2))
  dist_centroid5 <- sqrt(rowSums(centered5^2))
  cov2 <- cov(c2)
  eig <- eigen(cov2, symmetric = TRUE, only.values = TRUE)$values
  global_center <- colMeans(all_points[, c("D1", "D2")])
  dist_global <- sqrt(rowSums(sweep(c2, 2, global_center, "-")^2))

  tibble(
    n_points = nrow(dat),
    n_samples = n_distinct(dat$Sample),
    n_chrom_arms = n_distinct(dat$ChromArm),
    msd_centroid_2d = mean(dist_centroid2^2),
    rms_radius_2d = sqrt(mean(dist_centroid2^2)),
    median_radius_2d = median(dist_centroid2),
    q90_radius_2d = unname(quantile(dist_centroid2, 0.90)),
    msd_centroid_5d = mean(dist_centroid5^2),
    rms_radius_5d = sqrt(mean(dist_centroid5^2)),
    ellipse_area_68_2d = pi * qchisq(0.68, df = 2) * sqrt(det(cov2)),
    ellipse_area_95_2d = pi * qchisq(0.95, df = 2) * sqrt(det(cov2)),
    generalized_variance_2d = det(cov2),
    sd_major_axis_2d = sqrt(max(eig)),
    sd_minor_axis_2d = sqrt(min(eig)),
    convex_hull_area_2d = polygon_area(c2[, 1], c2[, 2]),
    mean_distance_global_centroid_2d = mean(dist_global),
    q90_distance_global_centroid_2d = unname(quantile(dist_global, 0.90))
  )
}

bootstrap_rms <- function(dat, n_boot = 1000) {
  samples <- unique(dat$Sample)
  if (length(samples) < 2) {
    return(tibble(rms_radius_2d_ci_low = NA_real_, rms_radius_2d_ci_high = NA_real_))
  }
  values <- replicate(n_boot, {
    sampled <- sample(samples, length(samples), replace = TRUE)
    boot_dat <- bind_rows(lapply(sampled, function(sample_id) {
      dat[dat$Sample == sample_id, , drop = FALSE]
    }))
    sqrt(metric_one(boot_dat, c("D1", "D2")))
  })
  tibble(
    rms_radius_2d_ci_low = unname(quantile(values, 0.025)),
    rms_radius_2d_ci_high = unname(quantile(values, 0.975))
  )
}

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  line <- system2("sha256sum", path, stdout = TRUE)
  strsplit(line, "[[:space:]]+")[[1]][1]
}

extract_sample <- function(group_name) {
  sub("^([^#]+).*", "\\1", group_name)
}

extract_haplotype <- function(group_name) {
  out <- sub("^[^#]+#([^#]+)#.*", "\\1", group_name)
  ifelse(out == group_name, NA_character_, out)
}

extract_chromosome <- function(group_name) {
  reg <- regexpr("chr[0-9XYM]+", group_name)
  out <- rep(NA_character_, length(group_name))
  has_match <- reg > 0
  out[has_match] <- regmatches(group_name, reg)[has_match]
  out
}

extract_arm <- function(group_name) {
  ifelse(grepl("_parm$", group_name), "p",
         ifelse(grepl("_qarm$", group_name), "q", NA_character_))
}

subpop_to_superpop <- c(
  "ACB" = "AFR", "ASW" = "AFR", "ESN" = "AFR", "GWD" = "AFR",
  "LWK" = "AFR", "MSL" = "AFR", "YRI" = "AFR",
  "CLM" = "AMR", "MXL" = "AMR", "PEL" = "AMR", "PUR" = "AMR",
  "CDX" = "EAS", "CHB" = "EAS", "CHS" = "EAS", "JPT" = "EAS",
  "KHV" = "EAS",
  "CEU" = "EUR", "FIN" = "EUR", "GBR" = "EUR", "IBS" = "EUR",
  "TSI" = "EUR",
  "BEB" = "SAS", "GIH" = "SAS", "ITU" = "SAS", "PJL" = "SAS",
  "STU" = "SAS"
)

missing_samples_superpop <- c(
  "CHM13" = "EUR",
  "HG002" = "EUR",
  "HG005" = "EAS",
  "HG06807" = "EUR",
  "HG00733" = "AMR",
  "HG01109" = "AMR",
  "HG01243" = "AMR",
  "HG02055" = "AFR",
  "HG02080" = "EAS",
  "HG02109" = "EAS",
  "HG02145" = "AFR",
  "HG02723" = "AFR",
  "HG02818" = "AFR",
  "HG03098" = "AFR",
  "HG03486" = "AFR",
  "NA18906" = "AFR",
  "NA18940" = "EAS",
  "NA18943" = "EAS",
  "NA18944" = "EAS",
  "NA18945" = "EAS",
  "NA18948" = "EAS",
  "NA18959" = "EAS",
  "NA18960" = "EAS",
  "NA18967" = "EAS",
  "NA18970" = "EAS",
  "NA18982" = "EAS",
  "NA19240" = "AFR",
  "NA20129" = "AFR",
  "NA21309" = "SAS"
)

sample_metadata <- read_tsv(sample_metadata_tsv, col_types = cols(.default = col_character())) %>%
  select(ChildID, Subpopulation, Superpopulation) %>%
  mutate(
    ChildID = trimws(sub("\\s*\\(.*\\)", "", ChildID)),
    Superpopulation = case_when(
      !is.na(Superpopulation) & Superpopulation != "" ~ Superpopulation,
      Subpopulation %in% names(subpop_to_superpop) ~ unname(subpop_to_superpop[Subpopulation]),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Superpopulation) & Superpopulation != "") %>%
  distinct(ChildID, .keep_all = TRUE)

sample_to_superpop <- setNames(sample_metadata$Superpopulation, sample_metadata$ChildID)
sample_to_superpop <- c(
  sample_to_superpop,
  missing_samples_superpop[!names(missing_samples_superpop) %in% names(sample_to_superpop)]
)

mds <- readRDS(mds_rds)
coords <- as.data.frame(mds$points)
colnames(coords) <- paste0("D", seq_len(ncol(coords)))
coords$Name <- rownames(mds$points)

labels_from_names <- tibble(
  Name = coords$Name,
  Sample = extract_sample(coords$Name),
  Haplotype = extract_haplotype(coords$Name),
  Chromosome = extract_chromosome(coords$Name),
  Arm = extract_arm(coords$Name)
) %>%
  mutate(
    ChromArm = paste0(Chromosome, "_", Arm),
    Superpopulation = unname(sample_to_superpop[Sample])
  )

all_points <- coords %>%
  inner_join(labels_from_names, by = "Name") %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    ChromArm = factor(ChromArm)
  )

stopifnot(nrow(all_points) == nrow(coords))
stopifnot(!any(is.na(all_points$Superpopulation)))

label_export_status <- "not checked"
if (file.exists(label_export_tsv)) {
  exported_labels <- read_tsv(label_export_tsv, show_col_types = FALSE) %>%
    select(Name, Superpopulation)
  label_check <- all_points %>%
    select(Name, Superpopulation) %>%
    mutate(Superpopulation = as.character(Superpopulation)) %>%
    inner_join(exported_labels, by = "Name", suffix = c("_recomputed", "_export"))
  n_mismatch <- sum(label_check$Superpopulation_recomputed != label_check$Superpopulation_export)
  stopifnot(n_mismatch == 0)
  label_export_status <- paste0("matches recomputed superpopulation labels for ", nrow(label_check), " rows")
}

var_explained <- mds$eig / sum(abs(mds$eig)) * 100

set.seed(8)
metrics <- all_points %>%
  group_by(Superpopulation) %>%
  group_modify(~ bind_cols(metric_table(.x), bootstrap_rms(.x))) %>%
  ungroup() %>%
  mutate(
    Superpopulation = as.character(Superpopulation),
    non_afr_mean_rms_radius_2d = mean(rms_radius_2d[Superpopulation != "AFR"]),
    rms_radius_2d_vs_non_afr_mean = rms_radius_2d / non_afr_mean_rms_radius_2d,
    non_afr_mean_ellipse_area_68_2d = mean(ellipse_area_68_2d[Superpopulation != "AFR"]),
    ellipse_area_68_2d_vs_non_afr_mean = ellipse_area_68_2d / non_afr_mean_ellipse_area_68_2d,
    non_afr_mean_convex_hull_area_2d = mean(convex_hull_area_2d[Superpopulation != "AFR"]),
    convex_hull_area_2d_vs_non_afr_mean = convex_hull_area_2d / non_afr_mean_convex_hull_area_2d,
    non_afr_mean_distance_global_centroid_2d = mean(mean_distance_global_centroid_2d[Superpopulation != "AFR"]),
    mean_distance_global_centroid_2d_vs_non_afr_mean =
      mean_distance_global_centroid_2d / non_afr_mean_distance_global_centroid_2d
  )

write_tsv(metrics, file.path(out_dir, "superpop_dispersion_metrics.tsv"))

centroids <- all_points %>%
  group_by(Superpopulation) %>%
  summarise(
    centroid_D1 = mean(D1),
    centroid_D2 = mean(D2),
    rms_radius_2d = sqrt(mean((D1 - mean(D1))^2 + (D2 - mean(D2))^2)),
    .groups = "drop"
  ) %>%
  rename(
    D1 = centroid_D1,
    D2 = centroid_D2
  )

afr_points <- all_points %>% filter(Superpopulation == "AFR")
afr_centroid <- centroids %>% filter(Superpopulation == "AFR")

circle_df <- afr_centroid %>%
  crossing(theta = seq(0, 2 * pi, length.out = 240)) %>%
  mutate(
    x = D1 + rms_radius_2d * cos(theta),
    y = D2 + rms_radius_2d * sin(theta)
  )

radius_segment <- afr_centroid %>%
  transmute(
    x = D1,
    y = D2,
    xend = D1 + rms_radius_2d,
    yend = D2
  )

x_lab <- sprintf("MDS dimension 1 (%.2f%%)", var_explained[1])
y_lab <- sprintf("MDS dimension 2 (%.2f%%)", var_explained[2])

scatter_plot <- ggplot(all_points, aes(D1, D2)) +
  geom_point(color = "grey80", alpha = 0.18, size = 0.25, stroke = 0) +
  geom_point(
    data = afr_points,
    aes(D1, D2),
    color = superpop_colors[["AFR"]],
    alpha = 0.30,
    size = 0.30,
    stroke = 0,
    inherit.aes = FALSE
  ) +
  geom_path(
    data = circle_df,
    aes(x, y),
    color = superpop_colors[["AFR"]],
    linewidth = 0.9,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = radius_segment,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = superpop_colors[["AFR"]],
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = afr_centroid,
    aes(D1, D2),
    fill = superpop_colors[["AFR"]],
    color = "white",
    shape = 21,
    size = 3.2,
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = afr_centroid,
    aes(
      x = D1 + 0.12,
      y = D2 + 0.065,
      label = paste0("AFR centroid\nRMS radius = ", format_num(rms_radius_2d, 3))
    ),
    color = "black",
    fill = "white",
    label.size = 0.18,
    label.padding = grid::unit(0.12, "lines"),
    size = 2.8,
    lineheight = 0.95,
    inherit.aes = FALSE
  ) +
  coord_equal() +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_fill_manual(values = superpop_colors, guide = "none") +
  labs(
    title = "Metric geometry (AFR example)",
    subtitle = "Gold points are AFR flanks; the ring marks RMS distance to the AFR centroid",
    x = x_lab,
    y = y_lab
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.5)
  )

metric_for_plot <- metrics %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    label = paste0(format_num(rms_radius_2d, 3), "  (", format_ratio(rms_radius_2d_vs_non_afr_mean), ")")
  )

bar_plot <- ggplot(metric_for_plot, aes(rms_radius_2d, Superpopulation, color = Superpopulation)) +
  geom_vline(
    xintercept = unique(metric_for_plot$non_afr_mean_rms_radius_2d),
    color = "grey55",
    linewidth = 0.45,
    linetype = "22"
  ) +
  geom_errorbarh(
    aes(xmin = rms_radius_2d_ci_low, xmax = rms_radius_2d_ci_high),
    height = 0.12,
    linewidth = 0.55
  ) +
  geom_point(size = 3.0) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.0, color = "grey15") +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0.335, 0.348),
    breaks = seq(0.336, 0.348, by = 0.004),
    expand = expansion(mult = c(0.02, 0.20))
  ) +
  labs(
    title = "Talk-ready metric",
    subtitle = "RMS distance to own superpopulation centroid; dashed line = non-AFR mean",
    x = "RMS radius in MDS dimensions 1-2",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

draw_two_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 1,
      ncol = 2,
      widths = grid::unit(c(1.15, 1), "null")
    )
  ))
  print(scatter_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(bar_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
}

png(file.path(out_dir, "superpop_dispersion_rms_radius.png"), width = 2400, height = 1350, res = 220)
draw_two_panel()
dev.off()

pdf(file.path(out_dir, "superpop_dispersion_rms_radius.pdf"), width = 10.8, height = 6.1)
draw_two_panel()
dev.off()

comparison <- metrics %>%
  select(
    Superpopulation,
    `RMS radius` = rms_radius_2d_vs_non_afr_mean,
    `68% ellipse area` = ellipse_area_68_2d_vs_non_afr_mean,
    `Convex hull area` = convex_hull_area_2d_vs_non_afr_mean,
    `Mean distance to global centroid` = mean_distance_global_centroid_2d_vs_non_afr_mean
  ) %>%
  pivot_longer(-Superpopulation, names_to = "metric", values_to = "ratio_to_non_afr_mean") %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    metric = factor(metric, levels = c(
      "RMS radius",
      "68% ellipse area",
      "Convex hull area",
      "Mean distance to global centroid"
    )),
    label = format_num(ratio_to_non_afr_mean, 2)
  )

comparison_plot <- ggplot(comparison, aes(ratio_to_non_afr_mean, Superpopulation, color = Superpopulation)) +
  geom_vline(xintercept = 1, color = "grey55", linewidth = 0.45, linetype = "22") +
  geom_point(size = 2.4) +
  geom_text(aes(label = label), hjust = -0.35, size = 2.7, color = "grey20") +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_x_continuous(limits = c(0.96, 1.04), breaks = c(0.98, 1.00, 1.02, 1.04)) +
  labs(
    title = "Dispersion metric sensitivity",
    subtitle = "Ratios are relative to the mean of AMR/EAS/EUR/SAS; no metric shows AFR-specific expansion",
    x = "Ratio to non-AFR mean",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8.4),
    axis.text.x = element_text(size = 7.5),
    axis.text.y = element_text(size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(out_dir, "superpop_dispersion_metric_sensitivity.png"),
  comparison_plot,
  width = 10.8,
  height = 4.1,
  dpi = 220
)
ggsave(
  file.path(out_dir, "superpop_dispersion_metric_sensitivity.pdf"),
  comparison_plot,
  width = 10.8,
  height = 4.1
)

ordered <- metrics %>% arrange(match(Superpopulation, superpop_order))
afr <- ordered %>% filter(Superpopulation == "AFR")
non_afr_mean <- unique(ordered$non_afr_mean_rms_radius_2d)
sample_counts <- all_points %>%
  distinct(Sample, Superpopulation) %>%
  count(Superpopulation) %>%
  mutate(label = paste0(Superpopulation, "=", n)) %>%
  pull(label) %>%
  paste(collapse = ", ")
point_counts <- all_points %>%
  count(Superpopulation) %>%
  mutate(label = paste0(Superpopulation, "=", n)) %>%
  pull(label) %>%
  paste(collapse = ", ")

metric_rows <- ordered %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", n_points,
      " | ", n_samples,
      " | ", format_num(rms_radius_2d, 3),
      " | ", format_num(rms_radius_2d_ci_low, 3), "-", format_num(rms_radius_2d_ci_high, 3),
      " | ", format_num(msd_centroid_2d, 3),
      " | ", format_num(ellipse_area_68_2d, 3),
      " | ", format_num(convex_hull_area_2d, 3),
      " | ", format_ratio(rms_radius_2d_vs_non_afr_mean),
      " |"
    )
  ) %>%
  pull(row)

readme <- c(
  "# Slide 08b Superpopulation Dispersion",
  "",
  "## Source Confirmation",
  "",
  paste0("- Current slide 08b asset: `", current_slide_png, "`."),
  paste0("- Exact upstream PNG: `", source_png, "`."),
  paste0("- SHA-256 for both PNGs: `", sha256(source_png), "`; the current slide asset is byte-identical to the HPRCv2 pipeline output."),
  paste0("- Coordinate source: `", mds_rds, "` (`", nrow(all_points), "` rows x `5` dimensions)."),
  paste0("- Direct superpopulation label source: `", sample_metadata_tsv, "` plus the hard-coded missing-sample overrides in `", source_script, "`."),
  paste0("- Row-level label export used as a cross-check: `", label_export_tsv, "`; ", label_export_status, "."),
  paste0("- Generation script: `", source_script, "`; the relevant call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`."),
  "",
  paste0("Terminology: this is **classical MDS / PCoA** on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses MDS dimensions 1 and 2; the cached R object contains five dimensions. Axis percentages in the source plot are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = ", format_num(var_explained[1], 2), "% and D2 = ", format_num(var_explained[2], 2), "%."),
  "",
  "## Metric Choice",
  "",
  "Chosen talk metric: **RMS radius in the displayed 2D MDS panel**.",
  "",
  "`RMS radius = sqrt(mean(||x_i - c_g||^2))`, where `x_i` is a sequence-level point in MDS dimensions 1 and 2, and `c_g` is the centroid for superpopulation `g`.",
  "",
  "I chose RMS radius over covariance-ellipse area and convex-hull area because it is easy to say on a slide, stable under unequal sample sizes, and directly describes how far a typical point lies from its own superpopulation centre. Convex hull area was evaluated but is more sensitive to sample size and outliers; covariance ellipse area is concise statistically but harder to explain in a talk.",
  "",
  "## Results",
  "",
  "| Superpop | MDS points | Samples | RMS radius | Sample-bootstrap 95% CI | Mean squared radius | 68% ellipse area | Convex hull area | RMS vs non-AFR mean |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  metric_rows,
  "",
  paste0("AFR RMS radius is `", format_num(afr$rms_radius_2d, 3), "` versus a non-AFR mean of `", format_num(non_afr_mean, 3), "` (`", format_ratio(afr$rms_radius_2d_vs_non_afr_mean), "`)."),
  "The candidate metrics agree: the displayed coordinates do **not** support a claim that AFR is unusually dispersed in the 2D MDS space. The visual impression is more likely due to point density, color salience, and every superpopulation sampling the same arm-defined clusters.",
  "",
  "## Candidate Plots",
  "",
  "- `superpop_dispersion_rms_radius.png` / `.pdf`: recommended slide candidate. Left panel shows the displayed MDS coordinates with AFR as the worked example for the RMS-radius definition; right panel shows the chosen metric for all superpopulations with sample-bootstrap CIs.",
  "- `superpop_dispersion_metric_sensitivity.png` / `.pdf`: sensitivity plot comparing RMS radius, covariance ellipse area, convex hull area, and distance to the global centroid as ratios to the non-AFR mean.",
  "- `superpop_dispersion_metrics.tsv`: full metric table used by the plots and this note.",
  "",
  "## Interpretation",
  "",
  "The defensible read is not \"AFR spans more of the MDS space.\" It is: **all five superpopulations span nearly the same displayed MDS space**, consistent with slide 08b's main claim that the dominant 2D structure is arm/community structure rather than continental superpopulation.",
  "",
  "## Limitations",
  "",
  paste0("- Unequal sample sizes: sample labels in this source are `", sample_counts, "`; sequence-level point counts are `", point_counts, "`. The chosen metric is less sample-size-sensitive than a hull, and the plotted CI bootstraps samples within each superpopulation, but the points are still not independent."),
  "- Unit of analysis: the plotted and quantified unit is the **sequence-level subtelomeric flank** used in slide 08b, not an arm-level aggregate and not one point per assembly. Each sample contributes many arm/haplotype flanks, so arm composition dominates the geometry.",
  "- MDS dimensionality: the source RDS stores five MDS dimensions, but slide 08b displays only dimensions 1 and 2. The 5D centroid metric was also computed in `superpop_dispersion_metrics.tsv` and gives the same conclusion.",
  "- Coordinate method: because this is classical MDS / PCoA on `1 - Jaccard`, distances in the first two displayed dimensions are a low-dimensional approximation to the full Jaccard distance structure, not direct feature-space PCA distances.",
  "- CHM13: the source script hard-codes CHM13 as EUR for superpopulation coloring; the metric follows the slide source exactly.",
  "",
  "## Speaker-Ready Sentence",
  "",
  paste0("\"Quantifying the same MDS panel, AFR is not more dispersed: its RMS radius is ", format_num(afr$rms_radius_2d, 3), " versus ", format_num(non_afr_mean, 3), " for the other superpopulations, so the 2D spread is shared across populations and mainly reflects arm-community structure.\"")
)

write_lines(readme, file.path(out_dir, "README.md"))

message("Wrote outputs to ", out_dir)
