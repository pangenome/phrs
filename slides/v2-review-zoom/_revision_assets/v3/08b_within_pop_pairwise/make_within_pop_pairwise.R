#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

out_dir <- "slides/v2-review-zoom/_revision_assets/v3/08b_within_pop_pairwise"
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

set.seed(80)

format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
format_ci <- function(lo, hi, digits = 3) paste0(format_num(lo, digits), "-", format_num(hi, digits))

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

bootstrap_sample_pair_mean <- function(dat, dims, n_boot = 500, pair_draws = 20000) {
  sample_ids <- unique(dat$Sample)
  if (length(sample_ids) < 2) {
    return(tibble(
      mean_ci_low = NA_real_,
      mean_ci_high = NA_real_,
      bootstrap_reps = 0L,
      bootstrap_pair_draws_per_rep = 0L
    ))
  }

  mat <- as.matrix(dat[, dims, drop = FALSE])
  idx_by_sample <- split(seq_len(nrow(dat)), dat$Sample)
  boot_values <- replicate(n_boot, {
    boot_samples <- sample(sample_ids, length(sample_ids), replace = TRUE)
    boot_idx <- unlist(idx_by_sample[boot_samples], use.names = FALSE)
    if (length(boot_idx) < 2) return(NA_real_)

    pos_i <- sample.int(length(boot_idx), pair_draws, replace = TRUE)
    pos_j <- sample.int(length(boot_idx), pair_draws, replace = TRUE)
    same_position <- pos_i == pos_j
    while (any(same_position)) {
      pos_j[same_position] <- sample.int(length(boot_idx), sum(same_position), replace = TRUE)
      same_position <- pos_i == pos_j
    }

    diff <- mat[boot_idx[pos_i], , drop = FALSE] - mat[boot_idx[pos_j], , drop = FALSE]
    mean(sqrt(rowSums(diff * diff)))
  })

  tibble(
    mean_ci_low = unname(quantile(boot_values, 0.025, na.rm = TRUE)),
    mean_ci_high = unname(quantile(boot_values, 0.975, na.rm = TRUE)),
    bootstrap_reps = n_boot,
    bootstrap_pair_draws_per_rep = pair_draws
  )
}

summarize_pairwise_distances <- function(dat, dims, dimensionality) {
  mat <- as.matrix(dat[, dims, drop = FALSE])
  distances <- as.numeric(dist(mat, method = "euclidean", diag = FALSE, upper = FALSE))
  qs <- quantile(distances, c(0.25, 0.50, 0.75), names = FALSE)
  ci <- bootstrap_sample_pair_mean(dat, dims)

  summary <- tibble(
    dimensionality = dimensionality,
    n_points = nrow(dat),
    n_samples = n_distinct(dat$Sample),
    n_chrom_arms = n_distinct(dat$ChromArm),
    n_pairs = length(distances),
    mean_pairwise_distance = mean(distances),
    median_pairwise_distance = qs[2],
    q25_pairwise_distance = qs[1],
    q75_pairwise_distance = qs[3],
    iqr_pairwise_distance = qs[3] - qs[1],
    min_pairwise_distance = min(distances),
    max_pairwise_distance = max(distances),
    sd_pairwise_distance = sd(distances)
  ) %>%
    bind_cols(ci)

  sampled <- tibble(
    dimensionality = dimensionality,
    pairwise_distance = sample(distances, min(length(distances), 12000), replace = FALSE)
  )

  list(summary = summary, sampled = sampled)
}

analysis_by_superpop <- lapply(superpop_order, function(superpop) {
  dat <- all_points %>% filter(Superpopulation == superpop)

  two_d <- summarize_pairwise_distances(dat, c("D1", "D2"), "2D displayed MDS")
  five_d <- summarize_pairwise_distances(dat, paste0("D", 1:5), "5D cached MDS")

  list(
    summary = bind_rows(two_d$summary, five_d$summary) %>%
      mutate(Superpopulation = superpop, .before = 1),
    sampled = bind_rows(two_d$sampled, five_d$sampled) %>%
      mutate(Superpopulation = superpop, .before = 1)
  )
})

pairwise_summary <- bind_rows(lapply(analysis_by_superpop, `[[`, "summary")) %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    dimensionality = factor(dimensionality, levels = c("2D displayed MDS", "5D cached MDS"))
  ) %>%
  arrange(dimensionality, Superpopulation) %>%
  mutate(
    Superpopulation = as.character(Superpopulation),
    dimensionality = as.character(dimensionality),
    mean_ci_method = "sample-bootstrap over samples; each replicate estimates the pair mean from 20,000 drawn within-population pairs"
  )

pairwise_distance_sample <- bind_rows(lapply(analysis_by_superpop, `[[`, "sampled")) %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    dimensionality = factor(dimensionality, levels = c("2D displayed MDS", "5D cached MDS"))
  ) %>%
  arrange(dimensionality, Superpopulation) %>%
  mutate(
    Superpopulation = as.character(Superpopulation),
    dimensionality = as.character(dimensionality)
  )

write_tsv(pairwise_summary, file.path(out_dir, "within_pop_pairwise_summary.tsv"))
write_tsv(pairwise_distance_sample, file.path(out_dir, "within_pop_pairwise_distance_sample.tsv"))

summary_2d <- pairwise_summary %>%
  filter(dimensionality == "2D displayed MDS") %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    mean_label = paste0(
      "mean ", format_num(mean_pairwise_distance, 3)
    )
  )

sample_2d <- pairwise_distance_sample %>%
  filter(dimensionality == "2D displayed MDS") %>%
  mutate(Superpopulation = factor(Superpopulation, levels = superpop_order))

x_upper_2d <- max(sample_2d$pairwise_distance, summary_2d$mean_pairwise_distance) + 0.18
summary_2d <- summary_2d %>%
  mutate(x_label = x_upper_2d - 0.015)

main_plot <- ggplot(sample_2d, aes(x = pairwise_distance, y = Superpopulation)) +
  geom_violin(
    aes(fill = Superpopulation),
    alpha = 0.18,
    width = 0.82,
    linewidth = 0.25,
    color = "grey35",
    trim = FALSE
  ) +
  geom_boxplot(
    aes(fill = Superpopulation),
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.55,
    linewidth = 0.35,
    color = "grey20"
  ) +
  geom_errorbarh(
    data = summary_2d,
    aes(
      xmin = mean_ci_low,
      xmax = mean_ci_high,
      y = Superpopulation,
      color = Superpopulation
    ),
    inherit.aes = FALSE,
    height = 0.06,
    linewidth = 0.9
  ) +
  geom_point(
    data = summary_2d,
    aes(x = mean_pairwise_distance, y = Superpopulation, color = Superpopulation),
    inherit.aes = FALSE,
    shape = 18,
    size = 3.6
  ) +
  geom_label(
    data = summary_2d,
    aes(
      x = x_label,
      y = Superpopulation,
      label = mean_label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3.15,
    color = "grey15",
    fill = "white",
    alpha = 0.92,
    label.size = 0.15,
    label.padding = grid::unit(0.12, "lines")
  ) +
  scale_fill_manual(values = superpop_colors, guide = "none") +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(0, x_upper_2d)) +
  labs(
    title = "Within-population pairwise variation in the MDS / PCoA panel",
    subtitle = sprintf(
      "Euclidean distances among unordered same-superpopulation PHR-flank pairs in displayed dimensions D1-D2 (%.2f%% + %.2f%%)",
      var_explained[1],
      var_explained[2]
    ),
    x = "Pairwise distance in displayed MDS dimensions 1-2",
    y = NULL,
    caption = "Violin/box use sampled pair distances for display; diamond is the exact mean; line is a sample-bootstrap 95% CI. Sequence-level flanks and pair distances are not independent."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.5, color = "grey25"),
    plot.caption = element_text(size = 8.2, color = "grey35", hjust = 0),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 18, 10, 12)
  )

ggsave(
  file.path(out_dir, "within_pop_pairwise_2d_distribution.png"),
  main_plot,
  width = 12.2,
  height = 6.6,
  dpi = 240
)
ggsave(
  file.path(out_dir, "within_pop_pairwise_2d_distribution.pdf"),
  main_plot,
  width = 12.2,
  height = 6.6
)

summary_long <- pairwise_summary %>%
  mutate(
    Superpopulation = factor(Superpopulation, levels = superpop_order),
    dimensionality = factor(dimensionality, levels = c("2D displayed MDS", "5D cached MDS")),
    mean_label = format_num(mean_pairwise_distance, 3),
    median_iqr_label = paste0(
      "median ", format_num(median_pairwise_distance, 3),
      "  IQR ", format_num(iqr_pairwise_distance, 3)
    )
  )

sensitivity_plot <- ggplot(summary_long, aes(x = mean_pairwise_distance, y = Superpopulation)) +
  geom_errorbarh(
    aes(xmin = mean_ci_low, xmax = mean_ci_high, color = Superpopulation),
    height = 0.10,
    linewidth = 0.8
  ) +
  geom_point(aes(color = Superpopulation), shape = 18, size = 3.3) +
  geom_text(
    aes(label = mean_label),
    nudge_x = 0.018,
    hjust = 0,
    size = 3.0,
    color = "grey15"
  ) +
  facet_wrap(~ dimensionality, nrow = 1, scales = "free_x") +
  scale_color_manual(values = superpop_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.20))) +
  labs(
    title = "Pairwise metric sensitivity to MDS dimensionality",
    subtitle = "Same within-superpopulation mean pair distance computed in the displayed 2D panel and in all five cached MDS dimensions",
    x = "Mean within-population pairwise distance",
    y = NULL,
    caption = "Intervals are sample-bootstrap 95% CIs estimated from drawn same-superpopulation pairs; exact medians and IQRs are in the TSV."
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5, color = "grey25"),
    plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.2),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(out_dir, "within_pop_pairwise_2d_vs_5d_sensitivity.png"),
  sensitivity_plot,
  width = 11.8,
  height = 5.0,
  dpi = 240
)
ggsave(
  file.path(out_dir, "within_pop_pairwise_2d_vs_5d_sensitivity.pdf"),
  sensitivity_plot,
  width = 11.8,
  height = 5.0
)

ordered_2d <- summary_2d %>% arrange(match(Superpopulation, superpop_order))
ordered_5d <- pairwise_summary %>%
  filter(dimensionality == "5D cached MDS") %>%
  arrange(match(Superpopulation, superpop_order))

metric_rows <- ordered_2d %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", comma(n_points),
      " | ", comma(n_samples),
      " | ", comma(n_pairs),
      " | ", format_num(mean_pairwise_distance, 3),
      " | ", format_ci(mean_ci_low, mean_ci_high, 3),
      " | ", format_num(median_pairwise_distance, 3),
      " | ", format_num(q25_pairwise_distance, 3), "-", format_num(q75_pairwise_distance, 3),
      " | ", format_num(iqr_pairwise_distance, 3),
      " |"
    )
  ) %>%
  pull(row)

sensitivity_rows <- ordered_5d %>%
  transmute(
    row = paste0(
      "| ", Superpopulation,
      " | ", format_num(mean_pairwise_distance, 3),
      " | ", format_ci(mean_ci_low, mean_ci_high, 3),
      " | ", format_num(median_pairwise_distance, 3),
      " | ", format_num(q25_pairwise_distance, 3), "-", format_num(q75_pairwise_distance, 3),
      " | ", format_num(iqr_pairwise_distance, 3),
      " |"
    )
  ) %>%
  pull(row)

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

source_png_sha <- sha256(source_png)
current_slide_png_sha <- sha256(current_slide_png)
stopifnot(identical(source_png_sha, current_slide_png_sha))

readme <- c(
  "# Slide 08b Within-Population Pairwise Variation",
  "",
  "## Metric",
  "",
  "The metric is the average Euclidean distance between every unordered pair of sequence-level PHR-flank points from the same superpopulation in the displayed MDS dimensions 1 and 2.",
  "",
  "This replaces the v2 dispersion framing: the candidate recommendation is based on same-superpopulation point-to-point distances only.",
  "",
  "## Source Confirmation",
  "",
  paste0("- Current slide 08b asset: `", current_slide_png, "`."),
  paste0("- Exact upstream PNG: `", source_png, "`."),
  paste0("- Upstream PNG SHA-256: `", source_png_sha, "`."),
  paste0("- Current slide PNG SHA-256: `", current_slide_png_sha, "`; the current slide asset is byte-identical to the HPRCv2 pipeline output."),
  paste0("- Coordinate source: `", mds_rds, "` (`", nrow(all_points), "` rows x `5` dimensions)."),
  paste0("- Superpopulation label source: `", sample_metadata_tsv, "` plus hard-coded missing-sample overrides in `", source_script, "`."),
  paste0("- Row-level label export cross-check: `", label_export_tsv, "`; ", label_export_status, "."),
  paste0("- Generation script: `", source_script, "`; the relevant call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`."),
  "",
  paste0("Terminology: this is **classical MDS / PCoA** on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses dimensions 1 and 2; the cached R object contains five dimensions. Axis percentages in the source plot are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = ", format_num(var_explained[1], 2), "% and D2 = ", format_num(var_explained[2], 2), "%."),
  "",
  "## Main 2D Results",
  "",
  "| Superpop | MDS points | Samples | Same-pop pairs | Mean distance | Sample-bootstrap 95% CI | Median distance | Q1-Q3 | IQR |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  metric_rows,
  "",
  "The 2D point-to-point summaries are close across superpopulations, so the displayed panel does not support a claim that one superpopulation has uniquely larger within-population spread.",
  "",
  "## 5D Sensitivity",
  "",
  "| Superpop | Mean distance | Sample-bootstrap 95% CI | Median distance | Q1-Q3 | IQR |",
  "|---|---:|---:|---:|---:|---:|",
  sensitivity_rows,
  "",
  "The five-dimensional check uses all cached MDS dimensions and preserves the same qualitative read as the displayed two-dimensional panel.",
  "",
  "## Outputs",
  "",
  "- `within_pop_pairwise_2d_distribution.png` / `.pdf`: main candidate plot for slide 08b; violin/box shows sampled same-superpopulation pair distances, diamond shows the exact mean, and the line shows the sample-bootstrap CI.",
  "- `within_pop_pairwise_2d_vs_5d_sensitivity.png` / `.pdf`: sensitivity check comparing the same metric in displayed 2D coordinates and cached 5D coordinates.",
  "- `within_pop_pairwise_summary.tsv`: exact per-superpopulation summaries for 2D and 5D distances, including n points, n samples, number of pairs, mean, median, quartiles, IQR, and CI columns.",
  "- `within_pop_pairwise_distance_sample.tsv`: sampled pair distances used for the violin/box display; this avoids writing tens of millions of pair rows while preserving an auditable plotted distribution sample.",
  "- `make_within_pop_pairwise.R`: reproducible generator for the outputs in this directory.",
  "",
  "## Limitations",
  "",
  paste0("- Unequal sample sizes: sample labels in this source are `", sample_counts, "`; sequence-level point counts are `", point_counts, "`. The number of same-population pairs scales with point count, so the summary table reports both n points and n pairs."),
  "- Non-independent PHR flanks: the plotted and quantified unit is a sequence-level subtelomeric flank, not one independent individual. Each sample contributes many arm/haplotype flanks, and pair distances share points.",
  "- Pairwise distances are also non-independent because each point appears in many pairs; the CI is a descriptive sample-bootstrap interval and should not be read as an independent-pair hypothesis test.",
  "- MDS dimensionality: the main candidate uses only the displayed D1-D2 panel. The 5D sensitivity check is included because the cached MDS object stores five coordinates, but both are low-dimensional summaries of the full Jaccard distance structure.",
  "- Coordinate method: because this is classical MDS / PCoA on `1 - Jaccard`, the reported values are distances in MDS coordinate space, not direct feature-space PCA distances.",
  "- CHM13: the source script hard-codes CHM13 as EUR for superpopulation coloring; this analysis follows the slide source exactly.",
  "",
  "## Speaker-Ready Sentence",
  "",
  "\"Using same-superpopulation point-to-point distances in the displayed MDS panel, the populations have similar within-population spread, so the visual structure is better read as shared arm/community geometry than as one population spreading more than the rest.\""
)

write_lines(readme, file.path(out_dir, "README.md"))

message("Wrote outputs to ", out_dir)
