#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ragg)
})

asset_dir <- getwd()
plot_dir <- file.path(asset_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input TSV: ", path)
  read.delim(path, check.names = FALSE)
}

read_matrix_tsv <- function(path) {
  x <- read_tsv(path)
  rn <- x[[1]]
  x <- x[-1]
  rownames(x) <- rn
  as.matrix(x)
}

save_plot_pair <- function(plot, stem, width = 1920, height = 1080, res = 144) {
  png_path <- file.path(plot_dir, paste0(stem, ".png"))
  pdf_path <- file.path(plot_dir, paste0(stem, ".pdf"))
  agg_png(png_path, width = width, height = height, units = "px", res = res, background = "white")
  print(plot)
  dev.off()
  ggsave(pdf_path, plot, width = width / res, height = height / res, units = "in", device = cairo_pdf)
}

seq_dist_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
community_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"

gm_summary_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_summary.tsv"
gm_per_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_cell.tsv"
gm_per_community_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_community_per_cell.tsv"
gm_cf_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_per_cell.tsv"
gm_cf_arm_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_arm.tsv"
gm_mantel_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_3d.tsv"
gm_3d_matrix_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_arm_3d_distance_matrix.tsv"
gm_3dg_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/gm12878_01.impute3.round4.clean.3dg.gz"

sperm_summary_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_summary.tsv"
sperm_per_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_cell.tsv"
sperm_per_community_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv"
sperm_cf_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_per_cell.tsv"
sperm_cf_arm_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_arm.tsv"
sperm_mantel_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_3d.tsv"
sperm_3d_matrix_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_arm_3d_distance_matrix.tsv"

format_p <- function(p) {
  if (p < 1e-3) {
    sprintf("%.1e", p)
  } else {
    sprintf("%.4f", p)
  }
}

make_mantel_proximity_plot <- function(dataset_label, tech_label, n_cells, d3_matrix_path,
                                       mantel_path, stem, within_color) {
  seq_dist <- read_matrix_tsv(seq_dist_path)
  d3_dist <- read_matrix_tsv(d3_matrix_path)
  comm <- read_tsv(community_path) %>%
    select(ChromArm, Community)
  comm_map <- setNames(comm$Community, comm$ChromArm)
  mantel <- read_tsv(mantel_path) %>% filter(test == "mantel_3d") %>% slice(1)

  common <- sort(intersect(rownames(d3_dist), rownames(seq_dist)))
  pairs <- t(combn(common, 2))
  df <- tibble(
    arm_a = pairs[, 1],
    arm_b = pairs[, 2],
    sequence_similarity = 1 - seq_dist[cbind(pairs[, 1], pairs[, 2])],
    proximity = -d3_dist[cbind(pairs[, 1], pairs[, 2])],
    pair_class = ifelse(
      !is.na(comm_map[pairs[, 1]]) & comm_map[pairs[, 1]] == comm_map[pairs[, 2]],
      "within",
      "between"
    )
  ) %>%
    filter(is.finite(sequence_similarity), is.finite(proximity)) %>%
    mutate(
      pair_class = factor(pair_class, levels = c("between", "within"))
    )

  label <- paste0(
    "Mantel rho=", sprintf("%.3f", mantel$rho), "\n",
    "p=", format_p(mantel$p_value), "; n=", nrow(df), " arm pairs\n",
    "Y is proximity = -mean 3D distance"
  )

  detail <- paste0(
    dataset_label, " | ", tech_label, " | ", n_cells, " cells | ",
    "source: ", basename(d3_matrix_path), " + ", basename(seq_dist_path)
  )

  p <- ggplot(df, aes(x = sequence_similarity, y = proximity)) +
    geom_point(
      data = df %>% filter(pair_class == "between"),
      color = "#c7ccd4",
      alpha = 0.40,
      size = 2.5
    ) +
    geom_point(
      data = df %>% filter(pair_class == "within"),
      color = within_color,
      fill = within_color,
      alpha = 0.88,
      size = 4.3,
      shape = 21,
      stroke = 0.45
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "#555555", linewidth = 0.9, linetype = "dashed") +
    annotate(
      "label",
      x = min(df$sequence_similarity) + 0.02,
      y = max(df$proximity) - 0.04,
      label = label,
      hjust = 0,
      vjust = 1,
      label.size = 0.35,
      fill = "#fff8db",
      color = "#111111",
      size = 5.1,
      lineheight = 0.96
    ) +
    annotate(
      "text",
      x = max(df$sequence_similarity),
      y = min(df$proximity) + 0.025,
      label = detail,
      hjust = 1,
      vjust = 0,
      color = "#696969",
      size = 3.35
    ) +
    scale_x_continuous(limits = c(0, 0.85), expand = expansion(mult = c(0.02, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.06))) +
    labs(
      title = paste0(dataset_label, ": sequence similarity vs 3D proximity"),
      subtitle = "Higher X = more shared subtelomeric sequence; higher Y = closer in 3D",
      x = "Sequence similarity (1 - Jaccard distance, arm-level)",
      y = "3D proximity\n(-mean distance; higher = closer)"
    ) +
    theme_minimal(base_family = "DejaVu Sans", base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 22, color = "#1a3a6b"),
      plot.subtitle = element_text(size = 14, color = "#444444"),
      axis.title.x = element_text(face = "bold", size = 16, color = "#222222"),
      axis.title.y = element_text(face = "bold", size = 12.5, color = "#222222", margin = margin(r = 12)),
      axis.text = element_text(size = 12, color = "#333333"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(24, 34, 30, 76)
    )

  save_plot_pair(p, stem, width = 1920, height = 1080, res = 144)

  manifest <- df %>%
    summarise(
      dataset = dataset_label,
      technology = tech_label,
      n_cells = n_cells,
      n_shared_arms = length(common),
      n_pairs = n(),
      n_within_pairs = sum(pair_class == "within"),
      n_between_pairs = sum(pair_class == "between"),
      mantel_rho = mantel$rho,
      mantel_p = mantel$p_value,
      source_3d_matrix = d3_matrix_path,
      source_sequence_matrix = seq_dist_path,
      source_community_map = community_path
    )
  write.table(
    manifest,
    file.path(plot_dir, paste0(stem, "_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

make_mantel_proximity_plot(
  "GM12878",
  "Dip-C",
  16,
  gm_3d_matrix_path,
  gm_mantel_path,
  "gm12878_mantel_proximity",
  "#c7443e"
)

make_mantel_proximity_plot(
  "Sperm",
  "single-cell Hi-C / 3DG",
  20,
  sperm_3d_matrix_path,
  sperm_mantel_path,
  "sperm_all20_mantel_proximity",
  "#b06920"
)

gm_summary <- read_tsv(gm_summary_path)
sperm_summary <- read_tsv(sperm_summary_path)
gm_cell <- read_tsv(gm_per_cell_path)
sperm_cell <- read_tsv(sperm_per_cell_path)
gm_comm <- read_tsv(gm_per_community_path)
sperm_comm <- read_tsv(sperm_per_community_path)
gm_cf <- read_tsv(gm_cf_cell_path)
sperm_cf <- read_tsv(sperm_cf_cell_path)
gm_cf_arm <- read_tsv(gm_cf_arm_path)
sperm_cf_arm <- read_tsv(sperm_cf_arm_path)

gm_fisher <- gm_summary %>% filter(test == "fisher_combined") %>% slice(1)
sperm_fisher <- sperm_summary %>% filter(test == "fisher_combined") %>% slice(1)

gm_s_all <- gm_comm %>% filter(community == "S_all")
sperm_s_all <- sperm_comm %>% filter(community == "S_all")

wb_plot_df <- bind_rows(
  gm_cell %>% transmute(cell_type = "GM12878 Dip-C", signal = "C-community sharing", cell_id = as.character(cell_id), ratio),
  gm_s_all %>% transmute(cell_type = "GM12878 Dip-C", signal = "S_all non-sharing", cell_id = as.character(cell_id), ratio),
  sperm_cell %>% transmute(cell_type = "Sperm scHi-C", signal = "C-community sharing", cell_id = as.character(cell_id), ratio),
  sperm_s_all %>% transmute(cell_type = "Sperm scHi-C", signal = "S_all non-sharing", cell_id = as.character(cell_id), ratio)
) %>%
  mutate(
    cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")),
    signal = factor(signal, levels = c("C-community sharing", "S_all non-sharing"))
  )

wb_summary <- tibble(
  cell_type = c("GM12878 Dip-C", "GM12878 Dip-C", "Sperm scHi-C", "Sperm scHi-C"),
  signal = c("C-community sharing", "S_all non-sharing", "C-community sharing", "S_all non-sharing"),
  display_ratio = c(gm_fisher$mean_ratio, 1.106, sperm_fisher$mean_ratio, 1.397),
  n_cells = c(nrow(gm_cell), nrow(gm_s_all), nrow(sperm_cell), nrow(sperm_s_all)),
  n_below_1 = c(sum(gm_cell$ratio < 1), sum(gm_s_all$ratio < 1), sum(sperm_cell$ratio < 1), sum(sperm_s_all$ratio < 1)),
  effect_label = c("6.9% closer", "11% farther", "60% closer", "40% farther"),
  p_label = c("Fisher p=2.4e-05", "0/16 cells < 1", "Fisher p=3.9e-51", "1/20 cells < 1")
) %>%
  mutate(
    cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")),
    signal = factor(signal, levels = c("C-community sharing", "S_all non-sharing")),
    y = ifelse(signal == "C-community sharing", 1.13, 1.74),
    label = paste0("W/B=", sprintf("%.3f", display_ratio), "\n", effect_label, "\n", p_label)
  )

wb_plot <- ggplot(wb_plot_df, aes(x = signal, y = ratio, fill = signal)) +
  geom_hline(yintercept = 1, linewidth = 0.9, linetype = "dashed", color = "#4d4d4d") +
  geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.78, color = "#2b2b2b", linewidth = 0.75) +
  geom_point(
    aes(color = signal),
    position = position_jitter(width = 0.09, height = 0, seed = 19),
    size = 3.6,
    alpha = 0.82,
    stroke = 0.15
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5.3, fill = "#ffffff", color = "#111111", stroke = 0.95) +
  geom_text(
    data = wb_summary,
    aes(x = signal, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold",
    lineheight = 0.94,
    color = "#1f1f1f"
  ) +
  facet_wrap(~ cell_type, nrow = 1) +
  scale_x_discrete(labels = c("C-community sharing" = "C-community\nsharing", "S_all non-sharing" = "S_all\nnon-sharing")) +
  scale_y_continuous(limits = c(0.05, 2.05), breaks = seq(0.25, 2.0, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c("C-community sharing" = "#2d6e9f", "S_all non-sharing" = "#d68627")) +
  scale_color_manual(values = c("C-community sharing" = "#1d4f75", "S_all non-sharing" = "#99510d")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "W/B distance ratio"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 20, color = "#1a3a6b"),
    axis.text.x = element_text(face = "bold", size = 16, color = "#1f1f1f", lineheight = 0.92),
    axis.text.y = element_text(size = 13, color = "#333333"),
    axis.title.y = element_text(size = 17, face = "bold", color = "#222222", margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(26, 40, 24, 42)
  )

save_plot_pair(wb_plot, "wb_negative_control_reduced_text")

arm_cor <- function(df) {
  ct <- suppressWarnings(cor.test(df$mean_jaccard, -df$mean_3d_distance, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), p_value = ct$p.value, n_pairs = nrow(df))
}

cf_plot_df <- bind_rows(
  gm_cf %>% transmute(cell_type = "GM12878 Dip-C", cell_id = as.character(cell_id), spearman_rho),
  sperm_cf %>% transmute(cell_type = "Sperm scHi-C", cell_id = as.character(cell_id), spearman_rho)
) %>%
  mutate(cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")))

cf_summary <- bind_rows(
  tibble(
    cell_type = "GM12878 Dip-C",
    median_rho = median(gm_cf$spearman_rho),
    n_positive = sum(gm_cf$spearman_rho > 0),
    n_cells = nrow(gm_cf)
  ) %>% bind_cols(arm_cor(gm_cf_arm) %>% rename(arm_level_rho = rho, arm_level_p = p_value, arm_level_n_pairs = n_pairs)),
  tibble(
    cell_type = "Sperm scHi-C",
    median_rho = median(sperm_cf$spearman_rho),
    n_positive = sum(sperm_cf$spearman_rho > 0),
    n_cells = nrow(sperm_cf)
  ) %>% bind_cols(arm_cor(sperm_cf_arm) %>% rename(arm_level_rho = rho, arm_level_p = p_value, arm_level_n_pairs = n_pairs))
) %>%
  mutate(cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")))

cf_ann <- cf_summary %>%
  mutate(
    y = c(0.238, 0.238),
    label = paste0("median rho=", sprintf("%.3f", median_rho), "\n", n_positive, "/", n_cells, " cells positive")
  )

arm_ann <- cf_summary %>%
  mutate(
    y = ifelse(cell_type == "GM12878 Dip-C", arm_level_rho + 0.018, arm_level_rho - 0.024),
    label = ifelse(
      cell_type == "GM12878 Dip-C",
      paste0("arm-level rho=", sprintf("%.3f", arm_level_rho), "\np=1.1e-18"),
      paste0("pooled arm-level rho=", sprintf("%.3f", arm_level_rho), "\np=0.197 caveat")
    )
  )

cf_plot <- ggplot(cf_plot_df, aes(x = cell_type, y = spearman_rho, fill = cell_type)) +
  geom_hline(yintercept = 0, linewidth = 0.9, linetype = "dashed", color = "#4d4d4d") +
  geom_boxplot(width = 0.32, outlier.shape = NA, alpha = 0.76, color = "#2b2b2b", linewidth = 0.8) +
  geom_point(
    aes(color = cell_type),
    position = position_jitter(width = 0.075, height = 0, seed = 23),
    size = 3.9,
    alpha = 0.84
  ) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 5.5, fill = "#ffffff", color = "#111111", stroke = 1.0) +
  geom_point(
    data = cf_summary,
    aes(x = cell_type, y = arm_level_rho),
    inherit.aes = FALSE,
    shape = 24,
    size = 6.0,
    fill = "#111111",
    color = "#111111"
  ) +
  geom_text(
    data = cf_ann,
    aes(x = cell_type, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.9,
    fontface = "bold",
    lineheight = 0.94,
    color = "#1f1f1f"
  ) +
  geom_text(
    data = arm_ann,
    aes(x = cell_type, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.4,
    lineheight = 0.95,
    color = "#333333"
  ) +
  scale_y_continuous(limits = c(-0.10, 0.38), breaks = seq(-0.10, 0.35, 0.05), expand = c(0, 0)) +
  scale_fill_manual(values = c("GM12878 Dip-C" = "#2d6e9f", "Sperm scHi-C" = "#5c7c45")) +
  scale_color_manual(values = c("GM12878 Dip-C" = "#1d4f75", "Sperm scHi-C" = "#3f5f2f")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Spearman rho\nsequence similarity vs 3D proximity"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 15) +
  theme(
    axis.text.x = element_text(face = "bold", size = 18, color = "#1f1f1f"),
    axis.text.y = element_text(size = 13, color = "#333333"),
    axis.title.y = element_text(size = 12.5, face = "bold", color = "#222222", margin = margin(r = 12)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(26, 40, 24, 82)
  )

save_plot_pair(cf_plot, "community_free_rho_distribution_reduced_text")

whole <- read.table(gzfile(gm_3dg_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(whole) <- c("chrom_hap", "coord", "x", "y", "z")
whole <- whole %>%
  mutate(
    chrom = sub("\\(.*\\)$", "", chrom_hap),
    hap = ifelse(grepl("\\(mat\\)", chrom_hap), "mat", ifelse(grepl("\\(pat\\)", chrom_hap), "pat", "unknown")),
    chrom_factor = factor(chrom, levels = c(as.character(1:22), "X", "Y"))
  )

whole_summary <- whole %>%
  summarise(
    dataset = "GM12878",
    technology = "Dip-C / 3DG",
    cell = "gm12878_01",
    particles = n(),
    chromosomes = paste(sort(unique(chrom)), collapse = ","),
    source_3dg = gm_3dg_path
  )
write.table(
  whole_summary,
  file.path(plot_dir, "gm12878_cell01_whole_genome_3dg_projection_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

chrom_cols <- setNames(hcl.colors(24, "Dark 3"), c(as.character(1:22), "X", "Y"))
whole_plot <- ggplot(whole, aes(x = x, y = y, color = chrom_factor)) +
  geom_point(alpha = 0.18, size = 0.38, stroke = 0) +
  coord_fixed() +
  scale_color_manual(values = chrom_cols, na.value = "#666666") +
  labs(
    title = "GM12878 Dip-C whole-genome 3D coordinate projection",
    subtitle = "One existing reconstructed cell (gm12878_01); 2D view of x/y coordinates from 3DG particles",
    x = "3DG x coordinate",
    y = "3DG y coordinate",
    color = "Chromosome"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 20, color = "#1a3a6b"),
    plot.subtitle = element_text(size = 12, color = "#444444"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 24, 20, 24)
  )

save_plot_pair(whole_plot, "gm12878_cell01_whole_genome_3dg_projection")
