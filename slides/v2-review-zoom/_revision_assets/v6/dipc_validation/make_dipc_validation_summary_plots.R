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

gm_summary_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_summary.tsv"
gm_per_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_cell.tsv"
gm_per_community_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_community_per_cell.tsv"
gm_cf_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_per_cell.tsv"
gm_cf_arm_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_arm.tsv"

sperm_summary_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_summary.tsv"
sperm_per_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_cell.tsv"
sperm_per_community_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv"
sperm_cf_cell_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_per_cell.tsv"
sperm_cf_arm_path <- "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_arm.tsv"

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

if (nrow(gm_s_all) != 16) stop("Expected 16 GM12878 S_all rows, found ", nrow(gm_s_all))
if (nrow(sperm_s_all) != 20) stop("Expected 20 sperm S_all rows, found ", nrow(sperm_s_all))

wb_plot_df <- bind_rows(
  gm_cell %>% transmute(cell_type = "GM12878 Dip-C", signal = "C-community", cell_id = as.character(cell_id), ratio),
  gm_s_all %>% transmute(cell_type = "GM12878 Dip-C", signal = "S_all negative control", cell_id = as.character(cell_id), ratio),
  sperm_cell %>% transmute(cell_type = "Sperm scHi-C", signal = "C-community", cell_id = as.character(cell_id), ratio),
  sperm_s_all %>% transmute(cell_type = "Sperm scHi-C", signal = "S_all negative control", cell_id = as.character(cell_id), ratio)
) %>%
  mutate(
    cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")),
    signal = factor(signal, levels = c("C-community", "S_all negative control"))
  )

sperm_s_all_pooled <- sum(sperm_s_all$mean_within * sperm_s_all$n_pairs) /
  sum(sperm_s_all$overall_between * sperm_s_all$n_pairs)

wb_summary <- tibble(
  cell_type = c("GM12878 Dip-C", "GM12878 Dip-C", "Sperm scHi-C", "Sperm scHi-C"),
  signal = c("C-community", "S_all negative control", "C-community", "S_all negative control"),
  display_ratio = c(gm_fisher$mean_ratio, 1.106, sperm_fisher$mean_ratio, 1.397),
  per_cell_mean_ratio = c(mean(gm_cell$ratio), mean(gm_s_all$ratio), mean(sperm_cell$ratio), mean(sperm_s_all$ratio)),
  n_cells = c(nrow(gm_cell), nrow(gm_s_all), nrow(sperm_cell), nrow(sperm_s_all)),
  n_below_1 = c(sum(gm_cell$ratio < 1), sum(gm_s_all$ratio < 1), sum(sperm_cell$ratio < 1), sum(sperm_s_all$ratio < 1)),
  effect_label = c("6.9% closer", "11% farther", "60% closer", "40% farther"),
  p_label = c("Fisher p=2.4e-05", "0/16 cells < 1", "Fisher p=3.9e-51", "1/20 cells < 1"),
  source_tsv = c(gm_per_cell_path, gm_per_community_path, sperm_per_cell_path, sperm_per_community_path)
) %>%
  mutate(
    cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")),
    signal = factor(signal, levels = c("C-community", "S_all negative control"))
  )

write.table(
  wb_summary,
  file.path(plot_dir, "wb_negative_control_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

wb_ann <- wb_summary %>%
  mutate(
    y = ifelse(signal == "C-community", 1.18, 1.78),
    label = paste0(
      "W/B=", sprintf("%.3f", display_ratio), "\n",
      effect_label, "\n",
      p_label
    )
  )

wb_plot <- ggplot(wb_plot_df, aes(x = signal, y = ratio, fill = signal)) +
  geom_hline(yintercept = 1, linewidth = 1.1, linetype = "dashed", color = "#4d4d4d") +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.78, color = "#2b2b2b", linewidth = 0.85) +
  geom_point(
    aes(color = signal),
    position = position_jitter(width = 0.11, height = 0, seed = 19),
    size = 4.8,
    alpha = 0.85,
    stroke = 0.2
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 6.5, fill = "#ffffff", color = "#111111", stroke = 1.1) +
  geom_text(
    data = wb_ann,
    aes(x = signal, y = y, label = label),
    inherit.aes = FALSE,
    size = 7.2,
    fontface = "bold",
    lineheight = 0.92,
    color = "#1f1f1f"
  ) +
  facet_wrap(~ cell_type, nrow = 1) +
  scale_x_discrete(labels = c("C-community" = "C-community\nsharing", "S_all negative control" = "S_all\nnon-sharing")) +
  scale_y_continuous(limits = c(0.05, 2.05), breaks = seq(0.25, 2.0, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c("C-community" = "#2d6e9f", "S_all negative control" = "#d68627")) +
  scale_color_manual(values = c("C-community" = "#1d4f75", "S_all negative control" = "#99510d")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "W/B ratio"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 23) +
  theme(
    strip.text = element_text(face = "bold", size = 28, color = "#1a3a6b"),
    axis.text.x = element_text(face = "bold", size = 21, color = "#1f1f1f", lineheight = 0.9),
    axis.text.y = element_text(size = 20, color = "#333333"),
    axis.title.y = element_text(size = 25, face = "bold", color = "#222222", margin = margin(r = 12)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(44, 56, 36, 54)
  )

agg_png(file.path(plot_dir, "wb_negative_control.png"), width = 1920, height = 1080, units = "px", res = 144, background = "white")
print(wb_plot)
dev.off()

cf_plot_df <- bind_rows(
  gm_cf %>% transmute(cell_type = "GM12878 Dip-C", cell_id = as.character(cell_id), spearman_rho),
  sperm_cf %>% transmute(cell_type = "Sperm scHi-C", cell_id = as.character(cell_id), spearman_rho)
) %>%
  mutate(cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")))

arm_cor <- function(df) {
  ct <- suppressWarnings(cor.test(df$mean_jaccard, -df$mean_3d_distance, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), p_value = ct$p.value, n_pairs = nrow(df))
}

cf_summary <- bind_rows(
  tibble(
    cell_type = "GM12878 Dip-C",
    median_rho = median(gm_cf$spearman_rho),
    n_positive = sum(gm_cf$spearman_rho > 0),
    n_cells = nrow(gm_cf),
    source_tsv = gm_cf_cell_path
  ) %>% bind_cols(arm_cor(gm_cf_arm) %>% rename(arm_level_rho = rho, arm_level_p = p_value, arm_level_n_pairs = n_pairs)),
  tibble(
    cell_type = "Sperm scHi-C",
    median_rho = median(sperm_cf$spearman_rho),
    n_positive = sum(sperm_cf$spearman_rho > 0),
    n_cells = nrow(sperm_cf),
    source_tsv = sperm_cf_cell_path
  ) %>% bind_cols(arm_cor(sperm_cf_arm) %>% rename(arm_level_rho = rho, arm_level_p = p_value, arm_level_n_pairs = n_pairs))
) %>%
  mutate(cell_type = factor(cell_type, levels = c("GM12878 Dip-C", "Sperm scHi-C")))

write.table(
  cf_summary,
  file.path(plot_dir, "community_free_rho_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cf_ann <- cf_summary %>%
  mutate(
    y = c(0.255, 0.255),
    label = paste0(
      "median rho=", sprintf("%.3f", median_rho), "\n",
      n_positive, "/", n_cells, " cells positive"
    )
  )

arm_ann <- cf_summary %>%
  mutate(
    y = ifelse(cell_type == "GM12878 Dip-C", arm_level_rho + 0.023, arm_level_rho - 0.03),
    label = ifelse(
      cell_type == "GM12878 Dip-C",
      paste0("arm-level rho=", sprintf("%.3f", arm_level_rho), "\np=1.1e-18"),
      paste0("arm-level pooled rho=", sprintf("%.3f", arm_level_rho), "\np=0.197 caveat")
    )
  )

cf_plot <- ggplot(cf_plot_df, aes(x = cell_type, y = spearman_rho, fill = cell_type)) +
  geom_hline(yintercept = 0, linewidth = 1.1, linetype = "dashed", color = "#4d4d4d") +
  geom_boxplot(width = 0.34, outlier.shape = NA, alpha = 0.76, color = "#2b2b2b", linewidth = 0.9) +
  geom_point(
    aes(color = cell_type),
    position = position_jitter(width = 0.08, height = 0, seed = 23),
    size = 5,
    alpha = 0.88
  ) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 7.2, fill = "#ffffff", color = "#111111", stroke = 1.2) +
  geom_point(
    data = cf_summary,
    aes(x = cell_type, y = arm_level_rho),
    inherit.aes = FALSE,
    shape = 24,
    size = 7.8,
    fill = "#111111",
    color = "#111111"
  ) +
  geom_text(
    data = cf_ann,
    aes(x = cell_type, y = y, label = label),
    inherit.aes = FALSE,
    size = 7.2,
    fontface = "bold",
    lineheight = 0.92,
    color = "#1f1f1f"
  ) +
  geom_text(
    data = arm_ann,
    aes(x = cell_type, y = y, label = label),
    inherit.aes = FALSE,
    size = 6.3,
    lineheight = 0.92,
    color = "#333333"
  ) +
  scale_y_continuous(limits = c(-0.105, 0.38), breaks = seq(-0.10, 0.35, 0.05), expand = c(0, 0)) +
  scale_fill_manual(values = c("GM12878 Dip-C" = "#2d6e9f", "Sperm scHi-C" = "#5c7c45")) +
  scale_color_manual(values = c("GM12878 Dip-C" = "#1d4f75", "Sperm scHi-C" = "#3f5f2f")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 23) +
  theme(
    axis.text.x = element_text(face = "bold", size = 26, color = "#1f1f1f"),
    axis.text.y = element_text(size = 20, color = "#333333"),
    axis.title.y = element_text(size = 25, face = "bold", color = "#222222", margin = margin(r = 12)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(44, 58, 38, 54)
  )

agg_png(file.path(plot_dir, "community_free_rho_distribution.png"), width = 1920, height = 1080, units = "px", res = 144, background = "white")
print(cf_plot)
dev.off()
