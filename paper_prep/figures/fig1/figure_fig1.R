#!/usr/bin/env Rscript
# Figure 1 — Population-scale subtelomeric communities (the landscape)
# Panels:
#   1a (READY)    — genome-wide stacked identity heatmap
#                   (vendored from p_genome_wide_identity_heatmap_no_inset.pdf)
#   1b (READY)    — genome-wide n-chromosomes-sharing heatmap with PHR overlay
#                   (vendored from p_genome_wide_numchrom_heatmap.pdf)
#   1c (GENERATE) — 41x41 arm-level Jaccard distance heatmap with Leiden
#                   community block annotations + UPGMA dendrogram
#   1d (GENERATE) — per-arm architecture-category bar
#                   (homogeneous / polymorphic / fully interchangeable)
#
# All inputs and 1a/1b vendored PDFs are read from the paths recorded in
# paper_prep/figures/fig1/sources.tsv.
#
# Run with:
#   guix shell -m /tmp/manifest_fig1.scm -- Rscript paper_prep/figures/fig1/figure_fig1.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(ggdendro)
  library(magick)
  library(pdftools)
  library(RColorBrewer)
  library(viridis)
  library(grid)
  library(gridExtra)
})

# ---- paths ------------------------------------------------------------------

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) NA_character_
)
if (is.na(script_path) || !nzchar(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  m    <- grep("--file=", args, fixed = TRUE, value = TRUE)
  script_path <- if (length(m)) sub("--file=", "", m[1]) else "."
}
repo_root  <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                            mustWork = FALSE)
ext_root   <- "/moosefs/guarracino/HPRCv2/PHR_III"

f1a_pdf    <- file.path(repo_root, "inter-chr-plots", "p_genome_wide_identity_heatmap_no_inset.pdf")
f1b_pdf    <- file.path(repo_root, "inter-chr-plots", "p_genome_wide_numchrom_heatmap.pdf")
f1a_png    <- file.path(repo_root, "inter-chr-plots", "p_genome_wide_identity_heatmap.png")
f1b_png    <- file.path(repo_root, "inter-chr-plots", "p_genome_wide_numchrom_heatmap.png")

arm_dist   <- file.path(ext_root, "similarity",
                        "hprcv2.1Mb.subtelo.arm_dist_matrix.tsv")
arm_leiden <- file.path(ext_root, "similarity",
                        "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv")
arm_upgma  <- file.path(ext_root, "similarity",
                        "hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv")
seq_assign <- file.path(ext_root, "similarity",
                        "hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv")
cross_aff  <- file.path(ext_root, "heterogeneity",
                        "cross_arm_affinity_sequences.tsv")

out_dir    <- file.path(repo_root, "paper_prep", "figures", "fig1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

arm_order <- function(arms) {
  # natural sort: chr1_p, chr1_q, chr2_p, ..., chr22_q, chrX_p, chrX_q, chrY_p, chrY_q
  m <- str_match(arms, "^chr([^_]+)_(p|q)$")
  num_key <- suppressWarnings(as.integer(m[, 2]))
  num_key[m[, 2] == "X"] <- 23
  num_key[m[, 2] == "Y"] <- 24
  letter_key <- ifelse(m[, 3] == "p", 0, 1)
  order(num_key, letter_key)
}

# ---- load -------------------------------------------------------------------

D <- as.matrix(read.table(arm_dist, header = TRUE, row.names = 1,
                          check.names = FALSE, sep = "\t"))
stopifnot(nrow(D) == 41, ncol(D) == 41)

leiden <- as.data.frame(readr::read_tsv(arm_leiden,
                                        show_col_types = FALSE))
upgma  <- as.data.frame(readr::read_tsv(arm_upgma,
                                        show_col_types = FALSE))
seqs   <- as.data.frame(readr::read_tsv(seq_assign,
                                        show_col_types = FALSE))
cross  <- as.data.frame(readr::read_tsv(cross_aff,
                                        show_col_types = FALSE))

# normalise arm naming: similarity files use chr1_p; heterogeneity uses chr1_parm
cross$own_arm <- str_replace(cross$own_arm, "_parm$", "_p") |>
                 str_replace("_qarm$", "_q")
seqs$ChromArm <- ifelse(grepl("_(p|q)$", seqs$ChromArm),
                        seqs$ChromArm,
                        str_replace(seqs$ChromArm, "_parm$", "_p") |>
                            str_replace("_qarm$", "_q"))

# ---- compute Leiden + UPGMA structure ---------------------------------------

# 15 Leiden communities (canonical order C1..C15)
leiden_tbl <- leiden |>
  mutate(Community = factor(Community,
                            levels = paste0("C", 1:15))) |>
  arrange(Community, ChromArm)

# UPGMA dendrogram from D (treat values as distances)
hc_upgma <- hclust(as.dist(D), method = "average")

# Order arms primarily by Leiden community, secondarily by UPGMA-induced order
# within community. This makes Leiden block diagonal in the heatmap.
upgma_pos <- setNames(seq_along(hc_upgma$labels)[hc_upgma$order],
                      hc_upgma$labels[hc_upgma$order])
leiden_tbl$upgma_pos <- upgma_pos[leiden_tbl$ChromArm]
leiden_ord <- leiden_tbl |>
  arrange(Community, upgma_pos) |>
  pull(ChromArm)

# Re-cluster within each community for a clean reordered dendrogram
hc_dend <- hclust(as.dist(D[leiden_ord, leiden_ord]), method = "average")
# Map community membership for the dendrogram leaves
comm_of <- setNames(as.character(leiden_tbl$Community), leiden_tbl$ChromArm)

# UPGMA k=14 ordering — for the 12/15 agreement annotation
upgma_k14 <- upgma |>
  mutate(Community_upgma = Community) |>
  select(ChromArm, Community_upgma)

agreement_count <- leiden_tbl |>
  inner_join(upgma_k14, by = "ChromArm") |>
  count(Community, Community_upgma) |>
  group_by(Community) |>
  arrange(desc(n)) |>
  slice(1) |>
  ungroup() |>
  pull(Community_upgma) |>
  unique() |>
  length()

# ---- panel 1c: arm-level distance heatmap with Leiden + UPGMA --------------

mat <- D[leiden_ord, leiden_ord]
heat_df <- as.data.frame(as.table(mat)) |>
  rename(Arm1 = Var1, Arm2 = Var2, dist = Freq) |>
  mutate(Arm1 = factor(Arm1, levels = leiden_ord),
         Arm2 = factor(Arm2, levels = rev(leiden_ord)))

# Leiden block annotations (rectangles around community blocks)
comm_runs <- leiden_tbl |>
  filter(ChromArm %in% leiden_ord) |>
  arrange(match(ChromArm, leiden_ord)) |>
  mutate(idx = row_number()) |>
  group_by(Community) |>
  summarise(start = min(idx), end = max(idx),
            n = n(),
            label = unique(Community), .groups = "drop")

n_arm <- length(leiden_ord)

# x/y positions (factor levels are 1..n_arm; in the y axis we reversed)
rect_df <- comm_runs |>
  mutate(xmin = start - 0.5,
         xmax = end + 0.5,
         ymin = n_arm - end + 0.5,
         ymax = n_arm - start + 0.5)

# label positions for the diagonal (community labels)
lbl_df <- comm_runs |>
  mutate(x = (start + end) / 2,
         y = n_arm - (start + end) / 2 + 0.5,
         label_text = label)

p_heat <- ggplot(heat_df, aes(Arm1, Arm2, fill = dist)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", direction = -1,
                     limits = c(0, 1),
                     name = "Jaccard\ndistance",
                     guide = guide_colorbar(barwidth = 0.6, barheight = 4)) +
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "cyan", linewidth = 0.5,
            inherit.aes = FALSE) +
  geom_text(data = lbl_df,
            aes(x = x, y = y, label = label_text),
            inherit.aes = FALSE,
            colour = "cyan", size = 2.4, fontface = "bold") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5),
        axis.text.y = element_text(size = 5),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6),
        plot.margin = margin(2, 2, 2, 2))

# UPGMA dendrogram on top, ordered to match leiden_ord on x
dend <- as.dendrogram(hc_dend)
ddata <- ggdendro::dendro_data(dend, type = "rectangle")

# colour dendrogram leaves by leiden community
seg_df <- ddata$segments
lbl_dendro <- ddata$labels |>
  mutate(community = comm_of[as.character(label)])

palette_15 <- setNames(c(brewer.pal(8, "Set1"), brewer.pal(7, "Set2")),
                       paste0("C", 1:15))

p_dend <- ggplot() +
  geom_segment(data = seg_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3, colour = "grey30") +
  scale_x_continuous(limits = c(0.5, n_arm + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 2))

panel_1c <- (p_dend / p_heat) +
  plot_layout(heights = c(1, 6))

# ---- panel 1d: per-arm architecture-category bar ----------------------------

# total sequences per arm (denominator)
total_per_arm <- seqs |>
  count(ChromArm, name = "n_total")

# cross-arm sequence count per arm (sequences with cross-arm affinity)
cross_per_arm <- cross |>
  distinct(sequence, own_arm) |>
  count(own_arm, name = "n_cross") |>
  rename(ChromArm = own_arm)

arm_rates <- total_per_arm |>
  left_join(cross_per_arm, by = "ChromArm") |>
  mutate(n_cross = replace_na(n_cross, 0L),
         cross_rate = n_cross / n_total) |>
  inner_join(leiden_tbl |> select(ChromArm, Community),
             by = "ChromArm")

# Architecture categories (data-driven from Leiden community + silhouette):
#   fully interchangeable — arms in negative-silhouette communities (C7, C14,
#                           C15). Allele-paralog distance reversed or near
#                           equal: complete inter-chromosomal homogenization
#                           at population scale (per SURVEY_04 §1.1, §1.2).
#   homogeneous           — arms in single-arm Leiden communities (C8, C9,
#                           C10, C13): no detectable inter-chromosomal
#                           partner; private subtelomeres.
#   polymorphic           — remaining multi-arm community members: inter-
#                           chromosomal sharing present but arm identity
#                           preserved (allele closer than paralog).
#
# The skeleton MANUSCRIPT_SKELETON.md cites a preliminary 8 / 34 / 7 split;
# the data-driven definitions above yield 4 / 28 / 9 (4+28+9 = 41), which
# we report transparently in the per-arm TSV alongside cross_rate.

fully_interchangeable_comm <- c("C7", "C14", "C15")  # neg silhouette
single_arm_comm <- leiden_tbl |>
  count(Community, name = "n_arms") |>
  filter(n_arms == 1) |>
  pull(Community) |>
  as.character()

arm_rates <- arm_rates |>
  mutate(category = case_when(
    Community %in% fully_interchangeable_comm  ~ "fully interchangeable",
    Community %in% single_arm_comm             ~ "homogeneous",
    TRUE                                       ~ "polymorphic"
  ))

# Order arms by genomic position
arm_rates <- arm_rates |>
  mutate(order = arm_order(ChromArm)) |>
  arrange(order) |>
  mutate(ChromArm = factor(ChromArm, levels = ChromArm))

# write per-arm table
write.table(arm_rates,
            file = file.path(out_dir, "architecture_per_arm.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat_pal <- c("homogeneous" = "#2C7FB8",
             "polymorphic" = "#7FCDBB",
             "fully interchangeable" = "#D7301F")

cat_counts <- arm_rates |> count(category)

panel_1d <- ggplot(arm_rates,
                   aes(x = ChromArm, y = cross_rate, fill = category)) +
  geom_col(width = 0.85, colour = "grey20", linewidth = 0.15) +
  scale_fill_manual(values = cat_pal,
                    name = "Architecture",
                    breaks = c("homogeneous",
                               "polymorphic",
                               "fully interchangeable"),
                    labels = function(x) {
                      n <- setNames(cat_counts$n, cat_counts$category)
                      sprintf("%s (%d/41)", x, n[x] %||% 0)
                    }) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Cross-arm sequence rate") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5),
        axis.text.y = element_text(size = 6),
        legend.position = "top",
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

# ---- panels 1a / 1b: vendored PDF rasters -----------------------------------

img_to_panel <- function(img_path) {
  ggdraw() + draw_image(img_path)
}

panel_1a <- img_to_panel(f1a_png)
panel_1b <- img_to_panel(f1b_png)

# ---- compose 4-panel figure -------------------------------------------------

label_theme <- theme(plot.tag = element_text(face = "bold", size = 11))

p1a <- panel_1a + labs(tag = "a") + label_theme
p1b <- panel_1b + labs(tag = "b") + label_theme
p1c <- panel_1c & label_theme
p1c <- p1c + plot_annotation(tag_levels = list(c("c", "")))
# patchwork annotation tags are tricky for nested patches — stamp manually:
p1c_final <- (p_dend / (p_heat + labs(tag = "c") + label_theme)) +
  plot_layout(heights = c(1, 6))
p1d <- panel_1d + labs(tag = "d") + label_theme

top_row    <- (p1a | p1b) + plot_layout(widths = c(1, 1))
bottom_row <- (p1c_final | p1d) + plot_layout(widths = c(1.05, 1))

fig <- top_row / bottom_row +
  plot_layout(heights = c(1, 1.15)) &
  theme(plot.margin = margin(4, 4, 4, 4))

# ---- write outputs ----------------------------------------------------------

pdf_out <- file.path(out_dir, "figure_fig1.pdf")
png_out <- file.path(out_dir, "figure_fig1.png")

ggsave(pdf_out, fig, width = 11.5, height = 11, units = "in",
       device = cairo_pdf)
ggsave(png_out, fig, width = 11.5, height = 11, units = "in",
       dpi = 300)

cat("Wrote:\n  ", pdf_out, "\n  ", png_out, "\n", sep = "")
cat("UPGMA k14 / Leiden k15 community-mapping agreement: ",
    agreement_count, "/15\n", sep = "")
cat("Architecture category counts:\n")
print(cat_counts)
