#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}

out_dir <- script_dir

dist_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_hic.dist_matrix.tsv"
seq_comm_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
region_bed_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_subtelomeric_regions.bed"
global_test_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_global_test.tsv"
seq_global_test_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_seq_global_test.tsv"

required <- c(dist_path, seq_comm_path, region_bed_path, global_test_path, seq_global_test_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
}

dist_df <- read.delim(dist_path, check.names = FALSE, row.names = 1)
dist_mat <- as.matrix(dist_df)
storage.mode(dist_mat) <- "numeric"
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

fit <- cmdscale(as.dist(dist_mat), eig = TRUE, k = 3)
points <- as.data.frame(fit$points)
colnames(points) <- c("D1", "D2", "D3")
points$arm <- rownames(points)

positive_eig <- fit$eig[fit$eig > 0]
variance <- fit$eig[1:3] / sum(positive_eig)

seq_comm <- read.delim(seq_comm_path, check.names = FALSE)
colnames(seq_comm)[1:3] <- c("arm", "community", "community_arms")
points <- merge(points, seq_comm, by = "arm", all.x = TRUE, sort = FALSE)
points$community[is.na(points$community)] <- "unassigned"
points$community_arms[is.na(points$community_arms)] <- points$arm[is.na(points$community_arms)]

multi_community <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C11", "C12")
points$display_community <- ifelse(
  points$community %in% multi_community,
  points$community,
  "singleton / no CHM13 partner"
)

community_levels <- c(multi_community, "singleton / no CHM13 partner")
points$display_community <- factor(points$display_community, levels = community_levels)

community_colors <- c(
  C1 = "#0072B2",
  C2 = "#D55E00",
  C3 = "#009E73",
  C4 = "#CC79A7",
  C5 = "#E69F00",
  C6 = "#56B4E9",
  C7 = "#C81D25",
  C11 = "#6A3D9A",
  C12 = "#555555",
  "singleton / no CHM13 partner" = "#B8B8B8"
)

make_projection <- function(df, y_col, projection_label, y_axis_label) {
  data.frame(
    arm = df$arm,
    D1 = df$D1,
    y = df[[y_col]],
    community = df$community,
    community_arms = df$community_arms,
    display_community = df$display_community,
    projection = projection_label,
    y_axis_label = y_axis_label,
    stringsAsFactors = FALSE
  )
}

plot_df <- rbind(
  make_projection(
    points,
    "D2",
    sprintf("D1-D2 projection (%.1f%% + %.1f%%)", 100 * variance[1], 100 * variance[2]),
    "MDS D2"
  ),
  make_projection(
    points,
    "D3",
    sprintf("D1-D3 projection (%.1f%% + %.1f%%)", 100 * variance[1], 100 * variance[3]),
    "MDS D3"
  )
)
plot_df$projection <- factor(plot_df$projection, levels = unique(plot_df$projection))
plot_df <- plot_df[order(plot_df$display_community == "singleton / no CHM13 partner"), ]

global_test <- read.delim(global_test_path, check.names = FALSE)
seq_global_test <- read.delim(seq_global_test_path, check.names = FALSE)
mantel <- global_test[global_test$test == "mantel", ]
seq_within <- seq_global_test[seq_global_test$test == "seq_within_vs_between", ]

metric_label <- sprintf(
  "n=38 arms; Mantel rho=%.3f, p<1e-4; sequence-community contact enrichment p=%.2g",
  mantel$U_statistic[1],
  seq_within$p_value[1]
)

wrap_text <- function(x, width = 150) {
  paste(strwrap(x, width = width), collapse = "\n")
}

caption_text <- paste(
  "Caption: CHM13 Hi-C contact-space MDS over PHR/subtelomeric arm regions from",
  "chm13_subtelomeric_regions.bed at 50 kb resolution, shown as D1-D2 and D1-D3",
  "projections of a 3D MDS. Colors mark sequence-defined subtelomeric communities;",
  "nearby points indicate similar bulk Hi-C contact profiles. This is contact-space",
  "MDS from bulk Hi-C, not a physical single-cell genome reconstruction."
)
caption_text <- wrap_text(caption_text, width = 118)
subtitle_text <- wrap_text(
  paste("3D MDS rendered as two readable projections; color = sequence community.", metric_label),
  width = 96
)

replacement_plot <- ggplot(
  plot_df,
  aes(x = D1, y = y, fill = display_community, color = display_community)
) +
  geom_hline(yintercept = 0, color = "#D7DCE2", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "#D7DCE2", linewidth = 0.6) +
  geom_point(shape = 21, size = 5.8, stroke = 1.0, alpha = 0.96) +
  geom_text(
    aes(label = arm),
    color = "#1E293B",
    size = 4.1,
    fontface = "bold",
    vjust = -1.05,
    check_overlap = TRUE
  ) +
  facet_wrap(~projection, nrow = 1) +
  scale_fill_manual(values = community_colors, drop = FALSE) +
  scale_color_manual(values = community_colors, drop = FALSE) +
  coord_equal(clip = "off") +
  labs(
    title = "CHM13 PHR Hi-C contact-space MDS",
    subtitle = subtitle_text,
    x = sprintf("MDS D1 (%.1f%% of positive eigenvalue variance)", 100 * variance[1]),
    y = "MDS projection coordinate",
    fill = "Sequence community",
    color = "Sequence community",
    caption = caption_text
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 22.5, color = "#183A68", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11.7, color = "#333333", margin = margin(b = 9), lineheight = 1.06),
    strip.text = element_text(face = "bold", size = 12.6, color = "#183A68"),
    axis.title = element_text(face = "bold", size = 13.8, color = "#111111"),
    axis.text = element_text(size = 11.5, color = "#333333"),
    panel.grid.major = element_line(color = "#E5EAF0", linewidth = 0.55),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.24, "in"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10.2),
    legend.key.width = unit(0.28, "in"),
    legend.key.height = unit(0.16, "in"),
    plot.caption = element_text(size = 8.2, color = "#2D2D2D", hjust = 0, lineheight = 1.05, margin = margin(t = 8)),
    plot.margin = margin(20, 42, 20, 42)
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4.8)),
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4.8))
  )

png_path <- file.path(out_dir, "best_replacement_chm13_phr_contact_mds.png")
pdf_path <- file.path(out_dir, "best_replacement_chm13_phr_contact_mds.pdf")

agg_png(png_path, width = 1920, height = 1080, units = "px", res = 144, background = "white")
print(replacement_plot)
dev.off()

ggsave(pdf_path, replacement_plot, width = 13.333, height = 7.5, units = "in", device = cairo_pdf)

coords_out <- points[, c("arm", "D1", "D2", "D3", "community", "community_arms", "display_community")]
write.table(
  coords_out,
  file.path(out_dir, "best_replacement_mds_coords.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

project_3d_oblique <- function(df) {
  # Oblique 3D projection for a slide-friendly view of the same D1/D2/D3 MDS.
  data.frame(
    arm = df$arm,
    px = 0.92 * df$D1 + 0.38 * df$D2,
    py = 0.92 * df$D3 + 0.42 * df$D2,
    community = df$community,
    community_arms = df$community_arms,
    display_community = df$display_community,
    stringsAsFactors = FALSE
  )
}

projection_3d <- project_3d_oblique(points)
axis_df <- rbind(
  data.frame(axis = "MDS D1", x = 0, y = 0, xend = 0.42 * 0.92, yend = 0),
  data.frame(axis = "MDS D2", x = 0, y = 0, xend = 0.42 * 0.38, yend = 0.42 * 0.42),
  data.frame(axis = "MDS D3", x = 0, y = 0, xend = 0, yend = 0.42 * 0.92)
)
axis_lab_df <- transform(axis_df, lx = xend * 1.10, ly = yend * 1.10)

plot_3d_view <- ggplot(
  projection_3d,
  aes(x = px, y = py, fill = display_community, color = display_community)
) +
  geom_segment(
    data = axis_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.10, "in"), type = "closed"),
    color = "#93A1B2",
    linewidth = 0.72
  ) +
  geom_text(
    data = axis_lab_df,
    aes(x = lx, y = ly, label = axis),
    inherit.aes = FALSE,
    color = "#475569",
    size = 4.0,
    fontface = "bold"
  ) +
  geom_point(shape = 21, size = 6.4, stroke = 1.05, alpha = 0.97) +
  geom_text(
    aes(label = arm),
    color = "#1E293B",
    size = 4.05,
    fontface = "bold",
    vjust = -1.0,
    check_overlap = TRUE
  ) +
  scale_fill_manual(values = community_colors, drop = FALSE) +
  scale_color_manual(values = community_colors, drop = FALSE) +
  coord_equal(clip = "off") +
  labs(
    title = "CHM13 PHR Hi-C contact-space MDS: 3D view",
    subtitle = wrap_text(
      paste(
        "The same D1-D2-D3 contact-space MDS as the projection slide, rendered in a",
        "single oblique 3D view. Colors match the sequence communities used throughout."
      ),
      width = 104
    ),
    x = NULL,
    y = NULL,
    fill = "Sequence community",
    color = "Sequence community",
    caption = wrap_text(
      paste(
        "Bulk Hi-C contact-space MDS over CHM13 PHR/subtelomeric arm regions at 50 kb.",
        "This is a visual embedding of contact profiles, not a physical single-cell genome reconstruction."
      ),
      width = 124
    )
  ) +
  theme_void(base_family = "DejaVu Sans", base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 24, color = "#183A68", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12.8, color = "#333333", margin = margin(b = 10), lineheight = 1.05),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10.2),
    legend.key.width = unit(0.28, "in"),
    legend.key.height = unit(0.16, "in"),
    plot.caption = element_text(size = 9.2, color = "#2D2D2D", hjust = 0, lineheight = 1.06, margin = margin(t = 7)),
    plot.margin = margin(24, 42, 20, 42)
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4.8)),
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4.8))
  )

png_3d_path <- file.path(out_dir, "chm13_phr_contact_mds_3d_view.png")
pdf_3d_path <- file.path(out_dir, "chm13_phr_contact_mds_3d_view.pdf")

agg_png(png_3d_path, width = 1920, height = 1080, units = "px", res = 144, background = "white")
print(plot_3d_view)
dev.off()

ggsave(pdf_3d_path, plot_3d_view, width = 13.333, height = 7.5, units = "in", device = cairo_pdf)

metrics <- data.frame(
  metric = c(
    "n_arms",
    "D1_positive_eigen_variance",
    "D2_positive_eigen_variance",
    "D3_positive_eigen_variance",
    "mantel_rho",
    "mantel_p_display",
    "sequence_community_contact_enrichment_p",
    "source_distance_matrix",
    "source_sequence_community_assignments",
    "source_region_bed"
  ),
  value = c(
    nrow(points),
    sprintf("%.6f", variance[1]),
    sprintf("%.6f", variance[2]),
    sprintf("%.6f", variance[3]),
    sprintf("%.6f", mantel$U_statistic[1]),
    "<1e-4",
    sprintf("%.8g", seq_within$p_value[1]),
    dist_path,
    seq_comm_path,
    region_bed_path
  )
)

write.table(
  metrics,
  file.path(out_dir, "best_replacement_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

candidate_inventory <- data.frame(
  candidate_id = c(
    "v8_current_chm13_whole_arm_contact_mds",
    "v8_optional_gm12878_cell01_whole_genome_3dg",
    "hg002_whole_arm_contact_mds",
    "submission_randiak_mds_3d",
    "submission_randiak_mds_3d_q_only",
    "gm12878_radial_community_overlay",
    "sperm_radial_community_overlay",
    "chm13_phr_contact_mds_50kb",
    "hg002_phr_contact_mds_50kb",
    "chm13_flanking_contact_mds_50kb",
    "cross_organism_overlay_summaries"
  ),
  local_path = c(
    "slides/v2-review-zoom/_revision_assets/v8/hic_dipc_clarity_split/inputs/chm13_hic_mds_3d_coords.png; /moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png",
    "slides/v2-review-zoom/_revision_assets/v8/hic_dipc_clarity_split/plots/gm12878_cell01_whole_genome_3dg_projection.png; /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/gm12878_01.impute3.round4.clean.3dg.gz",
    "/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/plots/MDS_3d_coords.png",
    "/moosefs/guarracino/HPRCv2/submission_Randiak/images/mds_3d.png",
    "/moosefs/guarracino/HPRCv2/submission_Randiak/images/mds_3d_q.png",
    "slides/v2-review-zoom/_revision_assets/v8/hic_dipc_clarity_split/inputs/gm12878_radial_community.png; /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf",
    "slides/v2-review-zoom/_revision_assets/v8/hic_dipc_clarity_split/inputs/sperm_all20_radial_community.png; /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf",
    "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_hic.dist_matrix.tsv; /moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_hic_mds_comparison.pdf",
    "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_hic.dist_matrix.tsv; /moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_hic_mds_comparison.pdf",
    "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/flanking/50000bp/chm13_hic_mds_comparison.pdf",
    "/moosefs/guarracino/HPRCv2/PHR_III/analysis/cross_organism/*.pdf; /moosefs/guarracino/HPRCv2/PHR_III/analysis/cross_organism/*.tsv"
  ),
  dataset = c(
    "CHM13",
    "GM12878",
    "HG002",
    "HG002",
    "HG002 q arms",
    "GM12878",
    "human sperm",
    "CHM13",
    "HG002",
    "CHM13",
    "human/Dip-C/sperm overlays"
  ),
  measurement = c(
    "Bulk Hi-C contact-space MDS over whole p-arm, centromere, q-arm regions",
    "Single-cell Dip-C physical 3DG x/y projection for one cell",
    "Bulk Hi-C contact-space MDS over diploid whole-arm/centromere regions",
    "Submission Randiak whole-arm/centromere contact-space MDS",
    "Submission Randiak q-arm-only contact-space MDS",
    "Normalized radial community position overlay",
    "Normalized radial community position overlay",
    "Bulk Hi-C contact-space MDS over CHM13 PHR/subtelomeric arm regions",
    "Bulk Hi-C contact-space MDS over HG002 PHR/subtelomeric arm regions",
    "Bulk Hi-C contact-space MDS over CHM13 flanking regions",
    "Cross-dataset overlay summaries, not coordinate structure plots"
  ),
  audit_decision = c(
    "Current slide source; visually hard to read because labels are tiny and the plot is mostly unused whitespace; not PHR-specific.",
    "Clearer physical-coordinate context, but not CHM13, not Hi-C, n=1 cell, and not a PHR/community validation plot.",
    "Not better for slide 10m.2: more labels and diploid haplotype text make it less readable.",
    "Not better: same whole-arm contact-space caveat and dense labels; older HG002 submission asset.",
    "Clearer subset but q-only HG002; too narrow and not a CHM13 replacement.",
    "Useful validation summary already handled by slide 11a; radial statistic, not a whole 3D structure or MDS replacement.",
    "Useful validation summary already handled by slide 11a; radial statistic, not a whole 3D structure or MDS replacement.",
    "Best replacement: CHM13, PHR-level, reproducible from local TSVs, stronger contrast, larger arm/community labels, and explicit contact-space caveat.",
    "Good backup but diploid/HG002 and less direct than the CHM13 slide target.",
    "Good methods backup for unique flanks, but less direct than PHR/community slide 10m.2.",
    "Useful for synthesis/robustness, but not a slide-ready 3D visualization."
  ),
  replacement_rank = c(6, 3, 7, 8, 5, 4, 4, 1, 2, 3, 9),
  caveat = c(
    "Contact-space MDS, not physical reconstruction; regions are whole p/c/q, not PHR or 500 kb windows.",
    "Physical 3DG projection, but one GM12878 cell only; cannot stand in for CHM13 Hi-C.",
    "Contact-space MDS, not physical reconstruction; whole-arm/centromere regions.",
    "Contact-space MDS, not physical reconstruction; HG002 submission context.",
    "Contact-space MDS, not physical reconstruction; q-arm subset only.",
    "Radial overlay, not 3D coordinate reconstruction.",
    "Radial overlay, not 3D coordinate reconstruction.",
    "Contact-space MDS from bulk Hi-C; not physical single-cell reconstruction.",
    "Contact-space MDS from bulk Hi-C; not physical single-cell reconstruction.",
    "Contact-space MDS from bulk Hi-C flanking regions; not physical single-cell reconstruction.",
    "Overlay summaries, not spatial coordinate plots."
  )
)

write.table(
  candidate_inventory,
  file.path(out_dir, "candidate_inventory.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote ", png_path)
message("Wrote ", pdf_path)
