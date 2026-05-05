#!/usr/bin/env Rscript
# Figure 4 — Pedigree-resolved exchanges + cross-species mouse generalisation
# 4a: WashU 3-gen patch landscape (PAN027 mat, PAN027 pat, PAN028 mat)
# 4b: CEPH1463 cross-assembler 11-feature parent-by-pair matrix
# 4c: RPE-1 t(X;10) sequence (Jaccard) and 3D (Hi-C contact) heatmap
# 4d: Mouse zygotene per-PHR-pair similarity vs Hi-C contact scatter

suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(data.table); library(scales); library(RColorBrewer)
})

OUT_DIR <- "paper_prep/figures/fig4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

melt_matrix <- function(m, varnames = c("a", "b"), value.name = "value") {
  dt <- as.data.table(m, keep.rownames = varnames[1])
  out <- melt(dt, id.vars = varnames[1], variable.name = varnames[2], value.name = value.name)
  out[[varnames[1]]] <- as.character(out[[varnames[1]]])
  out[[varnames[2]]] <- as.character(out[[varnames[2]]])
  out
}

# Inputs
PATCHES   <- "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv"
LEIDEN    <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
RPE_DIST  <- "/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/rpe1.dist_matrix.tsv"
RPE_CONT  <- "/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/hic_validation/50000bp/rpe1_self_async_cifi_contact_matrix.tsv"
MOUSE_TSV <- "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_correlation.tsv"

CHR_LEVELS <- paste0("chr", c(1:22, "X"))
ARM_LEVELS <- as.vector(rbind(paste0(CHR_LEVELS, "_p"), paste0(CHR_LEVELS, "_q")))

# Leiden community palette (15 communities + grey "unknown")
leiden <- fread(LEIDEN)
comm_colors <- c(
  C1  = "#a6cee3", C2  = "#1f78b4", C3  = "#b2df8a", C4  = "#33a02c",
  C5  = "#fb9a99", C6  = "#e31a1c", C7  = "#fdbf6f", C8  = "#ff7f00",
  C9  = "#cab2d6", C10 = "#6a3d9a", C11 = "#ffff99", C12 = "#b15928",
  C13 = "#8dd3c7", C14 = "#bebada", C15 = "#fb8072",
  unknown = "#d9d9d9"
)

theme_fig <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1),
      panel.spacing = unit(0.4, "lines"),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text  = element_text(size = base_size - 2),
      axis.text  = element_text(size = base_size - 2),
      axis.title = element_text(size = base_size - 1)
    )
}

# ---------- Panel 4a: WashU 3-gen patch landscape ----------
make_panel_4a <- function() {
  p <- fread(PATCHES)
  w <- p[ds == "WashU"]
  # Inter-chr only
  w[, ref_chr := sub("[pq]$", "", ref_chrarm)]
  w[, query_full_arm := paste0(query_chr, "_", sub("arm$", "", query_arm))]
  w[, query_full_arm := factor(query_full_arm, levels = ARM_LEVELS)]
  w_inter <- w[ref_chr != query_chr]
  w_inter[, status := factor(community_status,
                              levels = c("within_community", "cross_community", "unknown_community"),
                              labels = c("within community", "cross community", "unknown"))]
  # Use query community if known else 'unknown'
  w_inter[, comm_disp := ifelse(query_community %in% names(comm_colors) & query_community != "",
                                 query_community, "unknown")]
  w_inter[query_community == "unknown" | query_community == "", comm_disp := "unknown"]

  # Order labels chronologically (G2 mat, G2 pat, G3 mat)
  label_levels <- c(
    "PAN027 maternal (hap1) vs PAN010 (mother)",
    "PAN027 paternal (hap2) vs PAN011 (father)",
    "PAN028 maternal (hap1) vs PAN027 (mother)"
  )
  label_short <- c("PAN027 mat\n← PAN010", "PAN027 pat\n← PAN011", "PAN028 mat\n← PAN027")
  w_inter[, label_f := factor(label, levels = label_levels, labels = label_short)]

  totals <- w_inter[, .N, by = .(label_f, status)][, frac := N / sum(N), by = label_f]
  within_totals <- w_inter[, .(N = .N, within = sum(status == "within community")), by = label_f]
  within_totals[, label_pct := sprintf("%d/%d HQ inter-chr\n(%.0f%% within Leiden)",
                                        within, N, 100 * within / N)]

  # Plot: x = arm, y = patch position, point = patch midpoint, color = comm
  w_inter[, mid := (patch_start + patch_end) / 2]
  w_inter[, comm_disp := factor(comm_disp, levels = names(comm_colors))]

  ggplot(w_inter, aes(x = query_full_arm, y = mid / 1000, color = comm_disp)) +
    geom_point(data = w_inter[status == "within community"], size = 1.4, alpha = 0.85) +
    geom_point(data = w_inter[status != "within community"], size = 0.7, alpha = 0.4, color = "grey60") +
    scale_color_manual(values = comm_colors, na.value = "grey70",
                       breaks = c("C2","C3","C5","C6","C7","C11","C13","C15","unknown"),
                       labels = c("C2 (chr10p/chr18p)","C3 (chr3q/chr9q)","C5 (chr6p/chr9p)",
                                  "C6 (chr1q/chr13q)","C7 (acrocentric p)","C11 (chr5q/chr6q)",
                                  "C13 (chr8q/chr11q)","C15 (PAR1)","cross/unknown"),
                       name = "Leiden community\n(query arm)") +
    facet_wrap(~ label_f, ncol = 1, strip.position = "right") +
    geom_text(data = within_totals, aes(x = 1, y = 480, label = label_pct),
              hjust = 0, vjust = 1, size = 2.2, color = "black", inherit.aes = FALSE) +
    scale_x_discrete(drop = FALSE,
                     labels = function(x) sub("chr", "", x)) +
    scale_y_continuous("Position in terminal flank (kb)",
                       breaks = c(0, 250, 500), limits = c(0, 510), expand = c(0, 0)) +
    labs(title = "a  WashU 3-generation T2T pedigree",
         subtitle = "Inter-chromosomal patches in odgi untangle (HQ filter); 92% (494/538) sit in a Leiden community",
         x = "Query chromosome arm (terminal 500 kb)") +
    theme_fig(8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          legend.position = "right",
          legend.box.margin = margin(0, 0, 0, 0))
}

# ---------- Panel 4b: CEPH1463 cross-assembler 11-feature matrix ----------
make_panel_4b <- function() {
  # Hard-coded from SURVEY_14 §1.6 (hifiasm + verkko intersection, within-community)
  feats <- data.table(
    parent  = c("NA12877","NA12877","NA12877","NA12877","NA12878","NA12878","NA12878","NA12878","NA12889","NA12890","NA12892"),
    chr_pair= c("chr1/chr19","chr10/chr18","chr17/chr19","chr6/chr9","chr10/chr18","chr19/chr22","chr21/chr22","chr6/chr9","chr12/chr9","chr12/chr9","chr21/chr22"),
    community = c("C6","C2","C6","C5","C2","C6","C7","C5","C5","C5","C6"),
    hifi_best = c(0.982,0.998,0.984,0.975,0.972,0.897,0.992,0.947,0.957,0.971,0.994),
    verk_best = c(0.860,0.997,0.868,0.989,0.998,0.978,0.996,0.959,0.978,0.976,0.993),
    hifi_children = c(2,2,1,1,3,3,1,4,1,1,1),
    verk_children = c(1,2,1,1,2,4,2,1,1,1,1)
  )
  feats[, best := pmax(hifi_best, verk_best)]
  feats[, children_total := hifi_children + verk_children]
  feats[, label := sprintf("%.2f / %.2f\nh=%d v=%d", hifi_best, verk_best, hifi_children, verk_children)]
  feats[, comm_label := paste0(community, " · ", chr_pair)]

  parent_order <- c("NA12889","NA12890","NA12892","NA12877","NA12878")
  feats[, parent := factor(parent, levels = parent_order)]
  pair_order <- feats[order(community, chr_pair), unique(comm_label)]
  feats[, comm_label := factor(comm_label, levels = pair_order)]

  ggplot(feats, aes(x = comm_label, y = parent, fill = best)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label, color = best > 0.93),
              size = 1.85, lineheight = 0.8) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradientn(colours = brewer.pal(7, "YlGnBu"),
                         limits = c(0.85, 1.00),
                         oob = scales::squish, name = "best score\n(max hifi/verk)") +
    scale_y_discrete(limits = rev) +
    labs(title = "b  CEPH1463 cross-assembler-validated parent features (11)",
         subtitle = "Detected by hifiasm AND verkko in ≥1 child each, same Leiden community.\nchr10/chr18 (C2, Linardopoulou) detected in NA12877 paternal and NA12878 maternal independently.",
         x = "Leiden community · chromosome pair", y = "Parent (G1 grandparents below; G2 above)") +
    theme_fig(8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 7, family = "mono"),
          legend.position = "right",
          legend.key.height = unit(0.7, "cm"),
          legend.key.width = unit(0.3, "cm"))
}

# ---------- Panel 4c: RPE-1 t(X;10) sequence + 3D ----------
make_panel_4c <- function() {
  # Sequence: arm-level Jaccard distance matrix (rpe1.dist_matrix.tsv)
  d <- as.matrix(fread(RPE_DIST), rownames = 1)
  d_arms <- rownames(d)
  # Order arms canonically chr1..22, X with p then q
  order_arms <- intersect(ARM_LEVELS, d_arms)
  d <- d[order_arms, order_arms]
  # Convert distance to similarity (1 - d) for plotting
  s <- 1 - d
  diag(s) <- NA
  long_s <- melt_matrix(s, varnames = c("a", "b"), value.name = "sim")
  long_s[, a := factor(a, levels = order_arms)]
  long_s[, b := factor(b, levels = order_arms)]

  # 3D: per-haplotype contact matrix; collapse HAP1+HAP2 by max for arm-pair display
  cmat <- as.matrix(fread(RPE_CONT), rownames = 1)
  # Strip _HAP1/_HAP2 to get arm names
  raw <- rownames(cmat)
  arm_of <- sub("_HAP[12]_", "_", raw)
  arm_of <- sub("_HAP[12]$", "", arm_of)
  arm_of <- sub("_(p|q)$", "_\\1", arm_of)
  # Ensure consistent format chr#_p / chr#_q
  arm_of <- gsub("_HAP[12]", "", raw)
  # max-collapse to arm level
  unique_arms <- unique(arm_of)
  collapsed <- matrix(0, length(unique_arms), length(unique_arms),
                      dimnames = list(unique_arms, unique_arms))
  for (i in seq_along(raw)) {
    for (j in seq_along(raw)) {
      a <- arm_of[i]; b <- arm_of[j]
      v <- cmat[i, j]
      if (!is.na(v) && v > collapsed[a, b]) collapsed[a, b] <- v
    }
  }
  use_arms <- intersect(order_arms, unique_arms)
  c2 <- collapsed[use_arms, use_arms]
  diag(c2) <- NA
  long_c <- melt_matrix(c2, varnames = c("a", "b"), value.name = "contact")
  long_c[, a := factor(a, levels = use_arms)]
  long_c[, b := factor(b, levels = use_arms)]

  # Highlight box: chrX_q × chr10_q
  hl <- data.table(a = "chrX_q", b = "chr10_q")

  p1 <- ggplot(long_s, aes(x = a, y = b, fill = sim)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", na.value = "grey95",
                        limits = c(0, 0.5), oob = scales::squish, name = "Jaccard") +
    geom_tile(data = hl, fill = NA, color = "red", linewidth = 0.7,
              inherit.aes = FALSE, aes(x = a, y = b)) +
    geom_tile(data = hl[, .(a = b, b = a)], fill = NA, color = "red", linewidth = 0.7,
              inherit.aes = FALSE, aes(x = a, y = b)) +
    scale_x_discrete(labels = function(x) sub("chr", "", x)) +
    scale_y_discrete(limits = rev, labels = function(x) sub("chr", "", x)) +
    labs(subtitle = "Sequence (PHR Jaccard)",
         x = NULL, y = NULL) +
    coord_fixed() +
    theme_fig(7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
          axis.text.y = element_text(size = 4.5),
          legend.position = "right",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"))

  p2 <- ggplot(long_c, aes(x = a, y = b, fill = log10(contact + 1e-6))) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", na.value = "grey95",
                        name = "log10\ncontact") +
    geom_tile(data = hl, fill = NA, color = "red", linewidth = 0.7,
              inherit.aes = FALSE, aes(x = a, y = b)) +
    geom_tile(data = hl[, .(a = b, b = a)], fill = NA, color = "red", linewidth = 0.7,
              inherit.aes = FALSE, aes(x = a, y = b)) +
    scale_x_discrete(labels = function(x) sub("chr", "", x)) +
    scale_y_discrete(limits = rev, labels = function(x) sub("chr", "", x)) +
    labs(subtitle = "3D (RPE-1 async CiFi, 50 kb)",
         x = NULL, y = NULL) +
    coord_fixed() +
    theme_fig(7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
          axis.text.y = element_text(size = 4.5),
          legend.position = "right",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.25, "cm"))

  list(left = p1, right = p2)
}

# ---------- Panel 4d: Mouse zygotene PHR-pair scatter ----------
make_panel_4d <- function() {
  d <- fread(MOUSE_TSV)
  # Inter-chromosomal only (different chromosomes)
  d[, chr_a := sub("_[pq]$", "", arm_a)]
  d[, chr_b := sub("_[pq]$", "", arm_b)]
  d_inter <- d[chr_a != chr_b & !is.na(mean_jaccard) & !is.na(hic_contact)]
  d_inter <- d_inter[hic_contact > 0]
  rho <- cor(d_inter$mean_jaccard, d_inter$hic_contact, method = "spearman")
  # Compute p (Spearman)
  ct <- suppressWarnings(cor.test(d_inter$mean_jaccard, d_inter$hic_contact, method = "spearman"))
  p_val <- ct$p.value

  ggplot(d_inter, aes(x = mean_jaccard, y = hic_contact)) +
    geom_point(alpha = 0.45, size = 0.9, color = "#4575b4") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6,
                formula = y ~ x) +
    scale_x_continuous("Mean Jaccard similarity (PHR pair)",
                       limits = c(0, max(d_inter$mean_jaccard, na.rm = TRUE) * 1.02),
                       expand = c(0, 0)) +
    scale_y_log10("Hi-C contact (zygotene, 50 kb)",
                  labels = function(x) format(x, scientific = TRUE, digits = 1)) +
    annotate("text", x = 0.02, y = max(d_inter$hic_contact) * 0.9,
             hjust = 0, vjust = 1, size = 3,
             label = sprintf("ρ = %.3f\np = %.1e\nn = %d pairs",
                             rho, p_val, nrow(d_inter))) +
    labs(title = "d  Mouse zygotene cross-species generalisation",
         subtitle = "B6 + CAST T2T; per-PHR-pair similarity vs Hi-C contact (Zuo et al. 2021)") +
    theme_fig(8) +
    theme(legend.position = "none")
}

# ---------- Render ----------
message("Building panel 4a...")
p4a <- make_panel_4a()
ggsave(file.path(OUT_DIR, "_panel_4a.pdf"), p4a, width = 7.5, height = 5.0, device = cairo_pdf)

message("Building panel 4b...")
p4b <- make_panel_4b()
ggsave(file.path(OUT_DIR, "_panel_4b.pdf"), p4b, width = 6.5, height = 5.0, device = cairo_pdf)

message("Building panel 4c...")
p4c <- make_panel_4c()

# Compose 4c left+right inside one PDF page using grid.layout
cairo_pdf(file.path(OUT_DIR, "_panel_4c.pdf"), width = 7.5, height = 4.6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                            heights = unit(c(0.16, 0.84), "npc"),
                                            widths  = unit(c(0.5, 0.5), "npc"))))
# Title row spans both columns
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text("c  RPE-1 t(X;10) rediscovery: sequence → 3D",
          x = 0.02, y = 0.75, just = c("left","center"),
          gp = gpar(fontface = "bold", fontsize = 10))
grid.text("Self-discovered Leiden community C2 = {chr10_q, chrX_q}; chrX_HAP1 carries chr10q material",
          x = 0.02, y = 0.30, just = c("left","center"),
          gp = gpar(fontsize = 8, col = "grey30"))
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); print(p4c$left, newpage = FALSE); upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); print(p4c$right, newpage = FALSE); upViewport()
upViewport()
dev.off()

message("Building panel 4d...")
p4d <- make_panel_4d()
ggsave(file.path(OUT_DIR, "_panel_4d.pdf"), p4d, width = 6.5, height = 5.0, device = cairo_pdf)

# ---------- Compose 2x2 combined figure ----------
message("Composing combined figure...")

# Build 4c as a single grob (left + right share a row, with title above)
make_panel_4c_grob <- function(p4c) {
  # Use grobTree with a layout
  g_left  <- ggplotGrob(p4c$left)
  g_right <- ggplotGrob(p4c$right)
  # Title text
  title <- textGrob("c  RPE-1 t(X;10) rediscovery: sequence → 3D",
                    x = 0.02, y = 0.85, just = c("left", "center"),
                    gp = gpar(fontface = "bold", fontsize = 10))
  subtitle <- textGrob(
    "Self-discovered Leiden community C2 = {chr10_q, chrX_q}; chrX_HAP1 carries chr10q material",
    x = 0.02, y = 0.45, just = c("left", "center"),
    gp = gpar(fontsize = 8, col = "grey30"))
  body <- gTree(children = gList(
    rectGrob(gp = gpar(fill = NA, col = NA))
  ))
  # We'll just composite manually inside the combine pdf below.
  list(title = title, subtitle = subtitle, left = g_left, right = g_right)
}

g4c <- make_panel_4c_grob(p4c)
g4a <- ggplotGrob(p4a)
g4b <- ggplotGrob(p4b)
g4d <- ggplotGrob(p4d)

# Final figure size: 14 x 10 inches (landscape)
FIG_W <- 14
FIG_H <- 10

cairo_pdf(file.path(OUT_DIR, "figure_fig4.pdf"), width = FIG_W, height = FIG_H)
grid.newpage()
# Outer 2 rows x 2 cols, each row 50% height
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(0.5, 0.5), "npc"),
                                           widths  = unit(c(0.55, 0.45), "npc"))))

# Top-left: panel 4a
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(g4a)
upViewport()

# Top-right: panel 4b
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(g4b)
upViewport()

# Bottom-left: panel 4c (compound)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(layout = grid.layout(3, 2,
                                            heights = unit(c(0.10, 0.06, 0.84), "npc"),
                                            widths = unit(c(0.5, 0.5), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.draw(g4c$title); upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
grid.draw(g4c$subtitle); upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1)); grid.draw(g4c$left); upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2)); grid.draw(g4c$right); upViewport()
upViewport()
upViewport()

# Bottom-right: panel 4d
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(g4d)
upViewport()

upViewport()
dev.off()

# Render PNG of the combined figure (300 dpi-ish via cairo_ps then convert)
# We'll let an external tool produce PNG; here we also emit a 200dpi PNG via R cairo_ps
# Use png() with cairo backend
png(file.path(OUT_DIR, "figure_fig4.png"), width = FIG_W * 200, height = FIG_H * 200,
    res = 200, type = "cairo")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
                                           heights = unit(c(0.5, 0.5), "npc"),
                                           widths  = unit(c(0.55, 0.45), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); grid.draw(g4a); upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); grid.draw(g4b); upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
pushViewport(viewport(layout = grid.layout(3, 2,
                                            heights = unit(c(0.10, 0.06, 0.84), "npc"),
                                            widths = unit(c(0.5, 0.5), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2)); grid.draw(g4c$title); upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.draw(g4c$subtitle); upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1)); grid.draw(g4c$left); upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2)); grid.draw(g4c$right); upViewport()
upViewport()
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); grid.draw(g4d); upViewport()
upViewport()
dev.off()

message("Wrote: ", file.path(OUT_DIR, "figure_fig4.pdf"))
message("Wrote: ", file.path(OUT_DIR, "figure_fig4.png"))

# Clean up per-panel intermediate PDFs
for (f in c("_panel_4a.pdf", "_panel_4b.pdf", "_panel_4c.pdf", "_panel_4d.pdf")) {
  fp <- file.path(OUT_DIR, f)
  if (file.exists(fp)) file.remove(fp)
}
