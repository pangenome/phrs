#!/usr/bin/env Rscript
# Extended Data Figure 2 — Sequence-level (50-community) detail
# 4 panels in a 2x2 grid:
#   ED2a UMAP coloured by 50-community partition
#   ED2b Within-community Jaccard distance bimodality (8 arm-level communities)
#   ED2c Cross-arm affinity radial plot — 41 arms
#   ED2d Confusion matrix Arm-Leiden (15) vs Sequence-Leiden (50)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
  library(scales)
  library(RColorBrewer)
})

OUT_DIR <- "paper_prep/figures/ed2"

# ---- inputs ----
UMAP_RDS  <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.umap.rds"
SEQ_ASSIGN <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv"
ARM_ASSIGN <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
CROSS_ARM  <- "/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_affinity_sequences.tsv"

seq_a <- fread(SEQ_ASSIGN, sep = "\t", header = TRUE)
arm_a <- fread(ARM_ASSIGN, sep = "\t", header = TRUE)

# ---------- ED2a: UMAP coloured by 50-community ----------
um <- readRDS(UMAP_RDS)
um_dt <- data.table(Name = rownames(um), UMAP1 = um[, 1], UMAP2 = um[, 2])
um_dt <- merge(um_dt, seq_a[, .(Name, Community, ChromArm)], by = "Name")
# Order communities by size for legend ordering
order_lev <- seq_a[, .N, by = Community][order(-N)][, Community]
um_dt[, Community := factor(Community, levels = order_lev)]

# 50-community palette — interpolated from a perceptually broad set
make_palette <- function(n) {
  base <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"),
            brewer.pal(8, "Set2"), brewer.pal(12, "Set3"),
            brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))
  base <- unique(base)
  if (length(base) >= n) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}
pal50 <- make_palette(length(order_lev))
names(pal50) <- order_lev

# Top-8 mixed communities to label on the plot
top_label_communities <- seq_a[, .N, by = Community][order(-N)][1:8, Community]
centroids <- um_dt[Community %in% top_label_communities,
                   .(x = median(UMAP1), y = median(UMAP2), n = .N), by = Community]

p_a <- ggplot(um_dt, aes(UMAP1, UMAP2, colour = Community)) +
  geom_point(size = 0.35, alpha = 0.55, stroke = 0) +
  geom_label(data = centroids,
             aes(x, y, label = sprintf("%s n=%d", Community, n)),
             inherit.aes = FALSE, size = 2.4, fontface = "bold",
             colour = "grey10", fill = alpha("white", 0.85),
             label.padding = unit(0.10, "lines"), label.size = 0.15) +
  scale_colour_manual(values = pal50, guide = "none") +
  labs(title = "ED2a — UMAP of 15,668 PHRs by 50-community partition",
       subtitle = "Leiden k = 75, resolution = 0.8, modularity 0.97; top 8 communities labelled",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"))

# ---------- ED2b: within-community Jaccard distance bimodality ----------
target_comms <- c("C1","C2","C3","C5","C6","C7","C11","C12")
within <- list()
for (c in target_comms) {
  fp <- file.path(OUT_DIR, sprintf("within_community_jaccard_%s.tsv", c))
  if (file.exists(fp)) {
    x <- fread(fp, sep = "\t", header = TRUE)
    x[, community := c]
    within[[c]] <- x
  }
}
within_dt <- rbindlist(within)
# Convert similarity to distance
within_dt[, distance := 1 - jaccard]
# Subsample for plotting if too dense (keep <= 50000 per community)
set.seed(42)
within_dt <- within_dt[, .SD[sample(.N, min(.N, 50000))], by = community]

# Community member arms (for facet annotation)
comm_arms <- arm_a[, .(arms = paste(Arms, collapse = "; ")), by = Community]
setnames(comm_arms, "Community", "community")
within_dt <- merge(within_dt, comm_arms, by = "community", all.x = TRUE)
within_dt[, facet_label := sprintf("%s (%s)", community, arms)]
# Order facets by community number
within_dt[, comm_n := as.integer(sub("C", "", community))]
within_dt <- within_dt[order(comm_n)]
within_dt[, facet_label := factor(facet_label, levels = unique(facet_label))]

p_b <- ggplot(within_dt, aes(distance)) +
  geom_histogram(binwidth = 0.02, fill = "#3B7BB8", colour = "white", linewidth = 0.05) +
  facet_wrap(~ facet_label, ncol = 4, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0, 0))) +
  labs(title = "ED2b — within-community Jaccard distance distributions",
       subtitle = sprintf("8 arm-level communities (≤50,000 random pairs each); bimodality marks allele/paralog separation"),
       x = "Jaccard distance (1 − jaccard.similarity)", y = "pair count") +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        strip.text = element_text(size = 7, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.6, "lines"))

# ---------- ED2c: cross-arm affinity radial plot ----------
ca <- fread(CROSS_ARM, sep = "\t", header = TRUE)
# Restrict to true cross-arm (affinity > 1) and aggregate
ca_xa <- ca[cross_arm_affinity > 1]
edges <- ca_xa[, .N, by = .(own_arm, affinity_arm)]
setnames(edges, "N", "weight")
edges[, own_arm := sub("^chr", "", sub("arm$", "", own_arm))]
edges[, affinity_arm := sub("^chr", "", sub("arm$", "", affinity_arm))]

# Use arm-leiden universe (41 arms in the partition) — arm assignments file lists them.
arm_universe <- arm_a$ChromArm
arm_universe <- sub("^chr", "", arm_universe)
chrom_idx2 <- function(s) {
  c <- sub("_[pq]$", "", s)
  c <- ifelse(c == "X", 23L, ifelse(c == "Y", 24L, suppressWarnings(as.integer(c))))
  c
}
arm_idx2 <- function(s) ifelse(grepl("_p$", s), 0L, 1L)
arm_order_dt <- data.table(arm = arm_universe,
                           ord = chrom_idx2(arm_universe) * 2 + arm_idx2(arm_universe))
arm_order_dt <- arm_order_dt[order(ord)]
arm_levels <- arm_order_dt$arm
n_arms <- length(arm_levels)

# Layout: place arms around unit circle (rotate so chr1_p starts at 12 o'clock)
theta <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n_arms + 1)[-1]
# Clockwise: negate
theta <- pi / 2 - (theta - pi / 2)
arm_pos <- data.table(arm = arm_levels, theta = theta,
                      x = cos(theta), y = sin(theta))

# Per-arm community colour (from arm-leiden) — arms in the same community get
# the same hue so the visualisation reflects the partition.
arm_to_comm <- arm_a[, .(arm = sub("^chr", "", ChromArm), Community)]
arm_to_comm[, comm_n := as.integer(sub("C", "", Community))]
arm_pos <- merge(arm_pos, arm_to_comm[, .(arm, Community, comm_n)], by = "arm",
                 sort = FALSE, all.x = TRUE)
n_comms <- max(arm_to_comm$comm_n, na.rm = TRUE)
comm_pal <- colorRampPalette(brewer.pal(11, "Spectral"))(n_comms)

draw_ed2c <- function() {
  pushViewport(viewport(layout = grid.layout(1, 1)))
  pushViewport(viewport(x = 0.5, y = 0.5, width = 1, height = 1))

  grid.text("ED2c — cross-arm affinity (sequences absorbed into a foreign arm community)",
            x = 0.5, y = 0.97, gp = gpar(fontsize = 10.5, fontface = "bold"))
  grid.text(sprintf("%d arms (origin → affinity); chord width ∝ √(#sequences); colour = arm-Leiden community",
                    n_arms),
            x = 0.5, y = 0.93, gp = gpar(fontsize = 8, col = "grey25"))

  cx_ <- 0.50; cy_ <- 0.46; r_outer <- 0.34; r_inner <- 0.31
  xy <- function(th, rad) c(cx_ + rad * cos(th), cy_ + rad * sin(th))

  weight_max <- max(edges$weight)
  weight_to_alpha <- function(w) 0.20 + 0.80 * (sqrt(w) / sqrt(weight_max))
  weight_to_lwd   <- function(w) 0.4 + 3.0 * (sqrt(w) / sqrt(weight_max))

  edges_ord <- edges[order(weight)]
  for (i in seq_len(nrow(edges_ord))) {
    a <- arm_pos[arm == edges_ord$own_arm[i]]
    b <- arm_pos[arm == edges_ord$affinity_arm[i]]
    if (nrow(a) == 0 || nrow(b) == 0) next
    t <- seq(0, 1, length.out = 40)
    bx <- (1 - t) ^ 2 * (cx_ + r_inner * cos(a$theta)) +
          2 * (1 - t) * t * cx_ +
          t ^ 2 * (cx_ + r_inner * cos(b$theta))
    by <- (1 - t) ^ 2 * (cy_ + r_inner * sin(a$theta)) +
          2 * (1 - t) * t * cy_ +
          t ^ 2 * (cy_ + r_inner * sin(b$theta))
    col <- if (!is.na(b$comm_n)) comm_pal[b$comm_n] else "grey60"
    grid.lines(x = unit(bx, "npc"), y = unit(by, "npc"),
               gp = gpar(col = col, lwd = weight_to_lwd(edges_ord$weight[i]),
                         alpha = weight_to_alpha(edges_ord$weight[i])))
  }

  # Arm community-coloured arc segments + labels
  for (i in seq_len(nrow(arm_pos))) {
    a <- arm_pos[i]
    p1 <- xy(a$theta, r_outer);          p2 <- xy(a$theta, r_outer + 0.010)
    grid.lines(x = unit(c(p1[1], p2[1]), "npc"), y = unit(c(p1[2], p2[2]), "npc"),
               gp = gpar(col = "grey25", lwd = 0.4))
    # community-coloured dot
    p_col <- xy(a$theta, r_outer + 0.020)
    col_dot <- if (!is.na(a$comm_n)) comm_pal[a$comm_n] else "grey60"
    grid.circle(x = p_col[1], y = p_col[2], r = 0.008,
                gp = gpar(fill = col_dot, col = "grey20", lwd = 0.3))
    # label
    p3 <- xy(a$theta, r_outer + 0.034)
    rot <- (a$theta * 180 / pi) %% 360
    if (rot > 90 && rot < 270) {
      rot_label <- rot + 180; just <- c("right", "centre")
    } else {
      rot_label <- rot;       just <- c("left", "centre")
    }
    grid.text(a$arm, x = unit(p3[1], "npc"), y = unit(p3[2], "npc"),
              just = just, rot = rot_label, gp = gpar(fontsize = 6.5))
  }

  # Top edges legend (lower-left corner)
  top_edges <- edges_ord[order(-weight)][1:6]
  grid.text("Top 6 edges (origin → affinity, n)",
            x = 0.04, y = 0.16, just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold"))
  for (i in seq_len(nrow(top_edges))) {
    grid.text(sprintf("• %s → %s  (%d)",
                      top_edges$own_arm[i], top_edges$affinity_arm[i],
                      top_edges$weight[i]),
              x = 0.04, y = 0.14 - 0.022 * i, just = c("left", "top"),
              gp = gpar(fontsize = 7, col = "grey20"))
  }

  popViewport(2)
}

# ---------- ED2d: confusion matrix Arm-Leiden vs Seq-Leiden ----------
# Use seq-leiden assignments + arm assignments
arm_lookup <- arm_a[, .(ChromArm, ArmComm = Community)]
seq_a2 <- merge(seq_a[, .(Name, SeqComm = Community, ChromArm)],
                arm_lookup, by = "ChromArm", all.x = TRUE)
seq_a2 <- seq_a2[!is.na(ArmComm)]
conf <- seq_a2[, .N, by = .(ArmComm, SeqComm)]
# Order arm communities C1..C15 / seq C1..C50 by integer suffix
conf[, ArmComm_n := as.integer(sub("C", "", as.character(ArmComm)))]
conf[, SeqComm_n := as.integer(sub("C", "", as.character(SeqComm)))]
arm_levels_in <- paste0("C", sort(unique(conf$ArmComm_n)))
seq_levels_in <- paste0("C", sort(unique(conf$SeqComm_n)))
# Y axis: rows = arm communities, with C1 at the top (last factor level placed at top)
conf[, ArmComm := factor(as.character(ArmComm), levels = rev(arm_levels_in))]
conf[, SeqComm := factor(as.character(SeqComm), levels = seq_levels_in)]
# Compute ARI / NMI inline (small implementations)
ari_calc <- function(a, b) {
  ct <- table(a, b)
  n <- sum(ct)
  ai <- rowSums(ct); bj <- colSums(ct)
  index   <- sum(choose(ct, 2))
  expect  <- sum(choose(ai, 2)) * sum(choose(bj, 2)) / choose(n, 2)
  maxidx  <- 0.5 * (sum(choose(ai, 2)) + sum(choose(bj, 2)))
  (index - expect) / (maxidx - expect)
}
nmi_calc <- function(a, b) {
  ct <- table(a, b); n <- sum(ct)
  pij <- ct / n
  pi  <- rowSums(pij); pj <- colSums(pij)
  Hi  <- -sum(pi[pi > 0] * log(pi[pi > 0]))
  Hj  <- -sum(pj[pj > 0] * log(pj[pj > 0]))
  Iij <- 0
  for (r in seq_len(nrow(pij))) for (c in seq_len(ncol(pij))) {
    if (pij[r, c] > 0) Iij <- Iij + pij[r, c] * log(pij[r, c] / (pi[r] * pj[c]))
  }
  2 * Iij / (Hi + Hj)
}
ari <- ari_calc(seq_a2$ArmComm, seq_a2$SeqComm)
nmi <- nmi_calc(seq_a2$ArmComm, seq_a2$SeqComm)

p_d <- ggplot(conf, aes(SeqComm, ArmComm, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.05) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFEDA0", "#FEB24C", "#F03B20", "#7F0000"),
    trans = "sqrt", name = "# PHRs",
    labels = function(x) format(round(x), big.mark = ",")) +
  labs(title = sprintf("ED2d — Arm-Leiden (15) vs Sequence-Leiden (50) confusion (ARI %.2f, NMI %.2f)",
                       ari, nmi),
       subtitle = sprintf("n = %s PHRs; rows = arm-level communities, columns = sequence-level communities",
                          comma(nrow(seq_a2))),
       x = "Sequence-Leiden community", y = "Arm-Leiden community") +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        legend.key.height = unit(0.8, "lines"),
        legend.key.width = unit(0.4, "lines"))

# ---------- compose 2x2 ----------
pdf_path <- file.path(OUT_DIR, "figure_ed2.pdf")
png_path <- file.path(OUT_DIR, "figure_ed2.png")
W <- 14; H <- 11

cairo_pdf(pdf_path, width = W, height = H)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p_a, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); print(p_b, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); draw_ed2c(); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); print(p_d, newpage = FALSE); popViewport()
popViewport()
dev.off()
cat("wrote", pdf_path, "\n")

png(png_path, width = W, height = H, units = "in", res = 200)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p_a, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); print(p_b, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); draw_ed2c(); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); print(p_d, newpage = FALSE); popViewport()
popViewport()
dev.off()
cat("wrote", png_path, "\n")

# ---- metrics ----
cat("\n--- KEY METRICS ---\n")
cat(sprintf("ED2a: PHRs plotted = %d, communities = %d\n",
            nrow(um_dt), length(unique(um_dt$Community))))
cat(sprintf("ED2b: target communities = %s, total within-community pairs (subsampled): %d\n",
            paste(target_comms, collapse = ","), nrow(within_dt)))
cat(sprintf("ED2c: edges = %d, total absorbed sequences = %d, top edge = %s\n",
            nrow(edges), sum(edges$weight),
            paste0(edges[order(-weight)][1, .(own_arm, affinity_arm, weight)], collapse = "/")))
cat(sprintf("ED2d: rows (arm) = %d, cols (seq) = %d, ARI = %.4f, NMI = %.4f\n",
            length(unique(conf$ArmComm)), length(unique(conf$SeqComm)), ari, nmi))
