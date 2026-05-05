#!/usr/bin/env Rscript
# Extended Data Figure 1 — Pipeline and per-arm flank inventory
# 4 panels in a 2x2 grid:
#   ED1a Pipeline schematic (465 → 18,827 → 15,668 → 15/50 communities)
#   ED1b Per-arm flank counts (48 arms) with assembly QC overlay
#   ED1c PHR length distribution (median 105 kb, mean 144 kb)
#   ED1d Chr18_q (NA18982#1) chimera evidence schematic
#
# Outputs:
#   paper_prep/figures/ed1/figure_ed1.pdf
#   paper_prep/figures/ed1/figure_ed1.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
  library(scales)
  library(RColorBrewer)
})

OUT_DIR <- "paper_prep/figures/ed1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Inputs
LEN_TSV <- "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
CONTIG_TSV <- "/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/contig_classifications.tsv"

# ----- Load data -----
len <- fread(LEN_TSV, sep = "\t", header = TRUE, na.strings = c(".", "NA"))
# Parse arm from seq name suffix _chrN_[pq]arm
arm_match <- regmatches(len$seq, regexpr("_chr[0-9XY]+_[pq]arm$", len$seq))
arm_chr <- sub("^_(chr[0-9XY]+)_[pq]arm$", "\\1", arm_match)
arm_pq  <- sub("^_chr[0-9XY]+_([pq])arm$", "\\1", arm_match)
len[, arm_label := paste0(arm_chr, "_", arm_pq)]

# PHR length: end - start, only for signal-bearing rows (non-NA region_start)
len[, phr_length := as.numeric(region_end) - as.numeric(region_start)]

# ----- ED1b — per-arm flank counts -----
flank_counts <- len[, .N, by = arm_label]
setnames(flank_counts, "N", "n_flanks")

contigs <- fread(CONTIG_TSV, sep = "\t", header = TRUE)
# Per-arm counts of flagged (non-pass) and pass-with-mixed-strand-caveats classifications.
# Classify per arm using p_arm_chrs / q_arm_chrs columns. Build a per-(chr,arm) mixed-strand-caveat count.
contigs[, mixed_caveat := validation_status == "pass" & strand_confidence < 1.0]
# Aggregate per chr — assignment of arm comes from non-empty p_arm_chrs / q_arm_chrs.
flag_p <- contigs[p_arm_chrs != "" & p_arm_chrs != ".",
                  .(n_flag = sum(validation_status != "pass"),
                    n_caveat = sum(mixed_caveat),
                    n_total  = .N),
                  by = .(chr_arm = paste0(sub("^[^#]+#[0-9]+#", "", p_arm_chrs), "_p"))]
flag_q <- contigs[q_arm_chrs != "" & q_arm_chrs != ".",
                  .(n_flag = sum(validation_status != "pass"),
                    n_caveat = sum(mixed_caveat),
                    n_total  = .N),
                  by = .(chr_arm = paste0(sub("^[^#]+#[0-9]+#", "", q_arm_chrs), "_q"))]
flag_all <- rbind(flag_p, flag_q)[, .(n_flag = sum(n_flag), n_caveat = sum(n_caveat), n_total = sum(n_total)), by = chr_arm]
setnames(flag_all, "chr_arm", "arm_label")
flank_counts <- merge(flank_counts, flag_all, by = "arm_label", all.x = TRUE)
flank_counts[is.na(n_flag),   n_flag   := 0]
flank_counts[is.na(n_caveat), n_caveat := 0]
flank_counts[is.na(n_total),  n_total  := 0]

# Order arms naturally: chr1_p, chr1_q, chr2_p, ..., chr22_q, chrX_p, chrX_q, chrY_p, chrY_q
chrom_idx <- function(s) {
  c <- sub("_[pq]$", "", s)
  c <- sub("^chr", "", c)
  c <- ifelse(c == "X", 23L, ifelse(c == "Y", 24L, suppressWarnings(as.integer(c))))
  c
}
arm_idx <- function(s) ifelse(grepl("_p$", s), 0L, 1L)
flank_counts[, ord := chrom_idx(arm_label) * 2 + arm_idx(arm_label)]
flank_counts <- flank_counts[order(ord)]
flank_counts[, arm_label := factor(arm_label, levels = arm_label)]
# acrocentric / sex flag (under-represented)
flank_counts[, group := ifelse(grepl("^chr(13|14|15|21|22)_p$", arm_label), "acrocentric p",
                        ifelse(grepl("^chr(X|Y)_", arm_label), "sex chromosome", "autosome"))]

p_b <- ggplot(flank_counts, aes(x = arm_label, y = n_flanks, fill = group)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = median(flank_counts$n_flanks), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  geom_text(data = flank_counts[n_flag > 0],
            aes(x = arm_label, y = n_flanks + 25, label = sprintf("%d*", n_flag)),
            inherit.aes = FALSE, size = 2.4, colour = "grey20") +
  scale_fill_manual(values = c("autosome" = "#3B7BB8",
                               "acrocentric p" = "#E07B39",
                               "sex chromosome" = "#7B5EA7"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = "ED1b — per-arm flank inventory (48 arms)",
       subtitle = sprintf("Total %s flanks (median %d/arm); * = mixed-strand-caveat contigs (validation_status != pass)",
                          comma(sum(flank_counts$n_flanks)), as.integer(median(flank_counts$n_flanks))),
       x = NULL, y = "flanks per arm") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = c(0.85, 0.94),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# ----- ED1c — PHR length distribution -----
len_signal <- len[!is.na(phr_length) & phr_length > 0]
len_med <- median(len_signal$phr_length)
len_mean <- mean(len_signal$phr_length)
n_phr <- nrow(len_signal)

p_c <- ggplot(len_signal, aes(x = phr_length / 1e3)) +
  geom_histogram(binwidth = 10, fill = "#3B7BB8", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = len_med / 1e3, colour = "#C0392B", linewidth = 0.5) +
  geom_vline(xintercept = len_mean / 1e3, colour = "#27AE60", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = len_med / 1e3 + 5, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("median %.0f kb", len_med / 1e3), colour = "#C0392B", size = 2.6) +
  annotate("text", x = len_mean / 1e3 + 5, y = Inf, vjust = 3.0, hjust = 0,
           label = sprintf("mean %.0f kb", len_mean / 1e3), colour = "#27AE60", size = 2.6) +
  scale_x_continuous(limits = c(0, 510), expand = expansion(mult = c(0, 0))) +
  labs(title = "ED1c — PHR length distribution",
       subtitle = sprintf("n = %s signal-bearing flanks (%.1f%% of 18,827)",
                          comma(n_phr), 100 * n_phr / 18827),
       x = "PHR length (kb)", y = "flanks") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        panel.grid.minor = element_blank())

# ----- ED1a — pipeline schematic (drawn with grid primitives) -----
# Five stages: 465 assemblies -> 18,827 flanks -> 18,827 PAFs -> 15,668 PHRs -> 15 / 50 communities
draw_ed1a <- function() {
  pushViewport(viewport(layout = grid.layout(1, 1)))
  vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1)
  pushViewport(vp)

  # Title
  grid.text("ED1a — pipeline overview", x = 0.5, y = 0.95,
            gp = gpar(fontsize = 11, fontface = "bold"))

  # 5 stages stacked vertically with clear labels
  stages <- list(
    list(label = "465 assemblies",          sub = "HPRCv2 v1.1 — 233 individuals",                 colour = "#3B7BB8"),
    list(label = "18,827 flanks",            sub = "500 kb telomere-anchored windows (48 arms)",   colour = "#3B7BB8"),
    list(label = "18,827 PAFs",              sub = "wfmash v0.23 -p95 -t48 (all-vs-all)",          colour = "#7B5EA7"),
    list(label = "15,668 PHRs",              sub = "id≥0.95, ≥2 chrs, ≥3 kb region", colour = "#7B5EA7"),
    list(label = "15 / 50 communities",      sub = "Leiden arm + seq (silhouette 0.347 / 0.602)",  colour = "#E07B39")
  )
  ann <- list(
    "classify_contigs.py",
    "extract_telo_flanks.sh",
    "wfmash + impg sliding-window scan",
    "pggb + odgi similarity (Jaccard)",
    "Leiden (arm: r=1.16; seq: k=75 r=0.8)"
  )

  n <- length(stages)
  ytop <- 0.83
  ybot <- 0.18
  bh <- 0.075
  bw <- 0.40
  cx_box <- 0.30
  cx_ann <- 0.74
  ys <- seq(ytop, ybot, length.out = n)
  for (i in seq_len(n)) {
    grid.roundrect(x = cx_box, y = ys[i], width = bw, height = bh,
                   r = unit(0.012, "npc"),
                   gp = gpar(fill = stages[[i]]$colour, col = "grey30", lwd = 0.6, alpha = 0.92))
    grid.text(stages[[i]]$label, x = cx_box, y = ys[i] + 0.013,
              gp = gpar(col = "white", fontsize = 10, fontface = "bold"))
    grid.text(stages[[i]]$sub, x = cx_box, y = ys[i] - 0.018,
              gp = gpar(col = "white", fontsize = 7))
    grid.text(ann[[i]], x = cx_ann, y = ys[i],
              just = c("left", "centre"),
              gp = gpar(fontsize = 8, col = "grey25", fontface = "italic"))
    if (i < n) {
      grid.lines(x = unit(c(cx_box, cx_box), "npc"),
                 y = unit(c(ys[i] - bh / 2 - 0.002, ys[i + 1] + bh / 2 + 0.002), "npc"),
                 arrow = arrow(angle = 18, length = unit(0.020, "npc"), type = "closed"),
                 gp = gpar(col = "grey30", lwd = 0.8, fill = "grey30"))
    }
  }

  # Bottom annotations
  grid.text("83.2% of flanks (15,668 / 18,827) carry an inter-chromosomal match;",
            x = 0.5, y = 0.10, gp = gpar(fontsize = 8.5, col = "#27AE60", fontface = "bold"))
  grid.text("41 of 48 arms enter the community partition (7 zero-signal arms removed)",
            x = 0.5, y = 0.06, gp = gpar(fontsize = 7.5, col = "grey25"))

  popViewport(2)
}

# ----- ED1d — chr18_q chimera schematic -----
draw_ed1d <- function() {
  pushViewport(viewport(layout = grid.layout(1, 1)))
  pushViewport(viewport(x = 0.5, y = 0.5, width = 1, height = 1))

  grid.text("ED1d — chr18_q chimera (NA18982#1, JBKABS010000018.1)",
            x = 0.5, y = 0.95, gp = gpar(fontsize = 11, fontface = "bold"))

  # Contig coordinate system 0 .. 84.4 Mb. Render as horizontal track normalised to [0.08, 0.92].
  ctg_len <- 84.4   # Mb
  x0 <- 0.08; x1 <- 0.92
  y_track <- 0.62
  to_x <- function(mb) x0 + (mb / ctg_len) * (x1 - x0)

  # Background contig bar
  grid.rect(x = (x0 + x1) / 2, y = y_track, width = x1 - x0, height = 0.05,
            gp = gpar(col = "grey25", fill = "grey90", lwd = 0.6))

  # Chr18 region 0 .. ~83.37 Mb
  grid.rect(x = (to_x(0) + to_x(83.37)) / 2, y = y_track,
            width = to_x(83.37) - to_x(0), height = 0.05,
            gp = gpar(col = NA, fill = "#3B7BB8", alpha = 0.9))
  # NNN gap 83.37 .. 83.38 Mb (100 bp scaffold join)
  grid.rect(x = (to_x(83.37) + to_x(83.38)) / 2, y = y_track,
            width = to_x(83.38) - to_x(83.37), height = 0.05,
            gp = gpar(col = NA, fill = "#C0392B"))
  # ChrX PAR1 region 83.38 .. ~84.346 Mb (966 kb)
  grid.rect(x = (to_x(83.38) + to_x(84.346)) / 2, y = y_track,
            width = to_x(84.346) - to_x(83.38), height = 0.05,
            gp = gpar(col = NA, fill = "#7B5EA7", alpha = 0.9))
  # Terminal TTAGGG tract 2,826 bp at end (84.346 .. 84.349)
  grid.rect(x = (to_x(84.346) + to_x(84.349)) / 2, y = y_track,
            width = to_x(84.349) - to_x(84.346), height = 0.05,
            gp = gpar(col = NA, fill = "#27AE60"))

  # Track legend strip at top — two rows for clarity
  legend_y1 <- 0.88
  legend_y2 <- 0.83
  legend_items <- list(
    list(col = "#3B7BB8", txt = "chr18 segment",                row = 1),
    list(col = "#C0392B", txt = "NNN scaffold gap (100 bp)",   row = 1),
    list(col = "#7B5EA7", txt = "chrX PAR1 (966 kb, mapq 60)", row = 2),
    list(col = "#27AE60", txt = "(TTAGGG)n × 471 (2,826 bp)",  row = 2)
  )
  for (row in 1:2) {
    items <- Filter(function(it) it$row == row, legend_items)
    lx <- 0.10
    ly <- if (row == 1) legend_y1 else legend_y2
    for (it in items) {
      grid.rect(x = lx, y = ly, width = 0.020, height = 0.024,
                gp = gpar(fill = it$col, col = NA))
      grid.text(it$txt, x = lx + 0.018, y = ly, just = c("left", "centre"),
                gp = gpar(fontsize = 7))
      lx <- lx + 0.020 + 0.018 + nchar(it$txt) * 0.0080
    }
  }

  # Coordinate ticks (overview)
  ticks_mb <- c(0, 20, 40, 60, 80)
  for (mb in ticks_mb) {
    xx <- to_x(mb)
    grid.lines(x = unit(c(xx, xx), "npc"), y = unit(c(y_track - 0.025, y_track - 0.04), "npc"),
               gp = gpar(col = "grey25", lwd = 0.4))
    grid.text(sprintf("%d Mb", mb), x = xx, y = y_track - 0.06,
              just = c("centre", "top"), gp = gpar(fontsize = 6.5, col = "grey25"))
  }
  # Junction markers (called out separately)
  for (mb in c(83.37, 84.35)) {
    xx <- to_x(mb)
    grid.lines(x = unit(c(xx, xx), "npc"), y = unit(c(y_track + 0.025, y_track + 0.06), "npc"),
               gp = gpar(col = "grey25", lwd = 0.4))
    grid.text(sprintf("%.2f Mb", mb), x = xx, y = y_track + 0.075,
              just = c("centre", "bottom"), gp = gpar(fontsize = 6.5, col = "grey25"))
  }

  # Junction zoom panel (lower)
  # Show 82.5 .. 84.4 Mb at higher zoom
  zy <- 0.30
  zx0 <- 0.08; zx1 <- 0.92
  zoom_min <- 82.5; zoom_max <- 84.4
  to_zx <- function(mb) zx0 + ((mb - zoom_min) / (zoom_max - zoom_min)) * (zx1 - zx0)
  grid.text("zoom: junction region (82.5 — 84.4 Mb)", x = 0.5, y = zy + 0.18,
            gp = gpar(fontsize = 9, fontface = "italic"))
  grid.rect(x = (zx0 + zx1) / 2, y = zy, width = zx1 - zx0, height = 0.05,
            gp = gpar(col = "grey25", fill = "grey90", lwd = 0.6))
  grid.rect(x = (to_zx(82.5) + to_zx(83.37)) / 2, y = zy,
            width = to_zx(83.37) - to_zx(82.5), height = 0.05,
            gp = gpar(col = NA, fill = "#3B7BB8", alpha = 0.9))
  grid.rect(x = (to_zx(83.37) + to_zx(83.38)) / 2, y = zy,
            width = to_zx(83.38) - to_zx(83.37), height = 0.05,
            gp = gpar(col = NA, fill = "#C0392B"))
  grid.rect(x = (to_zx(83.38) + to_zx(84.346)) / 2, y = zy,
            width = to_zx(84.346) - to_zx(83.38), height = 0.05,
            gp = gpar(col = NA, fill = "#7B5EA7", alpha = 0.9))
  grid.rect(x = (to_zx(84.346) + to_zx(84.349)) / 2, y = zy,
            width = to_zx(84.349) - to_zx(84.346), height = 0.05,
            gp = gpar(col = NA, fill = "#27AE60"))

  # Annotated arrows in zoom
  jx <- to_zx(83.375)
  grid.lines(x = unit(c(jx, jx), "npc"),
             y = unit(c(zy + 0.04, zy + 0.10), "npc"),
             arrow = arrow(angle = 25, length = unit(0.012, "npc"), type = "closed"),
             gp = gpar(col = "#C0392B", fill = "#C0392B", lwd = 0.6))
  grid.text("100 bp NNN scaffold join", x = jx, y = zy + 0.13,
            gp = gpar(fontsize = 7, col = "#C0392B"))

  tx <- to_zx(84.347)
  grid.lines(x = unit(c(tx, tx), "npc"),
             y = unit(c(zy - 0.04, zy - 0.10), "npc"),
             arrow = arrow(angle = 25, length = unit(0.012, "npc"), type = "closed"),
             gp = gpar(col = "#27AE60", fill = "#27AE60", lwd = 0.6))
  grid.text("471 TTAGGG repeats", x = tx, y = zy - 0.13,
            gp = gpar(fontsize = 7, col = "#27AE60"))

  # Coordinate ticks zoom
  for (mb in c(82.5, 83.0, 83.37, 83.5, 84.0, 84.35)) {
    xx <- to_zx(mb)
    grid.lines(x = unit(c(xx, xx), "npc"), y = unit(c(zy - 0.025, zy - 0.04), "npc"),
               gp = gpar(col = "grey25", lwd = 0.4))
    grid.text(sprintf("%.2f", mb), x = xx, y = zy - 0.06,
              just = c("centre", "top"), gp = gpar(fontsize = 6.5, col = "grey25"))
  }

  # Footnote
  grid.text(paste0("Cross-aligner agreement: wfmash v0.23.0 + minimap2 v2.30 ",
                   "both report mapq 60 chrX PAR1 alignment for the 966 kb terminal block."),
            x = 0.5, y = 0.07, gp = gpar(fontsize = 7, col = "grey25"))
  grid.text("Sequence dropped from PHR set (15,669 → 15,668).",
            x = 0.5, y = 0.03, gp = gpar(fontsize = 7, col = "grey25", fontface = "italic"))

  popViewport(2)
}

# ----- Compose 2x2 layout -----
pdf_path <- file.path(OUT_DIR, "figure_ed1.pdf")
png_path <- file.path(OUT_DIR, "figure_ed1.png")

cairo_pdf(pdf_path, width = 12, height = 9)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
# ED1a (top-left)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw_ed1a()
popViewport()
# ED1b (top-right)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_b, newpage = FALSE)
popViewport()
# ED1c (bottom-left)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_c, newpage = FALSE)
popViewport()
# ED1d (bottom-right)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw_ed1d()
popViewport()
popViewport()
dev.off()
cat("wrote", pdf_path, "\n")

png(png_path, width = 12, height = 9, units = "in", res = 200)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); draw_ed1a(); popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); print(p_b, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); print(p_c, newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); draw_ed1d(); popViewport()
popViewport()
dev.off()
cat("wrote", png_path, "\n")

# ----- Print key metrics for sources.tsv -----
cat("\n--- KEY METRICS ---\n")
cat(sprintf("ED1b: total flanks = %d, n arms = %d, median per arm = %.0f, range = [%d, %d]\n",
            sum(flank_counts$n_flanks), nrow(flank_counts),
            median(flank_counts$n_flanks),
            min(flank_counts$n_flanks), max(flank_counts$n_flanks)))
cat(sprintf("ED1c: signal-bearing PHRs = %d (%.2f%%), median = %.0f bp, mean = %.0f bp\n",
            n_phr, 100 * n_phr / nrow(len), len_med, len_mean))
cat("\n--- per-arm flank counts (ED1b) ---\n")
print(flank_counts)
