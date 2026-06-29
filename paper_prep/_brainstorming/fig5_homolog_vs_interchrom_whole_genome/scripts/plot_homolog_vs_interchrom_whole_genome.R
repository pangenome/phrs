#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) args[[1]] else "."
out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "fig5_homolog_vs_interchrom_whole_genome"
summary_dir <- file.path(repo_root, "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries")
tract_path <- file.path(summary_dir, "homolog_vs_interchrom_top_tracts.tsv")
pair_path <- file.path(summary_dir, "homolog_vs_interchrom_pair_summary.tsv")
overall_path <- file.path(summary_dir, "homolog_vs_interchrom_overall.tsv")
chrom_size_path <- file.path(repo_root, "data/chm13.chrom.sizes")
annotation_path <- file.path(repo_root, "data/chm13-annotations.bed")
candidate_path <- file.path(
  repo_root,
  "paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity/config/candidate_windows.tsv"
)

required_paths <- c(tract_path, pair_path, overall_path, chrom_size_path, annotation_path, candidate_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))
}

read_tsv <- function(path) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
basis_order <- c("1:1", "4:4", "10:10")

tracts <- read_tsv(tract_path)
pairs <- read_tsv(pair_path)
overall <- read_tsv(overall_path)
chroms <- read.delim(chrom_size_path, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
names(chroms) <- c("chrom", "length")
annotations <- read.delim(annotation_path, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
names(annotations) <- c("chrom", "start", "end", "label")
candidates <- read_tsv(candidate_path)

for (col in c("query_start", "query_end", "bp", "windows", "mean_delta", "max_delta", "mean_inter_identity", "mean_same_identity")) {
  if (col %in% names(tracts)) tracts[[col]] <- suppressWarnings(as.numeric(tracts[[col]]))
}
for (col in c("windows", "bp", "mean_delta", "median_delta", "max_delta", "mean_inter_identity", "mean_same_identity")) {
  if (col %in% names(pairs)) pairs[[col]] <- suppressWarnings(as.numeric(pairs[[col]]))
}
for (col in names(overall)[names(overall) != "basis"]) {
  overall[[col]] <- suppressWarnings(as.numeric(overall[[col]]))
}
chroms$length <- suppressWarnings(as.numeric(chroms$length))
annotations$start <- suppressWarnings(as.numeric(annotations$start))
annotations$end <- suppressWarnings(as.numeric(annotations$end))
candidates$query_start <- suppressWarnings(as.numeric(candidates$query_start))
candidates$query_end <- suppressWarnings(as.numeric(candidates$query_end))

chroms <- chroms[chroms$chrom %in% chrom_order, ]
chroms <- chroms[match(chrom_order, chroms$chrom), ]
chroms <- chroms[!is.na(chroms$chrom), ]
chroms$y <- rev(seq_len(nrow(chroms)))
chroms$arm_mid <- chroms$length / 2

len_map <- setNames(chroms$length, chroms$chrom)
y_map <- setNames(chroms$y, chroms$chrom)
x_max <- max(chroms$length, na.rm = TRUE)

annotations <- annotations[annotations$chrom %in% chrom_order, ]
context_rows <- rbind(
  annotations[annotations$label %in% c("PAR1", "PAR2", "PHR-acro", "Centromere"), ],
  data.frame(
    chrom = sub("^.*#(chr[^.]+).*$", "\\1", candidates$query_name),
    start = candidates$query_start,
    end = candidates$query_end,
    label = candidates$event_id,
    stringsAsFactors = FALSE
  )
)
context_rows <- context_rows[context_rows$chrom %in% chrom_order, ]
context_rows$context_class <- ifelse(
  context_rows$label %in% c("PAR1", "PAR2"), "PAR",
  ifelse(context_rows$label == "PHR-acro", "acrocentric_PHR",
    ifelse(context_rows$label == "Centromere", "centromere", "candidate_window")
  )
)
context_rows <- context_rows[order(match(context_rows$chrom, chrom_order), context_rows$start, context_rows$end), ]
write.table(context_rows, file.path(out_dir, "context_intervals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

target_chroms <- sort(unique(tracts$target_chrom[tracts$target_chrom != ""]))
target_chroms <- target_chroms[target_chroms %in% chrom_order]
base_cols <- grDevices::hcl(seq(12, 372, length.out = length(target_chroms) + 1)[-1], c = 72, l = 48)
target_cols <- setNames(base_cols, target_chroms)
fixed_cols <- c(
  chr1 = "#7A5195", chr2 = "#BC5090", chr3 = "#00796B", chr4 = "#EF5675",
  chr5 = "#FFA600", chr6 = "#8A5A44", chr7 = "#4E79A7", chr8 = "#59A14F",
  chr9 = "#D95F02", chr10 = "#9C755F", chr11 = "#F28E2B", chr12 = "#BAB0AC",
  chr13 = "#2F4B7C", chr14 = "#665191", chr15 = "#A05195", chr16 = "#D45087",
  chr17 = "#F95D6A", chr18 = "#FF7C43", chr19 = "#1F77B4", chr20 = "#17BECF",
  chr21 = "#6A3D9A", chr22 = "#B15928", chrX = "#1B9E77", chrY = "#7570B3"
)
target_cols[names(fixed_cols)[names(fixed_cols) %in% names(target_cols)]] <- fixed_cols[names(fixed_cols) %in% names(target_cols)]

color_for_target <- function(target) {
  col <- target_cols[target]
  if (is.na(col)) "#444444" else unname(col)
}

fmt_bp <- function(bp) {
  ifelse(bp >= 1e6, sprintf("%.2f Mb", bp / 1e6), ifelse(bp >= 1e3, sprintf("%.0f kb", bp / 1e3), sprintf("%.0f bp", bp)))
}

fmt_mb_axis <- function(bp) {
  sprintf("%.0f", bp / 1e6)
}

short_event_label <- function(event_id) {
  out <- event_id
  out <- gsub("PAR1_XY_positive_control", "PAR1 X/Y", out, fixed = TRUE)
  out <- gsub("PAN027_chr9q_chr3q_PHR_candidate", "PAN027 chr9q->chr3q", out, fixed = TRUE)
  out <- gsub("PAN028_chr9q_chr3q_PHR_candidate", "PAN028 chr9q->chr3q", out, fixed = TRUE)
  out
}

make_top_pairs <- function(n = 12) {
  out <- data.frame()
  for (basis in basis_order) {
    rows <- pairs[pairs$basis == basis, ]
    rows <- rows[order(-rows$bp, rows$query_chrom, rows$target_chrom), ]
    rows <- head(rows, n)
    rows$rank <- seq_len(nrow(rows))
    out <- rbind(out, rows)
  }
  out
}

top_pairs <- make_top_pairs()
write.table(top_pairs, file.path(out_dir, "top_pair_summary_for_plot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw_title <- function() {
  par(mar = c(0.1, 4.8, 1.0, 0.8), xaxs = "i", yaxs = "i", family = "sans")
  plot.new()
  text(0.00, 0.84, "Fig5 IMPG homolog-vs-interchromosomal class winners across the whole genome",
    adj = 0, cex = 1.06, font = 2, col = "#202020"
  )
  text(0.00, 0.56,
    "Each row is a full-length CHM13 query chromosome. Colored tracts are 2 kb windows merged by target chromosome where best interchromosomal IMPG identity exceeds best same-chromosome/homolog identity.",
    adj = 0, cex = 0.54, col = "#4F4F4F"
  )
  text(0.00, 0.34,
    "Light context marks PAR1/PAR2, acrocentric PHR intervals, centromeres, and the chr9q-to-chr3q candidate windows. The plotted windows are from the corrected depth-filtered PAN027 paternal vs PAN011 class-winner summaries.",
    adj = 0, cex = 0.50, col = "#666666"
  )
}

draw_legend <- function() {
  par(mar = c(0.1, 4.8, 0.0, 0.8), xaxs = "i", yaxs = "i", family = "sans")
  plot.new()
  legend_targets <- c("chr3", "chr9", "chr13", "chr14", "chr15", "chr21", "chr22", "chrX", "chrY")
  legend_targets <- legend_targets[legend_targets %in% names(target_cols)]
  x <- 0.00
  text(x, 0.70, "winning target chromosome", adj = 0, cex = 0.54, font = 2, col = "#303030")
  x <- 0.18
  for (target in legend_targets) {
    rect(x, 0.61, x + 0.018, 0.78, col = color_for_target(target), border = NA)
    text(x + 0.023, 0.695, target, adj = 0, cex = 0.48, col = "#303030")
    x <- x + 0.075
  }
  rect(0.00, 0.22, 0.018, 0.39, col = "#D8DDE2", border = "#B9C0C7", lwd = 0.3)
  text(0.023, 0.305, "same/homolog wins or no retained comparison", adj = 0, cex = 0.48, col = "#303030")
  rect(0.250, 0.22, 0.268, 0.39, col = "#111111", border = NA)
  text(0.273, 0.305, "candidate window outline", adj = 0, cex = 0.48, col = "#303030")
  rect(0.450, 0.22, 0.468, 0.39, col = grDevices::adjustcolor("#6A6A6A", alpha.f = 0.38), border = NA)
  text(0.473, 0.305, "PAR / acrocentric / centromere context", adj = 0, cex = 0.48, col = "#303030")
}

draw_context <- function(y_offset = 0) {
  par(xpd = NA)
  cent <- context_rows[context_rows$context_class == "centromere", ]
  if (nrow(cent) > 0) {
    for (i in seq_len(nrow(cent))) {
      y <- y_map[cent$chrom[i]] + y_offset
      if (is.na(y)) next
      rect(cent$start[i], y - 0.135, cent$end[i], y + 0.135,
        col = grDevices::adjustcolor("#8A8F94", alpha.f = 0.30), border = NA
      )
    }
  }
  par_rows <- context_rows[context_rows$context_class == "PAR", ]
  if (nrow(par_rows) > 0) {
    for (i in seq_len(nrow(par_rows))) {
      y <- y_map[par_rows$chrom[i]] + y_offset
      if (is.na(y)) next
      rect(par_rows$start[i], y - 0.20, par_rows$end[i], y + 0.20,
        col = grDevices::adjustcolor("#4C78A8", alpha.f = 0.35), border = "#4C78A8", lwd = 0.25
      )
    }
  }
  acro <- context_rows[context_rows$context_class == "acrocentric_PHR", ]
  if (nrow(acro) > 0) {
    for (i in seq_len(nrow(acro))) {
      y <- y_map[acro$chrom[i]] + y_offset
      if (is.na(y)) next
      rect(acro$start[i], y - 0.20, acro$end[i], y + 0.20,
        col = grDevices::adjustcolor("#A6761D", alpha.f = 0.32), border = NA
      )
    }
  }
  cand <- context_rows[context_rows$context_class == "candidate_window", ]
  if (nrow(cand) > 0) {
    cand <- cand[order(cand$chrom, cand$start), ]
    cand$label_rank <- ave(cand$start, cand$chrom, FUN = seq_along)
    for (i in seq_len(nrow(cand))) {
      y <- y_map[cand$chrom[i]] + y_offset
      if (is.na(y)) next
      rect(cand$start[i], y - 0.37, cand$end[i], y + 0.37, col = NA, border = "#111111", lwd = 0.6)
      text((cand$start[i] + cand$end[i]) / 2, y + 0.54 + (cand$label_rank[i] - 1) * 0.20, short_event_label(cand$label[i]),
        cex = 0.31, font = 2, col = "#111111"
      )
    }
  }
  par(xpd = FALSE)
}

draw_basis_panel <- function(basis, bottom_axis = FALSE) {
  rows <- tracts[tracts$basis == basis & tracts$query_chrom %in% chrom_order, ]
  rows <- rows[rows$query_end > rows$query_start, ]
  rows <- rows[order(match(rows$query_chrom, chrom_order), rows$query_start, -rows$bp), ]
  par(mar = c(ifelse(bottom_axis, 2.7, 0.45), 4.8, 1.0, 0.8), xaxs = "i", yaxs = "i", family = "sans")
  plot(NA, xlim = c(0, x_max), ylim = c(0.35, nrow(chroms) + 0.85), axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7), col = "#E7E9EB", lwd = 0.35)
  for (i in seq_len(nrow(chroms))) {
    y <- chroms$y[i]
    len <- chroms$length[i]
    segments(0, y, len, y, col = "#D8DDE2", lwd = 6.2, lend = "butt")
    rect(0, y - 0.16, len, y + 0.16, border = "#B9C0C7", col = NA, lwd = 0.25)
  }
  draw_context()
  if (nrow(rows) > 0) {
    for (i in seq_len(nrow(rows))) {
      y <- y_map[rows$query_chrom[i]]
      len <- len_map[rows$query_chrom[i]]
      if (is.na(y) || is.na(len)) next
      x0 <- max(0, rows$query_start[i])
      x1 <- min(len, rows$query_end[i])
      if (x1 <= x0) next
      alpha <- max(0.58, min(0.98, 0.52 + rows$mean_delta[i] * 8))
      col <- color_for_target(rows$target_chrom[i])
      segments(x0, y, x1, y, col = grDevices::adjustcolor(col, alpha.f = alpha), lwd = 7.4, lend = "butt")
      if (rows$bp[i] >= 20000 || rows$target_chrom[i] %in% c("chr3", "chr9", "chrX", "chrY")) {
        segments((x0 + x1) / 2, y - 0.27, (x0 + x1) / 2, y + 0.27, col = col, lwd = 0.45)
      }
    }
  }
  label_rows <- rows[rows$bp >= 30000 | rows$target_chrom %in% c("chr3", "chr9", "chrX", "chrY"), ]
  if (nrow(label_rows) > 0) {
    labels <- aggregate(bp ~ query_chrom + target_chrom, label_rows, sum)
    mids <- aggregate(query_start ~ query_chrom + target_chrom, label_rows, mean)
    labels$mid <- mids$query_start
    labels <- labels[order(match(labels$query_chrom, chrom_order), -labels$bp), ]
    for (chrom in unique(labels$query_chrom)) {
      cc <- labels[labels$query_chrom == chrom, ]
      cc <- head(cc, 2)
      for (j in seq_len(nrow(cc))) {
        y <- y_map[chrom]
        if (is.na(y)) next
        col <- color_for_target(cc$target_chrom[j])
        text(cc$mid[j], y + 0.34 + (j - 1) * 0.20, cc$target_chrom[j],
          cex = 0.34, font = 2, col = col,
          adj = ifelse(cc$mid[j] > x_max * 0.9, 1, ifelse(cc$mid[j] < x_max * 0.05, 0, 0.5))
        )
      }
    }
  }
  axis(2, at = chroms$y, labels = chroms$chrom, las = 1, tick = FALSE, cex.axis = 0.48)
  axis(3, at = seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7), labels = FALSE, tcl = -0.16, col = "#B8BCC0")
  stats <- overall[overall$basis == basis, ]
  panel_title <- paste0("SweepGA/FastGA mapping ", basis)
  text(0, nrow(chroms) + 0.62, panel_title, adj = 0, cex = 0.67, font = 2, col = "#202020")
  if (nrow(stats) == 1) {
    text(x_max, nrow(chroms) + 0.62,
      paste0(
        fmt_bp(stats$inter_beats_same_bp), " interchromosomal-winner bp; ",
        format(stats$inter_beats_same_windows, big.mark = ","), " windows of ",
        format(stats$both_same_and_inter, big.mark = ","), " same+inter comparisons"
      ),
      adj = 1, cex = 0.45, col = "#555555"
    )
  }
  if (bottom_axis) {
    ticks <- seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7)
    axis(1, at = ticks, labels = fmt_mb_axis(ticks), cex.axis = 0.55)
    mtext("native query chromosome coordinate (Mb)", side = 1, line = 1.85, cex = 0.58)
  }
}

draw_pair_sidebar <- function() {
  par(mar = c(2.7, 0.2, 1.0, 0.3), xaxs = "i", yaxs = "i", family = "sans")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0, 0.98, "largest interchromosomal-winner pairs", adj = 0, cex = 0.62, font = 2, col = "#202020")
  y <- 0.93
  for (basis in basis_order) {
    rows <- top_pairs[top_pairs$basis == basis, ]
    rows <- head(rows, 5)
    text(0, y, paste0(basis, " mapping"), adj = 0, cex = 0.50, font = 2, col = "#303030")
    y <- y - 0.035
    for (i in seq_len(nrow(rows))) {
      col <- color_for_target(rows$target_chrom[i])
      rect(0.00, y - 0.010, 0.030, y + 0.010, col = col, border = NA)
      text(0.040, y, paste0(rows$query_chrom[i], " -> ", rows$target_chrom[i]), adj = 0, cex = 0.41, col = "#303030")
      text(0.98, y, fmt_bp(rows$bp[i]), adj = 1, cex = 0.41, col = "#555555")
      y <- y - 0.030
    }
    y <- y - 0.030
  }
  text(0, y - 0.01, "Candidate context", adj = 0, cex = 0.50, font = 2, col = "#303030")
  y <- y - 0.050
  cand <- context_rows[context_rows$context_class == "candidate_window", ]
  for (i in seq_len(nrow(cand))) {
    text(0.00, y, short_event_label(cand$label[i]), adj = 0, cex = 0.39, col = "#303030")
    text(0.98, y, paste0(cand$chrom[i], ":", sprintf("%.2f", cand$start[i] / 1e6), "-", sprintf("%.2f", cand$end[i] / 1e6), " Mb"),
      adj = 1, cex = 0.37, col = "#555555"
    )
    y <- y - 0.032
  }
  y <- y - 0.020
  text(0, y, "Controls", adj = 0, cex = 0.50, font = 2, col = "#303030")
  y <- y - 0.034
  text(0, y, "PAR1/PAR2 intervals are drawn on chrX/chrY.", adj = 0, cex = 0.38, col = "#555555")
  y <- y - 0.030
  text(0, y, "Acrocentric PHR intervals are drawn on chr13/14/15/21/22.", adj = 0, cex = 0.38, col = "#555555")
}

draw_page <- function() {
  layout(
    matrix(c(1, 1, 2, 2, 3, 6, 4, 6, 5, 6), nrow = 5, byrow = TRUE),
    heights = c(0.32, 0.22, 1, 1, 1),
    widths = c(0.76, 0.24)
  )
  draw_title()
  draw_legend()
  draw_basis_panel("1:1", bottom_axis = FALSE)
  draw_basis_panel("4:4", bottom_axis = FALSE)
  draw_basis_panel("10:10", bottom_axis = TRUE)
  draw_pair_sidebar()
}

pdf(file.path(out_dir, paste0(prefix, ".pdf")), width = 15.8, height = 12.4, onefile = FALSE, useDingbats = FALSE)
draw_page()
dev.off()

png(file.path(out_dir, paste0(prefix, ".png")), width = 3160, height = 2480, res = 200, type = "cairo")
draw_page()
dev.off()

svg(file.path(out_dir, paste0(prefix, ".svg")), width = 15.8, height = 12.4, onefile = FALSE)
draw_page()
dev.off()
