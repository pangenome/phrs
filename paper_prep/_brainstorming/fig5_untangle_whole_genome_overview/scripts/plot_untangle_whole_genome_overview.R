#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
prefix <- "fig5_untangle_whole_genome_overview"
segments_path <- file.path(panel_dir, "untangle_whole_genome_segments.tsv")
summary_path <- file.path(panel_dir, "untangle_whole_genome_summary.tsv")

segments <- read.delim(segments_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

segments$query_arm_order <- as.integer(segments$query_arm_order)
segments$local_query_start_0based <- as.numeric(segments$local_query_start_0based)
segments$local_query_end_0based_exclusive <- as.numeric(segments$local_query_end_0based_exclusive)
segments$query_length <- as.numeric(segments$query_length)
segments$segment_length_bp <- as.numeric(segments$segment_length_bp)
segments$interchromosomal <- as.character(segments$interchromosomal)

summary$query_arm_order <- as.integer(summary$query_arm_order)
summary$strict_total_bp <- as.numeric(summary$strict_total_bp)
summary$interchromosomal_bp <- as.numeric(summary$interchromosomal_bp)
summary$side_fragment_caveat_bp <- as.numeric(summary$side_fragment_caveat_bp)

arm_table <- unique(segments[, c("query_arm", "query_arm_order")])
arm_table <- arm_table[order(arm_table$query_arm_order), ]
arms <- arm_table$query_arm
n_arms <- length(arms)

pair_table <- unique(segments[, c("pair", "transmission")])
pair_table <- pair_table[order(pair_table$pair), ]
pair_order <- pair_table$pair
row_y <- rev(seq_along(pair_order))
names(row_y) <- pair_order
pair_short_labels <- setNames(
  c(
    "PAN027<-PAN010",
    "PAN027<-PAN011",
    "PAN028<-PAN027"
  )[match(pair_order, c("PAN027_vs_PAN010", "PAN027_vs_PAN011", "PAN028_vs_PAN027"))],
  pair_order
)
pair_short_labels[is.na(pair_short_labels)] <- pair_order[is.na(pair_short_labels)]

target_chroms <- sort(unique(segments$target_chrom[segments$interchromosomal == "1"]))
target_palette <- setNames(hcl(seq(15, 375, length.out = length(target_chroms) + 1)[-1], c = 72, l = 52), target_chroms)
target_palette[c("chr3", "chr9", "chrX", "chrY", "chr15", "chr16", "chr20")] <- c(
  "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#969696"
)
same_col <- "#D5D7D8"
same_border <- "#A6AAAD"

fmt_bp <- function(x) {
  ifelse(x >= 1e6, sprintf("%.2f Mb", x / 1e6), ifelse(x >= 1000, sprintf("%.1f kb", x / 1000), sprintf("%d bp", round(x))))
}

segment_fill <- function(row) {
  if (row[["interchromosomal"]] == "0") {
    return(same_col)
  }
  chrom <- row[["target_chrom"]]
  if (!is.na(target_palette[chrom])) {
    return(unname(target_palette[chrom]))
  }
  "#7F7F7F"
}

segment_x0 <- function(df) {
  (df$query_arm_order - 1) + (df$local_query_start_0based / pmax(df$query_length, 1))
}

segment_x1 <- function(df) {
  (df$query_arm_order - 1) + (df$local_query_end_0based_exclusive / pmax(df$query_length, 1))
}

draw_main_barcode <- function() {
  par(mar = c(4.7, 2.0, 3.2, 1.4), xaxs = "i", yaxs = "i", family = "sans")
  plot(NA, xlim = c(0, n_arms), ylim = c(0.35, length(pair_order) + 0.95), axes = FALSE, xlab = "", ylab = "")
  title("Fig5 strict untangle primary-path geometry across terminal windows", cex.main = 1.05, line = 1.7)
  mtext("One block per queried native terminal window/arm; coordinates are native assembly windows, not CHM13 projections", side = 3, line = 0.35, cex = 0.70)

  for (i in seq_len(n_arms)) {
    bg <- if (i %% 2 == 0) "#F7F7F7" else "#FFFFFF"
    rect(i - 1, 0.35, i, length(pair_order) + 0.95, col = bg, border = NA)
  }
  abline(v = seq(0, n_arms), col = "#E4E4E4", lwd = 0.35)

  for (pair in pair_order) {
    y <- row_y[pair]
    rows <- segments[segments$pair == pair, ]
    rect(0, y - 0.18, n_arms, y + 0.18, col = "#F1F1F1", border = "#CFCFCF", lwd = 0.5)
    text(0.08, y + 0.31, pair_short_labels[pair], adj = 0, cex = 0.62, font = 2, col = "#333333")
    rows <- rows[order(rows$query_arm_order, rows$local_query_start_0based, rows$local_query_end_0based_exclusive), ]
    for (j in seq_len(nrow(rows))) {
      x0 <- segment_x0(rows[j, ])
      x1 <- segment_x1(rows[j, ])
      y0 <- y - ifelse(rows$side_fragment_caveat[j] == "yes", 0.28, 0.16)
      y1 <- y + ifelse(rows$side_fragment_caveat[j] == "yes", -0.05, 0.16)
      col <- segment_fill(rows[j, ])
      border <- ifelse(rows$interchromosomal[j] == "0", same_border, "#333333")
      lty <- ifelse(rows$target_haplotype_number[j] %in% c("2"), 2, 1)
      if (rows$side_fragment_caveat[j] == "yes") {
        border <- "#333333"
        lty <- 3
      }
      rect(x0, y0, x1, y1, col = col, border = border, lwd = 0.35, lty = lty)
    }
  }

  event_rows <- segments[segments$event_id != "", ]
  event_keys <- unique(event_rows[, c("event_id", "event_label", "pair", "query_arm_order")])
  event_keys <- event_keys[order(event_keys$pair, event_keys$query_arm_order), ]
  label_offsets <- c(PAR1_XY_positive_control = 0.40, PAN027_chr9q_chr3q_PHR_candidate = 0.57, PAN028_chr9q_chr3q_PHR_candidate = 0.40)
  for (i in seq_len(nrow(event_keys))) {
    erows <- event_rows[event_rows$event_id == event_keys$event_id[i] & event_rows$pair == event_keys$pair[i], ]
    x0 <- min(segment_x0(erows))
    x1 <- max(segment_x1(erows))
    y <- row_y[event_keys$pair[i]]
    rect(x0, y - 0.32, x1, y + 0.32, border = "#111111", lwd = 1.0, col = NA)
    xm <- (x0 + x1) / 2
    label_y <- min(length(pair_order) + 0.88, y + label_offsets[event_keys$event_id[i]])
    segments(xm, y + 0.34, xm, label_y - 0.08, col = "#222222", lwd = 0.55)
    text(xm, label_y, event_keys$event_label[i], cex = 0.58, font = 2, adj = c(0.5, 0))
  }

  caveats <- segments[segments$side_fragment_caveat == "yes", ]
  if (nrow(caveats) > 0) {
    points((segment_x0(caveats) + segment_x1(caveats)) / 2, row_y[caveats$pair] - 0.38,
           pch = 24, bg = "#FFFFFF", col = "#222222", cex = 0.55)
  }

  axis(1, at = seq_len(n_arms) - 0.5, labels = arms, las = 2, tick = FALSE, cex.axis = 0.46)
  mtext("Genome-ordered queried arms; each arm block is scaled to its native terminal source window", side = 1, line = 3.55, cex = 0.66)

  legend_x <- 0.15
  legend_y <- 0.58
  legend_labels <- c("same chromosome", "chr3 donor", "chrY PAR", "other interchrom.", "side/caveat marker")
  legend_cols <- c(same_col, target_palette[["chr3"]], target_palette[["chrY"]], "#4B9CD3", "#FFFFFF")
  for (i in seq_along(legend_labels)) {
    x <- legend_x + (i - 1) * 4.2
    if (i == length(legend_labels)) {
      points(x + 0.10, legend_y, pch = 24, bg = "#FFFFFF", col = "#222222", cex = 0.65)
    } else {
      rect(x, legend_y - 0.06, x + 0.22, legend_y + 0.06, col = legend_cols[i], border = "#555555", lwd = 0.35)
    }
    text(x + 0.30, legend_y, legend_labels[i], adj = 0, cex = 0.53)
  }
}

draw_heatmap <- function() {
  par(mar = c(5.2, 5.0, 2.7, 1.0), xaxs = "i", yaxs = "i", family = "sans")
  inter <- summary[summary$interchromosomal_bp > 0, ]
  target_totals <- aggregate(interchromosomal_bp ~ target_arm, inter, sum)
  target_totals <- target_totals[order(-target_totals$interchromosomal_bp), ]
  must_keep <- c("chr3q", "chrYp", "chr15q", "chr16q", "chr20q")
  top_targets <- unique(c(must_keep[must_keep %in% target_totals$target_arm], head(target_totals$target_arm, 14)))
  top_targets <- top_targets[seq_len(min(length(top_targets), 16))]
  query_arms <- arms
  mat <- matrix(0, nrow = length(query_arms), ncol = length(top_targets), dimnames = list(query_arms, top_targets))
  for (i in seq_len(nrow(inter))) {
    q <- inter$query_arm[i]
    t <- inter$target_arm[i]
    if (q %in% query_arms && t %in% top_targets) {
      mat[q, t] <- mat[q, t] + inter$interchromosomal_bp[i]
    }
  }
  z <- log10(mat + 1)
  plot(NA, xlim = c(0, ncol(z)), ylim = c(0, nrow(z)), axes = FALSE, xlab = "", ylab = "")
  title("Query-arm vs target-arm support", cex.main = 0.86, line = 1.1)
  pal <- colorRampPalette(c("#FFFFFF", "#D7E8F5", "#6AA6C8", "#084B7A"))(80)
  max_z <- max(z)
  for (i in seq_len(nrow(z))) {
    for (j in seq_len(ncol(z))) {
      idx <- if (max_z == 0) 1 else max(1, ceiling((z[i, j] / max_z) * length(pal)))
      rect(j - 1, nrow(z) - i, j, nrow(z) - i + 1, col = pal[idx], border = "#EFEFEF", lwd = 0.25)
    }
  }
  axis(1, at = seq_len(ncol(z)) - 0.5, labels = colnames(z), las = 2, tick = FALSE, cex.axis = 0.54)
  axis(2, at = nrow(z) - seq_len(nrow(z)) + 0.5, labels = rownames(z), las = 1, tick = FALSE, cex.axis = 0.38)
  mtext("log10(interchromosomal strict bp + 1)", side = 3, line = 0.10, cex = 0.58)
}

draw_summary_bars <- function() {
  par(mar = c(5.2, 5.5, 2.7, 1.1), xaxs = "i", yaxs = "i", family = "sans")
  arm_totals <- aggregate(interchromosomal_bp ~ query_arm + query_arm_order, summary, sum)
  arm_totals <- arm_totals[order(-arm_totals$interchromosomal_bp), ]
  arm_totals <- arm_totals[arm_totals$interchromosomal_bp > 0, ]
  arm_totals <- head(arm_totals, 14)
  arm_totals <- arm_totals[order(arm_totals$interchromosomal_bp), ]
  x_max <- max(arm_totals$interchromosomal_bp) / 1000
  plot(NA, xlim = c(0, x_max * 1.18), ylim = c(0.5, nrow(arm_totals) + 0.5), axes = FALSE, xlab = "", ylab = "")
  title("Largest switch burdens by query arm", cex.main = 0.86, line = 1.1)
  for (i in seq_len(nrow(arm_totals))) {
    kb <- arm_totals$interchromosomal_bp[i] / 1000
    rect(0, i - 0.30, kb, i + 0.30, col = "#4B9CD3", border = "#2D5F7A", lwd = 0.35)
    text(kb + x_max * 0.025, i, fmt_bp(arm_totals$interchromosomal_bp[i]), adj = 0, cex = 0.52)
  }
  axis(2, at = seq_len(nrow(arm_totals)), labels = arm_totals$query_arm, las = 1, tick = FALSE, cex.axis = 0.62)
  axis(1, cex.axis = 0.58)
  mtext("interchromosomal strict primary-path support (kb)", side = 1, line = 3.1, cex = 0.62)
  note <- "Candidate windows: PAR1, PAN027 chr9q->chr3q, PAN028 chr9q->chr3q. Dashed/triangular marks retain side fragments as caveats."
  mtext(note, side = 3, line = 0.05, cex = 0.49)
}

draw_all <- function() {
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(1.25, 1.05))
  draw_main_barcode()
  draw_heatmap()
  draw_summary_bars()
}

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 15.5, height = 9.2, useDingbats = FALSE)
draw_all()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 15.5, height = 9.2, onefile = TRUE)
draw_all()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 3100, height = 1840, res = 200, type = "cairo")
draw_all()
dev.off()
