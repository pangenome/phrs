#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
prefix <- "fig5_whole_genome_length_scaled_tracks"

segments <- read.delim(file.path(panel_dir, "length_scaled_track_segments.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
chroms <- read.delim(file.path(panel_dir, "length_scaled_track_chromosomes.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
summary_df <- read.delim(file.path(panel_dir, "length_scaled_track_summary.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

num_cols <- c("page_order", "method_order", "query_chrom_length", "segment_start", "segment_end", "segment_length_bp", "display_support_bp", "support_fraction", "mean_identity")
for (col in intersect(num_cols, names(segments))) {
  segments[[col]] <- suppressWarnings(as.numeric(segments[[col]]))
}
chroms$page_order <- as.integer(chroms$page_order)
chroms$query_chrom_order <- as.integer(chroms$query_chrom_order)
chroms$query_chrom_length <- as.numeric(chroms$query_chrom_length)
summary_df$page_order <- as.integer(summary_df$page_order)
summary_df$method_order <- as.integer(summary_df$method_order)
summary_df$display_support_bp <- as.numeric(summary_df$display_support_bp)

chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
method_order <- c(
  "untangle_strict_primary_path",
  "wfmash_p95_qgrid2kb_1to1_ani",
  "sweepga_fastga_f32_qgrid2kb_1to1_ani"
)

target_chroms <- sort(unique(segments$display_target_chrom[segments$display_state == "interchromosomal"]))
target_chroms <- target_chroms[target_chroms != ""]
base_pal <- grDevices::hcl(seq(15, 375, length.out = length(target_chroms) + 1)[-1], c = 70, l = 52)
target_cols <- setNames(base_pal, target_chroms)
target_cols[c("chr3", "chr9", "chrX", "chrY", "chr15", "chr16", "chr20")] <- c(
  "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#A6761D", "#666666", "#969696"
)
same_col <- "#BFC5CA"
background_col <- "#F1F2F3"
track_border <- "#D2D5D7"
grid_col <- "#E7E8EA"
text_col <- "#252525"

fmt_mb <- function(bp) sprintf("%.0f", bp / 1e6)
fmt_bp <- function(bp) {
  ifelse(bp >= 1e6, sprintf("%.1f Mb", bp / 1e6), ifelse(bp >= 1e3, sprintf("%.0f kb", bp / 1e3), sprintf("%.0f bp", bp)))
}

short_event_label <- function(event_id) {
  if (event_id == "PAR1_XY_positive_control") {
    return("PAR1")
  }
  if (event_id == "PAN027_chr9q_chr3q_PHR_candidate") {
    return("PAN027 chr9q")
  }
  if (event_id == "PAN028_chr9q_chr3q_PHR_candidate") {
    return("PAN028 chr9q")
  }
  event_id
}

color_for_target <- function(target) {
  if (target == "same_chromosome") {
    return(same_col)
  }
  if (!is.na(target_cols[target])) {
    return(unname(target_cols[target]))
  }
  "#5B5B5B"
}

method_label <- function(method_id) {
  rows <- segments[segments$method_id == method_id, ]
  if (nrow(rows) == 0) {
    return(method_id)
  }
  rows$method_label[1]
}

short_method_label <- function(method_id) {
  if (method_id == "untangle_strict_primary_path") {
    return("Untangle strict")
  }
  if (method_id == "wfmash_p95_qgrid2kb_1to1_ani") {
    return("wfmash -p95, 2 kb query-grid, 1:1 ANI")
  }
  if (method_id == "sweepga_fastga_f32_qgrid2kb_1to1_ani") {
    return("SweepGA/FastGA -f32, 2 kb query-grid, 1:1 ANI")
  }
  method_label(method_id)
}

draw_header <- function(page_id, page_label, page_chroms) {
  par(mar = c(0.2, 4.9, 1.0, 0.8), xaxs = "i", yaxs = "i", family = "sans")
  plot.new()
  title <- paste0("Fig5 whole-genome target support: ", page_label)
  text(0.0, 0.82, title, adj = 0, font = 2, cex = 1.06, col = text_col)
  subtitle <- "Each chromosome is one native query-coordinate row; row length is proportional to actual query chromosome length. Colored ticks/blocks mark retained interchromosomal target support."
  text(0.0, 0.50, subtitle, adj = 0, cex = 0.58, col = "#555555")
  legend_items <- c("same chromosome", "chr3", "chr9", "chrX", "chrY", "other target")
  legend_cols <- c(same_col, color_for_target("chr3"), color_for_target("chr9"), color_for_target("chrX"), color_for_target("chrY"), "#5B5B5B")
  x <- 0.0
  for (i in seq_along(legend_items)) {
    rect(x, 0.15, x + 0.020, 0.27, col = legend_cols[i], border = "#555555", lwd = 0.35)
    text(x + 0.025, 0.21, legend_items[i], adj = 0, cex = 0.52, col = text_col)
    x <- x + 0.145
  }
}

draw_method_panel <- function(page_id, method_id, bottom_axis = FALSE) {
  page_chroms <- chroms[chroms$page_id == page_id, ]
  page_chroms <- page_chroms[match(chr_order, page_chroms$query_chrom), ]
  page_chroms <- page_chroms[!is.na(page_chroms$query_chrom), ]
  page_chroms <- page_chroms[page_chroms$query_chrom_length > 0, ]
  page_chroms$y <- rev(seq_len(nrow(page_chroms)))
  y_map <- setNames(page_chroms$y, page_chroms$query_chrom)
  len_map <- setNames(page_chroms$query_chrom_length, page_chroms$query_chrom)
  x_max <- max(page_chroms$query_chrom_length, na.rm = TRUE)
  rows <- segments[segments$page_id == page_id & segments$method_id == method_id, ]
  rows <- rows[rows$segment_end > rows$segment_start, ]
  rows <- rows[order(rows$display_state == "interchromosomal", rows$query_chrom, rows$segment_start), ]

  par(mar = c(ifelse(bottom_axis, 2.8, 0.45), 4.9, 1.15, 0.8), xaxs = "i", yaxs = "i", family = "sans")
  plot(NA, xlim = c(0, x_max), ylim = c(0.35, nrow(page_chroms) + 0.7), axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7), col = grid_col, lwd = 0.35)
  for (i in seq_len(nrow(page_chroms))) {
    y <- page_chroms$y[i]
    len <- page_chroms$query_chrom_length[i]
    segments(0, y, len, y, col = background_col, lwd = 7.2, lend = "butt")
    rect(0, y - 0.16, len, y + 0.16, border = track_border, col = NA, lwd = 0.25)
  }

  same_rows <- rows[rows$display_state == "same_chromosome", ]
  inter_rows <- rows[rows$display_state == "interchromosomal", ]
  if (nrow(same_rows) > 0) {
    for (i in seq_len(nrow(same_rows))) {
      y <- y_map[same_rows$query_chrom[i]]
      if (is.na(y)) next
      x0 <- max(0, same_rows$segment_start[i])
      x1 <- min(len_map[same_rows$query_chrom[i]], same_rows$segment_end[i])
      segments(x0, y, x1, y, col = adjustcolor(same_col, alpha.f = 0.62), lwd = 6.2, lend = "butt")
    }
  }
  if (nrow(inter_rows) > 0) {
    for (i in seq_len(nrow(inter_rows))) {
      y <- y_map[inter_rows$query_chrom[i]]
      if (is.na(y)) next
      x0 <- max(0, inter_rows$segment_start[i])
      x1 <- min(len_map[inter_rows$query_chrom[i]], inter_rows$segment_end[i])
      target <- inter_rows$display_target_chrom[i]
      col <- color_for_target(target)
      alpha <- ifelse(inter_rows$source_row_type[i] == "untangle_segment", 0.95, max(0.42, min(0.95, inter_rows$support_fraction[i])))
      segments(x0, y, x1, y, col = adjustcolor(col, alpha.f = alpha), lwd = 7.0, lend = "butt")
      xm <- (x0 + x1) / 2
      segments(xm, y - 0.25, xm, y + 0.25, col = col, lwd = 0.55)
    }
  }

  label_rows <- inter_rows[inter_rows$display_support_bp > 0, ]
  if (nrow(label_rows) > 0) {
    label_rows$mid <- (label_rows$segment_start + label_rows$segment_end) / 2
    label_rows$weight <- pmax(label_rows$display_support_bp, 1)
    groups <- unique(label_rows[, c("query_chrom", "display_target_chrom")])
    label_candidates <- data.frame()
    for (i in seq_len(nrow(groups))) {
      rr <- label_rows[label_rows$query_chrom == groups$query_chrom[i] & label_rows$display_target_chrom == groups$display_target_chrom[i], ]
      total <- sum(rr$display_support_bp, na.rm = TRUE)
      mid <- sum(rr$mid * rr$weight, na.rm = TRUE) / sum(rr$weight, na.rm = TRUE)
      label_candidates <- rbind(label_candidates, data.frame(
        query_chrom = groups$query_chrom[i],
        target = groups$display_target_chrom[i],
        total = total,
        mid = mid,
        stringsAsFactors = FALSE
      ))
    }
    for (chrom in unique(label_candidates$query_chrom)) {
      cc <- label_candidates[label_candidates$query_chrom == chrom, ]
      cc <- cc[order(-cc$total), ]
      cc <- head(cc, 2)
      for (i in seq_len(nrow(cc))) {
        y <- y_map[chrom]
        col <- color_for_target(cc$target[i])
        adj <- ifelse(cc$mid[i] > x_max * 0.90, 1, ifelse(cc$mid[i] < x_max * 0.08, 0, 0.5))
        text(cc$mid[i], y + 0.31 + (i - 1) * 0.22, cc$target[i], cex = 0.36, font = 2, col = col, adj = adj)
      }
    }
  }

  call_rows <- rows[rows$callout_event_id != "" & rows$segment_end > rows$segment_start, ]
  if (nrow(call_rows) > 0) {
    event_ids <- unique(call_rows$callout_event_id)
    for (event_id in event_ids) {
      rr_event <- call_rows[call_rows$callout_event_id == event_id, ]
      chroms_for_event <- unique(rr_event$query_chrom)
      for (chrom in chroms_for_event) {
        rr <- rr_event[rr_event$query_chrom == chrom, ]
        y <- y_map[chrom]
        if (is.na(y)) next
        x0 <- max(0, min(rr$segment_start, na.rm = TRUE))
        x1 <- min(len_map[chrom], max(rr$segment_end, na.rm = TRUE))
        xm <- (x0 + x1) / 2
        rect(x0, y - 0.37, x1, y + 0.37, border = "#111111", col = NA, lwd = 0.65)
        inter_rr <- rr[rr$display_state == "interchromosomal" & rr$display_support_bp > 0, ]
        target_suffix <- ""
        if (nrow(inter_rr) > 0) {
          target_totals <- aggregate(display_support_bp ~ display_target_chrom, inter_rr, sum)
          target_totals <- target_totals[order(-target_totals$display_support_bp), ]
          target_suffix <- paste0(" ", paste(head(target_totals$display_target_chrom, 2), collapse = "/"))
        }
        text(xm, y + 0.55, paste0(short_event_label(event_id), target_suffix), cex = 0.32, font = 2, col = "#111111")
      }
    }
  }

  axis(2, at = page_chroms$y, labels = page_chroms$query_chrom, las = 1, tick = FALSE, cex.axis = 0.49)
  axis(3, at = seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7), labels = FALSE, tcl = -0.16, col = "#B8BCC0")
  text(0, nrow(page_chroms) + 0.52, short_method_label(method_id), adj = 0, cex = 0.68, font = 2, col = text_col)
  inter_total <- sum(rows$display_support_bp[rows$display_state == "interchromosomal"], na.rm = TRUE)
  text(x_max, nrow(page_chroms) + 0.52, paste0("interchromosomal display support: ", fmt_bp(inter_total)), adj = 1, cex = 0.48, col = "#555555")
  if (bottom_axis) {
    axis(1, at = seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7), labels = fmt_mb(seq(0, ceiling(x_max / 5e7) * 5e7, by = 5e7)), cex.axis = 0.55)
    mtext("native query chromosome coordinate (Mb)", side = 1, line = 1.85, cex = 0.58)
  }
}

draw_page <- function(page_id) {
  page_chroms <- chroms[chroms$page_id == page_id, ]
  page_label <- page_chroms$page_label[1]
  layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(0.32, 1, 1, 1))
  draw_header(page_id, page_label, page_chroms)
  draw_method_panel(page_id, method_order[1], bottom_axis = FALSE)
  draw_method_panel(page_id, method_order[2], bottom_axis = FALSE)
  draw_method_panel(page_id, method_order[3], bottom_axis = TRUE)
}

page_ids <- unique(chroms$page_id[order(chroms$page_order)])

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 15.2, height = 12.0, onefile = TRUE, useDingbats = FALSE)
for (page_id in page_ids) {
  draw_page(page_id)
}
dev.off()

for (page_id in page_ids) {
  pdf(file.path(panel_dir, paste0(prefix, ".", page_id, ".pdf")), width = 15.2, height = 12.0, onefile = FALSE, useDingbats = FALSE)
  draw_page(page_id)
  dev.off()
  png(file.path(panel_dir, paste0(prefix, ".", page_id, ".png")), width = 3040, height = 2400, res = 200, type = "cairo")
  draw_page(page_id)
  dev.off()
  svg(file.path(panel_dir, paste0(prefix, ".", page_id, ".svg")), width = 15.2, height = 12.0, onefile = FALSE)
  draw_page(page_id)
  dev.off()
}
