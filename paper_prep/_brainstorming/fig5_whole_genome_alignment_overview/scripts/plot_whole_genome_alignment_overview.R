#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
prefix <- "fig5_whole_genome_alignment_overview"

binned <- read.delim(file.path(panel_dir, "whole_genome_binned_support.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
matrix_df <- read.delim(file.path(panel_dir, "whole_genome_support_matrix.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.delim(file.path(panel_dir, "whole_genome_method_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

binned$query_bin_start <- as.numeric(binned$query_bin_start)
binned$query_bin_end <- as.numeric(binned$query_bin_end)
binned$query_chrom_length <- as.numeric(binned$query_chrom_length)
binned$support_fraction <- as.numeric(binned$support_fraction)
binned$winner_support_bp <- as.numeric(binned$winner_support_bp)
binned$winner_mean_identity <- as.numeric(binned$winner_mean_identity)
binned$top_interchrom_support_bp <- as.numeric(binned$top_interchrom_support_bp)
matrix_df$support_bp <- as.numeric(matrix_df$support_bp)

chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
chrom_lengths <- aggregate(query_chrom_length ~ query_chrom, binned, max)
chrom_lengths <- chrom_lengths[match(chr_order[chr_order %in% chrom_lengths$query_chrom], chrom_lengths$query_chrom), ]
chrom_lengths$offset <- c(0, cumsum(head(chrom_lengths$query_chrom_length, -1)))
chrom_lengths$mid <- chrom_lengths$offset + chrom_lengths$query_chrom_length / 2
offset_map <- setNames(chrom_lengths$offset, chrom_lengths$query_chrom)
length_map <- setNames(chrom_lengths$query_chrom_length, chrom_lengths$query_chrom)
genome_end <- max(chrom_lengths$offset + chrom_lengths$query_chrom_length)

method_order <- manifest$method_id[manifest$status == "OK"]
method_order <- method_order[method_order %in% unique(binned$method_id)]
short_method_label <- function(method_id) {
  parts <- manifest[manifest$method_id == method_id, ]
  if (nrow(parts) == 0) {
    return(method_id)
  }
  p <- parts[1, ]
  if (p$method_family == "untangle") {
    return("untangle strict")
  }
  freq <- ifelse(p$fastga_frequency == "NA", "", paste0(" f", p$fastga_frequency))
  chop <- ifelse(p$chop_length_bp == "NA", "", paste0(" L", as.numeric(p$chop_length_bp) / 1000, "kb"))
  comp <- sub("_joint$", "", p$comparison_id)
  comp <- gsub("PAN027pat_vs_PAN011", "027pat/011", comp)
  comp <- gsub("PAN027mat_vs_PAN010", "027mat/010", comp)
  comp <- gsub("PAN028mat_vs_PAN027", "028mat/027", comp)
  paste0(gsub("SweepGA/FastGA", "SG/FG", p$method_label), freq, chop, " ", comp)
}
row_y <- rev(seq_along(method_order))
names(row_y) <- method_order

target_chroms <- sort(unique(binned$winner_target_chrom[binned$winner_target_chrom != "no_support"]))
target_chroms <- target_chroms[order(match(target_chroms, chr_order), target_chroms)]
pal <- grDevices::hcl(seq(10, 370, length.out = length(target_chroms) + 1)[-1], c = 70, l = 55)
target_cols <- setNames(pal, target_chroms)
target_cols[c("chr3", "chr9", "chrX", "chrY")] <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A")
target_cols["same_chromosome"] <- "#D4D7DA"
target_cols["no_support"] <- "#F2F2F2"

row_fill <- function(row) {
  fam <- row[["winner_target_family"]]
  chrom <- row[["winner_target_chrom"]]
  if (fam == "no_support") {
    return(target_cols[["no_support"]])
  }
  if (fam == "same_chromosome") {
    return(target_cols[["same_chromosome"]])
  }
  if (!is.na(target_cols[chrom])) {
    return(unname(target_cols[chrom]))
  }
  "#666666"
}

draw_tracks <- function() {
  par(mar = c(3.2, 9.8, 2.5, 1.0), xaxs = "i", yaxs = "i", family = "sans")
  label_space <- 0.19 * genome_end
  plot(NA, xlim = c(-label_space, genome_end), ylim = c(0.4, length(method_order) + 0.8), axes = FALSE, xlab = "", ylab = "")
  title("Fig5 whole-genome retained-target overview after strict / 1:1 ANI filtering", cex.main = 1.02, line = 1.1)
  mtext("1 Mb query-coordinate bins. Neutral = no retained support; gray = same query/target chromosome; colored bins are dominant retained target chromosomes with direct labels for recurrent interchromosomal support.", side = 3, line = 0.05, cex = 0.58)

  for (i in seq_len(nrow(chrom_lengths))) {
    x0 <- chrom_lengths$offset[i]
    x1 <- x0 + chrom_lengths$query_chrom_length[i]
    rect(x0, 0.4, x1, length(method_order) + 0.8, col = ifelse(i %% 2 == 0, "#FAFAFA", "#FFFFFF"), border = NA)
    segments(x0, 0.4, x0, length(method_order) + 0.8, col = "#E4E4E4", lwd = 0.4)
  }
  segments(genome_end, 0.4, genome_end, length(method_order) + 0.8, col = "#E4E4E4", lwd = 0.4)

  for (mid in method_order) {
    y <- row_y[mid]
    rows <- binned[binned$method_id == mid, ]
    rows <- rows[order(match(rows$query_chrom, chr_order), rows$query_bin_start), ]
    rect(0, y - 0.25, genome_end, y + 0.25, col = "#F5F5F5", border = "#CFCFCF", lwd = 0.25)
    for (j in seq_len(nrow(rows))) {
      x0 <- offset_map[rows$query_chrom[j]] + rows$query_bin_start[j]
      x1 <- offset_map[rows$query_chrom[j]] + pmin(rows$query_bin_end[j], length_map[rows$query_chrom[j]])
      col <- row_fill(rows[j, ])
      alpha_col <- adjustcolor(col, alpha.f = ifelse(rows$no_support[j] == "yes", 0.65, max(0.35, rows$support_fraction[j])))
      border <- ifelse(rows$callout_event_id[j] == "", NA, "#111111")
      rect(x0, y - 0.20, x1, y + 0.20, col = alpha_col, border = border, lwd = 0.25)
    }
    text(-0.012 * genome_end, y, short_method_label(mid), adj = 1, cex = 0.43)

    inter <- rows[rows$top_interchrom_target_chrom != "none" & rows$top_interchrom_support_bp > 0, ]
    if (nrow(inter) > 0) {
      top <- aggregate(top_interchrom_support_bp ~ top_interchrom_target_chrom, inter, sum)
      top <- top[order(-top$top_interchrom_support_bp), ]
      keep <- head(top$top_interchrom_target_chrom, 3)
      for (chrom in keep) {
        rr <- inter[inter$top_interchrom_target_chrom == chrom, ]
        rr <- rr[which.max(rr$top_interchrom_support_bp), ]
        x <- offset_map[rr$query_chrom] + (rr$query_bin_start + rr$query_bin_end) / 2
        text(x, y + 0.27, chrom, cex = 0.31, font = 2, col = "#222222")
      }
    }
  }
  axis(1, at = chrom_lengths$mid, labels = chrom_lengths$query_chrom, las = 2, tick = FALSE, cex.axis = 0.43)
  mtext("Query chromosome coordinates (native sample assemblies; chromosomes laid end-to-end)", side = 1, line = 2.35, cex = 0.56)
  legend_x <- 0.01 * genome_end
  legend_y <- 0.65
  legend_items <- c("no support", "same chrom.", "chr3", "chr9", "chrX", "chrY")
  legend_cols <- c(target_cols["no_support"], target_cols["same_chromosome"], target_cols["chr3"], target_cols["chr9"], target_cols["chrX"], target_cols["chrY"])
  for (i in seq_along(legend_items)) {
    x <- legend_x + (i - 1) * 0.075 * genome_end
    rect(x, legend_y - 0.07, x + 0.012 * genome_end, legend_y + 0.07, col = legend_cols[i], border = "#555555", lwd = 0.25)
    text(x + 0.015 * genome_end, legend_y, legend_items[i], adj = 0, cex = 0.42)
  }
}

draw_matrix <- function() {
  par(mar = c(4.9, 6.7, 2.2, 0.9), xaxs = "i", yaxs = "i", family = "sans")
  inter <- matrix_df[matrix_df$interchromosomal == "yes" & matrix_df$support_bp > 0, ]
  inter$setting <- paste(inter$method_label, ifelse(inter$fastga_frequency == "NA", "", paste0("f", inter$fastga_frequency)), ifelse(inter$chop_length_bp == "NA", "", paste0("L", as.numeric(inter$chop_length_bp) / 1000, "kb")))
  agg <- aggregate(support_bp ~ setting + target_chrom, inter, sum)
  target_totals <- aggregate(support_bp ~ target_chrom, agg, sum)
  target_totals <- target_totals[order(-target_totals$support_bp), ]
  targets <- unique(c("chr3", "chr9", "chrX", "chrY", head(target_totals$target_chrom, 12)))
  targets <- targets[targets %in% target_totals$target_chrom]
  settings <- unique(agg$setting)
  mat <- matrix(0, nrow = length(settings), ncol = length(targets), dimnames = list(settings, targets))
  for (i in seq_len(nrow(agg))) {
    if (agg$setting[i] %in% settings && agg$target_chrom[i] %in% targets) {
      mat[agg$setting[i], agg$target_chrom[i]] <- agg$support_bp[i]
    }
  }
  z <- log10(mat + 1)
  plot(NA, xlim = c(0, ncol(z)), ylim = c(0, nrow(z)), axes = FALSE, xlab = "", ylab = "")
  title("Compact interchromosomal support matrix by method/chop", cex.main = 0.86, line = 1.0)
  colors <- colorRampPalette(c("#FFFFFF", "#E7F0F4", "#8EBAD0", "#276A8C"))(80)
  max_z <- max(z)
  for (i in seq_len(nrow(z))) {
    for (j in seq_len(ncol(z))) {
      idx <- if (max_z == 0) 1 else max(1, ceiling((z[i, j] / max_z) * length(colors)))
      rect(j - 1, nrow(z) - i, j, nrow(z) - i + 1, col = colors[idx], border = "#EFEFEF", lwd = 0.25)
      if (mat[i, j] > 0) {
        label <- ifelse(mat[i, j] >= 1e6, sprintf("%.1fM", mat[i, j] / 1e6), sprintf("%.0fk", mat[i, j] / 1000))
        text(j - 0.5, nrow(z) - i + 0.5, label, cex = 0.32, col = "#1F1F1F")
      }
    }
  }
  axis(1, at = seq_len(ncol(z)) - 0.5, labels = colnames(z), las = 2, tick = FALSE, cex.axis = 0.55)
  axis(2, at = nrow(z) - seq_len(nrow(z)) + 0.5, labels = rownames(z), las = 1, tick = FALSE, cex.axis = 0.43)
  mtext("Target chromosome labels printed in cells; full query-arm x target-arm table is in whole_genome_support_matrix.tsv", side = 3, line = 0.05, cex = 0.50)
}

draw_callouts <- function() {
  par(mar = c(3.3, 9.8, 2.0, 1.0), xaxs = "i", yaxs = "i", family = "sans")
  call <- binned[binned$callout_event_id != "", ]
  if (nrow(call) == 0) {
    plot.new()
    text(0.5, 0.5, "No callout bins found")
    return()
  }
  events <- unique(call[, c("callout_event_id", "callout_label", "query_chrom", "callout_window_start", "callout_window_end")])
  events <- events[order(events$callout_event_id), ]
  row_keys <- method_order
  plot(NA, xlim = c(-0.55, nrow(events)), ylim = c(0.4, length(row_keys) + 0.8), axes = FALSE, xlab = "", ylab = "")
  title("Secondary callouts: PAR1 and chr9q/chr3q candidate windows", cex.main = 0.86, line = 1.0)
  for (e in seq_len(nrow(events))) {
    rect(e - 1, 0.4, e, length(row_keys) + 0.8, col = ifelse(e %% 2 == 0, "#FAFAFA", "#FFFFFF"), border = "#E2E2E2", lwd = 0.35)
    text(e - 0.5, length(row_keys) + 0.62, events$callout_label[e], cex = 0.52, font = 2)
    text(e - 0.5, length(row_keys) + 0.38, paste0(events$query_chrom[e], ":", format(as.numeric(events$callout_window_start[e]), big.mark = ","), "-", format(as.numeric(events$callout_window_end[e]), big.mark = ",")), cex = 0.42)
  }
  for (mid in row_keys) {
    y <- row_y[mid]
    text(-0.05, y, short_method_label(mid), adj = 1, cex = 0.43)
    for (e in seq_len(nrow(events))) {
      rr <- call[call$method_id == mid & call$callout_event_id == events$callout_event_id[e], ]
      if (nrow(rr) == 0) {
        rect(e - 0.92, y - 0.18, e - 0.08, y + 0.18, col = target_cols[["no_support"]], border = "#D0D0D0", lwd = 0.25)
        next
      }
      rr <- rr[order(rr$query_bin_start), ]
      x0 <- e - 0.92
      width <- 0.84 / nrow(rr)
      for (j in seq_len(nrow(rr))) {
        rect(x0 + (j - 1) * width, y - 0.18, x0 + j * width, y + 0.18, col = row_fill(rr[j, ]), border = NA)
      }
      inter_rr <- rr[rr$top_interchrom_target_chrom != "none" & rr$top_interchrom_support_bp > 0, ]
      if (nrow(inter_rr) > 0) {
        top <- aggregate(top_interchrom_support_bp ~ top_interchrom_target_chrom, inter_rr, sum)
        top <- top[order(-top$top_interchrom_support_bp), ]
        label <- paste(head(top$top_interchrom_target_chrom, 2), collapse = "/")
      } else {
        top <- aggregate(winner_support_bp ~ winner_target_chrom, rr, sum)
        top <- top[order(-top$winner_support_bp), ]
        label <- paste(head(top$winner_target_chrom, 2), collapse = "/")
      }
      text(e - 0.5, y + 0.26, label, cex = 0.30, font = 2)
    }
  }
  axis(1, at = seq_len(nrow(events)) - 0.5, labels = events$callout_event_id, las = 2, tick = FALSE, cex.axis = 0.42)
  mtext("Callouts reuse current query-grid genomic coordinates and remain subordinate to the whole-genome tracks.", side = 1, line = 2.4, cex = 0.52)
}

draw_all <- function() {
  layout(matrix(c(1, 2, 3), nrow = 3), heights = c(0.57, 0.22, 0.21))
  draw_tracks()
  draw_matrix()
  draw_callouts()
}

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 18, height = 15, onefile = TRUE)
draw_all()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 3600, height = 3000, res = 200)
draw_all()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 18, height = 15)
draw_all()
dev.off()
