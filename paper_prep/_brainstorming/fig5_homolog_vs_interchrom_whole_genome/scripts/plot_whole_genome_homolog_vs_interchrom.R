#!/usr/bin/env Rscript

out_dir <- "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome"
tract_path <- "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_top_tracts.tsv"
overall_path <- "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_overall.tsv"
fai_path <- "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa.fai"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tracts <- read.delim(tract_path, stringsAsFactors = FALSE, check.names = FALSE)
overall <- read.delim(overall_path, stringsAsFactors = FALSE, check.names = FALSE)
fai <- read.delim(fai_path, header = FALSE, stringsAsFactors = FALSE)
names(fai)[1:2] <- c("seq", "length")

chrom_from_name <- function(x) {
  sub(".*(chr([0-9]+|X|Y|M)).*", "\\1", x)
}

fai$chrom <- chrom_from_name(fai$seq)
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
fai <- fai[match(chrom_order, fai$chrom), ]
fai <- fai[!is.na(fai$chrom), ]
fai$y <- rev(seq_len(nrow(fai)))
length_by_chrom <- setNames(fai$length, fai$chrom)
y_by_chrom <- setNames(fai$y, fai$chrom)

acro <- c("chr13", "chr14", "chr15", "chr21", "chr22")
is_acro <- function(x) x %in% acro

classify_event <- function(query_chrom, target_chrom) {
  if ((query_chrom == "chrX" && target_chrom == "chrY") || (query_chrom == "chrY" && target_chrom == "chrX")) {
    return("PAR_XY")
  }
  if (query_chrom == "chr9" && target_chrom == "chr3") {
    return("chr9_chr3_candidate")
  }
  if (is_acro(query_chrom) && is_acro(target_chrom)) {
    return("acro_acro")
  }
  if (is_acro(query_chrom) || is_acro(target_chrom)) {
    return("acro_other")
  }
  "other_nonacro"
}

tracts$category <- mapply(classify_event, tracts$query_chrom, tracts$target_chrom)
tracts$mid_mb <- ((tracts$query_start + tracts$query_end) / 2) / 1e6
tracts$start_mb <- tracts$query_start / 1e6
tracts$end_mb <- tracts$query_end / 1e6
tracts$label <- sprintf("%s->%s %.0f kb", tracts$query_chrom, tracts$target_chrom, tracts$bp / 1000)

target_chroms <- chrom_order
target_cols <- c(
  chr1 = "#4E79A7", chr2 = "#A0CBE8", chr3 = "#D95F02", chr4 = "#FFBE7D",
  chr5 = "#59A14F", chr6 = "#8CD17D", chr7 = "#B6992D", chr8 = "#F1CE63",
  chr9 = "#1B9E77", chr10 = "#86BCB6", chr11 = "#E15759", chr12 = "#FF9D9A",
  chr13 = "#79706E", chr14 = "#BAB0AC", chr15 = "#9C755F", chr16 = "#D7B5A6",
  chr17 = "#B07AA1", chr18 = "#D4A6C8", chr19 = "#2F4B7C", chr20 = "#665191",
  chr21 = "#A05195", chr22 = "#D45087", chrX = "#7570B3", chrY = "#E7298A"
)

category_cols <- c(
  PAR_XY = "#E7298A",
  acro_acro = "#616161",
  acro_other = "#9E9E9E",
  chr9_chr3_candidate = "#D95F02",
  other_nonacro = "#2C7FB8"
)

fmt_mb <- function(x) sprintf("%.1f", x)
fmt_kb <- function(x) sprintf("%.0f kb", x / 1000)

summarize_categories <- function(rows) {
  if (nrow(rows) == 0) {
    return(data.frame(category = character(), bp = numeric(), stringsAsFactors = FALSE))
  }
  agg <- aggregate(bp ~ category, rows, sum)
  agg <- agg[order(-agg$bp), ]
  agg
}

draw_tracks <- function(basis) {
  rows <- tracts[tracts$basis == basis, ]
  rows <- rows[rows$query_chrom %in% fai$chrom, ]
  rows$y <- y_by_chrom[rows$query_chrom]
  rows$col <- target_cols[rows$target_chrom]
  rows$col[is.na(rows$col)] <- "#2C7FB8"
  rows$col[rows$category == "acro_acro"] <- "#616161"
  rows$col[rows$category == "acro_other"] <- "#9E9E9E"

  layout(matrix(c(1, 2), nrow = 1), widths = c(4.8, 1.5))
  par(mar = c(4.4, 5.3, 3.6, 0.7), xaxs = "i", yaxs = "i")
  x_max <- max(fai$length) / 1e6 * 1.02
  plot(NA, xlim = c(0, x_max), ylim = c(0.3, nrow(fai) + 0.8), axes = FALSE, xlab = "", ylab = "")
  title(sprintf("Whole-genome homolog-vs-interchrom winners (%s)", basis), line = 2.0, cex.main = 1.05)
  mtext("PAN027 paternal hap2 query; marks show exact tract midpoint where best interchromosomal IMPG similarity beats best same-chromosome/homolog match", side = 3, line = 0.55, cex = 0.62)

  axis(1, at = seq(0, floor(x_max / 50) * 50, by = 50), labels = paste0(seq(0, floor(x_max / 50) * 50, by = 50), " Mb"), cex.axis = 0.65)
  axis(2, at = fai$y, labels = fai$chrom, las = 2, tick = FALSE, cex.axis = 0.70)
  abline(v = seq(0, floor(x_max / 50) * 50, by = 50), col = "#EFEFEF", lwd = 0.5)

  for (i in seq_len(nrow(fai))) {
    segments(0, fai$y[i], fai$length[i] / 1e6, fai$y[i], col = "#D8D8D8", lwd = 4, lend = "butt")
  }

  if (nrow(rows) > 0) {
    ord <- order(rows$category != "acro_acro", rows$bp)
    for (i in ord) {
      rect(rows$start_mb[i], rows$y[i] - 0.10, rows$end_mb[i], rows$y[i] + 0.10, col = adjustcolor(rows$col[i], 0.70), border = NA)
      segments(rows$mid_mb[i], rows$y[i] - 0.32, rows$mid_mb[i], rows$y[i] + 0.32, col = rows$col[i], lwd = 1.1)
      points(rows$mid_mb[i], rows$y[i], pch = 21, bg = rows$col[i], col = "#222222", cex = ifelse(rows$bp[i] >= 20000, 0.55, 0.38), lwd = 0.25)
    }

    label_rows <- rows[
      (rows$category == "PAR_XY" & rows$bp >= 20000) |
        (rows$category == "chr9_chr3_candidate" & rows$bp >= 10000) |
        (rows$category == "other_nonacro" & rows$bp >= 20000) |
        (rows$category == "acro_other" & rows$bp >= 30000),
    ]
    non_control <- rows[rows$category == "other_nonacro" & rows$bp >= 20000, ]
    label_rows <- unique(rbind(label_rows, non_control))
    label_rows <- label_rows[order(label_rows$query_chrom, label_rows$mid_mb, -label_rows$bp), ]
    for (i in seq_len(nrow(label_rows))) {
      dy <- ifelse(i %% 2 == 0, 0.50, -0.56)
      text_y <- pmax(0.7, pmin(nrow(fai) + 0.45, label_rows$y[i] + dy))
      segments(label_rows$mid_mb[i], label_rows$y[i] + sign(dy) * 0.34, label_rows$mid_mb[i], text_y - sign(dy) * 0.10, col = "#777777", lwd = 0.35)
      text(label_rows$mid_mb[i] + 1.1, text_y, sprintf("%s %.0f kb", label_rows$target_chrom[i], label_rows$bp[i] / 1000), adj = 0, cex = 0.45, col = "#222222")
    }
  }

  mtext("Query coordinate, Mb. Chromosome bar lengths are actual query chromosome lengths.", side = 1, line = 3.1, cex = 0.70)

  par(mar = c(4.4, 0.4, 3.6, 1.1), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  title("Summary", line = 2.0, cex.main = 0.95)
  basis_overall <- overall[overall$basis == basis, ]
  y <- 0.96
  if (nrow(basis_overall) == 1) {
    lines <- c(
      sprintf("windows: %s", format(basis_overall$windows_total, big.mark = ",")),
      sprintf("inter > same: %s", format(basis_overall$inter_beats_same_windows, big.mark = ",")),
      sprintf("bp: %.2f Mb", basis_overall$inter_beats_same_bp / 1e6)
    )
    for (line in lines) {
      text(0, y, line, adj = 0, cex = 0.68)
      y <- y - 0.055
    }
  }
  y <- y - 0.025
  text(0, y, "Categories", adj = 0, font = 2, cex = 0.70)
  y <- y - 0.055
  cat_summary <- summarize_categories(rows)
  if (nrow(cat_summary) > 0) {
    for (i in seq_len(nrow(cat_summary))) {
      cat <- cat_summary$category[i]
      points(0.02, y, pch = 15, col = category_cols[[cat]], cex = 0.75)
      text(0.07, y, sprintf("%s  %.0f kb", gsub("_", " ", cat), cat_summary$bp[i] / 1000), adj = 0, cex = 0.58)
      y <- y - 0.050
    }
  }
  y <- y - 0.035
  text(0, y, "Top non-PAR/non-acro-acro tracts", adj = 0, font = 2, cex = 0.70)
  y <- y - 0.055
  top_other <- rows[!(rows$category %in% c("acro_acro", "PAR_XY")), ]
  top_other <- top_other[order(-top_other$bp), ]
  if (nrow(top_other) > 0) {
    for (i in seq_len(min(6, nrow(top_other)))) {
      row <- top_other[i, ]
      points(0.02, y, pch = 21, bg = row$col, col = "#222222", cex = 0.65, lwd = 0.25)
      text(
        0.07,
        y,
        sprintf("%s:%s-%s -> %s  %s",
                row$query_chrom,
                fmt_mb(row$query_start / 1e6),
                fmt_mb(row$query_end / 1e6),
                row$target_chrom,
                fmt_kb(row$bp)),
        adj = 0,
        cex = 0.50
      )
      y <- y - 0.050
    }
  }
  text(0, 0.06, "Color encodes winning target chromosome; gray marks acrocentric targets.", adj = 0, cex = 0.50, col = "#555555")
}

pdf(file.path(out_dir, "fig5_homolog_vs_interchrom_whole_genome.pdf"), width = 14.5, height = 9.2, useDingbats = FALSE)
for (basis in c("10:10", "4:4", "1:1")) {
  draw_tracks(basis)
}
dev.off()

png(file.path(out_dir, "fig5_homolog_vs_interchrom_whole_genome.10to10.png"), width = 2900, height = 1840, res = 200, type = "cairo")
draw_tracks("10:10")
dev.off()

svg(file.path(out_dir, "fig5_homolog_vs_interchrom_whole_genome.10to10.svg"), width = 14.5, height = 9.2, onefile = TRUE)
draw_tracks("10:10")
dev.off()
