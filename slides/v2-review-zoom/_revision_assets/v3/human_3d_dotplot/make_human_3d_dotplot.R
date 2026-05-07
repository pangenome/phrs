#!/usr/bin/env Rscript

out_dir <- "slides/v2-review-zoom/_revision_assets/v3/human_3d_dotplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

human_root <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human"

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE)
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("<2.2e-16")
  formatC(p, format = "e", digits = 2)
}

format_rho <- function(x) {
  if (is.na(x)) return("NA")
  formatC(x, format = "f", digits = 3)
}

chrom_from_arm <- function(x) sub("_[pq]$", "", x)

summarize_xy <- function(d, x_col, y_col, unit, dataset, sample, platform,
                         resolution, source_scope, source_path,
                         x_label, y_label, note) {
  d <- d[!is.na(d[[x_col]]) & !is.na(d[[y_col]]), ]
  ct <- suppressWarnings(stats::cor.test(d[[x_col]], d[[y_col]],
                                         method = "spearman", exact = FALSE))
  data.frame(
    dataset = dataset,
    sample = sample,
    platform = platform,
    resolution = resolution,
    unit = unit,
    source_scope = source_scope,
    n_points = nrow(d),
    n_nonzero_similarity = sum(d[[x_col]] > 0, na.rm = TRUE),
    n_nonzero_contact = sum(d[[y_col]] > 0, na.rm = TRUE),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value,
    x_column = x_col,
    y_column = y_col,
    x_label = x_label,
    y_label = y_label,
    source_path = source_path,
    note = note,
    stringsAsFactors = FALSE
  )
}

load_arm_pair <- function(path) {
  d <- read_tsv(path)
  d$chr_a <- chrom_from_arm(d$arm_a)
  d$chr_b <- chrom_from_arm(d$arm_b)
  d <- d[d$chr_a != d$chr_b, ]
  d
}

load_seq_pair <- function(path) {
  d <- read_tsv(path)
  d <- d[d$chr_a != d$chr_b, ]
  d
}

arm_sources <- list(
  list(
    dataset = "HG002 Pore-C 50 kb",
    sample = "HG002",
    platform = "Pore-C",
    resolution = "50 kb",
    path = file.path(human_root, "community_based/50000bp",
                     "hg002_porec_phr_pair_correlation.tsv"),
    color = "#007a78"
  ),
  list(
    dataset = "CHM13 Hi-C 50 kb",
    sample = "CHM13",
    platform = "Hi-C",
    resolution = "50 kb",
    path = file.path(human_root, "community_based/50000bp",
                     "chm13_phr_pair_correlation.tsv"),
    color = "#b64b2a"
  )
)

seq_sources <- list(
  list(
    dataset = "HG002 Pore-C 50 kb, sequence-level control",
    sample = "HG002",
    platform = "Pore-C",
    resolution = "50 kb",
    path = file.path(human_root, "exclusion_controls/no_strong/community_free",
                     "human_HG002_porec_hic_50000bp_seqlevel.tsv"),
    color = "#007a78"
  ),
  list(
    dataset = "CHM13 Hi-C 50 kb, sequence-level control",
    sample = "CHM13",
    platform = "Hi-C",
    resolution = "50 kb",
    path = file.path(human_root, "exclusion_controls/no_strong/community_free",
                     "human_CHM13_hic_50000bp_seqlevel.tsv"),
    color = "#b64b2a"
  )
)

arm_data <- lapply(arm_sources, function(src) {
  d <- load_arm_pair(src$path)
  list(source = src, data = d,
       summary = summarize_xy(
         d = d,
         x_col = "mean_jaccard",
         y_col = "hic_contact",
         unit = "arm pair",
         dataset = src$dataset,
         sample = src$sample,
         platform = src$platform,
         resolution = src$resolution,
         source_scope = "main community_based full arm set",
         source_path = src$path,
         x_label = "Mean PHR Jaccard similarity per chromosome-arm pair",
         y_label = "3D contact frequency",
         note = "Pointwise Spearman across arm pairs; not a Mantel test."
       ))
})

seq_data <- lapply(seq_sources, function(src) {
  d <- load_seq_pair(src$path)
  list(source = src, data = d,
       summary = summarize_xy(
         d = d,
         x_col = "jaccard",
         y_col = "hic_contact_norm",
         unit = "PHR sequence pair",
         dataset = src$dataset,
         sample = src$sample,
         platform = src$platform,
         resolution = src$resolution,
         source_scope = "community_free no_strong exclusion control",
         source_path = src$path,
         x_label = "PHR sequence-pair Jaccard similarity",
         y_label = "Balanced contact per valid bin pair",
         note = paste(
           "Pointwise Spearman across sequence pairs; not a Mantel test.",
           "Control excludes chr4_q, chr10_q, acrocentric p arms, and sex arms."
         )
       ))
})

summary_table <- do.call(rbind, c(lapply(arm_data, `[[`, "summary"),
                                  lapply(seq_data, `[[`, "summary")))
summary_path <- file.path(out_dir, "human_3d_dotplot_summary.tsv")
write.table(summary_table, summary_path, sep = "\t", row.names = FALSE,
            quote = FALSE)

plot_one <- function(item, x_col, y_col, x_label, y_label, display_unit) {
  d <- item$data
  src <- item$source
  summ <- item$summary
  x <- d[[x_col]]
  y <- d[[y_col]]
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]

  positive <- y[y > 0]
  floor_y <- if (length(positive)) min(positive) / 2 else 1e-12
  y_plot <- pmax(y, floor_y)

  plot(x, y_plot, log = "y", pch = 21,
       bg = grDevices::adjustcolor(src$color, alpha.f = 0.42),
       col = grDevices::adjustcolor("#1f1f1f", alpha.f = 0.32),
       lwd = 0.35, cex = 0.88,
       xlab = x_label,
       ylab = y_label,
       main = src$dataset,
       cex.main = 1.02, cex.lab = 0.92, cex.axis = 0.80,
       xaxs = "i")
  grid(col = "#e6e6e6", lwd = 0.7)

  if (length(unique(x)) > 2 && length(unique(y_plot)) > 2) {
    fit <- stats::lm(log10(y_plot) ~ x)
    xs <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
    lines(xs, 10 ^ stats::predict(fit, newdata = data.frame(x = xs)),
          col = "#111111", lwd = 1.35)
  }

  legend("topleft",
         legend = c(
           sprintf("n = %s %s", format(summ$n_points, big.mark = ","),
                   display_unit),
           sprintf("pointwise Spearman rho = %s",
                   format_rho(summ$spearman_rho)),
           sprintf("p = %s", format_p(summ$p_value)),
           sprintf("zero contacts at floor: %.1e", floor_y)
         ),
         bty = "n", cex = 0.72, text.col = "#222222")
}

render_plot_set <- function(items, png_name, pdf_name, title, subtitle,
                            x_col, y_col, x_label, y_label, display_unit) {
  png(file.path(out_dir, png_name), width = 1800, height = 1200, res = 180,
      type = "cairo")
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), mar = c(5.1, 5.1, 5.4, 1.2),
      oma = c(1.1, 0, 4.1, 0), family = "sans")
  for (item in items) {
    plot_one(item, x_col, y_col, x_label, y_label, display_unit)
  }
  mtext(title, outer = TRUE, side = 3, line = 2.3, font = 2, cex = 1.18)
  mtext(subtitle, outer = TRUE, side = 3, line = 1.0, cex = 0.78,
        col = "#444444")
  par(op)
  dev.off()

  grDevices::pdf(file.path(out_dir, pdf_name), width = 10, height = 6.67,
                 family = "Helvetica")
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), mar = c(5.1, 5.1, 5.4, 1.2),
      oma = c(1.1, 0, 4.1, 0), family = "sans")
  for (item in items) {
    plot_one(item, x_col, y_col, x_label, y_label, display_unit)
  }
  mtext(title, outer = TRUE, side = 3, line = 2.3, font = 2, cex = 1.18)
  mtext(subtitle, outer = TRUE, side = 3, line = 1.0, cex = 0.78,
        col = "#444444")
  par(op)
  dev.off()
}

render_plot_set(
  arm_data,
  png_name = "human_arm_pair_dotplot_candidate.png",
  pdf_name = "human_arm_pair_dotplot_candidate.pdf",
  title = "Human sequence similarity versus 3D contact at 50 kb",
  subtitle = paste(
    "Main Fig. 3/ED5 analog: each dot is an inter-chromosomal arm pair.",
    "Statistics are pointwise Spearman correlations, not Mantel tests."
  ),
  x_col = "mean_jaccard",
  y_col = "hic_contact",
  x_label = "Mean PHR Jaccard similarity per arm pair",
  y_label = "3D contact frequency (log scale; zero shown at floor)",
  display_unit = "arm pairs"
)

render_plot_set(
  seq_data,
  png_name = "human_phr_sequence_pair_control.png",
  pdf_name = "human_phr_sequence_pair_control.pdf",
  title = "Human PHR sequence-pair control versus 3D contact at 50 kb",
  subtitle = paste(
    "Community-free no_strong control: each dot is one PHR sequence pair.",
    "Use only with the exclusion label; this is not the main arm-pair/Mantel statistic."
  ),
  x_col = "jaccard",
  y_col = "hic_contact_norm",
  x_label = "PHR sequence-pair Jaccard similarity",
  y_label = "Balanced contact per valid bin pair (log scale; zero shown at floor)",
  display_unit = "PHR pairs"
)

message("Wrote human 3D dotplot assets to ", out_dir)
message("Summary: ", summary_path)
