#!/usr/bin/env Rscript

# Manuscript Fig 1b -- PHR length distribution by chromosome end, as a per-end
# heatstrip over the terminal 500 kbp in 10 kbp columns, q arms flipped so each
# row reads p telomere -> q telomere. Color = number of PHRs per 10 kbp bin, on
# a pseudo-log scale.
# Location-aware: reads data/ at the repo root, writes submission/fig/MainFigures.
# Ported verbatim (plot code) from the frozen slide generator
# slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

args <- commandArgs(trailingOnly = TRUE)
.cmd_args  <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
script_dir <- if (length(.this_file)) normalizePath(dirname(.this_file)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."))

data_path <- if (length(args) >= 1) args[[1]] else
  Sys.getenv("PHR_LENGTH_TSV", file.path(repo_root, "data/all-vs-all.1Mb.p95.id95.len.tsv"))
out_dir <- if (length(args) >= 2) args[[2]] else
  file.path(repo_root, "submission/fig/MainFigures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

search_window_kbp <- 500
binwidth_kbp <- 10
q_panel_gap_kbp <- 82
q_panel_start_kbp <- search_window_kbp + q_panel_gap_kbp
axis_end_kbp <- q_panel_start_kbp + search_window_kbp

fmt_int <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0(chroms, "p"), paste0(chroms, "q")))

read_lengths <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
  required <- c("seq", "region_start", "region_end")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
  df$arm_label <- sub(".*_chr([0-9XYM]+)_([pq])arm.*", "\\1\\2", df$seq)
  if (any(df$arm_label == df$seq)) stop("Could not parse arm labels from seq")
  valid <- df[df$region_start != ".", ]
  valid$region_start_int <- as.integer(valid$region_start)
  valid$region_end_int <- as.integer(valid$region_end)
  valid$length_kbp <- (valid$region_end_int - valid$region_start_int) / 1000
  valid$chrom <- sub("[pq]$", "", valid$arm_label)
  valid$arm <- sub("^.*([pq])$", "\\1", valid$arm_label)
  if (any(is.na(valid$length_kbp))) stop("NA interval lengths after parsing")
  if (any(valid$length_kbp < 0)) stop("Negative interval lengths")
  if (any(valid$length_kbp > search_window_kbp)) stop("Lengths above the ", search_window_kbp, " kbp window")
  valid
}

make_bin_labels <- function(binwidth) {
  bin_breaks <- seq(0, search_window_kbp, by = binwidth)
  bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")
  bin_labels[length(bin_labels)] <- paste0("[", search_window_kbp - binwidth, ",", search_window_kbp, "]")
  list(breaks = bin_breaks, labels = bin_labels)
}

assign_bins <- function(length_kbp, binwidth) {
  bins <- make_bin_labels(binwidth)
  length_bin <- cut(length_kbp, breaks = bins$breaks, include.lowest = TRUE,
                    right = FALSE, labels = bins$labels)
  length_bin[length_kbp == search_window_kbp] <- bins$labels[length(bins$labels)]
  factor(length_bin, levels = bins$labels)
}

make_arm_bins <- function(valid, binwidth) {
  bins <- make_bin_labels(binwidth)
  length_bin <- assign_bins(valid$length_kbp, binwidth)
  counts <- as.data.frame(
    table(arm_label = factor(valid$arm_label, levels = all_arms), bin = length_bin),
    stringsAsFactors = FALSE)
  names(counts)[names(counts) == "Freq"] <- "count"
  counts$bin_start_kbp <- rep(head(bins$breaks, -1), each = length(all_arms))
  counts$bin_end_kbp <- rep(tail(bins$breaks, -1), each = length(all_arms))
  counts$bin_mid_kbp <- (counts$bin_start_kbp + counts$bin_end_kbp) / 2
  counts$is_500_kbp_edge_bin <- counts$bin_end_kbp == search_window_kbp
  counts$chrom <- sub("[pq]$", "", counts$arm_label)
  counts$arm <- sub("^.*([pq])$", "\\1", counts$arm_label)
  counts$chrom_index <- match(counts$chrom, chroms)
  counts$panel <- ifelse(counts$arm == "p", "p arms", "q arms")
  arm_n <- as.integer(table(factor(valid$arm_label, levels = all_arms)))
  counts$arm_n <- arm_n[match(counts$arm_label, all_arms)]
  counts$within_arm_pct <- ifelse(counts$arm_n > 0, 100 * counts$count / counts$arm_n, NA_real_)
  counts$display_x_kbp <- ifelse(counts$arm == "q",
    q_panel_start_kbp + (search_window_kbp - counts$bin_mid_kbp), counts$bin_mid_kbp)
  counts
}

valid <- read_lengths(data_path)
arm_counts <- table(factor(valid$arm_label, levels = all_arms))
zero_signal_arms <- names(arm_counts)[as.integer(arm_counts) == 0]
arm_bins_10kbp <- make_arm_bins(valid, binwidth_kbp)

zero_labels <- unique(arm_bins_10kbp[arm_bins_10kbp$arm_n == 0, c("chrom", "arm", "chrom_index", "arm_label")])
zero_labels$label <- "No PHRs detected"
zero_labels$xmin <- ifelse(zero_labels$arm == "q", q_panel_start_kbp, 0)
zero_labels$xmax <- ifelse(zero_labels$arm == "q", axis_end_kbp, search_window_kbp)
zero_labels$ymin <- zero_labels$chrom_index - 0.41
zero_labels$ymax <- zero_labels$chrom_index + 0.41
zero_labels$x <- (zero_labels$xmin + zero_labels$xmax) / 2
zero_labels$hjust <- 0.5

panel_labels <- data.frame(
  x = c(search_window_kbp / 2, q_panel_start_kbp + search_window_kbp / 2),
  y = c(-1.15, -1.15), label = c("p arms", "q arms"), stringsAsFactors = FALSE)

x_breaks_p <- seq(0, search_window_kbp, by = 100)
x_breaks_q <- q_panel_start_kbp + seq(0, search_window_kbp, by = 100)
x_breaks <- c(x_breaks_p, x_breaks_q)
x_labels <- c(as.character(x_breaks_p), as.character(rev(x_breaks_p)))

theme_heatstrip_slide <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 8, colour = "#1A3A6B", margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 1, colour = "#333333", lineheight = 0.96, margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 5, colour = "#666666", hjust = 0, margin = margin(t = 8)),
      axis.title = element_text(size = base_size, colour = "#222222"),
      axis.text.x = element_text(size = base_size - 3, colour = "#222222"),
      axis.text.y = element_text(size = base_size - 3, colour = "#222222"),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 3, colour = "#222222"),
      legend.text = element_text(size = base_size - 4, colour = "#222222"),
      plot.margin = margin(26, 80, 10, 80))
}

p_heatstrip <- ggplot(arm_bins_10kbp, aes(x = display_x_kbp, y = chrom_index, fill = count)) +
  annotate("rect", xmin = search_window_kbp, xmax = q_panel_start_kbp, ymin = -Inf, ymax = Inf, fill = "white") +
  geom_tile(width = binwidth_kbp * 0.96, height = 0.82, colour = "#F0F2F4", linewidth = 0.06) +
  geom_rect(data = zero_labels, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#F4B6B0", colour = NA) +
  geom_vline(xintercept = c(search_window_kbp, q_panel_start_kbp), colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(0, axis_end_kbp), colour = "#CCD4DC", linewidth = 0.35) +
  geom_text(data = zero_labels, aes(x = x, y = chrom_index, label = label, hjust = hjust),
            inherit.aes = FALSE, size = 4.1, colour = "#7B241C") +
  geom_text(data = panel_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 8.2, colour = "#222222") +
  annotate("text", x = c(-14, axis_end_kbp + 14), y = length(chroms) + 1.55,
           label = c("← p telomere", "q telomere →"), hjust = c(1, 0), vjust = 0.5,
           size = 4.1, colour = "#555555") +
  scale_fill_gradientn("Number of\nPHRs",
    colours = c("#FFFFFF", "#DBEAF6", "#9ECAE1", "#4292C6", "#08519C"), na.value = "#F3F3F3",
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 20, 100, 400)) +
  scale_x_continuous("PHR length (kbp)",
    limits = c(0, axis_end_kbp), breaks = x_breaks, labels = x_labels,
    minor_breaks = NULL, oob = scales::oob_keep, expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_reverse("Chromosome", breaks = seq_along(chroms), labels = chroms,
    limits = c(length(chroms) + 0.78, -0.72), oob = scales::oob_keep,
    expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_colourbar(barwidth = unit(3.0, "in"), barheight = unit(0.13, "in"))) +
  coord_cartesian(clip = "off") +
  theme_heatstrip_slide()

ggsave(file.path(out_dir, "Fig1b_lengths.png"), p_heatstrip, width = 12.8, height = 7.2, dpi = 240, bg = "white")
ggsave(file.path(out_dir, "Fig1b_lengths.pdf"), p_heatstrip, width = 12.8, height = 7.2, bg = "white")

cat("wrote ", file.path(out_dir, "Fig1b_lengths.pdf"), "\n", sep = "")
cat(fmt_int(sum(arm_counts > 0)), "/48 ends with called intervals; zero-signal: ",
    paste(zero_signal_arms, collapse = ", "), "\n", sep = "")
