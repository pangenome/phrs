#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

data_path <- Sys.getenv(
  "PHR_LENGTH_TSV",
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
)

search_cap_kb <- 500
current_v7_binwidth_kb <- 25
recommended_binwidth_kb <- 10
sensitivity_binwidth_kb <- 5

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

out_dir <- script_dir()

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_num <- function(x, digits = 1) {
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

q_value <- function(x, p) {
  as.numeric(quantile(x, p, names = FALSE, type = 7))
}

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0(chroms, "p"), paste0(chroms, "q")))
chrom_factor_levels <- rev(chroms)

read_lengths <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
  required <- c("seq", "region_start", "region_end")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  df$arm_label <- sub(".*_chr([0-9XYM]+)_([pq])arm.*", "\\1\\2", df$seq)
  bad_arm <- df$arm_label == df$seq
  if (any(bad_arm)) {
    stop("Could not parse arm labels from seq for ", sum(bad_arm), " rows")
  }

  valid <- df[df$region_start != ".", ]
  valid$region_start_int <- as.integer(valid$region_start)
  valid$region_end_int <- as.integer(valid$region_end)
  valid$length_kb <- (valid$region_end_int - valid$region_start_int) / 1000
  valid$chrom <- sub("[pq]$", "", valid$arm_label)
  valid$arm <- sub("^.*([pq])$", "\\1", valid$arm_label)

  if (any(is.na(valid$length_kb))) {
    stop("Encountered NA interval lengths after parsing region_start/region_end")
  }
  if (any(valid$length_kb < 0)) {
    stop("Encountered negative interval lengths")
  }
  if (any(valid$length_kb > search_cap_kb)) {
    stop("Found reported lengths above the expected ", search_cap_kb, " kb cap")
  }

  list(raw = df, valid = valid)
}

make_bin_labels <- function(binwidth_kb) {
  bin_breaks <- seq(0, search_cap_kb, by = binwidth_kb)
  bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")
  bin_labels[length(bin_labels)] <- paste0("[", search_cap_kb - binwidth_kb, ",", search_cap_kb, "]")
  list(breaks = bin_breaks, labels = bin_labels)
}

assign_bins <- function(length_kb, binwidth_kb) {
  bins <- make_bin_labels(binwidth_kb)
  length_bin <- cut(
    length_kb,
    breaks = bins$breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = bins$labels
  )
  length_bin[length_kb == search_cap_kb] <- bins$labels[length(bins$labels)]
  factor(length_bin, levels = bins$labels)
}

make_histogram_bins <- function(valid, binwidth_kb) {
  bins <- make_bin_labels(binwidth_kb)
  length_bin <- assign_bins(valid$length_kb, binwidth_kb)
  bin_summary <- as.data.frame(table(length_bin), stringsAsFactors = FALSE)
  names(bin_summary) <- c("bin", "count")
  bin_summary$bin_start_kb <- head(bins$breaks, -1)
  bin_summary$bin_end_kb <- tail(bins$breaks, -1)
  bin_summary$bin_mid_kb <- (bin_summary$bin_start_kb + bin_summary$bin_end_kb) / 2
  bin_summary$is_cap_bin <- bin_summary$bin_end_kb == search_cap_kb
  bin_summary$fill <- ifelse(bin_summary$is_cap_bin, "cap bin", "measured bins")
  bin_summary
}

make_arm_bins <- function(valid, binwidth_kb) {
  bins <- make_bin_labels(binwidth_kb)
  length_bin <- assign_bins(valid$length_kb, binwidth_kb)
  counts <- as.data.frame(
    table(
      arm_label = factor(valid$arm_label, levels = all_arms),
      bin = length_bin
    ),
    stringsAsFactors = FALSE
  )
  names(counts)[names(counts) == "Freq"] <- "count"

  counts$bin_start_kb <- rep(head(bins$breaks, -1), each = length(all_arms))
  counts$bin_end_kb <- rep(tail(bins$breaks, -1), each = length(all_arms))
  counts$bin_mid_kb <- (counts$bin_start_kb + counts$bin_end_kb) / 2
  counts$is_cap_bin <- counts$bin_end_kb == search_cap_kb
  counts$chrom <- sub("[pq]$", "", counts$arm_label)
  counts$arm <- sub("^.*([pq])$", "\\1", counts$arm_label)
  counts$arm_side <- factor(
    counts$arm,
    levels = c("p", "q"),
    labels = c("p arms", "q arms")
  )
  counts$chrom_factor <- factor(counts$chrom, levels = chrom_factor_levels)

  arm_n <- as.integer(table(factor(valid$arm_label, levels = all_arms)))
  counts$arm_n <- arm_n[match(counts$arm_label, all_arms)]
  counts$within_arm_pct <- ifelse(
    counts$arm_n > 0,
    100 * counts$count / counts$arm_n,
    NA_real_
  )
  counts
}

dat <- read_lengths(data_path)
raw <- dat$raw
valid <- dat$valid
lengths <- valid$length_kb

arm_counts <- table(factor(valid$arm_label, levels = all_arms))
zero_signal_arms <- names(arm_counts)[as.integer(arm_counts) == 0]
cap_hits <- sum(lengths == search_cap_kb)
median_kb <- median(lengths)
p90_kb <- q_value(lengths, 0.90)
p95_kb <- q_value(lengths, 0.95)

summary_stats <- data.frame(
  metric = c(
    "source_tsv",
    "current_v7_bin_width_kb",
    "recommended_v8_bin_width_kb",
    "sensitivity_bin_width_kb",
    "search_window_kb",
    "non_empty_intervals",
    "arms_with_called_intervals",
    "total_chromosome_ends_checked",
    "zero_signal_arms",
    "min_kb",
    "p25_kb",
    "median_kb",
    "mean_kb",
    "p75_kb",
    "p90_kb",
    "p95_kb",
    "max_reported_kb",
    "reported_at_500kb_n",
    "reported_at_500kb_pct",
    "analysis_ceiling_note"
  ),
  value = c(
    data_path,
    current_v7_binwidth_kb,
    recommended_binwidth_kb,
    sensitivity_binwidth_kb,
    search_cap_kb,
    fmt_int(length(lengths)),
    fmt_int(sum(arm_counts > 0)),
    fmt_int(length(all_arms)),
    paste(zero_signal_arms, collapse = ", "),
    fmt_num(min(lengths), 0),
    fmt_num(q_value(lengths, 0.25), 0),
    fmt_num(median_kb, 0),
    fmt_num(mean(lengths), 1),
    fmt_num(q_value(lengths, 0.75), 0),
    fmt_num(p90_kb, 0),
    fmt_num(p95_kb, 0),
    fmt_num(max(lengths), 0),
    fmt_int(cap_hits),
    fmt_num(100 * cap_hits / length(lengths), 1),
    "analysis searched/measured terminal 500 kb windows; values at/near 500 kb are right-censored"
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_stats,
  file = file.path(out_dir, "length_distribution_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

arm_summary <- do.call(rbind, lapply(all_arms, function(arm_label) {
  x <- valid$length_kb[valid$arm_label == arm_label]
  data.frame(
    arm_label = arm_label,
    chrom = sub("[pq]$", "", arm_label),
    arm = sub("^.*([pq])$", "\\1", arm_label),
    n = length(x),
    median_kb = if (length(x) == 0) NA_real_ else median(x),
    mean_kb = if (length(x) == 0) NA_real_ else mean(x),
    p25_kb = if (length(x) == 0) NA_real_ else q_value(x, 0.25),
    p75_kb = if (length(x) == 0) NA_real_ else q_value(x, 0.75),
    p90_kb = if (length(x) == 0) NA_real_ else q_value(x, 0.90),
    p95_kb = if (length(x) == 0) NA_real_ else q_value(x, 0.95),
    max_kb = if (length(x) == 0) NA_real_ else max(x),
    reported_at_500kb_n = sum(x == search_cap_kb),
    stringsAsFactors = FALSE
  )
}))

write.table(
  arm_summary,
  file = file.path(out_dir, "arm_length_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

hist_10kb <- make_histogram_bins(valid, recommended_binwidth_kb)
hist_5kb <- make_histogram_bins(valid, sensitivity_binwidth_kb)
arm_bins_10kb <- make_arm_bins(valid, recommended_binwidth_kb)

write.table(
  hist_10kb[, c("bin", "bin_start_kb", "bin_end_kb", "count", "is_cap_bin")],
  file = file.path(out_dir, "histogram_bins_10kb.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hist_5kb[, c("bin", "bin_start_kb", "bin_end_kb", "count", "is_cap_bin")],
  file = file.path(out_dir, "histogram_bins_5kb.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  arm_bins_10kb[, c(
    "arm_label",
    "chrom",
    "arm",
    "bin",
    "bin_start_kb",
    "bin_end_kb",
    "count",
    "arm_n",
    "within_arm_pct",
    "is_cap_bin"
  )],
  file = file.path(out_dir, "arm_length_bins_10kb.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

theme_histogram_slide <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 7,
        colour = "#1A3A6B",
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        colour = "#333333",
        lineheight = 0.96,
        margin = margin(b = 14)
      ),
      plot.caption = element_text(
        size = base_size - 6,
        colour = "#666666",
        hjust = 0,
        margin = margin(t = 10)
      ),
      axis.title = element_text(size = base_size + 1, colour = "#222222"),
      axis.text = element_text(size = base_size - 2, colour = "#222222"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E5E9EF", linewidth = 0.35),
      legend.position = "none",
      plot.margin = margin(14, 24, 12, 22)
    )
}

make_histogram_plot <- function(bin_summary, binwidth_kb, title_suffix, subtitle_context) {
  max_count <- max(bin_summary$count)
  cap_bin_count <- bin_summary$count[bin_summary$is_cap_bin]
  callout_text <- paste(
    "terminal 500 kb window;",
    "edge values are right-censored",
    sep = "\n"
  )
  stats_text <- paste0(
    "n = ", fmt_int(length(lengths)), " called intervals",
    "\nmedian = ", fmt_num(median_kb, 0), " kb",
    "\n90th percentile = ", fmt_num(p90_kb, 0), " kb",
    "\nexact 500 kb calls = ", fmt_int(cap_hits),
    "\n", binwidth_kb, " kb cap-bin count = ", fmt_int(cap_bin_count)
  )

  ggplot(bin_summary, aes(x = bin_mid_kb, y = count, fill = fill)) +
    annotate(
      "rect",
      xmin = search_cap_kb - binwidth_kb,
      xmax = search_cap_kb,
      ymin = 0,
      ymax = Inf,
      fill = "#FFF2CC",
      alpha = 0.42
    ) +
    geom_col(
      width = binwidth_kb * 0.90,
      colour = if (binwidth_kb <= 5) "#F5F8FB" else "white",
      linewidth = if (binwidth_kb <= 5) 0.12 else 0.28
    ) +
    geom_vline(
      xintercept = median_kb,
      linetype = "dashed",
      colour = "#222222",
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = p90_kb,
      linetype = "dotted",
      colour = "#1A3A6B",
      linewidth = 0.62
    ) +
    geom_vline(
      xintercept = search_cap_kb,
      colour = "#D4820A",
      linewidth = 1.25
    ) +
    annotate(
      "label",
      x = search_cap_kb - 5,
      y = max_count * 0.94,
      label = callout_text,
      hjust = 1,
      vjust = 1,
      size = 4.85,
      lineheight = 0.92,
      label.size = 0.26,
      label.r = unit(2.5, "pt"),
      fill = "#FFF8E8",
      colour = "#573700"
    ) +
    annotate(
      "segment",
      x = search_cap_kb - 52,
      xend = search_cap_kb,
      y = max_count * 0.78,
      yend = max_count * 0.78,
      colour = "#D4820A",
      linewidth = 0.7,
      arrow = arrow(length = unit(0.12, "in"), type = "closed")
    ) +
    annotate(
      "label",
      x = 274,
      y = max_count * 0.64,
      label = stats_text,
      hjust = 0,
      vjust = 1,
      size = 4.65,
      lineheight = 0.94,
      label.size = 0.22,
      label.r = unit(2.5, "pt"),
      fill = "white",
      colour = "#222222"
    ) +
    annotate(
      "text",
      x = median_kb + 6,
      y = max_count * 0.48,
      label = paste0("median ", fmt_num(median_kb, 0), " kb"),
      hjust = 0,
      size = 4.25,
      colour = "#222222"
    ) +
    annotate(
      "text",
      x = p90_kb + 5,
      y = max_count * 0.35,
      label = paste0("90th ", fmt_num(p90_kb, 0), " kb"),
      hjust = 0,
      size = 4.0,
      colour = "#1A3A6B"
    ) +
    scale_fill_manual(values = c("measured bins" = "#4C78A8", "cap bin" = "#D4820A")) +
    scale_x_continuous(
      "Called PHR interval length (kb)",
      limits = c(0, search_cap_kb),
      breaks = seq(0, search_cap_kb, by = 50),
      expand = expansion(mult = c(0.01, 0.015))
    ) +
    scale_y_continuous(
      paste0("Intervals per ", binwidth_kb, " kb bin"),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = paste0("Called PHR lengths in the 500 kb window", title_suffix),
      subtitle = paste0(
        subtitle_context, "\n",
        "The right edge is a measurement ceiling, not evidence that sharing stops."
      ),
      caption = paste0(
        "Source: ", basename(data_path),
        " | bin width = ", binwidth_kb, " kb",
        " | ", fmt_int(sum(arm_counts > 0)), "/48 ends have called intervals",
        " | zero-signal ends: ", paste(zero_signal_arms, collapse = ", ")
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_histogram_slide()
}

p_hist_10kb <- make_histogram_plot(
  hist_10kb,
  recommended_binwidth_kb,
  "",
  "All non-empty inter-chromosomal PHR calls; 10 kb bins within the terminal 500 kb discovery window."
)

p_hist_5kb <- make_histogram_plot(
  hist_5kb,
  sensitivity_binwidth_kb,
  " - 5 kb sensitivity",
  "All non-empty inter-chromosomal PHR calls; 5 kb bins match the detection window step."
)

zero_labels <- unique(arm_bins_10kb[arm_bins_10kb$arm_n == 0, c("chrom", "arm", "arm_side", "chrom_factor", "arm_label")])
zero_labels$label <- "no called interval"
ceiling_labels <- data.frame(
  arm_side = factor(c("p arms", "q arms"), levels = c("p arms", "q arms")),
  chrom_factor = factor(c("1", "1"), levels = chrom_factor_levels),
  label = "500 kb ceiling",
  stringsAsFactors = FALSE
)

theme_heatstrip_slide <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 8,
        colour = "#1A3A6B",
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        colour = "#333333",
        lineheight = 0.96,
        margin = margin(b = 10)
      ),
      plot.caption = element_text(
        size = base_size - 5,
        colour = "#666666",
        hjust = 0,
        margin = margin(t = 8)
      ),
      axis.title = element_text(size = base_size, colour = "#222222"),
      axis.text.x = element_text(size = base_size - 3, colour = "#222222"),
      axis.text.y = element_text(size = base_size - 3, colour = "#222222"),
      strip.text = element_text(face = "bold", size = base_size + 1, colour = "#222222"),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.45, "in"),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 3, colour = "#222222"),
      legend.text = element_text(size = base_size - 4, colour = "#222222"),
      plot.margin = margin(12, 22, 10, 20)
    )
}

p_heatstrip <- ggplot(
  arm_bins_10kb,
  aes(x = bin_mid_kb, y = chrom_factor, fill = within_arm_pct)
) +
  geom_tile(width = recommended_binwidth_kb * 0.96, height = 0.82, colour = "#F0F2F4", linewidth = 0.06) +
  geom_tile(
    data = arm_bins_10kb[arm_bins_10kb$is_cap_bin, ],
    aes(x = bin_mid_kb, y = chrom_factor),
    inherit.aes = FALSE,
    fill = NA,
    colour = "#D4820A",
    linewidth = 0.18,
    width = recommended_binwidth_kb * 0.96,
    height = 0.82
  ) +
  geom_vline(xintercept = search_cap_kb, colour = "#D4820A", linewidth = 0.9) +
  geom_text(
    data = zero_labels,
    aes(x = 54, y = chrom_factor, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.35,
    colour = "#7A5A68"
  ) +
  geom_text(
    data = ceiling_labels,
    aes(x = 493, y = chrom_factor, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = -0.45,
    size = 3.0,
    colour = "#6B4300"
  ) +
  facet_grid(. ~ arm_side) +
  scale_fill_gradientn(
    "Within-end\n% of calls",
    colours = c("#FFFFFF", "#DBEAF6", "#9ECAE1", "#4292C6", "#08519C"),
    na.value = "#F3F3F3"
  ) +
  scale_x_continuous(
    "Called PHR interval length (kb), 10 kb columns",
    limits = c(0, search_cap_kb),
    breaks = seq(0, search_cap_kb, by = 100),
    minor_breaks = seq(0, search_cap_kb, by = 50),
    expand = expansion(mult = c(0.01, 0.015))
  ) +
  scale_y_discrete("Chromosome") +
  guides(fill = guide_colourbar(barwidth = unit(3.0, "in"), barheight = unit(0.13, "in"))) +
  labs(
    title = "PHR length distribution by chromosome end",
    subtitle = paste0(
      "Each row is one chromosome end, split into p-arm and q-arm panels. ",
      "Color is the within-end share of called intervals in each 10 kb bin.\n",
      "The analysis searched/measured terminal 500 kb windows; values at/near the right edge are right-censored."
    ),
    caption = paste0(
      "Source: ", basename(data_path),
      " | ", fmt_int(sum(arm_counts > 0)), "/48 ends have called intervals; zero-signal ends: ", paste(zero_signal_arms, collapse = ", "),
      " | bin width = ", recommended_binwidth_kb, " kb"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_heatstrip_slide()

assets <- data.frame(
  asset = c(
    "phr_length_histogram_10kb.png",
    "phr_length_histogram_10kb.pdf",
    "phr_length_histogram_5kb_sensitivity.png",
    "phr_length_histogram_5kb_sensitivity.pdf",
    "phr_length_arm_heatstrip_10kb.png",
    "phr_length_arm_heatstrip_10kb.pdf"
  ),
  bin_width_kb = c(10, 10, 5, 5, 10, 10),
  role = c(
    "recommended_simple_histogram",
    "recommended_simple_histogram",
    "finer_bin_sensitivity",
    "finer_bin_sensitivity",
    "recommended_per_end_old_style",
    "recommended_per_end_old_style"
  ),
  note = c(
    "Talk-fast overall distribution; fixes v7's 25 kb mega-bin issue.",
    "PDF companion for the talk-fast overall distribution.",
    "Shows exact 5 kb step sensitivity; busier than the 10 kb candidate.",
    "PDF companion for 5 kb sensitivity.",
    "All chromosome ends represented as p/q heatstrips with direct chromosome labels.",
    "PDF companion for the per-end heatstrip."
  ),
  stringsAsFactors = FALSE
)

write.table(
  assets,
  file = file.path(out_dir, "asset_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  file.path(out_dir, "phr_length_histogram_10kb.png"),
  p_hist_10kb,
  width = 12.8,
  height = 7.2,
  dpi = 240,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_histogram_10kb.pdf"),
  p_hist_10kb,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_histogram_5kb_sensitivity.png"),
  p_hist_5kb,
  width = 12.8,
  height = 7.2,
  dpi = 240,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_histogram_5kb_sensitivity.pdf"),
  p_hist_5kb,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_arm_heatstrip_10kb.png"),
  p_heatstrip,
  width = 12.8,
  height = 7.2,
  dpi = 240,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_arm_heatstrip_10kb.pdf"),
  p_heatstrip,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

message("Wrote slide 06a v8 length alternatives to: ", out_dir)
