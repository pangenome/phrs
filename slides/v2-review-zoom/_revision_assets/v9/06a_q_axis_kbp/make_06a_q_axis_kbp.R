#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

data_path <- Sys.getenv(
  "PHR_LENGTH_TSV",
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
)

search_window_kbp <- 500
binwidth_kbp <- 10
q_panel_gap_kbp <- 82
q_panel_start_kbp <- search_window_kbp + q_panel_gap_kbp
axis_end_kbp <- q_panel_start_kbp + search_window_kbp

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
  valid$length_kbp <- (valid$region_end_int - valid$region_start_int) / 1000
  valid$chrom <- sub("[pq]$", "", valid$arm_label)
  valid$arm <- sub("^.*([pq])$", "\\1", valid$arm_label)

  if (any(is.na(valid$length_kbp))) {
    stop("Encountered NA interval lengths after parsing region_start/region_end")
  }
  if (any(valid$length_kbp < 0)) {
    stop("Encountered negative interval lengths")
  }
  if (any(valid$length_kbp > search_window_kbp)) {
    stop("Found reported lengths above the expected ", search_window_kbp, " kbp window")
  }

  list(raw = df, valid = valid)
}

make_bin_labels <- function(binwidth) {
  bin_breaks <- seq(0, search_window_kbp, by = binwidth)
  bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")
  bin_labels[length(bin_labels)] <- paste0("[", search_window_kbp - binwidth, ",", search_window_kbp, "]")
  list(breaks = bin_breaks, labels = bin_labels)
}

assign_bins <- function(length_kbp, binwidth) {
  bins <- make_bin_labels(binwidth)
  length_bin <- cut(
    length_kbp,
    breaks = bins$breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = bins$labels
  )
  length_bin[length_kbp == search_window_kbp] <- bins$labels[length(bins$labels)]
  factor(length_bin, levels = bins$labels)
}

make_arm_bins <- function(valid, binwidth) {
  bins <- make_bin_labels(binwidth)
  length_bin <- assign_bins(valid$length_kbp, binwidth)
  counts <- as.data.frame(
    table(
      arm_label = factor(valid$arm_label, levels = all_arms),
      bin = length_bin
    ),
    stringsAsFactors = FALSE
  )
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
  counts$within_arm_pct <- ifelse(
    counts$arm_n > 0,
    100 * counts$count / counts$arm_n,
    NA_real_
  )
  counts$display_x_kbp <- ifelse(
    counts$arm == "q",
    q_panel_start_kbp + (search_window_kbp - counts$bin_mid_kbp),
    counts$bin_mid_kbp
  )
  counts
}

dat <- read_lengths(data_path)
valid <- dat$valid
lengths <- valid$length_kbp

arm_counts <- table(factor(valid$arm_label, levels = all_arms))
zero_signal_arms <- names(arm_counts)[as.integer(arm_counts) == 0]
reported_at_edge_n <- sum(lengths == search_window_kbp)
median_kbp <- median(lengths)
p90_kbp <- q_value(lengths, 0.90)
p95_kbp <- q_value(lengths, 0.95)

summary_stats <- data.frame(
  metric = c(
    "source_tsv",
    "bin_width_kbp",
    "search_window_kbp",
    "q_axis_orientation",
    "non_empty_intervals",
    "arms_with_called_intervals",
    "total_chromosome_ends_checked",
    "zero_signal_arms",
    "min_kbp",
    "p25_kbp",
    "median_kbp",
    "mean_kbp",
    "p75_kbp",
    "p90_kbp",
    "p95_kbp",
    "max_reported_kbp",
    "reported_at_500_kbp_n",
    "reported_at_500_kbp_pct",
    "terminal_window_note"
  ),
  value = c(
    data_path,
    binwidth_kbp,
    search_window_kbp,
    "q panel display is flipped: left edge is 500 kbp, right edge is 0 kbp",
    fmt_int(length(lengths)),
    fmt_int(sum(arm_counts > 0)),
    fmt_int(length(all_arms)),
    paste(zero_signal_arms, collapse = ", "),
    fmt_num(min(lengths), 0),
    fmt_num(q_value(lengths, 0.25), 0),
    fmt_num(median_kbp, 0),
    fmt_num(mean(lengths), 1),
    fmt_num(q_value(lengths, 0.75), 0),
    fmt_num(p90_kbp, 0),
    fmt_num(p95_kbp, 0),
    fmt_num(max(lengths), 0),
    fmt_int(reported_at_edge_n),
    fmt_num(100 * reported_at_edge_n / length(lengths), 1),
    "analysis measured terminal 500 kbp windows; >500 kbp was not measured"
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
  x <- valid$length_kbp[valid$arm_label == arm_label]
  data.frame(
    arm_label = arm_label,
    chrom = sub("[pq]$", "", arm_label),
    arm = sub("^.*([pq])$", "\\1", arm_label),
    n = length(x),
    median_kbp = if (length(x) == 0) NA_real_ else median(x),
    mean_kbp = if (length(x) == 0) NA_real_ else mean(x),
    p25_kbp = if (length(x) == 0) NA_real_ else q_value(x, 0.25),
    p75_kbp = if (length(x) == 0) NA_real_ else q_value(x, 0.75),
    p90_kbp = if (length(x) == 0) NA_real_ else q_value(x, 0.90),
    p95_kbp = if (length(x) == 0) NA_real_ else q_value(x, 0.95),
    max_reported_kbp = if (length(x) == 0) NA_real_ else max(x),
    reported_at_500_kbp_n = sum(x == search_window_kbp),
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

arm_bins_10kbp <- make_arm_bins(valid, binwidth_kbp)

write.table(
  arm_bins_10kbp[, c(
    "arm_label",
    "chrom",
    "arm",
    "panel",
    "bin",
    "bin_start_kbp",
    "bin_end_kbp",
    "bin_mid_kbp",
    "display_x_kbp",
    "count",
    "arm_n",
    "within_arm_pct",
    "is_500_kbp_edge_bin"
  )],
  file = file.path(out_dir, "arm_length_bins_10kbp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

zero_labels <- unique(arm_bins_10kbp[arm_bins_10kbp$arm_n == 0, c("chrom", "arm", "chrom_index", "arm_label")])
zero_labels$label <- "no called interval"
zero_labels$x <- ifelse(
  zero_labels$arm == "q",
  q_panel_start_kbp + search_window_kbp - 54,
  54
)
zero_labels$hjust <- ifelse(zero_labels$arm == "q", 1, 0)

edge_labels <- data.frame(
  x = c(search_window_kbp - 5, q_panel_start_kbp + 5),
  y = c(0.62, 0.62),
  label = ">500 kbp not measured",
  hjust = c(1, 0),
  stringsAsFactors = FALSE
)

panel_labels <- data.frame(
  x = c(search_window_kbp / 2, q_panel_start_kbp + search_window_kbp / 2),
  y = c(-0.34, -0.34),
  label = c("p arms", "q arms"),
  stringsAsFactors = FALSE
)

x_breaks_p <- seq(0, search_window_kbp, by = 100)
x_breaks_q <- q_panel_start_kbp + seq(0, search_window_kbp, by = 100)
x_breaks <- c(x_breaks_p, x_breaks_q)
x_labels <- c(as.character(x_breaks_p), as.character(rev(x_breaks_p)))

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
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 3, colour = "#222222"),
      legend.text = element_text(size = base_size - 4, colour = "#222222"),
      plot.margin = margin(12, 22, 10, 20)
    )
}

p_heatstrip <- ggplot(
  arm_bins_10kbp,
  aes(x = display_x_kbp, y = chrom_index, fill = within_arm_pct)
) +
  annotate(
    "rect",
    xmin = search_window_kbp,
    xmax = q_panel_start_kbp,
    ymin = -Inf,
    ymax = Inf,
    fill = "white"
  ) +
  geom_tile(width = binwidth_kbp * 0.96, height = 0.82, colour = "#F0F2F4", linewidth = 0.06) +
  geom_tile(
    data = arm_bins_10kbp[arm_bins_10kbp$is_500_kbp_edge_bin, ],
    aes(x = display_x_kbp, y = chrom_index),
    inherit.aes = FALSE,
    fill = NA,
    colour = "#D4820A",
    linewidth = 0.18,
    width = binwidth_kbp * 0.96,
    height = 0.82
  ) +
  geom_vline(xintercept = c(search_window_kbp, q_panel_start_kbp), colour = "#D4820A", linewidth = 0.9) +
  geom_vline(xintercept = c(0, axis_end_kbp), colour = "#CCD4DC", linewidth = 0.35) +
  geom_text(
    data = zero_labels,
    aes(x = x, y = chrom_index, label = label, hjust = hjust),
    inherit.aes = FALSE,
    size = 2.35,
    colour = "#7A5A68"
  ) +
  geom_text(
    data = edge_labels,
    aes(x = x, y = y, label = label, hjust = hjust),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.0,
    colour = "#6B4300"
  ) +
  geom_text(
    data = panel_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 8.2,
    colour = "#222222"
  ) +
  annotate(
    "text",
    x = c(0, axis_end_kbp),
    y = length(chroms) + 0.62,
    label = c("p telomere", "q telomere"),
    hjust = c(0, 1),
    size = 2.7,
    colour = "#555555"
  ) +
  scale_fill_gradientn(
    "Within-end\n% of calls",
    colours = c("#FFFFFF", "#DBEAF6", "#9ECAE1", "#4292C6", "#08519C"),
    na.value = "#F3F3F3"
  ) +
  scale_x_continuous(
    "Called PHR interval length (kbp), 10 kbp columns",
    limits = c(0, axis_end_kbp),
    breaks = x_breaks,
    labels = x_labels,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_reverse(
    "Chromosome",
    breaks = seq_along(chroms),
    labels = chroms,
    limits = c(length(chroms) + 0.78, -0.72),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(fill = guide_colourbar(barwidth = unit(3.0, "in"), barheight = unit(0.13, "in"))) +
  labs(
    title = "PHR length distribution by chromosome end",
    subtitle = paste0(
      "Each row is one chromosome end; q arms are flipped so the heatstrip reads p telomere to q telomere.\n",
      "Color is within-end share per 10 kbp bin. Terminal windows were 500 kbp; >500 kbp was not measured."
    ),
    caption = paste0(
      "Source: ", basename(data_path),
      " | ", fmt_int(sum(arm_counts > 0)), "/48 ends have called intervals; zero-signal ends: ", paste(zero_signal_arms, collapse = ", "),
      " | bin width = ", binwidth_kbp, " kbp"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_heatstrip_slide()

assets <- data.frame(
  asset = c(
    "phr_length_arm_heatstrip_10kbp.png",
    "phr_length_arm_heatstrip_10kbp.pdf"
  ),
  bin_width_kbp = c(binwidth_kbp, binwidth_kbp),
  role = c(
    "recommended_per_end_heatstrip",
    "recommended_per_end_heatstrip"
  ),
  note = c(
    "Slide-ready PNG with flipped q-arm orientation and >500 kbp not measured wording.",
    "PDF companion with flipped q-arm orientation and >500 kbp not measured wording."
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
  file.path(out_dir, "phr_length_arm_heatstrip_10kbp.png"),
  p_heatstrip,
  width = 12.8,
  height = 7.2,
  dpi = 240,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_arm_heatstrip_10kbp.pdf"),
  p_heatstrip,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

message("Wrote slide 06a v9 q-axis/kbp assets to: ", out_dir)
