#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

data_path <- Sys.getenv(
  "PHR_LENGTH_TSV",
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
)

search_cap_kb <- 500
binwidth_kb <- 25

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

dat <- read_lengths(data_path)
raw <- dat$raw
valid <- dat$valid
lengths <- valid$length_kb

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0(chroms, "p"), paste0(chroms, "q")))
arm_counts <- table(factor(valid$arm_label, levels = all_arms))
zero_signal_arms <- names(arm_counts)[as.integer(arm_counts) == 0]

bin_breaks <- seq(0, search_cap_kb, by = binwidth_kb)
bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")
bin_labels[length(bin_labels)] <- paste0("[", search_cap_kb - binwidth_kb, ",", search_cap_kb, "]")
valid$length_bin <- cut(
  valid$length_kb,
  breaks = bin_breaks,
  include.lowest = TRUE,
  right = FALSE,
  labels = bin_labels
)
valid$length_bin[valid$length_kb == search_cap_kb] <- bin_labels[length(bin_labels)]

bin_summary <- as.data.frame(table(valid$length_bin), stringsAsFactors = FALSE)
names(bin_summary) <- c("bin", "count")
bin_summary$bin_start_kb <- head(bin_breaks, -1)
bin_summary$bin_end_kb <- tail(bin_breaks, -1)
bin_summary$bin_mid_kb <- (bin_summary$bin_start_kb + bin_summary$bin_end_kb) / 2
bin_summary$is_cap_bin <- bin_summary$bin_end_kb == search_cap_kb
bin_summary$fill <- ifelse(bin_summary$is_cap_bin, "cap bin", "measured bins")

cap_hits <- sum(lengths == search_cap_kb)
stats <- data.frame(
  metric = c(
    "source_tsv",
    "non_empty_intervals",
    "arms_with_called_intervals",
    "total_arms_checked",
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
    fmt_int(length(lengths)),
    fmt_int(sum(arm_counts > 0)),
    fmt_int(length(all_arms)),
    paste(zero_signal_arms, collapse = ", "),
    fmt_num(min(lengths), 0),
    fmt_num(q_value(lengths, 0.25), 0),
    fmt_num(median(lengths), 0),
    fmt_num(mean(lengths), 1),
    fmt_num(q_value(lengths, 0.75), 0),
    fmt_num(q_value(lengths, 0.90), 0),
    fmt_num(q_value(lengths, 0.95), 0),
    fmt_num(max(lengths), 0),
    fmt_int(cap_hits),
    fmt_num(100 * cap_hits / length(lengths), 1),
    "analysis window ends at 500 kb; longer shared sequence is not measured"
  ),
  stringsAsFactors = FALSE
)

write.table(
  stats,
  file = file.path(out_dir, "length_distribution_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  bin_summary[, c("bin", "bin_start_kb", "bin_end_kb", "count", "is_cap_bin")],
  file = file.path(out_dir, "histogram_bins_25kb.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

median_kb <- median(lengths)
p90_kb <- q_value(lengths, 0.90)
p95_kb <- q_value(lengths, 0.95)
max_count <- max(bin_summary$count)

callout_text <- paste(
  "analysis window ends at 500 kb;",
  "longer shared sequence is not measured"
)

stats_text <- paste0(
  "n = ", fmt_int(length(lengths)), " called intervals",
  "\nmedian = ", fmt_num(median_kb, 0), " kb",
  "\n90th percentile = ", fmt_num(p90_kb, 0), " kb",
  "\n", fmt_int(cap_hits), " intervals reported at 500 kb"
)

theme_slide <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 8,
        colour = "#1a3a6b",
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        colour = "#333333",
        margin = margin(b = 16)
      ),
      plot.caption = element_text(
        size = base_size - 5,
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

p <- ggplot(bin_summary, aes(x = bin_mid_kb, y = count, fill = fill)) +
  annotate(
    "rect",
    xmin = search_cap_kb - 6,
    xmax = search_cap_kb,
    ymin = 0,
    ymax = Inf,
    fill = "#FFF2CC",
    alpha = 0.34
  ) +
  geom_col(width = binwidth_kb * 0.88, colour = "white", linewidth = 0.35) +
  geom_vline(
    xintercept = median_kb,
    linetype = "dashed",
    colour = "#222222",
    linewidth = 0.65
  ) +
  geom_vline(
    xintercept = search_cap_kb,
    colour = "#D4820A",
    linewidth = 1.25
  ) +
  annotate(
    "label",
    x = search_cap_kb - 7,
    y = max_count * 0.92,
    label = callout_text,
    hjust = 1,
    vjust = 1,
    size = 5.0,
    lineheight = 0.92,
    label.size = 0.28,
    label.r = unit(2.5, "pt"),
    fill = "#FFF8E8",
    colour = "#573700"
  ) +
  annotate(
    "segment",
    x = search_cap_kb - 55,
    xend = search_cap_kb,
    y = max_count * 0.78,
    yend = max_count * 0.78,
    colour = "#D4820A",
    linewidth = 0.7,
    arrow = arrow(length = unit(0.12, "in"), type = "closed")
  ) +
  annotate(
    "label",
    x = 284,
    y = max_count * 0.66,
    label = stats_text,
    hjust = 0,
    vjust = 1,
    size = 4.9,
    lineheight = 0.95,
    label.size = 0.22,
    label.r = unit(2.5, "pt"),
    fill = "white",
    colour = "#222222"
  ) +
  annotate(
    "text",
    x = median_kb + 6,
    y = max_count * 0.49,
    label = paste0("median ", fmt_num(median_kb, 0), " kb"),
    hjust = 0,
    size = 4.35,
    colour = "#222222"
  ) +
  scale_fill_manual(values = c("measured bins" = "#4C78A8", "cap bin" = "#D4820A")) +
  scale_x_continuous(
    "Called PHR interval length (kb)",
    limits = c(0, search_cap_kb),
    breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500),
    expand = expansion(mult = c(0.01, 0.015))
  ) +
  scale_y_continuous(
    "Intervals per 25 kb bin",
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Called PHR lengths in the 500 kb window",
    subtitle = paste0(
      "All non-empty inter-chromosomal PHR calls; 25 kb bins; no named-clade grouping.\n",
      "The right edge is a measurement ceiling, not evidence of absence beyond 500 kb."
    ),
    caption = paste0(
      "Source: ", basename(data_path),
      " | ", fmt_int(sum(arm_counts > 0)), "/48 arms have called intervals; zero-signal arms: ",
      paste(zero_signal_arms, collapse = ", "),
      " | bin width = ", binwidth_kb, " kb"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_slide()

ggsave(
  file.path(out_dir, "phr_length_histogram_restore.png"),
  p,
  width = 12.8,
  height = 7.2,
  dpi = 240,
  bg = "white"
)

ggsave(
  file.path(out_dir, "phr_length_histogram_restore.pdf"),
  p,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

message("Wrote restored slide 06a histogram assets to: ", out_dir)
