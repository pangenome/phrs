#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

data_path <- Sys.getenv(
  "PHR_LENGTH_TSV",
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
)
search_cap_kb <- 500

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
  ifelse(
    is.na(x),
    "NA",
    format(round(x, digits), nsmall = digits, trim = TRUE)
  )
}

q_value <- function(x, p) {
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, type = 7))
}

wrap_text <- function(x, width) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

read_lengths <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
  df$arm_label <- sub(".*_chr([0-9XYM]+)_([pq])arm.*", "\\1\\2", df$seq)
  bad <- df$arm_label == df$seq
  if (any(bad)) {
    stop("Could not parse arm labels from seq for ", sum(bad), " rows")
  }

  valid <- df[df$region_start != ".", ]
  valid$region_start_int <- as.integer(valid$region_start)
  valid$region_end_int <- as.integer(valid$region_end)
  valid$length_kb <- (valid$region_end_int - valid$region_start_int) / 1000

  if (any(valid$length_kb > search_cap_kb, na.rm = TRUE)) {
    stop("Found reported PHR lengths above the expected ", search_cap_kb, " kb search cap")
  }

  list(raw = df, valid = valid)
}

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0(chroms, "p"), paste0(chroms, "q")))

c7_arms <- c("13p", "14p", "15p", "21p", "22p")
c14_arms <- c("Xq", "Yq")
c15_arms <- c("Xp", "Yp")
c1_arms <- c("4q", "10q")
c2_arms <- c("10p", "18p")
named_arms <- c(c7_arms, c14_arms, c15_arms, c1_arms, c2_arms)

group_levels <- c(
  "Other signaled arms",
  "C2 10p-18p",
  "C1 DUX4/D4Z4",
  "C14 PAR2",
  "C15 PAR1",
  "C7 acrocentric p"
)

group_colors <- c(
  "Other signaled arms" = "#8E8E8E",
  "C2 10p-18p" = "#F28E2B",
  "C1 DUX4/D4Z4" = "#8C6D31",
  "C14 PAR2" = "#4C78A8",
  "C15 PAR1" = "#C44E52",
  "C7 acrocentric p" = "#2A9D8F"
)

group_arms <- list(
  "C7 acrocentric p" = c7_arms,
  "C15 PAR1" = c15_arms,
  "C14 PAR2" = c14_arms,
  "C1 DUX4/D4Z4" = c1_arms,
  "C2 10p-18p" = c2_arms
)

class_for_arm <- function(arm, n = NA_integer_) {
  if (arm %in% c7_arms) return("C7 acrocentric p")
  if (arm %in% c15_arms) return("C15 PAR1")
  if (arm %in% c14_arms) return("C14 PAR2")
  if (arm %in% c1_arms) return("C1 DUX4/D4Z4")
  if (arm %in% c2_arms) return("C2 10p-18p")
  if (!is.na(n) && n == 0L) return("No interchrom PHR")
  "Other signaled arms"
}

dat <- read_lengths(data_path)
raw <- dat$raw
valid <- dat$valid

arm_n <- table(factor(valid$arm_label, levels = all_arms))
zero_arms <- names(arm_n)[as.integer(arm_n) == 0]
other_arms <- setdiff(names(arm_n)[as.integer(arm_n) > 0], named_arms)
group_arms[["Other signaled arms"]] <- other_arms

valid$group <- factor(vapply(valid$arm_label, class_for_arm, character(1)), levels = group_levels)
plot_dat <- valid[!is.na(valid$group), c("arm_label", "group", "length_kb")]

summarize_group <- function(group_name) {
  arms <- group_arms[[group_name]]
  x <- plot_dat$length_kb[as.character(plot_dat$group) == group_name]
  cap_n <- sum(x == search_cap_kb)
  data.frame(
    group = group_name,
    arms = paste(arms, collapse = ", "),
    arm_count = length(arms),
    n = length(x),
    median_kb = if (length(x) == 0) NA_real_ else median(x),
    mean_kb = if (length(x) == 0) NA_real_ else mean(x),
    p25_kb = q_value(x, 0.25),
    p75_kb = q_value(x, 0.75),
    p90_kb = q_value(x, 0.90),
    p95_kb = q_value(x, 0.95),
    max_reported_kb = if (length(x) == 0) NA_real_ else max(x),
    reported_at_cap_n = cap_n,
    reported_at_cap_pct = if (length(x) == 0) NA_real_ else 100 * cap_n / length(x),
    stringsAsFactors = FALSE
  )
}

group_summary <- do.call(rbind, lapply(group_levels, summarize_group))
group_summary$group <- factor(group_summary$group, levels = group_levels)
group_summary$cap_label <- ifelse(
  group_summary$reported_at_cap_n == 0,
  "0 at cap",
  paste0(fmt_int(group_summary$reported_at_cap_n), " at cap")
)
group_summary$x_label <- paste0(
  group_levels,
  "\n",
  "n=", fmt_int(group_summary$n),
  "; ", group_summary$cap_label
)

summary_out <- group_summary
summary_out$x_label <- NULL
round_cols <- c(
  "median_kb",
  "mean_kb",
  "p25_kb",
  "p75_kb",
  "p90_kb",
  "p95_kb",
  "max_reported_kb",
  "reported_at_cap_pct"
)
summary_out[round_cols] <- lapply(summary_out[round_cols], function(x) round(x, 1))
write.table(
  summary_out,
  file = file.path(out_dir, "named_clade_violin_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

global_median <- median(plot_dat$length_kb)
global_mean <- mean(plot_dat$length_kb)
source_basename <- basename(data_path)
zero_note <- paste0(
  "No-signal arms in this TSV: ",
  paste(zero_arms, collapse = ", "),
  " (n=0; excluded from violins)"
)

axis_labels <- setNames(group_summary$x_label, group_levels)

cap_band <- data.frame(
  xmin = -Inf,
  xmax = Inf,
  ymin = search_cap_kb,
  ymax = search_cap_kb + 30
)

cap_text <- data.frame(
  x = 3.85,
  y = search_cap_kb + 20,
  label = ">500 kb unobserved: analysis/search stops at 500 kb"
)

median_text <- paste0(
  "Global median = ",
  round(global_median),
  " kb across ",
  fmt_int(nrow(plot_dat)),
  " non-empty intervals"
)

subtitle_text <- paste0(
  "Violin = observed intervals grouped by named clade/background; box = median/IQR; ",
  "triangles mark groups with rows reported at the cap. ",
  "The search stops at 500 kb, so above-cap lengths are unobserved."
)

theme_slide <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 6, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size, colour = "#3F3F3F", margin = margin(b = 12)),
      plot.caption = element_text(size = base_size - 3, colour = "#555555", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = base_size + 1, margin = margin(r = 10)),
      axis.text.x = element_text(size = base_size - 2, lineheight = 0.94, colour = "#202020"),
      axis.text.y = element_text(size = base_size - 1, colour = "#202020"),
      legend.position = "none",
      plot.margin = margin(10, 18, 10, 16)
    )
}

p_violin <- ggplot(plot_dat, aes(x = group, y = length_kb, fill = group)) +
  geom_rect(
    data = cap_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#F1E4D3",
    alpha = 0.68
  ) +
  geom_violin(
    width = 0.84,
    trim = TRUE,
    scale = "width",
    adjust = 0.75,
    colour = "#222222",
    linewidth = 0.32,
    alpha = 0.78
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "#111111",
    linewidth = 0.42,
    alpha = 0.9
  ) +
  geom_hline(
    yintercept = global_median,
    linetype = "dotted",
    colour = "#555555",
    linewidth = 0.42
  ) +
  geom_hline(
    yintercept = search_cap_kb,
    linetype = "longdash",
    colour = "#7A3E00",
    linewidth = 0.78
  ) +
  geom_point(
    data = group_summary[group_summary$reported_at_cap_n > 0, ],
    aes(x = group, y = search_cap_kb),
    inherit.aes = FALSE,
    shape = 24,
    fill = "white",
    colour = "#7A3E00",
    stroke = 0.65,
    size = 3.2
  ) +
  geom_label(
    data = cap_text,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    label.size = 0.18,
    fill = "white",
    colour = "#4C2D0A",
    size = 3.65,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 1.13,
    y = global_median + 15,
    label = median_text,
    hjust = 0,
    size = 3.05,
    colour = "#424242"
  ) +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = axis_labels, drop = FALSE) +
  scale_y_continuous(
    "Reported PHR interval length (kb; search-capped)",
    breaks = c(0, 100, 200, 300, 400, 500),
    limits = c(0, search_cap_kb + 32),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  labs(
    title = "PHR length distributions are right-censored at the 500 kb search cap",
    subtitle = wrap_text(subtitle_text, 118),
    caption = paste0(
      "Source: ", source_basename,
      "; unit = one non-empty inter-chromosomal PHR interval from region_end - region_start. ",
      zero_note
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_slide(11)

ggsave(
  file.path(out_dir, "candidate_06a_named_clade_violin_censor.png"),
  p_violin,
  width = 12.8,
  height = 7.2,
  dpi = 220,
  bg = "white"
)

ggsave(
  file.path(out_dir, "candidate_06a_named_clade_violin_censor.pdf"),
  p_violin,
  width = 12.8,
  height = 7.2,
  bg = "white"
)

message("Wrote violin/censor candidate assets")
message("  Source basename: ", source_basename)
message("  Non-empty intervals plotted: ", fmt_int(nrow(plot_dat)))
message("  Search cap: ", search_cap_kb, " kb; above-cap lengths are not measured")
