#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

data_path <- "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

out_dir <- script_dir()

read_lengths <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
  df$arm_label <- sub(".*_chr([0-9XYM]+)_([pq])arm.*", "\\1\\2", df$seq)
  bad <- df$arm_label == df$seq
  if (any(bad)) {
    stop("Could not parse arm labels from seq for ", sum(bad), " rows")
  }
  valid <- df[df$region_start != ".", ]
  valid$length_kb <- (as.integer(valid$region_end) - as.integer(valid$region_start)) / 1000
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

class_levels <- c(
  "C7 acrocentric p",
  "C15 PAR1",
  "C14 PAR2",
  "C1 DUX4/D4Z4",
  "C2 10p-18p",
  "Other signaled arms",
  "No interchrom PHR"
)

class_colors <- c(
  "C7 acrocentric p" = "#2A9D8F",
  "C15 PAR1" = "#C44E52",
  "C14 PAR2" = "#4C78A8",
  "C1 DUX4/D4Z4" = "#8C6D31",
  "C2 10p-18p" = "#F28E2B",
  "Other signaled arms" = "#9A9A9A",
  "No interchrom PHR" = "#D77A9A"
)

story_fills <- c(
  "C7 acrocentric p" = "#D7F0EA",
  "C15 PAR1" = "#FBE0E2",
  "C14 PAR2" = "#DDEBFF",
  "C1 DUX4/D4Z4" = "#F1E5CE",
  "C2 10p-18p" = "#FFE7C2",
  "Other signaled arms" = "#EFEFEF",
  "No interchrom PHR" = "#F7DCE5"
)

class_for_arm <- function(arm, n = 1L) {
  if (arm %in% c7_arms) return("C7 acrocentric p")
  if (arm %in% c15_arms) return("C15 PAR1")
  if (arm %in% c14_arms) return("C14 PAR2")
  if (arm %in% c1_arms) return("C1 DUX4/D4Z4")
  if (arm %in% c2_arms) return("C2 10p-18p")
  if (!is.na(n) && n == 0L) return("No interchrom PHR")
  "Other signaled arms"
}

q_value <- function(x, p) {
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, type = 7))
}

fmt_int <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_num <- function(x) {
  ifelse(is.na(x), "NA", format(round(x, 1), nsmall = 0, trim = TRUE))
}

wrap_text <- function(x, width) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

dat <- read_lengths(data_path)
raw <- dat$raw
valid <- dat$valid

arm_summary <- do.call(rbind, lapply(all_arms, function(arm) {
  x <- valid$length_kb[valid$arm_label == arm]
  data.frame(
    arm_label = arm,
    n = length(x),
    median_kb = if (length(x) == 0) NA_real_ else median(x),
    mean_kb = if (length(x) == 0) NA_real_ else mean(x),
    p25_kb = q_value(x, 0.25),
    p75_kb = q_value(x, 0.75),
    p90_kb = q_value(x, 0.90),
    p95_kb = q_value(x, 0.95),
    max_kb = if (length(x) == 0) NA_real_ else max(x),
    stringsAsFactors = FALSE
  )
}))

arm_summary$class <- mapply(class_for_arm, arm_summary$arm_label, arm_summary$n)
arm_summary$class <- factor(arm_summary$class, levels = class_levels)
arm_summary$chrom <- sub("[pq]$", "", arm_summary$arm_label)
arm_summary$arm <- sub("^.*([pq])$", "\\1", arm_summary$arm_label)

zero_arms <- arm_summary$arm_label[arm_summary$n == 0]
global_median <- median(valid$length_kb)
global_mean <- mean(valid$length_kb)

group_defs <- list(
  list(
    key = "C7 acrocentric p",
    label = "C7",
    arms = c7_arms,
    story = "rDNA-adjacent acrocentric p-arms are essentially full-window PHRs; this is the clearest length outlier.",
    bridge = "Say: the length plot already knows the acrocentric community."
  ),
  list(
    key = "C14 PAR2",
    label = "C14",
    arms = c14_arms,
    story = "PAR2 sits at the canonical pseudoautosomal q-end scale and gives the audience a familiar length reference.",
    bridge = "Say: ordinary subtelomeric PHRs are in the same size regime as PAR2."
  ),
  list(
    key = "C15 PAR1",
    label = "C15",
    arms = c15_arms,
    story = "PAR1 makes the sex-chromosome p-ends look like full-window shared sequence, not a generic subtelomeric tail.",
    bridge = "Say: both PAR clades are visible before community detection."
  ),
  list(
    key = "C1 DUX4/D4Z4",
    label = "C1",
    arms = c1_arms,
    story = "The 4q/10q DUX4-D4Z4 pair has a moderate center plus a real 500 kb tail, matching copy-number diversity.",
    bridge = "Say: the DUX4 clade is not a question mark; it is a named community."
  ),
  list(
    key = "C2 10p-18p",
    label = "C2",
    arms = c2_arms,
    story = "The Linardopoulou 10p/18p pair is a community callout more than a long-tail callout: short, tight, and specific.",
    bridge = "Say: not every named clade is long; some are specific exchange pairs."
  ),
  list(
    key = "No interchrom PHR",
    label = "No signal",
    arms = zero_arms,
    story = "These arms have no non-empty inter-chromosomal PHR intervals in the current TSV and therefore do not enter the 41-arm matrix.",
    bridge = "Say: the blank panels are biological absence, not a plotting failure."
  )
)

summarize_group <- function(def) {
  x <- valid$length_kb[valid$arm_label %in% def$arms]
  data.frame(
    group = def$key,
    label = def$label,
    arms = paste(def$arms, collapse = ", "),
    arm_count = length(def$arms),
    n = length(x),
    median_kb = if (length(x) == 0) NA_real_ else median(x),
    mean_kb = if (length(x) == 0) NA_real_ else mean(x),
    p90_kb = q_value(x, 0.90),
    p95_kb = q_value(x, 0.95),
    max_kb = if (length(x) == 0) NA_real_ else max(x),
    story = def$story,
    bridge = def$bridge,
    stringsAsFactors = FALSE
  )
}

clade_summary <- do.call(rbind, lapply(group_defs, summarize_group))

other_arms <- setdiff(arm_summary$arm_label[arm_summary$n > 0], named_arms)
other_x <- valid$length_kb[valid$arm_label %in% other_arms]
other_summary <- data.frame(
  group = "Other signaled arms",
  label = "Bulk",
  arms = paste(other_arms, collapse = ", "),
  arm_count = length(other_arms),
  n = length(other_x),
  median_kb = median(other_x),
  mean_kb = mean(other_x),
  p90_kb = q_value(other_x, 0.90),
  p95_kb = q_value(other_x, 0.95),
  max_kb = max(other_x),
  story = "The remaining signaled arms form the broad population-scale background.",
  bridge = "Say: this is the baseline the named clades stand out from.",
  stringsAsFactors = FALSE
)

write.table(
  arm_summary,
  file = file.path(out_dir, "arm_length_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  rbind(clade_summary, other_summary),
  file = file.path(out_dir, "clade_length_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

theme_slide <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size + 1, colour = "#444444", margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 3, colour = "#666666"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = unit(0.28, "in")
    )
}

# Candidate 06a alternative 1: ranked arm summary.
ranked <- arm_summary
ranked$p90_plot <- ifelse(is.na(ranked$p90_kb), 0, ranked$p90_kb)
ranked$rank_metric <- ifelse(ranked$n == 0, -1, ranked$p90_kb + ranked$median_kb / 1000)
ranked <- ranked[order(ranked$rank_metric, ranked$arm_label), ]
ranked$arm_rank <- factor(ranked$arm_label, levels = ranked$arm_label)

p_rank <- ggplot(ranked, aes(y = arm_rank)) +
  geom_vline(xintercept = global_median, linetype = "dashed", colour = "#555555", linewidth = 0.35) +
  geom_vline(xintercept = 334, linetype = "dotted", colour = class_colors[["C14 PAR2"]], linewidth = 0.45) +
  geom_segment(
    aes(x = 0, xend = p90_plot, yend = arm_rank, colour = class),
    linewidth = 2.15,
    alpha = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = ranked[ranked$n > 0, ],
    aes(x = median_kb),
    inherit.aes = TRUE,
    colour = "#111111",
    fill = "white",
    shape = 21,
    stroke = 0.4,
    size = 1.55
  ) +
  geom_text(
    data = ranked[ranked$n == 0, ],
    aes(x = 7, label = "n=0"),
    inherit.aes = TRUE,
    hjust = 0,
    size = 2.15,
    colour = "#6f3f51"
  ) +
  scale_colour_manual(values = class_colors, drop = FALSE) +
  scale_x_continuous(
    "PHR length per arm (kb)",
    limits = c(0, 540),
    breaks = c(0, 105, 200, 334, 500),
    labels = c("0", "105", "200", "334", "500")
  ) +
  labs(
    title = "Rank arms by the length scale of their inter-chromosomal PHRs",
    subtitle = paste0(
      "Line = 90th percentile per arm; black dot = median. ",
      "Dashed line = global median ", round(global_median), " kb; dotted line = PAR2 334 kb."
    ),
    caption = paste0(
      "Source: ", data_path, " | ",
      fmt_int(nrow(valid)), " non-empty intervals, ",
      sum(arm_summary$n > 0), "/48 arms with signal"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_slide(9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6.1),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(8, 18, 6, 8),
    legend.text = element_text(size = 7.3)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(out_dir, "candidate_06a_ranked_arm_summary.png"), p_rank, width = 12.8, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "candidate_06a_ranked_arm_summary.pdf"), p_rank, width = 12.8, height = 7.2, bg = "white")

# Candidate 06a alternative 2: focused clade and bulk facets.
valid$class <- factor(vapply(valid$arm_label, class_for_arm, character(1)), levels = class_levels)
focus_levels <- c(
  "Other signaled arms",
  "C2 10p-18p",
  "C1 DUX4/D4Z4",
  "C14 PAR2",
  "C15 PAR1",
  "C7 acrocentric p",
  "No interchrom PHR"
)
plot_valid <- valid[valid$class %in% focus_levels, c("length_kb", "class")]
plot_valid$class <- factor(as.character(plot_valid$class), levels = focus_levels)

facet_summary <- rbind(clade_summary, other_summary)
facet_summary <- facet_summary[match(focus_levels, facet_summary$group), ]
facet_summary$class <- factor(facet_summary$group, levels = focus_levels)
facet_summary$label_text <- ifelse(
  facet_summary$n == 0,
  paste0("n = 0\n", facet_summary$arms),
  paste0(
    "n = ", fmt_int(facet_summary$n),
    "\nmedian / p90 = ", fmt_num(facet_summary$median_kb), " / ", fmt_num(facet_summary$p90_kb), " kb",
    "\narms: ", facet_summary$arms
  )
)
facet_summary$label_text <- wrap_text(facet_summary$label_text, 36)

plot_valid$bin <- floor(plot_valid$length_kb / 10)
bin_counts <- aggregate(
  count ~ class + bin,
  data = data.frame(class = plot_valid$class, bin = plot_valid$bin, count = 1),
  FUN = sum
)
max_counts <- tapply(bin_counts$count, bin_counts$class, max)
facet_summary$y_pos <- as.numeric(max_counts[as.character(facet_summary$class)])
facet_summary$y_pos[is.na(facet_summary$y_pos)] <- 1
facet_summary$y_pos <- facet_summary$y_pos * 0.94
facet_summary$y_pos[facet_summary$n == 0] <- 50
facet_blank <- rbind(
  data.frame(class = facet_summary$class, x_blank = 510, y_blank = facet_summary$y_pos),
  data.frame(class = facet_summary$class, x_blank = 0, y_blank = 0)
)

p_focus <- ggplot(plot_valid, aes(x = length_kb, fill = class)) +
  geom_histogram(binwidth = 10, boundary = 0, colour = "white", linewidth = 0.18) +
  geom_blank(data = facet_blank, aes(x = x_blank, y = y_blank), inherit.aes = FALSE) +
  geom_vline(xintercept = global_median, linetype = "dashed", colour = "#555555", linewidth = 0.32) +
  geom_vline(xintercept = 334, linetype = "dotted", colour = class_colors[["C14 PAR2"]], linewidth = 0.38) +
  geom_label(
    data = facet_summary,
    aes(x = 16, y = y_pos, label = label_text),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.3,
    label.size = 0.12,
    fill = "white",
    alpha = 0.92
  ) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_x_continuous("PHR interval length (kb)", limits = c(0, 510), breaks = c(0, 105, 200, 334, 500)) +
  facet_wrap(~class, ncol = 3, scales = "free_y", drop = FALSE) +
  labs(
    title = "Focus the dense histograms on the named clades and the bulk background",
    subtitle = "This keeps the distribution-shape story while removing 48 tiny facets from the talk slide.",
    y = "Count",
    caption = paste0("Dashed = global median ", round(global_median), " kb; dotted = PAR2 334 kb. Source: ", basename(data_path))
  ) +
  theme_slide(10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9.5),
    panel.grid.major.y = element_line(colour = "#eeeeee"),
    plot.margin = margin(8, 10, 6, 8)
  )

ggsave(file.path(out_dir, "candidate_06a_focused_clade_facets.png"), p_focus, width = 12.8, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "candidate_06a_focused_clade_facets.pdf"), p_focus, width = 12.8, height = 7.2, bg = "white")

# Candidate 06b: stronger clade story matrix.
story <- clade_summary
story <- story[match(c("C7 acrocentric p", "C14 PAR2", "C15 PAR1", "C1 DUX4/D4Z4", "C2 10p-18p", "No interchrom PHR"), story$group), ]
story$signal <- ifelse(
  story$n == 0,
  "n = 0 in TSV\nabsent from 41-arm matrix",
  paste0(
    "n = ", fmt_int(story$n), " PHRs\n",
    "median / p90 = ", fmt_num(story$median_kb), " / ", fmt_num(story$p90_kb), " kb",
    ifelse(!is.na(story$max_kb) & story$max_kb == 500, "\nmax reaches 500 kb", "")
  )
)
story$title <- paste0(story$label, "\n", story$group)
story$fill <- story_fills[story$group]
story$story_wrapped <- wrap_text(story$story, 42)
story$bridge_wrapped <- wrap_text(story$bridge, 38)
story$arms_wrapped <- wrap_text(story$arms, 24)
story$signal_wrapped <- wrap_text(story$signal, 26)

cols <- data.frame(
  field = c("title", "arms_wrapped", "signal_wrapped", "story_wrapped", "bridge_wrapped"),
  header = c("Community", "Arms", "Length evidence", "Biological story", "One-sentence framing"),
  x0 = c(0.02, 0.22, 0.43, 0.67, 1.18),
  width = c(0.18, 0.19, 0.22, 0.49, 0.46),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) {
  data.frame(
    x = cols$x0[j] + cols$width[j] / 2,
    width = cols$width[j],
    y = rev(seq_len(nrow(story))),
    text = story[[cols$field[j]]],
    fill = story$fill,
    fontface = ifelse(j == 1, "bold", "plain"),
    stringsAsFactors = FALSE
  )
}))

hdr <- data.frame(
  x = cols$x0 + cols$width / 2,
  width = cols$width,
  y = nrow(story) + 1,
  text = cols$header,
  fill = "#E9EDF3",
  fontface = "bold",
  stringsAsFactors = FALSE
)

title_row <- data.frame(
  x = 0.84,
  width = 1.62,
  y = nrow(story) + 2.05,
  text = "Slide 06b should name the clades and say why each length pattern matters",
  fill = "#FFFFFF",
  fontface = "bold",
  stringsAsFactors = FALSE
)

note_row <- data.frame(
  x = 0.84,
  width = 1.62,
  y = 0.06,
  text = paste0(
    "Data check: current TSV zero-signal arms are ",
    paste(zero_arms, collapse = ", "),
    "; older notes list six and should be reconciled during final integration."
  ),
  fill = "#FFFFFF",
  fontface = "italic",
  stringsAsFactors = FALSE
)

p_story <- ggplot() +
  geom_tile(
    data = rbind(cells, hdr),
    aes(x = x, y = y, width = width, height = 0.94, fill = I(fill)),
    colour = "#C9CED6",
    linewidth = 0.25
  ) +
  geom_text(
    data = cells,
    aes(x = x - width / 2 + 0.012, y = y, label = text, fontface = fontface),
    hjust = 0,
    lineheight = 0.92,
    size = 2.72
  ) +
  geom_text(
    data = hdr,
    aes(x = x - width / 2 + 0.012, y = y, label = text, fontface = fontface),
    hjust = 0,
    size = 3.08
  ) +
  geom_text(
    data = title_row,
    aes(x = x - width / 2 + 0.012, y = y, label = text, fontface = fontface),
    hjust = 0,
    size = 3.85
  ) +
  geom_text(
    data = note_row,
    aes(x = x - width / 2 + 0.012, y = y, label = text, fontface = fontface),
    hjust = 0,
    size = 2.45,
    colour = "#555555"
  ) +
  coord_cartesian(xlim = c(0, 1.66), ylim = c(-0.28, nrow(story) + 2.42), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(out_dir, "candidate_06b_clade_story_matrix.png"), p_story, width = 12.8, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(out_dir, "candidate_06b_clade_story_matrix.pdf"), p_story, width = 12.8, height = 7.2, bg = "white")

message("Wrote slide 06 redesign assets to: ", out_dir)
