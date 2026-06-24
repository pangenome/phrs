#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
panel_dir <- if (length(args) >= 1) args[[1]] else "."
prefix <- "fig5_untangle_query_grid_style_panels"
segments_path <- file.path(panel_dir, "untangle_panel_segments.tsv")
summary_path <- file.path(panel_dir, "untangle_panel_summary.tsv")

segments <- read.delim(segments_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

event_order <- c(
  "PAR1_XY_positive_control",
  "PAN027_chr9q_chr3q_PHR_candidate",
  "PAN028_chr9q_chr3q_PHR_candidate"
)

event_labels <- c(
  PAR1_XY_positive_control = "PAR1 X/Y control",
  PAN027_chr9q_chr3q_PHR_candidate = "PAN027 chr9q -> chr3q",
  PAN028_chr9q_chr3q_PHR_candidate = "PAN028 chr9q -> chr3q"
)

target_palette <- c(
  chr3q_h1 = "#D95F02",
  chr3q_h2 = "#E6AB02",
  chr9q_h1 = "#1B9E77",
  chr9q_h2 = "#66A61E",
  chrXp = "#7570B3",
  chrYp = "#E7298A",
  side_chr15q = "#A6761D",
  side_chr16q = "#666666",
  side_chr20q = "#969696",
  other = "#8A8A8A"
)

role_border <- c(
  `same-chromosome context` = "#4C7F70",
  `PAR positive control` = "#A51E68",
  `primary donor` = "#8F3300",
  `side fragment` = "#5D5D5D",
  `low-confidence tail` = "#444444"
)

fmt_coord <- function(x) {
  ifelse(x >= 1e6, sprintf("%.3f Mb", x / 1e6), sprintf("%d kb", round(x / 1000)))
}

fmt_bp <- function(x) {
  ifelse(x >= 1000, sprintf("%.1f kb", x / 1000), sprintf("%d bp", x))
}

fmt_bp_summary <- function(x) {
  x <- gsub("bp", " bp", x, fixed = TRUE)
  x <- gsub(";", "; ", x, fixed = TRUE)
  x
}

hap_code <- function(haplotype_label) {
  ifelse(haplotype_label %in% c("haplotype1", "maternal"), "h1", "h2")
}

color_key <- function(arm, role, haplotype_label) {
  if (role %in% c("side fragment", "low-confidence tail")) {
    key <- paste0("side_", arm)
  } else if (arm %in% c("chr3q", "chr9q")) {
    key <- paste0(arm, "_", hap_code(haplotype_label))
  } else {
    key <- arm
  }
  key
}

color_for_target <- function(arm, role, haplotype_label) {
  key <- color_key(arm, role, haplotype_label)
  out <- target_palette[key]
  out[is.na(out)] <- target_palette[["other"]]
  unname(out)
}

role_lane <- function(role) {
  ifelse(role %in% c("side fragment", "low-confidence tail"), -0.28, 0)
}

draw_legend <- function() {
  par(xpd = NA)
  labels <- c("3q h1/mat", "3q h2/pat", "9q h1/mat", "9q h2/pat",
              "Xp", "Yp", "side/caveat")
  cols <- c(target_palette["chr3q_h1"], target_palette["chr3q_h2"],
            target_palette["chr9q_h1"], target_palette["chr9q_h2"],
            target_palette["chrXp"], target_palette["chrYp"], target_palette["side_chr15q"])
  x <- -0.02
  y <- 4.52
  gap <- 0.130
  for (i in seq_along(labels)) {
    x0 <- x + (i - 1) * gap
    rect(x0, y - 0.055, x0 + 0.018, y + 0.055, col = cols[i], border = NA)
    text(x0 + 0.024, y, labels[i], adj = 0, cex = 0.58)
  }
  segments(0.93, y, 0.975, y, col = "#555555", lwd = 1.1, lty = 3)
  text(0.982, y, "side fragments are caveat markers", adj = 0, cex = 0.62)
  par(xpd = FALSE)
}

draw_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(3.1, 8.8, 3.9, 4.1), xaxs = "i", yaxs = "i", family = "sans")
  plot(
    NA,
    xlim = c(-0.50, 1.52),
    ylim = c(0.25, 4.86),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  title("Fig5 strict untangle primary-path panels", cex.main = 1.04, line = 2.25)
  mtext(
    "Geometry: selected_segments.tsv strict nb=1 + SweepGA 1:1 primary path; native genomic query coordinates; not raw FASTA f16/f32 query-grid output",
    side = 3,
    line = 0.55,
    cex = 0.70
  )
  draw_legend()

  y_positions <- c(3.64, 2.36, 1.08)
  names(y_positions) <- event_order

  for (event_id in event_order) {
    y <- y_positions[event_id]
    srow <- summary[summary$event_id == event_id, ]
    rows <- segments[segments$event_id == event_id, ]
    rows <- rows[order(rows$native_query_start_0based, rows$native_query_end_0based_exclusive,
                       rows$event_role_order, rows$target_arm), ]
    w0 <- as.numeric(srow$query_window_start_0based[1])
    w1 <- as.numeric(srow$query_window_end_0based_exclusive[1])
    xscale <- function(x) (x - w0) / (w1 - w0)

    segments(-0.49, y + 0.56, 1.26, y + 0.56, col = "#D0D0D0", lwd = 0.8)
    text(-0.49, y + 0.36, event_labels[event_id], adj = 0, cex = 0.86, font = 2)
    text(-0.49, y + 0.16, srow$transmission[1], adj = 0, cex = 0.62, col = "#555555")
    text(-0.49, y - 0.04, srow$query_native_window_0based_half_open[1], adj = 0, cex = 0.55, col = "#555555")

    rect(0, y - 0.16, 1, y + 0.16, col = "#F3F3F3", border = "#C7C7C7", lwd = 0.7)

    for (j in seq_len(nrow(rows))) {
      q0 <- as.numeric(rows$native_query_start_0based[j])
      q1 <- as.numeric(rows$native_query_end_0based_exclusive[j])
      x0 <- max(0, min(1, xscale(q0)))
      x1 <- max(0, min(1, xscale(q1)))
      role <- rows$event_role[j]
      arm <- rows$target_arm[j]
      lane_y <- y + role_lane(role)
      h <- if (role %in% c("side fragment", "low-confidence tail")) 0.12 else 0.24
      fill <- color_for_target(arm, role, rows$target_haplotype_label[j])
      border <- role_border[role]
      if (is.na(border)) border <- "#333333"
      lty <- if (role %in% c("side fragment", "low-confidence tail")) 3 else 1
      rect(x0, lane_y - h / 2, x1, lane_y + h / 2, col = fill, border = border, lwd = 0.55, lty = lty)

      if (role %in% c("primary donor", "PAR positive control") && (x1 - x0) > 0.025) {
        text((x0 + x1) / 2, lane_y, arm, cex = 0.48, col = "white", font = 2)
      }
      if (role %in% c("side fragment", "low-confidence tail")) {
        segments((x0 + x1) / 2, lane_y - 0.07, (x0 + x1) / 2, y - 0.42,
                 col = "#555555", lwd = 0.45, lty = 3)
        text((x0 + x1) / 2, y - 0.48, arm, cex = 0.44, srt = 35, adj = 1, col = "#444444")
      }
    }

    ticks <- pretty(c(w0, w1), n = 4)
    ticks <- ticks[ticks >= w0 & ticks <= w1]
    for (tick in ticks) {
      x <- xscale(tick)
      segments(x, y - 0.36, x, y - 0.29, col = "#444444", lwd = 0.55)
      text(x, y - 0.43, fmt_coord(tick), cex = 0.46, adj = c(0.5, 1))
    }

    text(
      1.025,
      y + 0.10,
      sprintf(
        "%s rows; primary %s",
        srow$segment_rows[1],
        fmt_bp_summary(srow$primary_donor_arms_bp[1])
      ),
      adj = 0,
      cex = 0.50,
      col = "#333333"
    )
    text(
      1.025,
      y - 0.10,
      sprintf("side/caveat: %s", srow$side_fragment_caveat_labels[1]),
      adj = 0,
      cex = 0.47,
      col = "#555555"
    )
  }

  mtext("Native sample assembly query coordinates (0-based half-open; each row spans its selected 500 kb source window)", side = 1, line = 2.15, cex = 0.76)
}

pdf(file.path(panel_dir, paste0(prefix, ".pdf")), width = 13.5, height = 6.0, useDingbats = FALSE)
draw_panel()
dev.off()

svg(file.path(panel_dir, paste0(prefix, ".svg")), width = 13.5, height = 6.0, onefile = TRUE)
draw_panel()
dev.off()

png(file.path(panel_dir, paste0(prefix, ".png")), width = 2700, height = 1200, res = 200, type = "cairo")
draw_panel()
dev.off()
