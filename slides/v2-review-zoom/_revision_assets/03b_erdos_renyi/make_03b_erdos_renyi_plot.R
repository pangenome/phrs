# Static slide candidate adapted from github.com:ekg/erdos_renyi.
# Source idea: erdos_renyi_viz.R uses the ER thresholds 1/n and log(n)/n.
# This deck-specific version plots the actual HPRC point:
# n = 18,827 flanks, p ~ 0.12 wfmash pair-evaluation rate.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required")
}

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_obs <- 18827
p_obs <- 0.12
p_star <- log(n_obs) / n_obs
ratio <- p_obs / p_star

n_seq <- exp(seq(log(10), log(100000), length.out = 500))
thresholds <- rbind(
  data.frame(
    n = n_seq,
    p = 1 / n_seq,
    threshold = "Giant component threshold: 1/n"
  ),
  data.frame(
    n = n_seq,
    p = log(n_seq) / n_seq,
    threshold = "Connectivity threshold: log(n)/n"
  )
)

point_df <- data.frame(
  n = c(n_obs, n_obs),
  p = c(p_star, p_obs),
  label = c("ER threshold", "wfmash sampling")
)

format_scientific <- function(x) {
  formatC(x, format = "e", digits = 2)
}

plot <- ggplot(thresholds, aes(x = n, y = p, color = threshold, linetype = threshold)) +
  geom_line(linewidth = 1.15) +
  geom_segment(
    aes(x = n_obs, xend = n_obs, y = p_star, yend = p_obs),
    inherit.aes = FALSE,
    linewidth = 1.1,
    color = "#1F77B4",
    arrow = grid::arrow(ends = "both", length = grid::unit(0.11, "in"))
  ) +
  geom_hline(yintercept = p_obs, color = "#1F77B4", linewidth = 0.55, alpha = 0.45) +
  geom_point(
    data = point_df,
    aes(x = n, y = p),
    inherit.aes = FALSE,
    size = c(3.0, 5.4),
    color = c("#222222", "#1F77B4")
  ) +
  annotate(
    "label",
    x = 24000,
    y = 0.19,
    label = "wfmash sampling\np ~ 0.12",
    hjust = 0,
    label.size = 0.25,
    label.r = grid::unit(0.08, "lines"),
    fill = "white",
    color = "#1F77B4",
    size = 5.0
  ) +
  annotate(
    "label",
    x = 24000,
    y = 0.0065,
    label = sprintf("%.0fx above\nthreshold", ratio),
    hjust = 0,
    label.size = 0.25,
    label.r = grid::unit(0.08, "lines"),
    fill = "#F5FAFF",
    color = "#1F77B4",
    size = 5.0
  ) +
  annotate(
    "label",
    x = 5000,
    y = 0.00062,
    label = paste0("ER threshold at n = 18,827\np* = ", format_scientific(p_star)),
    hjust = 1,
    label.size = 0.25,
    label.r = grid::unit(0.08, "lines"),
    fill = "white",
    color = "#222222",
    size = 4.4
  ) +
  annotate(
    "text",
    x = 60,
    y = 0.078,
    label = "random graph becomes connected\nabove log(n) / n",
    hjust = 0,
    color = "#222222",
    size = 4.4,
    lineheight = 0.95
  ) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000, 100000),
    labels = c("10", "100", "1,000", "10,000", "100,000"),
    limits = c(10, 100000),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1"),
    limits = c(1e-5, 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_color_manual(
    values = c(
      "Connectivity threshold: log(n)/n" = "#222222",
      "Giant component threshold: 1/n" = "#8A8A8A"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Connectivity threshold: log(n)/n" = "solid",
      "Giant component threshold: 1/n" = "22"
    )
  ) +
  labs(
    title = "Erdos-Renyi connectivity check for the implicit pangenome graph",
    subtitle = "At n = 18,827 flanks, p* = log(n)/n = 5.23e-4; wfmash evaluates ~12% of pair space",
    x = "number of nodes / subtelomeric flanks (n)",
    y = "edge probability (p, log scale)",
    caption = "Threshold curves adapted from ekg/erdos_renyi@d9ec48f; HPRC point uses n = 18,827 and p ~ 0.12."
  ) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.45),
    plot.title = element_text(face = "bold", size = 24, color = "#111111"),
    plot.subtitle = element_text(size = 15, color = "#333333", margin = margin(t = 4, b = 10)),
    axis.title = element_text(size = 15, color = "#222222"),
    axis.text = element_text(size = 13, color = "#333333"),
    legend.position = c(0.69, 0.14),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "#DDDDDD"),
    legend.key.width = grid::unit(1.4, "lines"),
    plot.caption = element_text(size = 10, color = "#666666", hjust = 0)
  )

png_path <- file.path(out_dir, "erdos_renyi_connectivity_candidate.png")
pdf_path <- file.path(out_dir, "erdos_renyi_connectivity_candidate.pdf")

ggsave(png_path, plot, width = 13.333, height = 7.5, dpi = 180, bg = "white")
ggsave(pdf_path, plot, width = 13.333, height = 7.5, bg = "white")

cat("Wrote ", png_path, "\n", sep = "")
cat("Wrote ", pdf_path, "\n", sep = "")
cat(sprintf("n=%d p_star=%.8g p_obs=%.3g ratio=%.1f\n", n_obs, p_star, p_obs, ratio))
