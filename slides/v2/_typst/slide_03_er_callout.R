# slide_03_er_callout.R — ER connectivity threshold sanity check
# Output: slide_03_er_callout.pdf  (~3 in × 2 in, intended as inset on slide 3)
library(ggplot2)
n <- 18827
p_star  <- log(n) / n            # ER connectivity threshold ≈ 5.21e-4
p_obs   <- 0.12                  # wfmash k-mer-prefilter evaluation rate
ratio   <- p_obs / p_star        # ≈ 230×

df <- data.frame(
  label = c("ER threshold\np* = log(n)/n", "wfmash sampling\np ≈ 12%"),
  p     = c(p_star, p_obs),
  fill  = c("threshold", "observed")
)
df$label <- factor(df$label, levels = df$label)

ggplot(df, aes(label, p, fill = fill)) +
  geom_col(width = 0.55) +
  scale_y_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("1e-4","1e-3","1e-2","1e-1","1")
  ) +
  scale_fill_manual(values = c(threshold = "#888888", observed = "#1f77b4")) +
  geom_text(aes(label = sprintf("%.2g", p)), vjust = -0.5, size = 3.5) +
  annotate("text", x = 1.5, y = 0.4,
           label = sprintf("~%.0f× above threshold\n→ densely connected\n→ closure reaches\n   genome-wide", ratio),
           size = 3.2, hjust = 0.5) +
  labs(x = NULL, y = "edge probability  p   (log scale)",
       title = sprintf("n = %s flanks", format(n, big.mark = ","))) +
  guides(fill = "none") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave("slide_03_er_callout.pdf", width = 3.0, height = 2.0)
ggsave("slide_03_er_callout.png", width = 3.0, height = 2.0, dpi = 300)
