# slide_12_meiotic_stage_trajectory.R — Mantel ρ across the 4 meiotic stages
# Output: slides/v2/slide_12_stage_trajectory.pdf  (~4 in × 2.5 in, intended as inset
# next to panel 4d so the audience sees zygotene as a peak, not a single number)
suppressPackageStartupMessages({ library(ggplot2) })

# Source: end-to-end-report/report/08_mouse.md, "Mouse meiotic Hi-C validation (1Mb window)"
# section, Mantel rho column at 50 kb resolution, 4 stages, B6+CAST per-haplotype.
df <- data.frame(
  stage = factor(c("leptotene","zygotene","pachytene","diplotene"),
                 levels = c("leptotene","zygotene","pachytene","diplotene")),
  rho   = c(0.687, 0.718, 0.683, 0.577),
  is_bouquet = c(FALSE, TRUE, FALSE, FALSE)
)

ggplot(df, aes(x = stage, y = rho, group = 1)) +
  geom_line(linewidth = 0.7, colour = "#444444") +
  geom_point(aes(colour = is_bouquet, size = is_bouquet)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -1.0, size = 3.5) +
  scale_colour_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "#d62728"),
                      guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.4, `TRUE` = 4.0), guide = "none") +
  scale_y_continuous(limits = c(0.50, 0.78), breaks = seq(0.50, 0.75, 0.05)) +
  annotate("text", x = "zygotene", y = 0.76,
           label = "bouquet\n(telomeres clustered\nat nuclear envelope)",
           size = 3.0, hjust = 0.5, colour = "#d62728", fontface = "bold") +
  labs(x = "meiotic prophase stage (Zuo et al. 2021)",
       y = "Mantel ρ (similarity × Hi-C contact, 50 kb, 1 Mb window)",
       title = "Mouse meiotic Hi-C: zygotene peak") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

ggsave("slides/v2/slide_12_stage_trajectory.pdf", width = 4.0, height = 2.5)
ggsave("slides/v2/slide_12_stage_trajectory.png", width = 4.0, height = 2.5, dpi = 300)
