# slide_06_clade_callouts.R
# Renders a callout legend that names the outlier facets on the v1 length-
# distribution facet grid. Pure ggplot2; no SBATCH; runs in seconds.
library(ggplot2)

cal <- data.frame(
  row   = 5:1,                                        # C1 lands at top of table
  C     = c("C7", "C14", "C15", "C1", "C2"),
  arms  = c("13p, 14p, 15p, 21p, 22p",
            "Xq, Yq",
            "Xp, Yp  (+18q, n=1)",
            "4q, 10q",
            "10p, 18p"),
  clade = c("acrocentric short arms — fully homogenized (rDNA-adjacent)",
            "PAR2 — pseudoautosomal q-end (~334 kb scale)",
            "PAR1 — pseudoautosomal p-end",
            "4q–10q DUX4 / D4Z4 — long-tail, copy-number diverse",
            "10p–18p — Linardopoulou 2005 pair"),
  fill  = c("#E5D8FA",   # C7 acrocentric (matches slide_09 palette)
            "#CDEAD3",   # C14 PAR2
            "#CDEAD3",   # C15 PAR1
            "#FDE2C8",   # C1 DUX4
            "#FFF1AA"),  # C2 Linardopoulou
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x     = c(0.05, 0.22, 0.55),
  width = c(0.14, 0.30, 1.20),
  field = c("C", "arms", "clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j] / 2,
  width    = cols$width[j],
  y        = cal$row,
  text     = cal[[ cols$field[j] ]],
  fill     = cal$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width / 2,
  width    = cols$width,
  y        = 6,
  text     = c("C", "Arms (outlier facets)", "Named clade — why the tail is long"),
  fill     = "#EEEEEE",
  fontface = "bold",
  stringsAsFactors = FALSE
)

title <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 7,
  text     = "The fat-right-tail facets are the abstract's clades",
  fill     = "#FFFFFF",
  fontface = "bold",
  stringsAsFactors = FALSE
)

note <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 0,
  text     = "Pink facets (2p, 3p, 5p, 8q, 11q, 14q) = introvert arms, no cross-chrom hits.",
  fill     = "#FFFFFF",
  fontface = "italic",
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_tile(data = rbind(cells, hdr),
            aes(x = x, y = y, width = width, height = 0.95, fill = I(fill)),
            colour = "grey75") +
  geom_text(data = cells, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.0) +
  geom_text(data = hdr, aes(x = x - width/2 + 0.01, y = y, label = text,
                            fontface = fontface),
            hjust = 0, size = 3.2) +
  geom_text(data = title, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.6) +
  geom_text(data = note, aes(x = x - width/2 + 0.01, y = y, label = text,
                             fontface = fontface),
            hjust = 0, size = 2.6, colour = "grey35") +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(-0.5, 7.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_06_clade_callouts.pdf", p, width = 8.0, height = 3.2)
ggsave("slide_06_clade_callouts.png", p, width = 8.0, height = 3.2, dpi = 200)
