# slide_09_clade_legend.R
# Renders the abstract-anchored legend that replaces the v1 community list on
# the right margin of the PCA scatter. Runs in seconds, no data dependencies,
# pure ggplot2.
library(ggplot2)

leg <- data.frame(
  row   = 15:1,                       # so C1 lands at top of the rendered table
  C     = sprintf("C%d", 1:15),
  arms  = c("4q, 10q",
            "10p, 18p",
            "3q, 7p, 9q, 11p, 16q, 19p",
            "7q, 12q",
            "6p, 9p, 12p, 20q",
            "1q, 13q, 17q, 19q, 21q, 22q",
            "13p, 14p, 15p, 21p, 22p",
            "15q",
            "16p",
            "17p",
            "1p, 5q, 6q, 8p",
            "2q, 20p",
            "4p",
            "Xq, Yq",
            "18q (n=1), Xp, Yp"),
  clade = c("4q–10q DUX4 / D4Z4",
            "10p–18p (Linardopoulou 2005)",
            "f7501 duplicons (fixed + AFR-enriched)",
            "private 7q/12q pair",
            "RPL23A / WASH duplicons",
            "concerted q-arm clade (22/21/19/1/13/17q)",
            "acrocentric short arms (rDNA-adjacent)",
            "chr15_q (single arm)",
            "chr16_p (single arm)",
            "chr17_p (single arm)",
            "OR4F21 sharing (Linardopoulou block 5)",
            "2q/20p pair",
            "chr4_p (single arm)",
            "PAR2 (Xq/Yq)",
            "PAR1 (Xp/Yp)"),
  fill  = c("#FDE2C8","#FFF1AA","#FFFFFF","#FFFFFF","#FFFFFF",   # C1 DUX4, C2 10p18p
            "#D6E8FF","#E5D8FA","#FFFFFF","#FFFFFF","#FFFFFF",   # C6 q-arm, C7 acro
            "#FFFFFF","#FFFFFF","#FFFFFF","#CDEAD3","#CDEAD3"),  # C14 PAR2, C15 PAR1
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x      = c(0.05, 0.30, 0.70),
  width  = c(0.22, 0.38, 1.05),
  field  = c("C","arms","clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j]/2,
  width    = cols$width[j],
  y        = leg$row,
  text     = leg[[ cols$field[j] ]],
  fill     = leg$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width/2,
  width    = cols$width,
  y        = 16,
  text     = c("Community","Arms","Abstract clade / interpretation"),
  fill     = "#EEEEEE",
  fontface = "bold",
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
            hjust = 0, size = 3.3) +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(0.5, 16.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_09_clade_legend.pdf", p, width = 8.0, height = 5.5)
ggsave("slide_09_clade_legend.png", p, width = 8.0, height = 5.5, dpi = 200)
