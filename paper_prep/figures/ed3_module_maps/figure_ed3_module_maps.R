#!/usr/bin/env Rscript
# ED3 (module maps) - see README.md: legacy paper_prep/figures/ed3/ holds an unrelated
# superseded TAR1/telomere figure, so this figure lives in its own directory.
# Module maps of subtelomeric arm ends: 25 kb bins labelled by the partner arm that
# contributes the most aligned bp to them (HPRCv2 all-vs-all flank alignments).
# Rows are exemplar arms, one per Leiden community; fill is the partner's community,
# and tiles wide enough carry the partner arm label.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .args[grep("^--file=", .args)])
here <- if (length(.f)) dirname(normalizePath(.f)) else getwd()
d <- fread(file.path(here, "module_tiling.tsv"))
setnames(d, c("arm","bin","start","end","partner_arm","partner_community","bp"))
# exemplars: per community, the arm with the most labelled bins, then the widest partner
d[, arm_comm := sub("_.*", "", arm)]
ex <- d[, .(nbins = .N), by = arm][order(-nbins)][1:min(12, .N)]
dd <- d[arm %in% ex$arm]
dd[, arm := factor(arm, levels = rev(ex$arm))]
dd[, partner_arm := factor(partner_arm)]
# label only partner runs of >= 3 consecutive bins
setorder(dd, arm, bin)
dd[, run := cumsum(c(1, diff(as.integer(partner_arm)) != 0)), by = arm]
runs <- dd[, .(n = .N, x = mean(bin), lab = as.character(partner_arm[1])), by = .(arm, run)]
runs <- runs[n >= 3]
p <- ggplot(dd, aes(x = bin, y = arm, fill = partner_community)) +
  geom_tile(width = 1.02, height = 0.92) +
  geom_text(data = runs, aes(x = x, y = arm, label = lab), inherit.aes = FALSE,
            size = 2.6, colour = "#202124") +
  scale_fill_manual(values = c(C1="#1f77b4",C2="#ff7f0e",C3="#2ca02c",C4="#d62728",
    C5="#9467bd",C6="#8c564b",C7="#e377c2",C8="#7f7f7f",C9="#bcbd22",C10="#17becf",
    C11="#aec7e8",C12="#ffbb78",C13="#98df8a",C14="#ff9896",C15="#c5b0d5")) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20),
                     labels = function(x) x * 25, expand = c(0, 0)) +
  labs(x = "position in the telomeric 500 kb flank (kb, telomere at 0)",
       y = NULL, fill = "partner\ncommunity",
       title = "Subtelomeric module maps",
       subtitle = "each bin labelled by the partner arm contributing the most aligned bp; fill = that partner's community") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(size = 7.5, colour = "#5f6368"),
        legend.key.size = unit(0.7, "lines"))
ggsave(file.path(here, "figure_ed3_module_maps.pdf"), p, width = 9.5, height = 4.6)
ggsave(file.path(here, "figure_ed3_module_maps.png"), p, width = 9.5, height = 4.6, dpi = 150)
cat("wrote", file.path(here, "figure_ed3_module_maps.pdf"), "\n")
