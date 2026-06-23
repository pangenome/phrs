#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args) >= 1) args[[1]] else "/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity"
summary_path <- file.path(package_dir, "candidate_window_summary.tsv")
out_base <- file.path(package_dir, "figures", "fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity")

summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
summary <- summary[summary$comparison_id %in% c("PAN027pat_vs_PAN011_joint", "PAN028mat_vs_PAN027_joint"), ]
summary$event_label <- factor(summary$event_label, levels = c("PAR1", "PAN027 chr9q->chr3q", "PAN028 chr9q->chr3q"))

scaffold_levels <- c("raw", "0", "10k", "20k", "50k")
len_levels <- c("raw", "default", "1000", "5000", "10000")
score_levels <- c("raw", "ani", "log-length-ani")
map_levels <- c("many:many", "1:1", "4:many")
summary$scaffold_jump <- factor(summary$scaffold_jump, levels = scaffold_levels)
summary$min_aln_length <- factor(summary$min_aln_length, levels = len_levels)
summary$scoring <- factor(summary$scoring, levels = score_levels)
summary$num_mappings <- factor(summary$num_mappings, levels = map_levels)
summary <- summary[order(summary$event_label, summary$scaffold_jump, summary$min_aln_length, summary$scoring, summary$num_mappings), ]
summary$cell_label <- ifelse(
    summary$source_kind == "raw_baseline",
    "raw\nmany",
    paste0("j", summary$scaffold_jump, "\nl", summary$min_aln_length, "\n", ifelse(summary$scoring == "log-length-ani", "logANI", "ANI"), "\n", summary$num_mappings)
)

events <- levels(summary$event_label)
cells <- unique(summary$cell_label)
z <- matrix(0, nrow = length(events), ncol = length(cells), dimnames = list(events, cells))
status <- matrix("no_rows", nrow = length(events), ncol = length(cells), dimnames = list(events, cells))
rows <- matrix(0, nrow = length(events), ncol = length(cells), dimnames = list(events, cells))
expected <- matrix(0, nrow = length(events), ncol = length(cells), dimnames = list(events, cells))
for (i in seq_len(nrow(summary))) {
    e <- as.character(summary$event_label[i])
    c <- summary$cell_label[i]
    z[e, c] <- summary$chr3_union_bp[i]
    status[e, c] <- summary$status[i]
    rows[e, c] <- summary$row_count[i]
    expected[e, c] <- summary$expected_target_union_bp[i]
}

draw_panel <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(matrix(c(1, 2), nrow = 2), heights = c(8, 1.2))
    par(mar = c(9, 12, 4, 2), xpd = NA)
    max_bp <- max(z, na.rm = TRUE)
    n_breaks <- 8
    breaks <- seq(0, max_bp, length.out = n_breaks + 1)
    palette <- c("#f2f2f2", "#fee8c8", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#7f0000", "#4d0000")
    color_index <- matrix(findInterval(as.vector(z), breaks, all.inside = TRUE), nrow = nrow(z), dimnames = dimnames(z))
    color_matrix <- matrix(palette[as.vector(color_index)], nrow = nrow(z), dimnames = dimnames(z))
    color_matrix[z == 0] <- "#f7f7f7"
    image(
        x = seq_len(ncol(z)), y = seq_len(nrow(z)), z = t(color_index[nrow(z):1, ]),
        col = palette, axes = FALSE, xlab = "", ylab = "",
        main = "Fig5 raw FASTA f16 final SweepGA filter sensitivity: chr3 union bp"
    )
    rect(
        rep(seq_len(ncol(z)) - 0.5, each = nrow(z)),
        rep(seq_len(nrow(z)) - 0.5, times = ncol(z)),
        rep(seq_len(ncol(z)) + 0.5, each = nrow(z)),
        rep(seq_len(nrow(z)) + 0.5, times = ncol(z)),
        border = "#ffffff", col = as.vector(color_matrix[nrow(z):1, ])
    )
    axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.42, tick = FALSE)
    axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 1, cex.axis = 0.8, tick = FALSE)
    box()
    for (x in seq_len(ncol(z))) {
        for (y in seq_len(nrow(z))) {
            event <- rownames(z)[nrow(z) - y + 1]
            val <- z[event, x]
            lab <- if (val >= 1000) paste0(round(val / 1000), "k") else as.character(val)
            extra <- if (rows[event, x] > 0 && val == 0) "*" else ""
            text(x, y, paste0(lab, extra), cex = 0.42, col = ifelse(val > max_bp * 0.55, "white", "#111111"))
        }
    }
    mtext("* = candidate window has rows, but none on chr3", side = 1, line = 7.8, adj = 0, cex = 0.7)

    par(mar = c(2, 12, 1, 2))
    plot.new()
    legend_labels <- c("0", paste0(round(breaks[-1] / 1000), "k"))
    legend_x <- seq_along(palette)
    rect(legend_x - 0.45, 0.45, legend_x + 0.45, 0.75, col = palette, border = NA)
    text(legend_x, 0.3, legend_labels[-1], cex = 0.7)
    text(0.4, 0.9, "chr3 union bp", adj = 0, cex = 0.8)
    text(0.4, 0.05, "Rows encode PAR1, PAN027 chr9q->chr3q, PAN028 chr9q->chr3q. Columns encode scaffold jump, min alignment length, scoring, and mapping mode.", adj = 0, cex = 0.7)
}

dir.create(dirname(out_base), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(out_base, ".pdf"), width = 18, height = 6.5)
draw_panel()
dev.off()

svg(paste0(out_base, ".svg"), width = 18, height = 6.5)
draw_panel()
dev.off()

png(paste0(out_base, ".png"), width = 3600, height = 1300, res = 200)
draw_panel()
dev.off()
