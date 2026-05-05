#!/usr/bin/env Rscript
# Extended Data Figure 4 -- Gene enrichment, copy-weighted GO, pseudogene gradient
# Inputs: see paper_prep/figures/ed4/sources.tsv
# Output: figure_ed4.pdf, figure_ed4.png

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

REPO_ROOT  <- "/moosefs/erikg/phrs/.wg-worktrees/agent-697"
PLOTS_ROOT <- "/moosefs/guarracino/HPRCv2/PHR_III/plots"
OUT <- dirname(sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(OUT) || OUT == "") OUT <- "paper_prep/figures/ed4"

# ---------- Load data ----------
go_bp <- read_csv(file.path(REPO_ROOT, "phr_GO_BP_enrichment.csv"),
                  show_col_types = FALSE)
go_mf_coding <- read_csv(file.path(REPO_ROOT, "phr_coding_only_GO_BP_enrichment.csv"),
                         show_col_types = FALSE)
copy_vs_dedup <- read_csv(file.path(REPO_ROOT,
                                    "improved_copy_weighted_vs_deduplicated_comparison.csv"),
                          show_col_types = FALSE)
gene_copy <- read_csv(file.path(REPO_ROOT, "gene_copy_summary.csv"),
                      show_col_types = FALSE)
or4f_arm <- read_csv(file.path(PLOTS_ROOT, "or4f_pseudogene_fraction.csv"),
                     show_col_types = FALSE)

short_arm <- function(a) sub("_parm", "p", sub("_qarm", "q", a))

# ---------- Panel A: GO:BP top terms (vertical bar) ----------
panel_a <- function() {
  d <- go_bp %>%
    distinct(term_name, .keep_all = TRUE) %>%
    arrange(adjusted_p_value) %>%
    head(10) %>%
    mutate(term_short = ifelse(nchar(term_name) > 52,
                               paste0(substr(term_name, 1, 50), "..."),
                               term_name),
           neglogp = -log10(adjusted_p_value)) %>%
    arrange(neglogp)

  par(mar = c(3.6, 19.5, 2.4, 2.5))
  cols <- ifelse(grepl("snRNP|spliceos|trans splic|mRNA splic|RNA splic", d$term_name),
                 "#1f77b4", "#9aaecf")
  bp <- barplot(d$neglogp, names.arg = d$term_short, horiz = TRUE, las = 1,
                col = cols, border = NA,
                xlab = expression("-log"[10]*" adjusted p"),
                xlim = c(0, max(d$neglogp) * 1.18),
                cex.names = 0.6, cex.axis = 0.7, cex.lab = 0.8)
  mtext("a   GO:BP top terms (PHR-only gene set, n=23)", side = 3,
        line = 0.6, adj = 0, cex = 0.95, font = 2)
  abline(v = -log10(0.05), col = "grey60", lty = 3)
  text(-log10(0.05), 0.4, " p=0.05", adj = c(0, 0), cex = 0.6,
       col = "grey45", srt = 0)
  text(d$neglogp + 0.05, bp,
       sprintf("%d/%d", d$intersection_size, d$term_size),
       adj = 0, cex = 0.55, col = "grey25")
}

# ---------- Panel B: Copy-weighted vs deduplicated comparison ----------
panel_b <- function() {
  d <- copy_vs_dedup %>%
    arrange(copy_fold_enrichment) %>%
    mutate(name_short = ifelse(nchar(go_name) > 38,
                               paste0(substr(go_name, 1, 36), "..."),
                               go_name))

  par(mar = c(3.6, 13.5, 2.4, 4.0))
  cols <- c("BP" = "#d62728", "MF" = "#9467bd")[d$domain]
  bp <- barplot(d$copy_fold_enrichment, names.arg = d$name_short,
                horiz = TRUE, las = 1, col = cols, border = NA,
                xlab = "Copy-weighted fold enrichment",
                xlim = c(0, max(d$copy_fold_enrichment) * 1.30),
                cex.names = 0.65, cex.axis = 0.7, cex.lab = 0.8)
  mtext("b   Copy-weighted vs deduplicated GO enrichment",
        side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)
  text(d$copy_fold_enrichment + 6, bp,
       sprintf("%.0f×  (dedup p=%.3f)", d$copy_fold_enrichment, d$prev_pvalue),
       adj = 0, cex = 0.6, col = "grey25")
  legend("bottomright", legend = c("GO:BP", "GO:MF"),
         fill = c("#d62728", "#9467bd"), border = NA, cex = 0.65,
         bty = "n", inset = c(0.02, 0.02))
  mtext("Copy column = total gene copies in PHRs; dedup column = unique gene symbols",
        side = 1, line = 2.4, cex = 0.55, col = "grey25")
}

# ---------- Panel C: High-copy gene families ----------
panel_c <- function() {
  d <- gene_copy %>%
    filter(gene_biotype %in% c("protein_coding", "pseudogene")) %>%
    arrange(desc(total_copies)) %>%
    head(15) %>%
    arrange(total_copies) %>%
    mutate(label = sprintf("%s (%s)", gene_name,
                           sub("protein_coding", "coding", gene_biotype)))
  pal <- c("protein_coding" = "#2ca02c",
           "pseudogene" = "#bcbd22")
  cols <- pal[d$gene_biotype]

  par(mar = c(3.6, 11.0, 2.4, 1.0))
  bp <- barplot(d$total_copies, names.arg = d$label, horiz = TRUE, las = 1,
                col = cols, border = NA,
                xlab = "Copy count (across PHRs)",
                xlim = c(0, max(d$total_copies) * 1.18),
                cex.names = 0.62, cex.axis = 0.7, cex.lab = 0.8)
  mtext("c   High-copy gene families (top 15, coding + pseudogene)",
        side = 3, line = 0.6, adj = 0, cex = 0.85, font = 2)
  text(d$total_copies + 0.3, bp, d$total_copies,
       adj = 0, cex = 0.6, col = "grey25")
  legend("bottomright",
         legend = c("protein_coding", "pseudogene"),
         fill = pal, border = NA, cex = 0.65, bty = "n", inset = c(0.02, 0.02))
}

# ---------- Panel D: OR4F pseudogenisation gradient by arm ----------
panel_d <- function() {
  d <- or4f_arm %>%
    arrange(pseudo_frac) %>%
    mutate(short = short_arm(chr_arm),
           pct = pseudo_frac * 100,
           is_extreme = chr_arm %in% c("chr7_parm", "chr15_qarm"))

  par(mar = c(4.6, 4.4, 2.4, 1.0))
  cols <- ifelse(d$is_extreme, "#d62728", "#7f7f7f")
  bp <- barplot(d$pct, names.arg = d$short, horiz = FALSE, las = 2,
                col = cols, border = NA,
                ylab = "OR4F pseudogene fraction (%)",
                ylim = c(0, 110),
                cex.names = 0.65, cex.axis = 0.75, cex.lab = 0.85)
  mtext("d   OR4F pseudogenisation gradient by arm",
        side = 3, line = 0.6, adj = 0, cex = 0.95, font = 2)
  abline(h = 62.1, col = "#d62728", lty = 2, lwd = 1.4)
  text(par("usr")[2], 64.5, " mean = 62.1%", adj = c(1, 0),
       cex = 0.65, col = "#d62728")
  # Label extremes
  for (i in seq_along(d$is_extreme)) {
    if (d$is_extreme[i]) {
      text(bp[i], d$pct[i] + 4, sprintf("%.1f%%", d$pct[i]),
           cex = 0.7, col = "#d62728", font = 2)
    }
  }
  text(par("usr")[2] * 0.98, 100,
       sprintf("n = %d arms\nN = %d gene-copy entries",
               nrow(d), sum(d$total)),
       adj = c(1, 1), cex = 0.6, col = "grey25")
}

# ---------- Render ----------
render_to <- function(dev_open, dev_close) {
  dev_open()
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  panel_a()
  panel_b()
  panel_c()
  panel_d()
  dev_close()
}

pdf_path <- file.path(OUT, "figure_ed4.pdf")
png_path <- file.path(OUT, "figure_ed4.png")

render_to(function() pdf(pdf_path, width = 13, height = 9), dev.off)
render_to(function() png(png_path, width = 13 * 200, height = 9 * 200,
                         res = 200, type = "cairo"), dev.off)

cat("Wrote", pdf_path, "and", png_path, "\n")
