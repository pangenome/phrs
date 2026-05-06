#!/usr/bin/env Rscript
# Build a compact review-zoom candidate for the slide-14 OR4F/gene-family signal.
# The image is deliberately standalone and does not edit the deck source.

`%||%` <- function(x, y) if (length(x) == 0 || is.null(x) || is.na(x)) y else x

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1] %||% "")
if (script_arg == "") {
  script_arg <- "slides/v2-review-zoom/_revision_assets/14_gene_enrichment_or4f/make_or4f_gene_family_signal.R"
}
script_path <- normalizePath(script_arg, mustWork = FALSE)
out_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(out_dir, "../../../.."), mustWork = TRUE)

hprc_plots <- "/moosefs/guarracino/HPRCv2/PHR_III/plots"
hprc_enrich <- "/moosefs/guarracino/HPRCv2/PHR_III/enrichment"

or4f_path <- file.path(hprc_plots, "or4f_pseudogene_fraction.csv")
comm_gene_path <- file.path(hprc_enrich, "community_gene_enrichment.tsv")
comm_family_path <- file.path(hprc_enrich, "community_gene_families.tsv")
fisher_path <- file.path(hprc_enrich, "community_enrichment_fisher.tsv")
copy_summary_path <- file.path(repo_root, "paper_prep/_brainstorming/gene_copy_summary.csv")

required <- c(or4f_path, comm_gene_path, comm_family_path, fisher_path, copy_summary_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required source file(s):\n", paste(missing, collapse = "\n"))
}

short_arm <- function(a) {
  sub("^chr", "", sub("_parm", "p", sub("_qarm", "q", a)))
}

format_big <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

or4f <- read.csv(or4f_path, stringsAsFactors = FALSE)
or4f <- or4f[order(or4f$pseudo_frac), ]
or4f$short <- short_arm(or4f$chr_arm)
or4f$pct <- or4f$pseudo_frac * 100

or4f_total <- sum(or4f$total)
or4f_pseudo <- sum(or4f$n_pseudo)
or4f_coding <- sum(or4f$n_coding)
or4f_mean <- 100 * or4f_pseudo / or4f_total
or4f_low <- or4f[which.min(or4f$pseudo_frac), ]
or4f_high <- or4f[which.max(or4f$pseudo_frac), ]

comm_gene <- read.delim(comm_gene_path, stringsAsFactors = FALSE)
or4f_gene_rows <- comm_gene[grepl("^OR4F", comm_gene$gene_name), ]
or4f_unique_genes <- length(unique(or4f_gene_rows$gene_name))
or4f_unique_communities <- length(unique(or4f_gene_rows$community))
or4f_gene_max <- aggregate(total_arms ~ gene_name, or4f_gene_rows, max)
or4f_gene_max <- or4f_gene_max[order(-or4f_gene_max$total_arms, or4f_gene_max$gene_name), ]
or4f_top_spread <- paste(
  paste0(head(or4f_gene_max$gene_name, 2), " on ", head(or4f_gene_max$total_arms, 2), " arms"),
  collapse = "; "
)

hub_metric <- function(gene) {
  rows <- comm_gene[comm_gene$gene_name == gene, ]
  if (nrow(rows) == 0) return("not found")
  paste0(length(unique(rows$community)), " communities / ", max(rows$total_arms), " arms")
}

fisher <- read.delim(fisher_path, stringsAsFactors = FALSE)
or_c3 <- fisher[fisher$community == "C3" & fisher$gene_family == "OR (olfactory receptor)", ][1, ]
mtco_c7 <- fisher[fisher$community == "C7" & fisher$gene_family == "MTCO (mitochondrial pseudogene)", ][1, ]

copy_summary <- read.csv(copy_summary_path, stringsAsFactors = FALSE)
copy_total <- function(pattern) {
  sum(copy_summary$total_copies[grepl(pattern, copy_summary$gene_name)])
}
copy_named <- function(names) {
  sum(copy_summary$total_copies[copy_summary$gene_name %in% names])
}

or4f_coding_copies <- copy_total("^OR4F")
dux4_frg_copies <- copy_named(c("DUX4", "FRG2", "FRG2B"))
il9r_copies <- copy_total("^IL9R")

signal_table <- data.frame(
  signal = c(
    "OR4F decay gradient",
    "OR family community signal",
    "High-copy coding families",
    "D4Z4 / DUX4L community",
    "Subtelomeric duplicon backbone"
  ),
  slide_metric = c(
    paste0(format_big(or4f_total), " OR4F annotations across ", nrow(or4f),
           " arms; ", sprintf("%.1f", or4f_mean), "% pseudogene overall"),
    paste0(or4f_unique_genes, " OR4F genes in ", or4f_unique_communities,
           " communities; ", or4f_top_spread),
    paste0("Recovered ORA table: OR4F=", or4f_coding_copies,
           ", DUX4/FRG2/FRG2B=", dux4_frg_copies,
           ", IL9R/IL9RP=", il9r_copies, " copies"),
    "C1 chr4q/chr10q carries D4Z4; section 9 reports 22 DUX4L pseudogenes specific to C1",
    paste0("RPL23AP45 ", hub_metric("RPL23AP45"), "; SEPTIN14P22 ",
           hub_metric("SEPTIN14P22"), "; DDX11L16 ", hub_metric("DDX11L16"))
  ),
  caveat = c(
    paste0("Use as canonical slide metric; source is HPRCv2 ", basename(or4f_path)),
    paste0("Presence pattern, not BH-significant enrichment: C3 OR q=",
           sprintf("%.3f", or_c3$p_adjusted)),
    "Historical Erik copy-aware ORA; parked under paper_prep/_brainstorming, use only as context",
    "Use as linked biological comparator, not as an OR4F substitute",
    "Canonical community-enrichment qualitative signal; Fisher tests do not survive BH"
  ),
  stringsAsFactors = FALSE
)

write.table(signal_table,
            file.path(out_dir, "or4f_gene_family_signal_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

md <- c(
  "# Candidate slide table: OR4F and high-copy subtelomeric gene families",
  "",
  "| Signal | Slide-ready metric | Caveat |",
  "|---|---|---|",
  apply(signal_table, 1, function(row) {
    paste0("| ", row[["signal"]], " | ", row[["slide_metric"]], " | ", row[["caveat"]], " |")
  }),
  "",
  "Recommended one-line slide claim:",
  "",
  "> OR4F is the clean visual: 5,023 annotations across 16 subtelomeric arms form a 11.1% to 99.8% pseudogenization gradient, while the broader gene-family analysis shows the same PHRs are dominated by copy-rich duplicon backbones rather than statistically significant one-family enrichments after BH correction.",
  "",
  "Source note: the OR4F gradient comes from `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`; the old copy-aware ORA numbers are recovered from `paper_prep/_brainstorming/gene_copy_summary.csv` and should be treated as contextual because that directory is explicitly parked as noncanonical."
)
writeLines(md, file.path(out_dir, "or4f_gene_family_signal_table.md"))

wrap_to <- function(text, width = 45) paste(strwrap(text, width = width), collapse = "\n")

draw_table_panel <- function() {
  plot.new()
  par(xpd = NA)
  text(0.02, 0.97, "Gene-family signal to pair with OR4F", adj = c(0, 1),
       cex = 1.25, font = 2, col = "#222222")
  text(0.02, 0.91, "Use the OR4F gradient as the visible plot; keep enrichment wording precise.",
       adj = c(0, 1), cex = 0.78, col = "#555555")

  y <- 0.82
  row_gap <- 0.155
  colors <- c("#D1495B", "#4C78A8", "#2A9D8F", "#8E6C8A", "#6B6B6B")
  for (i in seq_len(nrow(signal_table))) {
    rect(0.02, y - 0.105, 0.98, y + 0.025, col = ifelse(i %% 2 == 1, "#F7F7F7", "#FFFFFF"),
         border = "#E1E1E1")
    text(0.05, y, signal_table$signal[i], adj = c(0, 1), cex = 0.74,
         font = 2, col = colors[i])
    text(0.05, y - 0.035, wrap_to(signal_table$slide_metric[i], 70), adj = c(0, 1),
         cex = 0.60, col = "#222222")
    text(0.05, y - 0.078, wrap_to(signal_table$caveat[i], 78), adj = c(0, 1),
         cex = 0.52, col = "#666666")
    y <- y - row_gap
  }

  text(0.02, 0.06,
       paste0("Fisher note: 116 community-family tests; none BH-significant. ",
              "Best OR row C3 q=", sprintf("%.3f", or_c3$p_adjusted),
              "; C7 MTCO q=", sprintf("%.3f", mtco_c7$p_adjusted), "."),
       adj = c(0, 0), cex = 0.58, col = "#444444")
}

draw_gradient_panel <- function() {
  par(mar = c(5.2, 4.8, 3.4, 1.0))
  cols <- rep("#B9B9B9", nrow(or4f))
  cols[or4f$chr_arm == or4f_low$chr_arm] <- "#2A9D8F"
  cols[or4f$chr_arm == or4f_high$chr_arm] <- "#D1495B"

  bp <- barplot(or4f$pct, names.arg = or4f$short, las = 2, col = cols,
                border = NA, ylim = c(0, 108), cex.names = 0.82,
                ylab = "OR4F pseudogene fraction (%)",
                xlab = "", cex.axis = 0.85, cex.lab = 0.92)
  abline(h = or4f_mean, col = "#D1495B", lty = 2, lwd = 2)
  mtext("OR4F pseudogenization gradient", side = 3, adj = 0, line = 1.2,
        font = 2, cex = 1.25)
  mtext(paste0(format_big(or4f_total), " annotations across ", nrow(or4f),
               " arms; ", format_big(or4f_pseudo), " pseudogene / ",
               format_big(or4f_coding), " coding"),
        side = 3, adj = 0, line = 0.05, cex = 0.72, col = "#555555")
  text(par("usr")[2] * 0.78, or4f_mean + 2.5,
       paste0("mean ", sprintf("%.1f", or4f_mean), "%"),
       adj = c(1, 0), cex = 0.72, col = "#D1495B")

  low_idx <- which(or4f$chr_arm == or4f_low$chr_arm)
  high_idx <- which(or4f$chr_arm == or4f_high$chr_arm)
  text(bp[low_idx], or4f$pct[low_idx] + 5,
       paste0("chr", short_arm(or4f_low$chr_arm), "\n", sprintf("%.1f%%", or4f_low$pct)),
       cex = 0.68, font = 2, col = "#2A9D8F")
  text(bp[high_idx], 101,
       paste0("chr", short_arm(or4f_high$chr_arm), "\n", sprintf("%.1f%%", or4f_high$pct)),
       cex = 0.68, font = 2, col = "#D1495B")

  mtext("Canonical visual source: HPRCv2 OR4F pseudogene-fraction table; not the cropped s14_or4f.png cache.",
        side = 1, line = 4.2, adj = 0, cex = 0.56, col = "#555555")
}

render <- function(path, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(path, width = 3200, height = 1800, res = 200, type = "cairo")
  } else {
    pdf(path, width = 16, height = 9, useDingbats = FALSE)
  }
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.15, 1.0))
  par(family = "sans", bg = "white")
  draw_gradient_panel()
  draw_table_panel()
  dev.off()
}

render(file.path(out_dir, "or4f_gene_family_signal.png"), "png")
render(file.path(out_dir, "or4f_gene_family_signal.pdf"), "pdf")

cat("Wrote OR4F/gene-family candidate assets to", out_dir, "\n")
