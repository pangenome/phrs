#!/usr/bin/env Rscript

# Functional component x chromosome-end heatmaps for the V7 CHM13
# copy-number-aware ontology analysis.
#
# The primary plotted unit is coordinate-distinct PHR copy burden. These plots
# do not count redundant GO/Reactome rows as independent biological events.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

.libPaths(c("~/R/library", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else
  "paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7/community_attribution"
out_dir <- if (length(args) >= 2) args[[2]] else
  file.path(base_dir, "functional_component_heatmap")

base_dir <- normalizePath(base_dir, mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

projection_path <- file.path(base_dir, "ARM_FUNCTIONAL_CLASS_PROJECTION.tsv")
class_summary_path <- file.path(base_dir, "FUNCTIONAL_CLASS_COMMUNITY_SUMMARY.tsv")
contributors_path <- file.path(dirname(base_dir), "EXACT_TERM_CONTRIBUTORS.tsv.gz")
term_results_path <- file.path(dirname(base_dir), "TERM_RESULTS.tsv.gz")
nj_newick_path <- "paper_prep/figures/nj_tree_arms/nj_tree.newick"
upgma_order_path <- "submission/fig/MainFigures/arm_order_tree.tsv"
chrom_sizes_path <- "/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/chrom.sizes"

required <- c(projection_path, class_summary_path, contributors_path, term_results_path)
missing_required <- required[!file.exists(required)]
if (length(missing_required) > 0) {
  stop("Missing required input(s): ", paste(missing_required, collapse = ", "))
}

read_tsv <- function(path) {
  fread(path, sep = "\t", quote = "", na.strings = c("", "NA"))
}

chrom_key <- function(chrom) {
  x <- sub("^chr", "", chrom)
  out <- suppressWarnings(as.integer(x))
  out[x == "X"] <- 23L
  out[x == "Y"] <- 24L
  out
}

all_chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0("chr", all_chroms, "_p"),
                            paste0("chr", all_chroms, "_q")))
all_arm_labels <- sub("^chr", "", gsub("_", "", all_arms))

format_arm_label <- function(x) sub("^chr", "", gsub("_", "", x))

community_number <- function(x) {
  out <- suppressWarnings(as.integer(sub("^C", "", x)))
  out[is.na(out)] <- 9999L
  out
}

projection <- read_tsv(projection_path)
class_summary <- read_tsv(class_summary_path)

needed_projection_cols <- c("row_type", "chrom_arm", "chrom_arm_label", "community",
                            "class_id", "display_label",
                            "class_unique_coordinate_copy_burden")
if (!all(needed_projection_cols %in% names(projection))) {
  stop("Projection table missing columns: ",
       paste(setdiff(needed_projection_cols, names(projection)), collapse = ", "))
}

class_order <- c(
  "DUX4_ZGA_TRANSCRIPTION_NUCLEAR_ENVELOPE_CELL_CYCLE",
  "WASH_ENDOSOMAL_ACTIN_EXOCYST",
  "DDX11_HELICASE_CHROMOSOME",
  "SEPTIN14_SEPTIN_CYTOKINESIS",
  "WBP1L_CXCL12_SIGNALING",
  "RPL23A_RIBOSOMAL_NUCLEOLAR"
)

class_labels <- c(
  DUX4_ZGA_TRANSCRIPTION_NUCLEAR_ENVELOPE_CELL_CYCLE = "DUX4 / ZGA",
  WASH_ENDOSOMAL_ACTIN_EXOCYST = "WASH / endosomal actin",
  DDX11_HELICASE_CHROMOSOME = "DDX11 / helicase",
  SEPTIN14_SEPTIN_CYTOKINESIS = "SEPTIN14 / cytokinesis",
  WBP1L_CXCL12_SIGNALING = "WBP1L / CXCL12",
  RPL23A_RIBOSOMAL_NUCLEOLAR = "RPL23A / nucleolar"
)

present_classes <- intersect(class_order, unique(projection$class_id))
if (length(present_classes) == 0) {
  stop("None of the expected display classes are present in ", projection_path)
}

community_by_arm <- unique(projection[, .(chrom_arm, community)])
community_by_arm <- community_by_arm[chrom_arm %in% all_arms]
community_by_arm <- community_by_arm[order(chrom_arm)]
community_by_arm <- community_by_arm[, .(community = community[1]), by = chrom_arm]
missing_arm_meta <- setdiff(all_arms, community_by_arm$chrom_arm)
if (length(missing_arm_meta) > 0) {
  community_by_arm <- rbind(
    community_by_arm,
    data.table(chrom_arm = missing_arm_meta, community = "NO_SIGNAL_ARM")
  )
}
community_by_arm[, chrom_arm_label := format_arm_label(chrom_arm)]

component_counts <- projection[
  row_type == "FUNCTIONAL_CLASS" & class_id %in% present_classes,
  .(copy_burden = sum(as.numeric(class_unique_coordinate_copy_burden), na.rm = TRUE)),
  by = .(class_id, display_label, chrom_arm)
]

full_label_map <- unique(projection[
  row_type == "FUNCTIONAL_CLASS" & class_id %in% present_classes &
    !is.na(display_label),
  .(class_id, full_display_label = display_label)
])

component_grid <- CJ(class_id = present_classes, chrom_arm = all_arms, unique = TRUE)
component_counts <- merge(component_grid, component_counts,
                          by = c("class_id", "chrom_arm"), all.x = TRUE)
component_counts[is.na(copy_burden), copy_burden := 0]
component_counts[, row_label := unname(class_labels[class_id])]
component_counts <- merge(component_counts, community_by_arm, by = "chrom_arm", all.x = TRUE)
component_counts <- merge(component_counts, full_label_map, by = "class_id", all.x = TRUE)

get_nj_order <- function() {
  if (!file.exists(nj_newick_path)) return(NULL)
  if (!requireNamespace("ape", quietly = TRUE)) return(NULL)
  tr <- ape::read.tree(nj_newick_path)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp, width = 4, height = 4)
  on.exit({
    invisible(dev.off())
    unlink(tmp)
  }, add = TRUE)
  plot(tr, show.tip.label = FALSE, plot = FALSE)
  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  tip_y <- lp$yy[seq_along(tr$tip.label)]
  tr$tip.label[order(tip_y)]
}

get_upgma_order <- function() {
  if (!file.exists(upgma_order_path)) return(NULL)
  x <- fread(upgma_order_path, sep = "\t", quote = "")
  if (!all(c("ChromArm", "position_left_to_right") %in% names(x))) return(NULL)
  x[order(position_left_to_right), ChromArm]
}

signal_arms <- community_by_arm[community != "NO_SIGNAL_ARM", chrom_arm]
no_signal_arms <- community_by_arm[community == "NO_SIGNAL_ARM", chrom_arm]

tree_signal_order <- get_nj_order()
tree_source <- "NJ Newick tip order from paper_prep/figures/nj_tree_arms/nj_tree.newick"
if (is.null(tree_signal_order)) {
  tree_signal_order <- get_upgma_order()
  tree_source <- "UPGMA tree order from submission/fig/MainFigures/arm_order_tree.tsv"
}
if (is.null(tree_signal_order)) {
  tree_signal_order <- signal_arms[order(community_number(community_by_arm[match(signal_arms, chrom_arm), community]),
                                         chrom_key(sub("_.*$", "", signal_arms)),
                                         signal_arms)]
  tree_source <- "Fallback community/genomic order; no NJ/UPGMA order available"
}
tree_signal_order <- tree_signal_order[tree_signal_order %in% signal_arms]
if (!setequal(tree_signal_order, signal_arms)) {
  missing_tree_arms <- setdiff(signal_arms, tree_signal_order)
  tree_signal_order <- c(tree_signal_order, missing_tree_arms)
}

nj_pos <- setNames(seq_along(tree_signal_order), tree_signal_order)
community_signal_order <- community_by_arm[community != "NO_SIGNAL_ARM"]
community_signal_order[, community_num := community_number(community)]
community_signal_order[, tree_pos := nj_pos[chrom_arm]]
community_signal_order[is.na(tree_pos), tree_pos := 9999L]
community_signal_order <- community_signal_order[
  order(community_num, tree_pos, chrom_key(sub("_.*$", "", chrom_arm)), chrom_arm),
  chrom_arm
]

orderings <- list(
  community_order = c(community_signal_order, no_signal_arms[order(chrom_key(sub("_.*$", "", no_signal_arms)), no_signal_arms)]),
  nj_order = c(tree_signal_order, no_signal_arms[order(chrom_key(sub("_.*$", "", no_signal_arms)), no_signal_arms)])
)

ordering_labels <- c(
  community_order = "Community order: C1-C15, tree order within each community",
  nj_order = "Tree order: NJ tip order, no-signal arms appended"
)

write_order_table <- function(order_name, order_vec) {
  dt <- data.table(
    ordering = order_name,
    position = seq_along(order_vec),
    chrom_arm = order_vec,
    chrom_arm_label = format_arm_label(order_vec)
  )
  dt <- merge(dt, community_by_arm[, .(chrom_arm, community)], by = "chrom_arm", all.x = TRUE)
  setorder(dt, position)
  dt[, ordering_source := ifelse(order_name == "nj_order", tree_source,
                                 ordering_labels[[order_name]])]
  dt
}

order_table <- rbindlist(Map(write_order_table, names(orderings), orderings), fill = TRUE)
fwrite(order_table, file.path(out_dir, "functional_component_arm_orderings.tsv"), sep = "\t")

matrix_export <- copy(component_counts)
matrix_export[, display_label := row_label]
matrix_export <- matrix_export[, .(
  class_id, display_label, full_display_label, chrom_arm, chrom_arm_label,
  community, copy_burden
)]
fwrite(matrix_export, file.path(out_dir, "functional_component_arm_matrix.tsv"), sep = "\t")

bin_copy <- function(x) {
  cut(x,
      breaks = c(-Inf, 0, 1, 2, 5, 10, Inf),
      labels = c("0", "1", "2", "3-5", "6-10", ">10"),
      right = TRUE)
}

fill_values <- c(
  "0" = "#f5f5f5",
  "1" = "#d8eef3",
  "2" = "#8fc6dd",
  "3-5" = "#f3b36d",
  "6-10" = "#d95f4f",
  ">10" = "#7f1d1d"
)

community_palette <- c(
  C1 = "#1b9e77", C2 = "#d95f02", C3 = "#7570b3", C4 = "#e7298a",
  C5 = "#66a61e", C6 = "#e6ab02", C7 = "#a6761d", C8 = "#1f78b4",
  C9 = "#b2df8a", C10 = "#fb9a99", C11 = "#6a3d9a", C12 = "#ff7f00",
  C13 = "#b15928", C14 = "#17becf", C15 = "#7f7f7f",
  NO_SIGNAL_ARM = "#bdbdbd"
)

make_blocks <- function(order_vec) {
  dt <- data.table(position = seq_along(order_vec), chrom_arm = order_vec)
  dt <- merge(dt, community_by_arm[, .(chrom_arm, community)], by = "chrom_arm", all.x = TRUE)
  setorder(dt, position)
  r <- rle(dt$community)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  data.table(
    community = r$values,
    start = starts,
    end = ends,
    mid = (starts + ends) / 2,
    n_arms = r$lengths
  )
}

plot_one <- function(order_name, order_vec, out_prefix) {
  plot_dt <- copy(component_counts)
  plot_dt <- plot_dt[chrom_arm %in% order_vec & class_id %in% present_classes]
  plot_dt[, x := match(chrom_arm, order_vec)]
  plot_dt[, y := match(class_id, rev(present_classes))]
  plot_dt[, x_label := format_arm_label(chrom_arm)]
  plot_dt[, fill_bin := bin_copy(copy_burden)]
  plot_dt[, text_label := fifelse(copy_burden > 0, as.character(copy_burden), "")]
  plot_dt[, text_color := fifelse(copy_burden > 5, "white", "#1f1f1f")]
  blocks <- make_blocks(order_vec)
  no_signal_block <- blocks[community == "NO_SIGNAL_ARM"]
  signal_blocks <- blocks[community != "NO_SIGNAL_ARM"]

  n_rows <- length(present_classes)
  x_breaks <- seq_along(order_vec)
  x_labels <- format_arm_label(order_vec)
  y_breaks <- seq_along(rev(present_classes))
  y_labels <- unname(class_labels[rev(present_classes)])

  # Invisible layer carrying one real row per copy bin, so every legend key
  # renders even for bins (3-5, 6-10) that no chromosome end actually reaches.
  legend_dummy <- data.table(
    x = 1L, y = 1L,
    fill_bin = factor(c("0", "1", "2", "3-5", "6-10", ">10"),
                      levels = c("0", "1", "2", "3-5", "6-10", ">10"))
  )

  # x-axis tick labels colored by p/q arm, matching the Fig 2b/2c palette.
  x_label_cols <- ifelse(grepl("_p$", order_vec), "#CC3B38", "#1F5EA8")

  p <- ggplot(plot_dt, aes(x = x, y = y)) +
    {if (nrow(no_signal_block) > 0) annotate(
      "rect",
      xmin = no_signal_block$start - 0.5,
      xmax = no_signal_block$end + 0.5,
      ymin = 0.5,
      ymax = n_rows + 0.5,
      fill = "#e8e8e8",
      alpha = 0.55
    )} +
    geom_tile(aes(fill = fill_bin), color = "#ffffff", linewidth = 0.35) +
    geom_tile(data = legend_dummy, aes(x = x, y = y, fill = fill_bin),
              alpha = 0, inherit.aes = FALSE,
              show.legend = c(fill = TRUE, colour = FALSE)) +
    geom_text(
      data = plot_dt[copy_burden > 0 & copy_burden <= 5],
      aes(label = text_label),
      color = "#1f1f1f", size = 3.9, fontface = "bold", show.legend = FALSE
    ) +
    geom_text(
      data = plot_dt[copy_burden > 5],
      aes(label = text_label),
      color = "white", size = 3.9, fontface = "bold", show.legend = FALSE
    ) +
    geom_point(
      data = data.table(x = 1, y = 1, arm = factor(c("p", "q"), levels = c("p", "q"))),
      aes(x = x, y = y, colour = arm), shape = 15,
      alpha = 0, inherit.aes = FALSE,
      show.legend = c(fill = FALSE, colour = TRUE)
    ) +
    scale_colour_manual(values = c(p = "#CC3B38", q = "#1F5EA8"),
                        name = NULL, labels = c("p arm", "q arm")) +
    scale_fill_manual(values = fill_values, drop = FALSE, name = "Gene copies in PHR") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.5, n_rows + 1.1), clip = "off") +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Chromosome end",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13,
                                 colour = x_label_cols),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 16, margin = margin(t = 6)),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.6, "cm"),
      legend.spacing.x = unit(1.4, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      plot.margin = margin(t = 10, r = 10, b = 8, l = 10)
    ) +
    guides(fill = guide_legend(nrow = 1, label.position = "right", order = 1,
                               override.aes = list(fill = unname(fill_values),
                                                   alpha = 1,
                                                   colour = "#cccccc")),
           colour = guide_legend(nrow = 1, label.position = "right", order = 2,
                                 override.aes = list(alpha = 1, size = 5)))

  if (nrow(blocks) > 0) {
    p <- p +
      geom_vline(
        xintercept = blocks[end < length(order_vec), end + 0.5],
        color = "#4d4d4d",
        linewidth = 0.25
      )
  }
  if (nrow(signal_blocks) > 0) {
    p <- p +
      annotate(
        "text",
        x = signal_blocks$mid,
        y = n_rows + 0.72,
        label = signal_blocks$community,
        size = ifelse(signal_blocks$community %in% c("C10", "C13"), 3.0, 3.7),
        fontface = "bold",
        color = community_palette[signal_blocks$community]
      )
  }
  if (nrow(no_signal_block) > 0) {
    p <- p +
      annotate(
        "text",
        x = no_signal_block$mid,
        y = n_rows + 0.72,
        label = "no PHR signal",
        size = 3.5,
        fontface = "bold",
        color = "#737373"
      )
  }

  pdf_path <- paste0(out_prefix, ".pdf")
  png_path <- paste0(out_prefix, ".png")
  ggsave(pdf_path, p, width = 15.0, height = 5.0, units = "in",
         device = cairo_pdf, bg = "white")
  ggsave(png_path, p, width = 15.0, height = 5.0, units = "in",
         dpi = 220, bg = "white")
  invisible(p)
}

plots <- list()
for (nm in names(orderings)) {
  prefix <- file.path(out_dir, paste0("functional_component_arm_heatmap.", nm))
  plots[[nm]] <- plot_one(nm, orderings[[nm]], prefix)
}

combined_pdf <- file.path(out_dir, "functional_component_arm_heatmap.review_pages.pdf")
pdf(combined_pdf, width = 12.8, height = 4.6, onefile = TRUE)
for (nm in names(plots)) print(plots[[nm]])
dev.off()

# Full exact supported-term matrix for auditing the "could we plot all 209?"
# question. This is not the primary figure because many exact terms are
# ancestor/direct duplicates over the same physical-copy patterns.
make_exact_term_audit <- function() {
  term_results <- fread(cmd = paste("zcat", shQuote(term_results_path)),
                        sep = "\t", quote = "", na.strings = c("", "NA"))
  primary_terms <- term_results[
    assignment == "midpoint" &
      primary_support == 1 &
      support_status == "PRIMARY_SUPPORTED",
    .(collection, relation, term_id, term_name, holm_p_global,
      bh_q_within_collection, term_phr_burden_a_T = a_T)
  ]
  if (nrow(primary_terms) == 0) return(NULL)

  contrib <- fread(cmd = paste("zcat", shQuote(contributors_path)),
                   sep = "\t", quote = "", na.strings = c("", "NA"))
  contrib <- contrib[
    v7_midpoint_partition == "PHR" &
      v7_midpoint_count_cell == "a_T" &
      v7_primary_support == 1
  ]
  if (nrow(contrib) == 0) return(NULL)

  sizes <- fread(chrom_sizes_path, header = FALSE, col.names = c("seqid", "size"))
  contrib <- merge(contrib, sizes, by = "seqid", all.x = TRUE)
  contrib[, midpoint := (as.numeric(start0) + as.numeric(end0)) / 2]
  contrib[, arm := fifelse(!is.na(size) & midpoint > size / 2, "q", "p")]
  contrib[, chrom_arm := paste0(seqid, "_", arm)]

  contrib <- unique(contrib[, .(collection, relation, term_id, term_name,
                                copy_id, chrom_arm, functional_source_symbol)])
  exact_counts <- contrib[
    chrom_arm %in% all_arms,
    .(copy_burden = uniqueN(copy_id),
      sources = paste(sort(unique(functional_source_symbol)), collapse = "|")),
    by = .(collection, relation, term_id, term_name, chrom_arm)
  ]
  exact_grid <- primary_terms[
    ,
    .(chrom_arm = all_arms),
    by = .(collection, relation, term_id, term_name, holm_p_global,
           bh_q_within_collection, term_phr_burden_a_T)
  ]
  exact_counts <- merge(exact_grid, exact_counts,
                        by = c("collection", "relation", "term_id",
                               "term_name", "chrom_arm"),
                        all.x = TRUE)
  exact_counts[is.na(copy_burden), copy_burden := 0]
  exact_counts[is.na(sources), sources := ""]
  exact_counts[, chrom_arm_label := format_arm_label(chrom_arm)]
  setorder(exact_counts, holm_p_global, collection, relation, term_id, chrom_arm)
  fwrite(exact_counts,
         file.path(out_dir, "exact_supported_term_arm_matrix.tsv.gz"),
         sep = "\t")

  wide <- dcast(exact_counts,
                collection + relation + term_id + term_name + holm_p_global +
                  bh_q_within_collection + term_phr_burden_a_T ~ chrom_arm,
                value.var = "copy_burden",
                fill = 0)
  arm_cols <- intersect(all_arms, names(wide))
  wide[, pattern_key := do.call(paste, c(.SD, sep = ",")), .SDcols = arm_cols]
  arm_count_matrix <- as.matrix(wide[, arm_cols, with = FALSE])
  wide[, nonzero_arm_counts := apply(arm_count_matrix, 1, function(vals) {
    vals <- as.integer(vals)
    nz <- which(vals > 0)
    if (length(nz) == 0) return("")
    paste(paste0(format_arm_label(arm_cols[nz]), ":", vals[nz]), collapse = "; ")
  })]
  pattern_summary <- wide[, .(
    n_exact_terms = .N,
    min_holm_p_global = min(holm_p_global, na.rm = TRUE),
    max_term_phr_burden_a_T = max(term_phr_burden_a_T, na.rm = TRUE),
    nonzero_arm_counts = nonzero_arm_counts[1],
    representative_terms = paste(head(paste0(collection, ":", relation, ":",
                                            term_name, " (", term_id, ")"), 8),
                                 collapse = " | ")
  ), by = pattern_key]
  pattern_summary[, pattern_id := sprintf("P%03d", .I)]
  setorder(pattern_summary, -n_exact_terms, min_holm_p_global)
  fwrite(pattern_summary,
         file.path(out_dir, "exact_supported_term_copy_pattern_summary.tsv"),
         sep = "\t")

  invisible(list(exact_counts = exact_counts, pattern_summary = pattern_summary))
}

exact_audit <- make_exact_term_audit()

manifest <- data.table(
  output = c(
    "functional_component_arm_heatmap.community_order.pdf",
    "functional_component_arm_heatmap.community_order.png",
    "functional_component_arm_heatmap.nj_order.pdf",
    "functional_component_arm_heatmap.nj_order.png",
    "functional_component_arm_heatmap.review_pages.pdf",
    "functional_component_arm_matrix.tsv",
    "functional_component_arm_orderings.tsv",
    "exact_supported_term_arm_matrix.tsv.gz",
    "exact_supported_term_copy_pattern_summary.tsv"
  ),
  description = c(
    "Six display-class heatmap, community-blocked chromosome-end order",
    "PNG companion to community-blocked heatmap",
    "Six display-class heatmap, NJ tree chromosome-end order",
    "PNG companion to NJ-order heatmap",
    "Two-page review PDF containing both heatmap orderings",
    "Long matrix behind the six display-class heatmaps",
    "Chromosome-end orderings and source notes",
    "Full 209 exact supported-term by chromosome-end copy-count matrix",
    "Audit showing how many exact terms share identical arm-copy patterns"
  )
)
fwrite(manifest, file.path(out_dir, "OUTPUTS.tsv"), sep = "\t")

cat("Wrote functional component heatmap outputs to ", out_dir, "\n", sep = "")
cat("Tree ordering source: ", tree_source, "\n", sep = "")
if (!is.null(exact_audit)) {
  cat("Exact-term audit patterns: ", nrow(exact_audit$pattern_summary), "\n", sep = "")
}
