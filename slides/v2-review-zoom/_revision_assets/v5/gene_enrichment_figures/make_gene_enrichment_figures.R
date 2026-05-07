#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 R package is required to render these figures.", call. = FALSE)
  }
  library(ggplot2)
  library(grid)
})

`%||%` <- function(x, y) if (length(x) == 0 || is.null(x) || is.na(x)) y else x

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
}

find_repo_root <- function(start_dir) {
  here <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "chm13.phrs.bed")) &&
        file.exists(file.path(here, "slides/v2-review-zoom"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find repo root from ", start_dir, call. = FALSE)
    }
    here <- parent
  }
}

out_dir <- script_dir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
repo_root <- find_repo_root(out_dir)

hprc_base <- Sys.getenv(
  "HPRCV2_PHR_III_DIR",
  "/moosefs/guarracino/HPRCv2/PHR_III"
)
hprc_enrich <- file.path(hprc_base, "enrichment")
hprc_plots <- file.path(hprc_base, "plots")

paths <- list(
  inventory_source = file.path(
    repo_root,
    "slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/source_inventory.tsv"
  ),
  inventory_candidates = file.path(
    repo_root,
    "slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/candidate_enrichment_signals.tsv"
  ),
  chm13_phrs = file.path(repo_root, "chm13.phrs.bed"),
  or4f = file.path(hprc_plots, "or4f_pseudogene_fraction.csv"),
  community_summary = file.path(hprc_enrich, "community_summary_table.tsv"),
  community_gene = file.path(hprc_enrich, "community_gene_enrichment.tsv"),
  community_family = file.path(hprc_enrich, "community_gene_families.tsv"),
  community_fisher = file.path(hprc_enrich, "community_enrichment_fisher.tsv"),
  community_specific = file.path(hprc_enrich, "community_specific_genes.tsv"),
  cross_community = file.path(hprc_enrich, "cross_community_genes.tsv")
)

missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  stop("Missing required source file(s):\n", paste(missing, collapse = "\n"),
       call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(path,
             sep = "\t",
             header = TRUE,
             quote = "",
             comment.char = "",
             check.names = FALSE,
             stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x,
              file = path,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              na = "")
}

fmt_int <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_num <- function(x, digits = 1) {
  formatC(x, format = "f", digits = digits)
}

fmt_q <- function(x) {
  formatC(x, format = "f", digits = 3)
}

wrap_one <- function(x, width = 42) {
  paste(strwrap(x, width = width), collapse = "\n")
}

wrap_vec <- function(x, width = 42) {
  vapply(x, wrap_one, character(1), width = width)
}

normalize_arm <- function(x) {
  x <- trimws(x)
  x <- sub("_parm$", "_p", x)
  x <- sub("_qarm$", "_q", x)
  x
}

short_arm <- function(x) {
  gsub("_", "", sub("^chr", "", normalize_arm(x)), fixed = TRUE)
}

split_arms <- function(x) {
  normalize_arm(unlist(strsplit(gsub("\\s+", "", x), ",", fixed = FALSE)))
}

source_label <- function(path) basename(path)

require_one <- function(df, label) {
  if (nrow(df) == 0) stop("Expected at least one row for ", label, call. = FALSE)
  df[1, , drop = FALSE]
}

candidate_caveat <- function(candidates, rank) {
  hit <- candidates[candidates$rank == rank, "caveat"]
  hit[1] %||% ""
}

candidate_source <- function(candidates, rank) {
  hit <- candidates[candidates$rank == rank, "source_path"]
  hit[1] %||% ""
}

or4f <- read.csv(paths$or4f, stringsAsFactors = FALSE)
summary_tbl <- read_tsv(paths$community_summary)
gene_tbl <- read_tsv(paths$community_gene)
family_tbl <- read_tsv(paths$community_family)
fisher_tbl <- read_tsv(paths$community_fisher)
specific_tbl <- read_tsv(paths$community_specific)
cross_tbl <- read_tsv(paths$cross_community)
source_inventory <- read_tsv(paths$inventory_source)
candidates <- read_tsv(paths$inventory_candidates)

bed <- read.delim(paths$chm13_phrs,
                  sep = "\t",
                  header = FALSE,
                  quote = "",
                  comment.char = "",
                  stringsAsFactors = FALSE)
names(bed)[1:4] <- c("chrom", "start0", "end0", "chrs_involved")
bed$arm <- paste0(bed$chrom, "_", ifelse(bed$start0 < 1000000, "p", "q"))
called_arms <- sort(unique(bed$arm))

community_interval_status <- do.call(rbind, lapply(seq_len(nrow(summary_tbl)), function(i) {
  comm <- summary_tbl$Community[i]
  arms <- split_arms(summary_tbl$Arms[i])
  missing_arms <- setdiff(arms, called_arms)
  called <- intersect(arms, called_arms)
  data.frame(
    community = comm,
    community_arms = paste(arms, collapse = ", "),
    community_arms_short = paste(short_arm(arms), collapse = " "),
    n_community_arms = length(arms),
    n_called_chm13_phr_arms = length(called),
    n_missing_chm13_phr_arms = length(missing_arms),
    missing_chm13_phr_arms = paste(missing_arms, collapse = ", "),
    interval_scope_label = if (length(missing_arms) > 0) {
      paste0("broader community signal; missing CHM13 called PHR: ",
             paste(missing_arms, collapse = ", "))
    } else {
      "all community arms have called CHM13 PHR rows"
    },
    stringsAsFactors = FALSE
  )
}))

write_tsv(community_interval_status,
          file.path(out_dir, "community_interval_status.tsv"))
community_interval_status <- read_tsv(file.path(out_dir, "community_interval_status.tsv"))

or4f$arm <- normalize_arm(or4f$chr_arm)
or4f$total <- as.numeric(or4f$total)
or4f$n_pseudo <- as.numeric(or4f$n_pseudo)
or4f$n_coding <- as.numeric(or4f$n_coding)
or4f$pseudo_frac <- as.numeric(or4f$pseudo_frac)
or4f_total <- sum(or4f$total)
or4f_pseudo <- sum(or4f$n_pseudo)
or4f_coding <- sum(or4f$n_coding)
or4f_pct <- 100 * or4f_pseudo / or4f_total
or4f_low <- or4f[which.min(or4f$pseudo_frac), , drop = FALSE]
or4f_high <- or4f[which.max(or4f$pseudo_frac), , drop = FALSE]

or4f_gene_rows <- gene_tbl[grepl("^OR4F", gene_tbl$gene_name), , drop = FALSE]
or4f_gene_max <- aggregate(total_arms ~ gene_name, or4f_gene_rows, max)
or4f_gene_max <- or4f_gene_max[order(-or4f_gene_max$total_arms, or4f_gene_max$gene_name), ]
or4f_top_arms <- max(or4f_gene_max$total_arms)
or4f_unique_genes <- length(unique(or4f_gene_rows$gene_name))
or4f_unique_communities <- length(unique(or4f_gene_rows$community))

or_c3 <- require_one(
  fisher_tbl[fisher_tbl$community == "C3" &
               fisher_tbl$gene_family == "OR (olfactory receptor)", , drop = FALSE],
  "C3 OR Fisher row"
)
mtco_c7_fisher <- require_one(
  fisher_tbl[fisher_tbl$community == "C7" &
               fisher_tbl$gene_family == "MTCO (mitochondrial pseudogene)", , drop = FALSE],
  "C7 MTCO Fisher row"
)

dux4l_c1 <- specific_tbl[
  specific_tbl$community == "C1" & grepl("^DUX4L", specific_tbl$gene_name),
  , drop = FALSE
]
if (nrow(dux4l_c1) == 0) stop("No C1 DUX4L rows found.", call. = FALSE)

cross_pick <- function(gene) {
  require_one(cross_tbl[cross_tbl$gene_name == gene, , drop = FALSE],
              paste0("cross-community gene ", gene))
}
rpl <- cross_pick("RPL23AP45")
septin <- cross_pick("SEPTIN14P22")
ddx <- cross_pick("DDX11L16")
fam <- cross_pick("FAM138D")
wash <- cross_pick("WASH6P")

family_pick <- function(comm, family_name) {
  require_one(family_tbl[family_tbl$community == comm &
                           family_tbl$gene_family == family_name, , drop = FALSE],
              paste(comm, family_name))
}
c5_ddx <- family_pick("C5", "DDX11L (DEAD-box helicase pseudogene)")
c5_wash <- family_pick("C5", "WASH (actin nucleation)")
c5_fam <- family_pick("C5", "FAM138 family")
iqsec3 <- require_one(gene_tbl[gene_tbl$gene_name == "IQSEC3", , drop = FALSE],
                      "IQSEC3 gene row")
gtpbp6 <- require_one(gene_tbl[gene_tbl$gene_name == "GTPBP6", , drop = FALSE],
                      "GTPBP6 gene row")
c15_summary <- require_one(summary_tbl[summary_tbl$Community == "C15", , drop = FALSE],
                           "C15 summary row")
c7_summary <- require_one(summary_tbl[summary_tbl$Community == "C7", , drop = FALSE],
                          "C7 summary row")
c7_mtco <- specific_tbl[
  specific_tbl$community == "C7" & grepl("^MTCO", specific_tbl$gene_name),
  , drop = FALSE
]
if (nrow(c7_mtco) == 0) stop("No C7 MTCO rows found.", call. = FALSE)

sig_count <- sum(toupper(as.character(fisher_tbl$significant)) == "TRUE")
min_q <- min(as.numeric(fisher_tbl$p_adjusted), na.rm = TRUE)

method_support <- data.frame(
  display_order = 1:4,
  lane = c(
    "Standard ORA",
    "Deduplicated coding-only ORA",
    "HPRCv2 copy-aware support",
    "Genome-wide copy-weighted ORA"
  ),
  counted_unit = c(
    "unique gene symbols",
    "unique coding symbols",
    "arms, samples, annotations, gene-family rows",
    "gene copies against a genome-wide background"
  ),
  source_scope = c(
    "parked CHM13/no-acro PHR gene-list ORA",
    "parked CHM13/no-acro coding-only ORA",
    "canonical HPRCv2 arm-level community enrichment",
    "parked exploratory copy-weighted hypergeometric workstream"
  ),
  plotted_in_v5 = c("no", "contrast only", "yes", "no"),
  support_summary = c(
    "Collapses repeated subtelomeric copies to one symbol.",
    "Inventory notes query size 9; useful only as historical contrast.",
    paste0(sum(summary_tbl$N_Arms), " community arms across ",
           nrow(summary_tbl), " communities; ", fmt_int(sum(summary_tbl$N_Genes)),
           " community gene records."),
    "True copy-weighted setup, but validation documents Type I error and FDR concerns."
  ),
  caveat = c(
    "Not copy-number-aware.",
    "Small parked query; not the HPRCv2 community source.",
    "Support-aware, not genome-wide copy-weighted ORA; Fisher rows do not survive BH.",
    "Do not present 598x, 928x, or 309x rows as final statistics without reanalysis."
  ),
  source_paths = c(
    "paper_prep/_brainstorming/phr_GO_*.csv",
    "paper_prep/_brainstorming/phr_coding_only_GO_*.csv",
    paste(paths$community_summary, paths$community_gene, paths$community_family,
          paths$community_fisher, sep = "; "),
    "paper_prep/_brainstorming/improved_copy_weighted_*.csv"
  ),
  stringsAsFactors = FALSE
)
write_tsv(method_support, file.path(out_dir, "method_comparison_support.tsv"))
method_support <- read_tsv(file.path(out_dir, "method_comparison_support.tsv"))

ranked_signal_support <- data.frame(
  rank = 1:8,
  signal = c(
    "OR4F pseudogenization gradient",
    "OR/OR4F community presence",
    "D4Z4/DUX4L C1 block",
    "Duplicon backbone genes",
    "C5 DDX11L/WASH/FAM138 module",
    "GTPBP6/IQSEC3 GTP anchors",
    "C15 PAR1 coding outlier",
    "C7 acrocentric MTCO signal"
  ),
  plotted_support_arms = c(
    nrow(or4f),
    or4f_top_arms,
    max(dux4l_c1$arms_in_comm),
    max(c(rpl$n_arms, septin$n_arms, ddx$n_arms, fam$n_arms, wash$n_arms)),
    max(c(c5_ddx$n_arms, c5_wash$n_arms, c5_fam$n_arms)),
    iqsec3$total_arms + gtpbp6$total_arms,
    c15_summary$N_Arms,
    c7_summary$N_Arms
  ),
  support_axis_label = c(
    paste0(nrow(or4f), " annotation-source arms"),
    paste0(or4f_top_arms, " arms for top OR4F genes"),
    paste0(max(dux4l_c1$arms_in_comm), " C1 arms"),
    paste0(max(c(rpl$n_arms, septin$n_arms, ddx$n_arms, fam$n_arms, wash$n_arms)),
           " arms for top backbone gene"),
    paste0(max(c(c5_ddx$n_arms, c5_wash$n_arms, c5_fam$n_arms)), " C5 arms"),
    paste0(iqsec3$total_arms + gtpbp6$total_arms, " canonical gene-support arms"),
    paste0(c15_summary$N_Arms, " PAR1 community arms"),
    paste0(c7_summary$N_Arms, " acrocentric p-arms")
  ),
  key_metric = c(
    paste0(fmt_int(or4f_total), " OR4F annotations; ",
           fmt_num(or4f_pct), "% pseudogene overall; ",
           fmt_num(100 * or4f_low$pseudo_frac), "% at ",
           short_arm(or4f_low$arm), " to ",
           fmt_num(100 * or4f_high$pseudo_frac), "% at ",
           short_arm(or4f_high$arm), "."),
    paste0(or4f_unique_genes, " OR4F genes across ",
           or4f_unique_communities, " communities; C3 OR q=",
           fmt_q(or_c3$p_adjusted), "."),
    paste0(nrow(dux4l_c1), " C1-specific DUX4L pseudogenes on chr4q/chr10q."),
    paste0("RPL23AP45 ", rpl$n_communities, " comm/", rpl$n_arms,
           " arms; SEPTIN14P22 ", septin$n_communities, "/", septin$n_arms,
           "; DDX11L16 ", ddx$n_communities, "/", ddx$n_arms, "."),
    paste0("C5 has DDX11L ", c5_ddx$n_unique_genes, " genes, WASH ",
           c5_wash$n_unique_genes, ", FAM138 ", c5_fam$n_unique_genes,
           "; IQSEC3 samples=", iqsec3$samples_in_comm, "."),
    paste0("GTPBP6 in C15 samples=", gtpbp6$samples_in_comm,
           "; IQSEC3 in C5 samples=", iqsec3$samples_in_comm,
           "; parked 309x GO statistic not used as final."),
    paste0(c15_summary$N_Specific, " C15-specific genes; ",
           fmt_num(c15_summary$Pct_Protein), "% protein-coding content."),
    paste0(nrow(c7_mtco), " C7-specific MTCO pseudogenes; MTCO Fisher q=",
           fmt_q(mtco_c7_fisher$p_adjusted), ".")
  ),
  claim_scope = c(
    "PHR annotation copies",
    "HPRCv2 arm/community support",
    "Community-specific genes",
    "Cross-community hub support",
    "HPRCv2 arm/community support",
    "Mixed canonical gene support",
    "Community-specific genes",
    "Community-specific genes"
  ),
  interval_note = c(
    "PHR sequence annotation signal, not a CHM13-only interval claim.",
    "Broader arm/community signal; gate any interval render through chm13.phrs.bed.",
    "Both C1 arms have called CHM13 PHR rows.",
    "Broader subtelomeric arm/community signal; not every listed arm is a rendered CHM13 interval.",
    "C5 includes chr6_p, which lacks a called CHM13 PHR row.",
    "GTPBP6 C15 includes chrY_p missing from CHM13 called PHR rows; IQSEC3 anchor is chr12_p.",
    "C15 includes chrY_p, which lacks a called CHM13 PHR row.",
    "C7 includes chr13_p, which lacks a called CHM13 PHR row."
  ),
  caveat = c(
    candidate_caveat(candidates, 1),
    candidate_caveat(candidates, 2),
    candidate_caveat(candidates, 3),
    candidate_caveat(candidates, 4),
    candidate_caveat(candidates, 5),
    candidate_caveat(candidates, 6),
    candidate_caveat(candidates, 7),
    candidate_caveat(candidates, 8)
  ),
  source_paths = c(
    candidate_source(candidates, 1),
    candidate_source(candidates, 2),
    candidate_source(candidates, 3),
    candidate_source(candidates, 4),
    candidate_source(candidates, 5),
    candidate_source(candidates, 6),
    candidate_source(candidates, 7),
    candidate_source(candidates, 8)
  ),
  stringsAsFactors = FALSE
)
write_tsv(ranked_signal_support, file.path(out_dir, "ranked_signal_support.tsv"))
ranked_signal_support <- read_tsv(file.path(out_dir, "ranked_signal_support.tsv"))

selected_communities <- c("C1", "C3", "C5", "C7", "C11", "C12", "C14", "C15", "C4")
family_specs <- data.frame(
  signal = c("OR", "RPL", "SEPTIN", "DDX11L", "WASH", "FAM138"),
  gene_family = c(
    "OR (olfactory receptor)",
    "RPL pseudogene",
    "SEPTIN (septin pseudogene)",
    "DDX11L (DEAD-box helicase pseudogene)",
    "WASH (actin nucleation)",
    "FAM138 family"
  ),
  stringsAsFactors = FALSE
)
special_signals <- c("DUX4L specific", "MTCO specific", "PAR1 coding anchors")
all_signals <- c(family_specs$signal, special_signals)

community_family_rows <- list()
row_i <- 1L
for (comm in selected_communities) {
  status <- require_one(community_interval_status[community_interval_status$community == comm, , drop = FALSE],
                        paste("interval status", comm))
  for (signal in all_signals) {
    support_arms <- 0
    support_value <- 0
    support_label <- ""
    source_table <- ""
    scope <- "HPRCv2 community-arm support"
    if (signal %in% family_specs$signal) {
      fam_name <- family_specs$gene_family[family_specs$signal == signal]
      fam_row <- family_tbl[family_tbl$community == comm & family_tbl$gene_family == fam_name,
                            , drop = FALSE]
      source_table <- source_label(paths$community_family)
      if (nrow(fam_row) > 0) {
        support_arms <- fam_row$n_arms[1]
        support_value <- fam_row$n_unique_genes[1]
        support_label <- paste0(support_arms, " arms")
      }
    } else if (signal == "DUX4L specific" && comm == "C1") {
      support_arms <- max(dux4l_c1$arms_in_comm)
      support_value <- nrow(dux4l_c1)
      support_label <- paste0(support_value, " genes; ", support_arms, " arms")
      source_table <- source_label(paths$community_specific)
      scope <- "community-specific genes"
    } else if (signal == "MTCO specific" && comm == "C7") {
      support_arms <- max(c7_mtco$arms_in_comm)
      support_value <- nrow(c7_mtco)
      support_label <- paste0(support_value, " genes; ", support_arms, " arms")
      source_table <- source_label(paths$community_specific)
      scope <- "community-specific genes"
    } else if (signal == "PAR1 coding anchors" && comm == "C15") {
      anchors <- gene_tbl[
        gene_tbl$community == "C15" &
          gene_tbl$gene_name %in% c("SHOX", "GTPBP6", "P2RY8", "PLCXD1", "PPP2R3B"),
        , drop = FALSE
      ]
      support_arms <- max(anchors$arms_in_comm)
      support_value <- nrow(anchors)
      support_label <- paste0(support_value, " coding; ", support_arms, " arms")
      source_table <- source_label(paths$community_gene)
      scope <- "community-specific coding anchors"
    }
    community_family_rows[[row_i]] <- data.frame(
      community = comm,
      signal = signal,
      support_arms = support_arms,
      support_value = support_value,
      tile_label = support_label,
      source_table = source_table,
      signal_scope = scope,
      n_missing_chm13_phr_arms = status$n_missing_chm13_phr_arms,
      missing_chm13_phr_arms = status$missing_chm13_phr_arms,
      interval_scope_label = status$interval_scope_label,
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}
community_family_support <- do.call(rbind, community_family_rows)
write_tsv(community_family_support, file.path(out_dir, "community_family_map_support.tsv"))
community_family_support <- read_tsv(file.path(out_dir, "community_family_map_support.tsv"))

open_device <- function(path, type) {
  if (type == "png") {
    png(path, width = 3200, height = 1800, res = 200, type = "cairo")
  } else {
    pdf(path, width = 16, height = 9, useDingbats = FALSE)
  }
}

close_device <- function() {
  invisible(dev.off())
}

render_method <- function(path, type) {
  open_device(path, type)
  on.exit(close_device(), add = TRUE)
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(
    "Copy-number-aware gene enrichment: what v5 can claim",
    x = unit(0.04, "npc"), y = unit(0.93, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 25, fontface = "bold", col = "#1F2933")
  )
  grid.text(
    "PHR/subtelomeric genes recur across arms and haplotypes. The slide assets use canonical HPRCv2 arm/community support, not parked weighted-ORA p-values.",
    x = unit(0.04, "npc"), y = unit(0.875, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 13.5, col = "#52616B")
  )

  card_colors <- c("#F3F4F6", "#F3F4F6", "#E8F5F0", "#FFF3E5")
  border_colors <- c("#CBD2D9", "#CBD2D9", "#2A9D8F", "#D9822B")
  status_colors <- c("#52616B", "#52616B", "#13795B", "#A65F00")
  x0 <- c(0.045, 0.285, 0.525, 0.765)
  y_top <- 0.78
  w <- 0.19
  h <- 0.48

  for (i in seq_len(nrow(method_support))) {
    pushViewport(viewport(x = unit(x0[i], "npc"), y = unit(y_top, "npc"),
                          width = unit(w, "npc"), height = unit(h, "npc"),
                          just = c("left", "top")))
    grid.roundrect(r = unit(0.04, "snpc"),
                   gp = gpar(fill = card_colors[i], col = border_colors[i], lwd = 1.4))
    grid.text(wrap_one(method_support$lane[i], 24),
              x = unit(0.06, "npc"), y = unit(0.90, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 13.2, fontface = "bold", col = "#1F2933",
                        lineheight = 1.02))
    grid.text(toupper(method_support$plotted_in_v5[i]),
              x = unit(0.06, "npc"), y = unit(0.74, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 10.5, fontface = "bold", col = status_colors[i]))
    grid.text(wrap_one(method_support$counted_unit[i], 26),
              x = unit(0.06, "npc"), y = unit(0.63, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 10.5, col = "#293845", lineheight = 1.05))
    grid.text(wrap_one(method_support$support_summary[i], 28),
              x = unit(0.06, "npc"), y = unit(0.43, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 9.5, col = "#3E4C59", lineheight = 1.05))
    grid.text(wrap_one(method_support$caveat[i], 30),
              x = unit(0.06, "npc"), y = unit(0.20, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 8.5, col = "#6B7280", lineheight = 1.05))
    popViewport()
    if (i < nrow(method_support)) {
      grid.lines(x = unit(c(x0[i] + w + 0.015, x0[i + 1] - 0.015), "npc"),
                 y = unit(c(0.54, 0.54), "npc"),
                 arrow = arrow(length = unit(0.08, "inches")),
                 gp = gpar(col = "#9AA5B1", lwd = 1.4))
    }
  }

  grid.roundrect(x = unit(0.5, "npc"), y = unit(0.15, "npc"),
                 width = unit(0.91, "npc"), height = unit(0.13, "npc"),
                 r = unit(0.02, "snpc"),
                 gp = gpar(fill = "#F8FAFC", col = "#D9E2EC", lwd = 1))
  grid.text(
    paste0("Use wording such as 'copy-aware support' or 'copy-number-aware candidate signal'. ",
           "HPRCv2 community-family Fisher tests: ", nrow(fisher_tbl),
           " rows, ", sig_count, " BH-significant; minimum q=", fmt_q(min_q), "."),
    x = unit(0.055, "npc"), y = unit(0.175, "npc"),
    just = c("left", "center"),
    gp = gpar(fontsize = 12, col = "#293845")
  )
  grid.text(
    "Figure source: local method_comparison_support.tsv generated from canonical inventory and HPRCv2 community tables.",
    x = unit(0.04, "npc"), y = unit(0.045, "npc"),
    just = c("left", "bottom"),
    gp = gpar(fontsize = 9.5, col = "#6B7280")
  )
}

ranked_plot <- function(ranked) {
  ranked$rank <- as.integer(ranked$rank)
  ranked$plotted_support_arms <- as.numeric(ranked$plotted_support_arms)
  ranked$row_label <- paste0(
    ranked$rank, ". ", ranked$signal, "\n",
    wrap_vec(ranked$key_metric, 54)
  )
  ranked$row_label <- factor(ranked$row_label, levels = rev(ranked$row_label))
  ranked$claim_scope <- factor(
    ranked$claim_scope,
    levels = c(
      "PHR annotation copies",
      "HPRCv2 arm/community support",
      "Cross-community hub support",
      "Community-specific genes",
      "Mixed canonical gene support"
    )
  )
  fill_values <- c(
    "PHR annotation copies" = "#4C78A8",
    "HPRCv2 arm/community support" = "#2A9D8F",
    "Cross-community hub support" = "#7B8CDE",
    "Community-specific genes" = "#D1495B",
    "Mixed canonical gene support" = "#D9822B"
  )
  ggplot(ranked, aes(x = row_label, y = plotted_support_arms, fill = claim_scope)) +
    geom_col(width = 0.64, color = "white", linewidth = 0.25) +
    geom_text(aes(label = support_axis_label),
              hjust = -0.05, size = 3.3, color = "#293845") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = fill_values, drop = FALSE, name = NULL) +
    scale_y_continuous(limits = c(0, max(ranked$plotted_support_arms) + 7),
                       breaks = seq(0, 24, by = 4),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title = "Ranked copy-aware candidate signals for review-zoom v5",
      subtitle = "Rows follow the upstream inventory rank; bar length is source/community arm support, not statistical significance.",
      x = NULL,
      y = "HPRCv2 community/source arms supporting signal",
      caption = paste0(
        "Canonical HPRCv2 support tables and OR4F plot table. ",
        "Fisher caveat: ", nrow(fisher_tbl), " family-community tests, ",
        sig_count, " BH-significant; best q=", fmt_q(min_q), "."
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 19, color = "#1F2933"),
      plot.subtitle = element_text(size = 11.5, color = "#52616B", margin = margin(b = 10)),
      plot.caption = element_text(size = 9.5, color = "#6B7280", hjust = 0),
      axis.text.y = element_text(size = 9.5, color = "#1F2933", lineheight = 0.92),
      axis.text.x = element_text(size = 10, color = "#52616B"),
      axis.title.x = element_text(size = 10.5, color = "#293845", margin = margin(t = 8)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9.5),
      plot.margin = margin(14, 72, 18, 18)
    )
}

community_plot <- function(community_family, interval_status) {
  selected_status <- interval_status[match(selected_communities, interval_status$community), ]
  x_labels <- vapply(seq_len(nrow(selected_status)), function(i) {
    arms <- unlist(strsplit(selected_status$community_arms_short[i], " ", fixed = TRUE))
    arms_wrapped <- paste(strwrap(paste(arms, collapse = " "), width = 14), collapse = "\n")
    star <- if (selected_status$n_missing_chm13_phr_arms[i] > 0) "*" else ""
    paste0(selected_status$community[i], star, "\n", arms_wrapped)
  }, character(1))
  names(x_labels) <- selected_status$community

  community_family$community <- factor(community_family$community,
                                       levels = selected_communities)
  community_family$signal <- factor(community_family$signal,
                                    levels = rev(all_signals))
  community_family$support_arms <- as.numeric(community_family$support_arms)
  community_family$tile_label[is.na(community_family$tile_label)] <- ""
  community_family$display_label <- gsub("; ", "\n", community_family$tile_label,
                                         fixed = TRUE)

  ggplot(community_family,
         aes(x = community, y = signal, fill = support_arms)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = display_label),
              size = 2.85, lineheight = 0.9, color = "#1F2933") +
    scale_x_discrete(labels = x_labels, position = "top") +
    scale_fill_gradientn(
      colors = c("#F4F6F7", "#D9EDEA", "#75B9A9", "#256D85"),
      values = c(0, 0.15, 0.55, 1),
      limits = c(0, max(community_family$support_arms)),
      breaks = c(0, 2, 4, 6),
      name = "Support\narms"
    ) +
    labs(
      title = "Community/family map: support level and interval scope are separate",
      subtitle = "Tiles show selected HPRCv2 arm/community support. Special rows report community-specific gene counts.",
      x = NULL,
      y = NULL,
      caption = paste0(
        "* Community contains at least one member arm without a called CHM13 PHR interval: ",
        "C5 chr6_p; C7 chr13_p; C14 chrY_q; C15 chrY_p. ",
        "Blank tiles mean no selected family/support row in the canonical HPRCv2 tables."
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18.5, color = "#1F2933"),
      plot.subtitle = element_text(size = 11.5, color = "#52616B", margin = margin(b = 10)),
      plot.caption = element_text(size = 9.3, color = "#6B7280", hjust = 0),
      axis.text.x = element_text(size = 8.9, color = "#1F2933", lineheight = 0.9),
      axis.text.y = element_text(size = 10.5, color = "#1F2933"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(14, 18, 18, 18)
    )
}

save_gg <- function(plot, basename, width = 16, height = 9) {
  ggsave(file.path(out_dir, paste0(basename, ".png")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 200, bg = "white")
  ggsave(file.path(out_dir, paste0(basename, ".pdf")),
         plot = plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
}

render_method(file.path(out_dir, "copy_aware_method_comparison.png"), "png")
render_method(file.path(out_dir, "copy_aware_method_comparison.pdf"), "pdf")
save_gg(ranked_plot(ranked_signal_support), "ranked_copy_aware_gene_signals")
save_gg(community_plot(community_family_support, community_interval_status),
        "community_family_signal_map")

cat("Wrote gene enrichment figures and support TSVs to ", out_dir, "\n", sep = "")
