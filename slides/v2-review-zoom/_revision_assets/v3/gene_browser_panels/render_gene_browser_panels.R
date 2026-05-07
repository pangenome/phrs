#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

if (requireNamespace("ragg", quietly = TRUE)) {
  png_device <- ragg::agg_png
} else {
  png_device <- "png"
}

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

find_repo_root <- function(start_dir) {
  here <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "chm13.phrs.bed")) &&
        file.exists(file.path(here, "phrs.genes.gff3"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find repo root containing chm13.phrs.bed and phrs.genes.gff3")
    }
    here <- parent
  }
}

out_dir <- script_dir()
repo_root <- find_repo_root(out_dir)
inventory_dir <- file.path(
  repo_root,
  "slides/v2-review-zoom/_revision_assets/v3/gene_browser_inventory"
)

paths <- list(
  target_loci = file.path(inventory_dir, "target_loci.tsv"),
  inventory_readme = file.path(inventory_dir, "README.md"),
  inventory_track_schema = file.path(inventory_dir, "track_schema.tsv"),
  chm13_phr_bed = file.path(repo_root, "chm13.phrs.bed"),
  phr_gene_gff = file.path(repo_root, "phrs.genes.gff3"),
  hprc_all_vs_all = Sys.getenv(
    "HPRC_ALL_VS_ALL_TSV",
    "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
  ),
  chm13_repeatmasker = Sys.getenv(
    "CHM13_REPEATMASKER_BED",
    "/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz"
  ),
  d4z4_dux4l_by_community = Sys.getenv(
    "D4Z4_DUX4L_BY_COMMUNITY_TSV",
    "/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv"
  ),
  or4f_pseudogene_fraction = Sys.getenv(
    "OR4F_PSEUDOGENE_FRACTION_CSV",
    "/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv"
  ),
  community_tar1_by_arm = Sys.getenv(
    "COMMUNITY_TAR1_BY_ARM_TSV",
    "/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv"
  )
)

for (path in paths) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
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

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_coord <- function(chrom, start0, end0) {
  paste0(chrom, ":", fmt_int(start0), "-", fmt_int(end0))
}

parse_bed_coord <- function(x) {
  pieces <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(pieces) != 2) {
    stop("Could not parse coordinate: ", x)
  }
  span <- strsplit(pieces[[2]], "-", fixed = TRUE)[[1]]
  data.frame(
    chrom = pieces[[1]],
    start0 = as.integer(span[[1]]),
    end0 = as.integer(span[[2]]),
    stringsAsFactors = FALSE
  )
}

parse_self_chr <- function(x) {
  pieces <- strsplit(x, ":", fixed = TRUE)[[1]]
  span <- strsplit(pieces[[2]], "-", fixed = TRUE)[[1]]
  data.frame(
    chrom = pieces[[1]],
    start0 = as.integer(span[[1]]),
    end0 = as.integer(span[[2]]),
    stringsAsFactors = FALSE
  )
}

parse_attrs <- function(attrs, key) {
  vapply(attrs, function(a) {
    fields <- strsplit(a, ";", fixed = TRUE)[[1]]
    hit <- fields[grepl(paste0("^", key, "="), fields)]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", key, "="), "", hit[[1]])
  }, character(1))
}

read_gff_genes <- function(path) {
  gff <- read.delim(path,
                    sep = "\t",
                    header = FALSE,
                    quote = "",
                    comment.char = "#",
                    stringsAsFactors = FALSE)
  names(gff) <- c(
    "chrom", "source", "type", "start1", "end1", "score",
    "strand", "phase", "attrs"
  )
  gff <- gff[gff$type == "gene", ]
  gff$start0 <- as.integer(gff$start1) - 1L
  gff$end0 <- as.integer(gff$end1)
  gff$gene_name <- parse_attrs(gff$attrs, "gene_name")
  missing_name <- is.na(gff$gene_name) | gff$gene_name == ""
  gff$gene_name[missing_name] <- parse_attrs(gff$attrs[missing_name], "gene")
  gff$gene_biotype <- parse_attrs(gff$attrs, "gene_biotype")
  gff
}

read_phr_bed <- function(path) {
  bed <- read.delim(path,
                    sep = "\t",
                    header = FALSE,
                    quote = "",
                    comment.char = "",
                    stringsAsFactors = FALSE)
  names(bed)[1:4] <- c("chrom", "start0", "end0", "chrs_involved")
  bed
}

read_tar1_repeatmasker <- function(path) {
  cmd <- paste(
    "zcat",
    shQuote(path),
    "| awk -F '\\t' '$4 == \"TAR1\" {print}'"
  )
  tar1 <- read.delim(pipe(cmd),
                     sep = "\t",
                     header = FALSE,
                     quote = "",
                     comment.char = "",
                     stringsAsFactors = FALSE)
  if (nrow(tar1) == 0) {
    return(data.frame())
  }
  names(tar1) <- c(
    "raw_chrom", "start0", "end0", "repeat_name", "score", "strand",
    "repeat_class", "repeat_family", "divergence", "rm_id"
  )
  tar1$chrom <- sub("^CHM13#0#", "", tar1$raw_chrom)
  tar1
}

target_loci <- read_tsv(paths$target_loci)
phr_bed <- read_phr_bed(paths$chm13_phr_bed)
genes_all <- read_gff_genes(paths$phr_gene_gff)
tar1_all <- read_tar1_repeatmasker(paths$chm13_repeatmasker)
all_vs_all <- read_tsv(paths$hprc_all_vs_all)
d4z4_summary <- read_tsv(paths$d4z4_dux4l_by_community)
or4f_fraction <- read.csv(paths$or4f_pseudogene_fraction, stringsAsFactors = FALSE)
tar1_summary <- read_tsv(paths$community_tar1_by_arm)

target_by_id <- setNames(seq_len(nrow(target_loci)), target_loci$target_id)

target_view <- function(target_id, label, community, community_text) {
  row <- target_loci[target_by_id[[target_id]], ]
  coord <- parse_bed_coord(row$view_coordinates)
  data.frame(
    target_id = target_id,
    view_id = label,
    row_label = paste0(
      label, "\n",
      fmt_coord(coord$chrom, coord$start0, coord$end0)
    ),
    chrom = coord$chrom,
    start0 = coord$start0,
    end0 = coord$end0,
    phr_start0 = coord$start0,
    phr_end0 = coord$end0,
    community = community,
    community_text = community_text,
    coord_system = "CHM13 v2.0; BED 0-based half-open",
    stringsAsFactors = FALSE
  )
}

parse_acro_view <- function(seq_name) {
  row <- all_vs_all[all_vs_all$seq == seq_name, ]
  if (nrow(row) != 1) {
    stop("Expected one HPRCv2 row for ", seq_name, ", saw ", nrow(row))
  }
  self <- parse_self_chr(row$self_chr)
  arm_label <- sub(".*_(chr[0-9XY]+)_([pq])arm$", "\\1\\2", row$seq)
  data.frame(
    target_id = "acrocentric_c7_p_arm_group",
    view_id = arm_label,
    row_label = paste0(
      arm_label, "\n",
      row$self_chr, " + offsets"
    ),
    chrom = self$chrom,
    start0 = self$start0,
    end0 = self$end0,
    phr_start0 = self$start0 + as.integer(row$region_start),
    phr_end0 = self$start0 + as.integer(row$region_end),
    community = "C7",
    community_text = "C7 acrocentric p-arm group",
    coord_system = "HPRCv2 CHM13#0 PanSN relative; offsets are 0-based half-open",
    stringsAsFactors = FALSE
  )
}

view_sets <- list(
  dux4_d4z4 = rbind(
    target_view(
      "dux4_d4z4_c1_chr4q",
      "chr4q",
      "C1",
      "C1 chr4q/chr10q DUX4-D4Z4"
    ),
    target_view(
      "dux4_d4z4_c1_chr10q",
      "chr10q",
      "C1",
      "C1 chr4q/chr10q DUX4-D4Z4"
    )
  ),
  or4f_c3 = target_view(
    "or4f_c3_chr3q",
    "chr3q",
    "C3",
    "C3 OR4F-rich subtelomere"
  ),
  or4f_c8 = target_view(
    "or4f_decay_c8_chr15q",
    "chr15q",
    "C8",
    "C8 OR4F pseudogene endpoint"
  ),
  tar1_c2 = target_view(
    "tar1_rich_c2_chr18p",
    "chr18p",
    "C2",
    "C2 TAR1-rich chr18p"
  ),
  acro_c7 = do.call(rbind, lapply(c(
    "CHM13#0#chr13:2544-502543_chr13_parm",
    "CHM13#0#chr14:2075-502074_chr14_parm",
    "CHM13#0#chr15:3258-503257_chr15_parm",
    "CHM13#0#chr21:2505-502504_chr21_parm",
    "CHM13#0#chr22:4138-504137_chr22_parm"
  ), parse_acro_view))
)

label_patterns <- list(
  dux4_d4z4 = c("^FRG2", "^DBET$", "^DUX4$", "^RPL23AP"),
  or4f_c3 = c(
    "^OR4F", "^WASH8P", "^DDX11L8", "^FAM138D",
    "^SEPTIN14P22", "^GTF2IP17"
  ),
  or4f_c8 = c("^OR4F", "^FAM138E", "^WASH3P", "^DDX11L9"),
  tar1_c2 = c("^TUBB8P8", "^TUBB8B", "^IL9RP4"),
  acro_c7 = c("^MTCO3P39")
)

array_patterns <- list(
  dux4_d4z4 = c("^DUX4L"),
  or4f_c3 = character(0),
  or4f_c8 = character(0),
  tar1_c2 = character(0),
  acro_c7 = character(0)
)

community_colors <- c(
  C1 = "#8B6BB1",
  C2 = "#D95F02",
  C3 = "#1B9E77",
  C7 = "#2A9D8F",
  C8 = "#B44E5A"
)

family_colors <- c(
  "DUX4/DUX4L" = "#7B61A8",
  "OR4F" = "#008C8C",
  "TAR1" = "#C65B28",
  "IL9R/IL9RP" = "#B156A0",
  "TUBB8" = "#2F6DB3",
  "MTCO" = "#667D2E",
  "FRG/DBET" = "#A15C38",
  "WASH/DDX/FAM" = "#545454",
  "SEPTIN/GTF2I/RPL23A" = "#8A8A8A",
  "Other displayed" = "#707070"
)

gene_family <- function(name) {
  ifelse(grepl("^DUX4", name), "DUX4/DUX4L",
  ifelse(grepl("^OR4F", name), "OR4F",
  ifelse(grepl("^IL9R", name), "IL9R/IL9RP",
  ifelse(grepl("^TUBB8", name), "TUBB8",
  ifelse(grepl("^MTCO", name), "MTCO",
  ifelse(grepl("^FRG|^DBET$", name), "FRG/DBET",
  ifelse(grepl("^WASH|^DDX11L|^FAM138", name), "WASH/DDX/FAM",
  ifelse(grepl("^SEPTIN|^GTF2IP|^RPL23A", name),
         "SEPTIN/GTF2I/RPL23A",
         "Other displayed"))))))))
}

label_match <- function(x, patterns) {
  if (length(patterns) == 0) {
    return(rep(FALSE, length(x)))
  }
  Reduce(`|`, lapply(patterns, grepl, x = x))
}

clip_to_view <- function(dat, views) {
  if (nrow(dat) == 0) return(dat)
  out <- list()
  for (i in seq_len(nrow(views))) {
    view <- views[i, ]
    hit <- dat[
      dat$chrom == view$chrom &
        dat$end0 > view$start0 &
        dat$start0 < view$end0,
    ]
    if (nrow(hit) == 0) next
    hit$view_id <- view$view_id
    hit$row_label <- view$row_label
    hit$window_start0 <- view$start0
    hit$window_end0 <- view$end0
    hit$xstart <- pmax(hit$start0, view$start0) - view$start0
    hit$xend <- pmin(hit$end0, view$end0) - view$start0
    hit$xstart_kb <- hit$xstart / 1000
    hit$xend_kb <- hit$xend / 1000
    out[[length(out) + 1L]] <- hit
  }
  if (length(out) == 0) return(dat[0, ])
  do.call(rbind, out)
}

build_panel_data <- function(panel_key, views) {
  row_levels <- rev(unique(views$row_label))
  views$row_label <- factor(views$row_label, levels = row_levels)
  views$width_kb <- (views$end0 - views$start0) / 1000

  phr <- views
  phr$xstart_kb <- (phr$phr_start0 - phr$start0) / 1000
  phr$xend_kb <- (phr$phr_end0 - phr$start0) / 1000

  genes <- clip_to_view(genes_all, views)
  if (nrow(genes) > 0) {
    genes$display_gene <- label_match(genes$gene_name, label_patterns[[panel_key]])
    genes$array_gene <- label_match(genes$gene_name, array_patterns[[panel_key]])
    genes_display <- genes[genes$display_gene, ]
    genes_array <- genes[genes$array_gene, ]
  } else {
    genes_display <- genes
    genes_array <- genes
  }

  if (nrow(genes_display) > 0) {
    genes_display$family <- gene_family(genes_display$gene_name)
    genes_display$arrow_x <- ifelse(
      genes_display$strand == "+",
      genes_display$xstart_kb,
      genes_display$xend_kb
    )
    genes_display$arrow_xend <- ifelse(
      genes_display$strand == "+",
      genes_display$xend_kb,
      genes_display$xstart_kb
    )
    genes_display$label_x <- (genes_display$xstart_kb + genes_display$xend_kb) / 2
    genes_display$gene_y <- 2.05
    genes_display$label_y <- 2.35 + 0.24 * (seq_len(nrow(genes_display)) %% 3)
  }

  if (nrow(genes_array) > 0) {
    genes_array$family <- "DUX4/DUX4L"
  }

  tar1 <- clip_to_view(tar1_all, views)
  if (nrow(tar1) > 0) {
    tar1$family <- "TAR1"
  }

  list(
    views = views,
    phr = phr,
    genes_display = genes_display,
    genes_array = genes_array,
    tar1 = tar1
  )
}

span_labels <- function(dat, label_prefix) {
  parts <- split(dat, dat$row_label, drop = TRUE)
  do.call(rbind, lapply(parts, function(x) {
    data.frame(
      row_label = unique(x$row_label),
      xstart_kb = min(x$xstart_kb),
      xend_kb = max(x$xend_kb),
      x = 0.5 * (min(x$xstart_kb) + max(x$xend_kb)),
      n = nrow(x),
      label = paste0(label_prefix, " x", nrow(x)),
      stringsAsFactors = FALSE
    )
  }))
}

metric_text <- function(panel_key) {
  if (panel_key == "dux4_d4z4") {
    chm13 <- d4z4_summary[d4z4_summary$join_key == "CHM13#0" &
                            d4z4_summary$ChromArm %in% c("chr4_q", "chr10_q"), ]
    min_n <- aggregate(n_dux4l ~ ChromArm, data = chm13, min)
    max_n <- aggregate(n_dux4l ~ ChromArm, data = chm13, max)
    parts <- paste0(sub("_", "", min_n$ChromArm), ": ",
                    min_n$n_dux4l, "-", max_n$n_dux4l, " DUX4L rows")
    return(paste(
      "D4Z4 proxy: clustered DUX4L/DUX4 gene array.",
      paste(parts, collapse = "; ")
    ))
  }
  if (panel_key == "or4f_c3") {
    row <- or4f_fraction[or4f_fraction$chr_arm == "chr3_qarm", ]
    return(paste0(
      "OR4F arm summary: chr3_qarm has ",
      row$total,
      " OR4F annotations; pseudogene fraction ",
      round(100 * row$pseudo_frac, 1),
      "%."
    ))
  }
  if (panel_key == "or4f_c8") {
    row <- or4f_fraction[or4f_fraction$chr_arm == "chr15_qarm", ]
    return(paste0(
      "OR4F decay endpoint: chr15_qarm has ",
      row$total,
      " OR4F annotations; ",
      row$n_pseudo,
      " are pseudogene annotations (",
      round(100 * row$pseudo_frac, 1),
      "%)."
    ))
  }
  if (panel_key == "tar1_c2") {
    row <- tar1_summary[tar1_summary$chr_arm == "chr18_parm", ]
    return(paste0(
      "TAR1 summary: chr18_parm is ",
      row$pct_with_tar1,
      "% TAR1-positive; mean ",
      row$mean_tar1_count,
      " TAR1 copies per sequence."
    ))
  }
  if (panel_key == "acro_c7") {
    rows <- tar1_summary[tar1_summary$chr_arm %in% c(
      "chr13_parm", "chr14_parm", "chr15_parm", "chr21_parm", "chr22_parm"
    ), ]
    pct <- paste0(sub("chr", "", sub("_parm", "p", rows$chr_arm)),
                  " ", rows$pct_with_tar1, "%")
    return(paste0(
      "C7 acrocentric p-arms use HPRCv2 PanSN offsets; TAR1 prevalence: ",
      paste(pct, collapse = "; "),
      "."
    ))
  }
  ""
}

panel_specs <- data.frame(
  panel_key = c("dux4_d4z4", "or4f_c3", "or4f_c8", "tar1_c2", "acro_c7"),
  output_prefix = c(
    "panel_01_dux4_d4z4_c1_chr4_chr10",
    "panel_02_or4f_c3_chr3q",
    "panel_03_or4f_decay_c8_chr15q",
    "panel_04_tar1_c2_chr18p",
    "panel_05_acrocentric_c7_p_arm_group"
  ),
  title = c(
    "C1 recovers the DUX4/D4Z4 subtelomeric system",
    "C3 chr3q OR4F-rich PHR interval",
    "C8 chr15q OR4F pseudogene endpoint",
    "C2 chr18p TAR1-rich PHR interval",
    "C7 acrocentric p-arm PHR small multiples"
  ),
  subtitle = c(
    "Same grammar for chr4q and chr10q: community band, PHR block, curated genes, repeat/proxy lane.",
    "Focused gene lane highlights OR4F and neighboring subtelomeric duplicon markers without raw clutter.",
    "Paired OR4F example showing a singleton community with mostly pseudogene OR4F annotation.",
    "TAR1 is drawn as a repeat lane, separate from nearby TUBB8B and IL9RP4 genes.",
    "HPRCv2 PanSN relative offsets show the shared rDNA-adjacent acrocentric p-arm community."
  ),
  width = c(13.333, 13.333, 13.333, 13.333, 13.333),
  height = c(7.5, 7.5, 7.5, 7.5, 7.5),
  stringsAsFactors = FALSE
)

make_panel_plot <- function(panel_key, title, subtitle) {
  dat <- build_panel_data(panel_key, view_sets[[panel_key]])
  row_levels <- levels(dat$views$row_label)
  max_x <- max((dat$views$end0 - dat$views$start0) / 1000)
  xpad <- max_x * 0.012

  community <- dat$views
  community$xstart_kb <- 0
  community$xend_kb <- (community$end0 - community$start0) / 1000

  metric <- metric_text(panel_key)
  metric_df <- data.frame(
    row_label = factor(row_levels[[1]], levels = row_levels),
    x = max_x * 0.995,
    y = 4.72,
    label = metric,
    stringsAsFactors = FALSE
  )

  empty_repeats <- dat$views
  empty_repeats$label <- ifelse(
    dat$views$view_id %in% unique(dat$tar1$view_id),
    "",
    ifelse(panel_key == "dux4_d4z4", "D4Z4 proxy from DUX4L gene array", "no focal repeat block")
  )
  empty_repeats$x <- 0.5 * (empty_repeats$end0 - empty_repeats$start0) / 1000

  p <- ggplot() +
    geom_rect(
      data = community,
      aes(xmin = xstart_kb, xmax = xend_kb,
          ymin = 4.18, ymax = 4.52, fill = community),
      color = NA,
      alpha = 0.90
    ) +
    geom_text(
      data = community,
      aes(x = 5, y = 4.35, label = community_text),
      hjust = 0,
      size = 3.8,
      color = "white",
      fontface = "bold"
    ) +
    geom_rect(
      data = dat$phr,
      aes(xmin = xstart_kb, xmax = xend_kb, ymin = 3.18, ymax = 3.50),
      fill = "#294C60",
      color = "#1B2F3A",
      linewidth = 0.35
    ) +
    geom_text(
      data = dat$phr,
      aes(x = (xstart_kb + xend_kb) / 2, y = 3.34, label = "PHR"),
      color = "white",
      size = 3.2,
      fontface = "bold"
    ) +
    geom_segment(
      data = dat$genes_display,
      aes(x = arrow_x, xend = arrow_xend, y = gene_y, yend = gene_y,
          color = family),
      linewidth = 2.6,
      lineend = "round",
      arrow = arrow(length = unit(0.075, "inches"), type = "closed")
    ) +
    geom_text(
      data = dat$genes_display,
      aes(x = label_x, y = label_y, label = gene_name, color = family),
      size = 3.05,
      fontface = "bold",
      check_overlap = TRUE
    ) +
    geom_rect(
      data = dat$genes_array,
      aes(xmin = xstart_kb, xmax = xend_kb, ymin = 0.80, ymax = 1.10),
      fill = family_colors[["DUX4/DUX4L"]],
      color = NA,
      alpha = 0.88
    ) +
    geom_rect(
      data = dat$tar1,
      aes(xmin = xstart_kb, xmax = xend_kb, ymin = 0.78, ymax = 1.12),
      fill = family_colors[["TAR1"]],
      color = "#783000",
      linewidth = 0.25,
      alpha = 0.90
    ) +
    geom_text(
      data = empty_repeats,
      aes(x = x, y = 0.55, label = label),
      color = "#5B5B5B",
      size = 2.8
    ) +
    geom_text(
      data = metric_df,
      aes(x = x, y = y, label = label),
      hjust = 1,
      color = "#333333",
      size = 3.1
    ) +
    facet_grid(row_label ~ ., switch = "y") +
    scale_fill_manual(values = community_colors, drop = FALSE) +
    scale_color_manual(values = family_colors, drop = FALSE) +
    scale_x_continuous(
      name = "Distance from panel window start (kb); row labels give absolute coordinates",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(limits = c(0.35, 4.95), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(xlim = c(-xpad, max_x + xpad), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      y = NULL,
      caption = "Track order: community, PHR interval, curated genes, repeat/proxy. GFF3 genes are converted to 0-based half-open intervals."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 20, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 11.5, color = "#3E3E3E", margin = margin(b = 10)),
      plot.caption = element_text(size = 8.5, color = "#555555", hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#E2E2E2", linewidth = 0.30),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 1,
        face = "bold",
        size = 10,
        lineheight = 0.95,
        margin = margin(r = 8)
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 10, color = "#333333"),
      axis.text.x = element_text(size = 9, color = "#333333"),
      plot.margin = margin(t = 12, r = 20, b = 10, l = 18)
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = guide_legend(order = 2, override.aes = list(linewidth = 3))
    )

  if (nrow(dat$genes_array) > 0) {
    array_labels <- span_labels(dat$genes_array, "DUX4L array")
    p <- p + geom_text(
      data = array_labels,
      aes(x = x, y = 1.30, label = label),
      color = family_colors[["DUX4/DUX4L"]],
      size = 3.1,
      fontface = "bold"
    )
  }

  if (nrow(dat$tar1) > 0) {
    tar1_labels <- span_labels(dat$tar1, "TAR1")
    tar1_labels$x <- tar1_labels$xstart_kb
    p <- p + geom_text(
      data = tar1_labels,
      aes(x = x, y = 1.30, label = label),
      hjust = 0,
      color = family_colors[["TAR1"]],
      size = 3.1,
      fontface = "bold",
      check_overlap = TRUE
    )
  }

  p
}

manifest_rows <- list()

for (i in seq_len(nrow(panel_specs))) {
  spec <- panel_specs[i, ]
  message("Rendering ", spec$output_prefix)
  p <- make_panel_plot(spec$panel_key, spec$title, spec$subtitle)
  png_path <- file.path(out_dir, paste0(spec$output_prefix, ".png"))
  pdf_path <- file.path(out_dir, paste0(spec$output_prefix, ".pdf"))
  ggsave(
    filename = png_path,
    plot = p,
    width = spec$width,
    height = spec$height,
    dpi = 220,
    device = png_device,
    bg = "white"
  )
  ggsave(
    filename = pdf_path,
    plot = p,
    width = spec$width,
    height = spec$height,
    device = cairo_pdf,
    bg = "white"
  )

  views <- view_sets[[spec$panel_key]]
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    panel = spec$output_prefix,
    png = basename(png_path),
    pdf = basename(pdf_path),
    target_ids = paste(unique(views$target_id), collapse = ";"),
    rows = paste(views$view_id, collapse = ";"),
    coordinate_systems = paste(unique(views$coord_system), collapse = "; "),
    track_order = "community_band > phr_interval > curated_gene_models > repeat_or_proxy_markers",
    primary_message = spec$title,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
write.table(
  manifest,
  file = file.path(out_dir, "panel_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

render_schema <- data.frame(
  track_order = 1:4,
  track_name = c(
    "community_band",
    "phr_interval",
    "curated_gene_models",
    "repeat_or_proxy_markers"
  ),
  color_rule = c(
    "Community-specific fill using fixed C1/C2/C3/C7/C8 palette",
    "Dark blue filled block for the selected PHR interval",
    "Family/biotype-aware directional glyphs with compact labels",
    "TAR1 repeat blocks or DUX4L/D4Z4 proxy ticks; empty lane is still retained"
  ),
  coordinate_rule = c(
    "No genomic span; row-level arm community label",
    "0-based half-open interval, shown as kb from window start",
    "GFF3 1-based closed converted to 0-based half-open before plotting",
    "RepeatMasker BED is already 0-based half-open; DUX4L proxy uses converted gene rows"
  ),
  stringsAsFactors = FALSE
)

write.table(
  render_schema,
  file = file.path(out_dir, "render_track_schema.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

input_paths <- unname(unlist(paths))
display_paths <- ifelse(
  startsWith(input_paths, paste0(repo_root, "/")),
  sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", repo_root), "/"),
      "",
      input_paths),
  input_paths
)

input_manifest <- data.frame(
  input = names(paths),
  path = display_paths,
  exists = file.exists(input_paths),
  stringsAsFactors = FALSE
)

write.table(
  input_manifest,
  file = file.path(out_dir, "input_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote ", nrow(panel_specs), " panel PNG/PDF pairs to ", out_dir)
