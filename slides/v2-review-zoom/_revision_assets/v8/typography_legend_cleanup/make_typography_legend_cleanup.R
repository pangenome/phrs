#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 R package is required to render these figures.", call. = FALSE)
  }
  library(ggplot2)
})

script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
}

find_repo_root <- function(start_dir) {
  here <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(here, "slides/v2-review-zoom"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find repo root from ", start_dir, call. = FALSE)
    }
    here <- parent
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

write_tsv <- function(x, path) {
  write.table(x,
              file = path,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              na = "")
}

wrap_one <- function(x, width = 34) {
  paste(strwrap(x, width = width), collapse = "\n")
}

fmt_q <- function(x) {
  formatC(x, format = "f", digits = 3)
}

repo_root <- find_repo_root(script_dir())
out_dir <- file.path(repo_root, "slides/v2-review-zoom/_revision_assets/v8/typography_legend_cleanup")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

v5_dir <- file.path(repo_root, "slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_figures")
paths <- list(
  ranked = file.path(v5_dir, "ranked_signal_support.tsv"),
  community_map = file.path(v5_dir, "community_family_map_support.tsv"),
  fisher = Sys.getenv(
    "HPRCV2_FISHER_TSV",
    "/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv"
  ),
  mouse_dir = Sys.getenv(
    "MOUSE_T2T_50000BP_DIR",
    "/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp"
  )
)

missing_core <- unlist(paths[c("ranked", "community_map")])
missing_core <- missing_core[!file.exists(missing_core)]
if (length(missing_core) > 0) {
  stop("Missing required copied v5 support table(s):\n",
       paste(missing_core, collapse = "\n"), call. = FALSE)
}

save_gg <- function(plot, basename, width = 16, height = 9) {
  ggsave(file.path(out_dir, paste0(basename, ".png")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 200, bg = "white")
  ggsave(file.path(out_dir, paste0(basename, ".pdf")),
         plot = plot, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf, bg = "white")
}

draw_callout <- function(text, x0, y0, x1, y1, cex = 1.0,
                         fill = "#fff7e6", border = "#a16207") {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = 1.8)
  text((x0 + x1) / 2, (y0 + y1) / 2, text,
       cex = cex, font = 2, col = "#1f2933")
}

extract_mantel <- function(path) {
  d <- read_tsv(path)
  as.numeric(d$U_statistic[d$test == "mantel"][1])
}

render_slide12 <- function() {
  pair_path <- file.path(paths$mouse_dir, "zuo2021_zygotene_phr_pair_correlation.tsv")
  if (!file.exists(pair_path)) {
    stop("Missing mouse zygotene pair table: ", pair_path, call. = FALSE)
  }

  d <- read_tsv(pair_path)
  d$chr_a <- sub("_[pq]$", "", d$arm_a)
  d$chr_b <- sub("_[pq]$", "", d$arm_b)
  d <- d[d$chr_a != d$chr_b & !is.na(d$mean_jaccard) &
           !is.na(d$hic_contact) & d$hic_contact > 0, ]
  rho <- stats::cor(d$mean_jaccard, d$hic_contact, method = "spearman")
  pval <- suppressWarnings(stats::cor.test(d$mean_jaccard, d$hic_contact,
                                           method = "spearman"))$p.value

  stages <- data.frame(
    stage = c("leptotene", "zygotene", "pachytene", "diplotene"),
    short_stage = c("lepto", "zygo", "pachy", "diplo"),
    stringsAsFactors = FALSE
  )
  stages$rho <- vapply(stages$stage, function(s) {
    extract_mantel(file.path(paths$mouse_dir, paste0("zuo2021_", s, "_global_test.tsv")))
  }, numeric(1))
  write_tsv(stages, file.path(out_dir, "slide12_stage_mantel_rho.tsv"))

  render_one <- function(path, type = c("png", "pdf")) {
    type <- match.arg(type)
    if (type == "png") {
      png(path, width = 3200, height = 1800, res = 200, type = "cairo")
    } else {
      pdf(path, width = 16, height = 9, useDingbats = FALSE)
    }
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)

    layout(matrix(c(1, 2, 1, 3), 2, 2, byrow = TRUE),
           widths = c(0.68, 0.32), heights = c(0.46, 0.54))

    par(mar = c(5.4, 5.8, 3.8, 1.2), family = "sans")
    plot(d$mean_jaccard, d$hic_contact, log = "y", pch = 21,
         bg = grDevices::adjustcolor("#4f84bd", 0.62), col = "#244d78",
         xlab = "Mean PHR Jaccard similarity",
         ylab = "Hi-C contact in zygotene (log scale)",
         main = "Zygotene pairs: sequence-similar subtelomeres contact more",
         cex.main = 1.35, cex.lab = 1.22, cex.axis = 1.08)
    fit <- stats::lm(log10(hic_contact) ~ mean_jaccard, data = d)
    xs <- seq(min(d$mean_jaccard), max(d$mean_jaccard), length.out = 100)
    lines(xs, 10 ^ stats::predict(fit, newdata = data.frame(mean_jaccard = xs)),
          col = "#111827", lwd = 2.3)
    legend("topleft",
           legend = sprintf("Spearman rho = %.3f\np = %.1e\nn = %d pairs",
                            rho, pval, nrow(d)),
           bty = "n", cex = 1.05, text.col = "#1f2933")

    par(mar = c(1.0, 1.0, 3.8, 1.2), family = "sans")
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    title("Talk-speed framing", cex.main = 1.30, font.main = 2, col.main = "#111827")
    text(0.03, 0.82,
         "Show the zygotene scatter\nbecause this is the bouquet\nstage: telomeres cluster at\nthe nuclear envelope.",
         adj = c(0, 1), cex = 1.05, col = "#1f2933")
    draw_callout("One sentence on slide:\nZygotene is the bouquet stage,\nwhere telomere clustering makes\n3D proximity most informative.",
                 0.04, 0.08, 0.96, 0.43, cex = 0.94)

    par(mar = c(5.4, 4.6, 3.2, 1.2), family = "sans")
    plot(seq_along(stages$stage), stages$rho, type = "b", pch = 21,
         bg = ifelse(stages$stage == "zygotene", "#d62728", "#2b7bbb"),
         col = "#1f2933", lwd = 2.0, cex = 1.9,
         xaxt = "n", xlim = c(0.85, 4.25), ylim = c(0.54, 0.75),
         xlab = "meiotic prophase stage",
         ylab = "Mantel rho",
         main = "Stage trajectory", cex.main = 1.20,
         cex.lab = 1.08, cex.axis = 0.98)
    axis(1, at = seq_along(stages$stage),
         labels = stages$short_stage, cex.axis = 1.00)
    text(seq_along(stages$stage), stages$rho + 0.012,
         sprintf("%.3f", stages$rho), cex = 0.95, col = "#111827")
    text(2, 0.742, "bouquet", cex = 1.02, font = 2, col = "#d62728")
  }

  render_one(file.path(out_dir, "slide12_mouse_zygotene_large_text.png"), "png")
  render_one(file.path(out_dir, "slide12_mouse_zygotene_large_text.pdf"), "pdf")
}

ranked_plot <- function(ranked) {
  ranked$rank <- as.integer(ranked$rank)
  ranked$plotted_support_arms <- as.numeric(ranked$plotted_support_arms)
  ranked$display_signal <- c(
    "OR4F pseudogenization gradient",
    "OR/OR4F community presence",
    "D4Z4/DUX4L C1 block",
    "Duplicon backbone genes",
    "C5 DDX11L/WASH/FAM138",
    "GTPBP6/IQSEC3 anchors",
    "C15 PAR1 coding outlier",
    "C7 acrocentric MTCO"
  )[match(ranked$rank, 1:8)]
  ranked$value_label <- c(
    "16 annotated arms",
    "14 arms",
    "2 C1 arms",
    "22 arms",
    "4 C5 arms",
    "3 support arms",
    "2 PAR1 arms",
    "5 acrocentric p-arms"
  )[match(ranked$rank, 1:8)]
  ranked <- ranked[order(ranked$plotted_support_arms, ranked$rank), ]
  ranked$display_signal <- factor(ranked$display_signal, levels = ranked$display_signal)

  ggplot(ranked, aes(x = display_signal, y = plotted_support_arms)) +
    geom_col(width = 0.58, fill = "#2f7f95", color = "white", linewidth = 0.4) +
    geom_text(aes(label = value_label), hjust = -0.08, size = 6.1,
              color = "#1f2933") +
    coord_flip(clip = "off") +
    scale_y_continuous(limits = c(0, 25), breaks = seq(0, 24, by = 4),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title = "Candidate gene-family signals by support count",
      subtitle = "Support units, not q-values. No bar is a BH-significant enrichment call.",
      x = NULL,
      y = "support count",
      caption = "Canonical Fisher screen: 116 family-community rows, 0 BH-significant rows; best q = 0.071."
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 29, color = "#111827"),
      plot.subtitle = element_text(size = 18, color = "#374151",
                                   margin = margin(b = 18)),
      plot.caption = element_text(size = 15, color = "#4b5563", hjust = 0),
      axis.text.y = element_text(size = 18, color = "#111827"),
      axis.text.x = element_text(size = 16, color = "#374151"),
      axis.title.x = element_text(size = 17, color = "#1f2933",
                                  margin = margin(t = 10)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(22, 150, 28, 30)
    )
}

make_simplified_map <- function(community_family) {
  selected_communities <- c("C1", "C3", "C5", "C7", "C11", "C12", "C15")
  community_labels <- c(
    C1 = "C1\n10q 4q",
    C3 = "C3\n11p 16q 19p\n3q 7p 9q",
    C5 = "C5*\n12p 20q 6p 9p",
    C7 = "C7*\n13p 14p 15p\n21p 22p",
    C11 = "C11\n1p 5q 6q 8p",
    C12 = "C12\n20p 2q",
    C15 = "C15*\nXp Yp"
  )

  family_rows <- data.frame(
    source_signal = c("OR", "RPL", "SEPTIN", "DDX11L"),
    display_signal = c("OR", "RPL", "SEPTIN", "DDX11L"),
    stringsAsFactors = FALSE
  )

  rows <- list()
  row_i <- 1L
  add_row <- function(comm, signal, support_arms, label) {
    rows[[row_i]] <<- data.frame(
      community = comm,
      community_label = community_labels[[comm]],
      signal = signal,
      support_arms = support_arms,
      label = label,
      stringsAsFactors = FALSE
    )
    row_i <<- row_i + 1L
  }

  for (comm in selected_communities) {
    for (i in seq_len(nrow(family_rows))) {
      src <- family_rows$source_signal[i]
      hit <- community_family[community_family$community == comm &
                                community_family$signal == src, , drop = FALSE]
      support <- if (nrow(hit) == 0) 0 else as.numeric(hit$support_arms[1])
      add_row(comm, family_rows$display_signal[i], support,
              if (support > 0) paste0(support, " arms") else "")
    }

    wash <- community_family[community_family$community == comm &
                               community_family$signal == "WASH", , drop = FALSE]
    fam <- community_family[community_family$community == comm &
                              community_family$signal == "FAM138", , drop = FALSE]
    wash_n <- if (nrow(wash) == 0) 0 else as.numeric(wash$support_arms[1])
    fam_n <- if (nrow(fam) == 0) 0 else as.numeric(fam$support_arms[1])
    support <- max(wash_n, fam_n)
    add_row(comm, "WASH/FAM138", support,
            if (support > 0) paste0(support, " arms") else "")

    anchor_signal <- switch(comm,
      C1 = "DUX4L specific",
      C7 = "MTCO specific",
      C15 = "PAR1 coding anchors",
      NA
    )
    if (is.na(anchor_signal)) {
      add_row(comm, "specific anchor", 0, "")
    } else {
      hit <- community_family[community_family$community == comm &
                                community_family$signal == anchor_signal, , drop = FALSE]
      support <- if (nrow(hit) == 0) 0 else as.numeric(hit$support_arms[1])
      label <- if (nrow(hit) == 0 || support == 0) "" else switch(comm,
        C1 = "DUX4L\n22 genes\n2 arms",
        C7 = "MTCO\n9 genes\n5 arms",
        C15 = "PAR1\n5 coding\n2 arms"
      )
      add_row(comm, "specific anchor", support, label)
    }
  }

  simplified <- do.call(rbind, rows)
  simplified$community <- factor(simplified$community, levels = selected_communities)
  simplified$community_label <- factor(
    simplified$community_label,
    levels = unname(community_labels[selected_communities])
  )
  simplified$signal <- factor(
    simplified$signal,
    levels = rev(c("OR", "RPL", "SEPTIN", "DDX11L", "WASH/FAM138", "specific anchor"))
  )
  simplified_out <- simplified
  simplified_out$label[simplified_out$label == ""] <- "."
  for (nm in names(simplified_out)) {
    if (is.character(simplified_out[[nm]]) || is.factor(simplified_out[[nm]])) {
      simplified_out[[nm]] <- gsub("\n", " / ", as.character(simplified_out[[nm]]),
                                   fixed = TRUE)
    }
  }
  write_tsv(simplified_out, file.path(out_dir, "slide14c_simplified_map_support.tsv"))
  simplified
}

community_plot <- function(simplified) {
  ggplot(simplified, aes(x = community_label, y = signal, fill = support_arms)) +
    geom_tile(color = "white", linewidth = 1.0) +
    geom_text(aes(label = label), size = 5.4, lineheight = 0.88,
              color = "#111827") +
    scale_fill_gradientn(
      colors = c("#f3f4f6", "#dbeee9", "#8fc8bc", "#2d7589"),
      values = c(0, 0.15, 0.55, 1),
      limits = c(0, max(simplified$support_arms)),
      guide = "none"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title = "Talk-ready community/family support map",
      subtitle = "Condensed to recurring families plus three anchors; tile text gives support directly.",
      x = NULL,
      y = NULL,
      caption = "* C5/C7/C15 include member arms without called CHM13 PHR intervals. Presence support; not BH-significant enrichment."
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 29, color = "#111827"),
      plot.subtitle = element_text(size = 17, color = "#374151",
                                   margin = margin(b = 14)),
      plot.caption = element_text(size = 14.5, color = "#4b5563", hjust = 0),
      axis.text.x = element_text(size = 16.5, color = "#111827", lineheight = 0.88),
      axis.text.y = element_text(size = 18, color = "#111827"),
      panel.grid = element_blank(),
      plot.margin = margin(20, 30, 28, 30)
    )
}

render_slide14_assets <- function() {
  ranked <- read_tsv(paths$ranked)
  write_tsv(ranked, file.path(out_dir, "slide14b_ranked_signal_support_source.tsv"))
  save_gg(ranked_plot(ranked), "slide14b_candidate_signals_talk_ready")

  community_family <- read_tsv(paths$community_map)
  simplified <- make_simplified_map(community_family)
  save_gg(community_plot(simplified), "slide14c_community_family_map_talk_ready")
}

render_slide13b <- function() {
  crop_script <- file.path(out_dir, "crop_png_top.py")
  source_png <- file.path(repo_root, "slides/v2-review-zoom/_typst/assets/s13_pedigree_bottom.png")
  output_png <- file.path(out_dir, "slide13b_pedigree_bottom_no_unused_legend.png")
  if (!file.exists(crop_script)) {
    stop("Missing deterministic PNG crop helper: ", crop_script, call. = FALSE)
  }
  if (!file.exists(source_png)) {
    stop("Missing slide 13b source PNG: ", source_png, call. = FALSE)
  }
  status <- system2(
    "python3",
    args = c(crop_script, source_png, output_png, "--height", "1495")
  )
  if (!identical(status, 0L)) {
    stop("slide 13b crop helper failed with exit status ", status, call. = FALSE)
  }
}

render_manifest <- function() {
  manifest <- data.frame(
    artifact = c(
      "slide12_mouse_zygotene_large_text.png",
      "slide12_mouse_zygotene_large_text.pdf",
      "slide12_stage_mantel_rho.tsv",
      "slide13b_pedigree_bottom_no_unused_legend.png",
      "crop_png_top.py",
      "slide13b_remove_unused_legend.typ",
      "slide14b_candidate_signals_talk_ready.png",
      "slide14b_candidate_signals_talk_ready.pdf",
      "slide14b_ranked_signal_support_source.tsv",
      "slide14c_community_family_map_talk_ready.png",
      "slide14c_community_family_map_talk_ready.pdf",
      "slide14c_simplified_map_support.tsv"
    ),
    purpose = c(
      "Slide 12 replacement asset with larger labels and preserved scatter plus stage trajectory.",
      "Vector/PDF companion for slide 12 replacement asset.",
      "Source stage Mantel rho values used in the slide 12 trajectory inset.",
      "Slide 13b materialized crop with the unused bottom legend removed.",
      "Deterministic PNG top-crop helper used because no pedigree regeneration recipe is available.",
      "Typst replacement snippet for slide 13b using the materialized no-legend crop.",
      "Slide 14b replacement asset with condensed labels, no legend, and doubled text.",
      "Vector/PDF companion for slide 14b replacement asset.",
      "Copied source support table from the v5 gene-enrichment figures.",
      "Slide 14c replacement asset with condensed map, direct labels, and no color legend.",
      "Vector/PDF companion for slide 14c replacement asset.",
      "Simplified support table used to render slide 14c."
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(out_dir, "asset_manifest.tsv"))
}

render_slide12()
render_slide13b()
render_slide14_assets()
render_manifest()

cat("Wrote typography/legend cleanup assets to ", out_dir, "\n", sep = "")
