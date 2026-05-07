#!/usr/bin/env Rscript

out_dir <- "slides/v2-review-zoom/_revision_assets/v4/10a_xaxis_orientation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mat_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv"
v2_comm_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_hic.arm-leiden.communities.tsv"
seq_comm_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv"
global_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv"
v3_asset <- "slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/candidate_10a_axis_box_fix.png"
v3_generator <- "slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/make_10a_axis_box_fix.R"
v3_audit <- "slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/matrix_order_audit.tsv"
v2_generator <- "slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R"
fig3_generator <- "paper_prep/figures/fig3/figure_fig3.R"

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE)
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

assert_true <- function(value, label) {
  if (!isTRUE(value)) fail("Assertion failed: ", label)
}

base_arm <- function(x) {
  sub("_(MATERNAL|PATERNAL|HAP1|HAP2)_", "_", x)
}

normalize_seq_arm <- function(x) {
  sub("arm$", "", x)
}

community_number <- function(x) {
  as.integer(sub("^[^0-9]*", "", x))
}

community_prefix <- function(x) {
  sub("[0-9].*$", "", x)
}

community_sort_key <- function(x) {
  sprintf("%s%05d", community_prefix(x), community_number(x))
}

chrom_rank <- function(arm) {
  chrom <- sub("_.*$", "", base_arm(arm))
  raw <- sub("^chr", "", chrom)
  out <- suppressWarnings(as.integer(raw))
  out[raw == "X"] <- 23L
  out[raw == "Y"] <- 24L
  out
}

arm_side_rank <- function(arm) {
  ifelse(grepl("_p$", arm), 1L, 2L)
}

hap_rank <- function(arm) {
  ifelse(grepl("_MATERNAL_", arm), 1L,
         ifelse(grepl("_PATERNAL_", arm), 2L,
                ifelse(grepl("_HAP1_", arm), 3L,
                       ifelse(grepl("_HAP2_", arm), 4L, 5L))))
}

format_sci <- function(x, digits = 2) {
  formatC(x, format = "e", digits = digits)
}

format_bool <- function(x) {
  ifelse(isTRUE(x), "TRUE", "FALSE")
}

detect_png_dims <- function(path) {
  if (!file.exists(path)) return("missing")
  if (nzchar(Sys.which("identify"))) {
    out <- tryCatch(
      system2("identify", c("-format", "%wx%h", path), stdout = TRUE,
              stderr = FALSE),
      error = function(e) character()
    )
    if (length(out) == 1 && nzchar(out)) return(out)
  }
  if (nzchar(Sys.which("file"))) {
    out <- tryCatch(system2("file", path, stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
    match <- regmatches(out, regexpr("[0-9]+ x [0-9]+", out))
    if (length(match) == 1 && nzchar(match)) return(gsub(" ", "", match))
  }
  "exists; dimensions not checked"
}

mat_df <- read_tsv(mat_path)
row_arms <- mat_df[[1]]
col_arms <- colnames(mat_df)[-1]
contact <- as.matrix(mat_df[, -1, drop = FALSE])
storage.mode(contact) <- "numeric"
rownames(contact) <- row_arms
colnames(contact) <- col_arms

assert_true(nrow(contact) == ncol(contact), "source matrix is square")
assert_true(identical(row_arms, col_arms),
            "source row names equal source column names")
assert_true(anyDuplicated(row_arms) == 0,
            "source row names are unique")
assert_true(anyDuplicated(col_arms) == 0,
            "source column names are unique")
assert_true(!anyNA(contact), "source matrix has no NA values")
sym_delta <- max(abs(contact - t(contact)))
assert_true(isTRUE(all.equal(contact, t(contact), tolerance = 1e-12)),
            "source matrix is symmetric")

seq_comm <- read_tsv(seq_comm_path)
seq_comm$base_arm <- normalize_seq_arm(seq_comm$arm)
seq_comm_by_base <- setNames(seq_comm$community, seq_comm$base_arm)

matrix_base_arms <- base_arm(row_arms)
seq_community <- unname(seq_comm_by_base[matrix_base_arms])
assert_true(!anyNA(seq_community),
            "all source matrix arm-haplotypes map to arm-level sequence communities")

v2_comm <- read_tsv(v2_comm_path)
v2_comm_by_arm <- setNames(v2_comm$community, v2_comm$arm)
v2_community <- unname(v2_comm_by_arm[row_arms])
assert_true(!anyNA(v2_community),
            "all source matrix arm-haplotypes map to v2 contact-community table")

global <- read_tsv(global_path)
wb <- global[global$test == "within_vs_between", , drop = FALSE]
assert_true(nrow(wb) == 1, "source global TSV has one within_vs_between row")

upper <- upper.tri(contact)
same_sequence_community <- outer(seq_community, seq_community, "==")
within_values <- contact[upper & same_sequence_community]
between_values <- contact[upper & !same_sequence_community]
computed_within <- mean(within_values)
computed_between <- mean(between_values)
computed_bw <- computed_between / computed_within
source_bw <- wb$between_mean / wb$within_mean

assert_true(length(within_values) == wb$n_within,
            "computed within pair count matches source global TSV")
assert_true(length(between_values) == wb$n_between,
            "computed between pair count matches source global TSV")
assert_true(isTRUE(all.equal(computed_within, wb$within_mean, tolerance = 1e-14)),
            "computed within mean matches source global TSV")
assert_true(isTRUE(all.equal(computed_between, wb$between_mean, tolerance = 1e-14)),
            "computed between mean matches source global TSV")
assert_true(isTRUE(all.equal(computed_bw, source_bw, tolerance = 1e-14)),
            "computed B/W matches source global TSV")

order_df <- data.frame(
  source_index = seq_along(row_arms),
  arm_haplotype = row_arms,
  base_arm = matrix_base_arms,
  sequence_community = seq_community,
  v2_contact_community = v2_community,
  stringsAsFactors = FALSE
)

ord <- with(order_df, order(
  community_sort_key(sequence_community),
  chrom_rank(arm_haplotype),
  arm_side_rank(arm_haplotype),
  base_arm,
  hap_rank(arm_haplotype),
  arm_haplotype
))

ordered_arms <- order_df$arm_haplotype[ord]
ordered_contact <- contact[ordered_arms, ordered_arms, drop = FALSE]
ordered_community <- order_df$sequence_community[ord]
ordered_base <- order_df$base_arm[ord]
ordered_community_by_arm <- setNames(ordered_community, ordered_arms)
n <- nrow(ordered_contact)

assert_true(identical(rownames(ordered_contact), colnames(ordered_contact)),
            "ordered row names equal ordered column names")
assert_true(identical(rownames(ordered_contact), ordered_arms),
            "ordered matrix row names equal ordered arm list")
assert_true(identical(colnames(ordered_contact), ordered_arms),
            "ordered matrix column names equal ordered arm list")

v3_implicit_display_col_order <- rev(ordered_arms)
corrected_display_row_order <- ordered_arms
corrected_display_col_order <- rev(v3_implicit_display_col_order)
assert_true(identical(corrected_display_col_order, ordered_arms),
            "corrected displayed X order is first-to-last ordered arms")
assert_true(!identical(v3_implicit_display_col_order, corrected_display_col_order),
            "v3 implicit displayed X order differs from corrected order")

display_contact <- ordered_contact[
  corrected_display_row_order,
  corrected_display_col_order,
  drop = FALSE
]

ordered_table <- data.frame(
  display_index = seq_along(ordered_arms),
  row_order_top_to_bottom = corrected_display_row_order,
  corrected_column_order_left_to_right = corrected_display_col_order,
  v3_implicit_column_order_left_to_right = v3_implicit_display_col_order,
  base_arm = ordered_base,
  sequence_community = ordered_community,
  source_row_index = match(corrected_display_row_order, row_arms),
  source_column_index = match(corrected_display_col_order, col_arms),
  stringsAsFactors = FALSE
)

write.table(ordered_table,
            file = file.path(out_dir, "ordered_arm_haplotypes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

make_blocks <- function(arms) {
  communities <- unname(ordered_community_by_arm[arms])
  rl <- rle(as.character(communities))
  ends <- cumsum(rl$lengths)
  starts <- c(1, head(ends, -1) + 1)
  data.frame(
    sequence_community = rl$values,
    start_index = starts,
    end_index = ends,
    n_arm_haplotypes = rl$lengths,
    stringsAsFactors = FALSE
  )
}

row_blocks <- make_blocks(corrected_display_row_order)
col_blocks <- make_blocks(corrected_display_col_order)
assert_true(identical(row_blocks$sequence_community, col_blocks$sequence_community),
            "row and corrected column blocks have matching communities")
assert_true(identical(row_blocks$start_index, col_blocks$start_index),
            "row and corrected column blocks have matching starts")
assert_true(identical(row_blocks$end_index, col_blocks$end_index),
            "row and corrected column blocks have matching ends")

blocks <- data.frame(
  sequence_community = row_blocks$sequence_community,
  row_start_index = row_blocks$start_index,
  row_end_index = row_blocks$end_index,
  column_start_index = col_blocks$start_index,
  column_end_index = col_blocks$end_index,
  n_arm_haplotypes = row_blocks$n_arm_haplotypes,
  x_min = col_blocks$start_index - 0.5,
  x_max = col_blocks$end_index + 0.5,
  y_min = n - row_blocks$end_index + 0.5,
  y_max = n - row_blocks$start_index + 1.5,
  stringsAsFactors = FALSE
)
blocks$members <- vapply(seq_len(nrow(blocks)), function(i) {
  paste(corrected_display_row_order[blocks$row_start_index[i]:blocks$row_end_index[i]],
        collapse = ";")
}, character(1))

write.table(blocks,
            file = file.path(out_dir, "sequence_community_boxes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

v3_mirror_basis <- paste(
  "v3 drew rasterImage(as.raster(t(color_matrix)[, n:1]));",
  "with rasterImage raster rows are already top-to-bottom and columns are left-to-right,",
  "so the symmetric matrix could look plausible while the displayed X axis mapped left-to-right to reverse(ordered_arms)."
)
corrected_x_axis_policy <- paste(
  "display columns left-to-right as ordered_arms[1:n], matching the sequence-community order;",
  "draw rasterImage(as.raster(color_matrix)) directly and compute x boxes in this displayed column coordinate system."
)
corrected_row_axis_policy <- paste(
  "display rows top-to-bottom as ordered_arms[1:n];",
  "compute y boxes with y_min = n - row_end + 0.5 and y_max = n - row_start + 1.5."
)

orientation_audit <- data.frame(
  display_index = seq_along(ordered_arms),
  row_order_top_to_bottom = corrected_display_row_order,
  row_sequence_community = unname(ordered_community_by_arm[corrected_display_row_order]),
  displayed_column_order_left_to_right = corrected_display_col_order,
  displayed_column_sequence_community = unname(ordered_community_by_arm[corrected_display_col_order]),
  v3_implicit_displayed_column_order_left_to_right = v3_implicit_display_col_order,
  v3_implicit_displayed_column_sequence_community = unname(ordered_community_by_arm[v3_implicit_display_col_order]),
  v3_x_axis_mirrored = "TRUE",
  orientation_summary = paste0(
    "v3_x_axis_mirrored=TRUE; corrected_x_axis_policy=",
    corrected_x_axis_policy
  ),
  corrected_x_axis_policy = corrected_x_axis_policy,
  corrected_row_axis_policy = corrected_row_axis_policy,
  v3_mirror_basis = v3_mirror_basis,
  stringsAsFactors = FALSE
)

write.table(orientation_audit,
            file = file.path(out_dir, "orientation_audit.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

same_v2_seq_partition <- function(a, b) {
  upper_mask <- upper.tri(matrix(FALSE, length(a), length(a)))
  same_a <- outer(a, a, "==")[upper_mask]
  same_b <- outer(b, b, "==")[upper_mask]
  sum(same_a == same_b)
}

pair_total <- choose(length(row_arms), 2)
v2_seq_pair_agree <- same_v2_seq_partition(v2_community, seq_community)

audit <- data.frame(
  check = c(
    "source_matrix_rows",
    "source_matrix_columns",
    "source_row_names_equal_column_names",
    "source_row_name_set_equal_column_name_set",
    "source_matrix_symmetric",
    "source_matrix_max_abs_symmetry_delta",
    "source_matrix_all_arms_have_sequence_community",
    "source_sequence_community_table",
    "v3_generator",
    "v3_asset_dimensions",
    "v3_audit",
    "v3_visual_x_axis_mirrored",
    "v3_implicit_displayed_column_order_left_to_right",
    "corrected_displayed_column_order_left_to_right",
    "corrected_x_axis_policy",
    "corrected_raster_transform",
    "original_fig3_panel_a_display_transform",
    "v2_redesign_generator_inspected",
    "v2_contact_community_table",
    "v2_contact_community_count",
    "sequence_community_count_in_matrix",
    "v2_contact_partition_pair_agreement_with_sequence_partition",
    "ordered_row_names_equal_ordered_column_names",
    "computed_within_mean_matches_source_tsv",
    "computed_between_mean_matches_source_tsv",
    "computed_bw_matches_source_tsv",
    "source_bw",
    "source_p_value"
  ),
  value = c(
    as.character(nrow(contact)),
    as.character(ncol(contact)),
    format_bool(identical(row_arms, col_arms)),
    format_bool(setequal(row_arms, col_arms)),
    "TRUE",
    format(sym_delta, scientific = TRUE, digits = 6),
    "TRUE",
    seq_comm_path,
    v3_generator,
    detect_png_dims(v3_asset),
    v3_audit,
    "TRUE",
    "reverse(ordered_arms); leftmost v3 displayed column was the last ordered arm-haplotype",
    "ordered_arms; leftmost v4 displayed column is the first ordered arm-haplotype",
    corrected_x_axis_policy,
    "rasterImage(as.raster(color_matrix)); no t(color_matrix)[, n:1] raster transform",
    "paper Fig. 3 panel A uses base image(..., t(vals_norm)[, n:1]); that transform is correct for image() but not for rasterImage()",
    v2_generator,
    v2_comm_path,
    as.character(length(unique(v2_community))),
    as.character(length(unique(seq_community))),
    paste0(v2_seq_pair_agree, "/", pair_total),
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    format(source_bw, scientific = FALSE, digits = 8),
    format_sci(wb$p_value, digits = 3)
  ),
  stringsAsFactors = FALSE
)

write.table(audit, file = file.path(out_dir, "matrix_order_audit.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

render_plot <- function(path, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(path, width = 1800, height = 1800, res = 180, type = "cairo")
  } else {
    pdf(path, width = 10, height = 10, onefile = TRUE, useDingbats = FALSE)
  }

  vals <- display_contact
  pos <- vals[vals > 0]
  floor_v <- if (length(pos)) stats::quantile(pos, 0.05) else 1e-6
  vals[vals == 0] <- floor_v
  log_vals <- log10(vals)
  zlim <- range(log_vals, finite = TRUE)
  norm <- (log_vals - zlim[1]) / diff(zlim)
  pal <- grDevices::colorRampPalette(
    c("white", "#fee0d2", "#fc9272", "#de2d26", "#7f0000"))(100)
  color_matrix <- matrix(
    pal[pmax(1, pmin(100, floor(norm * 99) + 1))],
    nrow = nrow(norm),
    ncol = ncol(norm),
    dimnames = dimnames(norm)
  )

  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)

  layout(matrix(c(1, 2), 1, 2), widths = c(0.82, 0.18))
  par(mar = c(7.4, 7.2, 6.4, 1.0), pty = "s", family = "sans")
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
       xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)

  rasterImage(as.raster(color_matrix), 0.5, 0.5, n + 0.5, n + 0.5,
              interpolate = FALSE)
  box(col = "#333333", lwd = 1.0)

  for (b in head(blocks$column_end_index, -1)) {
    abline(v = b + 0.5, col = grDevices::adjustcolor("#333333", alpha.f = 0.25),
           lwd = 0.55)
  }
  for (b in head(blocks$row_end_index, -1)) {
    abline(h = n - b + 0.5, col = grDevices::adjustcolor("#333333", alpha.f = 0.25),
           lwd = 0.55)
  }

  for (i in seq_len(nrow(blocks))) {
    rect(blocks$x_min[i], blocks$y_min[i], blocks$x_max[i], blocks$y_max[i],
         border = "#0b3c78", lwd = 1.35)
    if (blocks$n_arm_haplotypes[i] >= 4) {
      mid_x <- (blocks$column_start_index[i] + blocks$column_end_index[i]) / 2
      mid_y <- n - ((blocks$row_start_index[i] + blocks$row_end_index[i]) / 2) + 1
      text(mid_x, mid_y, blocks$sequence_community[i],
           cex = 0.62, col = "#0b3c78", font = 2)
    }
  }

  block_mid_x <- (blocks$column_start_index + blocks$column_end_index) / 2
  block_mid_y <- n - ((blocks$row_start_index + blocks$row_end_index) / 2) + 1
  axis(1, at = block_mid_x, labels = blocks$sequence_community,
       las = 2, cex.axis = 0.62, tick = FALSE, line = 0.1)
  axis(2, at = block_mid_y, labels = blocks$sequence_community,
       las = 1, cex.axis = 0.62, tick = FALSE, line = 0.1)
  axis(4, at = block_mid_y, labels = blocks$sequence_community,
       las = 1, cex.axis = 0.62, tick = FALSE, line = 0.1)

  first_col <- corrected_display_col_order[1]
  last_col <- corrected_display_col_order[n]
  first_row <- corrected_display_row_order[1]
  last_row <- corrected_display_row_order[n]
  first_col_comm <- ordered_community_by_arm[first_col]
  last_col_comm <- ordered_community_by_arm[last_col]
  first_row_comm <- ordered_community_by_arm[first_row]
  last_row_comm <- ordered_community_by_arm[last_row]

  title("HG002 Pore-C contacts ordered by sequence community",
        cex.main = 1.10, font.main = 2, line = 4.4)
  mtext("v4 orientation fix: X left = first ordered column; X right = last ordered column",
        side = 3, line = 2.5, cex = 0.78, col = "#333333")
  mtext("blue boxes = arm-level sequence communities C1-C15 in displayed row/column coordinates",
        side = 3, line = 1.4, cex = 0.76, col = "#0b3c78")
  mtext(sprintf("B/W = %.3f    p = %.1e    source within/between TSV checked",
                source_bw, wb$p_value),
        side = 1, line = 3.1, cex = 0.76, col = "#333333")
  mtext(sprintf("X left-to-right: FIRST %s (%s) -> LAST %s (%s)",
                first_col, first_col_comm, last_col, last_col_comm),
        side = 1, line = 4.5, cex = 0.67, col = "#333333")
  mtext("Corrected raster policy: draw color_matrix directly; v3's raster transpose mirrored X.",
        side = 1, line = 5.6, cex = 0.65, col = "#7a1f1f")
  mtext(sprintf("Y top-to-bottom: FIRST %s (%s) -> LAST %s (%s)",
                first_row, first_row_comm, last_row, last_row_comm),
        side = 2, line = 5.4, cex = 0.67, col = "#333333")

  par(mar = c(7.4, 1.0, 6.4, 4.0), pty = "m")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  y <- seq(0.22, 0.88, length.out = length(pal) + 1)
  for (i in seq_along(pal)) {
    rect(0.22, y[i], 0.38, y[i + 1], col = pal[i], border = NA)
  }
  rect(0.22, 0.22, 0.38, 0.88, border = "#555555", lwd = 0.8)
  text(0.44, 0.88, sprintf("%.0e", 10 ^ zlim[2]),
       adj = c(0, 0.5), cex = 0.70)
  text(0.44, 0.22, sprintf("%.0e", 10 ^ zlim[1]),
       adj = c(0, 0.5), cex = 0.70)
  text(0.78, 0.55, "contact value (log10)", srt = 90, cex = 0.72)
}

png_path <- file.path(out_dir, "candidate_10a_xaxis_orientation.png")
pdf_path <- file.path(out_dir, "candidate_10a_xaxis_orientation.pdf")
render_plot(png_path, "png")
render_plot(pdf_path, "pdf")

png_dims <- detect_png_dims(png_path)
pdf_dims <- "10x10 inches"

readme <- c(
  "# Review zoom v4 slide 10a X-axis orientation fix",
  "",
  "Task: `review-zoom-v4-slide10a-xaxis-orientation-fix`.",
  "",
  "## Outputs",
  "",
  paste0("- `candidate_10a_xaxis_orientation.png`: corrected square PNG candidate (", png_dims, ")."),
  paste0("- `candidate_10a_xaxis_orientation.pdf`: corrected square PDF candidate (", pdf_dims, ")."),
  "- `make_10a_xaxis_orientation.R`: generator with source/statistic assertions and display-axis audit.",
  "- `orientation_audit.tsv`: per-index row order, corrected displayed X order, v3 implicit displayed X order, and corrected X-axis policy.",
  "- `ordered_arm_haplotypes.tsv`: final displayed row/column order plus the v3 implicit X order.",
  "- `sequence_community_boxes.tsv`: exact displayed-coordinate community box coordinates.",
  "- `matrix_order_audit.tsv`: source, statistic, and orientation checks.",
  "",
  "## Source Files Inspected",
  "",
  paste0("- v3 slide 10a generator: `", v3_generator, "`."),
  paste0("- v3 slide 10a asset: `", v3_asset, "` (", detect_png_dims(v3_asset), ")."),
  paste0("- v3 audit: `", v3_audit, "`."),
  paste0("- Original manuscript Fig. 3 generator: `", fig3_generator, ":32-80`; panel A uses base `image(..., t(vals_norm)[, n:1])`."),
  paste0("- v2 redesign generator: `", v2_generator, ":28-93`; it also used the base-image transform shape with `rasterImage()`."),
  paste0("- Contact matrix: `", mat_path, "`."),
  paste0("- Sequence community table for boxes: `", seq_comm_path, "`."),
  paste0("- Source B/W and p-value TSV: `", global_path, "`."),
  "",
  "## Orientation Finding",
  "",
  "- The v3 validation asserted that analytical row and column names were identical after ordering. That is true, but it does not validate the displayed X axis.",
  "- The v3 renderer used `rasterImage(as.raster(t(color_matrix)[, n:1]))`. That transform is appropriate for base `image()` semantics in the original Fig. 3 panel, but `rasterImage()` already displays raster rows top-to-bottom and columns left-to-right.",
  "- Therefore v3's visual X axis was mirrored left/right: its implicit displayed column order was `reverse(ordered_arms)`. The symmetric contact matrix made the mistake easy to miss by visual block structure and by row/column-name equality checks.",
  paste0("- Corrected X-axis policy: ", corrected_x_axis_policy),
  paste0("- Corrected row-axis policy: ", corrected_row_axis_policy),
  "",
  "## Statistic Check",
  "",
  paste0("- Source within mean = ", format(wb$within_mean, digits = 12), " over ", wb$n_within, " within-community pairs."),
  paste0("- Source between mean = ", format(wb$between_mean, digits = 12), " over ", wb$n_between, " between-community pairs."),
  paste0("- B/W = between / within = ", format(source_bw, digits = 8), ", displayed as ", sprintf("%.3f", source_bw), "."),
  paste0("- p-value = ", format_sci(wb$p_value, digits = 3), ", displayed as ", sprintf("%.1e", wb$p_value), "."),
  "",
  "The generator recomputes the within and between means from the HG002 Pore-C contact matrix plus the expanded sequence-community table, then asserts that the recomputed values match `hg002_porec_global_test.tsv`. The B/W statistic, p-value, and sequence-community interpretation are preserved because those checks pass.",
  "",
  "## Validation",
  "",
  paste0("- v3 X-axis mirrored stated in `orientation_audit.tsv`: TRUE."),
  paste0("- Corrected displayed X axis starts at `", corrected_display_col_order[1], "` and ends at `", corrected_display_col_order[n], "`."),
  paste0("- Corrected PNG square: ", ifelse(grepl("^([0-9]+)x\\1$", png_dims), "TRUE", png_dims), "."),
  "- Community boxes are computed in displayed coordinates after the corrected raster policy.",
  "- The corrected asset is generated from source matrix values, not copied from the v3 PNG."
)

writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
message("Wrote: ", file.path(out_dir, "orientation_audit.tsv"))
message("Wrote: ", file.path(out_dir, "README.md"))
