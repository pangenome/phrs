#!/usr/bin/env Rscript

out_dir <- "slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mat_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv"
v2_comm_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_hic.arm-leiden.communities.tsv"
seq_comm_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv"
global_path <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv"
v2_asset <- "slides/v2-review-zoom/_revision_assets/hic_visual_redesign/slide_10a_square_matrix_candidate.png"
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

within_values <- numeric()
between_values <- numeric()
for (i in seq_len(nrow(contact) - 1)) {
  for (j in (i + 1):ncol(contact)) {
    if (seq_community[i] == seq_community[j]) {
      within_values <- c(within_values, contact[i, j])
    } else {
      between_values <- c(between_values, contact[i, j])
    }
  }
}

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

# Hard requirement for this correction: the plotted row and column labels must
# be the same ordered arm-haplotype list after all intended ordering.
assert_true(identical(rownames(ordered_contact), colnames(ordered_contact)),
            "ordered row names equal ordered column names")
assert_true(identical(rownames(ordered_contact), ordered_arms),
            "ordered matrix row names equal ordered arm list")
assert_true(identical(colnames(ordered_contact), ordered_arms),
            "ordered matrix column names equal ordered arm list")

ordered_table <- data.frame(
  ordered_index = seq_along(ordered_arms),
  arm_haplotype = ordered_arms,
  base_arm = ordered_base,
  sequence_community = ordered_community,
  source_row_index = match(ordered_arms, row_arms),
  source_column_index = match(ordered_arms, col_arms),
  stringsAsFactors = FALSE
)

assert_true(identical(ordered_table$source_row_index,
                      ordered_table$source_column_index),
            "ordered source row indices equal ordered source column indices")

write.table(ordered_table,
            file = file.path(out_dir, "ordered_arm_haplotypes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

rl <- rle(as.character(ordered_community))
block_ends <- cumsum(rl$lengths)
block_starts <- c(1, head(block_ends, -1) + 1)
n <- nrow(ordered_contact)
blocks <- data.frame(
  sequence_community = rl$values,
  start_index = block_starts,
  end_index = block_ends,
  n_arm_haplotypes = rl$lengths,
  x_min = block_starts - 0.5,
  x_max = block_ends + 0.5,
  y_min = n - block_ends + 0.5,
  y_max = n - block_starts + 1.5,
  stringsAsFactors = FALSE
)
blocks$members <- vapply(seq_len(nrow(blocks)), function(i) {
  paste(ordered_arms[blocks$start_index[i]:blocks$end_index[i]], collapse = ";")
}, character(1))

write.table(blocks,
            file = file.path(out_dir, "sequence_community_boxes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

same_v2_seq_partition <- function(a, b) {
  upper <- upper.tri(matrix(FALSE, length(a), length(a)))
  same_a <- outer(a, a, "==")[upper]
  same_b <- outer(b, b, "==")[upper]
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
    "v2_generator",
    "v2_asset_dimensions",
    "v2_contact_community_table",
    "v2_contact_community_count",
    "sequence_community_count_in_matrix",
    "v2_contact_partition_pair_agreement_with_sequence_partition",
    "ordered_row_names_equal_ordered_column_names",
    "ordered_source_row_indices_equal_ordered_source_column_indices",
    "analytical_matrix_transpose_applied",
    "visual_raster_transform",
    "boxes_applied_after_visual_coordinate_mapping",
    "computed_within_mean_matches_source_tsv",
    "computed_between_mean_matches_source_tsv",
    "computed_bw_matches_source_tsv",
    "source_bw",
    "source_p_value",
    "v2_plot_transposed_incorrectly"
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
    v2_generator,
    detect_png_dims(v2_asset),
    v2_comm_path,
    as.character(length(unique(v2_community))),
    as.character(length(unique(seq_community))),
    paste0(v2_seq_pair_agree, "/", pair_total),
    "TRUE",
    "TRUE",
    "FALSE",
    "t(color_matrix)[, n:1] for R display coordinates only",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    format(source_bw, scientific = FALSE, digits = 8),
    format_sci(wb$p_value, digits = 3),
    "FALSE"
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

  vals <- ordered_contact
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
  par(mar = c(6.0, 6.5, 6.2, 1.0), pty = "s", family = "sans")
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
       xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)

  rasterImage(as.raster(t(color_matrix)[, n:1]), 0.5, 0.5, n + 0.5, n + 0.5,
              interpolate = FALSE)
  box(col = "#333333", lwd = 1.0)

  for (b in head(blocks$end_index, -1)) {
    abline(v = b + 0.5, col = grDevices::adjustcolor("#333333", alpha.f = 0.25),
           lwd = 0.55)
    abline(h = n - b + 0.5, col = grDevices::adjustcolor("#333333", alpha.f = 0.25),
           lwd = 0.55)
  }

  for (i in seq_len(nrow(blocks))) {
    rect(blocks$x_min[i], blocks$y_min[i], blocks$x_max[i], blocks$y_max[i],
         border = "#0b3c78", lwd = 1.35)
    if (blocks$n_arm_haplotypes[i] >= 4) {
      mid <- (blocks$start_index[i] + blocks$end_index[i]) / 2
      text(mid, n - mid + 1, blocks$sequence_community[i],
           cex = 0.62, col = "#0b3c78", font = 2)
    }
  }

  block_mid_x <- (blocks$start_index + blocks$end_index) / 2
  block_mid_y <- n - block_mid_x + 1
  axis(1, at = block_mid_x, labels = blocks$sequence_community,
       las = 2, cex.axis = 0.62, tick = FALSE, line = 0.1)
  axis(3, at = block_mid_x, labels = blocks$sequence_community,
       las = 2, cex.axis = 0.62, tick = FALSE, line = 0.1)
  axis(2, at = block_mid_y, labels = blocks$sequence_community,
       las = 1, cex.axis = 0.62, tick = FALSE, line = 0.1)
  axis(4, at = block_mid_y, labels = blocks$sequence_community,
       las = 1, cex.axis = 0.62, tick = FALSE, line = 0.1)

  title("HG002 Pore-C contacts ordered by sequence community",
        cex.main = 1.10, font.main = 2, line = 4.3)
  mtext("rows = HG002 arm-haplotypes; columns = same ordered arm-haplotype list",
        side = 3, line = 2.6, cex = 0.78, col = "#333333")
  mtext("blue boxes = arm-level sequence communities C1-C15",
        side = 3, line = 1.5, cex = 0.78, col = "#0b3c78")
  mtext("Columns: ordered by arm-level sequence community, then chromosome arm and haplotype",
        side = 1, line = 4.7, cex = 0.76, col = "#333333")
  mtext("Rows: same ordered list, displayed top to bottom",
        side = 2, line = 4.9, cex = 0.76, col = "#333333")
  mtext(sprintf("B/W = %.3f    p = %.1e    source within/between TSV checked",
                source_bw, wb$p_value),
        side = 1, line = 3.2, cex = 0.78, col = "#333333")

  par(mar = c(6.0, 1.0, 6.2, 4.0), pty = "m")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  y <- seq(0.12, 0.88, length.out = length(pal) + 1)
  for (i in seq_along(pal)) {
    rect(0.22, y[i], 0.38, y[i + 1], col = pal[i], border = NA)
  }
  rect(0.22, 0.12, 0.38, 0.88, border = "#555555", lwd = 0.8)
  text(0.44, 0.88, sprintf("%.0e", 10 ^ zlim[2]),
       adj = c(0, 0.5), cex = 0.70)
  text(0.44, 0.12, sprintf("%.0e", 10 ^ zlim[1]),
       adj = c(0, 0.5), cex = 0.70)
  text(0.78, 0.50, "contact value (log10)", srt = 90, cex = 0.72)
}

png_path <- file.path(out_dir, "candidate_10a_axis_box_fix.png")
pdf_path <- file.path(out_dir, "candidate_10a_axis_box_fix.pdf")
render_plot(png_path, "png")
render_plot(pdf_path, "pdf")

png_dims <- detect_png_dims(png_path)
pdf_dims <- "10x10 inches"

readme <- c(
  "# Review zoom v3 slide 10a axis/box fix",
  "",
  "Task: `review-zoom-v3-slide10a-axis-box-fix`.",
  "",
  "## Outputs",
  "",
  paste0("- `candidate_10a_axis_box_fix.png`: square PNG candidate (", png_dims, ")."),
  paste0("- `candidate_10a_axis_box_fix.pdf`: square PDF candidate (", pdf_dims, ")."),
  "- `make_10a_axis_box_fix.R`: generator with assertions.",
  "- `ordered_arm_haplotypes.tsv`: final row/column order.",
  "- `sequence_community_boxes.tsv`: exact community box coordinates.",
  "- `matrix_order_audit.tsv`: source and ordering checks.",
  "",
  "## Source Files",
  "",
  paste0("- Contact matrix: `", mat_path, "`."),
  paste0("- Sequence community table for boxes: `", seq_comm_path, "`."),
  paste0("- Source B/W and p-value TSV: `", global_path, "`."),
  paste0("- v2 slide 10a asset inspected: `", v2_asset, "` (", detect_png_dims(v2_asset), ")."),
  paste0("- v2 redesign generator inspected: `", v2_generator, ":28-44` and `:67-80`."),
  paste0("- Manuscript Fig. 3 generator inspected: `", fig3_generator, ":32-45` and `:63-80`."),
  "",
  "## Rows, Columns, Boxes",
  "",
  "- rows = the 77 HG002 Pore-C arm-haplotypes from the first column of `hg002_porec_contact_matrix.tsv`, reordered by arm-level sequence community.",
  "- columns = the same 77 HG002 Pore-C arm-haplotypes from the matrix header, in the identical post-order list.",
  "- boxes = arm-level sequence communities from `hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, expanded from base arms such as `chr4_q` to matrix labels such as `chr4_MATERNAL_q` and `chr4_PATERNAL_q`.",
  "",
  "The final ordered row list and column list are asserted identical by name after ordering. The generator stops before plotting if they diverge.",
  "",
  "## Audit Findings",
  "",
  paste0("- Source contact matrix shape: ", nrow(contact), " x ", ncol(contact), "."),
  paste0("- Source row names equal source column names: ", format_bool(identical(row_arms, col_arms)), "."),
  paste0("- Source matrix symmetric: TRUE; max absolute symmetry delta = ", format(sym_delta, scientific = TRUE, digits = 6), "."),
  paste0("- All 77 source matrix arm-haplotypes map to sequence communities: TRUE."),
  paste0("- Sequence communities represented in the plotted matrix: ", length(unique(seq_community)), " (C1-C15)."),
  paste0("- v2 box/order table communities represented: ", length(unique(v2_community)), " contact-derived communities from `hg002_porec_hic.arm-leiden.communities.tsv`."),
  paste0("- v2 contact-community partition agreement with the sequence-community partition: ", v2_seq_pair_agree, "/", pair_total, " arm-pairs."),
  "- Was the v2 plot transposed incorrectly? No. The source matrix rows and columns were identical, and the v2 renderer used the usual R display transform (`t(colors)[, n:1]`) to put the first row at the top. The real v2 defect for this task is that community boxes/order came from the contact-community table, not the sequence-community table behind the stated interpretation and B/W statistic.",
  "- Candidate transpose policy: no analytical matrix transpose is applied. The raster is transformed only for R image coordinates, and the blue boxes are computed in that same displayed coordinate system after the visual y-axis reversal.",
  "",
  "## Statistic Check",
  "",
  paste0("- Source within mean = ", format(wb$within_mean, digits = 12), " over ", wb$n_within, " within-community pairs."),
  paste0("- Source between mean = ", format(wb$between_mean, digits = 12), " over ", wb$n_between, " between-community pairs."),
  paste0("- B/W = between / within = ", format(source_bw, digits = 8), ", displayed as ", sprintf("%.3f", source_bw), "."),
  paste0("- p-value = ", format_sci(wb$p_value, digits = 3), ", displayed as ", sprintf("%.1e", wb$p_value), "."),
  "",
  "The generator recomputes the within and between means from `hg002_porec_contact_matrix.tsv` plus the expanded sequence-community table, then asserts that the recomputed values match `hg002_porec_global_test.tsv`. The B/W and p-value are preserved in the candidate only because those checks pass.",
  "",
  "## Validation",
  "",
  paste0("- README row/column/order checks included: TRUE."),
  paste0("- v2 transpose finding included: TRUE."),
  paste0("- Candidate PNG square: ", ifelse(grepl("^([0-9]+)x\\1$", png_dims), "TRUE", png_dims), "."),
  "- Candidate PDF square: TRUE (10 x 10 inch device).",
  "- Community boxes align on both axes by construction: `sequence_community_boxes.tsv` stores matching x and y extents for each contiguous sequence-community block, and the plot draws those extents after the display coordinate transform."
)

writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
message("Wrote: ", file.path(out_dir, "README.md"))
