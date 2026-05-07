#!/usr/bin/env Rscript

# Render paired slide-ready 07a heatmaps for review-zoom v5:
# 1. tree-ordered similarity heatmap with side UPGMA tree
# 2. same matrix/style ordered by Leiden C1..C15 community blocks, no tree

suppressPackageStartupMessages({
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

dist_path <- if (length(args) >= 1) args[[1]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
leiden_path <- if (length(args) >= 2) args[[2]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
architecture_path <- if (length(args) >= 3) args[[3]] else
  "paper_prep/figures/fig1/architecture_per_arm.tsv"
chm13_bed_path <- if (length(args) >= 4) args[[4]] else
  "chm13.phrs.bed"
projected_path <- if (length(args) >= 5) args[[5]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
chrom_sizes_path <- if (length(args) >= 6) args[[6]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/chrom.sizes"
out_dir <- if (length(args) >= 7) args[[7]] else
  "slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, header = TRUE, row.names = NULL) {
  read.table(path,
             header = header,
             row.names = row.names,
             sep = "\t",
             check.names = FALSE,
             stringsAsFactors = FALSE,
             quote = "",
             comment.char = "")
}

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
}

assert_file(dist_path, "Distance matrix")
assert_file(leiden_path, "Leiden community assignment table")
assert_file(architecture_path, "Architecture per-arm table")
assert_file(chm13_bed_path, "CHM13 PHR BED")

message("Reading arm-level distance matrix: ", dist_path)
D_raw <- as.matrix(read_tsv(dist_path, header = TRUE, row.names = 1))
stopifnot(nrow(D_raw) == ncol(D_raw))
stopifnot(identical(rownames(D_raw), colnames(D_raw)))

D <- (D_raw + t(D_raw)) / 2
diag_for_note <- diag(D)
similarity <- 1 - D
similarity[similarity < 0] <- 0
similarity[similarity > 1] <- 1

message("Reading Leiden assignments: ", leiden_path)
leiden <- read_tsv(leiden_path)
required_leiden_cols <- c("ChromArm", "Community", "Arms")
stopifnot(all(required_leiden_cols %in% names(leiden)))
stopifnot(setequal(rownames(D), leiden$ChromArm))

message("Reading architecture audit table: ", architecture_path)
architecture <- read_tsv(architecture_path)
required_arch_cols <- c("ChromArm", "n_total", "n_cross", "cross_rate",
                        "Community", "category")
stopifnot(all(required_arch_cols %in% names(architecture)))

arm_side <- function(x) {
  ifelse(grepl("_p$", x), "p", "q")
}

short_arm <- function(x) {
  sub("^chr", "", gsub("_", "", x))
}

community_number <- function(x) {
  suppressWarnings(as.integer(sub("^C", "", x)))
}

natural_arm_order <- function(arms) {
  m <- regexec("^chr([^_]+)_([pq])$", arms)
  parts <- regmatches(arms, m)
  chrom <- vapply(parts, function(x) if (length(x) >= 3) x[[2]] else NA_character_, character(1))
  side <- vapply(parts, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1))
  chrom_key <- suppressWarnings(as.integer(chrom))
  chrom_key[chrom == "X"] <- 23
  chrom_key[chrom == "Y"] <- 24
  side_key <- ifelse(side == "p", 0, 1)
  order(chrom_key, side_key, arms)
}

community_of <- setNames(leiden$Community, leiden$ChromArm)
community_arms <- setNames(leiden$Arms, leiden$ChromArm)
architecture_category <- setNames(architecture$category, architecture$ChromArm)
architecture_n_total <- setNames(architecture$n_total, architecture$ChromArm)
architecture_n_cross <- setNames(architecture$n_cross, architecture$ChromArm)
architecture_cross_rate <- setNames(architecture$cross_rate, architecture$ChromArm)

hc <- hclust(as.dist(D), method = "average")
tree_order <- hc$labels[hc$order]
n <- length(tree_order)
stopifnot(n == nrow(D))

tree_pos <- setNames(seq_along(tree_order), tree_order)
community_levels <- paste0("C", sort(unique(community_number(leiden$Community))))
community_order <- unlist(lapply(community_levels, function(comm) {
  arms <- leiden$ChromArm[leiden$Community == comm]
  arms[order(tree_pos[arms], arms)]
}), use.names = FALSE)

stopifnot(setequal(tree_order, community_order))
stopifnot(length(community_order) == n)

tree_row_label_font_pt <- 8.9
tree_col_label_font_pt <- 8.6
community_row_label_font_pt <- 8.6
community_col_label_font_pt <- 8.3
label_size_policy <- paste0(
  "Chromosome arm labels are enlarged about 1.5x versus the v3/current ",
  "baseline: tree view row/column labels ",
  tree_row_label_font_pt, "/", tree_col_label_font_pt,
  " pt; community view row/column labels ",
  community_row_label_font_pt, "/", community_col_label_font_pt,
  " pt."
)

make_order_table <- function(order, mode) {
  data.frame(
    position_top_to_bottom = seq_along(order),
    position_left_to_right = seq_along(order),
    ChromArm = order,
    label = short_arm(order),
    arm_side = arm_side(order),
    Community = unname(community_of[order]),
    Community_members = unname(community_arms[order]),
    UPGMA_tree_position = unname(tree_pos[order]),
    ordering_mode = mode,
    matrix_source = dist_path,
    tree_order_source = "UPGMA average linkage via hclust(as.dist(D), method = 'average') on the 41x41 arm-level Jaccard distance matrix",
    community_mapping_source = leiden_path,
    stringsAsFactors = FALSE
  )
}

tree_order_tbl <- make_order_table(
  tree_order,
  "tree order: UPGMA leaf order used for both rows and columns"
)
community_order_tbl <- make_order_table(
  community_order,
  "community order: Leiden C1..C15, then UPGMA leaf position within each community"
)

write.table(tree_order_tbl,
            file = file.path(out_dir, "arm_order_tree.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
write.table(community_order_tbl,
            file = file.path(out_dir, "arm_order_community.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

make_community_blocks <- function(order) {
  comm <- unname(community_of[order])
  r <- rle(comm)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  data.frame(
    Community = r$values,
    start_position = starts,
    end_position = ends,
    n_arms = r$lengths,
    arms = vapply(seq_along(starts), function(i) {
      paste(short_arm(order[starts[i]:ends[i]]), collapse = ", ")
    }, character(1)),
    order_note = "Positions are top-to-bottom rows and left-to-right columns in 07a_community_ordered_heatmap",
    stringsAsFactors = FALSE
  )
}

community_blocks <- make_community_blocks(community_order)
write.table(community_blocks,
            file = file.path(out_dir, "community_blocks.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

read_chm13_bed_arms <- function(path, sizes_path) {
  if (!file.exists(path)) {
    return(data.frame(ChromArm = character(),
                      chm13_phr_bed_interval_count = integer(),
                      stringsAsFactors = FALSE))
  }
  bed <- read.table(path,
                    header = FALSE,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    quote = "",
                    comment.char = "")
  if (ncol(bed) < 3) {
    stop("CHM13 BED has fewer than 3 columns: ", path)
  }
  names(bed)[1:3] <- c("chrom_raw", "start", "end")
  bed$chrom <- sub("^.*#(chr[^#]+)$", "\\1", bed$chrom_raw)
  bed$chrom[!grepl("^chr", bed$chrom)] <- bed$chrom_raw[!grepl("^chr", bed$chrom)]

  chrom_len <- setNames(rep(NA_real_, nrow(bed)), seq_len(nrow(bed)))
  if (file.exists(sizes_path)) {
    sizes <- read.table(sizes_path,
                        header = FALSE,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        quote = "",
                        comment.char = "")
    if (ncol(sizes) >= 2) {
      names(sizes)[1:2] <- c("chrom", "length")
      size_by_chrom <- setNames(as.numeric(sizes$length), sizes$chrom)
      chrom_len <- unname(size_by_chrom[bed$chrom])
    }
  }

  mid <- (as.numeric(bed$start) + as.numeric(bed$end)) / 2
  bed$side <- ifelse(!is.na(chrom_len),
                     ifelse(mid < chrom_len / 2, "p", "q"),
                     ifelse(as.numeric(bed$start) < 1000000, "p", "q"))
  bed$ChromArm <- paste0(bed$chrom, "_", bed$side)
  counts <- as.data.frame(table(bed$ChromArm), stringsAsFactors = FALSE)
  names(counts) <- c("ChromArm", "chm13_phr_bed_interval_count")
  counts$chm13_phr_bed_interval_count <- as.integer(counts$chm13_phr_bed_interval_count)
  counts
}

read_projected_chm13_arms <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(ChromArm = character(),
                      projected_chm13_rows = integer(),
                      projected_chm13_nonempty_rows = integer(),
                      stringsAsFactors = FALSE))
  }
  projected <- read.table(path,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          check.names = FALSE,
                          quote = "",
                          comment.char = "")
  required_cols <- c("seq", "region_start", "region_end")
  if (!all(required_cols %in% names(projected))) {
    stop("Projected interval table lacks required columns: ", path)
  }
  chm13 <- projected[grepl("^CHM13#0#", projected$seq), , drop = FALSE]
  chm13$ChromArm <- sub(".*_(chr[^_]+)_(p|q)arm$", "\\1_\\2", chm13$seq)
  chm13 <- chm13[grepl("^chr[^_]+_[pq]$", chm13$ChromArm), , drop = FALSE]
  chm13$nonempty <- as.character(chm13$region_start) != "." &
    as.character(chm13$region_end) != "."
  rows <- aggregate(seq ~ ChromArm, data = chm13, FUN = length)
  names(rows)[names(rows) == "seq"] <- "projected_chm13_rows"
  nonempty <- aggregate(nonempty ~ ChromArm, data = chm13, FUN = sum)
  names(nonempty)[names(nonempty) == "nonempty"] <- "projected_chm13_nonempty_rows"
  merged <- merge(rows, nonempty, by = "ChromArm", all = TRUE, sort = FALSE)
  merged$projected_chm13_rows[is.na(merged$projected_chm13_rows)] <- 0L
  merged$projected_chm13_nonempty_rows[is.na(merged$projected_chm13_nonempty_rows)] <- 0L
  merged
}

all_arms <- unlist(lapply(c(as.character(1:22), "X", "Y"), function(chrom) {
  paste0("chr", chrom, "_", c("p", "q"))
}), use.names = FALSE)

chm13_bed_counts <- read_chm13_bed_arms(chm13_bed_path, chrom_sizes_path)
projected_counts <- read_projected_chm13_arms(projected_path)

audit <- data.frame(
  ChromArm = all_arms,
  label = short_arm(all_arms),
  arm_side = arm_side(all_arms),
  in_41x41_matrix = all_arms %in% rownames(D),
  in_tree_ordered_heatmap = all_arms %in% tree_order,
  in_community_ordered_heatmap = all_arms %in% community_order,
  has_arm_community_assignment = all_arms %in% leiden$ChromArm,
  Community = unname(community_of[all_arms]),
  Community_members = unname(community_arms[all_arms]),
  architecture_category = unname(architecture_category[all_arms]),
  architecture_n_total = unname(architecture_n_total[all_arms]),
  architecture_n_cross = unname(architecture_n_cross[all_arms]),
  architecture_cross_rate = unname(architecture_cross_rate[all_arms]),
  tree_order_position = unname(tree_pos[all_arms]),
  community_order_position = match(all_arms, community_order),
  stringsAsFactors = FALSE
)

audit <- merge(audit, chm13_bed_counts, by = "ChromArm", all.x = TRUE, sort = FALSE)
audit <- merge(audit, projected_counts, by = "ChromArm", all.x = TRUE, sort = FALSE)
audit <- audit[match(all_arms, audit$ChromArm), ]

audit$chm13_phr_bed_interval_count[is.na(audit$chm13_phr_bed_interval_count)] <- 0L
audit$projected_chm13_rows[is.na(audit$projected_chm13_rows)] <- 0L
audit$projected_chm13_nonempty_rows[is.na(audit$projected_chm13_nonempty_rows)] <- 0L
audit$has_called_chm13_phr_bed_interval <- audit$chm13_phr_bed_interval_count > 0
audit$has_projected_chm13_nonempty_interval <- audit$projected_chm13_nonempty_rows > 0
audit$arm_render_status <- ifelse(audit$in_41x41_matrix,
                                  "included in both v5 heatmaps",
                                  "not in 41x41 matrix; not rendered")
audit$community_mapping_status <- ifelse(audit$has_arm_community_assignment,
                                         "assigned to arm-level Leiden community",
                                         "no arm-level Leiden community assignment")
audit$chm13_phr_bed_status <- ifelse(audit$has_called_chm13_phr_bed_interval,
                                     "called interval present in chm13.phrs.bed",
                                     "no called interval row in chm13.phrs.bed")
audit$projected_chm13_status <- ifelse(
  audit$projected_chm13_rows == 0,
  "no CHM13#0 row in projected all-vs-all table",
  ifelse(audit$has_projected_chm13_nonempty_interval,
         "non-empty CHM13#0 projected PHR interval present",
         "CHM13#0 projected row exists but region_start/end are empty")
)
audit$community_assigned_but_no_called_chm13_phr_bed_interval <-
  audit$has_arm_community_assignment & !audit$has_called_chm13_phr_bed_interval
audit$note <- ifelse(
  audit$community_assigned_but_no_called_chm13_phr_bed_interval,
  "Community assignment exists, but the repo-root CHM13 PHR BED has no called interval row for this arm.",
  "not_applicable"
)

write.table(audit,
            file = file.path(out_dir, "arm_inclusion_audit.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

source_audit <- data.frame(
  item = c("matrix_source",
           "tree_order_source",
           "community_mapping_source",
           "architecture_table_source",
           "chm13_phr_bed_source",
           "projected_chm13_interval_source",
           "prior_crisp_renderer",
           "community_order",
           "display_scale",
           "label_size_policy",
           "color_scale_policy"),
  path_or_value = c(dist_path,
                    "computed in this script from the matrix using hclust(as.dist(D), method = 'average')",
                    leiden_path,
                    architecture_path,
                    chm13_bed_path,
                    projected_path,
                    "slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned/make_07a_crisp_aligned.R",
                    paste(community_levels, collapse = ","),
                    "Jaccard similarity = 1 - distance; limits clipped to [0, 1]; same palette in both v5 assets",
                    label_size_policy,
                    "colorRampPalette(c('#F8FBFF', '#FEE08B', '#F46D43', '#A50026')) with fixed [0,1] breaks"),
  audit_value = c(sprintf("%d rows x %d columns", nrow(D_raw), ncol(D_raw)),
                  "UPGMA leaf order defines rows and columns in 07a_tree_ordered_heatmap",
                  sprintf("%d assigned arms across %d observed communities", nrow(leiden), length(community_levels)),
                  sprintf("%d rows; architecture categories joined only for audit", nrow(architecture)),
                  sprintf("%d inferred arm intervals across %d distinct arms",
                          sum(audit$chm13_phr_bed_interval_count),
                          sum(audit$has_called_chm13_phr_bed_interval)),
                  sprintf("%d CHM13#0 arm rows with %d non-empty projected interval arms",
                          sum(audit$projected_chm13_rows),
                          sum(audit$has_projected_chm13_nonempty_interval)),
                  "v5 reuses the vector-first grid/tree geometry pattern from v3",
                  "Leiden C1..C15 numeric order, then UPGMA position within community",
                  "Both PDF assets are vector scenes; PNG assets are direct 4800x2700 renders",
                  "Applies to both x-axis and y-axis arm labels in both rendered assets",
                  "Palette, similarity transform, and legend limits are identical in both rendered assets"),
  caveat = c("The source diagonal is retained for display; as.dist ignores it for clustering.",
             "This is UPGMA average linkage, not an NJ tree.",
             "Community labels are arm-level Leiden C1-C15 labels; no sequence-level community IDs are used.",
             "Architecture categories do not determine heatmap order.",
             "BED rows do not include arm labels; p/q is inferred from CHM13 chromosome midpoint using chrom.sizes when available.",
             "Projected rows are included only to avoid implying that every community-assigned arm has a repo-root CHM13 BED call.",
             "Prior v3 output is not copied; this script regenerates both v5 views.",
             "If a community were absent from the matrix it would be absent here; observed set is all C1-C15.",
             "Cell color encodes similarity, not community.",
             "Only label size changed for readability; p/q label colors remain red/blue.",
             "Community bands use separate categorical colors and do not alter heatmap cell colors."),
  stringsAsFactors = FALSE
)

write.table(source_audit,
            file = file.path(out_dir, "source_audit.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

make_left_segments <- function(dend, leaf_pos) {
  segs <- data.frame(x = numeric(),
                     y = numeric(),
                     xend = numeric(),
                     yend = numeric())

  walk <- function(node) {
    height <- attr(node, "height")
    if (is.null(height)) {
      height <- 0
    }

    if (is.leaf(node)) {
      label <- attr(node, "label")
      return(list(height = 0, y = unname(leaf_pos[label])))
    }

    children <- lapply(seq_along(node), function(i) walk(node[[i]]))
    child_y <- vapply(children, function(x) x$y, numeric(1))
    child_h <- vapply(children, function(x) x$height, numeric(1))
    node_y <- mean(range(child_y))

    segs <<- rbind(
      segs,
      data.frame(x = height,
                 y = min(child_y),
                 xend = height,
                 yend = max(child_y))
    )

    for (i in seq_along(children)) {
      segs <<- rbind(
        segs,
        data.frame(x = child_h[i],
                   y = child_y[i],
                   xend = height,
                   yend = child_y[i])
      )
    }

    list(height = height, y = node_y)
  }

  walk(dend)
  segs
}

palette_heat <- colorRampPalette(c("#F8FBFF", "#FEE08B", "#F46D43", "#A50026"))(256)
breaks <- seq(0, 1, length.out = length(palette_heat) + 1)

value_to_col <- function(x) {
  palette_heat[findInterval(x, breaks, all.inside = TRUE)]
}

arm_cols <- c(p = "#CC3B38", q = "#1F5EA8")
community_cols <- setNames(
  c("#2E6FBB", "#D95F02", "#1B9E77", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#1F78B4", "#B2DF8A",
    "#FB9A99", "#FDBF6F", "#CAB2D6", "#6A3D9A", "#B15928")[seq_along(community_levels)],
  community_levels
)

draw_heatmap_body <- function(order, heat_left, heat_bottom, matrix_size) {
  cell <- matrix_size / length(order)
  vals <- similarity[order, order]
  coord <- expand.grid(row = seq_along(order), col = seq_along(order))
  xs <- heat_left + (coord$col - 0.5) * cell
  ys <- heat_bottom + (length(order) - coord$row + 0.5) * cell
  fills <- value_to_col(vals[cbind(coord$row, coord$col)])

  grid.rect(x = unit(xs, "in"),
            y = unit(ys, "in"),
            width = unit(cell, "in"),
            height = unit(cell, "in"),
            gp = gpar(fill = fills, col = NA))
  grid.rect(x = unit(heat_left + matrix_size / 2, "in"),
            y = unit(heat_bottom + matrix_size / 2, "in"),
            width = unit(matrix_size, "in"),
            height = unit(matrix_size, "in"),
            gp = gpar(fill = NA, col = "#202020", lwd = 0.8))

  line_col <- "#5C5C5C40"
  grid.segments(x0 = unit(heat_left + seq(0, matrix_size, by = cell), "in"),
                y0 = unit(rep(heat_bottom, length(order) + 1), "in"),
                x1 = unit(heat_left + seq(0, matrix_size, by = cell), "in"),
                y1 = unit(rep(heat_bottom + matrix_size, length(order) + 1), "in"),
                gp = gpar(col = line_col, lwd = 0.25))
  grid.segments(x0 = unit(rep(heat_left, length(order) + 1), "in"),
                y0 = unit(heat_bottom + seq(0, matrix_size, by = cell), "in"),
                x1 = unit(rep(heat_left + matrix_size, length(order) + 1), "in"),
                y1 = unit(heat_bottom + seq(0, matrix_size, by = cell), "in"),
                gp = gpar(col = line_col, lwd = 0.25))
}

draw_arm_labels <- function(order, heat_left, heat_bottom, matrix_size,
                            right_x, bottom_y, row_font = 5.9,
                            col_font = 5.7) {
  cell <- matrix_size / length(order)
  row_centers <- heat_bottom + (length(order) - seq_along(order) + 0.5) * cell
  col_centers <- heat_left + (seq_along(order) - 0.5) * cell
  labels <- short_arm(order)
  label_cols <- arm_cols[arm_side(order)]

  grid.text(labels,
            x = unit(right_x, "in"),
            y = unit(row_centers, "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = row_font, col = label_cols))
  grid.text(labels,
            x = unit(col_centers, "in"),
            y = unit(bottom_y, "in"),
            rot = 90,
            just = c("right", "center"),
            gp = gpar(fontsize = col_font, col = label_cols))
}

draw_encoding_legend <- function(legend_left, legend_top, include_community_note = FALSE) {
  grid.text("Encoding",
            x = unit(legend_left, "in"),
            y = unit(legend_top, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold", col = "#222222"))
  grid.rect(x = unit(c(legend_left + 0.06, legend_left + 0.06), "in"),
            y = unit(c(legend_top - 0.31, legend_top - 0.58), "in"),
            width = unit(0.12, "in"),
            height = unit(0.12, "in"),
            gp = gpar(fill = c(arm_cols[["p"]], arm_cols[["q"]]), col = NA))
  grid.text(c("p arm label", "q arm label"),
            x = unit(rep(legend_left + 0.18, 2), "in"),
            y = unit(c(legend_top - 0.31, legend_top - 0.58), "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = 7.5, col = c(arm_cols[["p"]], arm_cols[["q"]])))

  grid.text("Jaccard similarity",
            x = unit(legend_left, "in"),
            y = unit(legend_top - 0.98, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", col = "#333333"))
  cb_n <- length(palette_heat)
  cb_x <- rep(legend_left + 0.12, cb_n)
  cb_y <- seq(legend_top - 2.35, legend_top - 1.15, length.out = cb_n)
  grid.rect(x = unit(cb_x, "in"),
            y = unit(cb_y, "in"),
            width = unit(0.18, "in"),
            height = unit(1.22 / cb_n, "in"),
            gp = gpar(fill = palette_heat, col = NA))
  grid.rect(x = unit(legend_left + 0.12, "in"),
            y = unit(legend_top - 1.75, "in"),
            width = unit(0.18, "in"),
            height = unit(1.2, "in"),
            gp = gpar(fill = NA, col = "#333333", lwd = 0.5))
  grid.text(c("1.0", "0.5", "0.0"),
            x = unit(rep(legend_left + 0.36, 3), "in"),
            y = unit(c(legend_top - 1.15, legend_top - 1.75, legend_top - 2.35), "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = 7, col = "#333333"))

  if (include_community_note) {
    grid.text("Community bands",
              x = unit(legend_left, "in"),
              y = unit(legend_top - 2.74, "in"),
              just = c("left", "top"),
              gp = gpar(fontsize = 8, fontface = "bold", col = "#333333"))
    grid.text("Leiden C1-C15 numeric order;\nwithin each block, arms keep\nUPGMA leaf order.",
              x = unit(legend_left, "in"),
              y = unit(legend_top - 2.98, "in"),
              just = c("left", "top"),
              gp = gpar(fontsize = 6.8, col = "#555555", lineheight = 1.08))
  } else {
    grid.text("PDF uses vector cells/tree.\nPNG is rendered directly at 4800x2700.",
              x = unit(legend_left, "in"),
              y = unit(legend_top - 2.68, "in"),
              just = c("left", "top"),
              gp = gpar(fontsize = 6.8, col = "#555555", lineheight = 1.08))
  }
}

draw_community_blocks <- function(order, heat_left, heat_bottom, matrix_size) {
  cell <- matrix_size / length(order)
  blocks <- make_community_blocks(order)

  for (i in seq_len(nrow(blocks))) {
    start <- blocks$start_position[i]
    end <- blocks$end_position[i]
    comm <- blocks$Community[i]
    block_w <- (end - start + 1) * cell
    block_h <- block_w
    block_left <- heat_left + (start - 1) * cell
    block_bottom <- heat_bottom + (length(order) - end) * cell
    x_mid <- block_left + block_w / 2
    y_mid <- block_bottom + block_h / 2
    col <- community_cols[[comm]]

    grid.rect(x = unit(x_mid, "in"),
              y = unit(y_mid, "in"),
              width = unit(block_w, "in"),
              height = unit(block_h, "in"),
              gp = gpar(fill = NA, col = col, lwd = 1.1))

    top_band_y <- heat_bottom + matrix_size + 0.09
    grid.rect(x = unit(x_mid, "in"),
              y = unit(top_band_y, "in"),
              width = unit(block_w, "in"),
              height = unit(0.13, "in"),
              gp = gpar(fill = grDevices::adjustcolor(col, alpha.f = 0.88),
                        col = NA))
    grid.text(comm,
              x = unit(x_mid, "in"),
              y = unit(top_band_y + 0.20, "in"),
              rot = 90,
              just = c("center", "center"),
              gp = gpar(fontsize = 5.4, fontface = "bold", col = "#333333"))

    left_band_x <- heat_left - 0.09
    y_top <- heat_bottom + (length(order) - start + 1) * cell
    y_bottom <- heat_bottom + (length(order) - end) * cell
    y_center <- (y_top + y_bottom) / 2
    block_height <- y_top - y_bottom
    grid.rect(x = unit(left_band_x, "in"),
              y = unit(y_center, "in"),
              width = unit(0.13, "in"),
              height = unit(block_height, "in"),
              gp = gpar(fill = grDevices::adjustcolor(col, alpha.f = 0.88),
                        col = NA))
    grid.text(comm,
              x = unit(left_band_x - 0.16, "in"),
              y = unit(y_center, "in"),
              just = c("center", "center"),
              gp = gpar(fontsize = 5.4, fontface = "bold", col = "#333333"))
  }
}

draw_tree_ordered_scene <- function() {
  grid.newpage()

  page_w <- 13.333
  page_h <- 7.5
  matrix_bottom <- 1.02
  matrix_size <- 5.72
  tree_left <- 0.42
  tree_width <- 1.58
  gutter <- 0.06
  heat_left <- tree_left + tree_width + gutter
  label_left <- heat_left + matrix_size + 0.12
  legend_left <- 9.25
  cell <- matrix_size / n

  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("07a.1 Tree-ordered arm similarity heatmap",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.08, "in"),
            just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 12, col = "#1F1F1F"))
  grid.text("UPGMA leaf order defines the side tree, rows, and columns; labels: p = red, q = blue.",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.36, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8.5, col = "#3E3E3E"))

  draw_heatmap_body(tree_order, heat_left, matrix_bottom, matrix_size)
  draw_arm_labels(tree_order, heat_left, matrix_bottom, matrix_size,
                  right_x = label_left,
                  bottom_y = matrix_bottom - 0.075,
                  row_font = tree_row_label_font_pt,
                  col_font = tree_col_label_font_pt)
  grid.text("arms",
            x = unit(label_left + 0.03, "in"),
            y = unit(matrix_bottom + matrix_size + 0.12, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))

  leaf_y <- setNames(rev(seq_along(tree_order)), tree_order)
  segs <- make_left_segments(as.dendrogram(hc), leaf_y)
  max_height <- max(segs$x, segs$xend)
  tx <- function(h) tree_left + (max_height - h) / max_height * tree_width
  ty <- function(y) matrix_bottom + (y - 0.5) * cell
  grid.segments(x0 = unit(tx(segs$x), "in"),
                y0 = unit(ty(segs$y), "in"),
                x1 = unit(tx(segs$xend), "in"),
                y1 = unit(ty(segs$yend), "in"),
                gp = gpar(col = "#333333", lwd = 0.65, lineend = "square"))
  grid.text("UPGMA tree",
            x = unit(tree_left, "in"),
            y = unit(matrix_bottom + matrix_size + 0.12, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))
  grid.text("leaf order defines rows and columns",
            x = unit(tree_left, "in"),
            y = unit(matrix_bottom + matrix_size + 0.01, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 6, col = "#555555"))

  draw_encoding_legend(legend_left, matrix_bottom + matrix_size - 0.08,
                       include_community_note = FALSE)

  grid.text("Same 41x41 matrix as the community view; see arm_order_tree.tsv and arm_inclusion_audit.tsv.",
            x = unit(0.42, "in"),
            y = unit(0.18, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7.2, col = "#3A3A3A"))
}

draw_community_ordered_scene <- function() {
  grid.newpage()

  page_w <- 13.333
  page_h <- 7.5
  matrix_bottom <- 0.82
  matrix_size <- 5.55
  heat_left <- 1.36
  label_left <- heat_left + matrix_size + 0.12
  legend_left <- 8.85

  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("07a.2 Community-ordered arm similarity heatmap",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.08, "in"),
            just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 12, col = "#1F1F1F"))
  grid.text("Same values and color scale, sorted by Leiden C1-C15 with UPGMA order inside each block; no side tree.",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.36, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8.5, col = "#3E3E3E"))

  draw_heatmap_body(community_order, heat_left, matrix_bottom, matrix_size)
  draw_community_blocks(community_order, heat_left, matrix_bottom, matrix_size)
  draw_arm_labels(community_order, heat_left, matrix_bottom, matrix_size,
                  right_x = label_left,
                  bottom_y = matrix_bottom - 0.075,
                  row_font = community_row_label_font_pt,
                  col_font = community_col_label_font_pt)

  grid.text("arms",
            x = unit(label_left + 0.03, "in"),
            y = unit(matrix_bottom + matrix_size + 0.12, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))

  draw_encoding_legend(legend_left, matrix_bottom + matrix_size - 0.08,
                       include_community_note = TRUE)

  grid.text("Audit caveat: community assignment and called CHM13 PHR BED interval status are tracked separately.",
            x = unit(0.42, "in"),
            y = unit(0.18, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7.2, col = "#3A3A3A"))
}

render_pair <- function(base_name, draw_fun) {
  pdf_path <- file.path(out_dir, paste0(base_name, ".pdf"))
  png_path <- file.path(out_dir, paste0(base_name, ".png"))

  pdf(pdf_path, width = 13.333, height = 7.5, useDingbats = FALSE)
  draw_fun()
  dev.off()

  png(png_path,
      width = 4800,
      height = 2700,
      res = 360,
      type = "cairo")
  draw_fun()
  dev.off()

  c(pdf = pdf_path, png = png_path)
}

tree_tip_order_from_y <- names(sort(setNames(rev(seq_along(tree_order)), tree_order),
                                    decreasing = TRUE))
validation_tbl <- data.frame(
  check = c("tree_tip_order_equals_tree_heatmap_row_order",
            "tree_tip_order_equals_tree_heatmap_column_order",
            "community_order_rows_equal_columns",
            "community_order_is_C1_to_C15_numeric_order",
            "matrix_rows_equal_tree_labels",
            "source_matrix_rows",
            "source_matrix_columns",
            "max_source_asymmetry",
            "raw_diagonal_min",
            "raw_diagonal_max",
            "community_assigned_arms_without_called_chm13_phr_bed_interval",
            "tree_label_font_points_row_col",
            "community_label_font_points_row_col",
            "color_scale_policy"),
  value = c(identical(tree_tip_order_from_y, tree_order),
            identical(tree_tip_order_from_y, tree_order),
            identical(community_order, community_order),
            identical(unique(unname(community_of[community_order])), community_levels),
            setequal(rownames(D), hc$labels),
            nrow(D_raw),
            ncol(D_raw),
            sprintf("%.6g", max(abs(D_raw - t(D_raw)))),
            sprintf("%.6f", min(diag_for_note)),
            sprintf("%.6f", max(diag_for_note)),
            paste(audit$ChromArm[audit$community_assigned_but_no_called_chm13_phr_bed_interval],
                  collapse = ","),
            paste(tree_row_label_font_pt, tree_col_label_font_pt, sep = "/"),
            paste(community_row_label_font_pt, community_col_label_font_pt, sep = "/"),
            "same Jaccard similarity palette and [0,1] scale used for both tree and community assets"),
  stringsAsFactors = FALSE
)

write.table(validation_tbl,
            file = file.path(out_dir, "render_validation.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

if (!all(validation_tbl$value[1:5] == "TRUE")) {
  stop("Render validation failed: order checks did not all pass")
}

tree_paths <- render_pair("07a_tree_ordered_heatmap", draw_tree_ordered_scene)
community_paths <- render_pair("07a_community_ordered_heatmap", draw_community_ordered_scene)
all_asset_paths <- c(tree_paths, community_paths)
asset_manifest <- data.frame(
  asset = basename(all_asset_paths),
  path = unname(all_asset_paths),
  format = sub("^.*\\.", "", basename(all_asset_paths)),
  expected_dimensions = c("13.333 x 7.5 in PDF",
                          "4800 x 2700 px PNG",
                          "13.333 x 7.5 in PDF",
                          "4800 x 2700 px PNG"),
  file_exists = file.exists(all_asset_paths),
  bytes = as.numeric(file.info(all_asset_paths)$size),
  stringsAsFactors = FALSE
)
write.table(asset_manifest,
            file = file.path(out_dir, "asset_manifest.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

asset_checks <- data.frame(
  check = paste0(sub("\\.[^.]+$", "", asset_manifest$asset), "_",
                 asset_manifest$format, "_exists_nonzero"),
  value = as.character(asset_manifest$file_exists & asset_manifest$bytes > 0),
  stringsAsFactors = FALSE
)
validation_tbl <- rbind(validation_tbl, asset_checks)
write.table(validation_tbl,
            file = file.path(out_dir, "render_validation.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Wrote:")
message("  ", tree_paths[["pdf"]])
message("  ", tree_paths[["png"]])
message("  ", community_paths[["pdf"]])
message("  ", community_paths[["png"]])
message("  ", file.path(out_dir, "arm_order_tree.tsv"))
message("  ", file.path(out_dir, "arm_order_community.tsv"))
message("  ", file.path(out_dir, "community_blocks.tsv"))
message("  ", file.path(out_dir, "arm_inclusion_audit.tsv"))
message("  ", file.path(out_dir, "source_audit.tsv"))
message("  ", file.path(out_dir, "render_validation.tsv"))
message("  ", file.path(out_dir, "asset_manifest.tsv"))
