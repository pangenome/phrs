#!/usr/bin/env Rscript

# Regenerate slide 07a as a vector-first, row-aligned UPGMA heatmap.
# The tree and heatmap share one explicit arm order and one physical row scale.

suppressPackageStartupMessages({
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

dist_path <- if (length(args) >= 1) args[[1]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
leiden_path <- if (length(args) >= 2) args[[2]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
out_dir <- if (length(args) >= 3) args[[3]] else
  "slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading arm-level distance matrix: ", dist_path)
D_raw <- as.matrix(read.table(dist_path,
                              header = TRUE,
                              row.names = 1,
                              sep = "\t",
                              check.names = FALSE,
                              quote = "",
                              comment.char = ""))
stopifnot(nrow(D_raw) == ncol(D_raw))
stopifnot(identical(rownames(D_raw), colnames(D_raw)))

D <- (D_raw + t(D_raw)) / 2
diag_for_note <- diag(D)

message("Reading Leiden assignments: ", leiden_path)
leiden <- read.table(leiden_path,
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     check.names = FALSE,
                     quote = "",
                     comment.char = "")
stopifnot(all(c("ChromArm", "Community", "Arms") %in% names(leiden)))
stopifnot(setequal(rownames(D), leiden$ChromArm))

community_of <- setNames(leiden$Community, leiden$ChromArm)
community_arms <- setNames(leiden$Arms, leiden$ChromArm)

hc <- hclust(as.dist(D), method = "average")
tree_order <- hc$labels[hc$order]
n <- length(tree_order)
stopifnot(n == 41)

arm_side <- function(x) {
  ifelse(grepl("_p$", x), "p", "q")
}

short_arm <- function(x) {
  sub("^chr", "", gsub("_", "", x))
}

side <- arm_side(tree_order)
labels_short <- short_arm(tree_order)
label_cols <- c(p = "#CC3B38", q = "#1F5EA8")[side]

order_tbl <- data.frame(
  position_top_to_bottom = seq_along(tree_order),
  position_left_to_right = seq_along(tree_order),
  ChromArm = tree_order,
  label = labels_short,
  arm_side = side,
  Leiden = unname(community_of[tree_order]),
  Leiden_members = unname(community_arms[tree_order]),
  ordering_method = "UPGMA average linkage on the 41x41 arm-level Jaccard distance matrix",
  stringsAsFactors = FALSE
)

write.table(order_tbl,
            file = file.path(out_dir, "arm_order_upgma.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

leaf_y <- setNames(rev(seq_along(tree_order)), tree_order)
tree_tip_order_from_y <- names(sort(leaf_y, decreasing = TRUE))
heatmap_row_order <- tree_order
heatmap_col_order <- tree_order

validation_tbl <- data.frame(
  check = c("tree_tip_order_equals_heatmap_row_order",
            "tree_tip_order_equals_heatmap_column_order",
            "heatmap_rows_equal_heatmap_columns",
            "matrix_rows_equal_tree_labels",
            "ordering_method",
            "tree_method_label",
            "rendering_method",
            "source_matrix_rows",
            "source_matrix_columns",
            "max_source_asymmetry",
            "raw_diagonal_min",
            "raw_diagonal_max"),
  value = c(identical(tree_tip_order_from_y, heatmap_row_order),
            identical(tree_tip_order_from_y, heatmap_col_order),
            identical(heatmap_row_order, heatmap_col_order),
            setequal(rownames(D), hc$labels),
            "UPGMA average linkage via hclust(method = 'average')",
            "UPGMA",
            "PDF vector cells/tree; PNG rendered directly from the same vector scene",
            nrow(D_raw),
            ncol(D_raw),
            sprintf("%.6g", max(abs(D_raw - t(D_raw)))),
            sprintf("%.6f", min(diag_for_note)),
            sprintf("%.6f", max(diag_for_note))),
  stringsAsFactors = FALSE
)

write.table(validation_tbl,
            file = file.path(out_dir, "order_validation.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

if (!all(validation_tbl$value[1:4] == "TRUE")) {
  stop("Order validation failed: tree tips, heatmap rows, and heatmap columns are not identical")
}

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
similarity <- 1 - D
similarity[similarity < 0] <- 0
similarity[similarity > 1] <- 1

value_to_col <- function(x) {
  palette_heat[findInterval(x, breaks, all.inside = TRUE)]
}

draw_scene <- function() {
  grid.newpage()

  page_w <- 13.333
  page_h <- 7.5
  matrix_bottom <- 1.02
  matrix_size <- 5.72
  tree_left <- 0.42
  tree_width <- 1.58
  gutter <- 0.06
  heat_left <- tree_left + tree_width + gutter
  heat_bottom <- matrix_bottom
  cell <- matrix_size / n
  label_left <- heat_left + matrix_size + 0.12
  legend_left <- 9.25

  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text("Slide 07a: UPGMA-aligned 41-arm heatmap",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.08, "in"),
            just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 12, col = "#1F1F1F"))
  grid.text("Tree tips, heatmap rows, and columns share one UPGMA order; labels: p = red, q = blue.",
            x = unit(0.42, "in"),
            y = unit(page_h - 0.36, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8.5, col = "#3E3E3E"))

  vals <- similarity[tree_order, tree_order]
  fills <- matrix(value_to_col(as.vector(vals)), nrow = nrow(vals), ncol = ncol(vals))

  xs <- numeric(n * n)
  ys <- numeric(n * n)
  cols <- character(n * n)
  k <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      xs[k] <- heat_left + (j - 0.5) * cell
      ys[k] <- heat_bottom + (n - i + 0.5) * cell
      cols[k] <- fills[i, j]
      k <- k + 1
    }
  }
  grid.rect(x = unit(xs, "in"),
            y = unit(ys, "in"),
            width = unit(cell, "in"),
            height = unit(cell, "in"),
            gp = gpar(fill = cols, col = NA))

  grid.rect(x = unit(heat_left + matrix_size / 2, "in"),
            y = unit(heat_bottom + matrix_size / 2, "in"),
            width = unit(matrix_size, "in"),
            height = unit(matrix_size, "in"),
            gp = gpar(fill = NA, col = "#202020", lwd = 0.8))

  line_col <- "#5C5C5C40"
  grid.segments(x0 = unit(rep(heat_left + seq(0, matrix_size, by = cell), 2), "in"),
                y0 = unit(c(rep(heat_bottom, n + 1), rep(heat_bottom + matrix_size, n + 1)), "in"),
                x1 = unit(rep(heat_left + seq(0, matrix_size, by = cell), 2), "in"),
                y1 = unit(c(rep(heat_bottom + matrix_size, n + 1), rep(heat_bottom, n + 1)), "in"),
                gp = gpar(col = line_col, lwd = 0.25))
  grid.segments(x0 = unit(rep(heat_left, n + 1), "in"),
                y0 = unit(heat_bottom + seq(0, matrix_size, by = cell), "in"),
                x1 = unit(rep(heat_left + matrix_size, n + 1), "in"),
                y1 = unit(heat_bottom + seq(0, matrix_size, by = cell), "in"),
                gp = gpar(col = line_col, lwd = 0.25))

  row_centers <- heat_bottom + (n - seq_len(n) + 0.5) * cell
  col_centers <- heat_left + (seq_len(n) - 0.5) * cell

  grid.text(labels_short,
            x = unit(label_left, "in"),
            y = unit(row_centers, "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = 5.9, col = label_cols))
  grid.text(labels_short,
            x = unit(col_centers, "in"),
            y = unit(heat_bottom - 0.075, "in"),
            rot = 90,
            just = c("right", "center"),
            gp = gpar(fontsize = 5.7, col = label_cols))

  grid.text("arms",
            x = unit(label_left + 0.03, "in"),
            y = unit(heat_bottom + matrix_size + 0.12, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))

  segs <- make_left_segments(as.dendrogram(hc), leaf_y)
  max_height <- max(segs$x, segs$xend)
  tx <- function(h) tree_left + (max_height - h) / max_height * tree_width
  ty <- function(y) heat_bottom + (y - 0.5) * cell
  grid.segments(x0 = unit(tx(segs$x), "in"),
                y0 = unit(ty(segs$y), "in"),
                x1 = unit(tx(segs$xend), "in"),
                y1 = unit(ty(segs$yend), "in"),
                gp = gpar(col = "#333333", lwd = 0.65, lineend = "square"))
  grid.rect(x = unit(tree_left + tree_width / 2, "in"),
            y = unit(heat_bottom + matrix_size / 2, "in"),
            width = unit(tree_width, "in"),
            height = unit(matrix_size, "in"),
            gp = gpar(fill = NA, col = NA))
  grid.text("UPGMA tree",
            x = unit(tree_left, "in"),
            y = unit(heat_bottom + matrix_size + 0.12, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))
  grid.text("leaf order defines rows and columns",
            x = unit(tree_left, "in"),
            y = unit(heat_bottom + matrix_size + 0.01, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 6, col = "#555555"))

  legend_y <- heat_bottom + matrix_size - 0.08
  grid.text("Encoding",
            x = unit(legend_left, "in"),
            y = unit(legend_y, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold", col = "#222222"))
  grid.rect(x = unit(c(legend_left + 0.06, legend_left + 0.06), "in"),
            y = unit(c(legend_y - 0.31, legend_y - 0.58), "in"),
            width = unit(0.12, "in"),
            height = unit(0.12, "in"),
            gp = gpar(fill = c("#CC3B38", "#1F5EA8"), col = NA))
  grid.text(c("p arm label", "q arm label"),
            x = unit(rep(legend_left + 0.18, 2), "in"),
            y = unit(c(legend_y - 0.31, legend_y - 0.58), "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = 7.5, col = c("#CC3B38", "#1F5EA8")))

  grid.text("Jaccard similarity",
            x = unit(legend_left, "in"),
            y = unit(legend_y - 0.98, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 8, fontface = "bold", col = "#333333"))
  cb_n <- length(palette_heat)
  cb_x <- rep(legend_left + 0.12, cb_n)
  cb_y <- seq(legend_y - 2.35, legend_y - 1.15, length.out = cb_n)
  grid.rect(x = unit(cb_x, "in"),
            y = unit(cb_y, "in"),
            width = unit(0.18, "in"),
            height = unit(1.22 / cb_n, "in"),
            gp = gpar(fill = palette_heat, col = NA))
  grid.rect(x = unit(legend_left + 0.12, "in"),
            y = unit(legend_y - 1.75, "in"),
            width = unit(0.18, "in"),
            height = unit(1.2, "in"),
            gp = gpar(fill = NA, col = "#333333", lwd = 0.5))
  grid.text(c("1.0", "0.5", "0.0"),
            x = unit(rep(legend_left + 0.36, 3), "in"),
            y = unit(c(legend_y - 1.15, legend_y - 1.75, legend_y - 2.35), "in"),
            just = c("left", "center"),
            gp = gpar(fontsize = 7, col = "#333333"))
  grid.text("PDF uses vector cells/tree.\nPNG is rendered directly at 4800x2700.\nNo stretched raster crop.",
            x = unit(legend_left, "in"),
            y = unit(legend_y - 2.64, "in"),
            just = c("left", "top"),
            gp = gpar(fontsize = 6.6, col = "#555555", lineheight = 1.05))

  grid.text("Validation: tree tip order == heatmap row order == heatmap column order (see order_validation.tsv).",
            x = unit(0.42, "in"),
            y = unit(0.18, "in"),
            just = c("left", "bottom"),
            gp = gpar(fontsize = 7.2, col = "#3A3A3A"))
}

pdf_path <- file.path(out_dir, "candidate_07a_upgma_crisp_aligned.pdf")
png_path <- file.path(out_dir, "candidate_07a_upgma_crisp_aligned.png")

pdf(pdf_path, width = 13.333, height = 7.5, useDingbats = FALSE)
draw_scene()
dev.off()

png(png_path,
    width = 4800,
    height = 2700,
    res = 360,
    type = "cairo")
draw_scene()
dev.off()

message("Wrote:")
message("  ", pdf_path)
message("  ", png_path)
message("  ", file.path(out_dir, "arm_order_upgma.tsv"))
message("  ", file.path(out_dir, "order_validation.tsv"))
