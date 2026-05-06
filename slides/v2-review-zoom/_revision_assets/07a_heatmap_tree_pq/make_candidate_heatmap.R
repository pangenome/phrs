#!/usr/bin/env Rscript

# Candidate slide 07a heatmap:
# - same 41 x 41 arm-level Jaccard distance matrix as Fig 1c
# - rows and columns ordered by the UPGMA tree
# - left dendrogram is the ordering tree
# - p-arm labels are red and q-arm labels are blue

dist_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
leiden_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
out_dir <- "slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

D <- as.matrix(read.table(dist_path,
                          header = TRUE,
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE))
D <- D[rownames(D), rownames(D)]
D <- (D + t(D)) / 2

leiden <- read.table(leiden_path,
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = FALSE)
community_of <- setNames(leiden$Community, leiden$ChromArm)

hc <- hclust(as.dist(D), method = "average")
tree_order <- hc$labels[hc$order]
n <- length(tree_order)

arm_side <- function(x) {
  ifelse(grepl("_p$", x), "p", "q")
}

short_arm <- function(x) {
  sub("^chr", "", gsub("_", "", x))
}

arm_cols <- c(p = "#D13F3F", q = "#245FAD")
side <- arm_side(tree_order)
label_cols <- arm_cols[side]
labels_short <- short_arm(tree_order)

order_tbl <- data.frame(
  position = seq_along(tree_order),
  ChromArm = tree_order,
  label = labels_short,
  arm_side = side,
  Leiden = unname(community_of[tree_order]),
  stringsAsFactors = FALSE
)
write.table(order_tbl,
            file = file.path(out_dir, "candidate_upgma_tree_order.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

palette_heat <- colorRampPalette(c("#f7fbff", "#fddc7a", "#f36d33", "#b20d2a"))(256)
breaks <- seq(0, 1, length.out = length(palette_heat) + 1)

similarity <- 1 - D
similarity[similarity < 0] <- 0
similarity[similarity > 1] <- 1

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

plot_candidate <- function() {
  leaf_pos <- setNames(rev(seq_along(tree_order)), tree_order)
  row_y <- unname(leaf_pos[tree_order])
  segs <- make_left_segments(as.dendrogram(hc), leaf_pos)
  max_height <- max(segs$x, segs$xend)

  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.15, 5.15, 1.35))
  par(oma = c(0, 0, 1.5, 0), family = "sans")

  par(mar = c(5.7, 0.1, 0.6, 0.05))
  plot.new()
  plot.window(xlim = c(max_height, 0),
              ylim = c(0.5, n + 0.5),
              xaxs = "i",
              yaxs = "i")
  segments(segs$x, segs$y, segs$xend, segs$yend,
           col = "#3c3c3c",
           lwd = 0.75,
           lend = "square")

  par(mar = c(5.7, 0.25, 0.6, 3.0))
  plot.new()
  plot.window(xlim = c(0.5, n + 4.0),
              ylim = c(-4.4, n + 0.5),
              xaxs = "i",
              yaxs = "i")

  vals <- similarity[tree_order, tree_order]
  idx <- matrix(findInterval(as.vector(vals), breaks, all.inside = TRUE),
                nrow = nrow(vals),
                ncol = ncol(vals))
  img <- matrix(palette_heat[idx], nrow = nrow(idx), ncol = ncol(idx))

  rasterImage(as.raster(img), 0.5, 0.5, n + 0.5, n + 0.5)
  rect(0.5, 0.5, n + 0.5, n + 0.5, border = "#222222", lwd = 0.75)

  grid_col <- grDevices::adjustcolor("#555555", alpha.f = 0.23)
  abline(v = seq(0.5, n + 0.5, by = 1), col = grid_col, lwd = 0.35)
  abline(h = seq(0.5, n + 0.5, by = 1), col = grid_col, lwd = 0.35)

  comm_ids <- sort(unique(na.omit(community_of)))
  comm_cols <- setNames(grDevices::hcl.colors(length(comm_ids), "Set 3"), comm_ids)
  points(rep(n + 0.72, n),
         row_y,
         pch = 15,
         cex = 0.62,
         col = comm_cols[community_of[tree_order]],
         xpd = NA)
  text(n + 0.72, n + 1.4,
       labels = "Leiden",
       srt = 90,
       cex = 0.55,
       col = "#444444",
       xpd = NA)

  text(n + 1.08,
       row_y,
       labels = labels_short,
       adj = 0,
       cex = 0.58,
       col = label_cols,
       xpd = NA)

  text(seq_len(n),
       -0.35,
       labels = labels_short,
       srt = 90,
       adj = 1,
       cex = 0.58,
       col = label_cols,
       xpd = NA)

  text(n / 2,
       -4.0,
       labels = "Columns use the same UPGMA order as rows",
       cex = 0.68,
       col = "#333333",
       xpd = NA)

  par(mar = c(5.7, 0.2, 0.6, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  text(0.0, 0.95, "Label color", adj = 0, cex = 0.82, font = 2)
  points(0.08, 0.87, pch = 15, col = arm_cols["p"], cex = 1.1)
  text(0.17, 0.87, "p arm", adj = 0, cex = 0.76, col = arm_cols["p"])
  points(0.08, 0.81, pch = 15, col = arm_cols["q"], cex = 1.1)
  text(0.17, 0.81, "q arm", adj = 0, cex = 0.76, col = arm_cols["q"])

  text(0.0, 0.70, "Cell color", adj = 0, cex = 0.82, font = 2)
  cb <- matrix(rev(palette_heat), ncol = 1)
  rasterImage(as.raster(cb), 0.07, 0.22, 0.22, 0.66)
  rect(0.07, 0.22, 0.22, 0.66, border = "#333333", lwd = 0.6)
  text(0.28, 0.66, "1.0", adj = 0, cex = 0.68)
  text(0.28, 0.44, "0.5", adj = 0, cex = 0.68)
  text(0.28, 0.22, "0.0", adj = 0, cex = 0.68)
  text(0.0, 0.14,
       labels = "Jaccard similarity\n= 1 - distance",
       adj = 0,
       cex = 0.68,
       col = "#333333")
  text(0.0, 0.04,
       labels = "Right ticks mark\nLeiden community.",
       adj = 0,
       cex = 0.62,
       col = "#555555")

  mtext("Candidate 07a: tree-left UPGMA order on the 41-arm Jaccard matrix",
        outer = TRUE,
        line = 0.25,
        cex = 0.86,
        font = 2)
}

png(file.path(out_dir, "candidate_heatmap_upgma_tree_left_pq.png"),
    width = 2700,
    height = 2250,
    res = 300,
    type = "cairo")
plot_candidate()
dev.off()

pdf(file.path(out_dir, "candidate_heatmap_upgma_tree_left_pq.pdf"),
    width = 9.0,
    height = 7.5,
    useDingbats = FALSE)
plot_candidate()
dev.off()
