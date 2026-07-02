#!/usr/bin/env Rscript
# Fig2A variant: PGGB/odgi layout main component (component 8), nodes colored by
# arm-level Leiden k15 community. node_id = idx %/% 2 + 1 (odgi layout TSV encodes
# two handle-ends per node). node -> community from GFA P-line membership
# (dominant community among PHR paths crossing each node).
#
# Inputs (in this directory):
#   layout.tsv      odgi ...smooth.final.og.lay.tsv (idx,X,Y,component)
#   node_comm.tsv   node<TAB>community, from GFA P-lines x arm-leiden-k15 assignments

args <- commandArgs(trailingOnly = TRUE)
val <- function(flag, default) { h <- which(args == flag); if (length(h)==0 || h[length(h)]==length(args)) default else args[h[length(h)]+1] }

layout_tsv   <- val("--layout-tsv", "layout.tsv")
node_comm    <- val("--node-comm", "node_comm.tsv")
out_png      <- val("--out", "fig2a_by_community.png")
component_id <- as.integer(val("--component", "8"))
point_alpha  <- as.numeric(val("--point-alpha", "0.65"))
point_cex    <- as.numeric(val("--point-cex", "0.30"))
background   <- val("--background-color", "white")
border_color <- val("--border-color", "#b8c0cc")

stopifnot(file.exists(layout_tsv), file.exists(node_comm))

lay <- read.delim(layout_tsv, sep = "\t", stringsAsFactors = FALSE)
lay <- lay[lay$component == component_id, c("idx","X","Y")]
lay$node <- lay$idx %/% 2L + 1L

nc <- read.delim(node_comm, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("node","community"))
lay$community <- nc$community[match(lay$node, nc$node)]
lay$community[is.na(lay$community)] <- "NA"

# Canonical paper community palette (from submission/scripts/figures/
# make_fig2bc_jaccard_heatmaps.R; same in make_fig4b_porec_community.R),
# assigned to Leiden C1..C15 in numeric order. NA = gray.
pal <- c(C1="#2E6FBB", C2="#D95F02", C3="#1B9E77", C4="#7570B3", C5="#E7298A",
         C6="#66A61E", C7="#E6AB02", C8="#A6761D", C9="#1F78B4", C10="#B2DF8A",
         C11="#FB9A99", C12="#FDBF6F", C13="#CAB2D6", C14="#6A3D9A", C15="#B15928",
         `NA`="#D9D9D9")
levs <- names(pal)
lay$community <- factor(lay$community, levels = levs)

# draw NA first, then rarer communities, big ones on top
ord <- order(lay$community != "NA", table(lay$community)[as.character(lay$community)])
lay <- lay[ord, ]

plot_x <- lay$Y; plot_y <- lay$X   # rotate as in the published panel

png(out_png, width = 1920, height = 1080, res = 144, bg = background, type = "cairo")
par(mar = c(0.18,0.18,0.18,0.18), xaxs="i", yaxs="i", family="sans")
plot(plot_x, plot_y, pch = 16, cex = point_cex,
     col = grDevices::adjustcolor(pal[as.character(lay$community)], alpha.f = point_alpha),
     axes = FALSE, ann = FALSE, asp = 1)
usr <- par("usr"); rect(usr[1],usr[3],usr[2],usr[4], border = border_color, lwd = 2)

present <- levs[levs %in% as.character(unique(lay$community))]
present_named <- present[present != "NA"]
present_named <- present_named[order(as.integer(sub("^C","",present_named)))]
leg <- c(present_named, if ("NA" %in% present) "no community")
leg_col <- c(pal[present_named], if ("NA" %in% present) pal[["NA"]])
legend("topleft", legend = leg, col = leg_col, pch = 16, pt.cex = 1.1,
       ncol = 2, bty = "n", cex = 0.8, text.col = "#1f1f1f",
       title = "Community")
invisible(dev.off())

cat("wrote ", out_png, "\n", sep = "")
cat("nodes plotted: ", nrow(lay), "\n", sep = "")
print(sort(table(droplevels(lay$community)), decreasing = TRUE))
