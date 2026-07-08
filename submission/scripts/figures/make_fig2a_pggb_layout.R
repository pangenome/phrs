#!/usr/bin/env Rscript

# Manuscript Fig 2a -- PGGB/odgi 2D layout of the all-subtelomere PHR graph,
# main connected component, nodes colored by arm-level Leiden k15 community.
# Location-aware: reads gzipped inputs from data/ at the repo root, writes
# submission/fig/MainFigures/Fig2a_pggb_layout.png. Base R only, no moosefs.
#
# Inputs (data/):
#   fig2a_pggb_layout.og.lay.tsv.gz   odgi `...smooth.final.og.lay.tsv`
#                                     (columns idx,X,Y,component). Each graph node
#                                     contributes two handle-end rows, so
#                                     node_id = idx %/% 2 + 1.
#   fig2a_node_community.tsv.gz        node<TAB>community. Dominant arm-level Leiden
#                                     community among the PHR paths crossing each
#                                     node, derived once from the 5.9 GB pangenome
#                                     GFA P-lines:
#     awk -F'\t' '
#       FNR==NR{ if($1!="ChromArm" && $1!="") arm2comm[$1]=$2; next }   # arm-leiden-k15
#       $1=="P"{ if(match($2,/_chr[0-9XYxy]+_[pq]arm$/)){s=substr($2,RSTART+1);sub(/arm$/,"",s);a=s}else a="?";
#                c=arm2comm[a]; if(c=="")c="NA"; n=split($3,t,",");
#                for(i=1;i<=n;i++){x=t[i];sub(/[+-]$/,"",x);cnt[x SUBSEP c]++} }
#       END{ for(k in cnt){split(k,z,SUBSEP);v=cnt[k];if(v>b[z[1]]){b[z[1]]=v;bc[z[1]]=z[2]}}
#            for(nd in b) print nd"\t"bc[nd] }' \
#       data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv  <graph>.smooth.final.gfa
#
# The layout coordinates come from `odgi layout` (pggb v0.7.4, odgi v0.9.2); this
# script renders the layout NODES in base R (not `odgi draw`), as the published
# panel does. X/Y are swapped to rotate the tall main component into a 16:9 frame.

args <- commandArgs(trailingOnly = TRUE)
.cmd_args  <- commandArgs(trailingOnly = FALSE)
.this_file <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
script_dir <- if (length(.this_file)) normalizePath(dirname(.this_file)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."))

val <- function(flag, default) {
  h <- which(args == flag)
  if (length(h) && h[length(h)] < length(args)) args[h[length(h)] + 1] else default
}
layout_path    <- val("--layout",    file.path(repo_root, "data/fig2a_pggb_layout.og.lay.tsv.gz"))
node_comm_path <- val("--node-comm", file.path(repo_root, "data/fig2a_node_community.tsv.gz"))
out_dir        <- val("--out-dir",   file.path(repo_root, "submission/fig/MainFigures"))
out_name       <- val("--out-name",  "Fig2a_pggb_layout.png")

# --mono [color]: single-color, non-community layout with no legend (defaults to
# charcoal). Omit --mono for the community-colored manuscript panel.
mono <- any(args == "--mono")
mono_color <- if (mono) { v <- val("--mono", "#111111"); if (grepl("^--", v)) "#111111" else v } else NA_character_

point_alpha <- if (mono) 0.30 else 0.65
point_cex   <- if (mono) 0.12 else 0.42
background  <- "white"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(layout_path))

lay <- read.delim(gzfile(layout_path), sep = "\t", stringsAsFactors = FALSE)
# main connected component = the one with the most layout nodes (component 8)
comp_main <- as.integer(names(which.max(table(lay$component))))
lay <- lay[lay$component == comp_main, c("idx", "X", "Y")]
lay$node <- lay$idx %/% 2L + 1L

# Canonical paper community palette -- identical to
# submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R and
# make_fig4b_porec_community.R (Leiden C1..C15 in numeric order). NA = gray.
pal <- c(C1="#2E6FBB", C2="#D95F02", C3="#1B9E77", C4="#7570B3", C5="#E7298A",
         C6="#66A61E", C7="#E6AB02", C8="#A6761D", C9="#1F78B4", C10="#B2DF8A",
         C11="#FB9A99", C12="#FDBF6F", C13="#CAB2D6", C14="#6A3D9A", C15="#B15928",
         `NA`="#D9D9D9")
levs <- names(pal)

if (!mono) {
  stopifnot(file.exists(node_comm_path))
  nc <- read.delim(gzfile(node_comm_path), sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE, col.names = c("node", "community"))
  lay$community <- nc$community[match(lay$node, nc$node)]
  lay$community[is.na(lay$community)] <- "NA"
  lay$community <- factor(lay$community, levels = levs)
  # draw NA first, then rarer communities, largest on top
  lay <- lay[order(lay$community != "NA", table(lay$community)[as.character(lay$community)]), ]
  col_vec <- grDevices::adjustcolor(pal[as.character(lay$community)], alpha.f = point_alpha)
} else {
  col_vec <- grDevices::adjustcolor(mono_color, alpha.f = point_alpha)
}

plot_x <- lay$Y; plot_y <- lay$X   # rotate as in the published panel

png_path <- file.path(out_dir, out_name)
png(png_path, width = 1920, height = 820, res = 144, bg = background, type = "cairo")
par(mar = c(0.18, 0.18, 0.18, 0.18), xaxs = "i", yaxs = "i", family = "sans")
plot(plot_x, plot_y, pch = 16, cex = point_cex, col = col_vec,
     axes = FALSE, ann = FALSE, asp = 1)

# Surface the labelled communities by redrawing them last, on top of the larger
# communities that share their layout positions and otherwise occlude them
# (C15/PAR1 is tiny -> enlarge its points so the specks are visible).
if (!mono) {
  for (cc in c("C11", "C2", "C15")) {
    s <- which(lay$community == cc)
    points(plot_x[s], plot_y[s], pch = 16,
           cex = point_cex * if (cc == "C15") 1.8 else 1.0, col = pal[[cc]])
  }
}

if (!mono) {
  present <- levs[levs %in% as.character(unique(lay$community))]
  present_named <- present[present != "NA"]
  present_named <- present_named[order(as.integer(sub("^C", "", present_named)))]
  leg     <- c(present_named, if ("NA" %in% present) "no community")
  leg_col <- c(pal[present_named], if ("NA" %in% present) pal[["NA"]])
  usr <- par("usr")
  d <- 0.04 * (usr[4] - usr[3])   # physical inset from the top edge
  # shifted right of the left edge to leave room for the overpic panel letter "A"
  lgd <- legend(x = usr[1] + d + 0.06 * (usr[2] - usr[1]), y = usr[4] - d, xjust = 0, yjust = 1,
         legend = leg, col = leg_col, pch = 16, pt.cex = 1.7,
         ncol = 2, bty = "n", cex = 1.25, text.col = "#1f1f1f", title = "Community")

  # Gene/feature labels (cf. Fig 3): each is auto-placed in the nearest OPEN
  # pocket around a same-colour node -- a spot whose text box hits no graph
  # points and does not collide with the legend or another label.
  hx <- 0.0028 * (usr[2] - usr[1]); hy <- 0.0028 * (usr[4] - usr[3])
  cexL <- 1.05
  darken <- function(col, f = 0.7) {
    v <- grDevices::col2rgb(col)[, 1] * f
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
  }
  # anchor = a point in the community's own-colour-dominant region
  # pref = preferred direction for the label (radians; 0=right, pi/2=up, pi=left, -pi/2=down)
  anchors <- list(
    C1  = list(g = "DUX4 / D4Z4", at = c(152932,  88490), pref = -pi / 2),
    C11 = list(g = "OR4F",        at = c(180061, 153159), fixed = c(205000, 156500)),
    C2  = list(g = "TUBB8B",      at = c(218888,  92437)),
    C7  = list(g = "rDNA (acro)", at = c(300000, 108000), pref = -pi / 2),
    C14 = list(g = "PAR2",        at = c( 55000,  42000), pref = pi),
    C15 = list(g = "PAR1 (SHOX)", at = c(121000,  62000), pref = pi))
  box_overlap <- function(b, bs) {
    for (o in bs) if (b[1] < o[3] && b[3] > o[1] && b[2] < o[4] && b[4] > o[2]) return(TRUE)
    FALSE
  }
  boxes <- list(c(lgd$rect$left, lgd$rect$top - lgd$rect$h,
                  lgd$rect$left + lgd$rect$w, lgd$rect$top))   # legend box
  for (cc in names(anchors)) {
    info <- anchors[[cc]]; colc <- darken(pal[[cc]])
    sel <- which(lay$community == cc)
    j   <- sel[which.min((plot_x[sel] - info$at[1])^2 + (plot_y[sel] - info$at[2])^2)]
    ax <- plot_x[j]; ay <- plot_y[j]
    # measure the box for the BOLD text actually drawn, plus a small clearance ring
    pad <- 0.7 * strheight(info$g, cex = cexL, font = 2)
    w <- strwidth(info$g, cex = cexL, font = 2) + 2 * pad
    h <- strheight(info$g, cex = cexL, font = 2) + 2 * pad
    R <- 0.40 * (usr[2] - usr[1])
    loc <- which(abs(plot_x - ax) < R + w & abs(plot_y - ay) < R + h)
    lpx <- plot_x[loc]; lpy <- plot_y[loc]
    base <- seq(0, 2 * pi, length.out = 37)[-1]
    ang  <- if (is.null(info$pref)) base
            else base[order(abs(((base - info$pref + pi) %% (2 * pi)) - pi))]
    place <- c(ax, ay)
    if (!is.null(info$fixed)) {
      place <- info$fixed
    } else for (r in seq(0.55 * h + 0.2 * w, R, length.out = 34)) {
      hit <- FALSE
      for (a in ang) {
        cx <- ax + r * cos(a); cy <- ay + r * sin(a)
        b  <- c(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        if (b[1] < usr[1] || b[3] > usr[2] || b[2] < usr[3] || b[4] > usr[4]) next
        if (any(lpx > b[1] & lpx < b[3] & lpy > b[2] & lpy < b[4])) next
        if (box_overlap(b, boxes)) next
        place <- c(cx, cy); hit <- TRUE; break
      }
      if (hit) break
    }
    lx <- place[1]; ly <- place[2]
    boxes[[length(boxes) + 1]] <- c(lx - w / 2, ly - h / 2, lx + w / 2, ly + h / 2)
    segments(ax, ay, lx, ly, col = "black", lwd = 1.1)
    points(ax, ay, pch = 16, cex = 0.9, col = "black")
    for (dx in c(-hx, hx)) for (dy in c(-hy, hy))
      text(lx + dx, ly + dy, info$g, col = "white", font = 2, cex = cexL)
    text(lx, ly, info$g, col = "black", font = 2, cex = cexL)
  }
}
invisible(dev.off())

cat("wrote ", png_path, "  (", if (mono) paste0("mono ", mono_color) else "community",
    "; ", nrow(lay), " nodes, component ", comp_main, ")\n", sep = "")
