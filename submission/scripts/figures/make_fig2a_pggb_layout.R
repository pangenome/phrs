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

layout_path <- if (length(args) >= 1) args[[1]] else
  file.path(repo_root, "data/fig2a_pggb_layout.og.lay.tsv.gz")
node_comm_path <- if (length(args) >= 2) args[[2]] else
  file.path(repo_root, "data/fig2a_node_community.tsv.gz")
out_dir <- if (length(args) >= 3) args[[3]] else
  file.path(repo_root, "submission/fig/MainFigures")

point_alpha <- 0.65
point_cex   <- 0.30
background  <- "white"
border_col  <- "#b8c0cc"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(layout_path), file.exists(node_comm_path))

lay <- read.delim(gzfile(layout_path), sep = "\t", stringsAsFactors = FALSE)
# main connected component = the one with the most layout nodes (component 8)
comp_main <- as.integer(names(which.max(table(lay$component))))
lay <- lay[lay$component == comp_main, c("idx", "X", "Y")]
lay$node <- lay$idx %/% 2L + 1L

nc <- read.delim(gzfile(node_comm_path), sep = "\t", header = FALSE,
                 stringsAsFactors = FALSE, col.names = c("node", "community"))
lay$community <- nc$community[match(lay$node, nc$node)]
lay$community[is.na(lay$community)] <- "NA"

# Canonical paper community palette -- identical to
# submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R and
# make_fig4b_porec_community.R (Leiden C1..C15 in numeric order). NA = gray.
pal <- c(C1="#2E6FBB", C2="#D95F02", C3="#1B9E77", C4="#7570B3", C5="#E7298A",
         C6="#66A61E", C7="#E6AB02", C8="#A6761D", C9="#1F78B4", C10="#B2DF8A",
         C11="#FB9A99", C12="#FDBF6F", C13="#CAB2D6", C14="#6A3D9A", C15="#B15928",
         `NA`="#D9D9D9")
levs <- names(pal)
lay$community <- factor(lay$community, levels = levs)

# draw NA first, then rarer communities, largest on top
lay <- lay[order(lay$community != "NA", table(lay$community)[as.character(lay$community)]), ]

plot_x <- lay$Y; plot_y <- lay$X   # rotate as in the published panel

png_path <- file.path(out_dir, "Fig2a_pggb_layout.png")
png(png_path, width = 1920, height = 1080, res = 144, bg = background, type = "cairo")
par(mar = c(0.18, 0.18, 0.18, 0.18), xaxs = "i", yaxs = "i", family = "sans")
plot(plot_x, plot_y, pch = 16, cex = point_cex,
     col = grDevices::adjustcolor(pal[as.character(lay$community)], alpha.f = point_alpha),
     axes = FALSE, ann = FALSE, asp = 1)
usr <- par("usr"); rect(usr[1], usr[3], usr[2], usr[4], border = border_col, lwd = 2)

present <- levs[levs %in% as.character(unique(lay$community))]
present_named <- present[present != "NA"]
present_named <- present_named[order(as.integer(sub("^C", "", present_named)))]
leg     <- c(present_named, if ("NA" %in% present) "no community")
leg_col <- c(pal[present_named], if ("NA" %in% present) pal[["NA"]])
legend("topleft", legend = leg, col = leg_col, pch = 16, pt.cex = 1.1,
       ncol = 2, bty = "n", cex = 0.8, text.col = "#1f1f1f", title = "Community")
invisible(dev.off())

cat("wrote ", png_path, "\n", sep = "")
cat("main component: ", comp_main, "  nodes plotted: ", nrow(lay), "\n", sep = "")
print(sort(table(droplevels(lay$community)), decreasing = TRUE))
