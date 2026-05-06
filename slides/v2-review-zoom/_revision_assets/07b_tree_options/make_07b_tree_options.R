#!/usr/bin/env Rscript

# Generate larger candidate slide 07b tree assets from the arm-level distance
# matrix. This script intentionally writes only into the revision-assets folder;
# the Typst deck source is left untouched for the downstream renderer.

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)

matrix_path <- if (length(args) >= 1) args[[1]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
assignments_path <- if (length(args) >= 2) args[[2]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
out_dir <- if (length(args) >= 3) args[[3]] else
  "slides/v2-review-zoom/_revision_assets/07b_tree_options"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading matrix: ", matrix_path)
mat_df <- read.table(matrix_path,
                     header = TRUE,
                     sep = "\t",
                     check.names = FALSE,
                     row.names = 1,
                     quote = "",
                     comment.char = "")
D_raw <- as.matrix(mat_df)
stopifnot(nrow(D_raw) == ncol(D_raw))
stopifnot(identical(rownames(D_raw), colnames(D_raw)))

# The source matrix diagonal stores within-arm average distances, which are
# biologically useful but invalid as self-distances for NJ. Match the upstream
# NJ script by symmetrizing and zeroing the diagonal before ape::nj().
D <- (D_raw + t(D_raw)) / 2
diag(D) <- 0

message("Reading assignments: ", assignments_path)
assignments <- read.table(assignments_path,
                          header = TRUE,
                          sep = "\t",
                          check.names = FALSE,
                          quote = "",
                          comment.char = "",
                          stringsAsFactors = FALSE)
stopifnot(all(c("ChromArm", "Community", "Arms") %in% names(assignments)))
stopifnot(setequal(rownames(D), assignments$ChromArm))

comm_by_arm <- setNames(assignments$Community, assignments$ChromArm)
arms_by_comm <- setNames(assignments$Arms, assignments$ChromArm)

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0("chr", chroms, "_p"),
                            paste0("chr", chroms, "_q")))
present_arms <- rownames(D)
missing_arms <- setdiff(all_arms, present_arms)

missing_reason <- function(arm) {
  if (arm == "chr18_q") {
    return("zero inter-chromosomal PHR signal after removal of the NA18982#1 chr18_q/chrX PAR1 scaffold chimera")
  }
  "zero inter-chromosomal PHR signal at the >=95% identity PHR threshold"
}

arm_presence <- data.frame(
  chrom_arm = all_arms,
  in_41x41_matrix = all_arms %in% present_arms,
  leiden_community = ifelse(all_arms %in% names(comm_by_arm),
                            comm_by_arm[all_arms],
                            NA_character_),
  community_arms = ifelse(all_arms %in% names(arms_by_comm),
                          arms_by_comm[all_arms],
                          NA_character_),
  audit_status = ifelse(all_arms %in% present_arms,
                        "present",
                        "absent_by_construction"),
  reason = ifelse(all_arms %in% present_arms,
                  "present in arm-level distance matrix and NJ tree",
                  vapply(all_arms, missing_reason, character(1))),
  stringsAsFactors = FALSE
)

write.table(arm_presence,
            file = file.path(out_dir, "arm_presence.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

named_clades <- data.frame(
  community = c("C1", "C2", "C6", "C7", "C14", "C15"),
  short_label = c("C1 DUX4/D4Z4",
                  "C2 10p/18p",
                  "C6 concerted q-arm",
                  "C7 acro p",
                  "C14 PAR2",
                  "C15 PAR1"),
  abstract_label = c("4q-10q DUX4-containing homology",
                     "10p-18p homology",
                     "tightly linked 22q/21q/19q/1q/13q/17q clade",
                     "acrocentric short arms",
                     "Xq/Yq via PAR2",
                     "Xp/Yp via PAR1"),
  members = c("chr4_q, chr10_q",
              "chr10_p, chr18_p",
              "chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q",
              "chr13_p, chr14_p, chr15_p, chr21_p, chr22_p",
              "chrX_q, chrY_q",
              "chrX_p, chrY_p"),
  # Darkened versions of the slide 09 highlight rows, with C14/C15 kept in the
  # same PAR-green family but separated enough to distinguish PAR2 from PAR1.
  color = c("#C95F1B", "#A57900", "#2F6FB2", "#7B59C2", "#12884C", "#2FA36B"),
  slide09_highlight_family = c("#FDE2C8", "#FFF1AA", "#D6E8FF",
                               "#E5D8FA", "#CDEAD3", "#CDEAD3"),
  stringsAsFactors = FALSE
)

write.table(named_clades,
            file = file.path(out_dir, "clade_legend.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

clade_members <- setNames(strsplit(gsub(" ", "", named_clades$members), ","),
                          named_clades$community)
clade_short <- setNames(named_clades$short_label, named_clades$community)
clade_color <- setNames(named_clades$color, named_clades$community)
clade_abstract <- setNames(named_clades$abstract_label, named_clades$community)

message("Building NJ tree with ape::nj()")
set.seed(20260506)
tr_unrooted <- nj(as.dist(D))
write.tree(tr_unrooted, file = file.path(out_dir, "07b_unrooted_nj.newick"))

is_exact_clade <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  if (length(tips) < 2) return(NA)
  node <- tryCatch(getMRCA(tree, tips), error = function(e) NA_integer_)
  if (is.na(node)) return(FALSE)
  setequal(extract.clade(tree, node)$tip.label, tips)
}

acro_tips <- clade_members[["C7"]]
if (!isTRUE(is_exact_clade(tr_unrooted, acro_tips))) {
  stop("Cannot root on acrocentric p-arm MRCA because C7 is not monophyletic")
}
tr_rooted <- root(tr_unrooted,
                  node = getMRCA(tr_unrooted, acro_tips),
                  resolve.root = TRUE)
write.tree(tr_rooted, file = file.path(out_dir, "07b_rooted_acro.newick"))

clade_recovery <- data.frame(
  community = named_clades$community,
  short_label = named_clades$short_label,
  abstract_label = named_clades$abstract_label,
  members = named_clades$members,
  monophyletic_on_nj = vapply(clade_members,
                              function(tips) isTRUE(is_exact_clade(tr_unrooted, tips)),
                              logical(1)),
  stringsAsFactors = FALSE
)

write.table(clade_recovery,
            file = file.path(out_dir, "clade_recovery.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

matrix_audit <- data.frame(
  metric = c("rows",
             "columns",
             "present_arms",
             "missing_from_48_arm_set",
             "raw_diagonal_min",
             "raw_diagonal_max",
             "max_asymmetry_raw",
             "off_diagonal_min_after_zero_diag",
             "off_diagonal_max_after_zero_diag"),
  value = c(nrow(D_raw),
            ncol(D_raw),
            length(present_arms),
            paste(missing_arms, collapse = ", "),
            sprintf("%.6f", min(diag(D_raw))),
            sprintf("%.6f", max(diag(D_raw))),
            sprintf("%.6g", max(abs(D_raw - t(D_raw)))),
            sprintf("%.6f", min(D[upper.tri(D)])),
            sprintf("%.6f", max(D[upper.tri(D)]))),
  stringsAsFactors = FALSE
)

write.table(matrix_audit,
            file = file.path(out_dir, "matrix_audit.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

arm_short <- function(x) {
  gsub("_", "", sub("^chr", "", x))
}

tip_comm <- function(tips) {
  unname(comm_by_arm[tips])
}

tip_label <- function(tips) {
  paste0(arm_short(tips), " ", tip_comm(tips))
}

tip_col <- function(tips) {
  cc <- tip_comm(tips)
  ifelse(cc %in% named_clades$community, clade_color[cc], "#4A4A4A")
}

tip_face <- function(tips) {
  ifelse(tip_comm(tips) %in% named_clades$community, 2, 1)
}

legend_labels <- paste0(named_clades$community, " - ",
                        named_clades$abstract_label, "\n",
                        named_clades$members)

plot_rooted <- function() {
  plot_tr <- tr_rooted
  plot_tr$tip.label <- tip_label(tr_rooted$tip.label)

  op <- par(mar = c(0.45, 0.15, 1.05, 7.5), xpd = NA)
  on.exit(par(op), add = TRUE)

  plot(plot_tr,
       type = "phylogram",
       direction = "rightwards",
       use.edge.length = TRUE,
       show.tip.label = TRUE,
       tip.color = tip_col(tr_rooted$tip.label),
       font = tip_face(tr_rooted$tip.label),
       cex = 1.05,
       label.offset = 0.012,
       edge.width = 1.45,
       no.margin = FALSE)
  add.scale.bar(length = 0.1, cex = 0.75)
  title("Readable rooted NJ option: 41 PHR-positive subtelomeric arms",
        cex.main = 1.0,
        line = 0.2)
  mtext("Rooted at C7 acrocentric short-arm MRCA for display; branch lengths from the arm-level Jaccard distance matrix.",
        side = 3,
        line = -0.85,
        cex = 0.64,
        col = "#333333")

  for (community in named_clades$community) {
    node <- getMRCA(tr_rooted, clade_members[[community]])
    nodelabels(clade_short[[community]],
               node = node,
               frame = "rect",
               bg = adjustcolor(clade_color[[community]], alpha.f = 0.88),
               col = "white",
               cex = 0.58,
               adj = c(0.5, -0.25))
  }

}

plot_unrooted <- function() {
  plot_tr <- tr_unrooted
  plot_tr$tip.label <- tip_label(tr_unrooted$tip.label)

  op <- par(mar = c(5.0, 5.0, 3.0, 5.0), xpd = NA)
  on.exit(par(op), add = TRUE)

  plot(plot_tr,
       type = "unrooted",
       use.edge.length = TRUE,
       show.tip.label = TRUE,
       tip.color = tip_col(tr_unrooted$tip.label),
       font = tip_face(tr_unrooted$tip.label),
       cex = 0.64,
       label.offset = 0.006,
       edge.width = 1.3,
       lab4ut = "axial",
       no.margin = FALSE)
  title("Unrooted NJ option: no outgroup implied",
        cex.main = 1.0,
        line = 0.0)
  mtext("Same 41x41 arm-level Jaccard distance matrix; colored bold tips are the slide 09 / abstract clades.",
        side = 3,
        line = -0.9,
        cex = 0.62,
        col = "#333333")
}

plot_legend <- function() {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA)
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  y <- seq(0.86, 0.12, length.out = nrow(named_clades))
  text(0.02, 0.96, "Slide 09 / abstract clade map for slide 07b tree options",
       adj = 0, cex = 1.05, font = 2)
  for (i in seq_len(nrow(named_clades))) {
    rect(0.025, y[i] - 0.035, 0.055, y[i] + 0.035,
         col = named_clades$color[i], border = NA)
    text(0.075, y[i] + 0.014,
         paste0(named_clades$community[i], " - ",
                named_clades$abstract_label[i]),
         adj = 0, cex = 0.86, font = 2)
    text(0.075, y[i] - 0.022,
         named_clades$members[i],
         adj = 0, cex = 0.72, col = "#333333")
  }
}

rooted_pdf <- file.path(out_dir, "07b_rooted_acro_readable_large.pdf")
rooted_png <- file.path(out_dir, "07b_rooted_acro_readable_large.png")
unrooted_pdf <- file.path(out_dir, "07b_unrooted_nj_option.pdf")
unrooted_png <- file.path(out_dir, "07b_unrooted_nj_option.png")
legend_pdf <- file.path(out_dir, "07b_named_clade_legend.pdf")
legend_png <- file.path(out_dir, "07b_named_clade_legend.png")

pdf(rooted_pdf, width = 16, height = 9, useDingbats = FALSE)
plot_rooted()
dev.off()

png(rooted_png, width = 3840, height = 2160, res = 240, type = "cairo")
plot_rooted()
dev.off()

pdf(unrooted_pdf, width = 12, height = 12, useDingbats = FALSE)
plot_unrooted()
dev.off()

png(unrooted_png, width = 3200, height = 3200, res = 240, type = "cairo")
plot_unrooted()
dev.off()

pdf(legend_pdf, width = 9, height = 4.5, useDingbats = FALSE)
plot_legend()
dev.off()

png(legend_png, width = 2700, height = 1350, res = 240, type = "cairo")
plot_legend()
dev.off()

message("Wrote:")
message("  ", rooted_pdf)
message("  ", rooted_png)
message("  ", unrooted_pdf)
message("  ", unrooted_png)
message("  ", legend_pdf)
message("  ", legend_png)
