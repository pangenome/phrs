#!/usr/bin/env Rscript
# Neighbor-joining tree of HPRC v2 subtelomeric arm Jaccard distances.
#
# Input: 41x41 arm-level Jaccard distance matrix from
#   /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv
# Outputs (paper_prep/figures/nj_tree_arms/):
#   nj_tree.newick               - the unrooted NJ tree in Newick format
#   nj_tree_annotated.pdf/png    - rendered tree with abstract clades highlighted
#
# Annotated clades (from paper_prep/synthesis/ABSTRACT.md):
#   PAR1   : Xp/Yp
#   PAR2   : Xq/Yq
#   ACRO_p : acrocentric short arms (13p, 14p, 15p, 21p, 22p)
#   10p_18p: 10p-18p homology
#   TIGHT  : 22q, 21q, 19q, 1q, 13q, 17q tight clade
#   DUX4   : 4q-10q DUX4-containing homology

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ape)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
matrix_path <- if (length(args) >= 1) args[1] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
out_dir <- if (length(args) >= 2) args[2] else
  "paper_prep/figures/nj_tree_arms"
n_boot <- if (length(args) >= 3) as.integer(args[3]) else 1000

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load distance matrix --------------------------------------------------
mat_df <- read.table(matrix_path, header = TRUE, sep = "\t",
                     check.names = FALSE, row.names = 1)
D <- as.matrix(mat_df)
stopifnot(nrow(D) == ncol(D))
stopifnot(all(rownames(D) == colnames(D)))
# Symmetrize defensively (matrix should already be symmetric).
D <- (D + t(D)) / 2
diag(D) <- 0
cat(sprintf("Loaded %d x %d distance matrix\n", nrow(D), ncol(D)))

dist_obj <- as.dist(D)

# ---- neighbor-joining tree -------------------------------------------------
set.seed(20260506)
tr <- nj(dist_obj)
cat(sprintf("NJ tree: %d tips, %d internal nodes\n", Ntip(tr), Nnode(tr)))

# ---- abstract clade definitions -------------------------------------------
clades <- list(
  PAR1     = c("chrX_p", "chrY_p"),
  PAR2     = c("chrX_q", "chrY_q"),
  ACRO_p   = c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p"),
  P10_18   = c("chr10_p", "chr18_p"),
  TIGHT_q  = c("chr22_q", "chr21_q", "chr19_q", "chr1_q", "chr13_q", "chr17_q"),
  DUX4     = c("chr4_q", "chr10_q")
)
clade_labels <- c(
  PAR1    = "PAR1 (Xp/Yp)",
  PAR2    = "PAR2 (Xq/Yq)",
  ACRO_p  = "Acrocentric p-arms",
  P10_18  = "10p-18p",
  TIGHT_q = "Tight clade (22q/21q/19q/1q/13q/17q)",
  DUX4    = "DUX4 (4q/10q)"
)
clade_palette <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3", "#A65628"),
  names(clades)
)

# ---- root at acrocentric short-arm clade if monophyletic -------------------
acro_tips <- clades$ACRO_p
mrca_node <- tryCatch(getMRCA(tr, acro_tips), error = function(e) NA_integer_)
acro_monophyletic <- FALSE
if (!is.na(mrca_node)) {
  desc_tips <- tr$tip.label[ape::extract.clade(tr, mrca_node)$tip.label %in% tr$tip.label]
  desc_tips <- ape::extract.clade(tr, mrca_node)$tip.label
  acro_monophyletic <- setequal(desc_tips, acro_tips)
}
cat(sprintf("Acrocentric p-arms monophyletic on NJ tree: %s\n",
            ifelse(acro_monophyletic, "YES", "NO")))

if (acro_monophyletic) {
  tr_rooted <- root(tr, node = mrca_node, resolve.root = TRUE)
  rooting_note <- "rooted at MRCA of acrocentric short-arm clade"
} else {
  # Fall back to midpoint rooting via ape::root with outgroup heuristic.
  # If acrocentrics are not monophyletic, root on whichever single acrocentric
  # tip is most distal so the figure still has a stable orientation.
  tip_depths <- node.depth.edgelength(tr)[1:Ntip(tr)]
  names(tip_depths) <- tr$tip.label
  acro_present <- intersect(acro_tips, tr$tip.label)
  outgroup <- names(sort(tip_depths[acro_present], decreasing = TRUE))[1]
  tr_rooted <- root(tr, outgroup = outgroup, resolve.root = TRUE)
  rooting_note <- sprintf("rooted on outgroup %s (acrocentric p-arms not monophyletic)",
                          outgroup)
}
cat("Rooting: ", rooting_note, "\n", sep = "")

# ---- bootstrap support (perturbation-based) --------------------------------
# The input is a precomputed distance matrix derived from upstream alignment
# similarity, not a character matrix, so a Felsenstein column-bootstrap is not
# directly applicable. We use a perturbation-bootstrap: each replicate adds
# i.i.d. Gaussian noise to the off-diagonal distances (sigma scaled to a small
# fraction of the empirical IQR) and rebuilds the NJ tree. Reported support
# is the fraction of replicates in which each NJ clade in the reference tree
# is preserved.
boot_supports <- NULL
boot_done <- FALSE
boot_seconds <- NA_real_
boot_sigma <- NA_real_
if (n_boot > 0) {
  cat(sprintf("Running %d perturbation-bootstrap replicates...\n", n_boot))
  t0 <- Sys.time()
  off <- D[upper.tri(D)]
  # Scale sigma to a fraction of the spread of off-diagonal distances.
  # 25% of IQR gives a perturbation comparable to typical distance noise
  # in bootstrap-by-resampling on real character data.
  boot_sigma <- 0.25 * IQR(off)
  if (boot_sigma <= 0) boot_sigma <- 0.05
  cat(sprintf("Perturbation sigma = %.4f (25%% of off-diagonal IQR)\n",
              boot_sigma))
  n <- nrow(D)
  boot_trees <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    noise <- matrix(rnorm(n * n, sd = boot_sigma), n, n)
    noise <- (noise + t(noise)) / 2
    diag(noise) <- 0
    Db <- D + noise
    Db[Db < 0] <- 0
    diag(Db) <- 0
    rownames(Db) <- colnames(Db) <- rownames(D)
    bt <- tryCatch(nj(as.dist(Db)), error = function(e) NULL)
    boot_trees[[b]] <- bt
  }
  boot_trees <- Filter(function(x) inherits(x, "phylo"), boot_trees)
  n_valid <- length(boot_trees)
  if (n_valid > 0) {
    boot_supports <- prop.clades(tr_rooted, boot_trees, rooted = FALSE) /
      n_valid * 100
    boot_done <- TRUE
  }
  boot_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Bootstrap done (%d valid replicates, %.1f s)\n",
              n_valid, boot_seconds))
}

# ---- write Newick ----------------------------------------------------------
newick_path <- file.path(out_dir, "nj_tree.newick")
write.tree(tr_rooted, file = newick_path)
cat("Wrote ", newick_path, "\n", sep = "")

# ---- per-tip colours and clade membership ---------------------------------
tip_clade <- setNames(rep(NA_character_, Ntip(tr_rooted)), tr_rooted$tip.label)
for (cl in names(clades)) {
  hit <- intersect(clades[[cl]], names(tip_clade))
  tip_clade[hit] <- cl
}
tip_color <- ifelse(is.na(tip_clade), "#444444",
                    clade_palette[tip_clade])
tip_face  <- ifelse(is.na(tip_clade), 1, 2)  # bold for highlighted

# ---- recovery summary: is each abstract clade monophyletic? ---------------
recovery <- sapply(names(clades), function(cl) {
  tips <- intersect(clades[[cl]], tr_rooted$tip.label)
  if (length(tips) < 2) return(NA)
  mrca <- tryCatch(getMRCA(tr_rooted, tips), error = function(e) NA_integer_)
  if (is.na(mrca)) return(FALSE)
  desc <- ape::extract.clade(tr_rooted, mrca)$tip.label
  setequal(desc, tips)
})
cat("\nClade recovery (monophyletic on NJ tree?):\n")
for (cl in names(recovery)) {
  cat(sprintf("  %-8s %s  [%s]\n", cl,
              ifelse(isTRUE(recovery[[cl]]), "YES",
                     ifelse(isFALSE(recovery[[cl]]), "no", "n/a")),
              clade_labels[[cl]]))
}

# ---- plotting helper -------------------------------------------------------
plot_tree <- function(device_open, device_close) {
  device_open()
  op <- par(mar = c(2, 1, 4, 18), xpd = NA)
  plot(tr_rooted,
       type        = "phylogram",
       use.edge.length = TRUE,
       show.tip.label  = TRUE,
       tip.color    = tip_color[tr_rooted$tip.label],
       font         = tip_face[tr_rooted$tip.label],
       cex          = 0.85,
       label.offset = 0.005,
       edge.width   = 1.4,
       no.margin    = FALSE)
  add.scale.bar(length = 0.1, cex = 0.7)
  if (boot_done) {
    nodelabels(round(boot_supports), frame = "none", cex = 0.55,
               adj = c(1.1, -0.3), col = "#222222")
  }
  title(main = "NJ tree of HPRC v2 subtelomeric arm Jaccard distances",
        cex.main = 1.1)
  mtext(sprintf("%s%s",
                rooting_note,
                ifelse(boot_done,
                       sprintf("; bootstrap support at nodes (%d perturbation reps)",
                               n_valid),
                       "; bootstrap not computed")),
        side = 3, line = 0.2, cex = 0.75)
  legend("topright",
         inset      = c(-0.42, 0),
         legend     = clade_labels[names(clades)],
         fill       = clade_palette[names(clades)],
         border     = NA,
         bty        = "n",
         cex        = 0.78,
         title      = "Abstract clades",
         title.adj  = 0,
         xpd        = NA)
  par(op)
  device_close()
}

pdf_path <- file.path(out_dir, "nj_tree_annotated.pdf")
plot_tree(function() pdf(pdf_path, width = 11, height = 9), dev.off)
cat("Wrote ", pdf_path, "\n", sep = "")

png_path <- file.path(out_dir, "nj_tree_annotated.png")
plot_tree(function() png(png_path, width = 1600, height = 1300, res = 150),
          dev.off)
cat("Wrote ", png_path, "\n", sep = "")

# ---- side car summary for README -------------------------------------------
summary_path <- file.path(out_dir, ".nj_summary.tsv")
sumdf <- data.frame(
  clade        = names(clades),
  label        = clade_labels[names(clades)],
  monophyletic = unname(recovery),
  members      = sapply(clades, paste, collapse = ","),
  stringsAsFactors = FALSE
)
write.table(sumdf, file = summary_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(c(
  sprintf("rooting: %s", rooting_note),
  sprintf("bootstrap: %s",
          if (boot_done) sprintf("%d perturbation replicates (sigma=%.4f), %.1fs",
                                 n_valid,
                                 boot_sigma, boot_seconds)
          else sprintf("skipped (%d requested)", n_boot))
), file.path(out_dir, ".nj_meta.txt"))

# Per-clade min bootstrap support across the internal edges that subtend each
# named abstract clade. Saved alongside the recovery summary for the README.
if (boot_done) {
  internal_node_support <- setNames(boot_supports,
                                    seq.int(Ntip(tr_rooted) + 1,
                                            Ntip(tr_rooted) + Nnode(tr_rooted)))
  clade_support <- sapply(names(clades), function(cl) {
    tips <- intersect(clades[[cl]], tr_rooted$tip.label)
    if (length(tips) < 2) return(NA_real_)
    mrca <- tryCatch(getMRCA(tr_rooted, tips), error = function(e) NA_integer_)
    if (is.na(mrca)) return(NA_real_)
    as.numeric(internal_node_support[as.character(mrca)])
  })
  cat("\nBootstrap support at each named-clade MRCA:\n")
  for (cl in names(clade_support)) {
    cat(sprintf("  %-8s %s%%  [%s]\n", cl,
                ifelse(is.na(clade_support[[cl]]), "n/a",
                       sprintf("%.1f", clade_support[[cl]])),
                clade_labels[[cl]]))
  }
  sumdf$bootstrap_pct <- unname(clade_support)[match(sumdf$clade,
                                                     names(clade_support))]
  write.table(sumdf, file = summary_path, sep = "\t",
              quote = FALSE, row.names = FALSE)
}

cat("\nDone.\n")
