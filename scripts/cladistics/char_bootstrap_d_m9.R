#!/usr/bin/env Rscript
# D-M9: Character-level NJ bootstrap on PHR resampling.
#
# Closes the v5 deferred concern that "1000-replicate distance-matrix
# perturbation (Gaussian noise, sigma = 25% of off-diagonal IQR)" is a
# distance-matrix sensitivity analysis, NOT a phylogenetic bootstrap.
# This script does the actual character bootstrap: resample the 15,668
# PHRs (here 15,089 signal-bearing rows after filtering, see ANALYSIS_D_M9.md)
# with replacement, recompute the arm-level Jaccard distance matrix per
# replicate from the cached PHR-level cross-chromosome involvement, run
# ape::nj() (and an UPGMA tree) per replicate, summarise per-named-clade
# support over B replicates.
#
# Inputs:
#   --phr-tsv : per-PHR cross-chromosome involvement file
#               (default /home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/
#                        all-vs-all.p95.id95.len.tsv ; canonical path is
#                /moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv)
#   --ref-dist: 41/42x42 arm-level reference distance matrix to sanity-check
#               (default hic_validation/arm_dist_matrix.tsv)
#   --B       : number of bootstrap replicates (default 1000)
#   --seed    : RNG seed (default 20260518)
#   --out-dir : where to write per-replicate trees + summary
#
# Outputs:
#   <out-dir>/d_m9_support_summary.tsv
#       per-named-clade NJ + UPGMA bootstrap support over B replicates
#   <out-dir>/d_m9_jaccard_full.tsv
#       arm-level Jaccard distance recomputed on the FULL (non-resampled) PHR
#       set; sanity-check vs the saved reference matrix
#   <out-dir>/d_m9_meta.txt
#       B, runtime, named-clade tip lists, seed, paths
#   <out-dir>/d_m9_replicate_trees.RData (optional, only if --keep-trees)
#
# Run: Rscript scripts/cladistics/char_bootstrap_d_m9.R --B 1000

suppressPackageStartupMessages({
  library(ape)
})

# ---- args ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == paste0("--", name))
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
PHR_TSV  <- get_arg("phr-tsv",
                    "/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv")
REF_DIST <- get_arg("ref-dist",
                    "/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/hic_validation/arm_dist_matrix.tsv")
B        <- as.integer(get_arg("B", "1000"))
SEED     <- as.integer(get_arg("seed", "20260518"))
OUT_DIR  <- get_arg("out-dir", "/tmp/d_m9_bootstrap")
KEEP_TREES <- ("--keep-trees" %in% args)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[d-m9] PHR tsv     : %s\n", PHR_TSV))
cat(sprintf("[d-m9] Ref dist   : %s\n", REF_DIST))
cat(sprintf("[d-m9] B          : %d\n", B))
cat(sprintf("[d-m9] Seed       : %d\n", SEED))
cat(sprintf("[d-m9] Out dir    : %s\n", OUT_DIR))

# ---- load PHR table --------------------------------------------------------
t0_load <- Sys.time()
phr_df <- read.table(PHR_TSV, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "",
                     comment.char = "", check.names = FALSE)
cat(sprintf("[d-m9] PHR rows total                : %d\n", nrow(phr_df)))

# Keep only rows with non-empty cross-chromosome involvement.
# In the source file rows with no signal have chrs_involved == "." or empty.
phr_df <- phr_df[!(is.na(phr_df$chrs_involved) |
                     phr_df$chrs_involved %in% c("", ".")), ]
cat(sprintf("[d-m9] PHR rows after signal filter  : %d\n", nrow(phr_df)))

# Normalise arm tag: column "arm" is "parm"/"qarm", column "self_chr" e.g. chr10.
# Convert to "chr10_p" / "chr10_q" to match the saved matrix labels.
arm_pq <- ifelse(phr_df$arm == "parm", "p",
          ifelse(phr_df$arm == "qarm", "q", NA))
stopifnot(!any(is.na(arm_pq)))
phr_df$self_arm <- paste0(phr_df$self_chr, "_", arm_pq)
phr_df$pq <- arm_pq

# ---- build arm membership matrix ------------------------------------------
# For each PHR, the arms it contributes to are:
#   - its self arm (chrN_p or chrN_q)
#   - for each chromosome c in chrs_involved: arm c_<pq>, where <pq> is the
#     SAME p/q as the self arm. Subtelomeric homology is telomere-anchored
#     by construction (PHR flanks were extracted from telomere-anchored
#     contigs), so a P-arm flank's homology lands on P-arms of the involved
#     chromosomes and the same for Q-arms.
# This is the canonical interpretation used by the upstream arm-Jaccard
# pipeline; we sanity-check the full-data recomputation against the saved
# reference matrix below.
N_phr <- nrow(phr_df)

# Read reference matrix to get the canonical arm ordering.
ref_mat <- read.table(REF_DIST, header = TRUE, sep = "\t",
                      check.names = FALSE, row.names = 1)
ref_arms <- rownames(ref_mat)
stopifnot(all(ref_arms == colnames(ref_mat)))
N_arm <- length(ref_arms)
cat(sprintf("[d-m9] Reference arms              : %d\n", N_arm))

arm_index <- setNames(seq_along(ref_arms), ref_arms)

# Sparse-ish: instead of an N_phr x N_arm dense matrix (15k x 42 = ~5 MB doubles),
# store per-PHR a list of arm indices (1..N_arm) it contributes to.
# Then we build a per-PHR contribution table as a long-form data frame
# (row_id, arm_id) and tabulate into a count matrix per replicate.
make_phr_arm_long <- function() {
  out_phr <- integer(0)
  out_arm <- integer(0)
  for (i in seq_len(N_phr)) {
    pq <- phr_df$pq[i]
    self_a <- phr_df$self_arm[i]
    cs <- strsplit(phr_df$chrs_involved[i], ",", fixed = TRUE)[[1]]
    # Build candidate arm labels for the cross-chromosome list, with the
    # same p/q as the self arm.
    cand <- paste0(cs, "_", pq)
    cand <- unique(c(self_a, cand))
    cand <- cand[cand %in% ref_arms]
    out_phr <- c(out_phr, rep.int(i, length(cand)))
    out_arm <- c(out_arm, arm_index[cand])
  }
  list(phr = out_phr, arm = out_arm)
}
long_t0 <- Sys.time()
long <- make_phr_arm_long()
cat(sprintf("[d-m9] PHR-arm long pairs           : %d  (%.2fs)\n",
            length(long$phr),
            as.numeric(difftime(Sys.time(), long_t0, units = "secs"))))

# As a 0/1 sparse-like membership matrix (15k x 42, dense doubles ~5MB)
M <- matrix(0L, nrow = N_phr, ncol = N_arm,
            dimnames = list(NULL, ref_arms))
M[cbind(long$phr, long$arm)] <- 1L
arm_sizes_full <- colSums(M)
cat("[d-m9] Arm sizes (full, # of PHRs touching each arm):\n")
for (a in ref_arms) {
  cat(sprintf("  %-8s %5d\n", a, arm_sizes_full[a]))
}

# Per-PHR anchor arm index (length N_phr).
anchor_idx <- arm_index[phr_df$self_arm]
stopifnot(!any(is.na(anchor_idx)))

# Helper: arm-level distance matrix from per-PHR multiplicities w (length N_phr).
# Construction:
#   1. C[a, b] = sum_{i : anchor(i) = a} w_i * 1[b in fingerprint(i)]
#                = weighted count of PHRs anchored at arm a whose homology
#                  fingerprint (chrs_involved with same p/q tag) touches arm b.
#      Built as a sparse aggregation over the (anchor, fingerprint-arm) long
#      table. C is asymmetric.
#   2. Row-normalise to probability distributions: p_a[b] = C[a,b] / sum_c C[a,c].
#   3. Symmetric intersection-distance:
#         d(a, b) = 1 - sum_b min(p_a[b], p_b[b])
#      i.e. 1 - intersection of the two arms' cross-arm targeting
#      distributions; a proper metric on probability distributions and a
#      natural Jaccard-flavoured arm-arm distance that (a) recovers the
#      named clades on the full PHR set (see ANALYSIS_D_M9.md) and (b)
#      decomposes into per-PHR contributions, which is what makes the
#      character bootstrap (resample PHRs with replacement) meaningful.
arm_pair_dist <- function(w) {
  # Build C as N_arm x N_arm via sparse aggregation.
  # For each (phr_i, arm_j) entry in `long`, C[anchor_idx[phr_i], arm_j] += w[phr_i].
  C <- matrix(0, nrow = N_arm, ncol = N_arm,
              dimnames = list(ref_arms, ref_arms))
  contrib <- w[long$phr]
  # tapply group-sum is slow on long arrays; build a row-key directly.
  row_idx <- anchor_idx[long$phr]
  col_idx <- long$arm
  # Single linear key into a flat vector, then reshape.
  flat_key <- (col_idx - 1L) * N_arm + row_idx
  agg <- tabulate.weighted(flat_key, contrib, N_arm * N_arm)
  C[] <- agg
  rs <- rowSums(C)
  P  <- C / pmax(rs, 1)  # row-normalised; rows with rs=0 -> all 0 row
  # Intersection distance: d(a,b) = 1 - sum_i min(P[a,i], P[b,i]).
  # Vectorise: for each pair (a,b), sum of element-wise min over the columns.
  # Use a loop in C-like form; with N_arm = 42 this is 882 pairs × 42 cols.
  D <- matrix(0, N_arm, N_arm, dimnames = list(ref_arms, ref_arms))
  for (a in seq_len(N_arm - 1)) {
    for (b in (a + 1):N_arm) {
      inter <- sum(pmin(P[a, ], P[b, ]))
      D[a, b] <- D[b, a] <- 1 - inter
    }
  }
  diag(D) <- 0
  D
}

# tabulate weighted: equivalent to tapply(weights, key, sum) but faster.
tabulate.weighted <- function(key, weight, nbins) {
  out <- numeric(nbins)
  # rowsum on a 1-col matrix is fast enough.
  ord <- order(key)
  key_s <- key[ord]
  w_s <- weight[ord]
  rl <- rle(key_s)
  # cumulative sum of weights per run.
  idx_end <- cumsum(rl$lengths)
  idx_start <- c(1L, head(idx_end, -1L) + 1L)
  for (k in seq_along(rl$values)) {
    out[rl$values[k]] <- sum(w_s[idx_start[k]:idx_end[k]])
  }
  out
}

# Sanity-check: distance matrix on the full (non-resampled) PHR set.
w_full <- rep.int(1, N_phr)
D_full <- arm_pair_dist(w_full)
rownames(D_full) <- colnames(D_full) <- ref_arms

# Compare with saved reference. We expect tight agreement on which arm
# pairs are 'close' (subtelomeric homology) and which are 'far' (no shared
# sequence). Spearman correlation on the off-diagonal entries.
off_idx <- which(upper.tri(D_full), arr.ind = TRUE)
rho_off <- suppressWarnings(cor(
  as.matrix(ref_mat)[off_idx], D_full[off_idx],
  method = "spearman"))
cat(sprintf("[d-m9] Spearman(D_full, D_ref) off-diagonal: %.3f\n", rho_off))

# Write full-data Jaccard to out-dir.
write.table(round(D_full, 6), file = file.path(OUT_DIR, "d_m9_jaccard_full.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# ---- reference tree (full data) --------------------------------------------
set.seed(SEED)
tr_nj_ref <- nj(as.dist(D_full))
tr_up_ref <- as.phylo(hclust(as.dist(D_full), method = "average"))
cat("[d-m9] Built reference NJ + UPGMA trees on full PHR set\n")

# ---- named clade definitions ----------------------------------------------
clades <- list(
  PAR1     = c("chrX_p", "chrY_p"),
  PAR2     = c("chrX_q", "chrY_q"),
  ACRO_p   = c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p"),
  P10_18   = c("chr10_p", "chr18_p"),
  TIGHT_q  = c("chr22_q", "chr21_q", "chr19_q", "chr1_q", "chr13_q", "chr17_q"),
  DUX4     = c("chr4_q", "chr10_q")
)
# Filter clade tips to those actually present in the reference arm set.
clades <- lapply(clades, function(x) intersect(x, ref_arms))

# ---- clade-monophyly helpers ----------------------------------------------
# NJ tree is unrooted: clade "monophyly" = the tip set forms one side of a
# bipartition. Use prop.part-style check: for each internal edge, the descend-
# ant tip set vs the complement must match the target tip set (as a set).
# Implementation: for each internal node n in the (rerooted at tip 1) tree,
# get its descendants; check if either {desc} or {tips - desc} equals target.
get_all_bipartitions <- function(tr) {
  # Returns a list of sorted-tip-name vectors, one per internal edge.
  # Rerooting once at tip 1 makes descendants well-defined.
  tr2 <- root(tr, outgroup = tr$tip.label[1], resolve.root = TRUE)
  internal <- (Ntip(tr2) + 1):(Ntip(tr2) + Nnode(tr2))
  bp <- vector("list", length(internal))
  all_tips <- tr2$tip.label
  for (k in seq_along(internal)) {
    n <- internal[k]
    tips_n <- extract.clade(tr2, n)$tip.label
    bp[[k]] <- sort(tips_n)
  }
  bp
}

is_clade_in_unrooted <- function(tr, target_tips) {
  target_sorted <- sort(target_tips)
  all_tips <- tr$tip.label
  comp_sorted <- sort(setdiff(all_tips, target_sorted))
  bps <- get_all_bipartitions(tr)
  for (bp in bps) {
    if (identical(bp, target_sorted) || identical(bp, comp_sorted)) return(TRUE)
  }
  FALSE
}

is_clade_in_rooted <- function(tr, target_tips) {
  target_sorted <- sort(target_tips)
  if (length(target_tips) < 2) return(NA)
  mrca <- tryCatch(getMRCA(tr, target_tips), error = function(e) NA_integer_)
  if (is.na(mrca)) return(FALSE)
  desc <- sort(extract.clade(tr, mrca)$tip.label)
  identical(desc, target_sorted)
}

# Reference monophyly status on the full data.
nj_full_status <- sapply(clades, function(t) is_clade_in_unrooted(tr_nj_ref, t))
up_full_status <- sapply(clades, function(t) is_clade_in_rooted(tr_up_ref, t))
cat("[d-m9] Named-clade monophyly on FULL-data reference trees:\n")
for (cl in names(clades)) {
  cat(sprintf("  %-8s NJ=%s  UPGMA=%s  tips={%s}\n", cl,
              ifelse(nj_full_status[cl], "yes", "no"),
              ifelse(up_full_status[cl], "yes", "no"),
              paste(clades[[cl]], collapse = ",")))
}

# ---- bootstrap loop --------------------------------------------------------
set.seed(SEED)
boot_t0 <- Sys.time()

nj_support  <- setNames(integer(length(clades)), names(clades))
up_support  <- setNames(integer(length(clades)), names(clades))
valid_reps  <- 0L
rep_trees   <- if (KEEP_TREES) vector("list", B) else NULL

for (b in seq_len(B)) {
  draws <- sample.int(N_phr, N_phr, replace = TRUE)
  w     <- tabulate(draws, nbins = N_phr)
  Db    <- arm_pair_dist(w)
  rownames(Db) <- colnames(Db) <- ref_arms

  tr_nj_b <- tryCatch(nj(as.dist(Db)), error = function(e) NULL)
  tr_up_b <- tryCatch(as.phylo(hclust(as.dist(Db), method = "average")),
                      error = function(e) NULL)
  if (is.null(tr_nj_b) || is.null(tr_up_b)) next
  valid_reps <- valid_reps + 1L

  for (cl in names(clades)) {
    if (is_clade_in_unrooted(tr_nj_b, clades[[cl]])) {
      nj_support[cl] <- nj_support[cl] + 1L
    }
    if (isTRUE(is_clade_in_rooted(tr_up_b, clades[[cl]]))) {
      up_support[cl] <- up_support[cl] + 1L
    }
  }

  if (KEEP_TREES) rep_trees[[b]] <- list(nj = tr_nj_b, up = tr_up_b)

  if (b %% 100 == 0) {
    cat(sprintf("[d-m9] replicate %4d / %d   (%.1fs elapsed)\n",
                b, B,
                as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))))
  }
}

boot_secs <- as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))
cat(sprintf("[d-m9] Bootstrap done : %d valid / %d, %.1fs\n",
            valid_reps, B, boot_secs))

# ---- per-clade summary -----------------------------------------------------
nj_pct <- round(100 * nj_support / valid_reps, 1)
up_pct <- round(100 * up_support / valid_reps, 1)
v5_pct <- setNames(rep(100, length(clades)), names(clades))  # all 100% in v5

summary_df <- data.frame(
  clade        = names(clades),
  members      = sapply(clades, paste, collapse = ","),
  monophyletic_NJ_full   = unname(nj_full_status),
  monophyletic_UP_full   = unname(up_full_status),
  NJ_char_bootstrap_pct  = unname(nj_pct[names(clades)]),
  UP_char_bootstrap_pct  = unname(up_pct[names(clades)]),
  v5_sensitivity_pct     = unname(v5_pct[names(clades)]),
  stringsAsFactors = FALSE
)
write.table(summary_df, file = file.path(OUT_DIR, "d_m9_support_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- meta ------------------------------------------------------------------
load_secs <- as.numeric(difftime(Sys.time(), t0_load, units = "secs"))
writeLines(c(
  sprintf("d-m9 character bootstrap: B = %d (target 1000)", B),
  sprintf("  valid replicates       : %d", valid_reps),
  sprintf("  bootstrap runtime      : %.1f s", boot_secs),
  sprintf("  total runtime (incl I/O): %.1f s", load_secs),
  sprintf("  seed                   : %d", SEED),
  sprintf("  PHR table              : %s", PHR_TSV),
  sprintf("  PHR rows used (signal) : %d", N_phr),
  sprintf("  arms                   : %d (from %s)", N_arm, REF_DIST),
  sprintf("  Spearman(D_full,D_ref) : %.3f (off-diagonal)", rho_off),
  "",
  "Per-named-clade support (this analysis vs v5 sensitivity):"
), con = file.path(OUT_DIR, "d_m9_meta.txt"))

for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %-8s NJ=%5.1f%%  UPGMA=%5.1f%%   (v5 sensitivity=%.0f%%)\n",
              summary_df$clade[i],
              summary_df$NJ_char_bootstrap_pct[i],
              summary_df$UP_char_bootstrap_pct[i],
              summary_df$v5_sensitivity_pct[i]))
}

if (KEEP_TREES) {
  save(rep_trees, file = file.path(OUT_DIR, "d_m9_replicate_trees.RData"))
  cat(sprintf("[d-m9] Replicate trees saved to %s\n",
              file.path(OUT_DIR, "d_m9_replicate_trees.RData")))
}

cat("[d-m9] Done.\n")
