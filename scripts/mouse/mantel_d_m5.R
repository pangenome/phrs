#!/usr/bin/env Rscript
# D-M5: Arm-level Mantel test on mouse PHR Jaccard distance vs Zuo 2021
# meiotic Hi-C contact, per stage. Replaces the per-arm-pair Spearman with
# vegan::mantel (10,000 row+column permutations) plus an arm-block bootstrap
# 95% CI on rho.
#
# Cross-validates the Python implementation (scripts/mouse/mantel_d_m5.py).
#
# Usage:
#   Rscript scripts/mouse/mantel_d_m5.R \
#     --arm-dist /moosefs/.../subtelo_1Mb/similarity/mouse.dist_matrix.tsv \
#     --hic-dir  /moosefs/.../community_analysis_1Mb/50000bp \
#     --out-tsv  paper_prep/synthesis/ANALYSIS_D_M5_results.tsv
#
# Inputs match the Python script's expectations.

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
})

option_list <- list(
  make_option("--arm-dist", type = "character",
              help = "mouse arm-level Jaccard distance matrix TSV"),
  make_option("--invert-arm-dist", action = "store_true", default = TRUE,
              help = "convert distance to similarity via sim = 1 - dist before Mantel (default TRUE; matches upstream pipeline so rho is positive)"),
  make_option("--no-invert-arm-dist", action = "store_false",
              dest = "invert-arm-dist",
              help = "keep raw distance matrix (yields a negative rho)"),
  make_option("--hic-dir", type = "character",
              help = "directory with zuo2021_<stage>_oe_matrix.tsv files"),
  make_option("--hic-pattern", type = "character",
              default = "zuo2021_{stage}_oe_matrix.tsv",
              help = "filename pattern with {stage} placeholder"),
  make_option("--stages", type = "character",
              default = "leptotene,zygotene,pachytene,diplotene",
              help = "comma-separated stages"),
  make_option("--n-perm", type = "integer", default = 10000L),
  make_option("--n-boot", type = "integer", default = 10000L),
  make_option("--seed",   type = "integer", default = 42L),
  make_option("--out-tsv", type = "character",
              help = "output per-stage TSV")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$`arm-dist`), !is.null(opt$`hic-dir`),
          !is.null(opt$`out-tsv`))

log_msg <- function(...) message("[mantel_d_m5.R] ", ...)

load_labeled <- function(path) {
  log_msg("load_labeled path=", path)
  d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  comment.char = "", stringsAsFactors = FALSE)
  rn <- d[[1]]
  mat <- as.matrix(d[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rn
  if (!identical(rn, colnames(mat))) {
    stop("row/col labels disagree in ", path)
  }
  log_msg("  shape=(", nrow(mat), ",", ncol(mat), ") labels=", length(rn))
  mat
}

align_shared <- function(a, b) {
  shared <- intersect(rownames(a), rownames(b))
  log_msg("align_shared |A|=", nrow(a), " |B|=", nrow(b),
          " shared=", length(shared))
  list(shared = shared, a = a[shared, shared, drop = FALSE],
       b = b[shared, shared, drop = FALSE])
}

mantel_rho_dist <- function(a, b) {
  # vegan::mantel takes 'dist' objects. We treat both NxN matrices as dist
  # matrices (off-diagonal upper triangle). Use Spearman to match the
  # upstream pipeline's reported rho.
  as.numeric(vegan::mantel(as.dist(a), as.dist(b), method = "spearman",
                           permutations = 0)$statistic)
}

# --- Main ----------------------------------------------------------------
arm <- load_labeled(opt$`arm-dist`)
if (isTRUE(opt$`invert-arm-dist`)) {
  log_msg("converting arm distance -> similarity via sim = 1 - dist")
  arm <- 1.0 - arm
}
stages <- strsplit(opt$stages, ",")[[1]]
n_perm <- opt$`n-perm`
n_boot <- opt$`n-boot`
seed   <- opt$seed

rows <- list()
for (stage in stages) {
  hic_path <- file.path(opt$`hic-dir`,
                        sub("\\{stage\\}", stage, opt$`hic-pattern`))
  log_msg("--- stage=", stage, " hic_path=", hic_path)
  if (!file.exists(hic_path)) {
    log_msg("  MISSING — skipping")
    rows[[length(rows) + 1L]] <- data.frame(
      stage = stage, n_arms = NA_integer_, mantel_rho = NA_real_,
      perm_p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      per_pair_spearman = NA_real_, per_pair_n = NA_integer_,
      stringsAsFactors = FALSE)
    next
  }
  hic <- load_labeled(hic_path)
  al <- align_shared(arm, hic)
  n_shared <- length(al$shared)
  log_msg("  n_shared_arms=", n_shared)

  set.seed(seed)
  m <- vegan::mantel(as.dist(al$a), as.dist(al$b),
                     method = "spearman", permutations = n_perm)
  rho_obs <- as.numeric(m$statistic)
  perm_p  <- as.numeric(m$signif)
  log_msg("  mantel_rho=", sprintf("%.4f", rho_obs),
          " perm_p=", sprintf("%.6f", perm_p))

  # Arm-block bootstrap
  set.seed(seed + 1L)
  boot <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n_shared, n_shared, replace = TRUE)
    sa <- al$a[idx, idx, drop = FALSE]
    sb <- al$b[idx, idx, drop = FALSE]
    boot[i] <- tryCatch(mantel_rho_dist(sa, sb),
                        error = function(e) NA_real_)
  }
  valid <- boot[is.finite(boot)]
  ci_lo <- as.numeric(quantile(valid, 0.025))
  ci_hi <- as.numeric(quantile(valid, 0.975))
  log_msg("  bootstrap CI 95%: (", sprintf("%.4f", ci_lo), ", ",
          sprintf("%.4f", ci_hi), ") n_valid=", length(valid),
          "/", n_boot)

  # Per-pair Spearman (the v5 reading)
  va <- al$a[upper.tri(al$a)]
  vb <- al$b[upper.tri(al$b)]
  keep <- is.finite(va) & is.finite(vb)
  per_rho <- as.numeric(stats::cor(va[keep], vb[keep], method = "spearman"))
  log_msg("  per-pair spearman rho=", sprintf("%.4f", per_rho),
          " (n=", sum(keep), " pairs)")

  rows[[length(rows) + 1L]] <- data.frame(
    stage = stage, n_arms = n_shared,
    mantel_rho = rho_obs, perm_p = perm_p,
    ci_lo = ci_lo, ci_hi = ci_hi,
    per_pair_spearman = per_rho, per_pair_n = sum(keep),
    stringsAsFactors = FALSE)
}

out <- do.call(rbind, rows)
dir.create(dirname(opt$`out-tsv`), recursive = TRUE, showWarnings = FALSE)
write.table(out, opt$`out-tsv`, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
log_msg("wrote ", opt$`out-tsv`, " (", nrow(out), " rows)")
print(out)
