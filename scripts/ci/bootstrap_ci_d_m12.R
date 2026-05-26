#!/usr/bin/env Rscript
# D-M12: Bootstrap / Wilson / block-jackknife CIs for headline correlations.
#
# This is the R companion to scripts/ci/bootstrap_ci_d_m12.py. The Python
# version is the canonical implementation (numerics agree to 4 d.p.); this
# R script exists for downstream R-only consumers (paper_prep/figures/*.R)
# and as a cross-check. It re-implements the four CIs that have local data:
#   (a) Mantel rho bootstrap CI for HG002 Hi-C and CHM13 Hi-C
#   (a-F) Fisher z CI on the headline Mantel rho values (analytic stand-in)
#   (c) Wilson 95% CI for pedigree within-Leiden 494/538
#   (d) Hudson F_ST per superpopulation pair, block-jackknife over arms
#   (e-F) Fisher z CI on the mouse zygotene Spearman rho (stand-in only;
#         D-M5 arm-level Mantel CI is the canonical published CI).
# Deferred (no local data):
#   (b) Mantel trajectory under 5 exclusion sets x 5 resolutions.
#   (e) Per-pair bootstrap on the 344 mouse PHR pairs.
#
# Run:  Rscript scripts/ci/bootstrap_ci_d_m12.R

suppressPackageStartupMessages({
  # base only; no rare packages so this runs in any guix shell.
})

# ---- inputs -------------------------------------------------------------

ARM_DIST_TSV  <- "/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/arm_dist_matrix.tsv"
HG002_HIC_TSV <- "/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/hg002_contact_matrix.tsv"
CHM13_HIC_TSV <- "/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/hic_validation/res_50kb/chm13_contact_matrix.tsv"

PEDIGREE_NUM <- 494L
PEDIGREE_DEN <- 538L

SUPERPOPS <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# 10 arm/community pairs with significant Fisher p_adj < 0.05 (from
# end-to-end-report/report/04_heterogeneity.md "Population structure
# in cross-arm affinity"). Columns: comm, arm, cross+self per superpop.
FST_TABLE <- read.table(text = "
comm arm     cross_AFR cross_AMR cross_EAS cross_EUR cross_SAS self_AFR self_AMR self_EAS self_EUR self_SAS
C1   chr4_q   52 13 23 17 20  32 53 47 27 38
C3   chr16_q  60  5  2  1  6  54 76 76 61 64
C5   chr6_p   25 44 57 43 49  74 29 19 14 18
C15  chrX_p   71 62 52 43 48  18  0  0  0  0
C6   chr19_q  20 12  0  8 11  93 67 79 55 61
C3   chr9_q   45 17 19 18 12  66 62 58 43 59
C6   chr22_q  49 41 43 28 46  66 38 34 34 21
C15  chrY_p    5  1  0  1  0  13 15 20  9 18
C1   chr10_q   4  3  3  1 11  79 63 60 45 49
C11  chr6_q    1  1  0  4  4 107 78 79 55 64
", header = TRUE, stringsAsFactors = FALSE)

N_BOOT <- 10000L
N_PERM <- 10000L
RNG_SEED <- 20260518L

# ---- helpers ------------------------------------------------------------

load_matrix <- function(path) {
  raw <- read.table(path, header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE)
  M <- as.matrix(raw)
  arms_row <- rownames(M)
  arms_col <- colnames(M)
  stopifnot(setequal(arms_row, arms_col))
  M <- M[arms_row, arms_row]                      # square + row=col order
  M <- 0.5 * (M + t(M))                           # symmetrise round-off
  list(arms = arms_row, M = M)
}

upper_tri <- function(M) M[upper.tri(M)]

mantel_spearman <- function(sim, contact) {
  suppressWarnings(cor(upper_tri(sim), upper_tri(contact), method = "spearman"))
}

mantel_bootstrap_ci <- function(sim, contact, n_boot, rng_seed) {
  set.seed(rng_seed)
  n <- nrow(sim)
  out <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    s_b <- sim[idx, idx]
    c_b <- contact[idx, idx]
    out[b] <- suppressWarnings(
      cor(upper_tri(s_b), upper_tri(c_b), method = "spearman")
    )
  }
  list(lo = unname(quantile(out, 0.025, na.rm = TRUE)),
       hi = unname(quantile(out, 0.975, na.rm = TRUE)),
       dist = out)
}

mantel_permutation_p <- function(sim, contact, obs, n_perm, rng_seed) {
  set.seed(rng_seed + 1L)
  n <- nrow(sim)
  count <- 0L
  for (i in seq_len(n_perm)) {
    perm <- sample.int(n)
    sim_p <- sim[perm, perm]
    rho <- suppressWarnings(
      cor(upper_tri(sim_p), upper_tri(contact), method = "spearman")
    )
    if (!is.na(rho) && abs(rho) >= abs(obs)) count <- count + 1L
  }
  (count + 1) / (n_perm + 1)
}

wilson_ci <- function(k, n, alpha = 0.05) {
  if (n == 0) return(c(phat = NA_real_, lo = NA_real_, hi = NA_real_))
  p <- k / n
  z <- qnorm(1 - alpha / 2)
  denom <- 1 + z * z / n
  centre <- (p + z * z / (2 * n)) / denom
  half <- (z * sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
  c(phat = p, lo = centre - half, hi = centre + half)
}

hudson_fst <- function(ci, ni, cj, nj) {
  # Matches /moosefs/.../scripts/community/compute_fst_superpop.py:
  # F_ST = (HT - HS) / HT with HS = mean(2*pi*(1-pi), 2*pj*(1-pj)) and
  # HT = 2 * p_pool * (1 - p_pool), p_pool = (ci+cj)/(ni+nj).
  if (ni == 0 || nj == 0) return(NA_real_)
  pi <- ci / ni; pj <- cj / nj
  hs <- (2 * pi * (1 - pi) + 2 * pj * (1 - pj)) / 2
  p_pool <- (ci + cj) / (ni + nj)
  ht <- 2 * p_pool * (1 - p_pool)
  if (ht == 0) return(0)
  (ht - hs) / ht
}

block_jackknife_mean_ci <- function(values, alpha = 0.05) {
  values <- values[!is.na(values)]
  n <- length(values)
  if (n < 2) return(list(mean = NA_real_, lo = NA_real_, hi = NA_real_,
                         se = NA_real_, n = n))
  full <- mean(values)
  leave_one <- vapply(seq_len(n),
                      function(i) mean(values[-i]), numeric(1))
  pseudo <- n * full - (n - 1) * leave_one
  se <- sd(pseudo) / sqrt(n)
  tcrit <- qt(1 - alpha / 2, df = n - 1)
  list(mean = full, lo = full - tcrit * se, hi = full + tcrit * se,
       se = se, n = n)
}

fisher_z_ci <- function(rho, n, alpha = 0.05) {
  z <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  zcrit <- qnorm(1 - alpha / 2)
  c(rho = rho, lo = tanh(z - zcrit * se), hi = tanh(z + zcrit * se), n = n)
}

# ---- main ---------------------------------------------------------------

out_dir <- "scripts/ci"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# (c) Wilson CI for pedigree fraction --------------------------------------
ped <- wilson_ci(PEDIGREE_NUM, PEDIGREE_DEN)
cat(sprintf("\n[c] Pedigree within-Leiden: %d/%d = %.4f\n",
            PEDIGREE_NUM, PEDIGREE_DEN, ped["phat"]))
cat(sprintf("    Wilson 95%% CI: [%.4f, %.4f]  (width = %.4f)\n",
            ped["lo"], ped["hi"], ped["hi"] - ped["lo"]))

# (d) F_ST per superpop pair: block-jackknife -----------------------------
cat("\n[d] F_ST per superpopulation pair (Hudson estimator)\n")
cat("    Block-jackknife over 10 arm/community pairs (leave-one-out)\n\n")
fst_rows <- list()
for (i in seq_along(SUPERPOPS)) {
  if (i == length(SUPERPOPS)) next
  for (j in (i + 1):length(SUPERPOPS)) {
    p1 <- SUPERPOPS[i]; p2 <- SUPERPOPS[j]
    per_arm <- numeric(0)
    for (k in seq_len(nrow(FST_TABLE))) {
      row <- FST_TABLE[k, ]
      ci <- row[[paste0("cross_", p1)]]
      ni <- ci + row[[paste0("self_", p1)]]
      cj <- row[[paste0("cross_", p2)]]
      nj <- cj + row[[paste0("self_", p2)]]
      per_arm <- c(per_arm, hudson_fst(ci, ni, cj, nj))
    }
    bj <- block_jackknife_mean_ci(per_arm)
    cat(sprintf("    %s vs %s: F_ST = %+0.4f  95%% CI [%+0.4f, %+0.4f]  (n_arms=%d, se=%.4f)\n",
                p1, p2, bj$mean, bj$lo, bj$hi, bj$n, bj$se))
    fst_rows[[length(fst_rows) + 1L]] <- data.frame(
      pop1 = p1, pop2 = p2, mean_fst = bj$mean,
      ci_lo = bj$lo, ci_hi = bj$hi, se = bj$se, n_arms = bj$n,
      stringsAsFactors = FALSE)
  }
}
fst_tbl <- do.call(rbind, fst_rows)
write.table(fst_tbl,
            file.path(out_dir, "fst_block_jackknife_R.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# (a) Mantel bootstrap CI -------------------------------------------------
cat("\n[a] Mantel rho bootstrap CI (similarity vs Hi-C contact)\n")
cat(sprintf("    Bootstrap = arm-resampling with replacement, B = %d; perm p = %d.\n",
            N_BOOT, N_PERM))
cat("    DATA VERSION CAVEAT: local snapshot Feb 2026 != v5 rerun.\n\n")

ad <- load_matrix(ARM_DIST_TSV)
sim_full <- 1 - ad$M

reorder_mat <- function(M, src, target) {
  idx <- match(target, src)
  M[idx, idx]
}

mantel_out <- list()
for (entry in list(
      list(sample = "HG002", path = HG002_HIC_TSV),
      list(sample = "CHM13", path = CHM13_HIC_TSV))) {
  hc <- load_matrix(entry$path)
  common <- sort(intersect(ad$arms, hc$arms))
  sim <- reorder_mat(sim_full, ad$arms, common)
  contact <- reorder_mat(hc$M, hc$arms, common)
  obs <- mantel_spearman(sim, contact)
  boot <- mantel_bootstrap_ci(sim, contact, N_BOOT, RNG_SEED)
  perm <- mantel_permutation_p(sim, contact, obs, N_PERM, RNG_SEED)
  cat(sprintf("    %s: rho = %+0.4f  bootstrap 95%% CI [%+0.4f, %+0.4f]  (n_arms = %d, perm p = %.4e)\n",
              entry$sample, obs, boot$lo, boot$hi, length(common), perm))
  mantel_out[[length(mantel_out) + 1L]] <- data.frame(
    sample = entry$sample, n_arms = length(common),
    rho = obs, ci_lo = boot$lo, ci_hi = boot$hi, perm_p = perm,
    n_boot = N_BOOT, stringsAsFactors = FALSE)
}
mantel_tbl <- do.call(rbind, mantel_out)
write.table(mantel_tbl,
            file.path(out_dir, "mantel_bootstrap_ci_R.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# (a-Fisher) headline rho Fisher z stand-in -------------------------------
cat("\n[a-Fisher] Analytic Fisher z 95% CI for headline Mantel rho values\n")
headline <- data.frame(
  sample = c("CHM13",   "HG002",   "HG02559", "HG00658", "HG02148",
             "NA19036", "HG002",   "HG002"),
  tech   = c("Hi-C",    "Hi-C",    "Hi-C",    "Hi-C",    "Hi-C",
             "Hi-C",    "Pore-C",  "CiFi"),
  rho    = c(0.656,     0.657,     0.397,     0.276,     0.152,
             0.266,     0.486,     0.308),
  n_arms = c(38L,       41L,       37L,       37L,       37L,
             34L,       41L,       41L),
  stringsAsFactors = FALSE
)
for (k in seq_len(nrow(headline))) {
  fz <- fisher_z_ci(headline$rho[k], headline$n_arms[k])
  cat(sprintf("    %s %s: rho = %.3f  Fisher 95%% CI [%+0.3f, %+0.3f] (n_arms = %d)\n",
              headline$sample[k], headline$tech[k], fz["rho"], fz["lo"], fz["hi"],
              headline$n_arms[k]))
  headline$ci_lo[k] <- fz["lo"]
  headline$ci_hi[k] <- fz["hi"]
}
write.table(headline,
            file.path(out_dir, "mantel_fisher_z_ci_R.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# (e-Fisher) mouse Spearman stand-in --------------------------------------
cat("\n[e-Fisher] Mouse zygotene per-pair Spearman rho = 0.715  (n = 344)\n")
fz <- fisher_z_ci(0.715, 344)
cat(sprintf("    Fisher 95%% CI [%+0.3f, %+0.3f]\n", fz["lo"], fz["hi"]))
cat("    NOTE: pairs non-independent (D-M5). Use arm-level Mantel CI from D-M5 when published.\n")
