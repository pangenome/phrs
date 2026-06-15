#!/usr/bin/env Rscript
# DEBUG / inference (not a manuscript figure).
# Structure-aware significance for the mouse meiotic stage trajectory, treating
# the ~49 PHRs (NOT the ~1135 pairs) as the independent units.  Answers:
#   (1) Is the zygotene "bouquet" a real peak, and does it flip between
#       length-saturating and non-saturating PHRs?
#   (2) Are the arm-level (Mantel) correlations significant at all?
#
# TEST 1 -- PHR-node bootstrap (Snijders-Borgatti node bootstrap for dyadic data).
#   Resample the K PHRs with replacement; a pair (i,j) enters the resample with
#   multiplicity cnt_i * cnt_j (weighted Spearman by expansion).  Per replicate,
#   on the SAME resampled pair-set, compute per-pair Spearman rho of Jaccard vs
#   O/E contact for each stage, then the gap  D = rho_zygo - mean(other 3).
#   Done for masks all / nonsat / sat, and the flip  D_sat - D_nonsat, all from
#   ONE node resample so the flip CI is paired.  Pairs are restricted to those
#   present in all 4 stages so the contrast is paired across stages.
#
# TEST 2 -- Mantel permutation on the arm-level O/E matrix, per stage (all PHRs).
#   Permute arm labels of the contact matrix; p = frac(perm rho >= observed).
#
# Primary line: per-PHR-pair O/E, 50 kb.  Runs windows 1/2/4 Mb.
# INPUT: data/mouse_meiosis_sweep/seqlevel/<window>/mouse_<stage>_phr*_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): mouse_significance.{tsv,png,pdf}
# Env: BOOT (default 5000), PERM (default 5000), RES (default 50000), SAT_FRAC (0.95).
# Needs ggplot2.  Run: Rscript scripts/mouse/mouse_significance.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))

sweep_dir <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir   <- Sys.getenv("OUT_DIR", "/tmp")
windows   <- c("1Mb", "2Mb", "4Mb")
cap_bp    <- c("1Mb" = 1e6, "2Mb" = 2e6, "4Mb" = 4e6)
RES       <- as.integer(Sys.getenv("RES", "50000"))
B         <- as.integer(Sys.getenv("BOOT", "5000"))
P         <- as.integer(Sys.getenv("PERM", "5000"))
sat_frac  <- as.numeric(Sys.getenv("SAT_FRAC", "0.95"))
stages    <- c("leptotene", "zygotene", "pachytene", "diplotene")
set.seed(1)

sprho <- function(a, b) suppressWarnings(cor(rank(a), rank(b)))
# weighted Spearman via integer-multiplicity expansion (node bootstrap)
wrho <- function(a, b, w) {
  idx <- rep.int(seq_along(a), w)
  if (length(idx) < 3) return(NA_real_)
  suppressWarnings(cor(rank(a[idx]), rank(b[idx])))
}

find_file <- function(win, stage) {
  pat <- sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", stage, RES)
  f <- list.files(file.path(sweep_dir, win), pattern = pat, full.names = TRUE)
  if (length(f) != 1) stop("missing ", pat, " in ", win)
  f
}

# master table for a window: one row per inter-chrom pair present in ALL 4 stages,
# with jaccard, size_a/b, and O/E contact per stage; integer node ids a/b.
build_master <- function(win) {
  rd <- function(stage) {
    d <- read.delim(find_file(win, stage), sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
    k <- paste(pmin(d$seq_a, d$seq_b), pmax(d$seq_a, d$seq_b), sep = "||")
    d <- d[!duplicated(k), ]; rownames(d) <- k[!duplicated(k)]
    d
  }
  ds <- lapply(setNames(stages, stages), rd)
  keys <- Reduce(intersect, lapply(ds, rownames))
  base <- ds[[1]][keys, ]
  oe <- sapply(stages, function(s) ds[[s]][keys, "hic_contact_norm"])
  colnames(oe) <- stages
  nodes <- sort(unique(c(base$seq_a, base$seq_b)))
  list(a = match(base$seq_a, nodes), b = match(base$seq_b, nodes),
       jac = base$jaccard, size_a = base$size_a, size_b = base$size_b,
       oe = oe, K = length(nodes), seq_a = base$seq_a, seq_b = base$seq_b,
       arm_a = base$arm_a, arm_b = base$arm_b)
}

stage_rhos <- function(jac, oe, w, mask)                  # rho per stage on a mask
  setNames(vapply(stages, function(s) wrho(jac[mask], oe[mask, s], w[mask]),
                  numeric(1)), stages)
# per-stage peak/trough gap: gap_s = rho_s - mean(rho of the other 3 stages)
gaps_from_rhos <- function(r)
  setNames(vapply(stages, function(s) unname(r[s]) - mean(r[setdiff(stages, s)]),
                  numeric(1)), stages)

ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
boot_p <- function(v) {                                    # 2-sided bootstrap p vs 0
  v <- v[!is.na(v)]; if (!length(v)) return(NA)
  2 * min(mean(v <= 0), mean(v >= 0))
}

## ---- Mantel (arm-level O/E) ----
mantel_rho <- function(J, C) {
  m <- upper.tri(J) & !is.na(J) & !is.na(C)
  if (sum(m) < 3) return(NA_real_)
  sprho(J[m], C[m])
}
mantel_test <- function(M, stage) {
  agg <- aggregate(cbind(jac, oe) ~ ka + kb,
                   data = data.frame(ka = pmin(M$arm_a, M$arm_b),
                                     kb = pmax(M$arm_a, M$arm_b),
                                     jac = M$jac, oe = M$oe[, stage]), FUN = mean)
  arms <- sort(unique(c(agg$ka, agg$kb)))
  n <- length(arms)
  J <- matrix(NA_real_, n, n, dimnames = list(arms, arms)); C <- J
  ia <- match(agg$ka, arms); ib <- match(agg$kb, arms)
  for (r in seq_len(nrow(agg))) {
    J[ia[r], ib[r]] <- J[ib[r], ia[r]] <- agg$jac[r]
    C[ia[r], ib[r]] <- C[ib[r], ia[r]] <- agg$oe[r]
  }
  obs <- mantel_rho(J, C)
  perm <- replicate(P, { p <- sample(n); mantel_rho(J, C[p, p]) })
  list(obs = obs, p = (1 + sum(perm >= obs, na.rm = TRUE)) / (P + 1),
       n_arm = n, n_armpair = sum(upper.tri(J) & !is.na(J) & !is.na(C)))
}

## ---- run ----
filt_lv  <- c("all", "nonsat", "sat")
boot_rows <- list(); mant_rows <- list()
for (win in windows) {
  cat("== window", win, "(res", RES, ") ==\n")
  M <- build_master(win)
  cap <- cap_bp[win]
  w1  <- rep(1, length(M$jac))
  masks <- list(all    = rep(TRUE, length(M$jac)),
                nonsat = M$size_a <  sat_frac * cap & M$size_b <  sat_frac * cap,
                sat    = M$size_a >= sat_frac * cap & M$size_b >= sat_frac * cap)
  np <- sapply(masks, sum)

  # observed: per-stage rho and per-stage gap, for each mask
  obsR <- sapply(masks, function(mk) stage_rhos(M$jac, M$oe, w1, mk))      # stage x mask
  obsG <- sapply(masks, function(mk) gaps_from_rhos(stage_rhos(M$jac, M$oe, w1, mk)))

  # node bootstrap (one resample drives all masks -> paired flip), all stages
  bootG <- array(NA_real_, c(B, length(masks), length(stages)),
                 dimnames = list(NULL, names(masks), stages))
  for (bb in seq_len(B)) {
    cnt  <- tabulate(sample.int(M$K, M$K, replace = TRUE), nbins = M$K)
    w    <- cnt[M$a] * cnt[M$b]
    keep <- w > 0
    for (fi in seq_along(masks))
      bootG[bb, fi, ] <- gaps_from_rhos(stage_rhos(M$jac, M$oe, w, masks[[fi]] & keep))
  }

  for (fi in filt_lv) for (s in stages) {
    c95 <- ci(bootG[, fi, s])
    boot_rows[[length(boot_rows)+1]] <- data.frame(
      window = win, filter = fi, stage = s, n_pairs = np[[fi]],
      rho = obsR[s, fi], gap = obsG[s, fi],
      gap_lo = c95[[1]], gap_hi = c95[[2]], gap_p = boot_p(bootG[, fi, s]))
  }
  # length-confound flip per stage (paired): sat gap - nonsat gap
  for (s in stages) {
    flip <- bootG[, "sat", s] - bootG[, "nonsat", s]
    fci  <- ci(flip)
    boot_rows[[length(boot_rows)+1]] <- data.frame(
      window = win, filter = "sat-minus-nonsat", stage = s, n_pairs = NA,
      rho = NA, gap = obsG[s, "sat"] - obsG[s, "nonsat"],
      gap_lo = fci[[1]], gap_hi = fci[[2]], gap_p = boot_p(flip))
  }

  for (s in stages) {
    mt <- mantel_test(M, s)
    mant_rows[[length(mant_rows)+1]] <- data.frame(
      window = win, stage = s, mantel_rho = mt$obs, mantel_p = mt$p,
      n_arm = mt$n_arm, n_armpair = mt$n_armpair)
  }
}
boot <- do.call(rbind, boot_rows); mant <- do.call(rbind, mant_rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(boot, file.path(out_dir, "mouse_significance_bootstrap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mant, file.path(out_dir, "mouse_significance_mantel.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- forest plot: per-stage peak/trough gap, all stages ----
fp <- boot
fp$filter <- factor(fp$filter, levels = c("all","nonsat","sat","sat-minus-nonsat"))
fp$window <- factor(fp$window, levels = windows)
fp$stage  <- factor(fp$stage, levels = rev(stages))
fp$sig    <- !is.na(fp$gap_lo) & (fp$gap_lo > 0 | fp$gap_hi < 0)
g <- ggplot(fp, aes(gap, stage, colour = stage)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = gap_lo, xmax = gap_hi), orientation = "y",
                width = 0.25, linewidth = 0.8) +
  geom_point(aes(shape = sig), size = 2.8) +
  facet_grid(window ~ filter) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     name = "95% CI excludes 0") +
  scale_colour_manual(values = c(leptotene = "#1b9e77", zygotene = "#d62728",
                                 pachytene = "#7570b3", diplotene = "#e6ab02"),
                      guide = "none") +
  labs(x = expression("stage gap  " ~ rho[stage] - bar(rho)[others] ~
                       "  (per-pair O/E; 95% PHR-node bootstrap CI)"),
       y = NULL,
       title = sprintf("Per-stage peak/trough across the meiotic cycle -- PHR-node bootstrap (B=%d, res=%dkb)", B, RES/1000),
       subtitle = paste(
         "gap>0 = that stage is a peak, gap<0 = a trough, relative to the other three.  Filled point = CI excludes 0.",
         "Cols: all / non-saturating / saturating PHRs + the paired length-confound flip (sat-minus-nonsat).  Rows: flank window.",
         sep = "\n")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(out_dir, "mouse_significance.png"), g, width = 11, height = 8, dpi = 150)
ggsave(file.path(out_dir, "mouse_significance.pdf"), g, width = 11, height = 8)

## ---- console ----
cat("\n=== TEST 1: PHR-node bootstrap, per-stage gap = rho_stage - mean(other 3), per-pair O/E ===\n")
cat("    (gap>0 peak, gap<0 trough;  * = 95% CI excludes 0)\n")
for (i in seq_len(nrow(boot))) with(boot[i, ],
  cat(sprintf("  %-4s %-17s %-10s n=%-5s gap=%+.3f  CI[%+.3f, %+.3f]  p=%.4f%s\n",
      window, filter, stage, ifelse(is.na(n_pairs),"-",n_pairs), gap, gap_lo, gap_hi, gap_p,
      ifelse(!is.na(gap_lo) && (gap_lo > 0 | gap_hi < 0), "  *", ""))))
cat("\n=== TEST 2: Mantel permutation, arm-level O/E per stage (all PHRs) ===\n")
for (i in seq_len(nrow(mant))) with(mant[i, ],
  cat(sprintf("  %-4s %-10s mantel_rho=%+.3f  p=%.4f  (arms=%d, arm-pairs=%d)%s\n",
      window, stage, mantel_rho, mantel_p, n_arm, n_armpair,
      ifelse(mantel_p < 0.05, "  *", ""))))
cat("\nwrote ", file.path(out_dir, "mouse_significance{,_bootstrap,_mantel}.{png,pdf,tsv}"), "\n", sep = "")
