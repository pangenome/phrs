#!/usr/bin/env Rscript
# Step 2 of the proper bouquet test (run after compute_expected_trans.py).
#
# Forms the GENUINE per-PHR-pair trans-O/E and tests the MAGNITUDE hypothesis of
# Fig 4c (is the Jaccard <-> contact correlation strongest at zygotene?) at FINE
# resolution, with the PHR-node bootstrap (49 PHRs as the unit) and honest CIs.
#
#   O/E_mean(pair) = hic_contact_norm(pair) / expected_trans[chr_a, chr_b]
#
# (exact: the trans expected is constant per chromosome pair, so the mean over the
#  PHR x PHR block of observed/expected = mean-observed-density / expected).
# Reports, per window x resolution, the per-stage rho for the true O/E (and, for
# comparison, the raw density), and the zygotene gap = rho_zygo - mean(others)
# with a 95% PHR-node bootstrap CI + p.  A real bouquet = gap > 0 with CI above 0.
#
# INPUT:
#   data/mouse_meiosis_sweep/seqlevel/<window>/mouse_<stage>_phr*_<res>bp_seqlevel.tsv
#   <EXPECTED_DIR>/expected_trans_<stage>_<res>bp.tsv  (from compute_expected_trans.py)
# Output (OUT_DIR, default /tmp): mouse_true_oe_test.{tsv,png,pdf}
# Env: EXPECTED_DIR, BOOT (5000), SAT_FRAC (0.95). Needs ggplot2.
# Run: Rscript scripts/mouse/mouse_true_oe_test.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
sweep_dir  <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
exp_dir    <- Sys.getenv("EXPECTED_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/expected_trans"))
out_dir    <- Sys.getenv("OUT_DIR", "/tmp")
windows    <- c("1Mb", "2Mb", "4Mb")
res_list   <- c(5000, 10000, 20000, 50000)          # fine scales (+ 50 kb cross-check)
stages     <- c("leptotene", "zygotene", "pachytene", "diplotene")
B          <- as.integer(Sys.getenv("BOOT", "5000"))
set.seed(1)

wrho <- function(a, b, w) { idx <- rep.int(seq_along(a), w)
  if (length(idx) < 3) return(NA_real_); suppressWarnings(cor(rank(a[idx]), rank(b[idx]))) }
ci     <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
boot_p <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(NA); 2*min(mean(v<=0), mean(v>=0)) }
gap_of <- function(r) unname(r["zygotene"] - mean(r[setdiff(stages, "zygotene")]))

find_file <- function(win, stage, res) {
  pat <- sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", stage, res)
  list.files(file.path(sweep_dir, win), pattern = pat, full.names = TRUE)[1]
}

# symmetric chrom-pair -> expected map for one stage x resolution
load_expected <- function(stage, res) {
  fp <- file.path(exp_dir, sprintf("expected_trans_%s_%dbp.tsv", stage, res))
  if (!file.exists(fp)) stop("missing expected file: ", fp,
                             "\n  run compute_expected_trans.py first.")
  e <- read.delim(fp, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  k <- paste(pmin(e$region1, e$region2), pmax(e$region1, e$region2), sep = "|")
  setNames(e$expected, k)
}
expk <- function(ca, cb) paste(pmin(ca, cb), pmax(ca, cb), sep = "|")

# master: pairs present in all 4 stages, with density (norm) and true O/E per stage
build_master <- function(win, res) {
  emap <- lapply(setNames(stages, stages), load_expected, res = res)
  rd <- function(stage) {
    d <- read.delim(find_file(win, stage, res), sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
    k <- paste(pmin(d$seq_a, d$seq_b), pmax(d$seq_a, d$seq_b), sep = "||")
    d <- d[!duplicated(k), ]; rownames(d) <- k[!duplicated(k)]; d
  }
  ds   <- lapply(setNames(stages, stages), rd)
  keys <- Reduce(intersect, lapply(ds, rownames))
  base <- ds[[1]][keys, ]
  dens <- sapply(stages, function(s) ds[[s]][keys, "hic_contact_norm"]); colnames(dens) <- stages
  oe   <- sapply(stages, function(s) {
            ex <- emap[[s]][expk(base$chr_a, base$chr_b)]
            ds[[s]][keys, "hic_contact_norm"] / ex })
  colnames(oe) <- stages
  ok <- rowSums(is.na(oe)) == 0
  base <- base[ok, ]; dens <- dens[ok, , drop = FALSE]; oe <- oe[ok, , drop = FALSE]
  nodes <- sort(unique(c(base$seq_a, base$seq_b)))
  list(a = match(base$seq_a, nodes), b = match(base$seq_b, nodes),
       jac = base$jaccard, dens = dens, oe = oe, K = length(nodes), n = nrow(base))
}

rows <- list()
for (win in windows) for (res in res_list) {
  M <- tryCatch(build_master(win, res), error = function(e) { message(conditionMessage(e)); NULL })
  if (is.null(M) || M$n < 10) next
  w1 <- rep(1, M$n)
  for (metric in c("true-O/E", "density")) {
    Y <- if (metric == "true-O/E") M$oe else M$dens
    obsR <- sapply(stages, function(s) wrho(M$jac, Y[, s], w1))
    bt <- numeric(B)
    for (bb in seq_len(B)) {
      cnt <- tabulate(sample.int(M$K, M$K, replace = TRUE), nbins = M$K)
      w <- cnt[M$a] * cnt[M$b]; keep <- w > 0
      r <- sapply(stages, function(s) wrho(M$jac[keep], Y[keep, s], w[keep]))
      bt[bb] <- gap_of(r)
    }
    g <- gap_of(obsR); c95 <- ci(bt)
    rows[[length(rows)+1]] <- data.frame(window = win, res = res, metric = metric, n_pairs = M$n,
      rho_lepto = obsR["leptotene"], rho_zygo = obsR["zygotene"],
      rho_pachy = obsR["pachytene"], rho_diplo = obsR["diplotene"],
      zygo_gap = g, gap_lo = c95[[1]], gap_hi = c95[[2]], gap_p = boot_p(bt))
  }
}
if (!length(rows))
  stop("no results: expected_trans tables not found under ", exp_dir,
       "\n  run scripts/mouse/compute_expected_trans.py first (on a host with /moosefs + the coolers).")
tab <- do.call(rbind, rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(tab, file.path(out_dir, "mouse_true_oe_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tab$window <- factor(tab$window, levels = windows)
tab$reslab <- factor(paste0(tab$res/1000, " kb"), levels = paste0(res_list/1000, " kb"))
tab$sig    <- !is.na(tab$gap_lo) & (tab$gap_lo > 0 | tab$gap_hi < 0)
g <- ggplot(tab, aes(zygo_gap, reslab, colour = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = gap_lo, xmax = gap_hi), orientation = "y",
                width = 0.3, linewidth = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = sig), size = 2.8, position = position_dodge(width = 0.5)) +
  facet_wrap(~ window, ncol = 1) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "95% CI excludes 0") +
  scale_colour_manual(values = c(`true-O/E` = "#d62728", density = "#7f7f7f")) +
  labs(x = expression("zygotene gap  " ~ rho[zygo] - bar(rho)[others] ~ "  (95% PHR-node bootstrap CI)"),
       y = "Hi-C resolution",
       title = sprintf("Proper test: zygotene magnitude gap, true trans-O/E vs density (B=%d)", B),
       subtitle = "gap>0 with CI above 0 = a real bouquet peak.  true-O/E uses cooltools expected_trans; density = raw/hic_bins.") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(out_dir, "mouse_true_oe_test.png"), g, width = 9, height = 8, dpi = 150)
ggsave(file.path(out_dir, "mouse_true_oe_test.pdf"), g, width = 9, height = 8)

cat("=== Proper bouquet test: per-stage rho + zygotene gap, true trans-O/E vs density ===\n")
cat("    (zygo_gap>0 & CI above 0 = real bouquet peak; * = CI excludes 0)\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("  %-4s %-6s %-9s n=%-4d  rho L/Z/P/D = %+.3f/%+.3f/%+.3f/%+.3f  zygo_gap=%+.3f CI[%+.3f,%+.3f] p=%.4f%s\n",
      as.character(window), paste0(res/1000,"kb"), metric, n_pairs,
      rho_lepto, rho_zygo, rho_pachy, rho_diplo, zygo_gap, gap_lo, gap_hi, gap_p,
      ifelse(sig, "  *", ""))))
cat("\nwrote ", file.path(out_dir, "mouse_true_oe_test.{png,pdf,tsv}"), "\n", sep = "")
