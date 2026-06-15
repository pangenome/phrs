#!/usr/bin/env Rscript
# DEBUG / inference (not a manuscript figure).
# The WITHIN-PAIR / across-state design (Erik's framing): we do NOT compare PHRs
# to each other; we compare the SAME PHR pair across the four meiotic states.
# For a within-pair contrast, any per-pair-constant factor (length, hic_bins, a
# trans "expected") cancels exactly -- so normalization is irrelevant here, by
# construction.  Statistic per pair: the fraction of its contact at each state,
#   frac_s(pair) = contact_s / sum_s(contact)   (raw; identical under any per-pair norm).
#
# Two distinct questions, both with the PHR-node bootstrap (49 PHRs as the unit):
#  (1) BOUQUET (sequence-independent): do pairs put MORE than their even 1/4 share
#      of contact at zygotene?  test mean_pairs(frac_zygo) vs 0.25.
#  (2) CONCERTED (sequence-dependent): is that zygotene share predicted by
#      sequence similarity?  test rho(Jaccard, frac_zygo) vs 0.
#
# INPUT: data/mouse_meiosis_sweep/seqlevel/<window>/mouse_<stage>_phr*_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): mouse_within_pair.{tsv,png,pdf}
# Env: BOOT (4000).  Needs ggplot2.  Run: Rscript scripts/mouse/mouse_within_pair.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
sweep_dir  <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir    <- Sys.getenv("OUT_DIR", "/tmp")
windows    <- c("1Mb", "2Mb", "4Mb")
res_list   <- c(10000, 20000, 50000)
stages     <- c("leptotene", "zygotene", "pachytene", "diplotene")
B          <- as.integer(Sys.getenv("BOOT", "4000"))
set.seed(1)

ff <- function(win, st, res)
  list.files(file.path(sweep_dir, win),
             sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", st, res), full.names = TRUE)[1]
wmean <- function(x, w) sum(x * w) / sum(w)
wrho  <- function(a, b, w) { i <- rep.int(seq_along(a), w)
  if (length(i) < 3) return(NA_real_); suppressWarnings(cor(rank(a[i]), rank(b[i]))) }
ci    <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
p_vs  <- function(v, h0) { v <- v[!is.na(v)]; if (!length(v)) return(NA)
  2 * min(mean(v <= h0), mean(v >= h0)) }

master <- function(win, res) {
  rd <- function(st) { f <- ff(win, st, res); if (is.na(f)) return(NULL)
    d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_raw), ]
    k <- paste(pmin(d$seq_a, d$seq_b), pmax(d$seq_a, d$seq_b), sep = "|")
    d <- d[!duplicated(k), ]; rownames(d) <- k[!duplicated(k)]; d }
  ds <- lapply(setNames(stages, stages), rd)
  if (any(vapply(ds, is.null, logical(1)))) return(NULL)
  ky <- Reduce(intersect, lapply(ds, rownames)); if (length(ky) < 10) return(NULL)
  base <- ds[[1]][ky, ]
  raw  <- sapply(stages, function(s) ds[[s]][ky, "hic_contact_raw"])
  tot  <- rowSums(raw); ok <- tot > 0
  base <- base[ok, ]; raw <- raw[ok, , drop = FALSE]; tot <- tot[ok]
  nodes <- sort(unique(c(base$seq_a, base$seq_b)))
  list(a = match(base$seq_a, nodes), b = match(base$seq_b, nodes),
       jac = base$jaccard, frac = raw / tot, K = length(nodes), n = nrow(base))
}

rows <- list()
for (win in windows) for (res in res_list) {
  M <- master(win, res); if (is.null(M)) next
  w1 <- rep(1, M$n)
  # observed
  obs_share <- sapply(stages, function(s) wmean(M$frac[, s], w1))           # mean fraction per state
  obs_rho   <- sapply(stages, function(s) wrho(M$jac, M$frac[, s], w1))     # concerted
  # bootstrap
  bs <- array(NA_real_, c(B, length(stages), 2), dimnames = list(NULL, stages, c("share","rho")))
  for (bb in seq_len(B)) {
    cnt <- tabulate(sample.int(M$K, M$K, replace = TRUE), nbins = M$K)
    w <- cnt[M$a] * cnt[M$b]; k <- w > 0
    for (s in stages) {
      bs[bb, s, "share"] <- wmean(M$frac[k, s], w[k])
      bs[bb, s, "rho"]   <- wrho(M$jac[k], M$frac[k, s], w[k])
    }
  }
  for (s in stages) {
    cs <- ci(bs[, s, "share"]); cr <- ci(bs[, s, "rho"])
    rows[[length(rows)+1]] <- data.frame(window = win, res = res, stage = s, n_pairs = M$n,
      share = obs_share[[s]], share_lo = cs[[1]], share_hi = cs[[2]], share_p = p_vs(bs[, s, "share"], 0.25),
      rho = obs_rho[[s]], rho_lo = cr[[1]], rho_hi = cr[[2]], rho_p = p_vs(bs[, s, "rho"], 0))
  }
  cat(sprintf("  %s %dkb done (n=%d)\n", win, res/1000, M$n))
}
tab <- do.call(rbind, rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(tab, file.path(out_dir, "mouse_within_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tab$window <- factor(tab$window, levels = windows)
tab$reslab <- factor(paste0(tab$res/1000, " kb"), levels = paste0(res_list/1000, " kb"))
tab$stage  <- factor(tab$stage, levels = rev(stages))
tab$share_sig <- !is.na(tab$share_lo) & (tab$share_lo > 0.25 | tab$share_hi < 0.25)
pal <- c(leptotene = "#1b9e77", zygotene = "#d62728", pachytene = "#7570b3", diplotene = "#e6ab02")
g <- ggplot(tab, aes(share, stage, colour = stage)) +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = share_lo, xmax = share_hi), orientation = "y", width = 0.25, linewidth = 0.8) +
  geom_point(aes(shape = share_sig), size = 2.8) +
  facet_grid(window ~ reslab) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "differs from 1/4 (95% CI)") +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = "within-pair mean fraction of contact at the state (dashed = even 1/4 share; 95% PHR-node bootstrap CI)",
       y = NULL,
       title = "Within-pair / across-state bouquet test (normalization-irrelevant by design)",
       subtitle = paste(
         "Per PHR pair, share of its own contact at each meiotic state, averaged over pairs.  >1/4 = that state is enriched within the pair.",
         "This is the 'same PHR across states' comparison; length/expected cancel exactly.  Rows: flank window; cols: resolution.",
         sep = "\n")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(out_dir, "mouse_within_pair.png"), g, width = 11, height = 8, dpi = 150)
ggsave(file.path(out_dir, "mouse_within_pair.pdf"), g, width = 11, height = 8)

cat("\n=== WITHIN-PAIR across-state (normalization-irrelevant) ===\n")
cat("(1) BOUQUET: mean within-pair contact share per state; >0.25 = enriched.  * = CI excludes 0.25\n")
cat("(2) CONCERTED: rho(Jaccard, within-pair share); >0 = similar PHRs concentrate there.  # = CI excludes 0\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("  %-4s %-6s %-10s n=%-4d share=%.3f CI[%.3f,%.3f] p=%.3f%s  | rho=%+.3f CI[%+.3f,%+.3f] p=%.3f%s\n",
      as.character(window), paste0(res/1000,"kb"), as.character(stage), n_pairs,
      share, share_lo, share_hi, share_p, ifelse(share_lo>0.25|share_hi<0.25,"*"," "),
      rho, rho_lo, rho_hi, rho_p, ifelse(rho_lo>0|rho_hi<0,"#"," "))))
cat("\nwrote ", file.path(out_dir, "mouse_within_pair.{png,pdf,tsv}"), "\n", sep = "")
