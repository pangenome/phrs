#!/usr/bin/env Rscript
# DEBUG / inference (not a manuscript figure).
# Depth-corrected within-pair / across-state bouquet test (Erik's design), swept
# across every flank window x Hi-C resolution and RANKED by signal strength.
#
# Within-pair share is normalization-irrelevant (any per-pair factor cancels), but
# it does NOT correct per-STAGE library depth.  Here each stage is divided by its
# own cross-pair MEAN raw contact before forming the share -- a CONSERVATIVE
# correction that removes all per-stage scaling, including the real uniform bouquet
# elevation, so the surviving zygotene enrichment is a LOWER BOUND.
#   share_zygo(pair) = norm_zygo / sum_s(norm_s),  norm_s = raw_s / mean_pairs(raw_s)
# Test mean_pairs(share_zygo) vs the even 1/4, PHR-node bootstrap (49 PHRs).
# Strength = z = (share - 0.25) / sd(bootstrap); ranked descending.
#
# INPUT: data/mouse_meiosis_sweep/seqlevel/<window>/mouse_<stage>_phr*_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): mouse_within_pair_depth.{tsv,png,pdf}
# Env: BOOT (4000).  Needs ggplot2.  Run: Rscript scripts/mouse/mouse_within_pair_depth.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
sweep_dir  <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir    <- Sys.getenv("OUT_DIR", "/tmp")
windows    <- c("1Mb", "2Mb", "4Mb")
res_list   <- c(5000, 10000, 20000, 50000, 100000)
stages     <- c("leptotene", "zygotene", "pachytene", "diplotene")
B          <- as.integer(Sys.getenv("BOOT", "4000"))
set.seed(1)

wmean <- function(x, w) sum(x * w) / sum(w)
ff <- function(w, s, r)
  list.files(file.path(sweep_dir, w),
             sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", s, r), full.names = TRUE)[1]

master <- function(win, res) {
  rd <- function(s) { f <- ff(win, s, res); if (is.na(f)) return(NULL)
    d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_raw), ]
    k <- paste(pmin(d$seq_a, d$seq_b), pmax(d$seq_a, d$seq_b), sep = "|")
    d <- d[!duplicated(k), ]; rownames(d) <- k[!duplicated(k)]; d }
  ds <- lapply(setNames(stages, stages), rd)
  if (any(vapply(ds, is.null, logical(1)))) return(NULL)
  ky <- Reduce(intersect, lapply(ds, rownames)); if (length(ky) < 10) return(NULL)
  base <- ds[[1]][ky, ]
  raw <- sapply(stages, function(s) ds[[s]][ky, "hic_contact_raw"])
  nodes <- sort(unique(c(base$seq_a, base$seq_b)))
  list(a = match(base$seq_a, nodes), b = match(base$seq_b, nodes),
       raw = raw, K = length(nodes), n = length(ky))
}

rows <- list()
for (win in windows) for (res in res_list) {
  M <- master(win, res); if (is.null(M)) next
  norm <- sweep(M$raw, 2, colMeans(M$raw), "/")          # conservative per-stage depth correction
  tot <- rowSums(norm); ok <- tot > 0
  norm <- norm[ok, ]; tot <- tot[ok]; a <- M$a[ok]; b <- M$b[ok]; n <- sum(ok)
  fz <- norm[, "zygotene"] / tot
  obs <- wmean(fz, rep(1, n))
  bt <- vapply(seq_len(B), function(.) {
    cnt <- tabulate(sample.int(M$K, M$K, replace = TRUE), nbins = M$K)
    w <- cnt[a] * cnt[b]; k <- w > 0; wmean(fz[k], w[k])
  }, numeric(1))
  c95 <- quantile(bt, c(.025, .975), na.rm = TRUE); sdb <- sd(bt, na.rm = TRUE)
  rows[[length(rows)+1]] <- data.frame(
    window = win, res = res, n_pairs = n,
    zygo_share = obs, lo = c95[[1]], hi = c95[[2]],
    z = (obs - 0.25) / sdb, p = 2 * min(mean(bt <= 0.25), mean(bt >= 0.25)),
    depth_zygo_vs_diplo = colSums(M$raw)["zygotene"] / colSums(M$raw)["diplotene"])
}
tab <- do.call(rbind, rows)
tab <- tab[order(-tab$z), ]                                # rank by strength
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(tab, file.path(out_dir, "mouse_within_pair_depth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

tab$lab <- factor(sprintf("%s / %g kb", tab$window, tab$res/1000),
                  levels = rev(sprintf("%s / %g kb", tab$window, tab$res/1000)))
tab$sig <- tab$lo > 0.25 | tab$hi < 0.25
g <- ggplot(tab, aes(zygo_share, lab, colour = window)) +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.8) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "share > 1/4 (95% CI)") +
  scale_colour_brewer(palette = "Set1", name = "flank window") +
  labs(x = "depth-corrected within-pair zygotene contact share (dashed = even 1/4; 95% PHR-node bootstrap CI)",
       y = NULL,
       title = "Depth-corrected bouquet strength across window x resolution (ranked)",
       subtitle = paste(
         "Within-pair zygotene contact share after dividing each stage by its cross-pair mean (conservative -> lower bound).",
         "Ranked top = strongest (highest z = (share-0.25)/SE).  Filled = CI excludes the even 1/4 share.",
         sep = "\n")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(out_dir, "mouse_within_pair_depth.png"), g, width = 9, height = 7, dpi = 150)
ggsave(file.path(out_dir, "mouse_within_pair_depth.pdf"), g, width = 9, height = 7)

cat("=== Depth-corrected within-pair zygotene share, ranked by strength (z) ===\n")
cat("(share>0.25 = bouquet; conservative depth-correction = lower bound;  * = CI excludes 0.25)\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("  %2d. %-4s %-6s n=%-4d share=%.3f CI[%.3f,%.3f] z=%5.1f p=%.4f  (depth z/d=%.2f)%s\n",
      i, window, paste0(res/1000,"kb"), n_pairs, zygo_share, lo, hi, z, p,
      depth_zygo_vs_diplo, ifelse(sig, "  *", ""))))
cat("\nwrote ", file.path(out_dir, "mouse_within_pair_depth.{png,pdf,tsv}"), "\n", sep = "")
