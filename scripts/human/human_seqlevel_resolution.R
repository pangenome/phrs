#!/usr/bin/env Rscript
# Human single-sample sequence<->3D coupling across Hi-C/Pore-C/CiFi resolutions.
# Cross-PHR rho(Jaccard, contact density) per dataset x resolution, with the
# PHR-node bootstrap (PHRs, not pairs, as the unit) so the significance is valid.
# (The within-pair / across-state design does NOT apply to human: human is one
#  contact map per sample, no meiotic states to compare -- that test is mouse-only.)
# Confirms the coupling is resolution-invariant and solid at 50 kb for every panel.
#
# INPUT: data/human_seqlevel_sweep/human_<dataset>_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): human_seqlevel_resolution.{tsv,png,pdf}
# Env: BOOT (4000).  Needs ggplot2.  Run: Rscript scripts/human/human_seqlevel_resolution.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
data_dir   <- Sys.getenv("DATA_DIR", file.path(repo_root, "data/human_seqlevel_sweep"))
out_dir    <- Sys.getenv("OUT_DIR", "/tmp")
res_list   <- c(5000, 10000, 20000, 50000, 100000)
datasets   <- c("HG002_porec", "HG002_hic", "CHM13_hic", "HG002_cifi")
B          <- as.integer(Sys.getenv("BOOT", "4000"))
set.seed(1)

wrho <- function(a, b, w) { i <- rep.int(seq_along(a), w)
  if (length(i) < 3) return(NA_real_); suppressWarnings(cor(rank(a[i]), rank(b[i]))) }

rows <- list()
for (ds in datasets) for (res in res_list) {
  f <- file.path(data_dir, sprintf("human_%s_%dbp_seqlevel.tsv", ds, res))
  if (!file.exists(f)) next
  d <- read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_norm), ]
  nodes <- sort(unique(c(d$seq_a, d$seq_b)))
  a <- match(d$seq_a, nodes); b <- match(d$seq_b, nodes); K <- length(nodes)
  obs <- wrho(d$jaccard, d$hic_contact_norm, rep(1, nrow(d)))
  bt <- vapply(seq_len(B), function(.) {
    cnt <- tabulate(sample.int(K, K, replace = TRUE), nbins = K)
    w <- cnt[a] * cnt[b]; k <- w > 0
    wrho(d$jaccard[k], d$hic_contact_norm[k], w[k])
  }, numeric(1))
  c95 <- quantile(bt, c(.025, .975), na.rm = TRUE)
  rows[[length(rows)+1]] <- data.frame(dataset = ds, res = res, n = nrow(d),
    pct_zero = 100 * mean(d$hic_contact_norm <= 0),
    rho = obs, lo = c95[[1]], hi = c95[[2]],
    p = 2 * min(mean(bt <= 0, na.rm = TRUE), mean(bt >= 0, na.rm = TRUE)))
}
tab <- do.call(rbind, rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(tab, file.path(out_dir, "human_seqlevel_resolution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

tab$dataset <- factor(tab$dataset, levels = datasets)
tab$reslab  <- factor(paste0(tab$res/1000, " kb"), levels = paste0(res_list/1000, " kb"))
tab$sig     <- !is.na(tab$lo) & (tab$lo > 0 | tab$hi < 0)
g <- ggplot(tab, aes(rho, reslab, colour = dataset, group = dataset)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.8) +
  geom_point(aes(shape = sig), size = 3) +
  facet_wrap(~ dataset, ncol = 1, scales = "free_x") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "CI excludes 0") +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(x = expression(rho ~ "(PHR Jaccard vs 3D contact density; 95% PHR-node bootstrap CI)"),
       y = "resolution",
       title = "Human sequence<->3D coupling is resolution-invariant (single-sample)",
       subtitle = "Cross-PHR rho per dataset x resolution.  Solid at 50 kb for Pore-C / Hi-C / CHM13; CiFi weak (very sparse).") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "human_seqlevel_resolution.png"), g, width = 8, height = 9, dpi = 150)
ggsave(file.path(out_dir, "human_seqlevel_resolution.pdf"), g, width = 8, height = 9)

cat("=== Human cross-PHR rho(Jaccard, contact), PHR-node bootstrap ===\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("  %-11s %-6s n=%-4d zero=%2.0f%%  rho=%+.3f CI[%+.3f,%+.3f] p=%.4f%s\n",
      dataset, paste0(res/1000,"kb"), n, pct_zero, rho, lo, hi, p, ifelse(sig, "  *", ""))))
cat("\nwrote ", file.path(out_dir, "human_seqlevel_resolution.{png,pdf,tsv}"), "\n", sep = "")
