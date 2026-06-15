#!/usr/bin/env Rscript
# DEBUG / inference (not a manuscript figure).
# Length-FREE across-stage test, built to answer Erik's objection: PHR length is
# constant across meiotic stages, so length cannot by itself create an
# across-stage peak; and "normalize by length" only changes how pairs are
# weighted against each other, it does not change a single pair's stage profile.
#
# Statistic: for each inter-chrom PHR pair present in all 4 stages, the FRACTION
# of its contact at each stage,  frac_s(pair) = contact_s / sum_s'(contact_s').
# Because hic_bins (the region size / number of bin-pairs) is stage-constant, it
# cancels in this ratio:  raw_s/sum(raw) == norm_s/sum(norm).  So frac is
# identical under raw and under the per-bin "O/E" (= raw/hic_bins) normalization,
# and is fully length- and depth-free.  Then per stage we ask whether sequence
# similarity (Jaccard) predicts a pair's fractional contact at that stage:
#   rho_s = Spearman(Jaccard, frac_s)   over pairs.
# rho_zygo > 0 (and > other stages), significant, would mean similar PHRs
# preferentially gain contact at the zygotene bouquet -- a real, length-free
# concerted-bouquet signal.  Significance via the same PHR-node bootstrap
# (49 PHRs, not 1135 pairs, as the unit).
#
# INPUT: data/mouse_meiosis_sweep/seqlevel/<window>/mouse_<stage>_phr*_<res>bp_seqlevel.tsv
# Output (OUT_DIR, default /tmp): mouse_stage_enrichment.{tsv,png,pdf}
# Env: BOOT (5000), RES (50000), SAT_FRAC (0.95).  Needs ggplot2.
# Run: Rscript scripts/mouse/mouse_stage_enrichment.R

suppressPackageStartupMessages(library(ggplot2))

.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
sweep_dir  <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir    <- Sys.getenv("OUT_DIR", "/tmp")
windows    <- c("1Mb", "2Mb", "4Mb")
cap_bp     <- c("1Mb" = 1e6, "2Mb" = 2e6, "4Mb" = 4e6)
RES        <- as.integer(Sys.getenv("RES", "50000"))
B          <- as.integer(Sys.getenv("BOOT", "5000"))
sat_frac   <- as.numeric(Sys.getenv("SAT_FRAC", "0.95"))
stages     <- c("leptotene", "zygotene", "pachytene", "diplotene")
set.seed(1)

wrho <- function(a, b, w) {                               # weighted Spearman (node bootstrap)
  idx <- rep.int(seq_along(a), w)
  if (length(idx) < 3) return(NA_real_)
  suppressWarnings(cor(rank(a[idx]), rank(b[idx])))
}
ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
boot_p <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(NA); 2*min(mean(v<=0), mean(v>=0)) }

find_file <- function(win, stage) {
  pat <- sprintf("^mouse_%s_.*_%dbp_seqlevel\\.tsv$", stage, RES)
  list.files(file.path(sweep_dir, win), pattern = pat, full.names = TRUE)[1]
}
build_master <- function(win) {                            # pairs present in all 4 stages
  rd <- function(stage) {
    d <- read.delim(find_file(win, stage), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) & !is.na(d$hic_contact_raw), ]
    k <- paste(pmin(d$seq_a, d$seq_b), pmax(d$seq_a, d$seq_b), sep="||")
    d <- d[!duplicated(k), ]; rownames(d) <- k[!duplicated(k)]; d
  }
  ds <- lapply(setNames(stages, stages), rd)
  keys <- Reduce(intersect, lapply(ds, rownames))
  base <- ds[[1]][keys, ]
  raw  <- sapply(stages, function(s) ds[[s]][keys, "hic_contact_raw"]); colnames(raw) <- stages
  tot  <- rowSums(raw)
  keep <- tot > 0
  base <- base[keep, ]; raw <- raw[keep, , drop=FALSE]; tot <- tot[keep]
  frac <- raw / tot                                        # length- & depth-free stage profile
  nodes <- sort(unique(c(base$seq_a, base$seq_b)))
  list(a = match(base$seq_a, nodes), b = match(base$seq_b, nodes),
       jac = base$jaccard, size_a = base$size_a, size_b = base$size_b,
       frac = frac, K = length(nodes))
}

rows <- list()
for (win in windows) {
  M <- build_master(win); cap <- cap_bp[win]
  masks <- list(all    = rep(TRUE, length(M$jac)),
                nonsat = M$size_a <  sat_frac*cap & M$size_b <  sat_frac*cap,
                sat    = M$size_a >= sat_frac*cap & M$size_b >= sat_frac*cap)
  for (fi in names(masks)) {
    mk <- masks[[fi]]; w1 <- rep(1, length(M$jac))
    obs <- sapply(stages, function(s) wrho(M$jac[mk], M$frac[mk, s], w1[mk]))
    # node bootstrap
    bt <- matrix(NA_real_, B, length(stages), dimnames = list(NULL, stages))
    for (bb in seq_len(B)) {
      cnt <- tabulate(sample.int(M$K, M$K, replace=TRUE), nbins=M$K)
      w <- cnt[M$a]*cnt[M$b]; mkk <- mk & (w>0)
      bt[bb, ] <- sapply(stages, function(s) wrho(M$jac[mkk], M$frac[mkk, s], w[mkk]))
    }
    for (s in stages) {
      c95 <- ci(bt[, s])
      rows[[length(rows)+1]] <- data.frame(window=win, filter=fi, stage=s,
        n_pairs=sum(mk), rho=obs[[s]], lo=c95[[1]], hi=c95[[2]], p=boot_p(bt[, s]))
    }
  }
}
tab <- do.call(rbind, rows)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write.table(tab, file.path(out_dir, "mouse_stage_enrichment.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

tab$window <- factor(tab$window, levels=windows)
tab$filter <- factor(tab$filter, levels=c("all","nonsat","sat"))
tab$stage  <- factor(tab$stage, levels=rev(stages))
tab$sig    <- !is.na(tab$lo) & (tab$lo>0 | tab$hi<0)
g <- ggplot(tab, aes(rho, stage, colour=stage)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0.25, linewidth=0.8) +
  geom_point(aes(shape=sig), size=2.8) +
  facet_grid(window ~ filter) +
  scale_shape_manual(values=c(`FALSE`=1, `TRUE`=16), name="95% CI excludes 0") +
  scale_colour_manual(values=c(leptotene="#1b9e77", zygotene="#d62728",
                               pachytene="#7570b3", diplotene="#e6ab02"), guide="none") +
  labs(x = expression("rho( Jaccard , fraction of contact at that stage )  -- length-free; 95% PHR-node bootstrap CI"),
       y = NULL,
       title = sprintf("Length-FREE across-stage test: do similar PHRs gain contact at the bouquet? (res=%dkb, B=%d)", RES/1000, B),
       subtitle = paste(
         "Per pair, fraction of its contact at each stage (raw_s / sum_s); hic_bins cancels so this is identical under raw and per-bin 'O/E', and length/depth-free.",
         "rho>0 = sequence similarity predicts more contact at that stage.  Filled = CI excludes 0.  Cols: all/non-sat/sat PHRs.  Rows: flank window.",
         sep="\n")) +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), strip.text=element_text(face="bold"), legend.position="top")
ggsave(file.path(out_dir, "mouse_stage_enrichment.png"), g, width=11, height=8, dpi=150)
ggsave(file.path(out_dir, "mouse_stage_enrichment.pdf"), g, width=11, height=8)

cat("=== Length-free across-stage test: rho(Jaccard, fraction of contact at stage), per-pair, PHR-node bootstrap ===\n")
cat("    (rho>0 = similar PHRs gain contact at that stage; * = 95% CI excludes 0)\n")
for (i in seq_len(nrow(tab))) with(tab[i, ],
  cat(sprintf("  %-4s %-7s %-10s n=%-5d rho=%+.3f  CI[%+.3f, %+.3f]  p=%.4f%s\n",
      as.character(window), as.character(filter), as.character(stage), n_pairs, rho, lo, hi, p,
      ifelse(sig, "  *", ""))))
cat("\nwrote ", file.path(out_dir, "mouse_stage_enrichment.{png,pdf,tsv}"), "\n", sep="")
