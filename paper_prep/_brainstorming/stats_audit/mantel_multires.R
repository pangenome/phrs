#!/usr/bin/env Rscript
# Aggregate Mantel rho + B/W (within/between contact ratio) across
# 5 mcool resolutions × 8 Hi-C/Pore-C/CiFi datasets into a single SI table.
#
# Source: /moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/{res}bp/<sample>_global_test.tsv
# Each file has two rows: 'within_vs_between' (Mann-Whitney U + p) and 'mantel'
# (rho + p, n_between holds the # arms tested).
#
# Output: paper_prep/synthesis/stats_audit/mantel_multires_si_table.tsv

resolutions <- c(5000, 10000, 20000, 50000, 100000)
samples <- c("chm13", "hg002", "hg002_porec", "hg002_cifi",
             "hg00658", "hg02148", "hg02559", "na19036")
base <- "/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based"

rows <- list()
for (sample in samples) {
  for (res in resolutions) {
    fp <- file.path(base, paste0(format(res, scientific = FALSE), "bp"),
                    paste0(sample, "_global_test.tsv"))
    if (!file.exists(fp)) {
      warning("missing: ", fp)
      next
    }
    df <- read.table(fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    wb <- df[df$test == "within_vs_between", , drop = FALSE]
    mn <- df[df$test == "mantel", , drop = FALSE]
    bw_ratio <- if (nrow(wb)) wb$between_mean / wb$within_mean else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      sample        = sample,
      resolution_bp = res,
      within_mean   = if (nrow(wb)) wb$within_mean else NA_real_,
      between_mean  = if (nrow(wb)) wb$between_mean else NA_real_,
      bw_ratio      = bw_ratio,
      mw_U          = if (nrow(wb)) wb$U_statistic else NA_real_,
      mw_p          = if (nrow(wb)) wb$p_value     else NA_real_,
      n_within      = if (nrow(wb)) wb$n_within    else NA_integer_,
      n_between     = if (nrow(wb)) wb$n_between   else NA_integer_,
      mantel_rho    = if (nrow(mn)) mn$U_statistic else NA_real_,
      mantel_p      = if (nrow(mn)) mn$p_value     else NA_real_,
      n_arms        = if (nrow(mn)) mn$n_between   else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}

tbl <- do.call(rbind, rows)

# BH-FDR across the family of 40 within-vs-between Mann-Whitney tests
# (exclude exact-zero p-values from naive correction; treat Mantel p=0 as
# "<1e-4" per the 10,000-permutation lower bound the script uses).
tbl$mw_q_BH <- p.adjust(tbl$mw_p, method = "BH")
mantel_p_for_bh <- tbl$mantel_p
mantel_p_for_bh[mantel_p_for_bh == 0] <- 1e-4
tbl$mantel_p_floored <- mantel_p_for_bh
tbl$mantel_q_BH <- p.adjust(mantel_p_for_bh, method = "BH")

out <- "paper_prep/synthesis/stats_audit/mantel_multires_si_table.tsv"
write.table(tbl, out, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat(sprintf("wrote %s (%d rows)\n", out, nrow(tbl)))

# Wide summary: Mantel rho per (sample × resolution)
wide_rho <- reshape(
  tbl[, c("sample", "resolution_bp", "mantel_rho")],
  idvar = "sample", timevar = "resolution_bp", direction = "wide"
)
wide_rho <- wide_rho[match(samples, wide_rho$sample), ]
out_rho <- "paper_prep/synthesis/stats_audit/mantel_multires_rho_wide.tsv"
write.table(wide_rho, out_rho, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat(sprintf("wrote %s (%d × %d wide rho)\n",
            out_rho, nrow(wide_rho), ncol(wide_rho) - 1))

wide_bw <- reshape(
  tbl[, c("sample", "resolution_bp", "bw_ratio")],
  idvar = "sample", timevar = "resolution_bp", direction = "wide"
)
wide_bw <- wide_bw[match(samples, wide_bw$sample), ]
out_bw <- "paper_prep/synthesis/stats_audit/mantel_multires_bw_wide.tsv"
write.table(wide_bw, out_bw, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat(sprintf("wrote %s (%d × %d wide B/W)\n",
            out_bw, nrow(wide_bw), ncol(wide_bw) - 1))

cat("\nMantel rho across resolutions (rows = sample, cols = bp):\n")
print(wide_rho, row.names = FALSE)
cat("\nB/W ratio across resolutions:\n")
print(wide_bw, row.names = FALSE)
cat("\nN within-vs-between tests with q_BH < 0.05:",
    sum(tbl$mw_q_BH < 0.05, na.rm = TRUE), "/", sum(!is.na(tbl$mw_q_BH)), "\n")
cat("N Mantel tests with q_BH < 0.05:",
    sum(tbl$mantel_q_BH < 0.05, na.rm = TRUE), "/",
    sum(!is.na(tbl$mantel_q_BH)), "\n")
