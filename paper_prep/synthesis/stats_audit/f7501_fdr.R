#!/usr/bin/env Rscript
# f7501 per-arm × per-superpopulation Fisher's exact tests with
# BH-FDR correction and exact 95% confidence intervals on the OR.
#
# Inputs:  /moosefs/guarracino/HPRCv2/L78442/per_flank.1Mb/f7501_per_arm_summary.tsv
# Outputs: paper_prep/synthesis/stats_audit/f7501_per_arm_per_superpop_fisher.tsv
#          paper_prep/synthesis/stats_audit/f7501_per_arm_summary_with_q.tsv

suppressPackageStartupMessages({})

src <- "/moosefs/guarracino/HPRCv2/L78442/per_flank.1Mb/f7501_per_arm_summary.tsv"
out_long <- "paper_prep/synthesis/stats_audit/f7501_per_arm_per_superpop_fisher.tsv"
out_wide <- "paper_prep/synthesis/stats_audit/f7501_per_arm_summary_with_q.tsv"

# Superpopulation haplotype counts (SURVEY_01 §1.9 / §4): AFR=134, AMR=88,
# EAS=104, EUR=65, SAS=74; sum = 465 (1 missing CHM13 haplotype).
sp_n <- c(AFR = 134L, AMR = 88L, EAS = 104L, EUR = 65L, SAS = 74L)
sp_total <- sum(sp_n)  # 465

df <- read.table(src, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(df))) {
  arm <- df$arm[i]
  community <- df$leiden_community[i]
  n_carriers <- df$n_haplotypes[i]
  for (sp in names(sp_n)) {
    a <- df[[sp]][i]                      # carriers in sp
    b <- sp_n[[sp]] - a                   # non-carriers in sp
    c <- n_carriers - a                   # carriers in other sps
    d <- (sp_total - sp_n[[sp]]) - c      # non-carriers in other sps
    if (a < 0 || b < 0 || c < 0 || d < 0) {
      stop(sprintf("Negative cell at %s × %s: a=%d b=%d c=%d d=%d",
                   arm, sp, a, b, c, d))
    }
    m <- matrix(c(a, b, c, d), nrow = 2, byrow = FALSE,
                dimnames = list(c("carrier", "non_carrier"),
                                c(sp, paste0("not_", sp))))
    ft1 <- fisher.test(m, alternative = "greater")
    ft2 <- fisher.test(m)  # two-sided -> 95% OR CI
    # naive Wald OR for sanity (matches table's `best_enriched_OR` more closely)
    naive_or <- if (b == 0 || c == 0) Inf else (a * d) / (b * c)
    rows[[length(rows) + 1L]] <- data.frame(
      arm                = arm,
      leiden_community   = community,
      superpop           = sp,
      carriers_in_sp     = a,
      noncarriers_in_sp  = b,
      carriers_other_sp  = c,
      noncarriers_other  = d,
      n_carriers_total   = n_carriers,
      pct_carriers_in_sp = ifelse(sp_n[[sp]] > 0, 100 * a / sp_n[[sp]], NA_real_),
      naive_OR           = naive_or,
      conditional_MLE_OR = unname(ft2$estimate),
      OR_CI95_low        = ft2$conf.int[1],
      OR_CI95_high       = ft2$conf.int[2],
      p_one_sided        = ft1$p.value,
      p_two_sided        = ft2$p.value,
      stringsAsFactors   = FALSE
    )
  }
}

long <- do.call(rbind, rows)

# BH-FDR across the full per-arm × per-superpop family using the one-sided
# (greater) p-values, which is the test reported in SURVEY_01 §1.9.
long$q_BH <- p.adjust(long$p_one_sided, method = "BH")

# Sort by q for readability
long <- long[order(long$q_BH, long$p_one_sided), ]

write.table(long, out_long, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat(sprintf("wrote %s (%d rows; %d arms × %d superpops)\n",
            out_long, nrow(long), nrow(df), length(sp_n)))

# Also produce a "wide" file matching the original summary, augmented with
# the BH q-value AND 95% CI for the per-arm best-enriched test (the test
# already reported in the source table).
df$best_enriched_q_BH <- NA_real_
df$best_enriched_OR_CI95_low <- NA_real_
df$best_enriched_OR_CI95_high <- NA_real_
df$best_enriched_OR_conditional <- NA_real_
for (i in seq_len(nrow(df))) {
  sp <- df$best_enriched_sp[i]
  arm <- df$arm[i]
  hit <- long[long$arm == arm & long$superpop == sp, ]
  if (nrow(hit) == 1L) {
    df$best_enriched_q_BH[i]            <- hit$q_BH
    df$best_enriched_OR_CI95_low[i]     <- hit$OR_CI95_low
    df$best_enriched_OR_CI95_high[i]    <- hit$OR_CI95_high
    df$best_enriched_OR_conditional[i]  <- hit$conditional_MLE_OR
  }
}

write.table(df, out_wide, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat(sprintf("wrote %s (%d arms with BH q + 95%% CI on best-enriched test)\n",
            out_wide, nrow(df)))

# Print a compact summary
cat("\nTop 15 (arm × superpop) tests by q_BH:\n")
print(head(long[, c("arm", "superpop", "carriers_in_sp", "n_carriers_total",
                    "conditional_MLE_OR", "OR_CI95_low", "OR_CI95_high",
                    "p_one_sided", "q_BH")], 15), row.names = FALSE)

cat("\nN tests with q_BH < 0.05:", sum(long$q_BH < 0.05), "\n")
cat("N tests with q_BH < 0.10:", sum(long$q_BH < 0.10), "\n")
