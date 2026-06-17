#!/usr/bin/env Rscript

# Compute-light arm-level continuum diagnostics for manuscript revision C0a.
# Uses only base R so it can run on the head node without Guix/package setup.

options(stringsAsFactors = FALSE)

out_dir <- "paper_prep/manuscript_revision/C0_continuum"
dist_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
assign_path <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_arm <- function(x) {
  x <- sub("_parm$", "_p", x)
  x <- sub("_qarm$", "_q", x)
  x
}

display_arm <- function(x) {
  sub("^chr", "", gsub("_", "", x))
}

community_num <- function(x) {
  as.integer(sub("^C", "", x))
}

fmt <- function(x, digits = 4) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "NA"
  out
}

read_inputs <- function() {
  d <- read.delim(dist_path, row.names = 1, check.names = FALSE)
  d <- as.matrix(d)
  storage.mode(d) <- "numeric"

  rownames(d) <- normalize_arm(rownames(d))
  colnames(d) <- normalize_arm(colnames(d))
  if (!identical(rownames(d), colnames(d))) {
    stop("Distance matrix row/column arms do not match after normalization")
  }

  assignments <- read.delim(assign_path, check.names = FALSE)
  names(assignments) <- c("arm", "community")
  assignments$arm <- normalize_arm(assignments$arm)
  assignments <- assignments[match(rownames(d), assignments$arm), ]
  if (any(is.na(assignments$community))) {
    stop("Missing Leiden assignments for one or more matrix arms")
  }
  assignments$label <- display_arm(assignments$arm)
  assignments$community_num <- community_num(assignments$community)

  list(distance = d, assignments = assignments)
}

pair_table <- function(sim, assignments) {
  arms <- rownames(sim)
  idx <- which(upper.tri(sim), arr.ind = TRUE)
  pairs <- data.frame(
    arm_a = arms[idx[, 1]],
    arm_b = arms[idx[, 2]],
    similarity = sim[idx],
    distance = 1 - sim[idx],
    stringsAsFactors = FALSE
  )
  pairs$community_a <- assignments$community[match(pairs$arm_a, assignments$arm)]
  pairs$community_b <- assignments$community[match(pairs$arm_b, assignments$arm)]
  pairs$same_community <- pairs$community_a == pairs$community_b
  pairs$pair_label <- paste(display_arm(pairs$arm_a), display_arm(pairs$arm_b), sep = "-")
  pairs
}

summarize_values <- function(values) {
  values <- values[is.finite(values)]
  qs <- quantile(values, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
  data.frame(
    n = length(values),
    mean = mean(values),
    sd = if (length(values) > 1) sd(values) else NA_real_,
    min = qs[1],
    q05 = qs[2],
    q25 = qs[3],
    median = qs[4],
    q75 = qs[5],
    q95 = qs[6],
    max = qs[7],
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

top_pair_label <- function(tbl) {
  if (nrow(tbl) == 0) return(NA_character_)
  tbl <- tbl[order(-tbl$similarity, tbl$pair_label), ]
  paste0(tbl$pair_label[1], " (", fmt(tbl$similarity[1]), ")")
}

top_pairs <- function(tbl, n = 5) {
  if (nrow(tbl) == 0) return(NA_character_)
  tbl <- tbl[order(-tbl$similarity, tbl$pair_label), ]
  paste(paste0(head(tbl$pair_label, n), "=", fmt(head(tbl$similarity, n))), collapse = "; ")
}

system_summary <- function(name, arms, pairs, assignments, note) {
  arms <- normalize_arm(arms)
  sub_pairs <- pairs[pairs$arm_a %in% arms & pairs$arm_b %in% arms, ]
  values <- sub_pairs$similarity
  stats <- summarize_values(values)
  data.frame(
    system = name,
    community = paste(unique(assignments$community[match(arms, assignments$arm)]), collapse = ","),
    arms = paste(display_arm(arms), collapse = ", "),
    n_arms = length(arms),
    pair_count = nrow(sub_pairs),
    mean_similarity = stats$mean,
    median_similarity = stats$median,
    min_similarity = stats$min,
    max_similarity = stats$max,
    peak_pair = top_pair_label(sub_pairs),
    top_pairs = top_pairs(sub_pairs),
    note = note,
    stringsAsFactors = FALSE
  )
}

community_summary <- function(pairs, assignments) {
  out <- do.call(rbind, lapply(split(assignments, assignments$community), function(a) {
    arms <- a$arm
    sub_pairs <- pairs[pairs$arm_a %in% arms & pairs$arm_b %in% arms, ]
    stats <- summarize_values(sub_pairs$similarity)
    data.frame(
      community = a$community[1],
      arms = paste(display_arm(arms), collapse = ", "),
      n_arms = length(arms),
      pair_count = nrow(sub_pairs),
      mean_similarity = stats$mean,
      median_similarity = stats$median,
      min_similarity = stats$min,
      max_similarity = stats$max,
      peak_pair = top_pair_label(sub_pairs),
      top_pairs = top_pairs(sub_pairs, 3),
      stringsAsFactors = FALSE
    )
  }))
  out[order(community_num(out$community)), ]
}

exact_set_test <- function(sim, set_arms) {
  all_arms <- rownames(sim)
  set_idx <- match(set_arms, all_arms)
  k <- length(set_idx)
  observed <- mean(sim[set_idx, set_idx][upper.tri(sim[set_idx, set_idx])])
  means <- combn(seq_along(all_arms), k, FUN = function(ii) {
    mean(sim[ii, ii][upper.tri(sim[ii, ii])])
  })
  data.frame(
    set_size = k,
    n_all_arm_sets = length(means),
    observed_mean = observed,
    null_mean = mean(means),
    null_sd = sd(means),
    null_median = median(means),
    null_q95 = unname(quantile(means, 0.95)),
    exact_p_greater = (sum(means >= observed) + 1) / (length(means) + 1),
    percentile = mean(means <= observed),
    stringsAsFactors = FALSE
  )
}

plot_diagnostics <- function(sim, pairs, assignments, c6_arms, prefix) {
  ordered_assignments <- assignments[order(assignments$community_num, assignments$arm), ]
  ord <- match(ordered_assignments$arm, rownames(sim))
  sim_ord <- sim[ord, ord]

  draw <- function() {
    layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 1))
    par(mar = c(7, 7, 3, 1))
    image(
      x = seq_len(nrow(sim_ord)),
      y = seq_len(ncol(sim_ord)),
      z = t(sim_ord[nrow(sim_ord):1, ]),
      col = hcl.colors(80, "YlOrRd", rev = FALSE),
      axes = FALSE,
      xlab = "",
      ylab = "",
      zlim = c(0, 1),
      main = "Arm-level Jaccard similarity"
    )
    labs <- ordered_assignments$label
    axis(1, at = seq_along(labs), labels = labs, las = 2, cex.axis = 0.62)
    axis(2, at = seq_along(labs), labels = rev(labs), las = 2, cex.axis = 0.62)

    bounds <- cumsum(table(factor(ordered_assignments$community, levels = paste0("C", 1:15))))
    starts <- c(1, head(bounds, -1) + 1)
    for (i in seq_along(bounds)) {
      if (bounds[i] >= starts[i]) {
        rect(starts[i] - 0.5, nrow(sim_ord) - bounds[i] + 0.5,
             bounds[i] + 0.5, nrow(sim_ord) - starts[i] + 1.5,
             border = ifelse(names(bounds)[i] == "C6", "#006d77", "grey20"),
             lwd = ifelse(names(bounds)[i] == "C6", 2.5, 0.8))
      }
    }

    par(mar = c(5, 5, 3, 1))
    breaks <- seq(0, 1, by = 0.05)
    hist(
      pairs$similarity,
      breaks = breaks,
      col = "grey85",
      border = "white",
      xlab = "Arm-pair similarity (1 - Jaccard distance)",
      main = "Off-diagonal distribution"
    )
    hist(pairs$similarity[pairs$same_community], breaks = breaks, col = rgb(0.2, 0.5, 0.8, 0.45), border = NA, add = TRUE)
    c6_pairs <- pairs[pairs$arm_a %in% c6_arms & pairs$arm_b %in% c6_arms, ]
    rug(c6_pairs$similarity, col = "#006d77", lwd = 2)
    legend(
      "topright",
      bty = "n",
      cex = 0.82,
      fill = c("grey85", rgb(0.2, 0.5, 0.8, 0.45), NA),
      border = c("white", NA, NA),
      lty = c(NA, NA, 1),
      col = c(NA, NA, "#006d77"),
      legend = c("all off-diagonal", "within Leiden community", "C6 pair rug")
    )
  }

  png(paste0(prefix, ".png"), width = 2200, height = 1100, res = 180)
  draw()
  dev.off()
  pdf(paste0(prefix, ".pdf"), width = 12, height = 6.2)
  draw()
  dev.off()
}

inputs <- read_inputs()
dist <- inputs$distance
assignments <- inputs$assignments

sim <- 1 - dist
sim[sim < 0 & sim > -1e-10] <- 0
sim[sim > 1 & sim < 1 + 1e-10] <- 1
diag(sim) <- NA_real_

pairs <- pair_table(sim, assignments)
pairs <- pairs[order(pairs$community_a, pairs$community_b, pairs$arm_a, pairs$arm_b), ]

write_tsv(assignments[, c("arm", "label", "community")],
          file.path(out_dir, "arm_assignments_used.tsv"))
write_tsv(pairs, file.path(out_dir, "arm_pair_similarity_long.tsv"))

distribution_summary <- rbind(
  cbind(category = "all_off_diagonal", summarize_values(pairs$similarity)),
  cbind(category = "within_leiden_community", summarize_values(pairs$similarity[pairs$same_community])),
  cbind(category = "between_leiden_communities", summarize_values(pairs$similarity[!pairs$same_community]))
)
write_tsv(distribution_summary, file.path(out_dir, "similarity_distribution_summary.tsv"))

community_stats <- community_summary(pairs, assignments)
write_tsv(community_stats, file.path(out_dir, "community_similarity_summary.tsv"))

named_systems <- do.call(rbind, list(
  system_summary("C1_D4Z4_DUX4", c("chr4_q", "chr10_q"), pairs, assignments, "Known D4Z4/DUX4 4q-10q system"),
  system_summary("C2_10p_18p_TUBB8B", c("chr10_p", "chr18_p"), pairs, assignments, "Known 10p-18p tubulin/TUBB8B system"),
  system_summary("C6_q_arm_sextet", c("chr1_q", "chr13_q", "chr17_q", "chr19_q", "chr21_q", "chr22_q"), pairs, assignments, "Six q-arm community under continuum review"),
  system_summary("C7_acrocentric_p", c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p"), pairs, assignments, "Acrocentric short-arm community"),
  system_summary("C11_OR4F_core_5q_6q", c("chr5_q", "chr6_q"), pairs, assignments, "Browser-viewed OR4F 5q-6q core pair"),
  system_summary("C11_OR4F_full", c("chr1_p", "chr5_q", "chr6_q", "chr8_p"), pairs, assignments, "Full C11 OR4F-bearing arm-level community"),
  system_summary("C14_Xq_Yq", c("chrX_q", "chrY_q"), pairs, assignments, "Sex chromosome q-arm pair"),
  system_summary("C15_Xp_Yp", c("chrX_p", "chrY_p"), pairs, assignments, "Sex chromosome p-arm pair")
))
write_tsv(named_systems, file.path(out_dir, "named_system_peak_similarities.tsv"))

c6_arms <- c("chr1_q", "chr13_q", "chr17_q", "chr19_q", "chr21_q", "chr22_q")
c6_within <- pairs[pairs$arm_a %in% c6_arms & pairs$arm_b %in% c6_arms, ]
c6_to_other <- pairs[xor(pairs$arm_a %in% c6_arms, pairs$arm_b %in% c6_arms), ]
non_c6 <- pairs[!(pairs$arm_a %in% c6_arms | pairs$arm_b %in% c6_arms), ]
background_without_c6_within <- pairs[!(pairs$arm_a %in% c6_arms & pairs$arm_b %in% c6_arms), ]
between_communities <- pairs[!pairs$same_community, ]

c6_density <- rbind(
  cbind(region = "C6_within_sextet", summarize_values(c6_within$similarity)),
  cbind(region = "C6_to_non_C6", summarize_values(c6_to_other$similarity)),
  cbind(region = "non_C6_pairs_only", summarize_values(non_c6$similarity)),
  cbind(region = "off_diagonal_excluding_C6_within", summarize_values(background_without_c6_within$similarity)),
  cbind(region = "between_leiden_communities", summarize_values(between_communities$similarity))
)
c6_density$fold_vs_C6_within <- c6_density$mean[1] / c6_density$mean
write_tsv(c6_density, file.path(out_dir, "c6_neighborhood_density.tsv"))
write_tsv(c6_within[order(-c6_within$similarity), ],
          file.path(out_dir, "c6_within_pair_similarities.tsv"))

c6_test <- exact_set_test(sim, c6_arms)
wilcox_all <- wilcox.test(c6_within$similarity, background_without_c6_within$similarity, alternative = "greater")
wilcox_between <- wilcox.test(c6_within$similarity, between_communities$similarity, alternative = "greater")
c6_test$wilcox_p_vs_offdiag_excluding_C6_within <- wilcox_all$p.value
c6_test$wilcox_p_vs_between_communities <- wilcox_between$p.value
write_tsv(c6_test, file.path(out_dir, "c6_exact_set_test.tsv"))

plot_diagnostics(
  sim,
  pairs,
  assignments,
  c6_arms,
  file.path(out_dir, "arm_level_similarity_diagnostic")
)

report_path <- file.path(out_dir, "C0a_arm_level_report.md")
sink(report_path)
cat("# C0a Arm-Level Continuum Characterization\n\n")
cat("Date: 2026-06-17\n\n")
cat("## Scope\n\n")
cat("This is the compute-light C0a arm-level analysis. It uses the 41 x 41 arm-level Jaccard distance matrix and the arm-level Leiden k=15 assignments only. Similarity is computed as `1 - Jaccard distance` for off-diagonal arm pairs; the nonzero matrix diagonal is ignored. This report does not traverse or summarize the 15,668 x 15,668 sequence-level evidence, so it should be treated as arm-level support for the continuum framing rather than sequence-level proof.\n\n")
cat("## Inputs\n\n")
cat("- Distance matrix: `", dist_path, "`\n", sep = "")
cat("- Leiden assignments: `", assign_path, "`\n", sep = "")
cat("- Matrix arms: ", nrow(dist), "\n", sep = "")
cat("- Off-diagonal arm pairs: ", nrow(pairs), "\n\n", sep = "")

cat("## Arm-Level Similarity Distributions\n\n")
cat("| category | n | mean | median | q05 | q95 | max |\n")
cat("|---|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(distribution_summary))) {
  cat("| ", distribution_summary$category[i], " | ",
      distribution_summary$n[i], " | ",
      fmt(distribution_summary$mean[i]), " | ",
      fmt(distribution_summary$median[i]), " | ",
      fmt(distribution_summary$q05[i]), " | ",
      fmt(distribution_summary$q95[i]), " | ",
      fmt(distribution_summary$max[i]), " |\n", sep = "")
}
cat("\n")
cat("The arm-level matrix shows a two-tier pattern: within-Leiden arm pairs are much more similar on average than between-community pairs, but the off-diagonal background is not a hard zero class. That supports continuum language at arm level: named communities sit on a broad background of lower, variable similarity rather than forming perfectly isolated blocks.\n\n")

cat("## Named Systems\n\n")
cat("| system | arms | pair_count | mean | median | min | max | peak_pair |\n")
cat("|---|---|---:|---:|---:|---:|---:|---|\n")
for (i in seq_len(nrow(named_systems))) {
  cat("| ", named_systems$system[i], " | ",
      named_systems$arms[i], " | ",
      named_systems$pair_count[i], " | ",
      fmt(named_systems$mean_similarity[i]), " | ",
      fmt(named_systems$median_similarity[i]), " | ",
      fmt(named_systems$min_similarity[i]), " | ",
      fmt(named_systems$max_similarity[i]), " | ",
      named_systems$peak_pair[i], " |\n", sep = "")
}
cat("\n")

cat("## C6 / q-arm Sextet Density\n\n")
cat("C6 arms tested: ", paste(display_arm(c6_arms), collapse = ", "), ".\n\n", sep = "")
cat("| region | n | mean | median | q05 | q95 | max | fold_vs_C6_within |\n")
cat("|---|---:|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(c6_density))) {
  cat("| ", c6_density$region[i], " | ",
      c6_density$n[i], " | ",
      fmt(c6_density$mean[i]), " | ",
      fmt(c6_density$median[i]), " | ",
      fmt(c6_density$q05[i]), " | ",
      fmt(c6_density$q95[i]), " | ",
      fmt(c6_density$max[i]), " | ",
      fmt(c6_density$fold_vs_C6_within[i]), " |\n", sep = "")
}
cat("\n")
cat("Exact six-arm set test over all ", c6_test$n_all_arm_sets,
    " possible six-arm subsets: C6 observed mean = ", fmt(c6_test$observed_mean),
    ", null mean = ", fmt(c6_test$null_mean),
    ", null 95th percentile = ", fmt(c6_test$null_q95),
    ", percentile = ", fmt(100 * c6_test$percentile, 2),
    "%, exact greater-tail p = ", formatC(c6_test$exact_p_greater, format = "e", digits = 3),
    ". Wilcoxon greater-tail p versus all off-diagonal pairs excluding C6-within = ",
    formatC(c6_test$wilcox_p_vs_offdiag_excluding_C6_within, format = "e", digits = 3),
    "; versus between-community pairs = ",
    formatC(c6_test$wilcox_p_vs_between_communities, format = "e", digits = 3),
    ".\n\n", sep = "")

cat("Interpretation: the q-arm sextet block is visibly and quantitatively denser than the off-diagonal background in the 41-arm matrix. Because C6 was defined from this same arm-level similarity matrix, the exact set and Wilcoxon tests are descriptive diagnostics for heatmap-density decision-making, not independent discovery p-values. The result supports showing C6 as a dense local region while avoiding language that treats the sextet as a closed, sequence-level clade without the separate sequence-level evidence.\n\n")

cat("## Outputs\n\n")
cat("- `arm_assignments_used.tsv`: normalized arm labels and Leiden assignments used in this run.\n")
cat("- `arm_pair_similarity_long.tsv`: all 820 off-diagonal arm-pair similarities.\n")
cat("- `similarity_distribution_summary.tsv`: all / within-community / between-community similarity distributions.\n")
cat("- `community_similarity_summary.tsv`: per-community arm-level density and peak pairs.\n")
cat("- `named_system_peak_similarities.tsv`: named-system means and peak similarities.\n")
cat("- `c6_neighborhood_density.tsv`: C6 block, C6-neighbor, and background density summaries.\n")
cat("- `c6_within_pair_similarities.tsv`: the 15 C6 sextet arm pairs, sorted by similarity.\n")
cat("- `c6_exact_set_test.tsv`: exact six-arm-set density diagnostic plus Wilcoxon comparisons.\n")
cat("- `arm_level_similarity_diagnostic.png` and `.pdf`: community-ordered heatmap and similarity distribution diagnostic.\n")
cat("- `compute_arm_level_continuum.R`: reproducible script for all outputs.\n")
sink()

message("Wrote C0a arm-level continuum outputs to ", out_dir)
