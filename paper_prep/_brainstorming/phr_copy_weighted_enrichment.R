#!/usr/bin/env Rscript

# Copy-Number-Aware Enrichment Analysis for PHR Data
# Integration testing for copy-weighted ORA methodology

library(data.table)
library(dplyr)
library(readr)

cat("=== Copy-Number-Aware PHR Enrichment Analysis ===\n\n")

# Step 1: Load input data
cat("Step 1: Loading input data...\n")

# Gene copy summary (PHR gene families and their copy counts)
gene_copy_summary <- read_csv("gene_copy_summary.csv")
cat("Loaded", nrow(gene_copy_summary), "gene families with PHR copies\n")

# All gene copies with genomic locations
all_gene_copies <- read_csv("all_gene_copies_by_arm.csv")
cat("Loaded", nrow(all_gene_copies), "individual gene copies\n")

# Previous ORA results for comparison
phr_bp_results <- read_csv("phr_coding_only_GO_BP_enrichment.csv")
phr_mf_results <- read_csv("phr_coding_only_GO_MF_enrichment.csv")
cat("Loaded previous ORA results:", nrow(phr_bp_results), "BP terms,", nrow(phr_mf_results), "MF terms\n")

# Step 2: Build genome-wide copy count background
cat("\nStep 2: Building genome-wide copy count background...\n")

# Load full genome annotation to count genome-wide gene copies
# We'll need to process the GFF3 file to get all gene copies genome-wide
gff_file <- "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"

# For now, let's work with the gene families we have in PHRs and count their total copies
# This is a simplified approach - we'll use the total_copies column from gene_copy_summary

# Create a summary of copy counts in PHRs vs genome-wide
copy_background <- gene_copy_summary %>%
  select(gene_name, total_copies, gene_biotype) %>%
  rename(genome_wide_copies = total_copies)

# Count copies actually in PHRs for each gene family
phr_copies_per_family <- all_gene_copies %>%
  group_by(gene_family) %>%
  summarise(
    copies_in_phrs = n(),
    .groups = 'drop'
  )

# Merge background and PHR counts
copy_analysis <- copy_background %>%
  left_join(phr_copies_per_family, by = c("gene_name" = "gene_family")) %>%
  mutate(copies_in_phrs = ifelse(is.na(copies_in_phrs), 0, copies_in_phrs))

cat("Summary of copy counts:\n")
cat("Total gene families:", nrow(copy_analysis), "\n")
cat("Total genome-wide copies:", sum(copy_analysis$genome_wide_copies), "\n")
cat("Total PHR copies:", sum(copy_analysis$copies_in_phrs), "\n")

# Step 3: Prepare GO term mapping for copy-weighted enrichment
cat("\nStep 3: Preparing GO term mapping...\n")

# For this integration test, we'll focus on the gene families that were enriched
# in the previous analysis and see how copy-weighting affects the results

# Extract gene names from the PHR gene families
phr_gene_families <- copy_analysis$gene_name

# Create a mock GO term mapping for the gene families we have
# In a real analysis, this would come from a comprehensive GO annotation database
# For now, we'll simulate this based on the enriched terms from previous analysis

# Mock GO term assignments based on gene family names and known biology
go_mapping <- data.frame(
  gene_name = character(),
  go_id = character(),
  go_name = character(),
  go_domain = character(),
  stringsAsFactors = FALSE
)

# Add some known associations based on gene families
# DUX4, FRG2*, IL9R* families - transcription regulation and immune response
transcription_genes <- c("DUX4", "FRG2", "FRG2B")
immune_genes <- c("IL9R", "IL9RP1", "IL9RP3", "IL9RP4")
gtp_genes <- c("GTPBP6", "IQSEC3")

# Create GO mappings
for(gene in transcription_genes) {
  go_mapping <- rbind(go_mapping, data.frame(
    gene_name = gene,
    go_id = "GO:0006355",
    go_name = "regulation of transcription, DNA-templated",
    go_domain = "BP"
  ))
}

for(gene in immune_genes) {
  go_mapping <- rbind(go_mapping, data.frame(
    gene_name = gene,
    go_id = "GO:0006955",
    go_name = "immune response",
    go_domain = "BP"
  ))
  go_mapping <- rbind(go_mapping, data.frame(
    gene_name = gene,
    go_id = "GO:0005125",
    go_name = "cytokine activity",
    go_domain = "MF"
  ))
}

for(gene in gtp_genes) {
  go_mapping <- rbind(go_mapping, data.frame(
    gene_name = gene,
    go_id = "GO:0005525",
    go_name = "GTP binding",
    go_domain = "MF"
  ))
  go_mapping <- rbind(go_mapping, data.frame(
    gene_name = gene,
    go_id = "GO:0003924",
    go_name = "GTPase activity",
    go_domain = "MF"
  ))
}

cat("Created GO mappings for", length(unique(go_mapping$gene_name)), "gene families\n")
cat("Covering", length(unique(go_mapping$go_id)), "unique GO terms\n")

# Step 4: Copy-weighted hypergeometric test
cat("\nStep 4: Running copy-weighted hypergeometric tests...\n")

copy_weighted_results <- data.frame()

# Get unique GO terms
unique_go_terms <- unique(go_mapping[, c("go_id", "go_name", "go_domain")])

for(i in 1:nrow(unique_go_terms)) {
  go_term <- unique_go_terms$go_id[i]
  go_name <- unique_go_terms$go_name[i]
  go_domain <- unique_go_terms$go_domain[i]

  # Get genes annotated to this GO term
  go_genes <- go_mapping$gene_name[go_mapping$go_id == go_term]

  # Calculate copy-weighted parameters for phyper()
  # q = copies of GO-term genes drawn into PHRs
  q_copies <- copy_analysis %>%
    filter(gene_name %in% go_genes) %>%
    summarise(sum(copies_in_phrs)) %>%
    pull()

  # m = total copies of GO-term genes genome-wide
  m_copies <- copy_analysis %>%
    filter(gene_name %in% go_genes) %>%
    summarise(sum(genome_wide_copies)) %>%
    pull()

  # Total gene copies in PHRs
  k_total <- sum(copy_analysis$copies_in_phrs)

  # Total gene copies genome-wide
  N_total <- sum(copy_analysis$genome_wide_copies)

  # n = total gene copies genome-wide NOT in this GO term
  n_copies <- N_total - m_copies

  # Skip if no copies or invalid parameters
  if(q_copies == 0 || m_copies == 0 || n_copies <= 0) {
    next
  }

  # Run hypergeometric test
  # P(X >= q) where X ~ Hypergeometric(k, m, n)
  p_value <- phyper(q_copies - 1, m_copies, n_copies, k_total, lower.tail = FALSE)

  # Calculate enrichment metrics
  expected_copies <- k_total * (m_copies / N_total)
  fold_enrichment <- q_copies / expected_copies

  # Store results
  copy_weighted_results <- rbind(copy_weighted_results, data.frame(
    go_id = go_term,
    go_name = go_name,
    go_domain = go_domain,
    copies_in_phrs = q_copies,
    total_genome_copies = m_copies,
    expected_copies = expected_copies,
    fold_enrichment = fold_enrichment,
    p_value = p_value,
    total_phr_copies = k_total,
    total_genome_copies_all = N_total
  ))
}

# Apply multiple testing correction
copy_weighted_results$p_adjusted <- p.adjust(copy_weighted_results$p_value, method = "BH")

# Sort by p-value
copy_weighted_results <- copy_weighted_results[order(copy_weighted_results$p_value), ]

cat("Completed copy-weighted enrichment analysis\n")
cat("Found", nrow(copy_weighted_results), "testable GO terms\n")
cat("Significant results (p < 0.05):", sum(copy_weighted_results$p_value < 0.05), "\n")
cat("Significant after correction (p_adj < 0.05):", sum(copy_weighted_results$p_adjusted < 0.05), "\n")

# Step 5: Compare to previous results
cat("\nStep 5: Comparing to previous deduplicated ORA results...\n")

# Create comparison for overlapping terms
comparison_results <- data.frame()

# Check BP terms
for(i in 1:nrow(phr_bp_results)) {
  prev_term <- phr_bp_results$native[i]
  prev_pval <- phr_bp_results$p_value[i]
  prev_genes <- phr_bp_results$intersection_size[i]

  # Find matching term in copy-weighted results
  copy_match <- copy_weighted_results[copy_weighted_results$go_id == prev_term, ]

  if(nrow(copy_match) > 0) {
    comparison_results <- rbind(comparison_results, data.frame(
      go_id = prev_term,
      go_name = phr_bp_results$name[i],
      domain = "BP",
      prev_pvalue = prev_pval,
      prev_genes = prev_genes,
      copy_pvalue = copy_match$p_value[1],
      copy_copies = copy_match$copies_in_phrs[1],
      fold_change_pvalue = copy_match$p_value[1] / prev_pval
    ))
  }
}

# Check MF terms
for(i in 1:nrow(phr_mf_results)) {
  prev_term <- phr_mf_results$native[i]
  prev_pval <- phr_mf_results$p_value[i]
  prev_genes <- phr_mf_results$intersection_size[i]

  # Find matching term in copy-weighted results
  copy_match <- copy_weighted_results[copy_weighted_results$go_id == prev_term, ]

  if(nrow(copy_match) > 0) {
    comparison_results <- rbind(comparison_results, data.frame(
      go_id = prev_term,
      go_name = phr_mf_results$name[i],
      domain = "MF",
      prev_pvalue = prev_pval,
      prev_genes = prev_genes,
      copy_pvalue = copy_match$p_value[1],
      copy_copies = copy_match$copies_in_phrs[1],
      fold_change_pvalue = copy_match$p_value[1] / prev_pval
    ))
  }
}

cat("Found", nrow(comparison_results), "overlapping terms for comparison\n")

# Step 6: Save results
cat("\nStep 6: Saving results...\n")

write_csv(copy_weighted_results, "phr_copy_weighted_enrichment.csv")
cat("Saved copy-weighted enrichment results to phr_copy_weighted_enrichment.csv\n")

write_csv(comparison_results, "copy_weighted_vs_deduplicated_comparison.csv")
cat("Saved comparison results to copy_weighted_vs_deduplicated_comparison.csv\n")

write_csv(copy_analysis, "gene_copy_background_analysis.csv")
cat("Saved gene copy background analysis to gene_copy_background_analysis.csv\n")

# Step 7: Summary report
cat("\n=== SUMMARY REPORT ===\n")
cat("Copy-Number-Aware vs Deduplicated ORA Comparison:\n\n")

cat("Data Overview:\n")
cat("- Gene families in PHRs:", nrow(copy_analysis), "\n")
cat("- Total genome-wide copies:", sum(copy_analysis$genome_wide_copies), "\n")
cat("- Total PHR copies:", sum(copy_analysis$copies_in_phrs), "\n")
cat("- Average copies per family:", round(mean(copy_analysis$genome_wide_copies), 2), "\n\n")

cat("Enrichment Results:\n")
cat("- Copy-weighted testable terms:", nrow(copy_weighted_results), "\n")
cat("- Copy-weighted significant (p < 0.05):", sum(copy_weighted_results$p_value < 0.05), "\n")
cat("- Copy-weighted significant after correction:", sum(copy_weighted_results$p_adjusted < 0.05), "\n\n")

if(nrow(comparison_results) > 0) {
  cat("Direct Comparison (overlapping terms):\n")
  cat("- Overlapping terms found:", nrow(comparison_results), "\n")

  stronger_terms <- sum(comparison_results$fold_change_pvalue < 1, na.rm = TRUE)
  weaker_terms <- sum(comparison_results$fold_change_pvalue > 1, na.rm = TRUE)

  cat("- Terms with stronger significance (lower p-value):", stronger_terms, "\n")
  cat("- Terms with weaker significance (higher p-value):", weaker_terms, "\n")
}

cat("\n=== Analysis Complete ===\n")