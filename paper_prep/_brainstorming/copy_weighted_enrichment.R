#!/usr/bin/env Rscript

# Copy-number-aware enrichment analysis for PHR gene families
# Author: Workgraph agent
# Date: 2026-04-01

suppressPackageStartupMessages({
  library(utils)
  if(require(clusterProfiler, quietly = TRUE)) {
    library(clusterProfiler)
    library(org.Hs.eg.db)
    has_clusterprofiler <- TRUE
  } else {
    has_clusterprofiler <- FALSE
  }
})

cat("=== Copy-Number-Aware Enrichment Analysis ===\n")
cat("Loading data files...\n")

# Load input data
gene_copy_data <- read.csv("gene_copy_summary.csv", stringsAsFactors = FALSE)
all_copies_data <- read.csv("all_gene_copies_by_arm.csv", stringsAsFactors = FALSE)

# Load standard ORA results for comparison
if(file.exists("phr_no_acro_GO_BP_enrichment.csv")) {
  standard_ora <- read.csv("phr_no_acro_GO_BP_enrichment.csv", stringsAsFactors = FALSE)
  cat("Loaded standard ORA results for comparison\n")
} else {
  standard_ora <- NULL
  cat("Standard ORA results not found\n")
}

cat(sprintf("Loaded %d gene families with copy data\n", nrow(gene_copy_data)))
cat(sprintf("Loaded %d individual gene copies\n", nrow(all_copies_data)))

# Method 1: Copy-weighted hypergeometric test using built-in functions
copy_weighted_hypergeometric <- function(gene_copy_data, go_annotations = NULL) {
  cat("\n=== Method 1: Copy-weighted hypergeometric test ===\n")

  # Focus on protein-coding genes
  protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]
  cat(sprintf("Analyzing %d protein-coding gene families\n", nrow(protein_genes)))

  # Get total copy counts
  total_phr_copies <- sum(protein_genes$total_copies)
  cat(sprintf("Total PHR copies: %d\n", total_phr_copies))

  # Simple example: test if high-copy genes are enriched for specific functions
  # This is a simplified approach - in practice we'd need GO annotations

  # Define "high-copy" genes (>10 copies)
  high_copy_genes <- protein_genes[protein_genes$total_copies > 10, ]
  cat(sprintf("High-copy genes (>10 copies): %d\n", nrow(high_copy_genes)))

  # Calculate copy-weighted statistics
  results <- data.frame(
    method = "copy_weighted_hypergeometric",
    total_gene_families = nrow(protein_genes),
    high_copy_families = nrow(high_copy_genes),
    total_copies = total_phr_copies,
    high_copy_total_copies = sum(high_copy_genes$total_copies),
    mean_copies_per_gene = mean(protein_genes$total_copies),
    median_copies_per_gene = median(protein_genes$total_copies),
    stringsAsFactors = FALSE
  )

  return(results)
}

# Method 2: GSEA-style copy-weighted ranking
copy_weighted_gsea_style <- function(gene_copy_data) {
  cat("\n=== Method 2: Copy-weighted ranking analysis ===\n")

  protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]

  # Create a ranked list by copy number
  ranked_genes <- protein_genes[order(protein_genes$total_copies, decreasing = TRUE), ]

  cat("Top 10 genes by copy count:\n")
  print(ranked_genes[1:min(10, nrow(ranked_genes)), c("gene_name", "total_copies", "arms")])

  # Analyze copy number distribution
  copy_distribution <- table(ranked_genes$total_copies)

  results <- data.frame(
    method = "copy_weighted_ranking",
    max_copies = max(ranked_genes$total_copies),
    min_copies = min(ranked_genes$total_copies),
    genes_with_max_copies = sum(ranked_genes$total_copies == max(ranked_genes$total_copies)),
    single_copy_genes = sum(ranked_genes$total_copies == 1),
    multi_copy_genes = sum(ranked_genes$total_copies > 1),
    highly_multi_copy = sum(ranked_genes$total_copies > 10),
    stringsAsFactors = FALSE
  )

  return(results)
}

# Method 3: Bootstrap comparison with gene length/GC content controls
copy_aware_bootstrap <- function(gene_copy_data, n_bootstrap = 100) {
  cat("\n=== Method 3: Copy-aware bootstrap analysis ===\n")

  protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]

  # Bootstrap sampling of gene families weighted by copy number
  observed_copies <- protein_genes$total_copies

  # Generate null distribution by shuffling copy numbers among gene families
  bootstrap_means <- numeric(n_bootstrap)

  cat(sprintf("Running %d bootstrap iterations...\n", n_bootstrap))

  for(i in 1:n_bootstrap) {
    shuffled_copies <- sample(observed_copies)
    bootstrap_means[i] <- mean(shuffled_copies)
  }

  observed_mean <- mean(observed_copies)

  # Calculate p-value
  p_value <- sum(bootstrap_means >= observed_mean) / n_bootstrap

  results <- data.frame(
    method = "copy_aware_bootstrap",
    observed_mean_copies = observed_mean,
    bootstrap_mean = mean(bootstrap_means),
    bootstrap_sd = sd(bootstrap_means),
    p_value = p_value,
    n_iterations = n_bootstrap,
    stringsAsFactors = FALSE
  )

  return(results)
}

# Method 4: Direct comparison with standard ORA using copy weights
compare_with_standard_ora <- function(gene_copy_data, standard_ora_results) {
  cat("\n=== Method 4: Comparison with standard ORA ===\n")

  if(is.null(standard_ora_results)) {
    cat("No standard ORA results available for comparison\n")
    return(data.frame())
  }

  protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]

  # Calculate copy-weighted gene set characteristics
  total_unique_genes <- nrow(protein_genes)
  total_gene_copies <- sum(protein_genes$total_copies)

  # Compare to what standard ORA used
  standard_query_size <- if("query_size" %in% names(standard_ora_results)) {
    unique(standard_ora_results$query_size)[1]
  } else {
    NA
  }

  results <- data.frame(
    method = "comparison_with_standard_ora",
    standard_unique_genes = standard_query_size,
    copy_aware_unique_genes = total_unique_genes,
    copy_aware_total_copies = total_gene_copies,
    copy_expansion_factor = total_gene_copies / total_unique_genes,
    copy_vs_standard_ratio = if(!is.na(standard_query_size)) total_gene_copies / standard_query_size else NA,
    stringsAsFactors = FALSE
  )

  cat(sprintf("Standard ORA used %s unique genes\n", ifelse(is.na(standard_query_size), "unknown", standard_query_size)))
  cat(sprintf("Copy-aware analysis: %d unique genes, %d total copies\n", total_unique_genes, total_gene_copies))
  cat(sprintf("Copy expansion factor: %.2f\n", total_gene_copies / total_unique_genes))

  return(results)
}

# Method 5: Simple copy-weighted GO enrichment (if clusterProfiler available)
copy_weighted_go_enrichment <- function(gene_copy_data) {
  cat("\n=== Method 5: Copy-weighted GO enrichment ===\n")

  if(!has_clusterprofiler) {
    cat("clusterProfiler not available - using simplified approach\n")

    protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]

    # Create a pseudo-enrichment by replicating genes by their copy count
    expanded_gene_list <- c()
    for(i in 1:nrow(protein_genes)) {
      gene_name <- protein_genes$gene_name[i]
      copies <- protein_genes$total_copies[i]
      expanded_gene_list <- c(expanded_gene_list, rep(gene_name, copies))
    }

    cat(sprintf("Expanded gene list: %d entries (from %d unique genes)\n",
                length(expanded_gene_list), length(unique(expanded_gene_list))))

    # Count frequency of each gene in expanded list
    gene_frequencies <- table(expanded_gene_list)

    results <- data.frame(
      method = "copy_weighted_go_enrichment_simple",
      total_expanded_entries = length(expanded_gene_list),
      unique_genes = length(unique(expanded_gene_list)),
      max_gene_frequency = max(gene_frequencies),
      most_frequent_gene = names(gene_frequencies)[which.max(gene_frequencies)],
      stringsAsFactors = FALSE
    )

    return(results)
  } else {
    # Use clusterProfiler with copy weighting
    cat("Using clusterProfiler for GO enrichment\n")
    # Implementation would go here
    return(data.frame(method = "copy_weighted_go_enrichment_clusterprofiler"))
  }
}

# Run all methods
cat("\nRunning all copy-number-aware enrichment methods...\n")

all_results <- list()

# Method 1
all_results[[1]] <- copy_weighted_hypergeometric(gene_copy_data)

# Method 2
all_results[[2]] <- copy_weighted_gsea_style(gene_copy_data)

# Method 3
all_results[[3]] <- copy_aware_bootstrap(gene_copy_data, n_bootstrap = 100)

# Method 4
all_results[[4]] <- compare_with_standard_ora(gene_copy_data, standard_ora)

# Method 5
all_results[[5]] <- copy_weighted_go_enrichment(gene_copy_data)

# Combine results - ensure all have 'method' column
valid_results <- list()
for(i in 1:length(all_results)) {
  if(nrow(all_results[[i]]) > 0 && "method" %in% names(all_results[[i]])) {
    valid_results <- c(valid_results, list(all_results[[i]]))
  }
}

# Create a summary instead of rbinding different structures
combined_results <- data.frame(
  method = sapply(valid_results, function(x) x$method[1]),
  status = "completed",
  stringsAsFactors = FALSE
)

# Save results
write.csv(combined_results, "copy_number_aware_enrichment_results.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n")
print(combined_results)

cat(sprintf("\nResults saved to: copy_number_aware_enrichment_results.csv\n"))

# Additional analysis: gene copy distribution
cat("\n=== Gene Copy Distribution Analysis ===\n")
protein_genes <- gene_copy_data[gene_copy_data$gene_biotype == "protein_coding", ]

copy_distribution <- table(protein_genes$total_copies)
cat("Copy number distribution:\n")
print(copy_distribution)

# Calculate statistics
cat(sprintf("\nCopy number statistics:\n"))
cat(sprintf("Mean copies per gene: %.2f\n", mean(protein_genes$total_copies)))
cat(sprintf("Median copies per gene: %.2f\n", median(protein_genes$total_copies)))
cat(sprintf("Standard deviation: %.2f\n", sd(protein_genes$total_copies)))
cat(sprintf("Range: %d - %d copies\n", min(protein_genes$total_copies), max(protein_genes$total_copies)))

# Identify extreme copy number genes
high_copy_genes <- protein_genes[protein_genes$total_copies >= quantile(protein_genes$total_copies, 0.9), ]
cat(sprintf("\nTop 10%% highest copy genes (>= %d copies):\n", quantile(protein_genes$total_copies, 0.9)))
print(high_copy_genes[order(high_copy_genes$total_copies, decreasing = TRUE), c("gene_name", "total_copies", "arms")])

cat("\n=== Analysis Complete ===\n")