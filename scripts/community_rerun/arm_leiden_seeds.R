#!/usr/bin/env Rscript
# Arm-level Leiden as in detect_communities.R (similarity exp(-d / median nonzero d), resolution scan 0.1..3 by 0.01
# chosen by mean silhouette, then one more cluster_leiden at that resolution) repeated with seeds 1..N, on an arm
# distance matrix TSV. Per seed: scan optimum and the final partition. Summary: how often each final partition occurs,
# and the seeds whose final partition equals the reference assignments (e.g. the paper's).
# usage: arm_leiden_seeds.R ARM_DIST.tsv REFERENCE.tsv N OUT_PREFIX
suppressMessages({library(igraph); library(cluster)})
a <- commandArgs(TRUE); n_seeds <- as.integer(a[3])
m <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)); stopifnot(all(rownames(m) == colnames(m)))
ref <- read.delim(a[2]); ref <- setNames(as.character(ref[[2]]), ref[[1]])
stopifnot(setequal(names(ref), rownames(m)))
canon <- function(memb) {  # partition as sorted groups, independent of labels
  g <- split(names(memb), memb); g <- lapply(g, sort); paste(sort(sapply(g, paste, collapse = "+")), collapse = " | ")
}
ref_key <- canon(ref[rownames(m)])
sim <- exp(-m / median(m[m > 0])); diag(sim) <- 0
g <- graph_from_adjacency_matrix(sim, mode = "undirected", weighted = TRUE, diag = FALSE)
res_range <- seq(0.1, 3.0, by = 0.01)
rows <- vector("list", n_seeds)
for (s in seq_len(n_seeds)) {
  set.seed(s)
  sil <- sapply(res_range, function(r) {
    cm <- membership(cluster_leiden(g, resolution = r, objective_function = "modularity")); k <- max(cm)
    if (k >= 2 && k < nrow(m)) mean(silhouette(as.integer(cm), as.dist(m))[, 3]) else NA
  })
  best <- res_range[which.max(sil)]
  fin <- membership(cluster_leiden(g, resolution = best, objective_function = "modularity")); names(fin) <- rownames(m)
  rows[[s]] <- data.frame(seed = s, resolution = best, scan_silhouette = max(sil, na.rm = TRUE), k = max(fin),
                          final_silhouette = mean(silhouette(as.integer(fin), as.dist(m))[, 3]),
                          same_as_reference = canon(fin) == ref_key, partition = canon(fin))
}
d <- do.call(rbind, rows)
write.table(d, paste0(a[4], ".seeds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tab <- sort(table(d$partition), decreasing = TRUE)
message(sprintf("%s: %d seeds, %d distinct final partitions; reference partition in %d seeds (first: %s)", a[1], n_seeds,
                length(tab), sum(d$same_as_reference), paste(head(d$seed[d$same_as_reference], 5), collapse = ",")))
for (i in seq_along(tab)) {
  x <- d[d$partition == names(tab)[i], ][1, ]
  message(sprintf("  partition %d: %d seeds, k %d, final silhouette %.4f, reference %s", i, tab[i], x$k, x$final_silhouette, x$same_as_reference))
}
