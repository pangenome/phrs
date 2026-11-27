#!/usr/bin/env Rscript
# UPGMA on the arm-level distance matrix as in plot-similarity-subtelo.R (hclust average, k in 2..min(20, arms - 1)
# chosen by the mean silhouette). usage: upgma_arms.R ARM_DIST.tsv OUT.tsv
suppressMessages(library(cluster))
a <- commandArgs(TRUE)
m <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
stopifnot(all(rownames(m) == colnames(m)))
hc <- hclust(as.dist(m), method = "average")
ks <- 2:min(20, nrow(m) - 1)
sil <- sapply(ks, function(k) mean(silhouette(cutree(hc, k = k), as.dist(m))[, "sil_width"]))
k <- ks[which.max(sil)]
message(sprintf("%s: %d arms, UPGMA optimal k %d (silhouette %.3f)", a[1], nrow(m), k, max(sil)))
write.table(data.frame(ChromArm = rownames(m), Community = paste0("U", cutree(hc, k = k))), a[2], sep = "\t", quote = FALSE, row.names = FALSE)
