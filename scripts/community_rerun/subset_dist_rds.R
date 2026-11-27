#!/usr/bin/env Rscript
# Drop sequences from the cached Jaccard distance data frame (build_dist_matrix_rds.R output, square, row and column
# names = PHR sequence names). usage: subset_dist_rds.R IN.rds DROP.txt OUT.rds
a <- commandArgs(TRUE)
d <- readRDS(a[1]); drop <- readLines(a[2])
stopifnot("square with matching names" = nrow(d) == ncol(d) && all(rownames(d) == colnames(d)))
miss <- setdiff(drop, rownames(d))
if (length(miss)) stop(sprintf("FATAL: %d names to drop are not in the matrix, e.g. %s", length(miss), miss[1]))
keep <- setdiff(rownames(d), drop)
d2 <- d[keep, keep]
message(sprintf("%s: %d sequences; dropped %d; kept %d", a[1], nrow(d), length(drop), nrow(d2)))
saveRDS(d2, a[3])
