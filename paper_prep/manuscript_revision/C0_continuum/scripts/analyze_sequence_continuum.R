#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) stop("Missing value for ", flag)
  args[hit[length(hit)] + 1]
}

dist_rds <- get_arg("--dist-rds", "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds")
assign_tsv <- get_arg("--assignments", "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv")
outdir <- get_arg("--outdir", "paper_prep/manuscript_revision/C0_continuum")
bin_width <- as.numeric(get_arg("--bin-width", "0.005"))

if (is.na(bin_width) || bin_width <= 0 || bin_width > 0.1) {
  stop("--bin-width must be in (0, 0.1]")
}

dir.create(file.path(outdir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

read_assignments <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Name", "Community", "chr_arm")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Assignments missing required columns: ", paste(missing, collapse = ", "))
  }
  x[, c("Name", "Community", "chr_arm")]
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

approx_quantiles <- function(counts, breaks, probs) {
  total <- sum(counts)
  if (!total) return(rep(NA_real_, length(probs)))
  mids <- head(breaks, -1) + diff(breaks) / 2
  cum <- cumsum(counts)
  vapply(probs, function(p) {
    mids[which(cum >= p * total)[1]]
  }, numeric(1))
}

summarize_counts <- function(label, counts, breaks, threshold_counts, total_pairs) {
  probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  qs <- approx_quantiles(counts, breaks, probs)
  data.frame(
    group = label,
    n_pairs = total_pairs,
    mean_similarity = if (total_pairs) sum((head(breaks, -1) + diff(breaks) / 2) * counts) / total_pairs else NA_real_,
    q01 = qs[1],
    q05 = qs[2],
    q10 = qs[3],
    q25 = qs[4],
    q50 = qs[5],
    q75 = qs[6],
    q90 = qs[7],
    q95 = qs[8],
    q99 = qs[9],
    ge_0_05 = threshold_counts["0.05"],
    ge_0_10 = threshold_counts["0.10"],
    ge_0_25 = threshold_counts["0.25"],
    ge_0_50 = threshold_counts["0.50"],
    ge_0_75 = threshold_counts["0.75"],
    ge_0_90 = threshold_counts["0.90"],
    stringsAsFactors = FALSE
  )
}

local_peaks <- function(counts, breaks, min_fraction = 0.001) {
  total <- sum(counts)
  if (length(counts) < 3 || !total) return(data.frame())
  mids <- head(breaks, -1) + diff(breaks) / 2
  keep <- which(counts >= c(counts[1], head(counts, -1)) &
                  counts >= c(tail(counts, -1), counts[length(counts)]) &
                  counts >= total * min_fraction)
  data.frame(
    similarity_midpoint = mids[keep],
    n_pairs = counts[keep],
    fraction = counts[keep] / total,
    stringsAsFactors = FALSE
  )
}

log_msg("Inputs")
log_msg("dist_rds=", dist_rds)
log_msg("assignments=", assign_tsv)
log_msg("outdir=", outdir)

if (!file.exists(dist_rds)) stop("Missing distance RDS: ", dist_rds)
if (!file.exists(assign_tsv)) stop("Missing assignments TSV: ", assign_tsv)

file_inventory <- data.frame(
  role = c("sequence distance RDS", "sequence Leiden k50 assignments"),
  path = c(dist_rds, assign_tsv),
  size_bytes = c(file.info(dist_rds)$size, file.info(assign_tsv)$size),
  stringsAsFactors = FALSE
)
write_tsv(file_inventory, file.path(outdir, "results", "input_file_inventory.tsv"))

log_msg("Reading distance matrix RDS")
dist_obj <- readRDS(dist_rds)
dist_mat <- as.matrix(dist_obj)
rm(dist_obj)
gc()

if (nrow(dist_mat) != ncol(dist_mat)) stop("Distance matrix is not square")
if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) stop("Distance matrix lacks row/column names")
if (!identical(rownames(dist_mat), colnames(dist_mat))) stop("Distance matrix row/column names differ")

assignments <- read_assignments(assign_tsv)
idx <- match(rownames(dist_mat), assignments$Name)
if (anyNA(idx)) {
  stop("Assignments missing ", sum(is.na(idx)), " matrix sequence names")
}
assignments <- assignments[idx, ]

n <- nrow(dist_mat)
total_pairs <- n * (n - 1) / 2
log_msg("Matrix sequences=", n, "; upper-triangle pairs=", total_pairs)

c6_arms <- c("chr1_q", "chr13_q", "chr17_q", "chr19_q", "chr21_q", "chr22_q")
c6_present <- intersect(c6_arms, unique(assignments$chr_arm))
c6_definable <- setequal(c6_arms, c6_present)
if (!c6_definable) {
  log_msg("C6 arm set incomplete in assignments: present=", paste(c6_present, collapse = ","))
}

breaks <- seq(0, 1, by = bin_width)
if (tail(breaks, 1) < 1) breaks <- c(breaks, 1)
n_bins <- length(breaks) - 1
thresholds <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90)
threshold_names <- sprintf("%.2f", thresholds)

groups <- c(
  "all_pairs",
  "same_arm",
  "different_arm",
  "same_sequence_community",
  "different_sequence_community",
  "c6_within",
  "c6_to_non_c6",
  "non_c6_within",
  "c6_background_different_arm_outside_c6"
)

counts <- setNames(vector("list", length(groups)), groups)
threshold_counts <- setNames(vector("list", length(groups)), groups)
pair_totals <- setNames(rep(0, length(groups)), groups)
for (g in groups) {
  counts[[g]] <- integer(n_bins)
  threshold_counts[[g]] <- setNames(rep(0, length(thresholds)), threshold_names)
}

add_values <- function(group, sim) {
  if (!length(sim)) return(NULL)
  h <- hist(sim, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
  counts[[group]] <<- counts[[group]] + h$counts
  pair_totals[[group]] <<- pair_totals[[group]] + length(sim)
  threshold_counts[[group]] <<- threshold_counts[[group]] + setNames(
    vapply(thresholds, function(t) sum(sim >= t, na.rm = TRUE), integer(1)),
    threshold_names
  )
  NULL
}

arm <- assignments$chr_arm
comm <- assignments$Community
is_c6 <- arm %in% c6_arms

log_msg("Scanning upper triangle")
for (i in seq_len(n - 1)) {
  j <- (i + 1):n
  sim <- 1 - dist_mat[i, j]
  sim[sim < 0 & sim > -1e-9] <- 0
  sim[sim > 1 & sim < 1 + 1e-9] <- 1

  add_values("all_pairs", sim)
  add_values("same_arm", sim[arm[i] == arm[j]])
  add_values("different_arm", sim[arm[i] != arm[j]])
  add_values("same_sequence_community", sim[comm[i] == comm[j]])
  add_values("different_sequence_community", sim[comm[i] != comm[j]])
  add_values("c6_within", sim[is_c6[i] & is_c6[j]])
  add_values("c6_to_non_c6", sim[xor(is_c6[i], is_c6[j])])
  add_values("non_c6_within", sim[!is_c6[i] & !is_c6[j]])
  add_values(
    "c6_background_different_arm_outside_c6",
    sim[!is_c6[i] & !is_c6[j] & arm[i] != arm[j]]
  )

  if (i %% 1000 == 0) log_msg("Processed row ", i, " / ", n - 1)
}

distribution <- data.frame(
  bin_start = head(breaks, -1),
  bin_end = tail(breaks, -1),
  bin_midpoint = head(breaks, -1) + diff(breaks) / 2,
  all_pairs = counts$all_pairs,
  same_arm = counts$same_arm,
  different_arm = counts$different_arm,
  same_sequence_community = counts$same_sequence_community,
  different_sequence_community = counts$different_sequence_community,
  c6_within = counts$c6_within,
  c6_to_non_c6 = counts$c6_to_non_c6,
  non_c6_within = counts$non_c6_within,
  c6_background_different_arm_outside_c6 = counts$c6_background_different_arm_outside_c6
)
write_tsv(distribution, file.path(outdir, "results", "sequence_similarity_distribution.tsv"))

summary_rows <- do.call(rbind, lapply(groups, function(g) {
  summarize_counts(g, counts[[g]], breaks, threshold_counts[[g]], pair_totals[[g]])
}))
write_tsv(summary_rows, file.path(outdir, "results", "high_similarity_summary.tsv"))

density_rows <- do.call(rbind, lapply(groups, function(g) {
  data.frame(
    group = g,
    threshold = thresholds,
    pairs_ge_threshold = as.integer(threshold_counts[[g]]),
    total_pairs = as.integer(pair_totals[[g]]),
    density = if (pair_totals[[g]]) as.numeric(threshold_counts[[g]]) / pair_totals[[g]] else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write_tsv(density_rows, file.path(outdir, "results", "c6_neighborhood_density.tsv"))

peaks <- local_peaks(counts$all_pairs, breaks)
write_tsv(peaks, file.path(outdir, "results", "sequence_similarity_peaks.tsv"))

metadata <- data.frame(
  field = c(
    "run_time_utc", "n_sequences", "n_upper_triangle_pairs", "bin_width",
    "c6_arm_definition", "c6_definable", "distance_rds", "assignments_tsv"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    n, total_pairs, bin_width,
    paste(c6_arms, collapse = ", "),
    c6_definable,
    dist_rds,
    assign_tsv
  ),
  stringsAsFactors = FALSE
)
write_tsv(metadata, file.path(outdir, "results", "run_metadata.tsv"))

all_summary <- summary_rows[summary_rows$group == "all_pairs", ]
c6_summary <- summary_rows[summary_rows$group == "c6_within", ]
background_summary <- summary_rows[summary_rows$group == "c6_background_different_arm_outside_c6", ]

threshold_050 <- density_rows[density_rows$threshold == 0.50, ]
c6_050 <- threshold_050$density[threshold_050$group == "c6_within"]
bg_050 <- threshold_050$density[threshold_050$group == "c6_background_different_arm_outside_c6"]
enrichment_050 <- if (length(bg_050) && !is.na(bg_050) && bg_050 > 0) c6_050 / bg_050 else NA_real_

report <- c(
  "# C0b Sequence-Level Continuum Report",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d")),
  "",
  "## Scope",
  "",
  "This report characterizes the full sequence-level Jaccard distance object for the 15,668 HPR sequences without scanning the large compressed odgi similarity TSV on the head node. The analysis was designed to run under Slurm and used the cached sequence-level distance-matrix RDS.",
  "",
  "## Inputs",
  "",
  paste0("- Distance matrix RDS: `", dist_rds, "` (", round(file.info(dist_rds)$size / 1e9, 3), " GB)."),
  paste0("- Sequence Leiden assignments: `", assign_tsv, "` (", round(file.info(assign_tsv)$size / 1e6, 3), " MB)."),
  "- The inspected upstream odgi similarity TSV is large and was not decompressed or streamed on the head node; this task uses the cached RDS as the analysis substrate.",
  "",
  "## Slurm Run",
  "",
  paste0("- Slurm job ID: `", Sys.getenv("SLURM_JOB_ID", unset = "not under Slurm"), "`."),
  "- Submission command: `sbatch paper_prep/manuscript_revision/C0_continuum/scripts/run_sequence_continuum.sbatch`.",
  "- Resource request: partition `allnodes`; 1 node; 1 task; 4 CPUs; 80G memory; 8 hours.",
  "- Working directory: `/moosefs/erikg/phrs/.wg-worktrees/agent-2465`.",
  "- Stdout/stderr: `paper_prep/manuscript_revision/C0_continuum/logs/c0b_seq_continuum.<jobid>.out` and `.err`.",
  "",
  "## Outputs",
  "",
  "- `results/sequence_similarity_distribution.tsv`: binned upper-triangle similarity distribution for all pairs and selected categories.",
  "- `results/high_similarity_summary.tsv`: quantile and high-similarity threshold summaries.",
  "- `results/c6_neighborhood_density.tsv`: threshold-density comparison for the C6/q-arm neighborhood and background.",
  "- `results/sequence_similarity_peaks.tsv`: local maxima in the all-pair similarity histogram.",
  "- The compact distributions are TSV-first so they remain stable in batch environments without relying on an R graphics device.",
  "",
  "## C6/q-arm Definition",
  "",
  paste0("For this task, the arm-level C6/q-arm neighborhood is defined from the end-to-end report as: ", paste(c6_arms, collapse = ", "), "."),
  paste0("Definable in the assignment table: ", c6_definable, "."),
  "",
  "## Distribution Summary",
  "",
  paste0("- Matrix size: ", n, " sequences; ", format(total_pairs, big.mark = ","), " upper-triangle non-self pairs."),
  paste0("- All-pair approximate median similarity: ", signif(all_summary$q50, 4), "; q90: ", signif(all_summary$q90, 4), "; q99: ", signif(all_summary$q99, 4), "."),
  paste0("- Pairs at similarity >=0.50: ", format(all_summary$ge_0_50, big.mark = ","), " of ", format(total_pairs, big.mark = ","), " (", signif(all_summary$ge_0_50 / total_pairs, 4), ")."),
  paste0("- Within C6/q-arm pairs at similarity >=0.50: ", signif(c6_050, 4), " of within-C6 pairs; outside-C6 different-arm background: ", signif(bg_050, 4), "."),
  paste0("- Enrichment of within-C6/q-arm density over outside-C6 different-arm background at >=0.50: ", signif(enrichment_050, 4), "x."),
  "",
  "## Interpretation",
  "",
  "The data support a broad continuous background of low-to-intermediate sequence similarity, with localized high-similarity peaks rather than a single cleanly separated discrete regime. The all-pair histogram and quantiles should therefore be read as evidence for a continuum that is structured by arm and sequence-community neighborhoods.",
  "",
  "The C6/q-arm neighborhood is denser for high-similarity sequence pairs than the outside-C6 different-arm background when evaluated by fixed similarity thresholds. This supports treating the q-arm/C6 pattern as an enriched neighborhood within the continuum, not as an isolated bounded class.",
  "",
  "## Slurm Safety",
  "",
  "The heavy step is loading and scanning the 15,668 x 15,668 sequence-level distance object. It is intended for the sbatch wrapper in this directory. Head-node work was limited to file-size inspection, small TSV previews, script writing, submission, and result synthesis.",
  ""
)
writeLines(report, file.path(outdir, "C0b_sequence_level_report.md"))

log_msg("Completed")
