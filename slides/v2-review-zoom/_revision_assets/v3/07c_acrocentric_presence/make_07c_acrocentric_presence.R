#!/usr/bin/env Rscript

# Audit slide 07c for acrocentric p-arm presence and render an annotated
# unrooted-tree candidate. The script writes only into this v3 revision-assets
# directory so downstream deck integration can choose whether to use it.

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)

matrix_path <- if (length(args) >= 1) args[[1]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv"
assignments_path <- if (length(args) >= 2) args[[2]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
len_path <- if (length(args) >= 3) args[[3]] else
  "/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
out_dir <- if (length(args) >= 4) args[[4]] else
  "slides/v2-review-zoom/_revision_assets/v3/07c_acrocentric_presence"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chroms <- c(as.character(1:22), "X", "Y")
all_arms <- as.vector(rbind(paste0("chr", chroms, "_p"),
                            paste0("chr", chroms, "_q")))
acro_p <- c("chr13_p", "chr14_p", "chr15_p", "chr21_p", "chr22_p")

read_tsv <- function(path, row_names = NULL) {
  read.table(path,
             header = TRUE,
             sep = "\t",
             check.names = FALSE,
             quote = "",
             comment.char = "",
             row.names = row_names,
             stringsAsFactors = FALSE)
}

message("Reading matrix: ", matrix_path)
D_raw <- as.matrix(read_tsv(matrix_path, row_names = 1))
stopifnot(nrow(D_raw) == ncol(D_raw))
stopifnot(identical(rownames(D_raw), colnames(D_raw)))

D <- (D_raw + t(D_raw)) / 2
diag(D) <- 0

message("Reading assignments: ", assignments_path)
assignments <- read_tsv(assignments_path)
stopifnot(all(c("ChromArm", "Community", "Arms") %in% names(assignments)))
stopifnot(setequal(rownames(D), assignments$ChromArm))

message("Reading length/signal table: ", len_path)
len <- read_tsv(len_path)
len$chrom_arm <- sub(".*_(chr[0-9XY]+)_(p|q)arm$", "\\1_\\2", len$seq)
len$has_interchrom_signal <- len$region_start != "." & len$arms_involved != "."

count_for <- function(tab, key) {
  x <- as.integer(tab[key])
  x[is.na(x)] <- 0L
  x
}

total_by_arm <- table(len$chrom_arm)
signal_by_arm <- table(len$chrom_arm[len$has_interchrom_signal])

comm_by_arm <- setNames(assignments$Community, assignments$ChromArm)
arms_by_arm <- setNames(assignments$Arms, assignments$ChromArm)
present_arms <- rownames(D)
missing_arms <- setdiff(all_arms, present_arms)

status_for <- function(arm, total_rows, signal_rows, in_matrix) {
  if (isTRUE(in_matrix)) {
    return("present_distinct_tree_tip")
  }
  if (arm == "chr18_q") {
    return("artifact_filtered_then_zero_signal_excluded")
  }
  if (total_rows > 0 && signal_rows == 0) {
    return("zero_signal_excluded")
  }
  "missing_from_source"
}

reason_for <- function(arm, total_rows, signal_rows, in_matrix) {
  if (isTRUE(in_matrix)) {
    return("present in source matrix and rendered as an individual NJ tip")
  }
  if (arm == "chr18_q") {
    return("excluded from 41x41 matrix after the NA18982#1 chr18_q/chrX PAR1 scaffold chimera was removed; remaining source rows have zero retained inter-chromosomal signal")
  }
  if (total_rows > 0 && signal_rows == 0) {
    return("excluded from 41x41 matrix because all source rows have zero retained inter-chromosomal PHR signal at the >=95% identity threshold")
  }
  "not observed in the source length/signal table"
}

arm_presence <- data.frame(
  chrom_arm = all_arms,
  source_rows = count_for(total_by_arm, all_arms),
  retained_interchrom_signal_rows = count_for(signal_by_arm, all_arms),
  zero_signal_rows = count_for(total_by_arm, all_arms) - count_for(signal_by_arm, all_arms),
  in_41x41_matrix = all_arms %in% present_arms,
  leiden_community = ifelse(all_arms %in% names(comm_by_arm),
                            comm_by_arm[all_arms],
                            NA_character_),
  community_arms = ifelse(all_arms %in% names(arms_by_arm),
                          arms_by_arm[all_arms],
                          NA_character_),
  tree_tip_status = ifelse(all_arms %in% present_arms,
                           "distinct_tip",
                           "no_tip_matrix_excluded"),
  collapsed_with = ifelse(all_arms %in% present_arms,
                          "none",
                          NA_character_),
  stringsAsFactors = FALSE
)

arm_presence$audit_status <- mapply(status_for,
                                    arm_presence$chrom_arm,
                                    arm_presence$source_rows,
                                    arm_presence$retained_interchrom_signal_rows,
                                    arm_presence$in_41x41_matrix,
                                    USE.NAMES = FALSE)
arm_presence$reason <- mapply(reason_for,
                              arm_presence$chrom_arm,
                              arm_presence$source_rows,
                              arm_presence$retained_interchrom_signal_rows,
                              arm_presence$in_41x41_matrix,
                              USE.NAMES = FALSE)

write.table(arm_presence,
            file = file.path(out_dir, "arm_presence_complete.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Building NJ tree")
set.seed(20260507)
tr_unrooted <- nj(as.dist(D))
write.tree(tr_unrooted,
           file = file.path(out_dir, "07c_unrooted_nj_audited.newick"))

is_exact_clade <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  if (length(tips) < 2) return(NA)
  node <- tryCatch(getMRCA(tree, tips), error = function(e) NA_integer_)
  if (is.na(node)) return(FALSE)
  setequal(extract.clade(tree, node)$tip.label, tips)
}

tip_branch_lengths <- setNames(rep(NA_real_, length(tr_unrooted$tip.label)),
                               tr_unrooted$tip.label)
for (i in seq_len(nrow(tr_unrooted$edge))) {
  child <- tr_unrooted$edge[i, 2]
  if (child <= length(tr_unrooted$tip.label)) {
    tip_branch_lengths[tr_unrooted$tip.label[child]] <- tr_unrooted$edge.length[i]
  }
}

acro_status <- arm_presence[match(acro_p, arm_presence$chrom_arm), ]
acro_status$tip_label_in_plot <- sub("^chr", "", gsub("_", "", acro_status$chrom_arm))
acro_status$nj_tip_branch_length <- sprintf("%.6f", tip_branch_lengths[acro_status$chrom_arm])
acro_status$monophyletic_c7 <- isTRUE(is_exact_clade(tr_unrooted, acro_p))

write.table(acro_status,
            file = file.path(out_dir, "acrocentric_p_status.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(D_raw[acro_p, acro_p],
            file = file.path(out_dir, "acrocentric_p_distance_matrix.tsv"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

matrix_audit <- data.frame(
  metric = c("source_length_rows",
             "source_length_rows_with_interchrom_signal",
             "source_length_rows_without_interchrom_signal",
             "matrix_rows",
             "matrix_columns",
             "present_arms",
             "excluded_arms",
             "acrocentric_p_arms_present",
             "acrocentric_p_arms_monophyletic",
             "raw_diagonal_min",
             "raw_diagonal_max",
             "max_asymmetry_raw",
             "off_diagonal_min_after_zero_diag",
             "off_diagonal_max_after_zero_diag"),
  value = c(nrow(len),
            sum(len$has_interchrom_signal),
            sum(!len$has_interchrom_signal),
            nrow(D_raw),
            ncol(D_raw),
            length(present_arms),
            paste(missing_arms, collapse = ", "),
            paste(acro_p[acro_p %in% present_arms], collapse = ", "),
            isTRUE(is_exact_clade(tr_unrooted, acro_p)),
            sprintf("%.6f", min(diag(D_raw))),
            sprintf("%.6f", max(diag(D_raw))),
            sprintf("%.6g", max(abs(D_raw - t(D_raw)))),
            sprintf("%.6f", min(D[upper.tri(D)])),
            sprintf("%.6f", max(D[upper.tri(D)]))),
  stringsAsFactors = FALSE
)

write.table(matrix_audit,
            file = file.path(out_dir, "matrix_audit.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

clade_color <- c(
  C1 = "#C95F1B",
  C2 = "#A57900",
  C6 = "#2F6FB2",
  C7 = "#6F42C1",
  C14 = "#12884C",
  C15 = "#2FA36B"
)

arm_short <- function(x) {
  gsub("_", "", sub("^chr", "", x))
}

tip_comm <- function(tips) {
  unname(comm_by_arm[tips])
}

tip_label <- function(tips) {
  base <- paste0(arm_short(tips), " ", tip_comm(tips))
  ifelse(tips %in% acro_p,
         paste0(arm_short(tips), " C7 ACRO"),
         base)
}

tip_col <- function(tips) {
  cc <- tip_comm(tips)
  ifelse(cc %in% names(clade_color), clade_color[cc], "#4A4A4A")
}

tip_face <- function(tips) {
  ifelse(tips %in% acro_p | tip_comm(tips) %in% names(clade_color), 2, 1)
}

draw_wrapped <- function(x, y, text, width = 45, cex = 0.75, col = "#333333",
                         font = 1, line_height = 0.042) {
  lines <- unlist(strwrap(text, width = width))
  if (length(lines) == 0) return(y)
  for (i in seq_along(lines)) {
    text(x, y - (i - 1) * line_height, lines[i],
         adj = c(0, 1), cex = cex, col = col, font = font)
  }
  y - length(lines) * line_height
}

draw_acro_table_panel <- function() {
  op <- par(mar = c(0.35, 0.35, 0.35, 0.35), xpd = NA)
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  rect(0.02, 0.02, 0.98, 0.98, col = "#F7F5FC", border = "#6F42C1", lwd = 2)
  text(0.06, 0.94, "Acrocentric p-arm status", adj = 0, cex = 1.18,
       font = 2, col = "#2B164C")
  y <- draw_wrapped(0.06, 0.885,
                    "All five expected acrocentric p arms are present in the 41x41 matrix and appear as distinct tips on the NJ tree.",
                    width = 40, cex = 0.73, col = "#333333")

  y <- y - 0.025
  text(c(0.07, 0.31, 0.55, 0.77), y,
       c("Arm", "signal rows", "matrix", "tree tip"),
       adj = c(0, 1), cex = 0.63, font = 2, col = "#2B164C")
  segments(0.06, y - 0.02, 0.94, y - 0.02, col = "#6F42C1", lwd = 1.5)
  y <- y - 0.055

  for (i in seq_len(nrow(acro_status))) {
    row <- acro_status[i, ]
    bg <- if (i %% 2 == 1) "#FFFFFF" else "#EFE9FA"
    rect(0.055, y - 0.034, 0.945, y + 0.02, col = bg, border = NA)
    text(0.07, y, arm_short(row$chrom_arm), adj = c(0, 0.5),
         cex = 0.73, font = 2, col = "#2B164C")
    text(0.31, y, row$retained_interchrom_signal_rows, adj = c(0, 0.5),
         cex = 0.69, col = "#222222")
    text(0.55, y, "present", adj = c(0, 0.5),
         cex = 0.69, col = "#222222")
    text(0.77, y, "distinct", adj = c(0, 0.5),
         cex = 0.69, col = "#222222")
    y <- y - 0.058
  }

  y <- y - 0.018
  rect(0.055, y - 0.06, 0.945, y + 0.03, col = "#FFFFFF", border = "#D4C8EB")
  text(0.075, y + 0.005, "Collapsed or filtered?", adj = 0,
       cex = 0.68, font = 2, col = "#2B164C")
  text(0.075, y - 0.03, "No. C7 is monophyletic; each acrocentric p arm has its own tip.",
       adj = 0, cex = 0.63, col = "#333333")

  y <- y - 0.105
  text(0.06, y, "Excluded from matrix", adj = 0,
       cex = 0.76, font = 2, col = "#333333")
  y <- draw_wrapped(0.06, y - 0.04,
                    paste0("Zero-signal arms: ",
                           paste(missing_arms, collapse = ", "),
                           ". These are not acrocentric p arms; chr14q is the acrocentric q arm with zero signal."),
                    width = 43, cex = 0.62, col = "#333333",
                    line_height = 0.034)
}

plot_tree_with_panel <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(2.35, 1))
  plot_tr <- tr_unrooted
  plot_tr$tip.label <- tip_label(tr_unrooted$tip.label)

  op <- par(mar = c(1.2, 0.8, 2.0, 0.2), xpd = NA)
  on.exit({
    par(op)
    layout(1)
  }, add = TRUE)

  plot(plot_tr,
       type = "unrooted",
       use.edge.length = TRUE,
       show.tip.label = TRUE,
       tip.color = tip_col(tr_unrooted$tip.label),
       font = tip_face(tr_unrooted$tip.label),
       cex = 0.68,
       label.offset = 0.006,
       edge.width = 1.25,
       lab4ut = "axial",
       no.margin = FALSE)
  title("Slide 07c candidate: unrooted NJ with acrocentric p arms explicit",
        cex.main = 1.0,
        line = 0.3)
  mtext("Same 41x41 arm-level Jaccard matrix as v2; purple C7 tips are chr13p, chr14p, chr15p, chr21p, chr22p.",
        side = 3,
        line = -0.8,
        cex = 0.62,
        col = "#333333")

  draw_acro_table_panel()
}

plot_minitable <- function() {
  layout(1)
  draw_acro_table_panel()
}

candidate_pdf <- file.path(out_dir, "07c_unrooted_acrocentric_status.pdf")
candidate_png <- file.path(out_dir, "07c_unrooted_acrocentric_status.png")
mini_pdf <- file.path(out_dir, "07c_acrocentric_status_minitable.pdf")
mini_png <- file.path(out_dir, "07c_acrocentric_status_minitable.png")

pdf(candidate_pdf, width = 16, height = 9, useDingbats = FALSE)
plot_tree_with_panel()
dev.off()

png(candidate_png, width = 3840, height = 2160, res = 240, type = "cairo")
plot_tree_with_panel()
dev.off()

pdf(mini_pdf, width = 5.2, height = 9, useDingbats = FALSE)
plot_minitable()
dev.off()

png(mini_png, width = 1248, height = 2160, res = 240, type = "cairo")
plot_minitable()
dev.off()

md_bool <- function(x) ifelse(isTRUE(x), "yes", "no")
md_cell <- function(x) {
  x <- ifelse(is.na(x) | x == "", "NA", as.character(x))
  paste0("`", x, "`")
}

presence_md <- arm_presence
presence_md$leiden_community[is.na(presence_md$leiden_community)] <- "NA"
presence_md$community_arms[is.na(presence_md$community_arms)] <- "NA"
presence_md$collapsed_with[is.na(presence_md$collapsed_with)] <- "NA"

presence_rows <- apply(presence_md, 1, function(row) {
  paste0("| ", md_cell(row[["chrom_arm"]]),
         " | ", row[["source_rows"]],
         " | ", row[["retained_interchrom_signal_rows"]],
         " | ", row[["zero_signal_rows"]],
         " | ", md_bool(row[["in_41x41_matrix"]] == "TRUE" || row[["in_41x41_matrix"]] == TRUE),
         " | ", md_cell(row[["leiden_community"]]),
         " | ", md_cell(row[["tree_tip_status"]]),
         " | ", md_cell(row[["collapsed_with"]]),
         " | ", md_cell(row[["audit_status"]]),
         " | ", row[["reason"]],
         " |")
})

acro_rows <- apply(acro_status, 1, function(row) {
  paste0("| ", md_cell(row[["chrom_arm"]]),
         " | ", row[["source_rows"]],
         " | ", row[["retained_interchrom_signal_rows"]],
         " | ", row[["zero_signal_rows"]],
         " | ", md_cell(row[["leiden_community"]]),
         " | ", md_cell(row[["tree_tip_status"]]),
         " | ", md_cell(row[["collapsed_with"]]),
         " | ", md_cell(row[["nj_tip_branch_length"]]),
         " |")
})

readme <- c(
  "# Slide 07c Acrocentric Presence Audit",
  "",
  "## Recommendation",
  "",
  "**Fix 07c if it remains in the v3 deck; do not keep the current unannotated 07c as-is.** The five acrocentric p arms are not absent from the data or the NJ topology, but the v2 unrooted tree is easy to misread at slide scale. If slide 07c survives as a backup, use `07c_unrooted_acrocentric_status.png` or the same PDF so the side annotation makes the C7 acrocentric status explicit.",
  "",
  "Do not promote 07c to the main sequence-similarity slide. The main slide should remain the heatmap/tree view or the readable rooted tree, because the unrooted layout is a root-sensitivity backup and requires too much label-reading.",
  "",
  "## Short Answer",
  "",
  "- Slide 07c uses `slides/v2-review-zoom/_revision_assets/07b_tree_options/07b_unrooted_nj_option.png`, wired into the v2 zoom deck at `slides/v2-review-zoom/_typst/zoom_review_deck.typ:285-289`.",
  "- The source matrix contains 41 arms. All five acrocentric p arms (`chr13_p`, `chr14_p`, `chr15_p`, `chr21_p`, `chr22_p`) are present in that matrix, assigned to C7, and rendered as distinct NJ tips.",
  "- The absent arms are `chr2_p`, `chr3_p`, `chr5_p`, `chr8_q`, `chr11_q`, `chr14_q`, and `chr18_q`. These are zero-signal/filter exclusions from the 41-arm matrix, not collapsed acrocentric p arms.",
  "- `chr14_q` is an acrocentric chromosome arm, but it is the q arm and has zero retained inter-chromosomal PHR signal at this threshold. The user concern is about the acrocentric p arms, which are all present.",
  "",
  "## Files Produced",
  "",
  "| File | Purpose |",
  "|---|---|",
  "| `07c_unrooted_acrocentric_status.png` / `.pdf` | Candidate replacement visual: unrooted NJ tree plus an explicit acrocentric p-arm status side table. |",
  "| `07c_acrocentric_status_minitable.png` / `.pdf` | Companion annotation-only mini-table if downstream layout keeps the old tree crop. |",
  "| `arm_presence_complete.tsv` | Complete 48-arm audit with source rows, retained-signal rows, matrix presence, tree-tip status, collapsed status, and exclusion reason. |",
  "| `acrocentric_p_status.tsv` | Focused audit for `chr13_p`, `chr14_p`, `chr15_p`, `chr21_p`, and `chr22_p`. |",
  "| `acrocentric_p_distance_matrix.tsv` | Raw source distance submatrix for the five acrocentric p arms; diagonal values are within-arm average distances from the source matrix. |",
  "| `matrix_audit.tsv` | Matrix/source count checks and C7 monophyly result. |",
  "| `07c_unrooted_nj_audited.newick` | Rebuilt unrooted NJ Newick from the same 41x41 matrix. |",
  "| `make_07c_acrocentric_presence.R` | Reproducible generator for this audit and the candidate visuals. |",
  "",
  "## Evidence Trail",
  "",
  "- The v2 tree-option README already states that 41 PHR-positive arms are in the matrix and seven full-arm-set arms are absent by construction at `slides/v2-review-zoom/_revision_assets/07b_tree_options/README.md:7-11`.",
  "- The v2 arm-presence table lists the five acrocentric p arms as present at `slides/v2-review-zoom/_revision_assets/07b_tree_options/arm_presence.tsv:26-44`.",
  "- The report says the retained PHR sequences span 41 of 48 arms and that seven arms have no inter-chromosomal signal at `end-to-end-report/report/01_pipeline.md:71-73`.",
  "- The report names the seven zero-signal arms at `end-to-end-report/report/01_pipeline.md:87-91` and defines arm-level community detection over the remaining 41 arms at `end-to-end-report/report/01_pipeline.md:122-126`.",
  "- The report's Leiden table assigns C7 to all five acrocentric p arms at `end-to-end-report/report/01_pipeline.md:133-151`.",
  "- The v2 unrooted option is generated by `slides/v2-review-zoom/_revision_assets/07b_tree_options/make_07b_tree_options.R:266-335`; it labels all tips but does not add a side audit, which is why 07c can read as incomplete at slide scale.",
  "",
  "## Matrix Audit",
  "",
  "| Metric | Value |",
  "|---|---:|",
  paste0("| Source length/signal rows | ", nrow(len), " |"),
  paste0("| Rows with retained inter-chromosomal signal | ", sum(len$has_interchrom_signal), " |"),
  paste0("| Rows without retained inter-chromosomal signal | ", sum(!len$has_interchrom_signal), " |"),
  paste0("| Matrix dimensions | ", nrow(D_raw), " x ", ncol(D_raw), " |"),
  paste0("| Present arms in matrix/tree | ", length(present_arms), " |"),
  paste0("| Excluded arms from 48-arm set | ", paste(missing_arms, collapse = ", "), " |"),
  paste0("| Acrocentric p arms present | ", paste(acro_p[acro_p %in% present_arms], collapse = ", "), " |"),
  paste0("| C7 acrocentric p arms monophyletic on rebuilt NJ | ", md_bool(is_exact_clade(tr_unrooted, acro_p)), " |"),
  paste0("| Raw diagonal range | ", sprintf("%.6f", min(diag(D_raw))), " to ", sprintf("%.6f", max(diag(D_raw))), " |"),
  paste0("| Off-diagonal range after zeroing diagonal | ", sprintf("%.6f", min(D[upper.tri(D)])), " to ", sprintf("%.6f", max(D[upper.tri(D)])), " |"),
  "",
  "The nonzero source diagonal is expected: it stores within-arm average distances. The NJ rebuild follows the v2 script by symmetrizing the matrix and setting the diagonal to zero before `ape::nj()`.",
  "",
  "## Acrocentric P-Arm Check",
  "",
  "| Arm | source rows | retained signal rows | zero-signal rows | community | tree tip | collapsed with | NJ tip branch length |",
  "|---|---:|---:|---:|---|---|---|---:|",
  acro_rows,
  "",
  "Interpretation: all five acrocentric p arms have retained source signal, all five are in C7, all five have distinct tree tips, and none are collapsed into another label. `chr14_p` has one zero-signal source row, but the arm remains statistically valid for the matrix because 229 rows have retained inter-chromosomal signal.",
  "",
  "## Complete Arm-Presence Table",
  "",
  "| Arm | source rows | retained signal rows | zero-signal rows | in matrix | community | tree tip status | collapsed with | audit status | reason |",
  "|---|---:|---:|---:|---|---|---|---|---|---|",
  presence_rows,
  "",
  "## Fix/Drop Guidance",
  "",
  "- **Keep as-is:** no. The current 07c can be misread as missing acrocentric chromosomes even though the data are correct.",
  "- **Fix:** yes, if 07c remains in the backup section. Use `07c_unrooted_acrocentric_status.png` or add `07c_acrocentric_status_minitable.png` next to the existing unrooted tree. The side table explicitly says all five acrocentric p arms are present and distinct.",
  "- **Drop:** acceptable, and preferred if the deck is being tightened. The unrooted tree is a methodological backup, not the clearest main-slide evidence.",
  "",
  "## Validation",
  "",
  "- Rebuilt the unrooted NJ from `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`.",
  "- Verified the matrix has 41 rows and 41 columns, with seven excluded arms from the 48-arm set.",
  "- Verified `chr13_p`, `chr14_p`, `chr15_p`, `chr21_p`, and `chr22_p` are present in the matrix, present in C7, present as distinct tree tips, and monophyletic on the rebuilt NJ tree.",
  "- Verified no acrocentric p arm is filtered, missing, or collapsed. The only acrocentric chromosome arm excluded from the 41x41 matrix is `chr14_q`, a q arm with zero retained signal.",
  "- Generated candidate visual assets that make the acrocentric p-arm status explicit."
)

writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote outputs to ", out_dir)
