#!/usr/bin/env Rscript
# DEBUG / exploration figure (not a manuscript figure).
# Mouse meiotic-stage trajectory (lepto -> zygo -> pachy -> diplo) of the
# sequence-vs-3D Spearman rho, in a 5 x 3 facet grid:
#   rows = metric (3 contact normalisations x 2 aggregation levels)
#   cols = flank window (1 / 2 / 4 Mb)
# Produced for EACH of the 3 PHR length-class filters (one figure each):
#   all    = every inter-chromosomal PHR pair
#   nonsat = BOTH PHRs below SAT_FRAC of the window cap (length = real homology)
#   sat    = BOTH PHRs at/above SAT_FRAC of the window cap (length = window-capped)
# Splitting by length class shows the zygotene "bouquet peak" lives only in the
# saturating (length-confounded) PHRs; in the non-saturating set zygotene troughs.
#
# Aggregation:
#   per-PHR-pair = each inter-chromosomal PHR pair is a point
#   arm/Mantel   = pairs averaged to one point per chromosome-arm pair
#                  (arm-level Spearman over arm-pair entries = the Mantel rho)
# Contact normalisations:
#   O/E      = size-normalised observed/expected per bin-pair (hic_contact_norm; current Fig 4c)
#   raw      = total balanced contact, not size-normalised (hic_contact_raw)
#   sum-size = hic_contact_raw / (size_a + size_b)  -- the slide's "Sum-method" size norm, per pair
# (The old-figure "slide" arm contact is precomputed arm-level with no per-PHR
#  size, so it cannot be length-class filtered and is omitted from this grid.)
# Each facet overlays the 5 resolutions (colour); ring = each line's peak stage,
# number = n pairs.
#
# INPUT DATA lives in the repo under data/mouse_meiosis_sweep/seqlevel/ :
#   <window>/mouse_<stage>_phr[...]_<res>bp_seqlevel.tsv
# To (re)fetch it from the HPC (run on a host that mounts /moosefs, or via scp):
#   B=/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T
#   for w in 1Mb 2Mb 4Mb; do
#     mkdir -p data/mouse_meiosis_sweep/seqlevel/$w
#     scp "HOST:$B/seqlevel_correlation_$w/mouse_*phr*seqlevel.tsv" data/mouse_meiosis_sweep/seqlevel/$w/
#   done
#
# Outputs (OUT_DIR, default /tmp), one set per filter:
#   mouse_stage_resolution_grid_{all,nonsat,sat}.{png,pdf,tsv}
# Env: FILTER (default: all three) restricts to one of all|nonsat|sat;
#      SAT_FRAC (default 0.95) the saturation threshold.
# Needs ggplot2 + ggrepel.  Run: Rscript scripts/mouse/mouse_stage_resolution_grid.R

suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

## resolve repo root from this script's own location (scripts/mouse/)
.args  <- commandArgs(trailingOnly = FALSE)
.self  <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))

sweep_dir <- Sys.getenv("SWEEP_DIR", file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir   <- Sys.getenv("OUT_DIR", "/tmp")
windows   <- c("1Mb", "2Mb", "4Mb")
cap_bp    <- c("1Mb" = 1e6, "2Mb" = 2e6, "4Mb" = 4e6)
sat_frac  <- as.numeric(Sys.getenv("SAT_FRAC", "0.95"))
filters   <- { f <- Sys.getenv("FILTER", ""); if (f == "") c("all","nonsat","sat") else f }
res_list  <- c(5000, 10000, 20000, 50000, 100000)
stages    <- c("leptotene", "zygotene", "pachytene", "diplotene")
stage_lab <- c("lepto", "zygo", "pachy", "diplo")
reslab_lv <- paste0(res_list / 1000, " kb")
res_cols  <- setNames(c("#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#984ea3"), reslab_lv)
metric_lv <- c("per-PHR-pair (O/E)", "per-PHR-pair (raw)", "per-PHR-pair (sum-size)",
               "arm/Mantel (O/E)", "arm/Mantel (raw)")
filt_title <- c(
  all    = "ALL PHR pairs",
  nonsat = sprintf("NON-SATURATING PHRs only (both PHRs < %d%% of window cap)", round(100*sat_frac)),
  sat    = sprintf("SATURATING PHRs only (both PHRs >= %d%% of window cap)", round(100*sat_frac)))

sp <- function(a, b) suppressWarnings(cor(rank(a), rank(b)))

find_file <- function(win, stage, res) {
  pat <- sprintf("^mouse_%s_.*_%sbp_seqlevel\\.tsv$", stage, formatC(res, format = "d"))
  f <- list.files(file.path(sweep_dir, win), pattern = pat, full.names = TRUE)
  if (length(f) != 1) NA_character_ else f
}

compute <- function(f, cap, fm) {
  na <- c(pair_norm=NA, pair_raw=NA, pair_sumsize=NA, arm_norm=NA, arm_raw=NA, pair_n=NA, arm_n=NA)
  if (is.na(f)) return(na)
  d <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$chr_a != d$chr_b & !is.na(d$jaccard) &
         !is.na(d$hic_contact_norm) & !is.na(d$hic_contact_raw), ]
  if (fm == "nonsat") d <- d[d$size_a <  sat_frac * cap & d$size_b <  sat_frac * cap, ]
  if (fm == "sat")    d <- d[d$size_a >= sat_frac * cap & d$size_b >= sat_frac * cap, ]
  if (!nrow(d)) return(na)
  d$sumsize <- d$hic_contact_raw / (d$size_a + d$size_b)
  key <- apply(cbind(pmin(d$arm_a, d$arm_b), pmax(d$arm_a, d$arm_b)), 1, paste, collapse = "|")
  ag <- aggregate(cbind(jaccard, hic_contact_norm, hic_contact_raw) ~ key,
                  data = cbind(d, key), FUN = mean)
  c(pair_norm    = sp(d$jaccard, d$hic_contact_norm),
    pair_raw     = sp(d$jaccard, d$hic_contact_raw),
    pair_sumsize = sp(d$jaccard, d$sumsize),
    arm_norm     = sp(ag$jaccard, ag$hic_contact_norm),
    arm_raw      = sp(ag$jaccard, ag$hic_contact_raw),
    pair_n = nrow(d), arm_n = nrow(ag))
}

make_figure <- function(fm) {
  ## ---- gather ----
  rows <- list()
  for (w in windows) for (r in res_list) {
    vals <- sapply(stages, function(s) compute(find_file(w, s, r), cap_bp[w], fm))
    for (i in seq_along(stages))
      rows[[length(rows) + 1]] <- data.frame(window = w, res = r, stage = stages[i],
        t(vals[, i]), stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(tab, file.path(out_dir, sprintf("mouse_stage_resolution_grid_%s.tsv", fm)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## ---- long format ----
  mk <- function(lbl, rho, n) data.frame(window = tab$window, res = tab$res,
          stage = tab$stage, metric = lbl, rho = rho, n = n)
  long <- rbind(
    mk(metric_lv[1], tab$pair_norm,    tab$pair_n),
    mk(metric_lv[2], tab$pair_raw,     tab$pair_n),
    mk(metric_lv[3], tab$pair_sumsize, tab$pair_n),
    mk(metric_lv[4], tab$arm_norm,     tab$arm_n),
    mk(metric_lv[5], tab$arm_raw,      tab$arm_n))
  long$reslab <- factor(paste0(long$res / 1000, " kb"), levels = reslab_lv)
  long$x      <- as.integer(factor(long$stage, levels = stages))
  long$window <- factor(long$window, levels = windows)
  long$metric <- factor(long$metric, levels = metric_lv)
  long <- long[!is.na(long$rho), ]

  peaks <- do.call(rbind, by(long, list(long$metric, long$window, long$reslab),
                             function(d) if (nrow(d)) d[which.max(d$rho), ] else NULL))

  p <- ggplot(long, aes(x, rho, colour = reslab, group = reslab)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0) +
    geom_point(data = peaks, shape = 21, fill = NA, size = 4.6, stroke = 1.2) +
    geom_text_repel(aes(label = n), size = 3.4, fontface = "bold",
                    show.legend = FALSE, max.overlaps = Inf,
                    min.segment.length = 0, segment.size = 0.2,
                    box.padding = 0.22, point.padding = 0.18, seed = 1) +
    facet_grid(metric ~ window, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(stages), labels = stage_lab,
                       expand = expansion(mult = 0.08)) +
    scale_y_continuous(breaks = seq(-0.3, 0.95, 0.05), minor_breaks = NULL) +
    scale_colour_manual(values = res_cols, name = "resolution") +
    labs(x = NULL, y = expression(rho ~ "(sequence vs 3D)"),
         title = sprintf("Mouse stage-trajectory -- %s (DEBUG)", filt_title[fm]),
         subtitle = paste(
           "Spearman rho of PHR sequence similarity vs mouse meiotic Hi-C contact, per prophase stage (zygo = bouquet).  Colour = Hi-C bin resolution.",
           "Level -- per-PHR-pair: each pair a point;  arm/Mantel: averaged per arm pair = the Mantel matrix-correlation rho.",
           "Contact -- O/E: size-normalised observed/expected per bin-pair (Fig 4c);  raw: total balanced;  sum-size: raw/(size_a+size_b).  Cols (1/2/4 Mb): flank window;  ring = peak;  n = pairs.",
           sep = "\n")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          legend.key.width = unit(1.5, "cm"),
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 11),
          panel.grid.minor.x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2.4, size = 3.4)))

  ggsave(file.path(out_dir, sprintf("mouse_stage_resolution_grid_%s.png", fm)), p,
         width = 15, height = 15, dpi = 150, limitsize = FALSE)
  ggsave(file.path(out_dir, sprintf("mouse_stage_resolution_grid_%s.pdf", fm)), p,
         width = 15, height = 15, limitsize = FALSE)

  ## ---- console summary: per-PHR-pair O/E peak stage ----
  cat(sprintf("[%s] per-PHR-pair (O/E) peak stage per window x resolution (* = zygotene):\n", fm))
  for (w in windows) for (r in res_list) {
    v <- long[long$metric == "per-PHR-pair (O/E)" & long$window == w & long$res == r, ]
    v <- v[match(stages, v$stage), ]
    pk <- stages[which.max(v$rho)]
    cat(sprintf("    %-4s %-6s %-10s%s\n", w, paste0(r/1000,"kb"), pk,
                ifelse(length(pk) && pk == "zygotene", " *", "")))
  }
  invisible(NULL)
}

for (fm in filters) make_figure(fm)
cat("wrote ", file.path(out_dir, "mouse_stage_resolution_grid_{all,nonsat,sat}.{png,pdf,tsv}"), "\n", sep = "")
