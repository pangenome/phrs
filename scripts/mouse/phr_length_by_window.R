#!/usr/bin/env Rscript
# DEBUG figure (not a manuscript figure).
# Every mouse PHR drawn by its length, annotated by chromosome / arm / genome,
# stratified by flank window (1 / 2 / 4 Mb).  Shows directly how the window caps
# PHR length: a window-independent short set (<100 kb) plus a set that saturates
# the window.  This is the length confound behind raw/slide contact metrics.
#
# One lollipop per PHR (segment 0 -> length, point at length):
#   y = PHR (genome arm), grouped by chromosome then p/q then genome
#   x = PHR length (kb, log scale)
#   colour = chromosome ;  shape = genome (B6 / CAST)
#   facet column = flank window ;  dashed line = window cap
#
# PHR length is a property of the PHR x window only (independent of Hi-C stage /
# resolution), so one seqlevel file per window suffices.  seq id encodes
# genome#hap#arm:coords, e.g. B6#1#chr10_parm:8494-1000000:1-990000 ; size is
# the size_a / size_b column.
#
# INPUT lives in the repo: data/mouse_meiosis_sweep/seqlevel/<window>/*seqlevel.tsv
# Output (OUT_DIR, default /tmp): phr_length_by_window.{png,pdf,tsv}
# Needs ggplot2.  Run: Rscript scripts/mouse/phr_length_by_window.R

suppressPackageStartupMessages(library(ggplot2))

## resolve repo root from this script's own location (scripts/mouse/)
.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- if (length(.self)) normalizePath(dirname(.self)) else getwd()
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))

sweep_dir <- Sys.getenv("SWEEP_DIR",
                        file.path(repo_root, "data/mouse_meiosis_sweep/seqlevel"))
out_dir   <- Sys.getenv("OUT_DIR", "/tmp")
windows   <- c("1Mb", "2Mb", "4Mb")
cap_kb    <- c("1Mb" = 1000, "2Mb" = 2000, "4Mb" = 4000)

## complete PHR set per window = union over ALL stage/resolution seqlevel files
## (a PHR drops out of a single file when its pair has no contact that stage;
## length is intrinsic to the PHR, so any file that contains it gives its size).
phrs_one_window <- function(win) {
  fs <- list.files(file.path(sweep_dir, win),
                   pattern = "seqlevel\\.tsv$", full.names = TRUE)
  if (!length(fs)) stop("no seqlevel file under ", file.path(sweep_dir, win))
  long <- do.call(rbind, lapply(fs, function(f) {
    d <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    rbind(data.frame(seq = d$seq_a, chr = d$chr_a, arm = d$arm_a, size = d$size_a),
          data.frame(seq = d$seq_b, chr = d$chr_b, arm = d$arm_b, size = d$size_b))
  }))
  long <- long[!duplicated(long$seq), ]
  long$genome <- sub("#.*", "", long$seq)                # B6 / CAST
  long$pq     <- sub(".*_([pq])$", "\\1", long$arm)      # p / q
  long$window <- win
  long$len_kb <- long$size / 1000
  long
}

dat <- do.call(rbind, lapply(windows, phrs_one_window))
dat$window <- factor(dat$window, levels = windows)

## optional: keep only PHRs that do NOT saturate the window.  A PHR "saturates"
## when its homology fills the flank, so length just tracks the cap.  Truth test:
## even at the WIDEST window (4 Mb) it stays below SAT_FRAC of the cap, i.e. its
## length is set by real homology extent, not by the window.  FILTER=nonsat keeps
## those; FILTER=all (default) keeps everything.
filter_mode <- Sys.getenv("FILTER", "all")
sat_frac    <- as.numeric(Sys.getenv("SAT_FRAC", "0.95"))
out_tag     <- if (filter_mode == "nonsat") "_nonsat" else ""
dat$phr_key <- paste(dat$arm, dat$genome)                # stable id across windows
if (filter_mode == "nonsat") {
  w4   <- dat[dat$window == "4Mb", ]
  keep <- w4$phr_key[w4$len_kb < sat_frac * cap_kb["4Mb"]]
  dat  <- dat[dat$phr_key %in% keep, ]
  cat(sprintf("FILTER=nonsat: kept %d/49 PHRs (4Mb length < %.0f%% of 4Mb)\n",
              length(keep), 100 * sat_frac))
}

## stable PHR identity / y order: chromosome (numeric), then p before q, then genome
chr_num <- function(c) {
  n <- suppressWarnings(as.integer(sub("^chr", "", c)))
  ifelse(is.na(n), ifelse(sub("^chr", "", c) == "X", 100, 200), n)
}
dat$phr   <- paste0(dat$arm, " ", dat$genome)
dat$chr   <- factor(dat$chr, levels = unique(dat$chr[order(chr_num(dat$chr))]))
ord       <- order(chr_num(as.character(dat$chr)), dat$pq, dat$genome)
phr_levels <- unique(dat$phr[ord])
dat$phr   <- factor(dat$phr, levels = rev(phr_levels))   # top-to-bottom chr1..X

caps <- data.frame(window = factor(windows, levels = windows),
                   cap_kb = cap_kb[windows])

## 20-colour chromosome palette
chr_pal <- grDevices::hcl.colors(nlevels(dat$chr), "Spectral")

p <- ggplot(dat, aes(x = len_kb, y = phr)) +
  geom_vline(data = caps, aes(xintercept = cap_kb),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 5, xend = len_kb, yend = phr, colour = chr),
               linewidth = 0.5) +
  geom_point(aes(colour = chr, shape = genome), size = 2.1) +
  facet_grid(. ~ window) +
  scale_x_log10(breaks = c(10, 30, 100, 300, 1000, 2000, 4000),
                labels = c("10", "30", "100", "300", "1000", "2000", "4000")) +
  scale_colour_manual(values = chr_pal, name = "chromosome", guide = "none") +
  scale_shape_manual(values = c(B6 = 16, CAST = 17), name = "genome") +
  labs(x = "PHR length (kb, log scale)", y = NULL,
       title = if (filter_mode == "nonsat")
                 "Mouse non-saturating PHR length per flank window (DEBUG)"
               else "Mouse PHR length per flank window (DEBUG)",
       subtitle = if (filter_mode == "nonsat")
         paste0("Only PHRs below ", round(100*sat_frac),
                "% of the 4 Mb cap (length set by real homology, not the window); ",
                "dashed line = window cap.")
       else paste(
         "Every PHR (B6 hap1 + CAST hap2) drawn by length; dashed line = window cap.",
         "A window-independent short set (<100 kb) + a set saturating the window.",
         sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, paste0("phr_length_by_window", out_tag, ".png")), p,
       width = 11, height = 9, dpi = 150)
ggsave(file.path(out_dir, paste0("phr_length_by_window", out_tag, ".pdf")), p,
       width = 11, height = 9)

write.table(dat[order(dat$window, chr_num(as.character(dat$chr)), dat$pq, dat$genome),
                c("window", "genome", "chr", "arm", "pq", "len_kb")],
            file.path(out_dir, paste0("phr_length_by_window", out_tag, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

## console summary
cat("PHRs per window:\n")
print(table(dat$window))
cat("\nlength (kb) summary per window:\n")
for (w in windows) {
  v <- dat$len_kb[dat$window == w]
  cat(sprintf("  %-4s n=%2d  median=%.0f  <100kb=%d (%.0f%%)  >=95%%cap=%d (%.0f%%)\n",
              w, length(v), median(v),
              sum(v < 100), 100 * mean(v < 100),
              sum(v >= 0.95 * cap_kb[w]), 100 * mean(v >= 0.95 * cap_kb[w])))
}
cat("wrote ", file.path(out_dir, "phr_length_by_window.{png,pdf,tsv}"), "\n", sep = "")
