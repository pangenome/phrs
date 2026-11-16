#!/usr/bin/env Rscript
# SVbyEye synteny: DAD chr9 h2 vs chr1 h1
# args: SYN_DIR QC_DIR OUT_PREFIX
# SYN_DIR from synteny_chr1_chr9.sbatch; writes OUT_PREFIX.<v>.png (whole chromosomes + ends)
suppressPackageStartupMessages({library(SVbyEye); library(GenomicRanges); library(ggplot2); library(patchwork)})
a <- commandArgs(TRUE); d <- a[1]; qcd <- a[2]; out <- a[3]
TRACT <- list(v11 = c(133702042, 133722041), v131 = c(133702516, 133722515))
COPY9 <- list(v11 = c(133682643, 133744561), v131 = c(133683117, 133745035)); COPY1 <- c(8, 61926)

# flagger/nucflag intervals, v1.1 only
qc <- function(ctg, s, e) {
  fs <- c("flagger HiFi" = "flagger.PAN011.hifi.bed", "flagger ONT" = "flagger.PAN011.ONT.bed", "nucflag HiFi" = "nucflag.hifi.PAN011.bed")
  do.call(rbind, lapply(names(fs), function(k) {
    b <- read.table(file.path(qcd, fs[[k]]), sep = "\t", quote = "", fill = TRUE, col.names = paste0("V", 1:9))[, 1:4]
    b <- b[b$V1 == ctg & b$V3 > s & b$V2 < e, ]
    if (nrow(b)) data.frame(start = pmax(b$V2, s), end = pmin(b$V3, e), label = paste0(k, ": ", b$V4), track = k) else NULL
  }))
}
# one annotation row, track name on the y axis, no legend
ann <- function(p, ctg, df, space, level, lab, col) {
  if (is.null(df) || !nrow(df)) return(p)
  gr <- GRanges(ctg, IRanges(df$start, df$end), label = df$label)
  addAnnotation(ggplot.obj = p, annot.gr = gr, coordinate.space = space, shape = "rectangle", fill.by = "label",
                color.palette = setNames(rep(col, length(unique(df$label))), unique(df$label)),
                annotation.level = level, annotation.label = lab, offset.annotation = FALSE)
}
QCC <- c(Dup = "#E41A1C", Col = "#984EA3", Err = "#000000", MISJOIN = "#A65628")
# drop the legends of the annotation tracks
nolegend <- function(p) { for (s in p$scales$scales) { n <- if (is.character(s$name)) s$name else ""; if (any(grepl("^(fill|colour)", s$aesthetics)) && !(n %in% c("Identity", "Alignment\ndirection"))) s$guide <- "none" }; p }
breaks <- c(95, 99, 99.5, 99.9, 99.99)
for (v in c("v11", "v131")) {
  w <- readPaf(file.path(d, paste0("chr1_chr9.whole.", v, ".paf")), include.paf.tags = TRUE, restrict.paf.tags = "cg")
  w <- w[w$aln.len >= 10000, ]
  message(v, " whole: ", nrow(w), " alignments >= 10 kb")
  pw <- plotMiro(paf.table = w, color.by = "identity", perc.identity.breaks = breaks, genomic.scale = "Mbp") +
    ggtitle(paste0("A. Whole chromosomes (p arm at the left), alignments >= 10 kb (", nrow(w), ")"))
  L9w <- w$q.len[1]
  pw <- ann(pw, w$q.name[1], data.frame(start = COPY9[[v]][1], end = L9w, label = "w"), "query", 0.13, "wrong copy (chr9q end, 61.9 kb)", "#7B3294")
  pw <- ann(pw, w$t.name[1], data.frame(start = COPY1[1], end = COPY1[2], label = "c"), "target", 0.13, "true copy (chr1p start, 61.9 kb)", "#7B3294")
  e <- readPaf(file.path(d, paste0("chr1_chr9.", v, ".paf")), include.paf.tags = TRUE, restrict.paf.tags = "cg")
  e <- e[e$aln.len >= 5000, ]
  L <- setNames(read.table(file.path(d, paste0("lengths.", v, ".tsv")))$V2, read.table(file.path(d, paste0("lengths.", v, ".tsv")))$V1)
  w9 <- L[["chr9_window_start"]]; q <- e$q.name[1]; t <- e$t.name[1]
  pe <- plotMiro(paf.table = e, color.by = "identity", perc.identity.breaks = breaks, binsize = 2000, genomic.scale = "kbp") +
    ggtitle(paste0("B. Ends, 2 kb bins: chr9q h2 last 300 kb (", q, ") vs chr1p h1 first 300 kb (", t, ")"))
  r9 <- data.frame(start = c(TRACT[[v]][1], COPY9[[v]][1]) - w9 + 1, end = c(TRACT[[v]][2], COPY9[[v]][2]) - w9 + 1,
                   label = c("9q/3q exchange called in Fig5", "wrong in the DAD assembly: copy of chr1p h1 (61.9 kb)"))
  pe <- ann(pe, q, r9[1, ], "query", 0.4, "9q/3q exchange called in Fig5", "#E07B00")
  pe <- ann(pe, q, r9[2, ], "query", 0.25, "wrong in DAD assembly (copy of chr1p h1)", "#7B3294")
  pe <- ann(pe, t, data.frame(start = COPY1[1], end = COPY1[2], label = "c"), "target", 0.17, "true copy, supported by the reads", "#7B3294")
  if (v == "v11") {
    f9 <- qc("PAN011.chr9.haplotype2", L[["chr9_window_start"]], L[["chr9_end"]]); f1 <- qc("PAN011.chr1.haplotype1", 0, L[["span"]])
    for (k in seq_along(unique(f9$track))) { tk <- unique(f9$track)[k]; x <- f9[f9$track == tk, ]; x$start <- x$start - w9 + 1; x$end <- x$end - w9 + 1
      pe <- ann(pe, q, x, "query", 0.25 + 0.08 * k, paste(unique(x$label), collapse = ", "), QCC[[sub(".*: ", "", x$label[1])]]) }
    for (k in seq_along(unique(f1$track))) { tk <- unique(f1$track)[k]; x <- f1[f1$track == tk, ]
      pe <- ann(pe, t, x, "target", 0.17 + 0.08 * k, paste(unique(x$label), collapse = ", "), QCC[[sub(".*: ", "", x$label[1])]]) }
  } else pe <- pe + labs(caption = "flagger / nucflag: not available for v1.3.1")
  p <- (nolegend(pw) / nolegend(pe)) + plot_layout(guides = "collect") + plot_annotation(title = paste0("DAD (PAN011): chr9 h2 vs chr1 h1 in the DAD assembly, ", ifelse(v == "v11", "v1.1", "v1.3.1"), " (SVbyEye plotMiro, p arm at the left)"))
  ggsave(paste0(out, ".", v, ".png"), p, width = 16, height = 10, dpi = 150, bg = "white")
  message("wrote ", out, ".", v, ".png")
}
