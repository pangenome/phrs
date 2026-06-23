#!/usr/bin/env Rscript
# Meiotic prophase-I bouquet schematic for Figure 4 (companion to the mouse
# meiotic Hi-C panel).  Original drawing (base R); four stages matching our mouse
# Hi-C: leptotene -> zygotene -> pachytene -> diplotene.
#
# Literature-grounded cytology (telomeres = red dots = subtelomeric PHRs, on the
# nuclear envelope):
#   leptotene  -- chromosomes condense as thin UNPAIRED single threads; telomeres
#                 attach to the nuclear envelope, dispersed.
#   zygotene   -- BOUQUET: telomeres pulled to one pole by the centrosome/MTOC;
#                 homologues begin to synapse.
#   pachytene  -- full synapsis: homologues paired as bivalents (synaptonemal
#                 complex); bouquet resolved, telomeres dispersed again.
#   diplotene  -- SC disassembles; homologues separate but stay joined at CHIASMATA
#                 (sites of crossover).  Most bivalents have one chiasma (a CROSS),
#                 some two (a RING); telomeres at the bivalent arm tips.
# Refs: prophase-I/SC assembly (PMC8577265); bouquet/centrosome (Science abh3104,
#       PMC84352); chiasmata/bivalent shape + obligate CO, mouse ~20 bivalents,
#       ~1.3 CO/chr (Wikipedia chiasma/bivalent; PubMed 8565707, 17483430).
# Run: Rscript submission/scripts/figures/make_fig4_meiosis_stages.R
# Output (override OUT_DIR=...): meiosis_stages.{png,pdf}

out_dir <- Sys.getenv("OUT_DIR", "/tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

telo_col <- "#d62728"   # chromosome ends / subtelomeric PHRs
homo_col <- "#7d8590"   # homolog threads
chia_col <- "#3a3f45"   # chiasma marker
env_col  <- "#333333"   # nuclear envelope

bez <- function(p0, p1, p2, n = 80) {
  t <- seq(0, 1, length.out = n)
  cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
        (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2])
}
pol <- function(cx, cy, R, a) c(cx + R * cos(a * pi/180), cy + R * sin(a * pi/180))
offset_sides <- function(P, sep) {
  n <- nrow(P)
  d <- rbind(P[2, ] - P[1, ], P[3:n, ] - P[1:(n - 2), ], P[n, ] - P[n - 1, ])
  L <- sqrt(d[, 1]^2 + d[, 2]^2); L[L == 0] <- 1
  nrm <- cbind(-d[, 2]/L, d[, 1]/L)
  list(P + nrm * sep/2, P - nrm * sep/2)
}

## --- compact, locally-placed bivalents (leptotene / pachytene / diplotene) ---
# each bivalent sits over a home angle th on the envelope; telomeres on the rim
draw_local <- function(cx, cy, R, th, style) {
  E <- function(a) pol(cx, cy, R * 0.985, a)          # point on the envelope
  Iin <- function(a, rf) pol(cx, cy, R * rf, a)        # interior point
  tel <- function() {}
  if (style == "lepto") {                              # one unpaired thread
    p0 <- E(th - 9); p2 <- E(th + 9)
    lines(bez(p0, Iin(th, 0.55), p2), col = homo_col, lwd = 1.6)
    return(rbind(p0, p2))
  }
  if (style == "pachy") {                              # synapsed bivalent (paired + SC rungs)
    p0 <- E(th - 9); p2 <- E(th + 9)
    P <- bez(p0, Iin(th, 0.55), p2); s <- offset_sides(P, R * 0.055)
    lines(s[[1]], col = homo_col, lwd = 1.5); lines(s[[2]], col = homo_col, lwd = 1.5)
    ix <- round(seq(7, nrow(P) - 7, length.out = 5))
    segments(s[[1]][ix, 1], s[[1]][ix, 2], s[[2]][ix, 1], s[[2]][ix, 2], col = homo_col, lwd = 0.6)
    return(rbind(p0, p2))
  }
  if (style == "cross") {                              # 1-chiasma bivalent: two homologs cross
    ch <- Iin(th, 0.72)                                # the chiasma
    tA1 <- E(th - 13); tA2 <- E(th + 5); tB1 <- E(th - 5); tB2 <- E(th + 13)
    lines(bez(tA1, ch, tA2), col = homo_col, lwd = 1.5)
    lines(bez(tB1, ch, tB2), col = homo_col, lwd = 1.5)
    points(ch[1], ch[2], pch = 4, col = chia_col, cex = 0.7, lwd = 1.3)   # chiasma X
    return(rbind(tA1, tA2, tB1, tB2))
  }
  if (style == "ring") {                               # 2-chiasma bivalent: a loop
    c1 <- Iin(th - 4, 0.70); c2 <- Iin(th + 4, 0.70)   # two chiasmata
    lines(bez(c1, Iin(th, 0.55), c2), col = homo_col, lwd = 1.5)   # inner arc
    lines(bez(c1, Iin(th, 0.86), c2), col = homo_col, lwd = 1.5)   # outer arc
    tA1 <- E(th - 14); tA2 <- E(th - 7); tB1 <- E(th + 7); tB2 <- E(th + 14)
    segments(c1[1], c1[2], tA1[1], tA1[2], col = homo_col, lwd = 1.5)
    segments(c1[1], c1[2], tA2[1], tA2[2], col = homo_col, lwd = 1.5)
    segments(c2[1], c2[2], tB1[1], tB1[2], col = homo_col, lwd = 1.5)
    segments(c2[1], c2[2], tB2[1], tB2[2], col = homo_col, lwd = 1.5)
    points(rbind(c1, c2), pch = 4, col = chia_col, cex = 0.55, lwd = 1.1)
    return(rbind(tA1, tA2, tB1, tB2))
  }
}

## --- zygotene bouquet (telomeres pulled to one pole; loops hang down) ---
draw_bouquet <- function(cx, cy, R) {
  tel <- matrix(NA, 0, 2)
  bvs <- list(c(70,110,0,-0.60), c(75,105,-0.20,-0.38), c(80,100,0.20,-0.38),
              c(85,95,-0.33,-0.12), c(66,114,0.30,-0.12))
  for (b in bvs) {
    T1 <- pol(cx, cy, R*0.96, b[1]); T2 <- pol(cx, cy, R*0.96, b[2])
    C  <- c(cx + b[3]*R, cy + b[4]*R)
    P  <- bez(T1, C, T2); tt <- seq(0,1,length.out=nrow(P))
    sep <- R*0.15*sin(pi*tt); s <- offset_sides(P, sep)   # paired at the ends, open in the middle
    lines(s[[1]], col=homo_col, lwd=1.5); lines(s[[2]], col=homo_col, lwd=1.5)
    tel <- rbind(tel, T1, T2)
  }
  tel
}

draw_centrosome <- function(x, y, telos) {
  for (t in seq_len(nrow(telos))) segments(x, y, telos[t,1], telos[t,2], col="#c9ccd1", lwd=0.7)
  for (a in seq(0,330,30)) segments(x, y, x+0.42*cos(a*pi/180), y+0.42*sin(a*pi/180), col="#c9ccd1", lwd=0.8)
  points(x, y, pch=21, bg="#555b62", col="white", cex=1.7, lwd=0.7)
}

draw_stage <- function(cx, cy, R, label, kind) {
  highlight <- kind == "zygo"
  th <- seq(0, 2*pi, length.out=220)
  polygon(cx + R*cos(th), cy + R*sin(th), border=env_col, lwd=2, col="#fcfcfc")
  homes <- c(35, 107, 180, 253, 325)
  if (kind == "zygo") {
    telos <- draw_bouquet(cx, cy, R)
    draw_centrosome(cx, cy + R*1.20, telos)
  } else if (kind == "lepto") {
    telos <- do.call(rbind, lapply(homes, function(h) draw_local(cx, cy, R, h, "lepto")))
  } else if (kind == "pachy") {
    telos <- do.call(rbind, lapply(homes, function(h) draw_local(cx, cy, R, h, "pachy")))
  } else if (kind == "diplo") {
    sty <- c("cross","ring","cross","ring","cross")
    telos <- do.call(rbind, Map(function(h,s) draw_local(cx, cy, R, h, s), homes, sty))
  }
  points(telos, pch=21, bg=telo_col, col="white", cex=1.35, lwd=0.55)
  text(cx, cy - R*1.44, label, cex=1.4, font=2, col="#222222")
  if (highlight) {
    text(cx - R*0.62, cy + R*0.93, "bouquet", cex=1.2, font=2, col=telo_col, adj=c(1, 0.5))
    text(cx + R*0.26, cy + R*1.24, "centrosome", cex=1.2, col="#555b62", adj=c(0, 0.5))
  }
}

draw <- function() {
  par(mar=c(0,0,0,0), family="sans")
  plot(NA, xlim=c(1.5,37), ylim=c(0.5,11.9), asp=1, axes=FALSE, ann=FALSE, xaxs="i", yaxs="i")
  R <- 3.4; cy <- 7.2; xs <- c(5, 14.5, 24, 33.5)
  for (i in 1:3) arrows(xs[i]+R+0.4, cy, xs[i+1]-R-0.4, cy, length=0.12, lwd=2, col="#999999")
  draw_stage(xs[1], cy, R, "leptotene", "lepto")
  draw_stage(xs[2], cy, R, "zygotene",  "zygo")
  draw_stage(xs[3], cy, R, "pachytene", "pachy")
  draw_stage(xs[4], cy, R, "diplotene", "diplo")
  legend(x=19.25, y=1.15, xjust=0.5, yjust=0.5, horiz=TRUE, bty="n", pch=c(21,4),
         pt.bg=c(telo_col,NA), col=c("white",chia_col), pt.cex=c(1.5,1.0), pt.lwd=c(0.6,1.3),
         legend=c("telomere (subtelomeric PHRs)","chiasma"), cex=1.1, text.col="#222222")
}

png(file.path(out_dir, "meiosis_stages.png"), width=2400, height=771, res=200, type="cairo")
draw(); dev.off()
pdf(file.path(out_dir, "meiosis_stages.pdf"), width=12, height=3.85); draw(); dev.off()
cat("wrote ", file.path(out_dir, "meiosis_stages.{png,pdf}"), "\n", sep="")
