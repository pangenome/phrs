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

## --- interior-spanning chromosomes (leptotene / pachytene / diplotene) ---
# chromosomes are long threads that cross the nucleus interior between telomeres
# anchored on the envelope (as in the classic prophase-I cartoon).
rimpt <- function(cx, cy, R, a) pol(cx, cy, R * 0.985, a)

# leptotene: many dispersed telomeres; thin UNPAIRED threads criss-cross the centre
draw_lepto <- function(cx, cy, R) {
  set.seed(11)
  na <- 14
  angs <- seq(0, 360, length.out = na + 1)[1:na] + runif(na, -7, 7)
  T <- t(sapply(angs, function(a) rimpt(cx, cy, R, a)))
  ord <- sample(na)
  for (i in seq(1, na - 1, by = 2)) {
    p0 <- T[ord[i], ]; p2 <- T[ord[i + 1], ]
    ctrl <- c(cx + runif(1, -0.5, 0.5) * R, cy + runif(1, -0.5, 0.5) * R)
    lines(bez(p0, ctrl, p2, n = 90), col = homo_col, lwd = 1.25)
  }
  T
}

# pachytene: a few THICK synapsed bivalents spanning the interior (X pattern)
draw_pachy <- function(cx, cy, R) {
  set.seed(21)
  np <- 4
  a1 <- runif(1, 0, 40) + seq(0, 360, length.out = np + 1)[1:np]
  telos <- matrix(NA, 0, 2)
  for (k in seq_len(np)) {
    aa <- a1[k]; bb <- a1[k] + 180 + runif(1, -24, 24)
    p0 <- rimpt(cx, cy, R, aa); p2 <- rimpt(cx, cy, R, bb)
    ctrl <- c(cx + runif(1, -0.32, 0.32) * R, cy + runif(1, -0.32, 0.32) * R)
    lines(bez(p0, ctrl, p2, n = 110), col = homo_col, lwd = 4.3)
    telos <- rbind(telos, p0, p2)
  }
  telos
}

# diplotene: paired homologs spanning the interior, joined at chiasmata (X marks)
draw_diplo <- function(cx, cy, R) {
  set.seed(31)
  np <- 4
  a1 <- runif(1, 0, 40) + seq(0, 360, length.out = np + 1)[1:np]
  telos <- matrix(NA, 0, 2); chi <- matrix(NA, 0, 2)
  for (k in seq_len(np)) {
    aa <- a1[k]; bb <- a1[k] + 180 + runif(1, -22, 22)
    p0 <- rimpt(cx, cy, R, aa); p2 <- rimpt(cx, cy, R, bb)
    # push the control off-centre (own quadrant) so the arcs fan out, not bunch
    ctrl <- c(cx + sample(c(-1, 1), 1) * runif(1, 0.14, 0.40) * R,
              cy + sample(c(-1, 1), 1) * runif(1, 0.14, 0.40) * R)
    P <- bez(p0, ctrl, p2, n = 130); s <- offset_sides(P, R * 0.05)
    lines(s[[1]], col = homo_col, lwd = 2.6); lines(s[[2]], col = homo_col, lwd = 2.6)
    chi <- rbind(chi, P[round(nrow(P) * 0.5), ])
    telos <- rbind(telos, p0, p2)
  }
  points(chi, pch = 4, col = chia_col, cex = 0.85, lwd = 1.6)
  telos
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
  if (kind == "zygo") {
    telos <- draw_bouquet(cx, cy, R)
    draw_centrosome(cx, cy + R*1.20, telos)
  } else if (kind == "lepto") {
    telos <- draw_lepto(cx, cy, R)
  } else if (kind == "pachy") {
    telos <- draw_pachy(cx, cy, R)
  } else if (kind == "diplo") {
    telos <- draw_diplo(cx, cy, R)
  }
  points(telos, pch=21, bg=telo_col, col="white", cex=1.35, lwd=0.55)
  text(cx, cy - R*1.44, label, cex=2.6, font=2, col="#222222")
  if (highlight) {
    text(cx - R*0.62, cy + R*0.93, "bouquet", cex=2.1, font=2, col=telo_col, adj=c(1, 0.5))
    text(cx + R*0.26, cy + R*1.24, "centrosome", cex=2.1, col="#555b62", adj=c(0, 0.5))
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
  legend(x=25.6, y=0.85, xjust=0.5, yjust=0.5, horiz=TRUE, bty="n", pch=c(21,4),
         pt.bg=c(telo_col,NA), col=c("white",chia_col), pt.cex=c(1.5,1.0), pt.lwd=c(0.6,1.3),
         legend=c("telomere (subtelomeric PHRs)","chiasma"), cex=2.0, text.col="#222222")
}

png(file.path(out_dir, "Fig4d_meiosis_stages.png"), width=2400, height=771, res=200, type="cairo")
draw(); dev.off()
pdf(file.path(out_dir, "Fig4d_meiosis_stages.pdf"), width=12, height=3.85); draw(); dev.off()
cat("wrote ", file.path(out_dir, "Fig4d_meiosis_stages.{png,pdf}"), "\n", sep="")
