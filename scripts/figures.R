#!/usr/bin/env Rscript
# ============================================================
#  figures.R  --  render the main figures (base R, no external packages)
#  EJHG layout: max 4 main figures, each with at most 3 panels.
#    Fig. 1  genome-wide architecture + chr6 region        (2 panels)
#    Fig. 2  fine-mapping and candidate genes              (3 panels)
#    Fig. 3  credible-set effects across traits            (3 panels)
#    Fig. 4  no MHC signal + skeletal phenome profile      (3 panels)
#    Fig. S1 genetic relatedness (supplementary)           (2 panels)
#  input : figures/data/*.tsv, work/fig_manhattan.tsv, work/fig_region.tsv
#  output: figures/Fig*.png (300 dpi) + figures/Fig*.pdf (vector)
#  usage : Rscript figures.R
# ============================================================
D   <- "."
FD  <- file.path(D, "figures", "data")
OUT <- file.path(D, "figures")
WORK <- file.path(D, "work")

COL <- c(target = "#C0392B", autoimmune = "#2166AC", entrapment = "#E08214",
         neuro = "#7FB3D5", metabolic = "#4D9221", skeletal = "#762A83",
         other = "#7F7F7F")

SHT <- c("Bell's palsy" = "Bell's palsy",
         "Nerve, nerve root and plexus disorders" = "Nerve / plexus disorders",
         "Carpal tunnel syndrome" = "Carpal tunnel syndrome",
         "Epilepsy" = "Epilepsy",
         "Autoimmune diseases, non-thyroid" = "Autoimmune, non-thyroid",
         "Autoimmune hypothyroidism" = "Autoimmune hypothyroidism",
         "Type 1 diabetes" = "Type 1 diabetes",
         "Rheumatoid arthritis" = "Rheumatoid arthritis",
         "Psoriasis" = "Psoriasis",
         "Multiple sclerosis" = "Multiple sclerosis",
         "Type 2 diabetes" = "Type 2 diabetes")
short <- function(x) ifelse(x %in% names(SHT), SHT[x], x)

RENDER_PDF <- FALSE     # TRUE while drawing the vector version

render <- function(base, w, h, fn) {
  RENDER_PDF <<- FALSE
  png(file.path(OUT, paste0(base, ".png")), width = w, height = h, units = "in",
      res = 300, type = "cairo", bg = "white")
  fn(); dev.off()
  RENDER_PDF <<- TRUE
  pdf(file.path(OUT, paste0(base, ".pdf")), width = w, height = h, bg = "white")
  fn(); dev.off()
  RENDER_PDF <<- FALSE
  cat("  ", base, "\n")
}

rd <- function(dir, f, ...) read.delim(file.path(dir, f), header = TRUE,
                                       stringsAsFactors = FALSE, check.names = FALSE, ...)

# panel letter, correct also on log axes
panel_letter <- function(lab, cex = 1.0, dx = 0.012, dy = 0.985) {
  u <- par("usr")
  x <- if (isTRUE(par("xlog"))) 10^(u[1] + dx * (u[2] - u[1])) else u[1] + dx * (u[2] - u[1])
  y <- u[3] + dy * (u[4] - u[3])
  text(x, y, lab, adj = c(0, 1), cex = cex, font = 2, xpd = NA)
}

# gene track: boxes + optional hanging gene names
draw_genes <- function(g, xlimMb, y0, h, lab.cex = 0.42, show.lab = TRUE,
                       ylab.gap = 0.05) {
  s <- g$start / 1e6; e <- g$end / 1e6
  keep <- which(e >= xlimMb[1] & s <= xlimMb[2])
  s <- pmax(s, xlimMb[1]); e <- pmin(e, xlimMb[2])
  for (i in keep) {
    rect(s[i], y0 - h / 2, e[i], y0 + h / 2, col = "#BDBDBD",
         border = "#4D4D4D", lwd = 0.35)
    if (g$strand[i] == "+") segments(e[i], y0 - h / 2, e[i], y0 + h / 2, lwd = 0.9)
    else                    segments(s[i], y0 - h / 2, s[i], y0 + h / 2, lwd = 0.9)
    if (show.lab)
      text((s[i] + e[i]) / 2, y0 - h / 2 - ylab.gap, g$symbol[i], srt = 90,
           adj = c(1, 0.5), cex = lab.cex, col = "#1A1A1A", xpd = NA)
  }
  invisible(NULL)
}

# ============================================================
# Fig. 1 | genome-wide architecture and the chr6 loci
# ============================================================
fig1 <- function() {
  m <- rd(WORK, "fig_manhattan.tsv",
          colClasses = c("character", "integer", "numeric"))
  ch <- unique(m$chrom)
  num <- ch[grepl("^[0-9]+$", ch)]; num <- num[order(as.integer(num))]
  oth <- setdiff(ch, num)
  oth <- oth[order(match(oth, c("X", "Y", "MT", "XY")), oth)]
  ord <- c(num, oth)
  lab <- ifelse(ord == "23", "X", ord)          # FinnGen encodes X as 23
  len <- vapply(ord, function(cc) max(m$pos[m$chrom == cc]), numeric(1))
  off <- c(0, cumsum(len))[seq_along(ord)]
  names(off) <- ord
  m$g <- off[m$chrom] + m$pos
  tot <- sum(len)
  ymax <- max(m$mlogp, na.rm = TRUE)

  layout(matrix(1:2, nrow = 2), heights = c(1.3, 1))

  ## --- a: Manhattan ---
  par(mar = c(2.4, 5.5, 1.8, 0.7), las = 1, mgp = c(2.6, 0.55, 0), xpd = FALSE)
  ci <- match(m$chrom, ord)
  colv <- ifelse(ci %% 2 == 0, "#BDBDBD", "#4D4D4D")
  sig <- !is.na(m$mlogp) & m$mlogp > 7.3013
  # 1.09 M points as vectors makes a 43 MB PDF; 6-fold thinning is invisible at 300 dpi
# (keeps the vector file under the 10 MB limit most submission systems impose)
  gi <- which(!sig)
  if (RENDER_PDF) gi <- gi[seq(1, length(gi), by = 6)]
  plot(m$g[gi], m$mlogp[gi], pch = 16, cex = 0.16, col = colv[gi],
       xlim = c(0, tot), ylim = c(0, ymax * 1.22), xaxt = "n", xlab = "", ylab = "",
       bty = "l", cex.axis = 0.62)
  points(m$g[sig], m$mlogp[sig], pch = 16, cex = 0.55, col = COL[["target"]])
  abline(h = 7.3013, lty = 2, col = "#C0392B", lwd = 0.7)
  axis(1, at = off + len / 2, labels = lab, tick = FALSE, line = -0.75, cex.axis = 0.55)
  mtext(expression(-log[10](italic(P))), 2, line = 3.4, las = 0, cex = 0.75)
  text(tot * 0.995, 7.3013 + ymax * 0.03, expression(italic(P) == 5 %*% 10^-8),
       cex = 0.5, col = "#C0392B", adj = 1)

  x6 <- off[["6"]]
  LOC <- list(list(p = 44479182, t = "chr6:44.48 Mb, intergenic", y = 1.19),
              list(p = 45602126, t = "chr6:45.60 Mb, RUNX2", y = 1.08))
  for (k in LOC) {
    gx <- x6 + k$p
    pk <- max(m$mlogp[abs(m$g - gx) < 3e5], na.rm = TRUE)
    segments(gx, pk + 0.5, gx, ymax * k$y, lwd = 0.6, col = "#C0392B")
    text(gx, ymax * k$y, k$t, cex = 0.5, adj = 0, col = "#C0392B")
  }

  ## --- b: chr6 region ---
  par(mar = c(3.6, 5.5, 0.6, 0.7))
  r <- rd(WORK, "fig_region.tsv")
  r <- r[r$pos >= 44.0e6 & r$pos <= 45.0e6, ]
  r$mlogp <- -log10(pmax(r$pval, 1e-320))
  g <- rd(FD, "fig_genes.tsv")
  cs <- rd(FD, "fig_crediblesets.tsv")
  y2 <- max(r$mlogp, na.rm = TRUE)
  plot(r$pos / 1e6, r$mlogp, pch = 16, cex = 0.34, col = "#4D4D4D",
       xlim = c(44.0, 45.0), ylim = c(-5.6, y2 * 1.10), xaxt = "n", xlab = "",
       ylab = "", bty = "l", cex.axis = 0.62)
  axis(1, at = seq(44.0, 45.0, 0.2), labels = sprintf("%.1f", seq(44.0, 45.0, 0.2)),
       cex.axis = 0.62)
  mtext("Position on chromosome 6 (Mb)", 1, line = 2.3, cex = 0.72)
  mtext(expression(-log[10](italic(P))), 2, line = 3.4, las = 0, cex = 0.72)
  abline(h = 7.3013, lty = 2, col = "#C0392B", lwd = 0.7)
  abline(h = 0, col = "grey70", lwd = 0.4)
  c1 <- cs[cs$cs == "CS1", ]
  points(c1$pos / 1e6, r$mlogp[match(c1$pos, r$pos)], pch = 18, cex = 1.2,
         col = COL[["target"]])
  segments(min(c1$pos) / 1e6 - 0.005, y2 * 0.80, max(c1$pos) / 1e6 + 0.005, y2 * 0.80,
           col = COL[["target"]], lwd = 0.8)
  text(mean(c1$pos) / 1e6, y2 * 0.83, "four-variant credible set (CS1)",
       cex = 0.52, col = COL[["target"]])
  draw_genes(g, c(44.0, 45.0), y0 = -2.1, h = 0.45, show.lab = FALSE)
  gs <- g[g$symbol %in% c("MRPL14", "CAPN11", "SLC35B2", "CDC5L", "SUPT3H"), ]
  text((pmax(gs$start / 1e6, 44.0) + pmin(gs$end / 1e6, 45.0)) / 2, -2.55,
       gs$symbol, srt = 90, adj = c(1, 0.5), cex = 0.5, xpd = NA)
  mtext("16 protein-coding genes", 3, at = 44.0, adj = 0, line = -1.15, cex = 0.52,
        col = "grey30")
}

# ============================================================
# Fig. 2 | fine-mapping and candidate genes
# ============================================================
fig2 <- function() {
  r <- rd(WORK, "fig_region.tsv")
  r$mlogp <- -log10(pmax(r$pval, 1e-320))
  cs <- rd(FD, "fig_crediblesets.tsv")
  g <- rd(FD, "fig_genes.tsv")
  f3 <- rd(FD, "fig4_credibleset.tsv")     # data table shared with Fig. 3
  ymax <- max(r$mlogp, na.rm = TRUE)

  layout(matrix(1:4, nrow = 4), heights = c(2.30, 1.70, 1.42, 1.28))

  ## --- a: the three credible sets ---
  par(mar = c(0.7, 5.5, 1.5, 0.7), las = 1, mgp = c(2.6, 0.55, 0))
  r2 <- r[r$pos >= 42.9e6 & r$pos <= 46.0e6, ]
  plot(r2$pos / 1e6, r2$mlogp, pch = 16, cex = 0.3, col = "#BDBDBD",
       xlim = c(42.9, 46.0), ylim = c(0, ymax * 1.22), xaxt = "n", xlab = "",
       ylab = "", bty = "l", cex.axis = 0.62)
  abline(h = 7.3013, lty = 2, col = "#C0392B", lwd = 0.7)
  cs_col <- c(CS1 = COL[["target"]], CS2 = "#2166AC", CS3 = "#7F7F7F")
  for (k in names(cs_col)) {
    p <- cs$pos[cs$cs == k] / 1e6
    yv <- r$mlogp[match(cs$pos[cs$cs == k], r$pos)]
    points(p, yv, pch = 18, cex = 1.0, col = cs_col[[k]])
  }
  mtext(expression(-log[10](italic(P))), 2, line = 3.4, las = 0, cex = 0.72)
  legend("topleft", legend = c("Credible set 1 (4 variants, log10BF 18.0)",
                              "Credible set 2 (27 variants, RUNX2)",
                              "Credible set 3 (40 variants, low purity)"),
         pch = 18, col = unname(cs_col), bty = "n", cex = 0.48, inset = 0.005,
         x.intersp = 0.7, y.intersp = 1.1)
  legend("topright", legend = expression(italic(P) == 5 %*% 10^-8), lty = 2,
         col = "#C0392B", bty = "n", cex = 0.48)

  ## --- b: effect sizes at the five variants ---
  par(mar = c(3.2, 5.5, 1.5, 0.7))
  VAR <- c("6:44478351", "6:44478937", "6:44479182", "6:44479267", "6:45602126")
  b <- f3[f3$phenocode == "G6_BELLPA", ]
  b <- b[match(VAR, b$variant), ]
  xs2 <- 1:5
  plot(NA, xlim = c(0.4, 5.6), ylim = c(0, max(b$beta + 1.96 * b$se) * 1.30),
       xaxt = "n", xlab = "", ylab = "", bty = "l", cex.axis = 0.62)
  segments(xs2, b$beta - 1.96 * b$se, xs2, b$beta + 1.96 * b$se, lwd = 1.1,
           col = c(rep(COL[["target"]], 4), "#2166AC"))
  points(xs2, b$beta, pch = 19, cex = 1.0, col = c(rep(COL[["target"]], 4), "#2166AC"))
  axis(1, at = xs2, labels = sprintf("%.4f", c(44.478351, 44.478937, 44.479182,
                                               44.479267, 45.602126)), cex.axis = 0.58)
  mtext("Position on chromosome 6 (Mb)", 1, line = 2.0, cex = 0.68)
  mtext("Effect size (beta)", 2, line = 3.4, las = 0, cex = 0.68)
  mtext("Credible set 1", 3, at = 2.5, line = 0.2, cex = 0.58, col = COL[["target"]])
  mtext("CS2 (RUNX2)", 3, at = 5, line = 0.2, cex = 0.58, col = "#2166AC")

  ## --- c: blood cis-eQTL Z scores ---
  par(mar = c(0.4, 5.5, 1.9, 0.7))
  ge <- g[!is.na(g$eqtl_z), ]
  ge$z <- as.numeric(ge$eqtl_z)
  pp <- ge$eqtl_pos / 1e6
  zmax <- max(abs(ge$z)) * 1.42
  plot(NA, xlim = c(44.11, 45.38), ylim = c(0, zmax), xaxt = "n", xlab = "",
       ylab = "", bty = "l", cex.axis = 0.62)
  segments(pp, 0, pp, abs(ge$z), lwd = 2.6,
           col = ifelse(ge$z > 0, COL[["target"]], "#2166AC"))
  points(pp, abs(ge$z), pch = 19, cex = 0.5,
         col = ifelse(ge$z > 0, COL[["target"]], "#2166AC"))
  text(pp, abs(ge$z) + zmax * 0.02, ge$symbol, srt = 90, adj = c(0, 0.5),
       cex = 0.46, col = "grey15", xpd = NA)
  abline(h = 0, col = "grey60", lwd = 0.5)
  mtext("Blood cis-eQTL |Z| (eQTLGen)", 2, line = 3.4, las = 0, cex = 0.68)
  legend("topleft", legend = c("Z > 0", "Z < 0"), pch = 19,
         col = c(COL[["target"]], "#2166AC"), bty = "n", cex = 0.46)

  par(mar = c(3.2, 5.5, 1.0, 0.7))
  plot(NA, xlim = c(44.11, 45.38), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n")
  axis(1, at = seq(44.2, 45.4, 0.2), labels = sprintf("%.1f", seq(44.2, 45.4, 0.2)),
       cex.axis = 0.6)
  mtext("Position on chromosome 6 (Mb)", 1, line = 2.1, cex = 0.68)
  draw_genes(g, c(44.11, 45.38), y0 = 0.50, h = 0.26, lab.cex = 0.42, ylab.gap = 0.04)
  mtext("Protein-coding genes", 3, at = 44.11, adj = 0, line = 0.15, cex = 0.55,
        col = "grey30")
}

# ============================================================
# Fig. 3 | effects at the credible-set variants across traits
# ============================================================
fig3 <- function() {
  f <- rd(FD, "fig4_credibleset.tsv")
  tr <- c("Bell's palsy", "Nerve, nerve root and plexus disorders",
          "Carpal tunnel syndrome", "Epilepsy", "Autoimmune diseases, non-thyroid",
          "Autoimmune hypothyroidism", "Type 1 diabetes", "Rheumatoid arthritis",
          "Psoriasis", "Multiple sclerosis", "Type 2 diabetes")
  VAR <- c("6:44478351", "6:44478937", "6:44479182", "6:44479267")
  OFF <- c(-0.27, -0.09, 0.09, 0.27)
  VCOL <- c("#C0392B", "#E08214", "#7B3294", "#1A9850")
  layout(matrix(1:3, nrow = 3), heights = c(2.1, 1.9, 1.75))

  ## --- a: effect sizes ---
  par(mar = c(3.4, 8.6, 1.4, 0.7), las = 1, mgp = c(2.6, 0.55, 0))
  s1 <- f[f$cs == "CS1", ]
  xr <- range(c(s1$beta - 1.96 * s1$se, s1$beta + 1.96 * s1$se))
  plot(NA, xlim = xr, ylim = c(0.4, length(tr) + 0.6), yaxt = "n", xlab = "",
       ylab = "", bty = "l", cex.axis = 0.62)
  abline(v = 0, lty = 3, col = "grey45")
  for (i in seq_along(tr)) for (j in seq_along(VAR)) {
    k <- which(s1$trait == tr[i] & s1$variant == VAR[j])
    if (!length(k)) next
    y <- i + OFF[j]
    segments(s1$beta[k] - 1.96 * s1$se[k], y, s1$beta[k] + 1.96 * s1$se[k], y,
             col = VCOL[j], lwd = 1.0)
    points(s1$beta[k], y, pch = 19, cex = 0.65, col = VCOL[j])
  }
  axis(2, at = seq_along(tr), labels = short(tr), las = 1, cex.axis = 0.58, tick = FALSE)
  mtext("Effect size at credible set 1 (beta)", 1, line = 2.2, cex = 0.7)
  legend("topright", legend = sub("^6:", "chr6:", VAR), pch = 19, col = VCOL,
         bty = "n", cex = 0.52, ncol = 2, x.intersp = 0.6, y.intersp = 1.05)
  panel_letter("a")

  ## --- b: P values ---
  par(mar = c(3.4, 8.6, 0.5, 0.7))
  plot(NA, xlim = c(min(s1$p) * 0.4, 1.6), ylim = c(0.4, length(tr) + 0.6),
       yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "l", log = "x",
       cex.axis = 0.62)
  abline(v = 5e-8, lty = 2, col = "#C0392B", lwd = 0.7)
  abline(v = 0.05, lty = 3, col = "grey45", lwd = 0.7)
  for (i in seq_along(tr)) for (j in seq_along(VAR)) {
    k <- which(s1$trait == tr[i] & s1$variant == VAR[j])
    if (!length(k)) next
    points(s1$p[k], i + OFF[j], pch = 19, cex = 0.6, col = VCOL[j])
  }
  axis(1, at = 10^-c(0, 4, 8, 12, 16, 20, 24),
       labels = c("1", "1e-4", "1e-8", "1e-12", "1e-16", "1e-20", "1e-24"),
       cex.axis = 0.6)
  axis(2, at = seq_along(tr), labels = short(tr), las = 1, cex.axis = 0.58, tick = FALSE)
  mtext("P value", 1, line = 2.2, cex = 0.7)
  text(5e-8, length(tr) + 0.42, "5 x 10^-8", cex = 0.48, col = "#C0392B", adj = 0.5)
  text(0.05, length(tr) + 0.42, "0.05", cex = 0.48, col = "grey35", adj = 0)
  panel_letter("b")

  ## --- c: RUNX2 credible set ---
  par(mar = c(3.4, 8.6, 0.5, 0.7))
  s2 <- f[f$cs == "CS2", ]
  xr2 <- range(c(s2$beta - 1.96 * s2$se, s2$beta + 1.96 * s2$se))
  plot(NA, xlim = xr2, ylim = c(0.4, length(tr) + 0.6), yaxt = "n", xlab = "",
       ylab = "", bty = "l", cex.axis = 0.62)
  abline(v = 0, lty = 3, col = "grey45")
  for (i in seq_along(tr)) {
    k <- which(s2$trait == tr[i])
    if (!length(k)) next
    segments(s2$beta[k] - 1.96 * s2$se[k], i, s2$beta[k] + 1.96 * s2$se[k], i,
             col = "#2166AC", lwd = 1.0)
    points(s2$beta[k], i, pch = 19, cex = 0.65, col = "#2166AC")
  }
  axis(2, at = seq_along(tr), labels = short(tr), las = 1, cex.axis = 0.58, tick = FALSE)
  mtext("Effect size at RUNX2 credible set (beta)", 1, line = 2.2, cex = 0.7)
  panel_letter("c")
}

# ============================================================
# shared bar panel for Fig. 4a and 4b (identical category order)
# ============================================================
bar_traits <- function(d, val, ylab, fmt, ytop, show.labels = TRUE) {
  par(mar = c(if (show.labels) 8.6 else 1.0, 4.6, 1.5, 0.8), las = 1,
      mgp = c(3.0, 0.6, 0))
  h <- d[[val]]
  bp <- barplot(h, names.arg = rep("", nrow(d)), col = COL[d$category],
                border = "#333333", lwd = 0.4, ylim = c(0, ytop), ylab = "",
                cex.axis = 0.65, xpd = FALSE)
  mtext(ylab, 2, line = 2.7, las = 0, cex = 0.7)
  text(bp, ifelse(h > ytop * 0.14, h - ytop * 0.055, h + ytop * 0.035), fmt(h),
       cex = 0.42, col = ifelse(h > ytop * 0.14, "white", "#1A1A1A"))
  if (show.labels) {
    lab <- paste0(short(d$trait), " (", format(d$n_cases, big.mark = ","), ")")
    text(bp, -ytop * 0.035, lab, srt = 90, adj = c(1, 0.5), cex = 0.5, xpd = NA)
  }
  abline(h = 0, col = "grey30", lwd = 0.6)
  invisible(bp)
}

# ============================================================
# Fig. 4 | no MHC signal, and a skeletal phenome profile
# ============================================================
fig4 <- function() {
  a <- rd(FD, "fig3a_mhc.tsv")
  b <- rd(FD, "fig3b_hla.tsv")
  p <- rd(FD, "fig5_phenome.tsv"); p <- p[order(p$rank), ]

  layout(matrix(1:3, nrow = 3), heights = c(1.55, 1.75, 2.55))

  ## --- a: MHC proportion (shares the trait order with b) ---
  bar_traits(a, "pct_mhc", "Significant variants in MHC (%)",
             function(h) sprintf("%d/%d", a$n_mhc, a$n_gws), 100, show.labels = FALSE)
  panel_letter("a")

  ## --- b: HLA alleles ---
  bar_traits(b, "n_lt_5e8", "HLA alleles, P < 5 x 10^-8", function(h) h, 110)
  panel_letter("b")

  ## --- c: phenome-wide scan ---
  par(mar = c(4.4, 14.0, 1.8, 0.9), las = 1, mgp = c(2.9, 0.6, 0))
  bp <- barplot(rev(p$mlogp), horiz = TRUE, names.arg = rep("", nrow(p)),
                col = COL[rev(p$kind)], border = "#333333", lwd = 0.4,
                xlim = c(0, 27.5), xlab = "", cex.axis = 0.65)
  axis(2, at = bp, labels = rev(p$trait), las = 1, cex.axis = 0.6, tick = FALSE)
  abline(v = 7.3013, lty = 2, col = "#C0392B", lwd = 0.8)
  abline(v = 1.3010, lty = 3, col = "grey45", lwd = 0.7)
  text(7.3013 + 0.3, bp[1] + 0.35, "5 x 10^-8", cex = 0.5, col = "#C0392B", adj = 0)
  text(1.3010 + 0.3, bp[1] - 0.55, "0.05", cex = 0.5, col = "grey35", adj = 0)
  text(rev(p$mlogp) + 0.5, bp, sprintf("%.1f", rev(p$mlogp)), cex = 0.48, adj = 0)
  mtext(expression(-log[10](italic(P))), 1, line = 2.6, cex = 0.75)
  legend("bottomright", legend = c("Bell's palsy", "Skeletal or articular",
                                   "Autoimmune", "Entrapment neuropathy", "Metabolic"),
         fill = c(COL[["target"]], COL[["skeletal"]], COL[["autoimmune"]],
                  COL[["entrapment"]], COL[["metabolic"]]),
         border = "#333333", bty = "n", cex = 0.52)
  panel_letter("c")
}

# ============================================================
# Fig. S1 | genetic relatedness (supplementary)
# ============================================================
figS1 <- function() {
  a <- rd(FD, "fig6a_zcorr.tsv")
  a <- a[a$phenocode != "G6_BELLPA", ]
  a <- a[order(-a$r), ]
  b <- rd(FD, "fig6b_crosstrait_mr.tsv")
  b <- b[order(b$OR), ]
  layout(matrix(1:2, nrow = 2), heights = c(1, 1))

  par(mar = c(7.4, 9.4, 1.5, 0.8), las = 1, mgp = c(2.7, 0.6, 0))
  bp <- barplot(a$r, names.arg = rep("", nrow(a)), col = COL[a$category],
                border = "#333333", lwd = 0.4, ylim = c(0, max(a$r) * 1.26),
                ylab = "", cex.axis = 0.65, xpd = FALSE)
  mtext("z-score correlation with Bell's palsy", 2, line = 3.2, las = 0, cex = 0.7)
  text(bp, -0.008, short(a$trait), srt = 45, adj = c(1, 0.5), cex = 0.55, xpd = NA)
  text(bp, a$r + 0.009, sprintf("%+.3f", a$r), cex = 0.5)
  abline(h = 0, col = "grey30", lwd = 0.6)
  panel_letter("a")

  par(mar = c(3.8, 9.4, 1.0, 0.8))
  y <- seq_len(nrow(b))
  plot(NA, xlim = c(0.86, 3.35), ylim = c(0.4, nrow(b) + 0.6), yaxt = "n",
       xlab = "", ylab = "", bty = "l", log = "x", cex.axis = 0.65)
  abline(v = 1, lty = 3, col = "grey45", lwd = 0.8)
  segments(b$lo, y, b$hi, y, col = COL[b$category], lwd = 1.1)
  points(b$OR, y, pch = 19, cex = 0.7, col = COL[b$category])
  axis(2, at = y, labels = short(b$trait), las = 1, cex.axis = 0.58, tick = FALSE)
  axis(1, at = c(1, 1.5, 2, 2.5, 3), labels = c("1.0", "1.5", "2.0", "2.5", "3.0"),
       cex.axis = 0.65)
  mtext("Odds ratio for Bell's palsy (95% CI)", 1, line = 2.5, cex = 0.68)
  text(3.28, y, paste0("n = ", b$nsnp), cex = 0.48, adj = 1, col = "grey30")
  panel_letter("b")
}


# ============================================================
# Graphical abstract | 论文架构图
# EJHG 规格: 打印 9 x 5 cm, >= 300 dpi, 单文件提交
# 画布按 9:5 比例渲染(7.2 x 4.0 in), 缩放到 9 cm 宽时文字约 5 pt, 仍清晰
# ============================================================
gabstract <- function() {
  par(mar = c(0.15, 0.15, 0.15, 0.15), xpd = NA)
  plot(NA, xlim = c(0, 100), ylim = c(0, 56), axes = FALSE, xlab = "", ylab = "")

  text(49.5, 53.8, "Bell's palsy: genetic architecture in 5,062 Finnish cases",
       cex = 1.02, font = 2, col = "#111111")
  segments(1, 51.6, 99, 51.6, col = "#C0392B", lwd = 1.8)

  panel <- function(x0, x1, y0, y1, title, bullets, hcol, bcol, bcex = 0.64) {
    rect(x0, y0, x1, y1, col = "white", border = "#B8B8B8", lwd = 0.9)
    rect(x0, y1 - 4.4, x1, y1, col = hcol, border = NA)
    text((x0 + x1) / 2, y1 - 2.2, title, cex = 0.88, font = 2, col = "white")
    ys <- seq(y1 - 7.4, y0 + 2.4, length.out = length(bullets))
    for (i in seq_along(bullets))
      text(x0 + 1.7, ys[i], paste0("\u2022  ", bullets[i]), cex = bcex,
           adj = c(0, 0.5), col = bcol[i])
    invisible(NULL)
  }

  DATA <- c("FinnGen R13, endpoint G6_BELLPA",
            "5,062 cases / 431,627 controls",
            "10 comparator traits, one pipeline",
            "72,412-case autoimmune composite",
            "meta-analysis, GCST90018843")
  ANAL <- c("QC and meta-analysis",
            "SuSiE fine-mapping",
            "MHC SNP scan, chr6:25-35 Mb",
            "187 classical HLA alleles",
            "colocalisation, self-validated",
            "credible-set interrogation",
            "phenome scan, 2,754 phenotypes",
            "z-correlation, cross-trait MR")
  FIND <- c("two loci, both outside MHC",
            "chr6:44.48 Mb intergenic",
            "RUNX2, chr6:45.60 Mb",
            "0 of 49 variants in MHC",
            "0 of 187 HLA alleles",
            "not shared with any trait",
            "11 skeletal/articular traits")
  FCOL <- c("#C0392B", "#C0392B", "#C0392B", "#2166AC", "#2166AC", "#4D4D4D", "#762A83")

  panel(1, 30.5, 9.2, 49.4, "DATA", DATA, "#4D4D4D", rep("#1A1A1A", 5))
  panel(35.2, 64.8, 9.2, 49.4, "ANALYSES", ANAL, "#2166AC", rep("#1A1A1A", 8))
  panel(69.5, 99, 9.2, 49.4, "FINDINGS", FIND, "#C0392B", FCOL)

  arrows(31.4, 29, 34.3, 29, length = 0.09, lwd = 2.4, col = "#909090")
  arrows(65.6, 29, 68.5, 29, length = 0.09, lwd = 2.4, col = "#909090")

  rect(1, 0.8, 99, 7.0, col = "#F5F5F5", border = "#B8B8B8", lwd = 0.9)
  text(50, 4.4, "Facial canal and skeleton, not an autoimmune architecture",
       cex = 0.92, font = 2, col = "#C0392B")
  text(50, 2.1, "no MHC signal and no HLA association in 5,062 cases",
       cex = 0.68, col = "#4D4D4D")
}

# ============================================================
cat("Rendering figures...\n")
render("Fig1_genome_wide_architecture", 7.2, 7.4, fig1)
render("Fig2_finemapping_candidate_genes", 7.2, 8.8, fig2)
render("Fig3_credible_set_effects", 7.2, 7.9, fig3)
render("Fig4_no_MHC_signal_skeletal_phenotypes", 7.2, 10.0, fig4)
render("FigS1_genetic_relatedness", 7.2, 6.4, figS1)

# 图形摘要: 打印尺寸 9 x 5 cm, 以 600 dpi 出位图
RENDER_PDF <<- FALSE
png(file.path(OUT, "Graphical_abstract.png"), width = 7.2, height = 4.0,
    units = "in", res = 600, type = "cairo", bg = "white")
gabstract(); dev.off()
pdf(file.path(OUT, "Graphical_abstract.pdf"), width = 7.2, height = 4.0, bg = "white")
gabstract(); dev.off()
cat("   Graphical_abstract\n")
cat("Done.\n")
