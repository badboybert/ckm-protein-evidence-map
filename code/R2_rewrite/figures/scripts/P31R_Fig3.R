# Paper 3.1 rewrite (R2) -- Figure 3: renal worked example, UMOD-centred (decision S2-D24; plan FIGURE_MESSAGES_R2_2026-10-15.md v3
# section 3; figure spec figures/P31R_FIGURE_SPEC.md). This file REPLACES the R1-era generator (kept by the session scratchpad).
#
#   a  UMOD, the declared kidney control (S2-D25): per-signal forest of eGFR* (-ln eGFR) per SD of UMOD at the five conditionally
#      distinct protein signals, 503-person panel (filled, primary) with the 489-person panel beside (open), pointwise 95% intervals from
#      the approved frame figures/data/Figure_3/umod_signal_ratios_se.tsv (ratio and se are W6b outputs written by the S2-D25 re-run of
#      R/08_umod.R; ci = ratio -/+ 1.959964 se, DERIVED in the frame; the frame stores ln(eGFR) per SD, so the sign is flipped and the
#      bounds swapped here). Dashed line: the profiled common ratio theta of the covariance-aware second-order Q (results.tsv), a fit to
#      these same data ("fitted from these data; not an independent prediction"). The palindromic signal is ringed. Right-hand text keeps
#      three methods apart: the omnibus common-ratio test (Q, df, p; p without the palindromic signal), the conditional colocalization
#      (coloc.abf) PP.H4 at the two declared signals on both panels, and the inadmissible reference-LD fine-mapping fit (named as a
#      different method, not used as support; NO glyph: the plan's cross glyph collides with "not colocalized" elsewhere, visual-check
#      finding of the S2-D24 round). The sentence "Individual-ratio intervals are shown; the omnibus test additionally uses
#      between-estimate covariance" is printed: checked against the frame (se is each ratio's own delta-method SE, no between-signal
#      covariance) and gls_q.R (cov = "aware" keeps the off-diagonal covariance, which enters Q only). No per-signal number is printed
#      (the ledger holds none); the axis carries the scale.
#   b  kidney share among the 16 multi-signal candidates against the other 144 pool pairs, as an odds-ratio forest (log scale,
#      95% interval, type per row): primary (filled) and the sensitivity versions (open): (a) symmetric rule, (b) each locus counted once,
#      (c) cluster bootstrap of the Haldane OR (diamond; a different estimator, labelled), (d) inadmissible fits reclassified
#      (OR not finite: an interval from its lower bound with an arrow, NO point marker) and (d) the pre-declared version that keeps
#      ING1-CAD. No p in the header, no enrichment wording, no DPEP1.
#   c  DROPPED (S2-D25): the plan's optional UMOD regional-track panel is not drawn (it needed a download).
# Printed numbers: every one is read from the claims ledger at the point of printing (led(), spec section 3). Plotted series
# (panel b) are read from the two approved W6a/W6b results.tsv files through af3_read() (approved-manifest, path + md5) and are
# asserted equal to the ledger before drawing. Both panels are drawn on a millimetre canvas (1 data unit = 1 mm), so every label
# position is a physical position at the 170-mm print width.
# Run:  Rscript "V3/R2_rewrite/figures/scripts/P31R_Fig3.R"

.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), winslash = "/")
source(file.path(dirname(.this), "P31R_fig_common.R"))
source(file.path(P31R_FIGS, "af3", "P31R_af3_common.R"))

FIG <- "Figure_3"
af3_begin(FIG)

INK   <- "grey15"
GREY  <- "grey45"
DARK  <- PAL_NEUTRAL[["dark"]]
MID   <- PAL_NEUTRAL[["mid"]]
BODY  <- BODY_PT / ggplot2::.pt
LW    <- 0.32          # mm (about 0.9 pt): intervals, axis, box
LW_REF <- 0.25         # mm (about 0.7 pt): dashed reference line

# ---- sources (approved manifest only; fail closed) -----------------------------------------------------------------
SRC_A <- "R1_stage1/S1R1_R2_06a_renal/results.tsv"            # W6a: primary and sensitivities (a), (b), (c)
SRC_D <- "R1_stage1/S1R1_R2_06b_renal_d_umod/results.tsv"     # W6b: sensitivity (d), UMOD test and conditional PP.H4
SRC_SE <- "R2_rewrite/figures/data/Figure_3/umod_signal_ratios_se.tsv"   # per-signal ratio, se, 95% interval (S2-D25)
ra <- af3_read(SRC_A, "b", quote = "")
rb <- af3_read(SRC_D, "b", quote = "")
rv <- function(tab, id) {
  v <- tab$value[tab$quantity_id == id]
  if (length(v) != 1) stop("results.tsv: ", id, " found ", length(v), " times")
  as.numeric(v)
}
same <- function(x, id, idx = NULL, what = id) {      # plotted/stored value == ledger value (full precision; Inf == Inf)
  y <- p31r_ledger(id, idx)
  ok <- if (is.infinite(x)) identical(x, y) else abs(x - y) <= 1e-8 * max(1, abs(y))
  if (!isTRUE(ok)) stop("source value for ", what, " (", format(x, digits = 12), ") differs from ledger ", id, " (", format(y, digits = 12), ")")
  invisible(x)
}

# ================================================================== panel b: data
# row, label, status, marker, estimate / low / high (results.tsv id) and the ledger element each value is printed from
SPEC <- list(
  list(k = "primary", st = "primary",     mk = "circle",  est = c(SRC_A, "S6.primary.or_mle", "L862"),          lo = c(SRC_A, "S6.primary.or_lo", "L861"),         hi = c(SRC_A, "S6.primary.or_hi", "L860")),
  list(k = "a",       st = "sensitivity", mk = "circle",  est = c(SRC_A, "S6.b.a_pair.or_mle", "L849"),         lo = c(SRC_A, "S6.b.a_pair.or_lo", "L848"),        hi = c(SRC_A, "S6.b.a_pair.or_hi", "L847")),
  list(k = "b",       st = "sensitivity", mk = "circle",  est = c(SRC_A, "S6.b.locus_outcome.or_mle", "L853"),  lo = c(SRC_A, "S6.b.locus_outcome.or_lo", "L852"), hi = c(SRC_A, "S6.b.locus_outcome.or_hi", "L851")),
  list(k = "c",       st = "sensitivity", mk = "diamond", est = c(SRC_A, "S6.c.haldane_or", "L854"),            lo = c(SRC_A, "S6.c.pooled_lo", "L856"),           hi = c(SRC_A, "S6.c.pooled_hi", "L855")),
  list(k = "d",       st = "sensitivity", mk = "none",    est = c(SRC_D, "S6.d.or", NA),                        lo = c(SRC_D, "S6.d.lo", "L885"),                  hi = c(SRC_D, "S6.d.hi", NA)),
  list(k = "dkept",   st = "sensitivity", mk = "circle",  est = c(SRC_D, "S6.d.kept.or", "L872"),               lo = c(SRC_D, "S6.d.kept.lo", "L871"),             hi = c(SRC_D, "S6.d.kept.hi", "L870")))
tab <- function(src) if (src == SRC_A) ra else rb
FR <- data.table::rbindlist(lapply(SPEC, function(s) {
  g <- function(z) rv(tab(z[1]), z[2])
  e <- g(s$est); l <- g(s$lo); h <- g(s$hi)
  for (z in list(s$est, s$lo, s$hi)) if (!is.na(z[3])) same(g(z), z[3], what = z[2])
  data.table::data.table(vkey = s$k, status = s$st, marker = s$mk, or = e, ci_lo = l, ci_hi = h,
                         source_file = s$est[1], ids = paste(s$est[2], s$lo[2], s$hi[2], sep = ";"),
                         ledger = paste(s$est[3], s$lo[3], s$hi[3], sep = ";"))
}))
# the (d) row: the source itself says the OR and the upper bound are not finite (S6.d.or = S6.d.hi = Inf); the drawn form follows
stopifnot(!is.finite(FR[vkey == "d", or]), !is.finite(FR[vkey == "d", ci_hi]), is.finite(FR[vkey == "d", ci_lo]),
          all(is.finite(unlist(FR[vkey != "d", .(or, ci_lo, ci_hi)]))), all(FR[vkey != "d", ci_lo < or & or < ci_hi]))
# the two pool counts of the header are asserted against W6a's own count rows as well as the ledger
stopifnot(rv(ra, "S6.primary.a") == p31r_ledger("L264", 1), rv(ra, "S6.primary.a") + rv(ra, "S6.primary.b") == p31r_ledger("L264", 2),
          rv(ra, "S6.primary.c") == p31r_ledger("L264", 6), rv(ra, "S6.primary.c") + rv(ra, "S6.primary.d") == p31r_ledger("L264", 7))
n_recl <- as.numeric(rb$n[rb$quantity_id == "S6.d.reclassified"]); stopifnot(n_recl == p31r_ledger("L886"))

# ================================================================== panel a: stored values cross-checked against the ledger
same(rv(rb, "S6.umod.Q.aware_2nd.Q"), "L876"); same(rv(rb, "S6.umod.Q.aware_2nd.df"), "L877"); same(rv(rb, "S6.umod.Q.aware_2nd.p"), "L878")
same(rv(rb, "S6.umod.H4.A.503"), "L888"); same(rv(rb, "S6.umod.H4.B.503"), "L890"); same(rv(rb, "S6.umod.H4.A.489"), "L887"); same(rv(rb, "S6.umod.H4.B.489"), "L889")
same(rv(rb, "S6.d.UMOD_eGFR.s_max"), "L869")

# ================================================================== canvas helpers (1 data unit = 1 mm)
W <- 170
canvas <- function(H) ggplot() +
  scale_x_continuous(limits = c(0, W), expand = c(0, 0)) + scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))
tx <- function(x, y, label, hjust = 0, size = ANNOT, fontface = "plain", vjust = 0.5, colour = INK, lineheight = 0.95)
  annotate("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust, size = size, family = BASE_FONT,
           fontface = fontface, colour = colour, lineheight = lineheight)
tag <- function(letter, H) annotate("text", x = 0.6, y = H - 2.6, label = letter, hjust = 0, vjust = 0.5, size = TAG_PT / ggplot2::.pt,
                                     family = BASE_FONT, fontface = "bold", colour = "black")
num <- function(d) function(v) p31r_num(v, digits = d)

# ================================================================== panel a (per-signal forest + stored results), 86 mm
# Forest: eGFR* (= -ln eGFR) per SD of UMOD at the five conditionally distinct protein signals, 503-person panel (filled, primary) with the
# 489-person panel beside (open). Estimates and pointwise 95% intervals are read from the approved frame (ratio and se from the W6b step,
# interval = ratio +/- 1.959964 se, DERIVED in the frame); the frame stores ln(eGFR) per SD, so the sign is flipped here and the bounds swapped.
# No per-signal number is printed (the ledger holds none): the axis carries the scale. The common ratio is the profiled theta of the
# covariance-aware second-order Q (results.tsv), a fit to these same data.
se <- af3_read(SRC_SE, "a", quote = "")
stopifnot(nrow(se) == 10, setequal(se$panel, c(503, 489)), all(table(se$panel) == 5))
se <- se[order(-panel, signal)]
stopifnot(identical(as.integer(se$signal), rep(1:5, 2)), all(is.finite(c(se$ratio, se$se, se$ci_lo, se$ci_hi))), all(se$se > 0),
          all(se$ci_lo < se$ratio & se$ratio < se$ci_hi),
          all(abs((se$ci_hi - se$ci_lo) / 2 - 1.959964 * se$se) < 1e-9),                         # the DERIVED interval is ratio +/- 1.959964 se
          identical(se[panel == 503]$rsid, se[panel == 489]$rsid),
          identical(se[panel == 503]$rsid[1:2], c("rs12922822", "rs77924615")),                   # declared signals A and B (W6b results.tsv notes)
          identical(which(se[panel == 503]$palindromic == 1L), 5L), identical(which(se[panel == 489]$palindromic == 1L), 5L))
rba <- af3_read(SRC_D, "a", quote = "")                                                          # same file, read for panel a (registered as such for G12)
theta <- rv(rba, "S6.umod.Q.aware_2nd.theta")                                                      # 503 panel, primary cell (ln eGFR per SD)
theta489 <- rv(rba, "S6.umod.Q.p489.aware_2nd.theta")
stopifnot(theta < 0, abs(theta - theta489) < 2e-4, theta > min(se[panel == 503]$ratio), theta < max(se[panel == 503]$ratio))
se[, `:=`(est = -ratio, lo = -ci_hi, hi = -ci_lo)]                                               # eGFR* = -ln(eGFR): sign flipped, bounds swapped
stopifnot(se[, all(lo < est & est < hi & lo > 0)])
THETA <- -theta

HA <- 86
X0 <- 40; X1 <- 100; V0 <- -0.001; V1 <- 0.021
xv <- function(v) X0 + (v - V0) / (V1 - V0) * (X1 - X0)
yc <- HA - 21 - (0:4) * 8                       # signal row centres
OFF <- 1.7                                      # 503 panel above, 489 panel below the centre
se[, y := yc[signal] + ifelse(panel == 503, OFF, -OFF)]
Y_AXA <- yc[5] - 6.2
tk_a <- ticks("3a_x_egfrstar", c(0, 0.01, 0.02), c(V0, V1), "linear", "panel a, eGFR* per SD of UMOD axis")
tkx_a <- xv(c(0, 0.01, 0.02))
capsule <- function(x0, x1, ycen, r, n = 24) {
  th <- seq(pi / 2, -pi / 2, length.out = n); rt <- cbind(x1 + r * cos(th), ycen + r * sin(th))
  th2 <- seq(-pi / 2, -3 * pi / 2, length.out = n); lf <- cbind(x0 + r * cos(th2), ycen + r * sin(th2))
  m <- rbind(rt, lf, rt[1, ]); data.frame(x = m[, 1], y = m[, 2])
}
cap <- capsule(xv(min(se[signal == 5]$lo)) - 2.2, xv(max(se[signal == 5]$hi)) + 2.2, yc[5], OFF + 2.2)
lab_a <- c("rs12922822 (signal A)", "rs77924615 (signal B)", "rs28510439", "rs149109606", "rs62033284 (palindromic)")
XR2 <- 108                                      # right text column

pa <- canvas(HA) + tag("a", HA) +
  tx(7.5, HA - 2.6, "UMOD, the declared control: eGFR* per SD of UMOD at each of five protein signals", size = BODY, hjust = 0) +
  tx(2, HA - 9.2, "Lead variant", fontface = "bold") +
  tx(XR2, HA - 9.2, "Conditional colocalization with eGFR", fontface = "bold") +
  tx(XR2, HA - 12.4, "(coloc.abf) at the two declared signals", fontface = "bold") +
  tx(2, yc, lab_a) +
  # common ratio (fitted from these data), then the data
  annotate("segment", x = xv(THETA), xend = xv(THETA), y = yc[1] + OFF + 3.6, yend = yc[5] - OFF - 3.6, colour = GREY, linewidth = LW_REF, linetype = "dashed") +
  tx(xv(THETA), yc[1] + OFF + 5.8, "Common ratio, fitted from these data", hjust = 0.5, colour = GREY) +
  annotate("path", x = cap$x, y = cap$y, colour = INK, linewidth = LW_REF) +
  annotate("segment", x = xv(se$lo), xend = xv(se$hi), y = se$y, yend = se$y, colour = DARK, linewidth = LW) +
  annotate("point", x = xv(se[panel == 503]$est), y = se[panel == 503]$y, shape = 16, size = 2.4, colour = DARK) +
  annotate("point", x = xv(se[panel == 489]$est), y = se[panel == 489]$y, shape = 1, size = 2.4, stroke = 0.55, colour = DARK) +
  # axis
  annotate("segment", x = xv(0), xend = xv(0.02), y = Y_AXA, yend = Y_AXA, colour = INK, linewidth = LW) +
  annotate("segment", x = tkx_a, xend = tkx_a, y = Y_AXA, yend = Y_AXA - 1.3, colour = INK, linewidth = LW) +
  tx(tkx_a, Y_AXA - 4.0, tk_a, hjust = 0.5, size = BODY) +
  tx((X0 + X1) / 2, Y_AXA - 8.6, "eGFR* (−natural-log eGFR) per SD of UMOD", hjust = 0.5, size = BODY)

# right column rows A and B: the conditional PP.H4 of each declared signal, each panel named
ppA <- paste0("PP.H4 ", led("L888", NULL, num(2), "panel a, conditional PP.H4 signal A, 503", "PP.H4 {} \\("), " (", led("L936", NULL, num(0), "panel a", "{} panel"),
              " panel), ", led("L887", NULL, num(2), "panel a, conditional PP.H4 signal A, 489", "panel\\), {} \\(489"), " (", led("L937", NULL, num(0), "panel a", "{} panel"), " panel)")
ppB <- paste0("PP.H4 ", led("L890", NULL, num(2), "panel a, conditional PP.H4 signal B, 503", "PP.H4 {} \\("), " (", led("L936", NULL, num(0), "panel a", "{} panel"),
              " panel), ", led("L889", NULL, num(2), "panel a, conditional PP.H4 signal B, 489", "panel\\), {} \\(489"), " (", led("L937", NULL, num(0), "panel a", "{} panel"), " panel)")
pa <- pa + tx(XR2, yc[1], ppA) + tx(XR2, yc[2], ppB)

# right column: the omnibus test, then the fine-mapping fit (a different method, named on its own)
step <- 3.15
ln <- function(label, ytop, ...) tx(XR2, ytop, label, ...)
rule <- function(ytop) annotate("segment", x = XR2, xend = W - 1, y = ytop, yend = ytop, colour = GREY, linewidth = LW_REF)
T1 <- yc[3] + 2.8
q1 <- paste0("Q ", led("L876", NULL, num(2), "panel a, common-ratio test", "Q {},"), ", df ", led("L877", NULL, num(0), "panel a", "df {},"),
             ", p ", led("L878", NULL, num(2), "panel a", "p {} "), " (", led("L936", NULL, num(0), "panel a", "{} panel"), " panel).")
pa <- pa + tx(XR2, T1, "Test of a common ratio (five signals)", fontface = "bold") +
  ln(q1, T1 - step) + ln("No departure detected; this does not show", T1 - 2 * step) +
  ln("proportionality. Without the palindromic", T1 - 3 * step) +
  ln(paste0("signal, p ", led("L965", 3, num(3), "panel a, palindromic signal dropped", "p {}"), "."), T1 - 4 * step)
T2 <- T1 - 5.6 * step
f1 <- paste0("Inadmissible: s ", led("L869", NULL, num(3), "panel a, diagnostic maximum s", "s {} against"), " against a limit of ",
             led("L923", NULL, num(2), "panel a, admissibility limit", "limit of {};"), ";")
f2 <- "not used as support. Robustness to the"
pa <- pa + rule(T2 + 2.2) + tx(XR2, T2, "Fine-mapping fit with reference LD", fontface = "bold") +
  ln(f1, T2 - step) + ln(f2, T2 - 2 * step) +
  ln("reference LD was evaluated only on the", T2 - 3 * step) + ln("all-signal test.", T2 - 4 * step)

# key (inside the graphic, bottom left) and the interval note (bottom right)
ky <- c(13.4, 10.2, 7.0, 3.8)
pa <- pa +
  annotate("point", x = 3.6, y = ky[1], shape = 16, size = 2.4, colour = DARK) +
  tx(8, ky[1], paste0("Reference panel of ", led("L936", NULL, num(0), "panel a, key", "{} people"), " people (primary)")) +
  annotate("point", x = 3.6, y = ky[2], shape = 1, size = 2.4, stroke = 0.55, colour = DARK) +
  tx(8, ky[2], paste0("Nested subset of ", led("L937", NULL, num(0), "panel a, key", "{} people"), " people")) +
  annotate("segment", x = 1.6, xend = 5.6, y = ky[3], yend = ky[3], colour = GREY, linewidth = LW_REF, linetype = "dashed") +
  tx(8, ky[3], "Common ratio; not an independent prediction") +
  annotate("path", x = capsule(2.4, 4.8, ky[4], 1.3)$x, y = capsule(2.4, 4.8, ky[4], 1.3)$y, colour = INK, linewidth = LW_REF) +
  tx(8, ky[4], "Ringed: palindromic signal (A/T alleles)") +
  tx(XR2, ky[1], "Bars: pointwise 95% intervals. Individual-ratio") +
  tx(XR2, ky[2], "intervals are shown; the omnibus test") +
  tx(XR2, ky[3], "additionally uses between-estimate covariance.")

# ================================================================== panel b (odds-ratio forest), 80 mm
HB <- 80
X0 <- 66; X1 <- 126                       # forest x range, mm: OR 0.1 ... 200
LO10 <- log10(0.1); HI10 <- log10(200); PERDEC <- (X1 - X0) / (HI10 - LO10)
xm <- function(or) X0 + (log10(or) - LO10) * PERDEC
XC <- 130                                 # right text column
Y_HEAD <- HB - 2.6; Y_S1 <- HB - 8.4; Y_S2 <- HB - 12.0; Y_COL <- HB - 17.4
Y1 <- HB - 23.5; DY <- 7.4
yr <- Y1 - (seq_len(nrow(FR)) - 1) * DY
Y_AX <- tail(yr, 1) - 5.4
ticks_lab <- ticks("3b_x_or", c(0.1, 1, 10, 100), c(0.05, 300), "log10", "panel b, odds-ratio axis")
tick_x <- xm(c(0.1, 1, 10, 100))
n16 <- led("L264", 2, num(0), "panel b, header", "of {} candidates vs")
n11 <- led("L264", 1, num(0), "panel b, header", "^{} of ")
n57 <- led("L264", 6, num(0), "panel b, header", "vs {} of")
n144 <- led("L264", 7, num(0), "panel b, header", "of {} other")
sub1 <- paste0(n11, " of ", n16, " candidates vs ", n57, " of ", n144, " other pool pairs; exploratory and descriptive")
sub2 <- paste0(led("L886", NULL, num(0), "panel b, inadmissible-fit count", "^{} of 16"), " of ", led("L264", 2, num(0), "panel b, inadmissible-fit denominator", "of {} candidates rest"),
               " candidates rest on inadmissible fits")
fmt3 <- function(id, d, anchor, where) led(id, NULL, num(d), where, anchor)
row_txt <- function(k) {
  switch(k,
    primary = paste0(fmt3("L862", 2, "^{} \\(", "panel b, primary OR"), " (", fmt3("L861", 3, "\\({} to", "panel b, primary low"), " to ", fmt3("L860", 1, "to {}\\)$", "panel b, primary high"), ")"),
    a       = paste0(fmt3("L849", 2, "^{} \\(", "panel b, (a) OR"), " (", fmt3("L848", 3, "\\({} to", "panel b, (a) low"), " to ", fmt3("L847", 1, "to {}\\)$", "panel b, (a) high"), ")"),
    b       = paste0(fmt3("L853", 2, "^{} \\(", "panel b, (b) OR"), " (", fmt3("L852", 3, "\\({} to", "panel b, (b) low"), " to ", fmt3("L851", 1, "to {}\\)$", "panel b, (b) high"), ")"),
    c       = paste0(fmt3("L854", 2, "^{} \\(", "panel b, (c) Haldane OR"), " (", fmt3("L856", 2, "\\({} to", "panel b, (c) low"), " to ", fmt3("L855", 1, "to {}\\)$", "panel b, (c) high"), ")"),
    d       = paste0("OR not finite\ninterval from ", fmt3("L885", 3, "from {},", "panel b, (d) low"), ",\nno finite upper bound"),
    dkept   = paste0(fmt3("L872", 2, "^{} \\(", "panel b, (d) kept OR"), " (", fmt3("L871", 3, "\\({} to", "panel b, (d) kept low"), " to ", fmt3("L870", 1, "to {}\\)$", "panel b, (d) kept high"), ")"))
}
row_lab <- c(primary = "Primary: pair level,\nconditional-MLE OR",
             a = "(a) Symmetric acceptance rule",
             b = "(b) Each locus counted once",
             c = "(c) Cluster bootstrap of the Haldane OR\n(a different estimator)",
             d = "(d) Candidates with inadmissible\nfits reclassified",
             dkept = "(d) Pre-declared version,\nING1–CAD kept")
FR[, text := vapply(vkey, row_txt, "")]
FR[, label := row_lab[vkey]]
FR[, y := yr]
Fin <- FR[is.finite(ci_hi)]
pb <- canvas(HB) + tag("b", HB) +
  tx(7.5, Y_HEAD, "Kidney share among multi-signal candidates, as an odds ratio against the other pool pairs", size = BODY, hjust = 0) +
  tx(2, Y_S1, sub1) + tx(2, Y_S2, sub2) +
  tx(2, Y_COL, "Version", fontface = "bold") + tx(XC, Y_COL, "OR (95% interval)", fontface = "bold") +
  # null line, then the data
  annotate("segment", x = xm(1), xend = xm(1), y = Y_AX, yend = Y1 + DY * 0.55, colour = GREY, linewidth = LW_REF, linetype = "dashed") +
  annotate("segment", x = xm(Fin$ci_lo), xend = xm(Fin$ci_hi), y = Fin$y, yend = Fin$y, colour = DARK, linewidth = LW) +
  annotate("segment", x = xm(FR[vkey == "d", ci_lo]), xend = X1, y = FR[vkey == "d", y], yend = FR[vkey == "d", y], colour = DARK, linewidth = LW,
           arrow = arrow(length = unit(1.7, "mm"), angle = 25, type = "closed")) +
  annotate("point", x = xm(FR[marker == "circle" & status == "primary", or]), y = FR[marker == "circle" & status == "primary", y], shape = 16, size = 2.6, colour = DARK) +
  annotate("point", x = xm(FR[marker == "circle" & status == "sensitivity", or]), y = FR[marker == "circle" & status == "sensitivity", y], shape = 1, size = 2.6, stroke = 0.55, colour = DARK) +
  annotate("point", x = xm(FR[marker == "diamond", or]), y = FR[marker == "diamond", y], shape = 5, size = 2.8, stroke = 0.55, colour = DARK) +
  tx(2, FR$y, FR$label) + tx(XC, FR$y, FR$text) +
  # x axis, drawn on the same canvas (tick labels 9 pt = body)
  annotate("segment", x = X0, xend = X1, y = Y_AX, yend = Y_AX, colour = INK, linewidth = LW) +
  annotate("segment", x = tick_x, xend = tick_x, y = Y_AX, yend = Y_AX - 1.3, colour = INK, linewidth = LW) +
  tx(tick_x, Y_AX - 4.0, ticks_lab, hjust = 0.5, size = BODY) +
  tx((X0 + X1) / 2, Y_AX - 8.6, "Odds ratio (log scale)", hjust = 0.5, size = BODY)
# key, inside the graphic (bottom left)
ky <- c(12.6, 9.2, 5.8, 2.4)
pb <- pb +
  annotate("point", x = 3.6, y = ky[1], shape = 16, size = 2.4, colour = DARK) + tx(8, ky[1], "Primary (filled)") +
  annotate("point", x = 3.6, y = ky[2], shape = 1, size = 2.4, stroke = 0.55, colour = DARK) + tx(8, ky[2], "Sensitivity (open)") +
  annotate("point", x = 3.6, y = ky[3], shape = 5, size = 2.6, stroke = 0.55, colour = DARK) + tx(8, ky[3], "Haldane-corrected OR, a different estimator") +
  annotate("segment", x = 1.6, xend = 5.6, y = ky[4], yend = ky[4], colour = DARK, linewidth = LW, arrow = arrow(length = unit(1.5, "mm"), angle = 25, type = "closed")) +
  tx(8, ky[4], "No finite upper bound (OR not finite)")

# ================================================================== compose and save
P <- (pa / pb) + plot_layout(heights = c(HA, HB))
H_MM <- HA + HB

# frames = what was plotted/printed, written beside the figure (af3/source_data/Figure_3__<name>.csv)
fr_b <- data.table::copy(FR)[, .(row = seq_len(.N), vkey, status, marker, or, ci_lo, ci_hi, printed_text = gsub("\n", " | ", text), source_file, quantity_ids = ids, ledger_ids = ledger)]
fr_b[, ci_hi_note := ifelse(is.finite(ci_hi), "", "no finite upper bound: drawn as an arrow to the axis end; no point marker")]
fr_a <- data.table::data.table(
  item = c("Q", "df", "p", "p without palindromic signal", "PP.H4 signal A, 503 panel", "PP.H4 signal B, 503 panel", "PP.H4 signal A, 489 panel",
           "PP.H4 signal B, 489 panel", "diagnostic max s (fine-mapping fit)", "admissibility limit"),
  ledger_id = c("L876", "L877", "L878", "L965", "L888", "L890", "L887", "L889", "L869", "L923"),
  value = c(rv(rb, "S6.umod.Q.aware_2nd.Q"), rv(rb, "S6.umod.Q.aware_2nd.df"), rv(rb, "S6.umod.Q.aware_2nd.p"), p31r_ledger("L965", 3),
            rv(rb, "S6.umod.H4.A.503"), rv(rb, "S6.umod.H4.B.503"), rv(rb, "S6.umod.H4.A.489"), rv(rb, "S6.umod.H4.B.489"),
            rv(rb, "S6.d.UMOD_eGFR.s_max"), p31r_ledger("L923")),
  source = c(rep(SRC_D, 3), "R1_stage1/S1R1_R2_1b_verify/S6/adjudication/A_third/a3_umod_q.tsv (ledger L965 only; file not opened)", rep(SRC_D, 5), "ledger L923"))
fr_a[, panel_state := "values printed as text in panel a (test, conditional colocalization, fine-mapping fit); each is read from the ledger at print"]
# the forest as drawn: eGFR* = -ln eGFR, so estimate = -ratio and the bounds are swapped; the common ratio is -theta (503 panel)
fr_f <- data.table::copy(se)[, .(panel, signal, rsid, palindromic, plotted_est_egfrstar = est, plotted_lo = lo, plotted_hi = hi, stored_ratio = ratio, stored_se = se,
                                  stored_ci_lo = ci_lo, stored_ci_hi = ci_hi, common_ratio_plotted = THETA, source_file = SRC_SE,
                                  note = "pointwise 95% interval = ratio -/+ 1.959964 se (DERIVED in the frame); sign flipped for eGFR* = -ln(eGFR)")]

# merge duplicate registry entries (same id, idx, text, anchor): the gate sums counts within a (text, anchor) group, so this is cosmetic
af3_save(P, H_MM, panels = c("a", "b"),
         panel_map = c(a = "UMOD worked control: per-signal forest of eGFR* per SD of UMOD at five protein signals (503 panel, 489 beside) with pointwise 95% intervals, the common-ratio test, conditional colocalization and the inadmissible fine-mapping fit",
                       b = "Kidney share among the 16 multi-signal candidates vs the other 144 pool pairs: odds-ratio forest, primary and sensitivity versions"),
         frames = list(`3a_forest` = fr_f, `3a_stored_results` = fr_a, `3b_forest` = fr_b))
cat(sprintf("Figure_3: panel a forest (%d rows), panel b %d rows, canvas %d x %d mm
", nrow(se), nrow(FR), W, H_MM))
