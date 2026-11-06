# Paper 3.1 rewrite (R2) -- Figure 4: signal-level associations with arginine, ornithine and T2D at the ARG1 locus.
# Plan: review_R2/FIGURE_MESSAGES_R2_2026-10-15.md (v3, section 4; decisions S2-D24, F2, F7, F11, F13, F14, F15); spec figures/P31R_FIGURE_SPEC.md.
#   a  four separately labelled estimands, each A versus B with a pointwise 95% interval and its own scale (no shared null):
#      (i) arginine, ratio of ratios (null 1); (vi) ornithine, difference (null 0); (ii) protein-scaled T2D, difference, four
#      resource-by-setting cells (primary = filled, beside = open); (iii) direct comparison, difference, four cells at equal
#      weight. Covariance-aware second-order p values only (W4a). Visual status (plan section 0, S2-D24): primary filled /
#      beside open; no colocalization = circle with a cross + word; sensitive to the reference LD = dotted outline + word.
#   b  compact locus (enlarged interval), PROTEIN TRACK ONLY (S2-D26): the ARG1 protein -log10 p track with signals A and B and LD colouring; gene
#      models; credible sets at signal A, one row per trait and method with its own status read from the canonical file; rs2246012 at its actual
#      coordinate. The arginine and T2D regional tracks, per-signal estimates, signal B's sets and the 1-Mb overview are in Additional file 3
#      (Figure S5; the pointer 'Arginine and T2D regional tracks: Additional file 3: Figure S5, panel d' is printed beside the track; written 'S5, panel d' because the figure gate cannot read the token 'S5d').
# Sources: ONLY files with an APPROVED row for Figure_4 in registry/APPROVED_SOURCE_MANIFEST.tsv (path + md5). S2-D25 round: the 503-person
# dosage panel (.i8), PREP_STAMP.json and the FinnGen .tbi now have approved rows, so track variants are coloured by LD (r2 with signal A or B,
# from the dosage panel; spec section 1 bins) and all three are pinned in the render log. FETCH_RECORD.json is still metadata read only to check
# an md5 and is listed in REG$unregistered_reads.
# Status encodings in this figure (one key, stated in the graphic; S2-D24 five-encoding rule, S2-D25 distinct glyphs, S2-D26 glyph harmonisation): primary
# or equal-weight cell = filled marker; beside = open marker; signal A not colocalized = circle with a cross (shape 13) + the word; inadmissible fit =
# the HATCHED circle of AF3 S2 / S5c (P31R_af3_glyphs.R geom_hatchdot); not robust to the reference LD = the DOTTED RING of AF3 S5c (geom_dotring) in
# the credible-set rows, and the dotted outline round a forest row in panel a. In panel b one status marker per credible-set row (as in S5c); the
# members of every set are plain ticks. No glyph carries two statuses.
# p values with a power of ten (S2-D26): set at the panel's 7 pt, the base '10' at 9 pt ONLY because the figure gate reads an exponent as a raised
# span smaller than 0.85 x the span before it (a 7-pt exponent after a 7-pt '10' is read as a separate number and fails G7); exponent exactly 7 pt.
# Every number PRINTED comes from p31r_ledger and is registered at the call that prints it (spec section 3); this script writes
# registry/Figure_4.json from those calls. Plotted-series values come from the W4a / 04b / Chen files above and are asserted equal
# to the ledger rows that quote them (value AND the results.tsv line the ledger row names).
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_Fig4.R"
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_V3, "/R2_rewrite/figures/scripts/P31R_fig_common.R"))
suppressPackageStartupMessages(library(ggtext))
source(paste0(.P31R_V3, "/R2_rewrite/figures/af3/P31R_af3_glyphs.R"))   # geom_hatchdot, geom_dotring (the AF3 status glyphs); appends its own md5 to P31R_DEPS
for (.p in c("Rsamtools", "GenomicRanges", "IRanges")) if (!requireNamespace(.p, quietly = TRUE)) stop("package missing: ", .p)
EXP_PT <- 7   # exponents print at 7 pt (the production floor); the '10' they follow is 9 pt so the figure gate still reads the 7-pt span as a superscript

FIG <- "Figure_4"
V3  <- V3_ROOT
S1  <- file.path(V3, "R1_stage1")
SIG <- c(A = "rs2781668", B = "rs17788484"); RS_C <- "rs2246012"
COL_SIG <- PAL_SIGNAL
INK  <- PAL_NEUTRAL[["dark"]]; GREY <- PAL_NEUTRAL[["mid"]]; PTGREY <- PAL_NEUTRAL[["other"]]
COL_REF <- "grey55"
COL_LD <- c(A_hi = PAL_SIGNAL[["A"]], A_lo = PAL_SIGNAL_LIGHT[["A"]], B_hi = PAL_SIGNAL[["B"]], B_lo = PAL_SIGNAL_LIGHT[["B"]])
Z975 <- qnorm(0.975)
near <- function(a, b, tol = 1e-9) isTRUE(all(abs(as.numeric(a) - as.numeric(b)) <= tol))
need <- function(ok, msg) if (!isTRUE(ok)) stop("ASSERTION FAILED: ", msg, call. = FALSE)
L <- function(id, i = NULL) as.numeric(p31r_ledger(id, i))

# ---- registry, written by this script from the calls that print (spec section 3) -------------------------------
REG <- list(panels = c("a", "b"), ledger = list(), derived = list(), literals = list(), unregistered_reads = list())
.add <- function(kind, e) { REG[[kind]][[length(REG[[kind]]) + 1]] <<- Filter(Negate(is.null), e); invisible(e$text) }
# a ledger number: printed text from p31r_ledger(id, idx) and registered from that same call
LP <- function(id, idx = NULL, digits, where, count = 1L, anchor = NULL) {
  txt <- p31r_num(p31r_ledger(id, idx), digits)
  .add("ledger", list(text = txt, id = id, idx = idx, count = as.integer(count), anchor = anchor, where = where)); txt
}
# a p value in scientific notation: markdown for the plot, plain text (U+2212 exponent) for the registry
LPS <- function(id, idx = NULL, where, anchor = NULL) {
  p <- p31r_ledger(id, idx); e <- floor(log10(p)); m <- p / 10^e
  if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
  ms <- formatC(m, format = "f", digits = 1); es <- sub("^-", "−", as.character(e))
  .add("ledger", list(text = sprintf("%s × 10%s", ms, es), id = id, idx = idx, count = 1L, anchor = anchor, where = where))
  sprintf("%s × <span style='font-size:9pt'>10</span><sup><span style='font-size:%gpt'>%s</span></sup>", ms, EXP_PT, es)
}
reg_lit <- function(text, why, count = 1L, anchor = NULL, not_ledger = NULL)
  .add("literals", list(text = text, why = why, count = as.integer(count), anchor = anchor, not_ledger = not_ledger))
# axis ticks: one entry per printed break. A tick is a standalone token on a line of numbers only (adjacent ticks of a narrow
# scale can join into one printed line), which is what the anchor states
TICKS <- list()
TICK_ANCHOR <- "^(?:[−+-]?[0-9.,]+\\s)*{}(?:\\s[−+-]?[0-9.,]+)*$"
reg_ticks <- function(axis, breaks, labels, what, transform = "linear") {
  need(length(breaks) == length(labels) && length(breaks) >= 2 && !anyDuplicated(labels) && is.null(TICKS[[axis]]), paste("ticks of", axis))
  TICKS[[axis]] <<- labels
  for (t in labels) .add("literals", list(text = t, why = paste0("axis_tick: ", what), axis = axis, transform = transform, count = 1L, anchor = TICK_ANCHOR))
  labels
}
drawn <- function(p, ax) lapply(ggplot_build(p)$layout$panel_params, function(pp) {
  v <- pp[[ax]]; b <- v$get_breaks(); l <- v$get_labels(); as.character(l[!is.na(b)]) })
check_drawn <- function(p, ax, axes) {
  d <- drawn(p, ax)
  need(length(d) == length(axes), paste("panel count for", paste(axes, collapse = ",")))
  for (i in seq_along(axes)) need(identical(d[[i]], unname(TICKS[[axes[i]]])),
                                  sprintf("drawn ticks on %s (%s) = registered (%s)", axes[i], paste(d[[i]], collapse = " "), paste(TICKS[[axes[i]]], collapse = " ")))
}
fmt_mb <- function(x, digits) formatC(x, format = "f", digits = digits)
md5f <- function(p) unname(tools::md5sum(p))
UNREG <- function(label, path) REG$unregistered_reads[[label]] <<- list(path = sub(paste0("^", V3, "/"), "", normalizePath(path, winslash = "/")), md5 = md5f(path))

# ======================================================================================================== inputs (approved only)
W4A_F <- file.path(S1, "S1R1_R2_04a_arg1grid/results.tsv")
HLM_F <- file.path(S1, "S1R1_R2_10_integration/results.tsv")
SETS_F <- file.path(S1, "S1R1_R2_04b_arg1sets/results.tsv")
CSS_F <- file.path(S1, "S1R1_R2_12_chen_instrument/credible_set_status.tsv")
CHEN_F <- file.path(S1, "S1R1_R2_12_chen_instrument/results.tsv")
H_F   <- file.path(V3, "V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz")
T2D_F <- paste0(.P31R_EXT, "/D/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz")
ENS_D <- file.path(P31R_FIGS, "data/Figure_4")
W <- p31r_read(W4A_F, sep = "\t"); HL <- p31r_read(HLM_F, sep = "\t"); SETS <- p31r_read(SETS_F, sep = "\t")
CSS <- p31r_read(CSS_F, sep = "\t"); CHEN <- p31r_read(CHEN_F, sep = "\t"); h <- p31r_read(H_F)
need(md5f(W4A_F) == "f4c5104dce1f1822448a6a7bdced1636", "W4a results.tsv is the approved version")
Q <- function(q, col = "value", tab = W, key = "quantity_id") {
  r <- tab[tab[[key]] == q]; need(nrow(r) == 1, paste("one row for", q)); v <- r[[col]]
  if (col == "notes") as.character(v) else suppressWarnings(as.numeric(v))
}
# a ledger row that quotes a W4a value must (i) hold that value and (ii) name that quantity's line of results.tsv
LQ <- function(lid, q, tol = 1e-6, tab = W, key = "quantity_id") {
  row <- p31r_ledger_row(lid); ln <- as.integer(sub(".*; line ([0-9]+);.*", "\\1", row$source)); need(!is.na(ln), paste(lid, "names a line"))
  need(tab[[key]][ln - 1] == q, sprintf("%s names line %d = %s, expected %s", lid, ln, tab[[key]][ln - 1], q))
  need(near(Q(q, tab = tab, key = key), L(lid), tol), sprintf("%s value = %s", lid, q)); invisible(lid)
}

# ======================================================================================================== panel a: the contrasts
# Each cell: value, 95% interval, covariance-aware second-order p (W4a). Intervals: stored (ii: diff_lo/diff_hi, normal; i: Fieller
# ror_lo/ror_hi); where W4a stores a difference and its SE but no interval ((iii), (vi)) the interval is value +/- 1.959964 x SE, derived
# here and logged (no new analysis: arithmetic on stored values).
robust3ii <- function(q) { n <- Q(q, "notes"); need(grepl("-> (TRUE|FALSE)", n), paste("Rule 3(ii) verdict in notes of", q)); grepl("-> TRUE", n) }
cellT <- function(pre) list(d = Q(paste0(pre, ".diff")), se = Q(paste0(pre, ".diff"), "se"))
mk <- function(est, lo, hi, p, ok3ii, ...) data.table(est = est, lo = lo, hi = hi, p = p, dotted = !ok3ii, ...)
# (i) arginine
ri <- mk(Q("S4.arg.marg.ror"), Q("S4.arg.marg.ror_lo"), Q("S4.arg.marg.ror_hi"), Q("S4.arg.marg.p_aware_2nd"), robust3ii("S4.ld.boot.arg.marg.aware2nd.share_same"),
         label = "Arginine, marginal", kind = "primary", pid = "L797")
LQ("L797", "S4.arg.marg.p_aware_2nd"); LQ("L800", "S4.arg.marg.ror"); LQ("L802", "S4.arg.marg.ror_lo"); LQ("L801", "S4.arg.marg.ror_hi")
# (vi) ornithine (503 panel, marginal)
oz <- cellT("S4.orn.503")
rvi <- mk(oz$d, oz$d - Z975 * oz$se, oz$d + Z975 * oz$se, Q("S4.orn.503.p_aware_2nd"), robust3ii("S4.ld.boot.orn.marg.aware2nd.share_same"),
          label = "Ornithine, marginal", kind = "nocoloc", pid = "L807")
need(near(Q("S4.orn.503.kA") - Q("S4.orn.503.kB"), oz$d, 1e-6), "ornithine difference = kA - kB")
need(near(2 * pnorm(-abs(oz$d / oz$se)), rvi$p, 2e-7), "ornithine p from the stored SE")
LQ("L807", "S4.orn.503.p_aware_2nd"); LQ("L804", "S4.coloc_clause.ornithine.A.503.PPH4")
# (ii) protein-scaled T2D: four cells
c2 <- function(res, set, lab, kind, share, pid) {
  q <- sprintf("S4.grid.%s.%s", res, set)
  mk(Q(paste0(q, ".diff")), Q(paste0(q, ".diff_lo")), Q(paste0(q, ".diff_hi")), Q(paste0(q, ".p_aware_2nd")), robust3ii(share), label = lab, kind = kind, pid = pid)
}
rii <- rbind(c2("FinnGen", "cond", "FinnGen, conditional, primary", "primary", "S4.ld.boot.FinnGen.cond.aware2nd.share_same", "L814"),
             c2("FinnGen", "marg", "FinnGen, marginal, beside", "sens", "S4.ld.boot.FinnGen.marg.aware2nd.share_same", "L776"),
             c2("DIAMANTE", "marg", "DIAMANTE, marginal, primary", "primary", "S4.ld.boot.DIAMANTE.marg.aware2nd.share_same", "L808"),
             c2("DIAMANTE", "cond", "DIAMANTE, conditional, beside", "sens", "S4.ld.boot.DIAMANTE.cond.aware2nd.share_same", "L765"))
LQ("L814", "S4.grid.FinnGen.cond.p_aware_2nd"); LQ("L776", "S4.grid.FinnGen.marg.p_aware_2nd"); LQ("L808", "S4.grid.DIAMANTE.marg.p_aware_2nd"); LQ("L765", "S4.grid.DIAMANTE.cond.p_aware_2nd")
LQ("L809", "S4.grid.FinnGen.cond.diff"); LQ("L811", "S4.grid.FinnGen.cond.diff_lo"); LQ("L810", "S4.grid.FinnGen.cond.diff_hi")
LQ("L766", "S4.grid.DIAMANTE.marg.diff"); LQ("L768", "S4.grid.DIAMANTE.marg.diff_lo"); LQ("L767", "S4.grid.DIAMANTE.marg.diff_hi")
need(all(rii$est < 0 & rii$hi < 0.001), "(ii): all four differences negative; only FinnGen marginal's upper limit reaches 0")
# (iii) direct comparison: four cells, equal weight
c3 <- function(res, set, lab, share, pid) {
  q <- sprintf("S4.t3.%s.%s.503", res, set); z <- cellT(q)
  mk(z$d, z$d - Z975 * z$se, z$d + Z975 * z$se, Q(paste0(q, ".2nd.p_aware")), robust3ii(share), label = lab, kind = "primary", pid = pid)
}
riii <- rbind(c3("FinnGen", "cond", "FinnGen, conditional", "S4.ld.boot.t3.FinnGen.cond.aware2nd.share_same", "L818"),
              c3("FinnGen", "marg", "FinnGen, marginal", "S4.ld.boot.t3.FinnGen.marg.aware2nd.share_same", "L819"),
              c3("DIAMANTE", "marg", "DIAMANTE, marginal, beside", "S4.ld.boot.t3.DIAMANTE.marg.aware2nd.share_same", "L788"),
              c3("DIAMANTE", "cond", "DIAMANTE, conditional, beside", "S4.ld.boot.t3.DIAMANTE.cond.aware2nd.share_same", "L787"))
LQ("L818", "S4.t3.FinnGen.cond.503.2nd.p_aware"); LQ("L819", "S4.t3.FinnGen.marg.503.2nd.p_aware")
LQ("L788", "S4.t3.DIAMANTE.marg.503.2nd.p_aware"); LQ("L787", "S4.t3.DIAMANTE.cond.503.2nd.p_aware")
for (i in 1:4) need(near(2 * pnorm(-abs((riii$est[i]) / ((riii$hi[i] - riii$lo[i]) / (2 * Z975)))), riii$p[i], 2e-6), paste("(iii) p from the stored SE, cell", i))
need(all(riii$est > 0), "(iii): all four differences positive")
# sensitivity-to-reference-LD status: bootstrap shares quoted in the ledger agree with the file
LQ("L817", "S4.ld.boot.FinnGen.cond.aware2nd.share_same", 1e-6); LQ("L780", "S4.ld.boot.t3.FinnGen.cond.aware2nd.share_same", 1e-6)
need(near(L("L924"), 0.9), "Rule 3(ii) share threshold 0.90 = L924")
need(identical(rii$dotted, c(TRUE, TRUE, FALSE, FALSE)) && identical(riii$dotted, c(TRUE, FALSE, FALSE, FALSE)) && !ri$dotted && !rvi$dotted,
     "reference-LD sensitivity pattern read from the file: (ii) FinnGen cells; (iii) FinnGen conditional")
# Holm values: the stored Holm file agrees with the ledger rows (L898 FinnGen conditional, the primary-family comparison)
need(near(Q("S9.holm.arg1_t2d.p_holm", tab = HL), L("L898"), 1e-7) && near(Q("S9.holm.arg1_ornithine.p_holm", tab = HL), L("L897"), 1e-9), "Holm file = ledger L898 / L897")
need(near(Q("S9.holm.arg1_t2d.p_unadj", tab = HL), rii$p[1], 1e-9), "Holm 0.22 was computed on the FinnGen conditional p (the first (ii) row)")

# ---- drawing ---------------------------------------------------------------------------------------------------------
PITCH <- 2.7        # mm between 7-pt text lines (a 7-pt line box is 2.76 mm; the gate allows 15% overlap)
NOTE_W <- 70      # notes column (mm); forest column = 170 - NOTE_W
rt <- function(data, ...) geom_richtext(data = data, aes(...), hjust = 0, size = ANNOT, family = BASE_FONT, fill = NA, label.colour = NA,
                                         label.padding = unit(0, "pt"), label.margin = unit(0, "pt"), inherit.aes = FALSE)
PI <- "<i>p</i>"
mkp <- function(txt) paste0(PI, " ", txt)
forest <- function(d, xlim, breaks, trans, ttl, axid, what, px, tag = NULL) {
  n <- nrow(d); d <- copy(d)[, y := n:1]
  labs_x <- reg_ticks(axid, breaks, p31r_lab(breaks), what, if (trans == "log10") "log10" else "linear")
  null <- if (trans == "log10") 1 else 0
  PAD <- 0.045 * diff(if (trans == "log10") log10(xlim) else xlim)
  pw <- function(a, b) if (trans == "log10") 10^(log10(a) + b) else a + b
  p <- ggplot(d, aes(est, y)) +
    annotate("segment", x = null, xend = null, y = 0.5, yend = n + 0.5, linewidth = 0.3, colour = COL_REF, linetype = "dashed") +
    geom_rect(data = d[dotted == TRUE], aes(xmin = pw(lo, -PAD), xmax = pw(hi, PAD), ymin = y - 0.4, ymax = y + 0.4), inherit.aes = FALSE,
              fill = NA, colour = INK, linewidth = 0.3, linetype = "11") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.45, colour = INK, orientation = "y") +
    geom_point(data = d[kind == "primary"], shape = 21, size = 2.1, stroke = 0.45, fill = INK, colour = INK) +
    geom_point(data = d[kind == "sens"], shape = 21, size = 2.1, stroke = 0.6, fill = "white", colour = INK) +
    geom_point(data = d[kind == "nocoloc"], shape = 13, size = 2.6, stroke = 0.5, colour = INK) +
    rt(d[, px_ := px], x = px_, y = y, label = ptxt) +
    scale_x_continuous(trans = if (trans == "log10") "log10" else "identity", breaks = breaks, labels = labs_x, expand = expansion(0)) +
    scale_y_continuous(breaks = d$y, labels = d$label, limits = c(0.5, n + 0.5), expand = expansion(0)) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(x = NULL, y = NULL, title = ttl, tag = tag) +
    theme_p31r() +
    theme(plot.title = element_text(size = ANNOT_PT, face = "bold", hjust = 0, margin = margin(b = 1.5)), plot.title.position = "plot",
          axis.text.y = element_text(size = ANNOT_PT, colour = "black"), axis.text.x = element_text(size = ANNOT_PT, colour = "black"), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.margin = margin(0.5, 15, 0.5, 2, "mm"))
  check_drawn(p, "x", axid)
  p
}
# notes: one richtext line per call, y in fractions of the panel (top to bottom), x in mm of the notes column
notes <- function(lines, pitch_mm, panel_mm, dy = 0) {   # dy (mm, per line): the line carrying a 9-pt base of a power of ten sits higher, the rest lower
  n <- length(lines); d <- data.table(x = 5.5, y = panel_mm - (seq_len(n) - 0.5) * pitch_mm + dy, lab = lines)
  ggplot(d) + geom_richtext(aes(x = x, y = y, label = lab), hjust = 0, size = ANNOT, family = BASE_FONT, fill = NA, label.colour = NA,
                            label.padding = unit(0, "pt"), label.margin = unit(0, "pt")) +
    coord_cartesian(xlim = c(0, NOTE_W), ylim = c(0, panel_mm), expand = FALSE, clip = "off") + theme_void()
}
# ---- (i)
a_p_i  <- LP("L797", NULL, 2, "a (i): arginine ratio-of-ratios p, marginal covariance-aware second-order", anchor = "^p {}$")
ri[, ptxt := mkp(a_p_i)]
fi <- forest(ri, c(0.5, 2), c(0.5, 1, 2), "log10", "(i) Arginine: ratio of ratios, signal A / signal B", "4a_x_i", "a (i), x: ratio of ratios, log scale",
             px = 2.18, tag = "a")
n_i <- c(sprintf("Ratio of ratios %s (95%% CI %s to %s).",
                 LP("L800", NULL, 3, "a (i): ratio of ratios A/B", anchor = "^Ratio of ratios {} "),
                 LP("L802", NULL, 3, "a (i): lower limit", anchor = "CI {} to"),
                 LP("L801", NULL, 3, "a (i): upper limit", anchor = "to {}\\)\\.$")),
         "No difference detected; equivalence not established.")
# ---- (vi)
a_p_vi <- LPS("L807", NULL, "a (ii): ornithine A-minus-B p, covariance-aware second-order", anchor = "^p {}$")
rvi[, ptxt := mkp(a_p_vi)]
fvi <- forest(rvi, c(-1, 0.25), c(-1, -0.5, 0), "linear", "(ii) Ornithine: A − B, SD of natural-log ornithine per SD of ARG1", "4a_x_vi", "a (ii), x: difference per SD of ARG1",
              px = 0.31)
holm_orn <- LPS("L897", NULL, "a (ii): Holm-adjusted ornithine p", anchor = "^Holm-adjusted p {}")
n_vi <- c(sprintf("Holm-adjusted %s %s. Ornithine at A not", PI, holm_orn),
          sprintf("colocalized with ARG1 (PP.H4 %s): the ratio at A", LP("L804", NULL, 2, "a (ii): ornithine at signal A, PP.H4 with the protein, 503 panel", anchor = "PP\\.H4 {}\\)")),
          "is not tied to the protein signal.")
# ---- (ii)
rii[, ptxt := mkp(c(LP("L814", NULL, 3, "a (iii): FinnGen conditional p (primary)", anchor = "^p {}$"), LP("L776", NULL, 3, "a (iii): FinnGen marginal p (sensitivity)", anchor = "^p {}$"),
                    LP("L808", NULL, 3, "a (iii): DIAMANTE marginal p (primary)", anchor = "^p {}$"), LP("L765", NULL, 3, "a (iii): DIAMANTE conditional p (sensitivity)", anchor = "^p {}$")))]
fii <- forest(rii, c(-0.42, 0.06), c(-0.4, -0.2, 0), "linear", "(iii) Protein-scaled T2D: A − B difference in log-odds per SD of ARG1", "4a_x_ii",
              "a (iii), x: difference in T2D log-odds per SD of ARG1", px = 0.085)
# FinnGen marginal (sensitivity) row: its dotted outline is read from its own Rule 3(ii) verdict (stored bootstrap share 0.5942 < 0.90,
# W4a S4.ld.boot.FinnGen.marg.aware2nd.share_same, results.tsv line 334); that share has no ledger row, so the note describes the status
# in words and prints no number for it.
n_ii <- c(sprintf("Holm-adjusted %s %s, FinnGen conditional only.", PI, LP("L898", NULL, 2, "a (iii): Holm-adjusted p, FinnGen conditional only", anchor = "p {},")),
          sprintf("FinnGen conditional: bootstrap share %s against",
                  LP("L817", NULL, 2, "a (iii): bootstrap share, FinnGen conditional T2D", anchor = "share {} against")),
          sprintf("a %s rule (dotted outline); the FinnGen marginal",
                  LP("L924", NULL, 2, "a (iii)/(iv): reference-LD share rule", anchor = "a {} rule")),
          "row is also LD-sensitive.")
# ---- (iii)
riii[, ptxt := mkp(c(LP("L818", NULL, 3, "a (iv): FinnGen conditional direct-comparison p", anchor = "^p {}$"), LP("L819", NULL, 3, "a (iv): FinnGen marginal direct-comparison p", anchor = "^p {}$"),
                     LP("L788", NULL, 3, "a (iv): DIAMANTE marginal direct-comparison p", anchor = "^p {}$"), LP("L787", NULL, 3, "a (iv): DIAMANTE conditional direct-comparison p", anchor = "^p {}$")))]
fiii <- forest(riii, c(-0.1, 0.55), c(0, 0.2, 0.4), "linear", "(iv) Direct comparison: A − B, T2D log-odds per SD of natural-log arginine", "4a_x_iii",
               "a (iv), x: difference in T2D log-odds per SD of natural-log arginine", px = 0.60)
P05 <- "(?:< |below ){}"
n_iii <- c("Did not establish different scaling: the requirement",
           sprintf("(%s < %s under both FinnGen estimands) was not met.", PI, LP("L957", NULL, 2, "a (iv): requirement p < 0.05 under both FinnGen estimands", anchor = P05)),
           "Dotted outline: FinnGen conditional, bootstrap share",
           sprintf("%s against a %s rule. Both DIAMANTE %s values", LP("L780", NULL, 2, "a (iv): bootstrap share, FinnGen conditional direct comparison", anchor = "^{} against"),
                   LP("L924", NULL, 2, "a (iii)/(iv): reference-LD share rule", anchor = "a {} rule"), PI),
           sprintf("are below %s; the tests are not uniformly negative.", LP("L957", NULL, 2, "a (iv): DIAMANTE p values below 0.05", anchor = P05)))
# ---- key (inside the graphic). Line 1: the four status encodings of panel a; line 2: interval convention + the prominent pointer to
# the per-signal estimates (plan F14). Credible-set tick styles and the LD colours are keyed inside panel b.
KEY_MM <- 11.4
reg_lit("0.5", "scale: LD bin boundary of the track colours (r2 at least 0.5)", anchor = "≥ {},", not_ledger = "figure-spec LD bin boundary, not a result")
reg_lit("0.5", "scale: LD bin boundary of the track colours (light bin 0.2 to 0.5)", anchor = "–{},", not_ledger = "figure-spec LD bin boundary, not a result")
reg_lit("0.2", "scale: LD bin boundary of the track colours (light bin 0.2 to 0.5)", anchor = "light {}–", not_ledger = "figure-spec LD bin boundary, not a result")
reg_lit("0.2", "scale: LD bin boundary of the track colours (grey below 0.2)", anchor = "< {}$", not_ledger = "figure-spec LD bin boundary, not a result")
LDK <- data.table(x = c(2, 4.7, 12.4, 15.1), ld = c("A_hi", "A_lo", "B_hi", "B_lo"), y = KEY_MM - 8.9)
kd <- data.table(x = c(2, 52, 78, 127), y = KEY_MM - 3.3, shape = c("primary", "sens", "nocoloc", "dot"),
                 lab = c("filled, primary or direct-comparison cell", "open, beside (iii)", "ARG1–ornithine not colocalized at A", "sensitive to the reference LD"))
key <- ggplot() +
  geom_point(data = kd[shape == "primary"], aes(x, y), shape = 21, size = 2.1, stroke = 0.45, fill = INK, colour = INK) +
  geom_point(data = kd[shape == "sens"], aes(x, y), shape = 21, size = 2.1, stroke = 0.6, fill = "white", colour = INK) +
  geom_point(data = kd[shape == "nocoloc"], aes(x, y), shape = 13, size = 2.6, stroke = 0.5, colour = INK) +
  annotate("rect", xmin = kd[shape == "dot", x] - 2.2, xmax = kd[shape == "dot", x] + 2.2, ymin = KEY_MM - 3.3 - 1.1, ymax = KEY_MM - 3.3 + 1.1, fill = NA, colour = INK, linewidth = 0.3, linetype = "11") +
  geom_text(data = kd, aes(x = x + 3.2, y = y, label = lab), hjust = 0, size = ANNOT, family = BASE_FONT) +
  annotate("segment", x = 0.4, xend = 3.6, y = KEY_MM - 6.1, yend = KEY_MM - 6.1, linewidth = 0.45, colour = INK) +
  annotate("text", x = 5.2, y = KEY_MM - 6.1, hjust = 0, size = ANNOT, family = BASE_FONT,
           label = "pointwise 95% intervals, nominal p; dashed line, no difference") +
  annotate("text", x = 170, y = KEY_MM - 6.1, hjust = 1, size = ANNOT, family = BASE_FONT, fontface = "bold",
           label = "Per-signal estimates: Additional file 3: Figure S5") +
  geom_point(data = LDK, aes(x = x, y = y, colour = ld), size = 0.8 * 2.2, stroke = 0, show.legend = FALSE) + scale_colour_manual(values = COL_LD, guide = "none") +
  annotate("text", x = c(6.5, 16.9), y = KEY_MM - 8.9, label = c("A", "B"), hjust = 0, size = ANNOT, family = BASE_FONT) +
  annotate("text", x = 21, y = KEY_MM - 8.9, hjust = 0, size = ANNOT, family = BASE_FONT,
           label = "panel b, variant color: r² with signal A or B; dark ≥ 0.5, light 0.2–0.5, grey < 0.2") +
  coord_cartesian(xlim = c(0, 170), ylim = c(0, KEY_MM), expand = FALSE, clip = "off") + theme_void()
# ---- grid of panel a. Panel heights (mm) set the row proportions; row furniture (title, axis) is added by the plots
PH <- c(i = 2 * PITCH, vi = 3 * PITCH, ii = 6 * PITCH, iii = 7 * PITCH)
pitch <- c(i = PITCH, vi = PITCH, ii = PITCH, iii = PITCH)
nn <- list(i = notes(n_i, pitch[["i"]], PH[["i"]]), vi = notes(n_vi, pitch[["vi"]], PH[["vi"]], dy = c(0.3, -0.5, -0.5)), ii = notes(n_ii, pitch[["ii"]], PH[["ii"]]), iii = notes(n_iii, pitch[["iii"]], PH[["iii"]]))
natural <- function(p) { g <- patchwork::patchworkGrob(p)
  c(w = grid::convertWidth(sum(g$widths), "mm", valueOnly = TRUE), h = grid::convertHeight(sum(g$heights), "mm", valueOnly = TRUE)) }
# absolute panel sizes: the plots add their own axes and titles around them; the forest panel width is solved so the row is 170 mm wide
build_a <- function(Pw) wrap_plots(list(fi, nn$i, fvi, nn$vi, fii, nn$ii, fiii, nn$iii), ncol = 2, widths = unit(c(Pw, NOTE_W), "mm"), heights = unit(PH, "mm"))
na0 <- natural(build_a(40)); PW_A <- 170 - NOTE_W - (na0[["w"]] - 40 - NOTE_W)
panel_a <- build_a(PW_A); na <- natural(panel_a)
need(abs(na[["w"]] - 170) < 0.05, "panel a is 170 mm wide")

# ======================================================================================================== panel b: compact locus
mlp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
ZOOM <- c(131.48, 131.67)                                  # the enlarged interval (Mb, GRCh38): holds every credible-set member
hz <- h[pos38 >= ZOOM[1] * 1e6 & pos38 <= ZOOM[2] * 1e6]
need(nrow(hz) == 486, "486 harmonised variants in the enlarged interval")
need(md5f(p31r_use(T2D_F)) == "00f930552bdff89d40a921377d944f68", "FinnGen R10 T2D file is the approved version (md5 00f93055); pinned in the render log")
tf <- Rsamtools::TabixFile(T2D_F)                          # its .tbi index is read by Rsamtools; approved Figure_4 row (md5 7c7965b8), pinned below
need(md5f(p31r_use(paste0(T2D_F, ".tbi"))) == "7c7965b8e6794d20e897be42a8666895", "FinnGen T2D .tbi is the approved version")
chr <- if ("6" %in% Rsamtools::seqnamesTabix(tf)) "6" else "chr6"
raw <- Rsamtools::scanTabix(tf, param = GenomicRanges::GRanges(chr, IRanges::IRanges(ZOOM[1] * 1e6, ZOOM[2] * 1e6)))[[1]]
d <- fread(text = paste(raw, collapse = "\n"), header = FALSE, sep = "\t")
setnames(d, c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"))
m <- merge(hz, d[, .(pos38 = pos, ref, alt, bD = beta, sD = sebeta, af_alt, mlogpD = mlogp)], by = "pos38")
m <- m[(alt == a1 & ref == a2) | (alt == a2 & ref == a1)]
m[, flip := alt != a1][, afD := ifelse(flip, 1 - af_alt, af_alt)]
m <- m[palindromic != 1 | (pmin(af_ref, 1 - af_ref) <= 0.40 & pmin(afD, 1 - afD) <= 0.40 & abs(afD - af_ref) <= 0.10)]
m <- m[!duplicated(snp)][order(pos38)]
m[, `:=`(lpP = mlp(z__P__ARG1), lpM = mlp(z__M__GCST90200372), lpD = mlogpD)]
need(near(m$lpD, -log10(2 * pnorm(-abs(m$bD / m$sD))), 0.05), "FinnGen's -log10 p agrees with its beta/SE")
iA <- m[rsid_ref == SIG[["A"]], which = TRUE]; iB <- m[rsid_ref == SIG[["B"]], which = TRUE]
need(length(iA) == 1 && length(iB) == 1, "both signal variants in the T2D-joined set")
POSA <- m$pos38[iA] / 1e6; POSB <- m$pos38[iB] / 1e6
# ---- LD colouring (S2-D25): r2 of every plotted variant with signal A and with signal B, from the 503-person dosage panel (approved Figure_4 row).
# Bins of the figure spec: r2 >= 0.5 (dark signal hue) and 0.2-0.5 (light signal hue); below 0.2 stays grey. The variant takes the hue of the signal
# it is more strongly linked to. The reading is verified three ways: the stamp pins the (harmonised, dosage) pair, the dosage allele frequency
# equals af_ref, and r(A,B) and r2(rs2246012 with A, B) reproduce the ledger rows computed from PLINK.
I8_F <- file.path(V3, "V4_local/raw/l3-finemap/ARG1_6q23__ref_A1dosage_503xP.i8")
PREP_F <- file.path(V3, "V4_local/staging/l3-finemap/ARG1_6q23/PREP_STAMP.json")
stamp <- jsonlite::fromJSON(p31r_use(PREP_F))
i8key <- grep("ref_A1dosage_503xP[.]i8$", names(stamp$outputs), value = TRUE); need(length(i8key) == 1, "the stamp names the 503-person dosage file once")
need(md5f(p31r_use(H_F)) == stamp$outputs[["harmonised.tsv.gz"]] && md5f(p31r_use(I8_F)) == stamp$outputs[[i8key]] && md5f(I8_F) == "548e02d431c56888cd93afa2cf5d0fd6" &&
     stamp$n_snps == nrow(h) && stamp$ref_n == 503, "PREP_STAMP pins this harmonised file and this dosage matrix")
G <- readBin(I8_F, what = "integer", size = 1, signed = TRUE, n = 503 * nrow(h) + 10)
need(length(G) == 503 * nrow(h), "dosage matrix is 503 x 3067")
G <- matrix(G, nrow = 503, byrow = TRUE)
need(all(G %in% 0:2) && max(abs(colMeans(G) / 2 - h$af_ref)) < 1e-5, "dosage values 0-2 and allele frequency = af_ref (column order = harmonised row order)")
hA <- which(h$rsid_ref == SIG[["A"]]); hB <- which(h$rsid_ref == SIG[["B"]]); hC <- which(h$rsid_ref == RS_C)
need(length(hA) == 1 && length(hB) == 1 && length(hC) == 1 && !anyDuplicated(h$snp), "signal and instrument rows are unique in the harmonised panel")
rAll_A <- cor(G, G[, hA])[, 1]; rAll_B <- cor(G, G[, hB])[, 1]
need(round(cor(G[, hA], G[, hB]), 2) == round(L("L058"), 2), "dosage r(A,B) = ledger L058 at its printed precision")
need(near(rAll_A[hC]^2, L("L909"), 1e-5) && near(rAll_B[hC]^2, L("L910"), 1e-5), "dosage r2 of rs2246012 with A and B = ledger L909 / L910 (PLINK)")
rm(G)
ih <- match(m$snp, h$snp); need(!anyNA(ih), "every plotted variant is in the harmonised panel")
LD_HI <- 0.5; LD_LO <- 0.2
m[, `:=`(r2A = rAll_A[ih]^2, r2B = rAll_B[ih]^2)]
m[, `:=`(ldsig = fifelse(r2A >= r2B, "A", "B"), ldmax = pmax(r2A, r2B))]
m[, ld := fifelse(ldmax >= LD_HI, paste0(ldsig, "_hi"), fifelse(ldmax >= LD_LO, paste0(ldsig, "_lo"), "none"))]
need(m$ld[iA] == "A_hi" && m$ld[iB] == "B_hi", "each signal variant is in its own top LD bin")
# rs2246012: actual coordinate, read from the Chen results file (variant.C.pos37_pos38_a1_a2) and checked against the harmonised panel
cr <- CHEN[id == "variant.C.pos37_pos38_a1_a2"]; need(nrow(cr) == 1, "Chen coordinate row")
cp <- strsplit(cr$value, "_")[[1]]; POSC <- as.numeric(cp[2]) / 1e6
need(h[rsid_ref == RS_C, pos38] / 1e6 == POSC && h[rsid_ref == RS_C, pos37] == as.numeric(cp[1]), "rs2246012 coordinate: Chen results = harmonised panel")
need(POSC > POSA && POSC - POSA < 0.002 && POSA - POSB < 0.004, "rs2246012 lies ~1 kb from A and ~4 kb from B")
need(CHEN[id == "r2.C_A.503", value] == CHEN[id == "r2.C_A.503", value] && near(as.numeric(CHEN[id == "r2.C_A.503", value]), L("L909"), 1e-12) &&
     near(as.numeric(CHEN[id == "r2.C_B.503", value]), L("L910"), 1e-12), "Chen r2 (A, B; 503) = ledger L909 / L910")
need(which(CHEN$id == "r2.C_A.503") + 1 == as.integer(sub(".*; line ([0-9]+);.*", "\\1", p31r_ledger_row("L909")$source)) &&
     which(CHEN$id == "r2.C_B.503") + 1 == as.integer(sub(".*; line ([0-9]+);.*", "\\1", p31r_ledger_row("L910")$source)), "ledger L909 / L910 name the Chen r2.C_A / r2.C_B lines")
# ---- gene models (Ensembl release 116, GRCh38)
rec <- jsonlite::fromJSON(file.path(ENS_D, "FETCH_RECORD.json")); UNREG("FETCH_RECORD.json", file.path(ENS_D, "FETCH_RECORD.json"))
ENS_F <- p31r_use(file.path(ENS_D, rec$file))
need(md5f(ENS_F) == rec$md5 && rec$assembly == "GRCh38", "Ensembl response md5 / assembly")
ens <- as.data.table(jsonlite::fromJSON(ENS_F, simplifyDataFrame = TRUE)[, c("feature_type", "id", "Parent", "external_name", "biotype", "start", "end", "strand", "is_canonical", "assembly_name")])
need(all(ens$assembly_name == "GRCh38"), "every feature on GRCh38")
gn <- ens[feature_type == "gene" & biotype == "protein_coding", .(gene = id, symbol = external_name, g0 = start, g1 = end, strand)]
tx <- ens[feature_type == "transcript" & is_canonical %in% 1 & Parent %in% gn$gene, .(tx = id, gene = Parent, t0 = start, t1 = end)]
need(nrow(tx) == nrow(gn) && !anyDuplicated(tx$gene), "one canonical transcript per protein-coding gene")
pcx <- ens[feature_type == "transcript" & biotype == "protein_coding" & Parent %in% gn$gene, .(p0 = min(start), p1 = max(end)), by = .(gene = Parent)]
exn <- ens[feature_type == "exon" & Parent %in% tx$tx, .(tx = Parent, e0 = start, e1 = end)]
gn <- merge(merge(gn, tx, by = "gene"), pcx, by = "gene"); exn <- merge(exn, gn[, .(tx, symbol, strand)], by = "tx")
GZ <- gn[p1 > ZOOM[1] * 1e6 & p0 < ZOOM[2] * 1e6]
need(setequal(GZ$symbol, c("ARG1", "MED23", "ENPP3")), "the enlarged interval holds exactly ARG1, MED23 and ENPP3")
inP <- function(s, p) gn[symbol == s, p0 <= p & p <= p1]
need(all(sapply(c(POSA, POSB, POSC) * 1e6, function(p) inP("ARG1", p) && inP("MED23", p))), "A, B and rs2246012 lie inside the protein-coding spans of ARG1 and MED23")

# ---- credible sets at signal A: members from the 04b file (503 panel / in-sample), status from the canonical status file, per row
cs_members <- function(key) { r <- SETS[quantity_id == sprintf("S4.cs.%s.A.members", key)]; need(nrow(r) == 1, paste("members", key)); strsplit(r$value, ";")[[1]] }
CSK <- data.table(
  mkey  = c("2", "1_protein", "1_T2D", "3", "4", "5"),
  left  = c("SuSiE (in-sample), protein", "SuSiE-RSS, protein", "SuSiE-RSS, T2D", "ABF shared (FinnGen)", "ABF T2D (FinnGen)", "ABF shared (DIAMANTE)"),
  panel = c("UKB-PPP in-sample LD", "503", "503", "503", "503", "503"),
  lidx  = c("m2_insample", "m1_protein_503", "m1_T2D_503", "m3_503", "m4_503", "m5_503"))
CSK[, y := .N:1]
rows_cs <- list(); tick_cs <- list()
for (i in seq_len(nrow(CSK))) {
  k <- CSK[i]; r <- CSS[set_signal == "A" & method == k$mkey & panel == k$panel]
  need(nrow(r) == 3 && uniqueN(r$set_size) == 1 && uniqueN(r$admissibility_rule3i) == 1 && uniqueN(r$robustness_rule3ii) == 1, paste("one status per set:", k$mkey))
  mem <- cs_members(k$mkey); pos <- h$pos38[match(mem, h$rsid_ref)] / 1e6
  need(!anyNA(pos) && length(mem) == r$set_size[1] && length(mem) == L("L986", c("A", k$lidx)), paste("members = status-file size = ledger L986:", k$mkey))
  a <- r$admissibility_rule3i[1]; rb <- r$robustness_rule3ii[1]
  st <- if (a == "admissible") "admissible" else if (a == "inadmissible") "inadmissible" else if (grepl("^NOT robust", rb)) "notrobust" else stop("unclassified status ", k$mkey)
  need((k$mkey %in% c("1_protein", "1_T2D")) == (st == "inadmissible") && (k$mkey == "2") == (st == "admissible") && (k$mkey %in% c("3", "4", "5")) == (st == "notrobust"), paste("status pattern read from the file:", k$mkey))
  rows_cs[[i]] <- data.table(k, status = st, n = length(mem))
  tick_cs[[i]] <- data.table(y = k$y, x = pos, status = st)
}
rows_cs <- rbindlist(rows_cs); tick_cs <- rbindlist(tick_cs)
LQ("L820", "S4.s.max.ARG1_T2D", 1e-9, tab = SETS)
need(all(grepl("0.1482", CSS[method %in% c("1_T2D", "1_protein") & panel == "503", admissibility_text_in_source])), "the status file's inadmissibility text carries s = 0.1482 (= ledger L820 at 3 dp)")
rows_cs[, rtxt := NA_character_]
CSL <- function(status, n, lidx, mkey) {
  nn_ <- LP("L986", c("A", lidx), 0, paste0("b, credible set size, signal A, method ", mkey), anchor = paste0("^{} variants; ", switch(status, admissible = "admissible", inadmissible = "inadmissible", notrobust = "not robust")))
  switch(status,
    admissible = sprintf("%s variants; admissible; spans two genes", nn_),
    inadmissible = sprintf("%s variants; inadmissible fit", nn_),
    notrobust = sprintf("%s variants; not robust to the reference LD", nn_))
}
rows_cs[, rtxt := unlist(mapply(CSL, status, n, lidx, mkey, USE.NAMES = FALSE))]
# S2-D26: s is a PAIR-level statistic (L820: the maximum over every fit, both traits, all seven fits), not a per-fit value, and no ledger row holds the
# per-fit values; it is therefore printed ONCE, in a footnote under the credible-set rows, and never on a row.
cs_foot <- sprintf("Inadmissible fit: maximum s over all fits for the pair %s (limit %s).",
                   LP("L820", NULL, 3, "b, footnote: maximum s over all fits for the ARG1-T2D pair", anchor = "pair {} \\(limit"),
                   LP("L923", NULL, 2, "b, footnote: Rule 3(i) limit on s", anchor = "limit {}\\)\\.$"))

# ---- colouring of tracks: grey = r2 < 0.2 with both signals; light / dark signal hue = r2 0.2-0.5 / >= 0.5 with that signal; diamonds = the signals themselves
# S2-D26: the ARG1 protein track only. The variant set is the T2D-joined set (the same palindrome rule as the signal analysis and as AF3 S5d), so
# the track and S5d show the same variants.
TR <- m[, .(pos = pos38 / 1e6, z = lpP, trk = "P", ld)]
leads <- m[c(iA, iB), .(pos = pos38 / 1e6, z = lpP, trk = "P", sig = c("A", "B"))]
YLIM <- list(P = c(-9, 145))
for (k in names(YLIM)) need(all(TR[trk == k, z] >= YLIM[[k]][1] & TR[trk == k, z] <= YLIM[[k]][2]), paste("every point inside its y range:", k))
YAX <- c(P = "4b_y_protein")
YBR <- list(P = c(0, 50, 100))
LOG_TTL <- "−log<sub><span style='font-size:7pt'>10</span></sub> <i>p</i>"
reg_lit("10", "scale: base of the logarithm in the panel b y title (−log10 p)", anchor = "log {}")
LM <- 38; RM <- 170 - LM - 84        # left / middle (incl. tick labels) / right columns (mm)
trk <- function(k, ttl_y = FALSE, tag = NULL) {
  dd <- TR[trk == k]; ll <- leads[trk == k]
  labs_y <- reg_ticks(YAX[[k]], YBR[[k]], p31r_lab(YBR[[k]]), paste0("b, ", c(P = "ARG1 protein", M = "arginine", D = "T2D (FinnGen)")[[k]], " track, -log10 p"))
  p <- ggplot(dd, aes(pos, z)) +
    geom_point(data = dd[ld == "none"], colour = PTGREY, size = 0.55, stroke = 0) +
    geom_point(data = dd[grepl("_lo$", ld)], aes(colour = ld), size = 0.8, stroke = 0) +
    geom_point(data = dd[grepl("_hi$", ld)], aes(colour = ld), size = 0.8, stroke = 0) +
    scale_colour_manual(values = COL_LD, guide = "none") +
    geom_point(data = ll, aes(fill = sig), shape = 23, size = 2.1, stroke = 0.35, colour = "black", show.legend = FALSE) +
    scale_fill_manual(values = COL_SIG) +
    scale_x_continuous(limits = ZOOM, expand = expansion(0)) +
    scale_y_continuous(limits = YLIM[[k]], breaks = YBR[[k]], labels = labs_y, expand = expansion(0)) +
    labs(x = NULL, y = if (ttl_y) LOG_TTL else NULL) + theme_p31r() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_markdown(),
          plot.margin = margin(0.4, 1, 0.4, 0, "mm"))
  check_drawn(p, "y", YAX[[k]])
  p
}
tP <- trk("P", TRUE)
# text columns (x in mm; y in fractions of the row, 0 bottom, 1 top; lines centred on the row)
txtcol <- function(lines, wmm, bold = FALSE, cols = NULL, pitch_f = NULL, h_mm = 11, ys = NULL) {
  n <- length(lines); pf <- PITCH / h_mm
  y <- if (is.null(ys)) 0.5 + ((n - 1) / 2 - (seq_len(n) - 1)) * pf else ys
  d <- data.table(x = 0.8, y = y, lab = lines, col = if (is.null(cols)) "black" else cols)
  ggplot(d) + geom_richtext(aes(x = x, y = y, label = lab, colour = col), hjust = 0, size = ANNOT, family = BASE_FONT, fill = NA, label.colour = NA,
                            label.padding = unit(0, "pt"), label.margin = unit(0, "pt")) + scale_colour_identity() +
    coord_cartesian(xlim = c(0, wmm), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void()
}
H_TRK <- 15.0    # the protein track alone, taller than the three tracks it replaces (S2-D26)
lt_P <- txtcol("ARG1 protein", LM, h_mm = H_TRK)
r_txt <- paste("A and B: r =", LP("L058", NULL, 2, "b, r between the two signal variants (1000 Genomes EUR, n = 503)", anchor = "r = {}"))
r_ln <- c("", paste(r_txt, "(resampling interval"),
          sprintf("%s to %s; the sign is uncertain)", LP("L988", 1L, 2, "b, resampling interval of r, lower limit", anchor = "^{} to"), LP("L988", 2L, 2, "b, resampling interval of r, upper limit", anchor = "to {};")),
          "<b>Arginine and T2D regional tracks:</b>", "<b>Additional file 3: Figure S5, panel d</b>")
R_OFF <- c(2.3, 1.1, 0.1, -1.1, -2.1)
rP <- txtcol(r_ln, RM, h_mm = H_TRK, ys = 0.5 + R_OFF * PITCH / H_TRK) +
  geom_point(data = data.table(x = c(1.5, 22.5), y = 0.5 + R_OFF[1] * PITCH / H_TRK, s = c("A", "B")), aes(x, y, fill = s), shape = 23, size = 2.1, stroke = 0.35, colour = "black", show.legend = FALSE) +
  scale_fill_manual(values = COL_SIG) +
  geom_text(data = data.table(x = c(3.4, 24.4), y = 0.5 + R_OFF[1] * PITCH / H_TRK, lab = c(paste("A", SIG[["A"]]), paste("B", SIG[["B"]]))), aes(x, y, label = lab), hjust = 0, size = ANNOT, family = BASE_FONT)
# gene track
XMW <- 65.8                                                 # width (mm) of the middle column's plotting panel; asserted against the solved layout below
XM <- function(mm) ZOOM[1] + mm / XMW * diff(ZOOM)          # mm along the panel to Mb (positions of the tick key in the credible-set header row)
YG <- c(plus = 1.0, minus = -0.9)
gl <- copy(GZ)[, yy := fifelse(strand > 0, YG[["plus"]], YG[["minus"]])]
gl[, `:=`(c0 = pmax(p0, ZOOM[1] * 1e6), c1 = pmin(p1, ZOOM[2] * 1e6), k0 = pmax(t0, ZOOM[1] * 1e6), k1 = pmin(t1, ZOOM[2] * 1e6))]
ge <- merge(exn, gl[, .(tx, yy)], by = "tx")[e1 > ZOOM[1] * 1e6 & e0 < ZOOM[2] * 1e6]
gl[, `:=`(lx = fifelse(symbol == "ARG1", (t0 + t1) / 2e6, (c0 + c1) / 2e6), ly = fifelse(strand > 0, yy + 0.5, yy - 0.5), vj = fifelse(strand > 0, 0, 1))]
pG <- ggplot() +
  geom_segment(data = gl, aes(x = c0 / 1e6, xend = c1 / 1e6, y = yy, yend = yy), linewidth = 0.35, colour = "grey55") +
  geom_segment(data = gl[k1 > k0], aes(x = k0 / 1e6, xend = k1 / 1e6, y = yy, yend = yy), linewidth = 0.35, colour = "grey10") +
  geom_rect(data = ge, aes(xmin = pmax(e0, ZOOM[1] * 1e6) / 1e6, xmax = pmin(e1, ZOOM[2] * 1e6) / 1e6, ymin = yy - 0.28, ymax = yy + 0.28), fill = "grey10", colour = "grey10", linewidth = 0.15) +
  geom_text(data = gl, aes(x = lx, y = ly, label = symbol, vjust = vj), size = ANNOT, family = BASE_FONT, fontface = "italic") +
  scale_x_continuous(limits = ZOOM, expand = expansion(0)) + scale_y_continuous(limits = c(-2.8, 3.1), expand = expansion(0)) +
  labs(x = NULL, y = NULL) + theme_p31r() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(0.4, 1, 0.4, 0, "mm"))
lt_G <- txtcol("Genes (Ensembl)", LM, h_mm = 9)
rG <- txtcol(c("grey, protein-coding extent;", "black, canonical transcript"), RM, h_mm = 9)
# credible sets
NR <- nrow(CSK)
lw_cs <- 0.5
MKX <- XM(2.0)                                     # x (Mb) of the status marker of a credible-set row, at the left edge of the panel (S5c: one marker per row)
pCS <- ggplot() +
  geom_hline(yintercept = seq_len(NR), colour = "grey88", linewidth = 0.3) +
  geom_segment(data = tick_cs, aes(x = x, xend = x, y = y - 0.3, yend = y + 0.3), linewidth = 0.5, colour = INK) +        # members: plain ticks for every status
  geom_point(data = rows_cs[status == "admissible"], aes(x = MKX, y = y), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.4) +
  geom_hatchdot(data = rows_cs[status == "inadmissible"], aes(x = MKX, y = y), r_mm = 1.1, ink = INK) +
  geom_dotring(data = rows_cs[status == "notrobust"], aes(x = MKX, y = y), r_mm = 1.0, ink = INK) +
  # key of the three status markers, in the header row (each row also carries its status in words)
  annotate("point", x = XM(1.5), y = NR + 1.05, shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.4) +
  annotate("text", x = XM(3.7), y = NR + 1.05, label = "admissible", hjust = 0, size = ANNOT, family = BASE_FONT) +
  geom_hatchdot(data = data.frame(x = XM(19), y = NR + 1.05), aes(x = x, y = y), r_mm = 1.1, ink = INK) +
  annotate("text", x = XM(21.3), y = NR + 1.05, label = "inadmissible fit", hjust = 0, size = ANNOT, family = BASE_FONT) +
  geom_dotring(data = data.frame(x = XM(43), y = NR + 1.05), aes(x = x, y = y), r_mm = 1.0, ink = INK) +
  annotate("text", x = XM(45.2), y = NR + 1.05, label = "not robust", hjust = 0, size = ANNOT, family = BASE_FONT) +
  scale_x_continuous(limits = ZOOM, expand = expansion(0)) + scale_y_continuous(limits = c(0.45, NR + 1.55), expand = expansion(0)) +
  labs(x = NULL, y = NULL) + theme_p31r() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(0.4, 1, 0.4, 0, "mm"))
fy <- function(y) (y - 0.45) / (NR + 1.55 - 0.45)
lt_CS <- txtcol(c("Credible sets at signal A", rows_cs$left), LM, ys = fy(c(NR + 1.05, rows_cs$y)), h_mm = 20.5)
rt_CS <- txtcol(c("gene attribution at signal A unresolved", rows_cs$rtxt), RM, ys = fy(c(NR + 1.05, rows_cs$y)), h_mm = 20.5)
# rs2246012 row + x axis
XB <- c(131.50, 131.55, 131.60, 131.65)
labs_xb <- reg_ticks("4b_x", XB, fmt_mb(XB, 2), "b, x: chromosome 6 position (Mb, GRCh38)")
pRS <- ggplot() +
  annotate("segment", x = POSC, xend = ZOOM[2] + 0.0, y = 0.5, yend = 0.5, linewidth = 0.3, colour = INK) +
  geom_point(aes(x = POSC, y = 0.5), shape = 24, size = 2.4, stroke = 0.45, fill = INK, colour = INK) +
  annotate("text", x = POSC - 0.0045, y = 0.5, label = RS_C, hjust = 1, size = ANNOT, family = BASE_FONT) +
  scale_x_continuous(limits = ZOOM, breaks = XB, labels = labs_xb, expand = expansion(0)) + scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
  labs(x = "Chromosome 6 position (Mb, GRCh38)", y = NULL) + theme_p31r() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.x = element_text(size = ANNOT_PT),
        axis.title.x = element_text(margin = margin(t = 1)), plot.margin = margin(0.4, 1, 0.4, 0, "mm")) +
  coord_cartesian(clip = "off")
check_drawn(pRS, "x", "4b_x")
r_rs <- c("Published instrument; strongly tags the",
          sprintf("specified signal-A lead (r² %s with A,", LP("L909", NULL, 3, "b, rs2246012 r2 with signal A, 503 panel", anchor = "r² {} with")),
          sprintf("%s with B); tagging is not proof of", LP("L910", NULL, 4, "b, rs2246012 r2 with signal B, 503 panel", anchor = "^{} with B")),
          "an identical causal signal.")
rRS <- txtcol(r_rs, RM, h_mm = 10.8, ys = 0.5 + c(1.5, 0.5, -0.5, -1.5) * PITCH / 10.8)
lt_RS <- txtcol("Published arginine instrument", LM, h_mm = 10.8)
# tag b on the first left cell
lt_P <- lt_P + labs(tag = "b") + theme(plot.tag = element_text(size = TAG_PT, face = "bold", family = BASE_FONT))
FN_MM <- 3.4
BLANK <- ggplot() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()   # an empty cell (plot_spacer() changed the solved column widths)
lt_FN <- txtcol(cs_foot, LM, h_mm = FN_MM)       # one line under the credible-set rows; it runs on across the empty middle cell
HB <- c(P = H_TRK, G = 9, CS = 20.5, FN = FN_MM, RS = 10.8)
build_b <- function(Mw) wrap_plots(list(lt_P, tP, rP, lt_G, pG, rG, lt_CS, pCS, rt_CS, lt_FN, BLANK, BLANK, lt_RS, pRS, rRS),
                                    ncol = 3, widths = unit(c(LM, Mw, RM), "mm"), heights = unit(HB, "mm"))
nb0 <- natural(build_b(60)); MW_B <- 170 - LM - RM - (nb0[["w"]] - LM - 60 - RM)
panel_b <- build_b(MW_B); nb <- natural(panel_b)
need(abs(nb[["w"]] - 170) < 0.05, "panel b is 170 mm wide")
need(abs(MW_B - XMW) < 1.0, sprintf("middle panel is %.1f mm: set XMW to it (tick-key positions)", MW_B))

# ======================================================================================================== compose
HA_PANELS <- sum(PH)
HEIGHTS <- c(na[["h"]], KEY_MM, nb[["h"]] + 1)
Z0 <- theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = NA, colour = NA))
fig <- wrap_plots(list(wrap_elements(full = panel_a) + Z0, wrap_elements(full = key) + Z0, wrap_elements(full = panel_b) + Z0), ncol = 1, heights = unit(HEIGHTS, "mm")) + plot_annotation(theme = theme(plot.margin = margin(3.2, 0, 0, 0, "mm")))
nf <- natural(fig); H_FIG <- round(nf[["h"]], 1)
cat(sprintf("outer layout natural size %.1f x %.1f mm (sum of cells %.1f)
", nf[["w"]], nf[["h"]], sum(HEIGHTS)))
cat(sprintf("layout: panel a %.1f x %.1f mm (forest panel %.1f mm), key %.1f, panel b %.1f x %.1f mm (middle %.1f mm); image height %.1f mm
", na[["w"]], na[["h"]], PW_A, KEY_MM, nb[["w"]], nb[["h"]], MW_B, H_FIG))

# ======================================================================================================== facts (stdout)
cat(sprintf("F4 facts: %d variants in the enlarged interval joined to FinnGen R10 T2D (palindrome rule as the signal analysis); rs2246012 at %.6f Mb (A %.6f, B %.6f)\n",
            nrow(m), POSC, POSA, POSB))
print(rbind(cbind(f = "i", ri[, .(label, est = round(est, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 4), dotted)]),
            cbind(f = "vi", rvi[, .(label, est = round(est, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 4), dotted)]),
            cbind(f = "ii", rii[, .(label, est = round(est, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 4), dotted)]),
            cbind(f = "iii", riii[, .(label, est = round(est, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 4), dotted)])))
cat("LD classes of the plotted variants (r2 with A or B, 503 dosage panel): ", paste(names(table(m$ld)), as.integer(table(m$ld)), collapse = "; "), "
", sep = "")
cat("derived intervals (value +/- 1.959964 x stored SE): (vi) ornithine, (iii) all four direct-comparison cells; stored intervals: (i) Fieller, (ii) normal\n")
print(rows_cs[, .(mkey, left, status, n)])
cat(sprintf("  Ensembl release %s, %s\n", rec$release, rec$assembly))

# ======================================================================================================== plotted frames (what the panels draw)
# One frame per data table a panel plots, written beside the figure as af3/source_data/Figure_4__<name>.csv and listed in the registry under `frames`
# (the file naming and registry key of the AF3 units). Only drawn values: the LD field is the BINNED class of the track colours (A_hi, A_lo, B_hi,
# B_lo, none), never an r2; no individual-level value; no value of a trace that is not drawn (the arginine and T2D tracks are in Additional file 3).
fr_forest <- rbindlist(list(i = ri, vi = rvi, ii = rii, iii = riii), use.names = TRUE, idcol = "estimand")
fr_forest[, row := seq_len(.N), by = estimand]                                       # row = position in its forest, top to bottom
fr_forest <- fr_forest[, .(estimand, row, label, kind, est, lo, hi, p, ref_ld_sensitive = dotted, ledger_p = pid,
                           printed_p = trimws(gsub("<[^>]+>", "", ptxt)))]
need(nrow(TR) == nrow(m) && identical(TR$pos, m$pos38 / 1e6), "the plotted track is the T2D-joined variant set, in its order")
fr_track <- data.table(rsid = m$rsid_ref, pos_mb = TR$pos, mlogp = TR$z, ld = TR$ld)
fr_leads <- leads[, .(signal = sig, rsid = unname(SIG[sig]), pos_mb = pos, mlogp = z)]
fr_cs_rows <- rows_cs[, .(plot_row = y, method_key = mkey, method = left, ld_panel = panel, status, n_members = n, printed_text = rtxt)]
fr_cs_members <- merge(tick_cs, rows_cs[, .(y, method_key = mkey)], by = "y")[order(-y, x), .(method_key, pos_mb = x, status)]
fr_genes <- gl[, .(symbol, strand, plot_y = yy, protein_coding_start_mb = c0 / 1e6, protein_coding_end_mb = c1 / 1e6,
                   canonical_start_mb = k0 / 1e6, canonical_end_mb = k1 / 1e6, canonical_drawn = k1 > k0)]
fr_exons <- ge[, .(symbol, plot_y = yy, exon_start_mb = pmax(e0, ZOOM[1] * 1e6) / 1e6, exon_end_mb = pmin(e1, ZOOM[2] * 1e6) / 1e6)][order(symbol, exon_start_mb)]
fr_rs <- data.table(rsid = RS_C, pos_mb = POSC)
FRAMES <- list(`4a_forest` = fr_forest, `4b_protein_track` = fr_track, `4b_signal_leads` = fr_leads, `4b_credible_set_rows` = fr_cs_rows,
               `4b_credible_set_members` = fr_cs_members, `4b_genes` = fr_genes, `4b_exons` = fr_exons, `4b_rs2246012` = fr_rs)
need(!anyDuplicated(names(FRAMES)) && all(vapply(FRAMES, function(d) is.data.frame(d) && nrow(d) > 0, TRUE)), "every plotted frame is a non-empty table with its own name")
need(nrow(fr_forest) == nrow(ri) + nrow(rvi) + nrow(rii) + nrow(riii) && nrow(fr_cs_members) == nrow(tick_cs) && nrow(fr_cs_rows) == NR, "frame row counts equal the plotted rows")
need(!any(c("r2A", "r2B", "ldmax") %in% unlist(lapply(FRAMES, names))), "no r2 column in any frame (binned LD class only)")
dir.create(file.path(P31R_FIGS, "af3", "source_data"), showWarnings = FALSE)
for (.fn in names(FRAMES)) data.table::fwrite(FRAMES[[.fn]], file.path(P31R_FIGS, "af3", "source_data", sprintf("%s__%s.csv", FIG, .fn)))
REG$frames <- I(names(FRAMES))        # read by the figure gate: every listed frame must exist beside the figure
p31r_save(fig, FIG, w_mm = 170, h_mm = H_FIG)

# registry, written from the same calls that printed the numbers (spec section 3)
dir.create(file.path(P31R_FIGS, "registry"), showWarnings = FALSE)
writeLines(jsonlite::toJSON(REG, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(P31R_FIGS, "registry", paste0(FIG, ".json")), useBytes = TRUE)
message("registry written: ", length(REG$ledger), " ledger, ", length(REG$literals), " literal entries")
