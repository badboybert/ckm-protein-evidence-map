# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S5: ARG1 detail (page 1: per-signal estimates and the posterior of colocalization with ARG1;
# page 2, FigureS5_p2: the complete credible-set method grid at signals A and B, and the 1-Mb overview). Illustrates SM sections S-R9 and S-R13.
# One generator, two page units (FigureS5, FigureS5_p2); continuation pages instead of shrinking (decision F9, S2-D24).
#
#   a  PER-SIGNAL ESTIMATES (moved out of Figure 4 for inferential focus and space, F14; no value is suppressed: the plotted frame is written to
#      af3/source_data). Header (the reviewer's wording): "Per-signal estimates provide context; formal A-B comparisons are shown in Figure 4a."
#      Four estimands, each on its own axis in its own units (no shared axis, no shared null): arginine per SD of ARG1; ornithine per SD of ARG1;
#      T2D per SD of ARG1 (log-odds); T2D per SD of natural-log arginine (log-odds, direct comparison). Signal A and B points with pointwise 95% intervals
#      (Wald, estimate +/- 1.96 SE from the stored SE; not adjusted for multiplicity). The 503-person-panel reading is drawn; the 489-panel
#      readings are in the source-data frame. THE A-BASED REFERENCE FOR B (F7) is drawn in the units of each axis: a band holding signal A's
#      estimate and its 95% interval, behind the signal-B row of the same setting; it assumes the same ratio at both signals. No expected z value
#      is drawn on an effect axis. Visual status: filled = primary cell, open = beside (comparison) or sensitivity (DIAMANTE with the SE recovered
#      from p); crossed circle = "not colocalized with ARG1" (the ornithine association at signal A, whose clause fires).
#   b  POSTERIOR OF COLOCALIZATION with ARG1, from the earlier signal-specific grid (S1R1_arg1_coloc.tsv; ledger L057 and L265 pin it) and the East Asian
#      census/CKB-free rows. Default prior = circle, conservative prior = triangle. NOT the values of the ledger rows L051/L052: the whole-locus T2D values of this grid
#      differ from the whole-locus values printed in the Results (ledger L051); both are printed from their own sources.
#   c  (page 2) credible sets at signals A and B, one row per method, read cell by cell from credible_set_status.tsv. Admissible (SuSiE with in-sample LD) and
#      robust cells = filled; inadmissible SuSiE fit = crossed circle + the word; ABF set not robust to the reference LD = dotted outline + the word.
#   d  (page 2) the 1-Mb overview: NOT DRAWN while its sources are unapproved for this figure (explicit placeholder).
# Sources: ONLY through the approved manifest (AF3_S5). The table of 64 per-signal values is the Prep frame figures/data/Figure_S5/arg1_signal_values_W4a.tsv.
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS5.R"

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
source(paste0(.P31R_V3, "/R2_rewrite/figures/af3/P31R_af3_common.R"))
source(paste0(.P31R_V3, "/R2_rewrite/figures/af3/P31R_af3_glyphs.R"))
suppressPackageStartupMessages({ library(ggtext) })
EN <- "\u2013"
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
INK <- "grey15"; LINE_GREY <- "grey45"; N_ <- PAL_NEUTRAL
COL_SIG <- PAL_SIGNAL; COL_SIG_L <- PAL_SIGNAL_LIGHT
TXT <- function(...) geom_text(..., size = ANNOT, family = BASE_FONT)
fmt0 <- function(v) formatC(v, format = "f", digits = 0)
fmtd <- function(d) function(v) formatC(v, format = "f", digits = d)
# ticks: a narrow axis may join its labels into one printed line, so the anchor allows other numbers on the line (as Figure 2b's generator did)
TICK_ANCHOR <- "^(?:[\u2212+-]?[0-9.,]+\\s)*{}(?:\\s[\u2212+-]?[0-9.,]+)*$"
ticks_set <- function(axis, breaks, lims, why, transform = "linear", labels = p31r_lab(breaks)) {
  stopifnot(all(breaks > lims[1] & breaks < lims[2]), !anyDuplicated(breaks))
  labs <- labels
  for (t in labs) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis, transform = transform, count = 1L, anchor = TICK_ANCHOR)
  labs
}
mm <- function(x, w = 165) x / w
unit_reset <- function() { rm(list = ls(.p31r_read_set), envir = .p31r_read_set); invisible(NULL) }   # each page unit pins only the files IT read


# a dotted outline for "not robust to the reference LD" (visual-status rule 5): geom_point has no line type, so the ring is a path in data units
ring <- function(cx, cy, r_mm, ux, uy, n = 60) {
  th <- seq(0, 2 * pi, length.out = n)
  rbindlist(lapply(seq_along(cx), function(i) data.table(g = i, x = cx[i] + r_mm * ux * cos(th), y = cy[i] + r_mm * uy * sin(th))))
}
# =========================================================================================================== PAGE 1
UNIT <- "FigureS5"
af3_begin(UNIT)
fr <- af3_read("R2_rewrite/figures/data/Figure_S5/arg1_signal_values_W4a.tsv", "a", sep = "\t")
co <- af3_read("R1_stage1/S1R1_01_arg1_triangle/S1R1_arg1_coloc.tsv", "b", sep = "\t")
ecen <- af3_read("V4_local/staging/d2-eas-census/results_census.tsv", "b", sep = "\t")
eckf <- tryCatch(af3_read("V4_local/staging/d2-eas-census/results_ckbfree.tsv", "b", sep = "\t"), error = function(e) NULL)   # not approved for AF3_S5 (S6 only)

# ---- a: per-signal estimates ---------------------------------------------------------------------------------------------------------------------
stopifnot(nrow(fr) == 64L, !anyNA(fr$value), !anyNA(fr$se), all(abs(fr$ci_hi - (fr$value + qnorm(0.975) * fr$se)) < 1e-6))
f3 <- fr[panel == "503"]
stopifnot(nrow(f3) == 32L, setequal(unique(f3$signal), c("A", "B")), setequal(unique(f3$estimand_id), c("i", "ii", "iii", "vi")))
# stored values reproduce the ledger rows the Results print (spot checks of the four estimands)
fv <- function(e, r, s, sg) f3[estimand_id == e & resource == r & setting == s & signal == sg]$value
stopifnot(abs(fv("i", "CLSA arginine", "marginal", "A") - num("L053", 1)) < 1e-6, abs(fv("i", "CLSA arginine", "marginal", "B") - num("L053", 2)) < 1e-6,
          abs(fv("ii", "FinnGen", "conditional", "A") - num("L054", 1)) < 1e-6, abs(fv("ii", "FinnGen", "conditional", "B") - num("L054", 3)) < 1e-6,
          abs(f3[estimand_id == "ii" & resource == "FinnGen" & setting == "conditional" & signal == "A"]$se - num("L054", 2)) < 1e-6,
          abs(fv("ii", "DIAMANTE", "marginal", "A") - num("L055", 1)) < 1e-6, abs(fv("ii", "DIAMANTE", "marginal", "B") - num("L055", 2)) < 1e-6)
# the ornithine ledger values (L805 = A, marginal; B in the next row of the Results) and the clause value at A (L804, 503 panel)
orn_a <- f3[estimand_id == "vi" & setting == "marginal" & signal == "A"]
stopifnot(abs(orn_a$value - num("L805")) < 1e-6, abs(orn_a$coloc_clause_pph4_conditional_vs_ARG1_protein - num("L804")) < 1e-9)
# clause bookkeeping: values below the colocalization rule exist only for ornithine at A (the clause fires) and T2D at B in FinnGen (T2D has no association
# at B, so the clause does not fire there: pre-declaration Rule 3(iv); it is therefore NOT drawn crossed)
low <- f3[!is.na(coloc_clause_pph4_conditional_vs_ARG1_protein) & coloc_clause_pph4_conditional_vs_ARG1_protein < num("L236")]
stopifnot(setequal(unique(paste(low$estimand_id, low$signal, low$resource)), c("vi A CLSA ornithine", "ii B FinnGen")))
f3[, fires := estimand_id == "vi" & signal == "A"]
# status per TEST (S2-D25 visual check: each facet carries its own row labels, so label and marker cannot disagree):
#   (i) arginine, (vi) ornithine: marginal = primary, conditional = beside;
#   (ii) T2D per SD of ARG1: FinnGen conditional and DIAMANTE marginal = primary; FinnGen marginal and DIAMANTE conditional = beside (as Figure 4a);
#   (iii) direct comparison: BOTH FinnGen estimands are co-primary (the test requires p < 0.05 under both), DIAMANTE rows are beside;
#   DIAMANTE with the SE recovered from p = sensitivity in (ii) and (iii). Filled = primary or co-primary; open = beside or sensitivity.
f3[, stat := fcase(resource == "DIAMANTE_SEfromP", "sensitivity",
                   estimand_id %in% c("i", "vi"), fifelse(setting == "marginal", "primary", "beside"),
                   estimand_id == "ii", fifelse((resource == "FinnGen" & setting == "conditional") | (resource == "DIAMANTE" & setting == "marginal"), "primary", "beside"),
                   estimand_id == "iii", fifelse(resource == "FinnGen", "co-primary", "beside"))]
stopifnot(!anyNA(f3$stat))
f3[, `:=`(prim = stat %in% c("primary", "co-primary"), sens = resource == "DIAMANTE_SEfromP")]
cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
f3[, row := fifelse(resource %in% c("CLSA arginine", "CLSA ornithine"), paste0(cap1(setting), "
", stat),
             fifelse(sens, paste0("DIAMANTE, SE from p,
", setting, "; sensitivity"), paste0(resource, " ", setting, "
", stat)))]
f3[, rank := fcase(resource %in% c("CLSA arginine", "CLSA ornithine") & setting == "marginal", 1L, resource %in% c("CLSA arginine", "CLSA ornithine"), 2L,
                   resource == "FinnGen" & setting == "conditional", 1L, resource == "FinnGen", 2L,
                   resource == "DIAMANTE" & setting == "marginal", 3L, resource == "DIAMANTE", 4L,
                   setting == "marginal", 5L, default = 6L)]
ROWORD <- lapply(split(unique(f3[, .(estimand_id, row, rank)]), by = "estimand_id", keep.by = TRUE), function(d) d[order(rank)]$row)
stopifnot(identical(lengths(ROWORD)[c("i", "vi", "ii", "iii")], c(i = 2L, vi = 2L, ii = 6L, iii = 6L)),
          sum(f3[estimand_id == "iii"]$stat == "co-primary") == 4L, sum(f3[estimand_id == "ii"]$prim) == 4L)
f3[, band := signal == "B"]

mk_est <- function(e, ttl, xlab, axis_id, show_y, extra = NULL) {
  d <- copy(f3[estimand_id == e]); rows <- ROWORD[[e]]
  d[, yi := match(row, rows)]; d[, y := -yi + fifelse(signal == "A", 0.17, -0.17)]
  A <- d[signal == "A"][, .(row, yi, a_val = value, a_lo = ci_lo, a_hi = ci_hi)]
  rngx <- range(c(d$ci_lo, d$ci_hi)); pad <- diff(rngx) * 0.08; lims <- c(rngx[1] - pad, rngx[2] + pad)
  br <- scales::extended_breaks(n = 4)(lims); br <- br[br > lims[1] & br < lims[2]]
  p <- ggplot(d, aes(value, y)) +
    geom_rect(data = A, aes(xmin = a_lo, xmax = a_hi, ymin = -yi - 0.42, ymax = -yi + 0.42), inherit.aes = FALSE, fill = COL_SIG_L[["A"]], alpha = 0.45, colour = NA) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = LINE_GREY) +
    geom_segment(aes(x = ci_lo, xend = ci_hi, y = y, yend = y, colour = signal), linewidth = 0.45) +
    geom_point(data = d[prim == TRUE & fires == FALSE], aes(colour = signal, fill = signal), shape = 21, size = 1.9, stroke = 0.45) +
    geom_point(data = d[prim == FALSE & fires == FALSE], aes(colour = signal), shape = 21, fill = "white", size = 1.9, stroke = 0.6) +
    geom_point(data = d[fires == TRUE], aes(colour = signal), shape = 13, size = 2.3, stroke = 0.55) +
    scale_colour_manual(values = COL_SIG, guide = "none") + scale_fill_manual(values = COL_SIG, guide = "none") +
    scale_x_continuous(breaks = br, labels = ticks_set(paste0("S5_", axis_id, "_x"), br, lims, paste0("S5 ", axis_id, " x, ", xlab)), expand = expansion(0)) +
    scale_y_continuous(breaks = -seq_along(rows), labels = if (show_y) rows else rep("", length(rows)), expand = expansion(0)) +
    coord_cartesian(xlim = lims, ylim = c(-length(rows) - 0.5, if (e == "vi") 0.25 else -0.5), clip = "off") +
    labs(x = xlab, y = NULL, title = ttl) +
    theme_p31r() +
    theme(axis.ticks.y = element_blank(), plot.title = element_text(size = ANNOT_PT, face = "bold", hjust = 0, margin = margin(b = 2)), plot.title.position = "panel",
          axis.title.x = element_text(size = ANNOT_PT, margin = margin(t = 2)), axis.text.y = element_text(size = ANNOT_PT, lineheight = 0.9))
  if (!is.null(extra)) p <- p + extra(d, lims)
  p
}
# the clause label (ornithine at signal A): the printed PP.H4 is ledger L804, the 503-person panel value
orn_extra <- function(d, lims) {
  list(annotate("text", x = lims[1] + 0.02 * diff(lims), y = -0.1, hjust = 0, vjust = 0.5, size = ANNOT, lineheight = 0.9, family = BASE_FONT, colour = COL_SIG[["A"]],
                label = paste0("Ornithine at A: not colocalized\nwith ARG1 (PP.H4 = ", led("L804", NULL, fmtd(2), where = "S5a ornithine: signal A conditional PP.H4 (503 panel)", anchor = "^with ARG1 \\(PP\\.H4 = {}\\)$"), ")")))
}
pa1 <- mk_est("i",  "Arginine per SD of ARG1", "SD of natural-log arginine per SD of ARG1", "i", TRUE)
pa2 <- mk_est("vi", "Ornithine per SD of ARG1", "SD of natural-log ornithine per SD of ARG1", "vi", FALSE, orn_extra)
pa3 <- mk_est("ii", "T2D per SD of ARG1\n(protein-scaled)", "T2D log-odds per SD of ARG1", "ii", TRUE)
pa4 <- mk_est("iii", "T2D per SD of natural-log arginine\n(direct comparison)","T2D log-odds per SD of natural-log arginine", "iii", TRUE)
HDR <- paste0("Per-signal estimates provide context; formal A", EN, "B comparisons are shown in Figure 4a.")
p_hdr <- ggplot() + annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, fontface = "bold", label = HDR) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
# key strip for a (own strip; x in mm / 165); rsIDs are read from the per-signal frame
SIGA <- unique(f3[signal == "A"]$signal_variant); SIGB <- unique(f3[signal == "B"]$signal_variant)
stopifnot(length(SIGA) == 1L, length(SIGB) == 1L, grepl("^rs[0-9]+$", c(SIGA, SIGB)))
ka <- ggplot() +
  geom_point(aes(x = mm(3), y = 3.3), shape = 21, fill = COL_SIG[["A"]], colour = COL_SIG[["A"]], size = 1.9, stroke = 0.45) + TXT(aes(x = mm(5.5), y = 3.3, label = paste0("Signal A (", SIGA, ")")), hjust = 0) +
  geom_point(aes(x = mm(38), y = 3.3), shape = 21, fill = COL_SIG[["B"]], colour = COL_SIG[["B"]], size = 1.9, stroke = 0.45) + TXT(aes(x = mm(40.5), y = 3.3, label = paste0("Signal B (", SIGB, ")")), hjust = 0) +
  geom_point(aes(x = mm(75), y = 3.3), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.45) + TXT(aes(x = mm(77.5), y = 3.3, label = "Filled: primary or co-primary row"), hjust = 0) +
  geom_point(aes(x = mm(3), y = 2.55), shape = 21, fill = "white", colour = INK, size = 1.9, stroke = 0.6) + TXT(aes(x = mm(5.5), y = 2.55, label = "Open: beside or sensitivity row"), hjust = 0) +
  geom_point(aes(x = mm(75), y = 2.55), shape = 13, colour = INK, size = 2.3, stroke = 0.55) + TXT(aes(x = mm(77.5), y = 2.55, label = "Crossed: not colocalized with ARG1 (in the row's color)"), hjust = 0) +
  geom_rect(aes(xmin = mm(3), xmax = mm(11), ymin = 0.9, ymax = 1.95), fill = COL_SIG_L[["A"]], alpha = 0.45) +
  TXT(aes(x = mm(13), y = 1.7, label = "A-based reference for B: signal A's estimate and 95% CI, which assumes the same ratio at both signals."), hjust = 0) +
  TXT(aes(x = mm(13), y = 1.15, label = "Intervals are pointwise 95% CIs (Wald), not adjusted for multiplicity."), hjust = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, 3.8), expand = FALSE) + theme_void()

# ---- b: posterior of colocalization with ARG1 ----------------------------------------------------------------------------------------------------------
P12 <- c(def = num("L240", 3), cons = num("L241"))
stopifnot(P12[["def"]] == 1e-5, P12[["cons"]] == 1e-6)
ANA <- c(W = "marginal", A = "signal A (conditioned on B)", B = "signal B (conditioned on A)")
pick <- function(pr, an, p) { r <- co[pair == pr & analysis == an]; r <- r[abs(r$p12 - p) < p * 1e-6]; stopifnot(nrow(r) == 1L); r$pp_h4 }
stopifnot(all(abs(round(c(pick("P-D", ANA[["A"]], P12[["def"]]), pick("P-D", ANA[["B"]], P12[["def"]])), 3) - num("L057")) < 1e-9),
          all(abs(c(pick("P-D", ANA[["A"]], P12[["cons"]]), pick("P-D", ANA[["B"]], P12[["cons"]])) - num("L265")) < 1e-12))
en <- ecen[gene == "ARG1" & outcome == "T2D"]
stopifnot(nrow(en) == 1L, en$call == "NOT_SUPPORTED_H3", abs(en$PP.H3 - num("L266", 1)) < 1e-12, abs(en$PP.H4 - num("L266", 2)) < 1e-12)
ek_def <- num("L266", 5)                                                                   # Biobank Japan (CKB-free), default prior: the ledger value
ek <- eckf[gene == "ARG1" & outcome == "T2D"]; stopifnot(nrow(ek) == 1L, abs(ek$PP.H4 - ek_def) < 1e-12)
ek_con <- ek$`PP.H4_p12_1e-6`                                                              # the conservative-prior point (approved file, S2-D25)
# WHOLE-LOCUS values (S2-D25 visual check, "two surfaces carried different values for one quantity"): the whole-locus rows use the SAME values as the Results,
# read from the ledger: ARG1-T2D L051 (0.925 default / 0.553 conservative) and ARG1-arginine L052 (0.996 default). The earlier signal-specific grid
# (S1R1_arg1_coloc.tsv) holds a different whole-locus ARG1-T2D analysis (0.932 / 0.579); it is NOT plotted for the whole-locus T2D row, and the
# difference is recorded in the source-data frame `whole_locus_reconciliation`. The grid supplies every signal-specific value and the arginine
# whole-locus conservative-prior point, for which no ledger row exists (its default-prior value equals L052 at the printed precision).
W_T2D <- c(def = num("L051", 1), cons = num("L051", 2)); W_ARG_DEF <- num("L052")
stopifnot(abs(pick("P-M", ANA[["W"]], P12[["def"]]) - W_ARG_DEF) < 5e-4)
G_T2D <- c(def = pick("P-D", ANA[["W"]], P12[["def"]]), cons = pick("P-D", ANA[["W"]], P12[["cons"]]))
stopifnot(all(abs(G_T2D - W_T2D) < c(0.01, 0.03)), all(round(W_T2D, 3) == c(0.925, 0.553)))
RECON <- data.table(quantity = "ARG1-T2D whole-locus PP.H4", prior = c("default p12 = 1e-5", "conservative p12 = 1e-6"), plotted_value_ledger_L051 = unname(W_T2D),
                    earlier_grid_value_not_plotted = unname(G_T2D), difference = unname(G_T2D - W_T2D),
                    note = "different analyses of the same pair; the Results print the L051 values; this panel plots them")
ROWS <- c(whole = "Whole locus", A = "Signal A | B", B = "Signal B | A", eas = "AGEN (with CKB)", ekf = "Biobank Japan (CKB-free)")
mk <- function(tr, rw, pr, pp, src = "S1R1_arg1_coloc.tsv (signal-specific grid)") data.table(trait = tr, row = rw, prior = pr, pp = pp, src = src)
cd <- rbind(
  rbind(mk("Arginine", "whole", "def", W_ARG_DEF, "ledger L052 (as the Results)"), mk("Arginine", "whole", "cons", pick("P-M", ANA[["W"]], P12[["cons"]]), "S1R1_arg1_coloc.tsv (no ledger row for the conservative prior)")),
  mk("Arginine", "A", c("def", "cons"), c(pick("P-M", ANA[["A"]], P12[["def"]]), pick("P-M", ANA[["A"]], P12[["cons"]]))),
  mk("Arginine", "B", c("def", "cons"), c(pick("P-M", ANA[["B"]], P12[["def"]]), pick("P-M", ANA[["B"]], P12[["cons"]]))),
  mk("T2D", "whole", c("def", "cons"), unname(W_T2D), "ledger L051 (as the Results)"),
  mk("T2D", "A", c("def", "cons"), c(pick("P-D", ANA[["A"]], P12[["def"]]), pick("P-D", ANA[["A"]], P12[["cons"]]))),
  mk("T2D", "B", c("def", "cons"), c(pick("P-D", ANA[["B"]], P12[["def"]]), pick("P-D", ANA[["B"]], P12[["cons"]]))),
  mk("T2D", "eas", c("def", "cons"), c(en$PP.H4, en$`PP.H4_p12_1e-6`), "results_census.tsv"),
  mk("T2D", "ekf", c("def", "cons"), c(ek_def, ek_con), "results_ckbfree.tsv"))
cd <- cd[!is.na(pp)]
stopifnot(all(cd$pp >= 0 & cd$pp <= 1))
YD <- c(argH = 10.6, arg_whole = 9.6, arg_A = 8.6, arg_B = 7.6, t2dH = 6.3, t2d_whole = 5.3, t2d_A = 4.3, t2d_B = 3.3, easH = 2.0, t2d_eas = 1.0, t2d_ekf = 0.0)
cd[, ykey := paste0(fifelse(trait == "Arginine", "arg_", "t2d_"), row)]
stopifnot(all(cd$ykey %in% names(YD)))
cd[, prior := factor(prior, levels = c("def", "cons"), labels = c("Default prior", "Conservative prior"))]
DODGE <- c("Default prior" = -0.2, "Conservative prior" = 0.3)      # centre stays at +0.05 (row words); wide enough that the crossed circle and the triangle never overlap
cd[, yd := YD[ykey] + DODGE[as.character(prior)]]
cd[, rowlab := ROWS[row]]
COL_ROW <- c("Whole locus" = "black", "Signal A | B" = COL_SIG[["A"]], "Signal B | A" = COL_SIG[["B"]], "AGEN (with CKB)" = "black", "Biobank Japan (CKB-free)" = "black")
# not colocalized rows (default-prior PP.H4 below the rule): T2D signal B | A and both East Asian rows: crossed circle at the default prior + the word and PP.H4
NC <- c("t2d_B", "t2d_eas", "t2d_ekf")
stopifnot(all(cd[ykey %in% NC & prior == "Default prior"]$pp < num("L236")), all(cd[!ykey %in% NC & prior == "Default prior"]$pp >= num("L236")))
dl <- c(argH = "**Arginine (ARG1 protein with arginine)**", arg_whole = ROWS[["whole"]], arg_A = ROWS[["A"]], arg_B = ROWS[["B"]], t2dH = "**T2D (ARG1 protein with T2D)**",
        t2d_whole = ROWS[["whole"]], t2d_A = ROWS[["A"]], t2d_B = ROWS[["B"]], easH = "**T2D, East Asian ancestry**", t2d_eas = ROWS[["eas"]], t2d_ekf = ROWS[["ekf"]])
DB <- seq(0, 1, 0.2)
XLIM_B <- c(-0.05, 1.95)
DAX <- ticks_set("S5_b_x", DB, c(-0.04, 1.2), "S5 b x, PP.H4")
# row words (right of the axis): the status word travels with the crossed mark; PP.H4 and the most probable hypothesis are the ledger values
w_B <- paste0("not colocalized (PP.H4 ", led("L057", 2, fmtd(3), where = "S5b: signal B | A T2D PP.H4, default prior", anchor = "^not colocalized \\(PP\\.H4 {}\\)$"), ")")
w_E <- paste0("not colocalized (PP.H4 ", led("L266", 2, fmtd(3), where = "S5b: AGEN PP.H4 (default prior)", anchor = "^not colocalized \\(PP\\.H4 {}; PP\\.H3 0\\.99\\)$"), "; PP.H3 ",
              led("L266", 1, fmtd(2), where = "S5b: AGEN PP.H3 (default prior), the most probable hypothesis", anchor = "PP\\.H3 {}\\)$"), ")")
w_K <- paste0("not colocalized (PP.H4 ", led("L266", 5, fmtd(3), where = "S5b: Biobank Japan PP.H4 (default prior)", anchor = "^not colocalized \\(PP\\.H4 {}; PP\\.H1 0\\.74\\)$"), "; PP.H1 ",
              led("L266", 3, fmtd(2), where = "S5b: Biobank Japan PP.H1 (default prior), the most probable hypothesis", anchor = "PP\\.H1 {}\\)$"), ")")
dw <- data.table(y = c(YD[["t2d_B"]] + 0.05, YD[["t2d_eas"]] + 0.05, YD[["t2d_ekf"]] + 0.05), lab = c(w_B, w_E, w_K))
a_cons <- cd[ykey == "t2d_A" & prior == "Conservative prior"]
LAB674 <- led("L265", 1, fmtd(3), where = "S5b: T2D signal A | B PP.H4, conservative prior", anchor = "^{}$")
pb <- ggplot(cd, aes(pp, yd)) +
  annotate("segment", x = num("L236"), xend = num("L236"), y = YD[["t2d_ekf"]] - 0.5, yend = YD[["arg_whole"]] + 0.5, linewidth = 0.3, colour = LINE_GREY, linetype = "dashed") +
  annotate("text", x = num("L236") - 0.035, y = (YD[["arg_B"]] + YD[["t2dH"]]) / 2, label = "PP.H4 threshold", hjust = 1, size = ANNOT, family = BASE_FONT, colour = "grey30") +
  geom_point(data = cd[!(ykey %in% NC & prior == "Default prior")], aes(colour = rowlab, shape = prior), size = 1.5, stroke = 0) +
  geom_point(data = cd[ykey %in% NC & prior == "Default prior"], aes(colour = rowlab), shape = 13, size = 2.0, stroke = 0.4) +   # RULE: a crossed marker takes the colour of its row (signal colour for Signal B | A, black for the black rows), as in panel a
  annotate("text", x = a_cons$pp - 0.03, y = a_cons$yd, label = LAB674, hjust = 1, size = ANNOT, family = BASE_FONT, colour = COL_SIG[["A"]]) +
  geom_text(data = dw, aes(x = 1.05, y = y, label = lab), hjust = 0, size = ANNOT, family = BASE_FONT, inherit.aes = FALSE) +
  scale_shape_manual(values = c("Default prior" = 16, "Conservative prior" = 17), name = NULL) +
  scale_colour_manual(values = COL_ROW, guide = "none") +
  scale_x_continuous(breaks = DB, labels = DAX, expand = expansion(0), guide = guide_axis(cap = "upper")) +
  scale_y_continuous(breaks = unname(YD), labels = unname(dl), expand = expansion(0)) +
  coord_cartesian(xlim = XLIM_B, ylim = c(YD[["t2d_ekf"]] - 0.75, YD[["argH"]] + 0.5)) +
  labs(x = "PP.H4 with ARG1", y = NULL) +
  theme_p31r() +
  theme(axis.text.y = ggtext::element_markdown(size = ANNOT_PT), axis.ticks.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.02, 0.99), legend.justification = c(0, 1), legend.text = element_text(size = ANNOT_PT),
        legend.key.size = unit(8, "pt"), legend.spacing.y = unit(0, "pt"), legend.background = element_rect(fill = "white", colour = NA)) +
  guides(shape = guide_legend(override.aes = list(colour = "black", size = 1.5)))

topa <- (pa1 + labs(tag = "a")) + pa2 + plot_layout(widths = c(1, 1))
bota <- pa3 + pa4 + plot_layout(widths = c(1, 1))
page1 <- wrap_elements(full = p_hdr) / topa / bota / wrap_elements(full = ka) / (pb + labs(tag = "b")) +
  plot_layout(heights = c(5, 34, 74, 19, 62))
FR_A <- f3[, .(estimand_id, estimand, resource, setting, panel, signal, signal_variant, value, se, ci_lo, ci_hi, row, stat, prim, sens, fires,
               coloc_clause_pph4_conditional_vs_ARG1_protein)]
af3_save(page1, h_mm = 210, panels = c("a", "b"),
         panel_map = c(a = "per-signal estimates at signals A and B for four estimands, with the A-based reference for B",
                       b = "posterior of colocalization with ARG1: whole locus and signal-specific, both priors, and East Asian rows"),
         frames = list(per_signal_estimates = FR_A, posterior_colocalization = cd[, .(trait, row, prior, pp, src)], whole_locus_reconciliation = RECON))

# =========================================================================================================== PAGE 2
unit_reset()
UNIT <- "FigureS5_p2"
af3_begin(UNIT)
cs <- af3_read("R1_stage1/S1R1_R2_12_chen_instrument/credible_set_status.tsv", "c", sep = "\t")
stopifnot(!anyNA(cs$set_signal), setequal(unique(cs$set_signal), c("A", "B")), setequal(unique(cs$variant_role), c("A", "B", "C")))
# L986 holds the set size per method and signal; the file's set_size agrees cell by cell
L986 <- p31r_ledger("L986")
cs[, lkey := ifelse(method == "1_T2D", "m1_T2D_503", ifelse(method == "1_protein", "m1_protein_503", ifelse(method == "2", "m2_insample", paste0("m", method, "_", panel))))]
cs1 <- unique(cs[, .(set_signal, method, panel, lkey, set_size, adm = admissibility_rule3i, rob = robustness_rule3ii)])
stopifnot(nrow(cs1) == 18L, !anyDuplicated(cs1[, .(set_signal, lkey)]))
for (i in seq_len(nrow(cs1))) stopifnot(cs1$set_size[i] == L986[[cs1$set_signal[i]]][[cs1$lkey[i]]])
# members among the three named variants: rs2246012 (role C), the signal-A lead, the signal-B lead
mem <- cs[status == "included", .(inc = paste(sort(unique(sub("^rs", "rs", variant))), collapse = ", ")), by = .(set_signal, lkey)]
cs1 <- merge(cs1, mem, by = c("set_signal", "lkey"), all.x = TRUE)
cs1[is.na(inc), inc := ""]
cs1[, state := fcase(grepl("^inadmissible", adm), "inadmissible", grepl("^admissible", adm), "admissible",
                     grepl("^NOT robust", rob), "not robust", grepl("^robust", rob), "robust", default = NA_character_)]
stopifnot(!anyNA(cs1$state))
# the mark of a cell: admissible or robust = filled dot; inadmissible SuSiE fit = HATCHED circle (the crossed circle means "no colocalization" on page 1
# and in the main set); not robust to the reference LD = dotted ring. Each mark is built in millimetres (P31R_af3_glyphs.R).
ROWDEF <- data.table(lkey = c("m1_T2D_503", "m1_protein_503", "m2_insample", "m3_503", "m3_489", "m4_503", "m4_489", "m5_503", "m5_489"),
                     lab = c("SuSiE, T2D (FinnGen),\nreference panel (503)", "SuSiE, ARG1 protein (UKB-PPP),\nreference panel (503)", "SuSiE, ARG1 protein,\nin-sample LD (UKB-PPP)",
                             "Shared ABF set (FinnGen),\nconditional, reference panel (503)", "Shared ABF set (FinnGen),\nconditional, reference panel (489)",
                             "ABF set, T2D (FinnGen),\nconditional, reference panel (503)", "ABF set, T2D (FinnGen),\nconditional, reference panel (489)",
                             "Shared ABF set (DIAMANTE),\nconditional, reference panel (503)", "Shared ABF set (DIAMANTE),\nconditional, reference panel (489)"))
# row 1 is taller: its signal-B cell carries three lines (the set is empty; the wording must not read as a null)
ROWDEF[, h := c(1.45, rep(1, 8))]
ROWDEF[, top := -(cumsum(h) - h)]
cs1 <- merge(cs1, ROWDEF, by = "lkey")
cs1[, x := fifelse(set_signal == "A", 0.455, 0.745)]
# printed set sizes from the ledger (L986), one per cell
cs1[, sz_lab := mapply(function(s, k) led("L986", c(s, k), fmt0, where = paste0("S5c: credible-set size, signal ", s, ", ", k), anchor = "; {} variants?$"),
                       set_signal, lkey)]
cs1[, cell := paste0(sz_lab, fifelse(set_size == 1, " variant", " variants"))]
cs1[, inc2 := fifelse(set_size == 0, "no credible set returned;\nnot evidence of no association", fifelse(nzchar(inc), paste0("includes ", inc), "none of the three named variants"))]
stopifnot(sum(cs1$set_size == 0) == 1L, cs1[set_size == 0, set_signal == "B" & lkey == "m1_T2D_503" & state == "inadmissible"])
# reference panel numerals in the row labels: 503 drawn five times, 489 three times (ledger L936, L937)
led("L936", NULL, fmt0, where = "S5c row labels: 503-person reference panel", anchor = "reference panel \\({}\\)$", count = 5L)
led("L937", NULL, fmt0, where = "S5c row labels: 489-person reference panel", anchor = "reference panel \\({}\\)$", count = 3L)
YM <- 0.32                                                                          # first text line, below the row top (units)
pcg <- ggplot(cs1) +
  geom_hline(yintercept = ROWDEF$top - ROWDEF$h, colour = "grey85", linewidth = 0.3) +
  geom_point(data = cs1[state %in% c("admissible", "robust")], aes(x = x - 0.14, y = top - YM), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.4) +
  geom_hatchdot(data = cs1[state == "inadmissible"], aes(x = x - 0.14, y = top - YM), ink = INK) +
  geom_dotring(data = cs1[state == "not robust"], aes(x = x - 0.14, y = top - YM), r_mm = 1.0, ink = INK) +
  geom_text(aes(x = x - 0.125, y = top - YM, label = paste0(state, "; ", cell)), hjust = 0, size = ANNOT, family = BASE_FONT) +
  geom_text(aes(x = x - 0.125, y = top - 0.52, label = inc2), hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, colour = "grey30", lineheight = 1.0) +
  geom_text(data = ROWDEF, aes(x = 0.0, y = top - 0.5, label = lab), hjust = 0, size = ANNOT, family = BASE_FONT, lineheight = 1.0) +
  annotate("text", x = 0.455 - 0.14, y = 0.2, label = "Signal A set", hjust = 0, size = ANNOT, family = BASE_FONT, fontface = "bold") +
  annotate("text", x = 0.745 - 0.14, y = 0.2, label = "Signal B set", hjust = 0, size = ANNOT, family = BASE_FONT, fontface = "bold") +
  annotate("text", x = 0.0, y = 0.2, label = "Method, trait and reference LD", hjust = 0, size = ANNOT, family = BASE_FONT, fontface = "bold") +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(-sum(ROWDEF$h) - 0.02, 0.7), expand = FALSE) + theme_void()
kc2 <- ggplot() +
  geom_point(aes(x = mm(3), y = 1), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.4) + TXT(aes(x = mm(6), y = 1, label = "Admissible or robust"), hjust = 0) +
  geom_hatchdot(aes(x = mm(42), y = 1), data = data.frame(x = 1, y = 1), ink = INK) + TXT(aes(x = mm(45), y = 1, label = "Inadmissible SuSiE fit, hatched"), hjust = 0) +
  geom_dotring(aes(x = mm(86), y = 1), data = data.frame(x = 1, y = 1), r_mm = 1.0, ink = INK) + TXT(aes(x = mm(89), y = 1, label = "Not robust to the reference LD, dotted ring"), hjust = 0) +
  TXT(aes(x = mm(3), y = 0.35, label = "Shared ABF set: the set shared by the ARG1 protein and T2D associations. Credible-set membership and pairwise LD answer"), hjust = 0) +
  TXT(aes(x = mm(3), y = -0.15, label = "different questions, so no set is shown as an LD tag. Status is read from each method's own row of the credible-set table."), hjust = 0) +
  TXT(aes(x = mm(3), y = -0.65, label = "Figure 4b uses the same key: hatched (inadmissible fit) and a dotted ring (not robust to the reference LD)."), hjust = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.95, 1.6), expand = FALSE) + theme_void()

# ---- d: the 1-Mb overview of the ARG1 region (S2-D25: its sources now have approved AF3_S5 rows) --------------------------------------------------
# Association strength (-log10 p) for the ARG1 protein (UKB-PPP), arginine (CLSA) and T2D (FinnGen R10) over the colocalization window, gene models, and
# signals A and B (dashed lines; they lie a few kb apart, so at this scale they are one position). Points are NOT coloured by LD here: the 503-person
# dosage panel has no approved row for this figure. The interval enlarged in Figure 4b is shaded in the gene track.
ARG_F <- "V4_local/staging/l3-finemap/ARG1_6q23/harmonised.tsv.gz"
T2D_F <- paste0(.P31R_EXT, "/D/data/paper31_v4_raw/ukbfree_outcomes/finngen_R10_T2D.gz")
ENS_D <- "R2_rewrite/figures/data/Figure_4"
hh <- af3_read(ARG_F, "d")
invisible(af3_source(paste0(T2D_F, ".tbi"), "d"))
t2d_path <- af3_source(T2D_F, "d")
rec <- jsonlite::fromJSON(af3_source(file.path(ENS_D, "FETCH_RECORD.json"), "d"))
ENS_F <- af3_source(file.path(ENS_D, rec$file), "d")
stopifnot(unname(tools::md5sum(ENS_F)) == rec$md5, rec$assembly == "GRCh38", rec$release == 116)
WINA <- as.numeric(strsplit(sub("^6:", "", rec$window_grch38), "-")[[1]])
stopifnot(length(WINA) == 2L, all(hh$pos38 >= WINA[1] & hh$pos38 <= WINA[2]), !anyDuplicated(hh$snp))
for (.p in c("Rsamtools", "GenomicRanges", "IRanges")) if (!requireNamespace(.p, quietly = TRUE)) stop("package missing: ", .p)
tf <- Rsamtools::TabixFile(t2d_path)
chr <- if ("6" %in% Rsamtools::seqnamesTabix(tf)) "6" else "chr6"
raw <- Rsamtools::scanTabix(tf, param = GenomicRanges::GRanges(chr, IRanges::IRanges(WINA[1], WINA[2])))[[1]]
dd <- fread(text = paste(raw, collapse = "\n"), header = FALSE, sep = "\t")
setnames(dd, c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"))
mm_ <- merge(hh, dd[, .(pos38 = pos, ref, alt, bD = beta, sD = sebeta, af_alt, mlogpD = mlogp)], by = "pos38")
mm_ <- mm_[(alt == a1 & ref == a2) | (alt == a2 & ref == a1)]
mm_[, flip := alt != a1][, afD := ifelse(flip, 1 - af_alt, af_alt)]
mm_ <- mm_[palindromic != 1 | (pmin(af_ref, 1 - af_ref) <= 0.40 & pmin(afD, 1 - afD) <= 0.40 & abs(afD - af_ref) <= 0.10)]   # same palindrome rule as the signal analysis
mm_ <- mm_[!duplicated(snp)][order(pos38)]
mlp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
mm_[, `:=`(lpP = mlp(z__P__ARG1), lpM = mlp(z__M__GCST90200372), lpD = mlogpD)]
stopifnot(max(abs(mm_$lpD - -log10(2 * pnorm(-abs(mm_$bD / mm_$sD)))), na.rm = TRUE) < 0.05)
posA <- hh[rsid_ref == SIGA, pos38]; posB <- hh[rsid_ref == SIGB, pos38]
stopifnot(length(posA) == 1L, length(posB) == 1L, abs(posA - posB) < 10000)
ZOOM <- c(131.48, 131.67)       # layout constant: the interval Figure 4b enlarges (declared in P31R_Fig4.R; every credible-set member lies inside it)
# gene models (Ensembl 116, GRCh38): protein-coding extent (grey), canonical transcript (black), exons as ticks
ens <- as.data.table(jsonlite::fromJSON(ENS_F, simplifyDataFrame = TRUE)[, c("feature_type", "id", "Parent", "external_name", "biotype", "start", "end", "strand", "is_canonical", "assembly_name")])
stopifnot(all(ens$assembly_name == "GRCh38"))
gnA <- ens[feature_type == "gene" & biotype == "protein_coding", .(gene = id, symbol = external_name, strand)]
txA <- ens[feature_type == "transcript" & is_canonical %in% 1 & Parent %in% gnA$gene, .(tx = id, gene = Parent, t0 = start, t1 = end)]
stopifnot(nrow(txA) == nrow(gnA), !anyDuplicated(txA$gene))
pcxA <- ens[feature_type == "transcript" & biotype == "protein_coding" & Parent %in% gnA$gene, .(p0 = min(start), p1 = max(end)), by = .(gene = Parent)]
gnA <- merge(merge(gnA, txA, by = "gene"), pcxA, by = "gene")
exA <- merge(ens[feature_type == "exon" & Parent %in% gnA$tx, .(tx = Parent, e0 = start, e1 = end)], gnA[, .(tx, strand)], by = "tx")
stopifnot(all(c("ARG1", "MED23", "ENPP3") %in% gnA$symbol), all(c(posA, posB) >= gnA[symbol == "ARG1", p0] & c(posA, posB) <= gnA[symbol == "ARG1", p1]))
gnA <- gnA[p1 >= WINA[1] & p0 <= WINA[2]]
setorder(gnA, p0)
gnA[, `:=`(y = fifelse(strand > 0, 1, -1), lx = (pmax(p0, WINA[1]) + pmin(p1, WINA[2])) / 2 / 1e6)]
# label tiers for the minus-strand genes (OR2A4 and CTAGE9 lie 9 kb apart, MED23 and OR2A4 close): fixed by gene, asserted below
TIER <- c(EPB41L2 = 1L, MED23 = 1L, OR2A4 = 2L, CTAGE9 = 3L, CCN2 = 1L)   # S2-D26: labels sit two units off their line (plus strand above, minus strand below), never on it
gnA[, tier := fifelse(strand > 0, 1L, unname(TIER[symbol]))]
stopifnot(!anyNA(gnA$tier))
gnA[, ly := fifelse(strand > 0, 3.4, -3.4 - (tier - 1L) * 2.1)]
gnA[, hj := fifelse(p0 < WINA[1], 0, 0.5)]; gnA[hj == 0, lx := WINA[1] / 1e6 + 0.004]
gnA[, lab := fifelse(strand > 0, paste(symbol, "\u2192"), paste("\u2190", symbol))]
MBA <- 1e6; XA <- seq(131, 132, by = 0.25)
TRKD <- rbind(mm_[, .(trk = "P", pos = pos38 / MBA, y = lpP)], mm_[, .(trk = "M", pos = pos38 / MBA, y = lpM)], mm_[, .(trk = "T", pos = pos38 / MBA, y = lpD)])
TRK_LAB_D <- c(P = "ARG1 protein\n(UKB-PPP)", M = "Arginine\n(CLSA)", T = "T2D\n(FinnGen)")
LOG_TTL <- "\u2212log<sub><span style='font-size:7pt'>10</span></sub> <i>p</i>"
REG$literals[[length(REG$literals) + 1]] <- list(text = "10", why = "scale: base of the logarithm in the y title (-log10 p)", count = 1L, anchor = "log {}")
vlA <- geom_vline(xintercept = posA / MBA, colour = COL_SIG[["A"]], linetype = "dashed", linewidth = 0.35)
vlB <- geom_vline(xintercept = posB / MBA, colour = COL_SIG[["B"]], linetype = "dashed", linewidth = 0.35)
xlim_d <- WINA / MBA
trk_d <- function(tk, ytitle = FALSE) {
  d <- TRKD[trk == tk]; d[, strip_lab := factor(TRK_LAB_D[[tk]])]
  lim <- range(c(0, d$y)); EXP <- c(0.06, 0.1); edge <- lim + c(-EXP[1], EXP[2]) * diff(lim)
  brk <- scales::breaks_extended(n = 3)(lim); brk <- brk[brk >= edge[1] & brk <= edge[2] - 0.12 * diff(edge)]
  ggplot(d, aes(pos, y)) + vlA + vlB +
    geom_point(colour = PAL_NEUTRAL[["other"]], size = 0.55, stroke = 0) +
    facet_grid(strip_lab ~ ., switch = "y") +
    scale_x_continuous(limits = xlim_d, breaks = XA, expand = c(0, 0)) +
    scale_y_continuous(limits = lim, breaks = brk, labels = ticks_set(paste0("S5_d_y_", tk), brk, edge, paste0("S5 d y, track ", tk)), expand = expansion(mult = EXP)) +
    labs(x = NULL, y = if (ytitle) LOG_TTL else NULL) + theme_p31r() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_markdown(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95), plot.margin = margin(2, 12, 2, 4))
}
dP <- trk_d("P") + labs(tag = "d") + theme(plot.tag = element_text(size = TAG_PT, face = "bold", family = BASE_FONT)); dM <- trk_d("M", ytitle = TRUE); dT <- trk_d("T")
XA_LAB <- ticks_set("S5_d_x_chr6", XA, xlim_d, "S5 d x, chromosome 6 position (Mb, GRCh38)", labels = scales::label_number(accuracy = 0.01, style_negative = "minus")(XA))
segA <- gnA[, .(x = pmax(p0, WINA[1]) / MBA, xend = pmin(p1, WINA[2]) / MBA, y)]
canA <- gnA[, .(x = pmax(t0, WINA[1]) / MBA, xend = pmin(t1, WINA[2]) / MBA, y)]
exT <- merge(exA, gnA[, .(tx, y)], by = "tx")[, .(xc = (e0 + e1) / 2 / MBA, y)][xc >= xlim_d[1] & xc <= xlim_d[2]]
pdg2 <- ggplot() +
  annotate("rect", xmin = ZOOM[1], xmax = ZOOM[2], ymin = -8.4, ymax = 4.8, fill = "grey93", colour = NA) +
  annotate("text", x = ZOOM[2] + 0.004, y = 6.4, label = "interval enlarged in Figure 4b", hjust = 0, vjust = 0.5, size = ANNOT, family = BASE_FONT, colour = "grey30") +
  annotate("segment", x = posA / MBA, xend = posA / MBA, y = -2.1, yend = 2.1, colour = COL_SIG[["A"]], linetype = "dashed", linewidth = 0.35) +   # signal lines are stubs through the gene models only: they never cross a gene label or arrow
  annotate("segment", x = posB / MBA, xend = posB / MBA, y = -2.1, yend = 2.1, colour = COL_SIG[["B"]], linetype = "dashed", linewidth = 0.35) +
  geom_segment(data = segA, aes(x = x, xend = xend, y = y, yend = y), colour = "grey55", linewidth = 0.35) +
  geom_segment(data = canA, aes(x = x, xend = xend, y = y, yend = y), colour = INK, linewidth = 0.35) +
  geom_segment(data = exT, aes(x = xc, xend = xc, y = y - 0.35, yend = y + 0.35), colour = INK, linewidth = 0.3) +
  geom_text(data = gnA, aes(x = lx, y = ly, label = lab, hjust = hj), size = ANNOT, family = BASE_FONT, fontface = "italic") +
  annotate("text", x = posA / MBA - 0.006, y = 6.4, label = paste0("A ", SIGA, "; B ", SIGB), hjust = 1, vjust = 0.5, size = ANNOT, family = BASE_FONT, colour = INK) +
  scale_x_continuous(limits = xlim_d, breaks = XA, labels = XA_LAB, expand = c(0, 0)) + scale_y_continuous(limits = c(-9.4, 7.4), expand = c(0, 0)) +
  labs(x = "Chromosome 6 position (Mb, GRCh38)", y = NULL) + theme_p31r() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = margin(2, 12, 4, 4))
pd <- (dP / dM / dT / pdg2) + plot_layout(heights = c(1, 1, 1, 1.65))
page2 <- (wrap_elements(full = pcg + labs(tag = "c") + theme(plot.tag = element_text(size = TAG_PT, face = "bold", family = BASE_FONT))) / wrap_elements(full = kc2) /
          wrap_elements(full = pd)) +
  plot_layout(heights = c(84, 13, 96))
af3_save(page2, h_mm = 196, panels = c("c", "d"),
         panel_map = c(c = "credible sets at signals A and B for every method, status read per cell", d = "1-Mb overview of the ARG1 region: association tracks and gene models"),
         frames = list(credible_set_grid = cs1[, .(set_signal, method, panel, lkey, set_size, state, inc, inc2)], overview_tracks = TRKD, overview_genes = gnA[, .(symbol, strand, p0, p1, t0, t1)]))
