# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S6: East Asian sensitivities (the outcome-side scatter moved from the old Figure 5b, decision
# F17; the strict per-variant-sample-size subset and the CKB-free composite). Illustrates SM section S-R13. Two page units (continuation page, F9):
#   FigureS6     a  outcome-side (European instrument, East Asian outcome data): European against East Asian estimate for the BINARY-outcome evaluable pairs
#                   (38 of 52); filled = compatible, open = not significant, NO tier hue (F5); RING = palindromic instrument (A/T or C/G), derived from the
#                   instrument allele table (approved for this figure, S2-D25); 95% intervals of both estimates in light lines. Only the compatible pairs are
#                   labelled (white label boxes with leader lines, repelled), so no label stands for a point it does not name.
#                b  strict subset (outcome files that report per-variant sample sizes: the T2D and eGFR files): its T2D part, with the strict state
#   FigureS6_p2  c  the CKB-FREE COMPOSITE (S2-D25: census file approved for this panel): the 48-row ancestry-matched census with the 13 CKB-free re-run rows
#                   swapped in for the same gene x outcome pairs (T2D from AGEN and stroke from GIGASTROKE include CKB; every other outcome resource does
#                   not). All 48 pairs are drawn or named: PP.H4 at the default (circle) and conservative (triangle) prior; robust = filled dark,
#                   colocalized but not robust = light fill with a dark outline, not colocalized = crossed circle (the status glyph of the main figures; the
#                   hatched circle is reserved for an inadmissible fit), not run = a text cell with no marker. The composite reproduces L084 (8 and 4).
# HOLD H9 (applies here exactly as in the main figures): the East Asian eGFR file EGFR_EAS_2P2M has no magnitude on any published surface. The 14
# eGFR pairs of that file are counted (ledger L987, L967) and NOT plotted. The composite's eGFR rows come from a DIFFERENT file (TWB2+BBJ meta-analysis,
# accession TWB2_BBJ_eGFR_PMID39472450), which H9 does not cover; this is asserted below, and the plotted values there are posterior probabilities.
# Sources: ONLY through the approved manifest (AF3_S6).
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS6.R"

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
suppressPackageStartupMessages({ library(ggrepel) })
EN <- "\u2013"
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
TRUEish <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1")
INK <- "grey15"; LINE_GREY <- "grey45"; N_ <- PAL_NEUTRAL
TXT <- function(...) geom_text(..., size = ANNOT, family = BASE_FONT)
fmt0 <- function(v) formatC(v, format = "f", digits = 0)
# ticks on a narrow axis can be joined into one printed line ("-3 -1 -0.3 0 ..."): the anchor allows other numbers on the line (as Figure 2b's generator did)
ticks_set <- function(axis, breaks, lims, why, transform = "asinh") {
  stopifnot(all(breaks > lims[1] & breaks < lims[2]), !anyDuplicated(breaks))
  labs <- p31r_lab(breaks)
  for (t in labs) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis, transform = transform, count = 1L,
    anchor = "^(?:[\u2212+-]?[0-9.,]+\\s)*{}(?:\\s[\u2212+-]?[0-9.,]+)*$")
  labs
}
unit_reset <- function() { rm(list = ls(.p31r_read_set), envir = .p31r_read_set); invisible(NULL) }   # each page unit pins only the files IT read
mm <- function(x) x / 165

# =========================================================================================================== PAGE 1
UNIT <- "FigureS6"
af3_begin(UNIT)
tier <- af3_read("v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv", "a")
# level 2 (RC1b): the released tiered table masks the hold-H9 cells with a text marker; those rows are dropped below (H9), so read the three columns as numbers (marker -> NA)
for (.cc in c("eas_beta_mr", "eas_se_mr", "eas_q")) if (is.character(tier[[.cc]])) set(tier, j = .cc, value = suppressWarnings(as.numeric(tier[[.cc]])))
mr   <- af3_read("V4/source_data/tables_from_rds/phase2_mr_screen__mr_all_results.csv", "a")
eas  <- af3_read("V4/source_data/supp_tables/SuppTable_EAS_axis.csv", "b")
ins  <- af3_read("V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv", "a")   # palindromic flag (approved for panel a, S2-D25)

# ---- the outcome-side evaluable pairs: binary outcomes only (H9) ----------------------------------------------------------------------------
stopifnot(nrow(tier) == num("L010"), !anyDuplicated(tier[, .(gene, outcome)]), nrow(eas) == num("L010"), !anyDuplicated(eas[, .(gene, outcome)]))
ev <- tier[eas_outcome_state != "EAS_NOT_EVALUABLE"]
stopifnot(nrow(ev) == num("L080", 1), sum(ev$eas_outcome_state == "EAS_COMPATIBLE") == num("L080", 2))
HELD <- ev$eas_outcome_source == "EGFR_EAS_2P2M"
stopifnot(sum(HELD) == num("L987", "egfr_evaluable"), sum(HELD) == num("L967", "evaluable"), all(ev[HELD]$outcome == "eGFR"), !any(ev[!HELD]$outcome == "eGFR"))
bin <- ev[!HELD]                                                                       # H9: eGFR rows never enter a plotted frame
stopifnot(nrow(bin) == num("L987", "binary_evaluable"), !any(bin$eas_outcome_source == "EGFR_EAS_2P2M"),
          sum(bin$eas_outcome_state == "EAS_COMPATIBLE") == num("L080", 2) - num("L967", "compatible"))
# European side: the full-precision Wald ratio / IVW of the MR screen (the tiered table's mr_beta is rounded to 6 d.p.)
mr1 <- mr[method %in% c("Wald ratio", "Inverse variance weighted"), .(gene, outcome = outcome_name, b_eur = b, se_eur = se)]
stopifnot(!anyDuplicated(mr1[, .(gene, outcome)]))
bin <- merge(bin, mr1, by = c("gene", "outcome"))
stopifnot(nrow(bin) == num("L987", "binary_evaluable"), max(abs(bin$b_eur - bin$mr_beta)) < 1e-5)
# same-sign count of the binary pairs: 45 of 52 overall, none of the 7 opposite-sign pairs is an eGFR pair (L967: all 14 eGFR pairs same sign)
bin[, same := sign(b_eur) == sign(eas_beta_mr)]
stopifnot(sum(bin$same) == num("L987", "same_sign") - num("L967", "same_sign"), sum(!bin$same) == num("L987", "opposite_sign"))
bin[, state := factor(fifelse(eas_outcome_state == "EAS_COMPATIBLE", "Compatible", "Not significant"), levels = c("Compatible", "Not significant"))]
bin[, `:=`(lo_eas = eas_beta_mr - qnorm(0.975) * eas_se_mr, hi_eas = eas_beta_mr + qnorm(0.975) * eas_se_mr,
           lo_eur = b_eur - qnorm(0.975) * se_eur, hi_eur = b_eur + qnorm(0.975) * se_eur)]
# palindromic instruments (A/T or C/G), derived from the instrument alleles for every pair (the lead SNP of the pair joins the instrument table)
ins1 <- unique(ins[, .(gene, SNP, effect_allele, other_allele)])
pm <- merge(bin[, .(gene, outcome, lead_snps)], ins1, by.x = c("gene", "lead_snps"), by.y = c("gene", "SNP"))
stopifnot(nrow(pm) == nrow(bin))
pm[, pal := paste0(effect_allele, other_allele) %in% c("AT", "TA", "CG", "GC")]
bin <- merge(bin, pm[, .(gene, outcome, pal)], by = c("gene", "outcome"))
stopifnot(!anyNA(bin$pal), nrow(bin) == num("L987", "binary_evaluable"))

SIG <- 0.1
Tf <- function(x) asinh(x / SIG); Ti <- function(t) SIG * sinh(t)
tr <- scales::new_transform("asinh01", Tf, Ti)
allv <- c(bin$lo_eas, bin$hi_eas, bin$lo_eur, bin$hi_eur)
LT <- ceiling(max(abs(Tf(allv))) * 10) / 10 + 0.1                                          # symmetric limits, asinh units: one scale, so y = x is the diagonal
BRK <- c(-3, -1, -0.3, 0, 0.3, 1, 3); BRK <- BRK[abs(Tf(BRK)) < LT - 0.05]
quad <- data.table(xmin = c(-Inf, 0), xmax = c(0, Inf), ymin = c(0, -Inf), ymax = c(Inf, 0))
FILL <- c(Compatible = INK, "Not significant" = "white")
scat <- function(d, axis_id, ttl, nudge = NULL) {
  d <- copy(d)
  # S2-D26: a starting offset (asinh units) for named labels, so a label does not sit on an open point or run its leader under another label
  d[, `:=`(nx = 0, ny = 0)]
  if (!is.null(nudge)) for (k in names(nudge)) d[paste0(gene, "-", outcome) == k, `:=`(nx = nudge[[k]][1], ny = nudge[[k]][2])]
  ggplot(d, aes(b_eur, eas_beta_mr)) +
    geom_rect(data = quad, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = SHADE_OPPOSITE, colour = NA) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.35, colour = LINE_GREY) + geom_vline(xintercept = 0, linetype = "22", linewidth = 0.35, colour = LINE_GREY) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.35, colour = LINE_GREY) +
    geom_segment(aes(x = lo_eur, xend = hi_eur, y = eas_beta_mr, yend = eas_beta_mr), linewidth = 0.3, colour = "grey70") +
    geom_segment(aes(x = b_eur, xend = b_eur, y = lo_eas, yend = hi_eas), linewidth = 0.3, colour = "grey70") +
    geom_point(data = d[pal == TRUE], shape = 21, size = 3.6, stroke = 0.4, fill = NA, colour = INK) +
    geom_point(aes(fill = state), shape = 21, size = 1.8, colour = INK, stroke = 0.45) +
    geom_label_repel(data = d, aes(label = fifelse(state == "Compatible", paste0(gene, EN, fifelse(outcome == "Stroke", "stroke", outcome)), "")), size = ANNOT, family = BASE_FONT, seed = 11, nudge_x = d$nx, nudge_y = d$ny,
                     fill = alpha("white", 0.9), label.size = 0, label.padding = unit(0.6, "pt"), box.padding = 0.55, point.padding = 0.5, force = 4, force_pull = 0.6,
                     min.segment.length = 0, segment.size = 0.3, segment.colour = "grey25", max.overlaps = Inf, max.time = 4, max.iter = 100000) +
    annotate("text", x = Ti(-LT + 0.1), y = Ti(LT - 0.1), label = "asinh scale", hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, colour = "grey25") +
    scale_fill_manual(values = FILL, guide = "none") +
    scale_x_continuous(trans = tr, breaks = BRK, labels = ticks_set(paste0("S6_", axis_id, "_x"), BRK, Ti(c(-LT, LT)), paste0("S6 ", axis_id, " x, log-odds per SD of protein, asinh scale"))) +
    scale_y_continuous(trans = tr, breaks = BRK, labels = ticks_set(paste0("S6_", axis_id, "_y"), BRK, Ti(c(-LT, LT)), paste0("S6 ", axis_id, " y, log-odds per SD of protein, asinh scale"))) +
    coord_cartesian(xlim = Ti(c(-LT, LT)), ylim = Ti(c(-LT, LT)), expand = FALSE) +
    labs(x = "European estimate\nlog-odds per SD of protein", y = "East Asian estimate\nlog-odds per SD of protein", title = ttl) +
    theme_p31r() +
    theme(aspect.ratio = 1, plot.title = element_text(size = ANNOT_PT, face = "plain", hjust = 0, margin = margin(b = 3), lineheight = 1.05), plot.title.position = "plot",
          axis.title.x = element_text(lineheight = 0.9), axis.title.y = element_text(lineheight = 0.9))
}

# ---- a ---------------------------------------------------------------------------------------------------------------------------------------
tta <- paste0("Binary-outcome pairs: ", led("L987", "binary_evaluable", fmt0, where = "S6a title: binary evaluable pairs", anchor = "^Binary-outcome pairs: {} of 52"),
              " of ", led("L080", 1, fmt0, where = "S6a title: evaluable pairs", anchor = "of {} evaluable\\.$"), " evaluable.\n",
              "The ", led("L987", "egfr_evaluable", fmt0, where = "S6a title: eGFR pairs", anchor = "^The {} eGFR pairs"),
              " eGFR pairs are counted in Figure 5a,\nnot plotted.")
pa <- scat(bin, "a", tta, nudge = list("APOE-T2D" = c(2.2, -1.0)))

# ---- b: strict subset, its T2D part --------------------------------------------------------------------------------------------------------
st <- eas[eas_outcome_state_strict %in% c("EAS_COMPATIBLE", "EAS_INCONCLUSIVE")]
stopifnot(nrow(st) == num("L081", 1), sum(st$eas_outcome_state_strict == "EAS_COMPATIBLE") == num("L081", 2), setequal(unique(st$outcome), c("T2D", "eGFR")),
          all(st$sample_count_basis == "VARIANT_LEVEL"))
st_t2d <- st[outcome == "T2D"]
stopifnot(nrow(st_t2d) == nrow(st) - num("L967", "evaluable"), sum(st[outcome == "eGFR"]$eas_outcome_state_strict == "EAS_COMPATIBLE") == num("L967", "compatible"))
sb <- merge(bin[outcome == "T2D"], st_t2d[, .(gene, outcome, strict = eas_outcome_state_strict)], by = c("gene", "outcome"))
stopifnot(nrow(sb) == nrow(st_t2d), !anyNA(sb$strict))
sb[, state := factor(fifelse(strict == "EAS_COMPATIBLE", "Compatible", "Not significant"), levels = c("Compatible", "Not significant"))]
ttb <- paste0("Strict subset (per-variant sample sizes;\nT2D and eGFR files): ", led("L081", 1, fmt0, where = "S6b title: strict evaluable", anchor = "files\\): {} evaluable"),
              " evaluable, ", led("L081", 2, fmt0, where = "S6b title: strict compatible", anchor = "evaluable, {} compatible"), " compatible.\n",
              "Its T2D part is plotted; its ", led("L967", "evaluable", fmt0, where = "S6b title: eGFR pairs", anchor = "its {} eGFR pairs"), " eGFR pairs are not.")
pb <- scat(sb, "b", ttb, nudge = list("ARG1-T2D" = c(-1.0, 0.5)))

# ---- key (own strip, 7 pt) -----------------------------------------------------------------------------------------------------------------------
kab <- ggplot() +
  geom_point(aes(x = mm(3), y = 1), shape = 21, fill = INK, colour = INK, size = 1.8, stroke = 0.45) + TXT(aes(x = mm(5.5), y = 1, label = "Compatible"), hjust = 0) +
  geom_point(aes(x = mm(25), y = 1), shape = 21, fill = "white", colour = INK, size = 1.8, stroke = 0.45) + TXT(aes(x = mm(27.5), y = 1, label = "Not significant"), hjust = 0) +
  geom_point(aes(x = mm(52), y = 1), shape = 21, size = 3.6, stroke = 0.4, fill = NA, colour = INK) + TXT(aes(x = mm(56), y = 1, label = "Ring: palindromic instrument (A/T or C/G)"), hjust = 0) +
  TXT(aes(x = mm(108), y = 1, label = "Labeled: compatible pairs only"), hjust = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.4, 1.6), expand = FALSE) + theme_void()
top <- (pa + labs(tag = "a")) + (pb + labs(tag = "b")) + plot_layout(widths = c(1, 1))
fig <- top / wrap_elements(full = kab) + plot_layout(heights = c(90, 6))
frame_cols <- c("gene", "outcome", "evidence_tier", "eas_outcome_source", "eas_outcome_state", "b_eur", "se_eur", "eas_beta_mr", "eas_se_mr", "lo_eur", "hi_eur", "lo_eas", "hi_eas", "same", "state", "pal")
af3_save(fig, h_mm = 106, panels = c("a", "b"),
         panel_map = c(a = "outcome-side European against East Asian estimates, binary-outcome evaluable pairs; palindromic ring",
                       b = "strict per-variant-sample-size subset, T2D part"),
         frames = list(binary_outcome_pairs = bin[, ..frame_cols], strict_T2D = sb[, .(gene, outcome, eas_outcome_source, strict, state, b_eur, eas_beta_mr)]))

# =========================================================================================================== PAGE 2
unit_reset()
UNIT <- "FigureS6_p2"
af3_begin(UNIT)
cen <- af3_read("V4_local/staging/d2-eas-census/results_census.tsv", "c", sep = "\t")
ckf <- af3_read("V4_local/staging/d2-eas-census/results_ckbfree.tsv", "c", sep = "\t")
keyc <- function(d) paste(d$gene, d$outcome)
stopifnot(nrow(cen) == 48L, nrow(ckf) == 13L, !anyDuplicated(keyc(cen)), !anyDuplicated(keyc(ckf)), all(keyc(ckf) %in% keyc(cen)),
          all(!TRUEish(ckf$ckb_overlap)), all(ckf$outcome %in% c("T2D", "Stroke")))
# the composite: census rows, with the CKB-free re-run row swapped in for the same gene x outcome pair
cols <- c("gene", "outcome", "status", "resource", "accession", "ckb_overlap", "PP.H3", "PP.H4", "support", "call", "PP.H4_p12_1e-6", "support_p12_1e-6")
comp <- rbind(cen[!keyc(cen) %in% keyc(ckf), ..cols][, rerun := FALSE], ckf[, ..cols][, rerun := TRUE])
stopifnot(nrow(comp) == 48L, !anyDuplicated(keyc(comp)), sum(comp$rerun) == 13L)
comp[, `:=`(h4_def = suppressWarnings(as.numeric(PP.H4)), h4_con = suppressWarnings(as.numeric(`PP.H4_p12_1e-6`)), sup = TRUEish(support), rob = TRUEish(`support_p12_1e-6`))]
comp[, run := call != "NOT_TESTABLE_NO_EAS_GWAS"]
stopifnot(sum(comp$run) == 47L, !any(comp[run == FALSE, !is.na(h4_def)]), !anyNA(comp[run == TRUE, h4_def]), all(!TRUEish(comp[run == TRUE]$ckb_overlap)), all(!comp$rob | comp$sup))
# the composite reproduces the ledger: eight colocalized at the default prior, four robust to the conservative prior (L084)
stopifnot(sum(comp$sup) == num("L084", 1), sum(comp$rob) == num("L084", 2))
# H9: no row of the composite comes from EGFR_EAS_2P2M; its eGFR rows are the TWB2+BBJ meta-analysis, which H9 does not cover
stopifnot(!any(grepl("2P2M", comp$accession, ignore.case = TRUE)), all(comp[outcome == "eGFR" & run == TRUE]$accession == "TWB2_BBJ_eGFR_PMID39472450"),
          all(comp[run == TRUE & outcome %in% c("T2D", "Stroke")]$rerun))
comp[, status_f := factor(fifelse(!run, "Not run", fifelse(rob, "Robust to the conservative prior", fifelse(sup, "Colocalized, not robust", "Not colocalized"))),
                          levels = c("Robust to the conservative prior", "Colocalized, not robust", "Not colocalized", "Not run"))]
comp[, why := fcase(!run, "no East Asian outcome GWAS", sup, "", call == "NOT_SUPPORTED_H3", "distinct variants favored", call == "NOT_SUPPORTED_H1", "protein-only model favored",
                    call == "NOT_SUPPORTED_H4_below_rule", "PP.H4 below the rule", default = NA_character_)]
stopifnot(!anyNA(comp$why))
comp[, pair := paste0(gene, EN, fifelse(outcome == "Stroke", "stroke", outcome))]
comp[, grp := factor(fifelse(rerun, "Re-run with a CKB-free
outcome resource", "Outcome resource without CKB:
screen result, no re-run"),
                     levels = c("Re-run with a CKB-free
outcome resource", "Outcome resource without CKB:
screen result, no re-run"))]
comp[, ord := fifelse(is.na(h4_def), -1, h4_def)]
setorder(comp, grp, -ord)
comp[, y := rev(seq_len(.N)), by = grp]     # numeric row positions: the highest PP.H4 at the top of each facet (a discrete axis would take its order from the layers)
comp[, rowtxt := fifelse(sup, as.character(status_f), paste0(tolower(substr(as.character(status_f), 1, 1)), substring(as.character(status_f), 2), " (", why, ")"))]
cl <- melt(comp[run == TRUE, .(pair, y, grp, status_f, def = h4_def, con = h4_con)], id.vars = c("pair", "y", "grp", "status_f"), variable.name = "prior", value.name = "pp")
cl[, prior := factor(prior, levels = c("def", "con"), labels = c("Default prior", "Conservative prior"))]
XB <- c(0, 0.5, 1)
# R3 (USER): the protein-only wording in this panel, with this panel's own counts; each is a count of the rows plotted here, and the three reasons exhaust the row set
lit <- function(text, why, anchor = NULL, count = 1L) {
  REG$literals[[length(REG$literals) + 1]] <<- Filter(Negate(is.null), list(text = text, why = why, count = as.integer(count), anchor = anchor, not_ledger = "a count of the rows plotted in this panel"))
  text
}
nc_ <- comp[run == TRUE & sup == FALSE]
stopifnot(sum(nc_$why %in% c("protein-only model favored", "distinct variants favored", "PP.H4 below the rule")) == nrow(nc_))
subt <- paste0("Composite of the ancestry-matched screen and the re-runs: ", led("L084", 1, fmt0, where = "S6c subtitle: CKB-free composite colocalized", anchor = "^Composite of the ancestry-matched screen and the re-runs: {} colocalized"),
               " colocalized, ", led("L084", 2, fmt0, where = "S6c subtitle: CKB-free composite robust", anchor = "colocalized, {} robust"), " robust to the conservative prior.",
               "\nNot colocalized, ", lit(as.character(nrow(nc_)), "other: count of this panel's rows (not colocalized)", anchor = "^Not colocalized, {}:"),
               ": protein-only model favored (", lit(as.character(sum(nc_$why == "protein-only model favored")), "other: count of this panel's rows (protein-only model favored)", anchor = "protein-only model favored \\({}\\),"),
               "), distinct variants favored (", lit(as.character(sum(nc_$why == "distinct variants favored")), "other: count of this panel's rows (distinct variants favored)", anchor = "distinct variants favored \\({}\\),"),
               "), PP.H4 below the rule (", lit(as.character(sum(nc_$why == "PP.H4 below the rule")), "other: count of this panel's rows (PP.H4 below the rule)", anchor = "PP\\.H4 below the rule \\({}\\)\\.$"), ").")
pcc <- ggplot(cl, aes(pp, y)) +
  geom_vline(xintercept = num("L236"), linetype = "dashed", colour = LINE_GREY, linewidth = 0.3) +
  geom_segment(data = comp[run == TRUE & sup == TRUE], aes(x = h4_con, xend = h4_def, y = y, yend = y), inherit.aes = FALSE, colour = "grey70", linewidth = 0.4) +
  geom_point(data = cl[status_f == "Robust to the conservative prior"], aes(shape = prior), fill = N_[["dark"]], colour = N_[["dark"]], size = 1.8) +
  geom_point(data = cl[status_f == "Colocalized, not robust"], aes(shape = prior), fill = N_[["light"]], colour = N_[["dark"]], size = 1.8, stroke = 0.45) +
  geom_point(data = cl[status_f == "Not colocalized" & prior == "Default prior"], shape = 13, colour = N_[["dark"]], size = 1.9, stroke = 0.4) +
  geom_text(data = comp, aes(x = 1.06, y = y, label = rowtxt), inherit.aes = FALSE, hjust = 0, size = ANNOT, family = BASE_FONT) +
  geom_text(data = comp, aes(x = -0.07, y = y, label = pair), inherit.aes = FALSE, hjust = 1, size = ANNOT, family = BASE_FONT) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c("Default prior" = 21, "Conservative prior" = 24), name = NULL) +
  scale_x_continuous(breaks = XB, labels = ticks_set("S6_c_x", XB, c(-0.04, 1.9), "S6 c x, PP.H4", "linear"), expand = expansion(0), guide = guide_axis(cap = "upper"),
                     sec.axis = dup_axis(breaks = num("L236"), labels = "PP.H4 threshold", name = NULL)) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.6, 0.6))) +
  coord_cartesian(xlim = c(-0.04, 1.9), clip = "off") +
  labs(x = "PP.H4 with the outcome resource (circle, default prior; triangle, conservative prior)", y = NULL, title = subt) +
  theme_p31r() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(2, 2, 2, 24, "mm"), axis.text.x.top = element_text(size = ANNOT_PT, hjust = 1, margin = margin(b = 1)), axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(), plot.title = element_text(size = ANNOT_PT, face = "plain", hjust = 0, margin = margin(b = 3)), plot.title.position = "plot",
        strip.text.y = element_text(angle = -90, size = ANNOT_PT, lineheight = 0.95), strip.clip = "off", panel.spacing.y = unit(5, "pt"))
kc_ <- ggplot() +
  geom_point(aes(x = mm(3), y = 1), shape = 21, fill = N_[["dark"]], colour = N_[["dark"]], size = 1.9) + TXT(aes(x = mm(5.5), y = 1, label = "Robust"), hjust = 0) +
  geom_point(aes(x = mm(20), y = 1), shape = 21, fill = N_[["light"]], colour = N_[["dark"]], size = 1.9, stroke = 0.45) + TXT(aes(x = mm(22.5), y = 1, label = "Colocalized, not robust"), hjust = 0) +
  geom_point(aes(x = mm(58), y = 1), shape = 13, colour = N_[["dark"]], size = 2.0, stroke = 0.4) + TXT(aes(x = mm(60.5), y = 1, label = "Not colocalized"), hjust = 0) +
  TXT(aes(x = mm(86), y = 1, label = "Not run: no marker, the reason in words"), hjust = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.4, 1.6), expand = FALSE) + theme_void()
fig2 <- (pcc + labs(tag = "c")) / wrap_elements(full = kc_) + plot_layout(heights = c(156, 6))
af3_save(fig2, h_mm = 172, panels = "c",
         panel_map = c(c = "CKB-free composite of the ancestry-matched colocalization census (census rows with the 13 CKB-free re-runs swapped in)"),
         frames = list(ckbfree_composite = comp[, .(gene, outcome, resource, accession, rerun, run, call, h4_def, h4_con, sup, rob, status_f, why)]))
