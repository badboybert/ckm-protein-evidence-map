# Paper 3.1 rewrite (R2) -- Figure 5: transfer of the European pairs to East Asian populations (plan v3, section 5; S2-D24 F5/F16/F17).
#   a  two flows. European instrument, East Asian outcome data: all pairs -> evaluable | not evaluable; the evaluable pairs are
#      then split TWICE, as two separate stacked bars of equal length (by direction relative to the European estimate; by
#      East Asian significance). They are two partitions of the SAME pairs, not nested subsets: the 40 inconclusive pairs are
#      33 same-sign + 7 opposite-sign (ledger L987), so they are not a subset of the 45 same-sign pairs.
#      East Asian instrument and outcome data: eligible -> testable -> colocalized -> robust; non-colocalizing by hypothesis.
#   b  ancestry-matched paired forest (European and East Asian estimate with 95% CI per pair): the six binary-outcome pairs
#      that colocalized and FGFR4-T2D; FGF5-eGFR as PP.H4 only, the scale reason in words (estimates on different scales; not
#      drawn, and no direction is printed: neither the Results nor SM S-R13 states one; S2-D25 visual-check finding).
# NOT in this figure (plan v3 section 5, "Must not show"): the binary-outcome asinh scatter (moved to Additional file 3, Figure S6,
#   F17); any East Asian eGFR magnitude from EGFR_EAS_2P2M (hold H9: main and AF3 alike); post hoc tier colours (F5); SMPD3 as an
#   example; "replicated in East Asians"; an ancestry-specific reversal at TGFB1.
# Sources (APPROVED_SOURCE_MANIFEST.tsv, Figure_5 rows; nothing else is opened with p31r_read, so G12 can check every pin):
#   figures/data/Figure_5/eas_partition_members.tsv  (pair memberships, no magnitude columns: H9)
#   V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv  (ancestry-matched, panels a right and b)
#   the claims ledger (every printed number, read at the point of printing).
# Encodings follow figures/P31R_FIGURE_SPEC.md: East Asian colocalization is NEUTRAL shading (robust = dark; colocalized but not
# robust = light with a dark outline), never a tier hue; ancestry is SHAPE (circle European, diamond East Asian); the left flow in a
# carries no fill meaning (all bars neutral grey, meaning in words and position).
# The registry (figures/registry/Figure_5.json) is WRITTEN by this script AT THE POINT OF PRINTING (one entry per printed
# occurrence, ledger element idx, count drawn, anchor), as the earlier generator did.
# Run: Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_Fig5.R"

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

# ---- sources (manifest rows for Figure_5) -------------------------------------------------------------------------
F_MEM  <- file.path(P31R_FIGS, "data/Figure_5/eas_partition_members.tsv")
F_S23N <- file.path(V3_ROOT, "V4_local/source_data/figures/S23_four_quadrant_eur_eas/S23_native_census__source.csv")
F_INS  <- file.path(V3_ROOT, "V4/source_data/tables_from_rds/phase1_instruments__ukb_cis_instruments.csv")   # European instrument rsIDs under the panel b pair labels
mem <- p31r_read(F_MEM)
s23 <- p31r_read(F_S23N)
ins5 <- p31r_read(F_INS)

TRUEish <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1")
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))      # for assertions only (nothing printed)
N_ <- PAL_NEUTRAL                                                       # neutral greys, by name
OUT_LAB <- function(o) ifelse(o == "Stroke", "stroke", o)               # spec: outcome abbreviations; "stroke" lower case
EN <- "\u2013"

# ---- registry: written at the point of printing ------------------------------------------------------------------
REG <- list(ledger = list(), literals = list())
# led(): fetch ONE ledger element, print it with fmt, and record {text, id, idx, count, anchor} from THIS call.
led <- function(id, idx = NULL, fmt = function(v) as.character(v), where, anchor, count = 1L) {
  v <- p31r_ledger(id, idx)
  stopifnot(length(v) == 1, is.numeric(v))
  text <- fmt(v)
  REG$ledger[[length(REG$ledger) + 1]] <<- c(list(text = text, id = id), if (!is.null(idx)) list(idx = idx),
                                              list(count = as.integer(count), anchor = anchor, where = where))
  text
}
# ticks(): the labels a scale draws, from its own break vector (p31r_lab: U+2212 minus); every break must lie inside the limits
ticks <- function(axis, breaks, lims, transform, why) {
  stopifnot(all(breaks > lims[1] & breaks < lims[2]), !anyDuplicated(breaks))
  labs <- p31r_lab(breaks)
  for (t in labs) REG$literals[[length(REG$literals) + 1]] <<- list(text = t, why = paste0("axis_tick: ", why), axis = axis,
                                                                    transform = transform, count = 1L, anchor = "^{}$")
  labs
}
# "<label> <number>" printed as ONE line; the anchor is the label (no digits in it) followed by the number
lab_n <- function(label, id, idx = NULL, where, fmt = function(v) as.character(v)) {
  paste(label, led(id, idx, fmt, where = where, anchor = paste0("^", gsub("([.()])", "\\\\\\1", label), " {}$")))
}
n_lab <- function(label, id, idx = NULL, where) {   # "<number> <label>"
  paste(led(id, idx, where = where, anchor = paste0("^{} ", gsub("([.()])", "\\\\\\1", label), "$")), label)
}

# ================================================================================================ partition memberships (a left)
stopifnot(nrow(mem) == 81, !anyDuplicated(mem[, paste(gene, outcome)]),
          !any(c("beta_mr_eas", "se_mr_eas", "eas_beta_mr", "eas_se_mr", "p_eas") %in% names(mem)))   # H9: the frame holds no magnitude
ev   <- mem[eas_outcome_state != "EAS_NOT_EVALUABLE"]
st   <- ev[, eas_outcome_state]
sg   <- ev[, sign_vs_european]
stopifnot(all(sg %in% c("same sign", "opposite sign")), all(mem[eas_outcome_state == "EAS_NOT_EVALUABLE", sign_vs_european] == "not evaluable"))
tab <- list(evaluable = nrow(ev), same = sum(sg == "same sign"), opposite = sum(sg == "opposite sign"),
            comp = sum(st == "EAS_COMPATIBLE"), inc = sum(st == "EAS_INCONCLUSIVE"), dis = sum(st == "EAS_DISCORDANT"),
            comp_same = sum(st == "EAS_COMPATIBLE" & sg == "same sign"), comp_opp = sum(st == "EAS_COMPATIBLE" & sg == "opposite sign"),
            inc_same = sum(st == "EAS_INCONCLUSIVE" & sg == "same sign"), inc_opp = sum(st == "EAS_INCONCLUSIVE" & sg == "opposite sign"),
            dis_same = sum(st == "EAS_DISCORDANT" & sg == "same sign"), dis_opp = sum(st == "EAS_DISCORDANT" & sg == "opposite sign"),
            notev = nrow(mem[eas_outcome_state == "EAS_NOT_EVALUABLE"]), egfr = sum(ev$h9_hold_file == "TRUE" | ev$h9_hold_file == TRUE))
stopifnot(tab$comp + tab$inc + tab$dis == tab$evaluable, tab$same + tab$opposite == tab$evaluable,
          tab$evaluable + tab$notev == num("L010"))
# every count the figure prints is the frame's tally AND the ledger's value (L080 / L224 / L987 / L967 / L014)
stopifnot(identical(as.numeric(c(tab$evaluable, tab$comp, tab$inc, tab$dis, tab$notev)), as.numeric(unlist(p31r_ledger("L080")))),
          tab$same == num("L224", 1), tab$evaluable == num("L224", 2),
          tab$evaluable == num("L987", "evaluable"), tab$same == num("L987", "same_sign"), tab$opposite == num("L987", "opposite_sign"),
          tab$comp_same == num("L987", c("compatible", "same")), tab$comp_opp == num("L987", c("compatible", "opposite")),
          tab$inc_same == num("L987", c("inconclusive", "same")), tab$inc_opp == num("L987", c("inconclusive", "opposite")),
          tab$dis_same == num("L987", c("discordant", "same")), tab$dis_opp == num("L987", c("discordant", "opposite")),
          tab$notev == num("L987", "not_evaluable"), tab$egfr == num("L987", "egfr_evaluable"), tab$egfr == num("L967", "evaluable"))
# the three claims the figure draws in words: (i) none of the opposite-sign pairs is significant; (ii) the compatible pairs all
# keep the European sign; (iii) none of the three Convergent pairs was evaluable
stopifnot(tab$comp_opp == 0, tab$dis == 0, tab$inc_opp == tab$opposite,
          sum(mem$evidence_tier == "Convergent") == num("L014", "Convergent"),
          all(mem[evidence_tier == "Convergent", eas_outcome_state] == "EAS_NOT_EVALUABLE"))
# the 14 eGFR pairs: counted (L967), the held source is named; this frame and figure carry no eGFR magnitude
stopifnot(setequal(unique(ev[h9_hold_file == "TRUE", eas_outcome_source]), "EGFR_EAS_2P2M"),
          all(ev[h9_hold_file == "TRUE", sign_vs_european] == "same sign"))

# ================================================================================================ ancestry-matched (a right, b)
cen <- s23[set == "census_network"]
stopifnot(all(cen$eur_call == "H4"), all(TRUEish(cen$eur_support)))
tst <- cen[eas_call != "NOT_TESTABLE_NO_EAS_GWAS"]
col <- tst[eas_call == "SUPPORTED"]
stopifnot(all(TRUEish(tst$eas_support) == (tst$eas_call == "SUPPORTED")))
rob <- col[TRUEish(eas_support_p12_1e6)]
stopifnot(identical(as.numeric(c(nrow(cen), nrow(tst), nrow(col), nrow(rob))), as.numeric(unlist(p31r_ledger("L082")))))
# "robust" = the FULL rule (PP.H4 >= 0.80 AND PP.H4/PP.H3 >= 5) at p12 = 1e-6, not PP.H4 alone (ABO-T2D: PP.H4 0.801, ratio 4.1)
full_rule <- function(h4, h3) {
  h4 <- as.numeric(h4); h3 <- as.numeric(h3); s <- num("L241") / p31r_ledger("L240", 3)
  (s * h4 / (1 - (1 - s) * h4) >= num("L236")) & (s * h4 >= num("L237") * h3)
}
h4c <- as.numeric(col$eas_coloc_H4)
stopifnot(all(full_rule(col$eas_coloc_H4, col$eas_coloc_H3) == TRUEish(col$eas_support_p12_1e6)),
          all(abs(as.numeric(col$eas_coloc_H4_p12_1e6) - 0.1 * h4c / (1 - 0.9 * h4c)) <= 0.1 / (1 - 0.9 * h4c)^2 * 5e-4 + 5e-4))
pair <- function(d) paste0(d$protein, "\u2013", d$outcome)
l083 <- p31r_ledger("L083")
stopifnot(length(l083) == 2, setequal(pair(col), unlist(l083[[1]])), setequal(pair(rob), unlist(l083[[2]])))
n_h1 <- sum(tst$eas_call == "NOT_SUPPORTED_H1"); n_h3 <- sum(tst$eas_call == "NOT_SUPPORTED_H3")
below <- tst[eas_call == "NOT_SUPPORTED_H4_below_rule"]
stopifnot(n_h1 == num("L211", 1), n_h3 == num("L211", 2), nrow(below) == 1, below$protein == "DNER", below$outcome == "eGFR",
          abs(as.numeric(below$eas_coloc_H4) - num("L255")) < 5e-4,
          n_h1 + n_h3 + nrow(below) + nrow(col) == nrow(tst))
# rows of b: the colocalizing BINARY-outcome network pairs + FGFR4-T2D (off-network row: the instrument was identified only in
# this study's East Asian analysis; S2-D25 wording)
fg <- s23[set == "pilot_offnetwork" & protein == "FGFR4" & outcome == "T2D"]
stopifnot(nrow(fg) == 1, TRUEish(fg$eas_support), TRUEish(fg$eas_support_p12_1e6), full_rule(fg$eas_coloc_H4, fg$eas_coloc_H3),
          startsWith(fg$eas_instr_rule, "pilot_old"), all(startsWith(col$eas_instr_rule, "census_v4")))
cb <- rbind(col[outcome != "eGFR"], fg, fill = TRUE)
stopifnot(nrow(cb) == 7, nrow(col[outcome != "eGFR"]) == 6)
cb[, `:=`(eur = as.numeric(eur_beta_ha), eur_se = as.numeric(eur_se), eas = as.numeric(eas_beta_ha), eas_se = as.numeric(eas_se))]
cb[, robust := TRUEish(eas_support_p12_1e6)]
cb[, lab := paste0(protein, EN, OUT_LAB(outcome))]
# instrument rsIDs printed under each pair label (review R23): European = the UKB-PPP cis instrument (one per protein), East Asian = the pair's own instrument
stopifnot(all(table(ins5$gene)[cb$protein] == 1))
cb[, `:=`(eur_rs = ins5$SNP[match(protein, ins5$gene)], eas_rs = eas_instr)]
stopifnot(all(grepl("^rs[0-9]+$", c(cb$eur_rs, cb$eas_rs))))
# values the Results quote, checked against the plotted columns
ce <- cb[protein == "CELSR2"]; tg <- cb[protein == "TGFB1"]; fgc <- cb[protein == "FGFR4"]
stopifnot(abs(ce$eas - num("L212", 3)) < 1e-4, abs(ce$eur - num("L212", 4)) < 1e-4,
          abs(tg$eur - num("L085", 1)) < 1e-4, abs(tg$eas - num("L085", 2)) < 1e-4, tg$eur > 0, tg$eas < 0, tg$robust,
          abs(as.numeric(tg$eas_coloc_H4) - num("L214", 1)) < 5e-4,
          abs(as.numeric(fgc$eas) - num("L086", 1)) < 5e-5, 2 * pnorm(-abs(fgc$eur / fgc$eur_se)) > 0.25, abs(fgc$eur / fgc$eur_se) < 1.2)
stopifnot(sprintf("%.2f", 2 * pnorm(-abs(fgc$eur / fgc$eur_se))) == sprintf("%.2f", num("L215", 2)))
# colocalization states drawn: robust = CELSR2, TGFB1, FGFR4 here (+ FGF5-eGFR, a sign-only row); not robust = the other four
stopifnot(identical(sort(cb[robust == TRUE, protein]), c("CELSR2", "FGFR4", "TGFB1")),
          identical(sort(cb[robust == FALSE, protein]), c("ABO", "GALNT3", "IL6R", "SPON1")))
# TGFB1-CAD: the European colocalization is not robust (L214[3] = 0.753 < 0.80); the East Asian one is (dark fill)
stopifnot(num("L214", 3) < num("L236"))

# FGF5-eGFR: PP.H4 only, with the scale reason in words (F16). The scales differ (European -ln eGFR per SD of protein; East Asian
# per SD of rank-inverse-normal eGFR, TWB2+BBJ), so no estimate is drawn and no common scale is derived. No direction is printed
# (S2-D25: the Results and SM S-R13 carry none), so the eGFR* sign convention is not read and no sign is checked here.
f5 <- s23[set == "census_network" & protein == "FGF5" & outcome == "eGFR"]
stopifnot(nrow(f5) == 1, TRUEish(f5$eas_support), TRUEish(f5$eas_support_p12_1e6),
          abs(as.numeric(f5$eas_coloc_H4) - num("L213", 1)) < 5e-4, abs(as.numeric(f5$eas_coloc_H4_p12_1e6) - num("L213", 2)) < 5e-4,
          grepl("NA\\(scale\\)", f5$het_z), grepl("39472450", f5$eas_outcome_accession), nchar(f5$eas_outcome_resource) > 0)

# ============================================================================================================ shared drawing
REF <- list(colour = N_[["mid"]], linewidth = 0.3, linetype = "22")      # reference lines: dashed grey (0)
blank_axes <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(),
                    axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                    axis.line.x = element_blank(), axis.line.y = element_blank(), panel.grid = element_blank(),
                    panel.border = element_blank())
TXT <- function(d, ...) geom_text(data = d, aes(x = x, y = y, label = label, hjust = hj, vjust = vj, colour = I(tc)),
                                  size = ANNOT, family = BASE_FONT, lineheight = 0.95, ...)
tx  <- function(x, y, label, hj = 0, vj = 0.5, tc = N_[["dark"]]) data.table(x = x, y = y, label = label, hj = hj, vj = vj, tc = tc)
LEFT_MM <- 5                                                              # room for the panel tag; 1 data unit = 1 mm after it

# =============================================================================================================== a (mm canvas)
WA <- 164.7; HA <- 72.5; TOP_MM <- 3.5
K <- 1.5; X0 <- 10                                                         # 1.5 mm per pair: the evaluable pairs span 52 x 1.5 = 78 mm
W_ <- "white"
box <- function(x0, x1, y0, y1, fill = "white", colour = N_[["dark"]]) data.table(x0 = x0, x1 = x1, y0 = y0, y1 = y1, fill = fill, colour = colour)
seg <- function(x, xend, y, yend, arrow = FALSE, lt = "solid") data.table(x = x, xend = xend, y = y, yend = yend, arrow = arrow, lt = lt)

# left flow --------------------------------------------------------------------------------------------------------------
GREY <- N_[["other"]]
n_same <- tab$same; n_opp <- tab$opposite
xs <- function(n) X0 + n * K
L_boxes <- rbind(
  box(2, 32, 6, 12), box(2, 32, 17, 23), box(40, 84, 17, 26),
  # bar 1: by direction (same-sign pairs first, then opposite-sign)
  box(xs(0), xs(n_same), 32, 38, GREY), box(xs(n_same), xs(n_same + n_opp), 32, 38, GREY),
  # bar 2: by significance, in the SAME pair order (same-sign first): compatible (all same sign), inconclusive (same then opposite)
  box(xs(0), xs(tab$comp), 52, 58, GREY), box(xs(tab$comp), xs(tab$comp + tab$inc), 52, 58, GREY))
stopifnot(abs(xs(tab$comp + tab$inc) - xs(n_same + n_opp)) < 1e-9)           # the two bars are the same length: the same pairs
L_seg <- rbind(
  seg(17, 17, 12, 17, TRUE),                                                  # all pairs -> evaluable
  seg(17, 62, 14.5, 14.5), seg(62, 62, 14.5, 17, TRUE),                       # all pairs -> not evaluable
  seg(5, 5, 23, 55),                                                        # trunk: the same evaluable pairs are split twice
  seg(5, X0 - 0.6, 35, 35, TRUE), seg(5, X0 - 0.6, 55, 55, TRUE),
  # the 33 / 7 division inside the inconclusive segment, aligned under the opposite-sign segment of bar 1
  seg(xs(n_same), xs(n_same), 52, 58, lt = "13"))
L_txt <- rbind(
  tx(2, 2.5, "European instrument, East Asian outcome data"),
  tx(17, 9, lab_n("All pairs", "L010", NULL, "a: all pairs"), 0.5),
  tx(17, 20, lab_n("Evaluable", "L080", 1, "a: evaluable"), 0.5),
  tx(62, 20.6, c(paste0(lab_n("Not evaluable", "L080", 5, "a: not evaluable"), "\nincludes all ",
                        led("L014", "Convergent", where = "a: Convergent pairs, none evaluable", anchor = "^includes all {} Convergent pairs$"),
                        " Convergent pairs")), 0.5, 0.5),
  tx(X0, 29.2, "By direction relative to the European estimate"),
  tx(mean(c(xs(0), xs(n_same))), 35, lab_n("Same sign", "L224", 1, "a: evaluable, European sign kept"), 0.5),
  tx(xs(n_same + n_opp), 39.6, paste0(lab_n("Opposite sign", "L987", "opposite_sign", "a: opposite sign"), "\nnone significant"), 1, 1),
  tx(X0, 47.6, "By East Asian significance, the same pairs"),
  tx(mean(c(xs(0), xs(tab$comp))), 55, lab_n("Compatible", "L080", 2, "a: compatible"), 0.5),
  tx(mean(c(xs(tab$comp), xs(tab$comp + tab$inc))), 55, lab_n("Inconclusive", "L080", 3, "a: inconclusive"), 0.5),
  tx(xs(tab$comp + tab$inc) + 2, 55, lab_n("Discordant", "L080", 4, "a: discordant"), 0, 0.5),
  tx(mean(c(xs(tab$comp), xs(n_same))), 59.4, n_lab("same sign", "L987", c("inconclusive", "same"), "a: inconclusive, same sign"), 0.5, 1),
  tx(mean(c(xs(n_same), xs(n_same + n_opp))), 59.4, n_lab("opposite sign", "L987", c("inconclusive", "opposite"), "a: inconclusive, opposite sign"), 0.5, 1),
  tx(2, 65.4, paste0(n_lab("of the evaluable pairs are eGFR pairs in a single East Asian eGFR file", "L967", "evaluable",
                             "a: eGFR pairs in the held file"), "\n(standardized effects; compared by sign and p only)"), 0, 1))

# right flow -------------------------------------------------------------------------------------------------------------
RX0 <- 113; RX1 <- 153
R_boxes <- rbind(
  box(RX0, RX1, 6, 12), box(RX0, RX1, 17, 23), box(RX0, RX1, 28, 34, N_[["light"]], N_[["dark"]]), box(RX0, RX1, 39, 45, N_[["dark"]], N_[["dark"]]),
  box(RX0 + 3, WA - 1.2, 49.5, 55), box(RX0 + 3, WA - 1.2, 57, 62.5), box(RX0 + 3, WA - 1.2, 64.5, 72))
R_seg <- rbind(
  seg(133, 133, 12, 17, TRUE), seg(133, 133, 23, 28, TRUE), seg(133, 133, 34, 39, TRUE),
  seg(RX0, RX0 - 2.5, 20, 20), seg(RX0 - 2.5, RX0 - 2.5, 20, 68.25), seg(RX0 - 2.5, RX0 + 2.4, 52.25, 52.25, TRUE),
  seg(RX0 - 2.5, RX0 + 2.4, 59.75, 59.75, TRUE), seg(RX0 - 2.5, RX0 + 2.4, 68.25, 68.25, TRUE))
R_txt <- rbind(
  tx(RX0 - 2.5, 2.5, "East Asian instrument and outcome data"),
  tx(133, 9, lab_n("Eligible pairs", "L082", 1, "a: eligible"), 0.5),
  tx(133, 20, lab_n("Testable", "L082", 2, "a: testable"), 0.5),
  tx(133, 31, lab_n("Colocalized", "L082", 3, "a: colocalized"), 0.5),
  tx(133, 42, lab_n("Robust", "L082", 4, "a: robust"), 0.5, tc = W_),
  tx(RX0 + 3, 47.2, "Testable, not colocalized:"),
  tx(RX0 + 4.5, 52.25, lab_n("No outcome association", "L211", 1, "a: no outcome association")),
  tx(RX0 + 4.5, 59.75, lab_n("Distinct variants", "L211", 2, "a: distinct variants")),
  tx(RX0 + 4.5, 68.25, paste0("Below the rule: DNER", EN, "eGFR,\nPP.H4 ", led("L255", NULL, function(v) sprintf("%.3f", v), "a: below the rule", anchor = "^PP\\.H4 {}$"))))

Bx <- rbind(L_boxes, R_boxes); Sg <- rbind(L_seg, R_seg); Tx <- rbind(L_txt, R_txt)
pa <- ggplot() +
  geom_rect(data = Bx, aes(xmin = x0, xmax = x1, ymin = -y1, ymax = -y0, fill = I(fill), colour = I(colour)), linewidth = 0.3) +
  geom_segment(data = Sg[arrow == FALSE & lt == "solid"], aes(x = x, xend = xend, y = -y, yend = -yend), colour = N_[["dark"]], linewidth = 0.3) +
  geom_segment(data = Sg[arrow == FALSE & lt != "solid"], aes(x = x, xend = xend, y = -y, yend = -yend), colour = N_[["dark"]], linewidth = 0.4, linetype = "22") +
  geom_segment(data = Sg[arrow == TRUE], aes(x = x, xend = xend, y = -y, yend = -yend), colour = N_[["dark"]], linewidth = 0.3,
               arrow = arrow(length = unit(1.3, "mm"), angle = 22, type = "open")) +
  TXT(transform(Tx, y = -y)) +
  coord_cartesian(xlim = c(0, WA), ylim = c(-HA, 0), expand = FALSE, clip = "off") +
  labs(tag = "a") + theme_p31r() + blank_axes +
  theme(plot.margin = margin(TOP_MM, 0, 0, LEFT_MM, "mm"), plot.tag.position = c(0, 1), plot.tag.location = "plot")

# =============================================================================================================== b (mm rows)
# three aligned plots on one row axis: pair labels | forest | notes, under a key. Row positions are in mm (y down).
XL <- c(-0.52, 0.26)
TK <- c(-0.4, -0.2, 0, 0.2)
ROWS <- data.table(lab = c("ABO", "CELSR2", "GALNT3", "IL6R", "SPON1", "TGFB1", "FGFR4"),
                   y   = c(4, 10.5, 17, 23.5, 30, 38, 51))
cb[, y := ROWS$y[match(protein, ROWS$lab)]]
stopifnot(!anyNA(cb$y), nrow(cb) == nrow(ROWS))
Y_F5 <- 63; Y_S1 <- 46; Y_S2 <- 56.5; HB <- 68
DY <- 1.45
lo <- function(b, s) b - 1.96 * s; hi <- function(b, s) b + 1.96 * s
stopifnot(min(lo(cb$eur, cb$eur_se), lo(cb$eas, cb$eas_se)) > XL[1], max(hi(cb$eur, cb$eur_se), hi(cb$eas, cb$eas_se)) < XL[2])
FILL <- c(`TRUE` = N_[["dark"]], `FALSE` = N_[["light"]])
long <- rbind(cb[, .(y = y - DY, x = eur, xmin = lo(eur, eur_se), xmax = hi(eur, eur_se), robust, anc = "European")],
              cb[, .(y = y + DY, x = eas, xmin = lo(eas, eas_se), xmax = hi(eas, eas_se), robust, anc = "East Asian")])
long[, `:=`(fill = FILL[as.character(robust)], shape = c(European = 21, `East Asian` = 23)[anc])]
SEPS <- c(Y_S1, Y_S2)
sep_layer <- function() geom_hline(yintercept = SEPS, colour = REF$colour, linewidth = REF$linewidth, linetype = REF$linetype)
rowscale <- scale_y_reverse(limits = c(HB, 0), expand = c(0, 0))

# pair labels
lab_df <- rbind(data.table(kind = "pair label", x = 0, y = cb$y - 1.25, label = cb$lab, tc = N_[["dark"]]),
                data.table(kind = "instrument rsIDs, European / East Asian", x = 0, y = cb$y + 1.35, label = paste0(cb$eur_rs, " / ", cb$eas_rs), tc = "grey35"),
                data.table(kind = "pair label", x = 0, y = Y_F5, label = paste0("FGF5", EN, "eGFR"), tc = N_[["dark"]]))
pL <- ggplot() + geom_text(data = lab_df, aes(x = x, y = y, label = label, colour = I(tc)), hjust = 0, vjust = 0.5, size = ANNOT, family = BASE_FONT) +
  sep_layer() + rowscale + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_p31r() + blank_axes + theme(plot.margin = margin(0, 0, 0, LEFT_MM, "mm"))

# forest
f5_lines <- paste0("Estimates not shown\nPP.H4 ",
                   led("L213", 1, function(v) sprintf("%.3f", v), "b: FGF5-eGFR PP.H4 (both priors)", anchor = "^PP\\.H4 {}$"))
pF <- ggplot(long) + sep_layer() +
  annotate("segment", x = 0, xend = 0, y = 0, yend = Y_S2, colour = REF$colour, linewidth = REF$linewidth, linetype = REF$linetype) +
  geom_segment(aes(x = xmin, xend = xmax, y = y, yend = y), colour = N_[["mid"]], linewidth = 0.4) +
  geom_point(aes(x = x, y = y, shape = I(shape), fill = I(fill)), size = 2.3, stroke = 0.45, colour = N_[["dark"]]) +
  annotate("text", x = mean(XL), y = Y_F5, label = f5_lines, hjust = 0.5, vjust = 0.5, size = ANNOT, family = BASE_FONT, lineheight = 0.95) +
  rowscale +
  scale_x_continuous(breaks = TK, labels = ticks("5b_x", TK, XL, "linear", "b, effect on the log-odds scale per SD of protein"), expand = c(0, 0)) +
  coord_cartesian(xlim = XL, clip = "off") +
  labs(x = "Log-odds per SD of protein (95% CI)") +
  theme_p31r() + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
                       plot.margin = margin(0, 2, 0, 0, "mm"))

# notes
tg_note <- paste0("Instruments ", led("L085", 4, function(v) sprintf("%.0f", v), "b: TGFB1-CAD instruments kb apart",
                                      anchor = "^Instruments {} kb apart: instrument-dependent estimates;$"),
                  " kb apart: instrument-dependent estimates;\nnot evidence of ancestry-specific reversal.\nEuropean colocalization not robust (PP.H4 ",
                  led("L214", 3, function(v) sprintf("%.3f", v), "b: TGFB1-CAD European PP.H4, conservative prior", anchor = "^European colocalization not robust \\(PP\\.H4 {}\\)$"), ")")
fg_note <- paste0("From the East Asian screen; European p ", led("L215", 2, function(v) sprintf("%.2f", v), "b: FGFR4-T2D European p", anchor = "^From the East Asian screen; European p {} \\(non-rejection\\)$"),
                  " (non-rejection)")
f5_note <- "Scales differ: European \u2212(natural-log eGFR) per SD vs\nEast Asian per SD of rank-inverse-normal eGFR\n(Taiwan\u2013Japan meta-analysis); magnitudes not compared"
note_df <- data.table(x = 0, y = c(cb[protein == "TGFB1", y], cb[protein == "FGFR4", y], Y_F5), label = c(tg_note, fg_note, f5_note))
pR <- ggplot() + geom_text(data = note_df, aes(x = x, y = y, label = label), hjust = 0, vjust = 0.5, size = ANNOT, family = BASE_FONT, lineheight = 0.95) +
  sep_layer() + rowscale + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_p31r() + blank_axes + theme(plot.margin = margin(0, 0, 0, 2, "mm"))

# key (inside the graphic): ancestry = shape, fill = East Asian colocalization of the pair, bar = 95% CI
KY_H <- 9
ky_pts <- data.table(x = c(3, 31, 3, 74), y = c(2.2, 2.2, 6.4, 6.4), shape = c(21, 23, 21, 21),
                     fill = c(N_[["mid"]], N_[["mid"]], N_[["dark"]], N_[["light"]]))
ky_txt <- rbind(tx(5.5, 2.2, "European"), tx(33.5, 2.2, "East Asian"), tx(52, 2.2, "bars, 95% confidence intervals"),
                tx(5.5, 6.4, "East Asian colocalization robust to the conservative prior"),
                tx(76.5, 6.4, "colocalized, not robust"))
pK <- ggplot() + geom_point(data = ky_pts, aes(x = x, y = -y, shape = I(shape), fill = I(fill)), size = 2.3, stroke = 0.45, colour = N_[["dark"]]) +
  TXT(transform(ky_txt, y = -y)) +
  scale_x_continuous(limits = c(0, 164.7), expand = c(0, 0)) + scale_y_continuous(limits = c(-KY_H, 0), expand = c(0, 0)) +
  labs(tag = "b") + theme_p31r() + blank_axes + theme(plot.margin = margin(0, 0, 0, LEFT_MM, "mm"), plot.tag.position = c(0, 1), plot.tag.location = "plot")

AX_H <- 10.3                                                                   # tick labels + title under the forest
rowsb <- (pL | pF | pR) + plot_layout(widths = c(33, 57, 75))   # label column widened for the rsID line
# ONE patchwork (no nested wrap_elements): heights in absolute mm, so 1 data unit of the mm canvases = 1 mm on the page
# absolute heights in patchwork size the PANELS; margins (3.5 mm above a) and the forest axis (AX_H, measured on the PDF) add to them,
# and patchwork's default outer margin is removed so the canvas width is the page width less LEFT_MM
fig <- (pa / pK / rowsb) + plot_layout(heights = unit(c(HA, KY_H, HB), "mm")) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
H_MM <- TOP_MM + HA + KY_H + HB + AX_H
# ======================================================================================================== plotted frames (what the panels draw)
# One frame per data table a panel plots, written beside the figure as af3/source_data/Figure_5__<name>.csv and listed in the registry under `frames`
# (the file naming and registry key of the AF3 units). Panel a is drawn on a millimetre canvas (y down from the top of the panel; a pair is K mm wide in
# the bars). Panel b: the plotted estimates carry no East Asian eGFR magnitude (hold H9): the eGFR pair is not drawn as an estimate and is checked below.
stopifnot(nrow(long) == 2 * nrow(cb), !any(cb$outcome == "eGFR"), !any(grepl("2P2M", cb$eas_outcome_accession)))   # panel b estimates: no eGFR row, no held-source accession
fr_boxes <- copy(Bx)[, row := seq_len(.N)]
fr_segs <- copy(Sg)[, row := seq_len(.N)]
fr_lab <- copy(Tx)[, row := seq_len(.N)][, label := gsub("\n", " | ", label, fixed = TRUE)]
fr_forest <- copy(long)[, .(pair = rep(cb$lab, 2), ancestry = anc, est = x, ci_lo = xmin, ci_hi = xmax, ea_coloc_robust = robust, row_y_mm = y,
                            marker_shape = shape, marker_fill = fill, ea_outcome_accession = c(rep(NA_character_, nrow(cb)), cb$eas_outcome_accession),
                            instrument = c(cb$eur_rs, cb$eas_rs))]
fr_text <- rbind(lab_df[, .(kind, x, y, label)], note_df[, .(kind = "note", x, y, label)],
                 data.table(kind = "forest annotation", x = mean(XL), y = Y_F5, label = f5_lines))
fr_text[, label := gsub("\n", " | ", label, fixed = TRUE)]
FRAMES <- list(`5a_boxes` = fr_boxes, `5a_segments` = fr_segs, `5a_labels` = fr_lab, `5b_forest` = fr_forest, `5b_text` = fr_text)
stopifnot(!anyDuplicated(names(FRAMES)), all(vapply(FRAMES, function(d) is.data.frame(d) && nrow(d) > 0, TRUE)), nrow(fr_boxes) == nrow(Bx), nrow(fr_lab) == nrow(Tx))
dir.create(file.path(P31R_FIGS, "af3", "source_data"), showWarnings = FALSE)
for (.fn in names(FRAMES)) data.table::fwrite(FRAMES[[.fn]], file.path(P31R_FIGS, "af3", "source_data", sprintf("Figure_5__%s.csv", .fn)))
p31r_save(fig, "Figure_5", w_mm = 170, h_mm = H_MM)

# registry, written from what was printed
regf <- file.path(P31R_FIGS, "registry", "Figure_5.json")
dir.create(dirname(regf), showWarnings = FALSE)
writeLines(jsonlite::toJSON(list(panels = c("a", "b"), ledger = REG$ledger, derived = list(), literals = REG$literals, frames = I(names(FRAMES))),
                            auto_unbox = TRUE, pretty = TRUE), regf, useBytes = TRUE)
message("registry: ", length(REG$ledger), " ledger entries, ", length(REG$literals), " literals -> ", regf)
