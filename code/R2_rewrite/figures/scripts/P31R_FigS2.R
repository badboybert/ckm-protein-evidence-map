# Paper 3.1 rewrite (R2) -- Additional file 3, Figure S2: the 16 multi-signal candidates, coloc.abf and multi-signal PP.H4, with the Rule 3(i)
# LD-admissibility status of each pair's fit (moved here from the old Figure 3b; S2-D24). Illustrates SM sections S-R1 and S-R2. ONE panel.
#   tick  = coloc.abf PP.H4; dot = multi-signal (SuSiE-based, purity 0.5) PP.H4; dashed line = the colocalization rule.
#   admissibility (visual-status rule 3, "inadmissible fit"): a filled dot is an admissible fit; a HATCHED circle (three diagonal chords; the crossed
#   circle is reserved for "no colocalization", S2-D25 visual check) plus the word "inadmissible" is an inadmissible one (read per pair from the fit-set readings file, PRIMARY reading: seven fits per locus, build-shifted rows kept, the pair whose
#   protein s cannot be computed reclassified). The kidney / other split is the facet; no fill carries it.
# Sources: ONLY through the approved manifest (AF3_S2): the candidate edge audit and the fit-set readings.
# Run:  Rscript "<archive>/code/R2_rewrite/figures/scripts/P31R_FigS2.R"

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
suppressPackageStartupMessages(library(ggtext))
UNIT <- "FigureS2"
af3_begin(UNIT)
EN <- "\u2013"
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
KIDNEY <- c("eGFR", "CKD")
INK <- "grey15"; LINE_GREY <- "grey45"

aud <- af3_read("V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv", "", encoding = "UTF-8")
fs  <- af3_read("R1_stage1/S1R1_R2_10_integration/s6d_fitset_readings.tsv", "", sep = "\t", quote = "\"")
stopifnot(ncol(aud) >= 11)
setnames(aud, 1:10, c("gene", "outcome", "mr_dir", "abf_h4", "abf_h3", "susie_naive", "susie_h4", "converged", "fdr", "status"))
cand <- aud[status == "Drawn"]
cand[, `:=`(abf_h4 = as.numeric(abf_h4), abf_h3 = as.numeric(abf_h3), susie_h4 = as.numeric(susie_h4))]
stopifnot(nrow(cand) == num("L222", 1), !anyNA(cand$abf_h4), !anyNA(cand$susie_h4), !anyDuplicated(cand[, paste(gene, outcome)]),
          sum(cand$outcome %in% KIDNEY) == num("L225", 1), all(cand$susie_h4 >= num("L236")), all(cand$converged == "yes"), all(cand$fdr == "yes"),
          !any(cand$abf_h4 >= num("L236") & cand$abf_h4 >= num("L237") * cand$abf_h3))   # every one fails coloc.abf
# UMOD-eGFR and DPEP1-eGFR: ledger L207 = (SuSiE UMOD, SuSiE DPEP1, abf UMOD, abf DPEP1)
.u <- cand[gene == "UMOD" & outcome == "eGFR"]; .d <- cand[gene == "DPEP1" & outcome == "eGFR"]
stopifnot(abs(.u$susie_h4 - num("L207", 1)) < 1e-4, abs(.d$susie_h4 - num("L207", 2)) < 1e-4,
          abs(.u$abf_h4 / num("L207", 3) - 1) < 1e-3, abs(.d$abf_h4 / num("L207", 4) - 1) < 1e-3)
# exactly one candidate reaches 0.80 under coloc.abf (GUSB-eGFR) and it fails the PP.H4 / PP.H3 ratio instead
stopifnot(sum(cand$abf_h4 >= num("L236")) == 1L, cand[abf_h4 >= num("L236"), abf_h4 < num("L237") * abf_h3],
          identical(cand[abf_h4 >= num("L236"), paste(gene, outcome)], "GUSB eGFR"))

# ---- admissibility: the PRIMARY reading (seven fits per locus, build-shifted rows kept, non-computable ING1-CAD reclassified) ----------------
prim <- fs[fit_set == "seven" & build_shift_rows == "keep" & startsWith(ING1_version, "(d) estimate")]
stopifnot(nrow(prim) == 1L, prim$n_reclassified == num("L886"), abs(prim$NOTCH2_s_max - num("L866")) < 5e-4)
reclass <- strsplit(prim$reclassified, ";", fixed = TRUE)[[1]]
stopifnot(length(reclass) == num("L886"), !anyDuplicated(reclass))
cand[, key := paste0(gene, "_", outcome)]
stopifnot(all(reclass %in% cand$key))
cand[, inadm := key %in% reclass]
adm <- cand[inadm == FALSE, key]
stopifnot(setequal(adm, c("GUSB_eGFR", "NT5C3A_eGFR")), sum(cand$inadm) == num("L886"))          # the Results' statement, derived not typed
# NOTCH2-T2D is inadmissible only because s (primary reading) lies just above the limit; the other readings put it below
oth <- fs[!(fit_set == "seven" & build_shift_rows == "keep")]
stopifnot(all(oth[fit_set == "four" | build_shift_rows == "drop"]$NOTCH2_reclassified %in% c("FALSE", FALSE)), prim$NOTCH2_reclassified %in% c("TRUE", TRUE))
# ING1-CAD: s is not computable (protein window unavailable) and its outcome s (L865) exceeds the limit: inadmissible under every reading
stopifnot(all(grepl("ING1_CAD", fs[startsWith(ING1_version, "(d) estimate")]$reclassified)))
cand[, note := fifelse(key == "NOTCH2_T2D", "S", fifelse(key == "ING1_CAD", "I", ""))]

cand[, grp := factor(fifelse(outcome %in% KIDNEY, "Kidney", "Other"), levels = c("Kidney", "Other"))]
cand[, pair := paste0(gene, EN, outcome)]
setorder(cand, grp, abf_h4)
cand[, pair := factor(pair, levels = unique(pair))]
stopifnot(uniqueN(cand$pair) == nrow(cand))
cand[, state := factor(fifelse(inadm, "inadmissible", "admissible"), levels = c("admissible", "inadmissible"))]
# the word column right of the axis (visual-status rule 3: the word label travels with the mark)
S_TXT <- led("L866", NULL, function(v) formatC(v, format = "f", digits = 3), where = "S2 NOTCH2-T2D label: Rule 3(i) proxy s",
             anchor = "^inadmissible \\(s {}\\)$")
cand[, wlab := fifelse(key == "NOTCH2_T2D", paste0("inadmissible (s ", S_TXT, ")"),
                fifelse(key == "ING1_CAD", "inadmissible (s not computable)", fifelse(inadm, "inadmissible", "admissible")))]

X_B <- c(0, 0.5, 1); XLIM <- c(-0.04, 1.55)
XLAB <- ticks("S2_x_pph4", X_B, XLIM, "linear", "S2 x, PP.H4")
pb <- ggplot(cand, aes(y = pair)) +
  geom_vline(xintercept = num("L236"), linetype = "dashed", colour = LINE_GREY, linewidth = 0.3) +
  geom_segment(aes(x = abf_h4, xend = susie_h4, yend = pair), colour = "grey70", linewidth = 0.4) +
  geom_point(data = cand[inadm == FALSE], aes(x = susie_h4), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.35) +
  geom_hatchdot(data = cand[inadm == TRUE], aes(x = susie_h4), ink = INK) +
  # coloc.abf: a vertical tick (white halo, then ink), drawn over the dot so a near-coincident pair (GUSB-eGFR) shows both
  geom_errorbar(aes(xmin = abf_h4, xmax = abf_h4), width = 0.78, orientation = "y", colour = "white", linewidth = 0.9) +
  geom_errorbar(aes(xmin = abf_h4, xmax = abf_h4), width = 0.72, orientation = "y", colour = INK, linewidth = 0.45) +
  geom_text(aes(x = 1.07, label = wlab), hjust = 0, size = ANNOT, family = BASE_FONT) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = X_B, labels = XLAB, expand = expansion(0), guide = guide_axis(cap = "upper"),
                     sec.axis = dup_axis(breaks = num("L236"), labels = "Colocalization rule", name = NULL)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) +
  coord_cartesian(xlim = XLIM) +
  labs(x = "PP.H4 (tick, coloc.abf; dot, multi-signal)", y = NULL) +
  theme_p31r() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), strip.clip = "off", panel.spacing.y = unit(4, "pt"),
        axis.text.x.top = element_text(size = ANNOT_PT, hjust = 1, margin = margin(b = 1)), axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())

# header (own strip) and key (own strip), 7 pt
N_INADM <- led("L886", NULL, function(v) formatC(v, format = "f", digits = 0), where = "S2 header: candidates on inadmissible fits",
               anchor = "kept\\. {} of 16 candidates")
N_ALL <- led("L222", 1, function(v) formatC(v, format = "f", digits = 0), where = "S2 header: candidates", anchor = "kept\\. 14 of {} candidates")
ph <- ggplot() +
  annotate("text", x = 0, y = 2, hjust = 0, vjust = 1, size = ANNOT, family = BASE_FONT, lineheight = 1.05,
           label = paste0("Primary reading: seven fits per locus, build-shifted rows kept. ", N_INADM, " of ", N_ALL,
                          " candidates rest on inadmissible fits;\nonly GUSB", EN, "eGFR and NT5C3A", EN, "eGFR are admissible. Dot, multi-signal fit; tick, coloc.abf.")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2.2), expand = FALSE) + theme_void()
pk <- ggplot() +
  geom_point(aes(x = 0.000, y = 2), shape = 21, fill = INK, colour = INK, size = 1.9, stroke = 0.35) +
  geom_text(aes(x = 0.016, y = 2, label = "Admissible fit (LD-admissibility rule met)"), hjust = 0, size = ANNOT, family = BASE_FONT) +
  geom_hatchdot(aes(x = 0.000, y = 1), data = data.frame(x = 0, y = 1), ink = INK) +
  geom_text(aes(x = 0.016, y = 1, label = "Inadmissible fit (rule not met), hatched"), hjust = 0, size = ANNOT, family = BASE_FONT) +
  coord_cartesian(xlim = c(-0.01, 1), ylim = c(0.4, 2.6), expand = FALSE) + theme_void()

fig <- ph / pb / pk + plot_layout(heights = c(0.10, 1, 0.09))
af3_save(fig, h_mm = 125, panels = character(0), panel_map = character(0),
         frames = list(candidates_16 = cand[, .(gene, outcome, grp, abf_h4, abf_h3, susie_h4, inadm, state, note)]))
