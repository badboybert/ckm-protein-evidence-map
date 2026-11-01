# Paper 3.1 rewrite (R2) -- Figure 1: design and map overview (S2-D24 redesign, plan FIGURE_MESSAGES_R2 v3 section 1).
# One panel, 170 mm wide. Run:  Rscript "V3/R2_rewrite/figures/scripts/P31R_Fig1.R"   (any working directory)
#
# What changed against the R1/3a render (plan section 0 and 1, decisions F4, F6, F8; S2-D24):
#   REMOVED: the 84/10 protein-metabolite row, the ARG1 divergence line ("arginine follows .., T2D follows .."), the renal
#   Fisher p values and "as at UMOD", the "6 of 6 concordant" line and the post hoc 22/23 line, the 248/177/160 pool chain.
#   ADDED:   the corroboration row (L985, three pairs), the FGFR4-T2D provenance box, tier definitions inside the tier boxes
#   with the Convergent names beside them (F6; names are in the body, never in a title), one compact key inside the graphic.
#
# SOURCES. The approved source manifest (registry/APPROVED_SOURCE_MANIFEST.tsv) lists exactly two sources for this figure:
#   (1) the claims ledger (printed numbers, read with p31r_ledger at the point of printing, pinned by id); and
#   (2) the tiered pair table, copy B (md5 31576d12...), for the tier counts, the significance classes and the Convergent names.
# Nothing else is opened. The R1 script re-derived ~40 counts from raw V4 files that the manifest does not list (G12 would
# fail them), so those counts are now LEDGER-ONLY (led(..., lo = TRUE), PIPE status); the counts the tiered table can
# re-derive (81, 62, 3/45/33, 41, 80, DPEP1-AF, Convergent names, the eight outcomes) are still checked against the ledger.
# Every printed NAME that states a result (Convergent pairs, the pair that is not FDR-significant, the corroborated pairs and
# the one among them outside the 81) is built from data or the ledger and re-read from the DRAWN labels by the label audit.
# Layout is in millimetres on a 170 x H mm canvas (1 data unit = 1 mm); every text line is checked to fit its box.

.args <- commandArgs(trailingOnly = FALSE)
.me <- normalizePath(sub("^--file=", "", .args[grepl("^--file=", .args)]), winslash = "/")
source(file.path(dirname(.me), "P31R_fig_common.R"))
suppressPackageStartupMessages(library(systemfonts))
FIG <- "Figure_1"

# ---------------------------------------------------------------------------------------------------------------
# 0. Ledger access that registers what it prints
# ---------------------------------------------------------------------------------------------------------------
REG <- new.env(); REG$e <- list()
# led(): fetch ONE ledger element, format it, record {text, id, idx, anchor, pct} for the registry, return the text.
# Call it only where the returned string is drawn (each call = one drawn occurrence; counts are tallied from the calls).
# EVERY call carries an anchor (a regex with {} for the number, matched against the drawn line), so the registry is
# positional: the gate's G7 then fails a number printed in another number's place, not only a wrong value.
# A printed COUNT passes `check` (its value re-derived from the tiered table) or says lo = TRUE (LEDGER-ONLY: no approved
# file re-derives it for this figure). Only the analysis thresholds read from the ledger (THRESH) print without either.
THRESH <- c("L236", "L237", "L238", "L239")
led <- function(id, idx = NULL, anchor, fmt = p31r_num, pct = FALSE, check = NULL, lo = FALSE) {
  if (missing(anchor) || !grepl("{}", anchor, fixed = TRUE)) stop("led(", id, "): an anchor with {} is required")
  v <- p31r_ledger(id, idx)
  if (is.null(check) && !lo && !(id %in% THRESH)) stop("led(", id, "): a printed count needs check = <its re-derived value> or lo = TRUE")
  if (!is.null(check) && !isTRUE(all.equal(as.numeric(v), as.numeric(check)))) stop(sprintf("led(%s, %s): ledger %s != re-derived %s", id, format(idx), format(v), format(check)))
  s <- fmt(v)
  e <- list(text = s, id = id)
  if (!is.null(idx)) e$idx <- idx
  e$anchor <- anchor
  if (pct) e$pct <- TRUE
  REG$e[[length(REG$e) + 1]] <- e
  s
}
f2 <- function(v) sprintf("%.2f", v)                                   # thresholds printed to two decimals
fpct <- function(v) sprintf("%d%%", as.integer(round(100 * v)))         # a ledger fraction printed as a percentage
EN <- "–"; GE <- "≥"; TIMES <- "×"
`%||%` <- function(a, b) if (is.null(a)) b else a
en <- function(x) gsub("-", EN, x, fixed = TRUE)                        # PROTEIN-OUTCOME -> PROTEIN–OUTCOME

# ---------------------------------------------------------------------------------------------------------------
# 1. Re-derive what the tiered table can derive and assert it equals the ledger
# ---------------------------------------------------------------------------------------------------------------
F_TIER <- file.path(V3_ROOT, "v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv")
ti <- p31r_read(F_TIER)                                                  # canonical copy B (SOURCE_RECONCILIATION_tiered_table.md)
stopifnot(!anyDuplicated(ti[, .(gene, outcome)]))
ti_pairs <- paste(ti$gene, ti$outcome, sep = "-")
tiers <- table(ti$evidence_tier); sig <- table(ti$mr_sig_class)
CONV <- ti[evidence_tier == "Convergent"][order(gene)]
CONV_PAIRS <- paste(CONV$gene, CONV$outcome, sep = "-")
NOTFDR <- ti[mr_sig_class == "Suggestive"]
chk <- list(
  L010 = c(nrow(ti), p31r_ledger("L010")), L011 = c(uniqueN(ti$gene), p31r_ledger("L011")), L012 = c(sig[["Bonferroni"]], p31r_ledger("L012")),
  L014c = c(tiers[["Convergent"]], p31r_ledger("L014", "Convergent")), L014s = c(tiers[["Supported"]], p31r_ledger("L014", "Supported")),
  L014f = c(tiers[["Flagged"]], p31r_ledger("L014", "Flagged")),
  L959b = c(sig[["Bonferroni"]], p31r_ledger("L959", 1L)), L959f = c(sig[["Bonferroni"]] + sig[["FDR"]], p31r_ledger("L959", 5L)),
  L959n = c(nrow(ti), p31r_ledger("L959", 6L)),
  L002 = c(uniqueN(ti$outcome), p31r_ledger("L002")))                    # the 81 span all eight outcomes (Results map)
bad <- names(chk)[vapply(chk, function(v) !isTRUE(all.equal(as.numeric(v[[1]]), as.numeric(v[[2]]))), logical(1))]
if (length(bad)) stop("re-derived count != ledger: ", paste(sprintf("%s (%s vs %s)", bad, sapply(chk[bad], `[[`, 1), sapply(chk[bad], `[[`, 2)), collapse = "; "))
stopifnot(identical(p31r_ledger("L959", 4L), paste(NOTFDR$gene, NOTFDR$outcome, sep = "-")), nrow(NOTFDR) == 1,
          setequal(unlist(p31r_ledger("L015")), CONV_PAIRS))               # the Convergent names printed beside the tier definition
# the corroborated pairs (L985): all three are pairs selected in advance; two are among the 81, the third (FGFR4-T2D) is not
CORR <- unlist(p31r_ledger("L985", "pairs")); stopifnot(length(CORR) == p31r_ledger("L985", "n"), !anyDuplicated(CORR))
OFF_PAIR <- CORR[!CORR %in% ti_pairs]; stopifnot(length(OFF_PAIR) == 1, OFF_PAIR == "FGFR4-T2D")
message("Figure 1: ", length(chk) + 4, " tier-table checks equal the ledger")

# ---------------------------------------------------------------------------------------------------------------
# 2. Layout (mm from the TOP; converted to plot y = H - y at draw time)
# ---------------------------------------------------------------------------------------------------------------
PITCH <- 3.0; PADV <- 1.2; PADH <- 1.5
X0 <- 26; X1 <- 110; XC <- (X0 + X1) / 2           # main column
XR <- 163.5; XLANE <- 167.5                        # right edge of the right-hand boxes and rows; lane for the candidate connector
COL_SET <- c(fill = "white", border = "grey30")     # a set of combinations / associations / pairs (colour is kept for the tiers)
TIP <- 0.3                                          # a closed mitred head overshoots its end point by ~0.35 mm: stop short of a border
boxes <- list(); texts <- list(); segs <- list()
tw <- function(s, face) {                            # string width in mm (systemfonts, Arial, 7 pt)
  string_width(s, family = BASE_FONT, size = ANNOT_PT, res = 72, bold = face %in% c("bold", "bold.italic"),
               italic = face %in% c("italic", "bold.italic")) / 72 * 25.4
}
hbox <- function(n) n * PITCH + 2 * PADV
add_text <- function(x, ytop, lab, face = "plain", hjust = 0.5, col = "grey10", box = NULL) {
  texts[[length(texts) + 1]] <<- data.frame(x = x, y = ytop, label = lab, face = face, hjust = hjust, col = col, box = if (is.null(box)) NA else box)
}
add_box <- function(id, x0, x1, y0, lines, faces = c("bold", rep("plain", length(lines) - 1)), fill = COL_SET[["fill"]],
                    border = COL_SET[["border"]], lty = "solid", h = hbox(length(lines)), hjust = 0.5, lwd = 0.4, top = FALSE) {
  boxes[[id]] <<- data.frame(id = id, x0 = x0, x1 = x1, y0 = y0, y1 = y0 + h, fill = fill, border = border, lty = lty, lwd = lwd)
  w <- tw(lines, faces)
  if (any(w + 2 * PADH > x1 - x0 + 1e-9)) stop(sprintf("box %s: line '%s' needs %.1f mm, box is %.1f mm", id, lines[which.max(w)], max(w) + 2 * PADH, x1 - x0))
  y_first <- if (top) y0 + PADV + PITCH / 2 else y0 + (h - length(lines) * PITCH) / 2 + PITCH / 2
  x <- if (hjust == 0.5) (x0 + x1) / 2 else x0 + PADH
  for (i in seq_along(lines)) add_text(x, y_first + (i - 1) * PITCH, lines[i], faces[i], hjust, box = id)
  invisible(boxes[[id]])
}
arrow_v <- function(x, ya, yb, lty = "solid") segs[[length(segs) + 1]] <<- data.frame(x = x, y = ya, xend = x, yend = yb, lty = lty, arrow = TRUE)
seg <- function(x, y, xend, yend, lty = "solid", arrow = TRUE) segs[[length(segs) + 1]] <<- data.frame(x = x, y = y, xend = xend, yend = yend, lty = lty, arrow = arrow)

# --- inputs (two sources joined into the test family)
OUT_FG <- c("CAD", "HF", "AF", "T2D", "NAFLD"); OUT_MS <- "stroke"; OUT_CK <- c("eGFR", "CKD")   # Methods, outcomes sentence
stopifnot(setequal(c(OUT_FG, "Stroke", OUT_CK), unique(ti$outcome)))                           # the typed names are the 81's eight outcomes
s_prot1 <- led("L001", anchor = "^{} proteins with a cis-pQTL instrument$", lo = TRUE)
y <- 2
add_box("prot", X0, 96, y, c(sprintf("%s proteins with a cis-pQTL instrument", s_prot1), "UK Biobank Pharma Proteomics Project", "Olink Explore, European ancestry"))
add_box("out", 102, XR, y, c(sprintf("%s UK Biobank-free outcome GWAS", led("L002", anchor = "^{} UK Biobank-free outcome GWAS$", check = uniqueN(ti$outcome))),
                             sprintf("FinnGen: %s", paste(OUT_FG, collapse = ", ")),
                             sprintf("MEGASTROKE: %s; CKDGen: %s", OUT_MS, paste(OUT_CK, collapse = ", "))))
y_in_end <- y + hbox(3)
y_join <- y_in_end + 2.5; y <- y_in_end + 6
seg((X0 + 96) / 2, y_in_end, (X0 + 96) / 2, y_join, arrow = FALSE); seg((102 + XR) / 2, y_in_end, (102 + XR) / 2, y_join, arrow = FALSE)
seg((X0 + 96) / 2, y_join, (102 + XR) / 2, y_join, arrow = FALSE); arrow_v(XC, y_join, y)
# --- screen
y_fam <- y
add_box("fam", X0, X1, y, c(sprintf("%s protein%soutcome combinations", led("L003", anchor = "^{} protein–outcome combinations$", lo = TRUE), EN),
                            sprintf("%s proteins %s %s outcomes: the multiple-testing family", led("L001", anchor = "^{} proteins × [0-9]+ outcomes", lo = TRUE), TIMES, led("L002", anchor = "proteins × {} outcomes: the multiple-testing family$", check = uniqueN(ti$outcome)))))
y_scr0 <- y; y <- y + hbox(2); arrow_v(XC, y, y + 4); y <- y + 4
add_box("test", X0, X1, y, c(sprintf("%s testable", led("L004", anchor = "^{} testable$", lo = TRUE)), "instrument available in the outcome GWAS"))
y <- y + hbox(2); arrow_v(XC, y, y + 5)
add_text(XC + 2, y + 2.5, "two-sample Mendelian randomization", "italic", 0, "grey25")
y <- y + 5
add_box("nom", X0, X1, y, c(sprintf("%s nominal associations (p < %s)", led("L006", anchor = "^{} nominal associations", lo = TRUE), led("L238", anchor = "nominal associations \\(p < {}\\)$", fmt = f2)),
                            sprintf("%s study-wide (Bonferroni); %s at %s FDR", led("L007", anchor = "^{} study-wide \\(Bonferroni\\); ", lo = TRUE), led("L008", anchor = "\\(Bonferroni\\); {} at ", lo = TRUE),
                                    led("L239", anchor = "\\(Bonferroni\\); [0-9]+ at {} FDR$", fmt = fpct, pct = TRUE))))
y_nom <- y; y <- y + hbox(2); y_scr1 <- y
# --- colocalization step, with the multi-signal candidates beside it
Y_ARR <- 10
arrow_v(XC, y, y + Y_ARR)
add_text(XC + 2, y + Y_ARR / 2 - PITCH / 2, "Colocalization (coloc.abf):", "italic", 0, "grey25")
add_text(XC + 2, y + Y_ARR / 2 + PITCH / 2, sprintf("PP.H4 %s %s and PP.H4/PP.H3 %s %s", GE, led("L236", anchor = "^PP.H4 ≥ {} and PP.H4/PP.H3", fmt = f2), GE, led("L237", anchor = "and PP.H4/PP.H3 ≥ {}$")), "italic", 0, "grey25")
# the 16 multi-signal candidates: outside the 81, entered by the kidney row (dashed connector, right-hand lane)
cand <- c(sprintf("%s multi-signal candidates", led("L222", 1L, anchor = "^{} multi-signal candidates$", lo = TRUE)),
          "FDR-significant pairs that failed",
          sprintf("coloc.abf but reached PP.H4 %s %s", GE, led("L236", anchor = "coloc.abf but reached PP.H4 ≥ {}$", fmt = f2)),
          "in multi-signal fine-mapping;",
          sprintf("outside the %s graded pairs", led("L010", anchor = "^outside the {} graded pairs$", check = nrow(ti))))
X_CAND <- 114
add_box("cand", X_CAND, XR, y_nom - 1, cand, faces = c("bold", rep("plain", 4)), fill = "white", border = "grey40", lty = "22")
y_fdrline <- y_nom + (hbox(2) - 2 * PITCH) / 2 + PITCH / 2 + PITCH                       # centre of the nominal box's FDR line
seg(X1, y_fdrline, X_CAND - TIP, y_fdrline, lty = "22")                                  # the candidates come from (a subset of) the FDR set
y_cand_mid <- y_nom - 1 + hbox(length(cand)) / 2
# --- key, inside the graphic (S2-D24 visual status): two box styles and the italic convention
X_KEY <- 116; y_key <- y_fam + 1
boxes[["key_s"]] <- data.frame(id = "key_s", x0 = X_KEY, x1 = X_KEY + 5, y0 = y_key, y1 = y_key + 2.4, fill = "white", border = COL_SET[["border"]], lty = "solid", lwd = 0.4)
boxes[["key_d"]] <- data.frame(id = "key_d", x0 = X_KEY, x1 = X_KEY + 5, y0 = y_key + 4.2, y1 = y_key + 6.6, fill = "white", border = "grey40", lty = "22", lwd = 0.4)
add_text(X_KEY + 7, y_key + 1.2, "Solid: a counted set", "plain", 0, "grey10")
add_text(X_KEY + 7, y_key + 5.4, sprintf("Dashed: outside the %s graded pairs", led("L010", anchor = "^Dashed: outside the {} graded pairs$", check = nrow(ti))), "plain", 0, "grey10")
add_text(X_KEY, y_key + 9.6, "Italic: the method or rule on an arrow", "italic", 0, "grey25")
y <- y + Y_ARR
add_box("pairs", X0, X1, y, c(sprintf("%s colocalization-supported pairs over %s proteins", led("L010", anchor = "^{} colocalization-supported pairs over", check = nrow(ti)), led("L011", anchor = "colocalization-supported pairs over {} proteins$", check = uniqueN(ti$gene))),
                              sprintf("%s study-wide; %s at %s FDR (all but %s)", led("L959", 1L, anchor = "^{} study-wide; ", check = sig[["Bonferroni"]]), led("L959", 5L, anchor = "study-wide; {} at ", check = sig[["Bonferroni"]] + sig[["FDR"]]),
                                      led("L239", anchor = "study-wide; [0-9]+ at {} FDR \\(all but", fmt = fpct, pct = TRUE), en(p31r_ledger("L959", 4L)))))
y_col0 <- y_scr1; y <- y + hbox(2); y_col1 <- y
# --- grading: tier boxes carry the definitions; the Convergent names sit beside the Convergent definition (F6)
Y_GR <- Y_ARR + 1; y_bar <- y + Y_GR - 2.5
seg(XC, y, XC, y_bar, arrow = FALSE)                                                     # from the pairs box to the bar over the tiers
add_text(XC + 2, y + (Y_GR - 2.5) / 2 - PITCH / 2, "Graded on internal support and", "italic", 0, "grey25")
add_text(XC + 2, y + (Y_GR - 2.5) / 2 + PITCH / 2, "exposure-side corroboration (deCODE, ARIC)", "italic", 0, "grey25")
add_text(XR - PADH, y + (Y_GR - 2.5) / 2, "Figure 2a", "bold", 1, "grey10")
y <- y + Y_GR + 0.5
tint <- function(col, f) { m <- grDevices::col2rgb(col) / 255; grDevices::rgb(t(1 - (1 - m) * f)) }
TIERS <- c("Convergent", "Supported", "Flagged")
TW <- c(Convergent = 50, Supported = 38, Flagged = 38); GAPT <- 2
TXA <- X0 + c(0, cumsum(TW + GAPT))[1:3]; TXC <- TXA + TW / 2
tier_lines <- list(
  Convergent = c("Every internal criterion met; an", "exposure-side arm in the discovery", "direction; no weakness", paste(en(CONV_PAIRS), collapse = ", ")),
  Supported  = c("All remaining pairs:", "no weakness found"),
  Flagged    = c("At least one weakness in", "colocalization, corroboration", "or direction"))
H_TIER <- hbox(max(lengths(tier_lines)) + 1)
for (i in seq_along(TIERS)) {
  t_ <- TIERS[i]
  add_box(paste0("tier_", t_), TXA[i], TXA[i] + TW[[t_]], y, c(sprintf("%s %s", led("L014", t_, anchor = paste0("^{} ", t_, "$"), check = tiers[[t_]]), t_), tier_lines[[t_]]),
          faces = c("bold", rep("plain", length(tier_lines[[t_]]))), fill = tint(PAL_TIER[[t_]], 0.22), border = PAL_TIER[[t_]], lwd = 0.8, h = H_TIER, top = TRUE)
}
seg(TXC[1], y_bar, TXC[3], y_bar, arrow = FALSE)
for (xc in TXC) arrow_v(xc, y_bar, y + 0.0)
y_gr0 <- y_col1; y <- y + H_TIER; y_gr1 <- y
# --- follow-up analyses (each names the figure that shows it; each has its own starting set)
for (xc in TXC) seg(xc, y, xc, y + 2, arrow = FALSE)
seg(TXC[1], y + 2, TXC[3], y + 2, arrow = FALSE)
# FGFR4-T2D, one of the pairs selected in advance, is not one of the 81: it enters the external-check row from its own dashed box
X_OFF <- 112; Y_OFF <- 4.5; Y_OFF_ARR <- 4.5
off <- c(en(OFF_PAIR),
         "identified only in this study's East",
         sprintf("Asian analysis; outside the %s;", led("L010", anchor = "^Asian analysis; outside the {};$", check = nrow(ti))),
         "previously reported for T2D (Table 2)")                                          # S2-D25: the FGFR4-T2D association itself (plan v3 s1 note)
add_box("off", X_OFF, XR, y + Y_OFF, off, faces = c("bold", rep("plain", 3)), fill = "white", border = "grey40", lty = "22")
GAP_FU <- Y_OFF + hbox(length(off)) + Y_OFF_ARR
seg((X_OFF + XR) / 2, y + Y_OFF + hbox(length(off)), (X_OFF + XR) / 2, y + GAP_FU - TIP, lty = "22")
arrow_v(XC, y + 2, y + GAP_FU); y <- y + GAP_FU
y_fu0 <- y
corr_names <- paste0(paste(en(CORR[CORR != OFF_PAIR]), collapse = ", "), ", ", en(OFF_PAIR), " (in East Asians)")
FU <- list(
  list("External checks",
       c(sprintf("%s pairs corroborated under the rule (Methods):", led("L985", "n", anchor = "^{} pairs corroborated under", lo = TRUE)), corr_names), "Figure 2b"),
  list("Kidney",
       c(sprintf("%s of the %s multi-signal candidates are kidney pairs (descriptive);", led("L225", 1L, anchor = "^{} of the [0-9]+ multi-signal candidates", lo = TRUE), led("L222", 1L, anchor = "^[0-9]+ of the {} multi-signal candidates are", lo = TRUE)),
         "UMOD as a control comparison"), "Figure 3"),
  list("ARG1 locus", c("Two signals: arginine, ornithine and T2D", "compared separately"), "Figure 4"),
  list("East Asian",
       c(sprintf("%s of %s evaluable pairs kept the European sign (European instrument);", led("L224", 1L, anchor = "^{} of [0-9]+ evaluable pairs kept", lo = TRUE), led("L224", 2L, anchor = "^[0-9]+ of {} evaluable pairs kept", lo = TRUE)),
         sprintf("%s of %s testable pairs colocalized (East Asian instrument)", led("L082", 3L, anchor = "^{} of [0-9]+ testable pairs colocalized", lo = TRUE), led("L082", 2L, anchor = "^[0-9]+ of {} testable pairs colocalized", lo = TRUE))), "Figure 5"))
XN <- X0 + PADH; XT <- 60; FU_GAP <- 1.2; KIDNEY_ROW <- 2
XF <- XR - PADH                                                                          # right-aligned figure pointers
for (k in seq_along(FU)) {
  f <- FU[[k]]; id <- paste0("fu", k); fh <- hbox(max(length(f[[1]]), length(f[[2]])))
  boxes[[id]] <- data.frame(id = id, x0 = X0, x1 = XR, y0 = y, y1 = y + fh, fill = "white", border = COL_SET[["border"]], lty = "solid", lwd = 0.4)   # not the Flagged grey
  if (any(tw(f[[1]], "bold") > XT - XN - 1.5)) stop("follow-up name too wide: ", f[[1]][1])
  if (any(tw(f[[2]], "plain") > XF - tw(f[[3]], "bold") - 2 - XT)) stop("follow-up text too wide: ", f[[2]][which.max(tw(f[[2]], "plain"))], " (", round(max(tw(f[[2]], "plain")), 1), " mm > ", round(XF - tw(f[[3]], "bold") - 2 - XT, 1), ")")
  for (i in seq_along(f[[1]])) add_text(XN, y + fh / 2 + (i - (length(f[[1]]) + 1) / 2) * PITCH, f[[1]][i], "bold", 0, box = id)
  for (i in seq_along(f[[2]])) add_text(XT, y + fh / 2 + (i - (length(f[[2]]) + 1) / 2) * PITCH, f[[2]][i], "plain", 0, box = id)
  add_text(XF, y + fh / 2, f[[3]], "bold", 1, box = id)
  if (k == KIDNEY_ROW) y_kid <- y + fh / 2
  y <- y + fh + if (k < length(FU)) FU_GAP else 0
}
y_fu1 <- y
# the multi-signal candidates are reported by the kidney row: dashed connector down the right-hand lane
seg(XR, y_cand_mid, XLANE, y_cand_mid, lty = "22", arrow = FALSE)
seg(XLANE, y_cand_mid, XLANE, y_kid, lty = "22", arrow = FALSE)
seg(XLANE, y_kid, XR + TIP, y_kid, lty = "22")
H <- ceiling(y_fu1 + 2)                                                                  # page height follows the layout
stopifnot(H <= FIG_HMAX_MM)
# --- stage labels down the left margin, each with a bracket spanning its rows
stages <- data.frame(label = c("Inputs", "Screen", "Colocalization", "Grading", "Follow-up"),
                     y0 = c(2, y_scr0, y_col0, y_gr0, y_fu0), y1 = c(y_in_end, y_scr1, y_col1, y_gr1, y_fu1))
stopifnot(max(tw(stages$label, "bold")) < 21.5 - 1.5)

# ---------------------------------------------------------------------------------------------------------------
# 3. Label audit: re-read the DRAWN labels and check every result-bearing word against the data
# ---------------------------------------------------------------------------------------------------------------
TX <- do.call(rbind, texts)
LAB <- TX$label
one <- function(rx) {
  h <- regmatches(LAB, regexec(rx, LAB, perl = TRUE)); h <- h[lengths(h) > 0]
  if (length(h) != 1) stop("label audit: /", rx, "/ matched ", length(h), " drawn lines (want 1)")
  h[[1]][-1]
}
PR <- "([A-Z][A-Z0-9]*)–([A-Za-z0-9]+)"                                                   # a printed PROTEIN–OUTCOME pair
# the corroborated pairs: the two inside the 81, then the one outside it, tagged as East Asian
a <- one(paste0("^", PR, ", ", PR, ", ", PR, " \\(in East Asians\\)$"))
stopifnot(setequal(paste(a[c(1, 3, 5)], a[c(2, 4, 6)], sep = "-"), CORR), paste(a[5], a[6], sep = "-") == OFF_PAIR,
          all(paste(a[c(1, 3)], a[c(2, 4)], sep = "-") %in% ti_pairs))
a <- one(paste0("^", PR, "$")); stopifnot(paste(a[1], a[2], sep = "-") == OFF_PAIR)           # the FGFR4 box heading
a <- one(paste0("^", PR, ", ", PR, ", ", PR, "$"))                                         # the Convergent names beside the definition
stopifnot(setequal(paste(a[c(1, 3, 5)], a[c(2, 4, 6)], sep = "-"), CONV_PAIRS), setequal(CONV_PAIRS, unlist(p31r_ledger("L015"))))
a <- one(paste0("\\(all but ", PR, "\\)$"))
stopifnot(paste(a[1], a[2], sep = "-") == paste(NOTFDR$gene, NOTFDR$outcome, sep = "-"))
a <- list(FinnGen = strsplit(one("^FinnGen: (.+)$"), ", ")[[1]], MEGASTROKE = one("^MEGASTROKE: ([^;]+); CKDGen: .+$"), CKDGen = strsplit(one("^MEGASTROKE: .+; CKDGen: (.+)$"), ", ")[[1]])
stopifnot(!anyDuplicated(unlist(a)), setequal(unlist(a), c(OUT_FG, OUT_MS, OUT_CK)), length(unlist(a)) == p31r_ledger("L002"))
for (t_ in TIERS) one(sprintf("^[0-9]+ %s$", t_))
# words that must not appear on this figure (plan section 1, 'Must not show or imply'; S2-D24)
bad_w <- grep("replicat|validat|enrich|over-represent|diverg|follows|Fisher|concordant|Highest|screen: |\\b84\\b|as at UMOD|robust", LAB, ignore.case = TRUE, value = TRUE)
if (length(bad_w)) stop("label audit: removed wording is back: ", paste(bad_w, collapse = " | "))

# registry: merge identical calls; count = number of calls (= drawn occurrences). Cross-checks against the DRAWN labels
# (the generator's own text, not the PDF): (i) each entry's anchor, with its number, matches exactly `count` drawn lines;
# (ii) the drawn labels hold exactly as many numbers as led() calls, once gene-like tokens (T2D, PP.H4, ARG1), figure and
# table pointers are masked, so no number reaches the page without passing through led()
key <- vapply(REG$e, function(e) paste(e$text, e$id %||% "", paste(e$idx, collapse = "/"), e$anchor, isTRUE(e$pct), sep = "|"), character(1))
REGL <- lapply(unique(key), function(k) { e <- REG$e[[match(k, key)]]; n <- sum(key == k); if (n > 1) e$count <- n; e })
for (e in REGL) {
  rx <- gsub("{}", gsub("([.()+?*^$|\\[\\]])", "\\\\\\1", e$text), e$anchor, fixed = TRUE)
  hits <- sum(lengths(regmatches(LAB, gregexpr(rx, LAB, perl = TRUE))))
  if (hits != (e$count %||% 1L)) stop(sprintf("registry/label mismatch for '%s' (%s, anchor %s): drawn %d, registered %d", e$text, e$id, e$anchor, hits, e$count %||% 1L))
}
mlab <- gsub("Tables? [0-9]+", " ", gsub("\\b[A-Za-z]+[0-9][A-Za-z0-9]*\\b", " ", LAB[!grepl("^Figures? ", LAB)], perl = TRUE))
n_drawn <- sum(lengths(regmatches(mlab, gregexpr("(?<![A-Za-z0-9.,])[0-9](?:[0-9]{0,2}(?:,[0-9]{3})+|[0-9]*)(?:\\.[0-9]+)?%?", mlab, perl = TRUE))))
if (n_drawn != length(REG$e)) stop("drawn labels hold ", n_drawn, " numbers but led() printed ", length(REG$e), ": a number was typed")
message("Figure 1: label audit passed; ", length(REGL), " registry entries (", length(REG$e), " printed numbers)")

# ---------------------------------------------------------------------------------------------------------------
# 4. Draw
# ---------------------------------------------------------------------------------------------------------------
B <- do.call(rbind, boxes); SG <- do.call(rbind, segs)
# a DASHED arrow is drawn as a dashed shaft that stops inside the head plus a short SOLID segment carrying the head:
# linetype applies to a closed head's outline too, which drew ragged dashed heads (3a re-verification regression)
HEAD_IN <- 1.0
dh <- SG$arrow & SG$lty != "solid"
if (any(dh)) {
  d <- SG[dh, ]; L_ <- sqrt((d$xend - d$x)^2 + (d$yend - d$y)^2); ux <- (d$xend - d$x) / L_; uy <- (d$yend - d$y) / L_
  stopifnot(all(L_ > 2 * HEAD_IN))
  shaft <- transform(d, xend = xend - HEAD_IN * ux, yend = yend - HEAD_IN * uy, arrow = FALSE)
  head <- transform(d, x = xend - HEAD_IN * ux, y = yend - HEAD_IN * uy, lty = "solid", arrow = TRUE)
  SG <- rbind(SG[!dh, ], shaft, head)
}
B$ya <- H - B$y1; B$yb <- H - B$y0
TX$yp <- H - TX$y
SG$yp <- H - SG$y; SG$yendp <- H - SG$yend
ST <- stages; ST$ya <- H - ST$y1; ST$yb <- H - ST$y0; ST$ym <- (ST$ya + ST$yb) / 2

p <- ggplot() +
  geom_segment(data = ST, aes(x = 22.5, xend = 22.5, y = ya + 0.6, yend = yb - 0.6), colour = "grey70", linewidth = 0.6) +
  geom_text(data = ST, aes(x = 21, y = ym, label = label), hjust = 1, size = ANNOT, family = BASE_FONT, fontface = "bold", colour = "grey35") +
  geom_rect(data = B, aes(xmin = x0, xmax = x1, ymin = ya, ymax = yb, linetype = I(lty), linewidth = I(lwd)),
            fill = B$fill, colour = B$border) +
  geom_segment(data = SG[!SG$arrow, ], aes(x = x, y = yp, xend = xend, yend = yendp, linetype = I(lty)),
               colour = "grey35", linewidth = 0.4, lineend = "square") +
  geom_segment(data = SG[SG$arrow, ], aes(x = x, y = yp, xend = xend, yend = yendp),
               colour = "grey35", linewidth = 0.4, linetype = "solid",
               arrow = arrow(length = unit(1.4, "mm"), type = "closed", angle = 25), linejoin = "mitre") +
  geom_text(data = TX, aes(x = x, y = yp, label = label, fontface = face, hjust = hjust), colour = TX$col,
            size = ANNOT, family = BASE_FONT, vjust = 0.5, lineheight = 1) +
  scale_x_continuous(limits = c(0, 170), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 170), ylim = c(0, H), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_p31r() +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.length = unit(0, "pt"), panel.grid = element_blank(), panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "none", plot.margin = margin(0, 0, 0, 0))

p31r_save(p, FIG, w_mm = FIG_W_MM, h_mm = H)

# registry, written from the same led() calls that printed the numbers (after the save, so a failed render leaves the old one)
regd <- file.path(P31R_FIGS, "registry"); dir.create(regd, showWarnings = FALSE)
writeLines(enc2utf8(as.character(jsonlite::toJSON(list(panels = list(), ledger = REGL, derived = list(), literals = list()),
                                                  auto_unbox = TRUE, pretty = TRUE))),
           file.path(regd, paste0(FIG, ".json")), useBytes = TRUE)
message("registry written: ", length(REGL), " ledger entries, 0 derived, 0 literals")
