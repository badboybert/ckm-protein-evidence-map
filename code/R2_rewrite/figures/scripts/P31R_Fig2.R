# Paper 3.1 rewrite (R2) -- Figure 2: evidence dimensions and external checks (plan FIGURE_MESSAGES_R2 v3, section 2; decision S2-D24).
#   a  evidence map of the 81 colocalization-supported pairs: conservative-prior PP.H4 (x) by exposure-side corroboration (rows);
#      fill = evidence tier; notes: 16 core-complete = 3 + 11 + 2; the LD-admissibility mark (reference-LD fit inadmissible, tier
#      not rebuilt) on TGFB1-CAD and ARG1-T2D; 79 of 81 pairs not assessed; 35 of the 45 Supported pairs have no evaluable arm
#   b  the six pairs selected in advance, as a forest with columns: discovery and external estimates (95% CI, log-odds per SD of
#      protein), the resource's own p, the regional colocalization PP.H4 at the default and conservative priors ("external regional
#      check not run" for ARG1-T2D and FGFR4-T2D only), the corroboration-rule verdict with its reason, and a SCOPED literature status
#   c  TGFB1-CAD: conditional estimates at the two instruments (rs73045269, rs12983775; 32 kb apart) in FinnGen CHD and
#      CARDIoGRAMplusC4D, with the LD and protein-altering/assay qualifications (the post hoc 28-pair scatter is AF3 S1)
# MUST NOT show or imply (plan section 2): "replicated" or concordance as corroboration; RAD23B as robust at the Bonferroni cut; a
# "not run" cell drawn like an evaluated one or left blank; PCSK9 as a failed target; the TGFB1 allele direction (HOLD H4);
# independence of DIAMANTE / CARDIoGRAMplusC4D from FinnGen.
# Every number PRINTED comes from the claims ledger through p31r_ledger(id, idx) and is registered {text, id, idx} from THAT call
# (figures/registry/Figure_2.json is written from these calls). Plotted series: ONE named source each, pinned by md5 in the render log
# and listed in registry/APPROVED_SOURCE_MANIFEST.tsv (Figure_2 rows): tiered table copy B (panel a), W11 extchecks results.tsv
# (b estimates, p, palindromic flags), W3 extcoloc results.tsv (b PP.H4/PP.H3), W10 corroboration.tsv (b verdicts), outline v6 novelty
# rows N1-N7 (b literature codes), W2 tgfb1 results.tsv (c), W4b arg1sets results.tsv (the ARG1-T2D mark).
# Run:  Rscript "V3/R2_rewrite/figures/scripts/P31R_Fig2.R"

.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), winslash = "/")
source(file.path(dirname(.this), "P31R_fig_common.R"))
suppressPackageStartupMessages({ library(ggtext) })

FIG <- "Figure_2"
`%||%` <- function(a, b) if (is.null(a)) b else a
EN <- "–"

# ================================================================== registry (written at the point of printing)
REG <- list(ledger = list(), literals = list())
.reg_key <- function(e) paste(e$text, e$id %||% "", paste(unlist(e$idx %||% ""), collapse = "/"), e$anchor %||% "",
                              e$axis %||% "", e$why %||% "", sep = "|")
.reg_add <- function(kind, e) {
  for (i in seq_along(REG[[kind]])) if (.reg_key(REG[[kind]][[i]]) == .reg_key(e)) {
    REG[[kind]][[i]]$count <<- REG[[kind]][[i]]$count + 1L
    return(invisible(NULL))
  }
  e$count <- 1L
  REG[[kind]][[length(REG[[kind]]) + 1]] <<- e
  invisible(NULL)
}
reg_led <- function(text, id, idx = NULL, anchor = NULL, pct = FALSE) {
  e <- list(text = text, id = id)
  if (!is.null(idx)) e$idx <- if (is.numeric(idx)) as.integer(idx) else idx
  if (!is.null(anchor)) e$anchor <- anchor
  if (pct) e$pct <- TRUE
  .reg_add("ledger", e)
  text
}
L <- function(id, idx = NULL, digits = 0, anchor = NULL) reg_led(p31r_num(p31r_ledger(id, idx), digits = digits), id, idx, anchor)
# a percentage whose ledger value is already in percent (L314: ancestry shares 77, 13, 6)
LPCT <- function(id, idx, anchor) reg_led(paste0(p31r_num(p31r_ledger(id, idx), digits = 0), "%"), id, idx, anchor, pct = TRUE)
reg_ticks <- function(labels, axis, transform, what) {
  for (t in labels) .reg_add("literals", list(text = t, why = paste0("axis_tick: ", what), axis = axis, transform = transform, anchor = "^(?:[−+-]?[0-9.,]+\\s)*{}(?:\\s[−+-]?[0-9.,]+)*$"))
  labels
}
# p values (S2-D26): the whole p string is set at the panel's 7 pt, and its exponent at exactly 7 pt, raised. The base '10' of a power is the one 9-pt span:
# the figure gate reads an exponent only as a raised span smaller than 0.85 x the span before it (7 < 0.85 x 9), so a 7-pt exponent after a 7-pt '10' is
# read as a separate number and fails G7, and the gate rejects e-notation (G7f). ggtext's plain <sup> would give 5.6 pt, below the 7-pt floor.
plain_md <- function(x) gsub("<[^>]+>", "", gsub("\\*", "", x))
fmt_p <- function(p) {
  if (p >= 0.001) { s <- formatC(signif(p, 2), format = "fg", digits = 2); return(list(md = s, txt = s)) }
  e <- floor(log10(p)); m <- p / 10^e
  if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
  ms <- formatC(m, format = "f", digits = 1); es <- sub("^-", "−", as.character(e))
  list(md = sprintf("%s × <span style='font-size:9pt'>10</span><sup><span style='font-size:7pt'>%s</span></sup>", ms, es), txt = sprintf("%s × 10%s", ms, es))
}
P <- function(id, idx = NULL) { f <- fmt_p(p31r_ledger(id, idx)); reg_led(f$txt, id, idx); f }
num <- function(id, idx = NULL) as.numeric(p31r_ledger(id, idx))
TXT <- function(...) geom_text(..., size = ANNOT, family = BASE_FONT)

# text width at the annotation size, measured on the same cairo device and font the figure is saved with (mm)
.WDEV <- NULL
wmm <- function(txt, pt = ANNOT_PT, bold = FALSE) {
  if (is.null(.WDEV)) { cairo_pdf(tempfile(fileext = ".pdf"), width = 6, height = 6); .WDEV <<- dev.cur(); plot.new() }
  cur <- dev.cur(); dev.set(.WDEV); graphics::par(ps = pt)
  w <- max(vapply(strsplit(txt, "\n", fixed = TRUE)[[1]], function(s) graphics::strwidth(s, units = "inches", family = BASE_FONT, font = if (bold) 2 else 1), 0)) * 25.4
  if (cur != .WDEV && cur > 1) dev.set(cur)
  w
}
wrap_mm <- function(txt, w_mm) {            # greedy wrap on words; returns the string with "\n"
  ws <- strsplit(txt, " ", fixed = TRUE)[[1]]; out <- character(); cur <- ""
  for (w in ws) {
    cand <- if (nzchar(cur)) paste(cur, w) else w
    if (wmm(cand) <= w_mm || !nzchar(cur)) cur <- cand else { out <- c(out, cur); cur <- w }
  }
  paste(c(out, cur), collapse = "\n")
}
nlines <- function(s) length(strsplit(s, "\n", fixed = TRUE)[[1]])
LINE_MM <- 2.95                              # one 7-pt text line at lineheight 0.95, incl. leading

# ================================================================== data (fail-closed reads, pinned in the render log)
TIER_F  <- file.path(V3_ROOT, "v3.2/_engine_root/results/repro_phase1_hg38/_fig_results/phase4_annotation/table_s_pair_evidence_tiered.csv")   # canonical copy B
W11_F   <- file.path(V3_ROOT, "R1_stage1/S1R1_R2_11_extchecks/results.tsv")
W3_F    <- file.path(V3_ROOT, "R1_stage1/S1R1_R2_03_extcoloc/results.tsv")
W2_F    <- file.path(V3_ROOT, "R1_stage1/S1R1_R2_02_tgfb1/results.tsv")
W4B_F   <- file.path(V3_ROOT, "R1_stage1/S1R1_R2_04b_arg1sets/results.tsv")
CORR_F  <- file.path(V3_ROOT, "R1_stage1/S1R1_R2_10_integration/corroboration.tsv")
OUTL_F  <- file.path(P31R_ROOT, "review_R2/CLAIM_OUTLINE_STAGE3_2026-10-11_v6.md")
pe   <- p31r_read(TIER_F)
w11  <- p31r_read(W11_F, sep = "\t", colClasses = "character")
w3   <- p31r_read(W3_F, sep = "\t", colClasses = "character")
w2   <- p31r_read(W2_F, sep = "\t", colClasses = "character")
w4b  <- p31r_read(W4B_F, sep = "\t", colClasses = "character")
corr <- p31r_read(CORR_F, sep = "\t", colClasses = "character")
outl <- readLines(p31r_use(OUTL_F), encoding = "UTF-8", warn = FALSE)
qrow <- function(d, id) { r <- d[quantity_id == id]; if (nrow(r) != 1) stop("quantity id not found exactly once: ", id); r }
qv   <- function(d, id, col = "value") as.numeric(qrow(d, id)[[col]])

# ---------------------------------------------------------------- panel a invariants: what is plotted must reproduce the ledger
stopifnot(nrow(pe) == num("L010"), !anyDuplicated(pe[, .(gene, outcome)]), uniqueN(pe$gene) == num("L011"),
          uniqueN(pe$outcome) == num("L002"), !anyNA(pe$PP.H4_p12_1e6))
p6 <- pe$PP.H4_p12_1e6
TH_HI <- num("L236")   # 0.80: robust to the conservative prior = PP.H4 stays at the colocalization rule's minimum (Methods)
TH_LO <- 0.5           # one half (Results "fell below one half"); never printed, pinned by the L016 counts below
stopifnot(TH_HI == 0.8, sum(p6 >= TH_HI) == num("L016", "robust"), sum(p6 >= TH_LO & p6 < TH_HI) == num("L016", "atten"),
          sum(p6 < TH_LO) == num("L016", "fragile"))
tc <- table(pe$evidence_tier)
stopifnot(tc[["Convergent"]] == num("L014", "Convergent"), tc[["Supported"]] == num("L014", "Supported"), tc[["Flagged"]] == num("L014", "Flagged"))
stopifnot(identical(sort(paste0(pe[evidence_tier == "Convergent"]$gene, "-", pe[evidence_tier == "Convergent"]$outcome)), unlist(p31r_ledger("L015"))))
OUT_ORD <- c("T2D", "eGFR", "AF", "CAD", "HF", "Stroke", "CKD", "NAFLD")
OUT_LAB <- c(T2D = "T2D", eGFR = "eGFR", AF = "AF", CAD = "CAD", HF = "HF", Stroke = "stroke", CKD = "CKD", NAFLD = "NAFLD")
oc_n <- table(pe$outcome)
stopifnot(setequal(names(oc_n), OUT_ORD), all(diff(as.integer(oc_n[OUT_ORD])) <= 0))
for (o in OUT_ORD) stopifnot(oc_n[[o]] == num("L013", o))
xc <- table(pe$crossplatform_class)
ROWS <- data.table(
  cls = c("concordant", "discordant", "conflicting", "not assayed (panel gap)", "not evaluable (weak instrument)", "no usable variant",
          "assay discordant (within platform)"),
  lab = c("Same direction", "Opposite direction", "Cohorts disagree", "Not assayed", "Weak instrument", "No usable variant", "Assays disagree"),
  grp = c(rep("ev", 3), rep("ne", 4)),
  lid = c("L034", "L034", "L034", "L259", "L259", "L259", "L259"), lix = c(1, 2, 3, 3, 4, 5, 6))
stopifnot(setequal(names(xc), ROWS$cls))
for (i in seq_len(nrow(ROWS))) stopifnot(xc[[ROWS$cls[i]]] == num(ROWS$lid[i], ROWS$lix[i]))
stopifnot(sum(xc[ROWS[grp == "ev"]$cls]) == num("L259", 1), sum(xc[ROWS[grp == "ne"]$cls]) == num("L259", 2))
stopifnot(all(pe[PP.H4_p12_1e6 < TH_LO]$evidence_tier == "Flagged"))
L024 <- unlist(p31r_ledger("L024", "pairs"))
stopifnot(sum(startsWith(L024, "BRAP-")) == num("L024", "BRAP"), sum(startsWith(L024, "CCND2-")) == num("L024", "CCND2"),
          nrow(pe[gene == "BRAP"]) == num("L024", "BRAP"), nrow(pe[gene == "CCND2"]) == num("L024", "CCND2"),
          uniqueN(pe[gene %in% c("BRAP", "CCND2")]$crossplatform_class) == 1,
          uniqueN(pe[gene == "BRAP"]$lead_snps) == 1, uniqueN(pe[gene == "CCND2"]$lead_snps) == 1)
# Supported pairs with no evaluable arm (L226 elements 1-2): 35 of the 45; the lower block holds them
sup_ne <- pe[evidence_tier == "Supported" & crossplatform_class %in% ROWS[grp == "ne"]$cls]
stopifnot(nrow(sup_ne) == num("L226", 1), sum(pe$evidence_tier == "Supported") == num("L226", 2))
# core-complete (all core criteria): 16 = 3 Convergent + 11 Supported (no evaluable arm) + 2 Flagged (L962 core_complete)
cc <- pe[core_complete == TRUE | core_complete == "True" | core_complete == "TRUE" | core_complete == 1]
stopifnot(nrow(cc) == num("L962", c("core_complete", "core")), sum(cc$evidence_tier == "Convergent") == num("L962", c("core_complete", "convergent")),
          sum(cc$evidence_tier == "Supported") == num("L962", c("core_complete", "supported_no_arm")),
          sum(cc$evidence_tier == "Flagged") == num("L962", c("core_complete", "flagged")),
          all(cc[evidence_tier == "Supported"]$crossplatform_class %in% ROWS[grp == "ne"]$cls))
# LD admissibility (Rule 3(i)): computed for 2 of the 81 pairs; both are Supported and inadmissible (stored s above the limit)
S_LIM <- num("L923")
s_tg <- qv(w2, "S2.susie.protein.s_max"); s_ag <- qv(w4b, "S4.s.max.ARG1_T2D")
stopifnot(abs(s_tg - num("L726")) < 1e-9, abs(s_ag - num("L820")) < 1e-9, s_tg > S_LIM, s_ag > S_LIM, S_LIM == 0.1,
          num("L964", "assessed") == 2, num("L964", "not_assessed") == 79, num("L964", "n_pairs") == num("L010"),
          all(pe[(gene == "TGFB1" & outcome == "CAD") | (gene == "ARG1" & outcome == "T2D")]$evidence_tier == "Supported"))
INADM <- c(paste0("TGFB1", EN, "CAD"), paste0("ARG1", EN, "T2D"))

# ---------------------------------------------------------------- panel b invariants
PR <- data.table(pair = c("PAM-T2D", "RAD23B-T2D", "ARG1-T2D", "PCSK9-CAD", "TGFB1-CAD", "FGFR4-T2D (EAS)"),
                 k1 = c("PAM_T2D", "RAD23B_T2D", "ARG1_T2D", "PCSK9_CAD", "TGFB1_CAD", "FGFR4_T2D"),
                 ckey = c("PAM_T2D", "RAD23B_T2D", "ARG1_T2D", "PCSK9_CAD", "TGFB1_CAD", "FGFR4_T2D_EAS"),
                 lab = paste0(c("PAM", "RAD23B", "ARG1", "PCSK9", "TGFB1", "FGFR4"), EN, "T2D"),
                 resource = c("DIAMANTE", "DIAMANTE", "DIAMANTE", "CARDIoGRAMplusC4D", "CARDIoGRAMplusC4D", "TPMI, East Asian"),
                 ledger_ids = c("L040", "L041", "L042", "L043", "L044", "L045"), lidx = 1:6)
PR[4:5, lab := paste0(c("PCSK9", "TGFB1"), EN, "CAD")]
stopifnot(nrow(PR) == num("L268"))
est <- rbindlist(lapply(seq_len(nrow(PR)), function(i) {
  k <- PR$k1[i]; w <- qrow(w11, paste0("C1.", k, ".wald")); d <- qrow(w11, paste0("C1.", k, ".disc_wald"))
  lid <- PR$ledger_ids[i]
  # the W11 re-derivation holds the same Wald ratio, SE and p as the ledger's six rows (L040-L045); FGFR4's stored discovery is full precision
  stopifnot(abs(as.numeric(w$value) - num(lid, 1)) < 1e-9, abs(as.numeric(w$se) - num(lid, 2)) < 1e-9, abs(as.numeric(w$p) - num(lid, 3)) < 1e-12,
            abs(as.numeric(d$value) - num(lid, 5)) < 5e-5)
  pres <- qv(w11, paste0("C1.", k, ".p_resource"))
  stopifnot(abs(pres / num("L960", PR$lidx[i]) - 1) < 1e-9)       # the resource's own p (L960) classifies
  cl <- qrow(w11, paste0("C1.", k, ".class"))$value
  stopifnot(cl == "concordant", sign(as.numeric(w$value)) == sign(as.numeric(d$value)))
  data.table(pair = PR$pair[i], ext = as.numeric(w$value), ext_lo = as.numeric(w$ci_lo), ext_hi = as.numeric(w$ci_hi), ext_se = as.numeric(w$se),
             disc = as.numeric(d$value), disc_se = as.numeric(d$se), disc_lo = as.numeric(d$value) - 1.959964 * as.numeric(d$se),
             disc_hi = as.numeric(d$value) + 1.959964 * as.numeric(d$se), p_res = pres,
             pal = qrow(w11, paste0("C1.", k, ".palindromic"))$value, pal_note = qrow(w11, paste0("C1.", k, ".palindromic"))$notes)
}))
stopifnot(all(abs(est$ext_lo - (est$ext - 1.959964 * est$ext_se)) < 1e-6), all(abs(est$ext_hi - (est$ext + 1.959964 * est$ext_se)) < 1e-6))   # stored CI = Wald, first-order SE
# palindromic instruments: PAM-T2D oriented by allele frequency (frequencies agree), RAD23B-T2D by a forward-strand assumption
stopifnot(identical(est[pal == "TRUE"]$pair, c("PAM-T2D", "RAD23B-T2D")), grepl("frequencies agree", est[pair == "PAM-T2D"]$pal_note),
          grepl("forward-strand assumption", est[pair == "RAD23B-T2D"]$pal_note), all(est[pal == "FALSE"]$pal_note == "strand status: not palindromic"))
# corroboration verdicts, computed by the integration step from the key ids; three pairs corroborated (L985)
cm <- function(pk) { r <- corr[pkey == pk]; stopifnot(nrow(r) == 1); r }
vd <- rbindlist(lapply(seq_len(nrow(PR)), function(i) { r <- cm(PR$ckey[i]); data.table(pair = PR$pair[i], verdict = r$verdict, s3 = r$s3_branch, s2 = r$s2_branch) }))
CORR_PAIRS <- unlist(p31r_ledger("L985", "pairs"))
stopifnot(sum(vd$verdict == "CORROBORATED") == num("L985", "n"), setequal(sub(" \\(EAS\\)", "", vd[verdict == "CORROBORATED"]$pair), CORR_PAIRS),
          vd[pair == "RAD23B-T2D"]$s3 == "4", vd[pair == "PCSK9-CAD"]$s3 == "2 (with add-on)", vd[pair == "TGFB1-CAD"]$s2 == "(a)",
          all(vd[verdict != "CORROBORATED"]$pair %in% c("RAD23B-T2D", "PCSK9-CAD", "TGFB1-CAD")))
REASON <- c("RAD23B-T2D" = "shared variant unresolved", "PCSK9-CAD" = "conservative prior", "TGFB1-CAD" = "opposite conditional signs")
# external regional colocalization (W3): PAM, RAD23B, TGFB1, PCSK9 only; no S3 rows exist for ARG1-T2D or FGFR4-T2D, so those two cells read "not run"
stopifnot(!any(grepl("^S3\\.(ARG1|FGFR4)", w3$quantity_id)), all(c("PAM_T2D", "RAD23B_T2D", "TGFB1_CAD", "PCSK9_CAD") %in% sub("^S3\\.([A-Z0-9]+_[A-Z0-9]+)\\..*", "\\1", w3$quantity_id)))
s3 <- function(k, pr, h) qv(w3, sprintf("S3.%s.pmaf.%s.%s", k, pr, h))
for (z in list(c("L743", "PAM_T2D", "1e-05", "H4"), c("L744", "PAM_T2D", "1e-06", "H4"), c("L745", "PCSK9_CAD", "1e-06", "H3"), c("L746", "PCSK9_CAD", "1e-05", "H4"),
               c("L747", "PCSK9_CAD", "1e-06", "H4"), c("L748", "RAD23B_T2D", "1e-05", "H3"), c("L749", "RAD23B_T2D", "1e-06", "H3"), c("L750", "RAD23B_T2D", "1e-05", "H4"),
               c("L751", "RAD23B_T2D", "1e-06", "H4"), c("L752", "TGFB1_CAD", "1e-05", "H4"), c("L753", "TGFB1_CAD", "1e-06", "H4")))
  stopifnot(abs(s3(z[2], z[3], z[4]) - num(z[1])) < 1e-9)
stopifnot(num("L240", 3) == 1e-5, num("L241") == 1e-6)         # the two priors the cells are labelled with (never printed as numbers here)
# literature status: codes are generated from the novelty rows of outline v6 (N1-N7), asserted against the row text, each scoped to its claim
nl0 <- grep("^## Novelty ledger", outl); nl1 <- grep("^## Moved out", outl); stopifnot(length(nl0) == 1, length(nl1) == 1, nl1 > nl0)
nrow_txt <- function(k) { i <- nl0 + grep(paste0("^\\| ", k, " \\|"), outl[nl0:nl1]) - 1; stopifnot(length(i) == 1); outl[i] }
LIT <- data.table(pair = PR$pair, nid = c("N2", "N3", "N2", "N1", "N1", "N7"),
  code = c("pair reported before", "not found (main texts only; Table 2)", "pair reported before", "estimate published; regional test new",
           "estimate published; regional test new", "reported before"),
  must = c("PAM-T2D named in Yi", "not found in 15 main texts read; supplements NOT RUN", "ARG1-T2D in Yuan 2023: confirmatory",
           "single-SNP CAD estimates were published before.*regional test is the new element",
           "TGFB1 instrument is not stated in the main text.*reproduction of the TGFB1 estimate is not established", "FGFR4-T2D was reported before"))
for (i in seq_len(nrow(LIT))) stopifnot(grepl(LIT$must[i], nrow_txt(LIT$nid[i])))
# scope test: an absence code carries its scope; "new" appears only as "regional test new" (the N1 element); nothing reads "none found"
stopifnot(!any(grepl("none found", LIT$code)), all(!grepl("not found", LIT$code) | grepl("main texts only", LIT$code)),
          all(!grepl("\\bnew\\b", LIT$code) | grepl("regional test new", LIT$code)), !any(grepl("^new$|^not found$", LIT$code)))

# ---------------------------------------------------------------- panel c invariants (W2: conditional ratios, aware, second-order SE, panel 503)
RS <- c(finngen_chd = "FinnGen CHD", ieu_a_7 = "CARDIoGRAMplusC4D")
cdt <- rbindlist(lapply(names(RS), function(r) rbindlist(lapply(c("rs73045269", "rs12983775"), function(v) {
  x <- qrow(w2, sprintf("S2.ratio.%s.%s.value", r, v))
  data.table(res = r, snp = v, val = as.numeric(x$value), lo = as.numeric(x$ci_lo), hi = as.numeric(x$ci_hi), se = as.numeric(x$se),
             panel = x$panel, cov = x$cov, so = x$se_order)
}))))
stopifnot(all(cdt$panel == "503", all(cdt$cov == "aware"), all(cdt$so == "2nd")),
          abs(cdt[res == "finngen_chd" & snp == "rs73045269"]$val - num("L732")) < 1e-9, abs(cdt[res == "finngen_chd" & snp == "rs12983775"]$val - num("L731")) < 1e-9,
          abs(cdt[res == "ieu_a_7" & snp == "rs73045269"]$val - num("L734")) < 1e-9, abs(cdt[res == "ieu_a_7" & snp == "rs12983775"]$val - num("L733")) < 1e-9,
          all(abs(cdt$lo - (cdt$val - 1.959964 * cdt$se)) < 1e-4), all(abs(cdt$hi - (cdt$val + 1.959964 * cdt$se)) < 1e-4),
          cdt[snp == "rs73045269", all(val > 0)], cdt[snp == "rs12983775", all(val < 0)])
stopifnot(abs(qv(w2, "S2.het.finngen_chd.p_aware_2nd") - num("L728")) < 1e-18, abs(qv(w2, "S2.het.ieu_a_7.p_aware_2nd") - num("L729")) < 1e-20,
          abs(qv(w2, "S2.ld.rs1800470_rs12983775.r2") - num("L730")) < 1e-9)
# "sensitive to the reference LD": conditional colocalization at rs12983775 falls below the 0.90 bootstrap rule in both resources (L715, L716, L924)
stopifnot(abs(qv(w2, "S2.coloc_clause.boot.rs12983775.finngen_chd.share") - num("L715")) < 1e-9, abs(qv(w2, "S2.coloc_clause.boot.rs12983775.ieu_a_7.share") - num("L716")) < 1e-9,
          num("L715") < num("L924"), num("L716") < num("L924"), num("L085", 4) > 32, num("L085", 4) < 33)
cdt[, ylab := fifelse(snp == "rs73045269", "Discovery", "East Asian")]

# ================================================================== panel a: evidence map (as before, plus the S2-D24 notes and marks)
DX <- 0.024
lane_of <- function(x, lanes) {
  off <- rep(NA_real_, length(x)); used <- vector("list", length(lanes))
  for (i in order(x)) for (l in seq_along(lanes)) {
    if (all(abs(used[[l]] - x[i]) >= DX)) { used[[l]] <- c(used[[l]], x[i]); off[i] <- lanes[l]; break }
  }
  if (anyNA(off)) stop("lane overflow")
  off
}
a <- copy(pe)[, x := PP.H4_p12_1e6]
a[, pair := paste0(gene, EN, outcome)]
a[, locus := fifelse(gene %in% c("BRAP", "CCND2"), gene, "")]
a <- merge(a, ROWS[, .(crossplatform_class = cls, lab, grp)], by = "crossplatform_class")
LABPOS <- data.table(pair = paste0(c("ARG1", "IL6R", "PAM", "RAD23B", "PCSK9", "TGFB1"), EN, c("T2D", "CAD", "T2D", "T2D", "CAD", "CAD")),
                     tier = c(1, 1, 1, 2, 2, 1), hj = c(0.5, 1, 0.5, 1, 1, 1), dxl = c(0, 0.006, 0, 0.008, 0.012, 0.037))
stopifnot(all(LABPOS$pair %in% a$pair), all(a[evidence_tier == "Convergent"]$pair %in% LABPOS$pair), all(INADM %in% LABPOS$pair))
IL6R_X <- a[pair == paste0("IL6R", EN, "CAD")]$x
stopifnot(abs(IL6R_X - num("L204", 3)) < 5e-4, IL6R_X < TH_HI)
ON_HI <- a[abs(x - TH_HI) < DX / 2]
stopifnot(nrow(ON_HI) >= 1, all(ON_HI$x < TH_HI), paste0("FN1", EN, "CAD") %in% ON_HI$pair)
LABPOS[, txt := pair]
LABPOS[pair == paste0("IL6R", EN, "CAD"), txt := paste0(pair, " (", L("L204", 3, digits = 3), ")")]
LABPOS[pair %in% INADM, txt := paste0(pair, "‡")]               # LD-admissibility mark: a glyph plus a key entry, not a hue
a[, labelled := pair %in% LABPOS$pair]
a[, `:=`(off = NA_real_, y = NA_real_)]
LOC_UP <- 1.1
a[locus != "", loff := { o <- order(x); r <- numeric(.N); last <- -Inf
  for (j in o) { if (x[j] - last < DX) r[j] <- LOC_UP else last <- x[j] }; r }, by = locus]
stopifnot(a[locus != "" & loff > 0, .N] == 1, a[locus != "" & loff > 0]$pair == paste0("BRAP", EN, "Stroke"))
SUMM_LANES <- 1.6
TIERH <- 1.6                                   # lane units between label tiers (a 7-pt label box is about 1.2 lanes high at this panel height)
ytab <- data.table()
ycur <- 0
for (i in rev(seq_len(nrow(ROWS)))) {
  r <- ROWS[i]
  inrow <- a$crossplatform_class == r$cls
  sel <- inrow & a$locus == "" & !a$labelled
  has_lab <- any(inrow & a$labelled)
  lane_off <- lane_of(a$x[sel], if (has_lab) c(0, -1, -2, -3) else c(0, 1, -1, 2, -2))
  a$off[sel] <- lane_off
  lo <- min(lane_off); hi <- max(lane_off)
  if (has_lab) {
    stopifnot(all(diff(sort(a$x[inrow & a$labelled])) >= DX))
    hi <- hi + 1; a$off[inrow & a$labelled] <- hi
  }
  ntier <- if (has_lab) max(LABPOS$tier[LABPOS$pair %in% a$pair[inrow]]) else 0
  if (r$cls == "not assayed (panel gap)") {
    a$off[a$gene == "CCND2"] <- hi + 1.4; a$off[a$gene == "BRAP"] <- hi + 2.6; hi <- hi + 2.6
    ntier <- 1
  }
  pad_hi <- 1.1 + TIERH * ntier + if (i == 1) SUMM_LANES else 0
  if (r$cls == "assay discordant (within platform)") pad_hi <- pad_hi + 2.4   # room for the three-line definition of the not-evaluable rows INSIDE this one row band (S2-D25: no note crosses a divider)
  base <- ycur + 0.9 - lo
  a$y[inrow] <- base + a$off[inrow]
  ytab <- rbind(ytab, data.table(cls = r$cls, lab = r$lab, grp = r$grp, lid = r$lid, lix = r$lix, ymid = base + (lo + hi) / 2,
                                 ylo = ycur, yhi = base + hi + pad_hi))
  ycur <- base + hi + pad_hi
  if (i == 4) ycur <- ycur + 3
}
stopifnot(!anyNA(a$y))
a[locus != "", y0line := y]
a[locus != "", y := y + loff]
ytab$txt <- NA_character_
for (k in seq_len(nrow(ytab))) {
  n_pts <- sum(a$crossplatform_class == ytab$cls[k])
  stopifnot(n_pts == num(ytab$lid[k], ytab$lix[k]))
  ytab$txt[k] <- paste0(ytab$lab[k], " (", L(ytab$lid[k], ytab$lix[k], anchor = paste0(ytab$lab[k], " \\({}\\)")), ")")
}
GRP_LAB <- c(ev = paste0("Evaluable (", L("L259", 1), ")"), ne = paste0("Not evaluable (", L("L259", 2), ")"))
stopifnot(sum(a$grp == "ev") == num("L259", 1), sum(a$grp == "ne") == num("L259", 2))
a[, grpf := factor(GRP_LAB[grp], levels = GRP_LAB)]
ytab[, grpf := factor(GRP_LAB[grp], levels = GRP_LAB)]
TIER_LAB <- c(Convergent = paste0("Convergent (", L("L014", "Convergent", anchor = "Convergent \\({}\\)"), ")"),
              Supported = paste0("Supported (", L("L014", "Supported", anchor = "Supported \\({}\\)"), ")"),
              Flagged = paste0("Flagged (", L("L014", "Flagged"), ")"))
for (t in names(TIER_LAB)) stopifnot(grepl(paste0("(", sum(a$evidence_tier == t), ")"), TIER_LAB[[t]], fixed = TRUE))
a[, tierf := factor(evidence_tier, levels = names(PAL_TIER))]
lp <- merge(a[labelled == TRUE, .(pair, x, y, grpf)], LABPOS, by = "pair")
lp[, `:=`(xl = x + dxl, yl = y + 0.25 + TIERH * tier)]
seps <- ytab[, .(y = ylo[ylo > min(ylo)]), by = grpf]
stopifnot(nrow(seps) == nrow(ytab) - 2)
loc <- a[locus != "", .(x0 = min(x), x1 = max(x), y = y0line[1], n = .N), by = .(locus, grpf)]
stopifnot(all(loc[, .N, by = locus]$N == 1), all(a[locus != "", .N, by = .(locus, y0line)][, .N, by = locus]$N == 1))
stubs <- a[locus != "" & loff > 0, .(x, y0 = y0line, y1 = y, grpf)]
stopifnot(loc[locus == "BRAP"]$n == num("L024", "BRAP"), loc[locus == "CCND2"]$n == num("L024", "CCND2"))
loc[, txt := fifelse(locus == "BRAP", paste0("BRAP: ", L("L024", "BRAP"), " pairs"),
                     paste0("CCND2: ", L("L024", "CCND2", anchor = "CCND2: {} pairs"), " pairs"))]
loc[, `:=`(lx = fifelse(locus == "BRAP", x0 - 0.012, x0 - 0.014), ly = fifelse(locus == "BRAP", y + 1.0, y),
           lhj = fifelse(locus == "BRAP", 0, 1))]
stopifnot(loc[locus == "CCND2"]$x0 < TH_LO, loc[locus == "BRAP"]$x0 - 0.012 > TH_LO)
BAND <- data.table(x = c(0.25, 0.65, 0.90), lab = c(paste(L("L016", "fragile"), "pairs"), paste(L("L016", "atten"), "pairs"),
                                                     paste(L("L016", "robust"), "pairs")))
XB <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
XB_LAB <- reg_ticks(p31r_lab(XB), "2a_x_pph4", "linear", "panel a x, PP.H4")
SUMM1 <- paste0(L("L010", anchor = "{} pairs, "), " pairs, ", L("L011"), " proteins, ", L("L002", anchor = "{} outcomes"),
                " outcomes; all pass the colocalization rule at the default prior")
top_ev <- ytab[grp == "ev", max(yhi)]
summ_df <- data.table(grpf = factor(GRP_LAB[["ev"]], levels = GRP_LAB), x = 0.012, y = top_ev - 0.3, txt = SUMM1)
# S2-D24 notes and the exposure-side row definitions (moved out of the legend into the graphic), in free areas of panel a: the strip left of
# the data (no pair has PP.H4 below 0.3), the right of the "Opposite direction" row (no point right of 0.60), the left of the "Cohorts
# disagree" row, and the lower left under the last two rows (no point left of 0.78). Line breaks are explicit so each registered number
# sits in one printed line with its own anchor.
DEF_SAME <- paste(c("Same direction: every evaluable", "cohort had the discovery sign."), collapse = "\n")
NOTE_EV <- paste(c(
  "Core-complete (all core",
  paste0("criteria): ", L("L962", c("core_complete", "core"), anchor = "criteria\\): {}"), " = ",
         L("L014", "Convergent", anchor = "= {} Convergent"), " Convergent"),
  paste0("+ ", L("L962", c("core_complete", "supported_no_arm"), anchor = "\\+ {} Supported"), " Supported + ",
         L("L962", c("core_complete", "flagged"), anchor = "\\+ {} Flagged"), " Flagged.")), collapse = "\n")
# the LD-admissibility note is a separate block inside ONE row band
NOTE_LD <- paste(c(
  "‡ Reference-LD fit inadmissible;",
  paste0("tier not rebuilt. ", L("L964", "not_assessed", anchor = "{} of 81"), " of ", L("L010", anchor = "of {} pairs"), " pairs"),
  "not assessed for admissibility."), collapse = "\n")
DEF_OPP <- "none had the discovery sign"
DEF_COH <- "one of the two had the discovery sign"
DEF_NE <- wrap_mm(paste0("Not evaluable: not assayed, by neither cohort; weak instrument, F-statistic below ", L("L235", anchor = "below {}"),
                         "; no usable variant, lead variant missing; assays disagree, two assays of the protein differed in sign."), 84)
MIN_X <- min(a$x)
NOTE_W_MAX <- 36   # mm: the left strip beside the data is free to x = 0.30 (no pair left of 0.31; 1 unit = 125 mm, text starts at 0.012)
stopifnot(MIN_X > 0.30, max(vapply(strsplit(c(NOTE_EV, NOTE_LD), "\n"), function(v) max(vapply(v, wmm, 0)), 0)) < NOTE_W_MAX,
          max(vapply(strsplit(DEF_SAME, "\n")[[1]], wmm, 0)) < 40)   # DEF_SAME sits above the first data row, so it may be wider
# the plotted data must leave those areas free: no point right of 0.60 in "Opposite direction", none left of 0.62 in "Cohorts disagree",
# none left of 0.77 in the last two rows
stopifnot(max(a[crossplatform_class == "discordant"]$x) < 0.62, min(a[crossplatform_class == "conflicting"]$x) > 0.62,
          min(a[crossplatform_class %in% c("no usable variant", "assay discordant (within platform)")]$x) > 0.77,
          wmm(DEF_OPP) < 0.36 * 118, wmm(DEF_COH) < 0.60 * 118)
ne_row <- ytab[cls == "not evaluable (weak instrument)"]; na_row <- ytab[cls == "not assayed (panel gap)"]
op_row <- ytab[cls == "discordant"]; co_row <- ytab[cls == "conflicting"]
nu_row <- ytab[cls == "no usable variant"]; as_row <- ytab[cls == "assay discordant (within platform)"]
CCND2_Y <- a[gene == "CCND2"]$y[1]
conc_row <- ytab[cls == "concordant"]
note_df <- data.table(grpf = factor(GRP_LAB[c("ev", "ne", "ne")], levels = GRP_LAB), x = 0.012, txt = c(NOTE_EV, NOTE_LD, DEF_NE),
                      y = c(top_ev - 4.9, CCND2_Y - 0.9, (as_row$ylo + as_row$yhi) / 2), vj = c(1, 1, 0.5),
                      band_lo = c(conc_row$ylo, na_row$ylo, as_row$ylo), band_hi = c(conc_row$yhi, na_row$yhi, as_row$yhi))
# every note lies wholly inside ONE row band (lane units; one 7-pt line, box padding included, is 1.3 lanes at this panel height)
LPL <- 1.3
note_df[, n := vapply(strsplit(txt, "\n", fixed = TRUE), length, 0L)]
note_df[, `:=`(top = fifelse(vj == 1, y, y + n * LPL / 2), bot = fifelse(vj == 1, y - n * LPL, y - n * LPL / 2))]
stopifnot(note_df[, all(top <= band_hi - 0.1 & bot >= band_lo + 0.1)])
note_df[, c("n", "top", "bot", "band_lo", "band_hi") := NULL]
def_df <- data.table(grpf = factor(GRP_LAB[c("ev", "ev", "ev")], levels = GRP_LAB), x = c(1.0, 0.012, 0.012), y = c(op_row$ymid, co_row$ymid, top_ev - 1.8),
                     txt = c(DEF_OPP, DEF_COH, DEF_SAME), hj = c(1, 0, 0), vj = c(0.5, 0.5, 1))

TOP_EV <- max(a[grp == "ev"]$y)
VL <- rbind(data.table(grpf = factor(GRP_LAB[["ev"]], levels = GRP_LAB), x = c(TH_LO, TH_HI), y0 = -Inf, y1 = TOP_EV),
            data.table(grpf = factor(GRP_LAB[["ne"]], levels = GRP_LAB), x = c(TH_LO, TH_HI), y0 = -Inf, y1 = Inf))
stopifnot(all(lp[grpf == GRP_LAB[["ev"]]]$yl - 0.6 > TOP_EV + 0.45),
          max(lp[grpf == GRP_LAB[["ev"]]]$yl) + 0.6 < min(summ_df$y) - 1.0)
LBL <- function(d, fill) geom_label(data = d, aes(x = xl, y = yl, label = txt, hjust = hj), inherit.aes = FALSE, vjust = 0.5, size = ANNOT,
                                    family = BASE_FONT, fill = fill, border.colour = NA, label.padding = unit(0.5, "pt"), label.r = unit(0, "pt"))
XTITLE_A <- paste0("PP.H4 under the conservative prior; robust at ", L("L236", digits = 2), " or above")
blk <- ytab[, .(y = c(min(ylo), max(yhi))), by = grpf]
pa <- ggplot(a, aes(x, y)) +
  geom_blank(data = blk, aes(x = 0.5, y = y), inherit.aes = FALSE) +
  geom_hline(data = seps, aes(yintercept = y), colour = "grey80", linewidth = 0.25) +
  geom_segment(data = VL, aes(x = x, xend = x, y = y0, yend = y1), inherit.aes = FALSE, linetype = "22", linewidth = 0.35, colour = "grey45") +
  geom_segment(data = loc, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5, colour = "grey25") +
  geom_segment(data = stubs, aes(x = x, xend = x, y = y0, yend = y1), inherit.aes = FALSE, linewidth = 0.5, colour = "grey25") +
  geom_text(data = loc, aes(x = lx, y = ly, label = txt, hjust = lhj), inherit.aes = FALSE, vjust = 0.5, size = ANNOT, family = BASE_FONT) +
  geom_segment(data = lp, aes(x = x, xend = x, y = y, yend = yl - 0.42), inherit.aes = FALSE, linewidth = 0.25, colour = "grey30") +
  LBL(lp, "white") +
  TXT(data = summ_df, aes(x = x, y = y, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 1) +
  geom_label(data = def_df, aes(x = x, y = y, label = txt, hjust = hj, vjust = vj), inherit.aes = FALSE, size = ANNOT, family = BASE_FONT, colour = "grey25", lineheight = 0.95,
             fill = "white", border.colour = NA, label.padding = unit(0.5, "pt"), label.r = unit(0, "pt")) +
  geom_label(data = note_df, aes(x = x, y = y, label = txt, vjust = vj), inherit.aes = FALSE, hjust = 0, size = ANNOT, family = BASE_FONT, lineheight = 0.95,
             fill = "white", border.colour = NA, label.padding = unit(0.5, "pt"), label.r = unit(0, "pt")) +
  geom_point(aes(fill = tierf), shape = 21, size = 2.1, stroke = 0.3, colour = "grey15") +
  facet_grid(grpf ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = PAL_TIER, labels = TIER_LAB, name = "Evidence tier", drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.02), breaks = XB, labels = XB_LAB, expand = c(0, 0),
                     sec.axis = dup_axis(breaks = BAND$x, labels = BAND$lab, name = NULL)) +
  scale_y_continuous(breaks = ytab$ymid, labels = ytab$txt, expand = expansion(add = c(0.35, 0.35))) +
  labs(x = XTITLE_A, y = "Exposure-side corroboration (deCODE, ARIC)") +
  theme_p31r() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
        axis.ticks.y = element_blank(), axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(),
        strip.text.y = element_text(angle = -90), panel.spacing.y = unit(3, "pt"),
        legend.key.spacing.x = unit(6, "pt")) +
  guides(fill = guide_legend(override.aes = list(size = 2.6)))

# ================================================================== helpers for the two table-style panels (b, c): a 170-mm frame, 1 unit = 1 mm
LH <- 0.85                                                         # line height of the in-panel text (multiple of the font size)
LINE_MM <- 2.52                                                    # measured pitch of one 7-pt line at LH (R: LH x 1.2 x 7 pt)
mm_base <- function(h) ggplot() +
  coord_cartesian(xlim = c(0, 170), ylim = c(0, h), expand = FALSE, clip = "off") +
  theme_void(base_family = BASE_FONT) +
  theme(plot.margin = margin(0, 0, 0, 0), plot.tag.position = c(0, 1), plot.tag.location = "panel",
        plot.tag = element_text(size = TAG_PT, face = "bold", family = BASE_FONT, hjust = 0, vjust = 1))
TXT0 <- function(d, hj = 0, vj = 0.5, col = "black", fontface = "plain", size = ANNOT) geom_text(data = d, aes(x = x, y = y, label = txt), inherit.aes = FALSE,
  hjust = hj, vjust = vj, size = size, family = BASE_FONT, colour = col, fontface = fontface, lineheight = LH)
RTX <- function(d, hj = 0, vj = 0.5, col = "black") geom_richtext(data = d, aes(x = x, y = y, label = txt), inherit.aes = FALSE, hjust = hj, vjust = vj,
  size = ANNOT, family = BASE_FONT, colour = col, fill = NA, label.colour = NA, label.padding = unit(0, "pt"), label.margin = unit(0, "pt"), lineheight = LH)
AXT <- BODY_PT / ggplot2::.pt              # axis tick labels: 9 pt (locked house rule)
SEG <- function(d, lw = 0.35, col = "black", lt = "solid") geom_segment(data = d, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = lw, colour = col, linetype = lt)
PTM <- function(d, shape, fill, col, size = 2, stroke = 0.4) geom_point(data = d, aes(x = x, y = y), inherit.aes = FALSE, shape = shape, fill = fill, colour = col, size = size, stroke = stroke)
chk_w9 <- function(txt, w, what) { mm <- regmatches(txt, regexec("^(.* × )10(−[0-9]+)$", txt))[[1]]   # 7 pt text, 9-pt '10', 7-pt exponent
  m <- if (length(mm)) wmm(mm[2], pt = 7) + wmm("10", pt = 9) + wmm(mm[3], pt = 7) else wmm(txt, pt = 7); if (m > w) stop(what, " overflows its column: ", round(m, 1), " > ", w, " mm: ", txt); invisible(m) }
chk_w <- function(txt, w, what, bold = FALSE) { m <- max(vapply(strsplit(txt, "\n", fixed = TRUE)[[1]], wmm, 0, bold = bold)); if (m > w) stop(what, " overflows its column: ", round(m, 1), " > ", w, " mm: ", gsub("\n", "/", txt)); invisible(m) }

# ================================================================== panel b: the six pairs selected in advance
XC <- list(lab = 5, w_lab = 27.5, f0 = 34, fw = 28.5, p = 64, w_p = 16, h = 82, w_h = 29.5, v = 113, w_v = 28.5, l = 143, w_l = 26.5)
XL <- c(-0.75, 0.42); X2MM <- function(v) XC$f0 + (v - XL[1]) / diff(XL) * XC$fw
stopifnot(min(c(est$disc_lo, est$ext_lo)) > XL[1], max(c(est$disc_hi, est$ext_hi)) < XL[2])
CELL <- data.table(pair = PR$pair)
CELL[, p_md := mapply(function(i) P("L960", i)$md, PR$lidx)]
NOTRUN <- "not run; the rule needed concordance only"
CELL[, h4 := c(paste0(L("L743", digits = 2), " / ", L("L744", digits = 2)),
               paste0(L("L750", digits = 2), " / ", L("L751", digits = 4), "\nPP.H3 ", L("L748", digits = 2), " / ", L("L749", digits = 2)),
               wrap_mm(NOTRUN, XC$w_h),
               paste0(L("L746", digits = 2), " / ", L("L747", digits = 2), "\nPP.H3 ", L("L745", digits = 2), ", conservative"),
               paste0(L("L752", digits = 2), " / ", L("L753", digits = 2)),
               wrap_mm(NOTRUN, XC$w_h))]
CELL[, h4_run := !(pair %in% c("ARG1-T2D", "FGFR4-T2D (EAS)"))]
CELL[, verdict := vapply(seq_len(.N), function(i) if (vd$verdict[i] == "CORROBORATED") "Corroborated"
                         else paste0("Not corroborated:\n", wrap_mm(REASON[[pair[i]]], XC$w_v)), "")]
CELL[, vgrey := vd$verdict != "CORROBORATED"]
CELL[, lit := vapply(LIT$code, function(s) wrap_mm(s, XC$w_l), "")]
CELL[, lab1 := paste0(PR$lab, c("†", "§", "", "", "", ""))]
CELL[, lab2 := PR$resource]
stopifnot(identical(CELL$pair[PR$pair %in% est[pal == "TRUE"]$pair], c("PAM-T2D", "RAD23B-T2D")))
CELL[, fgfr := ifelse(pair == "FGFR4-T2D (EAS)", wrap_mm("sign and p only", XC$p - XC$f0 - 1.9), "")]
CELL[, nl := vapply(seq_len(.N), function(i) max(nlines(h4[i]), nlines(verdict[i]), nlines(lit[i]), if (nzchar(fgfr[i])) nlines(fgfr[i]) else 1L, 2L), 0)]
CELL[, rh := LINE_MM * nl + 1.3]
for (i in seq_len(nrow(CELL))) {
  chk_w(CELL$h4[i], XC$w_h, "regional column"); chk_w(CELL$verdict[i], XC$w_v, "verdict column"); chk_w(CELL$lit[i], XC$w_l, "literature column")
  if (wmm(CELL$lab1[i], pt = BODY_PT) > XC$w_lab) stop("pair label overflows its column at 9 pt: ", CELL$lab1[i]); chk_w(CELL$lab2[i], XC$w_lab, "resource label")
  if (nzchar(CELL$fgfr[i])) chk_w(CELL$fgfr[i], XC$p - XC$f0 - 1.5, "FGFR4 forest text")
  chk_w9(plain_md(CELL$p_md[i]), XC$w_p, "p column")
}
HDR_L <- c(p = "External\np", h = "Regional PP.H4,\ndefault / conservative", v = "Corroboration verdict\n(rule in Methods)", l = "Literature status\n(scoped; Table 2)")
chk_w(HDR_L[["h"]], XC$w_h, "header", TRUE); chk_w(HDR_L[["v"]], XC$w_v, "header", TRUE); chk_w(HDR_L[["l"]], XC$w_l, "header", TRUE)
HDR_H <- 2 * LINE_MM + 1.7
# key block (inside the graphic): marker key with the palindrome marks, then the ancestry statement (the ARG1 choice is stated once, in the Results: R13)
KEY_W <- 163
K1 <- list(c("sq", "discovery (FinnGen)"), c("dm", "external check"), c("t", "† palindromic, oriented by allele frequency"), c("t", "§ palindromic, forward-strand assumption"))
GAP <- 5
widths <- vapply(K1, function(z) (if (z[1] == "t") 0 else 3.4) + wmm(z[2]), 0)
if (sum(widths) + GAP * (length(K1) - 1) > KEY_W) stop("marker key line too long: ", round(sum(widths) + GAP * 3, 1))
K_TXT <- c(paste0("CARDIoGRAMplusC4D is mixed ancestry (", LPCT("L314", 3, "{} European"), " European, ",
                  LPCT("L314", 4, "{} South Asian"), " South Asian, ", LPCT("L314", 5, "{} East Asian"),
                  " East Asian)."))
K_TXT <- vapply(K_TXT, wrap_mm, "", w_mm = KEY_W)
stopifnot(all(vapply(K_TXT, function(s) max(vapply(strsplit(s, "\n")[[1]], wmm, 0)) <= KEY_W, TRUE)))
K_N <- 1L + sum(vapply(K_TXT, nlines, 0L))
AXB_H <- 6.0                                                       # axis line + tick labels (9 pt)
K_H <- LINE_MM * K_N + 1.4
Hb <- AXB_H + K_H + sum(CELL$rh) + HDR_H
CELL[, ytop := Hb - HDR_H - c(0, cumsum(rh)[-.N])]
CELL[, ymid := ytop - rh / 2][, ybot := ytop - rh]
hdr <- data.table(x = c(XC$h, XC$v, XC$l), y = Hb - 0.9, txt = unname(HDR_L[c("h", "v", "l")]))
hdr_md <- data.table(x = XC$p, y = Hb - 0.9, txt = "External<br>*p*")
hdr_ax <- data.table(x = XC$f0 + XC$fw - 1.5, y = Hb - 0.9, txt = "Log-odds per SD of\nprotein (95% CI)")
chk_w(hdr_ax$txt, 40, "forest header", TRUE)
EB <- merge(est, CELL[, .(pair, ymid)], by = "pair")[pair != "FGFR4-T2D (EAS)"]
OFF <- 1.25
disc_df <- EB[, .(x = X2MM(disc), y = ymid + OFF, xlo = X2MM(disc_lo), xhi = X2MM(disc_hi))]
ext_df  <- EB[, .(x = X2MM(ext), y = ymid - OFF, xlo = X2MM(ext_lo), xhi = X2MM(ext_hi))]
XLB <- c(-0.6, -0.3, 0, 0.3)
XLB_LAB <- reg_ticks(p31r_lab(XLB), "2b_x_logodds", "linear", "panel b x, log-odds per SD of protein")
y_ax <- K_H + AXB_H - 0.6
zero_seg <- data.table(x = X2MM(0), xend = X2MM(0), y = CELL[pair == "TGFB1-CAD", ybot], yend = CELL[pair != "FGFR4-T2D (EAS)", max(ytop)] - 0.3)
tick_seg <- data.table(x = X2MM(XLB), xend = X2MM(XLB), y = y_ax, yend = y_ax - 0.9)
ax_seg <- data.table(x = X2MM(XLB[1]), xend = X2MM(XLB[length(XLB)]), y = y_ax, yend = y_ax)
rowsep <- data.table(x = XC$lab, xend = 169.5, y = CELL$ybot[-nrow(CELL)], yend = CELL$ybot[-nrow(CELL)])
rule_top <- data.table(x = XC$lab, xend = 169.5, y = Hb - HDR_H, yend = Hb - HDR_H)
# key lines, bottom-up: the marker line sits lowest
yk <- LINE_MM * (seq_len(K_N) - 0.5) + 0.6
kx <- XC$lab; key_items <- list()
for (j in seq_along(K1)) {
  z <- K1[[j]]
  if (z[1] == "t") { key_items[[j]] <- list(kind = "t", x = kx, txt = z[2]); kx <- kx + wmm(z[2]) + GAP
  } else { key_items[[j]] <- list(kind = z[1], x = kx + 1.2, txt = z[2], xt = kx + 3.4); kx <- kx + 3.4 + wmm(z[2]) + GAP }
}
key_sq <- data.table(x = key_items[[1]]$x, y = yk[1]); key_dm <- data.table(x = key_items[[2]]$x, y = yk[1])
key_t <- rbindlist(lapply(key_items, function(k) data.table(x = k$xt %||% k$x, y = yk[1], txt = k$txt)))
kt <- data.table(txt = unlist(strsplit(K_TXT, "\n", fixed = TRUE)))
kt[, `:=`(x = XC$lab, y = yk[1] + LINE_MM * rev(seq_len(.N)))]
pc_df <- data.table(x = XC$p, y = CELL$ymid, txt = CELL$p_md)
# S8 step 3 (F1): pair labels at the theme's 9-pt axis-text size (AXT); the FGFR4 forest note is shortened to 'sign and p only' (ancestry wording is in the legend)
lab1_df <- data.table(x = XC$lab, y = CELL$ymid + 1.6, txt = CELL$lab1)
lab2_df <- data.table(x = XC$lab, y = CELL$ymid - 1.55, txt = CELL$lab2)
h4_df <- data.table(x = XC$h, y = CELL$ymid, txt = CELL$h4)
h4_run <- h4_df[CELL$h4_run == TRUE]; h4_nr <- h4_df[CELL$h4_run == FALSE]
vd_df <- data.table(x = XC$v, y = CELL$ymid, txt = CELL$verdict, grey = CELL$vgrey)
lit_df <- data.table(x = XC$l, y = CELL$ymid, txt = CELL$lit)
fg_df <- data.table(x = XC$f0 + 0.5, y = CELL[pair == "FGFR4-T2D (EAS)"]$ymid, txt = CELL[pair == "FGFR4-T2D (EAS)"]$fgfr)
pb <- mm_base(Hb) +
  SEG(rowsep, 0.25, "grey80") + SEG(rule_top, 0.35, "grey45") +
  SEG(zero_seg, 0.35, "grey45", "22") + SEG(ax_seg, 0.35) + SEG(tick_seg, 0.35) +
  geom_segment(data = disc_df, aes(x = xlo, xend = xhi, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.45, colour = PAL_NEUTRAL[["mid"]]) +
  geom_segment(data = ext_df, aes(x = xlo, xend = xhi, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.45, colour = PAL_NEUTRAL[["dark"]]) +
  PTM(disc_df, 22, PAL_NEUTRAL[["light"]], PAL_NEUTRAL[["mid"]]) + PTM(ext_df, 23, PAL_NEUTRAL[["dark"]], PAL_NEUTRAL[["dark"]]) +
  TXT0(data.table(x = X2MM(XLB), y = y_ax - 2.6, txt = XLB_LAB), hj = 0.5, size = AXT) +
  TXT0(lab1_df, size = AXT) + TXT0(lab2_df, col = "grey35") +
  TXT0(fg_df, col = "grey20") +
  RTX(pc_df) + RTX(hdr_md, vj = 1) + TXT0(hdr_ax, hj = 1, vj = 1, fontface = "bold") +
  TXT0(h4_run) + TXT0(h4_nr, col = "grey35", fontface = "italic") +
  TXT0(vd_df[grey == FALSE]) + TXT0(vd_df[grey == TRUE], col = "grey35") +
  TXT0(lit_df) +
  TXT0(hdr, vj = 1, fontface = "bold") +
  PTM(key_sq, 22, PAL_NEUTRAL[["light"]], PAL_NEUTRAL[["mid"]]) + PTM(key_dm, 23, PAL_NEUTRAL[["dark"]], PAL_NEUTRAL[["dark"]]) +
  TXT0(key_t) + TXT0(kt) +
  labs(tag = "b")

# ================================================================== panel c: TGFB1-CAD, estimates at the two instruments
CX <- list(lab = 5, sig = 33, f0 = 53, fw = 32, v = 88, t = 100)
XLC <- c(-0.6, 0.6); XC2MM <- function(v) CX$f0 + (v - XLC[1]) / diff(XLC) * CX$fw
stopifnot(min(cdt$lo) > XLC[1], max(cdt$hi) < XLC[2])
ROWH_C <- 3.9; AXC_H <- 6.0
cdt[, ord := match(paste(res, snp), paste(rep(names(RS), each = 2), rep(c("rs73045269", "rs12983775"), 2)))]
setorder(cdt, ord)
HC_HDR <- 3.4
pf <- P("L728"); pc2 <- P("L729")
# each p value ends its line: the exponent is a separate raised 7-pt span only while no 7-pt text follows it on the line (the PDF writer merges
# same-size neighbours into one span, which the gate then cannot read as an exponent)
TB <- c("Difference between instruments:",
        paste0("FinnGen CHD: *p* ", pf$md),
        paste0("CARDIoGRAMplusC4D: *p* ", pc2$md),
        paste0("Instruments ", L("L085", 4, anchor = "{} kb"), " kb apart. Opposite protein-scaled estimates"),
        "from different instruments do not show opposite causal",
        "effects. Conditional colocalization at the East Asian",
        "instrument is sensitive to the reference LD.")
TB_PLAIN <- plain_md(TB)
# a line's width: 7 pt for the text, plus the 9-pt '10' of each power of ten it carries
PV_TXT <- c(pf$txt, pc2$txt)
w_tb <- function(i) wmm(TB_PLAIN[i]) + sum(vapply(PV_TXT[vapply(PV_TXT, grepl, TRUE, x = TB_PLAIN[i], fixed = TRUE)], function(n) wmm("10", 9) - wmm("10"), 0))
stopifnot(all(vapply(seq_along(TB), w_tb, 0) <= 169 - CX$t))
NTB <- length(TB)
TB_D <- c(3.0, 3.5, 3.0, rep(LINE_MM, NTB - 4))             # centre-to-centre pitch of the lines (mm): the two 9-pt p lines need more than a 7-pt line
TB_Y0 <- 1.7                                                # first line centre below the panel top
Hc <- max(AXC_H + 4 * ROWH_C + HC_HDR, TB_Y0 + sum(TB_D) + LINE_MM / 2 + 0.3)
cy <- function(i) Hc - HC_HDR - ROWH_C * (i - 0.5)
cdt[, y := cy(ord)]
cax <- AXC_H - 0.6
XCB <- c(-0.5, 0, 0.5)
XCB_LAB <- reg_ticks(p31r_lab(XCB), "2c_x_logodds", "linear", "panel c x, log-odds per SD of protein")
cd_pts <- cdt[, .(x = XC2MM(val), y = y, snp = snp)]
cd_seg <- cdt[, .(x = XC2MM(lo), xend = XC2MM(hi), y = y, yend = y, snp = snp)]
c_ax <- data.table(x = XC2MM(XCB[1]), xend = XC2MM(XCB[3]), y = cax, yend = cax)
c_tk <- data.table(x = XC2MM(XCB), xend = XC2MM(XCB), y = cax, yend = cax - 0.9)
c_zero <- data.table(x = XC2MM(0), xend = XC2MM(0), y = cax, yend = Hc - HC_HDR)
sigcol <- c(rs73045269 = PAL_SIGNAL[["A"]], rs12983775 = PAL_SIGNAL[["B"]])
cval <- data.table(x = CX$v + 8, y = cdt$y, txt = vapply(seq_len(nrow(cdt)), function(i)
  L(c(finngen_chd_rs73045269 = "L732", finngen_chd_rs12983775 = "L731", ieu_a_7_rs73045269 = "L734", ieu_a_7_rs12983775 = "L733")[[paste0(cdt$res[i], "_", cdt$snp[i])]], digits = 3), ""))
c_res <- data.table(x = CX$lab, y = c(mean(cdt[res == "finngen_chd"]$y), mean(cdt[res == "ieu_a_7"]$y)), txt = unname(RS[c("finngen_chd", "ieu_a_7")]))
stopifnot(max(vapply(cdt$ylab, wmm, 0)) <= CX$f0 - CX$sig - 1)   # the instrument labels end before the forest
c_sig <- data.table(x = CX$sig, y = cdt$y, txt = cdt$ylab)
tb_df <- data.table(x = CX$t, y = Hc - TB_Y0 - c(0, cumsum(TB_D)), txt = TB)
stopifnot(min(tb_df$y) - LINE_MM / 2 > 0.2)
c_hdr <- data.table(x = CX$f0 - 0.5, y = Hc - 0.5, txt = "Log-odds per SD of protein (95% CI)")
pcpl <- mm_base(Hc) +
  SEG(c_zero, 0.35, "grey45", "22") + SEG(c_ax, 0.35) + SEG(c_tk, 0.35) +
  geom_segment(data = cd_seg, aes(x = x, xend = xend, y = y, yend = yend, colour = snp), inherit.aes = FALSE, linewidth = 0.55) +
  geom_point(data = cd_pts, aes(x = x, y = y, fill = snp), inherit.aes = FALSE, shape = 21, colour = "grey15", size = 2.3, stroke = 0.4) +
  scale_colour_manual(values = sigcol, guide = "none") + scale_fill_manual(values = sigcol, guide = "none") +
  TXT0(data.table(x = XC2MM(XCB), y = cax - 2.6, txt = XCB_LAB), hj = 0.5, size = AXT) +
  TXT0(c_res, vj = 0.5) + TXT0(c_sig) + TXT0(cval, hj = 1) +
  TXT0(data.table(x = CX$sig, y = Hc - 0.5, txt = "Instruments"), vj = 1, fontface = "italic") +
  TXT0(c_hdr, hj = 0, vj = 1, fontface = "italic") +
  RTX(tb_df, vj = 0.5) +
  labs(tag = "c")

# ================================================================== compose + save
HA <- as.numeric(Sys.getenv("P31R_FIG2_HA", "105.4"))     # total height of panel a: panel, legend, axis and title (b and c are fixed-size blocks below it)
H <- HA + Hb + Hc
message(sprintf("panel heights: a %.1f  b %.1f  c %.1f  total %.1f mm", HA, Hb, Hc, H))
stopifnot(H <= 225)
fig <- (pa + labs(tag = "a")) / wrap_elements(full = pb) / wrap_elements(full = pcpl) + plot_layout(heights = unit(c(1, Hb, Hc), c("null", "mm", "mm"))) & theme(plot.margin = margin(0, 0, 0, 0))
p31r_save(fig, FIG, w_mm = FIG_W_MM, h_mm = H)


# ================================================================== plotted frames (source data written at save time)
dfd <- file.path(P31R_FIGS, "data", "Figure_2"); dir.create(dfd, showWarnings = FALSE, recursive = TRUE)
fwrite(a[, .(gene, outcome, pair, evidence_tier, crossplatform_class, PP.H4_conservative = PP.H4_p12_1e6, core_complete, x_plotted = x, y_plotted = y)], file.path(dfd, "fig2a_plotted.tsv"), sep = "\t")
fwrite(merge(est, CELL[, .(pair, p_md, h4, verdict, lit)], by = "pair"), file.path(dfd, "fig2b_plotted.tsv"), sep = "\t")
fwrite(cdt[, .(res, snp, val, lo, hi, se)], file.path(dfd, "fig2c_plotted.tsv"), sep = "\t")
writeLines(jsonlite::toJSON(list(figure = FIG, note = "plotted frames written by P31R_Fig2.R at save time; every value is read from the pinned sources in render_log/Figure_2.tsv; discovery 95% CI = value +/- 1.959964 x se (arithmetic on stored values; the W11 discovery rows hold value and se only); external CI = stored W11 ci_lo/ci_hi (equal to value +/- 1.959964 x se, asserted)",
                                 h4_not_run_pairs = c("ARG1-T2D", "FGFR4-T2D (EAS)"), literature_codes = LIT[, .(pair, nid, code)]), auto_unbox = TRUE, pretty = TRUE),
           file.path(dfd, "fig2_plotted.provenance.json"), useBytes = TRUE)

# ================================================================== registry (written from the same calls that printed the numbers)
.norm <- function(t) gsub(" ", "", gsub("−", "-", t))
alltx <- c(vapply(REG$ledger, function(e) .norm(e$text), ""), vapply(REG$literals, function(e) .norm(e$text), ""))
for (kind in c("ledger", "literals")) for (e in REG[[kind]]) {
  if (sum(alltx == .norm(e$text)) > 1 && is.null(e$anchor)) stop("registry: shared text '", e$text, "' (", e$id %||% e$axis, ") has no anchor")
  stopifnot(is.integer(e$count), e$count >= 1L)
}
out <- list(panels = c("a", "b", "c"), ledger = REG$ledger, derived = list(), literals = REG$literals)
regd <- file.path(P31R_FIGS, "registry"); dir.create(regd, showWarnings = FALSE)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(regd, paste0(FIG, ".json")), useBytes = TRUE)
message("registry written: ", length(out$ledger), " ledger, ", length(out$literals), " literal entries")
if (!is.null(.WDEV)) try(dev.off(.WDEV), silent = TRUE)
