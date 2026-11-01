# W1 COMPARE_SEEN: run AFTER w1_main.R. Seen values are transcribed from the FULL pre-declaration
# (P31R_R2_STAGE1_PRESPEC.md, md5 3319dc39..., section 1 'Seen' block, lines 167-183; section 12 line 590-595),
# opened only after results.tsv and branches.tsv were written. Stage 1 values are read from results.tsv, never retyped.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
OUTDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_01_rad23b")
Sys.setenv(W1_LOG = file.path(OUTDIR, "run.log"))
source(file.path(OUTDIR, "w1_common.R"))
assert_md5(GOV$prespec_full[1], GOV$prespec_full[2], "GOV:prespec_full (opened for COMPARE_SEEN only)")
assert_md5(GOV$amend_full[1], GOV$amend_full[2], "GOV:amend_full (opened for COMPARE_SEEN only)")
# the outputs must be byte-identical to the freeze taken before the full text was opened
fr <- readLines(file.path(OUT, "PRE_COMPARE_md5.txt")); fr <- fr[grepl("^[0-9a-f]{32} \\*", fr)]
for (l in fr) { m <- sub(" .*$", "", l); f <- sub("^[0-9a-f]{32} \\*", "", l)
  if (!identical(md5(file.path(OUT, f)), m)) stop("output changed after the pre-compare freeze: ", f) }
# documents read by the implementer for METHOD (not inputs of any computation), recorded here with md5; no seen script or seen output was opened
for (m in list(c(paste0(.P31R_V3, "/paper3_deposit_request_v4/scripts/scripts/paper3.1_ckm/cascade/01_mr_screen.R"), "6c70f3195107279a9552a10fd72bdf19", "pipeline proxy and palindrome rules (orient_proxy, read_local_outcome): the reason a proxy was used"),
               c(paste0(.P31R_V3, "/v3.2/_engine_root/results/e_source/R/00_config.R"), "cc58fcf6be207052f792315589853307", "ENG_engine_config_00_config.R: primary_fdr_set (BH over the padded grid) and bonferroni_alpha"),
               c(paste0(.P31R_V3, "/R2_rewrite/supplement/P31R_SM_grading.md"), "9c37a37cb3ac40c4e232b53f36b5de46", "SM grading: the tier rules (Rule 1 'tier rules dated in SM grading')")))
  assert_md5(m[1], m[2], paste0("METHOD-READ: ", m[3]))
logf("seen scripts opened: NONE; seen outputs opened: NONE; seen corpus not opened at any point")
logf("COMPARE_SEEN: %d outputs byte-identical to the pre-compare freeze (PRE_COMPARE_md5.txt)", length(fr))

# NA audit on the RAW text (w1_main.R's own audit read the file with fread, which turns the string "NA" into a missing
# value, so its log line under-counted; the values are unaffected). Every NA value must be an expected absence.
raw <- strsplit(readLines(file.path(OUT, "results.tsv"))[-1], "\t", fixed = TRUE)
stopifnot(all(lengths(raw) %in% c(14L, 15L)))   # a trailing empty 'notes' cell is dropped by strsplit
ids <- vapply(raw, `[`, "", 1); vals <- vapply(raw, `[`, "", 8)
na_ids <- ids[vals %in% c("", "NA")]
expected_na <- grepl("^S1[.]ledger[.](PECAM1_CAD|ABO_Stroke|ABO_T2D)[.](eaf_ukbppp_queried|beta_exp_queried|se_exp_queried|p_exp_queried)$", na_ids) |
               grepl("^S1[.]ledger[.]ABO_Stroke[.](eaf_outcome_lead|beta_out_lead|se_out_lead|p_out_lead|outcome_lead_rsid_in_file)$", na_ids)
logf("NA audit (raw text of results.tsv, %d rows): %d rows with an NA/empty value, all expected absences = %s: %s", length(ids), length(na_ids), all(expected_na), paste(na_ids, collapse = ", "))
logf("  expected because: PECAM1 has no local UKB-PPP window; the ABO proxy rs576123 is absent from the ABO window; the ABO lead is absent from MEGASTROKE")
if (!all(expected_na)) stop("unexpected NA values: ", paste(na_ids[!expected_na], collapse = ", "))
R <- fread(file = file.path(OUT, "results.tsv"), sep = "\t", quote = "", colClasses = "character")
sv <- function(id, field = "value") { x <- R[quantity_id == id][[field]]; if (length(x) != 1L) stop("no unique result ", id); x }
half_unit <- function(s) {
  s <- sub("^[-+]", "", trimws(as.character(s)))
  if (grepl("[eE]", s)) { mant <- sub("[eE].*$", "", s); ex_ <- as.integer(sub("^.*[eE]", "", s))
    d <- if (grepl("[.]", mant)) nchar(sub("^[^.]*[.]", "", mant)) else 0L; return(0.5 * 10^(ex_ - d)) }
  d <- if (grepl("[.]", s)) nchar(sub("^[^.]*[.]", "", s)) else 0L
  0.5 * 10^(-d)
}
rows <- list()
add <- function(qid, seen, seen_label, s1, s1_label, explanation = "", numeric = TRUE) {
  if (numeric) {
    d <- abs(as.numeric(s1) - as.numeric(seen)); tol <- half_unit(seen) * (1 + 1e-9)
    flag <- if (d <= tol) "within printed precision" else "DIFFERENCE beyond printed precision"
    rows[[length(rows) + 1L]] <<- data.table(quantity_id = qid, seen_value = seen, seen_label = seen_label, stage1_value = s1, stage1_label = s1_label,
      abs_diff = trimws(formatC(d, digits = 3, format = "g")), explanation = paste(flag, explanation, sep = "; "))
  } else rows[[length(rows) + 1L]] <<- data.table(quantity_id = qid, seen_value = seen, seen_label = seen_label, stage1_value = s1, stage1_label = s1_label,
      abs_diff = "NA (label)", explanation = explanation)
}
weak <- "deterministic value from the same raw inputs; scoping used Python (pre-declaration section 10 table), Stage 1 used R data.table, PLINK --r and BEDMatrix, so agreement checks the tool stacks, not independent data"
L <- "full text section 1 Seen (l.167-183)"
add("S1.directlead.RAD23B_T2D.beta", "0.0307", L, sv("S1.directlead.RAD23B_T2D.beta"), "FinnGen R10 T2D at rs7035725, per G (forward strand)", weak)
add("S1.directlead.RAD23B_T2D.se", "0.0066", L, sv("S1.directlead.RAD23B_T2D.se"), "SE of the lead association", weak)
add("S1.directlead.RAD23B_T2D.p", "3.49e-6", L, sv("S1.directlead.RAD23B_T2D.p"), "lead p", weak)
add("S1.directlead.RAD23B_T2D.wald", "-0.4786", L, sv("S1.directlead.RAD23B_T2D.wald"), "forward-strand Wald at the lead", weak)
add("S1.directlead.RAD23B_T2D.wald_se", "0.1032", L, sv("S1.directlead.RAD23B_T2D.wald_se"), "first-order Wald SE at the lead", weak)
add("S1.rad23b.r2_plink", "0.9958", L, sv("S1.rad23b.r2_plink"), "PLINK --r square r2, 503 panel", weak)
add("S1.ledger.RAD23B_T2D.r_plink_A1A1", "0.998", L, sv("S1.ledger.RAD23B_T2D.r_plink_A1A1"), "signed r between lead A1 (C) and proxy A1 (G)", paste("seen 'r +0.998';", weak))
add("S1.ledger.RAD23B_T2D.phase_plink", "lead C in phase with proxy G", L, sv("S1.ledger.RAD23B_T2D.phase_plink"), "PLINK --r in-phase: lead C with proxy G, lead G with proxy T", "same phase", numeric = FALSE)
add("S1.ledger.RAD23B_T2D.p_out_queried", "3.15e-6", L, sv("S1.ledger.RAD23B_T2D.p_out_queried"), "FinnGen T2D p at the proxy rs7856112", weak)
add("S1.rad23b.cut.lead (threshold)", "3.2084e-6", L, as.character(0.05 / 15584), "0.05/15,584 (grid derived from the instrument table and mr_all_results)", "Bonferroni threshold")
add("S1.ledger.RAD23B_T2D.proxy_reason", "the lead is palindromic with EAF 0.55, not because it was missing", L, sv("S1.ledger.RAD23B_T2D.proxy_reason"), "computed reason", "same reason (Stage 1 adds the pipeline rule: EAF inside (0.42, 0.58))", numeric = FALSE)
add("S1.ledger.RAD23B_T2D.wald_se2", "0.1179", L, sv("S1.ledger.RAD23B_T2D.wald_se2"), "second-order SE of the proxy Wald", weak)
add("S1.ledger.RAD23B_T2D.wald_p2", "4.59e-5", L, sv("S1.ledger.RAD23B_T2D.wald_p2"), "second-order p of the proxy Wald", weak)
add("S1.rad23b.pp_h4_conservative", "0.8273", L, "NOT A W1 QUANTITY", "not computed", "conservative-prior PP.H4 is the TAB_tiered81 column PP.H4_p12_1e6, used by W1 only through crit_prior_stable; W1 did not recompute it", numeric = FALSE)
add("S1.rad23b.decode_F", "[value withheld, release rule F15]", L, sv("S1.rad23b.decode_F"), "deCODE F (table value; back-derived F in notes)", "read from TAB_decode_mr_results, so agreement is expected")
add("S1.rad23b.decode_p2", "0.00855", L, sv("S1.rad23b.decode_p2"), "2*Phi(-|decode_b/delta_se|), recomputed", "recomputed from the same table columns, so agreement is expected")
add("S1.ledger.ADM_AF.r2_plink", "1.0000", L, sv("S1.ledger.ADM_AF.r2_plink"), "PLINK r2 lead-proxy", weak)
add("S1.ledger.PECAM1_CAD.r2_plink", "1.0000", L, sv("S1.ledger.PECAM1_CAD.r2_plink"), "PLINK r2 lead-proxy", weak)
add("S1.ledger.ABO_T2D.r2_plink", "0.9771", L, sv("S1.ledger.ABO_T2D.r2_plink"), "PLINK r2 lead-proxy (same for ABO_Stroke)", weak)
add("S1.directlead.ADM_AF.p", "1.99e-6", L, sv("S1.directlead.ADM_AF.p"), "FinnGen AF p at the lead rs2403310", weak)
add("S1.directlead.PECAM1_CAD.p", "9.99e-7", L, sv("S1.directlead.PECAM1_CAD.p"), "FinnGen CHD p at the lead rs1962147", paste("seen adds 'class unchanged': Stage 1 lead p is below 0.05/15,584 as the proxy p is, so the first-order class is unchanged;", weak))
add("S1.directlead.ABO_T2D.p", "2.79e-18", L, sv("S1.directlead.ABO_T2D.p"), "FinnGen T2D p at the lead (row labelled rs920065566)", weak)
add("S1.ledger.ABO_T2D.outcome_lead_rsid_in_file", "lead present under rs920065566", L, sv("S1.ledger.ABO_T2D.outcome_lead_rsid_in_file"), "rsID on the FinnGen row at the lead position", "same; consistent with section 12 WP-B (ABO proxied although present by position)", numeric = FALSE)
add("S1.directlead.ABO_Stroke.wald", "ABO-stroke cannot be run, because the lead is absent from MEGASTROKE", "pre-declaration text (also in the redacted copy)", sv("S1.directlead.ABO_Stroke.wald"), "NOT RUN", "same: absent by GRCh38 position and by rsID over the whole streamed file", numeric = FALSE)
add("S1.se.core_bonf.1st", "41", L, sv("S1.se.core_bonf.1st"), "Bonferroni within the 81, first order", weak)
add("S1.se.core_bonf.2nd", "30", L, sv("S1.se.core_bonf.2nd"), "Bonferroni within the 81, second order", weak)
add("S1.se.core_fdr.1st", "80", L, sv("S1.se.core_fdr.1st"), "family BH within the 81, first order", weak)
add("S1.se.core_fdr.2nd", "71", L, sv("S1.se.core_fdr.2nd"), "family BH within the 81, second order", weak)
add("S1.se.family_bonf.1st", "94", L, sv("S1.se.family_bonf.1st"), "family Bonferroni, first order", weak)
add("S1.se.family_bonf.2nd", "72", paste(L, "and the assertion target (visible in the redacted copy)"), sv("S1.se.family_bonf.2nd"), "family Bonferroni, pipeline convention", "the target was visible before computing: weak evidence")
add("S1.se.family_bh.1st", "292", L, sv("S1.se.family_bh.1st"), "family BH, first order", weak)
add("S1.se.family_bh.2nd", "248", paste(L, "and the assertion target (visible in the redacted copy)"), sv("S1.se.family_bh.2nd"), "family BH, pipeline convention", "the target was visible before computing: weak evidence")
add("S1.se.movers_tier", "1 (only RAD23B changes tier)", L, paste(sv("S1.se.movers_tier"), sv("S1.se.movers_tier", "notes")), "tier movers by the SM grading rule", "same (RAD23B Convergent -> Supported under second-order SE; descriptive, Rule 2)", numeric = FALSE)
add("S1.se.movers_bonf", "11 (named lists asserted; names not transcribed in the pre-declaration)", "full text section 1 Analyses (l.164)", paste(sv("S1.se.movers_bonf"), "|", sv("S1.se.movers_bonf", "notes")), "Bonferroni movers within the 81",
    "count agrees (target visible in the redacted copy: weak evidence); the seen NAMES are not in the pre-declaration and the seen corpus was not opened, so the names are compared with nothing", numeric = FALSE)
add("S1.se.movers_fdr", "9 (named lists asserted; names not transcribed in the pre-declaration)", "full text section 1 Analyses (l.164)", paste(sv("S1.se.movers_fdr"), "|", sv("S1.se.movers_fdr", "notes")), "FDR movers within the 81",
    "count agrees (target visible: weak evidence); names compared with nothing (not transcribed; seen corpus not opened)", numeric = FALSE)
add("S1.rad23b.ext.resource_p_lead", "0.038", L, sv("S1.rad23b.ext.resource_p_lead"), "DIAMANTE printed p at rs7035725", "printed value read from the same file")
add("S1.rad23b.ext.resource_p_proxy", "0.042", L, sv("S1.rad23b.ext.resource_p_proxy"), "DIAMANTE printed p at rs7856112", "printed value read from the same file")
add("S1.rad23b.ext.wald1_p", "0.043", L, sv("S1.rad23b.ext.wald1_p"), "first-order rounded-Wald p", weak)
add("S1.rad23b.ext.wald1_lo", "0.036", L, sv("S1.rad23b.ext.wald1_lo"), "rounding interval, lower", weak)
add("S1.rad23b.ext.wald1_hi", "0.051", L, sv("S1.rad23b.ext.wald1_hi"), "rounding interval, upper", weak)
add("S1.rad23b.ext.wald2_p", "0.049", L, sv("S1.rad23b.ext.wald2_p"), "second-order rounded-Wald p", weak)
add("S1.rad23b.ext.wald2_lo", "0.042", L, sv("S1.rad23b.ext.wald2_lo"), "rounding interval, lower", weak)
add("S1.rad23b.ext.wald2_hi", "0.057", L, sv("S1.rad23b.ext.wald2_hi"), "rounding interval, upper", weak)
add("S1.rad23b.ext.bonf6", "fails 0.05/6 under every p source", L, paste(sv("S1.rad23b.ext.bonf6"), "|", sv("S1.rad23b.ext.bonf6", "notes")), "resource p and both Wald p (with intervals) against 0.0083", "same: not met under any p source", numeric = FALSE)
add("S1.orient.RAD23B_T2D", "UKB-PPP at rs7856112: T beta -0.0637, giving Wald -0.483, the same direction as the lead", L, paste(sv("S1.orient.RAD23B_T2D"), "|", sv("S1.orient.RAD23B_T2D", "notes")), "exposure-side orientation", "same direction (CONSISTENT); T beta -0.06374, phase-free Wald -0.4834", numeric = FALSE)
add("S1.ledger.RAD23B_T2D.beta_exp_queried", "-0.0637", L, sv("S1.ledger.RAD23B_T2D.beta_exp_queried"), "UKB-PPP beta at rs7856112 for T (the in-phase allele)", "read from the same file")
add("RNG", "n/a", "brief rule 8", "W1 draws no random numbers (no bootstrap is declared for section 1)", "", "no 'same RNG and seed' note needed", numeric = FALSE)
add("CLASS_FLIP", "n/a", "brief: s classifications", "W1 computes no s (no SuSiE in section 1)", "", "no CLASS_FLIP possible", numeric = FALSE)
CS <- rbindlist(rows)
fwrite(CS, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = FALSE)
nd <- sum(grepl("DIFFERENCE", CS$explanation))
logf("COMPARE_SEEN.tsv written: %d rows; numeric differences beyond printed precision: %d; md5 %s", nrow(CS), nd, md5(file.path(OUT, "COMPARE_SEEN.tsv")))
