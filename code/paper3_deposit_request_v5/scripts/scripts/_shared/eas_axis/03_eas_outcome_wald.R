# ============================================================================
# EAS AXIS 03: harmful-direction conversion + Wald + BH (contract SS8, SS10 item 3)
# ----------------------------------------------------------------------------
# Consumes 02's coverage table -- which already carries the located source row
# and the orientation decision -- and computes for every qc_evaluable pair the
# EAS Wald estimate under the same first-order formula as L1, in harmful
# orientation, with two-sided p and a BH q over the primary evaluable family
# (strict VARIANT_LEVEL family recomputed separately, pre-spec).
# 03 does NOT re-match and does NOT re-orient: matching and orientation are
# 02's single responsibility; this script only applies orientation_action.
# SECONDARY magnitude fields (z_diff/p_diff/q_diff) compare the two MR-SCALE
# estimates -- betaMR_EAS against the frozen EUR MR beta with its SE -- never
# an outcome-association beta against an MR beta (they are different scales).
# They never decide states (05's job).
# Includes the mandatory EUR regression test FIRST: the same harmful-direction
# function must reproduce the frozen L4 estimates for all 76 pairs (16 eGFR
# sign checks in particular) BEFORE any EAS estimate is written.
# ============================================================================
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
suppressPackageStartupMessages(library(data.table))
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"); setwd(B)
RCG     <- file.path("results", paste0("repro_", PHASE1_SET))
EAS_DIR <- file.path(RCG, "eas_axis")
cov <- fread(file.path(EAS_DIR, "eas_outcome_coverage.csv"))
stopifnot(nrow(cov) > 0L)

harmful <- function(beta, outcome) ifelse(outcome == "eGFR", -beta, beta)
wald    <- function(bY, bX)  bY / bX
wald_se <- function(seY, bX) seY / abs(bX)
# the Wald pair is exercised here so a regression in either is loud:
stopifnot(wald(1.2, -0.6) == -2, wald_se(0.3, -0.6) == 0.5)

# ── EUR regression test first: reproduce frozen L4 with the same function ────
t4  <- fread(file.path(RCG, "table_s_pair_evidence_tiered.csv"))
mr  <- as.data.table(readRDS(file.path(RCG, "mr_all_results.rds")))
chk <- merge(t4[, .(gene, outcome, mr_beta_frozen = mr_beta)],
             mr[, .(gene, outcome = outcome_name, b_raw = b, se_eur = se)],
             by = c("gene", "outcome"))
# the merge must recover every canonical pair -- that is the regression test, not the number
stopifnot(nrow(chk) == nrow(t4))
chk[, b_re := harmful(b_raw, outcome)]
# TOLERANCE 5e-7, NOT 1e-8: the L4 table stores mr_beta at six decimals, so a
# full-precision estimate legitimately differs from it by up to half an ULP of
# the sixth decimal (measured max on this generation: 4.9149e-7). 1e-8 flagged
# 74/76 correct pairs as failures (EASr2 finding 1). Sign equality stays exact.
bad <- chk[abs(b_re - mr_beta_frozen) > 5e-7 | sign(b_re) != sign(mr_beta_frozen)]
if (nrow(bad))
  stop("03: EUR regression test FAILED for ", nrow(bad), " pair(s) (e.g. ",
       bad$gene[1], "-", bad$outcome[1], ") -- the orientation function does not ",
       "reproduce frozen L4; writing no EAS estimate")
stopifnot(chk[outcome == "eGFR", .N] == 16L)
cat("03: EUR regression test OK -- 76/76 reproduce frozen L4 (16 eGFR signs verified)\n")

# ── instrument side: betaX/seX/F from the frozen table; F_lo from winners curse ──
ins <- as.data.table(readRDS(file.path("results", "phase1_hg38", "ukb_cis_instruments.rds")))
est <- copy(cov)
est[ins, on = c(gene = "gene", instrument_rsID = "SNP"),
    `:=`(betaX_UKB = i.beta, seX_UKB = i.se, F_original = i.F_stat)]
stopifnot(sum(is.na(est$betaX_UKB)) == 0L, all(is.finite(est$betaX_UKB)),
          all(est$betaX_UKB != 0))
wc <- fread(file.path(RCG, "winners_curse", "exposure_selection_corrected.csv"))
stopifnot(all(c("gene", "SNP", "F_lo") %in% names(wc)))
est[wc, on = c(gene = "gene", instrument_rsID = "SNP"), F_lo := i.F_lo]

# ── EAS estimates: apply 02's orientation, then harmful, then Wald ───────────
stopifnot(all(!est$qc_evaluable | est$orientation_action %in% c("AS_IS", "FLIP")))
est[, betaY_EAS_raw := fifelse(qc_evaluable & orientation_action == "FLIP",
                               -betaY_src, betaY_src)]          # aligned to UKB effect allele
est[qc_evaluable & orientation_action == "FLIP", eaf := 1 - eaf]
est[, betaY_EAS_harmful := fifelse(qc_evaluable, harmful(betaY_EAS_raw, outcome), NA_real_)]
est[, seY_EAS := fifelse(qc_evaluable, seY_src, NA_real_)]
est[, betaMR_EAS := fifelse(qc_evaluable, harmful(wald(betaY_EAS_raw, betaX_UKB), outcome), NA_real_)]
est[, seMR_EAS   := fifelse(qc_evaluable, wald_se(seY_EAS, betaX_UKB), NA_real_)]
est[, ci95_low_EAS  := betaMR_EAS - 1.96 * seMR_EAS]            # un-rounded, SS10.1
est[, ci95_high_EAS := betaMR_EAS + 1.96 * seMR_EAS]
est[, p_eas := fifelse(qc_evaluable, 2 * pnorm(-abs(betaMR_EAS / seMR_EAS)), NA_real_)]

# ── BH over the primary evaluable family; strict-basis family recomputed ─────
est[, q_eas := NA_real_]
fam <- which(est$qc_evaluable == TRUE)
est$q_eas[fam] <- p.adjust(est$p_eas[fam], method = "BH")
est[, q_eas_strict := NA_real_]
famS <- which(est$qc_evaluable == TRUE & est$sample_count_basis == "VARIANT_LEVEL")
est$q_eas_strict[famS] <- p.adjust(est$p_eas[famS], method = "BH")   # pre-spec: recomputed BH

# ── secondary magnitude fields, MR scale on BOTH sides; never state-deciding ─
est[chk, on = .(gene, outcome), `:=`(betaMR_EUR_harmful = i.b_re, seMR_EUR = i.se_eur)]
est[, magnitude_testable := qc_evaluable == TRUE & outcome != "CKD" &
      is.finite(betaMR_EUR_harmful) & is.finite(seMR_EUR) & is.finite(betaMR_EAS)]
est[, magnitude_untestable_reason := fifelse(magnitude_testable, NA_character_,
      fifelse(qc_evaluable != TRUE, "NOT_EVALUABLE",
      fifelse(outcome == "CKD", "PHENOTYPE_SCALE_NOT_COMPARABLE",
              "DIFFERENCE_STAT_MISSING")))]
est[, z_diff := fifelse(magnitude_testable,
      (betaMR_EAS - betaMR_EUR_harmful) / sqrt(seMR_EAS^2 + seMR_EUR^2), NA_real_)]
est[, p_diff := fifelse(is.finite(z_diff), 2 * pnorm(-abs(z_diff)), NA_real_)]
famD <- which(is.finite(est$p_diff))
est[, q_diff := NA_real_]; est$q_diff[famD] <- p.adjust(est$p_diff[famD], method = "BH")

fwrite(est, file.path(EAS_DIR, "eas_outcome_estimates.csv"))
cat(sprintf("03: %d/76 estimated | BH family=%d | strict family=%d | magnitude-testable=%d\n",
            sum(est$qc_evaluable), length(fam), length(famS), sum(est$magnitude_testable)))

# ── sensitivity sources (contract SS4.1): same estimator, SEPARATE families ──
# One BH family PER sensitivity source; never merged into the primary family
# and never state-deciding (05 reads the primary file only).
SEN_F <- file.path(EAS_DIR, "eas_outcome_coverage_sensitivity.csv")
if (!file.exists(SEN_F))
  stop("03: ", SEN_F, " absent -- run the current 02 first; the three registered ",
       "sensitivity sources must be analyzed, not merely localized")
sen <- fread(SEN_F)
sen[ins, on = c(gene = "gene", instrument_rsID = "SNP"),
    `:=`(betaX_UKB = i.beta, seX_UKB = i.se, F_original = i.F_stat)]
stopifnot(all(!sen$qc_evaluable | is.finite(sen$betaX_UKB)))
sen[, betaY_EAS_raw := fifelse(qc_evaluable & orientation_action == "FLIP",
                               -betaY_src, betaY_src)]
sen[qc_evaluable & orientation_action == "FLIP", eaf := 1 - eaf]
sen[, betaY_EAS_harmful := fifelse(qc_evaluable, harmful(betaY_EAS_raw, outcome), NA_real_)]
sen[, seY_EAS := fifelse(qc_evaluable, seY_src, NA_real_)]
sen[, betaMR_EAS := fifelse(qc_evaluable, harmful(wald(betaY_EAS_raw, betaX_UKB), outcome), NA_real_)]
sen[, seMR_EAS   := fifelse(qc_evaluable, wald_se(seY_EAS, betaX_UKB), NA_real_)]
sen[, ci95_low_EAS  := betaMR_EAS - 1.96 * seMR_EAS]
sen[, ci95_high_EAS := betaMR_EAS + 1.96 * seMR_EAS]
sen[, p_eas := fifelse(qc_evaluable, 2 * pnorm(-abs(betaMR_EAS / seMR_EAS)), NA_real_)]
sen[, q_eas := NA_real_]
for (sid in unique(sen$eas_source)) {
  fs <- which(sen$eas_source == sid & sen$qc_evaluable == TRUE)
  if (length(fs)) sen$q_eas[fs] <- p.adjust(sen$p_eas[fs], method = "BH")
}
fwrite(sen, file.path(EAS_DIR, "eas_outcome_estimates_sensitivity.csv"))
cat(sprintf("03 sensitivity: %d rows across %s; per-source BH families of %s\n",
            nrow(sen), paste(unique(sen$eas_source), collapse = "/"),
            paste(sen[qc_evaluable == TRUE, .N, by = eas_source]$N, collapse = "/")))
