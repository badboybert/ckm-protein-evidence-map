# A6 backfill — eGFR p<5e-8 concordance (phase2/16) on CORRECTED cascade, OFFLINE.
# Appends the 14 corrected eGFR pairs to the 38-row non-eGFR sens_5e8_concordance (-> 52),
# completing the fig4(d) 5e-8 panel. Local extractor override; Wuttke eGFR id from config.
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
Sys.setenv(CKM_OFFLINE = "1")
# PHASE1_SET, not PHASE1_TAG: the tag is defined in 00_config.R (line 1089) as a copy of the set,
# and config is sourced BELOW this line -- so this read was of an undefined name and the script
# stopped here in any fresh process. Naming the variable that exists at this point leaves the
# deliberate ordering (RESULTS_DIR / INSTRUMENTS_FILE set before config, so config adopts them)
# untouched.
RESULTS_DIR      <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), paste0("repro_", PHASE1_SET), "_fig_results")
INSTRUMENTS_FILE <- PHASE1_INSTRUMENTS
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
suppressMessages(library(TwoSampleMR))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_local_outcome_extract.R"))
# palindromic_fallthrough = TRUE for the same reason as 50_phase2_sensitivity_offline.R:16 -- the eGFR
# 5e-8 arm must read the outcome data the discovery screen read, not a stricter subset of it.
extract_outcome_data <- function(...) local_extract_outcome_data(..., palindromic_fallthrough = TRUE)
P2S_EGFR_FILE <- file.path(RESULTS_DIR, "phase4_annotation/table_s_pair_evidence.csv")
cat("=== A6 eGFR 5e-8 OFFLINE (16) ===\n")
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase2/16_sens_5e8_egfr.R"))
cat("=== eGFR 5e-8 done ===\n")
