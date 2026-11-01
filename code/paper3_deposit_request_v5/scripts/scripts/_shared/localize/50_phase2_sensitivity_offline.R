# A6 backfill — phase2 sensitivity (heterogeneity_Q + sens_5e8_concordance) on the CORRECTED
# cascade, OFFLINE. Reuses canonical phase2/11 + phase2/15 but overrides extract_outcome_data
# with the LOCAL extractor (no API). Writes to the A6 shadow tree (RESULTS_DIR=shadow).
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
RESULTS_DIR      <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), paste0("repro_", PHASE1_SET), "_fig_results")   # shadow as results root
INSTRUMENTS_FILE <- PHASE1_INSTRUMENTS  # measured 2026-08-17: 1,964 instruments / 1,948 proteins (the
                                        # comment here read "N_PROTEINS=1682", a previous generation's
                                        # count; the NAME resolves per PHASE1_SET, so only the note was stale
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
suppressMessages(library(TwoSampleMR))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_local_outcome_extract.R"))
# palindromic_fallthrough = TRUE: these sidecars must read the SAME outcome data the discovery screen
# read, or they are sensitivity analyses of an estimate that was never computed. With the flag off, an
# ambiguous palindromic DIRECT row counts as a direct match, removes itself from `missing` and so
# blocks the proxy path -- and harmonise_data then discards it, leaving the pair with fewer SNPs than
# discovery used. Measured: CNTNAP2-HF ran on 2 SNPs against discovery's 3, and its Cochran's Q went
# from p = 0.0170 (2 df) to p = 0.3179 (1 df), reversing the heterogeneity conclusion.
extract_outcome_data <- function(...) local_extract_outcome_data(..., palindromic_fallthrough = TRUE)   # global override -> offline, direct matches
P2S_STRONG_FILE <- file.path(RESULTS_DIR, "phase4_annotation/table_s_pair_evidence.csv")  # staged; measured
                                        # 2026-08-17: 76 pairs / 63 genes. The note read "corrected 52",
                                        # a previous generation's count -- the path itself is current.

cat("\n########## A6 phase2-sensitivity OFFLINE: heterogeneity_Q (11) ##########\n")
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase2/11_phase2s_sensitivity.R"))

cat("\n########## A6 phase2-sensitivity OFFLINE: sens_5e8 non-eGFR (15) ##########\n")
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase2/15_phase2_sens_5e8.R"))
cat("\n=== phase2-sensitivity offline (11+15) done ===\n")
