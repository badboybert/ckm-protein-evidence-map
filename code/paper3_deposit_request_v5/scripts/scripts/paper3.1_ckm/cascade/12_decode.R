# ============================================================================
# A3 (deCODE replication, localized): re-run deCODE cis-MR replication of the
# CORRECTED v7 strong-coloc hits, fully OFFLINE. Sets DECODE_LOCAL + the
# corrected inputs, then sources the now-localized decode/01 (outcome side via
# local_extract_outcome_data instead of OpenGWAS extract_outcome_data). deCODE
# pQTL (exposure) is already local (data/decode/pqtl). Genes without a deCODE
# pQTL file -> NO_DECODE_DATA (real data gap, not API).
# Run: Rscript scripts/paper3.1_ckm/cascade/12_decode.R
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
Sys.setenv(CKM_OFFLINE = "1")                      # blank tokens -> any stray network call fails loud
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

DECODE_LOCAL   <- TRUE
DECODE_COLOC   <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds")          # selected-generation v7 strong hits
DECODE_INSTR   <- PHASE1_INSTRUMENTS     # the selected Phase 1 set (R/_phase1_set.R)
DECODE_UKBMR   <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "mr_all_results.rds")            # corrected UKB MR (direction compare)
DECODE_OUT_DIR <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "decode_replication")
stopifnot(file.exists(DECODE_COLOC), file.exists(DECODE_INSTR), file.exists(DECODE_UKBMR))

cat(sprintf("=== A3 deCODE replication (localized, OFFLINE) on %s ===\n", PHASE1_TAG))
cat("  coloc :", DECODE_COLOC, "\n  instr :", DECODE_INSTR, "\n  ukbMR :", DECODE_UKBMR, "\n  out   :", DECODE_OUT_DIR, "\n\n")

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/decode/01_decode_replication_mr.R"))
cat("=== A3 deCODE DONE ===\n")
