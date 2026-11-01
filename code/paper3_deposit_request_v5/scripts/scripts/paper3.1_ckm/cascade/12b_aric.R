# ============================================================================
# ARIC cross-platform replication of the cis-pQTL effect — thin wrapper.
#
# Parallel arm to cascade/12_decode.R. Neither is primary and the two are never
# pooled: they are reported side by side, per pair, with the source named. See
# _shared/aric/01_aric_replication_mr.R for what this arm does and does not license.
#
# Run: Rscript scripts/paper3.1_ckm/cascade/12b_aric.R
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

ARIC_LOCAL <- TRUE                                  # outcomes come from the local extracts, not an API

ARIC_DIR     <- file.path(BASE_DIR, "data/aric")
ARIC_ZIP     <- file.path(ARIC_DIR, "EA.zip")
ARIC_SEQID   <- file.path(ARIC_DIR, "seqid.txt")
ARIC_INSTR   <- PHASE1_INSTRUMENTS                                                   # the selected Phase 1 set
ARIC_COLOC   <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds")
ARIC_OUT_DIR <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "aric_replication")

stopifnot(file.exists(ARIC_ZIP), file.exists(ARIC_SEQID),
          file.exists(ARIC_INSTR), file.exists(ARIC_COLOC))

# The strong set this arm replicates must be the one THIS generation produced. 06
# binds coloc_results_v7.rds to the two artifacts it combined via their
# coloc_inputs_key attributes; checking that equality here is what stops an
# interrupted 06 leaving an earlier, schema-compatible v7 in place while the ARIC
# estimates are computed from current inputs, with nothing in the output saying so.
.A_ATT <- function(f, a) { z <- attr(readRDS(f), a, exact = TRUE); if (is.null(z)) NA_character_ else as.character(z)[1] }
.A_S5 <- .A_ATT(ARIC_COLOC, "coloc_inputs_key_noneGFR")
.A_W5 <- .A_ATT(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7_stage1.rds"), "coloc_inputs_key")
.A_SE <- .A_ATT(ARIC_COLOC, "coloc_inputs_key_eGFR")
.A_WE <- .A_ATT(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_egfr_results.rds"), "coloc_inputs_key")
if (anyNA(c(.A_S5, .A_W5, .A_SE, .A_WE)))
  stop("12b_aric: a coloc artifact carries no coloc_inputs_key, so it cannot be shown to belong to this generation. Re-run 06_combine_coloc.R.")
if (!identical(.A_S5, .A_W5) || !identical(.A_SE, .A_WE))
  stop(sprintf("12b_aric: %s was not built from the coloc artifacts on disk now (non-eGFR %s vs %s; eGFR %s vs %s). Re-run 04 -> 05 -> 06.",
               ARIC_COLOC, substr(.A_S5, 1, 12), substr(.A_W5, 1, 12), substr(.A_SE, 1, 12), substr(.A_WE, 1, 12)))

cat(sprintf("=== ARIC replication on %s ===\n", PHASE1_TAG))
cat("  zip   :", ARIC_ZIP, "\n  instr :", ARIC_INSTR, "\n  coloc :", ARIC_COLOC, "\n  out   :", ARIC_OUT_DIR, "\n\n")

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/aric/01_aric_replication_mr.R"))
cat("=== ARIC replication DONE ===\n")
