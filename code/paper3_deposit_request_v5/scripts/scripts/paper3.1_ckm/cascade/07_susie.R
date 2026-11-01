# ============================================================================
# A2 (SuSiE, localized): run SuSiE fine-mapping coloc on the CORRECTED v7
# strong hits, fully OFFLINE. Thin wrapper that sets the corrected inputs then
# sources the now-localized phase3/15_susie_coloc.R (which no longer calls any
# API: local plink LD, local gene coords via the fixed R/coloc.R, no JWT
# preflight, no ld_matrix(pop="EUR") fallback).
# Run: Rscript scripts/paper3.1_ckm/cascade/07_susie.R
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
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))    # defines RESULTS_DIR (+ honours CKM_OFFLINE)

# corrected inputs (override phase3/15 defaults)
SUSIE_COLOC_INPUT <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds")        # selected-generation v7 strong hits
SUSIE_INSTR       <- PHASE1_INSTRUMENTS   # the selected Phase 1 set (R/_phase1_set.R)
SUSIE_OUT_DIR     <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "phase3_susie")
stopifnot(file.exists(SUSIE_COLOC_INPUT), file.exists(SUSIE_INSTR))

cat(sprintf("=== A2 SuSiE (localized, OFFLINE) on %s ===\n", PHASE1_TAG))
cat("  input :", SUSIE_COLOC_INPUT, "\n  instr :", SUSIE_INSTR, "\n  out   :", SUSIE_OUT_DIR, "\n\n")

# THE PRE-REGISTERED 10000 REFIT WAS ATTEMPTED AND THE SOFTWARE REFUSED IT (REBUILD_LEDGER §9ba
# registered it; §9bb registered the three outcomes before any result was seen; §9bf records what
# landed). At SUSIE_MAX_ITER = 10000 this trunk reached 30 of 60 pairs and then susie_rss threw on
# NCAN ~ CAD's exposure fit: "The estimated prior variance is unreasonably large ... mismatch
# between the summary statistics and the LD matrix." The engine's fail-loud guard aborted the whole
# run rather than swallow a data-foundation failure, which is what it is for, so the 10000 pass
# produced no final artifact. It did leave susie_checkpoint.rds -- 30 pairs at max_iter_used =
# 10000 -- which is not the canonical result and is cleared by FORCE_RERUN before the next run.
#
# WHAT THE ESCALATION DID *NOT* ESTABLISH, because the abort came first: NCAN's OUTCOME-side refit
# and ADM ~ AF's exposure refit were never reached, so their behaviour at 10000 is unknown. And
# susieR says that error is "usually caused by" summary-statistic/LD mismatch -- strong evidence
# consistent with mismatch, not proof of it. Neither may be written as a measured result.
#
# The ceiling is therefore left at the engine default of 1000 (phase3/15_susie_coloc.R:130). What
# is established is that this trunk COMPLETES at 1000 and ABORTS at 10000; nothing was measured
# about any other ceiling, and no claim of uniqueness is made. NCAN ~ CAD is not_evaluable either way -- via
# non-convergence at 1000, via the throw at 10000 -- so no classification depends on this line; the
# throw is reported as evidence about that locus, not as a result in the table.
#
# NOT A SILENT RETREAT: removing the assignment restores the engine default, it does not lower a
# ceiling below where it started. Re-raising it is the one thing §9bb forbids ("不得二次加碼").
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase3/15_susie_coloc.R"))
cat("=== A2 SuSiE DONE ===\n")
