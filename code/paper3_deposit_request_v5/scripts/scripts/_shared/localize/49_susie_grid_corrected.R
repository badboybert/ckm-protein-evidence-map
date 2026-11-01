# A6 backfill — SuSiE 4-setting grid (susie_grid) on the CORRECTED cascade, fully OFFLINE.
# Reuses canonical phase3/16_susie_grid.R via guarded overrides (offline JWT skip, local
# plink LD + local coords from the localized R/coloc.R). Writes to the A6 shadow tree.
# Smoke:  SG_SMOKE_RUN=1 Rscript .../49_susie_grid_corrected.R   (1 non-eGFR + 1 eGFR pair)
# Chunk:  SG_CHUNK_N=N SG_CHUNK_I=i ...   then SG_MERGE=1 ...
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
INSTRUMENTS_FILE <- PHASE1_INSTRUMENTS  # -> N_PROTEINS=1682 (previous generation)
PHASE1_DIR       <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), PHASE1_SET)
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
RC <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
SG_COLOC     <- file.path(RC, "coloc_results_v7.rds")          # corrected 52 strong

# THE APPROVED ARTIFACT IS PINNED HERE, IN THE GENERATION-SPECIFIC 49, not in the engine.
# The engine can check that the two halves of its inputs agree with each other, but it has no
# independent anchor: a foreign coloc_results_v7.rds together with its own matching SuSiE table
# passes every internal comparison. The realistic trigger is an interrupted 06 leaving an earlier
# schema-compatible file in place, after which the strong-pair set and the ABF posteriors come from
# that file while the SuSiE values are computed from current inputs -- and nothing in the saved
# result says so.
#
# THE ANCHOR IS THE INPUTS, NOT THE BYTES. This was a hardcoded md5 of the approved v7 file, which
# could not survive its own contract: 06_combine_coloc.R:148 stamps ckm_combined_at = Sys.time() into
# v7, so re-running 06 changes the md5 even when every number is identical, and the literal then had
# to be hand-edited -- a frozen value of exactly the class this rebuild removes.
# 06 already binds v7 to the two artifacts it combined, by their coloc_inputs_key attributes
# (06:85-95); those keys hash the GWAS, instrument, coordinate, chain and prior bytes of the session
# that produced them. Re-checking that same equality here needs no literal and still catches what
# this anchor exists for: an interrupted 06 leaving an earlier schema-compatible v7, whose stored
# keys then disagree with the upstream artifacts on disk.
.SG_V5_F  <- file.path(RC, "coloc_results_v7_stage1.rds")
.SG_EGF_F <- file.path(RC, "coloc_egfr_results.rds")
for (.f in c(SG_COLOC, .SG_V5_F, .SG_EGF_F))
  if (!file.exists(.f)) stop(sprintf("49: %s is missing -- run 06_combine_coloc.R for this generation first.", .f))
.SG_ATT <- function(f, a) { z <- attr(readRDS(f), a, exact = TRUE); if (is.null(z)) NA_character_ else as.character(z)[1] }
.SG_SEEN5 <- .SG_ATT(SG_COLOC, "coloc_inputs_key_noneGFR"); .SG_WANT5 <- .SG_ATT(.SG_V5_F,  "coloc_inputs_key")
.SG_SEENE <- .SG_ATT(SG_COLOC, "coloc_inputs_key_eGFR");    .SG_WANTE <- .SG_ATT(.SG_EGF_F, "coloc_inputs_key")
if (anyNA(c(.SG_SEEN5, .SG_WANT5, .SG_SEENE, .SG_WANTE)))
  stop("49: one of the coloc artifacts carries no coloc_inputs_key, so it cannot be shown to belong to this generation. Rerun 06_combine_coloc.R.")
# BOTH SIDES CARRY THE SAME STRING. 06 forms the eGFR key as "<md5 of 05_coloc_egfr.R>|<key>" and 05
# stamps its own output with that same composite, so these compare whole -- measured on disk, not
# assumed. The first attempt here stripped the prefix off one side only and could never match; a
# smoke test against the real artifacts caught it. The prefix is still checked separately, against 05
# as it is on disk now, which is what makes a changed eGFR computation visible.
.SG_EKEY <- sub("^.*[|]", "", .SG_SEENE); .SG_ESCR <- sub("[|].*$", "", .SG_SEENE)
.SG_05   <- unname(tools::md5sum(file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/05_coloc_egfr.R")))
if (!identical(.SG_SEEN5, .SG_WANT5) || !identical(.SG_SEENE, .SG_WANTE) ||
    is.na(.SG_05) || !identical(tolower(.SG_ESCR), tolower(.SG_05)))
  stop(sprintf("49: %s was not built from the coloc artifacts on disk now (non-eGFR %s vs %s; eGFR %s vs %s; 05 %s vs %s). Either 06_combine_coloc.R has not finished for this generation, or an earlier generation's file is in place. Re-run 04 -> 05 -> 06.",
               SG_COLOC, substr(.SG_SEEN5, 1, 12), substr(.SG_WANT5, 1, 12),
               substr(.SG_EKEY, 1, 12), substr(.SG_WANTE, 1, 12),
               substr(.SG_ESCR, 1, 12), substr(.SG_05, 1, 12)))
cat(sprintf("[anchor] coloc artifact bound to its inputs: non-eGFR %s | eGFR %s
", substr(.SG_SEEN5, 1, 12), substr(.SG_EKEY, 1, 12)))
SG_INSTR     <- INSTRUMENTS_FILE
SG_OUT       <- file.path(RC, "_fig_results", "phase3_susie")
SG_KEEP_EGFR <- FALSE   # non-eGFR grid (matches live susie_grid scope; 16's offline prep doesn't handle
                        # the Wuttke eGFR id — eGFR SuSiE is covered single-setting in susie_coloc_results)
if (nzchar(Sys.getenv("SG_SMOKE_RUN"))) SG_SMOKE <- TRUE
dir.create(SG_OUT, showWarnings = FALSE, recursive = TRUE)
cat("=== A6 SuSiE grid (corrected, OFFLINE) -> ", SG_OUT, " ===\n")
# Codex backfill-review fix: offline grid must use LOCAL LD only -- assert availability before sourcing
if (!exists("PLINK_BIN") || !file.exists(PLINK_BIN)) stop("49: PLINK_BIN not found for offline LD: ", if (exists("PLINK_BIN")) PLINK_BIN else "<unset>")
if (!exists("LD_EUR_DIR") || !file.exists(file.path(LD_EUR_DIR, "EUR.bim"))) stop("49: local 1KG EUR LD bfile (EUR.bim) not found under LD_EUR_DIR")
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase3/16_susie_grid.R"))
