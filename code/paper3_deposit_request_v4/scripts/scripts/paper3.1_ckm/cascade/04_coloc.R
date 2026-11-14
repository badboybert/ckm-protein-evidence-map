# ============================================================================
# A1: Corrected coloc (non-destructive). Reuses R/coloc.R functions but reads
# the CORRECTED suggestive set and writes to results/repro_corrected/ -- live
# results/phase3_coloc is NOT touched. eGFR (rsid_only) is skipped here and run
# by the separate eGFR coloc step. Outcome files are already R10/Wuttke (= the
# corrected versions), so coloc methodology is unchanged; only the input pair
# set expands (corrected suggestive ~1282 vs old ~842).
# Run: Rscript scripts/paper3.1_ckm/cascade/04_coloc.R
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
COLOC_SOURCE_ONLY <- TRUE                       # define functions, do NOT run 03's main()
# 2026-06-30 gate fix: force the 1682 (previous generation) instrument set BEFORE config is sourced, so R/coloc.R reads
# the pinned Phase 1 set, R/_phase1_set.R (NOT the phase1_instruments default). config uses if(!exists(...)),
# so setting these here wins. Makes a direct `Rscript 05...` foolproof against a 1635 stale-read.
PHASE1_DIR       <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), PHASE1_SET)
INSTRUMENTS_FILE <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"))

REPRO <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
dir.create(REPRO, showWarnings = FALSE, recursive = TRUE)

# redirect INPUT: corrected suggestive pairs
load_phase2_suggestive <- function() {
  sug <- file.path(REPRO, "mr_suggestive.rds")
  stopifnot(file.exists(sug))
  pairs <- readRDS(sug); setDT(pairs)
  # 2026-06-30 gate fix #3: eGFR is colocked SEPARATELY by 06 (rsid hg37 Wuttke). Since R/coloc.R
  # now has a Wuttke eGFR entry, main() would otherwise coloc eGFR here too -> double-coloc -> 08
  # max-PP.H4 dedup could pick the wrong row. Exclude eGFR here; 06 owns it.
  n0 <- nrow(pairs)
  # %in%, not !=: `NA != "eGFR"` is NA, and a NA row is dropped by `[`. A pair whose outcome_name is
  # missing or stale was therefore removed from the coloc universe silently -- it appears in neither
  # the results nor the status file, so nothing downstream can tell it was ever expected.
  # STRATIFIED BY outcome_id ALONE. Excluding on "name OR id says eGFR" while 05 selects on id only
  # left a gap: a row whose id is CAD/T2D but whose outcome_name reads "eGFR" was excluded from the
  # non-eGFR stratum by its name and from the eGFR stratum by its id, so it appeared in neither, and
  # 06 -- which reconstructs the universe the same way -- could not see that it was ever expected.
  # The id is canonical and the name is derived from it, so id-only makes the two strata exact
  # complements by construction. Deleting the name clause is the whole fix.
  pairs <- pairs[!(outcome_id %in% unname(OUTCOMES[["eGFR"]]))]   # Fix (audit Low): single-source the eGFR id from config (was hardcoded "Wuttke2019_eGFR_EA"), matching 21_coloc_split
  cat(sprintf("[corrected] suggestive pairs: %d (excluded %d eGFR -> handled by 06; %d unique proteins)\n",
              nrow(pairs), n0 - nrow(pairs), uniqueN(pairs$gene)))
  pairs
}

# redirect OUTPUT: coloc_results_v7_stage1.rds + checkpoint go to repro_corrected
PHASE3_DIR <- REPRO

cat("=== A1: CORRECTED coloc -> results/repro_corrected/ (live results/phase3_coloc untouched) ===\n")

# 2026-07-08 (multi-core is now the DEFAULT path — PI option 乙): coloc is per-locus INDEPENDENT, so by default
# this step parallelises across NCHUNK workers. 21_coloc_split partitions the EXACT same non-eGFR suggestive
# universe (mirrors load_phase2_suggestive above) into deterministic chunks; each 22_coloc_chunk_worker runs the
# SAME reviewed R/coloc.R main() on its chunk (union of chunks == one serial run); 23_merge_parallel_coloc
# unions them and ASSERTS exact coverage vs the split manifest (aborts on any missing/duplicate/errored chunk).
# COLOC_SERIAL=1 forces the single-process fallback (kept for debugging / an exact-serial audit).
.LOC <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize")
if (Sys.getenv("COLOC_SERIAL") == "1") {
  cat("[coloc] COLOC_SERIAL=1 -> single-process fallback (main())\n")
  .ckf <- file.path(REPRO, "coloc_checkpoint_v7.rds"); unlink(.ckf, force = TRUE)   # Fix (audit High): fresh serial run (no stale-checkpoint resume) so COLOC_SERIAL=1 is a CLEAN comparator
  if (file.exists(.ckf)) stop("04 coloc: could NOT remove stale coloc_checkpoint_v7.rds (locked/read-only) -- serial would silently RESUME it; unlink is not enough")   # verify the clear (unlink can no-op on a locked Windows file, per 21_split's guard)
  .stale5s <- c(file.path(REPRO, "coloc_results_v7_stage1.rds"), file.path(REPRO, "coloc_pair_status_v7.csv"))   # Fix (audit Med): serial must ALSO clear stale root v5/status (not just the checkpoint), like the parallel branch
  unlink(.stale5s, force = TRUE)
  if (any(file.exists(.stale5s))) stop("04 coloc (serial): could NOT remove stale root coloc v5/status (locked/read-only)")
  # THE SERIAL PATH MUST STAMP TOO. 23_merge attaches coloc_inputs_key on the parallel path and 06
  # refuses any result that does not name the inputs it was computed from -- so without this, the
  # documented COLOC_SERIAL=1 fallback produced a result that could never be published, which turns a
  # debugging aid into a dead end. main() cannot do it itself: R/_coloc_inputs_key.R sources
  # R/coloc.R, so calling it from inside the engine would recurse.
  # ...AND IT IS TAKEN BEFORE THE RUN, NOT ONLY AFTER. Computed once after main() returned, the key
  # named whatever bytes were on disk at that moment: replace the engine, an outcome GWAS, the tar,
  # the rsID map or the chain while main() is running, and the results carry a stamp for inputs they
  # were not computed from -- 06 recomputes the same new key and accepts them. Taken on both sides
  # and required to agree, that swap aborts instead of being invisible.
  source(file.path(BASE_DIR, "R/_coloc_inputs_key.R"))
  # load_phase2_suggestive() -- the SAME function main() called, not a second reading of the file.
  # `pairs` is local to it, so re-calling it is what guarantees the key describes the universe that
  # was actually iterated rather than one reconstructed with a slightly different filter.
  .pp <- as.data.table(load_phase2_suggestive())
  .sp <- unique(.pp[method %chin% c("Inverse variance weighted", "Wald ratio"), .(gene, outcome_id)])
  .kscope <- list(genes = sort(unique(.sp$gene)), outcome_ids = sort(unique(.sp$outcome_id)))
  .k_pre  <- coloc_inputs_key(.kscope)$key
  coloc_results <- main()
  .k_post <- coloc_inputs_key(.kscope)$key
  if (!identical(.k_pre, .k_post))
    stop(sprintf(paste0("04 coloc (serial): the coloc inputs CHANGED while main() was running ",
                        "(%s before, %s after). These results were computed from a mixture of both ",
                        "and cannot be stamped with either -- re-run on a stable tree."),
                 substr(.k_pre, 1, 12), substr(.k_post, 1, 12)))
  .v5f <- file.path(REPRO, "coloc_results_v7_stage1.rds")
  if (file.exists(.v5f)) {
    .v5 <- readRDS(.v5f)
    attr(.v5, "coloc_inputs_key") <- .k_pre
    saveRDS(.v5, .v5f)
    cat(sprintf("[coloc] serial result stamped with inputs key %s\n",
                substr(attr(.v5, "coloc_inputs_key"), 1, 12)))
  }
} else {
  cat("[coloc] DEFAULT multi-core path: 21 split -> N workers -> 23 merge\n")
  # THE ROOT OUTPUTS ARE NOT CLEARED HERE. They used to be unlinked before 21_coloc_split.R ran --
  # which is where the coordinate certificate and the input fingerprint are validated -- so a stale
  # or mismatched certificate destroyed the previous generation's only reproducible result and THEN
  # aborted. That abort condition was live in this tree when the deletion was removed, which is how
  # it was found; 03b has since been re-run and the certificate now matches, so it is no longer met
  # -- the deletion stays removed because the ordering is wrong in general, not because of one
  # transient state. 21 already moves any prior canonical root output into its own _superseded_<run>
  # directory after its validation passes, and the worker protocol does not require the root to be
  # empty while validating. Nothing needs clearing here.
  source(file.path(.LOC, "21_coloc_split.R"))                       # -> par_chunks/chunk_*.rds + _split_manifest.rds
  .PC   <- file.path(REPRO, "par_chunks")
  .man  <- readRDS(file.path(.PC, "_split_manifest.rds")); .NCH <- .man$nchunk
  .Rbin <- file.path(R.home("bin"), "Rscript.exe")
  .ld   <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/reports/Paper3_Program/run_logs/p3_coloc_workers")
  dir.create(.ld, showWarnings = FALSE, recursive = TRUE)
  # LAUNCHED DIRECTLY, NOT THROUGH A CLUSTER. makeCluster(N) + parLapply used to start N PSOCK nodes
  # whose entire job was to call system2() and block on it, so the machine ran 1 driver + N idle R
  # sessions + N workers -- 13 R processes to do the work of 6. Measured 2026-08-02: that layout put
  # the commit charge at 96% of the limit and every one of the 1,380 pairs came back
  # "region extraction python failed (exit 322)", because Windows could no longer create the python
  # process. The same six workers launched directly sat at 76-83% and ran clean. The nodes bought
  # nothing: system2(wait = FALSE) starts a worker without holding an R session open to watch it.
  #
  # Completion is read from what the worker itself writes, which is a stronger statement than an exit
  # code: _worker_complete.rds is its LAST action, after the pair loop and after its own abort checks,
  # so its presence means the worker reached the end rather than merely stopped.
  # EACH WORKER REPORTS ITS OWN EXIT CODE. system2(wait = FALSE) returns no pid on Windows, so the
  # first two versions of this loop inferred liveness -- from a marker file, then from the machine's
  # process table -- and both were weaker than the exit-code check they replaced: a marker cannot see
  # a worker that died before writing it, and a machine-wide Rscript count cannot tell this run's
  # workers from anyone else's, so an unrelated R session would make the driver wait forever. A
  # one-line batch file removes the inference entirely: cmd runs the worker, waits for it, and writes
  # its exit code where the driver can read it. That is exactly what parLapply returned, at the cost
  # of one cmd.exe per chunk instead of one idle R session per chunk.
  .exitf <- function(k) file.path(.PC, sprintf("out_%d", k), "_worker_exit.txt")
  for (k in seq_len(.NCH)) {
    unlink(.exitf(k))                                   # a stale code from a previous run is not this run's
    .bat <- file.path(.PC, sprintf("_run_worker_%d.bat", k))
    writeLines(c("@echo off",
                 sprintf('"%s" "%s" "%s" out_%d > "%s" 2> "%s"',
                         gsub("/", "\\\\", .Rbin),
                         gsub("/", "\\\\", file.path(.LOC, "22_coloc_chunk_worker.R")),
                         gsub("/", "\\\\", file.path(.PC, sprintf("chunk_%d.rds", k))), k,
                         gsub("/", "\\\\", file.path(.ld, sprintf("worker_%d.out.log", k))),
                         gsub("/", "\\\\", file.path(.ld, sprintf("worker_%d.err.log", k)))),
                 sprintf('echo %%ERRORLEVEL%% > "%s"', gsub("/", "\\\\", .exitf(k)))),
               .bat)
    system2("cmd", c("/c", shQuote(gsub("/", "\\\\", .bat))), wait = FALSE)
  }
  cat(sprintf("[coloc] %d workers launched -> waiting\n", .NCH))
  # WAIT FOR EVERY EXIT CODE, THEN JUDGE. The exit file is written by cmd after the worker has
  # terminated, so its presence means that worker is finished -- there is nothing to infer and no
  # window in which a half-written artifact reads as success. `_worker_complete.rds` is still
  # required afterwards, because a worker can exit 0 on a path that never reaches it.
  repeat {
    .have <- vapply(seq_len(.NCH), function(k) file.exists(.exitf(k)), logical(1))
    if (all(.have)) break
    Sys.sleep(15)
  }
  .codes <- vapply(seq_len(.NCH), function(k) {
    .v <- suppressWarnings(as.integer(trimws(readLines(.exitf(k), warn = FALSE)[1])))
    if (is.na(.v)) NA_integer_ else .v
  }, integer(1))
  .bad <- which(is.na(.codes) | .codes != 0L)
  if (length(.bad) > 0L)
    stop(sprintf(paste0("coloc worker(s) FAILED: chunk(s) %s exited %s -- see %s/worker_*.err.log; ",
                        "NOT merging a partial coloc"),
                 paste(.bad, collapse = ","), paste(.codes[.bad], collapse = ","), .ld))
  .nostamp <- which(!vapply(seq_len(.NCH), function(k)
    file.exists(file.path(.PC, sprintf("out_%d", k), "_worker_complete.rds")), logical(1)))
  if (length(.nostamp) > 0L)
    stop(sprintf(paste0("coloc worker(s) exited 0 but chunk(s) %s never wrote _worker_complete.rds ",
                        "-- they stopped before finishing; NOT merging a partial coloc"),
                 paste(.nostamp, collapse = ",")))
  unlink(file.path(.PC, sprintf("_run_worker_%d.bat", seq_len(.NCH))))
  cat(sprintf("[coloc] all %d workers exited 0 -> merging + coverage assert\n", .NCH))
  source(file.path(.LOC, "23_merge_parallel_coloc.R"))              # writes coloc_results_v7_stage1.rds; aborts unless FULL coverage
  coloc_results <- readRDS(file.path(REPRO, "coloc_results_v7_stage1.rds"))
}
cat("=== A1 DONE ===\n")
