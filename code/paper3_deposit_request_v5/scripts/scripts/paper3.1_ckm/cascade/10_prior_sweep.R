# ============================================================================
# A4 input: prior sweep (coloc p12 sensitivity) for the NON-eGFR strong hits on the
# corrected v7 cascade. Thin wrapper: sets the corrected inputs + offline flag,
# then sources the (now parametrized) phase3/10. eGFR pairs are handled by cascade/11.
# Smoke test:  PSW_LIMIT=2 Rscript .../10_prior_sweep.R
# Full run:    Rscript .../10_prior_sweep.R
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
Sys.setenv(CKM_OFFLINE = "1")
Sys.unsetenv(c("PSW_CHUNK_N","PSW_CHUNK_I","PSW_MERGE"))  # 2026-07-01 gate: corrected sweep is SERIAL-only — force-clear any leaked chunk/merge env (don't rely on shell hygiene)
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

PSW_COLOC        <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds")          # selected-generation v7 strong
PSW_INSTR        <- PHASE1_INSTRUMENTS     # the selected Phase 1 set (R/_phase1_set.R)
PSW_OUT          <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))                                  # corrected workspace
PSW_NONEGFR_ONLY <- TRUE
.WHO <- "10"
# EXISTENCE IS NOT CURRENCY. Unlike the deCODE, ARIC and eQTL wrappers, this one passed the combined
# coloc file by pathname alone, so an interrupted or resumed 06_combine_coloc -- an earlier strong
# set fused with current inputs -- would define the pair set this sweep is computed over, and the
# output would look like an ordinary complete table with no input key in it. The two stage
# artifacts each carry the key their own stratum was computed under; the combined file's two stamps
# are checked against them. Copied from cascade/12c_eqtl.R:35-49, same attributes, no new mechanism.
local({
  .rc <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
  .cr <- readRDS(PSW_COLOC)
  .stage <- c(coloc_inputs_key_noneGFR = file.path(.rc, "coloc_results_v7_stage1.rds"),
              coloc_inputs_key_eGFR    = file.path(.rc, "coloc_egfr_results.rds"))
  for (.a in names(.stage)) {
    .v <- attr(.cr, .a, exact = TRUE)
    if (is.null(.v) || !nzchar(as.character(.v)[1]))
      stop(sprintf("%s: %s carries no %s -- it predates the stamping or came from another generation. Re-run cascade/06_combine_coloc.R.", .WHO, basename(PSW_COLOC), .a))
    if (!file.exists(.stage[[.a]]))
      stop(sprintf("%s: cannot verify %s -- the stage artifact %s is absent.", .WHO, .a, basename(.stage[[.a]])))
    .sv <- attr(readRDS(.stage[[.a]]), "coloc_inputs_key", exact = TRUE)
    if (is.null(.sv) || !identical(as.character(.sv)[1], as.character(.v)[1]))
      stop(sprintf("%s: %s does not match %s -- the combined file was not built from the stage artifact now on disk. Re-run cascade/06_combine_coloc.R.", .WHO, .a, basename(.stage[[.a]])))
  }
  cat(sprintf("[%s] coloc artifact bound to its inputs\n", .WHO))
})
# AND THE LIVE INPUTS MUST NOT POSTDATE THE COLOC TABLE. The key comparison above proves the
# combined file and the stage artifacts agreed when they were written; it does not bind the sweep to
# the files it re-reads at run time. The sweep recomputes coloc.abf from R/coloc.R, the shared and
# generation-untagged cache/cis_rsid_map.tsv, and the outcome extracts. If any of those is replaced
# after the coloc table was built -- and cis_rsid_map.tsv was rebuilt on 2026-08-07, so this is a
# real operation here -- the stored keys still match, the sweep still completes, and the posterior
# movement it reports as prior sensitivity would partly be data or code movement. Nothing in
# prior_sweep_wide records which inputs produced it, so a reader cannot tell.
# Same test as the merge-base check: an input newer than the table it must be consistent with stops
# the run rather than silently recomputing under different data.
local({
  # THE FULL SET THE SWEEP RE-READS, not the two that came to mind. These are the dependency
  # classes R/_coloc_inputs_key.R already enumerates for exactly this purpose: the engine, the
  # config, the Phase-1 instrument table, the certified GRCh38 coordinates, the rsID map, the
  # downward liftover chain, the pQTL archives and the outcome GWAS. Naming a subset made the check
  # look like coverage while leaving the rest open -- rebuilding gene_hg38_coords.tsv or replacing a
  # same-schema outcome file would have passed it.
  # Compared by mtime, not by recomputing coloc_inputs_key(): the stored keys were computed under a
  # scope this wrapper does not have, so a live key would differ from them for a correct setup
  # (measured: 3a802cdd vs the stored 84be4820) and would abort every run.
  .live <- c("R/coloc.R"              = file.path(BASE_DIR, "R", "coloc.R"),
             "R/00_config.R"          = file.path(BASE_DIR, "R", "00_config.R"),
             "cache/cis_rsid_map.tsv" = file.path(BASE_DIR, "cache", "cis_rsid_map.tsv"),
             "phase1 instruments"     = PSW_INSTR,
             "gene_hg38_coords.tsv"   = if (exists("GENE_HG38_COORDS")) GENE_HG38_COORDS else NA_character_,
             "hg38ToHg19 chain"       = if (exists("LIFTOVER_CHAIN_38_TO_37")) LIFTOVER_CHAIN_38_TO_37 else NA_character_)
  .live <- .live[!is.na(.live)]
  # directories: take the newest member rather than the directory's own mtime, which does not move
  # when a file inside it is rewritten in place
  # NAMES SURVIVE ONLY IF THE LOOP RUNS OVER THEM. `for (x in c(a=1,b=2))` iterates over the VALUES
  # and drops the names, so names(x) is NULL inside the body and the named append below produced an
  # unnamed element that later subscripting could not find.
  .dirs <- c("pQTL archives" = file.path(DATA_DIR, "ukb_ppp", "european"),
             "outcome GWAS"  = OUTCOME_DIR)
  for (.n in names(.dirs)) {
    # top-level files only: OUTCOME_DIR's subdirectories hold localized extracts that the coloc
    # engine does not read, and stat-ing those instead of the full GWAS is what made the first
    # version of this guard watch the wrong files (objective review 2026-08-08, on 29).
    .fs <- list.files(.dirs[[.n]], full.names = TRUE)
    .fs <- .fs[!dir.exists(.fs)]
    if (length(.fs)) .live <- c(.live, setNames(.fs[which.max(file.mtime(.fs))], .n))
  }
  .live <- .live[file.exists(.live)]
  # THE REFERENCE IS THE EARLIEST STAGE, NOT THE COMBINED FILE. The combined table is written last
  # (03:19 on disk) while the strata it merges finished at 01:15 and 03:18, so comparing against it
  # accepted a dependency replaced after a stratum was computed but before the merge -- precisely
  # the window that matters, since the sweep recomputes from those inputs.
  .stages <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG),
                       c("coloc_results_v7_stage1.rds", "coloc_egfr_results.rds"))
  .stages <- .stages[file.exists(.stages)]
  .ct <- min(c(file.mtime(PSW_COLOC), file.mtime(.stages)))
  .KEY_ARM <- "coloc_inputs_key_noneGFR"
  .newer <- names(.live)[file.mtime(.live) > .ct]
  # MTIME IS THE SCREEN, CONTENT IS THE VERDICT (L8r19 #1). A newer mtime is not evidence that an
  # input CHANGED: editing a file and putting the original bytes back, or rewriting it with the same
  # content, moves the timestamp and tripped this stop on a generation whose inputs are provably
  # identical. When the screen trips, recompute the colocalisation inputs key over the non-eGFR
  # scope -- the same derivation 47_a6_stage_fig_inputs.R uses, from mr_suggestive -- and compare it
  # with the key the coloc artifact was stamped with. Equal keys mean the bytes this sweep will read
  # are the bytes the coloc table was computed from, whatever the timestamps say. Different keys
  # abort exactly as before, and now name both keys.
  # (The older comment above said the wrapper "does not have the scope" for a live key. It does:
  # the scope is derived from mr_suggestive, which is what made the measured mismatch disappear.)
  if (length(.newer)) {
    .verdict <- tryCatch({
      source(file.path(BASE_DIR, "R", "_coloc_inputs_key.R"))
      .sug <- data.table::as.data.table(readRDS(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG),
                                                          "mr_suggestive.rds")))
      .pr  <- .sug[method %in% c("Inverse variance weighted", "Wald ratio")]
      .eg  <- unname(OUTCOMES[["eGFR"]])
      .sel <- if (.KEY_ARM == "coloc_inputs_key_eGFR") .pr[outcome_id %in% .eg] else .pr[!(outcome_id %in% .eg)]
      .now <- coloc_inputs_key(list(genes = sort(unique(.sel$gene)),
                                    outcome_ids = sort(unique(.sel$outcome_id))))$key
      # THE eGFR ARM IS A COMPOSITE KEY, NOT THE BARE ONE (L8r20 #1). cascade/06 stamps it as
      # md5(05_coloc_egfr.R) | shared-key (06:87), and 47_a6_stage_fig_inputs.R compares the same
      # composite. Comparing the bare key against the composite can never be true, so this branch
      # blocked every run of 11 instead of clearing it.
      if (.KEY_ARM == "coloc_inputs_key_eGFR")
        .now <- paste(unname(tools::md5sum(file.path(BASE_DIR,
                        "scripts/paper3.1_ckm/cascade/05_coloc_egfr.R"))), .now, sep = "|")
      .stamp <- as.character(attr(readRDS(PSW_COLOC), .KEY_ARM, exact = TRUE))[1]
      list(now = .now, stamp = .stamp, ok = identical(.now, .stamp))
    }, error = function(e) list(now = NA_character_, stamp = NA_character_, ok = NA, err = conditionMessage(e)))
    if (isTRUE(.verdict$ok))
      cat(sprintf("[%s] %s postdate(s) the coloc table by mtime, but the inputs key still recomputes to %s -- content unchanged, proceeding\n",
                  .WHO, paste(.newer, collapse = ", "), substr(.verdict$stamp, 1, 12)))
    else
      stop(sprintf("%s: %s postdate(s) %s (%s) AND the inputs key no longer matches (live %s vs stamped %s%s). The sweep re-reads them, so its posterior movement would not be prior-only. Re-run cascade/06_combine_coloc.R before sweeping.",
                   .WHO, paste(.newer, collapse = ", "), basename(PSW_COLOC),
                   format(.ct, "%Y-%m-%d %H:%M"),
                   substr(ifelse(is.na(.verdict$now), "<uncomputable>", .verdict$now), 1, 12),
                   substr(ifelse(is.na(.verdict$stamp), "<absent>", .verdict$stamp), 1, 12),
                   if (!is.null(.verdict$err)) paste0("; ", .verdict$err) else ""))
  }
  cat(sprintf("[%s] live sweep inputs predate the coloc table\n", .WHO))
})

.lim <- suppressWarnings(as.integer(Sys.getenv("PSW_LIMIT", "")))
if (!is.na(.lim) && .lim > 0) {
  PSW_LIMIT <- .lim   # smoke-test cap
  # AND IT WRITES SOMEWHERE ELSE. The engine performs its canonical saves regardless of the cap, so
  # the documented PSW_LIMIT=2 invocation used to put two wide rows and six long rows into the
  # generation directory with nothing marking them partial; a later full cascade/11 would then merge
  # them into an ordinary-looking 18-row table. A capped run is a smoke test, so it goes to a smoke
  # directory and cannot reach production at all.
  PSW_OUT <- file.path(PSW_OUT, sprintf("_smoke_PSW_LIMIT_%d", .lim))
  dir.create(PSW_OUT, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("[%s] PSW_LIMIT=%d -> SMOKE run, output redirected to %s\n", .WHO, .lim, PSW_OUT))
}

cat(sprintf("=== A4 prior-sweep (non-eGFR) on %s ===\n", PHASE1_TAG))
# L8r17 #1: 10 refuses to start when the target table already holds eGFR rows -- that guard exists to
# stop an out-of-order 11-before-10 run, but it cannot tell that case apart from a COMPLETED 10 -> 11
# sequence, so re-running this generation as documented aborted before computing anything. This
# wrapper IS the documented entry point for a fresh non-eGFR sweep, and the 10 -> 11 order is its
# contract, so it clears the pair of outputs it is about to rebuild. A direct invocation of phase3/10
# still meets the guard untouched.
local({
  .psw <- if (exists("PSW_OUT")) PSW_OUT else file.path(RESULTS_DIR, "phase3_sensitivity")
  for (.f in file.path(.psw, c("prior_sweep_wide.rds", "prior_sweep_wide.csv",
                               "prior_sweep_long.rds"))) {
    if (file.exists(.f)) {
      if (!file.remove(.f))                      # never proceed on a failed delete: 10 would then
        stop(sprintf("cascade/10: could not remove %s before the sweep", .f))   # abort or, worse, merge into a stale table
      cat(sprintf("  [clear] %s (this wrapper rebuilds it; 10 -> 11 order is its contract)\n", basename(.f)))
    }
  }
})
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/paper3.1_ckm/phase3/10_phase3s_prior_sweep.R"))
