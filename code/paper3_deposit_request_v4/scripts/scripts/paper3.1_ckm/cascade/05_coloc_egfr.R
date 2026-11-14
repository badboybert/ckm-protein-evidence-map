# Phase 3-E: Re-run eGFR Coloc with Corrected Suggestive Pairs
# ============================================================================
# P3-F3 (cold audit) — corrected header: re-runs eGFR coloc on the CURRENT eGFR suggestive-pair set
# (Wuttke CKDGen, corrected cascade; count is data-driven, read at run time). Writes coloc_egfr_results.rds +
# chunk files; the FINAL combination into coloc_results_v7/final is done by 06_combine_coloc.R, NOT this script.
#
# Uses the LOCAL CKDGen file. Run via the corrected cascade (cascade/05), not standalone.
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
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
# liftover_to_hg37() lives in R/coloc.R. This script carries its own copy of the region/merge logic,
# but the GRCh38 -> GRCh37 projection must be the SAME function in both places: two copies of a
# coordinate conversion is exactly how the two branches came to disagree about which build the pQTL
# positions are in.
COLOC_SOURCE_ONLY <- TRUE   # define functions only. R/coloc.R:1213 calls main() on source otherwise,
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"))   # so this script ran an entire generic coloc before its own.
rm(COLOC_SOURCE_ONLY)

cat("=== Phase 3-E: Re-run eGFR Coloc (Corrected Pairs) ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
flush(stdout())

# ── Settings ─────────────────────────────────────────────────────────────────
EGFR_ID    <- unname(OUTCOMES[["eGFR"]])   # Fix (audit Med): single-source from config (was hardcoded "Wuttke2019_eGFR_EA"), matching 04/21; used throughout 05
# THE FILENAME COMES FROM THE SAME ENTRY AS THE SCALE. This spelled the filename out independently of
# GWAS_META, which holds both the filename and the phenotype SD estimated FROM that file. Two
# definitions of which bytes are the eGFR outcome meant one could be pointed at a different release
# while the other kept describing the old one -- and the mismatch that produces is a silently wrong
# posterior, not an error. One definition, so they cannot disagree.
EGFR_FILE  <- file.path(OUTCOME_DIR, GWAS_META[[EGFR_ID]]$file)
MAF_MIN    <- 0.01
INFO_MIN   <- 0.80
PP4_THRESH <- 0.80

stopifnot(file.exists(EGFR_FILE))
cat("eGFR file:", EGFR_FILE, "\n")
cat("File size:", round(file.info(EGFR_FILE)$size / 1e6, 0), "MB\n\n")

# THE LAUNCHER RUNS BEFORE ANY DATA IS READ, WHICH IS THE ONLY PLACE IT WORKS. It used to sit after
# load_rsid_map() and the suggestive read, so the driver spent ~8 minutes loading a 2.6M-row map it
# never uses and was still holding ~3.4 GB when the workers began loading theirs -- the peak of the
# whole run. 04_coloc.R's launcher is early for the same reason: load_rsid_map() lives inside
# main(), which its parallel path never calls, so its driver never reads the map at all. Everything
# below needs only config (RESULTS_DIR, PHASE1_TAG, BASE_DIR) and the environment.
# 2026-06-30 (gate obj #3): the CORRECTED cascade is SERIAL-only BY DEFAULT — force-clear any LEAKED chunk/
# merge env vars so a stale EGFR_CHUNK_N/EGFR_MERGE cannot silently glob/merge partial chunks. 2026-07-01: the
# dedicated parallel driver (cascade/04b_coloc_egfr_parallel.ps1) opts in with EGFR_PARALLEL=1; only THEN are
# the chunk/merge vars honoured, and the MERGE branch below hardens coverage so a partial glob still aborts.
if (Sys.getenv("EGFR_PARALLEL") != "1") Sys.unsetenv(c("EGFR_CHUNK_N","EGFR_CHUNK_I","EGFR_MERGE"))   # 2026-07-01 (Codex #2): EXACTLY "1" (a leaked EGFR_PARALLEL=0/false must NOT enable partial chunk/merge)
# --- optional parallelism: contiguous pair chunks, merged later (result == sequential) ---
.CN <- as.integer(Sys.getenv("EGFR_CHUNK_N", "0")); .CI <- as.integer(Sys.getenv("EGFR_CHUNK_I", "0"))
.MERGE <- Sys.getenv("EGFR_MERGE", "") %in% c("1","true","TRUE","yes")
.CHUNKED <- .CN > 0 && !.MERGE

# ── SELF-DRIVING PARALLEL, the way 04_coloc.R already works ────────────────────────────────────
# 2026-08-06. Until now this script only ever ran ONE chunk per process, and something outside had
# to split the work, start N processes, wait for them and ask for the merge. 04_coloc.R does all of
# that itself; 05 did not, and the dedicated driver it named (cascade/04b_coloc_egfr_parallel.ps1)
# starts workers with Start-Process, which REBUILD_LEDGER §9al records vanishing the whole process
# tree within 90 seconds, twice. So in practice the split was done by hand every time -- and on
# 2026-08-06 a hand-written runner called this script with no chunk variables at all, which the
# guard above quietly turns into a serial run: 1.34 pairs/min against 46 minutes for all 332 with
# five workers. A default that silently costs four hours is not a safe default.
#
# The opt-in guard above stays exactly as it is; this does not weaken it. Workers are started by
# this process with the variables set explicitly, so a LEAKED EGFR_CHUNK_N still cannot reach a
# worker or trigger a partial merge -- the driver owns the whole lifecycle, which is why 04 is
# allowed to be parallel by default and this was not.
#
# Launch mechanics are 04's, for the reasons recorded there: one cmd.exe per chunk writing its own
# %ERRORLEVEL%, rather than PSOCK nodes (which cost N idle R sessions and pushed commit to 96%,
# turning every pair into "python failed (exit 322)") and rather than Start-Process.
# EGFR_WORKERS IS READ AND VALIDATED ONLY ON THE PATH THAT USES IT. Validating it first made a
# typo'd value abort the serial, worker and merge modes as well -- three modes that predate this
# block and do not read the variable.
.SELF_DRIVE <- !.CHUNKED && !.MERGE && !(Sys.getenv("EGFR_SERIAL") %in% c("1","true","TRUE","yes"))
if (.SELF_DRIVE) {
  .WORKERS <- suppressWarnings(as.integer(Sys.getenv("EGFR_WORKERS", "5")))
  if (is.na(.WORKERS) || .WORKERS < 1L)
    stop(sprintf("05 eGFR: EGFR_WORKERS must be a positive integer, got '%s'",
                 Sys.getenv("EGFR_WORKERS")))
  if (.WORKERS == 1L) .SELF_DRIVE <- FALSE          # one worker is the serial path, not a driver
}
if (.SELF_DRIVE) {
  .drv_out <- Sys.getenv("COLOC_EGFR_OUT",
                         unset = file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG)))
  .drv_log <- file.path(BASE_DIR, "reports/Paper3_Program/run_logs/p3_egfr_workers")
  dir.create(.drv_log, showWarnings = FALSE, recursive = TRUE)
  .drv_self <- normalizePath(file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/05_coloc_egfr.R"),
                             winslash = "/", mustWork = TRUE)
  .drv_R    <- file.path(R.home("bin"), "Rscript.exe")
  .drv_exit <- function(k) file.path(.drv_log, sprintf("egfr_worker_%d_exit.txt", k))
  .w <- function(p) gsub("/", "\\\\", p)

  # EVERY CHUNK VARIABLE IS SET INSIDE THE BATCH FILE, NOT PASSED THROUGH system2(env=). Measured on
  # this machine 2026-08-06: system2(..., env = "X=1") exits 5 and the child sees nothing -- the
  # argument is a no-op here. The first version of this block used it for the merge pass, which
  # would have started the merge with neither EGFR_MERGE nor EGFR_PARALLEL set, and line 92 turns
  # that into a full SERIAL re-run of all 332 pairs. Silent, and four hours long.
  .mkbat <- function(path, vars, logstem) {
    writeLines(c("@echo off",
                 vars,
                 sprintf('"%s" "%s" > "%s" 2> "%s"', .w(.drv_R), .w(.drv_self),
                         .w(file.path(.drv_log, paste0(logstem, ".out.log"))),
                         .w(file.path(.drv_log, paste0(logstem, ".err.log"))))),
               path)
  }
  # CHUNK IDS ARE ZERO-BASED. The slicing below computes `.lo <- .CI * .sz + 1`, so .CI = 0 is the
  # FIRST chunk, and the pre-existing driver (04b_coloc_egfr_parallel.ps1:37) loops `0..N-1`. The
  # first version of this block used 1..N, which is the same mistake that was made by hand on
  # 2026-08-06: chunk 0 -- pairs 1-67 of 332 -- was never assigned to anyone, and chunk 5 got
  # `rows 336-332`, an empty slice. The merge's coverage assertion catches it, but only after the
  # whole run has been spent.
  .KS <- seq_len(.WORKERS) - 1L
  # NOTHING TO FREE HERE. An earlier version of this block sat after the rsID map and suggestive
  # reads and called rm() on them; from this position those objects do not exist yet and rm() would
  # abort the driver. The driver now never reads them at all, which is the point of moving it.
  # Any chunk output still on disk is a PREVIOUS run's. The merge below does check every chunk's
  # coloc_inputs_key, so a stale file from a different generation aborts rather than merges -- but a
  # stale file from an EARLIER RUN OF THIS SAME SCRIPT carries a matching key, and would be accepted
  # if its worker exited 0 without rewriting it. Matched by GLOB, not by expected id, the way 04b
  # does: a leftover from a run with a DIFFERENT worker count (2026-08-06 left a chunk05 behind) is
  # not named by this run's ids, and the merge counts files rather than checking their numbers, so a
  # stray would make it find N+1 chunks and abort.
  #
  # MOVED ASIDE, NOT DELETED. A worker deletes its own checkpoint when it finishes (line ~852), so a
  # completed chunk's result file is the ONLY surviving record of that chunk's work -- deleting it
  # costs a full recomputation of those pairs, and this happens immediately before the step most
  # likely to fail. Renaming is as good a guarantee of a clean directory and loses nothing. The
  # per-chunk checkpoints (egfr_chunk%02d_ckpt.rds) do not match these patterns and stay put, so
  # resume still works.
  .old <- c(Sys.glob(file.path(.drv_out, "coloc_egfr_chunk*.rds")),
            Sys.glob(file.path(.drv_out, "coloc_egfr_chunk*_status.csv")))
  if (length(.old)) {
    .quar <- file.path(.drv_out, format(Sys.time(), "_superseded_egfr_chunks_%Y%m%d_%H%M%S"))
    dir.create(.quar, showWarnings = FALSE, recursive = TRUE)
    if (!all(file.rename(.old, file.path(.quar, basename(.old)))))
      stop(sprintf("05 eGFR: could not move %d pre-existing chunk file(s) out of %s -- refusing to start a run that would mix them in",
                   length(.old), .drv_out))
    cat(sprintf("[eGFR] %d pre-existing chunk file(s) moved to %s\n", length(.old), basename(.quar)))
  }
  for (k in .KS) unlink(.drv_exit(k))                   # a stale code from a previous run is not this run's (04:136)
  for (k in .KS) {
    .bat <- file.path(.drv_log, sprintf("_run_egfr_worker_%d.bat", k))
    .mkbat(.bat,
           c("set EGFR_PARALLEL=1",
             sprintf("set EGFR_CHUNK_N=%d", .WORKERS),
             sprintf("set EGFR_CHUNK_I=%d", k),
             "set EGFR_MERGE="),
           sprintf("egfr_worker_%d", k))
    # cmd writes the worker's own exit code where the driver can read it; system2(wait = FALSE)
    # returns no pid on Windows, so there is nothing else to read (04:126-133).
    write(sprintf('echo %%ERRORLEVEL%% > "%s"', .w(.drv_exit(k))), .bat, append = TRUE)
    system2("cmd", c("/c", shQuote(.w(.bat))), wait = FALSE)
  }
  cat(sprintf("[eGFR] %d workers launched -> waiting (logs: %s)\n", .WORKERS, .drv_log))

  # Wait on what cmd writes AFTER the worker terminates, not on a process table: a machine-wide
  # Rscript count cannot tell this run's workers from anyone else's (04:129).
  repeat {
    if (all(vapply(.KS, function(k) file.exists(.drv_exit(k)), logical(1)))) break
    Sys.sleep(15)
  }
  .codes <- vapply(.KS, function(k) {
    .v <- suppressWarnings(as.integer(trimws(readLines(.drv_exit(k), warn = FALSE)[1])))
    if (is.na(.v)) NA_integer_ else .v
  }, integer(1))
  .bad <- .KS[is.na(.codes) | .codes != 0L]              # .KS[...] and not which(): report the CHUNK id, not its position
  if (length(.bad))
    stop(sprintf(paste0("05 eGFR: worker(s) FAILED -- chunk(s) %s exited %s. See %s/egfr_worker_*.err.log. ",
                        "NOT merging a partial eGFR coloc."),
                 paste(.bad, collapse = ","),
                 paste(.codes[is.na(.codes) | .codes != 0L], collapse = ","), .drv_log))
  # Every worker writes one result and one status file; the merge branch below re-checks coverage,
  # but a missing file here means a worker exited 0 without producing anything, which is worth
  # naming before the merge reports it as a count mismatch.
  .missing <- Filter(Negate(is.null), lapply(.KS, function(k) {
    .r <- file.path(.drv_out, sprintf("coloc_egfr_chunk%02d.rds", k))
    .s <- file.path(.drv_out, sprintf("coloc_egfr_chunk%02d_status.csv", k))
    if (!file.exists(.r) || !file.exists(.s)) k else NULL
  }))
  if (length(.missing))
    stop(sprintf(paste0("05 eGFR: worker(s) %s exited 0 but wrote no result/status file -- they ",
                        "stopped before finishing; NOT merging a partial eGFR coloc."),
                 paste(unlist(.missing), collapse = ",")))
  unlink(file.path(.drv_log, sprintf("_run_egfr_worker_%d.bat", .KS)))
  cat(sprintf("[eGFR] all %d workers exited 0 -> merging\n", .WORKERS))

  # Re-enter this same script in MERGE mode rather than duplicating the merge here: the merge branch
  # below is the reviewed one, and a second copy of it would drift. Same batch-file mechanism as the
  # workers, for the system2(env=) reason recorded above; wait = TRUE because there is nothing to
  # overlap with.
  .mbat <- file.path(.drv_log, "_run_egfr_merge.bat")
  .mkbat(.mbat, c("set EGFR_PARALLEL=1", sprintf("set EGFR_CHUNK_N=%d", .WORKERS), "set EGFR_CHUNK_I=",
                  "set EGFR_MERGE=1"), "egfr_merge")
  .mrc <- system2("cmd", c("/c", shQuote(.w(.mbat))), wait = TRUE)
  if (!identical(as.integer(.mrc), 0L))
    stop(sprintf("05 eGFR: the merge pass exited %s -- see %s/egfr_merge.err.log", .mrc, .drv_log))
  unlink(.mbat)
  cat("=== A2-E eGFR DONE (self-driven, ", .WORKERS, " workers) ===\n", sep = "")
  quit(save = "no", status = 0)
}

# ── Load rsID map ────────────────────────────────────────────────────────────
load_rsid_map <- function() {
  f <- file.path(BASE_DIR, "cache/cis_rsid_map.tsv")
  if (!file.exists(f)) { cat("[warn] cis_rsid_map.tsv not found\n"); return(NULL) }
  dt <- fread(f, colClasses = list(character = "chrom"))
  # KEYED ON THE ALLELES TOO -- same edit and same reason as R/coloc.R's load_rsid_map(): a
  # position-only key collapsed two variants at one coordinate to one rsID, and the indel was then
  # lost rather than mislabelled.
  dt[, lookup_key := paste(chrom, pos, toupper(allele0), toupper(allele1), sep = ":")]
  dt <- dt[order(as.integer(chrom), pos)]
  # ONE RSID -> ONE POSITION, BUT EVERY ALLELE PAIR AT THAT POSITION -- same edit and same reason
  # as R/coloc.R. `.SD[1]` collapsed the legacy-multiple-POSITIONS case it was written for AND the
  # distinct-variants case it was not: 1:1251285 carries G/A and G/T under one rsID (rs6657544).
  dt <- dt[, .SD[chrom == chrom[1L] & pos == pos[1L]], by = rsid]
  dt <- dt[, .SD[1], by = lookup_key]
  cat(sprintf("[rsID map] Loaded %s entries\n\n", formatC(nrow(dt), big.mark = ",")))
  dt
}

rsid_map <- load_rsid_map()

# ── Load suggestive eGFR pairs ───────────────────────────────────────────────
sug <- readRDS(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "mr_suggestive.rds"))   # CORRECTED input
setDT(sug)
pairs <- sug[outcome_id == EGFR_ID & method %in% c("Inverse variance weighted", "Wald ratio")]
pairs <- unique(pairs[, .(gene, outcome_id, outcome_name)])
cat(sprintf("eGFR suggestive pairs: %d\n\n", nrow(pairs)))

# THE eGFR CHECKPOINT IDENTITY. It used to be the literal "egfr_coloc|2026-07-09|<id>", which
# fingerprints nothing that can change a posterior: not this script, not the engine, the instruments,
# the cis window, the rsID map, the liftover chain, the outcome bytes or the suggestive input. A
# checkpoint from any of those could be resumed under any other. The non-eGFR path has carried a full
# fingerprint for exactly this reason; this reuses it, scoped to the eGFR genes and the one outcome,
# and prefixes the hash of THIS script because the shared key describes the generic engine and
# workers rather than the eGFR-specific code below. Built here, after `pairs` and before any chunk
# slicing, so every chunk and the merge all compute the same key from the same full pair set.
source(file.path(BASE_DIR, "R/_coloc_inputs_key.R"))
.EG_PROV_KEY <- paste(
  unname(tools::md5sum(file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/05_coloc_egfr.R"))),
  coloc_inputs_key(list(genes = sort(unique(pairs$gene)), outcome_ids = EGFR_ID))$key,
  sep = "|")
cat(sprintf("eGFR inputs key: %s\n\n", substr(.EG_PROV_KEY, 1, 12)))

# `unset=` and not nzchar(): reviewed 2026-08-01 as a possible hole (an exported but EMPTY
# COLOC_EGFR_OUT would give "" and write every result into the working directory). Not reachable on
# Windows -- SetEnvironmentVariable with an empty value DELETES the variable, so Sys.setenv(X="")
# leaves X absent from Sys.getenv() and `unset=` fires. Tested, not assumed.
.out_dir0 <- Sys.getenv("COLOC_EGFR_OUT", unset = file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG)))   # COLOC_EGFR_OUT redirects ALL eGFR coloc output to a scratch dir (non-destructive decoy-null capture); canonical untouched
if (.MERGE) {
  # 2026-07-01: HARDENED parallel merge. `pairs` above is the FULL eGFR suggestive set (chunk-subset happens
  # LATER at .CHUNKED, after this branch quits). Require EGFR_CHUNK_N chunk RESULT + STATUS files, union them,
  # and assert EXACT coverage (every expected pair has a status, none missing/extra/errored) so a partial or
  # stale glob ABORTS rather than silently shipping fewer pairs (the gate obj #3 risk).
  .NC <- as.integer(Sys.getenv("EGFR_CHUNK_N", "0"))
  if (.NC < 1) stop("06 eGFR MERGE: EGFR_CHUNK_N must be the worker chunk count")
  .expk <- paste(pairs$gene, EGFR_ID, sep = "_")
  .sf <- sort(Sys.glob(file.path(.out_dir0, "coloc_egfr_chunk*_status.csv")))
  if (length(.sf) != .NC) stop(sprintf("06 eGFR MERGE: found %d chunk STATUS files, expected %d -- a worker did not finish; refusing a PARTIAL merge", length(.sf), .NC))
  .psM <- rbindlist(lapply(.sf, fread), fill = TRUE)
  if (anyDuplicated(.psM$pair)) stop(sprintf("06 eGFR MERGE: %d duplicate/conflicting status key(s) across chunks (e.g. %s) -- chunks must be disjoint", sum(duplicated(.psM$pair)), paste(head(.psM$pair[duplicated(.psM$pair)], 5), collapse = ", ")))   # 2026-07-01 (Codex #3)
  .errM <- .psM[status == "error"]
  if (nrow(.errM) > 0) stop(sprintf("06 eGFR MERGE: %d pair(s) ERRORED across chunks: %s", nrow(.errM), paste(head(.errM$pair, 10), collapse = ", ")))
  .missM <- setdiff(.expk, .psM$pair); .extraM <- setdiff(.psM$pair, .expk)
  if (length(.missM) > 0)  stop(sprintf("06 eGFR MERGE: %d expected eGFR pair(s) have NO status (dropped by a partial chunk): %s", length(.missM), paste(head(.missM, 10), collapse = ", ")))
  if (length(.extraM) > 0) stop(sprintf("06 eGFR MERGE: %d UNEXPECTED pair(s) in status (stale/mismatched chunk): %s", length(.extraM), paste(head(.extraM, 10), collapse = ", ")))
  # 2026-07-01 (Codex #1): require EXACTLY EGFR_CHUNK_N result files too (each worker writes one, even if empty)
  # -- a missing/stale result.rds must NOT pass just because its status.csv survived.
  .rf <- sort(Sys.glob(file.path(.out_dir0, "coloc_egfr_chunk*.rds")))
  if (length(.rf) != .NC) stop(sprintf("06 eGFR MERGE: found %d chunk RESULT files, expected %d -- a worker's result is missing; refusing a PARTIAL merge", length(.rf), .NC))
  # EVERY CHUNK MUST CARRY THIS RUN'S KEY. rbindlist() drops attributes, so without checking them
  # first a chunk computed under different inputs would be merged in and the merged object would then
  # be stamped with the current key -- laundering the stale rows.
  .rr <- lapply(.rf, readRDS)
  if (any(!vapply(.rr, function(x) identical(attr(x, "coloc_inputs_key", exact = TRUE), .EG_PROV_KEY), logical(1))))
    stop("06 eGFR MERGE: at least one chunk result carries a missing or stale coloc_inputs_key -- it was computed under different inputs than this run; refusing to merge it")
  egfr_new <- rbindlist(.rr, fill = TRUE)
  attr(egfr_new, "coloc_inputs_key") <- .EG_PROV_KEY
  # 2026-07-01 (Codex #1): every status=ok pair MUST have a result row (a stale/empty result.rds would pass the
  # status checks yet silently drop rows).
  # A CLOSED STATUS VOCABULARY, A TWO-WAY BIJECTION, AND EACH ROW'S OWN ID. The three checks this
  # block claimed to perform were not all present. (a) The vocabulary was open, so `status="typo"`
  # counted as neither ok nor error and its pair passed with no result row. (b) Coverage was one-way:
  # a result row for a pair that is NOT status=ok was accepted. (c) The key was rebuilt with the
  # constant EGFR_ID rather than each row's own outcome_id, so a row carrying a different outcome was
  # relabelled as eGFR by the very check meant to catch it. 06 rejects all three before publication,
  # but only if this merge's output reaches it -- which is exactly what a laundered key defeats.
  .STAT_OK <- c("ok", "skip", "skip:no_instrument", "error")
  if (anyNA(.psM$status) || !all(.psM$status %chin% .STAT_OK))
    stop(sprintf("06 eGFR MERGE: missing or unrecognised pair status (%s) -- the vocabulary is %s",
                 paste(head(setdiff(unique(.psM$status), .STAT_OK), 5), collapse = ", "),
                 paste(.STAT_OK, collapse = "/")))
  .okk  <- .psM[status == "ok", pair]
  if (nrow(egfr_new) && (!"outcome_id" %in% names(egfr_new) || anyNA(egfr_new$outcome_id)))
    stop("06 eGFR MERGE: a chunk result row has no outcome_id, so it cannot be established as an eGFR result")
  .resk <- if (nrow(egfr_new)) paste(egfr_new$gene, egfr_new$outcome_id, sep = "_") else character(0)
  if (anyDuplicated(.resk))
    stop(sprintf("06 eGFR MERGE: %d duplicate (gene,outcome_id) result row(s) across chunks (must be disjoint): %s",
                 sum(duplicated(.resk)), paste(head(.resk[duplicated(.resk)], 5), collapse = ", ")))
  .mres <- setdiff(.okk, .resk)
  if (length(.mres) > 0) stop(sprintf("06 eGFR MERGE: %d pair(s) status=ok but with NO result row (stale/empty chunk result.rds): %s", length(.mres), paste(head(.mres, 10), collapse = ", ")))
  .xres <- setdiff(.resk, .okk)
  if (length(.xres) > 0) stop(sprintf("06 eGFR MERGE: %d result row(s) whose pair is NOT status=ok (stale/malformed chunk output): %s", length(.xres), paste(head(.xres, 10), collapse = ", ")))
  # 2026-07-01 (Codex #4): serial 05 refuses to save an empty eGFR coloc -> the parallel merge must match.
  if (nrow(egfr_new) == 0) stop("06 eGFR MERGE: 0 eGFR coloc results across ALL chunks (data/n_total_sum issue) -- refusing to save empty coloc_egfr_results.rds (matches serial guard)")
  saveRDS(egfr_new, file.path(.out_dir0, "coloc_egfr_results.rds"))
  fwrite(.psM, file.path(.out_dir0, "coloc_egfr_pair_status.csv"))
  # "strong" is the DUAL rule, here as everywhere else. Counting PP.H4 >= 0.8 alone reported a number
  # the engine's coloc_support and 06's final table both reject (H4 = 0.81 with H3 = 0.17 passes the
  # H4 threshold and fails the ratio), so the merge's own summary disagreed with the published count.
  cat(sprintf("MERGE OK: %d chunks; pair-status covers all %d expected eGFR pairs (%d ok / %d skip); %d rows, %d strong(dual rule H4>=.80 & H4/H3>=5)\n",
              .NC, length(.expk), sum(.psM$status == "ok"), sum(grepl("skip", .psM$status)), nrow(egfr_new), sum(is_coloc_pass(egfr_new$PP.H4, egfr_new$PP.H3), na.rm = TRUE)))
  quit(save = "no")
}
if (.CHUNKED) {
  .sz <- ceiling(nrow(pairs) / .CN); .lo <- .CI * .sz + 1; .hi <- min((.CI + 1) * .sz, nrow(pairs))
  pairs <- if (.lo <= nrow(pairs)) pairs[.lo:.hi] else pairs[0]
  cat(sprintf("CHUNK %d/%d: %d pairs (rows %d-%d)\n", .CI, .CN, nrow(pairs), .lo, .hi))
}

# ── Load instruments for gene coordinates ────────────────────────────────────
# read_phase1_instruments(), not a second readRDS of a path this script picks itself: config pinned
# and certified ONE snapshot when it was sourced, and the pair set above was selected against that
# generation. Re-opening the file re-reads whatever is on disk NOW, with no md5 comparison, so a
# Phase 1 replaced mid-run gave pair selection from generation A and window coordinates from B.
instruments <- as.data.table(read_phase1_instruments())

# ── Helper functions (same as R/03_coloc_egfr.R) ────────────────────────────
complement_allele <- function(a) chartr("ACGTacgt", "TGCAtgca", toupper(a))
is_palindromic    <- function(a1, a2) complement_allele(a1) == toupper(a2)

harmonize_alleles <- function(dt, ea_p, oa_p, ea_o, oa_o, beta_o) {
  ep <- toupper(dt[[ea_p]]); op <- toupper(dt[[oa_p]])
  eo <- toupper(dt[[ea_o]]); oo <- toupper(dt[[oa_o]])
  comp_ep <- complement_allele(ep); comp_op <- complement_allele(op)
  same      <- eo == ep  & oo == op
  flip_beta <- eo == op  & oo == ep
  same_str  <- eo == comp_ep & oo == comp_op
  flip_str  <- eo == comp_op & oo == comp_ep
  palin      <- is_palindromic(ep, op)
  # SNVs ONLY in the strand branches -- same reason and same edit as R/coloc.R: complement_allele()
  # is chartr(), which complements without reversing ("GGA" -> "CCT", not "TCC"), so a real
  # reverse-complemented indel is dropped while an unrelated multi-base pair can match by coincidence.
  # Reversing the string would not help: a point liftover does not re-anchor an indel on the target
  # reference. Exact and swapped matches still accept multi-base alleles.
  strand_snv <- nchar(ep) == 1L & nchar(op) == 1L & nchar(eo) == 1L & nchar(oo) == 1L
  keep       <- allele_compatible(ep, op, eo, oo)   # engine predicate (R/coloc.R, sourced above): the rsID merge below asks the same question, and two copies would drift
  needs_flip <- flip_beta | (strand_snv & !palin & flip_str)
  dt <- dt[keep]
  if (nrow(dt) > 0 && any(needs_flip[keep]))
    dt[needs_flip[keep], (beta_o) := -get(beta_o)]
  dt
}

get_pqtl_region <- function(gene, chrom, region_start, region_end) {
  tar_file <- resolve_ukbppp_tar(gene)   # 2026-07-01 (gate finding 1): canonical assay by instrument OID, not [1]
  if (is.null(tar_file)) return(NULL)

  py <- tempfile(fileext = ".py")
  writeLines(sprintf('
import tarfile, gzip, pandas as pd, sys
with tarfile.open(r"%s", "r:") as tf:
    for m in tf.getmembers():
        if "chr%s_" in m.name and m.name.endswith(".gz"):
            f = tf.extractfile(m)
            with gzip.open(f) as gz:
                df = pd.read_csv(gz, sep=" ")
            df = df[(df["GENPOS"] >= %d) & (df["GENPOS"] <= %d)]
            df.to_csv(sys.stdout, index=False)
            break
', gsub("\\\\", "/", tar_file), chrom, region_start, region_end), py)

  # stderr KEPT -- same edit and same reason as R/coloc.R: a discarded stderr turned a whole failed
  # run into a bare exit code (2026-08-02, exit 322 on all 1,380 pairs; six hours to establish
  # that python itself was fine and the machine had simply run out of commit).
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(.pyerr)
  unlink(py)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("eGFR region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))  # 2026-07-01 gate
  if (length(out) <= 1) return(NULL)

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("eGFR coloc region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep)
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("eGFR coloc region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  dt[, chrpos := paste0(chrom, ":", GENPOS)]
  if ("A1FREQ" %in% names(dt))
    dt <- dt[!is.na(A1FREQ) & A1FREQ >= MAF_MIN & A1FREQ <= 1 - MAF_MIN]

  if (!is.null(rsid_map)) {
    dt[, lookup_key := paste(chrpos, toupper(ALLELE0), toupper(ALLELE1), sep = ":")]
    dt <- merge(dt, rsid_map[, .(lookup_key, rsid)],
                by = "lookup_key", all.x = TRUE)
  } else {
    dt[, rsid := NA_character_]
  }
  dt
}

get_egfr_region <- function(chrom, region_start, region_end) {
  py <- tempfile(fileext = ".py")
  writeLines(sprintf('
import gzip, sys
chrom_target = "%s"
start = %d
end   = %d
with gzip.open(r"%s", "rt") as f:
    hdr = f.readline().strip()
    cols = hdr.split()
    ci = cols.index("Chr")
    pi = cols.index("Pos_b37")
    print("\\t".join(cols))
    for line in f:
        row = line.strip().split()
        if len(row) < len(cols):
            raise ValueError("short GWAS association row")
        try:
            chr_val = row[ci].lstrip("chr")
            if chr_val == chrom_target and start <= int(row[pi]) <= end:
                print("\\t".join(row))
        except (ValueError, IndexError):
            raise
', as.character(chrom), region_start, region_end,
   gsub("\\\\", "/", EGFR_FILE)), py)

  # stderr KEPT -- same edit and same reason as R/coloc.R: a discarded stderr turned a whole failed
  # run into a bare exit code (2026-08-02, exit 322 on all 1,380 pairs; six hours to establish
  # that python itself was fine and the machine had simply run out of commit).
  .pyerr <- tempfile(fileext = ".stderr")
  out <- system2(PYTHON_BIN, args = py, stdout = TRUE, stderr = .pyerr)
  .emsg <- if (file.exists(.pyerr)) paste(readLines(.pyerr, warn = FALSE), collapse = " | ") else ""
  unlink(.pyerr)
  unlink(py)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop(sprintf("eGFR region extraction python failed (exit %s) -- NOT a legit empty-region skip. python stderr: %s", attr(out, "status"), if (nzchar(.emsg)) .emsg else "<empty: no diagnostic from the process, which points at process creation rather than the script>"))  # 2026-07-01 gate
  if (length(out) <= 1) { cat("    [skip] eGFR region empty\n"); return(NULL) }

  dt <- suppressWarnings(tryCatch(fread(text = paste(out, collapse = "\n")), error = function(e) fail_loud_data(sprintf("eGFR coloc region fread PARSE failed on a NON-empty extraction: %s -- genuine parse/format failure, not a legit empty region", conditionMessage(e)))))  # 2026-07-01 (closure sweep)
  if (is.null(dt) || nrow(dt) == 0) fail_loud_data("eGFR coloc region parsed to 0 rows despite a non-empty extraction -- genuine parse/format failure")

  setnames(dt, trimws(names(dt)))
  if ("Chr" %in% names(dt))
    dt[, Chr := gsub("^chr", "", as.character(Chr))]

  dt[, rsid_out := RSID]
  dt[is.na(rsid_out) | !grepl("^rs", rsid_out), rsid_out := NA_character_]
  dt[, chrpos_out := paste0(Chr, ":", Pos_b37)]
  dt
}

rsid_pos_merge <- function(pqtl_sub, out_sub) {
  result <- list()
  matched_pqtl_var <- character(0)   # (chrpos:EA:OA) of the pQTL rows the primary merge already used

  has_p_rsid <- any(!is.na(pqtl_sub$rsid))
  has_o_rsid <- any(!is.na(out_sub$rsid_out))

  if (has_p_rsid && has_o_rsid) {
    # THE COORDINATE HAS TO AGREE TOO -- same edit and same reason as R/coloc.R, which this script
    # copies. rsID alone accepted an outcome record with that rsID at a different position, and the
    # pQTL position was then marked matched, so the correct record could not be reached by fallback.
    p_r <- pqtl_sub[!is.na(rsid) & !is.na(chrpos_pqtl), .(rsid, chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_r <- out_sub[!is.na(rsid_out) & !is.na(chrpos_out), .(rsid_out, chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    rsid_m <- merge(p_r, o_r, by.x = c("rsid", "chrpos_pqtl"), by.y = c("rsid_out", "chrpos_out"))
    if (nrow(rsid_m) > 0) {
      rsid_m[, snp_id := rsid]
      # THE PRIMARY MERGE CONSUMES A VARIANT, NOT A COORDINATE -- same edit as R/coloc.R, whose
      # matching this script copies. Excluding by position was wrong twice over: an rsID hit whose
      # alleles cannot harmonise closed the coordinate for the compatible record at it, and even a
      # harmonisable hit closed the coordinate for every OTHER pQTL variant sitting there (a
      # multiallelic site, or an SNV beside an indel). The (position, alleles) key excludes exactly
      # the pQTL row the primary merge used and nothing else.
      .ok <- allele_compatible(rsid_m$ea_pqtl, rsid_m$oa_pqtl, rsid_m$ea_out, rsid_m$oa_out)
      matched_pqtl_var <- rsid_m[!is.na(chrpos_pqtl) & !is.na(.ok) & .ok,
                                 paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":")]
      result[["rsid"]] <- rsid_m
      cat(sprintf("    rsID merge: %d SNPs\n", nrow(rsid_m)))
    }
  }

  has_p_pos <- any(!is.na(pqtl_sub$chrpos_pqtl))
  has_o_pos <- any(!is.na(out_sub$chrpos_out))

  if (has_p_pos && has_o_pos) {
    p_fb <- pqtl_sub[!is.na(chrpos_pqtl) &
                     !(paste(chrpos_pqtl, toupper(ea_pqtl), toupper(oa_pqtl), sep = ":") %in% matched_pqtl_var),
                     .(chrpos_pqtl, pqtl_beta, pqtl_varbeta, ea_pqtl, oa_pqtl)]
    o_fb <- out_sub[!is.na(chrpos_out), .(chrpos_out, out_beta, out_varbeta, ea_out, oa_out)]
    pos_m <- merge(p_fb, o_fb, by.x = "chrpos_pqtl", by.y = "chrpos_out")
    if (nrow(pos_m) > 0) {
      pos_m[, snp_id := chrpos_pqtl]
      result[["pos"]] <- pos_m
      cat(sprintf("    position fallback: +%d SNPs\n", nrow(pos_m)))
    }
  }

  if (length(result) == 0) { cat("    merge: 0 SNPs\n"); return(data.table()) }
  rbindlist(result, use.names = TRUE, fill = TRUE)
}

run_egfr_coloc_pair <- function(gene, chrom, gene_start, gene_end) {
  region_start <- max(1L, gene_start - CIS_WINDOW_KB * 1000L)
  region_end   <- gene_end + CIS_WINDOW_KB * 1000L

  cat(sprintf("  coloc: %s ~ eGFR | chr%s:%d-%d\n", gene, chrom, region_start, region_end))
  flush(stdout())

  pqtl_dt <- get_pqtl_region(gene, chrom, region_start, region_end)
  # 2026-07-01 (gate #14 finding 1): eGFR analogue of the 03 fix — a resolved tar returning NO pQTL data,
  # or missing required columns, is a DATA-FOUNDATION failure (not a legit low-SNP skip) -> error, not skip.
  if (is.null(pqtl_dt) || nrow(pqtl_dt) == 0)
    stop(sprintf("eGFR coloc: gene %s pQTL extraction returned NO data (tar resolved but empty) -- data-foundation failure", gene))

  req_cols <- c("ALLELE0", "ALLELE1", "BETA", "SE", "INFO", "N")
  if (length(setdiff(req_cols, names(pqtl_dt))) > 0)
    stop(sprintf("eGFR coloc: gene %s pQTL missing required columns: %s -- data-format failure", gene, paste(setdiff(req_cols, names(pqtl_dt)), collapse = ", ")))

  pqtl_dt <- pqtl_dt[!is.na(INFO) & INFO >= INFO_MIN & !is.na(BETA) & !is.na(SE)]
  if (nrow(pqtl_dt) < 30) {
    cat(sprintf("    [skip] pQTL after filter: %d SNPs\n", nrow(pqtl_dt))); return(NULL)
  }

  pqtl_dt[, pqtl_varbeta := SE^2]
  n_pqtl <- as.integer(pqtl_dt[["N"]][1])

  # PROJECTED TO GRCh37, like R/coloc.R's hg37 branch. pqtl_dt$chrpos is chrom:GENPOS from the
  # UKB-PPP archive and GENPOS is GRCh38, while get_egfr_region() below filters CKDGen's `Pos_b37`.
  # Comparing them directly does not merely fail to match -- it matches the WRONG variants: on the
  # delivered eGFR analysis, 25 SNPs across 8 genes were two different variants merged into one row.
  pqtl_37 <- liftover_to_hg37(pqtl_dt)
  cat(sprintf("    liftover hg38->hg37: %d/%d\n",
              if (is.null(pqtl_37)) 0L else sum(!is.na(pqtl_37$chrpos37)), nrow(pqtl_dt)))
  .psrc <- if (!is.null(pqtl_37)) pqtl_37 else pqtl_dt
  pqtl_sub <- data.table(
    rsid         = .psrc$rsid,
    chrpos_pqtl  = if (!is.null(pqtl_37)) pqtl_37$chrpos37 else rep(NA_character_, nrow(pqtl_dt)),
    pqtl_beta    = .psrc$BETA,
    pqtl_varbeta = .psrc$pqtl_varbeta,
    ea_pqtl      = toupper(.psrc$ALLELE1),
    oa_pqtl      = toupper(.psrc$ALLELE0)
  )
  cat(sprintf("    pQTL SNPs: %d (rsID: %d)\n", nrow(pqtl_sub), sum(!is.na(pqtl_sub$rsid))))

  # ...AND SO IS THE WINDOW. region_start/region_end are GRCh38, derived from the GRCh38 gene body,
  # while get_egfr_region() filters CKDGen's Pos_b37 -- the comment above says exactly that about the
  # variants and then the bounds were handed over unconverted anyway. For DNER, a current eGFR
  # suggestive pair (p = 5.4e-08), the GRCh38 window 2:228857519-230214586 corresponds to roughly
  # 2:229722235-231079301 in GRCh37: most of the true window is missed and an equally large wrong
  # flank is added, and the pair still returns enough SNPs to look valid. The two ENDPOINTS are not
  # lifted -- a window's endpoints can land in different chain blocks or in reverse order -- so the
  # bounds come from the projected positions of the pQTLs actually being colocalised.
  # AN OFF-CHROMOSOME LIFT IS A PROPERTY OF ONE VARIANT, NOT OF THE PAIR -- same edit and same reason
  # as R/coloc.R. That variant has no GRCh37 counterpart usable against this outcome, so it is
  # dropped by clearing its coordinate (which excludes it from the rsID merge and the position
  # fallback alike) and the pair proceeds. Only having NO variant left on this chromosome is fatal,
  # because then there is no window at all.
  .chr37  <- if (!is.null(pqtl_37)) sub(":.*", "", pqtl_37$chrpos37) else NA_character_
  .same37 <- .chr37 %in% as.character(chrom)
  .off37  <- !is.na(.chr37) & !.same37
  if (!any(.same37))
    fail_loud_data(sprintf(paste0("eGFR coloc %s: the pQTL liftover to GRCh37 leaves NO variant on ",
                                  "chr%s (mapped=%d, off-chromosome=%d, unmapped=%d), so there is no ",
                                  "window. Falling back to the GRCh38 bounds would query the wrong ",
                                  "region silently."),
                           gene, chrom, sum(.same37), sum(.off37), sum(is.na(.chr37))))
  if (any(.off37)) {
    cat(sprintf("    [liftover] %d pQTL(s) lifted off chr%s -- dropped for this outcome\n",
                sum(.off37), chrom))
    pqtl_37$chrpos37[.off37] <- NA_character_
  }
  .pos37 <- as.integer(sub(".*:", "", pqtl_37$chrpos37[.same37]))
  out_dt <- get_egfr_region(chrom, min(.pos37), max(.pos37))
  if (is.null(out_dt) || nrow(out_dt) < 20) {
    cat(sprintf("    [skip] eGFR SNPs=%d\n", if (is.null(out_dt)) 0L else nrow(out_dt)))
    return(NULL)
  }

  if ("n_total_sum" %in% names(out_dt)) {
    out_n_median <- as.integer(median(out_dt$n_total_sum, na.rm = TRUE))
    if (!is.finite(out_n_median) || out_n_median <= 0)   # 2026-07-01 Codex-phase3 #5: derived eGFR N must be finite/positive (else NA N silently into coloc.abf)
      stop(sprintf("06 eGFR coloc: derived eGFR N (median n_total_sum) = %s -- non-finite/non-positive, malformed N column", out_n_median))
  } else {
    stop("06 eGFR coloc: n_total_sum column absent — refusing stale fallback N=783978; check Wuttke eGFR file")  # 2026-06-30 fail-loud (gate fix #4)
  }

  out_sub <- out_dt[, .(
    rsid_out,
    chrpos_out,
    ea_out      = toupper(Allele1),
    oa_out      = toupper(Allele2),
    out_beta    = as.numeric(Effect),
    out_varbeta = as.numeric(StdErr)^2
  )]
  out_sub <- out_sub[!is.na(out_beta) & !is.na(out_varbeta)]

  common <- rsid_pos_merge(pqtl_sub, out_sub)
  cat(sprintf("    total merged: %d SNPs\n", nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] merged < 20\n"); return(NULL) }

  # HARMONISE FIRST, THEN DROP WHAT IS STILL AMBIGUOUS -- same edit and same reason as R/coloc.R.
  # A position matching two outcome records had both deleted before the alleles were looked at, even
  # when harmonisation would have kept exactly one. The drop stays for groups still multiple after it.
  common <- harmonize_alleles(common, "ea_pqtl", "oa_pqtl", "ea_out", "oa_out", "out_beta")
  cat(sprintf("    allele align: %d SNPs\n", nrow(common)))

  # A REUSED RSID IS NOT AN AMBIGUOUS VARIANT -- same edit and same reason as R/coloc.R. Only the
  # colliding rows are re-keyed, so unambiguous rsIDs keep their names for the rsID-based LD lookup.
  common[duplicated(snp_id) | duplicated(snp_id, fromLast = TRUE),
         snp_id := paste(chrpos_pqtl, ea_pqtl, oa_pqtl, sep = ":")]

  n_before <- nrow(common)
  common <- common[, if (.N == 1L) .SD, by = snp_id]
  if (nrow(common) < n_before)
    cat(sprintf("    dedup: %d -> %d\n", n_before, nrow(common)))
  if (nrow(common) < 20) { cat("    [skip] aligned < 20\n"); return(NULL) }

  ds1 <- list(beta     = common$pqtl_beta,
              varbeta  = common$pqtl_varbeta,
              N        = n_pqtl,
              sdY      = 1,
              type     = "quant",
              snp      = common$snp_id)

  # sdY = 1 asserts a variance-one phenotype. These are CKDGen effects on natural-log eGFR, so the
  # assertion was wrong by a factor of five and, since coloc's quantitative prior is 0.15 * sdY with
  # no role for N, so was every eGFR Bayes factor. Taken from the outcome's GWAS_META entry (R/coloc.R,
  # sourced above), which is where the estimate and its provenance live -- one definition, not a
  # literal repeated at each call site.
  ds2 <- list(beta     = common$out_beta,
              varbeta  = common$out_varbeta,
              N        = out_n_median,
              sdY      = GWAS_META[[EGFR_ID]]$sdY,
              type     = "quant",
              snp      = common$snp_id)

  # Non-destructive capture hook (eGFR stratum) for the coloc empirical-null.
  .capdir <- Sys.getenv("COLOC_CAPTURE_DIR", "")
  if (nzchar(.capdir)) {
    dir.create(.capdir, showWarnings = FALSE, recursive = TRUE)
    # THE OUTCOME'S OWN META, with only N localised -- the engine saves the whole meta object and this
    # copy hand-built a two-field subset, so the capture described the outcome less completely than the
    # generic path did. Once sdY moved into GWAS_META that omission became load-bearing: the empirical
    # null reads meta$sdY, an absent field assigns NULL, and coloc rejects a quantitative dataset that
    # has neither sdY nor MAF. A subset that has to be kept in step with its source is a second
    # definition; taking the source and overriding N is not.
    saveRDS(list(gene = gene, outcome_id = EGFR_ID, n_pqtl = n_pqtl,
                 meta = modifyList(GWAS_META[[EGFR_ID]], list(n = out_n_median)),
                 common = common[, .(snp_id, pqtl_beta, pqtl_varbeta, out_beta, out_varbeta)]),
            file = file.path(.capdir, sprintf("cap__%s__eGFR.rds", gene)))
  }
  cr <- tryCatch(
    # EXPLICIT, and identical to R/coloc.R:851. This called coloc.abf() with no priors and took the
    # package defaults -- numerically the same today, but _coloc_inputs_key.R reads the priors out of
    # R/coloc.R and records them as "the priors this run uses", so the fingerprint had a field named
    # priors whose value came from a file this arm does not execute. Same defect class as the key
    # hashing LIFTOVER_CHAIN while the engine read LIFTOVER_CHAIN_38_TO_37. Priors carry results
    # here: 17 of the previous generation's 52 pairs no longer qualified at p12 = 5e-6.
    coloc::coloc.abf(dataset1 = ds1, dataset2 = ds2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    # 2026-07-01 (gate): coloc.abf FAILURE is not a legit skip -> propagate to error status -> abort.
    error = function(e) stop(sprintf("eGFR coloc.abf failed (data problem, NOT a legit skip): %s", conditionMessage(e)))
  )
  if (is.null(cr)) return(NULL)

  pp <- cr$summary
  # 2026-07-01 (editor REQUIRE): non-finite PP on a REACHED eGFR pair is degenerate; it would be recorded
  # status=="ok" then silently filtered from strong (!is.na & >=0.80). Fail loud, don't skip. Mirror of 03.
  if (!is.finite(pp["PP.H4.abf"]) || !is.finite(pp["PP.H3.abf"]))
    fail_loud_data(sprintf("eGFR coloc %s: coloc.abf returned NON-FINITE PP (H4=%s H3=%s) on a reached pair -- degenerate, refusing a silent status=ok", gene, pp["PP.H4.abf"], pp["PP.H3.abf"]))
  cat(sprintf("  -> PP.H4=%.3f  PP.H3=%.3f  nSNPs=%d\n",
              pp["PP.H4.abf"], pp["PP.H3.abf"], as.integer(pp["nsnps"])))
  flush(stdout())

  data.table(
    gene               = gene,
    outcome_id         = EGFR_ID,
    n_snps_coloc       = as.integer(pp["nsnps"]),
    # ALL FIVE UNROUNDED -- same edit as R/coloc.R. Rounding H0-H2 meant the stored posteriors were
    # not the ones coloc.abf returned and need not sum to one.
    PP.H0              = as.numeric(pp["PP.H0.abf"]),
    PP.H1              = as.numeric(pp["PP.H1.abf"]),
    PP.H2              = as.numeric(pp["PP.H2.abf"]),
    # downstream recomputes the support
    # rule from these values, so rounding here could flip the call at the 0.8 threshold.
    PP.H3              = as.numeric(pp["PP.H3.abf"]),
    PP.H4              = as.numeric(pp["PP.H4.abf"]),
    coloc_support      = is_coloc_pass(pp["PP.H4.abf"], pp["PP.H3.abf"]),   # dual coloc rule (D1; centralized 2026-07-05, was H4-only)
    independent_signal = pp["PP.H3.abf"] >= PP4_THRESH,
    outcome_N          = out_n_median   # 2026-07-01 Codex-phase3 #5: record the derived eGFR N so downstream can audit the exact N used
  )
}

# ── Main ─────────────────────────────────────────────────────────────────────
out_dir <- Sys.getenv("COLOC_EGFR_OUT", unset = file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG)))   # CORRECTED workspace (live untouched); COLOC_EGFR_OUT redirects to scratch for non-destructive decoy-null capture
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

checkpoint_file <- if (.CHUNKED) file.path(out_dir, sprintf("egfr_chunk%02d_ckpt.rds", .CI)) else file.path(out_dir, "coloc_egfr_rerun_checkpoint.rds")
all_results <- list()
done_pairs  <- character(0)

pair_status <- list()   # 2026-06-30 gate obj #3: per-pair status (ok / skip:reason / error) -> no silent drop
# Fix (audit Med): a top-level SERIAL eGFR run (not a chunk worker, not merge) starts FRESH -- clear a possibly-
# stale checkpoint (old code but same pair keys) + stale root outputs, so it cannot silently resume them or leave
# them behind for a later 06 (matches 04_coloc.R's serial/parallel clearing).
if (!.CHUNKED && !.MERGE) {
  .egst <- c(checkpoint_file, file.path(.out_dir0, "coloc_egfr_results.rds"), file.path(.out_dir0, "coloc_egfr_pair_status.csv"))
  unlink(.egst, force = TRUE)
  if (any(file.exists(.egst))) stop("05 eGFR coloc (serial): could NOT clear stale checkpoint/root outputs (locked/read-only)")
}
# .EG_PROV_KEY is built once, near the pair set at the top of this script -- see the comment there.
# It used to be re-assigned here as a fixed literal, which is what made the checkpoint identity
# describe nothing that could change a posterior.
if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  if (is.null(cp$prov_key) || !identical(cp$prov_key, .EG_PROV_KEY)) {
    cat("[checkpoint] PROVENANCE mismatch (stale code/inputs) -> fresh start\n"); cp <- list(results = list(), done_pairs = character(0))
  }
  # 2026-07-01 (gate finding 3): discard a checkpoint whose done pairs aren't in the current eGFR run.
  if (!is.null(cp$done_pairs) && length(setdiff(cp$done_pairs, paste(pairs$gene, EGFR_ID, sep="_"))) > 0) {
    cat("[checkpoint] STALE (done pairs not in current run) -> fresh start\n"); cp <- list(results = list(), done_pairs = character(0))
  }
  all_results <- cp$results
  done_pairs  <- cp$done_pairs
  if (!is.null(cp$pair_status)) pair_status <- cp$pair_status
  cat(sprintf("Resuming: %d done, %d results\n\n", length(done_pairs), length(all_results)))
} else {
  cat("Fresh start\n\n")
}

stop_file <- file.path(out_dir, "STOP_EGFR_COLOC")
t0 <- proc.time()["elapsed"]

for (i in seq_len(nrow(pairs))) {
  # Stop signal
  if (file.exists(stop_file)) {
    saveRDS(list(results = all_results, done_pairs = done_pairs, pair_status = pair_status, prov_key = .EG_PROV_KEY), checkpoint_file)
    cat(sprintf("\n[PAUSED] at pair %d/%d. Checkpoint saved.\n", i, nrow(pairs)))
    file.remove(stop_file)
    # "Re-run to resume" was true for a chunk worker and false for the top-level serial run, which
    # clears the checkpoint before reading it (line 486) so that a serial result is always a clean
    # comparator rather than a partial resume of possibly different code. Following the instruction
    # therefore restarted every pair while saying it was resuming. The clearing is deliberate; the
    # message is what was wrong.
    stop(if (.CHUNKED) "Paused by STOP file. Re-run this chunk to resume from the checkpoint."
         else paste("Paused by STOP file. A top-level SERIAL re-run does NOT resume: it clears the",
                    "checkpoint on startup so the serial path stays a clean comparator. Re-running",
                    "starts every pair again; use the chunked path if you need resumable progress."))
  }

  g   <- pairs$gene[i]
  key <- paste(g, EGFR_ID, sep = "_")

  if (key %in% done_pairs) {
    cat(sprintf("[%d/%d] SKIP: %s ~ eGFR\n", i, nrow(pairs), g))
    next
  }

  # THE ROW COUNT IS TAKEN BEFORE `[1]` -- same defect and same fix as R/coloc.R: `[1]` on a zero-row
  # data.table returns one row of NA, so this skip was unreachable and the pair proceeded with NA
  # coordinates, surfacing as an error rather than as the skip it is.
  .gi_all <- instruments[gene == g]
  gene_info <- .gi_all[1]
  if (nrow(.gi_all) == 0) {
    cat(sprintf("[%d/%d] %s: no instrument data\n", i, nrow(pairs), g))
    pair_status[[key]] <- "skip:no_instrument"; done_pairs <- c(done_pairs, key); next
  }

  cat(sprintf("\n[%d/%d] %s ~ eGFR\n", i, nrow(pairs), g))
  flush(stdout())

  .cerr <- FALSE
  res <- tryCatch(
    run_egfr_coloc_pair(g,
                        gene_info$gene_chrom,
                        gene_info$gene_start,
                        gene_info$gene_end),
    error = function(e) {
      message(sprintf("  [ERROR] %s: %s", g, e$message)); .cerr <<- TRUE; NULL
    }
  )

  if (!is.null(res)) { all_results[[key]] <- res; pair_status[[key]] <- "ok" } else
    pair_status[[key]] <- if (.cerr) "error" else "skip"   # gate obj #3: error (bad) vs legit skip
  # An ERRORED pair is not done -- same reason as R/coloc.R: recording it made the resume skip it as
  # completed while the end-of-loop check aborted on its status, so the checkpoint had to be deleted
  # by hand. A legitimate NULL skip has .cerr FALSE and is still recorded as done.
  if (!.cerr) done_pairs <- c(done_pairs, key)

  # Checkpoint every 10 pairs
  if (i %% 10 == 0 || i == nrow(pairs)) {
    saveRDS(list(results = all_results, done_pairs = done_pairs, pair_status = pair_status, prov_key = .EG_PROV_KEY), checkpoint_file)
  }

  # Progress
  elapsed <- proc.time()["elapsed"] - t0
  rate <- sum(!(done_pairs %in% done_pairs[1:0])) / max(elapsed, 1)
  # Simple ETA
  done_this_run <- sum(seq_len(i) > 0) - sum(pairs$gene[1:i] %in%
    sub(paste0("_", EGFR_ID), "", done_pairs[done_pairs %in% paste(pairs$gene[1:i], EGFR_ID, sep="_")]))
}

# ── Combine and save ─────────────────────────────────────────────────────────
# 2026-06-30 gate obj #3: persist per-pair status + FAIL LOUD if any eGFR pair ERRORED (vs legit skip).
.eps <- data.table(pair = names(pair_status), status = unlist(pair_status))
# THE CANONICAL STATUS FILE IS THE WHOLE RUN'S, so only a whole run writes it. Every chunk wrote this
# same path, so six workers overwrote one another and what survived was whichever finished last -- an
# arbitrary subset wearing the name of the complete table. A successful merge replaces it, so a clean
# run ended correct; an interrupted one left a partial status looking canonical. The chunk's own
# status is written below under its own name, and the merge assembles the canonical one from those,
# so nothing needed this write. The error check below still runs on every path.
if (!.CHUNKED) fwrite(.eps, file.path(out_dir, "coloc_egfr_pair_status.csv"))
.eerr <- .eps[status == "error"]
if (nrow(.eerr) > 0) stop(sprintf("06 eGFR coloc: %d pair(s) ERRORED (not a legit skip): %s -- investigate", nrow(.eerr), paste(.eerr$pair, collapse = ", ")))

cat("\n=== eGFR Coloc Results ===\n")
cat(sprintf("Pairs attempted: %d\n", length(done_pairs)))
cat(sprintf("Pairs with results: %d\n", length(all_results)))

if (length(all_results) > 0) {
  egfr_new <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  egfr_new[, outcome_name := "eGFR"]

  # The dual rule, not H4 alone -- see the merge summary above. The H4-only and H4>=0.50 lines are
  # kept as diagnostics but no longer carry the word "strong", which names one specific criterion.
  cat(sprintf("strong coloc (PP.H4 >= 0.80 & PP.H4/PP.H3 >= 5): %d\n", sum(is_coloc_pass(egfr_new$PP.H4, egfr_new$PP.H3), na.rm = TRUE)))
  cat(sprintf("PP.H4 >= 0.80 (H4 alone, diagnostic): %d\n", sum(egfr_new$PP.H4 >= 0.80)))
  cat(sprintf("PP.H4 >= 0.50: %d\n", sum(egfr_new$PP.H4 >= 0.50)))
  cat(sprintf("PP.H3 >= 0.80 (independent): %d\n", sum(egfr_new$PP.H3 >= 0.80)))

  # Top hits
  cat("\n-- Top eGFR coloc hits --\n")
  print(egfr_new[order(-PP.H4)][1:min(15, .N),
    .(gene, PP.H4, PP.H3, PP.H1, n_snps_coloc, coloc_support)])
} else {
  egfr_new <- data.table()
  cat("No eGFR coloc results obtained\n")
}

# 2026-06-30 gate fix #4b: do NOT silently save an empty eGFR coloc. The per-pair tryCatch swallows
# the n_total_sum stop(), and 08 only checks file existence (not non-empty) -> an empty result would
# propagate into v6. Fail loud on the full (non-chunked) run.
if (!.CHUNKED && (is.null(egfr_new) || nrow(egfr_new) == 0)) {
  source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_guard.R"))
  rerun_abort("06 eGFR coloc produced 0 results (n_total_sum/data issue, or no eGFR suggestive pairs) — refusing to save empty coloc_egfr_results.rds")
}

# ── Save eGFR coloc results to corrected workspace (combine with non-eGFR v5 later; live untouched) ──
.out_f <- if (.CHUNKED) file.path(out_dir, sprintf("coloc_egfr_chunk%02d.rds", .CI)) else file.path(out_dir, "coloc_egfr_results.rds")
# THE KEY TRAVELS WITH THE RESULT, on chunks and on the serial output alike -- 06 recomputes it and
# compares, which is the only thing that binds these numbers to the inputs they were computed from.
attr(egfr_new, "coloc_inputs_key") <- .EG_PROV_KEY
saveRDS(egfr_new, .out_f)
cat(sprintf("\n  Saved eGFR coloc: %d rows -> %s\n", nrow(egfr_new), .out_f))
if (.CHUNKED) {   # 2026-07-01: persist this chunk's per-pair status so the parallel MERGE can prove FULL coverage (results alone can't show a pair was dropped)
  fwrite(data.table(pair = names(pair_status), status = unlist(pair_status, use.names = FALSE)),
         file.path(out_dir, sprintf("coloc_egfr_chunk%02d_status.csv", .CI)))
  cat(sprintf("  chunk %02d pair-status: %d pairs -> coloc_egfr_chunk%02d_status.csv\n", .CI, length(pair_status), .CI))
}

# ── Strong hits assessment ───────────────────────────────────────────────────
cat("\n── Cascade Assessment ──\n")
if (nrow(egfr_new) > 0) {
  strong <- egfr_new[is_coloc_pass(PP.H4, PP.H3)]   # dual coloc rule (D1; centralized 2026-07-05)
  if (nrow(strong) > 0) {
    cat(sprintf("*** %d eGFR strong coloc hits (PP.H4 >= 0.80 & PP.H4/PP.H3 >= 5) ***\n", nrow(strong)))
    cat("These are NEW strong hits requiring Phase 4 annotation:\n")
    print(strong[, .(gene, PP.H4, PP.H3, n_snps_coloc)])
    cat("\n*** Phase 4-5-6-7 cascade needed for these proteins ***\n")
  } else {
    cat("No eGFR strong coloc hits (PP.H4 >= 0.80 & PP.H4/PP.H3 >= 5)\n")
    cat("No cascade to Phase 4 needed for eGFR\n")
  }
}

# Clean up checkpoint
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

cat("\n=== Phase 3-E eGFR Coloc Re-run Complete ===\n")
