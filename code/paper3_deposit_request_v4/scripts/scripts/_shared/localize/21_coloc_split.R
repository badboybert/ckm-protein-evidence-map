# ============================================================================
# Parallel-coloc SPLITTER (REDESIGNED 2026-07-01 around the reviewed serial).
# Partitions the SERIAL wrapper's EXACT non-eGFR suggestive pair set (same filter
# as scripts/paper3.1_ckm/cascade/04_coloc.R load_phase2_suggestive) into NCHUNK
# deterministic chunks. Each worker (22) then runs the REVIEWED R/coloc.R main()
# on its chunk, so the coloc methodology + per-pair status + error-abort + prov_key
# checkpoint are the VETTED serial ones (union of chunks == one serial run).
#
# NO cross-run seeding (the old seeding laundered un-validated stale results and lost
# pair-status -> removed): each parallel run is FRESH; intra-run resume is handled by
# coloc's own per-chunk checkpoint. Old chunk outputs are CLEARED first so a
# 0-result chunk cannot leave a prior coloc_results_v7_stage1.rds behind (no stale merge).
# After this: 22_coloc_chunk_worker.R on each chunk; then 23_merge_parallel_coloc.R.
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
suppressPackageStartupMessages(library(data.table))
setwd(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"))
# HASHED BEFORE THEY ARE LOADED, and re-checked before the manifest is written. The worker already
# does this; the splitter did not, so it could build the pair universe from config A -- including
# A's eGFR id, which the merge later filters on -- while its input key hashed config B.
.SRC_FILES <- c("R/00_config.R", "R/_phase1_set.R", "R/_coloc_inputs_key.R", "R/coloc.R",
                "scripts/_shared/localize/21_coloc_split.R",
                "scripts/_shared/localize/22_coloc_chunk_worker.R",
                "scripts/_shared/localize/23_merge_parallel_coloc.R")
.SRC_MD5_AT_LOAD <- vapply(.SRC_FILES, function(p)
  if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_, character(1))
source("R/00_config.R")
source("R/_phase1_set.R"); source("R/_coloc_inputs_key.R")
REPRO   <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG)); PC <- file.path(REPRO, "par_chunks")
# 3 workers. Was 5 (REBUILD_LEDGER §9ak), measured against the 2,629,657-entry rsID map. The map
# was rebuilt on 2026-08-07 to 6,115,599 entries and each worker deduplicates it independently.
# Measured on this machine (31.7 GB): 5 workers -> 31.1 GB resident, 3.8 GB paged, no worker
# finished loading in 28 minutes; 4 workers -> 31.3 GB resident, 2.9 GB paged, same outcome.
# Each worker holds roughly 7 GB once the map is deduplicated, so 3 is what fits.
NCHUNK  <- 3L
EGFR_ID <- unname(OUTCOMES[["eGFR"]])   # config single source (= "Wuttke2019_eGFR_EA"); no hardcoded id

# CHECKED BEFORE ANYTHING IS CLEARED. The split used to accept any readable coordinate TSV while
# every worker requires that table's provenance certificate (R/coloc.R:339). With the certificate
# absent -- which is the state on disk -- 21 would delete all existing chunks, announce success, and
# every worker would then refuse to run. Whatever the workers demand, the split must demand first.
.cert <- paste0(GENE_HG38_COORDS, ".manifest.tsv")
# EXACTLY "1", which is what R/coloc.R:339 requires. nzchar() accepted "0", "false" and anything
# else non-empty as permission, so the one spelling a user would reach for to DISABLE the override
# enabled it -- and the workers, which do check for "1", would then refuse everything this split
# had just cleared the way for.
if (!identical(Sys.getenv("CKM_ALLOW_UNCERTIFIED_COORDS"), "1")) {
  # The SAME test the engine applies (R/coloc.R:352-372), not file.exists(). A stale, malformed,
  # wrong-builder or wrong-md5 certificate passed existence, so the split cleared every chunk and
  # moved the canonical outputs aside before each worker rejected the very same file.
  .why <- NULL
  if (!file.exists(.cert)) {
    .why <- paste0("there is no certificate (", basename(.cert), ") beside it")
  } else {
    .cm <- tryCatch(fread(.cert), error = function(e) NULL)
    .have <- digest::digest(readBin(GENE_HG38_COORDS, "raw", file.size(GENE_HG38_COORDS)),
                            algo = "md5", serialize = FALSE)
    if (is.null(.cm) || !all(c("file", "md5", "builder") %in% names(.cm)) || nrow(.cm) != 1L)
      .why <- "its certificate is unreadable or malformed"
    else if (!identical(tolower(.have), tolower(.cm$md5[1])))
      .why <- "its certificate describes a DIFFERENT file (md5 mismatch)"
    else if (!identical(.cm$builder[1], "03b_gene_coords_hg38.R"))
      .why <- paste0("its certificate names builder '", .cm$builder[1],
                     "', not 03b_gene_coords_hg38.R")
    else if (!"builder_md5" %in% names(.cm) || is.na(.cm$builder_md5[1]) ||
             !nzchar(.cm$builder_md5[1]))
      .why <- "its certificate does not record which version of 03b_gene_coords_hg38.R wrote it"
    else {
      # THE ENGINE'S OWN TEST, not a second copy of it. A nonempty builder_md5 is not the same as a
      # builder_md5 that matches the builder on disk, and the engine additionally validates
      # source_md5, the upstream builder, the GENCODE GTF, the UniProt mapping and the release
      # listing (R/coloc.R:339). Anything the workers will reject, this must reject BEFORE it clears
      # their chunks and moves the canonical outputs aside -- and calling their validator is the
      # only way that stays true when the engine's rules change.
      # FORCED TO TRUE, whatever was there. Testing exists() alone treated an existing FALSE as
      # "already set" and left it false, so sourcing the engine ran main() -- a full sequential
      # colocalisation, writing checkpoints and results, from inside a preflight.
      .prev <- if (exists("COLOC_SOURCE_ONLY", envir = globalenv()))
                 get("COLOC_SOURCE_ONLY", envir = globalenv()) else NULL
      assign("COLOC_SOURCE_ONLY", TRUE, envir = globalenv())
      .why <- tryCatch({
        sys.source(file.path(BASE_DIR, "R/coloc.R"), envir = globalenv())
        if (!exists(".assert_coords_certified", envir = globalenv()))
          "R/coloc.R no longer defines .assert_coords_certified(), so this preflight cannot mirror it"
        else { get(".assert_coords_certified", envir = globalenv())(); NULL }
      }, error = function(e) paste0("the engine refuses it: ", conditionMessage(e)))
      if (is.null(.prev)) rm("COLOC_SOURCE_ONLY", envir = globalenv()) else
        assign("COLOC_SOURCE_ONLY", .prev, envir = globalenv())
    }
  }
  if (!is.null(.why))
    stop("21 split: ", .why, ". Every worker applies this same test, so splitting would clear the ",
         "existing chunks and move the canonical outputs aside, and then nothing would be able to ",
         "run. Rebuild via 03b_gene_coords_hg38.R, or set CKM_ALLOW_UNCERTIFIED_COORDS=1 ",
         "deliberately -- which the workers also honour.")
}

# A LIVE WORKER MAKES RE-SPLITTING IMPOSSIBLE, rather than detectable afterwards. The chunk
# directories have fixed names, so a worker still running from the previous split keeps writing into
# out_k after this script recreates it -- its checkpoints land in the new run's directory and the
# next worker resumes them. Three separate merge-time checks were added to catch that; refusing to
# start is one check, and it removes the situation instead of describing it.
# A FUNCTION, because it is asked twice: once here so an obviously blocked split fails before doing
# any work, and again immediately before the directories are deleted. The worker claims its marker
# before it fingerprints anything, but the window between this check and that deletion is still a
# window, and it is the deletion that does the damage.
.assert_no_live_worker <- function() {
  .live <- list.files(PC, pattern = "^_worker_alive_", recursive = TRUE, full.names = TRUE)
  .live <- .live[vapply(.live, function(f) {
    pid <- suppressWarnings(as.integer(readLines(f, warn = FALSE)[1]))
    if (is.na(pid)) return(TRUE)           # unreadable marker: assume alive rather than guess
    # tasklist is the only reliable liveness test on this platform; a lookup failure is treated as
    # "alive" so a split never proceeds on the strength of a check that did not work.
    out <- suppressWarnings(tryCatch(
      system2("tasklist", c("/FI", shQuote(sprintf("PID eq %d", pid))), stdout = TRUE, stderr = FALSE),
      error = function(e) NULL))
    # A LOOKUP THAT DID NOT WORK MEANS ALIVE. The handler returned the string "ERROR" and the test
    # then asked whether it contained the pid -- it does not -- so a failed tasklist was read as
    # "dead", and the split proceeded to recreate directories a live worker was writing into. The
    # comment claimed the opposite of what the code did.
    if (is.null(out) || !length(out) ||
        (!is.null(attr(out, "status")) && !identical(as.integer(attr(out, "status")), 0L)))
      return(TRUE)
    any(grepl(as.character(pid), out, fixed = TRUE))
  }, logical(1))]
  if (length(.live))
    stop("21 split: ", length(.live), " coloc worker(s) are still running (",
         paste(basename(.live), collapse = ", "), "). Re-splitting now would recreate the ",
         "directories they are writing into, and their output would be merged as this split's. ",
         "Wait for them, or stop them and delete the markers.")
}
.assert_no_live_worker()

# A unique identity for THIS split. The directory names are fixed (chunk_k.rds, out_k), so clearing
# and recreating them gives a new run the same paths as the old one. A worker from a previous run
# that is still alive -- between checkpoint writes -- then writes into the new run's out_k, and if
# the pair set is unchanged the merge accepts it as one of its own. Every worker carries this id
# from the chunk file it was given, so such a worker now identifies itself as belonging elsewhere.
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
                 paste(sample(c(0:9, letters[1:6]), 8, replace = TRUE), collapse = ""))
cat(sprintf("split run id: %s\n", RUN_ID))

# EXACT serial non-eGFR universe: mirror 04_coloc.R (exclude BOTH outcome_name=="eGFR" AND the eGFR id),
# so the parallel pair set is byte-for-byte the serial one (incl. excluding any legacy Pattaro eGFR pairs).
sug <- readRDS(file.path(REPRO, "mr_suggestive.rds")); setDT(sug)
# The filter below drops on outcome_name AND outcome_id, so a row whose NAME is missing or stale
# is discarded even when its ID names a perfectly valid outcome the engine would have evaluated --
# and 23 applies the identical filter, so exact-coverage validation then passes against the already
# truncated universe and the pair is simply absent from the final table. Asserting the mapping first
# makes a label problem stop the run instead of quietly shrinking it.
# ... on the rows the METHOD FILTER keeps. Checked over the whole table, an auxiliary row carrying an
# ignored method and an outcome outside OUTCOMES aborted the split, although neither this script nor
# R/coloc.R::main() would ever have processed it. The rows that can actually be lost are the ones the
# eGFR filter below is about to test, and those are exactly these.
.truth <- data.table(outcome_name = names(OUTCOMES), outcome_id = unname(unlist(OUTCOMES)))
.seen  <- unique(sug[method %in% c("Inverse variance weighted", "Wald ratio"),
                     .(outcome_name, outcome_id)])
# ...AND IT ASSERTS THE FIELD THE FILTER USES. Requiring the (name, id) PAIR to match made this fatal
# on a row whose canonical id is perfectly valid but whose display name is stale or missing -- while
# the split below, 04, 23 and 06 all stratify by id alone, so serial 04 processes that row and the
# default parallel path aborts before splitting. The id is canonical and the name is derived from it,
# so the id is what has to be in OUTCOMES; a name that disagrees is a labelling problem, reported
# rather than fatal, because nothing downstream reads it to decide anything.
.badm  <- .seen[!.truth, on = .(outcome_id)]
if (nrow(.badm))
  stop("21 split: ", nrow(.badm), " outcome_id(s) in mr_suggestive.rds are not in the config's ",
       "OUTCOMES: ", paste(utils::head(unique(.badm$outcome_id), 5), collapse = ", "),
       ". The split stratifies by id, so an id the config does not know cannot be assigned to a ",
       "stratum or reconstructed by the merge's coverage check.")
.badn <- .seen[.truth, on = .(outcome_id)][!is.na(outcome_name) & outcome_name != i.outcome_name |
                                           is.na(outcome_name)]
if (nrow(.badn))
  warning("21 split: ", nrow(.badn), " row(s) carry an outcome_name that disagrees with the config's ",
          "name for their id (", paste(utils::head(.badn[, paste0(outcome_id, ": ", outcome_name,
          " vs ", i.outcome_name)], 3), collapse = "; "), "). The id is canonical and nothing ",
          "downstream reads the name to decide anything, so this does not stop the split.",
          immediate. = TRUE, call. = FALSE)

pu <- unique(sug[method %in% c("Inverse variance weighted","Wald ratio") &
                 !(outcome_id %in% EGFR_ID),
                 .(gene, outcome_id, outcome_name)])
pu[, key := paste(gene, outcome_id, sep = "_")]
# The pair key is two fields joined by "_", so it is only decodable, and only collision-free,
# while the gene side contains no "_". That holds for every gene in the current data and is
# asserted rather than assumed, because a future symbol containing one would make two different
# pairs share a key and the merge's coverage checks would compare the wrong things.
if (any(grepl("_", pu$gene, fixed = TRUE)))
  stop("21 split: ", sum(grepl("_", pu$gene, fixed = TRUE)), " gene symbol(s) contain '_' (",
       paste(head(unique(pu$gene[grepl("_", pu$gene, fixed = TRUE)]), 3), collapse = ", "),
       "), so gene_outcome keys are ambiguous. Change the separator before proceeding.")
# A blank gene or outcome makes a key like "_FinnGen_R10_I9_AF" -- not NA, not duplicated, and it
# passes every check while naming no protein. The FIELDS are tested, not only the key built from
# them.
if (anyNA(pu$gene) || !all(nzchar(trimws(pu$gene))) ||
    anyNA(pu$outcome_id) || !all(nzchar(trimws(pu$outcome_id))))
  stop("21 split: ", sum(is.na(pu$gene) | !nzchar(trimws(pu$gene))), " pair(s) have no gene and ",
       sum(is.na(pu$outcome_id) | !nzchar(trimws(pu$outcome_id))), " have no outcome id. A blank ",
       "field still produces a well-formed key, so every coverage check downstream would pass on a ",
       "pair that names nothing.")
if (anyNA(pu$key) || anyDuplicated(pu$key))
  stop("21 split: the pair set has missing or duplicated keys, which every downstream coverage ",
       "check assumes cannot happen.")
setorder(pu, key)
# Never more chunks than pairs. With six chunks and, say, four pairs, two workers receive an empty
# chunk and the engine cannot complete on nothing -- so the merge waits for output that will never
# arrive. Only reachable on a reduced pair set, which is exactly when someone is debugging.
if (nrow(pu) < NCHUNK) {
  cat(sprintf("only %d pair(s): using %d chunk(s) instead of %d\n", nrow(pu), nrow(pu), NCHUNK))
  # Zero pairs is not one empty chunk. `max(1L, 0)` contradicted the line above it and produced a
  # worker with nothing to do: the engine's empty pair_status has no pair or status column, its
  # status == "error" test cannot evaluate, and no completion stamp is written -- so a legitimately
  # empty analysis becomes an unmergeable one. If there is nothing to colocalise, that is the
  # result, and it is said here rather than discovered six processes later.
  if (nrow(pu) == 0L)
    stop("21 split: the non-eGFR suggestive set is empty, so there is nothing to colocalise. ",
         "That may be correct, but it is not something this pipeline can express as a split: ",
         "run the merge against an empty universe deliberately instead.")
  NCHUNK <- max(1L, nrow(pu))
}
pu[, chunk := (seq_len(.N) - 1L) %% NCHUNK + 1L]
cat(sprintf("full non-eGFR suggestive pairs: %d -> %d deterministic chunks\n", nrow(pu), NCHUNK))

# Manifest = SINGLE SOURCE of NCHUNK + the expected pair set, so 23 hardcodes neither and can assert
# EXACT coverage (no missing, no extra) against what was actually split.
# The COORDINATE TABLE is part of the split, not just the pair set. Each worker pins the table it
# happens to read at startup, and the merge combined their results without comparing those pins --
# so 03b installing a new table between two worker launches, or one failed chunk being rerun after
# a rebuild, produced one result set whose cis windows came from two different tables while every
# coverage check passed.
.coord_md5 <- if (file.exists(GENE_HG38_COORDS))
                unname(tools::md5sum(GENE_HG38_COORDS)) else NA_character_
if (is.na(.coord_md5))
  stop("cannot hash ", GENE_HG38_COORDS, ", so the workers cannot be bound to one coordinate table.")
cat(sprintf("  coordinate table pinned for all workers: %s\n", substr(.coord_md5, 1, 12)))
# ... and so is EVERYTHING ELSE the answer depends on. The coordinate table was the only input
# pinned here, so the engine, the configuration, the instrument table, the rsID map, the liftover
# chain, the outcome GWAS and the coloc priors could all be replaced between two chunks, or between
# the chunks and the merge, with every coverage check still passing.
# The key describes THIS run's data: the archives of the genes in the pair set and the outcomes it
# looks up -- not all 2,940 archives and all eight configured outcomes. The scope goes into the
# manifest so the workers and the merge fingerprint exactly the same material.
.scope <- list(genes = sort(unique(pu$gene)), outcome_ids = sort(unique(pu$outcome_id)))
# The sources must still be the ones this process loaded, or the key would name code that did not
# build this split.
.now_md5 <- vapply(.SRC_FILES, function(p)
  if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_, character(1))
.chg <- which(.now_md5 != .SRC_MD5_AT_LOAD | is.na(.now_md5) != is.na(.SRC_MD5_AT_LOAD))
if (length(.chg))
  stop("21 split: ", paste(basename(.SRC_FILES[.chg]), collapse = ", "), " changed on disk while ",
       "this split was being built. The pair universe came from the version loaded at startup, so ",
       "the fingerprint would name different code. Nothing was written.")
.ik <- coloc_inputs_key(.scope)
cat(sprintf("  inputs pinned for all workers: %s\n", substr(.ik$key, 1, 12)))
for (.n in names(.ik$detail$hashed))
  cat(sprintf("    %-14s %s\n", .n, substr(.ik$detail$hashed[[.n]], 1, 12)))

# ── NOTHING ABOVE THIS LINE WRITES OR DELETES ─────────────────────────────────────────────────
# The previous runnable state used to be destroyed first: the chunks were cleared and the canonical
# outputs moved aside, and only then did the script read the pair universe and build the key. An
# empty universe, an unresolvable archive or any error in the key therefore left the tree wiped and
# the manuscript's tables displaced, with nothing to run -- and the source-change check said
# "Nothing was written" while saying it. Every reason to refuse is now established first, so the
# destructive part is reached only once a complete split is known to be constructible.
.assert_no_live_worker()   # re-asked HERE: this is the moment the deletion happens
# FRESH: clear any stale chunk dirs/outputs/checkpoints so nothing from a prior run survives.
# 2026-07-01 (Codex redesign #2): force=TRUE + VERIFY the clear actually succeeded (Windows can leave
# locked/read-only files that unlink silently skips -> a stale chunk could then be resumed by coloc).
if (dir.exists(PC)) {
  unlink(PC, recursive = TRUE, force = TRUE)
  if (dir.exists(PC) && length(list.files(PC, recursive = TRUE)) > 0)
    stop(sprintf("21 split: could NOT clear stale par_chunks (%s) -- locked/read-only files remain; a stale chunk output/checkpoint could be resumed. Remove it manually and retry.", PC))
}
dir.create(PC, showWarnings = FALSE, recursive = TRUE)

# The root outputs the manuscript reads are NOT part of par_chunks, so clearing the chunks left the
# previous run's merged table in place looking current. If this split's merge later aborts, that
# older table is what downstream code finds. It is moved aside rather than deleted -- nothing is
# lost, but it can no longer be mistaken for the result of this run.
# OUTSIDE par_chunks. It used to live INSIDE it, and the next split's first act is to delete that
# directory recursively -- so "nothing is lost" held for exactly one further run and then silently
# stopped holding. Kept beside the canonical outputs instead, where nothing clears it.
.super <- file.path(REPRO, paste0("_superseded_", RUN_ID))
.moved <- character(0)
for (.f in c("coloc_results_v7_stage1.rds", "coloc_pair_status_v7.csv")) {
  .p <- file.path(REPRO, .f)
  if (file.exists(.p)) {
    dir.create(.super, showWarnings = FALSE, recursive = TRUE)
    if (!file.rename(.p, file.path(.super, .f)))
      # ...AND IT SAYS WHERE THE ALREADY-MOVED FILE IS. The result and its status file are a pair; if
      # the first move succeeds and the second fails, aborting here leaves them split across two
      # directories with nothing recording that, and the obvious recovery (re-run) then finds only
      # half of the previous generation at the canonical path.
      stop("21 split: could not move the previous ", .f, " aside. It would otherwise remain at the ",
           "canonical path and be read as this run's result if the merge does not complete.",
           if (length(.moved)) paste0(" ALREADY MOVED to ", .super, ": ",
                                      paste(.moved, collapse = ", "),
                                      " -- move them back before retrying.") else "")
    .moved <- c(.moved, .f)
    cat(sprintf("  previous %s moved to %s\n", .f, basename(.super)))
  }
}

for (k in 1:NCHUNK) {
  dir.create(file.path(PC, sprintf("out_%d", k)), showWarnings = FALSE, recursive = TRUE)
  ck <- pu[chunk == k]
  # The run id travels WITH the chunk, as an attribute, so a worker cannot be launched against this
  # run's directory while holding a different run's work: it reads the id from the pairs it was
  # actually given. An attribute rather than a column, so load_phase2_suggestive()'s data frame is
  # unchanged and the engine sees exactly what it saw before.
  .ckd <- ck[, .(gene, outcome_id, outcome_name, method = "Inverse variance weighted")]
  attr(.ckd, "run_id") <- RUN_ID
  saveRDS(.ckd, file.path(PC, sprintf("chunk_%d.rds", k)))
  cat(sprintf("  chunk %d: %d pairs\n", k, nrow(ck)))
}
saveRDS(list(nchunk = NCHUNK, expected_keys = pu$key, egfr_id = EGFR_ID,
             coord_md5 = .coord_md5, run_id = RUN_ID,
             inputs_key = .ik$key, inputs_detail = .ik$detail, key_scope = .scope,
             split_at = format(Sys.time(), tz = "UTC")),
        file.path(PC, "_split_manifest.rds"))
cat(sprintf("\nSplit done (NCHUNK=%d, %d pairs). Manifest written -> launch 22 workers, then 23_merge.\n",
            NCHUNK, nrow(pu)))
