# Parallel coloc worker: processes ONE chunk of the current generation's non-eGFR suggestive pairs.
# Identical methodology to 04_coloc.R / R/coloc.R::main()
# (coloc is per-locus independent -> union of chunks == one sequential run).
# Args: <chunk_rds> <out_subdir_name>
# HASHED BEFORE THEY ARE SOURCED, and re-checked afterwards. The key was built after the engine and
# the configuration had already been loaded, so it described whatever was on disk at that moment: a
# worker that loaded engine A and then had engine B installed under it stamped B's hash while main()
# still called A's functions -- including A's coloc priors, which the key parser reads back out of
# the file rather than out of memory.
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
.SRC_FILES <- c(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"),
                paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"),
                paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"),
                paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_coloc_inputs_key.R"))
.SRC_MD5_AT_LOAD <- vapply(.SRC_FILES, function(p)
  if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_, character(1))

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
COLOC_SOURCE_ONLY <- TRUE
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/coloc.R"))

.assert_src_unchanged <- function(when) {
  now <- vapply(.SRC_FILES, function(p)
    if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_, character(1))
  d <- which(now != .SRC_MD5_AT_LOAD | is.na(now) != is.na(.SRC_MD5_AT_LOAD))
  if (length(d))
    stop("22 worker: ", paste(basename(.SRC_FILES[d]), collapse = ", "), " changed on disk ", when,
         ". This process is still running the version it loaded at startup, so its fingerprint ",
         "would name code that did not produce these results.")
}
.assert_src_unchanged("between loading the engine and fingerprinting the inputs")
suppressPackageStartupMessages(library(data.table))

# gene coords for the FULL current-generation set (main() reads PHASE1_DIR/ukb_cis_instruments.rds;
# the canonical 1635 dir lacks the 47 backfilled proteins' coords -> they'd get NA region).
PHASE1_DIR <- file.path(RESULTS_DIR, PHASE1_SET)
stopifnot("need the pinned Phase 1 set (see R/_phase1_set.R): ukb_cis_instruments.rds" = file.exists(file.path(PHASE1_DIR, "ukb_cis_instruments.rds")))

args <- commandArgs(trailingOnly = TRUE)
chunk_file <- args[1]
out_sub    <- args[2]
REPRO <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
# out_sub BUILDS A PATH, so it is checked before it is used. Passed "..", PHASE3_DIR resolved to
# results/repro_corrected and main() wrote one chunk's results and status straight to the canonical
# manuscript paths; the worker then printed DONE.
if (!grepl("^out_[0-9]+$", out_sub))
  stop("22 worker: out_sub must be out_<k> (got ", sQuote(out_sub), "); it names a directory under ",
       "par_chunks and anything else can write outside it.")
if (!identical(basename(chunk_file), sprintf("chunk_%s.rds", sub("^out_", "", out_sub))))
  stop("22 worker: ", basename(chunk_file), " does not correspond to ", out_sub,
       "; the output directory and the work would describe different chunks.")
PHASE3_DIR <- file.path(REPRO, "par_chunks", out_sub)
dir.create(PHASE3_DIR, showWarnings = FALSE, recursive = TRUE)
# CLAIMED BEFORE ANY SLOW WORK, not after the fingerprint. The marker used to appear only once the
# input key had been built and the manifest validated -- seconds to minutes later -- so a splitter
# checking in that window saw no worker, and went on to recreate the directory this process was
# about to write into. 21 re-checks immediately before it deletes anything.
.ALIVE <- file.path(PHASE3_DIR, paste0("_worker_alive_", Sys.getpid()))
writeLines(as.character(Sys.getpid()), .ALIVE)
reg.finalizer(environment(), function(e) unlink(.ALIVE), onexit = TRUE)

# READ ONCE, AT THE START, and the pairs, the identity and the hash all come from THOSE BYTES.
# readRDS() followed by md5sum() is two reads, so a chunk replaced in between left the stamp
# describing bytes this worker never used; unserialising the raw vector removes the gap.
.CHUNK_RAW  <- readBin(chunk_file, "raw", file.size(chunk_file))
.CHUNK      <- unserialize(memDecompress(.CHUNK_RAW, type = "gzip"))
.CHUNK_MD5  <- digest::digest(.CHUNK_RAW, algo = "md5", serialize = FALSE)
.CHUNK_RUN  <- attr(.CHUNK, "run_id")
if (is.null(.CHUNK_RUN))
  stop("22 worker: ", chunk_file, " carries no run id. It was written by a splitter that predates ",
       "run identity, so this work cannot be tied to a split. Re-run 21_coloc_split.R.")
# The inputs are fingerprinted HERE too, for the same reason: taken at the end they would describe
# whatever a concurrent splitter had just installed rather than what this run actually used.
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_coloc_inputs_key.R"))
# The scope comes from the manifest, so this worker fingerprints the same material the split did.
# Computed independently it would cover only this chunk's genes and the keys could never match.
.MANI_F0 <- file.path(REPRO, "par_chunks", "_split_manifest.rds")
if (!file.exists(.MANI_F0))
  stop("22 worker: no split manifest at ", .MANI_F0, " -- run 21_coloc_split.R first.")
.STAMP_IK <- coloc_inputs_key(readRDS(.MANI_F0)$key_scope)

# ── CHECKED BEFORE main() RUNS, not after ─────────────────────────────────────────────────────
# The end-of-run check came too late to matter: by then this worker had already deposited
# checkpoints into the directory, and a worker belonging to the NEXT split would resume them. The
# engine's own checkpoint provenance records only outcome sample sizes and the coordinate hash
# (R/coloc.R:1019), so it cannot tell one split from another.
.MANI_F <- file.path(REPRO, "par_chunks", "_split_manifest.rds")
if (!file.exists(.MANI_F))
  stop("22 worker: no split manifest at ", .MANI_F, " -- run 21_coloc_split.R first.")
.MANI <- readRDS(.MANI_F)
if (!identical(.CHUNK_RUN, .MANI$run_id %||% NA_character_))
  stop("22 worker ", out_sub, ": this chunk belongs to split ", .CHUNK_RUN, " but the manifest in ",
       "par_chunks is for split ", .MANI$run_id %||% "<none>", ". A new split has taken over this ",
       "directory; nothing was computed.")
if (!identical(.STAMP_IK$key, .MANI$inputs_key %||% NA_character_))
  stop("22 worker ", out_sub, ": the inputs differ from those the split pinned (worker ",
       substr(.STAMP_IK$key, 1, 12), ", split ", substr(.MANI$inputs_key %||% "<none>", 1, 12),
       "):\n", paste(coloc_inputs_diff(.MANI$inputs_detail %||% list(), .STAMP_IK$detail),
                     collapse = "\n"), "\nNothing was computed.")

# A checkpoint in this directory has to belong to THIS split, or it is another run's work waiting to
# be resumed as if it were ours. The marker is written before the engine can create one.
.RUNMARK <- file.path(PHASE3_DIR, "_worker_run_id.txt")
.CKPT    <- file.path(PHASE3_DIR, "coloc_checkpoint_v7.rds")
.prev_run <- if (file.exists(.RUNMARK)) readLines(.RUNMARK, warn = FALSE)[1] else NA_character_
if (file.exists(.CKPT) && !identical(.prev_run, .CHUNK_RUN)) {
  unlink(.CKPT)
  if (file.exists(.CKPT))
    stop("22 worker ", out_sub, ": a checkpoint from split ",
         if (is.na(.prev_run)) "<unmarked>" else .prev_run,
         " is present and could not be removed. Resuming it would publish that split's results ",
         "under this one's provenance.")
  cat(sprintf("[worker %s] discarded a checkpoint left by split %s\n", out_sub,
              if (is.na(.prev_run)) "<unmarked>" else .prev_run))
}
writeLines(.CHUNK_RUN, .RUNMARK)

cat(sprintf("[worker %s] split %s, inputs %s\n", out_sub, .CHUNK_RUN, substr(.STAMP_IK$key, 1, 12)))

load_phase2_suggestive <- function() {
  p <- copy(.CHUNK); setDT(p)
  cat(sprintf("[worker %s] chunk pairs: %d\n", out_sub, nrow(p))); p
}

cat(sprintf("=== coloc worker %s (chunk=%s) ===\n", out_sub, chunk_file))
invisible(main())

# ── Completion stamp ────────────────────────────────────────────────────────────────────────
# Written LAST, after main() returns, and it is what the merge requires as proof this worker
# finished. Two things were previously unprovable at merge time:
#
#   - WHICH inputs this chunk was computed from. Only the coordinate table was recorded, via a
#     substring of the engine's provenance key, so a chunk re-run after the engine, the outcome
#     GWAS or the instrument table changed merged silently into results from before the change.
#
#   - WHETHER the status and result files being merged are the ones this run produced. The merge
#     read the checkpoint for its key and then read the CSV and the RDS separately, so a
#     correct-provenance checkpoint could bless stale files left behind by an earlier attempt --
#     the worker does not delete its previous outputs first, which makes that the ordinary case
#     for a re-run chunk, not an exotic one.
#
# Hashing the two files here binds them to this run: at merge time they either still hash to what
# this worker wrote, or they are not this worker's output.
.stamp_h   <- function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
# The inputs must not have moved under this worker while it ran. Taken only at the start, a change
# during the run would go unrecorded; taken only at the end, it would be recorded as what was used.
.ik_end <- coloc_inputs_key(.MANI$key_scope)
if (!identical(.ik_end$key, .STAMP_IK$key)) {
  # QUARANTINE FIRST. The engine's checkpoint provenance carries only its version literal, the
  # outcome sample sizes and the coordinate hash -- not the input key -- so a checkpoint written
  # while the inputs changed is indistinguishable from a clean one. Left in place, restoring the
  # original input and rerunning the same chunk reuses it: the completed pairs are skipped, both
  # fingerprints now agree, and the mixed results are stamped as clean.
  for (.q in c("coloc_checkpoint_v7.rds", "coloc_pair_status_v7.csv", "coloc_results_v7_stage1.rds"))
    if (file.exists(file.path(PHASE3_DIR, .q)))
      file.rename(file.path(PHASE3_DIR, .q),
                  file.path(PHASE3_DIR, paste0(.q, ".mixed_inputs_", Sys.getpid())))
  unlink(.RUNMARK)
}
if (!identical(.ik_end$key, .STAMP_IK$key))
  stop("22 worker ", out_sub, ": the inputs changed while this chunk was running (start ",
       substr(.STAMP_IK$key, 1, 12), ", end ", substr(.ik_end$key, 1, 12), "):\n",
       paste(coloc_inputs_diff(.STAMP_IK$detail, .ik_end$detail), collapse = "\n"),
       "\nThese results were computed against a mixture. Nothing was stamped; re-run this chunk.")
.assert_src_unchanged("while this chunk was running")
if (!identical(digest::digest(readBin(chunk_file, "raw", file.size(chunk_file)),
                              algo = "md5", serialize = FALSE), .CHUNK_MD5))
  stop("22 worker ", out_sub, ": ", chunk_file, " was replaced while this chunk was running. ",
       "A new split has taken over this directory; these results belong to the previous one.")
.stamp_csv <- file.path(PHASE3_DIR, "coloc_pair_status_v7.csv")
.stamp_rds <- file.path(PHASE3_DIR, "coloc_results_v7_stage1.rds")
saveRDS(list(out_sub      = out_sub,
             chunk_file   = normalizePath(chunk_file, winslash = "/"),
             chunk_md5    = .CHUNK_MD5,      # the bytes this worker actually read
             run_id       = .CHUNK_RUN,      # from that same read, not a later one
             inputs_key   = .STAMP_IK$key,
             inputs_detail= .STAMP_IK$detail,
             status_md5   = .stamp_h(.stamp_csv),
             results_md5  = .stamp_h(.stamp_rds),   # NA is legitimate: an all-skip chunk writes none
             finished_at  = format(Sys.time(), tz = "UTC"),
             r_version    = R.version.string),
        file.path(PHASE3_DIR, "_worker_complete.rds"))
cat(sprintf("=== worker %s DONE (inputs %s) ===\n", out_sub, substr(.STAMP_IK$key, 1, 12)))
