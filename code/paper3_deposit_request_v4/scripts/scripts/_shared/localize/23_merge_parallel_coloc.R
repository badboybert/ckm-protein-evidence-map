# ============================================================================
# Parallel-coloc MERGE (REDESIGNED 2026-07-01 around the reviewed serial).
# Unions the NCHUNK worker outputs (each produced by the REVIEWED R/coloc.R main())
# into results/repro_corrected/coloc_results_v7_stage1.rds + coloc_pair_status_v7.csv, and
# asserts EXACT coverage against the split manifest (== one serial coloc run).
# NCHUNK + the expected pair set come from the manifest (21 wrote it) -> nothing hardcoded.
# Run AFTER all workers finish.
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
options(width = 200)
suppressPackageStartupMessages(library(data.table))
setwd(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")); source("R/00_config.R")
source("R/_phase1_set.R"); source("R/_coloc_inputs_key.R")
REPRO <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG)); PC <- file.path(REPRO, "par_chunks")

manf <- file.path(PC, "_split_manifest.rds")
if (!file.exists(manf)) stop("23 merge: split manifest missing -- run 21_coloc_split.R first")
mani   <- readRDS(manf)
NCHUNK <- mani$nchunk

# Every worker must have used the coordinate table the split pinned, and it must still be the one on
# disk. Without this the merge would combine chunks computed against different tables -- different
# cis windows for the same pair -- and report success.
if (is.null(mani$coord_md5))
  stop("23 merge: the split manifest predates coordinate-table pinning, so the chunks cannot be ",
       "shown to share one set of cis windows. Re-run 21_coloc_split.R and the workers.")
.now <- if (file.exists(GENE_HG38_COORDS)) unname(tools::md5sum(GENE_HG38_COORDS)) else NA_character_
if (is.na(.now) || !identical(tolower(.now), tolower(mani$coord_md5)))
  stop("23 merge: the coordinate table changed after the split (split ",
       substr(mani$coord_md5, 1, 12), ", now ", if (is.na(.now)) "absent" else substr(.now, 1, 12),
       "). Results computed before and after would use different cis windows.")

# The coordinate table was the ONLY input compared. Everything else that decides a coloc result --
# the engine, the configuration, the instrument table, the rsID map, the liftover chain, the outcome
# GWAS files, the coloc priors -- could change between the split and here, or between two chunks,
# and nothing noticed. The same key is rebuilt now and compared with what the split recorded.
if (is.null(mani$inputs_key))
  stop("23 merge: the split manifest predates input pinning, so these chunks cannot be shown to ",
       "share one set of inputs. Re-run 21_coloc_split.R and the workers.")
.ik_now <- coloc_inputs_key(mani$key_scope)
if (!identical(.ik_now$key, mani$inputs_key))
  stop("23 merge: the inputs changed between the split and this merge (split ",
       substr(mani$inputs_key, 1, 12), ", now ", substr(.ik_now$key, 1, 12), "):\n",
       paste(coloc_inputs_diff(mani$inputs_detail, .ik_now$detail), collapse = "\n"),
       "\nRe-run 21_coloc_split.R and the workers.")
cat(sprintf("inputs key %s matches the split\n", substr(mani$inputs_key, 1, 12)))
# WHAT THIS DOES NOT COVER, stated rather than implied. Each worker fingerprints its inputs before
# main() and again after, and the merge compares both against the split. That detects a change that
# PERSISTS. An input replaced and restored while a worker was running leaves all three fingerprints
# equal, and no fingerprint taken at intervals can do otherwise. Two things narrow it: the engine
# re-reads and re-checks the coordinate table per pair (R/coloc.R), and 21 now refuses to run while
# any worker is alive, so the ordinary way this happened -- a new split under a running worker --
# cannot occur. The residual is a deliberate mid-run substitution, which nothing here claims to
# detect.

wk <- list(); ps <- list()
for (k in 1:NCHUNK) {
  # ... and each worker's own checkpoint has to agree, so a chunk rerun after a rebuild is caught
  # even when the table has since been restored to what the split pinned.
  #
  # THREE things were wrong with the first attempt at this, and each of them made it a no-op:
  #   - it looked for coloc_checkpoint.rds; the engine writes coloc_checkpoint_v7.rds, so the block
  #     never ran for any real worker;
  #   - a missing checkpoint, or one with no prov_key, was accepted -- the two states in which the
  #     worker's coordinate table is precisely unknown;
  #   - it asked whether the hash appeared ANYWHERE in the key, so a malformed key carrying it
  #     somewhere else passed.
  # All three are now fatal, and the COORDS= field is parsed and compared exactly.
  .ckf <- file.path(PC, sprintf("out_%d", k), "coloc_checkpoint_v7.rds")
  if (!file.exists(.ckf))
    stop(sprintf("23 merge: chunk %d has no checkpoint (%s), so the coordinate table it used ",
                 k, basename(.ckf)),
         "cannot be established. Re-run that chunk.")
  .ck <- unserialize(memDecompress(readBin(.ckf, "raw", file.size(.ckf)), type = "gzip"))
  if (is.null(.ck$prov_key) || !nzchar(.ck$prov_key))
    stop(sprintf("23 merge: chunk %d's checkpoint records no provenance key, so the coordinate ", k),
         "table it used cannot be established. Re-run that chunk.")
  # Bounded on the RIGHT as well as the left. COORDS= used to be the last field, so taking everything
  # after it was the whole value; the key now continues "|INSTR=...|WIN=..." (R/coloc.R:1105-1106),
  # and without the second sub() this compared "<md5>|INSTR=...|WIN=500" to a bare md5 and every
  # merge failed. Keys that still end at COORDS= are unaffected: there is no "|" left to cut.
  .got <- sub("\\|.*$", "", sub(".*\\|COORDS=", "", .ck$prov_key))
  if (!grepl("\\|COORDS=", .ck$prov_key) || !identical(tolower(.got), tolower(mani$coord_md5)))
    stop(sprintf("23 merge: chunk %d was computed against coordinate table %s, but the split ", k,
                 if (nzchar(.got)) substr(.got, 1, 12) else "<none recorded>"),
         "pinned ", substr(mani$coord_md5, 1, 12), ". Their cis windows are not the same. ",
         "Re-run that chunk.")
  # The worker's completion stamp, written after main() returns. pair-status alone was treated as
  # proof the worker finished, but the engine writes it BEFORE saving results (R/coloc.R:1087-1102),
  # so its presence proves only that the worker got that far.
  .stf <- file.path(PC, sprintf("out_%d", k), "_worker_complete.rds")
  if (!file.exists(.stf))
    stop(sprintf("23 merge: chunk %d has no completion stamp (%s), so it cannot be shown to have ",
                 k, basename(.stf)),
         "run to the end or to have used the same inputs as the other chunks. Re-run that chunk.")
  .st <- unserialize(memDecompress(readBin(.stf, "raw", file.size(.stf)), type = "gzip"))
  # Which SPLIT this chunk belongs to. The chunk directories have fixed names, so a worker left over
  # from a previous split writes into this run's out_k; with an unchanged pair set the merge had no
  # way to tell it apart from its own workers.
  if (is.null(mani$run_id))
    stop("23 merge: the split manifest predates run identity. Re-run 21_coloc_split.R and the ",
         "workers.")
  # WHICH chunk this directory processed. The stamp recorded chunk_file and chunk_md5 and nothing
  # compared them, so a complete permutation -- chunk_2 processed in out_1 and so on -- passed every
  # check while each directory's claim about its own work was false.
  .ckf_k <- file.path(PC, sprintf("chunk_%d.rds", k))
  if (!identical(basename(.st$chunk_file %||% ""), basename(.ckf_k)))
    stop(sprintf("23 merge: out_%d holds work for %s, not %s.", k,
                 basename(.st$chunk_file %||% "<none>"), basename(.ckf_k)))
  if (!identical(tolower(.st$chunk_md5 %||% ""),
                 tolower(digest::digest(readBin(.ckf_k, "raw", file.size(.ckf_k)),
                                        algo = "md5", serialize = FALSE))))
    stop(sprintf("23 merge: out_%d was computed from a different version of %s than the one on ", k,
                 basename(.ckf_k)), "disk now.")
  if (!identical(.st$run_id %||% NA_character_, mani$run_id))
    stop(sprintf("23 merge: chunk %d belongs to split %s, but this merge is for split %s. ", k,
                 .st$run_id %||% "<none>", mani$run_id),
         "A worker from an earlier run wrote into this run's directory. Re-run that chunk.")
  if (is.null(.st$inputs_key) || !identical(.st$inputs_key, mani$inputs_key))
    stop(sprintf("23 merge: chunk %d was computed from different inputs than the split pinned ", k),
         "(chunk ", substr(.st$inputs_key %||% "<none>", 1, 12), ", split ",
         substr(mani$inputs_key, 1, 12), "):\n",
         paste(coloc_inputs_diff(mani$inputs_detail, .st$inputs_detail %||% list()),
               collapse = "\n"),
         "\nRe-run that chunk.")

  # pair-status is ALWAYS written by coloc main() -> require it as proof the worker ran to completion.
  psf <- file.path(PC, sprintf("out_%d", k), "coloc_pair_status_v7.csv")
  if (!file.exists(psf))
    stop(sprintf("23 merge: chunk %d pair-status MISSING (%s) -- worker did not finish; refusing a PARTIAL merge. Re-run 22 for chunk %d.", k, psf, k))
  # The stamp names the bytes this worker wrote. Without this the merge read the checkpoint for its
  # provenance key and then read the CSV and RDS separately, so a correct-provenance checkpoint could
  # bless files left behind by an earlier attempt -- the worker does not clear its old outputs, which
  # makes that the ordinary case for a re-run chunk.
  # ONE read. Hashing the path and then calling fread() on it is two reads, so a file replaced in
  # between is verified in one state and merged in another -- the whole point of the stamp check.
  .csv_raw <- readBin(psf, "raw", file.size(psf))
  .csv_now <- digest::digest(.csv_raw, algo = "md5", serialize = FALSE)
  if (!identical(tolower(.csv_now), tolower(.st$status_md5 %||% "")))
    stop(sprintf("23 merge: chunk %d's pair-status is not the file its worker wrote (stamp %s, ", k,
                 substr(.st$status_md5 %||% "<none>", 1, 12)),
         "on disk ", substr(.csv_now, 1, 12), "). It is left over from another attempt, or it was ",
         "modified after the run. Re-run that chunk.")
  ps[[k]] <- fread(text = rawToChar(.csv_raw))   # the bytes just hashed, not a fresh read
  # The file is written by the engine and has exactly two columns, pair and status. An earlier
  # version of this block asserted a vocabulary invented rather than read -- it listed "skipped" and
  # "not_evaluable", which the engine never emits, and rejected "skip:no_instrument", which it emits
  # routinely (R/coloc.R:1064). The vocabulary below is the engine's: ok, skip, skip:<reason>, error.
  if (!identical(sort(names(ps[[k]])), c("pair", "status")))
    stop(sprintf("23 merge: chunk %d pair-status has columns [%s]; the engine writes exactly ", k,
                 paste(names(ps[[k]]), collapse = ", ")),
         "pair,status. This file was not written by the reviewed engine.")
  .bad_st <- ps[[k]][is.na(status) |
                     !(status %in% c("ok", "error") | grepl("^skip(:|$)", status)), unique(status)]
  if (length(.bad_st))
    stop(sprintf("23 merge: chunk %d pair-status contains unrecognised value(s): ", k),
         paste(sQuote(.bad_st), collapse = ", "),
         ". Only 'ok', 'error' or a 'skip'/'skip:<reason>' state may appear; anything else needs no ",
         "result row and would drop the pair silently.")
  # coloc_results_v7_stage1.rds is only written when a chunk has >=1 successful coloc row; a legit all-skip chunk
  # has NO v5 (that is fine -- coverage is proven by pair-status, not by the result file).
  vf <- file.path(PC, sprintf("out_%d", k), "coloc_results_v7_stage1.rds")
  # Same binding as the status file: the stamp says whether this worker wrote a results file at all,
  # and if so exactly which bytes. A stale RDS from an earlier attempt now cannot be picked up, and
  # neither can a missing one be passed off as an all-skip chunk.
  # Same for the results: read once, hash those bytes, unserialise those bytes.
  .rds_raw <- if (file.exists(vf)) readBin(vf, "raw", file.size(vf)) else NULL
  .rds_now <- if (is.null(.rds_raw)) NA_character_
              else digest::digest(.rds_raw, algo = "md5", serialize = FALSE)
  if (!identical(tolower(.rds_now %||% "NA"), tolower(.st$results_md5 %||% "NA")))
    stop(sprintf("23 merge: chunk %d's results file is not what its worker wrote (stamp %s, ", k,
                 if (is.na(.st$results_md5 %||% NA)) "<none written>"
                 else substr(.st$results_md5, 1, 12)),
         "on disk ", if (is.na(.rds_now)) "<absent>" else substr(.rds_now, 1, 12),
         "). Re-run that chunk.")
  if (file.exists(vf)) {
    x <- unserialize(memDecompress(.rds_raw, type = "gzip")); setDT(x)
    # EXACT SET EQUALITY between the result rows and the pairs this chunk's own status marks ok.
    # Counting statuses proved only that every pair was ATTEMPTED. Without this, a result file from
    # another run, a missing row for an ok pair, or a row for a pair marked skip all merged cleanly,
    # because nothing ever compared the two files that are supposed to describe the same work.
    # The status side is keyed by its own `pair` column. Building the key from gene and outcome_id
    # here was wrong: those columns exist in the RESULTS table, not in pair-status, which holds only
    # pair and status -- so this check aborted with "object 'gene' not found" on every chunk that
    # produced any result, which is every real chunk.
    .okk <- ps[[k]][status == "ok", pair]
    .gotk <- x[, paste(gene, outcome_id, sep = "_")]
    if (anyDuplicated(.gotk))
      stop(sprintf("23 merge: chunk %d has duplicate (gene, outcome) rows in its results", k))
    if (!setequal(.okk, .gotk))
      stop(sprintf("23 merge: chunk %d's results do not match its own pair-status. ", k),
           length(setdiff(.okk, .gotk)), " pair(s) marked ok have no result row (",
           paste(head(setdiff(.okk, .gotk), 3), collapse = ", "), "); ",
           length(setdiff(.gotk, .okk)), " result row(s) are not marked ok (",
           paste(head(setdiff(.gotk, .okk), 3), collapse = ", "),
           "). The two files describe different work; re-run that chunk.")
    wk[[k]] <- x; cat(sprintf("  chunk %d: %d coloc rows (== %d ok pairs)\n", k, nrow(x), length(.okk)))
  }
  else {
    # 2026-07-01 (Codex redesign #1): a chunk with NO result file must be a genuine ALL-SKIP chunk. coloc
    # writes pair-status BEFORE saving results, so if the status has any "ok" but the RDS is absent, the worker
    # crashed after status / before results -> those ok RESULTS are lost while coverage still passes. Abort.
    .nok <- sum(ps[[k]]$status == "ok")
    if (.nok > 0) stop(sprintf("23 merge: chunk %d has %d ok pair(s) in pair-status but NO coloc_results_v7_stage1.rds -- worker crashed after writing status, before saving results (results lost). Re-run 22 for chunk %d.", k, .nok, k))
    cat(sprintf("  chunk %d: 0 coloc rows (genuine all-skip; pair-status has no ok)\n", k))
  }
}

# An all-skip run is legitimate, and it used to produce a zero-COLUMN data.table: the merge wrote
# that, then evaluated v5[is_coloc_pass(PP.H4, PP.H3)] on columns that do not exist and died AFTER
# publishing. The empty case now carries the same schema as the non-empty one, so the summary below
# evaluates and the saved RDS is readable by everything downstream.
# Columns and types taken from the engine's own output, not written from memory.
.V5_EMPTY <- data.table::data.table(
  gene = character(), outcome_id = character(), n_snps_coloc = integer(),
  PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(), PP.H4 = numeric(),
  coloc_support = logical(), independent_signal = logical(), outcome_name = character())
v5 <- if (length(wk)) rbindlist(wk, fill = TRUE) else copy(.V5_EMPTY)
# THIS IS NO LONGER A ROUNDING QUESTION, SO IT IS NO LONGER A WARNING. The guard was written when the
# engine stored PP.H3/PP.H4 rounded to four decimals while deciding coloc_support from the unrounded
# posteriors: a raw 0.79996 stored as 0.8000 made the flag and the recomputed rule disagree on a
# perfectly legitimate row, so aborting would have failed the parallel path where the serial path
# succeeds. Both producers now store the posteriors unrounded, and is_coloc_pass() is the same
# deterministic function on the same numbers -- so the two answers cannot differ on a well-formed
# row. A disagreement now means the row was written by different code or is otherwise malformed,
# which is exactly what this merge must not publish.
if (nrow(v5) && all(c("coloc_support", "PP.H4", "PP.H3") %in% names(v5))) {
  .recalc <- is_coloc_pass(v5$PP.H4, v5$PP.H3)
  .disag  <- which(!is.na(v5$coloc_support) & v5$coloc_support != .recalc)
  if (length(.disag)) {
    .df <- v5[.disag, .(gene, outcome_id, PP.H4, PP.H3, engine_call = coloc_support,
                        recomputed_call = .recalc[.disag])]
    fwrite(.df, file.path(REPRO, "coloc_support_disagreement_pairs.csv"))
    stop("23 merge: ", length(.disag), " chunk row(s) whose stored coloc_support contradicts the dual ",
         "rule recomputed from the stored PP.H4/PP.H3. The posteriors are stored unrounded, so these ",
         "cannot disagree unless the row came from a different engine version or is corrupt: ",
         paste(utils::head(.df[, paste0(gene, "~", outcome_id, " (PP.H4=", PP.H4, ", PP.H3=", PP.H3,
                                        ", stored=", engine_call, ", recomputed=", recomputed_call, ")")], 3),
               collapse = "; "),
         ". Listed in coloc_support_disagreement_pairs.csv -- re-run the affected chunks.")
  }
}
# A non-empty merge must carry that schema too: a chunk written by an older engine would otherwise
# reach the summary and the downstream stages with columns missing.
if (nrow(v5) && length(setdiff(names(.V5_EMPTY), names(v5))))
  stop("23 merge: the merged results are missing column(s) the engine is expected to produce: ",
       paste(setdiff(names(.V5_EMPTY), names(v5)), collapse = ", "),
       ". A chunk was written by a different version of R/coloc.R.")
if (nrow(v5) && (!"outcome_name" %in% names(v5) || any(is.na(v5$outcome_name))))
  v5[, outcome_name := names(OUTCOMES)[match(outcome_id, OUTCOMES)]]
n0 <- nrow(v5)
# Fix (audit Blocker + Low): BIND the merge to the mr_suggestive the split saw -- run UNCONDITIONALLY (even if
# v5 is empty / all-skip), so a mr_suggestive rewrite between 21_split and here is caught regardless of whether
# any pair produced a coloc row. Mirror the serial main()/21_split universe (IVW/Wald, non-eGFR) first.
.sug  <- readRDS(file.path(REPRO, "mr_suggestive.rds")); setDT(.sug)
# STRATIFIED BY outcome_id ALONE, matching 04 and 06. Two things were wrong with the name clause:
# `!=` is not NA-safe, so a row with a missing outcome_name was dropped here while 04 kept it (the
# two universes then disagreed and this check failed the merge); and excluding on the NAME while 05
# selects the eGFR stratum on the ID meant a row whose id is non-eGFR but whose name reads "eGFR"
# belonged to neither stratum. The id is canonical; the name is derived from it.
.sugf <- .sug[method %in% c("Inverse variance weighted", "Wald ratio") &
              !(outcome_id %in% mani$egfr_id)]
.curkeys <- unique(paste(.sugf$gene, .sugf$outcome_id, sep = "_"))
if (!setequal(.curkeys, mani$expected_keys))
  stop(sprintf("23 merge: mr_suggestive CHANGED since 21_coloc_split (now %d non-eGFR IVW/Wald pairs vs manifest %d; +%d/-%d) -- re-run 21_split + 22 workers before merging",
               length(.curkeys), length(mani$expected_keys),
               length(setdiff(.curkeys, mani$expected_keys)), length(setdiff(mani$expected_keys, .curkeys))))
if (nrow(v5)) {
  # chunks are DISJOINT -> a duplicate (gene,outcome_id) row means overlapping/stale chunk output -> ABORT loudly.
  .dk <- paste(v5$gene, v5$outcome_id)
  if (anyDuplicated(.dk))
    stop(sprintf("23 merge: %d DUPLICATE (gene,outcome_id) coloc row(s) across chunks (must be disjoint -- overlapping/stale chunk): %s",
                 sum(duplicated(.dk)), paste(head(.dk[duplicated(.dk)], 5), collapse = ", ")))
  # re-order to MATCH the single-process (serial) processing order (mr_suggestive first-appearance) so the
  # DEFAULT parallel v5 is byte-identical to COLOC_SERIAL=1 (coloc is per-locus independent; only row ORDER differed).
  .serord <- unique(paste(.sugf$gene, .sugf$outcome_id))
  v5[, .ord := match(paste(gene, outcome_id), .serord)]
  if (anyNA(v5$.ord)) stop(sprintf("23 merge: %d coloc row(s) not found in the serial pair order -- v5/mr_suggestive universe mismatch", sum(is.na(v5$.ord))))
  v5 <- v5[order(.ord)][, .ord := NULL]
}

# Unified pair-status (parallel analogue of serial coloc's inline write) -> feeds ⟦99⟧ A2 + error-abort.
psall <- rbindlist(ps, fill = TRUE)
# Fix (audit Low): chunks are disjoint -> REJECT duplicate status keys instead of unique()-collapsing them, so a
# stale/overlapping chunk's conflicting status cannot be silently hidden (matches the eGFR merge's dup guard).
if (anyDuplicated(psall$pair))
  stop(sprintf("23 merge: %d DUPLICATE status key(s) across chunks (must be disjoint -- overlapping/stale chunk): %s",
               sum(duplicated(psall$pair)), paste(head(psall$pair[duplicated(psall$pair)], 5), collapse = ", ")))
.errp <- psall[status == "error"]
if (nrow(.errp) > 0)
  stop(sprintf("23 merge: %d coloc pair(s) ERRORED across chunks (not a legit skip): %s -- do not ship silent gaps", nrow(.errp), paste(head(.errp$pair, 10), collapse = ", ")))
# EXACT coverage vs the manifest: every expected pair has a status, and no unexpected pair sneaked in.
.miss  <- setdiff(mani$expected_keys, psall$pair)
.extra <- setdiff(psall$pair, mani$expected_keys)
if (length(.miss) > 0)
  stop(sprintf("23 merge: %d expected pair(s) have NO coloc status (silently dropped by a partial chunk): %s", length(.miss), paste(head(.miss, 10), collapse = ", ")))
if (length(.extra) > 0)
  stop(sprintf("23 merge: %d UNEXPECTED pair(s) in coloc status (not in the split set -- stale/mismatched chunk): %s", length(.extra), paste(head(.extra, 10), collapse = ", ")))

# Fix (audit Low): write pair-status in the SAME serial (mr_suggestive first-appearance) order as COLOC_SERIAL=1
# (v5 results are already reordered; the status CSV must match too for byte-equivalence). .curkeys = serial order.
psall <- psall[order(match(pair, .curkeys))]

# EVERYTHING that can fail is evaluated BEFORE anything is published, and the two canonical files
# are then written via temporary files and renamed. Previously the status CSV was overwritten first
# and the results RDS second, both before the summary expressions ran, so a failure in between --
# the all-skip schema defect was exactly such a failure -- left the canonical path holding one new
# file and one old one, with nothing recording that the pair was inconsistent.
.v5s <- v5[is_coloc_pass(PP.H4, PP.H3)]   # dual coloc rule (D1 2026-07-05 R2)
.by_outcome <- if (nrow(.v5s)) .v5s[, .N, by = outcome_name][order(-N)] else NULL
.n_ok <- sum(psall$status == "ok"); .n_skip <- sum(grepl("skip", psall$status))

.tmp_csv <- file.path(REPRO, "coloc_pair_status_v7.csv.tmp")
.tmp_rds <- file.path(REPRO, "coloc_results_v7_stage1.rds.tmp")
.dst_csv <- file.path(REPRO, "coloc_pair_status_v7.csv")
.dst_rds <- file.path(REPRO, "coloc_results_v7_stage1.rds")
fwrite(psall, .tmp_csv)
# THE KEY TRAVELS WITH THE RESULT. 06 could only establish currency by mtime, so a stale result file
# touched after mr_suggestive.rds passed every gate even though its numbers came from other GWAS,
# instrument, chain or coordinate bytes. The split manifest already pins those inputs, but it lives
# in par_chunks/ and is not bound to whatever file is later found at the canonical root -- carrying
# the key ON the object is what lets 06 recompute it and compare.
attr(v5, "coloc_inputs_key") <- mani$inputs_key
saveRDS(v5, .tmp_rds)
# The RESULTS go first and the status second, and a failed second rename PUTS THE FIRST BACK. Two
# renames cannot be one atomic act, but they can be all-or-nothing across the pair: publishing the
# status first and failing left a new status table beside an old result table, while the message
# said the previous results were untouched -- which had stopped being true one line earlier.
.prev_rds <- if (file.exists(.dst_rds))
  file.path(REPRO, paste0("coloc_results_v7_stage1.rds.rollback_", Sys.getpid())) else NA_character_
# THE ROLLBACK IS CHECKED, and the message says what actually happened. Every file.rename() here
# used to be called for effect and its return value discarded, while the error text asserted that the
# previous outputs were unchanged or had been restored -- so a rollback that itself failed reported
# success at restoring. A message describing a state the code did not verify is the failure this
# whole layer exists to prevent, one level up.
.moved <- if (!is.na(.prev_rds)) file.rename(.dst_rds, .prev_rds) else TRUE
if (!isTRUE(.moved))
  stop("23 merge: could not move the previous results aside (", .dst_rds, "), so a failed install ",
       "could not be rolled back. Nothing was installed.")
if (!file.rename(.tmp_rds, .dst_rds)) {
  .back <- if (!is.na(.prev_rds)) file.rename(.prev_rds, .dst_rds) else TRUE
  unlink(c(.tmp_csv, .tmp_rds))
  stop("23 merge: could not install the merged results. ",
       if (isTRUE(.back)) "The previous outputs are unchanged."
       else paste0("AND the previous results could not be put back -- they are at ", .prev_rds,
                   "; move them into place by hand before rerunning."))
}
if (!file.rename(.tmp_csv, .dst_csv)) {
  # Put the results back, so the two canonical files stay a matched pair.
  unlink(.dst_rds)
  .back <- if (!is.na(.prev_rds)) file.rename(.prev_rds, .dst_rds) else TRUE
  unlink(c(.tmp_csv, .tmp_rds))
  stop("23 merge: the results were installed but the pair-status could not be. ",
       if (isTRUE(.back)) "The results were rolled back; the canonical outputs are as they were before this merge."
       else paste0("The rollback ALSO failed: there is now no canonical result file, and the previous one is at ",
                   .prev_rds, ". Move it into place by hand before rerunning."))
}
if (!is.na(.prev_rds)) unlink(.prev_rds)

# Counts are derived from the manifest and the assembled table, never written as literals: a literal
# here names neither the pair count nor the row count, so it states a number belonging to nothing on screen.
cat(sprintf("\n=== non-eGFR coloc v5 (parallel) reassembled: %d row(s) from %d expected pair(s) (dedup -%d) ===\n",
            nrow(v5), length(mani$expected_keys), n0 - nrow(v5)))
cat(sprintf("pair-status covers all %d expected pairs (%d ok / %d skip)\n",
            length(mani$expected_keys), .n_ok, .n_skip))
cat(sprintf("strong (dual rule PP.H4>=0.80 & H4/H3>=5): %d pairs / %d proteins\n",
            nrow(.v5s), uniqueN(.v5s$gene)))
if (!is.null(.by_outcome)) print(.by_outcome)
cat("\nWrote coloc_results_v7_stage1.rds + coloc_pair_status_v7.csv. Next: 05_coloc_egfr + 06_combine.\n")
