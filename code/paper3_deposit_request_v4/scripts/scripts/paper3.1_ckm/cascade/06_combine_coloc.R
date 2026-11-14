# A1 combine: merge corrected non-eGFR coloc (v5) + eGFR coloc -> corrected v6.
# Run after A1 (05) + A1b (06) both done. Prints corrected coloc-supported summary
# vs the pre-correction manuscript. Writes results/repro_corrected/coloc_results_v7.rds.
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
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
REPRO <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))

v5_f  <- file.path(REPRO, "coloc_results_v7_stage1.rds")
egf_f <- file.path(REPRO, "coloc_egfr_results.rds")
stopifnot(file.exists(v5_f), file.exists(egf_f))

# THE FOUR CANONICAL MR FILES MUST COME FROM ONE RUN. 02_promote_mr.R installs them with four
# sequential copies; an interruption leaves some new and some from the previous generation, and
# nothing downstream noticed -- 04/05 read only mr_suggestive.rds, and the parallel fingerprint hashes
# that file's coloc-visible fields, which establishes consistency with whatever file is present rather
# than with the MR generation it came from. ckm_run_id, not ckm_entry_script_md5: two runs of the same
# script share a script hash. 02 has already shown that one run id implies matching Phase-1 and
# outcome stamps, so this one comparison covers the mixture.
.mr_names <- c("mr_all_results.rds", "mr_suggestive.rds",
               "mr_significant_bonferroni.rds", "mr_significant_fdr.rds")
.mr_run <- vapply(file.path(REPRO, .mr_names), function(f) {
  if (!file.exists(f)) return(NA_character_)
  z <- attr(readRDS(f), "ckm_run_id", exact = TRUE)
  if (is.null(z)) NA_character_ else as.character(z)[1]
}, character(1))
if (anyNA(.mr_run) || length(unique(.mr_run)) != 1L)
  stop("08 combine: the four canonical MR artifacts do not carry one ckm_run_id (",
       paste(sprintf("%s=%s", .mr_names, substr(ifelse(is.na(.mr_run), "MISSING", .mr_run), 1, 24)),
             collapse = "; "),
       ") -- 02_promote_mr.R was interrupted part way and the canonical set is mixed. Rerun it.")

# 2026-07-01 Codex-phase3 #4: refuse to combine STALE coloc. (a) 04/05 must have written their per-pair
# status files (proof of a COMPLETED clean run, not a leftover from an aborted one); (b) both coloc outputs
# must be NEWER than the current mr_suggestive.rds (else they were built for a previous suggestive set).
.ps5 <- file.path(REPRO, "coloc_pair_status_v7.csv"); .pse <- file.path(REPRO, "coloc_egfr_pair_status.csv")
if (!file.exists(.ps5) || !file.exists(.pse)) stop("08 combine: coloc pair-status file(s) missing -- 04/05 did not complete a clean run; refusing to combine possibly-stale v5/eGFR")
.sug_f <- file.path(REPRO, "mr_suggestive.rds")
if (!file.exists(.sug_f)) stop("08 combine: mr_suggestive.rds MISSING -- cannot verify coloc currency; run 01_mr + 02_promote first")   # 2026-07-01 reverify: missing suggestive must not silently skip the currency check
if (file.mtime(v5_f) <= file.mtime(.sug_f) || file.mtime(egf_f) <= file.mtime(.sug_f))   # <= : coloc must be STRICTLY newer (equal mtime = suspect stale)
  stop("08 combine: coloc v5/eGFR is NOT newer than the current mr_suggestive.rds -- stale coloc from a previous suggestive set; re-run 04/05 before combining")

v5  <- readRDS(v5_f);  setDT(v5)
egf <- readRDS(egf_f); setDT(egf)
# Fix (audit High): the v5 RESULT must cover EVERY status=ok pair in coloc_pair_status_v7.csv. 23_merge asserts
# this for the parallel path, but 06 must independently reject a stale/partial v5 it did not produce (the
# existence + mtime gates above would still pass a v5 that silently dropped ok pairs).
.st5 <- fread(.ps5); .ste <- fread(.pse)
# THE PAIR UNIVERSE, REBUILT FROM THE SAME INPUT 04 AND 05 READ. Every other check here compares the
# results against the status files, so a pair absent from BOTH is invisible and the published table is
# quietly short. Two ordinary paths produce that input: 04's eGFR filter used to drop a row whose
# outcome_name was missing, and each chunked eGFR worker writes its partial status to the canonical
# filename, so an aborted merge can leave one chunk's status looking like a complete run.
.sugd <- as.data.table(readRDS(.sug_f))
.prim <- .sugd[method %in% c("Inverse variance weighted", "Wald ratio")]   # the methods coloc iterates
.EGID <- unname(OUTCOMES[["eGFR"]])
# STRATIFIED BY outcome_id ALONE, matching 04 and 23. The non-eGFR set used to also exclude on
# outcome_name, while the eGFR set selected on id only -- so a row whose id is non-eGFR but whose
# name reads "eGFR" fell out of both, and because this reconstruction repeated the same rule, the
# exact-set checks below still passed on a universe that was already short one pair.
.exp5 <- unique(.prim[!(outcome_id %in% .EGID), paste(gene, outcome_id, sep = "_")])
.expe <- unique(.prim[outcome_id %in% .EGID, paste(gene, outcome_id, sep = "_")])
# EXACT SETS, not one-way coverage, and a closed status vocabulary. A one-way check passed a status
# table that ALSO held pairs the current suggestion set does not contain, and a duplicated or
# unrecognised status was invisible to every gate below.
if (anyDuplicated(.st5$pair) || anyDuplicated(.ste$pair))
  stop(sprintf("08 combine: duplicate pair-status rows (%d non-eGFR, %d eGFR) -- a partial chunk status was appended rather than replaced", sum(duplicated(.st5$pair)), sum(duplicated(.ste$pair))))
if (!setequal(.st5$pair, .exp5))
  stop(sprintf("08 combine: the non-eGFR pair-status set differs from the current suggestion universe (%d missing, %d unexpected; e.g. missing %s / unexpected %s)", length(setdiff(.exp5, .st5$pair)), length(setdiff(.st5$pair, .exp5)), paste(head(setdiff(.exp5, .st5$pair), 3), collapse = ", "), paste(head(setdiff(.st5$pair, .exp5), 3), collapse = ", ")))
if (!setequal(.ste$pair, .expe))
  stop(sprintf("08 combine: the eGFR pair-status set differs from the current suggestion universe (%d missing, %d unexpected; e.g. missing %s / unexpected %s)", length(setdiff(.expe, .ste$pair)), length(setdiff(.ste$pair, .expe)), paste(head(setdiff(.expe, .ste$pair), 3), collapse = ", "), paste(head(setdiff(.ste$pair, .expe), 3), collapse = ", ")))
.STATUS_OK <- c("ok", "skip", "skip:no_instrument", "error")
if (anyNA(.st5$status) || !all(.st5$status %chin% .STATUS_OK) ||
    anyNA(.ste$status) || !all(.ste$status %chin% .STATUS_OK))
  stop(sprintf("08 combine: missing or unrecognised pair status (%s) -- the status vocabulary is ok/skip/skip:no_instrument/error", paste(head(setdiff(unique(c(.st5$status, .ste$status)), .STATUS_OK), 5), collapse = ", ")))

# THE RESULTS MUST NAME THE INPUTS THEY WERE COMPUTED FROM. Everything above this line establishes
# that the pair BOOKKEEPING is complete and internally consistent; none of it can tell that the
# numbers came from the GWAS, instrument, coordinate and chain bytes this session is using. A stale
# result file whose pair set happens to match, touched after mr_suggestive.rds, satisfied every gate.
# 23_merge and 05 now stamp their output with the key; recomputing it here is what binds them.
source(file.path(BASE_DIR, "R/_coloc_inputs_key.R"))
.p5 <- .prim[!(outcome_id %in% .EGID)]
.pe <- .prim[outcome_id %in% .EGID]
.k5 <- coloc_inputs_key(list(genes = sort(unique(.p5$gene)),
                             outcome_ids = sort(unique(.p5$outcome_id))))$key
.ke <- paste(unname(tools::md5sum(file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/05_coloc_egfr.R"))),
             coloc_inputs_key(list(genes = sort(unique(.pe$gene)), outcome_ids = .EGID))$key, sep = "|")
.gk <- function(x) { z <- attr(x, "coloc_inputs_key", exact = TRUE); if (is.null(z)) "<none>" else as.character(z)[1] }
if (!identical(.gk(v5), .k5) || !identical(.gk(egf), .ke))
  stop(sprintf(paste0("08 combine: the coloc results do not name the inputs this session is using ",
                      "(non-eGFR result key %s vs current %s; eGFR result key %s vs current %s). They ",
                      "were computed from different GWAS/instrument/coordinate/chain bytes, or predate ",
                      "the stamping. Re-run 04 and 05."),
               substr(.gk(v5), 1, 12), substr(.k5, 1, 12), substr(.gk(egf), 1, 12), substr(.ke, 1, 12)))
# Fix (audit Med): a status CSV containing an ERROR row means an aborted coloc run -> refuse to combine (do not
# trust that a matching-ok result RDS is a clean, complete run).
if (nrow(.st5[status == "error"]) > 0) stop(sprintf("08 combine: %d non-eGFR pair(s) status=error in coloc_pair_status_v7.csv (aborted run): %s", nrow(.st5[status == "error"]), paste(head(.st5[status == "error", pair], 10), collapse = ", ")))
if (nrow(.ste[status == "error"]) > 0) stop(sprintf("08 combine: %d eGFR pair(s) status=error in coloc_egfr_pair_status.csv (aborted run): %s", nrow(.ste[status == "error"]), paste(head(.ste[status == "error", pair], 10), collapse = ", ")))
# v5 TWO-WAY bijection: result rows == status=ok set (no missing ok, no stale extra, no duplicate).
.ok5v  <- .st5[status == "ok", pair]
.res5v <- if (nrow(v5)) paste(v5$gene, v5$outcome_id, sep = "_") else character(0)
if (anyDuplicated(.res5v)) stop(sprintf("08 combine: %d DUPLICATE (gene,outcome_id) row(s) in coloc_results_v7_stage1.rds", sum(duplicated(.res5v))))
if (length(setdiff(.ok5v, .res5v)) > 0) stop(sprintf("08 combine: %d status=ok pair(s) ABSENT from coloc_results_v7_stage1.rds (stale/partial v5): %s", length(setdiff(.ok5v, .res5v)), paste(head(setdiff(.ok5v, .res5v), 10), collapse = ", ")))
if (length(setdiff(.res5v, .ok5v)) > 0) stop(sprintf("08 combine: %d row(s) in coloc_results_v7_stage1.rds are NOT status=ok (stale extra): %s", length(setdiff(.res5v, .ok5v)), paste(head(setdiff(.res5v, .ok5v), 10), collapse = ", ")))
# eGFR TWO-WAY bijection -- use egf's OWN outcome_id (Fix audit High: forced .EGID would let a wrong-id row pass).
.oke   <- .ste[status == "ok", pair]
# NOT BACK-FILLED. A missing outcome_id used to be replaced with the current eGFR id, which relabels
# a result of unknown provenance as this run's -- the opposite of what the bijection below is for.
if (!"outcome_id" %in% names(egf))
  stop("08 combine: coloc_egfr_results.rds has no outcome_id column, so the outcome its rows describe cannot be established; it will not be labelled as the current one")
.egcol <- egf$outcome_id
.rese  <- if (nrow(egf)) paste(egf$gene, .egcol, sep = "_") else character(0)
if (anyDuplicated(.rese)) stop(sprintf("08 combine: %d DUPLICATE row(s) in coloc_egfr_results.rds", sum(duplicated(.rese))))
if (length(setdiff(.oke, .rese)) > 0) stop(sprintf("08 combine: %d eGFR status=ok pair(s) ABSENT from coloc_egfr_results.rds (stale/partial): %s", length(setdiff(.oke, .rese)), paste(head(setdiff(.oke, .rese), 10), collapse = ", ")))
if (length(setdiff(.rese, .oke)) > 0) stop(sprintf("08 combine: %d row(s) in coloc_egfr_results.rds are NOT status=ok (stale/malformed eGFR rows): %s", length(setdiff(.rese, .oke)), paste(head(setdiff(.rese, .oke), 10), collapse = ", ")))
if (nrow(egf) == 0) stop("08 combine: eGFR coloc results are EMPTY -- 05 produced no eGFR pairs; refusing to save v6 without eGFR")
if (!"outcome_name" %in% names(v5))  v5[,  outcome_name := names(OUTCOMES)[match(outcome_id, OUTCOMES)]]
if (!"outcome_name" %in% names(egf)) egf[, outcome_name := "eGFR"]
# (the outcome_id backfill that stood here is gone: its absence is now fatal above, because relabelling
#  a result whose outcome cannot be established as the current one defeats the bijection check)

v6_raw <- rbind(v5, egf, fill = TRUE)
# Fix (audit Blocker): v5 (non-eGFR) and egf (eGFR) must be DISJOINT per (gene,outcome_id); a duplicate here is
# contamination (a stale eGFR row leaked into v5, or duplicate egf rows) and the max-PP.H4 collapse below would
# silently keep one -> ABORT loudly instead.
.dk6 <- paste(v6_raw$gene, v6_raw$outcome_id)
if (anyDuplicated(.dk6)) stop(sprintf("08 combine: %d DUPLICATE (gene,outcome_id) row(s) in union(v5,egf) -- v5/egf must be disjoint (double-coloc/stale contamination): %s", sum(duplicated(.dk6)), paste(head(.dk6[duplicated(.dk6)], 5), collapse = ", ")))
# 2026-06-30 (gate obj #6): assert a CLEAN union — dedup may only collapse identical gene+outcome rows,
# never lose a distinct pair or carry an unexpected extra (double-coloc / stale eGFR / duplicate chunks).
.exp_pairs <- unique(c(paste(v5$gene, v5$outcome_id), paste(egf$gene, egf$outcome_id)))
# ...AND NOTHING IS LEFT TO COLLAPSE. The guard above aborts on any duplicate (gene,outcome_id), so
# the max-PP.H4 dedup that stood here, and the line reporting how many rows it had collapsed, could
# never run. Code that cannot execute still reads as a recovery the script performs; the sort key is
# unchanged in effect because ties in (gene,outcome_id) no longer reach this line.
v6 <- v6_raw[order(gene, outcome_id)]
.got_pairs <- paste(v6$gene, v6$outcome_id)
if (!setequal(.got_pairs, .exp_pairs)) stop(sprintf("08 combine: v6 pair-set != union(v5,egf) — lost %d / extra %d pair(s)",
  length(setdiff(.exp_pairs, .got_pairs)), length(setdiff(.got_pairs, .exp_pairs))))
# Required-column guard BEFORE persisting (2026-07-05 R2-High: don't write a malformed v6 then fail after).
stopifnot(all(c("gene","outcome_name","PP.H4","PP.H3","n_snps_coloc") %in% names(v6)))   # R3: guard full CSV schema, not just PP columns
# THE UNION CARRIES BOTH KEYS. v5 and egf each arrive stamped and are each verified above, but rbind
# drops attributes, so the published table itself named nothing -- and 26_susie reads THIS file, not
# its two inputs. Retaining one input's attribute would certify half the table, so both are recorded,
# under names that say which stratum each belongs to.
attr(v6, "coloc_inputs_key_noneGFR") <- .k5
attr(v6, "coloc_inputs_key_eGFR")    <- .ke
attr(v6, "ckm_combined_at")          <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
saveRDS(v6, file.path(REPRO, "coloc_results_v7.rds"))

# Dual coloc rule (PP.H4>=0.80 AND PP.H4/PP.H3>=5); 2026-07-04 D1 (was H4-only -> admitted MFAP4-CKD at ratio 4.80).
strong <- v6[is_coloc_pass(PP.H4, PP.H3)]
# "pairs with a posterior", not "pairs tested": a skipped pair is attempted but has no row here, so
# nrow(v6) is the successful subset. The status files above hold the attempted count (ok + skip).
cat(sprintf("=== CORRECTED coloc v6: %d pairs with a posterior | strong (dual rule PP.H4>=0.80 & PP.H4/PP.H3>=5) = %d pairs, %d proteins ===\n",
            nrow(v6), nrow(strong), uniqueN(strong$gene)))
cat("--- corrected strong coloc by outcome ---\n")
print(strong[, .N, by = outcome_name][order(-N)])
fwrite(strong[order(outcome_name, -PP.H4), .(gene, outcome = outcome_name, outcome_id, PP.H4, PP.H3, n_snps_coloc)],   # P3 F-02: carry canonical outcome_id, not just short name
       file.path(REPRO, "corrected_strong_coloc_hits.csv"))
cat(sprintf("\nWrote coloc_results_v7.rds + corrected_strong_coloc_hits.csv\n"))
