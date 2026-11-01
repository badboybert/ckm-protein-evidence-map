# ============================================================================
# Stage 0b (dashed candidate layer): run coloc.susie on the candidate pool.
#
# Generalises phase3/27 to the full candidate pool from 28, REUSING the frozen
# engine per-pair functions without modifying them.
#
# MODE dispatch (env CAND_SUSIE_MODE = "egfr" | "nonegfr") — run BOTH, separately;
#   the two engines define same-named helpers with different signatures, so each mode
#   loads exactly ONE engine:
#     egfr    -> eval 26's head (precache-free; the phase3/27 precedent) -> run_susie_egfr_pair
#     nonegfr -> source config + coloc.R, then eval ONLY 15's function-def block
#                (get_ld_matrix_auto + run_susie_pair); we skip 15's head so its
#                strong-gene precache can't fail-loud on an unrelated strong pair.
#
# ERROR POLICY (reason-aware; Codex High #1): a NULL SuSiE result is NOT uniformly a benign skip.
#   The engines return NULL for BOTH genuinely benign per-locus cases (too few SNPs / no LD overlap / empty
#   coloc summary) AND genuine computation/data failures (coloc.susie threw / pQTL missing columns / no-data);
#   they THROW (fail_loud_data or any other error) for infrastructure failures. Each pair's engine output is
#   captured (run_one_capture) and classify_null() (_candidate_null_classify.R) decides:
#     NULL + known-benign reason -> record run_status="skip_benign" (error_msg=tag) and continue;
#     NULL + fatal/UNRECOGNISED  -> ABORT (default-fatal-on-unknown);  THROW -> ABORT.
#   ABORT saves the checkpoint WITHOUT the failing pair (retried on resume). We NEVER demote a computation
#   failure to a benign "uncomputable" edge. Shared infrastructure is prechecked once up front (fail loud early).
#
# Only 28-assigned testable candidates (susie_engine set) are run: MHC, undetermined-MHC, and no-coords
#   pairs carry susie_engine=NA and are skipped here; coloc-untested pairs were EXCLUDED by 28 (written to
#   candidate_pool_no_coloc.csv) and never enter this pool at all.
#
# NON-DESTRUCTIVE: writes ONLY candidate_susie_<mode>.csv (+ its own checkpoint).
#   Canonical susie_coloc_results.* are never read or written (we eval the engine head
#   strictly before its Main loop, the only place it writes).
#
# Run (both modes):
#   cd "<E>/20260331-CKM/Paper 3"
#   CAND_SUSIE_MODE=nonegfr Rscript scripts/paper3.1_ckm/phase3/29_candidate_susie.R
#   CAND_SUSIE_MODE=egfr    Rscript scripts/paper3.1_ckm/phase3/29_candidate_susie.R
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
suppressPackageStartupMessages({ library(data.table); library(coloc); library(ieugwasr); library(susieR) })
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")

MODE <- tolower(Sys.getenv("CAND_SUSIE_MODE", ""))
if (!MODE %in% c("egfr", "nonegfr")) stop("29: set CAND_SUSIE_MODE=egfr or nonegfr (run both, separately)")
# Load config HERE -- before the mode branch, before the LD/infra precheck that reads PLINK_BIN /
# LD_EUR_DIR / EGFR_FILE, and before assert_outcome_truth() on the pool. A first attempt (r24) put this
# after the branch and before the pool guard only; the constants at the precheck are used EARLIER, so the
# script still depended on the branch having sourced config as a side effect (r25 High). It worked, but
# "works because the engine happens to load my dependency" is what killed 31.
# NOTE the verification gap this exposed: the static check written for r23 asks "does this file source
# config?", not "does it source it BEFORE first use?" -- so it reported OK. parse() cannot see symbol
# resolution; that check could not see ORDER. Sourcing is idempotent, so the branch may still source it.
source(file.path(B, "R/00_config.R"))   # PLINK_BIN, LD_EUR_DIR, OUTCOMES, assert_outcome_truth()
# NOT EGFR_FILE: that is defined by evaluating 26's head in the eGFR branch below, not by config. A first
# version of this line claimed config provided it (r26 Low) -- and the static check written to verify this
# very ordering had EGFR_FILE in its config-symbol list too, so the check agreed with the false claim.
# The eGFR precheck's `exists("EGFR_FILE")` is therefore load-bearing, not defensive padding.
cat(sprintf("=== 29 candidate SuSiE | MODE=%s ===\n", MODE))
source(file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"))   # classify_null(): benign skip vs FATAL (default-fatal-on-unknown)

RC      <- file.path(B, file.path("results", paste0("repro_", PHASE1_TAG)))
POOL    <- file.path(RC, "phase3_susie/candidate_pool.csv")
OUT_DIR <- file.path(RC, "phase3_susie")
OUT_CSV <- file.path(OUT_DIR, sprintf("candidate_susie_%s.csv", MODE))
CKPT    <- file.path(OUT_DIR, sprintf("candidate_susie_%s_checkpoint.rds", MODE))
# This CSV is written once, hours later, at the very end. Left in place it survives a killed rerun as
# a COMPLETE, schema-correct file from the previous run, and 30 -- which checks shape, never age --
# audits the old posteriors with nothing in the output saying so. Removed here, before any validation,
# so no window remains in which a kill leaves the previous answer behind. Starting this script means
# its output is being replaced: an abort now costs a recompute, which is the cheaper error.
# unlink() reports failure by return value and R does not raise, so an unchecked call is a guard that
# can silently not be there -- on Windows a locked file is the realistic case.
unlink(OUT_CSV); if (file.exists(OUT_CSV)) stop(sprintf("29: could not remove the previous %s -- refusing to run with a stale result in place", OUT_CSV))
if (!file.exists(POOL)) stop(sprintf("29: candidate pool missing (%s) -- run 28 first", POOL))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TEMPLATE_COLS <- c("gene", "outcome_name", "outcome_id", "susie_engine", "PP.H4_abf", "PP.H4_susie", "PP.H3_susie",
                   "concordant_susie", "n_snps_merged", "n_snps_ld", "n_signals_exp", "n_signals_out",
                   "n_signal_pairs",
                   # PER-SIDE convergence, from the fit's own `$converged`. Until 2026-07-16 nothing in this
                   # study recorded convergence per fit except 31, and only per EDGE across both susie_rss
                   # calls -- so "did the primary fits converge?" was unanswerable for all 268 fits this
                   # script produces, and whether a non-convergence was the exposure or the outcome could
                   # not be told apart. Naming them here is what BACK-FILLS them on rows that lack them (the
                   # loop near the write); it is NOT what keeps them -- `setcolorder()` orders columns and
                   # does not subset, so an unnamed column would survive, just last. Said the other way
                   # round in an earlier comment, which was simply false.
                   "converged_exp", "converged_out",
                   "run_status", "error_msg")

# corrected inputs the engine head/functions expect + SuSiE out dir pointed at phase3_susie
SUSIE_COLOC_INPUT <- file.path(RC, "coloc_results_v7.rds")
SUSIE_INSTR       <- PHASE1_INSTRUMENTS
SUSIE_OUT_DIR     <- OUT_DIR
SUSIE_EGFR_COLOC  <- SUSIE_COLOC_INPUT
SUSIE_EGFR_INSTR  <- SUSIE_INSTR
SUSIE_EGFR_OUT    <- OUT_DIR

# ── Load exactly one engine ──────────────────────────────────────────────────
if (MODE == "egfr") {
  eng <- file.path(B, "scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R")
  src <- readLines(eng); cut <- grep("Main loop", src)[1]
  if (is.na(cut)) stop("29: no 'Main loop' marker in 26 -- engine layout changed")
  cat(sprintf("29: eval 26 head (lines 1-%d) -> run_susie_egfr_pair\n", cut - 1))
  eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")), envir = globalenv())
} else {
  source(file.path(B, "R/00_config.R"))
  COLOC_SOURCE_ONLY <- TRUE; suppressMessages(source(file.path(B, "R/coloc.R"))); rm(COLOC_SOURCE_ONLY)
  eng <- file.path(B, "scripts/paper3.1_ckm/phase3/15_susie_coloc.R")
  src <- readLines(eng)
  a <- grep("Helper: LD matrix", src)[1]; z <- grep("Main loop", src)[1] - 1L
  if (is.na(a) || is.na(z) || a > z) stop("29: cannot locate 15's function-def block -- engine layout changed")
  cat(sprintf("29: eval 15 function block (lines %d-%d) -> get_ld_matrix_auto + run_susie_pair (skips strong-precache)\n", a, z))
  eval(parse(text = paste(src[a:z], collapse = "\n")), envir = globalenv())
  rsid_map <- load_rsid_map()
}

# ── Shared-infra precheck (fail loud early; per-locus failures ABORT via the loop, never demote) ──
if (!(exists("PLINK_BIN") && !is.null(PLINK_BIN) && file.exists(PLINK_BIN)))
  stop("29: PLINK_BIN missing/not found -- LD infrastructure unavailable")
if (!exists("LD_EUR_DIR") || is.null(LD_EUR_DIR)) stop("29: LD_EUR_DIR undefined")
for (ext in c("bed", "bim", "fam"))
  if (!file.exists(file.path(LD_EUR_DIR, paste0("EUR.", ext))))
    stop(sprintf("29: 1KG EUR LD panel EUR.%s missing in %s -- infra failure", ext, LD_EUR_DIR))
if (MODE == "egfr" && !(exists("EGFR_FILE") && file.exists(EGFR_FILE)))
  stop("29: eGFR GWAS file missing -- infra failure for MODE=egfr")

# ── Candidates for this mode (MHC / undetermined-MHC / no-coords carry engine=NA and drop out on the filter
#    below; coloc-untested pairs were excluded from the pool entirely by 28) ──
pool <- fread(POOL)
# The pool's outcome_id is what run_susie_pair() dispatches on. 28 validates its own inputs against
# OUTCOMES, but this script must not assume that: a pool written by an older 28, or by a run under a
# re-pointed config, would dispatch the wrong outcome and every check here would still pass on internal
# consistency alone (r22 High). Validate at the boundary you read, not the one someone else read.
# `gene~outcome_name` is this script's KEY -- for pool_keys, for `done`, and for results[[key]]. A
# duplicate pool row is therefore silently skipped once the first copy is marked done, and only one record
# survives into the CSV. 28 and 30 both guard uniqueness, so a duplicate would eventually be caught -- but
# only AFTER 29 had written an incomplete mode file and a checkpoint keyed on it (r49 Medium). Check the
# key where the key is used.
if (anyDuplicated(pool, by = c("gene", "outcome_name")))
  stop(sprintf("29: candidate_pool.csv has %d duplicate (gene, outcome_name) key(s) -- this script keys its work, its checkpoint and its output on that pair, so a duplicate would be silently dropped rather than tested",
               sum(duplicated(pool, by = c("gene", "outcome_name")))))
assert_outcome_truth(pool, "candidate_pool.csv")
# ...and the engine label is what the NEXT line dispatches on, so it must be checked too -- checking
# outcome_id and then trusting susie_engine is the same "checked thing is not the deciding thing" shape
# that r23 and r26 found elsewhere (r27 High). The engine is derivable from the outcome; a stored label
# that disagrees means one of them is stale.
assert_engine_truth(pool, "candidate_pool.csv")
# A TESTABLE row with no engine would silently vanish at the filter below rather than being tested. 28
# leaves engine NA only for non-testable rows (MHC / undetermined-MHC / no-coords); if one of those ever
# reaches here as testable, the pool and this script disagree about what is runnable.
# RUNNABILITY IS DERIVED, NOT READ. `susie_engine` and `is_testable` are both stored labels, and this
# script dispatches on them -- so trusting either is the "checked thing is not the deciding thing" shape
# for the fourth time (r29 High). 28's own definition is the truth:
#   is_testable := !is.na(is_mhc) & is_mhc == FALSE & has_coords == TRUE
# A first version only checked that testable rows are not engine-less. It never checked the INVERSE, and
# that is the dangerous direction: an MHC or no-coords row carrying an engine that happens to match its
# outcome_id passes assert_engine_truth(), survives the mode filter, and gets SuSiE-run -- the one thing
# 28 exists to prevent. The post-filter MHC guard below only catches is_mhc == TRUE, not NA, and never
# looks at has_coords.
# derive_runnability() (R/00_config.R) derives `.exp_testable` from the FACTS -- the coordinate columns
# and is_mhc_gene() -- and asserts the stored is_mhc / has_coords / is_testable labels against them.
# A first version here derived `.exp_testable` from `is_mhc` and `has_coords` directly, which are stored
# labels themselves: a derivation whose inputs are labels is still trust, one level down (r33 High).
pool <- derive_runnability(pool, "candidate_pool.csv")
# Engine present IFF runnable -- both directions.
pool[, .has_eng := !is.na(susie_engine) & susie_engine != ""]
.ei <- pool[.has_eng != .exp_testable]
if (nrow(.ei))
  stop(sprintf("29: %d pool row(s) where susie_engine presence does not match runnability: %s -- a runnable row without an engine would vanish at the mode filter untested; a NON-runnable row WITH an engine would be SuSiE-run, which is exactly what 28 withholds the engine to prevent",
               nrow(.ei), paste(utils::head(.ei[, sprintf("%s~%s(testable=%s, engine=%s)", gene, outcome_name,
                                                          .exp_testable, susie_engine)], 8), collapse = ", ")))
pool <- pool[.exp_testable == TRUE & susie_engine == MODE][, c(".exp_testable", ".has_eng") := NULL]
# Same shape as the is_testable guard above: `if ("is_mhc" %in% names(pool) && ...)` would go quiet if
# the column vanished, i.e. exactly when the pool stopped being one 28 produced. 28 always writes it.
if (!"is_mhc" %in% names(pool))
  stop("29: candidate_pool.csv has no is_mhc column -- 28 always writes it, so this pool is not one 28 produced")
if (nrow(pool[is_mhc == TRUE]))
  stop("29: an MHC candidate was assigned an engine by 28 -- refusing to SuSiE-test MHC (28 bug)")
cat(sprintf("29: %s testable candidates = %d\n", MODE, nrow(pool)))

pool_keys <- sort(paste(pool$gene, pool$outcome_name, sep = "~"))
# checkpoint provenance: key on the FULL run-defining inputs, not just mode+keys, so a pool regenerated with
# the same (gene,outcome) keys but changed outcome_id / PP.H4_abf / coords / engine source cannot silently
# resume stale SuSiE results (Codex review High: provenance guard was too weak).
sel_cols <- intersect(c("gene", "outcome_name", "outcome_id", "PP.H4_abf",
                        "gene_chrom", "gene_start", "gene_end", "susie_engine"), names(pool))
file_sig <- function(f) if (file.exists(f)) list(size = file.size(f), mtime = as.numeric(file.mtime(f))) else NULL
eng_file <- if (MODE == "egfr") file.path(B, "scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R") else
                                file.path(B, "scripts/paper3.1_ckm/phase3/15_susie_coloc.R")
# CODE identity, by md5 (r12 High, 2026-07-16). Without this a resume is invisible downstream: the
# checkpoint replays rows computed under OLDER shared code, this script writes a CSV with a BRAND-NEW
# mtime, and 33's dependency gate -- which can only see mtimes -- passes it as if current code produced
# it. The gate made that failure LOOK checked, which is worse than leaving it unchecked.
# md5, not file_sig's size+mtime: content is the whole point here, these are small files, and the cheap
# proxy is exactly what the hole was made of. Must mirror 33's .dep for this mode: 29 sources R/coloc.R
# ONLY on the non-eGFR branch (the eGFR path evals 26's head, and 26 sources 00_config alone).
code_deps <- c(file.path(B, "R/00_config.R"),
               file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"),
               file.path(B, "scripts/paper3.1_ckm/phase3/29_candidate_susie.R"),
               file.path(B, "R/strand_align.R"),   # run-defining: the eval'd engine head sources it (cold re-review Low)
               eng_file,
               if (MODE == "nonegfr") file.path(B, "R/coloc.R"))
code_md5 <- vapply(code_deps, function(f) if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_, "")
if (anyNA(code_md5))
  stop(sprintf("29: cannot md5 run-defining code: %s -- refusing to key a checkpoint on an unknown code state",
               paste(basename(code_deps[is.na(code_md5)]), collapse = ", ")))
# SESSION dimension (r13 High, partial). The engines also read runtime DATA this signature does not
# cover -- the cis-rsid map, outcome GWAS, the UKB-PPP tars reached via resolve_ukbppp_tar(), and the LD
# reference panels. Those are GB-scale; hashing them per resume is not proportionate, and size+mtime over
# ~150 tars would be a cheap proxy of exactly the kind this whole exercise keeps finding insufficient.
# So state the boundary instead of implying coverage: this checkpoint is keyed on the POOL, the CODE
# (by md5) and the SESSION (R + package versions). It is NOT keyed on the reference data. A resume after
# the underlying data changed is not detected here -- and, mirroring phase3/33's "does NOT pin" section,
# saying so is the honest option. Tracked as REBUILD_LEDGER §2.1c(d).
sess <- list(r = R.version.string,
             pkg = vapply(c("data.table", "susieR", "coloc", "ieugwasr"),
                          function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), ""))
if (anyNA(sess$pkg))
  stop(sprintf("29: cannot read version for %s -- refusing to key a checkpoint on an unknown session",
               paste(names(sess$pkg)[is.na(sess$pkg)], collapse = ", ")))
prov <- list(mode = MODE, keys = pool_keys, sess = sess,
             rows = pool[order(gene, outcome_name), ..sel_cols],
             eng_sig = file_sig(eng_file), pool_sig = file_sig(POOL),
             code_md5 = code_md5)   # checkpoint provenance

results <- list(); done <- character(0); loaded_ckpt <- FALSE
if (nzchar(Sys.getenv("CKM_PALINDROME_SINK", "")) && file.exists(CKPT)) {
  file.remove(CKPT)   # palindrome capture demands a fresh run so EVERY pair calls strand_align_R -> complete sink (cold re-review High: a resumed pair skips capture)
  cat("29: palindrome capture on -- discarded checkpoint for a fresh, fully-captured run\n")
}
if (file.exists(CKPT)) {
  # THE RAW DATA IS PART OF THE RESUME DECISION, by mtime rather than by hash. The provenance above
  # deliberately leaves the outcome extracts, the ~150 UKB-PPP tars, the cis-rsID map and the LD
  # reference out of the signature because hashing GB-scale data on every resume is not
  # proportionate. That reasoning stands. What it left open is not hypothetical: a run aborts, one
  # of those inputs is corrected or replaced, the operator resumes exactly as this script intends,
  # and the finished CSV mixes fits computed from different raw data -- valid keys, valid schema,
  # current timestamp, and no record of which row used which. 30-32 accept it.
  #
  # mtime answers the one question that matters here -- was any of it REPLACED since these results
  # were computed -- and a stat is not a hash, so the proportionality objection does not apply. It
  # is paid only on a resume; a fresh run never reaches this branch.
  # OUTCOME_DIR, NOT LOCAL_EXTRACT_DIR. The first version of this guard stat-ed the localized
  # outcome RDS files, which these engines do not read: the non-eGFR path opens the full GWAS via
  # 15_susie_coloc.R:210 and R/coloc.R:719/781/847, and 26 reads OUTCOME_DIR/CKDGen_eGFR_EA_hg37.txt.gz
  # directly. All eight full outcome files sat outside the set, so replacing one and resuming was
  # exactly as undetected as before. Top level only -- the subdirectories under OUTCOME_DIR are not
  # read on this path.
  .odir <- list.files(OUTCOME_DIR, full.names = TRUE)
  .raw <- c(file.path(BASE_DIR, "cache", "cis_rsid_map.tsv"),
            .odir[!dir.exists(.odir)],
            list.files(file.path(DATA_DIR, "ukb_ppp", "european"), full.names = TRUE),
            list.files(LD_EUR_DIR, full.names = TRUE))
  .raw <- .raw[file.exists(.raw)]
  .ck_t <- file.mtime(CKPT)
  .newer <- .raw[file.mtime(.raw) > .ck_t]
  if (length(.newer))
    stop(sprintf("29: %d raw input file(s) postdate the checkpoint written %s -- e.g. %s. Resuming would skip the completed pairs and leave the output mixing fits computed from different data, with nothing in it recording which row used which. Delete %s for a fresh run.",
                 length(.newer), format(.ck_t, "%Y-%m-%d %H:%M"),
                 paste(basename(head(.newer, 3)), collapse = ", "), basename(CKPT)))
  cat(sprintf("29: %d raw input file(s) all predate the checkpoint\n", length(.raw)))

  cp <- readRDS(CKPT)
  if (is.null(cp$prov) || !identical(cp$prov$mode, MODE) || !identical(cp$prov$keys, pool_keys) ||
      !identical(cp$prov$rows, prov$rows) || !identical(cp$prov$eng_sig, prov$eng_sig) ||
      !identical(cp$prov$pool_sig, prov$pool_sig) || !identical(cp$prov$code_md5, prov$code_md5) ||
      !identical(cp$prov$sess, prov$sess))
    stop(sprintf("29: checkpoint %s is for a DIFFERENT candidate set/mode/engine-version/code state (pool rows, engine source, pool file, or a run-defining code md5 changed) -- delete it for a fresh run (provenance guard). NOTE: checkpoints written before 2026-07-16 carry no code_md5 and are rejected here by design.", CKPT))
  results <- cp$results; done <- cp$done; loaded_ckpt <- TRUE
  # Restore the RNG stream so a resumed run is stream-identical to an uninterrupted one (r13 Medium).
  # Checkpoints written before 2026-07-16 carry no rng -- but they are already rejected above for lacking
  # code_md5, so this can only be NULL on a checkpoint that failed the provenance guard, never here.
  if (is.null(cp$rng))
    stop(sprintf("29: checkpoint %s carries no RNG state -- it predates the 2026-07-16 resume fix; delete it for a fresh run", CKPT))
  assign(".Random.seed", cp$rng, envir = .GlobalEnv)
  cat(sprintf("29: resuming -- %d already done\n", length(done)))
}

run_one <- function(row) {
  if (MODE == "egfr")
    run_susie_egfr_pair(row$gene, row$gene_chrom, row$gene_start, row$gene_end, pp4_abf = row$PP.H4_abf)
  else
    run_susie_pair(row$gene, row$outcome_id, row$gene_chrom, row$gene_start, row$gene_end, rsid_map, pp4_abf = row$PP.H4_abf)
}
# Capture each pair's engine output (cat -> stdout, message -> stderr) so a NULL result can be classified
# benign-no-rescue vs genuine computation/data failure (editor middle-path). Live output is muffled and
# re-emitted to the run log afterwards, preserving the human-readable trace.
run_one_capture <- function(row) {
  msg_lines <- character(0)
  box <- new.env()
  out_lines <- capture.output({
    box$res <- withCallingHandlers(
      run_one(row),
      message = function(m) { msg_lines <<- c(msg_lines, sub("\n$", "", conditionMessage(m))); invokeRestart("muffleMessage") }
    )
  }, type = "output")
  list(res = box$res, log = c(out_lines, msg_lines))
}
# `converged_*` default to NA here, and the caller overwrites them when it knows better. This row is used
# for pairs where NO fit ran (the engine returned before fitting), so NA is the honest value -- FALSE would
# assert a fit that failed to converge, which is a different and false claim.
mk_na_row <- function(row, status)
  data.table(gene = row$gene, outcome_name = row$outcome_name, outcome_id = row$outcome_id,
             PP.H4_abf = row$PP.H4_abf, PP.H4_susie = NA_real_, PP.H3_susie = NA_real_,
             concordant_susie = NA, run_status = status, error_msg = NA_character_,
             converged_exp = NA, converged_out = NA)

# cold-review 29-M1: make delegated susieR calls reproducible across reruns.
# 2026-07-16 (r13 Medium): seeding ONCE here is not resume-stable -- a resume skips completed rows, so the
# remaining rows would draw from a different RNG stream than an uninterrupted run. Fixed by persisting
# .Random.seed in the checkpoint and restoring it on resume (see the CKPT block above), which makes a
# resumed run stream-identical to an uninterrupted one.
# NOT fixed by switching to 31's per-edge seed: that would change the uninterrupted-run stream too, and
# the pending re-run uses byte-identity of the retained rows against _REPROTEST_BASELINE_20260715/ as its
# determinism check. Preserve the semantics, fix the resume.
# Measured 2026-07-16: susie_rss (susieR 0.14.2) does not touch .Random.seed and returns identical fits
# under different seeds -- so on this path the whole question is currently inert. That is a property of the
# installed estimator, not a guarantee of the coloc.susie wrapper around it, so the guard stays.
#
# CONDITIONAL, and that is the whole point (r14 High): a first version of this fix restored cp$rng on
# resume and then ran an UNCONDITIONAL set.seed(1) forty lines later, silently clobbering the restore.
# Both halves were present; only their ORDER made the fix a no-op. Seed the stream ONLY when starting
# fresh -- on a resume the stream was already restored from the checkpoint above.
# Seed ONLY when no checkpoint was loaded. `!length(done)` was the first attempt and is the WRONG
# predicate (r16 High): a checkpoint saved by an abort on the FIRST pair carries done = character(0), so
# the restore above would be clobbered right here -- r14's exact failure mode, one layer down.
if (!loaded_ckpt) set.seed(1)
for (i in seq_len(nrow(pool))) {
  row <- pool[i]; key <- paste(row$gene, row$outcome_name, sep = "~")
  if (key %in% done) next
  cat(sprintf("\n[%d/%d] %s ~ %s (coloc.abf PP.H4=%.3f)\n", i, nrow(pool), row$gene, row$outcome_name, row$PP.H4_abf))
  flush(stdout())
  # RNG state BEFORE this pair consumes anything. The abort handlers below save the checkpoint with the
  # pair EXCLUDED from `done` -- it is retried on resume -- so they must persist the PRE-pair state.
  # Saving .Random.seed at abort time records the state AFTER the failing pair already advanced it, and
  # the retry would then start from a stream position an uninterrupted run never occupies (r16 High).
  # The success save at the bottom of the loop is different and correct: there the pair IS in `done`, so
  # the post-pair state is exactly the right thing to persist.
  .rng_pre <- .Random.seed

  cap <- tryCatch(run_one_capture(row), error = function(e) {
    saveRDS(list(prov = prov, results = results, done = done, rng = .rng_pre), CKPT)   # pair excluded -> retried, so persist the PRE-pair RNG
    stop(sprintf("29 ABORT at %s: %s\n  A candidate SuSiE failure that THREW = missing data / bad LD / code drift / infra, NOT a benign skip. Refusing to demote it to 'uncomputable'. Fix and resume from checkpoint (this pair is retried).",
                 key, conditionMessage(e)))
  })
  res <- cap$res
  if (length(cap$log)) cat(paste(cap$log, collapse = "\n"), "\n")   # re-emit captured engine output to the run log
  flush(stdout())

  # Belt-and-suspenders for the fail_loud_data() session-flag mechanism (00_config.R:33): a data-foundation
  # failure sets .CKM_DATA_FAILURE then stop()s, but the config header warns a re-raised classed condition can
  # be demoted to NULL by an inner sibling handler. Mirror the canonical loops (15:437 / 26:410): abort if set.
  if (exists(".CKM_DATA_FAILURE", envir = .GlobalEnv)) {
    saveRDS(list(prov = prov, results = results, done = done, rng = .rng_pre), CKPT)   # pair excluded -> retried, so persist the PRE-pair RNG
    stop(sprintf("29 ABORT at %s: engine flagged a data-foundation failure (.CKM_DATA_FAILURE): %s\n  NOT a benign skip. Fix and resume from checkpoint (this pair is retried).",
                 key, get(".CKM_DATA_FAILURE", envir = .GlobalEnv)))
  }

  if (is.null(res)) {
    # Reason-aware NULL handling (editor middle-path): a NULL is a benign no-rescue ONLY if its captured reason
    # matches a known-benign signature; a coloc.susie error / missing-column / UNRECOGNISED reason is FATAL and
    # ABORTS (never frozen as 'uncomputable'). classify_null() defaults unknown reasons to FATAL.
    cl <- classify_null(cap$log)
    if (cl$status == "FATAL") {
      saveRDS(list(prov = prov, results = results, done = done, rng = .rng_pre), CKPT)   # pair excluded -> retried, so persist the PRE-pair RNG
      stop(sprintf("29 ABORT at %s: SuSiE returned NULL for a NON-benign reason (%s): %s\n  A coloc.susie computation error / missing-column / unrecognised skip is NOT a benign no-rescue. Refusing to freeze it as 'uncomputable'. Fix and resume (this pair is retried).",
                   key, cl$tag, cl$detail))
    }
    rec <- mk_na_row(row, "skip_benign"); rec[, error_msg := cl$tag]
  } else if ("PP.H4_susie" %in% names(res) && nrow(as.data.table(res)) == 1L &&
             is.na(as.data.table(res)$PP.H4_susie)) {
    # A POST-FIT COLLAPSE NOW ARRIVES AS A ROW, NOT AS NULL. 15/26 were changed on 2026-07-16 so that a
    # coloc.susie failure or empty summary AFTER both susie_rss calls succeed still returns the fits'
    # convergence rather than a bare NULL. Right for susie_coloc_results -- and it breaks the contract THIS
    # branch relied on, exactly as it broke 31's: "non-NULL" no longer implies "usable", so such a row would
    # fall into the else-branch below and abort at its posterior-integrity check (`!is.finite(NA)` is TRUE).
    # The engine's NULL now means only "returned before fitting", so this shape has to be classified here.
    # Mirror the NULL branch above exactly -- same classify_null() input, same FATAL-aborts-and-checkpoints
    # rule -- and keep the convergence the row now carries, which is the whole reason the row exists.
    .r <- as.data.table(res)
    cl <- classify_null(cap$log)
    if (cl$status == "FATAL") {
      saveRDS(list(prov = prov, results = results, done = done, rng = .rng_pre), CKPT)
      stop(sprintf("29 ABORT at %s: SuSiE collapsed after fitting for a NON-benign reason (%s): %s\n  A coloc.susie computation error / missing-column / unrecognised skip is NOT a benign no-rescue. Refusing to freeze it as 'uncomputable'. Fix and resume (this pair is retried).",
                   key, cl$tag, cl$detail))
    }
    # CARRY EVERYTHING THE ROW KNOWS, NOT JUST THE CONVERGENCE. mk_na_row() is built for the pre-fit case,
    # where nothing is known; on THIS path the fits ran, so the engine's row already carries the SNP and
    # credible-set counts and they are facts. Rebuilding from mk_na_row and copying only the convergence
    # writes those known facts as NA -- the same "the recording does not reach the artifact" bug this whole
    # change exists to fix, one field-set over. It also defeats susie_fit_ran() (R/00_config.R), which uses
    # exactly these columns as the evidence a fit ran: blank them and the honesty check goes inert on the
    # very rows it was extended to cover.
    rec <- mk_na_row(row, "skip_benign")
    for (.c in intersect(c("n_snps_merged", "n_snps_ld", "n_signals_exp", "n_signals_out", "n_signal_pairs"),
                         names(.r)))
      rec[, (.c) := .r[[.c]]]
    rec[, `:=`(error_msg = cl$tag,
               converged_exp = isTRUE(.r$converged_exp), converged_out = isTRUE(.r$converged_out))]
  } else {
    r <- as.data.table(res); if (!"outcome_name" %in% names(r)) r[, outcome_name := row$outcome_name]
    # A non-NULL return is NOT proof of a usable result (r15 High, 2026-07-16). Without this, an engine or
    # schema change that drops PP.H4_susie sails straight through: the fixed-schema block near the write
    # back-fills the column as NA, the row still reads run_status=="ok", and 30/32 then read that NA as
    # "ran fine, did not reach 0.80" -- a COMPUTATION failure reported as biological evidence. That is
    # exactly the type-(b)-must-not-look-benign rule this layer was built to enforce, arriving through the
    # ok path instead of the NULL path. Measured 2026-07-16 on the frozen output: 129 ok rows, 0 with a
    # missing or NA PP.H4_susie -- so this costs nothing today.
    .need <- c("PP.H4_susie", "PP.H3_susie")
    .miss <- setdiff(.need, names(r))
    # Check NA over EVERY required column, not just the first. A first version listed both and tested only
    # `anyNA(r$PP.H4_susie)` -- so a return carrying PP.H3_susie = NA was still stamped ok (r16 High).
    # The incomplete-fix pattern, occurring inside the fix for the incomplete-fix pattern.
    # Presence is not integrity: `Inf >= 0.80` is TRUE, so a malformed posterior would be COUNTED as a
    # rescue by this script's own summary and read as one by 30/31 (r39 High). assert_pp_integrity()
    # (R/00_config.R) is the single source for this -- three rounds of writing it per-site produced three
    # checks that did not line up with the six thresholds they were meant to protect.
    .nacol <- if (length(.miss)) character(0) else
              .need[vapply(.need, function(.c) any(!is.finite(r[[.c]]) | r[[.c]] < 0 | r[[.c]] > 1), TRUE)]
    if (length(.miss) || length(.nacol)) {
      saveRDS(list(prov = prov, results = results, done = done, rng = .rng_pre), CKPT)   # pair excluded -> PRE-pair RNG
      stop(sprintf("29 ABORT at %s: engine returned a non-NULL result that is NOT usable -- %s. This is a computation/schema failure, not a benign no-rescue, and must never be frozen as run_status=='ok'.",
                   key, if (length(.miss)) paste("missing column(s):", paste(.miss, collapse = ", "))
                        else paste("missing, non-finite or out-of-range value in required column(s):",
                                   paste(.nacol, collapse = ", "))))
    }
    r[, run_status := "ok"]; r[, error_msg := NA_character_]; rec <- r
  }
  results[[key]] <- rec; done <- c(done, key)
  if (i %% 5 == 0 || i == nrow(pool)) { saveRDS(list(prov = prov, results = results, done = done, rng = .Random.seed), CKPT)
    cat(sprintf("  [checkpoint: %d done]\n", length(done))) }
  flush(stdout())
}

# ── Write output with a FIXED schema even when a mode has zero candidates (cold-review 29-M2) ──
final <- if (length(results)) rbindlist(results, fill = TRUE, use.names = TRUE) else data.table()
for (col in TEMPLATE_COLS)
  if (!col %in% names(final))
    final[, (col) := if (col %in% c("gene", "outcome_name", "outcome_id", "susie_engine", "run_status", "error_msg")) NA_character_ else NA]
# Batch-level integrity before the write. The per-edge assertion above runs inside run_one_capture's
# caller; this is the guarantee the FILE carries, and it is what 30/31/32 will trust (r39).
assert_pp_integrity(final, c("PP.H4_susie", "PP.H3_susie"), sprintf("candidate_susie_%s.csv (about to be written)", MODE))
assert_status_premise(final, c("PP.H4_susie", "PP.H3_susie"), sprintf("candidate_susie_%s.csv (about to be written)", MODE))
final[, susie_engine := MODE]   # stamp engine/mode so 30 can verify file-to-engine membership (Codex review Low)
setcolorder(final, intersect(TEMPLATE_COLS, names(final)))
fwrite(final, OUT_CSV)

n_ok   <- sum(final$run_status == "ok", na.rm = TRUE)
n_resc <- sum(final$PP.H4_susie >= 0.80, na.rm = TRUE)
cat(sprintf("\n29 [%s]: %d run | ok=%d | skip_benign=%d | rescued(PP.H4_susie>=0.80)=%d\n",
            MODE, nrow(final), n_ok, sum(final$run_status == "skip_benign", na.rm = TRUE), n_resc))
cat(sprintf("29: wrote %s\n", OUT_CSV))
if (n_resc > 0) {
  cat("\n== rescued candidates (PP.H4_susie >= 0.80) ==\n")
  print(final[PP.H4_susie >= 0.80, .(gene, outcome_name, PP.H4_abf, PP.H4_susie, PP.H3_susie,
                                     n_signals_exp, n_signals_out)][order(-PP.H4_susie)])
}
