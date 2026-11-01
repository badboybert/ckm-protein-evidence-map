# ============================================================================
# Stage 0c (dashed candidate layer): FREEZE the candidate-edge audit table.
#
# Combines the candidate pool (28) with the two SuSiE runs (29 egfr + nonegfr) and
# applies the LOCKED classification rule (REQUIREMENTS.md 2.2):
#   dashed_edge / susie_rescued      <- run_status=="ok" AND coloc_status=="coloc_fail" AND PP.H4_susie>=0.80
#                                       (a real SuSiE multi-signal rescue). ONLY these become dashed edges.
#   not_evaluable / MHC_excluded     <- is_mhc==TRUE (never a dashed edge).
#   excluded / undetermined_mhc      <- is_mhc==NA.
#   excluded / no_coords             <- missing instrument coords.
#   excluded / susie_not_supportive  <- run_status=="ok" AND PP.H4_susie<0.80.
#   excluded / no_susie_rescue:<tag> <- run_status=="skip_benign" (benign no-rescue: too few SNPs / empty
#                                       credible set / empty coloc summary). 29 is reason-aware and ABORTS on
#                                       any real computation/data failure, so a coloc.susie error never lands here.
#   excluded / susie_review:<status> <- defensive catch-all. It USED to also catch an "ok" row with no
#                                      PP.H4_susie; the invariant gates now stop before that reaches
#                                      classification, so that arm should be unreachable (r43). Any row
#                                      that still lands here is a rule gap worth investigating, not a
#                                      benign exclusion. (e.g. an "ok"
#                                       row lacking PP.H4_susie -- an invariant/schema drift). Surfaced as a
#                                       distinct, non-benign reason so such a case is visible, never hidden.
#   EVERY candidate appears exactly once.
#
# Guards (cold-review 30): required SuSiE columns asserted; every TESTABLE candidate
#   MUST have a SuSiE row (a missing run is an incomplete/stale 29 run -> stop, not a
#   silent exclusion); dashed requires run_status=="ok"+coloc_fail (a stale/failed row
#   with a stray PP.H4_susie>=0.80 cannot become dashed); pool asserted unique.
#
# NON-DESTRUCTIVE: writes ONE new file candidate_edge_audit.csv.
# Run: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/30_candidate_edge_audit.R
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
Sys.setenv(CKM_OFFLINE = "1")
suppressPackageStartupMessages(library(data.table))
B  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
# 30 PRODUCES the audit every downstream script trusts, so it is the last place that should be taking
# labels on faith -- and until 2026-07-16 it sourced no config at all, so it could not call the truth
# guards even in principle (r33 High). Load config BEFORE anything is read or normalised.
source(file.path(B, "R/00_config.R"))   # assert_outcome_truth/assert_engine_truth/derive_runnability/is_coloc_pass
# MOVED BELOW THE CONFIG SOURCE. These two lines stood ABOVE it and used PHASE1_TAG, which
# R/00_config.R is what defines -- so a clean run died here with "object 'PHASE1_TAG' not found"
# before reading anything. 28 has the order right (source at :19/:23, RC at :25).
RC <- file.path(B, file.path("results", paste0("repro_", PHASE1_TAG)))
D  <- file.path(RC, "phase3_susie")
source(file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"))   # CAND_BENIGN_TAGS: single source (Codex round-2 Low: no duplication)

POOL <- file.path(D, "candidate_pool.csv")
S_EG <- file.path(D, "candidate_susie_egfr.csv")
S_NE <- file.path(D, "candidate_susie_nonegfr.csv")
OUT  <- file.path(D, "candidate_edge_audit.csv")
# Same reason as 29's and 31's; 31 and 32 both read this file and neither can tell its age.
unlink(OUT); if (file.exists(OUT)) stop(sprintf("30: could not remove the previous %s -- refusing to run with a stale result in place", OUT))
for (.f in c(POOL, S_EG, S_NE))
  if (!file.exists(.f)) stop(sprintf("30: missing input %s -- run 28 then 29 (both modes) first", .f))

pool <- fread(POOL)
# Truth-check the pool BEFORE the engine normalisation below rewrites non-{egfr,nonegfr} values to NA --
# checking after that is checking the repair, not the artifact (the shape r32 found in 31). 30 gained the
# ability to call these in r33 but did not call them; sourcing a guard is not invoking it (r34 Blocker).
assert_outcome_truth(pool, "candidate_pool.csv")
assert_engine_truth(pool, "candidate_pool.csv")
if (anyDuplicated(pool, by = c("gene", "outcome_name")))   # cold-review 30-M1
  stop("30: candidate pool has duplicate (gene,outcome_name) keys -- refusing (28 should have asserted uniqueness)")
# fread reads an empty CSV field as "" (not NA). 28 writes susie_engine=NA for non-testable rows
# (MHC / undetermined-MHC / no-coords), which round-trips to "". Normalise non-{egfr,nonegfr} back to NA so
# is.na(susie_engine) is correct everywhere below (esp. the missing_run guard) and the frozen susie_engine is
# clean. RUNTIME bug caught on first execution 2026-07-12 (static review missed the NA/CSV round-trip).
pool[!susie_engine %in% c("egfr", "nonegfr"), susie_engine := NA_character_]

# ── SuSiE runs: assert schema, keep outcome_id for provenance, require disjoint modes ──
# PP.H3_susie belongs here: an `ok` row missing it is unusable SuSiE output, and 30's job is to REFUSE
# a stale 29 file -- 29's own producer-side assertion is not a substitute for the consumer checking
# (r37 High). Without it, such a row falls to the `susie_review:ok` catch-all and is frozen as a benign
# `excluded`, which is a computation failure recorded as a biological one.
# error_msg belongs here: it CARRIES the benign-skip reason this file gates on and classifies with, so a
# stale 29 output missing the column was accepted at the read boundary and only failed later, or worse,
# made the benign-reason check conditional on the column it is meant to check (r47 Medium).
REQ_SUSIE <- c("gene", "outcome_name", "outcome_id", "susie_engine", "PP.H4_susie", "PP.H3_susie",
               "run_status", "error_msg",
               # Required since 2026-07-16: 29 records per-side convergence, and 30 must refuse a
               # candidate_susie_*.csv that predates it -- otherwise 30 (and 31 after it) consume a
               # half-stale file and 33's mtime gate only notices after the bad inputs are already spent.
               "converged_exp", "converged_out")
read_susie <- function(f, expect_mode) {
  s <- fread(f)
  miss <- setdiff(REQ_SUSIE, names(s))
  if (length(miss)) stop(sprintf("30: %s missing required column(s): %s -- malformed/wrong-version 29 output", basename(f), paste(miss, collapse = ", ")))  # cold-review 30-H2
  # file-to-engine membership: susie_engine is REQUIRED (in REQ_SUSIE), so a stale/pre-engine-column file fails
  # the check above; here assert each mode file contains ONLY its own engine's rows (Codex round-2 Low hardening).
  if (nrow(s)) {
    bad <- s[is.na(susie_engine) | susie_engine != expect_mode]
    if (nrow(bad)) stop(sprintf("30: %s contains %d row(s) with susie_engine != '%s' -- swapped/stale mode file", basename(f), nrow(bad), expect_mode))
  }
  # PRESENCE was checked above; CORRECTNESS here -- a file tainted by an earlier fill has the columns and
  # NA values (susie_fit_ran / assert_convergence_honest, R/00_config.R). Before the rename, while the
  # engine's own column names are intact.
  assert_convergence_honest(s, basename(f))
  keep <- intersect(c(REQ_SUSIE, "susie_engine", "PP.H3_susie", "n_signals_exp", "n_signals_out", "n_signal_pairs",
                      "n_snps_merged", "n_snps_ld", "error_msg"), names(s))
  s <- s[, ..keep]
  # Truth-check BEFORE the rename: after `outcome_id` becomes `outcome_id_susie` the shared guards no
  # longer recognise the column, so this is the only window in which they can run at all (r40 High).
  assert_outcome_truth(s, basename(f))
  assert_engine_truth(s, basename(f))
  # A row is proved to EXIST by a sentinel stamped here, not by any value it happens to carry. r40 used
  # `!is.na(run_status)` for that, which is a VALUE from the row -- so a stale row with real PP and a
  # missing run_status read as "no row at all", escaped every downstream check, and could be frozen into
  # the audit under not_evaluable/excluded carrying its PP (r41 High). Existence is structural; ask the
  # structure.
  s[, .has_susie_row := TRUE]
  setnames(s, "outcome_id", "outcome_id_susie")
  if ("susie_engine" %in% names(s)) setnames(s, "susie_engine", "susie_engine_file")
  s
}
susie <- rbindlist(list(read_susie(S_EG, "egfr"), read_susie(S_NE, "nonegfr")), fill = TRUE, use.names = TRUE)
.dup <- susie[, .N, by = .(gene, outcome_name)][N > 1]
if (nrow(.dup))
  stop(sprintf("30: (gene,outcome) present in BOTH SuSiE mode files -- engine partition violated for %s",
               paste(.dup[, paste(gene, outcome_name)], collapse = ", ")))

# orphan guard: every SuSiE row MUST correspond to a current-pool candidate. merge(all.x) would silently
# drop stale rows, hiding a mismatched/old 29 output (Codex review Medium). Refuse to freeze if any exist.
orphan <- susie[!pool, on = .(gene, outcome_name)]
if (nrow(orphan))
  stop(sprintf("30: %d SuSiE row(s) not in the candidate pool (stale/mismatched 29 output): %s -- rerun 29 against the current pool before freezing",
               nrow(orphan), paste(head(orphan[, paste0(gene, "~", outcome_name)], 12), collapse = ", ")))

aud <- merge(pool, susie, by = c("gene", "outcome_name"), all.x = TRUE)
# all.x leaves the sentinel NA for pool rows with no SuSiE row -- that IS the "no row" case, and it is the
# only thing that means it. Normalise it to FALSE so every downstream test reads a boolean, never an NA.
aud[is.na(.has_susie_row), .has_susie_row := FALSE]

# outcome_id provenance: SuSiE file's id must match the pool's frozen id (cold-review M)
if ("outcome_id_susie" %in% names(aud)) {
  # NA-AWARE, BUT ONLY WHERE A SuSiE RESULT EXISTS. `!is.na(x) & !is.na(y) & x != y` reads a one-sided NA
  # as agreement, so a stale 29 row with the right (gene, outcome_name), real PP and run_status, but NO
  # outcome_id would pass and then be frozen into the audit UNDER THE POOL'S id -- the next line drops the
  # SuSiE id and the disagreement with it (r40 High).
  # The scope matters and a first version got it wrong: the 17 non-runnable candidates have no SuSiE row
  # at all, so `merge(all.x = TRUE)` leaves their outcome_id_susie legitimately NA. Applying the NA-aware
  # predicate to every row turned "this pair was never run" into "these two disagree" and stopped 30 on
  # the frozen data. The gate is the STRUCTURAL sentinel `.has_susie_row`, not `!is.na(run_status)` --
  # that was the r40 attempt, and r41 found it asks a payload value whether its own row exists.
  .ne <- function(x, y) xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
  mism <- aud[.has_susie_row == TRUE & .ne(outcome_id, outcome_id_susie)]
  if (nrow(mism))
    stop(sprintf("30: SuSiE outcome_id disagrees with the pool's (a one-sided NA counts) for %s -- stale/mismatched SuSiE file",
                 paste(mism[, sprintf("%s~%s(pool=%s, susie=%s)", gene, outcome_name,
                                      outcome_id, outcome_id_susie)], collapse = ", ")))
  aud[, outcome_id_susie := NULL]
}

# engine provenance: a SuSiE row's engine must match the pool's assigned engine (Codex review Low)
if ("susie_engine_file" %in% names(aud)) {
  emism <- aud[!is.na(susie_engine_file) & !is.na(susie_engine) & susie_engine_file != susie_engine]
  if (nrow(emism))
    stop(sprintf("30: SuSiE engine != pool susie_engine for %s -- wrong-engine/stale SuSiE file",
                 paste(emism[, paste(gene, outcome_name)], collapse = ", ")))
  aud[, susie_engine_file := NULL]
}

# RUNNABILITY IS DERIVED, NOT READ -- and this file matters most, because it PRODUCES the audit every
# downstream script trusts (r32 High). "28 assigned an engine" was the old proxy for "testable", so a
# runnable edge whose engine label went missing had `!is.na(susie_engine)` FALSE, skipped this stop
# entirely, and fell through to the `susie_review:no_run` catch-all below -- silently classified
# `excluded` instead of halting the freeze. 28's rule is the truth:
#     is_testable := !is.na(is_mhc) & is_mhc == FALSE & has_coords == TRUE
# Frozen data 2026-07-16: 0 runnable-without-engine, 0 non-runnable-with-engine, 0 rows reasoned
# `susie_review:no_run` -- so no edge has ever been lost this way, and this costs nothing today.
# Facts, not labels -- see derive_runnability() in R/00_config.R. This file PRODUCES the audit, so it is
# the last place that should take is_mhc / has_coords on faith (r33 High).
aud <- derive_runnability(aud, "candidate_edge_audit (pool-derived)")
# And the coloc-fail arm of the dashed rule is is_coloc_fail() -- evaluated AND not passing -- not the
# stored coloc_status label, and not the negation of the pass rule either (r37/r47). The
# downstream checks in 31/32 are not a substitute for the producer checking before it writes (r33 High).
assert_coloc_status_truth(aud, "candidate_edge_audit (pool-derived)")
# ...and CREATE the derived column the classification below reads. assert_coloc_status_truth() only
# ASSERTS -- it does not mutate -- so the r34 fix classified on a column that never existed and this
# script would have died at the dashed rule (r35 Blocker). It slipped through because the verification
# simulated 30's classification in a test harness that computed .coloc_fail_rd ITSELF: the test proved
# the LOGIC, not the SCRIPT. Mirrors 31.
# is_coloc_fail(), not !is_coloc_pass(): "unevaluable" must not read as "failed" (r37 High). Any row
# whose ABF probabilities are missing or out of range is caught by the evaluability stop below rather
# than silently counted as a coloc negative.
.une_aud <- aud[!is_coloc_evaluable(PP.H4_abf, PP.H3_abf)]
if (nrow(.une_aud))
  stop(sprintf("30: %d candidate(s) have missing or out-of-range ABF probabilities: %s -- they were never evaluated by coloc; classifying them as coloc-fail candidates would treat absent evidence as negative evidence",
               nrow(.une_aud), paste(head(.une_aud[, sprintf("%s~%s(H4=%s, H3=%s)", gene, outcome_name,
                                          format(PP.H4_abf), format(PP.H3_abf))], 8), collapse = ", ")))
aud[, .coloc_fail_rd := is_coloc_fail(PP.H4_abf, PP.H3_abf)]
# 28 admits ONLY coloc-fail pairs to the candidate pool -- that is the pool's definition. Checking that
# the stored coloc_status matches the ABF rule (above) does not check that every row actually FAILS: a
# coloc-PASSING row with a correctly-stored "coloc_pass" label passes that assertion, cannot become
# dashed, and is then frozen into this audit as `not_evaluable` or `excluded` -- putting a SOLID edge
# into the candidate-edge audit (r36 High). Frozen data: 151/151 are coloc-fail.
.cp_pool <- aud[.coloc_fail_rd != TRUE]
if (nrow(.cp_pool))
  stop(sprintf("30: %d candidate pool row(s) PASS the dual coloc.abf rule and do not belong in the candidate layer at all: %s -- 28 admits only coloc-fail pairs; this pool is stale or malformed",
               nrow(.cp_pool), paste(head(.cp_pool[, sprintf("%s~%s(H4=%s, H3=%s)", gene, outcome_name,
                                          format(PP.H4_abf), format(PP.H3_abf))], 8), collapse = ", ")))
.eng_mismatch <- aud[(!is.na(susie_engine) & susie_engine != "") != .exp_testable]
if (nrow(.eng_mismatch))
  stop(sprintf("30: %d candidate(s) where susie_engine presence does not match derived runnability: %s -- the engine label is stale, and every classification below keys off it",
               nrow(.eng_mismatch), paste(head(.eng_mismatch[, sprintf("%s~%s(testable=%s, engine=%s)", gene,
                                                outcome_name, .exp_testable, susie_engine)], 8), collapse = ", ")))
# Every RUNNABLE candidate MUST have a SuSiE row. A missing run is an incomplete/stale 29 run, NOT a
# biological negative -> stop (cold-review 30-H1).
# An `ok` SuSiE row with missing/malformed PP is unusable output, and the classification fallback would
# freeze it as `excluded / susie_review:ok` -- a computation failure recorded as a biological one. 29
# asserts this producer-side; 30's job is to REFUSE a stale 29 file, so it checks too (r37 High).
assert_pp_integrity(aud, c("PP.H4_susie", "PP.H3_susie"), "candidate_edge_audit (SuSiE-joined)")
assert_status_premise(aud[.has_susie_row == TRUE], c("PP.H4_susie", "PP.H3_susie"),
                      "candidate_edge_audit (rows that HAVE a SuSiE result)")
# A SuSiE row for a NON-RUNNABLE candidate is not an orphan by the pool-membership test -- the pair IS in
# the pool -- so it merges in and is frozen into the audit carrying PP and a run_status, for a pair 28
# says must never have been SuSiE-run at all (r39 High). Runnability, not pool membership, is the test.
.ran_unrunnable <- aud[.exp_testable != TRUE & .has_susie_row == TRUE]
if (nrow(.ran_unrunnable))
  stop(sprintf("30: %d NON-runnable candidate(s) carry a SuSiE result: %s -- 28 withholds the engine from MHC / undetermined-MHC / no-coords pairs precisely so they are never tested; a result here means a stale 29 output or a pool that changed under it",
               nrow(.ran_unrunnable), paste(head(.ran_unrunnable[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
missing_run <- aud[.exp_testable == TRUE & .has_susie_row == FALSE]
if (nrow(missing_run))
  stop(sprintf("30: %d testable candidate(s) have NO SuSiE row (29 run incomplete/stale): %s -- rerun 29 (both modes) to completion before freezing",
               nrow(missing_run), paste(head(missing_run[, paste0(gene, "~", outcome_name)], 12), collapse = ", ")))

# ── freeze-time integrity gate (editor middle-path) ──────────────────────────────────────────
# 29 is reason-aware: it ABORTS on any coloc.susie/missing-column/UNRECOGNISED NULL, so a COMPLETED 29
# run yields only run_status in {ok, skip_benign}, and every skip_benign carries a known-benign error_msg
# tag. Re-assert both here so a stale/pre-reason-aware 29 output (old "skip_no_susie") or a benign tag that
# drifted out of the shared vocabulary can NEVER be frozen -- a genuine computation failure must not hide as
# a benign "uncomputable" exclusion in the one artifact the map offers as proof of no moving goalposts.
bad_status <- aud[!is.na(run_status) & !run_status %in% c("ok", "skip_benign")]
if (nrow(bad_status))
  stop(sprintf("30: %d SuSiE row(s) have run_status not in {ok, skip_benign} (stale/pre-reason-aware 29 output): %s -- rerun 29 with the current reason-aware code before freezing",
               nrow(bad_status), paste(head(unique(bad_status$run_status), 5), collapse = ", ")))
if ("error_msg" %in% names(aud)) {
  bad_reason <- aud[run_status == "skip_benign" & (is.na(error_msg) | !error_msg %in% CAND_BENIGN_TAGS)]
  if (nrow(bad_reason))
    stop(sprintf("30: %d benign-skip row(s) carry an unknown/non-benign error_msg (a real failure could be hiding): %s -- 29 must ABORT+fix, not freeze",
                 nrow(bad_reason), paste(head(bad_reason[, paste0(gene, "~", outcome_name, "[", error_msg, "]")], 8), collapse = ", ")))
}

# ── risk-aligned MR direction (eGFR flipped to eGFR*, matching 12_rebuild:188-190) ──
# FINITE, not merely non-NA: `Inf < 0` is FALSE, so an infinite effect size becomes "risk" and rides
# into the audit as this edge's direction -- which 32 then reads as the figure's edge colour (r44 High).
.bad_b <- aud[!is.finite(b)]
if (nrow(.bad_b))
  stop(sprintf("30: %d candidate(s) carry a missing or non-finite MR effect size: %s -- mr_direction is derived from it and becomes the edge colour downstream",
               nrow(.bad_b), paste(head(.bad_b[, sprintf("%s~%s(b=%s)", gene, outcome_name, format(b))], 8),
                                   collapse = ", ")))
aud[, mr_beta_star := fifelse(outcome_name == "eGFR", -b, b)]
aud[, mr_direction := fifelse(is.na(mr_beta_star), NA_character_, fifelse(mr_beta_star < 0, "protective", "risk"))]

# ── classification (priority order) ──────────────────────────────────────────
aud[, classification := NA_character_]; aud[, reason := NA_character_]
# CLASSIFY FROM THE DERIVED FIELDS, not the stored labels. derive_runnability() above already asserted
# the labels against the facts, so on any passing run these are identical -- but reading the label here
# would mean the helper derived a truth and then the decision was made from the thing it was checking
# (r34 High). `.is_mhc_rd` NA = undetermined MHC, which 28 treats as not runnable; same semantics.
aud[.is_mhc_rd == TRUE,                        `:=`(classification = "not_evaluable", reason = "MHC_excluded")]
aud[is.na(classification) & is.na(.is_mhc_rd), `:=`(classification = "excluded", reason = "undetermined_mhc")]
aud[is.na(classification) & .has_coords_rd == FALSE, `:=`(classification = "excluded", reason = "no_coords")]
# dashed REQUIRES a successful run on a genuine coloc-fail candidate (cold-review 30-H3)
aud[is.na(classification) & run_status == "ok" & .coloc_fail_rd == TRUE &
      !is.na(PP.H4_susie) & PP.H4_susie >= 0.80, `:=`(classification = "dashed_edge", reason = "susie_rescued")]
aud[is.na(classification) & run_status == "ok" & !is.na(PP.H4_susie) & PP.H4_susie < 0.80,
    `:=`(classification = "excluded", reason = "susie_not_supportive")]
# remaining unclassified = genuine benign no-rescue (skip_benign; 29 already ABORTED on any real failure).
# Label it explicitly with the benign tag -- never the old conflating "susie_uncomputable:skip_no_susie".
aud[is.na(classification) & run_status == "skip_benign",
    `:=`(classification = "excluded", reason = paste0("no_susie_rescue:", fifelse(is.na(error_msg) | !nzchar(error_msg), "benign", error_msg)))]
aud[is.na(classification), `:=`(classification = "excluded",
    reason = fifelse(!is.na(run_status), paste0("susie_review:", run_status), "susie_review:no_run"))]
if (aud[is.na(classification), .N])
  stop(sprintf("30: %d candidate(s) unclassified -- rule gap, refusing to freeze", aud[is.na(classification), .N]))

# THE TWO PER-PAIR ANCHORS (UMOD~eGFR, DPEP1~eGFR) STOOD HERE AND WERE REMOVED 2026-08-04 (PI).
#
# They pinned the exact classification each pair had in the 1,682-protein generation -- UMOD must be
# dashed_edge, DPEP1 must be excluded / susie_not_supportive / PP.H4_susie < 0.80 -- and stopped the
# run on any departure, deliberately: 'the tripwire is inverted, not removed, so the new expected
# state is still guarded'. On this generation DPEP1~eGFR came back dashed_edge at PP.H4_susie 0.9541
# with both fits converged, so the tripwire fired and demanded the investigation it was written for.
#
# The investigation found the premise wrong, not the result. The previous generation's gene
# coordinate table held GRCh37 positions while the UKB-PPP pQTL archive is GRCh38: lifting this
# generation's DPEP1 coordinates (chr16:89,613,308-89,638,456) to GRCh37 returns 89,679,716 --
# exactly the start the 1,682 generation used. Its cis window was therefore centred 66 kb from the
# gene, and the 0.0824 that the anchor encodes was fitted over that displaced window. 1,585 of the
# 1,589 shared genes have different coordinates between the two generations (median shift 274 kb,
# max 2.88 Mb); 50 of the 61 current strong genes and 148 of the 182 candidate-pool genes are among
# them. A per-pair expectation from that generation cannot be a tripwire for this one -- it would
# fire on every rebuild, which makes it a roadblock rather than a surprise.
#
# Re-pinning the new values was considered and rejected: freezing one run's answer is the same
# mistake with the opposite sign. What remains is the structural check above -- an unclassified
# candidate still stops the run (:272) -- because that guards a rule gap, not a remembered number.
#
# Neither pair is unguarded scientifically: both appear in the audit output and in the dashed layer
# 32 finalises, where their classification is visible and reviewable.

col_order <- intersect(c("gene", "outcome_name", "outcome_id", "susie_engine", "coloc_status",
                         "PP.H4_abf", "PP.H3_abf", "b", "mr_beta_star", "mr_direction", "pval",
                         "is_bonferroni", "is_mhc", "has_coords", "PP.H4_susie", "PP.H3_susie",
                         "n_signals_exp", "n_signals_out", "run_status", "error_msg", "classification", "reason"), names(aud))
# `col_order` is an explicit projection, so the derived helpers (.has_susie_row, .exp_testable,
# .is_mhc_rd, .has_coords_rd, .coloc_fail_rd) never reach the frozen CSV -- verified r34/r41.
aud <- aud[, ..col_order][order(factor(classification, levels = c("dashed_edge", "not_evaluable", "excluded")),
                                -PP.H4_susie, -PP.H4_abf)]
# The sort above puts every dashed row first, so a write interrupted part-way leaves a syntactically
# valid CSV holding the COMPLETE dashed set and only some of the later excluded/not-evaluable rows --
# 31 checks the dashed subset and 32 never looks, so the loss is invisible. Rename is atomic within a
# directory: the audit either exists whole or does not exist.
.tmp <- paste0(OUT, ".tmp")
fwrite(aud, .tmp)
if (!file.rename(.tmp, OUT)) stop(sprintf("30: could not move %s into place as %s", .tmp, OUT))

cat(sprintf("\n30: candidate-edge audit frozen -> %s (%d candidates)\n", OUT, nrow(aud)))
cat("\n== classification ==\n"); print(aud[, .N, by = classification][order(-N)])
cat("\n== reason ==\n");         print(aud[, .N, by = reason][order(-N)])
cat(sprintf("\n== DASHED EDGES (n=%d) ==\n", aud[classification == "dashed_edge", .N]))
print(aud[classification == "dashed_edge",
          .(gene, outcome_name, mr_direction, PP.H4_abf, PP.H4_susie, PP.H3_susie, is_bonferroni)][order(-PP.H4_susie)])
