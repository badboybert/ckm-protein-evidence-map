# ============================================================================
# 31 (Stage 0 · sensitivity): Editor B1 robustness re-run of the 17 dashed edges.
#
# Re-runs coloc.susie for each dashed edge through the FROZEN engine (15/26_susie_coloc*.R) with ONLY two
# susie_rss settings PATCHED IN MEMORY at load time (A": no on-disk copy; single source of truth = the frozen
# engine; a guard asserts EXACTLY 2 patch sites, so any source drift fails loud instead of mis-patching):
#   min_abs_corr 0 -> 0.5  (SuSiE credible-set PURITY filter ON)
#   check_LD    FALSE->TRUE (INERT -- see below; kept only to match the frozen engine's exact text)
#
# `check_LD` IS NOT A susieR PARAMETER. Verified against susieR 0.14.2 (the version that produced these
# results) 2026-07-16: it is absent from formals(susie_rss); no function in the susieR namespace mentions it;
# susie_rss(check_LD=FALSE) and susie_rss(check_LD=TRUE) return identical objects; and so does
# susie_rss(banana_kumquat_9000=TRUE) -- it is swallowed by `...` on its way to susie_suff_stat, which also
# has `...`. It has always been a no-op. It stays in the fixed-string patch below ONLY because the frozen
# engines contain it verbatim and the A" assertion matches on the exact string; removing it would force a
# re-run of frozen engines for no numeric gain.
# SO: this script changes exactly ONE thing that susieR acts on -- min_abs_corr 0 -> 0.5, a POST-HOC
# credible-set purity filter. susie_suff_stat consumes it strictly AFTER the IBSS loop and after the
# non-convergence warning, so CONVERGENCE IS INVARIANT to everything patched here: an edge that does not
# converge in the strict run did not converge in the primary run either. Confirmed on the frozen data --
# all 9 survivors have strict_PP.H4_susie - orig_PP.H4_susie EXACTLY 0 (same fit; only the credible-set
# extraction differs), while non-survivors collapse.
# This script therefore tests PURITY. It does NOT test external-LD consistency; nothing here does. The
# diagnostic for that would be susieR::kriging_rss() (the one the warning's own URL points at), which has
# never been run. See REBUILD_LEDGER 2.5.
#
# ALL warnings from the fit are captured (n_warn / warning_msg). `ld_warn` is NOT "check_LD fired" -- that
# could not happen. It is a regex count of LD-related WORDS over that same text, and on the three held edges
# it matches "LD" inside susieR's remediation advice ("Please check consistency between summary statistics
# and LD matrix"). It means "the message mentions LD". Nothing gates on it; 32 checks non-convergence
# against the convergence text itself (r53 Low, r51 High).
# A dashed edge SURVIVES iff it still reaches strict PP.H4_susie >= 0.80. This tests whether each rescue
# rests on PURE credible sets or on impure ones that the primary run's disabled purity filter let through
# (editor B1). Over-splitting against an external LD panel is the SUSPECTED cause of impurity and is the
# reason the filter matters -- but this script does not test that attribution, and nothing here does.
#
# MULTI-CORE (PSOCK): edges run in parallel; each worker gets an ISOLATED tempdir (PLINK LD scratch,
#   so concurrent plink runs never collide) and each edge a DETERMINISTIC per-edge seed, so results
#   are reproducible regardless of scheduling. MODE dispatch (env CAND_SUSIE_MODE = egfr|nonegfr) — run
#   BOTH. Worker count via env CAND_SENS_CORES (default 4; PSOCK workers each load their own engine +
#   rsid_map, so keep this modest for memory).
#
# NON-DESTRUCTIVE (no CANONICAL write): the only persistent NEW output is candidate_susie_sensitivity_<mode>.csv.
#   Frozen engines (15/26), all canonical results, and the Stage-0 candidate_* files are never written/overwritten.
#   We eval ONLY the (in-memory-patched) frozen engine's function block / head-before-Main-loop (the Main loop is the only canonical
#   writer) — mirrors 29. Side effects are limited to transient scratch (PLINK/python temp files, auto-cleaned)
#   and a harmless dir.create() of the already-existing phase3_susie dir in the eGFR head; no canonical artifact
#   is touched.
#
# Run (both modes):
#   cd "<E>/20260331-CKM/Paper 3"
#   CAND_SUSIE_MODE=nonegfr Rscript scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R
#   CAND_SUSIE_MODE=egfr    Rscript scripts/paper3.1_ckm/phase3/31_candidate_susie_sensitivity.R
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
suppressPackageStartupMessages({ library(data.table); library(parallel) })
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")

MODE <- tolower(Sys.getenv("CAND_SUSIE_MODE", ""))
if (!MODE %in% c("egfr", "nonegfr")) stop("31: set CAND_SUSIE_MODE=egfr or nonegfr (run both, separately)")
NCORES <- suppressWarnings(as.integer(Sys.getenv("CAND_SENS_CORES", "4")))
if (is.na(NCORES) || NCORES < 1L) NCORES <- 1L
cat(sprintf("=== 31 sensitivity (strict SuSiE: credible-set purity ON, min_abs_corr=0.5) | MODE=%s | cores=%d ===\n", MODE, NCORES))

# The MAIN process needs config too. Until 2026-07-16 the only source() in this file was inside
# worker_init(), i.e. workers only -- so the assert_outcome_truth() call below died with "could not find
# function" and the whole chain stopped at 31 (r23 Blocker). It slipped through because the fix was
# verified with parse(), which checks SYNTAX and never resolves a symbol. Sourcing here is idempotent:
# worker_init() still sources it inside each worker, where this process's environment does not reach.
source(file.path(B, "R/00_config.R"))   # assert_outcome_truth(), is_coloc_pass(), OUTCOMES
# CAND_BENIGN_TAGS is used by the write-boundary check at the BOTTOM of this script, in the MAIN process.
# worker_init() sources this file too, but that reaches the workers only -- so the main process would have
# run every edge and then died at the check, discarding the whole run (r44 Blocker). This is r23 exactly:
# same script, same shape -- a guard whose symbol exists only inside the workers. Sourcing is idempotent.
source(file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"))   # CAND_BENIGN_TAGS

RC   <- file.path(B, file.path("results", paste0("repro_", PHASE1_TAG)))
D    <- file.path(RC, "phase3_susie")
AUD  <- file.path(D, "candidate_edge_audit.csv")
POOL <- file.path(D, "candidate_pool.csv")
ORIG <- file.path(D, sprintf("candidate_susie_%s.csv", MODE))
OUT  <- file.path(D, sprintf("candidate_susie_sensitivity_%s.csv", MODE))
# Same reason as 29's, one layer down: this file is written once at the end, and 32 validates its
# schema, keys, row count, outcome and engine but not which run produced it -- so a survivor of a
# killed rerun changes the DRAWN/HELD/DROPPED verdict invisibly. Removed before any validation so
# the window is closed, and the removal is checked because unlink() only returns its failure.
unlink(OUT); if (file.exists(OUT)) stop(sprintf("31: could not remove the previous %s -- refusing to run with a stale result in place", OUT))
# SINGLE SOURCE for this script's output COLUMN ORDER, and for the header of the zero-edge early write --
# used by both, so they cannot drift apart (a hand-copied empty header is exactly how 32's
# `sens[, ..scols]` would break on a mode that legitimately has no dashed edges; r12 Medium, 2026-07-16).
# IT IS NOT A PROJECTION AND IT ENFORCES NOTHING. `setcolorder()` orders columns; it does not subset, and
# unlike 29 this script has no back-fill loop. What actually guarantees the columns is that every exit of
# run_edge() emits the same set -- ok, post-fit collapse, error and pre-fit NULL alike. Naming a column
# here does not create it; dropping one from run_edge would just silently shorten the file.
SENS_COLS <- c("gene", "outcome_name", "susie_engine", "outcome_id", "PP.H4_abf",
  "orig_PP.H4_susie", "orig_n_exp", "orig_n_out", "strict_PP.H4_susie", "strict_PP.H3_susie",
  "strict_n_exp", "strict_n_out", "strict_n_pairs", "strict_n_snps_ld",
  # PER-SIDE convergence, straight from the fit. `n_warn` is counted per EDGE, across BOTH susie_rss calls,
  # so it cannot say whether the exposure (UKB-PPP) or the outcome (e.g. CKDGen eGFR) was the one that
  # failed to converge -- and "the three held edges are all eGFR" was read as an eGFR-side finding on
  # exactly that unsupported inference. The engines now return `$converged` per side; keeping only n_warn
  # would throw the answer away at the last step. Those are two different questions and two different
  # papers.
  "strict_converged_exp", "strict_converged_out",
  "n_warn", "ld_warn",
  "run_status", "survives", "collapse_reason", "warning_msg", "error_msg")
ENG <- if (MODE == "egfr") file.path(B, "scripts/paper3.1_ckm/phase3/26_susie_coloc_egfr.R") else
                           file.path(B, "scripts/paper3.1_ckm/phase3/15_susie_coloc.R")
for (.f in c(AUD, POOL, ORIG, ENG))
  if (!file.exists(.f)) stop(sprintf("31: missing input %s", .f))

# ── the dashed edges for this mode + pool coords + original (loose-setting) PP.H4_susie ──
.aud_all <- fread(AUD)
assert_outcome_truth(.aud_all, "candidate_edge_audit.csv")   # outcome_id drives engine dispatch (r22 High)
dash <- .aud_all[classification == "dashed_edge",
                   .(gene, outcome_name, aud_outcome_id = outcome_id, aud_engine = susie_engine)]
pool <- fread(POOL)
# The pool is treated as a KEYED SET here -- `.runnable` is built from it and compared to the SuSiE cover
# with fsetdiff(), which is a SET operation and therefore blind to duplicates: a duplicated runnable row
# would pass that comparison invisibly. 29 and 30 both guard this; 31 did not, because the r49 fix was
# applied where the review pointed and not to the other reader of the same file (r50 Medium).
# Check the key wherever the key is used.
if (anyDuplicated(pool, by = c("gene", "outcome_name")))
  stop(sprintf("31: candidate_pool.csv has %d duplicate (gene, outcome_name) key(s) -- this script treats the pool as a keyed set (the runnable cover, the audit join, the ORIG join), and set operations cannot see a duplicate",
               sum(duplicated(pool, by = c("gene", "outcome_name")))))
# This is the boundary that actually feeds the engine: run_susie_pair() dispatches on the POOL's
# outcome_id, not the audit's. Validating the audit alone left the live path unguarded (r23 Medium).
assert_outcome_truth(pool, "candidate_pool.csv")
# VALIDATE BEFORE NORMALISING. The normalise below rewrites any non-{egfr,nonegfr} engine to NA -- which
# is right for the ""-from-CSV-round-trip case it was written for, but it also erases a GARBAGE engine
# into the same NA that means "legitimately not testable". Checking after that would be checking the
# repair, not the artifact (r28 Medium). Engine truth first, on what the file actually says.
assert_engine_truth(pool, "candidate_pool.csv")
pool[!susie_engine %in% c("egfr", "nonegfr"), susie_engine := NA_character_]   # NA/CSV round-trip, same as 30
# Every dashed edge went through SuSiE by definition, so an engine-less one is impossible-by-construction
# -- NA is legitimate only for `not_evaluable` (frozen data: dashed 17/0 NA, not_evaluable 17/17). Without
# this, a stale pre-stamp audit row and a stale pool row could both be NA and "agree" (r28). Reject before
# the agreement check, so both-NA can never be read as consensus.
# The audit's dashed set is compared with `fsetdiff()` below and used with `%in%` -- both SET operations,
# both blind to a duplicate row. Assert the key before any of them, not after the joins downstream
# (r51 Medium). Same shape as the pool guard above; same reason.
if (anyDuplicated(dash, by = c("gene", "outcome_name")))
  stop(sprintf("31: the audit's dashed_edge set has %d duplicate (gene, outcome_name) key(s) -- it is about to be compared as a SET, and a set comparison cannot see a duplicate",
               sum(duplicated(dash, by = c("gene", "outcome_name")))))
.d_na <- dash[is.na(aud_engine) | aud_engine == ""]
if (nrow(.d_na))
  stop(sprintf("31: %d dashed audit edge(s) carry no susie_engine: %s -- a dashed edge ran through SuSiE by definition, so this audit is stale or was not produced by 30",
               nrow(.d_na), paste(utils::head(.d_na[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
# THE DASHED SET ITSELF IS DERIVABLE, and until now it was taken on trust from a stored label -- the
# fifth and deepest instance of "the checked thing is not the deciding thing" (r30 High). 30's rule is
# explicit (30:6, 30:138):
#     dashed_edge  <==  run_status == "ok" & coloc_status == "coloc_fail" & PP.H4_susie >= 0.80
# all of which lives in the pool + 29's loose output. Re-derive it and cross-check.
# The direction that matters is the one no downstream check can recover: a stale audit that OMITS a
# genuinely rescued edge produces a sensitivity file and an audit that are perfectly consistent with each
# other, so 32 sees no drift and the edge is simply gone. An extra non-rescued row would at least surface
# as a strict failure; a missing one surfaces as nothing.
# Verified on the frozen data 2026-07-16 before adding this: derived 17, audit 17, identical set, 0 in
# either direction. Costs nothing today.
# CHECK EACH FILE BEFORE THE FILL, on the values not just the presence. A half-stale candidate_susie
# file (one mode re-run under the convergence recording, the other not) would otherwise reach here, get its
# convergence column created by fill=TRUE for the stale mode, and be consumed as if complete. Same guard 32
# uses on the sensitivity files, applied to the loose inputs 31 reads.
for (.lf in file.path(D, sprintf("candidate_susie_%s.csv", c("egfr", "nonegfr"))))
  assert_convergence_honest(fread(.lf), basename(.lf))
.su_all <- rbindlist(lapply(file.path(D, sprintf("candidate_susie_%s.csv", c("egfr", "nonegfr"))), fread),
                     fill = TRUE)
assert_outcome_truth(.su_all, "candidate_susie_{egfr,nonegfr}.csv")
assert_engine_truth(.su_all, "candidate_susie_{egfr,nonegfr}.csv")
# The dashed set is derived from `PP.H4_susie >= 0.80` below, and `Inf >= 0.80` is TRUE -- so integrity
# has to hold BEFORE the derivation, not after it (r39 High).
assert_pp_integrity(.su_all, c("PP.H4_susie", "PP.H3_susie"), "candidate_susie_{egfr,nonegfr}.csv")
assert_status_premise(.su_all, c("PP.H4_susie", "PP.H3_susie"), "candidate_susie_{egfr,nonegfr}.csv")
# COMPLETENESS FIRST. The derivation below inner-joins the pool to the SuSiE output, so an edge missing
# from BOTH the loose file and the audit is missing from `.derived` too -- and the cross-check then passes
# on a consensus of absences. That is the same "both sides absent" hole r28 found in 32, reproduced inside
# the guard added to fix r30 (r31 High). Prove the SuSiE output is a complete, duplicate-free cover of the
# runnable pool BEFORE deriving anything from it. Frozen data: runnable 134, SuSiE 134 unique, 0 missing,
# 0 orphan, 0 duplicate.
pool <- derive_runnability(pool, "candidate_pool.csv")   # facts, not labels (r33)
.runnable <- pool[.exp_testable == TRUE, .(gene, outcome_name)]
.sk <- unique(.su_all[, .(gene, outcome_name)])
if (nrow(.sk) != nrow(.su_all))
  stop(sprintf("31: candidate_susie_{egfr,nonegfr}.csv carry %d duplicate (gene,outcome_name) key(s) -- ambiguous loose provenance",
               nrow(.su_all) - nrow(.sk)))
.miss_su <- fsetdiff(.runnable, .sk); .orph_su <- fsetdiff(.sk, .runnable)
if (nrow(.miss_su) || nrow(.orph_su))
  stop(sprintf("31: the loose SuSiE output is not a complete cover of the runnable pool. Runnable but never tested (%d): %s | tested but not runnable (%d): %s -- an edge absent from BOTH the SuSiE file and the audit would otherwise be invisible to every check below.",
               nrow(.miss_su), paste(utils::head(.miss_su[, paste0(gene, "~", outcome_name)], 8), collapse = ", "),
               nrow(.orph_su), paste(utils::head(.orph_su[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))

# `coloc_status` is a STORED LABEL, and the dashed rule's coloc-fail arm is derivable from the raw ABF
# probabilities via is_coloc_fail() -- evaluated AND not passing; config's single source for that rule,
# and NOT the negation of is_coloc_pass(), which reads an NA-ABF row as failed (r37/r47). Deriving the dashed set from
# the label was trusting a label inside the fix for trusting a label (r31 High).
# ASSERT ON THE FULL POOL, BEFORE THE JOIN. A first version asserted after the inner join to .su_all, so
# it only ever covered the 134 runnable rows with loose output while the comment claimed "151/151 agree"
# -- the check was narrower than the claim, and a stale coloc-passing NON-runnable row would still have
# carried its label into the audit's not_evaluable/excluded provenance (r32 High).
# Frozen data 2026-07-16: 151/151 agree, and 0 of the 17 dashed actually pass.
assert_coloc_status_truth(pool, "candidate_pool.csv")   # full pool, before any join (r32/r33)
# is_coloc_fail(), not the negation of the pass rule: 31 is explicitly a stale-artifact guard and must
# not assume 30 was current, so it needs its own evaluability floor (r38 High).
pool[, .coloc_fail_rd := is_coloc_fail(PP.H4_abf, PP.H3_abf)]
.dj <- merge(pool[, .(gene, outcome_name, .coloc_fail_rd)],
             .su_all[, .(gene, outcome_name, run_status, PP.H4_susie)], by = c("gene", "outcome_name"))
.derived <- .dj[run_status == "ok" & .coloc_fail_rd == TRUE &
                !is.na(PP.H4_susie) & PP.H4_susie >= 0.80, .(gene, outcome_name)]
.only_der <- fsetdiff(.derived, dash[, .(gene, outcome_name)])
.only_aud <- fsetdiff(dash[, .(gene, outcome_name)], .derived)
if (nrow(.only_der) || nrow(.only_aud))
  stop(sprintf("31: the audit's dashed_edge set disagrees with 30's own rule re-derived from the pool + 29's output (run_status=='ok' & coloc_status=='coloc_fail' & PP.H4_susie>=0.80). Missing from the audit (%d): %s | in the audit but not rescued (%d): %s -- the audit is stale. A missing edge cannot be recovered downstream: 32 would see a self-consistent audit/sensitivity pair and never know. Refusing.",
               nrow(.only_der), paste(utils::head(.only_der[, paste0(gene, "~", outcome_name)], 8), collapse = ", "),
               nrow(.only_aud), paste(utils::head(.only_aud[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))

edges_all <- merge(dash, pool, by = c("gene", "outcome_name"))
.p_na <- edges_all[is.na(susie_engine) | susie_engine == ""]
if (nrow(.p_na))
  stop(sprintf("31: %d dashed edge(s) have no susie_engine in the POOL: %s -- they would vanish at the mode filter below without ever being tested",
               nrow(.p_na), paste(utils::head(.p_na[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
# no dashed edge may silently vanish or reroute on audit/pool drift (Codex review Low): every audit dashed edge
# must have a pool row, and its outcome_id / susie_engine must AGREE between the two artifacts.
if (nrow(edges_all) != nrow(dash))
  stop(sprintf("31: %d of %d dashed audit edge(s) have NO pool row -- audit/pool drift, refusing", nrow(dash) - nrow(edges_all), nrow(dash)))
# A one-sided NA IS a mismatch. The `!is.na(x) & !is.na(y) & x != y` shape treats "audit has an id, pool
# has NA" as agreement, which is exactly the drift this guard exists to catch (r23 Medium).
# NOTE the parenthetical that used to sit here -- "both NA is agreement, one NA is not" -- was true of the
# pool at large and FALSE of the rows reaching this line, which are dashed edges only. It taught the very
# bug r28 found in 32. Both-NA cannot occur here: the precheck above rejects an engine-less dashed row on
# either side. `.ne` therefore only ever fires on a genuine non-missing disagreement.
.ne <- function(x, y) xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
mm <- edges_all[.ne(aud_outcome_id, outcome_id) | .ne(aud_engine, susie_engine)]
if (nrow(mm)) stop(sprintf("31: audit vs pool outcome_id/susie_engine mismatch for %s -- provenance drift", paste(mm[, paste(gene, outcome_name)], collapse = ", ")))
edges <- edges_all[susie_engine == MODE]
# post-filter completeness (Codex review Low): edges selected for this mode MUST equal what the AUDIT lists for
# it. Closes the gap where a pool engine blank/normalized-to-NA passes the non-NA-only mismatch guard above and
# then silently drops the edge at this filter.
n_mode_aud <- nrow(dash[aud_engine == MODE])
if (nrow(edges) != n_mode_aud)
  stop(sprintf("31: mode %s selected %d edge(s) but the audit lists %d dashed edge(s) for it -- an edge dropped at the engine filter (blank/mismatched pool engine?); refusing", MODE, nrow(edges), n_mode_aud))
edges[, c("aud_outcome_id", "aud_engine") := NULL]
# ORIG is joined on the NAME only, so a stale row -- or one from the other engine -- carrying the same
# pair name would be copied in as this edge's loose-setting provenance (r26 Medium). Validate it like
# every other read boundary, and require it to be THIS mode's file, before taking any value from it.
assert_convergence_honest(fread(ORIG), basename(ORIG))   # ORIG is loose provenance copied into the output; must not be half-stale
.orig_all <- fread(ORIG)
assert_outcome_truth(.orig_all, sprintf("candidate_susie_%s.csv", MODE))
assert_engine_truth(.orig_all, sprintf("candidate_susie_%s.csv", MODE))
# NA is NOT acceptable here, unlike in the pool: 29 stamps EVERY row it emits with its MODE, so an
# engine-less row in this file can only mean a stale pre-stamp version -- and it would still feed
# orig_PP.H4_susie through the name-only join below (r27 Medium). A zero-row file has no row to excuse.
.wrong_mode <- .orig_all[is.na(susie_engine) | susie_engine != MODE]
if (nrow(.wrong_mode))
  stop(sprintf("31: %s contains %d row(s) whose susie_engine is missing or not '%s' -- wrong/stale mode file; refusing to copy loose-setting provenance from it",
               basename(ORIG), nrow(.wrong_mode), MODE))
# The loose values copied below are this edge's evidence for having been rescued at all, so require them
# to actually say that -- reading orig_PP.H4_susie without checking the run_status and threshold that made
# it a dashed edge is trusting the label again, one level down (r30 High).
.bad_orig <- .orig_all[paste(gene, outcome_name) %in% dash[, paste(gene, outcome_name)] &
                       (run_status != "ok" | is.na(PP.H4_susie) | PP.H4_susie < 0.80)]
if (nrow(.bad_orig))
  stop(sprintf("31: %d dashed edge(s) whose loose SuSiE row does not meet the rescue rule it was classified by (run_status=='ok' & PP.H4_susie>=0.80): %s -- the audit and 29's output disagree; refusing to re-run them under strict settings as though they had been rescued",
               nrow(.bad_orig), paste(utils::head(.bad_orig[, sprintf("%s~%s(status=%s, PP.H4=%s)", gene, outcome_name,
                                                                      run_status, format(PP.H4_susie))], 8), collapse = ", ")))
orig  <- .orig_all[, .(gene, outcome_name, orig_PP.H4_susie = PP.H4_susie,
                         orig_n_exp = n_signals_exp, orig_n_out = n_signals_out)]
edges <- merge(edges, orig, by = c("gene", "outcome_name"), all.x = TRUE)
if (nrow(edges[is.na(orig_PP.H4_susie)])) stop("31: a dashed edge lacks its original PP.H4_susie -- incomplete join to the 29 output")
if (anyDuplicated(edges, by = c("gene", "outcome_name"))) stop("31: dashed edges not unique after joins")
cat(sprintf("31 [%s]: %d dashed edge(s) (of %d total dashed) to re-run with strict settings\n", MODE, nrow(edges), nrow(dash)))
# Write an EMPTY but schema-correct CSV instead of nothing: 32 requires BOTH mode files (its input
# existence check stops on a missing one), so quitting silently turns "this mode legitimately has no
# dashed edges" into a broken chain (r12 Medium). Mirrors 29's fixed-schema empty write.
if (nrow(edges) == 0L) {
  cat("31: no dashed edges for this mode -- writing an empty, schema-correct CSV so 32 can still run\n")
  fwrite(as.data.table(setNames(rep(list(character(0)), length(SENS_COLS)), SENS_COLS)), OUT)
  cat(sprintf("31: wrote %s (0 rows)\n", OUT))
  quit(save = "no")
}

# ── worker init: replicate 29's engine environment; PATCH the frozen engine in memory (A"): read 15/26,
#    flip ONLY the two susie_rss settings via a guarded gsub (assert EXACTLY 2 sites), then eval its function
#    block / head-before-Main-loop. No on-disk strict copy; single source of truth = the frozen 15/26. ──
worker_init <- function(B, MODE, ENG) {
  Sys.setenv(CKM_OFFLINE = "1")
  suppressPackageStartupMessages({ library(data.table); library(coloc); library(ieugwasr); library(susieR) })
  # DEFINED INSIDE worker_init, and it has to be. clusterCall() serialises this function to PSOCK workers --
  # separate R processes. The rule is NOT "the worker has nothing the main process defined" (this script
  # exports `B`, `MODE`, `ENG`, `worker_init` and `run_edge` below); it is: **nothing exists in a worker
  # unless it was exported, passed as an argument, sourced, library()'d, or created there.** A helper
  # declared at top level in this file is none of those, so it is a free variable in the worker and fails
  # with "could not find function".
  # Learned by running: the first version of this guard sat above worker_init, every worker died on it, and
  # four static review rounds had passed it -- parse() cannot see symbol resolution across a process
  # boundary. Same shape as the two earlier findings in this file (the main process not sourcing 00_config;
  # CAND_BENIGN_TAGS living only in worker_init).
  #
  # Asserts the strict settings reached the code that will actually be EVALUATED -- see the scope note at
  # the patch site below. Requires BOTH: exactly 2 strict sites present, and ZERO loose sites left behind.
  .assert_patched <- function(block, eng, what) {
    n_strict <- sum(grepl("min_abs_corr = 0.5, verbose = FALSE, check_LD = TRUE", block, fixed = TRUE))
    n_loose  <- sum(grepl("min_abs_corr = 0, verbose = FALSE, check_LD = FALSE", block, fixed = TRUE))
    if (n_strict != 2L || n_loose != 0L)
      stop(sprintf("31 worker: the %s of %s that is about to be evaluated carries %d strict susie_rss site(s) and %d un-patched loose site(s) -- expected 2 and 0. The patch did not reach the code that runs; this would silently produce a LOOSE fit and report it as the STRICT sensitivity.",
                   what, eng, n_strict, n_loose))
    invisible(TRUE)
  }
  source(file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"))   # classify_null(): benign collapse vs coloc.susie error
  # PLINK LD scratch is isolated BY CONSTRUCTION: ieugwasr:::ld_matrix_local() writes via tempfile() into this
  # worker's own tempdir(); each PSOCK worker is a separate R process with a unique tempdir(), and tempfile()
  # yields a fresh name per call -> no cross-worker or intra-worker collision (verified 2026-07-12). No TMPDIR
  # override is needed (and it would be a no-op: R fixes tempdir() at startup).
  # THE GENERATION IS RESOLVED HERE, IN THE WORKER, BEFORE IT IS USED. The next line builds RC from
  # PHASE1_TAG and reads PHASE1_INSTRUMENTS, but nothing in a PSOCK worker defines either: the config
  # source further down sits inside a branch, and only B / MODE / ENG are exported. Every worker died
  # with "object 'PHASE1_TAG' not found" -- the same shape as the note above about helpers declared
  # outside worker_init, and the same shape as 30's ordering bug found in this run.
  source(file.path(B, "R/_phase1_set.R"))
  source(file.path(B, "R/00_config.R"))
  g  <- globalenv(); RC <- file.path(B, file.path("results", paste0("repro_", PHASE1_TAG)))
  assign("SUSIE_COLOC_INPUT", file.path(RC, "coloc_results_v7.rds"), envir = g)
  assign("SUSIE_INSTR",       PHASE1_INSTRUMENTS, envir = g)
  assign("SUSIE_OUT_DIR",     file.path(RC, "phase3_susie"), envir = g)
  assign("SUSIE_EGFR_COLOC",  get("SUSIE_COLOC_INPUT", g), envir = g)
  assign("SUSIE_EGFR_INSTR",  get("SUSIE_INSTR", g), envir = g)
  assign("SUSIE_EGFR_OUT",    get("SUSIE_OUT_DIR", g), envir = g)
  # A"): load the FROZEN engine and patch ONLY the two susie_rss settings IN MEMORY (no on-disk copy).
  src  <- readLines(ENG)
  .pat <- "min_abs_corr = 0, verbose = FALSE, check_LD = FALSE"
  .rep <- "min_abs_corr = 0.5, verbose = FALSE, check_LD = TRUE"
  .n   <- sum(grepl(.pat, src, fixed = TRUE))
  if (.n != 2L) stop(sprintf("31 worker: expected EXACTLY 2 susie_rss setting sites to patch in %s, found %d -- engine layout changed; refusing to run a mis-patched sensitivity", basename(ENG), .n))
  src  <- gsub(.pat, .rep, src, fixed = TRUE)
  # The count above is taken over the WHOLE engine file, but only a SLICE of it is ever evaluated (below:
  # the function block, not the main loop). Those are different scopes, and the guard was checking the one
  # that does not run: an engine reorganised so that a susie_rss call inside the evaluated slice no longer
  # matched `.pat` -- while two matches survived elsewhere in the file -- would pass `.n == 2L`, run the
  # LOOSE settings, and report the result as the STRICT sensitivity. The whole 17 -> 6/3/8 split rests on
  # this patch actually reaching the code that runs, so assert it on the code that runs. Both sites sit
  # inside the evaluated slice in both engines -- established by the assertion itself, at run time. An
  # earlier version of this comment quoted the line numbers (331/338, slices [113,384] and [1,379]); a same
  # -day edit to the engines moved every one of them. That is exactly why the assertion greps for the text
  # instead of counting lines, and why this comment no longer names any. Asserted AFTER the slice is taken.
  if (MODE == "nonegfr") {
    source(file.path(B, "R/00_config.R"))
    assign("COLOC_SOURCE_ONLY", TRUE, envir = g)
    suppressMessages(source(file.path(B, "R/coloc.R")))
    if (exists("COLOC_SOURCE_ONLY", envir = g, inherits = FALSE)) rm("COLOC_SOURCE_ONLY", envir = g)
    a <- grep("Helper: LD matrix", src)[1]; z <- grep("Main loop", src)[1] - 1L
    if (is.na(a) || is.na(z) || a > z) stop("31 worker: cannot locate 15's function block")
    .assert_patched(src[a:z], basename(ENG), "15 function block")
    eval(parse(text = paste(src[a:z], collapse = "\n")), envir = g)
    assign("rsid_map", load_rsid_map(), envir = g)
  } else {
    cut <- grep("Main loop", src)[1]
    if (is.na(cut)) stop("31 worker: no 'Main loop' marker in 26")
    .assert_patched(src[seq_len(cut - 1)], basename(ENG), "26 pre-main-loop block")
    eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")), envir = g)
  }
  if (!(exists("PLINK_BIN", envir = g) && file.exists(get("PLINK_BIN", g)))) stop("31 worker: PLINK_BIN missing/not found")
  TRUE
}

# ── per-edge run: deterministic seed; capture stdout/messages/ALL warnings; NULL -> classify (benign collapse
#    vs coloc.susie error), throw -> an error ROW. An ENGINE RUNTIME failure never aborts the batch.
#    An INVARIANT failure does: the usability assertion below stop()s, because a non-NULL result missing
#    PP.H4_susie is a schema/computation failure that must not be recorded as "this edge just didn't reach
#    0.80" (r15). "Never aborts the batch" was written before that stop existed and was left standing --
#    a comment teaching the code's old behaviour (r30 Low). ──
run_edge <- function(row, MODE) {
  set.seed(1L + sum(utf8ToInt(paste0(row$gene, "~", row$outcome_name))))   # deterministic per edge -> reproducible
  warn_msgs <- character(0); msg_lines <- character(0); box <- new.env()
  # capture the engine's stdout (cat, incl. "[skip] ..." reasons), messages (incl. "[coloc.susie error]"), and
  # ALL warnings as TEXT (not just an LD-regex count -- Codex review Medium: SuSiE convergence / numerical
  # warnings must not be silently muffled, or a survivor could hide a non-convergence warning).
  out_lines <- capture.output({
    box$res <- withCallingHandlers(
      tryCatch(
        if (MODE == "egfr")
          run_susie_egfr_pair(row$gene, row$gene_chrom, row$gene_start, row$gene_end, pp4_abf = row$PP.H4_abf)
        else
          run_susie_pair(row$gene, row$outcome_id, row$gene_chrom, row$gene_start, row$gene_end, rsid_map, pp4_abf = row$PP.H4_abf),
        error = function(e) structure(list(msg = conditionMessage(e)), class = "edge_error")),
      warning = function(w) { warn_msgs[[length(warn_msgs) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") },
      message = function(m) { msg_lines[[length(msg_lines) + 1L]] <<- sub("\n$", "", conditionMessage(m)); invokeRestart("muffleMessage") })
  }, type = "output")
  res      <- box$res
  ld_warn  <- sum(grepl("LD|correlation|inconsist|kriging|check_ld|refin", warn_msgs, ignore.case = TRUE))
  warn_txt <- if (length(warn_msgs)) paste(unique(trimws(warn_msgs)), collapse = " || ") else NA_character_
  base <- data.table(gene = row$gene, outcome_name = row$outcome_name, susie_engine = MODE, outcome_id = row$outcome_id,
                     PP.H4_abf = row$PP.H4_abf, orig_PP.H4_susie = row$orig_PP.H4_susie,
                     orig_n_exp = row$orig_n_exp, orig_n_out = row$orig_n_out,
                     n_warn = length(warn_msgs), ld_warn = ld_warn, warning_msg = warn_txt)
  if (inherits(res, "edge_error"))
    return(cbind(base, data.table(strict_PP.H4_susie = NA_real_, run_status = "error",
                                  # NA: the engine threw, and the throw could be anywhere -- whether the
                                  # fits ran at all is unknown, so their convergence is unknown too.
                                  strict_converged_exp = NA, strict_converged_out = NA,
                                  collapse_reason = NA_character_, error_msg = res$msg)))
  if (is.null(res)) {
    # distinguish a BENIGN strict collapse (purity filtered out the credible sets / empty summary / too few SNPs
    # = the informative "did not survive") from a genuine coloc.susie computation ERROR (must NOT be counted as a
    # collapse, or the survivor count is biased). Reuse the Stage-0 classifier (Codex review Medium):
    #   benign -> run_status="collapsed"; fatal / coloc.susie-error / unrecognised -> run_status="error".
    cl <- classify_null(c(out_lines, msg_lines))
    st <- if (identical(cl$status, "benign")) "collapsed" else "error"
    return(cbind(base, data.table(strict_PP.H4_susie = NA_real_, run_status = st,
                                  # NA: since 2026-07-16 a NULL from the engine means it returned BEFORE
                                  # fitting (no pQTL data, <20 SNPs, no LD). No fit ran; convergence is
                                  # genuinely unknown. A post-fit collapse no longer arrives as NULL -- see
                                  # the branch immediately below.
                                  strict_converged_exp = NA, strict_converged_out = NA,
                                  collapse_reason = cl$tag, error_msg = if (st == "error") cl$detail else NA_character_)))
  }
  r <- as.data.table(res)

  # A POST-FIT COLLAPSE NOW ARRIVES AS A ROW, NOT AS NULL -- and this branch has to exist because of that.
  # 15/26 were changed on 2026-07-16 so that a coloc.susie failure or empty summary AFTER both susie_rss
  # calls succeed still returns the fits' convergence instead of collapsing to a bare NULL (their own
  # no-fit fallback has no convergence columns, so a fit whose convergence was KNOWN was being written as
  # if it were not). That change is correct for susie_coloc_results -- and it silently broke the contract
  # THIS function relies on: the assertion below was written when "non-NULL" implied "usable", and would
  # now abort the entire run on every benignly-collapsed edge. That is not hypothetical: 4 of the 17
  # dashed candidates are `dropped:collapsed` / `empty_coloc_summary`.
  # So: detect the shape, classify it exactly as the NULL path does, and keep EVERYTHING the row now
  # carries. The row has MORE information than before, not less; the only bug would be to throw it away --
  # which is exactly what the first version of this branch did, three lines under that sentence. It kept
  # the convergence and rebuilt the rest from `base`, so the credible-set and SNP counts the fits had
  # already produced were written as NA. They are facts on this path: the fits ran, and only coloc.susie
  # failed after them. Blanking them also defeats susie_fit_ran() (R/00_config.R), which reads those very
  # columns as the evidence a fit ran.
  if ("PP.H4_susie" %in% names(r) && nrow(r) == 1L && is.na(r$PP.H4_susie)) {
    cl <- classify_null(c(out_lines, msg_lines))
    st <- if (identical(cl$status, "benign")) "collapsed" else "error"
    .known <- data.table(
      strict_PP.H4_susie = NA_real_,
      strict_PP.H3_susie = if ("PP.H3_susie" %in% names(r)) r$PP.H3_susie else NA_real_,
      strict_n_exp       = if ("n_signals_exp" %in% names(r)) r$n_signals_exp else NA_integer_,
      strict_n_out       = if ("n_signals_out" %in% names(r)) r$n_signals_out else NA_integer_,
      strict_n_pairs     = if ("n_signal_pairs" %in% names(r)) r$n_signal_pairs else NA_integer_,
      strict_n_snps_ld   = if ("n_snps_ld" %in% names(r)) r$n_snps_ld else NA_integer_,
      strict_converged_exp = isTRUE(r$converged_exp), strict_converged_out = isTRUE(r$converged_out),
      run_status = st, collapse_reason = cl$tag,
      error_msg = if (st == "error") cl$detail else NA_character_)
    return(cbind(base, .known))
  }
  # A non-NULL return is NOT proof of a usable result (r15 High, 2026-07-16). `skip_absent = TRUE` below is
  # the silent path: if the engine stops returning PP.H4_susie, setnames quietly skips it, `keep` picks up
  # fewer columns, this row is still stamped run_status=="ok", the fixed-schema block near the write
  # back-fills strict_PP.H4_susie as NA -- and 32 then reads that NA as `dropped:strict_fail`, i.e. a
  # COMPUTATION failure presented as "this edge simply did not reach 0.80". Assert BEFORE the rename, while
  # the engine's own column names are still intact. Measured 2026-07-16 on the frozen output: 13 ok rows,
  # 0 with a missing or NA strict_PP.H4_susie -- free today.
  .need <- c("PP.H4_susie", "PP.H3_susie")
  .miss <- setdiff(.need, names(r))
  # NA over EVERY required column, not just the first -- same fix as 29 (r16 High).
  # RANGE, not just presence. `Inf >= 0.80` and `1.2 >= 0.80` are both TRUE, so a malformed PP would be
  # judged to have survived the strict threshold (r38 High). A probability outside [0,1] is a computation
  # failure, not a high one.
  # Same rule as everywhere else -- see assert_pp_integrity() in R/00_config.R. Inlined here only because
  # this runs inside a parallel worker, where config is sourced but the row is a single engine return
  # rather than a table with a run_status column.
  .nacol <- if (length(.miss)) character(0) else
            .need[vapply(.need, function(.c) any(!is.finite(r[[.c]]) | r[[.c]] < 0 | r[[.c]] > 1), TRUE)]
  if (length(.miss) || length(.nacol))
    stop(sprintf("31 worker [%s ~ %s]: engine returned a non-NULL result that is NOT usable -- %s. This is a computation/schema failure, not a benign collapse, and must never be stamped run_status=='ok'.",
                 row$gene, row$outcome_name,
                 if (length(.miss)) paste("missing column(s):", paste(.miss, collapse = ", "))
                 else paste("missing, non-finite or out-of-range value in required column(s):",
                        paste(.nacol, collapse = ", "))))
  setnames(r, c("PP.H4_susie", "PP.H3_susie", "n_signals_exp", "n_signals_out", "n_snps_ld", "n_signal_pairs",
                "converged_exp", "converged_out"),
              c("strict_PP.H4_susie", "strict_PP.H3_susie", "strict_n_exp", "strict_n_out", "strict_n_snps_ld", "strict_n_pairs",
                "strict_converged_exp", "strict_converged_out"),
           skip_absent = TRUE)
  keep <- intersect(c("strict_PP.H4_susie", "strict_PP.H3_susie", "strict_n_exp", "strict_n_out", "strict_n_snps_ld", "strict_n_pairs",
                      "strict_converged_exp", "strict_converged_out"), names(r))
  cbind(base, r[, ..keep], data.table(run_status = "ok", collapse_reason = NA_character_, error_msg = NA_character_))
}

# ── run (PSOCK parallel; falls back to serial when 1 core / 1 edge) ──
edge_list <- split(edges, seq_len(nrow(edges)))
nw <- max(1L, min(NCORES, nrow(edges)))
cat(sprintf("31 [%s]: launching %d worker(s)\n", MODE, nw))
if (nw > 1L) {
  cl <- makeCluster(nw, type = "PSOCK")
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  clusterExport(cl, c("B", "MODE", "ENG", "worker_init", "run_edge"), envir = environment())
  # Propagate the palindrome-capture env into each PSOCK worker EXPLICITLY: workers do not inherit the master's
  # options(ckm.palindrome_sink), and OS env inheritance across makePSOCKcluster is platform-dependent, so
  # strand_align_R (which reads getOption then Sys.getenv) would otherwise capture nothing from 31's parallel fits
  # (cold re-review Medium). No-op when capture is off (empty env).
  .psink <- Sys.getenv("CKM_PALINDROME_SINK", ""); .prid <- Sys.getenv("CKM_RUN_ID", "")
  clusterExport(cl, c(".psink", ".prid"), envir = environment())
  invisible(clusterEvalQ(cl, if (nzchar(.psink)) Sys.setenv(CKM_PALINDROME_SINK = .psink, CKM_RUN_ID = .prid)))
  invisible(clusterCall(cl, function(B, MODE, ENG) worker_init(B, MODE, ENG), B, MODE, ENG))
  results <- parLapplyLB(cl, edge_list, function(r) run_edge(r, MODE))
} else {
  worker_init(B, MODE, ENG)
  results <- lapply(edge_list, function(r) run_edge(r, MODE))
}

final <- rbindlist(results, fill = TRUE, use.names = TRUE)
# ensure the full schema exists even if EVERY edge in a mode collapsed/errored, so the summary print (which
# references strict_* cols) cannot fail after the CSV is written (Codex High).
for (col in c("strict_PP.H4_susie", "strict_PP.H3_susie", "strict_n_exp", "strict_n_out", "strict_n_snps_ld", "strict_n_pairs", "ld_warn", "n_warn"))
  if (!col %in% names(final)) final[, (col) := NA_real_]
for (col in c("collapse_reason", "error_msg", "warning_msg"))
  if (!col %in% names(final)) final[, (col) := NA_character_]
final[, survives := run_status == "ok" & !is.na(strict_PP.H4_susie) & strict_PP.H4_susie >= 0.80]
setcolorder(final, intersect(SENS_COLS, names(final)))   # SENS_COLS: single source, defined near OUT
assert_pp_integrity(final, c("strict_PP.H4_susie", "strict_PP.H3_susie"),
                    sprintf("candidate_susie_sensitivity_%s.csv (about to be written)", MODE))
assert_status_premise(final, c("strict_PP.H4_susie", "strict_PP.H3_susie"),
                      sprintf("candidate_susie_sensitivity_%s.csv (about to be written)", MODE))
# A `collapsed` row asserts a BENIGN collapse -- CAND_BENIGN_TAGS is the single source for which those
# are. Checked here at the write boundary as well as in 32, because a producer that emits an
# unrecognised collapse reason has already lost the distinction 32 would be trying to recover (r43).
.cr_bad <- final[run_status == "collapsed" & (is.na(collapse_reason) | !collapse_reason %in% CAND_BENIGN_TAGS)]
if (nrow(.cr_bad))
  stop(sprintf("31: %d row(s) are run_status=='collapsed' but carry a missing or non-benign collapse_reason (expected one of %s): %s -- refusing to write a computation failure as a benign collapse",
               nrow(.cr_bad), paste(CAND_BENIGN_TAGS, collapse = "/"),
               paste(utils::head(.cr_bad[, sprintf("%s~%s(reason=%s)", gene, outcome_name,
                                                   as.character(collapse_reason))], 8), collapse = ", ")))
fwrite(final, OUT)

cat(sprintf("\n31 [%s]: %d edge(s) | SURVIVED(strict PP.H4>=0.80)=%d | collapsed(benign)=%d | ERROR=%d | with-warnings=%d\n",
    MODE, nrow(final), final[survives == TRUE, .N], final[run_status == "collapsed", .N],
    final[run_status == "error", .N], final[!is.na(n_warn) & n_warn > 0, .N]))
if (final[run_status == "error", .N] > 0)
  cat(sprintf("31: [!] %d edge(s) hit a coloc.susie ERROR (NOT a benign collapse) -- see run_status=='error' rows / error_msg\n", final[run_status == "error", .N]))
cat(sprintf("31: wrote %s\n", OUT))
print(final[order(-strict_PP.H4_susie),
            .(gene, outcome_name, orig = orig_PP.H4_susie, strict = strict_PP.H4_susie,
              n = paste0(strict_n_exp, "/", strict_n_out), n_warn, run_status, survives, collapse_reason)])
