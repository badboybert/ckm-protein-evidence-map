# ============================================================================
# 32 (Stage 0 · finalize): apply the LOCKED dashed-edge drawing rule (journal-editor final sign-off
# 2026-07-12) to the B1 sensitivity re-run, DERIVING which of the 16 candidate dashed edges are
# DRAWN / HELD / DROPPED (currently 7/1/8 — see the supersede note below). Turns the originally hand-applied
# draw rule into a reproducible, asserted program output.
#
# ⚠️ 2026-07-18 — COUNTS IN THE COMMENTS BELOW ARE SUPERSEDED (retained as derivation rationale). The
# 2026-07-12 sign-off (17 candidates -> draw 6 / 3 held / 8 dropped, verified at SuSiE max_iter=100) was
# reconciled after a uniform max_iter=1000 rebuild: DPEP1-eGFR LEFT the candidate set (converged
# PP.H4_SuSiE=0.08 < 0.80), ATP5IF1-eGFR became DRAWN, ZFYVE19-eGFR stays HELD (does not converge even at
# 1000 iterations). CURRENT frozen output = 16 candidates -> 7 DRAWN / 1 HELD / 8 DROPPED (8 strict
# survivors), in dashed_edges_final.csv. The inline "17 / 6-3-8 / 3 held / 9 survivors / 100 iterations"
# figures below describe the pre-rebuild state; the CODE is data-driven and produces the current 16/7/1/8
# from the max_iter=1000 sensitivity CSVs (candidate_susie_sensitivity_{egfr,nonegfr}.csv = 6 + 10 = 16).
#
# RULE (an edge is DRAWN iff ALL three hold):
#   (1) MR FDR q<0.05                           -- the manuscript's PRIMARY rule: BH once over the full
#       N_PROTEINS x length(OUTCOMES) grid with untested cells at p=1, via primary_fdr_set()
#       (R/00_config.R). 28 builds the pool from the same call; this RE-ASSERTS it independently.
#       Until 2026-07-15 both read mr_significant_fdr.rds instead -- the saved by-method artifact (three-
#       surface block in R/00_config.R), so this rule was enforced under a definition it did not state. Comparing the two definitions on the frozen data PREDICTS that all
#       17 candidate dashed edges are FDR-significant under both, so this gate's verdict should not move
#       and only the pool should (151 -> 142, the 9 leavers already `excluded`/`not_evaluable`). That is
#       what the gated re-run must CONFIRM -- this script still reads the pre-fix audit until then.
#   (2) strict PP.H4_susie >= 0.80              -- RE-DERIVED here from run_status + strict_PP.H4_susie
#       (SuSiE purity ON, min_abs_corr=0.5), then cross-checked against 31's `survives`. Reading that
#       column directly would classify the rule by its own output -- a tautology (r20).
#   (3) ZERO warnings of ANY kind             -- 31's n_warn == 0. NOTE the rule is deliberately
#       ALL-warning, not LD/convergence-specific: 31 captures "ALL warnings as TEXT (not just an LD-regex
#       count -- ... SuSiE convergence / numerical warnings must not be silently muffled, or a survivor
#       could hide a non-convergence warning)". So an edge carrying ANY warning is held. This is stricter
#       than the LOCKED prose ("LD/convergence warnings"), and the held:non_converged label names the
#       observed cause of all 3 held edges, not the gate's scope. Do not narrow the code to match the
#       label -- that would re-open exactly the muffling 31 was built to prevent (r13 Low, 2026-07-16).
#   held:non_converged = survives(2) but n_warn>0 (fails 3) -> evaluated-but-not-drawn.
#   dropped:collapsed  = no pure credible set (run_status collapsed).
#   dropped:strict_fail= ran ok but strict PP.H4 < 0.80.
#
# COUNTS ARE DERIVED + PRINTED, NEVER HARD-CODED. Assertions check the RULE + a clean partition + zero
# coloc.susie errors + the FDR gate (an external cross-check) -- NOT "== 6" (which would be a hard-coded
# expectation that breaks on legitimate data change; and asserting the rule against its own definition would
# be tautological, the failure mode the editor flagged for `concordant_susie`).
# As of the 2026-07-12 sign-off this yields DRAWN=6 / HELD=3 / DROPPED=8 (documentation, not an assertion).
#
# NON-DESTRUCTIVE: reads the sensitivity CSVs + candidate_edge_audit.csv + mr_all_results.rds (rule 1 is
#   re-derived from it via primary_fdr_set(); it previously read mr_significant_fdr.rds, the saved
#   by-method artifact, which is not the rule this file states -- 2026-07-15); writes ONE
#   new file dashed_edges_final.csv (the authoritative drawn/held/dropped set for Stage 1). The Stage-0 audit
#   (loose-setting classification of every candidate in the pool) is NOT modified. The pool size is
#   whatever 28 emits -- do not restate it here; narrowing the FDR gate to the primary rule is expected
#   to change it, and a count frozen in a comment is how a stale claim outlives the thing it described.
# Run: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/32_dashed_finalize.R
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
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
B  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
# ROOTED AT THE SELECTED GENERATION, not at a literal. This read and wrote results/repro_corrected
# -- the PREVIOUS generation -- with no PHASE1_TAG anywhere, and 32 is the sole producer of
# dashed_edges_final.csv, which ten L7/L8 scripts read. So it would have silently handed the new
# generation the old generation's dashed layer. The config has to be sourced FIRST, because that is
# what defines PHASE1_TAG and RESULTS_DIR.
source(file.path(B, "R/00_config.R"))   # primary_fdr_set() -- OWNER of the primary FDR rule (see rule (1) above); also PHASE1_TAG / RESULTS_DIR
RC <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))
D  <- file.path(RC, "phase3_susie")
# CAND_BENIGN_TAGS: the single source for "which collapses are benign". 32 turns run_status=="collapsed"
# into `dropped:collapsed` -- a BIOLOGICAL reading -- so it must check that the collapse really was one of
# the benign kinds, not merely that the label says "collapsed" (r43 High).
source(file.path(B, "scripts/paper3.1_ckm/phase3/_candidate_null_classify.R"))
AUD  <- file.path(D, "candidate_edge_audit.csv")
SENS <- file.path(D, c("candidate_susie_sensitivity_nonegfr.csv", "candidate_susie_sensitivity_egfr.csv"))
MR_ALL <- file.path(RC, "mr_all_results.rds")   # rule (1) re-derives the PRIMARY FDR set from this
OUT  <- file.path(D, "dashed_edges_final.csv")
for (.f in c(AUD, SENS, MR_ALL)) if (!file.exists(.f)) stop(sprintf("32: missing input %s", .f))

# ── the 17 candidate dashed edges (loose-setting) + their audit provenance ──
.aud_all <- fread(AUD)
assert_outcome_truth(.aud_all, "candidate_edge_audit.csv")   # one helper, every read boundary (r22)
assert_engine_truth(.aud_all, "candidate_edge_audit.csv")   # engine is a function of the outcome (r27)
dash <- .aud_all[classification == "dashed_edge",
                 .(gene, outcome_name, outcome_id, susie_engine, PP.H4_abf, PP.H3_abf, mr_direction,
                   is_bonferroni)]
if (anyDuplicated(dash, by = c("gene", "outcome_name"))) stop("32: audit dashed edges not unique")
# A dashed edge is BY DEFINITION one that FAILED single-variant coloc -- that is the whole premise of the
# layer, and the caption says so. The audit's classification is a stored label; the raw ABF probabilities
# are right here. Recheck them rather than finalize a coloc-PASSING edge as exploratory-dashed on the
# strength of a stale label (r31 High). Frozen data: 0 of the 17 satisfy is_coloc_pass.
# State the PREMISE positively: every dashed edge must be `is_coloc_fail` -- evaluated, and not passing.
# A pass-only stop leaves the other way to violate the premise wide open: a stale pre-r37 row with NA or
# out-of-range ABF is neither passing nor evaluable, sails through, and can be drawn (r38 High).
.bad_prem <- dash[!is_coloc_fail(PP.H4_abf, PP.H3_abf)]
if (nrow(.bad_prem)) {
  .bad_prem[, .why := fifelse(!is_coloc_evaluable(PP.H4_abf, PP.H3_abf),
                              "UNEVALUABLE ABF -- never tested", "PASSES coloc -- this is a solid edge")]
  stop(sprintf("32: %d audit dashed edge(s) do not satisfy the dashed premise (evaluated by coloc.abf AND not passing): %s -- the audit's classification is stale. Refusing to finalize.",
               nrow(.bad_prem), paste(.bad_prem[, sprintf("%s~%s(H4=%s, H3=%s: %s)", gene, outcome_name,
                                                          format(PP.H4_abf), format(PP.H3_abf), .why)],
                                      collapse = ", ")))
}

# ── the strict re-run result (credible-set purity ON, min_abs_corr=0.5; `check_LD` in 31's patch is a
#    no-op -- not a susieR parameter. See REBUILD_LEDGER 2.5 and 31's header.) ──
# CHECK EACH FILE BEFORE THE FILL, NOT THE FILLED TABLE. `rbindlist(..., fill = TRUE)` CREATES a column
# that only one file has and fills the other with NA -- so a check placed after it passes whenever ONE mode
# was re-run and the other is stale, which is the likeliest way to get here. The guard I first wrote sat
# after this line and was defeated by the exact mechanism it existed to catch: "a guard that runs after the
# damage is not a guard" is already recorded in this project as an earlier finding, and I repeated it.
.sens_each <- lapply(SENS, fread)
for (.i in seq_along(SENS)) {
  .m <- setdiff(c("strict_converged_exp", "strict_converged_out"), names(.sens_each[[.i]]))
  if (length(.m))
    stop(sprintf("32: %s lacks %s -- it predates 31's per-side convergence recording. rbindlist(fill=TRUE) would create the column from the OTHER mode's file and fill this one with NA, so dashed_edges_final.csv would carry 'convergence unknown' for fits that ran, silently. Re-run 31 for BOTH modes.",
                 basename(SENS[.i]), paste(.m, collapse = ", ")))
}
sens <- rbindlist(.sens_each, fill = TRUE)
# PRESENCE IS NOT CORRECTNESS. A file that was ALREADY tainted by an earlier fill=TRUE HAS both columns --
# every value NA -- and a presence check waves it through. Check the premise instead: a row that shows
# evidence the fit RAN (a posterior, or an LD SNP count) cannot honestly say "convergence unknown". This is
# the failure mode this layer has a name for -- verifying a label exists rather than that its premise holds
# -- committed inside the guard written to fix that very class.
# The first version of this check said "a posterior, or an LD SNP count" in its comment and then tested
# only `run_status == "ok" & !is.na(strict_PP.H4_susie)` -- a comment describing the right check above the
# wrong one. It missed every strict post-fit collapse (`run_status == "collapsed"`,
# collapse_reason == "empty_coloc_summary"), which is exactly where the fit ran and the convergence is known.
# susie_fit_ran() in R/00_config.R is the single definition; 31 renames the columns, so point it at the
# strict_* names.
assert_convergence_honest(sens, "32: the sensitivity output",
                          conv = c("strict_converged_exp", "strict_converged_out"),
                          posterior = "strict_PP.H4_susie", n_exp = "strict_n_exp",
                          n_out = "strict_n_out", n_ld = "strict_n_snps_ld")
# The audit is truth-checked above, but draw_status is decided by the STRICT values, and those come from
# HERE -- so this is the table that had to be validated and was not (r26 High). A stale sensitivity row
# with the right (gene, outcome_name) but an old outcome_id or the wrong engine passes the row-count, the
# duplicate check and the join, and then draws/holds/drops on a different GWAS. Same shape as r23's
# finding in 31: the guard was on the table that is read, not on the table that decides.
assert_outcome_truth(sens, "candidate_susie_sensitivity_*.csv")
assert_engine_truth(sens, "candidate_susie_sensitivity_*.csv")
if (anyDuplicated(sens, by = c("gene", "outcome_name"))) stop("32: sensitivity rows not unique")
if (nrow(sens) != nrow(dash))
  stop(sprintf("32: sensitivity rows=%d != dashed audit rows=%d -- re-run 31 for both modes on the current dashed set", nrow(sens), nrow(dash)))
# Carry outcome_id/susie_engine through so they can be checked against the audit, instead of dropping
# them here and joining on the name alone (r26 High).
# strict_PP.H3_susie was NOT carried until r39, so the round-38 claim that "32 repeats the strict PP
# integrity check" was false for H3 -- the column never arrived to be checked.
# UNLIKE 29/31's setcolorder, `sens[, ..scols]` below IS a projection: a column not named here does not
# reach dashed_edges_final.csv at all. That is how the per-side convergence was lost at the last step --
# 15/26/16 record it, 29/31 carry it, and this line dropped it, leaving reporting to infer convergence from
# `n_warn` again. `n_warn` is counted per EDGE across both susie_rss calls and cannot say WHICH side failed;
# that is the whole reason `$converged` is recorded per side. Eighth layer of one change: every layer was
# "the recording does not reach the next step".
scols <- c("gene", "outcome_name", "outcome_id", "susie_engine", "orig_PP.H4_susie",
           "strict_PP.H4_susie", "strict_PP.H3_susie",
           "strict_converged_exp", "strict_converged_out",
           "run_status", "survives", "n_warn", "ld_warn", "collapse_reason", "warning_msg")
# REQUIRE the columns rather than discovering they are absent: `sens[, ..scols]` errors on a missing
# column, but the error would name a column and not the reason, and the reason is what the next reader needs.
.sc_miss <- setdiff(scols, names(sens))   # the per-file check above covers convergence; this covers the rest of scols
if (length(.sc_miss))
  stop(sprintf("32: the sensitivity output lacks %s -- `sens[, ..scols]` is a projection, so a column missing here does not reach dashed_edges_final.csv at all. Re-run 31 (both modes).",
               paste(.sc_miss, collapse = ", ")))
m <- merge(dash, sens[, ..scols], by = c("gene", "outcome_name"), all.x = TRUE,
           suffixes = c("", ".sens"))
if (nrow(m) != nrow(dash) || anyNA(m$run_status))
  stop("32: dashed <-> sensitivity join is incomplete -- every dashed edge must have a sensitivity row")
# The audit and the sensitivity run must agree on WHICH GWAS each edge was tested against. NA-aware: a
# one-sided NA is a disagreement, not a match (the shape fixed in 31 at r23).
.ne <- function(x, y) xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
# `.ne` treats both-NA as agreement -- correct for fields that may legitimately be absent, WRONG here:
# every row in `m` is a dashed edge, and a dashed edge went through SuSiE by definition, so an engine-less
# one cannot exist. Without this, two stale pre-stamp rows with matching outcome_id and no engine at all
# would "agree" and be finalized with no engine provenance (r28 High). Reject NA on both sides FIRST, so
# the agreement check below can never certify a consensus of absences.
# Strict PP RANGE, not just presence -- `Inf >= 0.80` is TRUE, so a malformed value would be read as
# having cleared the threshold and drawn. 31 asserts this producer-side; 32 is the consumer and the last
# gate before dashed_edges_final.csv, so it checks too (r38 High).
assert_pp_integrity(m, c("strict_PP.H4_susie", "strict_PP.H3_susie"),
                    "the sensitivity rows joined to the dashed audit (32's last gate before dashed_edges_final.csv)")
assert_status_premise(m, c("strict_PP.H4_susie", "strict_PP.H3_susie"),
                      "the sensitivity rows joined to the dashed audit (32's last gate)")
# `collapsed` becomes `dropped:collapsed`, which the audit and the caption read as "no pure credible set"
# -- a benign, biological outcome. assert_status_premise() only proves such a row carries no posterior; it
# cannot tell a benign collapse from a computation failure wearing the same label. CAND_BENIGN_TAGS is the
# single source for that distinction, and neither 31's write nor this consumer asked it (r43 High).
.cr_bad <- m[run_status == "collapsed" &
             (is.na(collapse_reason) | !collapse_reason %in% CAND_BENIGN_TAGS)]
if (nrow(.cr_bad))
  stop(sprintf("32: %d edge(s) are run_status=='collapsed' but their collapse_reason is missing or not a benign tag (%s): %s -- `dropped:collapsed` asserts a BENIGN collapse. An unrecognised reason may be a computation failure, and freezing it as a biological negative is exactly the type-(b)-must-not-look-benign rule this layer exists to enforce.",
               nrow(.cr_bad), paste(CAND_BENIGN_TAGS, collapse = "/"),
               paste(.cr_bad[, sprintf("%s~%s(reason=%s)", gene, outcome_name,
                                       as.character(collapse_reason))], collapse = ", ")))
.cr_ok <- m[run_status == "ok" & !is.na(collapse_reason) & nzchar(trimws(as.character(collapse_reason)))]
if (nrow(.cr_ok))
  stop(sprintf("32: %d edge(s) are run_status=='ok' yet carry a collapse_reason: %s -- a row cannot both have succeeded and report why it collapsed",
               nrow(.cr_ok), paste(.cr_ok[, sprintf("%s~%s(reason=%s)", gene, outcome_name,
                                                    as.character(collapse_reason))], collapse = ", ")))

.na_eng <- m[is.na(susie_engine) | susie_engine == "" | is.na(susie_engine.sens) | susie_engine.sens == ""]
if (nrow(.na_eng))
  stop(sprintf("32: %d dashed edge(s) carry no susie_engine in the audit and/or the sensitivity run: %s -- a dashed edge ran through SuSiE by definition; an engine-less row means a stale or pre-stamp artifact. Refusing to finalize.",
               nrow(.na_eng), paste(.na_eng[, sprintf("%s~%s(audit=%s, sens=%s)", gene, outcome_name,
                                                      susie_engine, susie_engine.sens)], collapse = ", ")))
.oc <- m[.ne(outcome_id, outcome_id.sens) | .ne(susie_engine, susie_engine.sens)]
if (nrow(.oc))
  stop(sprintf("32: %d edge(s) where the audit and the sensitivity run disagree on outcome_id or susie_engine: %s -- the two describe a different GWAS or a different engine for the same pair. Refusing to finalize.",
               nrow(.oc), paste(.oc[, sprintf("%s~%s(audit id=%s/eng=%s, sens id=%s/eng=%s)", gene, outcome_name,
                                              outcome_id, susie_engine, outcome_id.sens, susie_engine.sens)],
                                collapse = ", ")))
# The comment above this block named the "wrong engine" case from the start; until r27 only outcome_id was
# actually compared, so the case it named was the one case not caught. Both sides are also checked against
# the outcome-derived engine by assert_engine_truth() above -- agreement between two stale labels is not
# truth, which is failure mode (v) all over again.
m[, c("outcome_id.sens", "susie_engine.sens") := NULL]
# a coloc.susie ERROR (not a benign collapse) must be resolved before finalizing (would bias the drawn set)
if (m[run_status == "error", .N])
  stop(sprintf("32: %d edge(s) have run_status=='error' (coloc.susie computation failure, NOT a benign collapse): %s -- investigate before finalizing",
               m[run_status == "error", .N], paste(m[run_status == "error", paste0(gene, "~", outcome_name)], collapse = ", ")))

# ── rule (1) FDR: re-assert every candidate dashed edge is FDR-significant (external cross-check) ──
# 2026-07-15: this read `mr_significant_fdr.rds` -- BH within each method family over analysable rows
# only, before the primary dedup (cascade/01_mr_screen.R's `fdr_sig` line; it is also NOT the same
# computation as the sensitivity the manuscript prints -- three surfaces, see the block above
# primary_fdr_set() in R/00_config.R) -- while the rule at the top of this file calls
# it "MR FDR q<0.05", i.e. the manuscript's PRIMARY definition (BH once over the full grid, untested
# cells padded at p=1). The gate was asserting a rule it did not describe. primary_fdr_set()
# (R/00_config.R) now OWNS that definition -- though not yet as its only implementation; the reporting
# layer still re-implements it and cascade/14 still reads the saved by-method artifact (REBUILD_LEDGER §2.1c).
# 28 builds the pool from the same call, so this cross-check stays external to the pool while testing the
# same rule the file claims.
# Verified on the frozen data before the change: all 17 candidate dashed edges are FDR-significant under
# BOTH definitions, which PREDICTS this gate's verdict will not move -- only the rule it enforces now
# matches the rule it states. The gated re-run is what confirms that; until it runs, this file still
# reads the pre-fix audit.
# Keyed two-column join, not paste()-strings: rule (1) is the gate this whole file exists to enforce, and
# a string key would make it silently wrong under a collision. The Bonferroni sets were converted in r24
# and this one was left -- the same instance-vs-adjacent-instance miss, in the more load-bearing of the
# two (r25 Medium).
# `mr_direction` rides from the audit into dashed_edges_final.csv and becomes the EDGE COLOUR on the
# published figure (red = higher protein, higher risk; blue = protective). Nothing checked it: 30 derives
# it from `b` after the eGFR flip, but 32 selected the label and dropped `b`, so a stale audit gives a
# drawn edge the WRONG colour and every gate here still passes (r43 High).
# `mr` is already read for rule (1) below; derive the direction from it and cross-check.
# 30's rule (30:221): beta_star = -b for eGFR (risk-aligned), else b; direction = protective if < 0.
.mr_pm <- as.data.table(readRDS(MR_ALL))[method %in% c("Wald ratio", "Inverse variance weighted"),
                                         .(gene, outcome_name, .b = b)]
# Uniqueness BEFORE the merge. primary_fdr_set()'s duplicate-primary guard runs further down, which is
# after this merge has already had its chance to fan `m` out -- and every direction check below would then
# be operating on multiplied rows (r51 Medium). A guard that runs after the damage is not a guard.
if (anyDuplicated(.mr_pm, by = c("gene", "outcome_name")))
  stop(sprintf("32: mr_all_results.rds has %d (gene, outcome_name) pair(s) with more than one primary-method row -- merging that into the dashed table would multiply its rows before any check sees them",
               sum(duplicated(.mr_pm, by = c("gene", "outcome_name")))))
m <- merge(m, .mr_pm, by = c("gene", "outcome_name"), all.x = TRUE)
# FINITE, not merely non-NA. `Inf < 0` is FALSE, so an infinite effect size silently becomes "risk" and
# colours an edge on the figure; the r43 comment claimed this required a finite b while the code only
# rejected NA (r44 High).
if (m[!is.finite(.b), .N])
  stop(sprintf("32: %d dashed edge(s) have a missing or non-finite primary MR effect size in mr_all_results.rds: %s -- the edge direction cannot be derived from it, and that direction is the edge colour",
               m[!is.finite(.b), .N], paste(m[!is.finite(.b), sprintf("%s~%s(b=%s)", gene, outcome_name,
                                                format(.b))], collapse = ", ")))
m[, .beta_star := fifelse(outcome_name == "eGFR", -.b, .b)]
m[, .dir_rd := fifelse(.beta_star < 0, "protective", "risk")]
.dd <- m[is.na(mr_direction) | mr_direction != .dir_rd]
if (nrow(.dd))
  stop(sprintf("32: %d dashed edge(s) where the audit's mr_direction disagrees with the direction derived from mr_all_results.rds (eGFR flipped): %s -- this label is the edge colour on the figure. Refusing to finalize.",
               nrow(.dd), paste(.dd[, sprintf("%s~%s(audit=%s, derived=%s, b=%s)", gene, outcome_name,
                                              as.character(mr_direction), .dir_rd, format(.b))], collapse = ", ")))
m[, mr_direction := .dir_rd][, c(".b", ".beta_star", ".dir_rd") := NULL]

.fdr_set <- unique(primary_fdr_set(as.data.table(readRDS(MR_ALL)))[, .(gene, outcome_name)])[, .fdr := TRUE]
m <- merge(m, .fdr_set, by = c("gene", "outcome_name"), all.x = TRUE)
m[, is_fdr := !is.na(.fdr)][, .fdr := NULL]
if (!all(m$is_fdr))
  stop(sprintf("32: %d dashed edge(s) NOT in the FDR-significant set (rule 1 violated): %s",
               m[!is_fdr, .N], paste(m[!is_fdr, paste0(gene, "~", outcome_name)], collapse = ", ")))

# An `ok` row with no strict_PP.H4_susie is a CONTRADICTION, and fcase below would resolve it silently to
# `dropped:strict_fail` via survives==FALSE -- turning a computation/schema failure into the claim that the
# edge simply did not reach 0.80. 29 and 31 now assert usability before stamping `ok`; this is the
# belt-and-suspenders stop at the point where the mislabelling would actually happen (r15 High).
# Measured 2026-07-16: all 4 current dropped:strict_fail edges carry real values; only the 4
# dropped:collapsed rows are NA, which is correct (collapsed = no pure credible set).
if (m[run_status == "ok" & is.na(strict_PP.H4_susie), .N])
  stop(sprintf("32: %d edge(s) have run_status=='ok' but NO strict_PP.H4_susie: %s -- an ok row must carry a value; this is a 29/31 schema or engine failure, NOT a strict-threshold failure. Refusing to finalize.",
               m[run_status == "ok" & is.na(strict_PP.H4_susie), .N],
               paste(m[run_status == "ok" & is.na(strict_PP.H4_susie), paste0(gene, "~", outcome_name)], collapse = ", ")))

# ── RE-DERIVE rule arm (2) rather than trusting 31's `survives` column ───────────────────────
# This file's header claims it applies the locked rule INDEPENDENTLY. Until 2026-07-16 that was true of
# arm (1) -- which re-derives the primary FDR set from mr_all_results.rds -- but false of arm (2), which
# simply read `survives` out of 31's CSV. `survives` is 31's own encoding of arm (2)
# (`run_status == "ok" & !is.na(strict_PP.H4_susie) & strict_PP.H4_susie >= 0.80`), so classifying from it
# was checking the rule against its own output: a TAUTOLOGY, the exact failure the editor flagged for
# `concordant_susie`. A stale, threshold-drifted or malformed `survives` would have been obeyed silently,
# and no guard here would have noticed (r20 High).
# Re-derive from the raw columns, then CROSS-CHECK against the upstream flag. The cross-check is the part
# with evidential value: it can only fail if 31 and 32 disagree about what the rule says.
# Verified on the frozen data before the change: 17 rows, 9 survivors either way, identical row by row,
# and 0 rows with survives==TRUE but strict_PP.H4_susie < 0.80. So this costs nothing today.
# `is_bonferroni` rides in from the audit and lands in dashed_edges_final.csv unexamined. The draw rule
# never consults Bonferroni, so NOTHING here would notice a stale label -- the same provenance-label
# failure shape just removed from 28 (r22 Medium). `mr` is already read above for rule (1); derive the
# set and cross-check, then emit the DERIVED value. Frozen-data check 2026-07-16: saved 75 pairs,
# derived 75, setequal TRUE.
# Keyed two-column set + join, not paste()-strings: a space in a gene or outcome name would make the
# string key ambiguous and the membership test silently wrong (r24). Latent on this vocabulary; fixed
# because these lines were open.
.bonf_rd <- unique(as.data.table(readRDS(MR_ALL))[
  method %in% c("Wald ratio", "Inverse variance weighted") & pval < BONFERRONI_ALPHA,
  .(gene, outcome_name)])[, .bonf := TRUE]
m <- merge(m, .bonf_rd, by = c("gene", "outcome_name"), all.x = TRUE)
m[, is_bonferroni_rd := !is.na(.bonf)][, .bonf := NULL]
# as_logical_strict(), not as.logical(): a present-but-unparseable label coerces to NA, the comparison
# goes NA, the filter drops it, and this cross-check reports agreement about a label it never read (r48 --
# the same bug r47 found in derive_runnability, repeated here because the guard was written as a LOCAL
# helper inside that function instead of a shared one).
m[, .bonf_lab := as_logical_strict(is_bonferroni, "candidate_edge_audit (dashed)", "is_bonferroni")]
.bd <- m[is.na(.bonf_lab) | .bonf_lab != is_bonferroni_rd]
if (nrow(.bd))
  stop(sprintf("32: %d edge(s) where the audit's is_bonferroni disagrees with the set derived from mr_all_results.rds at BONFERRONI_ALPHA=%.4e: %s -- the audit label is stale. Refusing to carry it into dashed_edges_final.csv.",
               nrow(.bd), BONFERRONI_ALPHA,
               paste(.bd[, sprintf("%s~%s(audit=%s, derived=%s)", gene, outcome_name, is_bonferroni, is_bonferroni_rd)],
                     collapse = ", ")))
m[, is_bonferroni := is_bonferroni_rd][, is_bonferroni_rd := NULL]

# Rule arm (3) draws on `n_warn == 0`, and `n_warn` is a stored COUNT sitting next to the warning text
# it counts. Trusting it while ignoring `warning_msg`/`ld_warn` is the tenth instance of "the checked
# thing is not the deciding thing" (r36 High): a stale row with n_warn=0 and a non-empty warning_msg, or
# ld_warn > n_warn, would be DRAWN under a rule that says zero warnings of any kind.
# Frozen data 2026-07-16: 17 rows, exactly two shapes -- 14 with 0/0/no-message and 3 with 1/1/message
# (the three held non-convergent edges). 0 incoherent rows, so this costs nothing today.
.has_wmsg <- !is.na(m$warning_msg) & nzchar(trimws(as.character(m$warning_msg)))
# `ld_warn` is REQUIRED to be sane, not "sane if present". The first version wrote
# `!is.na(ld_warn) & (...)`, so an NA ld_warn passed -- one-sided NA read as agreement, the fourth
# instance of that shape and the first inside the fix written for another instance of it (r37 High).
# 31 always writes ld_warn; frozen data has 0 NA.
# `ld_warn` is a COUNT of a regex hit in `warning_msg`, and both sit in the same row -- so check the
# count against the text it counts, the same way n_warn is checked against the presence of text (r52
# Medium). 31 computes it as grepl("LD|correlation|inconsist|kriging|check_ld|refin", ...); re-derive a
# PRESENCE check with the same vocabulary rather than trusting a stale count.
# NOTE what ld_warn does and does not mean: on the three held edges it fires on the word "LD" inside
# susieR's remediation advice ("Please check consistency between summary statistics and LD matrix"),
# NOT on an LD warning as such. It means "the message mentions LD". `held:non_converged` is therefore
# checked against the convergence text (is_nonconvergence_warning), never against this.
.LD_WARN_PAT <- "LD|correlation|inconsist|kriging|check_ld|refin"
.has_ld_text <- !is.na(m$warning_msg) & grepl(.LD_WARN_PAT, m$warning_msg, ignore.case = TRUE)
.wi <- m[!is.finite(n_warn) | n_warn < 0 | n_warn != as.integer(n_warn) |
         !is.finite(ld_warn) | ld_warn < 0 | ld_warn != as.integer(ld_warn) | ld_warn > n_warn |
         (n_warn == 0) != !.has_wmsg |
         (ld_warn > 0) != .has_ld_text]
if (nrow(.wi))
  stop(sprintf("32: %d edge(s) whose warning bookkeeping is incoherent -- n_warn must be a non-negative whole number, ld_warn must not exceed it, and n_warn==0 must mean no warning text: %s. Rule (3) draws on n_warn alone, so an incoherent row could be drawn under a rule that requires zero warnings.",
               nrow(.wi), paste(.wi[, sprintf("%s~%s(n_warn=%s, ld_warn=%s, msg=%s)", gene, outcome_name,
                                              n_warn, ld_warn,
                                              fifelse(is.na(warning_msg) | !nzchar(trimws(as.character(warning_msg))),
                                                      "<none>", "<present>"))], collapse = ", ")))

m[, survives_rd := run_status == "ok" & !is.na(strict_PP.H4_susie) & strict_PP.H4_susie >= 0.80]
m[, .surv_lab := as_logical_strict(survives, "candidate_susie_sensitivity_*", "survives")]
.sd <- m[is.na(.surv_lab) | .surv_lab != survives_rd]
if (nrow(.sd))
  stop(sprintf("32: %d edge(s) where 31's `survives` disagrees with rule arm (2) re-derived here (run_status=='ok' & strict_PP.H4_susie >= 0.80): %s -- 31 and 32 do not agree on the locked rule. Refusing to finalize.",
               nrow(.sd), paste(.sd[, sprintf("%s~%s(survives=%s, run_status=%s, strict=%s)",
                                              gene, outcome_name, survives, run_status,
                                              format(strict_PP.H4_susie))], collapse = ", ")))

# ── DERIVE draw_status from the locked rule (counts fall out of the data; nothing hard-coded) ──
# `held:non_converged` ASSERTS non-convergence, and the LOCKED caption repeats that assertion. Verify it
# against the warning text before the label is written, rather than inferring it from "was warned at all"
# (r51 High). Frozen data: all 3 held edges say "IBSS algorithm did not converge in 100 iterations!".
.held <- m[survives_rd == TRUE & n_warn > 0]
if (nrow(.held)) {
  .not_nc <- .held[!is_nonconvergence_warning(warning_msg)]
  if (nrow(.not_nc))
    stop(sprintf("32: %d surviving edge(s) would be labelled `held:non_converged` but their warning does not say the fit failed to converge: %s -- the label is a claim the caption repeats verbatim. Held for a different reason needs a different status, and that is an editorial decision, not a silent relabel.",
                 nrow(.not_nc), paste(.not_nc[, sprintf("%s~%s(msg=%s)", gene, outcome_name,
                                          substr(gsub("[\r\n]+", " ", warning_msg), 1, 60))], collapse = ", ")))
}
m[, draw_status := fcase(
  survives_rd == TRUE & n_warn == 0,           "drawn",
  survives_rd == TRUE & n_warn >  0,           "held:non_converged",
  run_status == "collapsed",                   "dropped:collapsed",
  run_status == "ok" & survives_rd == FALSE,   "dropped:strict_fail",
  default = NA_character_)]

# ── structural assertions (partition completeness; NOT a hard-coded count, NOT a tautology) ──
if (m[is.na(draw_status), .N])
  stop(sprintf("32: %d edge(s) unclassified by the rule -- rule gap, refusing: %s",
               m[is.na(draw_status), .N], paste(m[is.na(draw_status), paste0(gene, "~", outcome_name, "(", run_status, ")")], collapse = ", ")))
if (nrow(m) != nrow(dash)) stop("32: partition lost/gained rows -- refusing")

lev <- c("drawn", "held:non_converged", "dropped:collapsed", "dropped:strict_fail")
setorder(m, -strict_PP.H4_susie, -orig_PP.H4_susie, na.last = TRUE)
m <- m[order(factor(draw_status, levels = lev))]
m[, c(".bonf_lab", ".surv_lab") := NULL]   # temporaries, never into dashed_edges_final.csv
fwrite(m, OUT)

# ── report (DERIVED counts, printed) ──
n_drawn <- m[draw_status == "drawn", .N]; n_held <- m[draw_status == "held:non_converged", .N]; n_drop <- m[draw_status %like% "^dropped", .N]
cat(sprintf("\n32: dashed-edge finalization -> %s\n", OUT))
cat("\n== draw_status (DERIVED from the locked rule) ==\n"); print(m[, .N, by = draw_status][order(-N)])
cat(sprintf("\nDRAWN=%d | HELD(non-converged)=%d | DROPPED=%d   (of %d evaluated candidate dashed edges)\n", n_drawn, n_held, n_drop, nrow(m)))
if (!(n_drawn == 7L && n_held == 1L && n_drop == 8L))
  cat(sprintf("32: [note] counts differ from the 2026-07-18 editor-ruled baseline (7/1/8) -- re-check with the editor.\n"))
cat("\n== DRAWN (the dashed edges to draw in Stage 1) ==\n")
print(m[draw_status == "drawn", .(gene, outcome_name, coloc_abf = PP.H4_abf, strict_susie = strict_PP.H4_susie, mr_direction, is_bonferroni)])
cat("\n== HELD: non-converged, evaluated-but-not-drawn (in-sample LD unobtainable for CKDGen; held for non-convergence) ==\n")
print(m[draw_status == "held:non_converged", .(gene, outcome_name, strict_susie = strict_PP.H4_susie, n_warn, ld_warn)])
cat("\n== DROPPED under strict settings ==\n")
print(m[draw_status %like% "^dropped", .(gene, outcome_name, orig = orig_PP.H4_susie, strict = strict_PP.H4_susie, draw_status)])
