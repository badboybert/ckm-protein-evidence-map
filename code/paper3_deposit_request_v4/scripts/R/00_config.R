# Paper 3: Proteome-wide MR for CKM Drug Targets
# Configuration: paths, parameters, API tokens

# ── Base Paths ──────────────────────────────────────────────────────────────────
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
BASE_DIR    <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
DATA_DIR    <- file.path(BASE_DIR, "data")
UKB_DIR     <- file.path(DATA_DIR, "ukb_ppp/european")
DECODE_DIR  <- file.path(DATA_DIR, "decode")
OUTCOME_DIR <- file.path(DATA_DIR, "outcome_gwas")
# Result/phase dirs are overridable (pre-set before sourcing) so the corrected cascade
# can repoint them at results/repro_corrected/ without editing every downstream script.
if (!exists("RESULTS_DIR")) RESULTS_DIR <- file.path(BASE_DIR, "results")
if (!exists("FIGURES_DIR")) FIGURES_DIR <- file.path(BASE_DIR, "figures")
SCRIPTS_DIR <- file.path(BASE_DIR, "scripts")

# PHASE1_DIR FOLLOWS THE SELECTED PHASE-1 SET. It defaulted to results/phase1_instruments
# independently of CKM_PHASE1_SET, so selecting the rebuilt generation moved PHASE1_INSTRUMENTS and
# left this pointing at the 1,635-gene table -- and localize/02 and /03 read THIS one. The documented
# commands for those two scripts would therefore have gone on building outcome extracts and the
# proxy cache from the previous generation, omitting 427 rebuilt proteins and keeping 26 the rebuild
# does not have, while every log line said the rebuilt set had been selected. One selector.
# (the selector is sourced once, at the top of this file -- sourcing it again re-ran its rebuild
#  guard and its selection message, and could fail on the second call after succeeding on the
#  first if the rebuild state changed in between)
if (!exists("PHASE1_DIR")) PHASE1_DIR <- dirname(PHASE1_INSTRUMENTS)
if (!exists("PHASE2_DIR")) PHASE2_DIR <- file.path(RESULTS_DIR, "phase2_mr_screen")
if (!exists("PHASE3_DIR")) PHASE3_DIR <- file.path(RESULTS_DIR, "phase3_coloc")
if (!exists("PHASE4_DIR")) PHASE4_DIR <- file.path(RESULTS_DIR, "phase4_annotation")
# ONE ROOT PER SESSION. All four are `if (!exists())`, which is what lets a caller repoint them --
# but it also means that sourcing this file, changing RESULTS_DIR, and sourcing it again leaves
# RESULTS_DIR at the new root while PHASE2/3/4_DIR still point into the old one. The result is a run
# split across two generations with every individual assignment looking correct. A caller repointing
# them deliberately (04 sets PHASE3_DIR to the generation directory) keeps them under RESULTS_DIR, so
# this rejects the split without rejecting the intended override.
local({
  .r <- normalizePath(RESULTS_DIR, winslash = "/", mustWork = FALSE)
  for (.n in c("PHASE2_DIR", "PHASE3_DIR", "PHASE4_DIR")) {
    .p <- normalizePath(get(.n), winslash = "/", mustWork = FALSE)
    if (!startsWith(.p, .r))
      stop(sprintf(paste0("config: %s is %s, which is not under RESULTS_DIR (%s). This session would ",
                          "read one generation's phase directories while writing another's -- source ",
                          "config once per process, or repoint both together."), .n, .p, .r))
  }
})

PYTHON_BIN  <- paste0(.P31R_EXT, "/C/Users/User/AppData/Local/Programs/Python/Python310/python.exe")
RSCRIPT_BIN <- file.path(R.home("bin"), "Rscript")

# 2026-07-01 (gate finding 1): resolve a gene's UKB-PPP tar by the INSTRUMENT's OID. Multi-assay proteins
# (e.g. IDO1, TNF) have several panel tars; picking [1] can make coloc/SuSiE/prior test a DIFFERENT assay
# than the MR exposure. Uses the canonical (assay-deduped, single-OID/gene) instrument set that
# PHASE1_INSTRUMENTS names -- whichever generation that is, not a fixed 1,682.
# 2026-07-01 (closure B1): fail-loud for DATA-FOUNDATION failures that must NOT be demoted to a skip/NA by
# an outer per-pair `tryCatch(error=...)`. A re-raised *classed* condition does NOT reliably survive R's
# nested-tryCatch dispatch (a `stop(e)` re-raise inside a data_failure handler gets caught by a sibling
# `error=` handler), so we set a SESSION FLAG instead: the producer loop re-checks `.CKM_DATA_FAILURE`
# right after each per-pair tryCatch and aborts. stop() still fires (demoted to the per-pair skip), but the
# flag survives it. Fresh Rscript per step => no cross-script leakage.
fail_loud_data <- function(msg) { assign(".CKM_DATA_FAILURE", msg, envir = .GlobalEnv); stop(msg) }

resolve_ukbppp_tar <- function(gene) {
  .g <- gene
  tars <- list.files(file.path(DATA_DIR, "ukb_ppp/european"), pattern = paste0("^", tolower(.g), "_"), full.names = TRUE)
  # A ZERO-MATCH IS NOT A MISSING ASSAY. The filename carries the symbol that was current when
  # UKB-PPP released the data, and 15 genes in this instrument set have been renamed since --
  # RIGI's archive is ddx58_..._oid21226_..., WARS1's is wars_..._oid21084_.... This used to stop
  # here, before the OID map below (which identifies the archive unambiguously) was ever built,
  # aborting colocalisation on 17 of the current suggestive pairs.
  #
  # THE SINGLE-MATCH BRANCH IS CHECKED TOO. It used to `return(tars[1])` right here, before the map
  # was consulted -- so the comment promising resolution "by the instrument's OID" was false for the
  # most common path, and a gene whose one symbol-matched archive carried a DIFFERENT OID than Phase 1
  # selected was accepted in silence: MR on one assay, colocalisation on another. The map is built
  # first now, so all three branches (one match, none, several) answer to the same pinned OID.
  # THE CACHE IS TIED TO THE INSTRUMENT SET IT CAME FROM. It was built once per session and kept
  # across a generation switch, so after switching TNF still resolved OID20074 from the previous map
  # instead of this set's OID20848, and a gene absent from the stale map fell through to "first tar"
  # -- CXCL8 resolved OID20153 rather than its instrument's OID21430. That opens the WRONG assay's
  # archive: MR on one assay, colocalisation on another.
  # ...and it is built from THE BYTES THE THRESHOLD WAS COMPUTED FROM. This used to re-open
  # PHASE1_INSTRUMENTS by path, hours after config was sourced, and stamp whatever it found with the
  # md5 recorded at source time. A rebuild landing in between therefore put a NEW table's OID map
  # beside the OLD table's denominator, wearing the old table's hash -- so the very check meant to
  # catch a generation switch certified the mismatch instead. One read, one table, one identity.
  if (!exists(".UKBPPP_OIDMAP", envir = .GlobalEnv) ||
      !identical(attr(get(".UKBPPP_OIDMAP", envir = .GlobalEnv), "ckm_instruments_md5"),
                 INSTRUMENTS_MD5)) {
    .i <- as.data.frame(read_phase1_instruments())
    .map <- unique(.i[, c("gene", "oid")])
    attr(.map, "ckm_instruments_md5") <- INSTRUMENTS_MD5
    assign(".UKBPPP_OIDMAP", .map, envir = .GlobalEnv)
  }
  .mm <- get(".UKBPPP_OIDMAP", envir = .GlobalEnv)
  .oid <- .mm$oid[.mm$gene == .g]
  .oid_in <- function(p) grepl(paste0("_oid", gsub("[^0-9]", "", .oid), "_"), tolower(basename(p)))
  # ONE SYMBOL MATCH: the common path, and now the OID is verified on it rather than assumed. A gene
  # with no OID in the instrument set keeps the previous permissive behaviour -- that is a different
  # problem and is not made fatal here.
  if (length(tars) == 1L) {
    if (length(.oid) == 1L && !.oid_in(tars))
      stop(sprintf("resolve_ukbppp_tar: gene %s resolves to %s, which does not carry its instrument's OID %s -- MR would use one assay and colocalisation another; data-foundation failure, NOT a legit skip", .g, basename(tars), .oid))
    return(tars[1])
  }
  # THE RENAMED-SYMBOL FALLBACK. Requires exactly one instrument OID for the gene and exactly one
  # archive carrying that OID as a whole `_oid<digits>_` token -- the underscores stop OID21226 from
  # matching an oid212260 archive. Zero or more than one is fatal: "first tar" is precisely how this
  # function once opened one assay's archive while MR had used another.
  if (length(tars) == 0L) {
    if (length(.oid) != 1L)
      stop(sprintf("resolve_ukbppp_tar: no tar matches gene symbol %s and the instrument set gives %d OID(s) for it -- refusing to guess an archive", .g, length(.oid)))
    tars <- list.files(file.path(DATA_DIR, "ukb_ppp/european"),
                       pattern = paste0("_oid", gsub("[^0-9]", "", .oid), "_.*[.]tar$"), full.names = TRUE)
    if (length(tars) != 1L)
      stop(sprintf("resolve_ukbppp_tar: gene %s (OID %s) matched %d archives by OID; expected exactly one -- data-foundation failure, NOT a legit skip", .g, .oid, length(tars)))
    return(tars[1])
  }
  # ONE INSTRUMENT OID, ONE MATCH, TOKEN-BOUNDED -- the same three requirements as the renamed-symbol
  # fallback above, which this path did not share. `return(tars[1])` on a gene with no mapped OID and
  # `tars[.hit[1]]` on several matches are both "pick the first", which is exactly how this function
  # once opened one assay's archive while MR had used another; and an unbounded "oid123" also matches
  # oid1234. TNF has four assay archives, so the multi-match path is not hypothetical.
  if (length(.oid) != 1L)
    stop(sprintf("resolve_ukbppp_tar: gene %s matches %d archives and the instrument set gives %d OID(s) for it -- refusing to pick one", .g, length(tars), length(.oid)))
  .hit <- grep(paste0("_oid", gsub("[^0-9]", "", .oid), "_"), tolower(basename(tars)), fixed = FALSE)
  if (length(.hit) != 1L)
    stop(sprintf("resolve_ukbppp_tar: gene %s (OID %s) matched %d of its %d archives by OID; expected exactly one -- data-foundation failure, NOT a legit skip", .g, .oid, length(.hit), length(tars)))
  tars[.hit]
}

# ── Analysis Parameters ─────────────────────────────────────────────────────────
P_INSTRUMENT  <- 1e-6   # Instrument p-value threshold (initial screen)
P_GWAS        <- 5e-8   # Genome-wide significance (sensitivity analysis)
CIS_WINDOW_KB <- 500    # Cis window: ±500 kb from gene boundary
# The window a COLOC REGION is built with, as distinct from the window INSTRUMENTS were selected in.
# One name used to carry both, and that conflation is what broke the pre-specified +/-1 Mb window
# sensitivity: 55_coloc_1mb_sensitivity.R reassigned CIS_WINDOW_KB to widen the coloc region, which
# also made the session claim its instruments were selected at 1 Mb -- so the gene-coordinate
# certificate (which records the window 03b checked the Phase-1 cis boundaries against) correctly
# refused the table for 56 of 74 pairs. The certificate was right; the variable was wrong.
#
# These are genuinely different quantities. The certificate attests that Phase 1 selected instruments
# within CIS_WINDOW_KB of each gene body, which stays true whatever region a later analysis integrates
# over. The gene coordinates themselves carry no window at all. R/coloc.R's .PROV_KEY already treated
# the region width as a run-level parameter that changes every pair's result.
#
# Default equal, so nothing changes unless a script sets it deliberately. Only the coloc-region
# builders in R/coloc.R read it; the certificate guard (R/coloc.R:531, R/00_config.R:1203), 03b's
# Phase-1 boundary check, localize/58 and phase1/14 all still read CIS_WINDOW_KB, because those are
# statements about instrument selection and must not move with a sensitivity analysis.
COLOC_REGION_KB <- CIS_WINDOW_KB
INFO_MIN      <- 0.80   # Imputation quality
MAF_MIN       <- 0.01   # Minimum minor allele frequency
F_STAT_MIN    <- 10     # Minimum F-statistic (weak instrument filter)
# NOT the primary instrument-selection rule. Phase 1 selects by DISTANCE-based greedy clumping in a
# 1,000 kb window (08_phase1_extract_hg38.py); no LD reference panel is applied there, and the
# manuscript Methods say so. LD_R2 parameterises the LD re-clumping SENSITIVITY arm only --
# 58_ld_reclump_sensitivity.R passes it to PLINK --clump-r2, and 54 labels its rows with it.
# Changing it does not change which instruments the analysis uses; it changes that sidecar.
LD_R2         <- 0.001  # LD clumping r² threshold (SENSITIVITY ARM ONLY -- see above)
LD_KB         <- 1000   # LD clumping window (kb; also the distance-clump window Phase 1 does use)

# ── Colocalization pass rule (D1 dual rule; SINGLE SOURCE OF TRUTH) ───────────────
# A protein-disease coloc.abf pair PASSES iff PP.H4 >= 0.80 AND PP.H4/PP.H3 >= 5.
# Centralized 2026-07-05 after the rule had drifted across ~30 duplicated call sites
# (a stray H4-only combiner re-admitted MFAP4-CKD). Call is_coloc_pass() everywhere
# instead of re-inlining the predicate. NA in either PP -> FALSE (not NA), so it is
# safe in base/data.table/dplyr subsetting. This is the coloc.abf PP.H4/PP.H3 rule
# ONLY -- do NOT use it for SuSiE concordance (PP.H4_susie) or prior-stability
# (PP.H4_p12_*), which are different criteria on different columns.
COLOC_PP4_MIN  <- 0.80   # min PP.H4
COLOC_H4H3_MIN <- 5      # min PP.H4 / PP.H3 ratio

# "UNEVALUABLE" IS NOT "FAILED" (r37 High). is_coloc_pass() returns FALSE for NA -- correct for its own
# job, which is selecting the 52: a pair with no probabilities must not enter a PASSING set, and FALSE
# keeps it out. But the candidate layer uses the NEGATION of that predicate to mean "this pair was TESTED
# and FAILED", and under that reading a row with missing or malformed ABF probabilities is read as
# evidence of failure when it is the absence of evidence. Such a pair could then be rescued by SuSiE and
# drawn as a dashed edge on the strength of a coloc test that never produced a number.
# So: ask evaluability SEPARATELY, and never by negating the pass rule.
# Frozen data 2026-07-16: coloc_results_v7 has 1,084 rows with 0 NA / 0 non-finite / 0 out-of-range in
# either PP -- so this costs nothing today, and 28's pool inherits that cleanliness.
is_coloc_evaluable <- function(pp_h4, pp_h3) {
  is.finite(pp_h4) & is.finite(pp_h3) &
    pp_h4 >= 0 & pp_h4 <= 1 & pp_h3 >= 0 & pp_h3 <= 1
}

# The candidate layer's "coloc-fail" arm, stated positively: TESTED, and did not pass.
is_coloc_fail <- function(pp_h4, pp_h3) {
  is_coloc_evaluable(pp_h4, pp_h3) & !is_coloc_pass(pp_h4, pp_h3)
}

# ── posterior probabilities must BE probabilities, wherever one is compared to a threshold ────
# `Inf >= 0.80` is TRUE. So is `1.2 >= 0.80`. Every rescue and survival decision in the candidate layer
# is a `>=` against a stored posterior, so a malformed value is not caught by the rule -- it SATISFIES
# the rule, and a computation failure is drawn as a high-posterior edge.
# This is centralised because doing it per-site failed three rounds running (r37, r38, r39): a sweep found
# SIX threshold comparisons (29 x2, 30, 31 loose + strict, 32) against THREE separately-written integrity
# checks that did not line up. Call this wherever a PP column is about to be trusted; do not re-implement.
# `run_status == "ok"` rows only -- a collapsed/skipped row legitimately carries NA.
# Frozen data 2026-07-16: 0 violations across the loose SuSiE output (134 rows), the sensitivity output
# (17), and coloc_results_v7 (1,084).
assert_pp_integrity <- function(d, cols, what, status_col = "run_status", ok_value = "ok") {
  d <- data.table::as.data.table(d)
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop(sprintf("assert_pp_integrity: %s lacks %s -- cannot verify the posteriors it is about to be trusted for",
                 what, paste(miss, collapse = ", ")))
  x <- if (status_col %in% names(d)) d[get(status_col) == ok_value] else d
  if (!nrow(x)) return(invisible(TRUE))
  bad <- Reduce(`|`, lapply(cols, function(.c) {
    v <- x[[.c]]
    !is.finite(v) | v < 0 | v > 1
  }))
  if (any(bad))
    stop(sprintf("assert_pp_integrity: %s has %d row(s) whose %s is not a probability in [0,1]: %s -- `Inf >= 0.80` and `1.2 >= 0.80` are both TRUE, so such a value would SATISFY the rescue/survival rule and be drawn. This is a computation failure, not a high posterior.",
                 what, sum(bad), paste(cols, collapse = "/"),
                 paste(utils::head(x[bad, sprintf("%s~%s(%s)", gene, outcome_name,
                        paste(vapply(cols, function(.c) sprintf("%s=%s", .c, format(.SD[[.c]])), ""), collapse = ", ")),
                        by = seq_len(sum(bad)), .SDcols = cols]$V1, 8), collapse = ", ")))
  invisible(TRUE)
}
is_coloc_pass <- function(pp_h4, pp_h3) {
  # EVALUABLE FIRST. `Inf >= 0.80` is TRUE and `Inf / (pp_h3 + 1e-12) >= 5` is TRUE, so a malformed PP.H4
  # SATISFIED this rule and was counted as one of the 52 (r52 High). r37 protected the FAIL side from the
  # NA case, said explicitly that is_coloc_pass()'s NA-is-FALSE behaviour was "correct for its own job",
  # and never asked the same question about Inf -- the pass side was protected from absent evidence and
  # left open to impossible evidence.
  # `!is.na()` is now subsumed: is_coloc_evaluable() requires finite and in [0,1], which excludes NA.
  # Frozen data 2026-07-16: coloc_results_v7 has 1,084 rows, 0 non-finite / 0 out-of-range, and this rule
  # still selects exactly 52 -- so the SEMANTICS changed and the 52 did not.
  # MULTIPLIED, NOT DIVIDED. `pp_h4 / (pp_h3 + 1e-12) >= 5` excluded part of its own stated boundary:
  # at PP.H4 = 0.80 and PP.H3 = 0.16 the true ratio is exactly 5 and the documented rule passes, but
  # the epsilon put the computed value fractionally below 5 and it failed. The multiplied form is
  # exact and needs no epsilon: at PP.H3 = 0 any evaluable PP.H4 >= 0.80 passes, which is the right
  # reading of an infinite ratio, and at PP.H4 = PP.H3 = 0 the 0.80 arm already fails.
  # Checked against the previous generation's stored tables: 38 non-eGFR, 14 eGFR, 52 combined under
  # BOTH forms -- the boundary is corrected without moving that generation's selected set.
  is_coloc_evaluable(pp_h4, pp_h3) &
    pp_h4 >= COLOC_PP4_MIN &
    pp_h4 >= COLOC_H4H3_MIN * pp_h3
}

# ── The manuscript's PRIMARY FDR set ──────────────────────────────────────────
# BH applied ONCE over the full N_PROTEINS x length(OUTCOMES) grid, with untested cells padded at
# p = 1, over exactly one primary-method row per (gene, outcome). This is the definition the
# manuscript reports (make_draft_manuscript_redesigned.py:84-91) and it uses the same denominator as
# BONFERRONI_ALPHA below -- the two corrections finally agree on what was tested.
#
# IT IS NOT WHAT `mr_significant_fdr.rds` CONTAINS. That artifact is BH WITHIN each method family over
# analysable rows only. The pair is deliberate and documented (cert r7 #2, 2026-07-05; documented at
# make_analysis_reports.R:50-51, PRINTED side by side at :174). Do not "reconcile" them by changing
# cascade/01_mr_screen.R (the `fdr_sig` line) -- that collapses [2] onto [1] and makes
# make_analysis_reports.R:174 print the same number twice while calling one a sensitivity analysis.#
# THREE FDR SURFACES EXIST, and two of them are both called "234". Name them, never merge them
# (journal-editor C-2 + objective review r16, 2026-07-16):
#   ⚠ THE COUNTS BELOW ARE FROZEN-DATA OBSERVATIONS, NOT DEFINITIONS. They were measured on the
#   1,682-gene generation (grid 1,682 x 8 = 13,456). The GRCh38 rebuild makes the grid
#   1,948 x 8 = 15,584 and every count below has to be re-measured; do not carry 221/234 forward,
#   and do not let a reporting layer reimplement "13,456" as a constant. What survives the rebuild
#   is the DEFINITION of each surface, not its size.
#   [1] PRIMARY, 221   -- BH ONCE over N_PROTEINS_CANONICAL x length(OUTCOMES), untested cells
#                         padded p=1, one primary-method row per pair. primary_fdr_set(). The headline.
#   [2] SAVED BY-METHOD ARTIFACT, 234 -- mr_significant_fdr.rds: BH WITHIN each method family, ANALYSABLE
#                         BEFORE the primary dedup. What 01_mr_screen.R's `fdr_sig` line writes.
#   [3] MANUSCRIPT SENSITIVITY, 234 -- what the generator computes fresh and PRINTS: POOLED BH over the
#                         11,534 analysable rows, AFTER the primary dedup. A DIFFERENT PROCEDURE from [2].
# [2] and [3] select the IDENTICAL SET on the frozen data -- verified 2026-07-16 by setequal(), not merely
# by both being 234. That agreement is an OBSERVED FACT ABOUT THIS DATA, not a guarantee: it holds because
# 0 pairs carry both Wald and IVW (so the dedup is a no-op) and no pair falls in the 1.9% gap between the
# 11,322-row Wald family and the 11,534-row pooled denominator. Different data can separate them.
# The manuscript is internally consistent -- it describes [3] and computes [3]. The conflation to avoid is
# calling the ARTIFACT [2] "the sensitivity value the manuscript reports": it is a different computation
# that currently agrees.
#
# This function exists because of two real defects (2026-07-15):
#   1. Consumers read the saved by-method artifact [2] while the rule they state is the primary one [1].
#      FIXED here: phase3/28_candidate_edge_enumerate.R, phase3/32_dashed_finalize.R.
#      NOT FIXED, deliberately, and tracked in REBUILD_LEDGER §2.1c: cascade/14_assign_tiers.R, which
#      feeds the artifact to assign_evidence_tiers(fdr_rds=) in R/theme_paper3.R. There it only relabels
#      mr_sig_class "Suggestive" -> "FDR" (theme_paper3.R:408); the tier criterion crit_bonferroni
#      (:455) tests == "Bonferroni" and never consults FDR. Measured on the frozen 52 (2026-07-16):
#      30 Bonferroni / 19 FDR / 3 Suggestive, and 19/19 of the FDR-labelled pairs are FDR-significant
#      under BOTH definitions -- 0 would move. So the defect is real but provably inert, and fixing it
#      means changing assign_evidence_tiers()'s signature in theme_paper3.R, which carries a byte
#      contract. That is a PI call, not something this function should decide by being called.
#   2. The primary rule had NO single source. This function is now the source for the two phase3
#      consumers above -- but the reporting layer STILL re-implements it (make_analysis_reports.R:52,
#      make_draft_manuscript_redesigned.py:84, fig1/fig2/volcano/discovery_forest). Do not read
#      "single source" as done: it is the owner of the definition, not yet its only implementation.
#      Consolidation is tracked in REBUILD_LEDGER §2.1c and belongs to the reporting-layer audit.
# Measured difference on the frozen data: [2] has 234 rows, this one ([1]) 221 -- 13 pairs,
# all with p in 0.00083-0.00101 (the BH boundary, which is exactly where a larger denominator bites).
# None of the 13 is among the 52 colocalization-supported pairs, so tier labels do not move.
# ── outcome_name <-> outcome_id, asserted against the SINGLE SOURCE (OUTCOMES) ────────────────
# Checking that a table's mapping is internally 1:1 is CONSISTENCY, not TRUTH: two artifacts stale in the
# SAME way are perfectly 1:1 and perfectly wrong, and every engine in phase3 dispatches on outcome_id
# (objective review r21/r22, 2026-07-16). Demonstrated: a constructed T2D = FinnGen_R9_T2D passes a 1:1
# check and is caught here. Call this at EVERY artifact read boundary -- one helper, not a per-site
# re-implementation, because per-site copies are what let this diverge in the first place.
# `required = TRUE` by default, and that default is the point. A first version silently returned when the
# columns were absent, described as a no-op for tables that "carry only outcome_name". But EVERY current
# caller passes a table that has both -- so the no-op was never exercised, and all it could ever do was go
# quiet at the exact moment the table stopped looking like what the caller thought it was. That is the
# same shape this guard exists to catch. A caller that genuinely wants the permissive form must now say so.
assert_outcome_truth <- function(d, what, required = TRUE) {
  d <- data.table::as.data.table(d)
  miss <- setdiff(c("outcome_name", "outcome_id"), names(d))
  if (length(miss)) {
    if (!required) return(invisible(FALSE))
    stop(sprintf("assert_outcome_truth: %s lacks %s -- cannot verify the outcome mapping it is being trusted for. Pass required=FALSE only if this table genuinely has no outcome_id and nothing downstream dispatches on one.",
                 what, paste(miss, collapse = ", ")))
  }
  truth <- data.table::data.table(outcome_name = names(OUTCOMES), outcome_id = unname(unlist(OUTCOMES)))
  got <- unique(d[, .(outcome_name, outcome_id)])
  bad <- got[!truth, on = .(outcome_name, outcome_id)]
  if (nrow(bad))
    stop(sprintf("assert_outcome_truth: %s carries %d (outcome_name, outcome_id) pair(s) absent from config's OUTCOMES: %s -- the artifact predates a config change, or config was re-pointed. Engines dispatch on outcome_id; refusing.",
                 what, nrow(bad), paste(bad[, paste0(outcome_name, "=", outcome_id)], collapse = ", ")))
  invisible(TRUE)
}

# ── which SuSiE engine a row MUST use, derived from its outcome ───────────────────────────────
# `susie_engine` is a stored label, and rounds 23/26/27 each found the same failure: the table that was
# checked was not the table that decided. Checking outcome_id and then dispatching on an UNCHECKED engine
# label is that shape again -- a pool row with a valid CAD outcome_id but susie_engine="egfr" passes
# assert_outcome_truth() and is then run through the eGFR engine (r27 High).
# The engine is not independent information: it is a function of the outcome. Derive it, never trust it.
# NA is legitimate ONLY for non-testable rows -- 28 deliberately leaves MHC / undetermined-MHC /
# no-coords rows engine-less so 29 never SuSiE-runs them. This function therefore SKIPS engine-less rows,
# and the caller that knows whether a row is testable must reject NA itself.
# THAT DELEGATION HAS ALREADY FAILED ONCE (r28): 32 accepted "both sides NA" as agreement, so a pair of
# stale pre-stamp rows with matching outcome_id and no engine at all could be finalized carrying no engine
# provenance. If you call this on a table whose rows are testable BY CONSTRUCTION -- anything selected as
# `dashed_edge`, anything 29/31 emitted -- you must also reject NA. Frozen data 2026-07-16: NA engine
# appears on exactly the 17 `not_evaluable` audit rows and nowhere else (dashed 17/0, excluded 117/0,
# sensitivity 17/0).
expected_susie_engine <- function(outcome_id) {
  data.table::fifelse(outcome_id == OUTCOMES[["eGFR"]], "egfr", "nonegfr", na = NA_character_)
}

# Assert stored engine == derived engine, for rows that HAVE a stored engine.
assert_engine_truth <- function(d, what, required = TRUE) {
  d <- data.table::as.data.table(d)
  miss <- setdiff(c("outcome_id", "susie_engine"), names(d))
  if (length(miss)) {
    if (!required) return(invisible(FALSE))
    stop(sprintf("assert_engine_truth: %s lacks %s -- cannot verify the engine label it is being trusted for. Pass required=FALSE only if nothing downstream dispatches on an engine.",
                 what, paste(miss, collapse = ", ")))
  }
  x <- d[!is.na(susie_engine) & susie_engine != ""]
  if (!nrow(x)) return(invisible(TRUE))
  bad <- x[susie_engine != expected_susie_engine(outcome_id)]
  if (nrow(bad))
    stop(sprintf("assert_engine_truth: %s carries %d row(s) whose susie_engine contradicts its outcome_id: %s -- the engine label is a function of the outcome, so a disagreement means one of them is stale. Refusing.",
                 what, nrow(bad),
                 paste(utils::head(bad[, sprintf("%s~%s(id=%s, engine=%s, expected=%s)", gene, outcome_name,
                                                 outcome_id, susie_engine, expected_susie_engine(outcome_id))], 8),
                       collapse = ", ")))
  invisible(TRUE)
}

# ── runnability, derived all the way down to the underlying facts ─────────────────────────────
# 28's rule is `is_testable := !is.na(is_mhc) & is_mhc == FALSE & has_coords == TRUE`. Rounds 29-32
# replaced trust in `is_testable`/`susie_engine` with a derivation from `is_mhc` and `has_coords` -- and
# r33 pointed out those are STORED LABELS TOO. A derivation whose inputs are labels is still trust, one
# level down; the chain has to bottom out in facts. So:
#   has_coords  <-  the coordinate columns themselves
#   is_mhc      <-  is_mhc_gene(), config's single source for the MHC rule
# Frozen data 2026-07-16: stored has_coords matches the coordinate columns on all 151 rows; stored is_mhc
# matches is_mhc_gene() on all 151, with 0 undetermined. Both guards therefore cost nothing today.
# NOTE is_mhc_gene() may legitimately return a non-scalar / error for an unknown gene -- that is
# "undetermined MHC", which 28 treats as NOT testable (it withholds the engine). Preserve that: NA here
# means undetermined, and undetermined is not runnable.
# Read a stored logical label, and REFUSE anything that is not one.
# Two failure modes, and the second one bit AFTER the first was fixed:
#   `as.logical("bogus")` is NA; `NA != TRUE` is NA; a data.table filter does not select NA; and
#   `is.na("bogus")` is FALSE -- so every `as.logical(label) != derived | is.na(label)` check in this
#   layer silently ACCEPTED a label it could not parse (r47).
#   Then: `as.logical(2)` is TRUE. So is `as.logical(-1)` and `as.logical(0.5)`. A first version of this
#   guard only rejected values that coerce to NA, so numeric garbage sailed through AS a valid TRUE
#   whenever the derived value happened to be TRUE (r49). Rejecting "unreadable" is not the same as
#   accepting only "a logical", and R's coercion is lenient in exactly the direction that hides this.
# So: parse EXPLICITLY against a closed vocabulary; do not delegate to as.logical() at all.
# NA passes through as NA -- "legitimately missing" is the caller's to allow (is_mhc = undetermined MHC)
# or to reject; only "present but not a logical" stops here.
# This lived as a local helper inside derive_runnability() when first written, which is exactly why 32 did
# not have it and repeated the r47 bug two rounds later (r48). Reusable things go here.
LOGICAL_TRUE_TOKENS  <- c("TRUE", "true", "T", "True")
LOGICAL_FALSE_TOKENS <- c("FALSE", "false", "F", "False")
as_logical_strict <- function(x, what, nm) {
  if (is.logical(x)) return(x)                       # already one; NA stays NA
  bad <- rep(FALSE, length(x))
  if (is.character(x)) {
    ok <- is.na(x) | x %in% c(LOGICAL_TRUE_TOKENS, LOGICAL_FALSE_TOKENS)
    bad <- !ok
    v <- ifelse(is.na(x), NA, x %in% LOGICAL_TRUE_TOKENS)
  } else if (is.numeric(x)) {
    # 0/1 only. `as.logical(2)` is TRUE and `as.logical(-1)` is TRUE; a flag column holding 2 is a defect,
    # not a TRUE, and reading it as one is how a malformed label becomes agreement.
    ok <- is.na(x) | x %in% c(0, 1)
    bad <- !ok
    v <- ifelse(is.na(x), NA, x == 1)
  } else {
    stop(sprintf("as_logical_strict: %s's %s is a %s -- a logical label must be logical, character or 0/1 numeric",
                 what, nm, paste(class(x), collapse = "/")))
  }
  if (any(bad))
    stop(sprintf("as_logical_strict: %s has %d row(s) whose %s is present but is not a logical value (e.g. %s) -- accepted tokens are %s / %s or 0 / 1. R would coerce some of these silently (as.logical(2) is TRUE, as.logical(\"bogus\") is NA) and the result would be read as agreement.",
                 what, sum(bad), nm, paste(utils::head(unique(as.character(x[bad])), 3), collapse = ", "),
                 paste(LOGICAL_TRUE_TOKENS, collapse = "/"), paste(LOGICAL_FALSE_TOKENS, collapse = "/")))
  as.logical(v)
}

derive_runnability <- function(d, what) {
  d <- data.table::as.data.table(d)
  # Registered HERE, before any temporary exists and before any stop can fire: this function mutates the
  # CALLER's data.table by reference, so a stop between creating a temporary and returning would leave it
  # behind in their table (r48 Low). suppressWarnings: the columns may not exist yet when this runs.
  on.exit(suppressWarnings(d[, c(".hc_lab", ".mhc_lab", ".tst_lab") := NULL]), add = TRUE)
  # `is_testable` is deliberately NOT required: this helper's job is to derive runnability from FACTS,
  # so it must work on a table that has never carried the label. Where the label IS present it is checked
  # (below), but its absence is not an error. A first version omitted it from `need` and then used it
  # unconditionally anyway -- contract and implementation disagreeing (r35 Medium).
  need <- c("gene", "is_mhc", "has_coords", "gene_chrom", "gene_start", "gene_end")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(sprintf("derive_runnability: %s lacks %s -- cannot derive runnability from facts, and deriving it from the stored labels instead is what this function exists to stop",
                 what, paste(miss, collapse = ", ")))
  d[, .has_coords_rd := !(is.na(gene_chrom) | is.na(gene_start) | is.na(gene_end))]
  d[, .is_mhc_rd := vapply(gene, function(g) {
        m <- tryCatch(is_mhc_gene(g), error = function(e) NA)
        if (length(m) != 1L) NA else as.logical(m)
      }, TRUE)]
  # COERCE ONCE, AND REJECT WHAT WOULD NOT PARSE. `as.logical("bogus")` is NA, `NA != TRUE` is NA, and a
  # data.table filter does not select NA -- so a label this guard could not even READ slipped through and
  # the guard reported agreement about it. `is.na(original)` does not save it either: "bogus" is not NA
  # (r47 Medium). Coerce into a temporary, and treat "could not parse" as its own failure, distinct from
  # "legitimately missing".
  d[, .hc_lab := as_logical_strict(has_coords, what, "has_coords")]
  bad_c <- d[is.na(.hc_lab) | .hc_lab != .has_coords_rd]
  if (nrow(bad_c))
    stop(sprintf("derive_runnability: %s has %d row(s) whose stored has_coords contradicts its coordinate columns: %s -- a stale FALSE makes a real edge vanish untested; a stale TRUE sends missing coordinates into the engine",
                 what, nrow(bad_c), paste(utils::head(bad_c[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
  # NA-AWARE. A first version used `!is.na(x) & !is.na(y) & x != y`, which reads a ONE-SIDED NA as
  # agreement -- so a stored is_mhc = NA against a derived TRUE/FALSE (or the reverse) passed silently
  # (r34 High). NA here means "undetermined MHC", which 28 treats as NOT runnable; a stored NA where the
  # rule now has an answer, or a stored answer where the rule now has none, is exactly the drift this
  # checks for. Same predicate shape as 31/32 use for outcome_id/engine.
  .ne <- function(x, y) xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
  # unreadable != legitimately-NA: an NA is_mhc means "undetermined MHC" and is allowed through to the
  # .ne() comparison below; a "bogus" is_mhc is not.
  d[, .mhc_lab := as_logical_strict(is_mhc, what, "is_mhc")]
  bad_m <- d[.ne(.mhc_lab, .is_mhc_rd)]
  if (nrow(bad_m))
    stop(sprintf("derive_runnability: %s has %d row(s) whose stored is_mhc contradicts is_mhc_gene() (a one-sided NA counts): %s -- config's MHC rule is the single source; a stale label here decides what gets SuSiE-run",
                 what, nrow(bad_m), paste(utils::head(bad_m[, sprintf("%s~%s(stored=%s, derived=%s)", gene,
                                            outcome_name, is_mhc, .is_mhc_rd)], 8), collapse = ", ")))
  d[, .exp_testable := (!is.na(.is_mhc_rd) & .is_mhc_rd == FALSE & .has_coords_rd == TRUE)]
  if ("is_testable" %in% names(d)) {
    d[, .tst_lab := as_logical_strict(is_testable, what, "is_testable")]
    bad_t <- d[is.na(.tst_lab) | .tst_lab != .exp_testable]
    if (nrow(bad_t))
      stop(sprintf("derive_runnability: %s has %d row(s) whose stored is_testable contradicts the fact-derived rule: %s",
                   what, nrow(bad_t), paste(utils::head(bad_t[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
  }
  d[]
}

# coloc-fail, derived from the raw ABF probabilities rather than the stored `coloc_status` label.
assert_coloc_status_truth <- function(d, what) {
  d <- data.table::as.data.table(d)
  need <- c("coloc_status", "PP.H4_abf", "PP.H3_abf")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(sprintf("assert_coloc_status_truth: %s lacks %s -- cannot check the label against the rule it encodes",
                 what, paste(miss, collapse = ", ")))
  # UNEVALUABLE FIRST. This helper compared the label to `!is_coloc_pass(...)`, which reads an NA-ABF row
  # as "failed" -- the very conflation r37 introduced is_coloc_fail() to remove. The call sites were fixed
  # and the shared helper was not (r38 High), which is the whole reason the rule lives in one place.
  une <- d[!is_coloc_evaluable(PP.H4_abf, PP.H3_abf)]
  if (nrow(une))
    stop(sprintf("assert_coloc_status_truth: %s has %d row(s) with missing or out-of-range ABF probabilities: %s -- these were never evaluated by coloc, and no coloc_status label can be true of them",
                 what, nrow(une), paste(utils::head(une[, sprintf("%s~%s(H4=%s, H3=%s)", gene, outcome_name,
                                                                  format(PP.H4_abf), format(PP.H3_abf))], 8),
                                        collapse = ", ")))
  # An NA label makes the comparison NA, and a data.table filter does not select NA -- so a MISSING
  # coloc_status slipped through and this helper reported agreement about a label that was not there
  # (r41 Medium). Reject the missing/unknown label first, then compare something that cannot be NA.
  ok_lab <- c("coloc_fail", "coloc_pass")
  bad_lab <- d[is.na(coloc_status) | !coloc_status %in% ok_lab]
  if (nrow(bad_lab))
    stop(sprintf("assert_coloc_status_truth: %s has %d row(s) with a missing or unrecognised coloc_status (expected one of %s): %s -- an absent label cannot agree with the ABF rule, and an NA comparison is silently dropped by the filter below",
                 what, nrow(bad_lab), paste(ok_lab, collapse = "/"),
                 paste(utils::head(bad_lab[, sprintf("%s~%s(status=%s)", gene, outcome_name,
                                                     as.character(coloc_status))], 8), collapse = ", ")))
  bad <- d[(coloc_status == "coloc_fail") != is_coloc_fail(PP.H4_abf, PP.H3_abf)]
  if (nrow(bad))
    stop(sprintf("assert_coloc_status_truth: %s has %d row(s) whose stored coloc_status contradicts is_coloc_fail(PP.H4_abf, PP.H3_abf): %s -- a coloc-PASSING pair is a SOLID edge and must never be frozen as an exploratory dashed candidate",
                 what, nrow(bad), paste(utils::head(bad[, sprintf("%s~%s(status=%s, H4=%s, H3=%s)", gene, outcome_name,
                                                                  coloc_status, format(PP.H4_abf), format(PP.H3_abf))], 8),
                                        collapse = ", ")))
  invisible(TRUE)
}


# ── a status label's PREMISE, both directions ────────────────────────────────────────────────
# assert_pp_integrity() answers "an `ok` row must carry usable probabilities". It never asked the
# inverse -- "a NON-ok row must carry NONE" -- so a stale `skip_benign` row with a real PP.H4_susie was
# frozen as `no_susie_rescue`, and a stale `collapsed` sensitivity row with a real strict PP was
# finalized as `dropped:collapsed`. Both are the status label being TRUSTED rather than its premise being
# CHECKED (r42 High): a row that reports no usable result must not be sitting on one.
# Frozen data 2026-07-16: loose has 5 skip_benign rows, strict has 4 collapsed -- 0 of the 9 carry a
# non-NA posterior. Free today.
assert_status_premise <- function(d, cols, what, status_col = "run_status", ok_value = "ok") {
  d <- data.table::as.data.table(d)
  miss <- setdiff(c(cols, status_col), names(d))
  if (length(miss))
    stop(sprintf("assert_status_premise: %s lacks %s", what, paste(miss, collapse = ", ")))
  # A row that exists must say what happened to it. NA is not a status.
  no_stat <- d[is.na(get(status_col)) | !nzchar(trimws(as.character(get(status_col))))]
  if (nrow(no_stat))
    stop(sprintf("assert_status_premise: %s has %d row(s) with no %s at all: %s -- a row that exists must report an outcome; an NA status slips past every `!is.na(status) & status != ...` check written downstream",
                 what, nrow(no_stat), status_col,
                 paste(utils::head(no_stat[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
  non_ok <- d[get(status_col) != ok_value]
  if (!nrow(non_ok)) return(invisible(TRUE))
  carries <- Reduce(`|`, lapply(cols, function(.c) !is.na(non_ok[[.c]])))
  if (any(carries))
    stop(sprintf("assert_status_premise: %s has %d row(s) whose %s is not '%s' yet which carry a %s: %s -- the status says there is no usable result and the row is holding one, so one of the two is stale",
                 what, sum(carries), status_col, ok_value, paste(cols, collapse = "/"),
                 paste(utils::head(non_ok[carries, sprintf("%s~%s(%s=%s)", gene, outcome_name,
                        status_col, get(status_col))], 8), collapse = ", ")))
  invisible(TRUE)
}

# ── does a captured warning actually say the fit did not converge? ───────────────────────────
# `held:non_converged` is a CLAIM about why an edge was withheld, and the LOCKED caption repeats it:
# "crossed the SuSiE threshold only with a non-convergent fit". It was assigned on `n_warn > 0` alone --
# any warning at all -- so an edge warned for some other reason would be held correctly and LABELLED
# falsely, and that false reason would reach the caption (r51 High).
# NOT checkable via `ld_warn`: that is `grepl("LD|correlation|inconsist|kriging|check_ld|refin", ...)`,
# and on the three held edges it matches the word "LD" inside susieR's REMEDIATION ADVICE ("Please check
# consistency between summary statistics and LD matrix"), not an LD warning as such. `ld_warn` means "the
# message mentions LD", which is a different claim from either "it did not converge" or "the LD was
# inconsistent". Ask the convergence text directly.
# Frozen data 2026-07-16: all 3 held edges carry exactly "IBSS algorithm did not converge in 100
# iterations!" -- so the label is true today, and this makes it checkable rather than assumed.
SUSIE_NONCONVERGENCE_PATTERNS <- c("did not converge", "IBSS algorithm did not converge",
                                   "failed to converge", "not converged")
is_nonconvergence_warning <- function(msg) {
  !is.na(msg) & Reduce(`|`, lapply(SUSIE_NONCONVERGENCE_PATTERNS,
                                   function(p) grepl(p, msg, fixed = TRUE, useBytes = FALSE)))
}


# ── did a SuSiE fit actually run for this row? ───────────────────────────────────────────────
# ONE definition, because four guards needed it and three of them got it wrong in three different ways.
# The question every convergence guard has to ask is not "is the column there" but "does this row claim not
# to know something it demonstrably knows". A row whose fit ran and whose convergence is NA is FALSE, not
# missing -- that is what a half-stale artifact looks like once the columns exist.
# `PP.H4_susie` ALONE IS THE WRONG TEST, and this is the trap: a post-fit collapse (coloc.susie failed or
# returned an empty summary AFTER both susie_rss calls succeeded) legitimately has PP.H4_susie = NA while
# S1/S2 existed and their convergence is known. Testing the posterior misses exactly the rows the
# convergence recording was extended to cover.
# n_signals_exp / n_signals_out are the direct evidence: they are length(S1$sets$cs) / length(S2$sets$cs),
# so non-NA means the fit objects existed. The others are OR-ed because being conservative here means
# "more rows must carry honest convergence", which is the safe direction to be wrong in.
susie_fit_ran <- function(d, posterior = "PP.H4_susie", n_exp = "n_signals_exp", n_out = "n_signals_out",
                          n_ld = "n_snps_ld") {
  d <- data.table::as.data.table(d)
  ev <- lapply(c(posterior, n_exp, n_out, n_ld), function(cn) if (cn %in% names(d)) !is.na(d[[cn]]) else rep(FALSE, nrow(d)))
  Reduce(`|`, ev)
}

# Reject rows that claim "convergence unknown" about a fit that ran. `what` names the artifact for the message.
assert_convergence_honest <- function(d, what, conv = c("converged_exp", "converged_out"), ...) {
  d <- data.table::as.data.table(d)
  miss <- setdiff(conv, names(d))
  if (length(miss))
    stop(sprintf("%s lacks %s -- it predates the per-fit convergence recording. Re-run the engine that writes it; a presence-only check upstream would let a half-stale table through.",
                 what, paste(miss, collapse = ", ")))
  ran <- susie_fit_ran(d, ...)
  bad <- ran & (is.na(d[[conv[1]]]) | is.na(d[[conv[2]]]))
  if (any(bad)) {
    lab <- if (all(c("gene", "outcome_name") %in% names(d))) d[bad, paste0(gene, "~", outcome_name)] else which(bad)
    stop(sprintf("%s has %d row(s) whose fit demonstrably ran but whose convergence is NA (e.g. %s) -- 'unknown' is false there, not missing. This is a half-stale artifact from a partial re-run. Re-run the engine from scratch (delete its checkpoint).",
                 what, sum(bad), paste(utils::head(lab, 6), collapse = ", ")))
  }
  invisible(TRUE)
}

FDR_ALPHA_PRIMARY_METHODS <- c("Wald ratio", "Inverse variance weighted")
primary_fdr_set <- function(mr_all, n_proteins = NULL, alpha = NULL, protein_universe = NULL) {
  stopifnot("primary_fdr_set: mr_all must be a data.frame/data.table" = is.data.frame(mr_all))
  d <- data.table::as.data.table(mr_all)
  need <- c("gene", "outcome_name", "method", "pval")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(sprintf("primary_fdr_set: mr_all lacks %s", paste(miss, collapse = ", ")))
  a  <- if (is.null(alpha)) FDR_ALPHA else alpha
  # `alpha` was the one parameter this function did not validate while validating every other one
  # aggressively (r9 Medium). A vector, NA, Inf, a character, or a value outside [0,1] reaches the final
  # `fdr_q < a` and yields an empty, nonsensical, or silently-everything set. Latent today (both callers
  # use the default) -- which is exactly when it is cheap to close.
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 0 || a > 1)
    stop(sprintf("primary_fdr_set: alpha must be a single finite number in (0, 1], got %s",
                 paste(utils::capture.output(str(a)), collapse = " ")))
  # Default to N_PROTEINS_CANONICAL, NOT the derived N_PROTEINS. Two reasons, both load-bearing:
  #  1. BONFERRONI_ALPHA (below) is `bonferroni_alpha(N_PROTEINS_CANONICAL)`. The entire argument for
  #     this function is that the FDR and Bonferroni corrections should agree on what was tested --
  #     defaulting to a DIFFERENT count than Bonferroni uses would reintroduce the split it fixes.
  #  2. The derived N_PROTEINS is `tryCatch(..., error = function(e) NA_integer_)`: on any read failure
  #     it is NA, and `grid <- NA * 8` is NA, so `if (grid < n_an)` would raise "missing value where
  #     TRUE/FALSE needed" instead of the guard message -- and a caller that pre-set its own N_PROTEINS
  #     before sourcing config would silently change the denominator.
  np <- if (is.null(n_proteins)) N_PROTEINS_CANONICAL else n_proteins
  # is.finite() is doing real work here: `is.na(Inf)` is FALSE and `Inf != as.integer(Inf)` warns and
  # returns NA, so without it an infinite denominator slips past into `grid` and the fail-loud path never
  # fires. NA, NaN and Inf must all land on the same stop.
  if (!is.numeric(np) || length(np) != 1L || !is.finite(np) || np <= 0 || np != as.integer(np))
    stop(sprintf("primary_fdr_set: n_proteins must be a single positive whole number, got %s -- the full-grid denominator cannot be established and BH would silently fall back to an analysable-only correction",
                 paste(utils::capture.output(str(np)), collapse = " ")))
  np <- as.integer(np)
  # `d[...]` returns a NEW data.table, so the := calls below cannot write back into the caller's
  # mr_all by reference. Verified by objective review 2026-07-15; do not replace with `d[, ...]`.
  pm <- d[method %in% FDR_ALPHA_PRIMARY_METHODS]
  if (!nrow(pm))
    stop("primary_fdr_set: no Wald ratio / IVW rows -- an empty set would report agreement vacuously")
  # EXACTLY ONE primary row per pair, enforced -- not "prefer Wald over IVW".
  # The previous shape ranked Wald=1 / IVW=2, rejected ties only at the BEST rank, then kept the first
  # after sorting. A pair carrying BOTH Wald and IVW therefore had n_best==1, passed the check, and had
  # its IVW estimate silently discarded by sort order -- an estimate nobody chose, from a pair whose
  # instrument set is self-contradictory (Wald is the single-SNP ratio, IVW the multi-SNP meta; a pair
  # should carry one XOR the other). Duplicate rows also double-count inside the BH family, moving both
  # the rank ordering and the effective denominator. Objective review 2026-07-16 (r5 High).
  # Measured the same day on mr_all_results.rds: 11,534 primary rows, 0 pairs with >1 -- so failing loud
  # costs nothing today and closes the hole before data change opens it.
  # phase3/28 has a weaker sibling guard (it stops only when duplicate rows DISAGREE on b/pval); this one
  # is deliberately stricter, because identical duplicates are harmless to a point estimate but not to BH.
  dup <- pm[, .N, by = .(gene, outcome_name)][N > 1L]
  if (nrow(dup))
    stop(sprintf("primary_fdr_set: %d (gene,outcome) pair(s) carry >1 primary-method row (Wald and/or IVW) -- the primary row is ambiguous and BH would double-count them: %s",
                 nrow(dup), paste(utils::head(dup[, paste(gene, outcome_name)], 12L), collapse = ", ")))
  # The BH family and the grid denominator must describe the same universe. Without these, a stale or
  # non-CKM row joins the numerator while `grid` stays N_PROTEINS_CANONICAL x length(OUTCOMES) -- the
  # correction would then be computed over a universe that does not exist (r5 Medium #3).
  if (anyNA(pm$gene) || anyNA(pm$outcome_name))
    stop("primary_fdr_set: primary rows carry NA gene/outcome_name -- such a row cannot be placed on the grid it is being corrected against")
  oob <- setdiff(unique(as.character(pm$outcome_name)), names(OUTCOMES))
  if (length(oob))
    stop(sprintf("primary_fdr_set: outcome_name(s) outside the CKM universe (%s): %s -- they would enter the BH family while the denominator stays %d x %d",
                 paste(names(OUTCOMES), collapse = "/"), paste(oob, collapse = ", "), np, length(OUTCOMES)))
  if (data.table::uniqueN(pm$gene) > np)
    stop(sprintf("primary_fdr_set: %d distinct genes exceed n_proteins (%d) -- the grid denominator understates what was actually tested, so BH would be anti-conservative",
                 data.table::uniqueN(pm$gene), np))
  # MEMBERSHIP, not just count. The count check above only catches "too many"; a stale or non-CKM gene
  # that REPLACES a missing canonical one keeps the count under the cap and sails through, putting a row
  # in the BH numerator that has no cell on the denominator's grid (r6 High, 2026-07-16).
  # Defaults to PROTEIN_UNIVERSE (derived from INSTRUMENTS_FILE at config load). If it cannot be
  # established we STOP: falling back to the count check would be silently reverting to the exact
  # weaker test this guard exists to replace.
  u <- if (is.null(protein_universe)) PROTEIN_UNIVERSE else protein_universe
  if (is.null(u) || !length(u))
    stop(sprintf("primary_fdr_set: the protein universe is unavailable (INSTRUMENTS_FILE = %s) -- cannot verify the BH numerator lives on the %d x %d grid; pass protein_universe= explicitly or restore the instrument set",
                 INSTRUMENTS_FILE, np, length(OUTCOMES)))
  u <- unique(as.character(u))
  if (anyNA(u))
    stop("primary_fdr_set: the protein universe contains NA -- an NA gene would match nothing and silently widen what counts as 'in universe'")
  # Membership is only half the claim. The universe and the DENOMINATOR must be the same universe:
  # if INSTRUMENTS_FILE drifts to 1,683 genes, an out-of-contract gene is "in universe" while `grid`
  # stays np x length(OUTCOMES) -- and the count guard above misses it unless the analysed rows exceed
  # np (objective review r7, 2026-07-16). config only WARNS on that drift (see N_PROTEINS below); here
  # it must stop, because BH is what consumes the mismatch.
  if (length(u) != np)
    stop(sprintf("primary_fdr_set: the protein universe has %d genes but the BH denominator uses n_proteins=%d -- membership would be checked against one universe and corrected against another (INSTRUMENTS_FILE = %s)",
                 length(u), np, INSTRUMENTS_FILE))
  oog <- setdiff(unique(as.character(pm$gene)), u)
  if (length(oog))
    stop(sprintf("primary_fdr_set: %d gene(s) outside the %d-protein universe: %s -- they would enter the BH family while the denominator stays %d x %d",
                 length(oog), length(u), paste(utils::head(oog, 12L), collapse = ", "), np, length(OUTCOMES)))
  if (!all(is.finite(pm$pval) & pm$pval >= 0 & pm$pval <= 1))
    stop("primary_fdr_set: primary rows carry a non-finite or out-of-range pval -- BH would be silently wrong")
  # The FDR set is what 28 builds the candidate pool from and what 32 re-asserts rule (1) against, and it
  # carries outcome_id downstream into engine dispatch. Callers were validating coloc and Bonferroni
  # against OUTCOMES but never `mr` itself -- the one table this set is actually made of (r22 High).
  # Belongs HERE, not at each call site: this function is the boundary the FDR set leaves through.
  assert_outcome_truth(pm, "the primary-method rows of mr_all (primary_fdr_set input)")
  n_an <- nrow(pm); grid <- np * length(OUTCOMES)
  # If the grid were smaller than the analysable count the padding would vanish and this would silently
  # become the analysable-only BH it exists to replace -- i.e. it would no-op without saying so.
  if (grid < n_an)
    stop(sprintf("primary_fdr_set: grid (%d x %d = %d) < analysable primary tests (%d) -- the p=1 padding would vanish and this would silently revert to an analysable-only correction",
                 np, length(OUTCOMES), grid, n_an))
  pm[, fdr_q := stats::p.adjust(c(pval, rep(1, grid - n_an)), method = "BH")[seq_len(n_an)]]
  pm[fdr_q < a]
}

# ── Outcome GWAS IDs (IEU OpenGWAS, GRCh37) ────────────────────────────────────
OUTCOMES <- c(
  "CAD"    = "FinnGen_R10_I9_CHD",
  "Stroke" = "MEGASTROKE_EUR_GCST006906",
  "HF"     = "FinnGen_R10_I9_HEARTFAIL",
  "T2D"    = "FinnGen_R10_T2D",
  "eGFR"   = "Wuttke2019_eGFR_EA",
  "CKD"    = "Wuttke2019_CKD_EA",
  "NAFLD"  = "FinnGen_R10_NAFLD",
  "AF"     = "FinnGen_R10_I9_AF"
)

# ── 2026-06-30 UKB-free overhaul: authoritative outcome N/ncase (SINGLE SOURCE OF TRUTH) ──
# Data = FinnGen R10 (CAD/HF/T2D/NAFLD/AF) / Wuttke 2019 EA (eGFR/CKD) / MEGASTROKE EUR (Stroke).
# Keys above RENAMED to trait-based (2026-06-30, grep-gated): no more finn-b-*/ieu-a-1102/ebi-a-GCST006414
# — those OpenGWAS accessions re-fetch UKB-containing (Nielsen AF) / Mixed-ancestry (Pattaro CKD) data. results regen at rerun.
# Prefer reading these over hardcoding ncase/N; the stopifnot() catches drift/typo at source-load time.
OUTCOME_N_TOTAL <- c(CAD=412181L, HF=412181L, T2D=400197L, NAFLD=412181L,
                     eGFR=567460L, CKD=480698L, Stroke=446696L, AF=261395L)
OUTCOME_N_CASE  <- c(CAD=46959L, HF=29672L, T2D=65085L, NAFLD=2568L,
                     eGFR=NA_integer_, CKD=41395L, Stroke=40585L, AF=50743L)
OUTCOME_UKB     <- c(CAD=FALSE, HF=FALSE, T2D=FALSE, NAFLD=FALSE,
                     eGFR=FALSE, CKD=FALSE, Stroke=FALSE, AF=FALSE)  # 8 disease outcomes all UKB-free
# CHECKED AGAINST THE RELEASE MANIFEST, not against four hand-picked cells. The old assertion tested
# T2D total, CAD cases, CKD cases and AF total -- so a wrong HF total or NAFLD case count passed the
# very guard that says it catches drift at source time. The five FinnGen outcomes are now DERIVED
# from data/outcome_gwas/finngen_R10_manifest.tsv (num_cases + num_controls), which is the release's
# own record and the thing the literals above were transcribed from; a transcription error can no
# longer survive. Wuttke and MEGASTROKE publish no equivalent manifest, so those six stay literal and
# are asserted explicitly below with their source named.
local({
  .fgm <- file.path(OUTCOME_DIR, "finngen_R10_manifest.tsv")
  if (!file.exists(.fgm))
    stop("config: finngen_R10_manifest.tsv is missing, so the five FinnGen sample sizes cannot be ",
         "checked against the release that produced them.")
  .m <- data.table::fread(.fgm, select = c("phenocode", "num_cases", "num_controls"))
  .want <- c(CAD = "I9_CHD", HF = "I9_HEARTFAIL", T2D = "T2D", NAFLD = "NAFLD", AF = "I9_AF")
  for (.k in names(.want)) {
    .r <- .m[phenocode == .want[[.k]]]
    if (nrow(.r) != 1L)
      stop(sprintf("config: the FinnGen R10 manifest has %d row(s) for %s (%s); its sample size cannot be established", nrow(.r), .k, .want[[.k]]))
    .tot <- as.integer(.r$num_cases[1]) + as.integer(.r$num_controls[1])
    if (!identical(OUTCOME_N_TOTAL[[.k]], .tot) || !identical(OUTCOME_N_CASE[[.k]], as.integer(.r$num_cases[1])))
      stop(sprintf(paste0("config: %s is declared N=%d ncase=%d but FinnGen R10's own manifest says ",
                          "N=%d ncase=%d (%s). Every Wald ratio and every case-proportion in the ",
                          "analysis is computed from these."),
                   .k, OUTCOME_N_TOTAL[[.k]], OUTCOME_N_CASE[[.k]], .tot, .r$num_cases[1], .want[[.k]]))
  }
})
stopifnot(
  # Wuttke 2019 EA and MEGASTROKE EUR: no release manifest exists, so these stay literal. Sources are
  # named so the next reader can check them against the papers rather than against another copy here.
  OUTCOME_N_TOTAL[["eGFR"]]  == 567460L,  # Wuttke 2019 EA eGFR
  OUTCOME_N_TOTAL[["CKD"]]   == 480698L,  # Wuttke 2019 EA CKD
  OUTCOME_N_CASE[["CKD"]]    == 41395L,   # Wuttke 2019 EA, NOT Pattaro 7480
  is.na(OUTCOME_N_CASE[["eGFR"]]),        # quantitative: no case count
  OUTCOME_N_TOTAL[["Stroke"]] == 446696L, # MEGASTROKE EUR
  OUTCOME_N_CASE[["Stroke"]]  == 40585L,
  !any(OUTCOME_UKB)                      # all 8 disease outcomes UKB-free (only UKB-PPP exposure + SBP-intermediate carry UKB)
)

# Backup outcome IDs — DEPRECATED for the corrected UKB-free cascade (2026-06-30). The old backups
# (ieu-a-7 / ieu-b-4873 / ebi-a-GCST006867 / finn-b-I9_AF) were UKB-overlapping or old-FinnGen and are a
# stale-fallback TRAP. Set equal to the UKB-free primaries so any fallback cannot reintroduce UKB overlap.
OUTCOMES_BACKUP <- OUTCOMES[intersect(c("CAD","HF","T2D","AF"), names(OUTCOMES))]

# ── IEU OpenGWAS JWT Token ──────────────────────────────────────────────────────
OPENGWAS_JWT <- Sys.getenv("OPENGWAS_JWT", "")   # NEW-2 fix: read from env, NOT hardcoded (literal scrubbed). Only needed for the one-time OpenGWAS download; the localized pipeline reads local files and needs no token. The previously-committed token is COMPROMISED — rotate it (self-expires ~2026-06-18).
Sys.setenv(OPENGWAS_JWT = OPENGWAS_JWT)   # exports whatever env already holds (blank if unset)

# ── Significance Thresholds for Phase 2 ────────────────────────────────────────
# Bonferroni denominator = (#unique cis-pQTL instrument genes) × (#CKM outcomes).
# DERIVED from the instrument set actually in use — NEVER hardcode the protein count.
# A stale literal here silently mis-draws every
# significance line. TO SWITCH GENERATIONS, SET CKM_PHASE1_SET (or pre-set PHASE1_SET) before
# sourcing -- this line used to say "override INSTRUMENTS_FILE", which the assertion twenty lines
# below rejects: overriding it to anything other than PHASE1_INSTRUMENTS aborts, so following the
# instruction left the run with no threshold at all. INSTRUMENTS_FILE follows the selection, it
# does not drive it. Scripts that load a proteome-wide MR table should prefer
#   bonferroni_alpha(data.table::uniqueN(<that_table>$gene))
# so the threshold always tracks the data being plotted.
N_CKM_OUTCOMES   <- length(OUTCOMES)                                          # 8 primary CKM outcomes (derived)
bonferroni_alpha <- function(n_genes, n_outcomes = N_CKM_OUTCOMES) {
  # VALIDATED, because this is exported and the failure is silent: bonferroni_alpha(0) returned Inf
  # (every comparison passes), NA returned NA (none pass, and `pval < NA` is FALSE everywhere), and a
  # vector returned a vector that then recycled against p-values row by row.
  for (.a in list(list(n_genes, "n_genes"), list(n_outcomes, "n_outcomes")))
    if (!is.numeric(.a[[1]]) || length(.a[[1]]) != 1L || !is.finite(.a[[1]]) || .a[[1]] < 1 ||
        .a[[1]] != as.integer(.a[[1]]))
      stop(sprintf("bonferroni_alpha: %s must be a single whole number >= 1, got %s", .a[[2]],
                   paste(utils::capture.output(utils::str(.a[[1]])), collapse = " ")))
  0.05 / (n_genes * n_outcomes)
}
# The instrument set the THRESHOLD is computed from, which must be the one the ANALYSIS reads.
if (!exists("INSTRUMENTS_FILE")) INSTRUMENTS_FILE <- PHASE1_INSTRUMENTS
# ONE INSTRUMENT SET PER RUN, asserted. PHASE1_INSTRUMENTS, PHASE1_DIR and INSTRUMENTS_FILE are each
# retained by their own exists() guard, so a caller that pre-sets one of them -- or a second source()
# in the same session after changing one -- left them naming different generations with nothing
# noticing. Measured: selecting phase1_hg38 while pre-setting INSTRUMENTS_FILE to the previous RDS
# put the MR screen on 1,948 genes while BONFERRONI_ALPHA was computed from 1,682, a threshold 15.8%
# looser than the analysis it was applied to. Pre-setting PHASE1_DIR instead split the MR screen from
# the outcome extractor. The three describe one thing and must agree about it.
if (!identical(normalizePath(INSTRUMENTS_FILE, mustWork = FALSE),
               normalizePath(PHASE1_INSTRUMENTS, mustWork = FALSE)))
  stop(sprintf(paste0("config: INSTRUMENTS_FILE (%s) and PHASE1_INSTRUMENTS (%s) name different ",
                      "instrument sets. The Bonferroni denominator comes from the first and the ",
                      "analysis reads the second, so they would describe different studies."),
               INSTRUMENTS_FILE, PHASE1_INSTRUMENTS))
if (!identical(normalizePath(PHASE1_DIR, mustWork = FALSE),
               normalizePath(dirname(PHASE1_INSTRUMENTS), mustWork = FALSE)))
  stop(sprintf(paste0("config: PHASE1_DIR (%s) is not the directory of PHASE1_INSTRUMENTS (%s). ",
                      "cascade/01 reads the second and localize/02 and /03 read the first, so the ",
                      "MR screen and the outcome extraction would run on different generations."),
               PHASE1_DIR, dirname(PHASE1_INSTRUMENTS)))
# N_PROTEINS_CANONICAL: the canonical instrument-GENE count. SINGLE SOURCE for (a) the "expected
# count" integrity assertions (01_mr / 02_promote) and (b) the Bonferroni denominator below.
#
# COUNTED FROM THE INSTRUMENT TABLE IN USE, never written here. It was the literal 1682L, and the
# rebuild made that stale -- but replacing it with 1948L would only move the staleness to the next
# rebuild. A count that tracks the data has no business being typed by hand: this is the same rule
# the comment twenty lines above already states ("NEVER hardcode the protein count").
#
# It also removes a real integrity hazard. The literal made the guards point the wrong way -- a
# correctly rebuilt 1,948-gene table ABORTED at 01_mr's assertion while a stale table that happened
# to hold 1,682 passed silently -- and it put a human at the exact moment the multiple-testing
# denominator is chosen. Choosing among 1,635 / 1,682 / 1,948 with the results already visible is
# indistinguishable from threshold shopping to a reader, however honest the intent. Derived, there
# is nothing to choose: the denominator is however many proteins were put to test, which is what the
# manuscript's own rule says it is.
#
# Fail-loud rather than NA: a missing denominator makes `pval < NA` FALSE everywhere and flags
# nothing at all, which is why a literal was used here in the first place.
# NA and "" are not proteins. unique() keeps both and length() would count them into the
# denominator, and a blank would also enter PROTEIN_UNIVERSE and pass its membership check.
# ONE READ, and the hash is of the bytes that were parsed. Hashing the path in a second filesystem
# read leaves a window in which the file can change, so the threshold could come from one file and
# the recorded identity from another while every path check still passed.
# THE BYTES ARE PINNED FOR THE SESSION, not just the set NAME. _phase1_set.R refuses a changed
# CKM_PHASE1_SET on a second source(), but an in-place rebuild of the SAME name replaces the file
# under that name: the second source then re-reads, re-certifies against the rebuilt table's new
# (valid) certificate, and silently moves N_PROTEINS_CANONICAL, BONFERRONI_ALPHA and
# PROTEIN_UNIVERSE mid-session, while everything already computed still describes the first table.
.PREV_INSTR_MD5  <- if (exists("INSTRUMENTS_MD5")) INSTRUMENTS_MD5 else NULL
.INSTR_BYTES     <- tryCatch(readBin(INSTRUMENTS_FILE, "raw", file.info(INSTRUMENTS_FILE)$size),
                             error = function(e) NULL)
# COMPARED BEFORE IT IS PUBLISHED. Assigning INSTRUMENTS_MD5 first and checking afterwards left
# the refusal half-applied: the session kept the first table's thresholds, accessor and universe but
# already carried the SECOND table's md5, so anything running after a caught error analysed
# generation A while stamping it B -- and a further source() then saw the new hash as "previous",
# agreed with itself, and moved the thresholds.
.NEW_INSTR_MD5   <- if (is.null(.INSTR_BYTES)) {
                      NA_character_
                    } else {
                      digest::digest(.INSTR_BYTES, algo = "md5", serialize = FALSE)
                    }
if (!is.null(.PREV_INSTR_MD5) && !is.na(.PREV_INSTR_MD5) &&
    !identical(.PREV_INSTR_MD5, .NEW_INSTR_MD5))
  stop(sprintf(paste0("config: %s changed while this session was running (was %s, now %s). Every ",
                      "threshold, gene universe and provenance stamp computed before this point ",
                      "describes the first table; re-sourcing would replace them and leave the two ",
                      "halves of the session describing different instrument sets. Start a fresh ",
                      "R process."),
               INSTRUMENTS_FILE, substr(.PREV_INSTR_MD5, 1, 12), substr(.NEW_INSTR_MD5, 1, 12)))
INSTRUMENTS_MD5  <- .NEW_INSTR_MD5
rm(.PREV_INSTR_MD5, .NEW_INSTR_MD5)
# "unknown" lets the header say how the bytes are compressed. It was hardcoded "gzip", which is
# saveRDS()'s default but not its only setting: an otherwise identical table written with
# compress="xz"/"bzip2"/FALSE is read fine by readRDS() and failed here, taking every threshold with
# it. What is being read is an RDS, not specifically a gzip stream.
.INSTR_DT        <- tryCatch(unserialize(memDecompress(.INSTR_BYTES, "unknown")),
                             error = function(e) NULL)
# REFUSED, NOT DROPPED. Silently discarding NA/blank genes kept those instruments in the table the
# accessor hands out while removing their protein from the denominator: the row is analysed and the
# multiple-testing grid does not know it exists. Both present tables have none, so this changes
# nothing today and cannot be reached quietly tomorrow.
.INSTR_GENES     <- tryCatch({ .g <- as.character(.INSTR_DT[["gene"]])
                               # PADDING IS A DIFFERENT IDENTIFIER, and it was only used to
                               # detect emptiness: "GENE" and " GENE " were counted as two
                               # proteins in the denominator while 15/16 trim exposure names and
                               # collapse them back into one. Refused rather than trimmed here --
                               # silently rewriting an identifier makes the table and the count
                               # disagree about what the table says.
                               .bad <- is.na(.g) | !nzchar(trimws(.g)) | .g != trimws(.g)
                               if (any(.bad))
                                 stop(sprintf("%d instrument row(s) have an unusable gene identifier (blank, NA, or whitespace-padded)",
                                              sum(.bad)))
                               unique(.g) }, error = function(e) NULL)
if (is.null(.INSTR_GENES) || !length(.INSTR_GENES))
  stop(sprintf(paste0("config: cannot count the multiple-testing denominator -- INSTRUMENTS_FILE ",
                      "(%s) is unreadable or has no `gene` column. Every significance threshold ",
                      "rests on it, and a missing denominator would silently flag nothing."),
               INSTRUMENTS_FILE))
# -- THE WHOLE TABLE ITS BUILDER WROTE, not merely a readable one ------------------------------
# The count below is whatever this file happens to hold, and nothing asked whether that is the
# table 09 produced. A TRUNCATION is the dangerous shape, because it fails in the loosening
# direction and looks entirely normal: a one-row RDS yields N_PROTEINS_CANONICAL = 1 and a
# Bonferroni threshold of 0.00625 -- 1,900x looser than the real one -- and every guard that
# compares the analysis against "the canonical count" then agrees with it, because they all read
# this same number. A same-size substitution is likewise invisible to a count.
# 09 writes a single-row certificate beside the RDS naming the md5 and row count of the bytes it
# installed, atomically, after the table itself. Comparing the bytes we PARSED against it settles
# identity and completeness in one check -- and it is the same certificate idiom this file already
# applies to the coordinate table and 09 applies to the LD panel, not a new mechanism.
# The previous generation predates the certificate, so its absence is not automatically an error --
# it is the same deliberate act, under the same flag, that selecting a previous generation already
# requires.
# file.exists() FIRST: fread() does not fail on a missing path, it treats the string as literal
# input and returns an empty data.table with a warning, so "no certificate" arrived here looking
# like "a certificate with the wrong number of rows" -- the previous generation, which predates the
# certificate, then aborted with a message about a mismatch that did not exist and could not be
# reached by the flag that is supposed to allow it.
.IMANI_F <- paste0(INSTRUMENTS_FILE, ".manifest.tsv")
.IMANI   <- if (file.exists(.IMANI_F))
              tryCatch(data.table::fread(.IMANI_F), error = function(e) NULL) else NULL
if (is.null(.IMANI)) {
  # TWO DIFFERENT REFUSALS, said separately. Combining them told a caller who HAD set the flag to
  # set the flag.
  if (!ckm_env_yes("CKM_ALLOW_STALE_PHASE1"))
    stop(sprintf(paste0("config: %s has no .manifest.tsv certificate, so there is nothing to say ",
                        "whether it is the whole table its builder wrote. A truncated table would ",
                        "set the multiple-testing denominator to its own row count and loosen every ",
                        "threshold silently. Rebuild it with 09_phase1_harmonise_hg38.R, or set ",
                        "CKM_ALLOW_STALE_PHASE1=1 to use a recorded previous generation on ",
                        "purpose."), INSTRUMENTS_FILE))
  # ...and the flag authorises a recorded previous generation, not any uncertified table that
  # happens to be readable -- the same restriction the selector applies.
  if (!ckm_is_known_previous(INSTRUMENTS_FILE))
    stop(sprintf(paste0("config: CKM_ALLOW_STALE_PHASE1 is set, but %s has no certificate AND is ",
                        "not one of the recorded previous generations (%s). Nothing states what ",
                        "this table is, so no threshold derived from it could be described in the ",
                        "manuscript. Add its md5 to CKM_PREVIOUS_GENERATION_MD5 if it really is ",
                        "one."),
                 INSTRUMENTS_FILE, paste(names(CKM_PREVIOUS_GENERATION_MD5), collapse = ", ")))
  cat("** config:", basename(INSTRUMENTS_FILE), "has no certificate; proceeding because",
      "CKM_ALLOW_STALE_PHASE1 is set. The denominator is UNVERIFIED. **\n")
} else if (nrow(.IMANI) == 1L &&
           (!all(c("source", "source_md5") %in% names(.IMANI)) ||
            !file.exists(file.path(dirname(INSTRUMENTS_FILE), .IMANI$source[1])) ||
            !identical(tolower(as.character(.IMANI$source_md5[1])),
                       tolower(unname(tools::md5sum(file.path(dirname(INSTRUMENTS_FILE),
                                                              .IMANI$source[1]))))))) {
  # ...AND IT MUST DESCRIBE THE EXTRACTION THAT IS ON DISK NOW. The certificate authenticates the
  # RDS against itself, which a NEW extraction cannot disturb: a harmonised table is legitimately a
  # strict subset of its extraction (1,964 of 2,053 rows here), so when a later extraction lands
  # carrying the same rows plus new ones, the old RDS still passes the selector's subset-and-beta
  # test and its own certificate, while its certificate names a source CSV that no longer exists in
  # that form. The result is a rebuild that looks complete and an analysis still reading the
  # generation before it.
  # ...and the three ways this can fail are three different sentences, because "the certificate
  # cannot authenticate the extraction" and "the extraction changed" call for different actions.
  # The check used to be conditional on the source columns AND the source file existing, which made
  # it optional exactly when it mattered: a certificate with no source fields, or an RDS whose
  # extraction had been deleted, fell through to the rows+md5 branch and was accepted as canonical.
  .src_f <- if (all(c("source", "source_md5") %in% names(.IMANI)))
              file.path(dirname(INSTRUMENTS_FILE), .IMANI$source[1]) else NA_character_
  stop(if (is.na(.src_f))
         sprintf(paste0("config: the certificate beside %s has no source/source_md5 column, so it ",
                        "cannot say which extraction this table was harmonised from. Re-run ",
                        "09_phase1_harmonise_hg38.R, which writes both."),
                 basename(INSTRUMENTS_FILE))
       else if (!file.exists(.src_f))
         sprintf(paste0("config: %s names extraction %s, which is not on disk. The instrument set ",
                        "cannot be checked against the extraction it claims to summarise, and ",
                        "localize/02 and /03 read that extraction's directory. Restore it or ",
                        "re-run the Phase-1 extraction and 09_phase1_harmonise_hg38.R."),
                 basename(INSTRUMENTS_FILE), .IMANI$source[1])
       else
         sprintf(paste0("config: %s was harmonised from %s with md5 %s, but that file now hashes ",
                        "to %s. The extraction was rebuilt and this table was not: re-run ",
                        "09_phase1_harmonise_hg38.R so the instruments describe the extraction on ",
                        "disk."),
                 basename(INSTRUMENTS_FILE), .IMANI$source[1],
                 substr(as.character(.IMANI$source_md5[1]), 1, 12),
                 substr(unname(tools::md5sum(.src_f)), 1, 12)))
} else if (nrow(.IMANI) != 1L || !all(c("rows", "md5") %in% names(.IMANI)) ||
           !identical(tolower(as.character(.IMANI$md5[1])), tolower(INSTRUMENTS_MD5)) ||
           !identical(as.integer(.IMANI$rows[1]), nrow(.INSTR_DT))) {
  stop(sprintf(paste0("config: %s does not match its certificate (certificate: %s rows / md5 %s; ",
                      "file: %d rows / md5 %s). The multiple-testing denominator is counted from ",
                      "this table, so a partial or substituted one would move every significance ",
                      "threshold without changing anything a reader can see."),
               INSTRUMENTS_FILE,
               if ("rows" %in% names(.IMANI)) paste(.IMANI$rows[1]) else "?",
               if ("md5" %in% names(.IMANI)) substr(as.character(.IMANI$md5[1]), 1, 12) else "?",
               nrow(.INSTR_DT), substr(INSTRUMENTS_MD5, 1, 12)))
# ...AND THE BUILD IT DECLARES, plus the code that made it. Everything above authenticates the BYTES
# of the selected table; none of it looks at which genome build those bytes are in, which is the one
# property this whole generation is being rebuilt over. 03b_gene_coords_hg38.R does enforce these,
# but only when it runs -- the MR screen and every other consumer read the instrument table long
# before that. Four string comparisons and two md5 of small script files; no extra table is read.
} else if (!all(c("gene_coord_build", "chrpos_build", "builder_md5", "extractor_code_md5") %in% names(.IMANI)) ||
           !identical(as.character(.IMANI$gene_coord_build[1]), "GRCh38") ||
           !identical(as.character(.IMANI$chrpos_build[1]), "GRCh38") ||
           !identical(tolower(as.character(.IMANI$builder_md5[1])),
                      tolower(unname(tools::md5sum(file.path(BASE_DIR, "scripts/_shared/localize/09_phase1_harmonise_hg38.R"))))) ||
           !identical(tolower(as.character(.IMANI$extractor_code_md5[1])),
                      tolower(unname(tools::md5sum(file.path(BASE_DIR, "scripts/_shared/localize/08_phase1_extract_hg38.py")))))) {
  stop(sprintf(paste0("config: %s declares gene_coord_build=%s / chrpos_build=%s and names producers ",
                      "%s / %s, which are absent, not GRCh38, or not the 09/08 on disk. The table's ",
                      "bytes authenticate; the BUILD they are in does not, and applying one build's ",
                      "boundaries to another build's positions is the defect this generation exists ",
                      "to remove."),
               INSTRUMENTS_FILE,
               if ("gene_coord_build" %in% names(.IMANI)) as.character(.IMANI$gene_coord_build[1]) else "<absent>",
               if ("chrpos_build" %in% names(.IMANI)) as.character(.IMANI$chrpos_build[1]) else "<absent>",
               if ("builder_md5" %in% names(.IMANI)) substr(as.character(.IMANI$builder_md5[1]), 1, 12) else "<absent>",
               if ("extractor_code_md5" %in% names(.IMANI)) substr(as.character(.IMANI$extractor_code_md5[1]), 1, 12) else "<absent>"))
}
N_PROTEINS_CANONICAL <- length(.INSTR_GENES)
# THE CERTIFIED TABLE ITSELF, for consumers that would otherwise re-open the path. Only the
# PATHNAME was exported, so a script that sources config and then calls readRDS(PHASE1_INSTRUMENTS)
# gets whatever the path holds at that later moment, while its Bonferroni threshold and gene
# universe still come from the bytes config captured. The certificate above authenticates those
# bytes; it is not consulted by the later read, and it cannot be: a rebuild that finishes AFTER
# config was sourced installs a valid table with a valid new certificate, and the running process
# then analyses generation B against generation A's threshold with every check still passing. The
# measured cost of that split is 3.208419e-06 vs 3.715815e-06, 15.8%, with three stored MR results
# in between. copy() because this is a data.table and consumers mutate with `:=`; a shared object
# would let one script's edits reach the next reader in the same process.
# ...and the snapshot itself is HELD IN THE FUNCTION'S OWN ENVIRONMENT, not left in .GlobalEnv.
# data.table mutates by reference: `data.table::set(.INSTR_DT, ...)` changed what every later
# caller received while INSTRUMENTS_MD5, N_PROTEINS_CANONICAL and BONFERRONI_ALPHA went on
# describing the original bytes -- the certificate proves what was READ, and a mutable global
# undoes that between the read and the use.
read_phase1_instruments <- local({
  .snapshot <- .INSTR_DT
  function() data.table::copy(.snapshot)
})
rm(.INSTR_DT)
# DERIVED from whatever INSTRUMENTS_FILE is in effect -- informational + a consistency tripwire ONLY.
# ONE read, TWO derivations: the gene COUNT (the tripwire below) and the gene UNIVERSE.
# PROTEIN_UNIVERSE exists because primary_fdr_set() must prove its BH numerator lives on the grid its
# denominator describes, and a COUNT cannot do that: a stale/non-CKM gene that replaces a missing
# canonical one keeps uniqueN(gene) under the cap and passes silently (objective review r6, 2026-07-16).
# NULL when INSTRUMENTS_FILE is unreadable -- primary_fdr_set() STOPS on that rather than degrading to
# the count check it was just shown to be insufficient.
#
# WHAT THIS IS NOT (r8/r10, 2026-07-16): it is the universe INSTRUMENTS_FILE held when config was
# SOURCED (read once, stored -- re-pointing INSTRUMENTS_FILE afterwards changes nothing; override it
# before sourcing, or pass protein_universe= explicitly). It is NOT a pinned canonical identity.
# A SAME-SIZE substitution -- a different 1,682-gene set -- passes every guard: right count, right size,
# and membership checked against the substituted set itself.
# An earlier version of this note argued the rest of the chain catches it, because mr_all_results.rds is
# built FROM these instruments so a one-sided swap would leave mr's genes outside the new universe.
# THAT WAS WRONG (r10): only 1,504 of the 1,682 genes have a primary MR row -- mr_1682_no_result_genes.csv
# lists the other 178 -- so a swap touching only unanalysed genes passes membership AND size while
# changing what the denominator denotes.
# So state the guarantee exactly: these guards prove the BH NUMERATOR's genes are in the universe and
# that the universe is the right SIZE. They do NOT prove the universe's IDENTITY -- and nothing else
# does either. An earlier version of this note said identity was "pinned by phase3/33's md5 of
# ukb_cis_instruments.rds"; that is not pinning (r11). 33 RECORDS the md5 of whatever file is present
# when it runs; it never compares against an expected canonical value. A substitution made before the
# re-run is simply blessed with a new md5 and recorded as if it were always so.
# So: canonical universe identity is UNPINNED, everywhere, today. Whether to pin it (a locked sorted gene
# list or an expected-md5 manifest, at the cost of breaking every legitimate rebuild until it is updated
# by hand) is a PI call tracked in REBUILD_LEDGER §2.1c(c).
# The SAME read as the count above -- one readRDS, not two. The universe and the count cannot
# disagree by construction, so the old drift warning ("INSTRUMENTS_FILE has N genes but canonical is
# M") now compares a number with itself and is gone; what it used to detect was the literal going
# stale, which is the thing that no longer exists.
PROTEIN_UNIVERSE <- .INSTR_GENES
N_PROTEINS       <- length(.INSTR_GENES)
# THE BYTES THIS THRESHOLD WAS COMPUTED FROM. Config stored a gene vector and nothing else, so a
# same-size replacement of the RDS was undetectable and a different-size one left the earlier
# threshold sitting in memory while later code re-read the file. The two real tables differ
# materially -- 1,948 vs 1,682 genes, only 1,589 shared -- so "same size" is not a remote scenario
# once a third generation exists. Consumers that re-read PHASE1_INSTRUMENTS can compare against this.
if (is.na(INSTRUMENTS_MD5))
  stop("config: cannot hash ", INSTRUMENTS_FILE, " -- the instrument set the significance ",
       "thresholds rest on cannot be identified.")
# 2026-07-01 (Codex Bonferroni #1): the Bonferroni denominator uses the CANONICAL count, NOT the derived
# N_PROTEINS. A derived NA (missing/unreadable INSTRUMENTS_FILE) must NOT silently blank the significance
# threshold (NA denom -> `pval < NA` -> nothing flagged). Canonical is fixed + 01_mr-guarded -> never NA, never 1635.
BONFERRONI_ALPHA <- bonferroni_alpha(N_PROTEINS_CANONICAL)
FDR_ALPHA        <- 0.05                  # Benjamini-Hochberg FDR
SUGGESTIVE_P     <- 0.05                  # Enter Phase 3 threshold (IVW p < 0.05)

# ── Local 1KG EUR LD Reference (Phase 3-C SuSiE) ──────────────────────────────
PLINK_BIN  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/ld_reference/plink.exe")
LD_EUR_DIR <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/ld_reference/1kg_eur")

# ── Artifact paths that FOLLOW THE SELECTED GENERATION ────────────────────────────────────────
# These were spelled "..._1682" wherever they appeared, in the producer (cascade/01) and in every
# consumer. That was true of one generation and became a lie the moment Phase 1 was rebuilt: the
# GRCh38 rebuild writes 1,948 genes, and a reader who opens mr_all_results_1682.rds after the
# rebuild has no way to know which set is inside. Worse, a producer and a consumer that disagree
# about the suffix silently pair a new MR table with an old proxy cache.
# DERIVED FROM PHASE1_SET, so the name cannot outlive the thing it names, and so the previous
# generation's artifacts keep their own names instead of being overwritten (PI decision,
# 2026-07-31). The extracts and the proxy cache are shared by cascade/01, cascade/09_steiger,
# localize/41b, 53 and 54, which is why they are defined once here rather than in each script.
PHASE1_TAG        <- PHASE1_SET
LOCAL_EXTRACT_DIR <- file.path(OUTCOME_DIR, paste0("local_extracts_", PHASE1_TAG))
LD_PROXY_PREFIX   <- file.path(LD_EUR_DIR,  paste0("instr_ld_proxies_", PHASE1_TAG))
INSTR_SNPLIST     <- file.path(LD_EUR_DIR,  paste0("instr_", PHASE1_TAG, ".snplist"))
# A helper rather than a list of constants, because the MR outputs live under RESULTS_DIR, which is
# itself overridable (the corrected cascade repoints it at a shadow tree) and so cannot be baked in
# here. `mr_tag("all_results")` -> "mr_all_results_phase1_hg38.rds".
mr_artifact <- function(stem, ext = "rds")
  sprintf("mr_%s_%s.%s", stem, PHASE1_TAG, ext)

# ── Local liftover + gene-coord references (replace coloc's Ensembl REST / web pyliftover) ──
LIFTOVER_CHAIN   <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/ld_reference/hg19ToHg38.over.chain.gz")  # vendored UCSC chain; use single-arg LiftOver() = offline
# The other direction, for the two GRCh37 outcomes (CKDGen CKD / eGFR). The pQTL side is GRCh38, so
# a position comparison against their `Pos_b37` needs the pQTL positions projected DOWN, not the
# outcome projected up -- projecting the outcome would move millions of rows instead of thousands.
LIFTOVER_CHAIN_38_TO_37 <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/ld_reference/hg38ToHg19.over.chain.gz")
GENE_HG38_COORDS <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/ld_reference/gene_hg38_coords.tsv")      # local gene->GRCh38 coords; replaces get_gene_hg38_coords() Ensembl REST

# ── MHC / HLA region: PRE-SPECIFIED fine-mapping exclusion ────────────────────────
# 2026-07-02 (editor-endorsed): SuSiE fine-mapping is unreliable under the MHC's extreme long-range LD, and
# coloc.abf's single-causal-variant assumption is ALSO violated there -> a pre-specified exclusion. A cis pair
# whose gene lies in the extended MHC is reported SuSiE "not_evaluable (MHC_excluded)", and its coloc.abf is
# flagged for cautious interpretation (symmetric: we do not trust abf in the MHC while excluding SuSiE).
# Extended MHC (GRCh38): chr6:25-34 Mb. is_mhc_gene() returns NA if coords are unavailable so a caller must
# fail loud rather than silently treat an undetermined pair as non-MHC.
MHC_CHR <- "6"; MHC_START <- 25000000L; MHC_END <- 34000000L
# This was a SECOND, uncertified reader of the coordinate table, with a cache that was never
# invalidated. A previous revision "fixed" it by adding a comment saying it defers to the certified
# reader. It deferred to nothing, and a comment asserting a behaviour the code does not have is worse
# than the untouched defect, because it stops anyone from looking again. It now actually defers.
#
# Why this function in particular: MHC membership decides whether a pair is reported as
# not_evaluable for colocalisation. Reading a retired table would let this function and the coloc
# layer disagree about which genes are in the MHC, and 99_completeness_audit.R cross-checks one
# against the other -- so the disagreement would surface as an unexplained abort rather than as the
# stale-table problem it actually is.
is_mhc_gene <- function(g) {
  if (exists(".assert_coords_certified", envir = .GlobalEnv)) {
    cur <- get(".assert_coords_certified", envir = .GlobalEnv)()
    # The RUN PIN too, not only the certificate. get_gene_hg38_coords() refuses a table that changed
    # mid-run; this function validated the new table and simply adopted it, so a rebuild during an
    # MHC loop produced decisions spanning two versions -- and 99_completeness_audit.R compares this
    # function against the coloc layer, so the disagreement would surface as an unexplained abort.
    # ESTABLISHES the pin as well as enforcing it. Both SuSiE entry points call this function before
    # their first coordinate lookup, so leaving the pin to get_gene_hg38_coords() meant a rebuild
    # during that opening loop produced MHC decisions from two table versions with nothing set to
    # notice.
    .pin <- get0(".COLOC_RUN_COORD_MD5", envir = .GlobalEnv)
    if (!is.null(cur) && is.null(.pin)) {
      assign(".COLOC_RUN_COORD_MD5", cur$md5, envir = .GlobalEnv)
      .pin <- cur$md5
    }
    if (!is.null(cur) && !is.null(.pin) && !identical(.pin, cur$md5))
      stop("is_mhc_gene(): the coordinate table changed while this run was in progress (pinned ",
           substr(.pin, 1, 12), ", now ", substr(cur$md5, 1, 12), "). MHC decisions made before and ",
           "after would come from different tables.")
    if (!is.null(cur) &&
        !identical(get0(".GENE_COORDS_CACHE_MD5", envir = .GlobalEnv), cur$md5)) {
      assign(".GENE_COORDS_CACHE",
             data.table::fread(text = rawToChar(cur$bytes)), envir = .GlobalEnv)
      assign(".GENE_COORDS_CACHE_MD5", cur$md5, envir = .GlobalEnv)
    }
  } else if (!exists(".GENE_COORDS_CACHE_MD5", envir = .GlobalEnv) ||
             !identical(get0(".GENE_COORDS_CACHE_MD5", envir = .GlobalEnv),
                        paste0(tolower(unname(tools::md5sum(GENE_HG38_COORDS))), "|",
                               tolower(unname(tools::md5sum(paste0(GENE_HG38_COORDS, ".manifest.tsv"))))))) {
    # KEYED ON THE CERTIFICATION STAMP, not on the cache's existence. Only the two paths that
    # validate a table set .GENE_COORDS_CACHE_MD5, whereas .GENE_COORDS_CACHE is an ordinary global
    # anything can create: a pre-existing one skipped the certificate entirely and was then trusted
    # for every MHC decision in the session (a hand-made row claiming chr6:26,000,000 was accepted
    # without a single file being read).
    # Reached only by tools that never load the coloc engine. Uncertified -- and it says so out loud
    # rather than silently, because this is the path that can disagree with the coloc layer.
    # The SAME contract the coloc engine applies, not a weaker one. Checking only that the sidecar's
    # md5 matches the table accepts a self-consistent but stale 03b output -- one that
    # get_gene_hg38_coords() would refuse -- and then decides MHC membership from it, which is what
    # 99_completeness_audit.R cross-checks against the coloc layer.
    .cert <- paste0(GENE_HG38_COORDS, ".manifest.tsv")
    .m <- if (file.exists(.cert)) tryCatch(data.table::fread(.cert), error = function(e) NULL)
    .b3  <- file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/03b_gene_coords_hg38.R")
    .b93 <- file.path(BASE_DIR, "scripts/_shared/localize/93_build_gene_coords_hgnc.R")
    .gtf <- file.path(BASE_DIR, "data/annotation/gencode.v50.basic.annotation.gtf.gz")
    # Read ONCE into memory and parsed from those bytes, exactly as the coloc reader does. Hashing
    # the path and then calling fread() on it again is two reads: 03b's rename landing between them
    # validates the OLD table and caches the NEW one, which is never revalidated -- and that is
    # precisely the installation intermediate state both readers are supposed to refuse.
    .bytes <- if (file.exists(GENE_HG38_COORDS))
                readBin(GENE_HG38_COORDS, "raw", file.info(GENE_HG38_COORDS)$size) else NULL
    .tblmd5 <- if (!is.null(.bytes)) digest::digest(.bytes, algo = "md5", serialize = FALSE) else NA
    .h <- function(f) if (file.exists(f)) tolower(unname(tools::md5sum(f))) else NA_character_
    .same <- function(a, b) !is.na(b) && identical(tolower(a), b)
    # FATAL, matching the coloc engine. A warning here is not a smaller version of the same
    # protection: this cache is never revalidated, so one warning at the start of a session is all
    # anyone gets while every subsequent MHC decision comes from whatever was loaded then. MHC
    # membership decides whether a pair is reported as not_evaluable, and 99_completeness_audit.R
    # cross-checks this function against the coloc layer -- so a disagreement surfaces as an
    # unexplained abort rather than as the stale table it actually is.
    # The SAME contract the coloc reader applies, including source/source_md5 and a single-row
    # certificate. Anything weaker means a table coloc refuses can still decide MHC membership --
    # and 99_completeness_audit.R cross-checks the two against each other.
    if (is.null(.m) || nrow(.m) != 1L ||
        # REQUIRED HERE TOO, matching R/coloc.R's .assert_coords_certified(). The four build/window/
        # manifest fields were checked only `if present` below, so a certificate that simply omitted
        # them passed this branch while the engine refused the same table -- two guards over one file
        # disagreeing about what it must declare, with the weaker one deciding MHC membership.
        !all(c("md5", "builder", "builder_md5", "upstream_builder_md5", "gtf_md5",
               "source", "source_md5", "phase1_gene_coord_build", "phase1_chrpos_build",
               "phase1_manifest_md5", "cis_window_kb") %in% names(.m)) ||
        is.na(.tblmd5) || !identical(tolower(.tblmd5), tolower(.m$md5[1])) ||
        !.same(.m$source_md5[1],
               .h(file.path(BASE_DIR, "data/annotation", basename(.m$source[1])))) ||
        !identical(.m$builder[1], "03b_gene_coords_hg38.R") ||
        !.same(.m$builder_md5[1], .h(.b3)) ||
        !.same(.m$upstream_builder_md5[1], .h(.b93)) ||
        !.same(.m$gtf_md5[1], .h(.gtf)) ||
        !all(c("idmapping_md5", "release_listing_md5", "instruments_md5", "file", "rows")
             %in% names(.m)) ||
        # The three the coloc reader compares and this one did not: without them a certificate that
        # coloc refuses can still decide MHC membership, and 99_completeness_audit.R then reports
        # the disagreement as an unexplained abort.
        # The generation this session pinned, matching the coloc reader. Hashing the path afresh
        # accepts a table certified against a generation the session is not running.
        !.same(.m$instruments_md5[1],
               if (exists("INSTRUMENTS_MD5")) tolower(INSTRUMENTS_MD5) else .h(PHASE1_INSTRUMENTS)) ||
        # The same contract the coloc reader applies to the recorded Phase-1 provenance: a build
        # relabelled after installation, a changed manifest, or a different cis window must not leave
        # MHC membership being decided from a table certified against something else.
        !all(c(.m$phase1_gene_coord_build[1], .m$phase1_chrpos_build[1]) == "GRCh38") ||
        !.same(.m$phase1_manifest_md5[1], .h(paste0(PHASE1_INSTRUMENTS, ".manifest.tsv"))) ||
        # Digits first, like `rows` below and like the engine: as.integer("500.9") is 500, so a
        # certificate stating a window this run does not use was accepted as though it matched.
        (exists("CIS_WINDOW_KB") &&
         (!grepl("^[[:space:]]*[0-9]+[[:space:]]*$", as.character(.m$cis_window_kb[1])) ||
          !identical(as.integer(.m$cis_window_kb[1]), as.integer(CIS_WINDOW_KB)))) ||
        !identical(.m$file[1], basename(GENE_HG38_COORDS)) ||
        !grepl("^[[:space:]]*[0-9]+[[:space:]]*$", as.character(.m$rows[1])) ||
        !identical(as.integer(.m$rows[1]),
                   sum(.bytes == as.raw(10L)) - 1L) ||
        !.same(.m$idmapping_md5[1],
               .h(file.path(BASE_DIR, "data/annotation/HUMAN_9606_idmapping.dat.gz"))) ||
        !.same(.m$release_listing_md5[1],
               .h(file.path(BASE_DIR, "data/ukb_ppp/eur_folder_contents.json"))))
      {
        .msg <- paste0("is_mhc_gene(): ", basename(GENE_HG38_COORDS), " does not match its ",
                       "certificate, so MHC decisions would come from an unverified table. Rebuild ",
                       "it with cascade/03b_gene_coords_hg38.R, or source R/coloc.R first so the ",
                       "certified reader is used.")
        # The override is HONOURED, not merely advertised. The previous message named
        # CKM_ALLOW_UNCERTIFIED_COORDS=1 and then stopped unconditionally, so a historical rerun the
        # coloc engine allows was refused here with instructions that did not work.
        # ...AND HONOURED ON THE SAME TERMS AS THE ENGINE. This used ckm_env_yes(), which accepts
        # true/t/yes, while R/coloc.R:599 requires the literal "1" -- so `=yes` let the MHC reader
        # through while colocalisation still refused the same table, and the two guards over one file
        # disagreed about what switches them off. The stricter of the two wins: an override that
        # silences the certificate must be unambiguous, not spelled four ways.
        if (!identical(Sys.getenv("CKM_ALLOW_UNCERTIFIED_COORDS"), "1"))
          stop(.msg, " Set CKM_ALLOW_UNCERTIFIED_COORDS=1 to proceed deliberately.")
        warning(.msg, " Proceeding because CKM_ALLOW_UNCERTIFIED_COORDS=1.", call. = FALSE,
                immediate. = TRUE)
      }
    # Parsed from the bytes that were validated, not from a fresh read of the path.
    assign(".GENE_COORDS_CACHE",
           if (!is.null(.bytes)) data.table::fread(text = rawToChar(.bytes)) else NULL,
           envir = .GlobalEnv)
    # ...and stamped, so this branch is not re-run on every lookup -- and so a later call cannot
    # mistake an uninstalled cache for a certified one. Not stamped when there was nothing to
    # validate: the check must run again rather than bless a NULL.
    # THE PAIR, not the table alone: everything validated above comes from the certificate, so a
    # cache keyed on the table's bytes goes stale the moment the certificate does.
    if (!is.null(.bytes))
      assign(".GENE_COORDS_CACHE_MD5",
             paste0(tolower(.tblmd5), "|",
                    tolower(unname(tools::md5sum(paste0(GENE_HG38_COORDS, ".manifest.tsv"))))),
             envir = .GlobalEnv)
  }
  gc <- get(".GENE_COORDS_CACHE", envir = .GlobalEnv)
  if (is.null(gc) || !all(c("gene","chrom","start","end") %in% names(gc))) return(NA)
  r <- gc[gc$gene == g]
  if (nrow(r) == 0) return(NA)   # 2026-07-02 (Codex #caveat): gene absent from coords = UNDETERMINED (NA), not silently non-MHC
  as.character(r$chrom[1]) %in% c(MHC_CHR, paste0("chr", MHC_CHR)) && r$start[1] <= MHC_END && r$end[1] >= MHC_START
}

# ── deCODE Download Token ───────────────────────────────────────────────────────
DECODE_TOKEN <- Sys.getenv("DECODE_TOKEN", "")   # NEW-2 fix: from env, NOT hardcoded (scrubbed). Previously-committed token is COMPROMISED — rotate it.

# ── Offline mode (local-first safety net) ──────────────────────────────────────
# Set env CKM_OFFLINE=1 to blank the API CREDENTIALS this project holds (OpenGWAS, deCODE), so a
# stray AUTHENTICATED call fails loud instead of silently phoning home. It is not a network sandbox:
# unauthenticated endpoints, other credentials and ordinary HTTP clients still work. The message it
# used to print claimed any network call would fail, which is not something this file can enforce.
# Affirmative values only -- see ckm_env_yes(). nzchar() is TRUE for "0" and "false", so the two
# spellings that most obviously mean "not offline" both blanked the credentials.
if (ckm_env_yes("CKM_OFFLINE")) {
  OPENGWAS_JWT <- ""; Sys.setenv(OPENGWAS_JWT = "")
  # ...IN THE ENVIRONMENT AS WELL, like the JWT above. Only the R variable was blanked, so any
  # downloader reading Sys.getenv("DECODE_TOKEN") -- or any subprocess inheriting it -- still had a
  # working credential while this line claimed all of them were gone.
  DECODE_TOKEN <- ""; Sys.setenv(DECODE_TOKEN = "")
  cat("** CKM_OFFLINE set: OpenGWAS and deCODE credentials blanked, so any AUTHENTICATED call",
      "fails loudly. Unauthenticated network access is not blocked from here. **\n")
}

cat("Config loaded. Base dir:", BASE_DIR, "\n")
