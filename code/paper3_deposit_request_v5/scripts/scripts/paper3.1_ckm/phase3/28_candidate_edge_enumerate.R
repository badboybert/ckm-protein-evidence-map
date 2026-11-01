# ============================================================================
# Stage 0a (dashed candidate layer): ENUMERATE the candidate-edge pool.
#
# Candidate rule (REQUIREMENTS.md 2.2): candidate = MR full-grid FDR q<0.05 AND
#   coloc.abf dual-rule FAIL = is_coloc_fail() = EVALUATED (finite PP.H4/PP.H3 in [0,1]) AND not
#   passing (PP.H4>=0.80 & PP.H4/PP.H3>=5). NOT the negation of is_coloc_pass(): that returns FALSE
#   for NA, which would read a pair with no ABF probabilities at all as a coloc FAILURE (r37/r47).
#   A pair with NO coloc row is NOT a coloc FAILURE (never tested) -> it is NOT a
#   candidate; such pairs are written to a side file and EXCLUDED from the pool.
#   MHC-region genes (and undetermined-MHC, and no-coords pairs) get susie_engine=NA
#   so 29 never SuSiE-runs them; 30 classifies them without a SuSiE result.
#
# NON-DESTRUCTIVE: reads only frozen corrected artifacts; writes candidate_pool.csv
#   (+ candidate_pool_no_coloc.csv, ALWAYS written -- empty if every FDR-sig pair has a coloc row. It
#    used to be written only when non-empty and never cleared, so a stale one could survive a re-run.)
#
# Run: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase3/28_candidate_edge_enumerate.R
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
suppressPackageStartupMessages(library(data.table))
B <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
source(file.path(B, "R/00_config.R"))   # is_coloc_pass(), is_mhc_gene(), OUTCOMES

RC          <- file.path(B, file.path("results", paste0("repro_", PHASE1_TAG)))
F_COLOC     <- file.path(RC, "coloc_results_v7.rds")
# 2026-07-15: `F_FDR <- .../mr_significant_fdr.rds` was REMOVED, not merely left unread. That artifact is
# the saved by-method one (three-surface block in R/00_config.R) and this script was using it
# while calling the result "MR FDR-significant"; the pool's gate now derives the manuscript's PRIMARY
# rule from F_MR via primary_fdr_set() (R/00_config.R). Requiring a file nothing reads would fail a
# rebuild whose actual inputs are complete, and would tell the next reader this script depends on a
# lineage it does not.
F_BONF      <- file.path(RC, "mr_significant_bonferroni.rds")
F_MR        <- file.path(RC, "mr_all_results.rds")
F_INSTR     <- PHASE1_INSTRUMENTS
OUT         <- file.path(RC, "phase3_susie/candidate_pool.csv")
OUT_NOCOLOC <- file.path(RC, "phase3_susie/candidate_pool_no_coloc.csv")

for (.f in c(F_COLOC, F_BONF, F_MR, F_INSTR))
  if (!file.exists(.f)) stop(sprintf("28: required corrected input missing: %s", .f))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)   # cold-review L1: ensure out dir exists

# CLEAR BOTH OUTPUTS BEFORE ANYTHING CAN ABORT (objective review 2026-08-06, finding 11). The pool
# was only ever replaced by the final fwrite, and the no-coloc sidecar is written earlier, so an
# interruption between the two writes left the PREVIOUS pool on disk beside this run's sidecar --
# and 29 accepts the pool on existence and schema alone. An abandoned run must leave nothing that
# reads as complete.
for (.stale in c(OUT, OUT_NOCOLOC)) if (file.exists(.stale)) {
  if (!file.remove(.stale) || file.exists(.stale))
    stop(sprintf("28: could not remove the previous run's %s -- it would survive a failure here and be read as this run's pool", basename(.stale)))
}

egfr_id <- if (!is.null(OUTCOMES[["eGFR"]])) OUTCOMES[["eGFR"]] else
           stop("28: OUTCOMES$eGFR undefined -- cannot dispatch eGFR engine")

# ── coloc.abf: dual-rule pass (52 solid) vs fail (candidate universe) ────────
# THE RAW OBJECT IS KEPT, because as.data.table() drops the attributes and the provenance check
# below reads them off the stamped result 06 wrote.
coloc_raw <- readRDS(F_COLOC)
coloc <- as.data.table(coloc_raw)
stopifnot(all(c("gene", "outcome_name", "outcome_id", "PP.H4", "PP.H3") %in% names(coloc)))
# A coloc row with no usable probabilities was never TESTED, and the candidate rule below reads
# "not passing" as "tested and failed" -- so an unevaluable row would enter the dashed layer as though
# coloc had returned a negative (r37 High). Reject it here, at the only place that admits pairs.
.une <- coloc[!is_coloc_evaluable(PP.H4, PP.H3)]
if (nrow(.une))
  stop(sprintf("28: %d coloc row(s) have missing or out-of-range ABF probabilities: %s -- these were not evaluated, and the candidate rule would read them as coloc FAILURES. Investigate before enumerating candidates.",
               nrow(.une), paste(utils::head(.une[, sprintf("%s~%s(H4=%s, H3=%s)", gene, outcome_name,
                                             format(PP.H4), format(PP.H3))], 8), collapse = ", ")))
coloc[, coloc_pass := is_coloc_pass(PP.H4, PP.H3)]
.d <- coloc[, .N, by = .(gene, outcome_name)][N > 1]
if (nrow(.d)) {
  chk <- coloc[.d, on = .(gene, outcome_name)][, .(nH4 = uniqueN(round(PP.H4, 6)),
                                                   nH3 = uniqueN(round(PP.H3, 6)),
                                                   nid = uniqueN(outcome_id)), by = .(gene, outcome_name)]
  if (nrow(chk[nH4 > 1 | nH3 > 1 | nid > 1]))
    stop(sprintf("28: conflicting duplicate coloc rows for %s", paste(chk[nH4 > 1 | nH3 > 1 | nid > 1, paste(gene, outcome_name)], collapse = ", ")))
  coloc <- unique(coloc, by = c("gene", "outcome_name"))
}
if (nrow(coloc[, .(n = uniqueN(outcome_id)), by = outcome_name][n > 1]) ||
    nrow(coloc[, .(n = uniqueN(outcome_name)), by = outcome_id][n > 1]))
  stop("28: outcome_name<->outcome_id not 1:1 in coloc_results_v7 (engine dispatch keys on outcome_id)")
# 1:1 is INTERNAL CONSISTENCY, not TRUTH (r21 High, 2026-07-16). If MR and coloc are stale in the SAME
# way -- both carrying a superseded outcome_id for a name -- the check above passes and 29 dispatches its
# engine on an outcome_id config no longer recognises. Assert against OUTCOMES, which is the single source
# for that mapping. Verified 2026-07-16: mr_all_results, coloc_results_v7 and mr_significant_bonferroni
# all already match OUTCOMES exactly, so this costs nothing today.
# assert_outcome_truth() lives in R/00_config.R -- one helper called at every read boundary. A local
# copy here was the first attempt and is exactly how the rule diverges (r22).
assert_outcome_truth(coloc, "coloc_results_v7")
n_solid <- coloc[coloc_pass == TRUE, .N]
cat(sprintf("28: coloc.abf pairs=%d | solid(pass)=%d | fail=%d\n", nrow(coloc), n_solid, coloc[coloc_pass == FALSE, .N]))
# IDENTITY, NOT A COUNT. This demanded exactly 52 solid pairs -- the PREVIOUS generation's frozen
# number -- so on this generation (74) it aborted immediately, and on any other it would have
# aborted too. What it was actually trying to catch is a stale or foreign coloc_results_v7, and 06
# already stamps that file with the keys of the inputs its posteriors were computed from. Requiring
# those (as 26:61 does) tests the thing the count was standing in for. Only PRESENCE is required:
# 06 verified each stratum against the live inputs, and re-deriving them here adds no guarantee.
# PRESENCE IS NOT CURRENCY (objective review 2026-08-06, finding 10). Requiring only that the two
# stamps EXIST accepts a combined file whose stamps came from an earlier successful 06 while its
# CONTENT came from an interrupted or resumed one -- a previous strong set fused with current MR and
# instrument data, with nothing in the output naming the mixture. The two stage artifacts each carry
# the key their own stratum was computed under, so equality is checkable for free, off attributes
# already on disk; no input is re-hashed.
.stage <- c(coloc_inputs_key_noneGFR = file.path(RC, "coloc_results_v7_stage1.rds"),
            coloc_inputs_key_eGFR    = file.path(RC, "coloc_egfr_results.rds"))
for (.a in names(.stage)) {
  .v <- attr(coloc_raw, .a, exact = TRUE)
  if (is.null(.v) || !nzchar(as.character(.v)[1]))
    stop(sprintf("28: %s carries no %s, so the inputs its posteriors were computed from cannot be established. It predates the stamping, or it came from another generation. Re-run cascade/06_combine_coloc.R.", basename(F_COLOC), .a))
  if (!file.exists(.stage[[.a]]))
    stop(sprintf("28: cannot verify %s: the stage artifact %s is absent, so the combined file's stamp cannot be checked against the stratum it claims to come from.", .a, basename(.stage[[.a]])))
  .sv <- attr(readRDS(.stage[[.a]]), "coloc_inputs_key", exact = TRUE)
  if (is.null(.sv) || !identical(as.character(.sv)[1], as.character(.v)[1]))
    stop(sprintf("28: %s = %s does not match %s (%s). The combined file was not built from the stage artifact now on disk -- an interrupted or resumed 06 fuses a previous strong set with current inputs. Re-run cascade/06_combine_coloc.R.",
                 .a, substr(as.character(.v)[1], 1, 12), basename(.stage[[.a]]),
                 if (is.null(.sv)) "unstamped" else substr(as.character(.sv)[1], 1, 12)))
}
cat(sprintf("28: coloc provenance verified against both stage artifacts (%s / %s)
",
            substr(as.character(attr(coloc_raw, "coloc_inputs_key_noneGFR", exact = TRUE))[1], 1, 12),
            substr(as.character(attr(coloc_raw, "coloc_inputs_key_eGFR", exact = TRUE))[1], 1, 12)))

# ── FDR-significant MR set (carry outcome_id for a provenance cross-check; cold-review M) ──
# 2026-07-15: this read `mr_significant_fdr.rds` and called the result "MR FDR-significant". That
# artifact is BH within each method family over analysable rows only, before the primary dedup
# (cascade/01_mr_screen.R's `fdr_sig` line) -- NOT the primary rule, and NOT the same computation as the
# sensitivity value the manuscript prints either (three distinct surfaces; see the block above
# primary_fdr_set() in R/00_config.R -- do not restate it here, it would drift). The manuscript's PRIMARY FDR is BH once over the full grid with
# untested cells padded at p=1, and `make_analysis_reports.R:174` reports the two side by side,
# labelled. This script was using the saved by-method rule [2] while describing it as the primary one [1].
# `primary_fdr_set()` (R/00_config.R) now OWNS that definition and is what this script applies. It is
# NOT yet the definition's only implementation: the reporting layer still re-implements it
# (make_analysis_reports.R:52, make_draft_manuscript_redesigned.py:84, fig1/fig2/volcano/discovery_forest),
# and cascade/14_assign_tiers.R still reads the saved by-method artifact [2]. Both are tracked, with the proof
# that 14's copy is inert, in REBUILD_LEDGER §2.1c. Do not read "owner" as "sole" (objective review r5).
#
# EXPECTED effect, NOT a result: comparing the two definitions on the frozen data (2026-07-15) predicts
# this pool goes 151 -> 142, that the 9 leavers are pairs the current audit already classifies
# `excluded` (8) / `not_evaluable` (1), and that all 17 candidate dashed edges stay FDR-significant so
# 32's rule-1 gate and the 6/3/8 set do not move. That comparison was made against the CURRENT artifacts
# by an independent path -- it is what this script's re-run must CONFIRM, not something it has shown.
# Until the chain is re-run, candidate_pool.csv on disk is still the old 151, gated by [2].
# The chain is 28 -> 29 (BOTH modes) -> 30 -> 31 (BOTH modes) -> 32 -> 33. Neither 29 nor 31 is
# optional: narrowing the pool orphans the leavers' rows in the existing candidate_susie_*.csv and 30's
# orphan guard hard-stops; and 31 reads the audit 30 rewrites, so its sensitivity CSVs would then
# predate their own input (33's dependency gate refuses to certify that).
mr <- as.data.table(readRDS(F_MR))          # read once; the primary FDR and the pm block below share it

# ── SECOND-ORDER (DELTA) SE ON THIS GATE (PI decision 2026-08-06) ────────────
# The plan adopts z2 = z_Y / sqrt(1 + z_Y^2/F) and requires it applied CONSISTENTLY to the entry
# gates, not only to the headline table. THIS POOL IS AN ENTRY GATE -- it is *defined* by the FDR
# set -- so leaving it on first-order p would admit pairs the adopted rule excludes. That is why
# "the set only shrinks, so nothing needs re-running" does NOT carry here: for L1/L2 a shrinking
# set is a filter over work already done, but a pool that is smaller is a pool that was defined
# differently, and everything downstream of it (29 -> 30 -> 31 -> 32) is keyed on that definition.
# Measured on this generation: primary FDR 292 -> 248 (44 leave, 0 join); 36 of the previous
# pool's 213 pairs leave with them.
#
# THE COPY IS DELIBERATE. Only the FDR gate moves to second-order; `mr` keeps its first-order
# `pval` because the Bonferroni cross-check below compares the derived set against
# mr_significant_bonferroni.rds, a first-order artifact -- rewriting `pval` in place would make
# that assertion fire on a difference this block introduced rather than on a stale file.
#
# SINGLE-INSTRUMENT ONLY. The identity rests on se_Wald = se_Y/|b_X|, hence z_MR = z_Y, which holds
# for one instrument. Multi-instrument rows have no agreed second-order rule yet (see
# IMPLEMENTATION_SPEC "多工具規則", options (i)/(ii)); they stay on first-order here, which IS
# option (i). The counts are printed rather than implied, and a gene that reached here without an
# instrument F would also fall through to first-order, so it is counted separately instead of
# hiding inside the multi-instrument tally.
instr  <- as.data.table(readRDS(F_INSTR))
.ftab  <- instr[, .(F_stat = F_stat[1]), by = gene]
mr_fdr <- merge(copy(mr), .ftab, by = "gene", all.x = TRUE)
.single <- mr_fdr$method %in% "Wald ratio"
.usable <- .single & is.finite(mr_fdr$se) & is.finite(mr_fdr$b) & is.finite(mr_fdr$F_stat) & mr_fdr$F_stat > 0
mr_fdr[.usable, pval := 2 * pnorm(-abs(b / sqrt(se^2 + b^2 / F_stat)))]
cat(sprintf("28: second-order SE on %d single-instrument row(s); first-order kept for %d multi-instrument row(s) [option (i)] and %d single-instrument row(s) lacking a usable F\n",
            sum(.usable), sum(!.single), sum(.single & !.usable)))
if (!sum(.usable))
  stop("28: the second-order correction matched no row -- the gate would silently stay first-order")
mr_fdr[, F_stat := NULL]
fdr <- primary_fdr_set(mr_fdr)
stopifnot(all(c("gene", "outcome_name", "outcome_id") %in% names(fdr)))
fdr_keys <- unique(fdr[, .(gene, outcome_name, outcome_id_fdr = outcome_id)])
if (anyDuplicated(fdr_keys, by = c("gene", "outcome_name")))
  stop("28: FDR set maps a (gene,outcome_name) to >1 outcome_id -- ambiguous")
cat(sprintf("28: MR FDR-significant pairs = %d\n", nrow(fdr_keys)))

# is_bonferroni annotates dashed-edge provenance, and reading it from a saved artifact means a stale
# file would mislabel it silently. `mr` is already in memory, so DERIVE the set and use the artifact only
# as an external cross-check -- the same shape as rule (1) in 32, and for the same reason: a guard whose
# input is downstream of what it verifies proves nothing (r21 Medium, 2026-07-16).
# BONFERRONI_ALPHA is bonferroni_alpha(N_PROTEINS_CANONICAL) -- the same full-grid denominator
# primary_fdr_set() uses. Verified 2026-07-16: saved 75 pairs, derived 75, setequal TRUE.
bf <- as.data.table(readRDS(F_BONF))
assert_outcome_truth(bf, "mr_significant_bonferroni")
# Two-column keyed sets, not paste()-strings. `paste(gene, outcome_name)` collapses a composite key into
# a string that a space in either field would make ambiguous, and the comparison below would then run over
# the SAME lossy representation -- it could not detect its own collision. The frozen vocabulary has no such
# names, so this is latent, not live; it is fixed because these lines were open anyway (r24).
.bonf_saved   <- unique(bf[, .(gene, outcome_name)])
.bonf_derived <- unique(mr[method %in% c("Wald ratio", "Inverse variance weighted") & pval < BONFERRONI_ALPHA,
                           .(gene, outcome_name)])
.only_saved   <- fsetdiff(.bonf_saved, .bonf_derived)
.only_derived <- fsetdiff(.bonf_derived, .bonf_saved)
if (nrow(.only_saved) || nrow(.only_derived))
  stop(sprintf("28: mr_significant_bonferroni.rds (%d pairs) is not set-equal to the Bonferroni set derived from mr_all_results.rds at BONFERRONI_ALPHA=%.4e (%d pairs). saved-only: %s | derived-only: %s -- one of them is stale; refusing to annotate is_bonferroni from a file that disagrees with the data.",
               nrow(.bonf_saved), BONFERRONI_ALPHA, nrow(.bonf_derived),
               paste(utils::head(.only_saved[, paste0(gene, "~", outcome_name)], 8), collapse = ", "),
               paste(utils::head(.only_derived[, paste0(gene, "~", outcome_name)], 8), collapse = ", ")))
# Build from the DERIVED set, not from `bf`. The cross-check above makes the two identical on a passing
# run, so this is inert today -- but the comment claimed the artifact was "only" a cross-check while the
# code still took its value from it (r22 Low). Say what you do; do what you say.
bonf_keys <- copy(.bonf_derived)

pm <- mr[method %in% c("Wald ratio", "Inverse variance weighted")]   # `mr` read above
# conflict guard mirroring 12_rebuild_pair_evidence.R:180 (P4-F02): if a pair has >1 primary-method row
# disagreeing on b/pval, DO NOT resolve by row order -- an arbitrary Wald-vs-IVW pick could flip the frozen
# risk/protective direction downstream in 30. One primary row (Wald XOR IVW) per pair is the norm, so this is
# a tripwire, not an expected stop (Codex review Medium: silent first-row MR pick).
.pmc <- pm[, .(nb = uniqueN(b), np = uniqueN(pval), .N), by = .(gene, outcome_name)][N > 1 & (nb > 1 | np > 1)]
if (nrow(.pmc))
  stop(sprintf("28: conflicting duplicate primary MR rows for %s -- ambiguous, refusing first-row pick (mirror 12_rebuild P4-F02)", paste(.pmc[, paste(gene, outcome_name)], collapse = ", ")))
pm <- pm[!duplicated(paste(gene, outcome_name))]   # rows agree (or single) -> first is safe; method preference moot
pm <- pm[, intersect(c("gene", "outcome_name", "b", "pval", "se"), names(pm)), with = FALSE]

# THE POOL MUST CARRY THE p THAT ADMITTED THE ROW (objective review 2026-08-06, finding 7).
# `pm` comes from the untouched first-order `mr`, so the emitted `pval` was the SUPERSEDED
# statistic while membership was decided on the second-order one -- e.g. PMVK~AF enters on
# 8.238e-08 and the pool reported 1.892e-31. Membership is unaffected; the number printed beside it
# was not the number that put it there, and nothing said so. Both are now emitted under names that
# say which is which, and the gate's own statistic leads.
setnames(pm, "pval", "pval_first_order")
pm <- merge(pm, unique(mr_fdr[method %in% c("Wald ratio", "Inverse variance weighted"),
                              .(gene, outcome_name, pval = pval)])[!duplicated(paste(gene, outcome_name))],
            by = c("gene", "outcome_name"), all.x = TRUE)
if (anyNA(pm$pval))
  stop(sprintf("28: %d pair(s) have no gate p-value after the join -- the pool would report a p that did not admit them", sum(is.na(pm$pval))))

gcoord <- unique(instr[, .(gene, gene_chrom, gene_start, gene_end)])   # `instr` already read for the F table above
if (anyDuplicated(gcoord$gene)) stop("28: a gene has >1 coordinate row in instruments -- ambiguous region window")

# ── candidate = FDR-significant joined to coloc ──────────────────────────────
cand <- merge(fdr_keys, coloc[, .(gene, outcome_name, outcome_id, PP.H4_abf = PP.H4,
                                   PP.H3_abf = PP.H3, coloc_pass)],
              by = c("gene", "outcome_name"), all.x = TRUE)
mism <- cand[!is.na(outcome_id) & !is.na(outcome_id_fdr) & outcome_id != outcome_id_fdr]
if (nrow(mism)) stop(sprintf("28: FDR outcome_id != coloc outcome_id for %s -- provenance mismatch", paste(mism[, paste(gene, outcome_name)], collapse = ", ")))

# FDR-significant pairs with NO coloc row are NOT coloc-failures -> side file, EXCLUDE from pool (cold-review 28-H1)
nocoloc <- cand[is.na(coloc_pass)]
# 2026-07-15: this used to write the file ONLY when non-empty, and never removed it otherwise -- so a
# stale sidecar from an earlier run could survive a complete re-run and be read as this run's exclusion
# audit. That matters right now: narrowing the FDR gate to the primary rule changes which pairs reach
# this branch, so a leftover file would describe a pool that no longer exists. Always write it, empty or
# not: an empty file is a result ("no FDR-significant pair lacked a coloc row"), a missing or stale one
# is an unanswered question.
fwrite(nocoloc[, .(gene, outcome_name, outcome_id = outcome_id_fdr)], OUT_NOCOLOC)
cat(if (nrow(nocoloc))
      sprintf("28: [WARN] %d FDR-significant pair(s) have NO coloc.abf row -- NOT candidates (a missing test != a failed test). Written to %s and EXCLUDED from the pool.\n", nrow(nocoloc), basename(OUT_NOCOLOC))
    else
      sprintf("28: every FDR-significant pair has a coloc.abf row; %s written empty (0 rows).\n", basename(OUT_NOCOLOC)))
cand <- cand[coloc_pass == FALSE]          # pool = strictly coloc-FAILED
cand[, `:=`(coloc_status = "coloc_fail", outcome_id_fdr = NULL)]

# ── annotate ─────────────────────────────────────────────────────────────────
cand <- merge(cand, pm, by = c("gene", "outcome_name"), all.x = TRUE)
cand <- merge(cand, gcoord, by = "gene", all.x = TRUE)
cand[, is_mhc := vapply(gene, function(g) { m <- is_mhc_gene(g); if (length(m) != 1) NA else as.logical(m) }, logical(1))]  # TRUE/FALSE/NA preserved
# Keyed join, not paste()-strings. r24 converted the Bonferroni CHECK and left the EMISSION on string
# keys -- the check was then keyed while the label it protects was not (r25 Medium).
cand <- merge(cand, copy(bonf_keys)[, .bonf := TRUE], by = c("gene", "outcome_name"), all.x = TRUE)
cand[, is_bonferroni := !is.na(.bonf)][, .bonf := NULL]
cand[, is_egfr    := (!is.na(outcome_id) & outcome_id == egfr_id)]
cand[, has_coords := !(is.na(gene_chrom) | is.na(gene_start) | is.na(gene_end))]
# testable = coloc-fail, definitively non-MHC, full coords. Only these get a SuSiE engine (cold-review 28-H2/29-H1):
# MHC(TRUE), undetermined-MHC(NA), and no-coords pairs get engine=NA so 29 never runs them.
cand[, is_testable  := (!is.na(is_mhc) & is_mhc == FALSE & has_coords == TRUE)]
cand[, susie_engine := fifelse(is_testable, fifelse(is_egfr, "egfr", "nonegfr"), NA_character_)]

# row-grain assertion after all joins (cold-review 28-M2)
if (anyDuplicated(cand, by = c("gene", "outcome_name")))
  stop("28: candidate pool not unique by (gene,outcome_name) after MR/coord joins -- refusing to freeze a fanned-out pool")

setorder(cand, susie_engine, -PP.H4_abf, gene, outcome_name, na.last = TRUE)
fwrite(cand, OUT)

cat(sprintf("\n28: candidate pool = %d coloc-fail pairs | testable(SuSiE)=%d [%d eGFR, %d non-eGFR] | MHC=%d MHC-NA=%d no-coords=%d\n",
            nrow(cand), cand[is_testable == TRUE, .N],
            cand[is_testable == TRUE & is_egfr == TRUE, .N], cand[is_testable == TRUE & is_egfr == FALSE, .N],
            cand[is_mhc == TRUE, .N], cand[is.na(is_mhc), .N], cand[has_coords == FALSE, .N]))
cat(sprintf("28: wrote %s\n", OUT))
cat("\n== coloc-fail candidates by outcome ==\n"); print(cand[, .N, by = outcome_name][order(-N)])
cat("\n== top 15 by coloc.abf PP.H4 (closest-to-pass first) ==\n")
print(head(cand[order(-PP.H4_abf), .(gene, outcome_name, PP.H4_abf, PP.H3_abf, is_mhc, susie_engine)], 15))
