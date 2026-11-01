# ============================================================================
# Analysis B (must-add): colocalization specificity/sensitivity ROBUSTNESS TABLES  [rev3]
# ----------------------------------------------------------------------------
# LOCATION: moved 2026-07-15 (PI) from `reporting/coloc_robustness_table.R` into phase3 as
# `34_coloc_robustness_table.R` -- it is a phase-3 colocalization analysis, not a reporting-layer
# projection, and belongs with 28-33 (candidate enumerate/susie/audit/sensitivity/finalize/provenance).
# Paths below are ABSOLUTE, so the move changes no behaviour. RUN ORDER (L8): after
# build_table5_prioritisation.py, BEFORE build_supptables_from_sources.py (whose manifest
# requires this script's outputs) and verify_supptables.py.
# ----------------------------------------------------------------------------
# Pure ASSEMBLY (no new statistics) of frozen tables into TWO supplementary tables that
# make the two-sided coloc calibration explicit:
#   CORE (the ABF-coloc pairs): ABF PP.H4 vs SuSiE PP.H4 (concordance, MHC-flagged) and
#     PP.H4 across 3 coloc priors (p12 = 1e-6/5e-6/1e-5) + prior stability + signal structure.
#   DASHED (the SuSiE-rescued candidates): ABF PP.H4 (mostly ~0) vs SuSiE PP.H4 (orig) vs
#     SuSiE PP.H4 under the STRICT re-fit + survival + raw run_status/warnings + draw_status.
# Deterministic frozen sources (no assign_evidence_tiers ambient-global dependency):
#   table_s_pair_evidence_tiered.csv  (canonical ROOT; tiers, ABF, stable)
#   phase3_susie/susie_coloc_results.csv  (ABF vs SuSiE, mhc_excluded, signal counts)
#   prior_sweep_wide.csv  (3-prior sweep + PP.H4_original)
#   phase3_susie/dashed_edges_final.csv  (the dashed candidates)
# Integrity-gated by STRUCTURE, not by this generation's counts (2026-09-22): the three core
# sources must have identical, collision-free key sets; the joins must be one-to-one and must not
# change the row count; the tier values must be exactly the three known levels; core and dashed
# must not overlap; draw_status must be one of the four known levels AND must agree with the raw
# run fields per status; cross-source ABF identity to 4 dp. The row counts, the tier distribution
# and the draw_status distribution used to be pinned here as literals (76 / 25 / 3-41-32 / 15-1-5-4,
# the 2026-08-16 hg38 anchor). Every one of them was REDUNDANT with a structural check already in
# this file -- key uniqueness plus setequal already forces the three sources to the same row count,
# and the per-status field invariants already test what the draw_status counts were standing in for
# -- so all they added was an abort whenever the generation legitimately moved, which is what they
# did on 2026-09-22 (81 pairs / 62 genes / 3-45-33 / 26 dashed). The Convergent IDENTITIES below
# stay: those are a manuscript claim, not a count.
# Writes 3 CSVs (OUT must pre-exist; script never creates directories).
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
suppressPackageStartupMessages(library(tidyverse))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))   # L8r1 finding 1: PHASE1_TAG was undefined standalone
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))     # COLOC_H4H3_MIN for the reported (non-gating) SuSiE ratio flag -- the single-variant rule's constant, not a copy of it
RC  <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), paste0("repro_", PHASE1_SET))   # L8r2 #1: _phase1_set.R defines PHASE1_SET (PHASE1_TAG was undefined standalone)
OUT <- Sys.getenv("SUPP_OUT", paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/reports/Paper3.1_CKM_Evidence_Map/supp_tables"))
stopifnot("output dir does not exist (refusing to create it)" = dir.exists(OUT))
TIER_DISPLAY <- c(Convergent = "Convergent", Supported = "Supported", Flagged = "Flagged")  # MUST match R/theme_paper3.R
key   <- function(g, o) paste(g, o, sep = "\x1f")                # unit-separator -> collision-free keys (#6)
same4 <- function(a, b) isTRUE(all(round(a, 4) == round(b, 4)))  # identity after 4-dp rounding, not a loose tolerance (#7)

# ---- frozen sources ----
tiered <- read_csv(file.path(RC, "table_s_pair_evidence_tiered.csv"), show_col_types = FALSE) %>%   # obj7A: canonical ROOT (was _fig_results); avoids mixing staged tiers with root SuSiE/prior sources
  transmute(gene, outcome, pp_h4_abf = pp_h4, tier = as.character(evidence_tier),
            t_stable = stable, t_pp_h4_1e6 = PP.H4_p12_1e6,
            # L8r17 #2: the CANONICAL four-class SuSiE category (supportive / inconclusive /
            # discordant / not_evaluable). Re-deriving it here from concordant_susie counted
            # non-converged fits (ADM-AF, NCAN-CAD) as discordant although the pipeline calls
            # them not_evaluable, so this file disagreed with every other consumer.
            susie_class = as.character(susie_class))
susie <- read_csv(file.path(RC, "phase3_susie/susie_coloc_results.csv"), show_col_types = FALSE) %>%
  transmute(gene, outcome = outcome_name, s_pp_h4_abf = PP.H4_abf, pp_h4_susie = PP.H4_susie,
            pp_h3_susie = PP.H3_susie, susie_concordant_abf = concordant_susie, mhc_excluded,
            n_signals_exp, n_signals_out, n_signal_pairs)
prior <- read_csv(file.path(RC, "prior_sweep_wide.csv"), show_col_types = FALSE) %>%
  transmute(gene, outcome = outcome_name, prior_pp_h4_1e6 = PP.H4_p12_1e6, prior_pp_h4_5e6 = PP.H4_p12_5e6,
            prior_pp_h4_1e5 = PP.H4_p12_1e5, p_pp_h4_orig = PP.H4_original, prior_stable = stable)

# coverage / uniqueness / bidirectional key equality (#6)
stopifnot("no core pairs at all" = nrow(tiered) > 0L,
          "tiered keys not unique" = !any(duplicated(key(tiered$gene, tiered$outcome))),
          "susie keys not unique"  = !any(duplicated(key(susie$gene, susie$outcome))),
          "prior keys not unique"  = !any(duplicated(key(prior$gene, prior$outcome))),
          "tiered vs susie key set differs" = setequal(key(tiered$gene, tiered$outcome), key(susie$gene, susie$outcome)),
          "tiered vs prior key set differs" = setequal(key(tiered$gene, tiered$outcome), key(prior$gene, prior$outcome)))
.tc <- tiered %>% count(tier)
stopifnot("tier values are not exactly the three known levels" =
            setequal(.tc$tier, c("Convergent","Supported","Flagged")))

core_tbl <- tiered %>%
  left_join(susie, by = c("gene","outcome"), relationship = "one-to-one") %>%
  left_join(prior, by = c("gene","outcome"), relationship = "one-to-one")
stopifnot("core join changed row count" = nrow(core_tbl) == nrow(tiered))
# cross-source ABF + stable + prior identity (#7): the three sources must describe the same coloc
stopifnot("ABF PP.H4 disagrees tiered vs susie"     = same4(core_tbl$pp_h4_abf, core_tbl$s_pp_h4_abf),
          "ABF PP.H4 disagrees tiered vs prior"      = same4(core_tbl$pp_h4_abf, core_tbl$p_pp_h4_orig),
          "prior-stable flag disagrees tiered vs prior" = identical(core_tbl$t_stable, core_tbl$prior_stable),
          "prior PP.H4_1e6 disagrees tiered vs prior" = same4(core_tbl$t_pp_h4_1e6, core_tbl$prior_pp_h4_1e6))

core_out <- core_tbl %>% transmute(
  gene, outcome, tier = unname(TIER_DISPLAY[tier]),
  pp_h4_abf, pp_h4_susie, pp_h3_susie, susie_concordant_abf, susie_class, mhc_excluded,
  n_signals_exp, n_signals_out, n_signal_pairs,
  prior_pp_h4_1e6, prior_pp_h4_5e6, prior_pp_h4_1e5, prior_stable,
  note = case_when(
    mhc_excluded %in% TRUE                                      ~ "MHC-excluded (coloc call not SuSiE-driven)",
    is.na(pp_h4_susie)                                          ~ "SuSiE not evaluable",
    susie_concordant_abf %in% TRUE & prior_stable %in% TRUE     ~ "robust: ABF+SuSiE concordant, prior-robust",
    !(prior_stable %in% TRUE) & prior_pp_h4_1e6 >= 0.50        ~ "prior-attenuated coloc (0.50-0.80)",
    !(prior_stable %in% TRUE) & prior_pp_h4_1e6 <  0.50        ~ "prior-fragile coloc (<0.50)",
    !(susie_concordant_abf %in% TRUE)                           ~ "ABF/SuSiE discordant",
    TRUE                                                        ~ "supported"))

# ---- 16 dashed candidates ----
DRAW_LEVELS <- c("drawn", "held:non_converged", "dropped:collapsed", "dropped:strict_fail")
dashed <- read_csv(file.path(RC, "phase3_susie/dashed_edges_final.csv"), show_col_types = FALSE)
.cnt <- table(factor(dashed$draw_status, levels = DRAW_LEVELS))
# 2026-07-18: dashed counts corrected 17->16 (6/3/8 -> 7/1/8). The 2026-07-12 6/3/8 sign-off used max_iter=100
# SuSiE fits; the max_iter=1000 rebuild converged them, revealing ATP5IF1 as a clean drawn edge and DPEP1's
# 0.95 as a NON-CONVERGED artifact (converged PP.H4=0.08 -> excluded, no longer a dashed candidate). Verified
# against dashed_edges_final.csv; manuscript prose reconciled to 7/1/8 (DPEP1 excluded, not held). LOCKED.md
# carries a dated 2026-07-18 correction note; its verbatim body still awaits the editor's formal re-sign.
stopifnot("no dashed candidates at all" = nrow(dashed) > 0L,
          "dashed keys not unique" = !any(duplicated(key(dashed$gene, dashed$outcome_name))),
          "dashed overlaps the core set" = length(intersect(key(dashed$gene, dashed$outcome_name), key(tiered$gene, tiered$outcome))) == 0L,
          "dashed has an unexpected draw_status" = all(dashed$draw_status %in% DRAW_LEVELS))
# per-status field invariants (#8): draw_status must be consistent with the raw run fields
dr <- dashed$draw_status
stopifnot("drawn invariant (ok + survives + converged)" =
            all(dashed$run_status[dr=="drawn"]=="ok" & dashed$survives[dr=="drawn"]==TRUE & dashed$n_warn[dr=="drawn"]==0L),
          "held:non_converged invariant (ok + survives + exactly one warning)" =
            all(dashed$run_status[dr=="held:non_converged"]=="ok" & dashed$survives[dr=="held:non_converged"]==TRUE & dashed$n_warn[dr=="held:non_converged"]==1L),
          "dropped:collapsed invariant (collapsed + not survive)" =
            all(dashed$run_status[dr=="dropped:collapsed"]=="collapsed" & dashed$survives[dr=="dropped:collapsed"]==FALSE),
          "dropped:strict_fail invariant (ok + not survive)" =
            all(dashed$run_status[dr=="dropped:strict_fail"]=="ok" & dashed$survives[dr=="dropped:strict_fail"]==FALSE))

# RATIO COLUMN (PI ruling 2026-08-01, condition of record): the dashed layer deliberately applies
# PP.H4_SuSiE >= 0.80 ALONE, with no PP.H4/PP.H3 ratio arm -- SuSiE's H3 is per-credible-set-pair and
# the single-variant ratio does not carry over. That decision is only defensible if the ratio it
# declines to apply is VISIBLE for every row, including drawn edges below the single-variant rule.
# Reported, never gating.
dashed_out <- dashed %>% transmute(
  gene, outcome = outcome_name, draw_status,
  pp_h4_abf = PP.H4_abf, pp_h4_susie_orig = orig_PP.H4_susie, pp_h4_susie_strict = strict_PP.H4_susie,
  pp_h3_susie_strict = strict_PP.H3_susie,
  susie_h4_h3_ratio_strict = if_else(strict_PP.H3_susie > 0,
                                     strict_PP.H4_susie / strict_PP.H3_susie, Inf),
  susie_ratio_below_single_variant_rule = if_else(strict_PP.H3_susie > 0,
                                     strict_PP.H4_susie < COLOC_H4H3_MIN * strict_PP.H3_susie, FALSE),
  survives, run_status, n_warn, ld_warn,
  note = paste0(draw_status, if_else(PP.H4_abf < 0.02,
                if_else(draw_status == "drawn", " | ABF-negative (PP.H4<0.02; SuSiE-only rescue)", " | ABF-negative (PP.H4<0.02)"),   # obj7A4#4: only DRAWN rows are rescues; held/dropped are not
                "")))

# obj7A: the summary/counts read these fields; assert non-missing so na.rm cannot silently drop a row
stopifnot("dashed pp_h4_abf / survives must be non-missing" =
            !any(is.na(dashed_out$pp_h4_abf)) && !any(is.na(dashed_out$survives)))

# ---- summary ----
summ <- tibble(
  metric = c("core_pairs", "core_robust", "core_prior_stable", "core_mhc_excluded",
             "core_susie_not_evaluable", "core_abf_susie_discordant", "core_susie_inconclusive",
             "dashed_edges", "dashed_drawn", "dashed_held", "dashed_dropped",
             "dashed_ABF_negative_lt0.02", "dashed_survives_strict_refit"),
  value  = c(nrow(core_out),
             sum(core_out$note == "robust: ABF+SuSiE concordant, prior-robust"),
             sum(core_out$prior_stable %in% TRUE),
             sum(core_out$mhc_excluded %in% TRUE),
             sum(core_out$susie_class == "not_evaluable"),
             sum(core_out$susie_class == "discordant"), sum(core_out$susie_class == "inconclusive"),   # L8r15 #1 + L8r17 #2: not the note field (its
             # prior-band branches fire first and masked this), and not raw concordant_susie=FALSE
             # (that lumps non-converged fits in with genuine discordance). The canonical class.
             nrow(dashed_out),
             .cnt[["drawn"]], .cnt[["held:non_converged"]],
             .cnt[["dropped:collapsed"]] + .cnt[["dropped:strict_fail"]],
             sum(dashed_out$pp_h4_abf < 0.02, na.rm = TRUE),
             sum(dashed_out$survives == TRUE, na.rm = TRUE)))

write_csv(core_out,   file.path(OUT, "SuppTable_coloc_robustness_core.csv"))
write_csv(dashed_out, file.path(OUT, "SuppTable_coloc_robustness_dashed.csv"))
write_csv(summ,       file.path(OUT, "coloc_robustness_summary.csv"))   # audit sidecar (NOT a SuppTable_ -> stays out of the six-table verify manifest)

cat("=== Analysis B: colocalization specificity/sensitivity robustness ===\n")
cat(sprintf("core rows: %d | dashed rows: %d\n", nrow(core_out), nrow(dashed_out)))
print(as.data.frame(summ))
cat("\nWrote SuppTable_coloc_robustness_{core,dashed}.csv + coloc_robustness_summary.csv ->", OUT, "\n")
