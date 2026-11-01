# Phase 2-Sens Extension: p < 5e-8 sensitivity for eGFR strong hits
# ============================================================================
# The original 15_phase2_sens_5e8.R excludes eGFR and (under the corrected wrapper) derives the non-eGFR
# strong set from P2S_STRONG_FILE. This script runs eGFR-specific p<5e-8 MR and appends the eGFR pairs
# (the current table has 16 eGFR pairs; 14 retain a 5e-8 instrument) to the concordance table.
#
# Usage: cd "<E>/20260331-CKM/Paper 3" && Rscript scripts/paper3.1_ckm/phase2/16_sens_5e8_egfr.R
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
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_local_outcome_extract.R"))   # local_extract_outcome_data + .LO_META

cat("=== Phase 2-Sens Extension: eGFR p < 5e-8 Sensitivity ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# NOT INVOKED DIRECTLY. 50c_sens5e8_egfr_offline.R sets RESULTS_DIR to this generation's
# _fig_results BEFORE sourcing R/00_config.R, so config adopts it and RESULTS_DIR/PHASE2_DIR
# already resolve under the generation here. A previous edit re-derived a generation directory
# from RESULTS_DIR, which under 50c produced .../_fig_results/repro_<tag> -- a path that does not
# exist -- and turned a working read into a hard stop. Leave the resolution to the caller.
out_dir <- file.path(RESULTS_DIR, "phase2_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load instruments and filter to p < 5e-8 ─────────────────────────────────
# THE TABLE CONFIG CERTIFIED, not a second read of the same path -- this script applies config's
# BONFERRONI_ALPHA to what it loads here, so the two must be the same bytes.
instr <- read_phase1_instruments()
setDT(instr)

instr_5e8 <- instr[pval < 5e-8]
cat(sprintf("Instruments p<5e-8: %d SNPs, %d proteins\n\n",
            nrow(instr_5e8), uniqueN(instr_5e8$gene)))

# ── eGFR strong hits (16 in the current table; derived from P2S_EGFR_FILE, else fail-loud) ─────
# DERIVED WHEN NOT SUPPLIED, rather than demanded. The refusal below is right about WHAT it refuses
# -- the stale hardcoded 13-gene list -- but it also made the Usage line above unrunnable: a direct
# Rscript launch cannot set an in-process variable, so this script always stopped. The pair table of
# the generation being run is the same thing a caller would have passed, so it is derived from
# PHASE1_TAG here; if that file is absent the refusal still stands, because then there genuinely is
# nothing to derive the eGFR strong set from.
if (!exists("P2S_EGFR_FILE")) {
  .p2s_default <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3"), "results", paste0("repro_", PHASE1_TAG),
                            "_fig_results", "phase4_annotation", "table_s_pair_evidence.csv")
  if (file.exists(.p2s_default)) {
    P2S_EGFR_FILE <- .p2s_default
    cat(sprintf("16_sens_5e8_egfr: P2S_EGFR_FILE derived from this generation -> %s\n", P2S_EGFR_FILE))
  }
}
if (exists("P2S_EGFR_FILE")) {            # A6: derive eGFR strong-pair genes from corrected table_s (16 in the current tableER)
  .ts <- data.table::fread(P2S_EGFR_FILE)
  egfr_hits <- data.table(gene = sort(unique(.ts[outcome == "eGFR", gene])), outcome_name = "eGFR")
} else {
  # P2-F2 (cold audit 2026-07-08): the hardcoded 13-gene eGFR fallback is STALE vs the corrected 14-eGFR set
  # (missing DNER-eGFR). A DIRECT run without P2S_EGFR_FILE (its Usage header permits it) would silently drop
  # DNER-eGFR from the p<5e-8 sensitivity table -> fail loud instead.
  stop("16_sens_5e8_egfr: P2S_EGFR_FILE not set -- refusing the stale hardcoded 13-gene eGFR fallback (cold-audit P2-F2); the corrected wrapper must set P2S_EGFR_FILE = the current table_s_pair_evidence.csv (14 eGFR pairs incl DNER).")
}
cat(sprintf("eGFR strong hits to check: %d\n", nrow(egfr_hits)))

# Check instrument availability
genes_with_5e8 <- unique(instr_5e8$gene)
egfr_hits[, has_instrument_5e8 := gene %in% genes_with_5e8]
cat(sprintf("With p<5e-8 instrument: %d / %d\n\n",
    sum(egfr_hits$has_instrument_5e8), nrow(egfr_hits)))

# ── Format exposure data (only genes in eGFR hits with p<5e-8 instruments) ──
egfr_genes_5e8 <- egfr_hits[has_instrument_5e8 == TRUE, gene]
instr_sub <- instr_5e8[gene %in% egfr_genes_5e8]

if (nrow(instr_sub) == 0) {
  cat("No eGFR hit genes have p<5e-8 instruments. Nothing to do.\n")
  # Still create empty concordance rows
} else {
  exp_dat <- format_data(
    dat               = as.data.frame(instr_sub),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize",
    phenotype_col     = "exposure"
  )
  cat(sprintf("Exposure data formatted: %d rows, %d genes\n\n",
              nrow(exp_dat), uniqueN(instr_sub$gene)))

  # ── Fetch eGFR outcome data ──────────────────────────────────────────────
  # NO OPENGWAS FALLBACK. "ebi-a-GCST003372" is an API accession; with a local read it is not a
  # weaker default but an id .LO_META cannot resolve, so a config that failed to define eGFR
  # would surface as an obscure extractor error instead of the missing setting it is.
  if (!exists("OUTCOMES") || is.null(OUTCOMES[["eGFR"]]))
    stop("16: OUTCOMES[[\"eGFR\"]] is not set by R/00_config.R -- cannot choose the eGFR outcome")
  EGFR_ID <- OUTCOMES[["eGFR"]]
  cat(sprintf("Fetching eGFR outcome data (%s)...\n", EGFR_ID))
  # LOCAL READ, NOT THE API. extract_outcome_data() is the OpenGWAS call: this generation is
  # offline by rule, so invoked directly this script could never complete. The local extractor
  # takes the same first two arguments; rsq/align_alleles/palindromes/maf_threshold were API
  # proxy-service parameters with no local counterpart -- proxies come from the LD cache.
  out_dat <- tryCatch(
    local_extract_outcome_data(
      snps     = unique(exp_dat$SNP),
      outcomes = EGFR_ID,
      proxies  = TRUE,
      # This script calls the extractor DIRECTLY rather than through the name the offline wrappers
      # override, so 50c's fallthrough-enabled closure never reaches it (p2fix_r2) -- the flag has to
      # be set here or the eGFR arm keeps the pre-alignment rule while 11 and 15 no longer do.
      palindromic_fallthrough = TRUE
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )

  if (!is.null(out_dat) && nrow(out_dat) > 0) {
    if ("n" %in% names(out_dat))                    out_dat$n                    <- suppressWarnings(as.numeric(out_dat$n))
    if ("samplesize.outcome" %in% names(out_dat))   out_dat$samplesize.outcome   <- suppressWarnings(as.numeric(out_dat$samplesize.outcome))
    cat(sprintf("eGFR outcome rows: %d\n\n", nrow(out_dat)))

    # ── Harmonise & MR ───────────────────────────────────────────────────
    # DIRECT AND PROXY ROWS HARMONISED SEPARATELY, mirroring cascade/01_mr_screen.R:758-765 and
    # phase2/11 and /15. A proxy row carries the PROXY variant's frequency under the instrument's
    # allele labels, which action = 2 discards as a palindromic frequency conflict; discovery sets
    # eaf.outcome to NA and sends proxies through action = 1.
    .od      <- as.data.frame(out_dat)
    .is_prox <- if ("proxy.outcome" %in% names(.od)) !is.na(.od$proxy.outcome) & .od$proxy.outcome else rep(FALSE, nrow(.od))
    .hz      <- function(d, act) if (!nrow(d)) NULL else harmonise_data(exp_dat, d, action = act)
    .direct  <- .hz(.od[!.is_prox, , drop = FALSE], 2L)
    .pd      <- .od[.is_prox, , drop = FALSE]
    if (nrow(.pd)) .pd$eaf.outcome <- NA_real_
    .proxy   <- .hz(.pd, 1L)
    harm <- if (is.null(.direct)) .proxy else if (is.null(.proxy)) .direct else
              as.data.frame(data.table::rbindlist(list(.direct, .proxy), fill = TRUE))
    if (is.null(harm) || nrow(harm) == 0)
      stop(sprintf("16_sens_5e8_egfr: harmonisation produced no rows from %d outcome row(s) -- refusing to continue", nrow(.od)))
    harm <- harm[harm$mr_keep, ]
    cat(sprintf("Harmonised: %d rows (direct %d, proxy %d)\n\n",
                nrow(harm), sum(!.is_prox), sum(.is_prox)))

    mr_res <- mr(harm, method_list = c("mr_wald_ratio", "mr_ivw"))
    mr_dt  <- as.data.table(mr_res)
    mr_dt[, gene := trimws(gsub(" \\(.*\\).*", "", exposure))]
    mr_dt[, outcome_name := "eGFR"]
  } else {
    # Codex backfill-review fix: if there ARE eGFR p<5e-8 instruments to test, empty extraction is a silent
    # coverage failure -> fail loud rather than emit an all-NA concordance file.
    if (exists("exp_dat") && nrow(exp_dat) > 0)
      stop(sprintf("16_sens_5e8_egfr: %d eGFR p<5e-8 instrument rows but eGFR outcome extraction returned 0 -- refusing all-NA concordance", nrow(exp_dat)))
    cat("No eGFR p<5e-8 instruments to test; outcome extraction empty (expected).\n")
    mr_dt <- data.table()
  }
}

# ── Build concordance table for the eGFR hits (16 in the current 76-pair table) ─────────────────
# Load original Phase 2 MR results for comparison
orig <- readRDS(file.path(PHASE2_DIR, "mr_all_results.rds"))
# GENERATION-CHECKED, for the same reason as in 15: this table is compared against estimates built
# from config's current instrument set.
.orig_stamp <- attr(orig, "ckm_phase1_md5")
if (is.null(.orig_stamp) || is.na(.orig_stamp) || !nzchar(.orig_stamp))
  stop("16: mr_all_results.rds carries no instrument-set stamp, so it cannot be shown to describe ",
       "the same Phase-1 generation as this run (md5 ", substr(INSTRUMENTS_MD5, 1, 12),
       "). Regenerate it with the current selection.")
if (!identical(tolower(.orig_stamp), tolower(INSTRUMENTS_MD5)))
  stop(sprintf(paste0("16: mr_all_results.rds was built from instrument set %s but this run uses ",
                      "%s."), substr(.orig_stamp, 1, 12), substr(INSTRUMENTS_MD5, 1, 12)))
setDT(orig)
orig_subset <- orig[method %in% c("Wald ratio", "Inverse variance weighted") &
                    outcome_name == "eGFR",
                    .(gene, outcome_name, b_orig = b, se_orig = se, pval_orig = pval)]

egfr_hits <- merge(egfr_hits, orig_subset, by = c("gene", "outcome_name"), all.x = TRUE)
# Codex backfill-review fix: a strong pair must have an original primary MR to compare against -- fail loud, never
# resolve CONCORDANT off a missing b_orig.
if (anyNA(egfr_hits$b_orig)) stop("16_sens_5e8_egfr: eGFR strong pair(s) lack original primary MR (b_orig NA): ", paste(egfr_hits[is.na(b_orig), paste(gene, outcome_name, sep = "-")], collapse = ", "))

# Merge sensitivity results
if (exists("mr_dt") && nrow(mr_dt) > 0) {
  sens_sub <- mr_dt[method %in% c("Wald ratio", "Inverse variance weighted"),
                    .(gene, outcome_name, b_sens = b, se_sens = se,
                      pval_sens = pval, method_sens = method)]
  egfr_hits <- merge(egfr_hits, sens_sub, by = c("gene", "outcome_name"), all.x = TRUE)
} else {
  egfr_hits[, `:=`(b_sens = NA_real_, se_sens = NA_real_,
                    pval_sens = NA_real_, method_sens = NA_character_)]
}

# Classify concordance
egfr_hits[, concordant := fcase(
  !has_instrument_5e8,                             "NO_INSTRUMENT",
  is.na(pval_sens),                                "NOT_TESTED",
  sign(b_sens) == sign(b_orig) & pval_sens < SUGGESTIVE_P, "CONCORDANT",
  sign(b_sens) == sign(b_orig),                    "SAME_DIRECTION",
  default = "DISCORDANT"
)]

# Resolve NOT_TESTED where instrument is same single SNP
egfr_hits[, resolution_note := NA_character_]
api_fail_idx <- egfr_hits[has_instrument_5e8 == TRUE & is.na(b_sens), which = TRUE]
if (length(api_fail_idx) > 0) {
  cat(sprintf("\nResolving %d API-failure NOT_TESTED:\n", length(api_fail_idx)))
  # For single-SNP instruments that pass both thresholds, sensitivity = original
  for (idx in api_fail_idx) {
    g <- egfr_hits$gene[idx]
    n_orig <- nrow(instr[gene == g])
    n_5e8  <- nrow(instr_5e8[gene == g])
    .fin <- all(is.finite(c(egfr_hits$b_orig[idx], egfr_hits$se_orig[idx], egfr_hits$pval_orig[idx])))
    if (n_orig == 1 && n_5e8 == 1 && .fin) {
      .conc <- if (egfr_hits$pval_orig[idx] < SUGGESTIVE_P) "CONCORDANT" else "SAME_DIRECTION"   # reapply the same rule as above, from config
      egfr_hits[idx, `:=`(
        b_sens = b_orig, se_sens = se_orig, pval_sens = pval_orig,
        method_sens = "Wald ratio",
        concordant = .conc,
        resolution_note = sprintf("Single-SNP p<5e-8; sens=orig; reclassified %s", .conc)
      )]
      cat(sprintf("  %s: resolved to %s (single SNP)\n", g, .conc))
    }
  }
}

cat("\n=== eGFR p<5e-8 Concordance ===\n")
print(egfr_hits[, .(gene, outcome_name, has_instrument_5e8,
                     b_orig = round(b_orig, 3), b_sens = round(b_sens, 3),
                     pval_orig = formatC(pval_orig, format = "e", digits = 2),
                     pval_sens = formatC(pval_sens, format = "e", digits = 2),
                     concordant)])

cat("\nSummary:\n")
print(egfr_hits[, .N, by = concordant])

# ── Merge with the existing non-eGFR concordance (data-driven; currently 38 pairs) ──────────────
existing <- fread(file.path(out_dir, "sens_5e8_concordance.csv"))
# GENERATION-CHECKED BEFORE ANYTHING IS WRITTEN. This file is merged with rows computed from
# config's current instrument set and then overwritten in place; the check that used to be the only
# one ran later, on the RDS, by which point the mixed CSV was already on disk.
# EVERY ROW, not row 1. The stamp is a column, so checking `[1]` authenticated one row and then
# line ~220 rewrote the whole column with this run's hash: a CSV whose first row was current and
# whose remaining rows came from another generation passed and was published as one set.
if (!"ckm_phase1_md5" %in% names(existing))
  stop("16: sens_5e8_concordance.csv carries no instrument-set stamp, so it cannot be shown to ",
       "describe the same Phase-1 generation as this run (md5 ", substr(INSTRUMENTS_MD5, 1, 12),
       "). Re-run 15 under this selection.")
.cs <- tolower(as.character(existing$ckm_phase1_md5))
if (any(is.na(.cs) | !nzchar(.cs)))
  stop(sprintf(paste0("16: %d of %d rows in sens_5e8_concordance.csv carry no instrument-set ",
                      "stamp. Re-run 15 under this selection."),
               sum(is.na(.cs) | !nzchar(.cs)), length(.cs)))
if (!all(.cs == tolower(INSTRUMENTS_MD5)))
  stop(sprintf(paste0("16: %d of %d rows in sens_5e8_concordance.csv were written from instrument ",
                      "set(s) %s but this run uses %s. Merging would publish one table describing ",
                      "two generations."),
               sum(.cs != tolower(INSTRUMENTS_MD5)), length(.cs),
               paste(substr(unique(.cs[.cs != tolower(INSTRUMENTS_MD5)]), 1, 12), collapse = ", "),
               substr(INSTRUMENTS_MD5, 1, 12)))
# ...AND THE COMPANION RDS IS VALIDATED NOW, BEFORE EITHER FILE IS WRITTEN. The RDS check used to
# live in the block at the bottom, after this CSV had already been overwritten: a refusal there left
# the concordance CSV at the new generation and mr_sens_5e8.rds at the old one. The two files are
# one analysis and are validated together.
if (exists("mr_dt") && nrow(mr_dt) > 0) {
  existing_mr <- readRDS(file.path(out_dir, "mr_sens_5e8.rds"))
  # THE SAME GENERATION 15 WROTE. Every flag at the bottom is recomputed over the merged table with
  # THIS run's BONFERRONI_ALPHA, so an existing table from another selection would be re-flagged
  # with a threshold that does not belong to it.
  .prev_stamp <- attr(existing_mr, "ckm_phase1_md5")
  if (is.null(.prev_stamp) || is.na(.prev_stamp) || !nzchar(.prev_stamp))
    stop("16: ", file.path(out_dir, "mr_sens_5e8.rds"), " carries no instrument-set stamp. It ",
         "predates the stamp or was written by an older 15; regenerate it by re-running 15 under ",
         "this selection.")
  if (!identical(tolower(.prev_stamp), tolower(INSTRUMENTS_MD5)))
    stop(sprintf(paste0("16: mr_sens_5e8.rds was written by 15 from instrument set %s, but this ",
                        "run uses %s. Appending would put two generations in one table and flag ",
                        "all of it with this run's threshold."),
                 substr(.prev_stamp, 1, 12), substr(INSTRUMENTS_MD5, 1, 12)))
}

# P2-F4 (cold audit): overwrite the pre-eGFR backup each run + fail loud on copy failure (was: no overwrite -> a
# stale backup from an earlier non-eGFR set could persist; unchecked return could silently lose the backup).
if (!file.copy(file.path(out_dir, "sens_5e8_concordance.csv"),
               file.path(out_dir, "sens_5e8_concordance.csv.bak_pre_egfr"), overwrite = TRUE))
  stop("16: failed to write sens_5e8_concordance.csv.bak_pre_egfr")

# THE ROWS THIS RUN COMPUTED WIN. duplicated() keeps the FIRST occurrence, and `existing` came
# first, so once the CSV already held eGFR rows -- which it does after any previous run of this
# script -- the freshly computed ones were the copies that got discarded, and the published table
# was the older run's.
combined <- rbind(existing[!egfr_hits, on = c("gene", "outcome_name")], egfr_hits, fill = TRUE)
combined <- combined[!duplicated(combined, by = c("gene", "outcome_name")), ]
combined[, ckm_phase1_md5 := unname(INSTRUMENTS_MD5)]
# Codex backfill-review fix: validate the combined concordance covers EXACTLY the corrected strong pairs (a stale
# non-eGFR file would otherwise silently yield a partial/mismatched combined table).
if (exists("P2S_EGFR_FILE")) {
  .pe <- fread(P2S_EGFR_FILE); .exp <- unique(.pe[, paste(gene, outcome, sep = "_")])
  .got <- combined[, paste(gene, outcome_name, sep = "_")]
  if (!setequal(.got, .exp)) stop(sprintf("16_sens_5e8_egfr: combined concordance keys != corrected strong pairs (missing: %s ; extra: %s)", paste(setdiff(.exp, .got), collapse = ","), paste(setdiff(.got, .exp), collapse = ",")))
}
fwrite(combined, file.path(out_dir, "sens_5e8_concordance.csv"))
cat(sprintf("\nCombined concordance: %d rows after dedup (non-eGFR + %d eGFR)\n",
    nrow(combined), nrow(egfr_hits)))

# Also update mr_sens_5e8.rds if we have new MR results
if (exists("mr_dt") && nrow(mr_dt) > 0) {
  # existing_mr was read and generation-checked above, before the concordance CSV was written.
  saveRDS(existing_mr, file.path(out_dir, "mr_sens_5e8.rds.bak_pre_egfr"))
  # REPLACED, NOT APPENDED. The header invites running this script directly, and it appended its
  # eGFR rows to a file that may already contain them: each re-run added another copy, and because
  # the flags below are recomputed by p.adjust() over the whole table, the number of copies changed
  # the BH ranks -- FDR significance depended on how many times 16 had been run. The rows this run
  # computed replace the same (gene, outcome, method) rows rather than joining them.
  data.table::setDT(existing_mr); data.table::setDT(mr_dt)
  .k <- c("gene", "outcome_name", "method")
  if (!all(.k %in% names(existing_mr)) || !all(.k %in% names(mr_dt)))
    stop("16: mr_sens_5e8.rds and the new eGFR results do not share the (",
         paste(.k, collapse = ", "), ") key, so re-running cannot replace rows instead of ",
         "duplicating them.")
  .dropped <- nrow(existing_mr) - nrow(existing_mr[!mr_dt, on = .k])
  if (.dropped)
    cat(sprintf("16: replacing %d row(s) this run recomputed (previous run of this script)\n",
                .dropped))
  combined_mr <- rbind(existing_mr[!mr_dt, on = .k], mr_dt, fill = TRUE)
  # Recompute significance flags for appended rows
  # P2-N1 / P6-N-02 (cold audit): RECOMPUTE flags over the FULL combined table, not only the NA/appended subset.
  # The prior `is.na(fdr_sig)` fill left the existing (non-eGFR) rows' BH-FDR stale -- computed before the eGFR
  # rows were appended -- mis-flagging 11 pairs. FDR is recomputed by method to match the main screen
  # (cascade/01: p.adjust(pval, "BH"), by=method).
  BONF <- BONFERRONI_ALPHA
  combined_mr[, bonf_sig   := pval < BONF]
  combined_mr[, fdr_sig    := p.adjust(pval, method = "BH") < FDR_ALPHA, by = method]
  combined_mr[, suggestive := pval < SUGGESTIVE_P]
  # RE-STAMPED: rbind() drops the attributes, so the merged table would go back to disk unstamped
  # and the check above would refuse the very file this script just wrote.
  attr(combined_mr, "ckm_phase1_md5")  <- unname(INSTRUMENTS_MD5)
  attr(combined_mr, "ckm_phase1_file") <- PHASE1_INSTRUMENTS
  saveRDS(combined_mr, file.path(out_dir, "mr_sens_5e8.rds"))
  cat(sprintf("MR results updated: %d rows\n", nrow(combined_mr)))
}

# ── Overall summary ──────────────────────────────────────────────────────────
cat("\n=== Overall p<5e-8 Concordance ===\n")
print(combined[, .N, by = concordant])

cat("\n=== Done ===\n")
