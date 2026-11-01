# Phase 2-Sens: MR with p < 5e-8 instruments (sensitivity analysis)
# =========================================================
# Compares Phase 2 main results (p<1e-6) vs stricter threshold (p<5e-8)
# Outputs concordance stats and table for paper
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/phase2_sens_5e8.R"

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

# ── BEFORE ANY WORK: the table this run will be compared against ─────────────────────────────
# GENERATION-CHECKED, AND CHECKED FIRST. This comparison table is read row-for-row against the 5e-8
# estimates, so a generation difference would be reported as threshold sensitivity. The check used
# to sit two thirds of the way down, AFTER mr_sens_5e8.rds had already been written -- so a run that
# was going to be refused still published its results file. Reading it here costs one file read and
# makes the refusal happen before anything exists on disk.
orig <- readRDS(file.path(PHASE2_DIR, "mr_all_results.rds"))
.orig_stamp <- attr(orig, "ckm_phase1_md5")
if (is.null(.orig_stamp) || is.na(.orig_stamp) || !nzchar(.orig_stamp))
  stop("15: mr_all_results.rds carries no instrument-set stamp, so it cannot be shown to describe ",
       "the same Phase-1 generation as this run (md5 ", substr(INSTRUMENTS_MD5, 1, 12),
       "). Regenerate it with the current selection.")
if (!identical(tolower(.orig_stamp), tolower(INSTRUMENTS_MD5)))
  stop(sprintf(paste0("15: mr_all_results.rds was built from instrument set %s but this run uses ",
                      "%s. The concordance table would report a generation difference as a ",
                      "threshold effect."),
               substr(.orig_stamp, 1, 12), substr(INSTRUMENTS_MD5, 1, 12)))

cat("=== Phase 2-Sens: p < 5e-8 Sensitivity Analysis ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

out_dir <- file.path(RESULTS_DIR, "phase2_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load instruments and filter ────────────────────────────────────────────
# THE TABLE CONFIG CERTIFIED, not a second read of the same path -- BONFERRONI_ALPHA below is
# derived from config's snapshot, so re-opening the path here is what lets the two disagree.
instr <- read_phase1_instruments()
setDT(instr)

instr_5e8 <- instr[pval < 5e-8]
cat(sprintf("Instruments p<1e-6 (original): %d SNPs, %d proteins\n",
            nrow(instr), uniqueN(instr$gene)))
cat(sprintf("Instruments p<5e-8 (sensitivity): %d SNPs, %d proteins\n",
            nrow(instr_5e8), uniqueN(instr_5e8$gene)))
cat(sprintf("Proteins losing ALL instruments: %d\n\n",
            uniqueN(instr[!gene %in% instr_5e8$gene, gene])))

# ── Format exposure data ───────────────────────────────────────────────────
exp_dat <- format_data(
  dat               = as.data.frame(instr_5e8),
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
cat(sprintf("Exposure data formatted: %d rows\n\n", nrow(exp_dat)))

# ── Fetch outcome data (one at a time) ─────────────────────────────────────
cat("Fetching outcome data from IEU OpenGWAS...\n")
outcome_ids <- unname(OUTCOMES[names(OUTCOMES) != "eGFR"])
snp_list    <- unique(exp_dat$SNP)

out_list <- list()
for (oid in outcome_ids) {
  cat(sprintf("  Fetching %s ...", oid))
  res <- tryCatch(
    extract_outcome_data(
      snps          = snp_list,
      outcomes      = oid,
      proxies       = TRUE,
      rsq           = 0.8,
      align_alleles = 1,
      palindromes   = 1,
      maf_threshold = 0.3
    ),
    error = function(e) { message(" ERROR: ", e$message); NULL }
  )
  if (!is.null(res) && nrow(res) > 0) {
    if ("n"                    %in% names(res)) res$n                    <- suppressWarnings(as.numeric(res$n))
    if ("ncase"                %in% names(res)) res$ncase                <- suppressWarnings(as.numeric(res$ncase))
    if ("ncontrol"             %in% names(res)) res$ncontrol             <- suppressWarnings(as.numeric(res$ncontrol))
    if ("samplesize.outcome"   %in% names(res)) res$samplesize.outcome   <- suppressWarnings(as.numeric(res$samplesize.outcome))
    out_list[[oid]] <- res
    cat(sprintf(" %d rows\n", nrow(res)))
  } else {
    cat(" 0 rows (skipped)\n")
  }
}

# Codex backfill-review fix: require EXACT non-eGFR outcome coverage -- a single skipped/empty outcome would
# otherwise be silently dropped, yielding a partial concordance table.
.miss_oid <- setdiff(outcome_ids, names(out_list))
if (length(.miss_oid)) stop("15_phase2_sens_5e8: outcome extraction incomplete -- missing/empty outcome(s): ", paste(.miss_oid, collapse = ", "))
out_dat <- rbindlist(out_list, fill = TRUE)
cat(sprintf("\nTotal outcome rows: %d\n\n", nrow(out_dat)))

# ── Harmonise & run MR ─────────────────────────────────────────────────────
# DIRECT AND PROXY ROWS HARMONISED SEPARATELY, mirroring cascade/01_mr_screen.R:758-765 and
# phase2/11. action = 2 drops a proxy row as a palindromic frequency conflict because the row carries
# the PROXY variant's frequency under the instrument's allele labels; discovery sets eaf.outcome to NA
# and sends proxies through action = 1. Without this, the p<5e-8 arm reported Wald ratios built on a
# different instrument set than the estimate it is testing.
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
  stop(sprintf("15_phase2_sens_5e8: harmonisation produced no rows from %d outcome row(s) -- refusing to continue", nrow(.od)))
harm <- harm[harm$mr_keep, ]
cat(sprintf("Harmonised (mr_keep=TRUE): %d rows (direct %d, proxy %d)\n\n",
            nrow(harm), sum(!.is_prox), sum(.is_prox)))

mr_res <- mr(harm, method_list = c("mr_wald_ratio", "mr_ivw"))
mr_dt  <- as.data.table(mr_res)

# Add gene name
mr_dt[, gene := trimws(gsub(" \\(.*\\).*", "", exposure))]

# Add outcome name
outcome_map <- c(
  "FinnGen_R10_I9_CHD"     = "CAD",
  "MEGASTROKE_EUR_GCST006906" = "Stroke",
  "FinnGen_R10_I9_HEARTFAIL" = "HF",
  "FinnGen_R10_T2D"        = "T2D",
  "ebi-a-GCST003372"  = "eGFR",
  "Wuttke2019_CKD_EA"        = "CKD",
  "FinnGen_R10_NAFLD"      = "NAFLD",
  "FinnGen_R10_I9_AF"  = "AF"
)
mr_dt[, outcome_name := outcome_map[id.outcome]]

# Multiple testing – use same Bonferroni threshold as main Phase 2 for direct comparability
n_tests   <- nrow(mr_dt[method %in% c("Wald ratio", "Inverse variance weighted")])
bonf_thr  <- BONFERRONI_ALPHA  # config-derived = 0.05 / (N_PROTEINS * N_CKM_OUTCOMES); NOT 0.05/n_tests
mr_dt[method %in% c("Wald ratio", "Inverse variance weighted"),
      bonf_sig := pval < bonf_thr]
# NO FDR FLAG HERE. This script deliberately omits eGFR (line 53: eGFR has no OpenGWAS outcome and
# is fetched from the local extract by 16), so the rows present at this point are SEVEN outcomes of
# an eight-outcome family. A BH correction over an incomplete family is not a provisional version of
# the right answer, it is a different answer: measured on the previous generation's real sensitivity
# table, dropping its 14 eGFR rows moves 11 non-eGFR verdicts (118 significant here, 129 after 16
# appends eGFR and recomputes). Whether a pair was FDR-significant would then depend on whether the
# second script had been run yet. 16 computes it once, over the complete family, by method at
# FDR_ALPHA -- the family the headline screen uses (cascade/01_mr_screen.R:407).
mr_dt[, fdr_sig := NA]
mr_dt[method %in% c("Wald ratio", "Inverse variance weighted"),
      suggestive := pval < SUGGESTIVE_P]

# STAMPED WITH THE GENERATION IT WAS BUILT FROM. 16 reloads this file, appends its own rows and
# recomputes every flag with ITS config; unstamped, a 15 run under one selection and a 16 run under
# another merged into one table with no way to tell, and the reader sees a single consistent-looking
# sensitivity analysis.
attr(mr_dt, "ckm_phase1_md5")  <- unname(INSTRUMENTS_MD5)
attr(mr_dt, "ckm_phase1_file") <- PHASE1_INSTRUMENTS
saveRDS(mr_dt, file.path(out_dir, "mr_sens_5e8.rds"))
cat(sprintf("Results saved: %d rows\n", nrow(mr_dt)))

# ── Summary ────────────────────────────────────────────────────────────────
main_res <- mr_dt[method %in% c("Wald ratio", "Inverse variance weighted")]
cat(sprintf("\nBonferroni threshold (p < %.2e, main Phase 2 BONFERRONI_ALPHA): \n", bonf_thr))
cat(sprintf("  Bonferroni significant: %d pairs\n",  sum(main_res$bonf_sig, na.rm = TRUE)))
cat("  FDR significant: not computed here -- the family is complete only once 16_sens_5e8_egfr.R\n")
cat("                   appends eGFR; it computes fdr_sig over all", length(OUTCOMES), "outcomes.\n")
cat(sprintf("  Suggestive (p<%g): %d pairs\n", SUGGESTIVE_P, sum(main_res$suggestive, na.rm = TRUE)))

# ── Concordance with strong coloc hits ─────────────────────────────────────
cat("\n=== Concordance with the strong coloc hits ===\n")

# Original strong hits (live default below; A6 corrected wrapper sets P2S_STRONG_FILE =
# corrected table_s_pair_evidence.csv so the concordance reference tracks the 52-pair cascade)
if (exists("P2S_STRONG_FILE")) {
  .sh <- data.table::fread(P2S_STRONG_FILE)
  strong_hits <- unique(.sh[outcome != "eGFR", .(gene, outcome_name = outcome)])  # 15 covers non-eGFR (eGFR via 16)
  cat(sprintf("Strong-hit reference from %s: %d non-eGFR pairs\n", basename(P2S_STRONG_FILE), nrow(strong_hits)))
} else {
  # P2-F1 (cold audit 2026-07-08): the hardcoded 33-pair non-eGFR fallback is STALE vs the corrected 38-pair set
  # (missing ANXA4-AF/DPEP1-AF/DTX3-AF/PCSK9-AF/TGFB2-AF/PDE5A-CAD/FGR-HF/IL6-T2D/MINDY1-T2D/TRIM40-T2D; 5 stale
  # extras). A DIRECT run without P2S_STRONG_FILE (which its own header permits) would silently overwrite
  # sens_5e8_concordance.csv with the wrong strong-pair set -> fail loud instead of shipping a stale sensitivity table.
  stop("15_phase2_sens_5e8: P2S_STRONG_FILE not set -- refusing the stale hardcoded 33-pair non-eGFR fallback (cold-audit P2-F1). The corrected wrapper (localize/50) must set P2S_STRONG_FILE = the current table_s_pair_evidence.csv so the concordance reference tracks the live 52-pair set.")
}

# Check which can be tested (gene has instrument in p<5e-8 set)
genes_with_5e8 <- unique(instr_5e8$gene)
strong_hits[, has_instrument_5e8 := gene %in% genes_with_5e8]

# Merge with sensitivity MR results
sens_subset <- main_res[, .(gene, outcome_name, b_sens = b, se_sens = se,
                             pval_sens = pval, method_sens = method)]
strong_hits <- merge(strong_hits, sens_subset, by = c("gene", "outcome_name"), all.x = TRUE)

# Load original MR results for comparison (read and generation-checked at the top of the run)
setDT(orig)
orig_subset <- orig[method %in% c("Wald ratio", "Inverse variance weighted"),
                    .(gene, outcome_name, b_orig = b, se_orig = se, pval_orig = pval, method_orig = method)]
strong_hits <- merge(strong_hits, orig_subset, by = c("gene", "outcome_name"), all.x = TRUE)
# Codex backfill-review fix: a strong pair must have an original primary MR to compare against -- fail loud (mirror 16)
if (anyNA(strong_hits$b_orig)) stop("15_phase2_sens_5e8: strong pair(s) lack original primary MR (b_orig NA): ", paste(strong_hits[is.na(b_orig), paste(gene, outcome_name, sep = "-")], collapse = ", "))

# Classify concordance
strong_hits[, concordant := fcase(
  !has_instrument_5e8,                          "NO_INSTRUMENT",
  is.na(pval_sens),                             "NOT_TESTED",
  sign(b_sens) == sign(b_orig) & pval_sens < SUGGESTIVE_P, "CONCORDANT",
  sign(b_sens) == sign(b_orig),                 "SAME_DIRECTION",
  default = "DISCORDANT"
)]

cat("\nConcordance table:\n")
print(strong_hits[, .(gene, outcome_name, has_instrument_5e8,
                       b_orig = round(b_orig, 3), b_sens = round(b_sens, 3),
                       pval_orig = formatC(pval_orig, format = "e", digits = 2),
                       pval_sens = formatC(pval_sens, format = "e", digits = 2),
                       concordant)])

cat("\nSummary:\n")
print(strong_hits[, .N, by = concordant])

# Resolve NOT_TESTED via the SINGLE-INSTRUMENT IDENTITY check (Codex + editor Ruling #1; mirrors 16_sens_5e8_egfr).
# b_sens=b_orig is exact ONLY when the gene has exactly one instrument at 1e-6 AND one at 5e-8: since instr_5e8 is a
# SUBSET of instr, both being size 1 guarantees the SAME single SNP -> identical Wald ratio. A multi-instrument gene
# that is single-SNP Wald for THIS outcome (e.g. NOTCH2-T2D: 2 instruments) has n_orig>1 and a multi-SNP 5e-8 IVW ->
# NOT an identity -> correctly left NOT_TESTED. (method_orig=="Wald ratio" alone was insufficient -- the closed hole.)
strong_hits[, resolution_note := NA_character_]
api_fail_idx <- strong_hits[has_instrument_5e8 == TRUE & is.na(b_sens), which = TRUE]
if (length(api_fail_idx) > 0) {
  cat(sprintf("\nResolving %d NOT_TESTED via single-instrument identity:\n", length(api_fail_idx)))
  for (idx in api_fail_idx) {
    g <- strong_hits$gene[idx]; on <- strong_hits$outcome_name[idx]
    n_orig <- nrow(instr[gene == g]); n_5e8 <- nrow(instr_5e8[gene == g])
    .fin <- all(is.finite(c(strong_hits$b_orig[idx], strong_hits$se_orig[idx], strong_hits$pval_orig[idx])))
    if (n_orig == 1 && n_5e8 == 1 && .fin) {
      .conc <- if (strong_hits$pval_orig[idx] < SUGGESTIVE_P) "CONCORDANT" else "SAME_DIRECTION"   # reapply the same rule as above, from config
      strong_hits[idx, `:=`(
        b_sens = b_orig, se_sens = se_orig, pval_sens = pval_orig,
        method_sens = "Wald ratio", concordant = .conc,
        resolution_note = sprintf("Single-instrument (n_orig==1 & n_5e8==1); same SNP -> sens=orig; reclassified %s", .conc))]
      cat(sprintf("  %s-%s: resolved %s (single instrument)\n", g, on, .conc))
    } else {
      cat(sprintf("  %s-%s: left NOT_TESTED (n_orig=%d, n_5e8=%d, finite_orig=%s)\n", g, on, n_orig, n_5e8, .fin))
    }
  }
  cat("\nUpdated concordance summary:\n"); print(strong_hits[, .N, by = concordant])
}

# Save concordance table
# STAMPED, like the RDS beside it. 16 reads this CSV, merges its eGFR rows in and overwrites it --
# and it did so before reaching any generation check, so a CSV written under one instrument set was
# published combined with another set's results even when the later check aborted the run. A CSV
# cannot carry attributes, so the stamp is a column.
strong_hits[, ckm_phase1_md5 := unname(INSTRUMENTS_MD5)]
fwrite(strong_hits, file.path(out_dir, "sens_5e8_concordance.csv"))
cat("\nConcordance CSV saved:", file.path(out_dir, "sens_5e8_concordance.csv"), "\n")

cat("\n=== Done ===\n")
