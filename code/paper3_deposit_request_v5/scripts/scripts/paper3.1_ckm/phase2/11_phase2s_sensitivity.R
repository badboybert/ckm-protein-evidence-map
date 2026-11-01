# Phase 2-S: Multi-SNP Sensitivity Analyses
# =========================================================
# For the 40 proteins with >=2 cis-pQTL instruments:
#   (a) Cochran's Q heterogeneity (IVW)
#   (b) Leave-one-out (LOO)
#   (c) MR-PRESSO (outlier detection)
#
# Requires IEU OpenGWAS API (JWT token in 00_config.R).
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/phase2s_sensitivity.R"

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
  library(MRPRESSO)
})

source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

cat("=== Phase 2-S: Multi-SNP Sensitivity Analyses ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

out_dir <- file.path(RESULTS_DIR, "phase2_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
instr   <- readRDS(file.path(PHASE1_DIR, "ukb_cis_instruments.rds"))
mr_all  <- readRDS(file.path(PHASE2_DIR, "mr_all_results.rds"))
setDT(instr); setDT(mr_all)

# ── Identify 40 multi-SNP proteins ───────────────────────────────────────────
snp_counts <- instr[, .N, by = gene]
multi_genes <- snp_counts[N >= 2, gene]
cat(sprintf("Multi-SNP proteins (>=2 SNPs): %d\n", length(multi_genes)))  # expected: 40

multi_instr <- instr[gene %in% multi_genes]
cat(sprintf("  SNP distribution: %s\n",
    paste(sort(table(snp_counts[N >= 2, N])), collapse=" SNPs=")))

# ── Format exposure data ───────────────────────────────────────────────────────
exp_dat <- format_data(
  dat               = as.data.frame(multi_instr),
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

# ── Fetch outcome data ─────────────────────────────────────────────────────────
cat("Fetching outcome data from IEU OpenGWAS (one outcome at a time)...\n")
outcome_ids <- unname(OUTCOMES[names(OUTCOMES) != "eGFR"])  # skip eGFR (no coordinates)
cat("Outcomes:", paste(outcome_ids, collapse=", "), "\n\n")

# Fetch one GWAS at a time to avoid type-coercion errors when binding chunks
snp_list <- unique(exp_dat$SNP)
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
    # Force numeric types on columns that may be mixed
    if ("n" %in% names(res))           res$n           <- suppressWarnings(as.numeric(res$n))
    if ("ncase"     %in% names(res))   res$ncase       <- suppressWarnings(as.numeric(res$ncase))
    if ("ncontrol"  %in% names(res))   res$ncontrol    <- suppressWarnings(as.numeric(res$ncontrol))
    if ("samplesize.outcome" %in% names(res)) {
      res$samplesize.outcome <- suppressWarnings(as.numeric(res$samplesize.outcome))
    }
    out_list[[oid]] <- res
    cat(sprintf(" %d rows\n", nrow(res)))
  } else {
    cat(" 0 rows (skipped)\n")
  }
}

# Codex backfill-review fix: require EXACT outcome coverage (mirror 15) -- a skipped/empty outcome must not yield a partial heterogeneity_Q
.miss_oid <- setdiff(outcome_ids, names(out_list))
if (length(.miss_oid)) stop("11_phase2s_sensitivity: outcome extraction incomplete -- missing/empty outcome(s): ", paste(.miss_oid, collapse = ", "))
out_dat <- rbindlist(out_list, fill = TRUE)
cat(sprintf("\nTotal outcome data: %d rows across %d GWASs\n\n", nrow(out_dat), length(out_list)))
cat(sprintf("Outcome data fetched: %d rows\n\n", nrow(out_dat)))

# ── Harmonise ─────────────────────────────────────────────────────────────────
# DIRECT AND PROXY ROWS HARMONISED SEPARATELY, mirroring cascade/01_mr_screen.R:758-765. A proxy row
# carries the PROXY variant's frequency under the instrument's allele labels; action = 2 reads that as
# a palindromic frequency conflict and drops the instrument, so this arm lost SNPs the discovery
# estimate kept. Discovery sets eaf.outcome to NA -- no frequency for the instrument's own alleles
# exists in the outcome study -- and sends proxies through action = 1, which still refuses
# incompatible alleles but does not adjudicate by frequency. Direct rows keep action = 2 unchanged.
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
  stop(sprintf("11_phase2s_sensitivity: harmonisation produced no rows from %d outcome row(s) -- refusing to continue and leave stale sidecars in place", nrow(.od)))
harm <- harm[harm$mr_keep, ]
cat(sprintf("Harmonised (mr_keep=TRUE): %d rows (direct %d, proxy %d)\n\n",
            nrow(harm), sum(!.is_prox), sum(.is_prox)))

# ── (a) Cochran's Q heterogeneity ─────────────────────────────────────────────
cat("Running Cochran's Q heterogeneity...\n")
het_res <- tryCatch(
  mr_heterogeneity(harm, method_list = c("mr_ivw", "mr_egger_regression")),
  error = function(e) { message("Q error: ", e$message); NULL }
)
# Codex backfill-review fix: 50 promises heterogeneity_Q.rds -- fail loud rather than silently skip the write
if (is.null(het_res) || nrow(het_res) == 0)
  stop(sprintf("11_phase2s_sensitivity: mr_heterogeneity produced no Q results (harmonised rows=%d) -- refusing to continue without heterogeneity_Q.rds", nrow(harm)))
if (!is.null(het_res)) {
  saveRDS(as.data.table(het_res), file.path(out_dir, "heterogeneity_Q.rds"))
  cat(sprintf("Q results: %d rows\n", nrow(het_res)))
  ivw_het <- het_res[het_res$method == "Inverse variance weighted", ]
  cat(sprintf("  IVW Cochran's Q: %d pairs\n", nrow(ivw_het)))
  sig_het <- ivw_het[!is.na(ivw_het$Q_pval) & ivw_het$Q_pval < 0.05, ]
  cat(sprintf("  Significant Q (p<0.05): %d pairs\n", nrow(sig_het)))
  if (nrow(sig_het) > 0) {
    cat("  Heterogeneous proteins:\n")
    print(sig_het[order(sig_het$Q_pval), c("exposure","outcome","Q","Q_df","Q_pval")])
  }
}

# ── (b) Leave-one-out (LOO) ────────────────────────────────────────────────────
cat("\nRunning Leave-one-out analysis...\n")
loo_res <- tryCatch(
  mr_leaveoneout(harm, method = mr_ivw),
  error = function(e) { message("LOO error: ", e$message); NULL }
)
# Same fail-loud as Q above (phase2audit F-03). Without it the handler's NULL merely skipped the
# write, the script carried on into MR-PRESSO and exited 0, and a PREVIOUS loo_results.rds stayed in
# place -- with no run id, Phase-1 md5 or outcome stamp in the sidecar, nothing downstream could tell
# that file from this run's output.
if (is.null(loo_res) || nrow(loo_res) == 0)
  stop(sprintf("11_phase2s_sensitivity: mr_leaveoneout produced no results (harmonised rows=%d) -- refusing to continue and leave a stale loo_results.rds in place", nrow(harm)))
if (!is.null(loo_res)) {
  loo_dt <- as.data.table(loo_res)
  saveRDS(loo_dt, file.path(out_dir, "loo_results.rds"))
  cat(sprintf("LOO results: %d rows (%d exposure-outcome pairs)\n",
              nrow(loo_dt), uniqueN(paste(loo_dt$exposure, loo_dt$outcome))))

  # Flag pairs where any single SNP drives the result
  # (when b_loo changes sign or becomes non-significant upon SNP removal)
  loo_full <- loo_dt[loo_dt$SNP == "All", .(exposure, outcome, b_full = b)]
  loo_each <- loo_dt[loo_dt$SNP != "All"]
  loo_check <- merge(loo_each, loo_full, by = c("exposure", "outcome"))
  loo_check[, sign_flip := sign(b) != sign(b_full)]
  flip_pairs <- loo_check[sign_flip == TRUE, .(exposure, outcome, SNP_removed = SNP)]
  if (nrow(flip_pairs) > 0) {
    cat("\n  LOO sign-flip warnings (SNP removal changes direction):\n")
    print(flip_pairs)
  } else {
    cat("  LOO: No sign-flips detected (robust)\n")
  }
}

# ── (c) MR-PRESSO ─────────────────────────────────────────────────────────────
cat("\nRunning MR-PRESSO (outlier detection)...\n")
# Only run for proteins with >=3 SNPs in harmonised data
harm_dt   <- as.data.table(harm)
snp_in_harm <- harm_dt[, .N, by = .(exposure, outcome)]
presso_pairs <- snp_in_harm[N >= 3]
cat(sprintf("  Pairs eligible for MR-PRESSO (>=3 SNPs after harmonisation): %d\n",
            nrow(presso_pairs)))

# 2026-07-05 Task-1 final-gate sensitivity addition: seed the MR-PRESSO NbDistribution Monte-Carlo
# (previously unseeded -> non-reproducible global/outlier p) + build an ALWAYS-persisted ledger below.
set.seed(2026)
presso_results <- list()
presso_ledger  <- list()   # one row per PRESSO-eligible (>=3 SNP) pair, ok OR error -> mrpresso_results.rds
for (i in seq_len(nrow(presso_pairs))) {
  exp_i <- presso_pairs$exposure[i]
  out_i <- presso_pairs$outcome[i]
  n_i   <- presso_pairs$N[i]
  harm_sub <- as.data.frame(harm_dt[exposure == exp_i & outcome == out_i])

  tryCatch({
    res <- mr_presso(
      BetaOutcome   = "beta.outcome",
      BetaExposure  = "beta.exposure",
      SdOutcome     = "se.outcome",
      SdExposure    = "se.exposure",
      OUTLIERtest   = TRUE,
      DISTORTIONtest = TRUE,
      data          = harm_sub,
      NbDistribution = 1000,
      SignifThreshold = 0.05
    )
    # 2026-07-05 Task-1 (review fix): MRPRESSO returns the Global Test Pvalue as a STRING (e.g. "<0.001")
    # when it falls below the NbDistribution resolution; %.3f on a string would throw and mislabel a SUCCESS
    # as error. Normalize to numeric for the ledger/compare, keep the raw bound in note, log with %s.
    global_p_raw  <- res$`MR-PRESSO results`$`Global Test`$Pvalue
    global_p_num  <- suppressWarnings(as.numeric(sub("^<\\s*", "", as.character(global_p_raw))))
    global_p_note <- if (grepl("^<", as.character(global_p_raw))) paste0("global_p bound ", global_p_raw) else NA_character_
    n_out_i       <- sum(res$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, na.rm = TRUE)
    presso_results[[paste(exp_i, out_i)]] <- list(
      exposure = exp_i, outcome = out_i, global_p = global_p_num, n_outliers = n_out_i)
    cat(sprintf("  %s ~ %s: global_p=%s, outliers=%d\n",
                gsub(" .*", "", exp_i), out_i, format(global_p_raw), n_out_i))
    # 2026-07-05 Task-1: ledger row (ok)
    presso_ledger[[paste(exp_i, out_i)]] <- data.table(
      exposure = exp_i, outcome = out_i, n_snps = n_i,
      presso_global_p = global_p_num,
      n_outliers      = as.integer(n_out_i),
      status = "ok", note = global_p_note)
  }, error = function(e) {
    message("  PRESSO error (", gsub(" .*", "", exp_i), " ~ ", out_i, "): ", e$message)
    # 2026-07-05 Task-1: ledger row (error) -- <<- writes the loop-frame list from inside the handler
    presso_ledger[[paste(exp_i, out_i)]] <<- data.table(
      exposure = exp_i, outcome = out_i, n_snps = n_i,
      presso_global_p = NA_real_, n_outliers = NA_integer_,
      status = paste0("error:", e$message), note = NA_character_)
  })
}

# 2026-07-05 Task-1 final-gate sensitivity addition: ALWAYS persist mrpresso_results.rds as a ledger with one
# row per PRESSO-eligible (>=3 SNP) pair (status ok / error:<msg>), so downstream can distinguish "ran, nothing
# flagged / all errored" from "never ran". Previously the file was written ONLY when >=1 PRESSO succeeded.
if (length(presso_ledger) > 0) {
  presso_ledger_dt <- rbindlist(presso_ledger, fill = TRUE)
  presso_ledger_dt[, gene       := trimws(gsub("\\(.*\\)", "", exposure))]        # mirror summary-block derivation
  presso_ledger_dt[, outcome_id := trimws(sub(".*\\|\\|\\s*id:", "", outcome))]
  if (!"note" %in% names(presso_ledger_dt)) presso_ledger_dt[, note := NA_character_]
  setcolorder(presso_ledger_dt, c("gene", "exposure", "outcome", "outcome_id", "n_snps",
                                  "presso_global_p", "n_outliers", "status", "note"))
} else {
  # zero PRESSO-eligible pairs (no harmonised pair reached >=3 SNPs): 1-row sentinel, still an artifact
  presso_ledger_dt <- data.table(
    gene = NA_character_, exposure = NA_character_, outcome = NA_character_, outcome_id = NA_character_,
    n_snps = NA_integer_, presso_global_p = NA_real_, n_outliers = NA_integer_,
    status = "no_eligible_pairs", note = "no harmonised pair had >=3 SNPs -> MR-PRESSO not applicable")
}
saveRDS(presso_ledger_dt, file.path(out_dir, "mrpresso_results.rds"))
cat(sprintf("\nMR-PRESSO ledger: %d row(s) -> mrpresso_results.rds\n", nrow(presso_ledger_dt)))
cat(sprintf("  ok=%d | error=%d | global p<0.05=%d\n",
            sum(presso_ledger_dt$status == "ok", na.rm = TRUE),
            sum(startsWith(presso_ledger_dt$status, "error"), na.rm = TRUE),
            sum(presso_ledger_dt$presso_global_p < 0.05, na.rm = TRUE)))

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n=== Phase 2-S Summary ===\n")

# Join Q to MR results for 40 proteins
mr_multi <- mr_all[gene %in% multi_genes & method == "Inverse variance weighted"]
cat(sprintf("IVW results for 40 proteins: %d rows\n", nrow(mr_multi)))

if (!is.null(het_res)) {
  het_ivw <- as.data.table(het_res)[method == "Inverse variance weighted",
                                     .(exposure, outcome, Q, Q_df, Q_pval)]
  # Match back to genes (exposure is like "IDO1 (Q16695)" format)
  het_ivw[, gene := trimws(gsub("\\(.*\\)", "", exposure))]
  # outcome column is TwoSampleMR label "Name || id:ID"; extract pure ID for join
  het_ivw[, outcome_id := trimws(sub(".*\\|\\|\\s*id:", "", outcome))]
  mr_multi_het <- merge(
    mr_multi,
    het_ivw[, .(gene, outcome_id, Q, Q_df, Q_pval)],
    by = c("gene", "outcome_id"), all.x = TRUE
  )
  cat("\nMulti-SNP IVW results with heterogeneity stats:\n")
  # Run fix: cosmetic summary print (heterogeneity_Q.rds already saved above) -- print only columns that exist, so a
  # merge Q.x/Q.y name collision or an empty multi-SNP set cannot crash the script AFTER the real output is written.
  mshcols <- intersect(c("gene","outcome_name","n_snps","b","se","pval","Q","Q_df","Q_pval"), names(mr_multi_het))
  if (nrow(mr_multi_het) && "b" %in% names(mr_multi_het)) print(mr_multi_het[!is.na(b), ..mshcols])
}

cat("\nOutputs saved to:", out_dir, "\n")
cat("  heterogeneity_Q.rds\n  loo_results.rds\n  mrpresso_results.rds\n")
