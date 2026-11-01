# Phase 3-S: Steiger Filtering for the strong coloc hits (dual rule PP.H4>=0.80 & PP.H4/PP.H3>=5)
# =========================================================
# Tests causal direction: protein -> disease (expected)
# vs disease -> protein (reverse causation, unexpected)
#
# Method (2026-07-01 CORRECTED, gate #13 finding 4): Steiger directionality from the SNP->EXPOSURE and
# SNP->OUTCOME associations of the SAME instrument SNPs, summed across all instruments.
#   r2_exp = sum get_r_from_bsen(pQTL beta, se, N)^2
#   r2_out = sum of SNP->disease variance from the LOCAL OUTCOME EXTRACT:
#            get_r_from_lor(beta_out, eaf, ncase, ncontrol, prev) [case-control] /
#            get_r_from_bsen(beta_out, se_out, N) [quant].
# The OLD version derived r2_out from the MR causal-estimate p-value (NOT the SNP->outcome association;
# invalid for multi-SNP IVW) — that mislabelled steiger_dir and could distort tiering.
# steiger_dir = TRUE (r2_exp > r2_out) = protein->disease (expected). Flag FALSE.
#
# Run: Rscript "<E>/20260331-CKM/Paper 3/scripts/phase3s_steiger.R"

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
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

Sys.setenv(CKM_OFFLINE = "1")   # Steiger is fully local (no API); blank tokens so any stray call fails loud
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/00_config.R"))

cat("=== Phase 3-S: Steiger Filtering ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

out_dir <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG))   # CORRECTED workspace (live untouched)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# THE PREVIOUS RESULT GOES FIRST, BEFORE ANY INPUT IS READ. The outputs are written only after the
# fitting loop, so a run that stops anywhere before that -- a missing population prevalence, a
# degenerate pair, an absent or truncated outcome extract, an unreadable coloc/instrument/MR RDS --
# used to leave the PREVIOUS run's steiger_results untouched at the canonical names. Nothing
# downstream would notice: 13 and 99 check the file's structure, not its run id, its prevalence, or
# the code that made it, so a table built with the old sample-fraction transformation would be
# consumed as this run's result. Absence is a hard error for every consumer; a stale table is not.
#
# MOVED ABOVE THE INPUT READS (objective review 2026-08-08). It used to sit just before the loop --
# after the eight outcome extracts and the coloc/instrument/MR artifacts had already been read --
# so every one of those reads was a way to abort with the stale result still in place, which is
# exactly the set of failures the removal exists to cover.
for (.old in file.path(out_dir, c("steiger_results.rds", "steiger_results.csv"))) {
  if (file.exists(.old) && !file.remove(.old))
    stop(sprintf("09 Steiger: could not remove the previous %s -- refusing to run with a stale result in place that a failed run would leave behind", basename(.old)))
}

# Sourced AFTER the invalidation, not before: out_dir needs 00_config.R (RESULTS_DIR, PHASE1_TAG),
# which is the hard bound on how early the removal can go, but the prevalence table is not needed
# to know the path -- so a missing or malformed prevalence source is no longer a way to abort with
# the stale result still canonical.
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_outcome_prevalence.R"))   # OUTCOME_POP_PREVALENCE (L4-only; kept out of the shared config so an L4 edit cannot re-date the L1-L3 provenance chain)

# ── Outcome metadata (N, type, ncase) ─────────────────────────────────────────
# 2026-07-01: type/name are study design (config); N/ncase are DERIVED from the outcome extracts (17's
# stamped samplesize.outcome / ncase.outcome) — SINGLE source of truth, no hardcoded literals that could
# drift from the MR/coloc extraction or go stale. In particular eGFR N is now the derived Wuttke n_total_sum
# (consistent with 06/26), NOT a hardcoded 567460 (which disagreed with 03's 133413 and 06's derived median).
OUTCOME_META <- data.table(
  outcome_id   = c("FinnGen_R10_I9_CHD", "FinnGen_R10_I9_HEARTFAIL", "FinnGen_R10_T2D",
                   "FinnGen_R10_NAFLD",  "FinnGen_R10_I9_AF",   "MEGASTROKE_EUR_GCST006906",
                   "Wuttke2019_CKD_EA",    "Wuttke2019_eGFR_EA"),
  outcome_name = c("CAD", "HF", "T2D", "NAFLD", "AF", "Stroke", "CKD", "eGFR"),
  type         = c("cc","cc","cc","cc","cc","cc","cc","quant")
)
.ext_dir <- LOCAL_EXTRACT_DIR
OUTCOME_META[, `:=`(n_total = NA_integer_, n_case = NA_integer_)]
for (.i in seq_len(nrow(OUTCOME_META))) {
  .ef <- file.path(.ext_dir, paste0(OUTCOME_META$outcome_id[.i], ".rds"))
  if (!file.exists(.ef)) stop(sprintf("09 Steiger: outcome extract missing (%s) — cannot DERIVE N/ncase; run 17 (MR extraction) first", .ef))
  .oe <- as.data.table(readRDS(.ef))
  OUTCOME_META$n_total[.i] <- as.integer(round(median(.oe$samplesize.outcome, na.rm = TRUE)))
  if ("ncase.outcome" %in% names(.oe) && OUTCOME_META$type[.i] == "cc")
    OUTCOME_META$n_case[.i] <- as.integer(round(median(.oe$ncase.outcome, na.rm = TRUE)))
}
OUTCOME_META[, n_control := n_total - n_case]
cat("09 Steiger: OUTCOME_META N/ncase DERIVED from extracts (single source, no hardcoded N):\n")
print(OUTCOME_META[, .(outcome_name, type, n_total, n_case)])

# ── Load strong hits ───────────────────────────────────────────────────────────
coloc  <- readRDS(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds"))   # CORRECTED combined coloc
instr  <- readRDS(PHASE1_INSTRUMENTS)  # the certified instruments for THIS generation (PHASE1_TAG); reading the live tree instead would skip proteins new to the current set
mr_all <- readRDS(file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "mr_all_results.rds"))   # CORRECTED MR
# THE MR TABLE'S OWN RUN, compared with the extracts below. 09 combines the two; without this it can
# read run A's canonical MR table beside run B's extracts whenever Phase 1 and the outcome files are
# unchanged -- which is the normal case for a rerun that stopped before promoting.
.MR_RUN_ID <- attr(mr_all, "ckm_run_id")
setDT(coloc); setDT(instr); setDT(mr_all)

# 2026-07-01: cached loader for the local outcome extracts (SNP -> outcome associations, produced by 17)
.EXTRACT_DIR <- LOCAL_EXTRACT_DIR
.EXT_CACHE   <- new.env()
get_outcome_extract <- function(oid) {
  if (!is.null(.EXT_CACHE[[oid]])) return(.EXT_CACHE[[oid]])
  f <- file.path(.EXTRACT_DIR, paste0(oid, ".rds"))
  if (!file.exists(f)) return(NULL)
  .x <- readRDS(f)
  # WHICH RUN THIS EXTRACT IS FROM. 01_mr_screen.R writes the eight extracts one at a time, so an
  # interrupted or partly rerun 01 leaves this directory holding a mixture: some files from the
  # current instrument set and outcome GWAS, some from an earlier one. The directory name carries
  # only the generation tag, which does not change when an outcome FILE is replaced, and this
  # function used to read whatever was there. Steiger would then compare instruments from one run
  # against outcome associations from another and report a direction. 01 stamps every extract with
  # the instrument set and the outcome file it was built from; this compares both against what is
  # on disk now.
  # AND THAT IT IS THIS PHENOTYPE. Everything below establishes that the file is internally
  # consistent and from the right run -- none of it establishes that the CHD extract copied to
  # FinnGen_R10_T2D.rds is not being read as T2D. The extract carries the outcome id it was built
  # for and the source file it was built from; the filename is not evidence of either.
  if (!("id.outcome" %in% names(.x)) || !identical(unique(as.character(.x$id.outcome)), oid))
    stop(sprintf(paste0("09 Steiger: %s holds id.outcome '%s', not '%s'. The associations in it are ",
                        "for a different phenotype than the file it is named after."),
                 basename(f), paste(unique(as.character(.x$id.outcome)), collapse = ","), oid))
  .rid <- attr(.x, "ckm_run_id")
  if (is.null(.rid))
    stop(sprintf(paste0("09 Steiger: %s carries no ckm_run_id, so it cannot be shown to come from ",
                        "the same run as the MR table beside it."), basename(f)))
  if (!is.null(.MR_RUN_ID) && !identical(.rid, .MR_RUN_ID))
    stop(sprintf(paste0("09 Steiger: %s comes from run %s, but the MR table beside it comes from ",
                        "run %s. Steiger would combine two runs."),
                 basename(f), .rid, .MR_RUN_ID))
  .p1 <- attr(.x, "ckm_phase1_md5"); .of <- attr(.x, "ckm_outcome_file")
  .om <- attr(.x, "ckm_outcome_md5"); .pr <- attr(.x, "ckm_entry_script_md5")
  # AND THE CODE THAT MADE IT. The three hashes above name the INPUTS, so two runs of different
  # extraction or proxy code over the same Phase 1 and the same outcome files are
  # indistinguishable, and an extract from the earlier one is accepted beside a later run's
  # results. Compared against 01_mr_screen.R as it is on disk now.
  .p01 <- file.path(BASE_DIR, "scripts/paper3.1_ckm/cascade/01_mr_screen.R")
  if (is.null(.pr) || !file.exists(.p01) ||
      !identical(tolower(.pr), tolower(unname(tools::md5sum(.p01)))))
    stop(sprintf(paste0("09 Steiger: %s was produced by a different version of ",
                        "01_mr_screen.R than the one on disk (extract %s, current %s). ",
                        "Rerun cascade/01_mr_screen.R."),
                 basename(f), if (is.null(.pr)) "unstamped" else substr(.pr, 1, 12),
                 if (file.exists(.p01)) substr(unname(tools::md5sum(.p01)), 1, 12) else "absent"))
  if (is.null(.p1) || is.null(.of) || is.null(.om))
    stop(sprintf(paste0("09 Steiger: %s carries no instrument/outcome stamp, so it cannot be shown ",
                        "to belong to the run that produced the MR table. Rerun ",
                        "cascade/01_mr_screen.R."), basename(f)))
  if (!identical(tolower(.p1), tolower(unname(INSTRUMENTS_MD5))))
    stop(sprintf(paste0("09 Steiger: %s was extracted against a DIFFERENT instrument set ",
                        "(extract %s, active %s). Rerun cascade/01_mr_screen.R."),
                 basename(f), substr(.p1, 1, 12), substr(unname(INSTRUMENTS_MD5), 1, 12)))
  .src <- file.path(OUTCOME_DIR, .of)
  if (!file.exists(.src))
    stop(sprintf("09 Steiger: %s names outcome file %s, which is not on disk", basename(f), .of))
  if (!identical(tolower(unname(tools::md5sum(.src))), tolower(.om)))
    stop(sprintf(paste0("09 Steiger: %s was extracted from a different version of %s ",
                        "(extract %s, on disk %s). Rerun cascade/01_mr_screen.R."),
                 basename(f), .of, substr(.om, 1, 12),
                 substr(unname(tools::md5sum(.src)), 1, 12)))
  x <- as.data.table(.x); .EXT_CACHE[[oid]] <- x; x
}

strong <- coloc[is_coloc_pass(PP.H4, PP.H3)]   # dual coloc rule (D1 2026-07-04; was H4-only -> 53)
cat(sprintf("Strong coloc hits: %d pairs, %d proteins\n",
            nrow(strong), uniqueN(strong$gene)))

# ── Helper: r² from p-value and N ─────────────────────────────────────────────
r2_from_pn <- function(p, n) {
  z <- qnorm(p / 2, lower.tail = FALSE)
  (z / sqrt(n))^2
}

# ── Steiger for each pair ──────────────────────────────────────────────────────
results <- list()

for (i in seq_len(nrow(strong))) {
  g   <- strong$gene[i]
  oid <- strong$outcome_id[i]

  # ── Exposure: get instrument for this gene ──────────────────────────────────
  # 2026-07-01 Codex-phase3 #3: these are STRONG coloc pairs -- each MUST have its instrument / MR row /
  # outcome meta / extract / SNP overlap (it passed MR + coloc). A missing one = STALE/INCOMPLETE input
  # (e.g. stale mr_all_results.rds), NOT a legit skip -> fail loud instead of silently dropping the strong pair.
  instr_g <- instr[gene == g]
  if (nrow(instr_g) == 0)
    stop(sprintf("09 Steiger: strong pair %s ~ %s has NO instrument -- stale/incomplete instrument set (a strong coloc hit must carry its instrument)", g, oid))

  # ── MR result for this pair ─────────────────────────────────────────────────
  mr_row <- mr_all[gene == g & outcome_id == oid &
                   method %in% c("Wald ratio", "Inverse variance weighted")]
  if (nrow(mr_row) == 0)
    stop(sprintf("09 Steiger: strong pair %s ~ %s has NO MR result -- stale/incomplete mr_all_results.rds (a strong coloc hit passed suggestive MR)", g, oid))
  mr_row <- mr_row[1]  # take first row if multiple (shouldn't happen)

  # ── Outcome metadata ────────────────────────────────────────────────────────
  out_meta <- OUTCOME_META[outcome_id == oid]
  if (nrow(out_meta) == 0)
    stop(sprintf("09 Steiger: strong pair %s ~ %s has NO outcome meta for %s -- OUTCOME_META/config incomplete", g, oid, oid))
  n_out <- out_meta$n_total[1]

  # ── r²_exp (SNP->protein) and r²_out (SAME SNPs -> disease), summed across instruments ──────────
  F_exp <- max(instr_g$F_stat, na.rm = TRUE)   # kept for reporting (lead F)
  p_mr  <- mr_row$pval[1]                       # kept for reporting only (NOT used for r2_out anymore)
  oe <- get_outcome_extract(oid)
  if (is.null(oe))
    stop(sprintf("09 Steiger: strong pair %s ~ %s -- NO outcome extract for %s (17 must have produced %s). The path is derived from the generation, so a message naming another generation's directory would send you to a folder that does not exist.", g, oid, oid, file.path(.EXTRACT_DIR, paste0(oid, ".rds"))))
  mo <- merge(instr_g[, .(SNP, b_exp = beta, se_exp = se, n_exp = samplesize)],
              oe[, .(SNP, b_out = beta.outcome, se_out = se.outcome, eaf_out = eaf.outcome,
                     n_out = samplesize.outcome,
                     # PROXY PROVENANCE, CARRIED (2026-08-07 editor ruling). The extractor relabels a
                     # proxy row with the TARGET instrument's rsID, so merging on SNP alone attributes
                     # the proxy's outcome association to the instrument and drops every trace that it
                     # happened. 53_correlated_ivw_recompute.R and 62_proxy_provenance_fetch.R already
                     # carry these columns; this was the one step in the pipeline that did not.
                     px = proxy.outcome, px_snp = proxy_snp.outcome, px_r2 = proxy_r2.outcome)],
              by = "SNP")
  if (nrow(mo) == 0)
    stop(sprintf("09 Steiger: strong pair %s ~ %s -- NO instrument SNP overlaps the outcome GWAS (a strong coloc hit that passed MR must overlap; stale extract/instrument?)", g, oid))
  .re_vec <- get_r_from_bsen(mo$b_exp, mo$se_exp, mo$n_exp)^2                            # SNP -> exposure (per SNP)
  if (out_meta$type[1] == "cc") {                                                       # SNP -> disease
    # get_r_from_lor()'s `prevalence` is the POPULATION prevalence of the disease. The sampling
    # fraction is already carried by ncase/ncontrol, so passing n_case/n_total here fed the sample
    # case fraction into a slot that means something else, and every binary outcome's r2_out,
    # steiger_z and steiger_p was computed on the wrong liability transformation -- silently, since
    # nothing in the output records which prevalence was used. Case-control GWAS are deliberately
    # case-enriched, so the two numbers are not close.
    #
    # There is no defensible default here: it is a per-disease epidemiological input. Rather than
    # substitute a plausible-looking number, refuse. OUTCOME_POP_PREVALENCE must be supplied with a
    # cited source per outcome before binary Steiger results can be produced.
    # Two sources, kept apart on purpose (R/_outcome_prevalence.R): a registry prevalence for the
    # FinnGen endpoints, and -- for the two studies where no population prevalence exists -- the
    # sample's own case fraction, which makes the de-ascertainment an exact identity. DERIVED here,
    # never a literal: the identity holds only when K equals ncase/(ncase+ncontrol) exactly.
    if (oid %in% OUTCOME_NO_DEASCERTAINMENT) {
      .prev <- out_meta$n_case[1] / out_meta$n_total[1]
      # Named for what it IS, not for what it resembles: the case fraction implied by the
      # extract-wide counts this call passes to get_r_from_lor, which is what makes the
      # de-ascertainment an identity. It is NOT the study's published sample fraction --
      # CKDGen reports 41,395 of 480,698 (0.0861) while the analysed extract medians give
      # 41,395 of 410,882 (0.1007) -- and labelling it as the latter would misdescribe it.
      .prev_src <- "extract_implied_case_fraction_no_deascertainment"
    } else {
      .prev <- OUTCOME_POP_PREVALENCE[[oid]]
      .prev_src <- "registry_endpoint_prevalence"
    }
    if (is.null(.prev) || !is.finite(.prev) || .prev <= 0 || .prev >= 1)
      stop(sprintf("09 Steiger: no prevalence configured for binary outcome %s (%s). Either add %s to OUTCOME_POP_PREVALENCE in R/_outcome_prevalence.R with its source, or -- if that study has no definable population prevalence -- add it to OUTCOME_NO_DEASCERTAINMENT in the same file, which supplies the sample case fraction and applies no de-ascertainment. NOT R/00_config.R: that file is sourced by every layer and its md5 is hashed into the L3 provenance sidecar, so an L4-only table there re-dates the shared config and breaks the L1-L3 chain.", oid, strong$outcome_name[i], oid))
    # ncase/ncontrol here stay the extract-wide values ON PURPOSE: .prev is derived as
    # n_case/n_total from the same object, and the no-de-ascertainment identity holds only if the
    # two agree exactly. Changing these would silently break it.
    .ro_vec <- get_r_from_lor(mo$b_out, mo$eaf_out, out_meta$n_case[1], out_meta$n_control[1], .prev)^2
    # The Fisher-z denominator is built from the SNPs the test is on, but the case count is NOT
    # available per SNP and must not be treated as if it were. 01_mr_screen.R writes the study-wide
    # case count to every row; only the analysed-sample size varies by variant (CKD: 168,302-444,971
    # against a 410,882 median), and CKDGen publishes no per-variant case counts. A first version
    # here subtracted the study-wide case count from the per-SNP N, which silently assumed every
    # participant missing at that variant was a control.
    #
    # The assumption made instead, stated: the case fraction is constant across variants at the
    # STUDY-WIDE value, so a variant analysed in N subjects contributes n_eff = 4*phi*(1-phi)*N.
    # That uses the per-variant N, which is real, without inventing a per-variant case split. It is
    # an explicitly approximate missingness model -- variant participation can change cohort
    # composition in a meta-analysis -- and the CKD z is to be read as approximate accordingly.
    #
    # phi MUST come from the study's own totals, not from the extract. Dividing the study-wide case
    # count by the median per-variant N mixes two granularities and gives 41,395/410,882 = 0.1007
    # where CKDGen publishes 41,395 cases and 439,303 controls, i.e. 0.0861. That constructed
    # fraction has one legitimate use, immediately below: .prev must equal ncase/(ncase+ncontrol) as
    # passed to get_r_from_lor or the no-de-ascertainment identity breaks. It is not a study-wide
    # case fraction and must not be reused as one.
    .phi <- OUTCOME_STUDY_CASE_FRACTION[[oid]]
    if (is.null(.phi) || !is.finite(.phi) || .phi <= 0 || .phi >= 1)
      stop(sprintf("09 Steiger: no study-wide case fraction configured for %s (%s). The Fisher-z denominator needs the study's own published case:total ratio, which is NOT the extract-implied one -- add %s to OUTCOME_STUDY_CASE_FRACTION in R/_outcome_prevalence.R with its source.", oid, strong$outcome_name[i], oid))
    .neff_vec <- 4 * .phi * (1 - .phi) * mo$n_out
    n_eff_out <- median(.neff_vec[is.finite(.neff_vec) & .neff_vec > 3], na.rm = TRUE)
    if (!is.finite(n_eff_out) || n_eff_out <= 3)
      stop(sprintf("09 Steiger: %s ~ %s -- no usable per-SNP effective N over %d SNP; the Fisher-z denominator would be undefined", g, oid, nrow(mo)))
    .phi_used <- .phi; .neff_model <- "4*phi*(1-phi)*n_per_snp"
    # WHERE phi CAME FROM, not only its value (L4 plan item 1, finding 3). phi and the
    # population prevalence are different quantities with different sources -- the first is
    # this study's own case:total ratio, the second a registry estimate -- and a reader given
    # only the number cannot tell which denominator the Fisher-z used.
    .phi_src <- "study_case_fraction (OUTCOME_N_CASE / OUTCOME_N_TOTAL, config-derived)"
  } else {
    .prev <- NA_real_; .prev_src <- "not_applicable_continuous_outcome"
    .phi_used <- NA_real_; .neff_model <- "median_n_per_snp"
    .phi_src <- "not_applicable_continuous_outcome"
    .ro_vec <- get_r_from_bsen(mo$b_out, mo$se_out, mo$n_out)^2
    # The SAME SNPs that r2_out was computed on, not the extract-wide total. r2_out already uses the
    # per-SNP mo$n_out; taking the whole-extract N for the Fisher-z denominator made the two halves
    # of the test disagree about the sample. On this generation eGFR's per-SNP N spans 381,345 to
    # 567,450 around a 525,153 median, so 14 of 16 eGFR pairs and both CKD pairs carried a z computed
    # at an N none of their SNPs had.
    n_eff_out <- as.integer(median(mo$n_out, na.rm = TRUE))
    if (!is.finite(n_eff_out) || n_eff_out <= 3)
      stop(sprintf("09 Steiger: %s ~ %s -- no usable per-SNP outcome N (median over %d SNP is %s); the Fisher-z denominator would be undefined", g, oid, nrow(mo), n_eff_out))
  }
  n_exp <- as.integer(median(mo$n_exp, na.rm = TRUE))
  # 2026-07-01 (editor REQUIRE; Codex G1 fix): a DEGENERATE Steiger must NOT be written as an affirmative
  # steiger_dir/flag_reverse (it becomes a published "reverse causation (Steiger borderline)" label, and
  # unlike prior D2 / SuSiE D3 there was no 99 value-guard behind it). TWO degeneracy shapes:
  #  (a) a side with NO usable (non-NA) r term -> sum(na.rm=TRUE) would SILENTLY make its r2 == 0 and, if the
  #      other side is positive, escape the both-zero test. Require >=1 usable term per side.
  #  (b) non-finite r2, or both r2 == 0. (A TRUE one-sided ~0 with usable terms is still allowed.)
  if (sum(!is.na(.re_vec)) == 0 || sum(!is.na(.ro_vec)) == 0)
    fail_loud_data(sprintf("Steiger DEGENERATE for %s ~ %s: NO usable r term on the %s side (all NA over %d SNP) -- cannot assign a trustworthy direction", g, oid, if (sum(!is.na(.re_vec)) == 0) "exposure" else "outcome", nrow(mo)))
  r2e <- sum(.re_vec, na.rm = TRUE)
  r2o <- sum(.ro_vec, na.rm = TRUE)
  # ATTENUATION-CORRECTED OUTCOME r2, reported ALONGSIDE r2_out and never in place of it (editor
  # ruling 2026-08-07). For a proxy in LD r2_LD with the target, r2_out(proxy) ~= r2_LD *
  # r2_out(target), so the proxy UNDERSTATES the outcome-side variance explained -- and Steiger
  # declares forward when r2_exp > r2_out, so the bias makes forward EASIER to declare. It does not
  # err toward caution. Exact on the linear scale, first order under the liability transform, and
  # eaf_out is the proxy's: this is a bounded adjustment, not a reconstruction of the target's
  # statistics. Only the outcome side is touched; r2_exp uses the target's own UKB-PPP numbers.
  # A PROXY WITH NO USABLE r2 IS A BROKEN JOIN, NOT A PERFECT PROXY. The substitution below silently
  # replaced a missing or out-of-range px_r2 with 1 -- the identity correction -- and because
  # `NA < 1` is not selected by which(), the row also fell out of the proxy index and the output
  # claimed proxy_used=FALSE. Correction and provenance disappeared together, leaving a row
  # indistinguishable from a direct-SNP analysis. px_r2 comes from the outcome extract's own
  # metadata, so an unusable value on a row flagged as a proxy is a data defect, not a case to
  # absorb. The ifelse is kept for px==FALSE rows, whose px_r2 is NA by definition.
  .px_bad <- which(mo$px %in% TRUE & (is.na(mo$px_r2) | !is.finite(mo$px_r2) |
                                      mo$px_r2 <= 0 | mo$px_r2 > 1 |
                                      is.na(mo$px_snp) | !nzchar(as.character(mo$px_snp))))
  if (length(.px_bad))
    fail_loud_data(sprintf("Steiger %s ~ %s: %d proxy row(s) carry an unusable px_r2 or no proxy SNP (r2=%s, snp=%s) -- the LD metadata did not join, so neither the attenuation correction nor the proxy provenance can be trusted for %s",
                           g, oid, length(.px_bad),
                           paste(mo$px_r2[.px_bad], collapse = ", "),
                           paste(mo$px_snp[.px_bad], collapse = ", "),
                           paste(mo$SNP[.px_bad], collapse = "; ")))
  .rl <- ifelse(is.na(mo$px_r2) | !is.finite(mo$px_r2) | mo$px_r2 <= 0 | mo$px_r2 > 1, 1, mo$px_r2)
  r2o_corr <- sum(.ro_vec / .rl, na.rm = TRUE)
  # A PERFECT-LD PROXY IS STILL A PROXY (editor ruling 2026-08-08, deciding against the reading
  # that r2 = 1 puts a pair on the "no proxy" side). The ruling says "a proxy-based pair", not "a
  # pair with proxy r2 < 1": the outcome association was measured at ANOTHER variant and transported
  # through an LD reference, and the point correction being the identity removes neither that
  # provenance nor the uncertainty of relying on the surrogate relationship. The two clauses have
  # different roles -- the direction-change clause responds to the numerical correction and cannot
  # fire at r2 = 1; the factor-of-two margin clause protects a direction call resting on a proxy and
  # can. So ONE index drives both the margin clause and the provenance fields.
  #
  # This index previously read `mo$px %in% TRUE & mo$px_r2 < 1`, which excluded ADM~AF
  # (rs2403310 <- rs2015148) and PECAM1~CAD (rs1962147 <- rs8070460), both r2 = 1, from BOTH uses:
  # they were published proxy_used = FALSE with blank proxy_snps, i.e. the table said the instrument
  # itself carried the outcome association. The impact on their steiger_class is measured AFTER this
  # scope change, never before -- measuring first would reverse the pre-result ordering the margin
  # rule exists to preserve.
  .px_i <- which(mo$px %in% TRUE)
  if (!is.finite(r2e) || !is.finite(r2o) || (r2e == 0 && r2o == 0))
    fail_loud_data(sprintf("Steiger DEGENERATE for %s ~ %s (r2_exp=%.3g, r2_out=%.3g, nSNP=%d) -- cannot assign a trustworthy causal direction; refusing to publish a degenerate reverse-causation label", g, oid, r2e, r2o, nrow(mo)))
  steiger_dir <- r2e > r2o                                                              # TRUE = protein->disease
  # Steiger z-test on the summed r (Fisher-z on the two correlations)
  .re <- sqrt(min(r2e, 0.999)); .ro <- sqrt(min(r2o, 0.999))
  z1 <- 0.5 * log((1 + .re) / (1 - .re)); z2 <- 0.5 * log((1 + .ro) / (1 - .ro))
  steiger_z   <- (z1 - z2) / sqrt(1/max(n_exp - 3, 1) + 1/max(n_eff_out - 3, 1))
  steiger_p   <- 2 * pnorm(-abs(steiger_z))
  .roc <- sqrt(min(r2o_corr, 0.999))
  z2c  <- 0.5 * log((1 + .roc) / (1 - .roc))
  steiger_z_corr <- (z1 - z2c) / sqrt(1/max(n_exp - 3, 1) + 1/max(n_eff_out - 3, 1))
  steiger_p_corr <- 2 * pnorm(-abs(steiger_z_corr))
  # WRITTEN BEFORE ANY CORRECTED VALUE HAD BEEN COMPUTED (editor ruling 2026-08-07). A proxy-based
  # pair whose direction changes under the correction, or whose two r2 fall within a factor of two of
  # each other, is direction-indeterminate rather than forward-confirmed. Steiger is one of the five
  # Convergent criteria, so this rule could not be written after seeing which pairs it catches.
  # Pairs with no proxy are left alone: the correction is the identity for them, and applying the
  # margin clause to them would move pairs the ruling does not cover.
  .dir_corr <- r2e > r2o_corr
  .ratio <- if (min(r2e, r2o_corr) > 0) max(r2e, r2o_corr) / min(r2e, r2o_corr) else Inf
  steiger_class <- if (!length(.px_i)) {
                     if (steiger_dir) "forward" else "reverse"
                   } else if (!identical(steiger_dir, .dir_corr) || .ratio < 2) {
                     "direction_indeterminate"
                   } else if (steiger_dir) "forward" else "reverse"

  results[[i]] <- data.table(
    gene          = g,
    outcome_id    = oid,
    outcome_name  = strong$outcome_name[i],
    # WHICH K WAS USED AND WHERE IT CAME FROM. Without these two a Steiger table is indistinguishable
    # from one computed under any other K -- including the sample fractions this script used until
    # 2026-08-05. The script's own comment complained that "nothing in the output records which
    # prevalence was used"; filling the table did not close that, writing these does.
    prev_used     = .prev,
    prev_source   = .prev_src,
    # AND WHICH FISHER-Z DENOMINATOR. n_eff_out is what the z is actually divided by; phi_used and
    # n_eff_model say how it was arrived at. Without them a binary row and a continuous row are
    # indistinguishable in the table even though their denominators come from different constructions,
    # and a changed phi leaves no trace at all.
    phi_used      = .phi_used,
    phi_source    = .phi_src,
    n_eff_out     = as.numeric(n_eff_out),
    n_eff_model   = .neff_model,
    PP.H4         = strong$PP.H4[i],
    method        = mr_row$method[1],
    n_snps        = mr_row$n_snps[1],
    mr_b          = mr_row$b[1],
    mr_p          = p_mr,
    lead_F        = F_exp,
    r2_exp        = round(r2e, 6),
    r2_out        = round(r2o, 6),
    steiger_z     = round(steiger_z, 3),
    steiger_p     = round(steiger_p, 4),
    steiger_dir   = steiger_dir,   # TRUE = forward (protein->disease)
    flag_reverse  = !steiger_dir,  # TRUE = potential reverse causation
    # WHICH SNPs THE OUTCOME SIDE ACTUALLY CAME FROM. r2_out above is computed on the rows the
    # outcome extract supplied; for these SNPs those rows are the proxy's, not the instrument's.
    proxy_used    = length(.px_i) > 0L,
    proxy_snps    = if (length(.px_i))
                      paste(sprintf("%s<-%s(r2=%.4f)", mo$SNP[.px_i], mo$px_snp[.px_i],
                                    mo$px_r2[.px_i]), collapse = "; ") else NA_character_,
    proxy_r2_min  = if (length(.px_i)) min(mo$px_r2[.px_i], na.rm = TRUE) else NA_real_,
    r2_out_corrected    = round(r2o_corr, 6),
    steiger_z_corrected = round(steiger_z_corr, 3),
    steiger_p_corrected = round(steiger_p_corr, 4),
    steiger_class = steiger_class
  )
}

results_dt <- rbindlist(results, fill = TRUE)
results_dt <- results_dt[order(results_dt$PP.H4, decreasing = TRUE)]

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(results_dt, file.path(out_dir, "steiger_results.rds"))
fwrite(results_dt,  file.path(out_dir, "steiger_results.csv"))
cat(sprintf("\nSteiger results saved: %d pairs\n", nrow(results_dt)))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Steiger Filtering Summary ===\n")
cat(sprintf("Total pairs tested:            %d\n", nrow(results_dt)))
cat(sprintf("steiger_dir = TRUE (forward):  %d  (protein -> disease, expected)\n",
    sum(results_dt$steiger_dir, na.rm = TRUE)))
cat(sprintf("steiger_dir = FALSE (reverse): %d  (flag for reverse causation)\n",
    sum(!results_dt$steiger_dir, na.rm = TRUE)))
cat(sprintf("Steiger p < 0.05:              %d\n",
    sum(results_dt$steiger_p < 0.05, na.rm = TRUE)))

flagged <- results_dt[flag_reverse == TRUE | (steiger_p < 0.05 & !steiger_dir)]
if (nrow(flagged) > 0) {
  cat("\n⚠️  Flagged pairs (possible reverse causation):\n")
  print(flagged[, .(gene, outcome_name, PP.H4, mr_b, mr_p, r2_exp, r2_out,
                     steiger_z, steiger_p, steiger_dir)])
} else {
  cat("\n✓ All pairs pass Steiger direction test (protein -> disease)\n")
}

cat("\nFull results:\n")
print(results_dt[, .(gene, outcome_name, PP.H4, r2_exp, r2_out,
                      steiger_z, steiger_p, steiger_dir)])
