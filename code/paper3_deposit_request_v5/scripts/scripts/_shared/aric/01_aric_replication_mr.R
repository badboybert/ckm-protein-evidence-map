#!/usr/bin/env Rscript
# ============================================================================
# 01_aric_replication_mr.R -- cross-platform replication of the cis-pQTL effect
# using ARIC (Atherosclerosis Risk in Communities, European American, n = 7,213),
# as a parallel arm to the deCODE analysis in `_shared/decode/01`.
#
# WHAT THIS ARM IS. It re-estimates the EXPOSURE effect (the cis-pQTL) for the same
# instrument SNP in an independent cohort on a different assay platform, and reuses
# the same outcome GWAS bit-for-bit. It is therefore a cross-platform check on the
# instrument's effect -- an assay/epitope artefact check plus an exposure estimate
# free of the discovery cohort's selection -- and NOT an independent replication of
# the protein-disease association. Naming it otherwise is what got the deCODE arm
# mislabelled in the manuscript.
#
# WHY ARIC EXISTS ALONGSIDE deCODE, not instead of it. deCODE's download endpoint
# terminates transfers after ~30 s and does not honour Range, leaving 29 of the 76
# pairs unretrievable; ARIC's release is a 352 MB S3 object with Range support and
# covers every gene this study needs. Both arms are reported; neither is primary and
# they are never pooled -- pooling would be a meta-analysis across cohorts, and the
# estimates carry different precision (n 7,213 vs 35,287).
#
# ⚠️ THE F >= 10 GATE IS A WITHIN-COHORT POWER THRESHOLD, NOT AN INSTRUMENT-VALIDITY
# THRESHOLD. ARIC is ~4.9x smaller than deCODE, so the same true effect clears it
# less often. That attenuation is expected and is declared in advance; the threshold
# is NOT relaxed to compensate, because relaxing a criterion after seeing that a new
# source fails it more often is exactly the move this project is trying not to make.
# A pair that is concordant in direction but misses F >= 10 here is recorded as
# NOT EVALUABLE IN ARIC, never as a replication failure.
#
# ⚠️ MULTI-APTAMER RULE IS AN INTERSECTION-UNION TEST. Several genes carry more than
# one SomaScan aptamer. Requiring every aptamer to agree in direction, then taking the
# LEAST favourable value on each statistic (max delta_p, min F, min |beta|), is a
# level-alpha test under any correlation structure and needs no multiplicity
# correction. Taking min(p) across aptamers -- which `decode/01:568` currently does --
# is an uncorrected multiple test; this script does not copy that, and the same fix
# is owed to the deCODE arm.
#
# Inputs are local and read-only. Writes only to ARIC_OUT_DIR.
# Run via the wrapper: scripts/paper3.1_ckm/cascade/12b_aric.R
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

if (!exists("ARIC_LOCAL")) ARIC_LOCAL <- TRUE
if (isTRUE(ARIC_LOCAL)) source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/scripts/_shared/localize/_local_outcome_extract.R"))

BASE <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3")
if (!exists("ARIC_DIR"))     ARIC_DIR     <- file.path(BASE, "data/aric")
if (!exists("ARIC_ZIP"))     ARIC_ZIP     <- file.path(ARIC_DIR, "EA.zip")
if (!exists("ARIC_SEQID"))   ARIC_SEQID   <- file.path(ARIC_DIR, "seqid.txt")
if (!exists("ARIC_N"))       ARIC_N       <- 7213L        # asserted against OBS_CT below, never assumed
if (!exists("ARIC_INSTR"))   ARIC_INSTR   <- file.path(PHASE1_DIR, "ukb_cis_instruments.rds")
if (!exists("ARIC_COLOC"))   ARIC_COLOC   <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "coloc_results_v7.rds")
if (!exists("ARIC_OUT_DIR")) ARIC_OUT_DIR <- file.path(RESULTS_DIR, paste0("repro_", PHASE1_TAG), "aric_replication")
dir.create(ARIC_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== ARIC replication MR ===\nDate:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("zip :", ARIC_ZIP, "\nout :", ARIC_OUT_DIR, "\n\n")
for (.f in c(ARIC_ZIP, ARIC_SEQID, ARIC_INSTR, ARIC_COLOC))
  if (!file.exists(.f)) stop(sprintf("01_aric: missing required input %s", .f))

fail_loud <- function(msg) stop(sprintf("01_aric FATAL: %s", msg), call. = FALSE)

# A partial run must never be mistaken for a complete one: clear stale outputs up
# front so a crash leaves nothing behind rather than last run's file.
for (.s in file.path(ARIC_OUT_DIR, c("aric_replication_summary.csv", "aric_mr_results.rds",
                                     "aric_seqid_map.csv")))
  if (file.exists(.s)) {
    # file.remove() returns FALSE on a locked Windows file and this ignored it, so an interrupted
    # rerun could leave the PREVIOUS generation's artifact in place, complete-looking and read as
    # this run's (objective review 2026-08-06, finding 11). An artifact that cannot be cleared is
    # exactly the state the clearing exists to prevent, so it stops the run.
    if (!file.remove(.s) || file.exists(.s))
      fail_loud(sprintf("could not remove the previous run's %s -- it would survive a failure in this run and be read as this run's result", basename(.s)))
  }

# ── 1. pairs under test: the strong set of the selected generation ───────────
coloc_all <- as.data.table(readRDS(ARIC_COLOC))
strong <- coloc_all[is_coloc_pass(PP.H4, PP.H3)]
if (!nrow(strong)) fail_loud("no strong pairs in coloc_results_v7.rds")
strong_genes <- sort(unique(strong$gene))
cat(sprintf("strong pairs: %d | genes: %d\n", nrow(strong), length(strong_genes)))

instruments <- as.data.table(readRDS(ARIC_INSTR))
if (!all(c("SNP", "gene", "effect_allele", "other_allele", "beta", "se") %in% names(instruments)))
  fail_loud("instrument table lacks SNP/gene/effect_allele/other_allele/beta/se")

# ── 2. gene <-> SeqId <-> UniProt map, PERSISTED ─────────────────────────────
# ARIC's per-aptamer files are named by SeqId alone, so this map is the only thing
# tying a file to a gene. That makes it a single point of failure: it is written out
# as a product and fingerprinted rather than re-derived silently on each run.
sq <- fread(ARIC_SEQID)
if (!all(c("seqid_in_sample", "uniprot_id", "entrezgenesymbol") %in% names(sq)))
  fail_loud("seqid.txt lacks seqid_in_sample/uniprot_id/entrezgenesymbol")

# ⚠️ THE KEY IS UniProt, NOT THE SYMBOL. ARIC records NECTIN2 under its legacy symbol
# PVRL2 (SeqId_6245_4, Q92692). Keying on `entrezgenesymbol` silently drops it -- the
# same failure that made the deCODE arm report "deCODE never assayed this protein" for
# a protein deCODE had assayed. Being produced by a script does not make a number
# UniProt-keyed; every one of this project's coverage errors has been a key error, not
# a data error. Symbol matching is kept only as a cross-check, and any gene the two
# routes disagree on is written out rather than silently resolved.
prot <- fread(file.path(PHASE1_DIR, "ukb_cis_instruments.csv"))
if (!all(c("gene", "uniprot") %in% names(prot))) fail_loud("instrument csv lacks gene/uniprot")
gene_acc <- unique(prot[gene %chin% strong_genes & nzchar(uniprot), .(gene, uniprot)])
.miss_acc <- setdiff(strong_genes, gene_acc$gene)
if (length(.miss_acc))
  cat(sprintf("[warn] %d strong gene(s) carry no uniprot in the instrument table, so only the symbol route can reach them: %s\n",
              length(.miss_acc), paste(.miss_acc, collapse = ", ")))

gmap <- unique(merge(gene_acc, sq[nzchar(uniprot_id), .(seqid = seqid_in_sample, uniprot = uniprot_id,
                                                        aric_symbol = trimws(entrezgenesymbol))],
                     by = "uniprot")[, .(gene, seqid, uniprot, aric_symbol)])
gmap[, resolved_by := fifelse(aric_symbol == gene, "uniprot+symbol", "uniprot_only")]

.by_symbol <- sort(intersect(strong_genes, unique(trimws(sq$entrezgenesymbol))))
.key_disagree <- sort(setdiff(unique(gmap$gene), .by_symbol))
cat(sprintf("key check -- UniProt route: %d gene(s); bare-symbol route: %d; found ONLY by UniProt: %s\n",
            uniqueN(gmap$gene), length(.by_symbol),
            if (length(.key_disagree)) paste(.key_disagree, collapse = ", ") else "none"))
if (length(setdiff(.by_symbol, unique(gmap$gene))))
  cat(sprintf("[warn] found ONLY by symbol (UniProt route missed them): %s\n",
              paste(setdiff(.by_symbol, unique(gmap$gene)), collapse = ", ")))
fwrite(gmap[order(gene, seqid)], file.path(ARIC_OUT_DIR, "aric_seqid_map.csv"))

ARIC_NOT_IN_PANEL <- sort(setdiff(strong_genes, gmap$gene))
cat(sprintf("genes with >=1 ARIC aptamer: %d | not on ARIC's panel: %d\n",
            uniqueN(gmap$gene), length(ARIC_NOT_IN_PANEL)))
if (length(ARIC_NOT_IN_PANEL))
  cat(sprintf("  not on panel: %s\n", paste(ARIC_NOT_IN_PANEL, collapse = ", ")))

# ── 3. which aptamer files the archive actually holds ────────────────────────
zl <- unzip(ARIC_ZIP, list = TRUE)
setDT(zl)
zl[, base := basename(Name)]
zl[, seqid := sub("[.]PHENO1[.]glm[.]linear$", "", base)]
have <- zl[grepl("[.]PHENO1[.]glm[.]linear$", base), .(seqid, Name, Length)]
gmap[, present := seqid %chin% have$seqid]
ARIC_NOT_RETRIEVED <- sort(setdiff(unique(gmap[present == FALSE, gene]), unique(gmap[present == TRUE, gene])))
if (length(ARIC_NOT_RETRIEVED))
  cat(sprintf("on panel but no file in the archive: %s\n", paste(ARIC_NOT_RETRIEVED, collapse = ", ")))

# ── 4. read one aptamer's cis file ───────────────────────────────────────────
read_aptamer <- function(sid) {
  # NB: the parameter must NOT be called `seqid` -- that is a column of `have`, and the
  # column would shadow it inside the data.table i-expression.
  nm <- have[seqid == sid, Name]
  if (!length(nm)) return(NULL)
  con <- unz(ARIC_ZIP, nm[1])
  d <- tryCatch(fread(text = readLines(con), showProgress = FALSE),
                error = function(e) fail_loud(sprintf("could not read %s from the archive: %s", nm[1], conditionMessage(e))))
  req <- c("ID", "REF", "ALT", "A1", "BETA", "SE", "P", "OBS_CT")
  if (!all(req %in% names(d)))
    fail_loud(sprintf("%s lacks required PLINK2 columns (%s)", nm[1], paste(setdiff(req, names(d)), collapse = ", ")))
  d
}

# ── 5. per pair ──────────────────────────────────────────────────────────────
# The outcome side is deliberately identical to the discovery MR: same GWAS, same
# extraction. Only the exposure beta changes. Stating that here because the whole
# interpretation of this arm depends on it.
out_rows <- list(); mr_rows <- list(); k <- 0L

for (i in seq_len(nrow(strong))) {
  gene_i    <- strong$gene[i]
  outcome_i <- strong$outcome_id[i]
  outname_i <- strong$outcome_name[i]
  k <- k + 1L
  base_row <- data.table(gene = gene_i, outcome_id = outcome_i, outcome_name = outname_i,
                         aric_b = NA_real_, aric_se = NA_real_, aric_pval_iut = NA_real_,
                         aric_b_min = NA_real_, aric_b_max = NA_real_,
                         delta_se = NA_real_, delta_pval_iut = NA_real_, aric_F_min = NA_real_,
                         assay_Q = NA_real_, assay_Q_df = NA_real_, assay_Q_pval = NA_real_,
                         n_assays_evaluable = 0L, evaluable_class = "none",
                         bx_aric = NA_real_, bx_ukb = NA_real_, direction_concordant = NA,
                         weak_instrument = NA, high_confidence = FALSE,
                         n_assays = 0L, assay_concordance = NA_character_,
                         n_aptamers_mapped = NA_integer_, n_aptamers_absent_from_archive = NA_integer_,
                         n_aptamers_skipped = NA_integer_, skip_reasons = NA_character_,
                         status = NA_character_)

  if (gene_i %chin% ARIC_NOT_IN_PANEL)   { base_row[, status := "NOT_ON_ARIC_PANEL"]; out_rows[[k]] <- base_row; next }
  if (gene_i %chin% ARIC_NOT_RETRIEVED)  { base_row[, status := "NOT_RETRIEVED"];     out_rows[[k]] <- base_row; next }

  instr_i <- instruments[gene == gene_i]
  if (!nrow(instr_i)) { base_row[, status := "NO_INSTRUMENT"]; out_rows[[k]] <- base_row; next }
  target_rsids <- unique(instr_i$SNP)

  # proxies = TRUE MATCHES THE DISCOVERY ARM. `local_extract_outcome_data` defaults proxies = FALSE
  # so each caller opts in; the main cascade is proxy-included and this arm was not, so for the five
  # of 76 pairs whose discovery estimate used an LD proxy the two arms were NOT reading the same
  # outcome data -- contradicting this script's own claim that only the exposure beta is swapped.
  # It also mislabelled our configuration as a property of the third-party data: ABO~T2D and
  # ABO~Stroke were recorded OUTCOME_MISSING while discovery had p = 2.6e-18 and 5.6e-4 on the same
  # cells. Same class of error as NOT_IN_RELEASE vs NOT_RETRIEVED.
  # palindromic_fallthrough = TRUE IS THE OTHER HALF. proxies = TRUE alone changes nothing for a
  # palindromic instrument, because a palindromic direct row still counts as a direct match, so the
  # SNP never enters `missing` and no proxy is sought -- and harmonise_data(action = 2) then drops
  # it here. Measured: PECAM1-CAD (rs1962147, G/C) and RAD23B-T2D (rs7035725, G/C) reached
  # NO_OVERLAP_AFTER_HARMONISE with proxies = TRUE set. The discovery screen excludes those rows
  # before computing `missing`; this opts into the same rule.
  out_dat <- tryCatch(local_extract_outcome_data(snps = target_rsids, outcomes = outcome_i,
                                                 proxies = TRUE, palindromic_fallthrough = TRUE),
                      error = function(e) fail_loud(sprintf("outcome extraction THREW for %s ~ %s: %s",
                                                            gene_i, outcome_i, conditionMessage(e))))
  if (is.null(out_dat) || !nrow(out_dat)) { base_row[, status := "OUTCOME_MISSING"]; out_rows[[k]] <- base_row; next }

  seqs <- gmap[gene == gene_i & present == TRUE, unique(seqid)]
  # A MAPPED APTAMER THAT IS NOT IN THE ARCHIVE NEVER REACHES THE LOOP, so the skip accounting below
  # cannot see it (round 3, finding 5). A gene is called "not retrieved" only when NONE of its
  # aptamers is present, so a partially-populated archive yields n_aptamers_skipped = 0 and an
  # output shaped exactly like a correct single-aptamer evaluation. Example on this data: silently
  # losing ADAM15's F = 1.06 aptamer would leave only its F = 18.84 one and reverse the
  # weak-instrument verdict. Counted here, before the loop, from the map rather than the archive.
  .absent <- gmap[gene == gene_i & present == FALSE, unique(seqid)]
  per_assay <- list()
  # WHY AN APTAMER PRODUCED NOTHING IS PART OF THE RESULT (objective review 2026-08-06, findings 5
  # and 8). Every `next` below used to be silent, so "this protein has no rsID in ARIC", "the target
  # was present but harmonisation discarded it", and "the fit returned nothing" all arrived at the
  # same label, NO_OVERLAP_AFTER_HARMONISE -- and the third is a code/schema failure wearing the
  # costume of the first two. Measured now: 5 pairs are no-rsID (ABO~Stroke, ABO~T2D, ADM~AF,
  # NCAN~CAD, TINAGL1~CAD) and 2 are harmonisation drops (PECAM1~CAD, RAD23B~T2D).
  skip_why <- if (length(.absent)) rep("aptamer_absent_from_archive", length(.absent)) else character(0)

  for (sq_id in seqs) {
    ap <- read_aptamer(sq_id)
    if (is.null(ap)) { skip_why <- c(skip_why, "no_aptamer_file"); next }
    if (any(ap$OBS_CT != ARIC_N, na.rm = TRUE) && abs(median(ap$OBS_CT, na.rm = TRUE) - ARIC_N) > 1)
      cat(sprintf("  [note] %s median OBS_CT = %s (declared ARIC_N = %d)\n",
                  sq_id, format(median(ap$OBS_CT, na.rm = TRUE)), ARIC_N))
    hit <- ap[ID %chin% target_rsids]
    if (!nrow(hit)) { skip_why <- c(skip_why, "no_target_rsid"); next }

    # PLINK2: A1 is the tested (effect) allele; the other allele is whichever of
    # REF/ALT is not A1. Deriving it rather than assuming REF avoids a silent sign
    # flip on variants where A1 == REF.
    hit[, other := fifelse(A1 == ALT, REF, ALT)]
    ex <- data.frame(SNP = hit$ID, beta.exposure = hit$BETA, se.exposure = hit$SE,
                     effect_allele.exposure = hit$A1, other_allele.exposure = hit$other,
                     eaf.exposure = if ("A1_FREQ" %in% names(hit)) hit$A1_FREQ else NA_real_,
                     pval.exposure = hit$P, exposure = gene_i, id.exposure = sq_id,
                     stringsAsFactors = FALSE)
    # DIRECT AND PROXY ROWS ARE HARMONISED SEPARATELY, and the proxy rows use action = 1.
    #
    # THE REASON IS NOT "assume palindromes are unambiguous". It is that a proxy row's strand was
    # ALREADY resolved inside the orienter, by the same allele-frequency rule harmonise_data(2)
    # would apply, under a coupling guarantee that r2 >= 0.8 supplies. Running the frequency test a
    # second time here re-adjudicates a question that has been answered, and for a palindrome it
    # answers it "ambiguous" and discards the row -- which is exactly how a proxy fetched to rescue
    # an unresolvable direct palindrome gets thrown away again one step later.
    # Direct rows are unchanged: action = 2, frequency-based palindrome rejection, as before.
    .oh <- function(d, act) {
      if (is.null(d) || !nrow(d)) return(NULL)
      tryCatch(harmonise_data(ex, d, action = act),
               error = function(e) fail_loud(sprintf("harmonise THREW for %s / %s ~ %s (action=%d): %s",
                                                     gene_i, sq_id, outcome_i, act, conditionMessage(e))))
    }
    # `%in% TRUE` not `== TRUE`: an NA in proxy.outcome would make `==` return NA, and NA in a row
    # subscript yields a row of NAs rather than dropping it -- the row would then be harmonised
    # twice, once in each arm.
    .is_prox <- if ("proxy.outcome" %in% names(out_dat)) out_dat$proxy.outcome %in% TRUE
                else rep(FALSE, nrow(out_dat))
    harm <- rbind(.oh(out_dat[!.is_prox, , drop = FALSE], 2L),
                  .oh(out_dat[ .is_prox, , drop = FALSE], 1L))
    if (is.null(harm) || !nrow(harm) || !sum(harm$mr_keep)) { skip_why <- c(skip_why, "harmonise_dropped"); next }
    hk <- as.data.table(harm)[mr_keep == TRUE]
    # Derived from the rows that SURVIVED harmonisation, not from what was fetched: a proxy that was
    # fetched and then dropped did not contribute to this estimate and must not be labelled as if it
    # had. Minimum r2 across the contributing proxies is the least favourable one.
    .pk <- if ("proxy.outcome" %in% names(hk)) hk$proxy.outcome %in% TRUE else rep(FALSE, nrow(hk))
    .pu_lab <- if (!any(.pk)) "direct" else if (all(.pk)) "proxy" else "mixed"
    .pu_r2  <- if (any(.pk) && "proxy_r2.outcome" %in% names(hk)) min(hk$proxy_r2.outcome[.pk], na.rm = TRUE) else NA_real_

    pm <- if (nrow(hk) == 1) "mr_wald_ratio" else "mr_ivw"
    r  <- tryCatch(mr(as.data.frame(hk), method_list = pm),
                   error = function(e) fail_loud(sprintf("mr() THREW for %s / %s ~ %s: %s",
                                                         gene_i, sq_id, outcome_i, conditionMessage(e))))
    if (is.null(r) || !nrow(r) || !is.finite(r$b[1]) || !is.finite(r$se[1])) { skip_why <- c(skip_why, "mr_no_estimate"); next }

    Fv <- min((hk$beta.exposure / hk$se.exposure)^2, na.rm = TRUE)
    b_mr <- r$b[1]; se_mr <- r$se[1]
    # second-order (delta) SE: the first-order Wald SE ignores exposure uncertainty,
    # which at ARIC's n is not negligible.
    d_se <- sqrt(se_mr^2 + b_mr^2 / Fv)
    d_p  <- 2 * pnorm(-abs(b_mr / d_se))

    # THE DIRECTION VERDICT COMPARES EXPOSURE EFFECTS, NOT MR RATIOS. Comparing an MR
    # ratio against an instrument beta mixes dimensions; and because the two ratios share
    # the outcome numerator, sign(b_MR,ARIC) == sign(b_MR,UKB) is algebraically identical
    # to sign(beta_X,ARIC) == sign(beta_X,UKB) anyway. Comparing the exposure betas
    # directly makes that explicit and keeps the outcome GWAS out of the verdict entirely,
    # which is what lets this axis be immune to any exposure-outcome sample overlap.
    .u <- instr_i[SNP %chin% hk$SNP]
    .aln <- merge(hk[, .(SNP, ea_a = effect_allele.exposure, oa_a = other_allele.exposure,
                         b_a = beta.exposure)],
                  .u[, .(SNP, ea_u = effect_allele, oa_u = other_allele, b_u = beta)], by = "SNP")
    .aln[, same := toupper(ea_a) == toupper(ea_u) & toupper(oa_a) == toupper(oa_u)]
    .aln[, flip := toupper(ea_a) == toupper(oa_u) & toupper(oa_a) == toupper(ea_u)]
    .aln <- .aln[same | flip]
    .aln[, b_u_aligned := fifelse(flip, -b_u, b_u)]
    bx_dir <- if (nrow(.aln)) all(sign(.aln$b_a) == sign(.aln$b_u_aligned)) else NA

    per_assay[[sq_id]] <- data.table(seqid = sq_id, aric_b = b_mr, aric_se = se_mr,
                                     aric_pval = r$pval[1], delta_se = d_se, delta_pval = d_p,
                                     aric_F_stat = Fv, n_snp = nrow(hk),
                                     bx_aric = .aln$b_a[1], bx_ukb = .aln$b_u_aligned[1],
                                     bx_dir_concordant = bx_dir,
                                     # WHICH ROWS THIS ESTIMATE RESTS ON. A verdict built on an LD
                                     # proxy plus a frequency-based strand rule is not the same
                                     # evidence as one built on the instrument's own outcome row,
                                     # and it must be visible per pair rather than inferred from a
                                     # methods paragraph. "mixed" is possible for a multi-SNP fit.
                                     outcome_proxy_used = .pu_lab, outcome_proxy_r2 = .pu_r2)
    mr_rows[[length(mr_rows) + 1L]] <- cbind(data.table(gene = gene_i, outcome_id = outcome_i,
                                                        outcome_name = outname_i), per_assay[[sq_id]])
  }

  if (!length(per_assay)) {
    # A fit that returned nothing is NOT an overlap problem: it is a failure of this pipeline and
    # must not be reported as a property of ARIC. Ordered most-diagnostic first.
    .st0 <- if ("mr_no_estimate" %chin% skip_why) "ARIC_MR_NO_ESTIMATE"
            else if ("harmonise_dropped" %chin% skip_why) "NO_OVERLAP_AFTER_HARMONISE"
            else if ("no_target_rsid" %chin% skip_why) "NO_TARGET_RSID_IN_ARIC"
            else if ("no_aptamer_file" %chin% skip_why) "APTAMER_FILE_UNREADABLE"
            else "NO_APTAMER_EVALUATED"
    # n_assays COUNTS APTAMERS ATTEMPTED (review finding 7). This branch left the base row's
    # n_assays at 0, so a pair where two aptamers were read and both failed was indistinguishable
    # from a pair where no aptamer exists at all. Current triggers include ADM~AF, which has two
    # aptamers, both read and both checked.
    base_row[, `:=`(status = .st0, n_aptamers_mapped = length(seqs) + length(.absent),
                    n_aptamers_absent_from_archive = length(.absent), n_aptamers_skipped = length(skip_why),
                    n_assays = length(seqs) + length(.absent), n_assays_evaluable = 0L,
                    evaluable_class = "none",
                    skip_reasons = paste(sort(unique(skip_why)), collapse = "|"))]
    out_rows[[k]] <- base_row; next
  }

  pa <- rbindlist(per_assay)

  # intersection-union across aptamers: every one must agree in direction, and each
  # statistic takes its least favourable value.
  dirs <- sign(pa$aric_b)
  conc <- if (length(unique(dirs)) == 1L) "ROBUST" else "DISCORDANT"
  # EVERY statistic takes its least favourable value independently -- the SE is the LARGEST, not
  # whichever SE happened to sit beside the smallest |beta| (objective review 2026-08-06, finding 3).
  # Pairing them re-imported a single assay's optimism through the back door: for FLT1~CAD that
  # emitted SE 0.2309 where the least favourable is 1.2555, across 8 multi-assay pairs.
  # NO SINGLE ESTIMATE IS EMITTED FOR A MULTI-APTAMER PAIR (statistical + adversarial editor,
  # 2026-08-07, both vetoing the previous rule; the deCODE arm is treated identically so the two
  # arms remain comparable).
  #
  # Round 1 took least-|beta| and largest SE INDEPENDENTLY, emitting a beta from one aptamer beside
  # an SE from another -- for FLT1~CAD, b = 1.0185 (own SE 0.2309) beside se = 1.2555 from the other
  # aptamer. That pair belongs to no estimator, so any interval drawn from it is invalid, and it
  # affected all eight current multi-aptamer pairs. Round 2 replaced it with "the least favourable
  # aptamer, chosen by largest delta p, supplies beta and SE together" -- MY OWN INVENTION, which
  # both editors vetoed: delta p and the reported beta/se are monotone functions of the same z, so
  # choosing the aptamer BY its delta p and then reporting that aptamer's interval is
  # result-dependent selection. It destroys the interval's coverage, and in the direction that makes
  # the pair look worse, which is a different estimand rather than conservatism.
  #
  # So the estimate fields are NA and the observed RANGE is emitted instead. Consumers needing a
  # plottable interval draw one line PER APTAMER from the assay-level output, and must ERROR on
  # n_assays > 1 rather than silently take row 1. The IUT statistics (largest p, largest delta p,
  # smallest F) remain independent least-favourable values: taking a least-favourable TEST
  # STATISTIC is what the rule licenses; assembling an ESTIMATE out of parts is not.
  # THE SUPPRESSION APPLIES TO MULTI-APTAMER PAIRS ONLY (review finding 3, 2026-08-07). Unlike the
  # deCODE arm this script has no separate single-assay branch, so blanking the estimate fields
  # unconditionally erased the estimate of EVERY one-aptamer pair -- 40 of them on the current
  # archive. The reason those fields are NA for several aptamers is that no single number
  # represents several estimators; with one aptamer there is exactly one estimator and its own
  # beta, se and delta se are the honest answer.
  .b_rng <- range(pa$aric_b, na.rm = TRUE)
  .one <- nrow(pa) == 1L
  agg <- data.table(aric_b  = if (.one) pa$aric_b[1]  else NA_real_,
                    aric_se = if (.one) pa$aric_se[1] else NA_real_,
                    aric_b_min = .b_rng[1], aric_b_max = .b_rng[2],
                    aric_pval_iut = max(pa$aric_pval),
                    delta_se = if (.one) pa$delta_se[1] else NA_real_,
                    delta_pval_iut = max(pa$delta_pval), aric_F_min = min(pa$aric_F_stat))
  # Cochran's Q across aptamers -- a LOWER BOUND on heterogeneity, because the aptamers are
  # estimated against the SAME outcome GWAS and their errors are therefore positively correlated.
  # A large Q is informative; a small Q is not evidence of agreement.
  # WEIGHTED BY THE DELTA SE, NOT THE FIRST-ORDER SE (review finding 4). The first-order MR SE
  # treats the exposure estimate as known, so exposure uncertainty was being charged to
  # between-aptamer heterogeneity. Both arms compute a second-order (delta) SE precisely because
  # that assumption is false; Q must use it. Measured on deCODE FLT1~CAD, where one assay has
  # F = 2.36: first-order gives Q = 4.624, df = 1, p = 0.0315; the delta SE gives Q = 0.569,
  # p = 0.451. The p-value remains UNCALIBRATED regardless, because the aptamers share one outcome
  # GWAS and their errors are positively correlated -- Q is reported as a LOWER BOUND on
  # heterogeneity and its p-value must not be read as a test.
  .q <- local({
    b <- pa$aric_b; se <- pa$delta_se
    ok <- is.finite(b) & is.finite(se) & se > 0
    if (sum(ok) < 2L) c(Q = NA_real_, df = NA_real_, p = NA_real_) else {
      w <- 1 / se[ok]^2; bb <- sum(w * b[ok]) / sum(w)
      Qv <- sum(w * (b[ok] - bb)^2); dfv <- sum(ok) - 1L
      c(Q = Qv, df = dfv, p = stats::pchisq(Qv, dfv, lower.tail = FALSE))
    }
  })

  # exposure-side direction, aligned per aptamer above; ALL aptamers must agree (IUT)
  dir_conc <- if (all(!is.na(pa$bx_dir_concordant))) all(pa$bx_dir_concordant) else NA
  weak     <- agg$aric_F_min < 10
  # `n_assays` COUNTS APTAMERS ATTEMPTED, not those that produced an estimate; evaluable_class is
  # three-valued over that denominator. high_confidence is an INTERSECTION and is forced FALSE
  # unless EVERY attempted aptamer was evaluable -- a FALSE arising that way means the pair was not
  # evaluable, and must never be reported as a failure to replicate.
  .n_att <- nrow(pa) + length(.absent) + length(skip_why)
  .n_ev  <- nrow(pa)
  .ev_class <- if (.n_ev == 0L) "none" else if (.n_ev == .n_att) "all" else "some"
  hi_conf  <- .ev_class == "all" && conc == "ROBUST" && isTRUE(dir_conc) && !weak &&
              isTRUE(agg$delta_pval_iut < 0.05)

  # Direction disagreement is only a replication FAILURE when both sides are
  # non-null; a significant estimate against a near-zero noise estimate is an
  # absence of signal in one cohort, not a conflict. Below F >= 10 this arm reports
  # not-evaluable rather than failure.
  st <- if (!isTRUE(dir_conc) && !weak && isTRUE(agg$delta_pval_iut < 0.05)) "NOT_REPLICATED"
        else if (weak) "NOT_EVALUABLE_IN_ARIC_WEAK_F"
        else if (isTRUE(dir_conc)) "REPLICATED"
        else "NOT_EVALUABLE_DIRECTION_UNRESOLVED"

  base_row[, `:=`(aric_b = agg$aric_b, aric_se = agg$aric_se,
                  aric_b_min = agg$aric_b_min, aric_b_max = agg$aric_b_max,
                  aric_pval_iut = agg$aric_pval_iut,
                  delta_se = agg$delta_se, delta_pval_iut = agg$delta_pval_iut,
                  aric_F_min = agg$aric_F_min,
                  assay_Q = unname(.q[[1]]), assay_Q_df = unname(.q[[2]]),
                  assay_Q_pval = unname(.q[[3]]),
                  n_assays_evaluable = .n_ev, evaluable_class = .ev_class,
                  # `aric_b` is an MR RATIO; the old `ukb_b` beside it was the UKB EXPOSURE
                  # coefficient -- different units, so the pair could not be read as two exposure
                  # effects nor as two MR estimates (finding 4). The comparable pair is the aligned
                  # exposure beta on each side, which the direction verdict already uses; both are
                  # emitted under names that say which is which, least favourable across aptamers.
                  # Same rule as the estimate: with ONE aptamer these are that aptamer's own
                  # exposure betas; with several, naming one re-imports the selection the editors
                  # vetoed, so they are NA and the per-aptamer values live in the assay-level output.
                  bx_aric = if (.one) pa$bx_aric[1] else NA_real_,
                  bx_ukb  = if (.one) pa$bx_ukb[1]  else NA_real_,
                  # PROXY PROVENANCE SURVIVES AGGREGATION. Unlike bx_*, this is not an estimate
                  # whose selection across aptamers would re-import a vetoed choice: it is a
                  # property of the evidence, and the least favourable value across aptamers is the
                  # honest one. Any aptamer resting on a proxy makes the pair "mixed" at best.
                  outcome_proxy_used = if (all(pa$outcome_proxy_used == "direct")) "direct"
                                       else if (all(pa$outcome_proxy_used == "proxy")) "proxy" else "mixed",
                  outcome_proxy_r2 = if (all(is.na(pa$outcome_proxy_r2))) NA_real_ else min(pa$outcome_proxy_r2, na.rm = TRUE),
                  direction_concordant = dir_conc, weak_instrument = weak,
                  high_confidence = hi_conf, n_assays = .n_att,
                  assay_concordance = conc, status = st,
                  # A FAILED APTAMER IS RECORDED EVEN WHEN ANOTHER ONE SUCCEEDS (round 2, finding 4).
                  # `skip_why` was written only in the all-failed branch, so a pair with one good
                  # aptamer and one `mr_no_estimate` emitted REPLICATED / n_assays=1 /
                  # skip_reasons=NA -- indistinguishable from a legitimate single-aptamer evaluation,
                  # with the fit failure gone from both the summary and the assay-level RDS.
                  n_aptamers_mapped = length(seqs) + length(.absent),
                  n_aptamers_absent_from_archive = length(.absent),
                  n_aptamers_skipped = length(skip_why),
                  skip_reasons = if (length(skip_why)) paste(sort(unique(skip_why)), collapse = "|") else NA_character_)]
  out_rows[[k]] <- base_row
}

all_summary <- rbindlist(out_rows, fill = TRUE)
if (nrow(all_summary) != nrow(strong))
  fail_loud(sprintf("summary has %d rows for %d strong pairs -- a pair was dropped silently",
                    nrow(all_summary), nrow(strong)))

# evaluable_pair MIRRORS THE deCODE ARM (decode/01:957-960), so a downstream consumer can read one
# field name across both arms instead of knowing that one folds F into `status` and the other does
# not. DERIVED AND ADDITIVE: `status` is untouched, so nothing that already reads this file changes
# meaning. This arm's `status` does gate on F, so the column is expected to agree exactly with
# status %in% {REPLICATED, NOT_REPLICATED}; the check below is what makes that expectation testable
# rather than assumed -- if the two definitions are NOT equivalent, that must surface here and not
# be absorbed by quietly taking whichever number is larger.
all_summary[, evaluable_pair := evaluable_class == "all" &
                                status %chin% c("REPLICATED", "NOT_REPLICATED") &
                                !is.na(direction_concordant) &
                                is.finite(aric_F_min) & aric_F_min >= 10]
local({
  .by_status <- sum(all_summary$status %chin% c("REPLICATED", "NOT_REPLICATED"))
  .by_column <- sum(all_summary$evaluable_pair %in% TRUE)
  if (.by_status != .by_column)
    fail_loud(sprintf("ARIC: evaluable_pair (%d) disagrees with status-based evaluability (%d) -- the two definitions are not equivalent; resolve before using either number.",
                      .by_column, .by_status))
})

cat("\n-- status over the strong set --\n")
print(all_summary[, .N, by = status][order(-N)])
cat(sprintf("\ntestable: %d | evaluable_pair: %d | high-confidence: %d\n",
            sum(all_summary$status %chin% c("REPLICATED", "NOT_REPLICATED")),
            sum(all_summary$evaluable_pair %in% TRUE),
            sum(all_summary$high_confidence %in% TRUE)))

# WRITE THEN RENAME (review finding 16). Writing straight to the final path leaves a parseable
# PREFIX of the rows if the run is interrupted, and a short summary is indistinguishable from one
# where the missing pairs genuinely had no ARIC coverage -- there is no completion marker to check.
# The deCODE arm already does this; the two arms now fail the same way.
local({
  .fin <- file.path(ARIC_OUT_DIR, "aric_replication_summary.csv")
  .tmp <- paste0(.fin, ".", Sys.getpid(), ".tmp")
  fwrite(all_summary, .tmp)
  if (!file.exists(.tmp)) stop("ARIC: the temporary summary was not written")
  if (!file.rename(.tmp, .fin)) { unlink(.tmp); stop("ARIC: could not move the summary into place") }
})
if (length(mr_rows)) saveRDS(rbindlist(mr_rows, fill = TRUE), file.path(ARIC_OUT_DIR, "aric_mr_results.rds"))
cat(sprintf("\nwrote %s\n=== ARIC replication done ===\n",
            file.path(ARIC_OUT_DIR, "aric_replication_summary.csv")))
