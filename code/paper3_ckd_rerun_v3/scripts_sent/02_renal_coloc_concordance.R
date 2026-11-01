# ============================================================================
# Reconstruct kidney_cystatin_bun_concordance.csv (hand-run orphan; Stage-B audit PB-AUD-003) from
# scripted sources. Same verified single-SNP renal-triad Wald as 01_renal_triad_mr.R, restricted to the
# renal coloc-supported pairs (outcome in {eGFR, CKD}) + CKD Wald + SNP/tier.
#   row set     <- table_s_pair_evidence_tiered.csv rows with outcome in (eGFR, CKD)  -> 15 pairs (frozen 52-set)
#   egfrcrea/ckd<- cascade Wald (mr_all_results.rds; eGFR=Wuttke2019, CKD=Wuttke2019_CKD_EA) reused
#   cys/bun     <- single-SNP Wald vs CKDGen_eGFRcys_Stanzick2021 / CKDGen_BUN_EA_Wuttke2019 (blank if
#                  lead SNP absent from that GWAS; coloc pairs are kept regardless, unlike the mirage set)
# Outputs the columns CONSUMED downstream (build_supptables_from_sources.py reads gene/outcome/
# evidence_tier/SNP/egfrcrea_*/cys_*/bun_*/ckd_*); the original file's vestigial cols (mr_beta,
# mr_direction, prim_b, exp_cys, and the DEPRECATED stored cys_conc/bun_conc that mis-score CKD pairs and
# are ignored downstream in favour of the sign rule) are intentionally dropped.
# Output (env KIDNEY_OUT overrides; default = canonical):
#   results/eas_replication/kidney_cystatin_bun_concordance.csv
# Run: Rscript scripts/_shared/renal_triad/02_renal_coloc_concordance.R
# ============================================================================
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/R/_phase1_set.R"))
suppressMessages(library(data.table))
RC  <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), paste0("repro_", PHASE1_SET))   # L7 補遺 2026-08-16: _phase1_set.R 定義的是 PHASE1_SET（PHASE1_TAG 在 00_config 才有，獨立執行時未定義）
KM  <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/outcome_gwas/ckdgen_kidney_markers")
EAS <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results/eas_replication")
OUT <- Sys.getenv("KIDNEY_OUT", file.path(EAS, "kidney_cystatin_bun_concordance.csv"))

instr <- as.data.table(readRDS(PHASE1_INSTRUMENTS))
instr <- instr[, .(gene, SNP, ea = toupper(effect_allele), oa = toupper(other_allele), beta_exp = beta)][!duplicated(gene)]
mr <- as.data.table(readRDS(file.path(RC, "mr_all_results.rds")))
crea <- mr[outcome_name == "eGFR" & method == "Wald ratio", .(gene, egfrcrea_b = b, egfrcrea_p = pval)]
ckd  <- mr[outcome_name == "CKD"  & method == "Wald ratio", .(gene, ckd_b = b, ckd_p = pval)]

tier <- fread(file.path(RC, "table_s_pair_evidence_tiered.csv"))
ren <- tier[outcome %in% c("eGFR", "CKD"), .(gene, outcome, evidence_tier)]   # 15 renal coloc pairs (frozen 52-set)

wald_vs <- function(genes_dt, sumfile, a1, a2, eff, sef, rs, pv) {   # single-SNP Wald ratio; harmonised
  s <- fread(sumfile, showProgress = FALSE)
  s <- s[get(rs) %in% genes_dt$SNP,
         .(SNP = get(rs), o1 = toupper(get(a1)), o2 = toupper(get(a2)),
           oeff = as.numeric(get(eff)), osef = as.numeric(get(sef)), op = as.numeric(get(pv)))][!duplicated(SNP)]
  m <- merge(genes_dt[, .(gene, SNP, ea, oa, beta_exp)], s, by = "SNP", all.x = TRUE)
  m[, align := fifelse(ea == o1 & oa == o2, 1, fifelse(ea == o2 & oa == o1, -1, NA_real_))]
  # p = the source GWAS two-sided P-value for the SNP-outcome association (editor Finding 5, 2026-07-18):
  # dividing the SNP effect by the constant beta_exp to form the Wald ratio leaves the two-sided p unchanged,
  # so we take the GWAS's own P-value column instead of recomputing it from the ROUNDED Effect/StdErr, which
  # produced artifacts (DNER's Effect rounds to 0 -> recomputed p=1.0 vs the true source p=0.96; C2CD2L/HCLS1
  # collided on an identical recomputed p). NA when the SNP is absent/unharmonisable (align NA), preserving
  # the "blank if lead SNP absent from that GWAS" semantics. The sign-based concordance is unaffected (it uses
  # b's sign, not p), so the 13/14 cystatin and 13/15 BUN counts do not change.
  m[, .(gene, b = (oeff * align) / beta_exp, p = fifelse(is.na(align), NA_real_, op))]
}

gset <- merge(unique(ren[, .(gene)]), instr, by = "gene")
cys <- wald_vs(gset, file.path(KM, "CKDGen_eGFRcys_Stanzick2021.gz"), "Allele1", "Allele2", "Effect", "StdErr", "RSID", "P.value")
bun <- wald_vs(gset, file.path(KM, "CKDGen_BUN_EA_Wuttke2019.txt.gz"), "Allele1", "Allele2", "Effect", "StdErr", "RSID", "P-value")

out <- merge(ren, instr[, .(gene, SNP)], by = "gene")
out <- merge(out, crea, by = "gene", all.x = TRUE)
out <- merge(out, ckd,  by = "gene", all.x = TRUE)
out <- merge(out, cys[, .(gene, cys_b = b, cys_p = p)], by = "gene", all.x = TRUE)
out <- merge(out, bun[, .(gene, bun_b = b, bun_p = p)], by = "gene", all.x = TRUE)
out <- out[, .(gene, outcome, evidence_tier, SNP, egfrcrea_b, egfrcrea_p, ckd_b, ckd_p, cys_b, cys_p, bun_b, bun_p)]
setorder(out, gene, outcome)
fwrite(out, OUT)
cat(sprintf("wrote %s: %d rows (%d eGFR + %d CKD)\n", OUT, nrow(out), sum(out$outcome=="eGFR"), sum(out$outcome=="CKD")))
