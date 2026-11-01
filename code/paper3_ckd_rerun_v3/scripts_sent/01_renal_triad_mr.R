# ============================================================================
# Reconstruct the renal-triad MR orphan egfr_mirage_classification.csv (hand-run; Stage-B audit
# PB-AUD-003 / GM-F1) from scripted sources. NO new method — single-SNP Wald, which is exactly what
# the cascade already does for eGFR (verified: all 71 proteins are single-instrument; crea == cascade
# eGFR Wald in mr_all_results.rds).
#   crea  <- cascade eGFR(Wuttke2019) Wald  (mr_all_results.rds; reused, not recomputed)
#   cys   <- single-SNP Wald vs eGFRcys  (CKDGen_eGFRcys_Stanzick2021.gz)
#   bun   <- single-SNP Wald vs BUN      (CKDGen_BUN_EA_Wuttke2019.txt.gz)
#   set   <- proteins with crea_p < 1e-3 (the GATE-1 n=71 set)
#   class <- concordant (cys sig, same sign) / mirage (cys p>=0.05) / other (cys sig, opposite sign)
# Output (env MIRAGE_OUT overrides; default = canonical):
#   results/repro_corrected/pam_conditional/egfr_mirage_classification.csv
# Run: Rscript scripts/_shared/renal_triad/01_renal_triad_mr.R
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
# _phase1_set.R defines PHASE1_SET; PHASE1_TAG belongs to 00_config.R, which this script does
# not source -- the reference stopped a direct launch at object-not-found.
RC <- file.path(paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/results"), paste0("repro_", PHASE1_SET))
KM <- paste0(.P31R_EXT, "/E/20260331-CKM/Paper 3/data/outcome_gwas/ckdgen_kidney_markers")
OUT <- Sys.getenv("MIRAGE_OUT", file.path(RC, "pam_conditional", "egfr_mirage_classification.csv"))

instr <- as.data.table(readRDS(PHASE1_INSTRUMENTS))
# eaf CARRIED, because the palindrome test below needs it. Without a frequency this script could only
# align on allele letters, and for a palindromic SNP the letters match in BOTH orientations.
.eafcol <- intersect(c("eaf.exposure", "eaf", "A1FREQ", "eaf_exposure"), names(instr))[1]
if (is.na(.eafcol))
  stop("01_renal_triad_mr: the instrument table carries no effect-allele frequency column ",
       "(looked for eaf.exposure / eaf / A1FREQ) -- the palindrome test cannot run without one")
instr <- instr[, .(gene, SNP, ea = toupper(effect_allele), oa = toupper(other_allele), beta_exp = beta,
                   eaf_exp = as.numeric(get(.eafcol)))]
instr <- instr[!duplicated(gene)]                       # 1 instrument per gene (verified)

mr <- as.data.table(readRDS(file.path(RC, "mr_all_results.rds")))
crea <- mr[outcome_name == "eGFR" & method == "Wald ratio", .(gene, crea_b = b, crea_p = pval)]
# SELECTED AT THE LOOSEST THRESHOLD, CLASSIFIED ONCE, PUBLISHED TWICE. The cystatin-scrub figure's
# panel (b) reports the mirage fraction at FOUR crea thresholds (Bonferroni / 1e-4 / 1e-3 / 0.05),
# but this script only computed cystatin/BUN Wald for the Gate-1 set (crea_p < 1e-3) -- so the
# loosest threshold's bar (n~248) described a computation that never ran, and all four bars were
# hardcoded in the figure script. Widening the SELECTION to p < 0.05 makes every threshold's set a
# subset of what is actually computed; the canonical Gate-1 artifact below keeps its exact old
# semantics because the renal figure and its Wilcoxon test consume it.
GATE1_P <- 1e-3
SWEEP_P <- 0.05
sel <- crea[crea_p < SWEEP_P]
dt <- merge(sel, instr, by = "gene")

# single-SNP Wald vs an outcome sumstat, harmonised to the pQTL effect allele
wald_vs <- function(dt, sumfile, a1, a2, eff, sef, rs, frq = "Freq1") {
  s <- fread(sumfile, showProgress = FALSE)
  if (!frq %in% names(s))
    stop(sprintf("01_renal_triad_mr: %s has no %s column -- strand for a palindromic SNP cannot be resolved without a frequency",
                 basename(sumfile), frq))
  s <- s[get(rs) %in% dt$SNP,
         .(SNP = get(rs), o1 = toupper(get(a1)), o2 = toupper(get(a2)),
           oeff = as.numeric(get(eff)), osef = as.numeric(get(sef)),
           ofrq = as.numeric(get(frq)))]
  s <- s[!duplicated(SNP)]
  m <- merge(dt[, .(gene, SNP, ea, oa, beta_exp, eaf_exp)], s, by = "SNP", all.x = TRUE)
  m[, align := fifelse(ea == o1 & oa == o2, 1, fifelse(ea == o2 & oa == o1, -1, NA_real_))]
  # AN UNRESOLVABLE PALINDROME IS NOT AN ALIGNED ROW. Alignment above uses allele LETTERS only, and for
  # a palindromic SNP (A/T, C/G) the letters match in both orientations -- so the sign of `b` was being
  # decided by which orientation the outcome file happened to print. The frequency is the only
  # discriminator, and it discriminates nothing when both sides sit at ~0.5. Same rule and same
  # thresholds as _local_outcome_extract.R:251-252, which the discovery screen already applies.
  # MEASURED TRIGGER on current data: TWF2 / rs7622851, C/G, exposure eaf 0.497754, outcome Freq1
  # 0.5019. It was emitted as cys_b = +0.04232 and called `concordant`; the opposite -- equally
  # supported -- orientation gives -0.0423 and `other`. Dropping it takes the renal arm from 90 rows
  # (53 concordant / 35 mirage / 2 other) to 89 (52 / 35 / 2), and that count change is the point:
  # the row was never evidence for either class.
  .cmpl <- c(A = "T", T = "A", C = "G", G = "C")
  m[, .palindromic := !is.na(ea) & !is.na(oa) & ea %chin% names(.cmpl) & unname(.cmpl[ea]) == oa]
  m[, .unresolvable := .palindromic &
       (is.na(eaf_exp) | (eaf_exp > 0.42 & eaf_exp < 0.58) |
        is.na(ofrq)    | (ofrq    > 0.42 & ofrq    < 0.58))]
  if (m[.unresolvable == TRUE, .N])
    cat(sprintf("[renal_triad] %s: dropping %d palindromic SNP(s) whose strand the frequencies cannot resolve: %s\n",
                basename(sumfile), m[.unresolvable == TRUE, .N],
                paste0(m[.unresolvable == TRUE, paste0(gene, "/", SNP, " ", ea, "/", oa,
                                                       " eaf=", signif(eaf_exp, 4), " frq=", signif(ofrq, 4))],
                       collapse = "; ")))
  m[.unresolvable == TRUE, align := NA_real_]
  m[, b := (oeff * align) / beta_exp]
  m[, se := abs(osef / beta_exp)]
  m[, p := 2 * pnorm(-abs(b / se))]
  setkey(m, gene)
  m[, .(gene, b, se, p)]
}

cys <- wald_vs(dt, file.path(KM, "CKDGen_eGFRcys_Stanzick2021.gz"), "Allele1", "Allele2", "Effect", "StdErr", "RSID")
bun <- wald_vs(dt, file.path(KM, "CKDGen_BUN_EA_Wuttke2019.txt.gz"), "Allele1", "Allele2", "Effect", "StdErr", "RSID")

out <- merge(dt[, .(gene, crea_b, crea_p)], cys[, .(gene, cys_b = b, cys_p = p)], by = "gene")
out <- merge(out, bun[, .(gene, bun_b = b, bun_p = p)], by = "gene")
# the mirage set requires eGFRcys to be testable (lead SNP present in the cystatin GWAS); proteins whose
# lead SNP is absent from eGFRcys cannot be classified mirage/concordant and are excluded (-> n=71).
out <- out[!is.na(cys_b)]
out[, mr_egfr_p := crea_p]
out[, class := fifelse(cys_p >= 0.05, "mirage",
                fifelse(sign(cys_b) == sign(crea_b), "concordant", "other"))]
setcolorder(out, c("gene", "crea_b", "crea_p", "cys_b", "cys_p", "bun_b", "bun_p", "mr_egfr_p", "class"))
setorder(out, gene)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
# THE CANONICAL FILE KEEPS ITS GATE-1 SEMANTICS EXACTLY. Consumers (the renal figure, its Wilcoxon
# comparison, cystatin-scrub panel (c)) were built against the crea_p < 1e-3 set; publishing the
# widened selection there would silently change every count they draw.
out_gate1 <- out[crea_p < GATE1_P]
fwrite(out_gate1, OUT)
cat(sprintf("wrote %s: %d rows (Gate-1, crea_p < %g)\n", OUT, nrow(out_gate1), GATE1_P))
print(out_gate1[, .N, by = class])

# THE THRESHOLD SWEEP, for cystatin-scrub panel (b) to DERIVE from instead of hardcoding. Bonferroni
# is derived from the analysed universe (proteins x outcomes), not written as a literal -- the figure
# used to say 3.7e-6, the PREVIOUS generation's alpha. Each threshold's denominator is the
# classifiable set at that threshold (lead SNP present in the cystatin GWAS), and the fraction is
# mirage / classifiable -- the same rule at every threshold, on actually computed rows.
.bonf <- 0.05 / (uniqueN(instr$gene) * uniqueN(mr$outcome_name))
.thr <- data.table(threshold_label = c(sprintf("Bonferroni (%.3g)", .bonf), "p<1e-4", "p<1e-3", "p<0.05"),
                   threshold = c(.bonf, 1e-4, 1e-3, 0.05))
sweep <- .thr[, {
  .set <- out[crea_p < threshold]
  .(n_classifiable = nrow(.set),
    n_mirage       = .set[class == "mirage", .N],
    pct_mirage     = if (nrow(.set)) round(100 * .set[class == "mirage", .N] / nrow(.set), 1) else NA_real_)
}, by = .(threshold_label, threshold)]
SWEEP_OUT <- file.path(dirname(OUT), "egfr_mirage_threshold_sweep.csv")
fwrite(sweep, SWEEP_OUT)
cat(sprintf("wrote %s:\n", SWEEP_OUT)); print(sweep)
