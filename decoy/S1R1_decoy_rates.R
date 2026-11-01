# Decoy-locus stress test: pass rates with cluster-bootstrap intervals, per stratum and per outcome, and the CKD
# selection-matched illustration. Read-only. Own code (R), independent of the Python null_uncertainty.py of the R1 check.
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
suppressPackageStartupMessages(library(data.table))
V3 <- .P31R_V3
OD <- file.path(V3, "R1_stage1/S1R1_03_decoy_method")
F <- c(noneGFR_draws = "v3.2/_engine_root/results/v4_local/coloc_null/null_draws.csv",
       noneGFR_obs   = "v3.2/_engine_root/results/v4_local/coloc_null/null_observed.csv",
       eGFR_draws    = "v3.2/_engine_root/results/repro_phase1_hg38/coloc_null/null_draws_egfr.csv",
       eGFR_obs      = "v3.2/_engine_root/results/repro_phase1_hg38/coloc_null/null_observed_egfr.csv",
       ckd_draws     = "V4_local/source_data/coloc_null/ckd_stratum_draws.csv",
       ckd_annot     = "R1_stage1/S1R1_03_decoy_method/S1R1_F2_draws_annotated_from_skeptic1.csv")  # R1 skeptic1 annotation (evidence zip): sel = |aligned z_out| >= 1.96 at the recipient instrument row
L <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); L <<- c(L, s) }
for (k in names(F)) say("md5 %-14s %s  %s", k, unname(tools::md5sum(file.path(V3, F[[k]]))), F[[k]])
set.seed(20261002); NB <- 2000

boot_rate <- function(d, two_way = TRUE) {
  tg <- factor(paste(d$gene, d$outcome)); dn <- factor(d$donor)
  it <- as.integer(tg); id <- as.integer(dn); nt <- nlevels(tg); nd <- nlevels(dn)
  p <- as.numeric(d$pass)
  one <- function() {
    wt <- tabulate(sample.int(nt, nt, TRUE), nt)
    w <- wt[it]; if (two_way) w <- w * tabulate(sample.int(nd, nd, TRUE), nd)[id]
    sum(w * p) / sum(w)
  }
  quantile(replicate(NB, one()), c(.025, .975), names = FALSE)
}
summ <- function(d, obs, label) {
  rate <- mean(d$pass); tci <- boot_rate(d, FALSE); wci <- boot_rate(d, TRUE)
  data.table(stratum = label, targets_tested = sum(!is.na(obs$pass)), observed_passes = sum(obs$pass, na.rm = TRUE),
             draws = nrow(d), decoy_passes = sum(d$pass), rate = rate,
             target_boot_lo = tci[1], target_boot_hi = tci[2], twoway_boot_lo = wci[1], twoway_boot_hi = wci[2],
             median_H4 = median(d$H4), q95_H4 = quantile(d$H4, .95, names = FALSE))
}
dn <- fread(file.path(V3, F[["noneGFR_draws"]])); on <- fread(file.path(V3, F[["noneGFR_obs"]]))
de <- fread(file.path(V3, F[["eGFR_draws"]]));   oe <- fread(file.path(V3, F[["eGFR_obs"]]))
say("non-eGFR outcomes: %s", paste(sort(unique(dn$outcome)), collapse = ", "))
res <- rbind(summ(dn, on, "non-eGFR (seven outcomes, CKD included)"), summ(de, oe, "eGFR"))
for (o in sort(unique(dn$outcome))) res <- rbind(res, summ(dn[outcome == o], on[outcome == o], o))

# CKD selection-matched: |aligned z_out| >= 1.96 at the recipient's MR instrument row (sel), as the real MR screen selected
ck <- fread(file.path(V3, F[["ckd_annot"]]))
ckd0 <- fread(file.path(V3, F[["ckd_draws"]]))
stopifnot(nrow(ck) == nrow(ckd0), identical(ck$gene, ckd0$gene), identical(ck$donor, ckd0$donor), identical(as.logical(ck$pass), as.logical(ckd0$pass)))
ck[, `:=`(pass = as.logical(pass), pass_mine = as.logical(pass_mine), inst_kept = as.logical(inst_kept), sel = as.logical(sel), sel_p = as.logical(sel_p))]
say("CKD instrument absent from the target capture whatever the draw: %d draws; removed BY TRUNCATION (k < n_t, instrument index > k): %d (%.1f%%)",
    ck[!inst_kept & !(k < n_t), .N], ck[!inst_kept & k < n_t, .N], 100 * ck[!inst_kept & k < n_t, .N] / nrow(ck))
say("CKD idx (instrument row index) NA in %d draws; inst_kept FALSE with idx <= k: %d", ck[is.na(idx), .N], ck[!inst_kept & !is.na(idx) & idx <= k, .N])
say("CKD annotated draws: %d; pass column agrees with recomputed rule in %d of %d rows", nrow(ck), sum(ck$pass == ck$pass_mine), nrow(ck))
say("CKD draws that truncate the target vector: %d (%.1f%%); that cut a LOCATED instrument row (idx > k): %d (%.1f%%); targets whose instrument is located in their capture: %d of %d",
    sum(ck$k < ck$n_t), 100 * mean(ck$k < ck$n_t), ck[!is.na(idx) & idx > k, .N], 100 * ck[!is.na(idx) & idx > k, .N] / nrow(ck),
    uniqueN(ck[!is.na(idx), gene]), uniqueN(ck$gene))
say("CKD selection-matched draws come from %d targets (of %d with a located instrument)", uniqueN(ck[sel_p == TRUE, gene]), uniqueN(ck[!is.na(idx), gene]))
say("CKD decoy passes with the located instrument row cut away: %d of %d passes", ck[pass == TRUE & !is.na(idx) & idx > k, .N], ck[pass == TRUE, .N])
res <- rbind(res, summ(ck, on[outcome == unique(ck$outcome)], "CKD (all decoy draws)"),
             summ(ck[sel_p == TRUE], on[outcome == unique(ck$outcome)], "CKD, selection-matched, p < 0.05 at the located instrument (sel_p; 154 contributing targets)"),
             summ(ck[sel == TRUE], on[outcome == unique(ck$outcome)], "CKD, selection-matched at the instrument rsID only (|z| >= 1.96; same 154 targets)"))
num <- names(res)[vapply(res, is.numeric, TRUE)]
fwrite(res, file.path(OD, "S1R1_decoy_rates.tsv"), sep = "\t")
say("== pass rates (%%) with 95%% bootstrap intervals (target-cluster; two-way target x donor), B = %d", NB)
for (i in seq_len(nrow(res))) with(res[i], say("  %-62s %6d draws %4d passes  %.3f%% [target %.3f-%.3f; two-way %.3f-%.3f]  observed %d of %d",
  stratum, draws, decoy_passes, 100 * rate, 100 * target_boot_lo, 100 * target_boot_hi, 100 * twoway_boot_lo, 100 * twoway_boot_hi,
  observed_passes, targets_tested))
writeLines(L, file.path(OD, "S1R1_decoy_rates.log"))
