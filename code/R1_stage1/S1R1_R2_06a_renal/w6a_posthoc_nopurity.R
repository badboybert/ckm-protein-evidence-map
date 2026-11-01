# W6a post hoc sensitivity (S2-D45, USER 2026-10-21; Stage 6.6 checker finding M3): the kidney comparison with the multi-signal candidates
# defined WITHOUT the credible-set purity filter. DECLARED POST HOC: written after the primary and the checker's note were seen.
#
# Definition. The primary kept a candidate when its SuSiE PP.H4 stayed >= 0.80 after a refit with credible-set purity filtering
# (min |r| 0.5) with converged fits (Additional file 2, SM candidates). The map's own multi-signal class uses no purity filter.
# Here a candidate is an association whose FIRST fit, without the purity filter, gave PP.H4 >= 0.80 with a converged fit:
#   - the 16 Drawn rows of the edge audit (kept by the primary);
#   - the 9 Dropped rows (4 fell below 0.80 under the purity filter; 5 had no pure credible set). Their audit flag
#     'SuSiE converged (no warning)' = no belongs to the purity refit (the collapse), not to the first fit;
#   - NOT ZFYVE19-eGFR (Held): it reached 0.80 only with a non-convergent fit, which this definition still excludes.
# Pool, kidney flag (eGFR or CKD) and estimator are the primary's: the 160-pair pool written by w6a_compute.R, the conditional-MLE
# odds ratio with its exact two-sided CI (exact_cond_or.R) and the fisher.test p. The primary is reproduced first (control).
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
suppressPackageStartupMessages(library(data.table))
HERE <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal")
IN <- list(pool = file.path(HERE, "w6a_pool160_genomic_clusters.tsv"),
           audit = paste0(.P31R_V3, "/V4/source_data/supp_tables/SuppTable_candidate_edge_audit.csv"),
           ecor = file.path(HERE, "exact_cond_or.R"),
           prim = file.path(HERE, "results.tsv"))
MD5 <- c(pool = "c1e1426892d7b54443399befe2ac98ae", audit = "aa6baa01970fd92f009d1f2cf00d1741", ecor = "1217326d2a6880b174f91d377c5b21b9")
md5f <- function(p) unname(tools::md5sum(p))
for (k in names(MD5)) if (!identical(md5f(IN[[k]]), MD5[[k]])) stop("INPUT md5 MISMATCH: ", k, " ", IN[[k]])
source(IN$ecor)
pool <- fread(IN$pool)
stopifnot(nrow(pool) == 160L, sum(pool$cand) == 16L)
aud <- fread(IN$audit, encoding = "UTF-8")
stopifnot(nrow(aud) == 26L)
pk <- function(g, o) paste(g, o, sep = "|")
aud[, key := pk(Gene, Endpoint)]
pool[, key := pk(gene, outcome_name)]
stopifnot(all(aud$key %in% pool$key))                                    # every audited association is in the 160-pair pool
stopifnot(setequal(aud[Status == "Drawn", key], pool[cand == TRUE, key])) # the primary's 16 are the audit's Drawn rows
naive_ok <- suppressWarnings(as.numeric(aud$`SuSiE PP.H4 (naive)`)) >= 0.80
stopifnot(all(naive_ok))                                                 # all 26 reached 0.80 in the first fit (the audit's entry criterion)
stopifnot(sum(aud$Status == "Dropped") == 9L, sum(aud$Status == "Held (not drawn)") == 1L, aud[Status == "Held (not drawn)", key] == "ZFYVE19|eGFR")
nf <- aud[Status %in% c("Drawn", "Dropped"), key]
stopifnot(length(nf) == 25L)
pool[, cand_nf := key %in% nf]

tab4 <- function(cand, kid) matrix(c(sum(cand & kid), sum(!cand & kid), sum(cand & !kid), sum(!cand & !kid)), 2,
                                   dimnames = list(c("cand", "other"), c("kidney", "nonkidney")))
fit <- function(tab) {
  e <- exact_cond_or(tab)
  list(a = tab[1, 1], b = tab[1, 2], c = tab[2, 1], d = tab[2, 2], or = e$estimate, lo = e$conf.int[1], hi = e$conf.int[2],
       p = fisher.test(tab)$p.value)
}
# control: reproduce the primary from its own outputs
prim <- fit(tab4(pool$cand, pool$kid))
R <- fread(IN$prim, sep = "\t", na.strings = c("NA", ""))
pv <- function(id, col = "value") as.numeric(R[quantity_id == id][[col]])
stopifnot(prim$a == pv("S6.primary.a"), prim$b == pv("S6.primary.b"), prim$c == pv("S6.primary.c"), prim$d == pv("S6.primary.d"))
stopifnot(abs(prim$or - pv("S6.primary.or_mle")) < 1e-6, abs(prim$lo - pv("S6.primary.or_mle", "ci_lo")) < 1e-6,
          abs(prim$hi - pv("S6.primary.or_mle", "ci_hi")) < 1e-6, abs(prim$p - pv("S6.primary.p")) < 1e-9)
cat(sprintf("CONTROL reproduced: %d/%d vs %d/%d, OR %.6f (%.6f to %.6f), p %.10g\n", prim$a, prim$a + prim$b, prim$c, prim$c + prim$d, prim$or, prim$lo, prim$hi, prim$p))
# the post hoc sensitivity
s <- fit(tab4(pool$cand_nf, pool$kid))
cat(sprintf("NO-FILTER: %d/%d vs %d/%d, OR %.6f (%.6f to %.6f), p %.10g\n", s$a, s$a + s$b, s$c, s$c + s$d, s$or, s$lo, s$hi, s$p))
out <- data.table(quantity_id = c("S6.posthoc_nopurity.a", "S6.posthoc_nopurity.b", "S6.posthoc_nopurity.c", "S6.posthoc_nopurity.d",
                                  "S6.posthoc_nopurity.or_mle", "S6.posthoc_nopurity.p"),
                  prespec_ref = "post hoc (S2-D45, 2026-10-21)",
                  estimand = "kidney (eGFR or CKD) share among multi-signal candidates defined without the purity filter (first fit PP.H4 >= 0.80, converged; 25 pairs) vs the other pool pairs; pool = 160",
                  value = c(s$a, s$b, s$c, s$d, s$or, s$p), ci_lo = c(NA, NA, NA, NA, s$lo, NA), ci_hi = c(NA, NA, NA, NA, s$hi, NA),
                  p = c(NA, NA, NA, NA, s$p, s$p),
                  notes = c("candidate kidney pairs", "candidate non-kidney pairs", "other kidney pairs", "other non-kidney pairs",
                            "conditional-MLE odds ratio, exact two-sided 95% CI (exact_cond_or.R)", "fisher.test two-sided p"))
fwrite(out, file.path(HERE, "w6a_posthoc_nopurity.tsv"), sep = "\t")
cat("wrote w6a_posthoc_nopurity.tsv; inputs md5:", paste(names(MD5), MD5, collapse = "; "), "\n")
