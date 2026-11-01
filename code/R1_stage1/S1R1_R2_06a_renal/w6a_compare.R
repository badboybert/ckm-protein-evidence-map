# W6a step 2 -- COMPARE_SEEN.tsv.  Runs only after results.tsv / branches.tsv exist; asserts they are the
# pre-compare versions (md5 pinned below from run.log, 2026-10-08T10:47:35+0800), then reads the FULL
# pre-declaration and the seen wpC outputs and writes one row per compared quantity.  Changes no Stage 1 value.
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
suppressMessages({ library(data.table); library(jsonlite) })
OUT <- Sys.getenv("W6A_OUT", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))   # override to re-run into another folder (e.g. verify_L3/rerun)
RESDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
LOGF <- file.path(OUT, "run.log")
lg <- function(...) { line <- paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] ", paste0(..., collapse = "")); cat(line, "\n", sep = "", file = LOGF, append = TRUE); cat(line, "\n") }
md5f <- function(p) unname(as.character(tools::md5sum(p)))
cat("\n==================== RUN ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " (w6a_compare.R) ====================\n", sep = "", file = LOGF, append = TRUE)

# FIXER F1 (FIXLOG.md): results.tsv / branches.tsv were regenerated AFTER the full text had been opened, to fix verification-panel
# defect F1 (root tolerance of the exact CI). The original pre-compare versions are preserved in fix_F1_prefix/ and must still carry
# the original pins; fix_F1_diffgate.R (run by w6a_run_all.R before this step) proves the regenerated files differ from them only in
# the Fisher OR/CI cells. PRECOMPARE below pins the post-fix versions.
SRC <- Sys.getenv("W6A_SRC", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))
PRECOMPARE_ORIG <- c(results = "724f202391469efa60d670c22b07eeb3", branches = "7da7a8002a67f7ed8eea9ae3223a1bbb")
orig <- c(results = md5f(file.path(SRC, "fix_F1_prefix", "results.tsv")), branches = md5f(file.path(SRC, "fix_F1_prefix", "branches.tsv")))
if (!identical(orig, PRECOMPARE_ORIG)) stop("fix_F1_prefix/ no longer holds the original pre-compare results.tsv / branches.tsv")
PRECOMPARE <- c(results = "2a85c748f20495b74227cb4780a5f759", branches = "bfe38e2f73bf8f430547bf698f951132")   # post-fix F1
got <- c(results = md5f(file.path(OUT, "results.tsv")), branches = md5f(file.path(OUT, "branches.tsv")))
lg("COMPARE: results.tsv md5 ", got[["results"]], " / branches.tsv md5 ", got[["branches"]], " (pre-compare pins ", PRECOMPARE[["results"]], " / ", PRECOMPARE[["branches"]], ")")
if (!identical(got, PRECOMPARE)) stop("results.tsv or branches.tsv differs from the pre-compare version -- a value changed after the full text was opened")
FULL <- paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md")
stopifnot(md5f(FULL) == "3319dc39740d9f5b3b8822133d9108b7")
ft <- readLines(FULL, encoding = "UTF-8")
RES <- fromJSON(file.path(RESDIR, "RESOURCES.json"), simplifyVector = FALSE); rmap <- setNames(RES, vapply(RES, function(e) e$id, ""))
SEEN <- c("SEEN_wpC_pool160_reconstructed_csv", "SEEN_wpC_pool160_genomic_clusters_csv", "SEEN_wpC_renal_fisher_R_csv", "SEEN_wpC_renal_fisher_python_csv",
          "SEEN_wpC_drawn16_failure_modes_csv", "SEEN_wpC_susie_candidate_fit_ledger_from_logs_csv", "SEEN_wpC_pool_meta_json")
SP <- list()
for (id in SEEN) { e <- rmap[[id]]; m <- md5f(e$path); lg("SEEN OUTPUT read for COMPARE only: ", id, " md5 ", m, " (RESOURCES ", e$md5, ")"); stopifnot(identical(m, e$md5)); SP[[id]] <- e$path }

R <- fread(file.path(OUT, "results.tsv"), sep = "\t", colClasses = "character")
rv <- function(id) { v <- R[quantity_id == id, value]; if (length(v) != 1) stop("missing results id ", id); v }
rn <- function(id) as.numeric(rv(id))
CMP <- list()
add <- function(qid, seen_value, seen_label, stage1_value, stage1_label, explanation, dp = NA, rel = FALSE) {
  sv <- suppressWarnings(as.numeric(seen_value)); tv <- suppressWarnings(as.numeric(stage1_value))
  ad <- if (!is.na(sv) && !is.na(tv)) abs(sv - tv) else NA
  beyond <- if (is.na(ad)) NA else if (!is.na(dp)) ad > 0.5 * 10^(-dp) + 1e-12 else ad > 1e-9 * max(1, abs(sv))
  flag <- if (is.na(beyond)) "" else if (beyond) "DIFFERS beyond printed precision: " else "agrees at printed precision. "
  CMP[[length(CMP) + 1]] <<- data.table(quantity_id = qid, seen_value = as.character(seen_value), seen_label = seen_label,
    stage1_value = as.character(stage1_value), stage1_label = stage1_label, abs_diff = if (is.na(ad)) "" else sprintf("%.6g", ad),
    explanation = paste0(flag, explanation))
}
L <- function(n) sprintf("full text P31R_R2_STAGE1_PRESPEC.md l.%d (§6 Seen)", n)
chkline <- function(n, pat) stopifnot(grepl(pat, ft[n], fixed = TRUE))
chkline(485, "Pair test"); chkline(486, "Locus–outcome"); chkline(487, "Excluding GUSB"); chkline(488, "(a):"); chkline(489, "Next-lowest")
chkline(491, "Fisher variants"); chkline(493, "numpy, seed 20261007"); chkline(494, "R, same seed"); chkline(495, "R, 10 seeds"); chkline(496, "numpy, 12 seeds")
chkline(497, "26 of 160"); chkline(498, "9 / 6 / 1"); chkline(576, "seen: 1.003")
WEAK <- "Deterministic value; the count/target was visible in the redacted copy or DEF_SM_CANDIDATES before compute, so agreement is weak evidence."
DET <- "Deterministic value not visible before compute; computed independently in R from the W6a rebuild (p: fisher.test; OR and exact CI: fisher.test's equations solved at tol 1e-13, exact_cond_or.R, FIXLOG F1)."
FTX <- paste("OR / exact CI. The seen value is R fisher.test, whose uniroot runs at its default tolerance (about 1.2e-4 on psi or 1/psi), which mostly moves the UPPER bound in the 3rd decimal.",
             "Stage 1 first printed fisher.test's value too; verification panel L2 flagged it (defect F1), and the fixer re-solved fisher.test's own equations at tol 1e-13",
             "(agrees with scipy odds_ratio(kind='conditional') to 1e-13 relative; FIXLOG.md). A difference here is that root tolerance, not a different estimand or table.",
             "The fix was made after the full text had been opened, on the verifier's evidence; it moves Stage 1 away from the seen value, not toward it.")

# ---- text Seen values (full pre-declaration §6 / §11)
add("S6.primary.p", "0.0329", L(485), rv("S6.primary.p"), "results.tsv S6.primary.p", paste("p was visible to 2 s.f. (0.033) in DEF_SM_CANDIDATES.", WEAK), dp = 4)
add("S6.primary.or_mle", "3.33", L(485), rv("S6.primary.or_mle"), "results.tsv S6.primary.or_mle", FTX, dp = 2)
add("S6.primary.or_lo", "1.003", L(485), rv("S6.primary.or_lo"), "results.tsv S6.primary.or_lo", FTX, dp = 3)
add("S6.primary.or_hi", "12.897", L(485), rv("S6.primary.or_hi"), "results.tsv S6.primary.or_hi", FTX, dp = 3)
add("S6.primary.or_lo (D6)", "1.003", "full text l.576 (§11 D6 'seen: 1.003')", rv("S6.primary.or_lo"), "results.tsv S6.primary.or_lo", paste("D6 compares at 3 dp: Stage 1 1.003 > 1.000, same side as seen.", DET), dp = 3)
add("S6.b.locus_outcome.p", "0.0798", L(486), rv("S6.b.locus_outcome.p"), "results.tsv", paste("p was visible as 0.080 in DEF_SM_CANDIDATES.", WEAK), dp = 4)
add("S6.b.locus_outcome.or_mle", "3.08", L(486), rv("S6.b.locus_outcome.or_mle"), "results.tsv", FTX, dp = 2)
add("S6.b.locus_outcome.or_lo", "0.863", L(486), rv("S6.b.locus_outcome.or_lo"), "results.tsv", FTX, dp = 3)
add("S6.b.locus_outcome.or_hi", "12.48", L(486), rv("S6.b.locus_outcome.or_hi"), "results.tsv", FTX, dp = 2)
add("S6.exclGUSB.*", "p 0.0554, OR 3.03 (0.888-11.91)", L(487), "NOT COMPUTED", "-", "Not a declared §6 analysis ('Excluding GUSB' is a scoping reviewer variant); Stage 1 computes (a), which reclassifies GUSB instead of excluding it.")
add("S6.a.p", "0.0573", L(488), rv("S6.a.p"), "results.tsv S6.a.p", DET, dp = 4)
add("S6.a.or_mle", "2.98", L(488), rv("S6.a.or_mle"), "results.tsv", FTX, dp = 2)
add("S6.a.or_lo", "0.874", L(488), rv("S6.a.or_lo"), "results.tsv", FTX, dp = 3)
add("S6.a.or_hi", "11.70", L(488), rv("S6.a.or_hi"), "results.tsv", FTX, dp = 2)
add("S6.a.locus_outcome", "8/13 vs 45/121, p 0.134", L(488), "NOT COMPUTED", "-",
    "Stage 1 read §6(b) 'for pairs, for locus-outcome groups and for (a)' as three tables, with (a) at the pair level (the only (a) CI the wording rule uses). The full text's Seen block also records a locus-outcome table under (a); if §6(b) intends it, it is missing from Stage 1. Not computed after reading the full text (compute-first rule); flagged for the verifier/user. It does not enter any branch.")
add("S6.a.susie_ratio.APOBR~HF", "5.45", L(489), rv("S6.a.susie_ratio.APOBR~HF"), "results.tsv", "LOG_30 values (7 s.f.) are the same input for both; deterministic.", dp = 2)
add("S6.a.susie_ratio.NT5C3A~eGFR", "5.78", L(489), rv("S6.a.susie_ratio.NT5C3A~eGFR"), "results.tsv", "as above", dp = 2)
add("S6.a.susie_ratio.GUSB~eGFR", "4.80", L(489), rv("S6.a.susie_ratio.GUSB~eGFR"), "results.tsv", "as above", dp = 2)
fm <- fread(file.path(OUT, "w6a_drawn16_failure_modes.tsv"))
add("S6.fail.GUSB_coloc_ratio", "4.9956", L(489), sprintf("%.10g", fm[gene == "GUSB", ratio]), "w6a_drawn16_failure_modes.tsv (TAB_coloc_final full precision)", "deterministic from TAB_coloc_final", dp = 4)
add("S6.pool.clusters", "113", L(490), rv("S6.pool.clusters"), "results.tsv", WEAK, dp = 0)
add("S6.pool.mixed", "14", L(490), rv("S6.pool.mixed"), "results.tsv", WEAK, dp = 0)
add("S6.G1.p", "0.0195", L(491), rv("S6.G1.p"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", DET), dp = 4)
add("S6.G1.or_mle", "4.91", L(491), rv("S6.G1.or_mle"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 2)
add("S6.G1.or_lo", "1.20", L(491), rv("S6.G1.or_lo"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 2)
add("S6.G1.or_hi", "29.1", L(491), rv("S6.G1.or_hi"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 1)
add("S6.G3.p", "0.046", L(491), rv("S6.G3.p"), "results.tsv", "as above", dp = 3)
add("S6.G2.p", "0.234", L(491), rv("S6.G2.p"), "results.tsv", "as above", dp = 3)
add("S6.G4.table", "10/13 vs 43/100", L(491), sprintf("%s/%d vs %s/%d", rv("S6.G4.a"), as.integer(rv("S6.G4.a")) + as.integer(rv("S6.G4.b")), rv("S6.G4.c"), as.integer(rv("S6.G4.c")) + as.integer(rv("S6.G4.d"))), "results.tsv S6.G4.a-d", "identical counts; G4 in Stage 1 = G1 under the (a) rule, which reclassifies exactly GUSB-eGFR (asserted)")
add("S6.G4.p", "0.0357", L(491), rv("S6.G4.p"), "results.tsv", "as above", dp = 4)
add("S6.G4.or_mle", "4.36", L(491), rv("S6.G4.or_mle"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 2)
add("S6.G4.or_lo", "1.04", L(491), rv("S6.G4.or_lo"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 2)
add("S6.G4.or_hi", "26.16", L(491), rv("S6.G4.or_hi"), "results.tsv", paste("SUPPLEMENT-ONLY (§12).", FTX), dp = 2)
add("S6.c.haldane_or", "3.18", L(492), rv("S6.c.haldane_or"), "results.tsv", "asserted exactly (3.1818) by the pre-declaration; weak evidence", dp = 2)

# ---- bootstrap: Monte Carlo quantities; explain generator / indexing differences
BOOTX <- paste("Monte Carlo quantity: bounds carry MC error and are not asserted (§6(c)).",
               "Scoping numpy runs used PCG64 (different generator): agreement expected only within MC error.")
add("S6.c.seed20261007.lo (vs numpy seed 20261007)", "1.01", L(493), rv("S6.c.seed20261007.lo"), "results.tsv", BOOTX, dp = 2)
add("S6.c.seed20261007.hi (vs numpy seed 20261007)", "12.3", L(493), rv("S6.c.seed20261007.hi"), "results.tsv", BOOTX, dp = 1)
add("S6.c.seed20261007.share_le1 (vs numpy)", "0.024", L(493), rv("S6.c.seed20261007.share_le1"), "results.tsv", BOOTX, dp = 3)

# Diagnostic: the scoping R script (seen/prespec_review3/renal_mc.R, opened at the compare phase and logged) used the SAME
# generator and seeds (Mersenne-Twister, 20261007+), but indexed clusters by the Python chain order (chromosome sorted
# as a STRING: 1,10,11,...) and took quantiles on the OR scale.  Stage 1 indexes clusters by numeric chromosome order
# and takes quantiles on the log-OR scale.  Rerun seed 20261007 both ways with the scoping script's draw code.
gs <- fread(SP$SEEN_wpC_pool160_genomic_clusters_csv)
gs[, `:=`(kid = kid %in% c(TRUE, "True"), cand = cand %in% c(TRUE, "True"))]
Kseen <- gs[, .(a = sum(cand & kid), b = sum(cand & !kid), c = sum(!cand & kid), d = sum(!cand & !kid)), by = locus][order(locus), as.matrix(.SD), .SDcols = c("a", "b", "c", "d")]
gm <- fread(file.path(OUT, "w6a_pool160_genomic_clusters.tsv"))
Kmine <- gm[, .(a = sum(cand & kid), b = sum(cand & !kid), c = sum(!cand & kid), d = sum(!cand & !kid)), by = locus][order(locus), as.matrix(.SD), .SDcols = c("a", "b", "c", "d")]
scoping_draws <- function(K, seed) { set.seed(seed); replicate(20000, { i <- sample.int(nrow(K), nrow(K), TRUE); x <- colSums(K[i, , drop = FALSE]); (x[1] + .5) * (x[4] + .5) / ((x[2] + .5) * (x[3] + .5)) }) }
bs_seenorder <- scoping_draws(Kseen, 20261007); bs_mineorder <- scoping_draws(Kmine, 20261007)
d_seen <- sprintf("%.4f-%.4f, %.4f", quantile(bs_seenorder, .025), quantile(bs_seenorder, .975), mean(bs_seenorder <= 1))
d_mine <- sprintf("%.4f-%.4f, %.4f", quantile(bs_mineorder, .025), quantile(bs_mineorder, .975), mean(bs_mineorder <= 1))
mylog <- readRDS(file.path(OUT, "boot_ckpt", "seed_20261007.rds"))
lg("DIAGNOSTIC (compare phase; not a Stage 1 result): scoping draw code, seed 20261007, scoping cluster order -> ", d_seen,
   " ; same code with the Stage 1 cluster order -> ", d_mine, " ; Stage 1 draws (log scale) identical OR values as draws with Stage 1 order: ",
   isTRUE(all.equal(exp(mylog), unname(bs_mineorder))))
add("S6.c.seed20261007 (vs scoping R, same seed)", "1.04-12.30, 2.14%", L(494),
    sprintf("%.4f-%.4f, %.5f", rn("S6.c.seed20261007.lo"), rn("S6.c.seed20261007.hi"), rn("S6.c.seed20261007.share_le1")), "results.tsv S6.c.seed20261007.{lo,hi,share_le1}",
    paste0("DIFFERS beyond printed precision, EXPLAINED: same RNG and seed (Mersenne-Twister, 20261007) but a different cluster INDEX ORDER, so the same random integers pick different clusters; agreement was expected only within MC error, not exactly. ",
           "Diagnostic (compare phase, not a result): the scoping draw code on the scoping cluster order gives ", d_seen, " (reproduces the seen 1.04-12.30, 2.14%); the same code on the Stage 1 order gives ", d_mine,
           "; the Stage 1 draws equal the latter draw-for-draw = ", isTRUE(all.equal(exp(mylog), unname(bs_mineorder))), ". Residual bound differences (<0.001) are log- vs OR-scale quantile interpolation."))
sl <- vapply(20261007:20261016, function(s) rn(sprintf("S6.c.seed%d.lo", s)), 0); sh <- vapply(20261007:20261016, function(s) rn(sprintf("S6.c.seed%d.hi", s)), 0)
ss <- vapply(20261007:20261016, function(s) rn(sprintf("S6.c.seed%d.share_le1", s)), 0)
add("S6.c.seedrange (vs scoping R, 10 seeds)", "lower 1.014-1.052; upper 12.20-12.81; share 2.07-2.37%", L(495),
    sprintf("lower %.3f-%.3f; upper %.2f-%.2f; share %.2f-%.2f%%", min(sl), max(sl), min(sh), max(sh), 100 * min(ss), 100 * max(ss)), "results.tsv S6.c.seed<k>.*",
    "Same generator and seeds, different cluster index order (see the seed-20261007 row): per-seed values are different MC realisations, so the ranges overlap but are not identical. Expected; not a defect.")
add("S6.c.pooled_lo (vs numpy 12 seeds)", "lower bound 1.000-1.036", L(496), rv("S6.c.pooled_lo"), "results.tsv S6.c.pooled_lo (200,000 pooled draws)",
    "Stage 1 pooled lower bound lies inside the scoping numpy per-seed range; different generator; MC quantity, not asserted.")
add("S6.c.pooled_ci", "not seen (Stage 1 pooled estimate is new)", "-", sprintf("%.4f-%.4f", rn("S6.c.pooled_lo"), rn("S6.c.pooled_hi")), "results.tsv", "first computed in Stage 1; no seen counterpart")

# ---- fits and failure modes
fl <- fread(file.path(OUT, "w6a_fit_ledger.tsv"))
add("S6.fitledger.n_rescued_naive", "26 (14 eGFR / 12 non-eGFR)", L(497), sprintf("%s (%s)", rv("S6.fitledger.n_rescued_naive"), R[quantity_id == "S6.fitledger.n_rescued_naive", notes]), "results.tsv", "deterministic log parse; also equals the logs' own tallies (asserted)")
n10_either <- sum(fl$cs_exp == 10 | fl$cs_out == 10); n10_exp <- sum(fl$cs_exp == 10); n10_out <- sum(fl$cs_out == 10)
add("S6.fitledger.n_10_signals (compare-phase count from w6a_fit_ledger.tsv)", "99 of 160", L(497), sprintf("either trait 10 CS: %d; protein 10 CS: %d; outcome 10 CS: %d", n10_either, n10_exp, n10_out),
    "w6a_fit_ledger.tsv cs_exp/cs_out (LOG_29), counted at the compare phase", paste0("Not a declared Stage 1 quantity; counted from the existing Stage 1 ledger only to compare. The seen 99 equals the count under: ", paste(c(if (n10_either == 99) "either trait", if (n10_exp == 99) "protein (exposure) fit with 10 credible sets (= L cap)", if (n10_out == 99) "outcome"), collapse = ", "), if (!any(c(n10_either, n10_exp, n10_out) == 99)) "no definition tried (unexplained)" else "", "."))
add("S6.fitledger.core_10_signals", "51 of 81 core fits", L(497), "NOT COMPUTED", "-", "The 81 core fits are outside W6a's §6 scope (no declared Stage 1 quantity).")
add("S6.fail.n9/n6/n1", "9 / 6 / 1", L(498), sprintf("%s / %s / %s", rv("S6.fail.n9"), rv("S6.fail.n6"), rv("S6.fail.n1")), "results.tsv", "the 9/6/1 counts were in the redacted copy (assertion target): weak evidence. Class membership compared below against the seen file.")
add("S6.fail.kidney_in_9 / kidney_in_other", "7 of the 9; 4 of the other 7", L(498), sprintf("%s; %s", rv("S6.fail.kidney_in_9"), rv("S6.fail.kidney_in_other")), "results.tsv", "deterministic")

# ---- seen files: pool, clusters, Fisher tables, failure modes, fit ledger, BH threshold
ps <- fread(SP$SEEN_wpC_pool160_reconstructed_csv)
pm <- fread(file.path(OUT, "w6a_pool177_rebuild.tsv"))[mhc == FALSE]
kk <- function(d, o = "outcome_name") paste(d$gene, d[[o]], sep = "~")
add("S6.pool.rows (pair set vs seen pool160_reconstructed.csv)", sprintf("%d rows", nrow(ps)), "seen/wpC/pool160_reconstructed.csv", sprintf("%d rows; in both %d; only seen %d; only Stage 1 %d", nrow(pm), length(intersect(kk(ps), kk(pm))), length(setdiff(kk(ps), kk(pm))), length(setdiff(kk(pm), kk(ps)))),
    "w6a_pool177_rebuild.tsv (mhc==FALSE)", if (setequal(kk(ps), kk(pm))) "identical pair set" else "pair sets differ: see counts")
mm <- merge(ps, pm, by = c("gene", "outcome_name"), suffixes = c(".seen", ".s1"))
for (cc in c("b", "se", "se2", "p2", "PP.H3", "PP.H4", "start", "end")) {
  dd <- abs(mm[[paste0(cc, ".seen")]] - mm[[paste0(cc, ".s1")]]); rel <- dd / pmax(abs(mm[[paste0(cc, ".seen")]]), 1e-300)
  add(sprintf("S6.pool.col.%s (160 common rows)", cc), "seen/wpC/pool160_reconstructed.csv", "seen file, full precision", sprintf("max abs diff %.3g; max rel diff %.3g; n=%d", max(dd), max(rel), nrow(mm)), "w6a_pool177_rebuild.tsv",
      if (max(rel) < 1e-9) "identical to floating-point precision" else "DIFFERS: unexplained")
}
for (cc in c("cand", "kid", "mhc")) {
  a <- mm[[paste0(cc, ".seen")]] %in% c(TRUE, "True"); b <- mm[[paste0(cc, ".s1")]] %in% c(TRUE, "True")
  add(sprintf("S6.pool.col.%s (160 common rows)", cc), sprintf("%d TRUE", sum(a)), "seen/wpC/pool160_reconstructed.csv", sprintf("%d TRUE; %d mismatching rows", sum(b), sum(a != b)), "w6a_pool177_rebuild.tsv", if (all(a == b)) "identical" else "DIFFERS")
}
# cluster partition: identical iff the seen-locus x Stage1-locus cross-table is one-to-one
cm <- merge(gs[, .(gene, outcome_name, locus_seen = locus, kid_seen = kid, cand_seen = cand)], gm[, .(gene, outcome_name, locus_s1 = locus, kid, cand)], by = c("gene", "outcome_name"))
one2one <- nrow(cm) == 160 && cm[, uniqueN(locus_s1), by = locus_seen][, all(V1 == 1)] && cm[, uniqueN(locus_seen), by = locus_s1][, all(V1 == 1)]
add("S6.pool.clusters (partition vs seen pool160_genomic_clusters.csv)", sprintf("%d clusters", uniqueN(gs$locus)), "seen/wpC/pool160_genomic_clusters.csv",
    sprintf("%d clusters; partition identical (one-to-one) = %s; kid/cand flags identical = %s", uniqueN(gm$locus), one2one, all(cm$kid_seen == cm$kid & cm$cand_seen == cm$cand)), "w6a_pool160_genomic_clusters.tsv",
    "Locus NUMBERS differ by design (scoping sorted chromosome as a string; Stage 1 numerically); the partition is what the estimator uses.")
fr <- fread(SP$SEEN_wpC_renal_fisher_R_csv, colClasses = "character")
dp_of <- function(x) ifelse(grepl(".", x, fixed = TRUE), nchar(sub("^[^.]*[.]", "", x)), 0L)
fmap <- list(A = "S6.primary", C = "S6.a", D = "S6.b.locus_outcome", G1 = "S6.G1", G2 = "S6.G2", G3 = "S6.G3", G4 = "S6.G4")
for (lab in names(fmap)) {
  row <- fr[grepl(paste0("^", lab, " "), analysis)]; stopifnot(nrow(row) == 1); pre <- fmap[[lab]]
  add(paste0(pre, ".p (vs renal_fisher_R.csv ", lab, ")"), row$p, "seen/wpC/renal_fisher_R.csv (printed 5 s.f.)", rv(paste0(pre, ".p")), "results.tsv", "R fisher.test both; deterministic", dp = dp_of(row$p))
  add(paste0(pre, ".or_mle (vs renal_fisher_R.csv ", lab, ")"), row$cond_MLE_OR, "seen/wpC/renal_fisher_R.csv (3 dp)", rv(paste0(pre, ".or_mle")), "results.tsv", FTX, dp = dp_of(row$cond_MLE_OR))
  add(paste0(pre, ".or_lo (vs renal_fisher_R.csv ", lab, ")"), row$lo95, "seen/wpC/renal_fisher_R.csv (3 dp)", rv(paste0(pre, ".or_lo")), "results.tsv", FTX, dp = dp_of(row$lo95))
  add(paste0(pre, ".or_hi (vs renal_fisher_R.csv ", lab, ")"), row$hi95, "seen/wpC/renal_fisher_R.csv (3 dp)", rv(paste0(pre, ".or_hi")), "results.tsv", FTX, dp = dp_of(row$hi95))
}
for (lab in c("B", "D2")) add(sprintf("renal_fisher_R.csv %s", lab), fr[grepl(paste0("^", lab, " "), analysis), sprintf("%s: p %s OR %s (%s-%s)", analysis, p, cond_MLE_OR, lo95, hi95)], "seen/wpC/renal_fisher_R.csv", "NOT COMPUTED", "-",
                              if (lab == "B") "GUSB-excluded variant: not declared in §6." else "Locus-outcome under (a): see row S6.a.locus_outcome (interpretation of §6(b) flagged).")
fp <- fread(SP$SEEN_wpC_renal_fisher_python_csv)
pmap <- list("A pair" = "S6.primary", "C GUSB" = "S6.a", "D locus" = "S6.b.locus_outcome", "G1 " = "S6.G1", "G2 " = "S6.G2", "G3 " = "S6.G3", "G4 " = "S6.G4")
for (lab in names(pmap)) {
  row <- fp[startsWith(analysis, lab)]; stopifnot(nrow(row) == 1); pre <- pmap[[lab]]
  s1 <- sprintf("%s/%s/%s/%s p=%s", rv(paste0(pre, ".a")), rv(paste0(pre, ".b")), rv(paste0(pre, ".c")), rv(paste0(pre, ".d")), rv(paste0(pre, ".p")))
  sv <- sprintf("%d/%d/%d/%d p=%.10g", row$cand_kid, row$cand_non, row$other_kid, row$other_non, row$fisher_p)
  same <- all(c(row$cand_kid, row$cand_non, row$other_kid, row$other_non) == as.integer(c(rv(paste0(pre, ".a")), rv(paste0(pre, ".b")), rv(paste0(pre, ".c")), rv(paste0(pre, ".d"))))) &&
    abs(row$fisher_p - rn(paste0(pre, ".p"))) < 1e-9
  add(paste0(pre, ".table+p (vs renal_fisher_python.csv)"), sv, "seen/wpC/renal_fisher_python.csv (scipy)", s1, "results.tsv", if (same) "identical counts; p identical to 1e-9 (scipy vs R two-sided exact)" else "DIFFERS: unexplained")
}
meta <- fromJSON(SP$SEEN_wpC_pool_meta_json)
add("S6.pool.bh_threshold", sprintf("%.16g", meta$thr_bh_second_order), "seen/wpC/pool_meta.json", rv("S6.pool.bh_threshold"), "results.tsv (printed to 10 s.f.)", "deterministic", dp = NA)
d16 <- fread(SP$SEEN_wpC_drawn16_failure_modes_csv)
d16[, cls := fifelse(startsWith(abf_failure_mode, "ratio"), "ratio", fifelse(startsWith(abf_failure_mode, "H3"), "H3", "H4"))]
fm[, cls := fifelse(startsWith(mode, "ratio"), "ratio", fifelse(startsWith(mode, "H3"), "H3", "H4"))]
fc <- merge(d16[, .(gene, outcome_name = outcome, cls_seen = cls, sym_seen = susie_passes_symmetric %in% c(TRUE, "True"))], fm[, .(gene, outcome_name, cls)], by = c("gene", "outcome_name"))
sus_s1 <- fl[candidate == TRUE, .(gene, outcome_name, sym = susie_sym_pass %in% c(TRUE, "TRUE"))]
fc <- merge(fc, sus_s1, by = c("gene", "outcome_name"))
add("S6.fail.class_membership (16 pairs vs drawn16_failure_modes.csv)", sprintf("%d rows", nrow(d16)), "seen/wpC/drawn16_failure_modes.csv", sprintf("%d matched; class identical %d/%d; SuSiE symmetric-rule flag identical %d/%d", nrow(fc), sum(fc$cls_seen == fc$cls), nrow(fc), sum(fc$sym_seen == fc$sym), nrow(fc)),
    "w6a_drawn16_failure_modes.tsv; w6a_fit_ledger.tsv", if (nrow(fc) == 16 && all(fc$cls_seen == fc$cls) && all(fc$sym_seen == fc$sym)) "identical" else "DIFFERS")
fs <- fread(SP$SEEN_wpC_susie_candidate_fit_ledger_from_logs_csv)
setnames(fs, setdiff(names(fs), c("gene", "outcome")), paste0("seen_", setdiff(names(fs), c("gene", "outcome"))))
fj <- merge(fs, fl, by.x = c("gene", "outcome"), by.y = c("gene", "outcome_name"))
cmpf <- list(window = list(fs = "seen_region", s1 = function(d) sprintf("chr%s:%d-%d", d$window_chr, as.integer(d$window_start), as.integer(d$window_end))),
             merged = list(fs = "seen_merged", s1 = "n_merged"), ld_aligned = list(fs = "seen_ld_aligned", s1 = "n_ld_aligned"), ld_matrix = list(fs = "seen_ld_matrix", s1 = "ld_matrix_dim"),
             outcome_N = list(fs = "seen_outcome_N", s1 = "outcome_N_regional"), cs_exp = list(fs = "seen_cs_exp", s1 = "cs_exp"), cs_out = list(fs = "seen_cs_out", s1 = "cs_out"),
             n_pairs = list(fs = "seen_n_pairs", s1 = "n_signal_pairs"), susie_h4 = list(fs = "seen_susie_h4", s1 = "susie_h4_log3"), susie_h3 = list(fs = "seen_susie_h3", s1 = "susie_h3_log3"))
for (nm in names(cmpf)) {
  a <- fj[[cmpf[[nm]]$fs]]; b <- if (is.function(cmpf[[nm]]$s1)) cmpf[[nm]]$s1(fj) else fj[[cmpf[[nm]]$s1]]
  eq <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b)) | (!is.na(suppressWarnings(as.numeric(a))) & !is.na(suppressWarnings(as.numeric(b))) & abs(suppressWarnings(as.numeric(a)) - suppressWarnings(as.numeric(b))) < 1e-9)
  add(sprintf("S6.fitledger.%s (160 fits vs susie_candidate_fit_ledger_from_logs.csv)", nm), sprintf("%d rows", nrow(fs)), "seen/wpC/susie_candidate_fit_ledger_from_logs.csv",
      sprintf("matched %d; identical %d; differing %d%s", nrow(fj), sum(eq), sum(!eq), if (any(!eq)) paste0(" (e.g. ", paste(head(sprintf("%s~%s seen=%s s1=%s", fj$gene[!eq], fj$outcome[!eq], a[!eq], b[!eq]), 3), collapse = "; "), ")") else ""),
      "w6a_fit_ledger.tsv", if (all(eq)) "identical (log parse)" else "DIFFERS: see examples; field-definition difference or unexplained")
}
add("CLASS_FLIP", "n/a", "-", "n/a", "-", "W6a computes no SuSiE s; the §6 scoping s classification (UMOD s 0.29) belongs to W6b.")
add("RNG note (rule 8)", "scoping: numpy PCG64 seed 20261007 (renal_enrichment.py; 12-seed numpy runs) and R Mersenne-Twister seeds 20261007-16 (renal_mc.R)", "seen scripts",
    "Stage 1: R Mersenne-Twister/Rejection, seeds 20261007-20261016, 20,000 draws each, pooled", "w6a_compute.R",
    "Same RNG and seed as the scoping R run, but agreement is NOT expected draw-for-draw: the cluster index order differs (string vs numeric chromosome sort). Diagnostic above reproduces the seen R value exactly under the scoping order.")

out <- rbindlist(CMP)
fwrite(out, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = FALSE)
lg("WROTE COMPARE_SEEN.tsv rows=", nrow(out), " md5=", md5f(file.path(OUT, "COMPARE_SEEN.tsv")), " | rows flagged DIFFERS: ", sum(grepl("DIFFERS", out$explanation)),
   " | results.tsv md5 still ", md5f(file.path(OUT, "results.tsv")))
