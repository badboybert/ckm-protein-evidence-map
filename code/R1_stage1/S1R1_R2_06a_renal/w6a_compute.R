# W6a -- pre-declaration §6 renal: pool rebuild, primary Fisher test, (a), (b), (c), G1-G4,
# 9/6/1 failure modes, SuSiE fit ledger.  Stage 1 of Paper 3.1 R2.
# Computes from the REDACTED pre-declaration only.  Writes ONLY into this folder.
# Run:  Rscript w6a_run_all.R   (this file is step 1; w6a_compare.R is step 2, run after results.tsv exists)
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
suppressMessages({ library(data.table); library(digest); library(jsonlite) })

OUT <- Sys.getenv("W6A_OUT", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))   # override to re-run into another folder (e.g. verify_L3/rerun)
RESDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
STAGING <- paste0(.P31R_V3, "/V4_local/staging/l3-finemap")
stopifnot(dir.exists(OUT))
LOGF <- file.path(OUT, "run.log")
T0 <- Sys.time()

lg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] ", paste0(..., collapse = ""))
  cat(line, "\n", sep = "", file = LOGF, append = TRUE); cat(line, "\n", sep = ""); flush.console()
}
md5f <- function(p) unname(as.character(tools::md5sum(p)))

cat("
==================== RUN ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " (w6a_compute.R) ====================
", sep = "", file = LOGF, append = TRUE)   # append-only: the log keeps every run, so timestamps prove ordering
lg("=== W6a compute start === R ", R.version.string, " | data.table ", as.character(packageVersion("data.table")),
   " | digest ", as.character(packageVersion("digest")), " | jsonlite ", as.character(packageVersion("jsonlite")))
lg("RNG for (c): ", paste(RNGkind(), collapse = "/"), " ; seeds 20261007..20261016 (pre-declaration §6(c)); set.seed(seed, kind='Mersenne-Twister', normal.kind='Inversion', sample.kind='Rejection')")

# ---------------------------------------------------------------- 1. governing md5s (STOP on mismatch)
GOV <- list(
  BRIEFS          = list(p = file.path(RESDIR, "STAGE1_BRIEFS.md"), m = "40d0d52f7055542d774c03f1edf04e00"),
  PRESPEC_RED     = list(p = file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"), m = "a281c3464619b1b67c28e13e7392c41c"),
  AMEND1_RED      = list(p = file.path(RESDIR, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"), m = "f265d600973c43239bb8377e0cd39d33"),
  PRESPEC_FULL    = list(p = paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), m = "3319dc39740d9f5b3b8822133d9108b7"),
  AMEND1_FULL     = list(p = paste0(.P31R_V3, "/R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), m = "377508832ea9bd384133e869c456c909"),
  GLS_Q_R         = list(p = file.path(RESDIR, "gls_q.R"), m = "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  RESOURCES_JSON  = list(p = file.path(RESDIR, "RESOURCES.json"), m = "6d21238155a62fb58f71e05a01638dd0"))
for (nm in names(GOV)) {
  got <- md5f(GOV[[nm]]$p)
  lg("GOVERNING md5 ", nm, " = ", got, " (expected ", GOV[[nm]]$m, ") ", GOV[[nm]]$p)
  if (!identical(got, GOV[[nm]]$m)) stop("GOVERNING md5 MISMATCH for ", nm, " -- stopping (rule 1)")
}
lg("NOTE: the full pre-declaration and full amendment are md5-asserted only; their CONTENT is not read by this script (compute-first rule 10).")

# gls_q.R: not used by W6a (no GLS in §6 renal); md5 asserted above; self-test run for the log (§0 rule 6)
st <- tryCatch(system2(file.path(R.home("bin"), "Rscript.exe"), c(shQuote(GOV$GLS_Q_R$p), "--self-test"), stdout = TRUE, stderr = TRUE),
               error = function(e) paste("self-test failed to launch:", conditionMessage(e)))
for (s in st) lg("gls_q.R --self-test | ", s)

# ---------------------------------------------------------------- 2. resource map
RES <- fromJSON(GOV$RESOURCES_JSON$p, simplifyVector = FALSE)
rmap <- setNames(RES, vapply(RES, function(e) e$id, ""))
INPUT_IDS <- c("TAB_mr_all_results", "TAB_instruments", "TAB_coloc_final", "TAB_edge_audit",
               "LOG_28_candidate_pool", "LOG_29_cand_egfr", "LOG_29_cand_nonegfr", "LOG_30_candidate_audit",
               "LOG_31_candidate_egfr", "LOG_31_candidate_nonegfr", "DEF_SM_CANDIDATES", "GRID_ens116_genes",
               "SEEN_wpC_renal_enrichment_py")
P <- list()
for (id in INPUT_IDS) {
  e <- rmap[[id]]; if (is.null(e)) stop("resource id missing from RESOURCES.json: ", id)
  got <- md5f(e$path)
  lg("INPUT ", id, " role=", e$role, " md5=", got, " (RESOURCES ", e$md5, ") bytes=", file.size(e$path), " ", e$path)
  if (!identical(got, e$md5)) stop("INPUT md5 MISMATCH for ", id, " -- stopping (rule 1)")
  P[[id]] <- e$path
}
lg("SEEN SCRIPT OPENED (method reading, before any computation): SEEN_wpC_renal_enrichment_py ",
   "(seen/wpC/renal_enrichment.py, md5 762c4a9817a515668452e0fb97f9f69c) at 2026-10-08T10:31:31+08:00. ",
   "Reason: the W6a brief instructs reading it to find the estimator's inputs (the V4 source tables, LOG_28, LOG_30) ",
   "and the pinned §6(c) estimator; it also supplies the G1-G4 definitions, which the redacted copy (§12) does not carry. ",
   "No other seen script and no seen output was opened before results.tsv was written.")

# ---------------------------------------------------------------- 3. staging snapshot (start)
snap <- function(tag) {
  f <- list.files(STAGING, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  d <- data.table(file = sub(paste0("^", STAGING, "/"), "", f), bytes = file.size(f), md5 = vapply(f, md5f, ""))
  setorder(d, file)
  fwrite(d, file.path(OUT, sprintf("staging_snapshot_%s.tsv", tag)), sep = "\t")
  agg <- digest(paste(d$file, d$md5, collapse = "\n"), algo = "md5", serialize = FALSE)
  lg("STAGING SNAPSHOT ", tag, ": ", nrow(d), " files under V4_local/staging/l3-finemap, aggregate md5(file+md5 list) = ", agg)
  agg
}
SNAP_START <- snap("start")

# ---------------------------------------------------------------- 4. result / assert machinery
ROWS <- list(); NASSERT_FAIL <- 0L; ASSERT_LOG <- list()
fmt <- function(x) if (is.null(x) || length(x) == 0) "" else if (is.numeric(x)) ifelse(is.na(x), "NA", sprintf("%.10g", x)) else as.character(x)
addrow <- function(qid, ref, estimand, value, resource = "", panel = "none (no LD in W6a)", cov = "NA", se_order = "NA",
                   se = NA, ci_lo = NA, ci_hi = NA, p = NA, n = NA, unit = "", notes = "") {
  ROWS[[length(ROWS) + 1]] <<- data.table(quantity_id = qid, prespec_ref = ref, estimand = estimand, resource = resource,
    panel = panel, cov = cov, se_order = se_order, value = fmt(value), se = fmt(se), ci_lo = fmt(ci_lo), ci_hi = fmt(ci_hi),
    p = fmt(p), n = fmt(n), unit = unit, notes = notes)
}
chk <- function(id, ok, got, target, derivation) {
  ok <- isTRUE(ok)
  lg(if (ok) "ASSERT OK   " else "ASSERT FAIL ", id, " | got=", paste(got, collapse = ","), " target=", paste(target, collapse = ","), " | ", derivation)
  ASSERT_LOG[[length(ASSERT_LOG) + 1]] <<- data.table(assert_id = id, ok = ok, got = paste(got, collapse = ","),
                                                     target = paste(target, collapse = ","), derivation = derivation)
  if (!ok) {
    NASSERT_FAIL <<- NASSERT_FAIL + 1L
    addrow(paste0("ASSERT_FAILED.", id), "§0 rule 9", "assertion target", paste(got, collapse = ","),
           notes = paste0("ASSERT_FAILED: target ", paste(target, collapse = ","), "; derivation: ", derivation))
  }
  invisible(ok)
}
rdlines <- function(p) readLines(p, encoding = "UTF-8", warn = FALSE)
# FIXER F1 (FIXLOG.md): conditional-MLE OR and exact CI from fisher.test's own equations, roots solved at a tight tolerance
# (fisher.test's default uniroot tolerance moved the upper bounds in the 3rd decimal). Defined before the mutation injection point.
ECOR <- file.path(Sys.getenv("W6A_SRC", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal")), "exact_cond_or.R")
ECOR_MD5 <- "1217326d2a6880b174f91d377c5b21b9"
lg("SHARED SCRIPT exact_cond_or.R md5 ", md5f(ECOR), " (expected ", ECOR_MD5, ") ", ECOR)
if (!identical(md5f(ECOR), ECOR_MD5)) stop("exact_cond_or.R md5 MISMATCH -- stopping (rule 1)")
source(ECOR)
EXACT_TOL <- 1e-13; EXACT_CONF <- 0.95

# ---------------------------------------------------------------- 5. read inputs (assert row counts: fread + spaces)
nlines_data <- function(p) length(rdlines(p)) - 1L
MR  <- fread(P$TAB_mr_all_results, na.strings = c("NA", ""))
INS <- fread(P$TAB_instruments, na.strings = c("NA", ""))
CO  <- fread(P$TAB_coloc_final, na.strings = c("NA", ""))
CA  <- fread(P$TAB_edge_audit, na.strings = c("NA", ""), encoding = "UTF-8")
GE  <- fread(P$GRID_ens116_genes, na.strings = c("NA", ""))
for (nm in c("MR", "INS", "CO", "CA", "GE")) {
  pth <- switch(nm, MR = P$TAB_mr_all_results, INS = P$TAB_instruments, CO = P$TAB_coloc_final, CA = P$TAB_edge_audit, GE = P$GRID_ens116_genes)
  chk(paste0("read.", nm, ".rows"), nrow(get(nm)) == nlines_data(pth) && nrow(get(nm)) > 0, nrow(get(nm)), nlines_data(pth),
      "fread row count equals file line count minus header (guards the fread-on-a-path-with-spaces failure)")
}
# placeholder tokens / numeric parsing of every column used as a number
num_ok <- function(d, cols) all(vapply(cols, function(cc) is.numeric(d[[cc]]), TRUE))
chk("parse.MR.numeric", num_ok(MR, c("b", "se", "pval", "n_snps")), "class check", "numeric", "MR b/se/pval/n_snps parsed as numeric (no placeholder tokens)")
chk("parse.INS.numeric", num_ok(INS, c("beta", "se", "F_stat", "gene_start", "gene_end")), "class check", "numeric", "INS beta/se/F/gene coords numeric")
chk("parse.CO.numeric", num_ok(CO, c("PP.H3", "PP.H4")), "class check", "numeric", "coloc PP.H3/PP.H4 numeric")
for (cc in c("gene", "outcome_name", "method", "b", "se")) lg("NULL RATE MR$", cc, " = ", sum(is.na(MR[[cc]])), "/", nrow(MR))
for (cc in c("gene", "beta", "se", "gene_chrom", "gene_start", "gene_end")) lg("NULL RATE INS$", cc, " = ", sum(is.na(INS[[cc]])), "/", nrow(INS))
for (cc in c("gene", "outcome_name", "PP.H3", "PP.H4")) lg("NULL RATE CO$", cc, " = ", sum(is.na(CO[[cc]])), "/", nrow(CO))
L28 <- paste(rdlines(P$LOG_28_candidate_pool), collapse = "\n")
num_from <- function(txt, re) as.integer(regmatches(txt, regexec(re, txt))[[1]][-1])

# ---------------------------------------------------------------- 6. family-wide FDR with second-order SEs (DEF_SM_CANDIDATES + LOG_28)
FAM <- 1948L * 8L   # 15,584 = proteins x outcomes (definition text); derived, then checked against the INS gene count below
chk("family.size", FAM == 15584L && uniqueN(INS$gene) == 1948L && uniqueN(MR$outcome_name) == 8L,
    c(FAM, uniqueN(INS$gene), uniqueN(MR$outcome_name)), c(15584, 1948, 8), "family = unique INS genes x unique MR outcomes")
prim <- MR[method %in% c("Wald ratio", "Inverse variance weighted")]
chk("prim.nodup", !anyDuplicated(prim[, .(gene, outcome_name)]), sum(duplicated(prim[, .(gene, outcome_name)])), 0,
    "one primary test per gene x outcome (Wald ratio for 1 SNP, IVW for >=2)")
nW <- prim[method == "Wald ratio", .N]; nIVW <- prim[method == "Inverse variance weighted", .N]
L28_wald <- num_from(L28, "second-order SE on (\\d+) single-instrument row")
L28_multi <- num_from(L28, "first-order kept for (\\d+) multi-instrument row")
chk("se.wald_rows", nW == L28_wald, nW, L28_wald, "Wald-ratio rows in TAB_mr_all_results vs LOG_28 '#second-order rows'")
chk("se.multi_rows", MR[method != "Wald ratio", .N] == L28_multi, MR[method != "Wald ratio", .N], L28_multi,
    "non-Wald rows (IVW+Egger+WM) vs LOG_28 '#first-order multi-instrument rows'")
ni <- INS[, .N, by = gene]
one <- INS[gene %in% ni[N == 1L, gene]]
prim[, `:=`(bx = NA_real_, sx = NA_real_, bx_alt = NA_real_, sx_alt = NA_real_, inst_case = NA_character_)]
w1 <- prim$method == "Wald ratio" & prim$gene %in% one$gene
prim[w1, `:=`(bx = one$beta[match(gene, one$gene)], sx = one$se[match(gene, one$gene)], inst_case = "single-instrument gene")]
wm <- which(prim$method == "Wald ratio" & !(prim$gene %in% one$gene))
lg("JOIN Wald rows -> INS instrument: ", sum(w1), " via single-instrument genes; ", length(wm),
   " Wald rows belong to genes with >1 INS instrument (", paste(unique(prim$gene[wm]), collapse = ","), "); ",
   "the instrument used is not recorded in TAB_mr_all_results, so the second-order term is bounded by the min and max |bx|/sx over that gene's instruments")
for (i in wm) {
  g <- INS[gene == prim$gene[i]]
  Fz <- (g$beta / g$se)^2
  lo <- which.min(Fz); hi <- which.max(Fz)
  prim[i, `:=`(bx = g$beta[lo], sx = g$se[lo], bx_alt = g$beta[hi], sx_alt = g$se[hi], inst_case = "multi-instrument gene: weakest-F instrument (primary, conservative); strongest-F instrument (alt)")]
}
chk("se.wald_join_yield", sum(!is.na(prim[method == "Wald ratio", bx])) == nW, sum(!is.na(prim[method == "Wald ratio", bx])), nW,
    "every Wald row gets an instrument bx/sx")
# second-order (delta method, first two terms): var2 = se1^2 + b^2 (sx/bx)^2
prim[, se2 := se]
prim[method == "Wald ratio", se2 := sqrt(se^2 + b^2 * (sx / bx)^2)]
prim[, se2_alt := se2]
prim[!is.na(bx_alt), se2_alt := sqrt(se^2 + b^2 * (sx_alt / bx_alt)^2)]
prim[, p2 := 2 * pnorm(-abs(b / se2))]
prim[, p2_alt := 2 * pnorm(-abs(b / se2_alt))]
bh <- function(p) {   # BH over the full family; untestable pairs padded with p = 1 (definition: the family is 15,584 tests)
  pp <- c(p, rep(1, FAM - length(p)))
  q <- p.adjust(pp, "BH")[seq_along(p)]
  q2 <- p.adjust(p, "BH", n = FAM)
  stopifnot(isTRUE(all.equal(q, q2)))
  q
}
chk("bh.padding", FAM - nrow(prim) == 15584L - nrow(prim) && nrow(prim) == nW + nIVW, FAM - nrow(prim), 15584 - nW - nIVW,
    "untestable pairs padded as p = 1 to reach 15,584")
prim[, q2 := bh(p2)]; prim[, q2_alt := bh(p2_alt)]
fdr <- prim[q2 <= 0.05]
fdr_alt <- prim[q2_alt <= 0.05]
thr <- max(fdr$p2)
same_fdr <- setequal(paste(fdr$gene, fdr$outcome_name), paste(fdr_alt$gene, fdr_alt$outcome_name))
chk("fdr.multi_instrument_invariance", same_fdr, nrow(fdr_alt), nrow(fdr),
    "FDR set identical whether the 21 multi-instrument-gene Wald rows take the weakest- or strongest-F instrument's second-order term")
lg("multi-instrument-gene Wald rows: min p2 (primary) = ", sprintf("%.3g", min(prim$p2[wm])), "; any FDR-significant: ",
   any(prim$q2[wm] <= 0.05), " ; BH p threshold (largest p2 passing) = ", sprintf("%.6g", thr))
L28_fdr <- num_from(L28, "MR FDR-significant pairs = (\\d+)")
chk("fdr.rows", nrow(fdr) == L28_fdr, nrow(fdr), L28_fdr, "BH(second-order, family 15,584) <= 0.05 vs LOG_28 'MR FDR-significant pairs'")
# first-order count for context (not a declared quantity)
prim[, q1 := bh(2 * pnorm(-abs(b / se)))]
lg("context: first-order BH count over the same family = ", sum(prim$q1 <= 0.05))

# ---------------------------------------------------------------- 7. colocalization rule and pool
CO[, ratio := PP.H4 / PP.H3]
CO[, core := PP.H4 >= 0.80 & ratio >= 5]   # PP.H3 = 0 gives ratio = Inf (passes), as R and the definition imply
L28_pass <- num_from(L28, "solid\\(pass\\)=(\\d+)")
L28_pairs <- num_from(L28, "coloc.abf pairs=(\\d+)")
chk("coloc.pairs", nrow(CO) == L28_pairs && !anyDuplicated(CO[, .(gene, outcome_name)]), nrow(CO), L28_pairs, "TAB_coloc_final rows vs LOG_28 coloc.abf pairs")
chk("coloc.core_pass", sum(CO$core) == L28_pass, sum(CO$core), L28_pass, "H4>=0.80 & H4/H3>=5 count vs LOG_28 solid(pass)")
fdr <- merge(fdr, CO[, .(gene, outcome_name, PP.H3, PP.H4, ratio, core, independent_signal)], by = c("gene", "outcome_name"), all.x = TRUE)
chk("join.fdr_coloc_yield", sum(!is.na(fdr$PP.H4)) == nrow(fdr), sum(!is.na(fdr$PP.H4)), nrow(fdr),
    "every FDR-significant pair has a coloc.abf row (LOG_28: candidate_pool_no_coloc.csv empty)")
pool <- fdr[core == FALSE]
L28_pool <- num_from(L28, "candidate pool = (\\d+) coloc-fail pairs \\| testable\\(SuSiE\\)=(\\d+) \\[(\\d+) eGFR, (\\d+) non-eGFR\\] \\| MHC=(\\d+)")
chk("pool.colocfail", nrow(pool) == L28_pool[1], nrow(pool), L28_pool[1], "FDR pairs failing the coloc rule vs LOG_28")
# by-outcome table in LOG_28
l28l <- rdlines(P$LOG_28_candidate_pool)
s <- grep("== coloc-fail candidates by outcome ==", l28l, fixed = TRUE); e <- grep("== top 15 by coloc.abf PP.H4", l28l, fixed = TRUE)
tb <- regmatches(l28l[s:e], regexec("^\\s*\\d+:\\s+(\\S+)\\s+(\\d+)\\s*$", l28l[s:e])); tb <- do.call(rbind, tb[lengths(tb) == 3])
L28_byout <- data.table(outcome_name = tb[, 2], N_log = as.integer(tb[, 3]))
mybyout <- pool[, .N, by = outcome_name]
cmp <- merge(L28_byout, mybyout, by = "outcome_name", all = TRUE)
chk("pool.colocfail_by_outcome", nrow(cmp) == nrow(L28_byout) && all(cmp$N == cmp$N_log, na.rm = FALSE) && !anyNA(cmp),
    paste(cmp$outcome_name, cmp$N, sep = "=", collapse = ";"), paste(cmp$outcome_name, cmp$N_log, sep = "=", collapse = ";"),
    "coloc-fail pairs by outcome vs LOG_28 table")
# top-15 by PP.H4 in LOG_28 (gene, outcome, PP.H4_abf to 7 digits)
tt <- regmatches(l28l[e:length(l28l)], regexec("^\\s*\\d+:\\s+(\\S+)\\s+(\\S+)\\s+([0-9.e+-]+)\\s+([0-9.e+-]+)\\s+(TRUE|FALSE)\\s+(\\S+)\\s*$", l28l[e:length(l28l)]))
tt <- do.call(rbind, tt[lengths(tt) == 7])
my15 <- pool[order(-PP.H4)][1:15]
ok15 <- nrow(tt) == 15 && all(my15$gene == tt[, 2]) && all(my15$outcome_name == tt[, 3]) &&
  all(abs(signif(my15$PP.H4, 7) - as.numeric(tt[, 4])) < 1e-9 * 10)
chk("pool.top15_identity", ok15, paste(my15$gene, collapse = ","), paste(tt[, 2], collapse = ","), "top 15 coloc-fail pairs by PP.H4 (gene, outcome, PP.H4 to log precision) vs LOG_28")

# gene coordinates (INS, one row per gene) and build check
gc0 <- INS[, .(nchr = uniqueN(gene_chrom), nst = uniqueN(gene_start), nen = uniqueN(gene_end)), by = gene]
chk("ins.coords_unique", all(gc0$nchr == 1 & gc0$nst == 1 & gc0$nen == 1), sum(gc0$nchr > 1 | gc0$nst > 1 | gc0$nen > 1), 0, "one coordinate set per gene in INS")
g1 <- unique(INS[, .(gene, chr = as.character(gene_chrom), start = as.numeric(gene_start), end = as.numeric(gene_end))])
bchk <- merge(g1, unique(GE[, .(gene = symbol, start38, end38)]), by = "gene")
lg("BUILD CHECK INS gene coords vs GRID_ens116_genes (GRCh38): ", nrow(bchk), " shared symbols; exact start match ",
   sum(bchk$start == bchk$start38), ", |diff|<=1kb ", sum(abs(bchk$start - bchk$start38) <= 1000), "; examples: ",
   paste(head(sprintf("%s %d/%d", bchk$gene, as.integer(bchk$start), as.integer(bchk$start38)), 6), collapse = "; "))
chk("build.ins_grch38", nrow(bchk) >= 5 && mean(abs(bchk$start - bchk$start38) <= 1000) >= 0.9,
    sprintf("%d/%d within 1 kb", sum(abs(bchk$start - bchk$start38) <= 1000), nrow(bchk)), ">=90% within 1 kb",
    "INS gene_start agrees with Ensembl 116 GRCh38 starts -> MHC window applied on GRCh38 as the definition states")
pool <- merge(pool, g1, by = "gene", all.x = TRUE)
chk("join.pool_coords_yield", sum(!is.na(pool$start)) == nrow(pool), sum(!is.na(pool$start)), nrow(pool), "every coloc-fail pair has gene coordinates")
pool[, mhc := chr == "6" & start <= 34e6 & end >= 25e6]
chk("pool.mhc", sum(pool$mhc) == L28_pool[5], sum(pool$mhc), L28_pool[5], "genes overlapping chr6 25-34 Mb (GRCh38) vs LOG_28 MHC")
tst <- pool[mhc == FALSE]
KID <- c("eGFR", "CKD")
tst[, kid := outcome_name %in% KID]
chk("pool.rows", nrow(tst) == L28_pool[2] && nrow(tst) == 160L, nrow(tst), c(L28_pool[2], 160), "testable pool vs LOG_28 and the §6(c) assertion target 160")
chk("pool.egfr_split", sum(tst$outcome_name == "eGFR") == L28_pool[3] && sum(tst$outcome_name != "eGFR") == L28_pool[4],
    c(sum(tst$outcome_name == "eGFR"), sum(tst$outcome_name != "eGFR")), L28_pool[3:4], "eGFR / non-eGFR split vs LOG_28")

# the 160 pairs the SuSiE stage actually ran (LOG_29 headers): an independent list from the pipeline record
parse_hdr <- function(p) {
  x <- rdlines(p); m <- regmatches(x, regexec("^\\[(\\d+)/(\\d+)\\] (\\S+) ~ (\\S+) \\(coloc\\.abf PP\\.H4=([0-9.]+)\\)", x))
  m <- do.call(rbind, m[lengths(m) == 6]); data.table(k = as.integer(m[, 2]), N = as.integer(m[, 3]), gene = m[, 4], outcome_name = m[, 5], abf_h4_log = as.numeric(m[, 6]))
}
h29 <- rbind(cbind(mode = "egfr", parse_hdr(P$LOG_29_cand_egfr)), cbind(mode = "nonegfr", parse_hdr(P$LOG_29_cand_nonegfr)))
pk <- function(d) paste(d$gene, d$outcome_name, sep = "~")
chk("pool.vs_log29_list", nrow(h29) == 160L && setequal(pk(h29), pk(tst)) && !anyDuplicated(pk(h29)),
    sprintf("log29=%d, in_both=%d", nrow(h29), length(intersect(pk(h29), pk(tst)))), "160 identical",
    "rebuilt 160 pairs vs the 61+99 pairs LOG_29 ran SuSiE on")
d29 <- merge(h29, tst, by = c("gene", "outcome_name"))[, abs(round(PP.H4, 3) - abf_h4_log)]
chk("pool.vs_log29_abf_h4", length(d29) == 160L && all(d29 < 1e-9), sprintf("n=%d max|diff|=%.2g", length(d29), max(d29)), "n=160, 0 at 3 dp", "coloc.abf PP.H4 printed in LOG_29 headers equals TAB_coloc_final rounded to 3 dp")

# ---------------------------------------------------------------- 8. candidates: two independent routes
drawn <- CA[Status == "Drawn", .(gene = Gene, outcome_name = Endpoint)]
# route 2: the logs -- LOG_30 dashed edges (26), LOG_31 strict survivors, minus fits with a convergence warning
parse_l30 <- function(p) {
  x <- rdlines(p)
  m <- regmatches(x, regexec("^\\s*\\d+:\\s+(\\S+)\\s+(\\S+)\\s+(risk|protective)\\s+([0-9.e+-]+)\\s+([0-9.e+-]+)\\s+([0-9.e+-]+)\\s*$", x))
  m <- do.call(rbind, m[lengths(m) == 7])
  data.table(gene = m[, 2], outcome_name = m[, 3], mr_direction = m[, 4], h4_abf_l30 = as.numeric(m[, 5]), h4_susie_l30 = as.numeric(m[, 6]), h3_susie_l30 = as.numeric(m[, 7]))
}
l30 <- parse_l30(P$LOG_30_candidate_audit)
l30txt <- paste(rdlines(P$LOG_30_candidate_audit), collapse = "\n")
chk("log30.rows", nrow(l30) == 26L && grepl("dashed_edge\\s+26", l30txt), nrow(l30), 26, "LOG_30 dashed-edge table rows")
parse_l31 <- function(p) {
  x <- rdlines(p)
  m <- regmatches(x, regexec("^\\s*\\d+:\\s+(\\S+)\\s+(\\S+)\\s+([0-9.e+-]+)\\s+([0-9.e+-]+|NA)\\s+(\\d+/\\d+)\\s+(\\d+)\\s+(ok|collapsed|error|ERROR)(\\s+(TRUE|FALSE))?\\s*$", x))
  m <- do.call(rbind, m[lengths(m) == 10])
  d <- data.table(gene = m[, 2], outcome_name = m[, 3], h4_orig_l31 = as.numeric(m[, 4]), h4_strict_l31 = suppressWarnings(as.numeric(m[, 5])),
                  cs_strict_l31 = m[, 6], n_warn_l31 = as.integer(m[, 7]), run_status_l31 = m[, 8], survives_l31 = m[, 10])
  # the non-eGFR log wraps the 'survives' column onto a continuation block: recover it
  if (any(d$survives_l31 == "")) {
    s2 <- regmatches(x, regexec("^\\s*(\\d+):\\s+(TRUE|FALSE)\\s+(<NA>|\\S+)\\s*$", x)); s2 <- do.call(rbind, s2[lengths(s2) == 4])
    stopifnot(nrow(s2) == nrow(d)); d$survives_l31 <- s2[order(as.integer(s2[, 2])), 3]
  }
  d
}
l31 <- rbind(parse_l31(P$LOG_31_candidate_egfr), parse_l31(P$LOG_31_candidate_nonegfr))
chk("log31.rows", nrow(l31) == 26L, nrow(l31), 26, "LOG_31 strict re-run rows (14 eGFR + 12 non-eGFR)")
chk("log31.survives_parsed", all(l31$survives_l31 %in% c("TRUE", "FALSE")), paste(unique(l31$survives_l31), collapse = ","), "TRUE/FALSE",
    "survives column recovered for every LOG_31 row")
chk("log31.survives_consistent", all((l31$survives_l31 == "TRUE") == (!is.na(l31$h4_strict_l31) & l31$h4_strict_l31 >= 0.80)),
    sum((l31$survives_l31 == "TRUE") != (!is.na(l31$h4_strict_l31) & l31$h4_strict_l31 >= 0.80)), 0, "LOG_31 survives == strict PP.H4 >= 0.80")
fromlogs <- l31[survives_l31 == "TRUE" & n_warn_l31 == 0L, .(gene, outcome_name)]
held <- l31[survives_l31 == "TRUE" & n_warn_l31 > 0L, .(gene, outcome_name)]
chk("cand.routes_agree", nrow(drawn) == 16L && setequal(pk(drawn), pk(fromlogs)), sprintf("edge_audit=%d logs=%d", nrow(drawn), nrow(fromlogs)),
    "16 identical", "candidates: TAB_edge_audit Status==Drawn vs (LOG_31 strict PP.H4>=0.80 with a converged fit, n_warn==0)")
lg("held (strict PP.H4>=0.80 only with a convergence warning): ", paste(pk(held), collapse = ","))
chk("cand.in_pool", all(pk(drawn) %in% pk(tst)) && setequal(pk(l30), pk(l31)) && all(pk(l30) %in% pk(tst)),
    sum(pk(drawn) %in% pk(tst)), 16, "candidates and the 26 dashed edges are inside the rebuilt 160")
tst[, cand := pk(tst) %in% pk(drawn)]
chk("pool.candidates", sum(tst$cand) == 16L, sum(tst$cand), 16, "§6(c) assertion target: 16 candidates")

# ---------------------------------------------------------------- 9. tables, Fisher tests, ORs
tab4 <- function(cand, kid) c(a = sum(cand & kid), b = sum(cand & !kid), c = sum(!cand & kid), d = sum(!cand & !kid))
fish <- function(t) {
  m <- matrix(c(t["a"], t["c"], t["b"], t["d"]), nrow = 2)   # rows: candidate / other; cols: kidney / non-kidney
  f <- fisher.test(m)                                         # p value only (an exact sum; no root finding)
  e <- exact_cond_or(m, conf.level = EXACT_CONF, tol = EXACT_TOL)   # FIXER F1: MLE and exact CI, roots at tol 1e-13
  ft <- c(unname(f$estimate), f$conf.int)
  c(t, p = f$p.value, or_mle = e$estimate, or_lo = e$conf.int[1], or_hi = e$conf.int[2],
    resid_max = max(abs(e$resid)), ft_mle = ft[1], ft_lo = ft[2], ft_hi = ft[3], ft_reldiff = max(abs(c(e$estimate, e$conf.int) - ft) / ft),
    or_sample = unname(t["a"] * t["d"] / (t["b"] * t["c"])), or_haldane = unname((t["a"] + .5) * (t["d"] + .5) / ((t["b"] + .5) * (t["c"] + .5))))
}
FTNOTE <- function(v) sprintf(" | R fisher.test (uniroot at its default tolerance) gave %.10g; superseded by exact_cond_or.R (FIXLOG F1)", v)
emit_fisher <- function(prefix, ref, estimand, r, notes) {
  chk(paste0("fisher.roots.", prefix), r[["resid_max"]] < 1e-10, sprintf("%.2g", r[["resid_max"]]), "<1e-10",
      "exact_cond_or: max |residual| of the defining equations E[X]=x, P(X>=x)=0.025, P(X<=x)=0.025 at the reported MLE and bounds (FIXLOG F1)")
  chk(paste0("fisher.matches_fishertest.", prefix), r[["ft_reldiff"]] < 0.01, sprintf("%.2g", r[["ft_reldiff"]]), "<0.01",
      "same estimand as fisher.test: MLE and bounds within 1% of fisher.test's default-tolerance values (FIXLOG F1)")
  for (k in c("a", "b", "c", "d")) addrow(paste0(prefix, ".", k), ref, estimand, r[[k]], resource = "rebuilt 160-pair pool", unit = "count",
                                          n = r[["a"]] + r[["b"]] + r[["c"]] + r[["d"]],
                                          notes = paste0(notes, if (k == "a") " | a=candidate&kidney, b=candidate&non-kidney, c=other&kidney, d=other&non-kidney" else ""))
  addrow(paste0(prefix, ".p"), ref, estimand, r[["p"]], p = r[["p"]], unit = "p (Fisher exact, two-sided, R fisher.test)", notes = notes)
  addrow(paste0(prefix, ".or_mle"), ref, estimand, r[["or_mle"]], ci_lo = r[["or_lo"]], ci_hi = r[["or_hi"]], p = r[["p"]],
         unit = "OR (conditional MLE)", notes = paste0(notes, " | exact 95% CI: fisher.test's equations solved at tol 1e-13 (exact_cond_or.R)", FTNOTE(r[["ft_mle"]])))
  addrow(paste0(prefix, ".or_lo"), ref, estimand, r[["or_lo"]], unit = "OR exact 95% CI lower", notes = paste0(notes, FTNOTE(r[["ft_lo"]])))
  addrow(paste0(prefix, ".or_hi"), ref, estimand, r[["or_hi"]], unit = "OR exact 95% CI upper", notes = paste0(notes, FTNOTE(r[["ft_hi"]])))
  addrow(paste0(prefix, ".or_sample"), ref, estimand, r[["or_sample"]], unit = "OR (sample, ad/bc)", notes = "descriptive")
}
EST_PAIR <- "kidney (eGFR or CKD) share among multi-signal candidates vs other testable associations; unit = gene-outcome pair; pool = 160 (FDR second-order, failed coloc, MHC excluded)"

r_pair <- fish(tab4(tst$cand, tst$kid))
chk("primary.table", all(r_pair[c("a", "b", "c", "d")] == c(11, 5, 57, 87)), paste(r_pair[c("a", "b", "c", "d")], collapse = "/"), "11/5/57/87",
    "pre-declaration §6 'Primary, kept: 11/16 vs 57/144'")
emit_fisher("S6.primary", "§6 Primary", EST_PAIR, r_pair, "primary pair-level Fisher test")
addrow("S6.primary.or_haldane", "§6 Primary / §6(c)", EST_PAIR, r_pair[["or_haldane"]], unit = "OR (Haldane +0.5 on every cell)", notes = "observed table")
emit_fisher("S6.b.pair", "§6(b)", EST_PAIR, r_pair, "(b) pair level: conditional-MLE OR with exact 95% CI (identical to the primary table)")

# (a) H4/H3 >= 5 on the SuSiE route too (coloc.abf route already requires it)
sus <- merge(l30, l31, by = c("gene", "outcome_name"))
sus[, ratio_susie := ifelse(h3_susie_l30 == 0, Inf, h4_susie_l30 / h3_susie_l30)]
sus[, susie_sym_pass := h4_susie_l30 >= 0.80 & ratio_susie >= 5]
chk("a.strict_equals_naive_h4", all(abs(sus[pk(sus) %in% pk(drawn), h4_strict_l31 - h4_susie_l30]) < 5e-7),
    sprintf("max|diff|=%.2g", max(abs(sus[pk(sus) %in% pk(drawn), h4_strict_l31 - h4_susie_l30]))), "<5e-7",
    "for the 16 candidates the purity-filtered PP.H4 (LOG_31) equals the naive PP.H4 (LOG_30), so LOG_30's PP.H3 is taken for the purity fit")
fail_a <- sus[pk(sus) %in% pk(drawn) & susie_sym_pass == FALSE]
tst[, cand_a := cand & !(pk(tst) %in% pk(fail_a))]
r_a <- fish(tab4(tst$cand_a, tst$kid))
addrow("S6.a.reclassified", "§6(a)", "candidates failing H4/H3>=5 on the SuSiE route (reclassified as other testable associations)",
       paste(sprintf("%s (SuSiE PP.H4 %.7g, PP.H3 %.7g, ratio %.4f)", pk(fail_a), fail_a$h4_susie_l30, fail_a$h3_susie_l30, fail_a$ratio_susie), collapse = "; "),
       resource = "LOG_30_candidate_audit (SuSiE PP.H4/PP.H3, 7 significant digits) + TAB_coloc_final", n = nrow(fail_a), unit = "pairs",
       notes = "rule: SuSiE PP.H4 >= 0.80 and PP.H4/PP.H3 >= 5 (PP.H3 = 0 passes); a failing candidate stays in the 160 as an 'other' pair")
for (i in seq_len(nrow(sus[pk(sus) %in% pk(drawn)]))) {
  x <- sus[pk(sus) %in% pk(drawn)][i]
  addrow(sprintf("S6.a.susie_ratio.%s~%s", x$gene, x$outcome_name), "§6(a)", "SuSiE PP.H4/PP.H3 for a drawn candidate", x$ratio_susie,
         resource = "LOG_30_candidate_audit", unit = "ratio", notes = sprintf("PP.H4 %.7g, PP.H3 %.7g (log precision 7 sig. digits); passes=%s", x$h4_susie_l30, x$h3_susie_l30, x$susie_sym_pass))
}
emit_fisher("S6.a", "§6(a)", paste0(EST_PAIR, "; candidates restricted by the H4/H3>=5 rule on both routes"), r_a, "(a) symmetric acceptance rule")
emit_fisher("S6.b.a_pair", "§6(b) for (a)", paste0(EST_PAIR, "; (a) rule"), r_a, "(b) conditional-MLE OR with exact CI for (a)")

# clusters: 1-Mb single linkage of gene starts
chain <- function(d, by) {
  d <- copy(d); d[, chrn := as.integer(chr)]
  setorderv(d, c(by, "chrn", "start"))
  grpkey <- if (length(by)) do.call(paste, c(as.list(d[, ..by]), list(d$chrn))) else as.character(d$chrn)
  newc <- c(TRUE, grpkey[-1] != grpkey[-nrow(d)] | diff(d$start) > 1e6)
  d[, locus := cumsum(newc)]
  d
}
chk("chain.chr_numeric", !anyNA(suppressWarnings(as.integer(tst$chr))), sum(is.na(suppressWarnings(as.integer(tst$chr)))), 0, "all pool chromosomes are autosomes 1-22")
tst[, grp := ifelse(kid, "kidney", outcome_name)]
lo_ <- chain(tst, "grp")
lc <- lo_[, .(kid = kid[1], cand = any(cand), n = .N), by = locus]
chk("locus_outcome.kid_homog", all(lo_[, uniqueN(kid), by = locus]$V1 == 1), 0, 0, "locus-outcome groups are homogeneous in kidney status by construction")
r_lo <- fish(tab4(lc$cand, lc$kid))
addrow("S6.pool.locus_outcome_groups", "§6(b)", "locus-outcome groups (within outcome group, eGFR+CKD merged; 1-Mb single linkage of gene starts)", nrow(lc), unit = "groups")
emit_fisher("S6.b.locus_outcome", "§6(b)", "kidney share among candidate locus-outcome groups vs other groups (a group is a candidate if it holds any candidate)", r_lo,
            "(b) locus-outcome grouping: conditional-MLE OR with exact CI")

gcl <- chain(tst, character(0))
setorder(gcl, locus, chrn, start, gene, outcome_name)
fwrite(gcl[, .(gene, outcome_name, chr, start, locus, kid, cand, cand_a)], file.path(OUT, "w6a_pool160_genomic_clusters.tsv"), sep = "\t")
# intermediate export (added at the compare phase, 2026-10-08, to let COMPARE_SEEN diff se2/p2; it changes no computed value)
p177 <- merge(pool[, .(gene, outcome_name, method, b, se, se2, p2, q2, inst_case, PP.H3, PP.H4, coloc_ratio = ratio, chr, start, end, mhc)],
              tst[, .(gene, outcome_name, kid, cand, cand_a)], by = c("gene", "outcome_name"), all.x = TRUE)
setorder(p177, -mhc, gene, outcome_name)
fwrite(p177, file.path(OUT, "w6a_pool177_rebuild.tsv"), sep = "\t")
cl <- gcl[, .(n = .N, anykid = any(kid), allkid = all(kid), cand = any(cand), cand_a = any(cand_a),
              genes = paste(sort(unique(gene)), collapse = ","), outs = paste(sort(unique(outcome_name)), collapse = ",")), by = locus]
nclu <- nrow(cl); nmix <- sum(cl$anykid & !cl$allkid)
chk("pool.clusters", nclu == 113L, nclu, 113, "§6(c) assertion target: 113 genomic clusters (1-Mb single linkage of gene starts, chained across all outcomes)")
chk("pool.mixed", nmix == 14L, nmix, 14, "§6(c) assertion target: 14 mixed kidney/non-kidney clusters")
hal <- r_pair[["or_haldane"]]
chk("c.haldane_or", round(hal, 4) == 3.1818, sprintf("%.6f", hal), "3.1818", "§6(c) assertion target: Haldane +0.5 OR on the observed pair table")
chk("sm.pair_counts", all(r_pair[c("a", "b", "c", "d")] == c(11, 5, 57, 87)) && round(r_pair[["p"]], 3) == 0.033,
    sprintf("%d/%d vs %d/%d p=%.4f", r_pair[["a"]], r_pair[["a"]] + r_pair[["b"]], r_pair[["c"]], r_pair[["c"]] + r_pair[["d"]], r_pair[["p"]]),
    "11/16 vs 57/144, p = 0.033", "DEF_SM_CANDIDATES manuscript counts (weak evidence: reproduction target)")
chk("sm.locus_counts", r_lo[["a"]] == 9 && r_lo[["a"]] + r_lo[["b"]] == 14 && r_lo[["c"]] == 44 && r_lo[["c"]] + r_lo[["d"]] == 120 && round(r_lo[["p"]], 3) == 0.080,
    sprintf("%d/%d vs %d/%d p=%.4f", r_lo[["a"]], r_lo[["a"]] + r_lo[["b"]], r_lo[["c"]], r_lo[["c"]] + r_lo[["d"]], r_lo[["p"]]),
    "9/14 vs 44/120, p = 0.080", "DEF_SM_CANDIDATES manuscript locus counts (weak evidence: reproduction target)")
chk("sm.kidney_counts", sum(pool$outcome_name %in% KID) == 80L && sum(tst$kid) == 68L, c(sum(pool$outcome_name %in% KID), sum(tst$kid)), c(80, 68),
    "DEF_SM_CANDIDATES: 80 kidney of 177; 68 kidney of 160")

addrow("S6.pool.fdr_rows", "§6 pool rules (DEF_SM_CANDIDATES, LOG_28)", "BH q<=0.05 over 15,584 tests, second-order SE on Wald ratios, first-order on IVW; untestable padded p=1", nrow(fdr),
       resource = "TAB_mr_all_results + TAB_instruments", se_order = "2nd (Wald); 1st (IVW)", unit = "pairs", notes = sprintf("BH p threshold %.6g; 21 multi-instrument-gene Wald rows bounded both ways, FDR set invariant=%s", thr, same_fdr))
addrow("S6.pool.bh_threshold", "§6 pool rules", "largest second-order p passing BH at 0.05", thr, unit = "p", se_order = "2nd (Wald); 1st (IVW)")
addrow("S6.pool.colocfail", "§6 pool rules", "FDR pairs failing the coloc rule (not H4>=0.80 with H4/H3>=5)", nrow(pool), resource = "TAB_coloc_final", unit = "pairs")
addrow("S6.pool.mhc", "§6 pool rules", "coloc-fail pairs whose gene overlaps chr6 25-34 Mb GRCh38 (excluded)", sum(pool$mhc), resource = "TAB_instruments gene coords", unit = "pairs")
addrow("S6.pool.rows", "§6(c) assert", "testable pool (MHC excluded)", nrow(tst), unit = "pairs", notes = sprintf("eGFR %d, non-eGFR %d, kidney (eGFR+CKD) %d; identical to the 160 pairs LOG_29 ran", sum(tst$outcome_name == "eGFR"), sum(tst$outcome_name != "eGFR"), sum(tst$kid)))
addrow("S6.pool.candidates", "§6(c) assert", "multi-signal candidates (Drawn)", sum(tst$cand), resource = "TAB_edge_audit; cross-checked against LOG_31", unit = "pairs",
       notes = paste("candidates:", paste(pk(tst[cand == TRUE]), collapse = ",")))
addrow("S6.pool.clusters", "§6(c) assert", "genomic clusters: 1-Mb single linkage of gene starts chained across all outcomes", nclu, unit = "clusters")
addrow("S6.pool.mixed", "§6(c) assert", "clusters holding both kidney and non-kidney pairs", nmix, unit = "clusters")

# ---------------------------------------------------------------- 10. (c) genomic-cluster bootstrap of the pair-level OR
K <- as.matrix(gcl[, .(a = sum(cand & kid), b = sum(cand & !kid), c = sum(!cand & kid), d = sum(!cand & !kid)), by = locus][order(locus), .(a, b, c, d)])
stopifnot(nrow(K) == nclu, all(colSums(K) == r_pair[c("a", "b", "c", "d")]))
CK <- file.path(OUT, "boot_ckpt"); dir.create(CK, showWarnings = FALSE)
SEEDS <- 20261007:20261016; B <- 20000L
boots <- list()
for (sd in SEEDS) {
  f <- file.path(CK, sprintf("seed_%d.rds", sd))
  if (file.exists(f)) {
    v <- readRDS(f); if (length(v) == B && !anyNA(v)) { boots[[as.character(sd)]] <- v; lg("(c) seed ", sd, ": resumed from checkpoint ", basename(f), " md5 ", md5f(f)); next }
  }
  t1 <- Sys.time()
  set.seed(sd, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  idx <- matrix(sample.int(nclu, nclu * B, replace = TRUE), nrow = B, byrow = TRUE)   # row = one draw, in generation order
  A <- matrix(K[, "a"][idx], nrow = B); Bm <- matrix(K[, "b"][idx], nrow = B); Cm <- matrix(K[, "c"][idx], nrow = B); Dm <- matrix(K[, "d"][idx], nrow = B)
  v <- log((rowSums(A) + .5) * (rowSums(Dm) + .5) / ((rowSums(Bm) + .5) * (rowSums(Cm) + .5)))
  stopifnot(length(v) == B, !anyNA(v))
  saveRDS(v, f)
  boots[[as.character(sd)]] <- v
  lg(sprintf("(c) seed %d: %d draws in %.1fs; 95%% percentile CI %.4f-%.4f; share OR<=1 %.5f; checkpoint %s md5 %s", sd, B, as.numeric(difftime(Sys.time(), t1, units = "secs")),
             exp(quantile(v, .025, type = 7)), exp(quantile(v, .975, type = 7)), mean(v <= 0), basename(f), md5f(f)))
}
allb <- unlist(boots, use.names = FALSE)
chk("c.draw_count", length(allb) == length(SEEDS) * B && !anyNA(allb), length(allb), length(SEEDS) * B, "pooled bootstrap draws, no NA")
q <- exp(quantile(allb, c(.025, .975), type = 7))
addrow("S6.c.haldane_or", "§6(c)", EST_PAIR, hal, unit = "OR (Haldane +0.5)", notes = "observed table; asserted 3.1818")
addrow("S6.c.pooled_lo", "§6(c)", "pair-level Haldane OR, genomic-cluster bootstrap (113 clusters resampled with replacement)", q[[1]], n = length(allb), unit = "OR, 2.5th percentile",
       notes = "200,000 draws pooled (10 seeds x 20,000); percentile on log OR, R quantile type 7 (= numpy linear); Monte Carlo error not asserted")
addrow("S6.c.pooled_hi", "§6(c)", "pair-level Haldane OR, genomic-cluster bootstrap", q[[2]], n = length(allb), unit = "OR, 97.5th percentile", notes = "as pooled_lo")
addrow("S6.c.pooled_ci", "§6(c)", "pair-level Haldane OR with pooled genomic-cluster bootstrap 95% CI", hal, ci_lo = q[[1]], ci_hi = q[[2]], n = length(allb), unit = "OR")
addrow("S6.c.share_le1", "§6(c)", "share of bootstrap draws with Haldane OR <= 1", mean(allb <= 0), n = length(allb), unit = "share")
for (sd in SEEDS) {
  v <- boots[[as.character(sd)]]; qq <- exp(quantile(v, c(.025, .975), type = 7))
  addrow(sprintf("S6.c.seed%d.lo", sd), "§6(c)", "per-seed 2.5th percentile", qq[[1]], n = B, unit = "OR")
  addrow(sprintf("S6.c.seed%d.hi", sd), "§6(c)", "per-seed 97.5th percentile", qq[[2]], n = B, unit = "OR")
  addrow(sprintf("S6.c.seed%d.share_le1", sd), "§6(c)", "per-seed share OR <= 1", mean(v <= 0), n = B, unit = "share")
}
sl <- vapply(boots, function(v) exp(quantile(v, .025, type = 7)), 0); sh <- vapply(boots, function(v) exp(quantile(v, .975, type = 7)), 0)
addrow("S6.c.seedrange_lo", "§6(c)", "range of per-seed lower bounds", sprintf("%.4f-%.4f", min(sl), max(sl)), unit = "OR")
addrow("S6.c.seedrange_hi", "§6(c)", "range of per-seed upper bounds", sprintf("%.4f-%.4f", min(sh), max(sh)), unit = "OR")

# ---------------------------------------------------------------- 11. G1-G4 (cluster-level Fisher variants; supplement only)
GNOTE <- "SUPPLEMENT-ONLY (§12 fate; §6: mixed-cluster Fisher variants go to the SM). Definitions from the pinned estimator renal_enrichment.py (the redacted §12 does not carry them), recomputed on the W6a rebuild"
gdefs <- list(
  G1 = list(d = cl, kid = cl$anykid, cand = cl$cand, e = "genomic clusters; kidney = any kidney pair; candidate = any candidate"),
  G2 = list(d = cl, kid = cl$allkid, cand = cl$cand, e = "genomic clusters; kidney = all pairs kidney (mixed clusters counted non-kidney); candidate = any candidate"),
  G3 = list(d = cl[!(anykid & !allkid)], kid = cl[!(anykid & !allkid)]$anykid, cand = cl[!(anykid & !allkid)]$cand, e = "genomic clusters with the mixed clusters dropped"),
  G4 = list(d = cl, kid = cl$anykid, cand = cl$cand_a, e = "G1 with the (a)-rule failures reclassified as non-candidates (the estimator hard-codes GUSB-eGFR)"))
for (g in names(gdefs)) {
  r <- fish(tab4(gdefs[[g]]$cand, gdefs[[g]]$kid))
  addrow(sprintf("S6.%s.n_clusters", g), "§12 (G-variants)", gdefs[[g]]$e, length(gdefs[[g]]$cand), unit = "clusters", notes = GNOTE)
  emit_fisher(sprintf("S6.%s", g), "§12 (G-variants)", gdefs[[g]]$e, r, GNOTE)
}
chk("G4.reclass_is_GUSB", identical(pk(fail_a), "GUSB~eGFR"), paste(pk(fail_a), collapse = ","), "GUSB~eGFR",
    "the (a)-rule reclassification set equals the estimator's hard-coded G4 reclassification (GUSB-eGFR)")

# ---------------------------------------------------------------- 12. failure modes of the 16 candidates (coloc.abf route)
fm <- tst[cand == TRUE, .(gene, outcome_name, kid, PP.H4, PP.H3, ratio)]
fm[, mode := fifelse(PP.H4 >= 0.80 & ratio < 5, "ratio-only (coloc.abf PP.H4>=0.80, PP.H4/PP.H3<5)",
             fifelse(PP.H3 > PP.H4, "H3-dominant (PP.H4<0.80, PP.H3>PP.H4: distinct variants favoured)",
                     "H4-leaning (PP.H4<0.80, PP.H4>=PP.H3)"))]
chk("fail.partition", nrow(fm) == 16L && !anyNA(fm$mode) && all(fm$PP.H4 < 0.80 | fm$ratio < 5), nrow(fm), 16, "every candidate failed coloc.abf, and falls in exactly one mode")
fmt_ <- fm[, .(n = .N, n_kid = sum(kid), pairs = paste(sprintf("%s~%s", gene, outcome_name), collapse = ",")), by = mode][order(-n)]
chk("fail.961", identical(as.integer(fmt_$n), c(9L, 6L, 1L)), paste(fmt_$n, collapse = "/"), "9/6/1", "pre-declaration §6 'the 9 / 6 / 1 failure-mode classification'")
fwrite(fm[order(mode, gene)], file.path(OUT, "w6a_drawn16_failure_modes.tsv"), sep = "\t")
nm_ <- c("n9", "n6", "n1")
for (i in seq_len(min(3, nrow(fmt_)))) addrow(paste0("S6.fail.", nm_[i]), "§6 'Also built: 9/6/1'", fmt_$mode[i], fmt_$n[i], resource = "TAB_coloc_final (full precision)",
                                            unit = "candidates", notes = paste0("kidney ", fmt_$n_kid[i], "; ", fmt_$pairs[i], " | id assigned by rank of class size"))
addrow("S6.fail.kidney_in_9", "§6 'Also built'", paste("kidney candidates in the largest class:", fmt_$mode[1]), fmt_$n_kid[1], unit = "candidates")
addrow("S6.fail.kidney_in_other", "§6 'Also built'", "kidney candidates in the other classes", sum(fmt_$n_kid[-1]), unit = "candidates")

# ---------------------------------------------------------------- 13. SuSiE fit ledger from LOG_29/30/31 (log precision)
parse_l29 <- function(p, mode) {
  x <- rdlines(p); hs <- grep("^\\[(\\d+)/(\\d+)\\] (\\S+) ~ (\\S+) \\(coloc", x)
  stopifnot(length(hs) > 0); ends <- c(hs[-1] - 1, length(x))
  g <- function(blk, re, n = 1) { m <- regmatches(blk, regexec(re, blk, perl = TRUE)); m <- m[lengths(m) > 0]; if (!length(m)) rep(NA_character_, n) else m[[1]][2:(n + 1)] }
  rbindlist(lapply(seq_along(hs), function(i) {
    blk <- x[hs[i]:ends[i]]
    h <- regmatches(blk[1], regexec("^\\[(\\d+)/(\\d+)\\] (\\S+) ~ (\\S+) \\(coloc\\.abf PP\\.H4=([0-9.]+)\\)", blk[1]))[[1]]
    win <- g(blk, "susie: \\S+ ~ (\\S+) \\| chr(\\w+):(\\d+)-(\\d+)", 4)
    cs <- g(blk, "SuSiE (?:CS|credible sets): exp=(\\d+), out=(\\d+)", 2)
    co <- g(blk, "coloc\\.susie: (\\d+) (?:signal )?pairs \\| (?:best )?PP\\.H4=([0-9.]+) PP\\.H3=([0-9.]+)", 3)
    data.table(mode = mode, k = as.integer(h[2]), gene = h[4], outcome_name = h[5], abf_h4_log3 = as.numeric(h[6]),
               outcome_file = win[1], window_chr = win[2], window_start = as.numeric(win[3]), window_end = as.numeric(win[4]),
               n_merged = as.integer(g(blk, "merged: (\\d+) SNPs")),
               n_ld_aligned = as.integer(g(blk, "LD-aligned: (\\d+) SNPs")),
               ld_matrix_dim = as.integer(g(blk, "\\[LD local\\] (?:Returned )?(\\d+) x \\d+ matrix")),
               outcome_N_regional = as.integer(g(blk, "outcome N: (\\d+) \\(regional")),
               strand_align = g(blk, "\\[strand_align [^]]+\\] (.*)$"),
               cs_exp = as.integer(cs[1]), cs_out = as.integer(cs[2]), n_signal_pairs = as.integer(co[1]),
               susie_h4_log3 = as.numeric(co[2]), susie_h3_log3 = as.numeric(co[3]))
  }))
}
fl <- rbind(parse_l29(P$LOG_29_cand_egfr, "egfr"), parse_l29(P$LOG_29_cand_nonegfr, "nonegfr"))
chk("fitledger.rows", nrow(fl) == 160L && setequal(pk(fl), pk(tst)), nrow(fl), 160, "one LOG_29 fit record per pool pair")
chk("fitledger.complete", !anyNA(fl[, .(window_start, window_end, cs_exp, cs_out, n_signal_pairs, susie_h4_log3, susie_h3_log3, n_ld_aligned)]),
    sum(is.na(fl[, .(window_start, window_end, cs_exp, cs_out, n_signal_pairs, susie_h4_log3, susie_h3_log3, n_ld_aligned)])), 0,
    "no NA in window / CS / coloc.susie fields of the 160 parsed records")
fl <- merge(fl, sus[, .(gene, outcome_name, h4_abf_l30, h4_susie_l30, h3_susie_l30, ratio_susie, susie_sym_pass, h4_strict_l31, cs_strict_l31, n_warn_l31, run_status_l31, survives_l31)],
            by = c("gene", "outcome_name"), all.x = TRUE)
fl[, rescued_naive := susie_h4_log3 >= 0.80]
fl[, candidate := pk(fl) %in% pk(drawn)]
fl[, precision := "LOG_29: PP.H4/PP.H3 3 decimals, counts exact; LOG_30: 7 significant digits; LOG_31 strict PP.H4: as printed (7-10 digits)"]
setorder(fl, mode, k)
fwrite(fl, file.path(OUT, "w6a_fit_ledger.tsv"), sep = "\t")
r29 <- c(egfr = num_from(paste(rdlines(P$LOG_29_cand_egfr), collapse = "\n"), "rescued\\(PP\\.H4_susie>=0\\.80\\)=(\\d+)"),
         nonegfr = num_from(paste(rdlines(P$LOG_29_cand_nonegfr), collapse = "\n"), "rescued\\(PP\\.H4_susie>=0\\.80\\)=(\\d+)"))
chk("fitledger.rescued", sum(fl[mode == "egfr"]$rescued_naive) == r29[["egfr"]] && sum(fl[mode == "nonegfr"]$rescued_naive) == r29[["nonegfr"]] &&
      setequal(pk(fl[rescued_naive == TRUE]), pk(l30)),
    c(sum(fl[mode == "egfr"]$rescued_naive), sum(fl[mode == "nonegfr"]$rescued_naive)), r29, "parsed 3-dp PP.H4>=0.80 count per log equals the log's own 'rescued' tally, and equals LOG_30's 26")
nconv <- c(egfr = length(grep("IBSS algorithm did not converge", rdlines(P$LOG_29_cand_egfr))), nonegfr = length(grep("IBSS algorithm did not converge", rdlines(P$LOG_29_cand_nonegfr))))
FLR <- "LOG_29/30/31 (pipeline SuSiE fits, E: engine; parsed, not refitted)"
addrow("S6.fitledger.n_fits", "§6 'Also built: SuSiE fit ledger'", "first-pass SuSiE+coloc.susie fits parsed from LOG_29", nrow(fl), resource = FLR, unit = "fits", notes = "per-fit table: w6a_fit_ledger.tsv")
addrow("S6.fitledger.n_rescued_naive", "§6 fit ledger", "fits with naive SuSiE PP.H4 >= 0.80", sum(fl$rescued_naive), resource = FLR, unit = "fits", notes = sprintf("eGFR %d, non-eGFR %d (log precision 3 dp)", r29[["egfr"]], r29[["nonegfr"]]))
addrow("S6.fitledger.n_strict_survive", "§6 fit ledger", "purity-filtered (min_abs_corr 0.5) refits with PP.H4 >= 0.80", sum(l31$survives_l31 == "TRUE"), resource = FLR, unit = "fits")
addrow("S6.fitledger.n_strict_collapsed", "§6 fit ledger", "purity-filtered refits collapsing to no pure credible set", sum(l31$run_status_l31 == "collapsed"), resource = FLR, unit = "fits")
addrow("S6.fitledger.n_strict_below", "§6 fit ledger", "purity-filtered refits with PP.H4 < 0.80", sum(l31$run_status_l31 == "ok" & l31$survives_l31 == "FALSE"), resource = FLR, unit = "fits")
addrow("S6.fitledger.held_nonconverged", "§6 fit ledger", "strict PP.H4 >= 0.80 only with a convergence warning (held, not drawn)", paste(pk(held), collapse = ","), resource = FLR, n = nrow(held), unit = "pairs")
addrow("S6.fitledger.ibss_warnings_log29", "§6 fit ledger", "IBSS non-convergence warnings printed in LOG_29 (not attributable to fits in the log)", sum(nconv), resource = FLR, unit = "warnings",
       notes = sprintf("eGFR log %d, non-eGFR log %d; R prints the warnings at the end of the run, without the fit they came from", nconv[["egfr"]], nconv[["nonegfr"]]))
for (i in seq_len(nrow(fl[!is.na(h4_susie_l30)]))) {
  x <- fl[!is.na(h4_susie_l30)][i]; pre <- sprintf("S6.fitledger.%s~%s", x$gene, x$outcome_name)
  addrow(paste0(pre, ".window"), "§6 fit ledger", "SuSiE window (LOG_29)", sprintf("chr%s:%d-%d", x$window_chr, as.integer(x$window_start), as.integer(x$window_end)), resource = FLR, unit = "GRCh38 bp")
  addrow(paste0(pre, ".pph4_naive"), "§6 fit ledger", "naive SuSiE coloc PP.H4", x$h4_susie_l30, resource = FLR, unit = "PP", notes = "LOG_30, 7 significant digits")
  addrow(paste0(pre, ".pph3_naive"), "§6 fit ledger", "naive SuSiE coloc PP.H3", x$h3_susie_l30, resource = FLR, unit = "PP", notes = "LOG_30, 7 significant digits")
  addrow(paste0(pre, ".pph4_strict"), "§6 fit ledger", "purity-filtered SuSiE coloc PP.H4", if (is.na(x$h4_strict_l31)) "not estimable (collapsed: no pure credible set)" else x$h4_strict_l31, resource = FLR, unit = "PP", notes = paste0("LOG_31 as printed; status ", x$run_status_l31, "; n_warn ", x$n_warn_l31))
  addrow(paste0(pre, ".cs"), "§6 fit ledger", "credible sets exp/out naive (LOG_29) | strict (LOG_31)", sprintf("%d/%d | %s", x$cs_exp, x$cs_out, x$cs_strict_l31), resource = FLR, unit = "count")
  addrow(paste0(pre, ".status"), "§6 fit ledger", "candidate status", if (x$candidate) "Drawn (candidate)" else if (pk(x) %in% pk(held)) "Held (non-convergent)" else "Dropped", resource = FLR, unit = "label")
}

# ---------------------------------------------------------------- 14. write results.tsv
res <- rbindlist(ROWS)
chk("results.unique_ids", !anyDuplicated(res$quantity_id), sum(duplicated(res$quantity_id)), 0, "quantity ids unique")
chk("results.no_empty_values", all(nzchar(res$value) & res$value != "NA"), sum(!nzchar(res$value) | res$value == "NA"), 0, "no empty/NA value cells")
fwrite(rbindlist(ROWS), file.path(OUT, "results.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(ASSERT_LOG), file.path(OUT, "w6a_asserts.tsv"), sep = "\t", quote = FALSE)
lg("WROTE results.tsv rows=", nrow(rbindlist(ROWS)), " md5=", md5f(file.path(OUT, "results.tsv")), " | asserts: ", length(ASSERT_LOG), " run, ", NASSERT_FAIL, " failed")

# ---------------------------------------------------------------- 15. branches.tsv (partial: W6b completes §6 renal wording)
r3 <- function(x) sprintf("%.3f", x)
inc1 <- function(lo, hi) lo <= 1 && hi >= 1
q_lo <- q[[1]]; q_hi <- q[[2]]
grouping_hit <- inc1(r_lo[["or_lo"]], r_lo[["or_hi"]]) || inc1(q_lo, q_hi)
accept_hit <- inc1(r_a[["or_lo"]], r_a[["or_hi"]])
pairOR_le1 <- r_pair[["or_mle"]] <= 1
lst <- c(if (grouping_hit) "grouping", if (accept_hit) "the acceptance rule")
br <- list(
  data.table(rule_ref = "§6 Renal OR wording: selector 'pair OR <= 1'",
             inputs = "S6.primary.or_mle; S6.primary.or_lo; S6.primary.or_hi",
             branch = if (pairOR_le1) "pair OR <= 1" else "pair OR > 1 (otherwise-branch: add 'sensitive to {list}')",
             filled_template = if (pairOR_le1) sprintf("not more frequent among kidney pairs (OR %s, 95%% CI %s-%s)", r3(r_pair[["or_mle"]]), r3(r_pair[["or_lo"]]), r3(r_pair[["or_hi"]]))
                               else sprintf("{D6 base phrase} among kidney pairs (11/16 vs 57/144; OR %s, exact 95%% CI %s-%s; p %s)", r3(r_pair[["or_mle"]]), r3(r_pair[["or_lo"]]), r3(r_pair[["or_hi"]]), sprintf("%.2g", r_pair[["p"]])),
             qualifiers = "PARTIAL: W6b fills the final §6 renal sentence (D6 and the LD-guard qualifier need §6(d)); OR = conditional MLE with exact CI per §6(b)"),
  data.table(rule_ref = "§6 qualifier 'grouping' (locus-outcome exact CI or (c) includes 1)",
             inputs = "S6.b.locus_outcome.or_lo; S6.b.locus_outcome.or_hi; S6.c.pooled_lo; S6.c.pooled_hi",
             branch = if (grouping_hit) "includes 1 -> 'grouping' named" else "neither includes 1 -> 'grouping' not named",
             filled_template = sprintf("locus-outcome exact CI %s-%s (includes 1: %s); (c) pooled bootstrap CI %s-%s (includes 1: %s)", r3(r_lo[["or_lo"]]), r3(r_lo[["or_hi"]]),
                                       inc1(r_lo[["or_lo"]], r_lo[["or_hi"]]), r3(q_lo), r3(q_hi), inc1(q_lo, q_hi)),
             qualifiers = "travels with the renal OR claim (Principle 0)"),
  data.table(rule_ref = "§6 qualifier 'the acceptance rule' ((a) exact CI includes 1)",
             inputs = "S6.b.a_pair.or_lo; S6.b.a_pair.or_hi",
             branch = if (accept_hit) "includes 1 -> 'the acceptance rule' named" else "excludes 1 -> not named",
             filled_template = sprintf("(a) OR %s, exact 95%% CI %s-%s", r3(r_a[["or_mle"]]), r3(r_a[["or_lo"]]), r3(r_a[["or_hi"]])),
             qualifiers = "travels with the renal OR claim (Principle 0)"),
  data.table(rule_ref = "§6 qualifier 'the LD guard' ((d) CI includes 1)",
             inputs = "S6.d.lo; S6.d.hi (W6b)",
             branch = "PENDING_W6b", filled_template = "computed by W6b (§6(d))", qualifiers = "W6b"),
  data.table(rule_ref = "§11 D6 renal base phrase",
             inputs = "S6.primary.or_lo; S6.d.lo (W6b)",
             branch = if (round(r_pair[["or_lo"]], 3) > 1.000) "primary lower bound > 1.000 at 3 dp -> depends on (d): PENDING_W6b" else "primary lower bound <= 1.000 at 3 dp -> 'proportionally more frequent' (resolved without (d))",
             filled_template = if (round(r_pair[["or_lo"]], 3) > 1.000) sprintf("primary exact lower bound %s > 1.000; 'over-represented at the pair level' only if (d)'s exact lower bound also > 1.000, else 'proportionally more frequent'", r3(r_pair[["or_lo"]]))
                               else sprintf("proportionally more frequent (primary exact lower bound %s)", r3(r_pair[["or_lo"]])),
             qualifiers = sprintf("partial list from W6a: sensitive to {%s}%s", paste(lst, collapse = ", "), "; + 'the LD guard' if W6b's (d) CI includes 1")))
fwrite(rbindlist(br), file.path(OUT, "branches.tsv"), sep = "\t", quote = FALSE)
lg("WROTE branches.tsv rows=", length(br), " md5=", md5f(file.path(OUT, "branches.tsv")))

# ---------------------------------------------------------------- 16. end: staging snapshot
SNAP_END <- snap("end")
chk("staging.unchanged", identical(SNAP_START, SNAP_END), SNAP_END, SNAP_START, "nothing written under V4_local/staging/l3-finemap during compute")
lg("=== W6a compute end; elapsed ", sprintf("%.1f", as.numeric(difftime(Sys.time(), T0, units = "secs"))), " s; asserts failed: ", NASSERT_FAIL,
   " ; results.tsv md5 ", md5f(file.path(OUT, "results.tsv")), " ; branches.tsv md5 ", md5f(file.path(OUT, "branches.tsv")))
writeLines(c(md5f(file.path(OUT, "results.tsv")), md5f(file.path(OUT, "branches.tsv"))), file.path(OUT, "w6a_compute_md5.txt"))
