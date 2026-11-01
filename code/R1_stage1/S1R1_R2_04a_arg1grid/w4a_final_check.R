#!/usr/bin/env Rscript
# w4a_final_check.R -- W4a self-verification of the shipped outputs (row counts, id coverage, no NA where a value is
# expected, ranges, frozen prefix, staging untouched). Each coverage check is mutated in memory to prove it can fail.
#   Rscript w4a_final_check.R
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
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04a_arg1grid")
WORK <- file.path(OUT, "work")
source(file.path(OUT, "w4a_funcs.R"))
log_init(file.path(OUT, "run.log"))
logf("==== W4a final check start ====")
fail <- character(0)
ok <- function(cond, msg) { if (!isTRUE(cond)) fail <<- c(fail, msg); logf("CHECK %s: %s", if (isTRUE(cond)) "PASS" else "FAIL", msg) }
R <- fread(file = file.path(OUT, "results.tsv"), sep = "\t", colClasses = "character", na.strings = NULL)
B <- fread(file = file.path(OUT, "branches.tsv"), sep = "\t", colClasses = "character", na.strings = NULL)
C <- fread(file = file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", colClasses = "character", na.strings = NULL)
ok(identical(names(R), c("quantity_id", "prespec_ref", "estimand", "resource", "panel", "cov", "se_order", "value", "se", "ci_lo", "ci_hi", "p", "n", "unit", "notes")), "results.tsv columns as Rule 12")
ok(identical(names(B), c("rule_ref", "inputs", "branch", "filled_template", "qualifiers")), "branches.tsv columns as Rule 12")
ok(identical(names(C), c("quantity_id", "seen_value", "seen_label", "stage1_value", "stage1_label", "abs_diff", "explanation")), "COMPARE_SEEN.tsv columns as Rule 10")
ok(!anyDuplicated(R$quantity_id), sprintf("results.tsv quantity ids unique (%d rows)", nrow(R)))
fz <- readLines(file.path(WORK, "results.precompare_frozen.tsv")); cur <- readLines(file.path(OUT, "results.tsv"))
ok(identical(cur[seq_along(fz)], fz) && md5f(file.path(WORK, "results.precompare_frozen.tsv")) == readLines(file.path(WORK, "FREEZE.txt"))[1],
   sprintf("computed rows byte-identical to the pre-compare freeze (%d computed + %d appended)", length(fz) - 1, length(cur) - length(fz)))
app <- R[(length(fz)):nrow(R)]
ok(all(grepl("READ_FROM_SEEN", app$notes)), "every appended row is labelled READ_FROM_SEEN")
comp <- R[seq_len(length(fz) - 1)]
ok(!any(comp$value %in% c("", "NA", "NaN", "Inf", "-Inf")), "no NA/empty value in any computed row")
num <- suppressWarnings(as.numeric(comp$value))
nonnum <- comp$quantity_id[is.na(num)]
ok(all(grepl("^S4.ld.rgrid.crossing", nonnum)), sprintf("only r-grid 'none in grid' rows are non-numeric (%d: %s)", length(nonnum), paste(nonnum, collapse = ",")))
pr <- comp[grepl("(p_aware_2nd|p_ign|p_aware_1st|p_cojo|p_aware|p_lo|p_hi|share_lt05|share_same|share|power_point|power_pred|PPH4|pA|pB)$", comp$quantity_id)]
pv <- as.numeric(pr$value)
ok(all(is.finite(pv) & pv >= 0 & pv <= 1), sprintf("all %d p / share / power / PP.H4 values in [0, 1]", nrow(pr)))
# required id coverage (brief W4a)
req <- c(
  as.vector(outer(outer(c("FinnGen", "DIAMANTE"), c("marg", "cond"), paste, sep = "."),
                  c("ratioA", "ratioB", "diff", "diff_se", "diff_lo", "diff_hi", "p_aware_2nd", "p_ign", "p_aware_1st", "p_cojo", "exp_zB", "obs_zB", "power_point", "power_pred"),
                  function(a, b) paste0("S4.grid.", a, ".", b))),
  as.vector(outer(c("marg", "cond"), c("ratioA", "ratioB", "p_aware_2nd", "p_ign", "ror", "ror_lo", "ror_hi"), function(a, b) paste0("S4.arg.", a, ".", b))),
  "S4.arg.pA", "S4.arg.pB", "S4.ld.r.503", "S4.ld.r.489", "S4.ld.boot.r_lo", "S4.ld.boot.r_hi", "S4.ld.rgrid.crossing.aware", "S4.ld.rgrid.crossing.ign",
  as.vector(outer(c("FinnGen.cond.aware2nd", "DIAMANTE.marg.aware2nd", "arg.marg.aware2nd", "orn.marg.aware2nd", "t3.FinnGen.cond.aware2nd", "t3.FinnGen.marg.aware2nd"),
                  c("share_lt05", "p_lo", "p_hi"), function(a, b) paste0("S4.ld.boot.", a, ".", b))),
  as.vector(outer(outer(c("FinnGen", "DIAMANTE"), c("marg", "cond"), paste, sep = "."), outer(c("503", "489"), c("2nd", "1st"), paste, sep = "."),
                  function(a, b) paste0("S4.t3.", a, ".", b, ".p_aware"))),
  "S4.t3.boot.cond.share", "S4.t3.boot.marg.share",
  as.vector(outer(c("503", "489"), c("kA", "kB", "p_aware_2nd", "p_ign"), function(a, b) paste0("S4.orn.", a, ".", b))),
  as.vector(outer(outer(c("arginine", "ornithine", "T2D"), c("A", "B"), paste, sep = "."), outer(c("503", "489"), c("z", "PPH4"), paste, sep = "."),
                  function(a, b) paste0("S4.coloc_clause.", a, ".", b))))
cover <- function(ids) setdiff(req, ids)
miss <- cover(comp$quantity_id)
ok(!length(miss), sprintf("all %d brief-required quantity ids present (missing: %s)", length(req), if (length(miss)) paste(miss, collapse = ",") else "none"))
mut <- cover(setdiff(comp$quantity_id, "S4.t3.FinnGen.cond.503.2nd.p_aware"))
ok(identical(mut, "S4.t3.FinnGen.cond.503.2nd.p_aware"), "MUTATION: dropping S4.t3.FinnGen.cond.503.2nd.p_aware is caught by the coverage check (and only it)")
reqb <- c("4 Arginine opening clause", "4 Arginine A vs B", "4 T2D", "Amendment 1 A2", "4 Combined (test (iii))", "4 Ornithine",
          "Rule 3(ii) arginine A vs B (marginal, aware, 503)", "Rule 3(ii) T2D A vs B (FinnGen R10 conditional, aware)",
          "Rule 3(ii) A2: T2D A vs B DIAMANTE marginal aware (conditional share printed beside)", "Rule 3(ii) test (iii) conditional (FinnGen aware)",
          "Rule 3(ii) test (iii) marginal (FinnGen aware)", "Rule 3(ii) ornithine A vs B (marginal, aware, 503)",
          paste("Rule 3(iv)", as.vector(outer(c("arginine", "ornithine", "T2D"), c("A", "B"), function(a, b) paste(a, "signal", b)))), "Rule 4 family entry")
mb <- setdiff(reqb, B$rule_ref)
ok(!length(mb), sprintf("all %d required branch rows present (missing: %s)", length(reqb), if (length(mb)) paste(mb, collapse = ",") else "none"))
ok(identical(setdiff(reqb, setdiff(B$rule_ref, "4 T2D")), "4 T2D"), "MUTATION: dropping the '4 T2D' branch row is caught")
ok(sum(B$rule_ref == "Rule 4 family entry") == 4, "four Rule 4 family entries (arginine, T2D, test (iii), ornithine)")
ok(all(grepl("PENDING_HOLM", B$qualifiers[B$rule_ref %in% c("4 Arginine A vs B", "4 T2D", "4 Combined (test (iii))", "4 Ornithine")])), "family statements carry 'nominal: PENDING_HOLM'")
ok(nrow(C) > 0 && all(nzchar(C$explanation)), sprintf("COMPARE_SEEN.tsv has %d rows, each with an explanation", nrow(C)))
# staging untouched across every step since the main run's start snapshot
s0 <- fread(file = file.path(WORK, "staging_snapshot_start.tsv"), sep = "\t", colClasses = "character")
f <- list.files(paste0(.P31R_V3, "/V4_local/staging/l3-finemap"), recursive = TRUE, full.names = TRUE, all.files = TRUE)
s1 <- data.frame(file = f, md5 = unname(tools::md5sum(f)), stringsAsFactors = FALSE); s1 <- s1[order(s1$file), ]
ok(identical(s0$file, s1$file) && identical(s0$md5, s1$md5), sprintf("V4_local/staging/l3-finemap unchanged since the start snapshot (%d files)", nrow(s1)))
logf("final check: %d failures%s", length(fail), if (length(fail)) paste0(": ", paste(fail, collapse = " | ")) else "")
logf("output md5: results.tsv %s; branches.tsv %s; COMPARE_SEEN.tsv %s", md5f(file.path(OUT, "results.tsv")), md5f(file.path(OUT, "branches.tsv")), md5f(file.path(OUT, "COMPARE_SEEN.tsv")))
logf("==== W4a final check end ====")
if (length(fail)) quit(status = 1)
