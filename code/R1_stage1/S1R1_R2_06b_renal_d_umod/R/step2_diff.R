# step2_diff.R -- Stage 2 step 2 scope report (replaces the FIXER F1 scope gate, whose pre-fix baseline fix_F1_prefix/ differs
# structurally from any seven-fit output). Compares the regenerated results.tsv / branches.tsv with the pre-step-2 copies kept in
# _pre_step2_20261010/ (branches: after the Stage 1b normaliser, same as the regenerated file after 10b).
# FIX ROUND 2026-10-10: the gate is now ID-LIMITED and content-pinned (R/step2_gate_lib.R, ../step2_expected_changes.tsv). The ONLY
# differences from the baseline allowed are the (file, id, column) cells and added rows listed in step2_expected_changes.tsv (each with a
# one-line reason), and each listed cell must carry exactly the recorded content (md5). Any other difference, a listed cell with other
# content, a removed row or a listed change that did not happen prints a PROBLEM line naming the id and the column and the script
# exits non-zero. (Before: only a column-level check; +1e-3 on an unrelated value, a per-fit s or an s_max passed.)
# step2_diff.tsv still lists every differing cell with old and new text. Test overrides (mutation tests run on COPIES): environment
# variables W6B_NEW (dir with results.tsv/branches.tsv), W6B_EXPECTED (table), W6B_DIFFOUT (report path; default = the real one).
# [portability] declared transform (rel_portability.py): archive-relative roots. Environment P31R_ARCHIVE overrides; by default the archive root is
# derived from this file's place (<archive>/code/<path>). P31R_EXTERNAL points at third-party inputs that are not in the release.
if (!exists(".P31R_ARCHIVE")) {
  .p31r_f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .P31R_ARCHIVE <- if (nzchar(Sys.getenv("P31R_ARCHIVE"))) Sys.getenv("P31R_ARCHIVE") else {
    if (!length(.p31r_f)) stop("run with Rscript <file>, or set the environment variable P31R_ARCHIVE to the archive root")
    normalizePath(file.path(dirname(.p31r_f[1]), "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
  }
  .P31R_V3 <- file.path(.P31R_ARCHIVE, "code")
  .P31R_EXT <- Sys.getenv("P31R_EXTERNAL", file.path(.P31R_ARCHIVE, "external_inputs"))
  .P31R_CKM <- file.path(.P31R_EXT, "CKM papers")
}
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod/R/step2_gate_lib.R"))
D0 <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"); D <- Sys.getenv("W6B_NEW", D0); PRE <- file.path(D0, "_pre_step2_20261010")
EXPF <- Sys.getenv("W6B_EXPECTED", file.path(D0, "step2_expected_changes.tsv")); DIFFOUT <- Sys.getenv("W6B_DIFFOUT", file.path(D0, "step2_diff.tsv"))
rd <- function(p) fread(p, sep = "\t", colClasses = "character", na.strings = NULL, quote = "\"")
nr <- rd(file.path(D, "results.tsv")); orr <- rd(file.path(PRE, "results.tsv"))
nb <- rd(file.path(D, "branches.tsv")); ob <- rd(file.path(PRE, "branches.tsv"))
stopifnot(identical(names(nr), names(orr)), identical(names(nb), names(ob)))
out <- list(); add <- function(file, row, col, old, new) out[[length(out) + 1]] <<- data.table(file, row, column = col, old, new)
removed <- setdiff(orr$quantity_id, nr$quantity_id); added <- setdiff(nr$quantity_id, orr$quantity_id)
ok_new <- grepl("_4fit$|[.]s[.](protein|outcome)[.](p2np[.]503|susie[.]489|p1fr[.]489)$", added)
cat(sprintf("results: %d old rows, %d new rows; removed %d; added %d (of which seven-fit additions %d)\n", nrow(orr), nrow(nr), length(removed), length(added), sum(ok_new)))
m <- merge(orr, nr, by = "quantity_id", suffixes = c(".o", ".n"), sort = FALSE)
OKCOL <- c("value", "ci_lo", "ci_hi", "p", "notes", "estimand", "n", "unit", "prespec_ref")
bad <- 0L
for (cc in setdiff(names(orr), "quantity_id")) {
  a <- m[[paste0(cc, ".o")]]; b <- m[[paste0(cc, ".n")]]; d <- which(a != b)
  for (i in d) { add("results.tsv", m$quantity_id[i], cc, a[i], b[i]); if (!(cc %in% OKCOL)) bad <- bad + 1L }
}
if (nrow(nb) != nrow(ob)) stop("step2_diff: branches row count ", nrow(ob), " -> ", nrow(nb))
for (cc in names(ob)) for (i in which(ob[[cc]] != nb[[cc]])) add("branches.tsv", sprintf("row %d (%s)", i, ob$rule_ref[i]), cc, ob[[cc]][i], nb[[cc]][i])
G <- rbindlist(out); if (!nrow(G)) G <- data.table(file = character(), row = character(), column = character(), old = character(), new = character())
fwrite(G, DIFFOUT, sep = "\t", quote = TRUE)
cat(sprintf("step2_diff: results cells changed %d (in a column outside the allowed list: %d); branches cells changed %d (branch column: %d; rule_ref column: %d)\n",
            sum(G$file == "results.tsv"), bad, sum(G$file == "branches.tsv"), sum(G$file == "branches.tsv" & G$column == "branch"), sum(G$file == "branches.tsv" & G$column == "rule_ref")))
PROB <- g_run(D, PRE, EXPF)
cat(sprintf("step2_diff gate (id-limited, content-pinned; table %s): %d expected-change entries, %d problems
", basename(EXPF), nrow(g_rd(EXPF)), length(PROB)))
if (length(PROB)) { for (p in PROB) message("PROBLEM ", p); stop("step2_diff: ", length(PROB), " cell(s) outside the expected-change list or with other content (first: ", PROB[1], ")") }
