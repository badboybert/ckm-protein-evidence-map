#!/usr/bin/env Rscript
# W7 mutation tests: each run of w7_main.R with W7_MUTATE=<id> injects ONE in-memory defect (no input file is touched) and must trip
# the check that names it. Outputs of mutated runs go to mutations/<id>/. Writes mutation_tests.tsv in this folder.
# Rerun:  Rscript "<archive>/code/R1_stage1/S1R1_R2_07_ancestry/w7_mutations.R"
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
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_07_ancestry")
MAIN <- file.path(OUT, "w7_main.R")
RS <- file.path(R.home("bin"), "Rscript")
tests <- list(
  list(id = "folder_digest", kind = "halt", expect = "CKB_REGIONS_FOLDER digest mismatch", check = "CKB folder digest assert"),
  list(id = "eas_sign", kind = "halt", expect = "n_dir_disagree == 0", check = "direction_class recomputation assert"),
  list(id = "ckbfree_drop", kind = "halt", expect = "nrow(jj) == nrow(ckf)", check = "census x CKB-free join-yield assert"),
  list(id = "dia_drop_arg1", kind = "assert_row", expect = "ASSERT.S8.all28_same_sign", check = "s8 assertion target: 27 same sign (ARG1 becomes not evaluable)"),
  list(id = "dia_flip_arg1", kind = "assert_row", expect = "ASSERT.S8.all28_same_sign", check = "s8 assertion target: 27 same sign (ARG1 check sign flipped)")
)
res <- rbindlist(lapply(tests, function(t) {
  Sys.setenv(W7_MUTATE = t$id)
  out <- suppressWarnings(system2(RS, shQuote(MAIN), stdout = TRUE, stderr = TRUE))
  Sys.unsetenv("W7_MUTATE")
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  if (t$kind == "halt") {
    err <- grep("^Error", out, value = TRUE)
    tripped <- status != 0 && length(err) > 0 && any(grepl(t$expect, err, fixed = TRUE))
    seen <- paste(err, collapse = " | ")
  } else {
    rf <- file.path(OUT, "mutations", t$id, "results.tsv")
    r <- if (file.exists(rf)) fread(rf, sep = "\t", colClasses = "character") else data.table(quantity_id = character(), notes = character())
    row <- r[quantity_id == t$expect]
    tripped <- status == 0 && nrow(row) == 1 && startsWith(row$notes, "ASSERT_FAILED")
    seen <- if (nrow(row)) paste0(row$value, " | ", row$notes) else "row missing"
  }
  cat(sprintf("%-15s %-10s exit %s tripped %s | %s\n", t$id, t$kind, status, tripped, substr(seen, 1, 160)))
  data.table(mutation = t$id, check_named = t$check, kind = t$kind, exit_status = status, tripped = tripped, evidence = seen)
}))
fwrite(res, file.path(OUT, "mutation_tests.tsv"), sep = "\t")
cat("mutation tests tripped:", sum(res$tripped), "of", nrow(res), "\n")
if (!all(res$tripped)) quit(status = 1)
