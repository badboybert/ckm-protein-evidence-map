# w4b_run_all.R -- reruns every W4b step in order, each in its own R process.
#   Rscript w4b_run_all.R            resume (SuSiE fits and bootstrap chunks already on disk are reused)
#   Rscript w4b_run_all.R --fresh    delete the checkpoints and data/ derived files first, then recompute everything
# To re-run into another folder (verifier L3), set W4B_OUT to that folder; inputs are read in place, md5-asserted.
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
here <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets")
out <- Sys.getenv("W4B_OUT", here)
if ("--fresh" %in% commandArgs(trailingOnly = TRUE)) {
  for (d in c("susie", "boot", "data", "parts", "plink")) unlink(file.path(out, d), recursive = TRUE)
  cat("removed checkpoints and derived data under", out, "\n")
}
steps <- c("w4b_00_setup.R", "w4b_01_data.R", "w4b_02_susie.R", "w4b_02b_fitset7.R", "w4b_03_omnibus.R", "w4b_04_credsets.R", "w4b_05_rules.R",
           "w4b_06_compare_seen.R", "w4b_06b_step2_gate.R", "w4b_06c_gate_mutations.R", "w4b_07_finish.R")
for (s in steps) {
  t0 <- Sys.time()
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(here, s)))
  cat(sprintf("%s exit %s in %.0fs\n", s, st, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (st != 0) stop("step failed: ", s)
}
