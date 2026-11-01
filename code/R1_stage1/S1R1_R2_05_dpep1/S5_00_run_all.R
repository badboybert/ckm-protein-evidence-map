# S5_00_run_all.R -- the one command that reruns W5 (§5 DPEP1) end to end:
#   "Rscript" "<archive>/code/R1_stage1/S1R1_R2_05_dpep1/S5_00_run_all.R"
# Each step runs in its own Rscript process (no environment bleed). The bootstrap (S5_03) is checkpointed: a killed run
# resumes from work/boot_draws.tsv. Delete work/ and ld/ for a run from scratch. Never edit a step while it is running
# (Rscript streams the file; an edit mid-run desynchronises the parser -- this happened once in development, see README).
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_05_dpep1/S5_common.R"))
rs_exe <- file.path(R.home("bin"), "Rscript.exe")
t0 <- Sys.time()
say("==== S5_00_run_all start ====")
steps <- c("S5_01_prepare.R", "S5_02_analyses.R", "S5_03_bootstrap.R", "S5_04_assemble.R", "S5_05_selfcheck.R")
ca <- commandArgs(trailingOnly = TRUE)
if (length(ca) >= 2 && ca[1] == "--from") {           # resume from a step whose predecessors' outputs exist in work/
  stopifnot(ca[2] %in% steps); steps <- steps[match(ca[2], steps):length(steps)]
  say("resume: running %s only", paste(steps, collapse = ", "))
}
for (st in steps) {
  t1 <- Sys.time()
  out <- system2(rs_exe, q(file.path(OUT, st)), stdout = file.path(WORK, paste0(sub("\\.R$", "", st), ".stdout.txt")), stderr = file.path(WORK, paste0(sub("\\.R$", "", st), ".stderr.txt")))
  say("step %s exit %s in %.0f s", st, out, as.numeric(difftime(Sys.time(), t1, units = "secs")))
  if (!identical(as.integer(out), 0L)) stop("step failed: ", st, " (see work/*.stderr.txt)")
}
snap_staging("end")
a <- fread(file.path(WORK, "staging_md5_start.tsv")); b <- fread(file.path(WORK, "staging_md5_end.tsv"))
same <- identical(a$file, b$file) && identical(a$md5, b$md5)
say("STAGING start vs end: %d vs %d files; identical file list and md5s: %s", nrow(a), nrow(b), same)
if (!same) stop("V4_local/staging/l3-finemap changed during the run")
say("==== S5_00_run_all done in %.0f s ====", as.numeric(difftime(Sys.time(), t0, units = "secs")))
