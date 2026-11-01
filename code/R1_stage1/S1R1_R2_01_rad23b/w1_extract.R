# W1 extraction runner. Usage: Rscript w1_extract.R targets | <JOB> | all
# 'targets' builds extracts/targets.tsv (PLINK --make-just-bim + liftover); a JOB streams one large file
# into a small checkpointed extract (resumable: a finished extract is reused). w1_main.R calls the same
# functions, so running this first only pre-fills the checkpoints (it was run in parallel to save time).
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
args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_01_rad23b")
job <- if (length(args)) args[1] else "all"
Sys.setenv(W1_LOG = file.path(OUTDIR, "extracts", paste0("extract_", job, ".log")))
source(file.path(OUTDIR, "w1_common.R"))
cat("", file = LOG)
assert_governing(); load_resources()
if (job == "targets") { staging_snapshot("pre_extract"); build_targets(); quit(save = "no") }
tg <- fread(file = file.path(OUT, "extracts", "targets.tsv"), colClasses = c(chr = "character"))
jobs <- extraction_jobs(tg)
todo <- if (job == "all") names(jobs) else job
for (j in todo) run_job(j, jobs)
logf("w1_extract.R %s finished", job)
