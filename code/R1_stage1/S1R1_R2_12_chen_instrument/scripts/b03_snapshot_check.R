# b03: re-hash every pinned input after all scripts ran and compare with staging_snapshot_start.tsv
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
PK <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_12_chen_instrument")
s <- read.delim(file.path(PK, "staging_snapshot_start.tsv"), stringsAsFactors = FALSE)
now <- unname(tools::md5sum(s$path)); ok <- identical(now, s$md5)
cat(sprintf("[%s] snapshot recheck after b02: %d inputs, %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), nrow(s), if (ok) "ALL md5 UNCHANGED" else "MISMATCH"), file = file.path(PK, "run.log"), append = TRUE)
stopifnot(ok); cat("OK\n")
