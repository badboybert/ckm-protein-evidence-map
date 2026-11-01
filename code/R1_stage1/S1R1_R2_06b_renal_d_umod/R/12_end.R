# 12_end.R -- end-of-run checks: governing md5s unchanged, new inputs unchanged, W6a results.tsv unchanged since read,
# staging snapshot (end) identical to the start snapshot (proves nothing was written under V4_local/staging/l3-finemap/).
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
source(file.path(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod"), "R", "00_common.R"))
say("==== W6b 12_end ====")
assert_governing()
nr <- fread(NEWREG, colClasses = "character")
chg <- nr[vapply(path, md5, "") != md5]
log_assert("END.new_inputs_unchanged", 0, nrow(chg), nrow(chg) == 0, "md5 of every registered new input re-computed at end")
w6a <- file.path(WORK, "w6a_results_read.txt")
if (file.exists(w6a)) {
  r <- readLines(w6a); now <- md5(r[1])
  log_assert("END.w6a_results_unchanged", r[2], now, now == r[2], paste("W6a results.tsv md5 at read vs at end:", r[1]))
}
root <- file.path(ROOT, "V3/V4_local/staging/l3-finemap")
f <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
endt <- data.table(path = sub(paste0("^", root, "/"), "", f), bytes = file.size(f), md5 = unname(tools::md5sum(f)))
setorder(endt, path)
fwrite(endt, file.path(OUT, "staging_snapshot_end.tsv"), sep = "\t")
st <- fread(file.path(OUT, "staging_snapshot_start.tsv"))
same <- nrow(st) == nrow(endt) && all(st$path == endt$path) && all(st$md5 == endt$md5)
log_assert("END.staging_unchanged", sprintf("%d files identical", nrow(st)), sprintf("%d files at end; identical %s", nrow(endt), same), same,
           "md5 of every file under V4_local/staging/l3-finemap/ at start vs end")
say("12_end done")
