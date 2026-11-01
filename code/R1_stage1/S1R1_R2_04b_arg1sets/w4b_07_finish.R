# w4b_07_finish.R -- Stage 1 W4b, step 7: staging snapshot at end (must equal the start), no genotype file in the folder,
# output manifest with md5s.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_common.R"))
logf("==== W4b step 7 (finish) start")
fs <- sort(list.files(file.path(V3, "V4_local/staging/l3-finemap"), recursive = TRUE, full.names = TRUE, all.files = TRUE))
s1 <- data.table(file = sub(paste0(V3, "/"), "", fs, fixed = TRUE), bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
f1 <- file.path(OUT, "staging_snapshot_end.tsv"); fwrite(s1, f1, sep = "\t")
s0 <- fread(file.path(OUT, "staging_snapshot_start.tsv"))
same <- nrow(s0) == nrow(s1) && all(s0$file == s1$file) && all(s0$md5 == s1$md5) && all(s0$bytes == s1$bytes)
logf("staging snapshot (end): %d files, md5 of table %s; identical to start (%d files): %s", nrow(s1), md5(f1), nrow(s0), same)
if (!same) logf("ASSERT_FAILED staging changed between start and end snapshots")
geno <- list.files(OUT, pattern = "\\.(bed|i8|bgen|pgen)$", recursive = TRUE)
logf("genotype-format files inside the W4b folder: %d %s", length(geno), paste(geno, collapse = ","))
if (length(geno)) logf("ASSERT_FAILED genotype file copied into the folder")
outs <- sort(list.files(OUT, recursive = TRUE, full.names = TRUE))
outs <- outs[!grepl("/(boot/boot_[0-9]|susie/[A-Za-z0-9_]+__)", outs)]   # checkpoint units listed by count only
man <- data.table(file = sub(paste0(OUT, "/"), "", outs, fixed = TRUE), bytes = file.size(outs), md5 = unname(tools::md5sum(outs)))
fwrite(man[file != "OUTPUT_MANIFEST.tsv"], file.path(OUT, "OUTPUT_MANIFEST.tsv"), sep = "\t")
logf("checkpoint units: %d SuSiE fits (susie/*.rds), %d bootstrap chunks (boot/boot_*.rds)", length(list.files(file.path(OUT, "susie"), "__.*\\.rds$")),
     length(list.files(file.path(OUT, "boot"), "^boot_[0-9].*\\.rds$")))
logf("==== W4b step 7 done")
