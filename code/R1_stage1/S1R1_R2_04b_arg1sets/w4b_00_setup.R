# w4b_00_setup.R -- Stage 1 W4b, step 0: governing md5s, tool versions, gls_q.R self-test, staging snapshot (start).
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
logf("==== W4b step 0 (setup) start; R %s; wd %s", R.version.string, getwd())
for (k in names(GOV)) assert_md5(GOV[[k]]$path, GOV[[k]]$md5, k)
logf("NOTE: the full pre-declaration texts are md5-asserted here but NOT opened before results.tsv is written (rule 10)")

# tool versions
pv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
for (p in c("data.table", "susieR", "coloc", "Rsamtools", "readxl", "jsonlite")) logf("package %s %s", p, pv(p))
if (pv("coloc") != "5.2.3") logf("WARNING coloc version %s != 5.2.3", pv("coloc"))
dir.create(file.path(OUT, "plink"), showWarnings = FALSE)
v <- run_plink("--version", "version")
logf("PLINK --version: %s", paste(v, collapse = " "))

# gls_q.R self-test (rule 6): run the shared file as a script
gq <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(GOV$GLS_Q_R$path), "--self-test"), stdout = TRUE, stderr = TRUE)
for (l in gq) logf("gls_q.R --self-test | %s", l)
if (!any(grepl("SELF-TEST PASS", gq))) stop("gls_q.R self-test did not pass")

# staging snapshot (rule 2): md5 of every file under V4_local/staging/l3-finemap/
snap <- function(tag) {
  fs <- sort(list.files(file.path(V3, "V4_local/staging/l3-finemap"), recursive = TRUE, full.names = TRUE, all.files = TRUE))
  dt <- data.table(file = sub(paste0(V3, "/"), "", fs, fixed = TRUE), bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
  f <- file.path(OUT, paste0("staging_snapshot_", tag, ".tsv"))
  fwrite(dt, f, sep = "\t")
  logf("staging snapshot (%s): %d files, %.1f MB; md5 of the snapshot table %s -> %s", tag, nrow(dt), sum(dt$bytes) / 1e6,
       md5(f), f)
  dt
}
s0 <- snap("start")
logf("==== W4b step 0 done")
