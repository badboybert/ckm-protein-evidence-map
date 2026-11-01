# 01_setup.R -- governing md5s, staging snapshot (start), gls_q.R self-test, tool versions.
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
say("==== W6b 01_setup ====")
assert_governing()

# staging snapshot (rule 2): md5 of every file under V4_local/staging/l3-finemap/
snap <- function(tag) {
  root <- file.path(ROOT, "V3/V4_local/staging/l3-finemap")
  f <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  dt <- data.table(path = sub(paste0("^", root, "/"), "", f), bytes = file.size(f), md5 = unname(tools::md5sum(f)))
  setorder(dt, path)
  out <- file.path(OUT, sprintf("staging_snapshot_%s.tsv", tag))
  fwrite(dt, out, sep = "\t")
  say("staging snapshot %s: %d files, md5 of listing %s", tag, nrow(dt), md5(out))
  dt
}
if (!file.exists(file.path(OUT, "staging_snapshot_start.tsv"))) snap("start") else say("staging snapshot start exists (kept from first run)")

# tool versions
say("R: %s", R.version.string)
for (p in c("data.table", "susieR", "coloc", "Rsamtools", "rtracklayer", "jsonlite")) say("pkg %s %s", p, as.character(packageVersion(p)))
pv <- system2(plink_exe(), "--version", stdout = TRUE)
say("PLINK19_paper6 --version: %s", paste(pv, collapse = " "))
say("RNG: %s; seeds: 20261007 for every bootstrap (rule 8)", paste(RNGkind(), collapse = "/"))

# shared GLS estimator self-test (rule 6)
st <- system2(file.path(R.home("bin"), "Rscript.exe"), c(q(GOV$gls_q[1]), "--self-test"), stdout = TRUE, stderr = TRUE)
for (l in st) say("gls_q.R --self-test | %s", l)
if (!any(grepl("SELF-TEST PASS", st))) stop("gls_q.R self-test did not pass")
say("01_setup done")
