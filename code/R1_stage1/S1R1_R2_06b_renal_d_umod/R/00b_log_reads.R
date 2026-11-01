# 00b_log_reads.R -- records in run.log every seen script opened (with the reason) and every engine script read (rule: seen scripts).
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
e <- res_entry("R1_arg1_triangle_R")
stopifnot(md5(e$path) == e$md5)
say("SEEN SCRIPT OPENED (method question): %s = %s md5 %s (role %s).", e$id, e$path, e$md5, e$role)
say("  Reason: Rule 5 fixes the SE order as that of the R1 implementation of the original rule (delta method, second order); neither the redacted pre-declaration nor R1_ARG1_rule (definition) gives the formula.")
say("  Read only lines 55-84 and 113-116 (single-SNP conditioning at(), ratio(), ratio_se1(), het()); no R1 output file was opened.")
say("  Taken: ratio SE = abs(k) * sqrt((se_num/b_num)^2 + (se_den/b_den)^2) (2nd order); 1st order = se_num/abs(b_den); R1 het() ignores the A-B covariance (= p_ign here).")
for (id in c("ENG_l3fm__prep.py", "ENG_l3fm__panel2.py", "ENG_l3fm__susie.R", "ENG_l3fm__conditional.py", "ENG_l3fm__common.py")) {
  x <- res_entry(id); stopifnot(md5(x$path) == x$md5)
  say("ENGINE SCRIPT READ (role %s; read to understand, not run): %s md5 %s", x$role, x$path, x$md5)
}
rb <- file.path(ROOT, "V3/V4_local/builders/l3fm__robustness.py")
say("ENGINE SCRIPT READ (docstring only, p1fr rule): %s md5 %s (RESOURCES ENG_l3fm__robustness.py)", rb, md5(rb))
say("DEFINITIONS READ: R1_ARG1_rule %s; DEF_SM_CANDIDATES; LOG_28; UMOD/DPEP1 PREP_STAMP, UMOD PANEL2_STAMP, UMOD COND_STAMP (definitions)", res_entry("R1_ARG1_rule")$md5)
