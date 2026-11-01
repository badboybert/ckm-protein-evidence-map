# W6a -- the one command that reruns everything:
#   "Rscript" "<archive>/code/R1_stage1/S1R1_R2_06a_renal/w6a_run_all.R"
# Order: (1) w6a_compute.R -> results.tsv, branches.tsv (+ intermediates, bootstrap checkpoints, asserts)
#        (2) w6a_mutation_tests.R -> every named assertion trips on its mutation
#        (2b) fix_F1_diffgate.R (FIXER F1) -> results/branches differ from fix_F1_prefix/ only in the Fisher OR/CI cells; its own mutations trip
#        (3) w6a_compare.R -> COMPARE_SEEN.tsv (asserts results/branches equal the post-fix md5 pins first, and the prefix copies the original pins)
# To rerun into another folder (e.g. a verifier's verify_L3/rerun):  set W6A_OUT=<folder> (must exist) and W6A_SRC=<this folder>.
# Delete <OUT>/boot_ckpt to recompute the bootstrap instead of resuming from its per-seed checkpoints.
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
SRC <- Sys.getenv("W6A_SRC", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))
OUT <- Sys.getenv("W6A_OUT", SRC)
rs <- file.path(R.home("bin"), "Rscript.exe")
Sys.setenv(W6A_OUT = OUT, W6A_SRC = SRC)   # inherited by the child Rscript processes
cat(sprintf("\n[%s] w6a_run_all.R: SRC=%s OUT=%s ; PLINK not used (W6a computes no LD quantity: §0 rules 3-4 not applicable)\n",
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), SRC, OUT), file = file.path(OUT, "run.log"), append = TRUE)
F1_RE <- "^S6\\.(primary|b\\.pair|a|b\\.a_pair|b\\.locus_outcome|G[1-4])\\.or_(mle|lo|hi)$"   # FIXER F1 scope (FIXLOG.md)
for (f in c("w6a_compute.R", "w6a_mutation_tests.R", "fix_F1_diffgate.R", "w6a_compare.R")) {
  t0 <- Sys.time()
  rc <- if (f == "fix_F1_diffgate.R") {   # FIXER F1: the regenerated results/branches may differ from fix_F1_prefix/ only in the F1 cells
    Sys.setenv(FIX_F1_MUTATE = "1")
    system2(rs, c(shQuote(file.path(SRC, f)), shQuote(OUT), shQuote(file.path(SRC, "fix_F1_prefix")), shQuote(F1_RE)), stdout = file.path(OUT, "fix_F1_diffgate.out"), stderr = file.path(OUT, "fix_F1_diffgate.out"))
  } else system2(rs, shQuote(file.path(SRC, f)))
  cat(sprintf("[%s] w6a_run_all.R: %s exit %d in %.1f s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), f, rc, as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      file = file.path(OUT, "run.log"), append = TRUE)
  if (rc != 0) stop(f, " failed with exit ", rc)
}
# final staging snapshot: nothing under V4_local/staging/l3-finemap changed during the whole run (start snapshot written by w6a_compute.R)
ST <- paste0(.P31R_V3, "/V4_local/staging/l3-finemap")
f <- list.files(ST, recursive = TRUE, full.names = TRUE, all.files = TRUE)
now <- data.frame(file = sub(paste0("^", ST, "/"), "", f), md5 = unname(tools::md5sum(f)))
st0 <- utils::read.delim(file.path(OUT, "staging_snapshot_start.tsv"))
mm <- merge(st0, now, by = "file", all = TRUE)   # join by file name: base order() and data.table's C-locale order differ
same <- nrow(st0) == nrow(now) && nrow(mm) == nrow(now) && !anyNA(mm$md5.x) && !anyNA(mm$md5.y) && all(mm$md5.x == mm$md5.y)
cat(sprintf("[%s] w6a_run_all.R: FINAL staging snapshot %d files; identical to the compute-start snapshot = %s
", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), nrow(now), same),
    file = file.path(OUT, "run.log"), append = TRUE)
if (!same) stop("staging changed")
