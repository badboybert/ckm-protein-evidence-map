# run_all.R -- reruns every W6b step in order. Every step is resumable (finished units are skipped).
#   Rscript "<archive>/code/R1_stage1/S1R1_R2_06b_renal_d_umod/R/run_all.R"
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
D <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06b_renal_d_umod/R")
rs <- file.path(R.home("bin"), "Rscript.exe")
step <- function(f, args = character(0)) {
  cat(sprintf("[run_all] %s %s\n", f, paste(args, collapse = " ")))
  st <- system2(rs, c(shQuote(file.path(D, f), type = "cmd"), args))
  if (st != 0) stop("step failed: ", f)
}
step("01_setup.R"); step("00b_log_reads.R"); step("02_candidates.R"); step("03_windows.R"); step("04_harmonise.R")
step("05_susie.R", "4")                    # 4 parallel workers; checkpoints per candidate x set
step("05b_fitset7.R", "4")                # Stage 2 step 2: the 3 extra fits (p2np/503, susie/489, p1fr/489) for the seven-fit s_max; checkpoint per candidate
step("07_fisher_d.R"); step("08_umod.R"); step("09_collect.R")
step("10_branches.R")                      # reads W6a results.tsv (md5 recorded)
# Stage 2 step 2: 10_branches.R writes value-bearing `inputs`; the Stage 1b key-id normaliser is re-applied (10b), then the scope report
# replaces the FIXER F1 scope gate (F1 compared with fix_F1_prefix/, whose rows differ structurally from the seven-fit output; the
# F1 gate and its 4/4 mutations are kept in W6a/fix_F1_diffgate.R and fix_F1_diffgate*.tsv from the Stage 1 run).
py <- "python"
st <- system2(py, c("-X", "utf8", "-I", shQuote(file.path(D, "10b_normalise_inputs.py"), type = "cmd")))
if (st != 0) stop("step failed: 10b_normalise_inputs.py (exit ", st, ")")
step("step2_diff.R")
step("11_compare_seen.R")                 # opens the full pre-declaration and the seen outputs ONLY after results/branches exist
step("13_mutation_tests.R"); step("14_join_yields.R")
step("12_end.R")
