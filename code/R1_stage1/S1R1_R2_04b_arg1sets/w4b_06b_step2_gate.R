# w4b_06b_step2_gate.R -- ID-LIMITED, content-pinned frozen-copy gate for W4b (fix round 2026-10-10). Compares results.tsv / branches.tsv with
# the Stage 1 FREEZE copy fix1_prefix/ (md5 of its results.tsv asserted here). The ONLY differences allowed are the cells and added rows listed
# in step2_expected_changes.tsv (FIXLOG F1 + Stage 2 step 2 seven-fit; one reason each), each with exactly the recorded content (md5). s_max
# (S4.s.max.ARG1_*), the new-fit s rows (S4.s.<trait>.<p2np.503|susie.489|p1fr.489>), the four-fit sensitivity rows and the genename /
# credible-set notes are all covered. Anything else is printed as PROBLEM <file> id=.. column=.. and the script exits non-zero (stop).
# W4B_OUT (as in w4b_common.R) selects the output folder that is checked and whose run.log is appended; default = this folder.
# Test overrides (mutation tests run on COPIES): W4B_GATE_NEW (dir with results.tsv/branches.tsv), W4B_GATE_EXPECTED (table); with an
# override nothing is written to run.log.
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
HERE <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets")
source(file.path(HERE, "w4b_step2_gate_lib.R"))
OUTD <- Sys.getenv("W4B_OUT", HERE)   # output folder under test (same convention as w4b_common.R); scripts, expected table and fix1_prefix are inputs read from HERE
NEW <- Sys.getenv("W4B_GATE_NEW", OUTD); EXPF <- Sys.getenv("W4B_GATE_EXPECTED", file.path(HERE, "step2_expected_changes.tsv")); PRE <- file.path(HERE, "fix1_prefix")
real <- !nzchar(Sys.getenv("W4B_GATE_NEW"))
say <- function(...) { l <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...)); cat(l, "\n", sep = ""); if (real) cat(l, "\n", sep = "", file = file.path(OUTD, "run.log"), append = TRUE) }
want <- "72c729526796514c5d975856989cab4f"; got <- unname(tools::md5sum(file.path(PRE, "results.tsv")))
if (!identical(got, want)) stop("fix1_prefix/results.tsv md5 ", got, " != ", want)
PROB <- g_run(NEW, PRE, EXPF)
say("STEP2 GATE (id-limited, content-pinned; table %s): %d expected-change entries; %d problems", basename(EXPF), nrow(g_rd(EXPF)), length(PROB))
if (length(PROB)) { for (p in PROB) message("PROBLEM ", p); stop("w4b_06b: ", length(PROB), " cell(s) outside the expected-change list or with other content (first: ", PROB[1], ")") }
