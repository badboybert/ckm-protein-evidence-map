# w4b_06c_gate_mutations.R -- mutation tests of the id-limited gate w4b_06b_step2_gate.R (fix round 2026-10-10). Each mutation is applied to a COPY of
# results.tsv / branches.tsv (and of the expected-changes table for the legitimate-change case) in a temp folder; the real gate runs on the copy as a
# separate Rscript process. A mutation is TRIPPED only if the exit status is non-zero AND the output names the mutated id and column (and the kind of
# problem). Controls (unmutated; a legitimate change recorded in the table) are TRIPPED when they exit 0 as required. Writes mutation_tests.tsv.
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
OUTD <- Sys.getenv("W4B_OUT", HERE)   # output folder under test; mutation copies are made from it and mutation_tests.tsv is written to it
GATE <- file.path(HERE, "w4b_06b_step2_gate.R"); EXPT <- file.path(HERE, "step2_expected_changes.tsv"); PREDIR <- file.path(HERE, "fix1_prefix")
M <- list(); rec <- function(id, check, mutation, tripped, detail) {
  M[[length(M) + 1]] <<- data.table(mutation_id = id, check_named = check, mutation = mutation, result = if (tripped) "TRIPPED" else "NOT TRIPPED", detail = detail)
  cat(sprintf("[%s] MUTATION %s -> %s (%s)\n", format(Sys.time(), "%H:%M:%S"), id, if (tripped) "TRIPPED" else "NOT TRIPPED", detail)) }
gate_copy <- function(tag, mut = function(R, B, E) list(R = R, B = B, E = E)) {
  td <- file.path(tempdir(), paste0("gate_", tag)); dir.create(td, showWarnings = FALSE, recursive = TRUE)
  R <- g_rd(file.path(OUTD, "results.tsv")); B <- g_rd(file.path(OUTD, "branches.tsv")); E <- g_rd(EXPT)
  m <- mut(R, B, E)
  fwrite(m$R, file.path(td, "results.tsv"), sep = "\t", quote = TRUE); fwrite(m$B, file.path(td, "branches.tsv"), sep = "\t", quote = TRUE); fwrite(m$E, file.path(td, "exp.tsv"), sep = "\t", quote = TRUE)
  Sys.setenv(W4B_GATE_NEW = td, W4B_GATE_EXPECTED = file.path(td, "exp.tsv"))
  o <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript.exe"), shQuote(GATE), stdout = TRUE, stderr = TRUE))
  Sys.unsetenv(c("W4B_GATE_NEW", "W4B_GATE_EXPECTED"))
  st <- attr(o, "status"); list(rc = if (is.null(st)) 0L else st, out = paste(o, collapse = "\n")) }
only1 <- function(DT, id) { n <- sum(DT$quantity_id == id); stopifnot(n == 1L); which(DT$quantity_id == id) }   # a mutation target must be unique
plus <- function(x) sprintf("%.10g", as.numeric(x) + 1e-3)
CK <- "w4b_06b_step2_gate.R id-limited gate (exit non-zero, message names id and column)"
gate_rec <- function(id, mutation, mut, must_have) {
  g <- gate_copy(id, mut); ok <- g$rc != 0L && all(vapply(must_have, function(x) grepl(x, g$out, fixed = TRUE), TRUE))
  rec(id, CK, mutation, ok, sprintf("exit %d; message must contain: %s", g$rc, paste(must_have, collapse = " | "))) }
ctl <- g_run(OUTD, PREDIR, EXPT); g0 <- gate_copy("G0")
rec("G0", "w4b_06b gate (control: unmutated output must pass)", "none (copy of the real results.tsv/branches.tsv/table)", g0$rc == 0L && length(ctl) == 0L && grepl("0 problems", g0$out, fixed = TRUE),
    sprintf("exit %d; %d problems (a control is TRIPPED when it passes as required)", g0$rc, length(ctl)))
gate_rec("G1", "+1e-3 on the unrelated cell S4.omni.arginine.5sig.aware.Q value", function(R, B, E) { i <- only1(R, "S4.omni.arginine.5sig.aware.Q"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.omni.arginine.5sig.aware.Q column=value", "UNEXPECTED CHANGE"))
gate_rec("G2a", "a pre-existing s_max row changed in a column not on the list (S4.s.max.ARG1_T2D se := 0.001)", function(R, B, E) { i <- only1(R, "S4.s.max.ARG1_T2D"); R[i, se := "0.001"]; list(R = R, B = B, E = E) },
         c("id=S4.s.max.ARG1_T2D column=se", "UNEXPECTED CHANGE"))
gate_rec("G2b", "+1e-3 on a listed s_max (S4.s.max.ARG1_T2D value): listed id, wrong content", function(R, B, E) { i <- only1(R, "S4.s.max.ARG1_T2D"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.s.max.ARG1_T2D column=value", "differs from the recorded expected content"))
gate_rec("G2c", "+1e-3 on the four-fit sensitivity row S4.s.max_4fit.ARG1_T2D (value)", function(R, B, E) { i <- only1(R, "S4.s.max_4fit.ARG1_T2D"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.s.max_4fit.ARG1_T2D column=*", "differs from the recorded expected row"))
gate_rec("G3", "+1e-3 on a new-fit s (S4.s.P__ARG1.susie.489 value)", function(R, B, E) { i <- only1(R, "S4.s.P__ARG1.susie.489"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.s.P__ARG1.susie.489 column=*", "differs from the recorded expected row"))
gate_rec("G3b", "+1e-3 on a pre-existing per-fit s (S4.s.P__ARG1.p1fr.503 value)", function(R, B, E) { i <- only1(R, "S4.s.P__ARG1.p1fr.503"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.s.P__ARG1.p1fr.503 column=value", "UNEXPECTED CHANGE"))
gate_rec("G9a", "a listed genename note edited (S4.genename.ARG1_arginine.glist.total notes + ' x')", function(R, B, E) { i <- only1(R, "S4.genename.ARG1_arginine.glist.total"); R[i, notes := paste0(notes, " x")]; list(R = R, B = B, E = E) },
         c("id=S4.genename.ARG1_arginine.glist.total column=notes", "differs from the recorded expected content"))
gate_rec("G9b", "an unlisted genename cell edited (value of S4.genename.ARG1_arginine.glist.total)", function(R, B, E) { i <- only1(R, "S4.genename.ARG1_arginine.glist.total"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S4.genename.ARG1_arginine.glist.total column=value", "UNEXPECTED CHANGE"))
gate_rec("G5", "an unlisted new row appended", function(R, B, E) { x <- R[1]; x[, quantity_id := "S4.zz.unlisted_new_row"]; list(R = rbind(R, x), B = B, E = E) },
         c("id=S4.zz.unlisted_new_row column=*", "UNEXPECTED NEW ROW"))
gate_rec("G6", "an unrelated branches.tsv cell changed (row 10, qualifiers)", function(R, B, E) { stopifnot(!any(E$file == "branches.tsv" & E$key == sprintf("row 10 (%s)", B$rule_ref[10]))); B[10, qualifiers := paste0(qualifiers, " x")]; list(R = R, B = B, E = E) },
         c("branches.tsv id=row 10", "column=qualifiers", "UNEXPECTED CHANGE"))
gate_rec("G7", "a results.tsv row removed (S4.omni.arginine.5sig.aware.Q)", function(R, B, E) { i <- only1(R, "S4.omni.arginine.5sig.aware.Q"); list(R = R[-i], B = B, E = E) },
         c("id=S4.omni.arginine.5sig.aware.Q column=*", "row REMOVED"))
gate_rec("G8", "a listed change reverted to the frozen value (S4.s.max.ARG1_T2D value): EXPECTED CHANGE ABSENT", function(R, B, E) {
  i <- only1(R, "S4.s.max.ARG1_T2D"); O <- g_rd(file.path(PREDIR, "results.tsv")); R[i, value := O[quantity_id == "S4.s.max.ARG1_T2D"]$value]; list(R = R, B = B, E = E) },
  c("id=S4.s.max.ARG1_T2D column=value", "EXPECTED CHANGE ABSENT"))
legit <- function(R, B, E, upd) { i <- only1(R, "S4.s.max.ARG1_T2D"); R[i, estimand := paste0(estimand, " [legit]")]
  j <- which(E$file == "results.tsv" & E$key == "S4.s.max.ARG1_T2D" & E$column == "estimand"); stopifnot(length(j) == 1L)
  if (upd) E[j, new_md5 := g_md5(R$estimand[i])]; list(R = R, B = B, E = E) }
g4 <- gate_copy("G4", function(R, B, E) legit(R, B, E, TRUE))
rec("G4", "w4b_06b gate (control: a legitimate allowed change recorded in the table must pass)", "S4.s.max.ARG1_T2D estimand text changed AND its table md5 updated", g4$rc == 0L, sprintf("exit %d (must be 0)", g4$rc))
gate_rec("G4b", "the same legitimate-looking change WITHOUT updating the table", function(R, B, E) legit(R, B, E, FALSE), c("id=S4.s.max.ARG1_T2D column=estimand", "differs from the recorded expected content"))
MT <- rbindlist(M); fwrite(MT, file.path(OUTD, "mutation_tests.tsv"), sep = "\t")
cat(sprintf("[%s] W4b gate mutation tests: %d of %d tripped (G0 and G4 are controls: TRIPPED = passed as required)\n", format(Sys.time(), "%H:%M:%S"), sum(MT$result == "TRIPPED"), nrow(MT)))
if (any(MT$result != "TRIPPED")) stop("a gate mutation did not trip the check it names")
