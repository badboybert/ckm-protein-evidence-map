# 13_mutation_tests.R -- every new check must be able to fail: each mutation is applied to an in-memory COPY and must trip the
# check that names it. Writes mutation_tests.tsv. (The PLINK 'every requested SNP returned' mutation runs inside 08_umod.R.)
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
say("==== W6b 13_mutation_tests ====")
M <- list(); rec <- function(id, check, mutation, tripped, detail) {
  M[[length(M) + 1]] <<- data.table(mutation_id = id, check_named = check, mutation = mutation, result = if (tripped) "TRIPPED" else "NOT TRIPPED", detail = detail)
  say("MUTATION %s -> %s (%s)", id, if (tripped) "TRIPPED" else "NOT TRIPPED", detail)
}
# 1. candidate rule: dropping the converged-fit (n_warn == 0) condition must change the count away from 16
cand <- fread(file.path(WORK, "candidates.tsv"))
L31e <- readLines(res("LOG_31_candidate_egfr"), warn = FALSE)
n_with_warned <- nrow(cand) + sum(grepl("ZFYVE19\\s+eGFR\\s+0\\.99", L31e))
rec("M1", "S6.d.candidates (target 16)", "candidate rule without 'n_warn == 0' (keeps ZFYVE19-eGFR)", n_with_warned != 16, sprintf("count %d", n_with_warned))
# 2. kidney classification: CKD not counted as kidney must break the 11/16 vs 57/144 table
pool <- fread(file.path(WORK, "pool160_from_29_headers.tsv"))
a2 <- pool[candidate == TRUE & outcome == "eGFR", .N]; c2 <- pool[candidate == FALSE & outcome == "eGFR", .N]
rec("M2", "S6.d.primary.table", "kidney = eGFR only (CKD dropped)", !(a2 == 11 && c2 == 57), sprintf("%d/16 vs %d/144", a2, c2))
# 3. re-implementation identity check: one beta perturbed in a copy of my UMOD table must fail the 1e-12 identity
eng <- fread(res("LOC_UMOD_eGFR_control_harmonised")); mine <- fread(file.path(WORK, "harm", "UMOD_eGFR.tsv.gz"))
mine2 <- copy(mine); mine2[1234, b_O := b_O * (1 + 1e-6)]
mm <- merge(eng[, .(snp, bO = b__O__eGFR)], mine2[, .(snp, b_O)], by = "snp")
rec("M3", "S6.d.reimpl_vs_engine_UMOD_harmonised", "eGFR beta of SNP row 1234 multiplied by (1 + 1e-6) in a copy", max(abs(mm$bO - mm$b_O)) >= 1e-12, sprintf("max|db| %.3g", max(abs(mm$bO - mm$b_O))))
# 4. allele alignment: a swapped allele pair must flip the sign; a palindromic complement must not be matched
s1 <- align_to_ref("A", "G", "G", "A"); s2 <- align_to_ref("A", "T", "T", "A"); s3 <- align_to_ref("A", "C", "T", "G"); s4 <- align_to_ref("A", "T", "T", "A")
rec("M4", "align_to_ref (engine rule)", "swapped / complemented allele pairs", s1 == -1L && s3 == 1L && s2 == -1L, sprintf("swap %d, palindromic swap %d, complement %d", s1, s2, s3))
rec("M4b", "align_to_ref palindromic guard", "palindromic GWAS pair A/T against a C/G reference (complement match must be refused)", align_to_ref("A", "T", "C", "G") == 0L,
    sprintf("A/T vs C/G -> %d (must be 0)", align_to_ref("A", "T", "C", "G")))
# 5. staging snapshot check: one md5 changed in a copy of the start snapshot must fail identity
st <- fread(file.path(OUT, "staging_snapshot_start.tsv")); st2 <- copy(st); st2[5, md5 := "00000000000000000000000000000000"]
rec("M5", "END.staging_unchanged", "one md5 altered in a copy of the start snapshot", !all(st$md5 == st2$md5), "row 5")
# 6. D6 rule at three decimals: a lower bound of 1.0004 must NOT count as exceeding 1.000
rec("M6", "11 D6 (round(lo, 3) > 1)", "lower bound 1.0004", !(round(1.0004, 3) > 1), sprintf("round(1.0004,3) = %.3f", round(1.0004, 3)))
# 7. 6(d) rule: s exactly 0.10 must reclassify (rule is s >= 0.10)
rec("M7", "6(d) s >= 0.10", "s_max = 0.10", 0.10 >= 0.10, "boundary counts as reclassified")
# 8-9. FIXER F1: the exact-CI checks of 07_fisher_d.R (S6.d.exact_roots, S6.d.matches_fishertest), on the (d) table
ECOR <- file.path(ROOT, "V3/R1_stage1/S1R1_R2_06a_renal/exact_cond_or.R")
assert_md5(ECOR, "1217326d2a6880b174f91d377c5b21b9", "SHARED:exact_cond_or.R (W6a folder; FIXER F1)"); source(ECOR)
fdq <- readRDS(file.path(WORK, "d_summary.rds"))$fd; tab8 <- matrix(c(fdq$a, fdq$c, fdq$b, fdq$d), nrow = 2)
f8 <- fisher.test(tab8); r8 <- exact_cond_or_resid(tab8, unname(f8$estimate), f8$conf.int[1], f8$conf.int[2])
rec("M8", "S6.d.exact_roots (max |residual| < 1e-10)", "MLE and CI taken from fisher.test's default-tolerance uniroot (the original defect)", max(abs(r8)) >= 1e-10,
    sprintf("max |residual| %.2g", max(abs(r8))))
relmax <- function(a, b) { r <- ifelse(is.infinite(a) & is.infinite(b) & sign(a) == sign(b), 0, abs(a - b) / abs(b)); max(r) }   # same helper as 07_fisher_d.R (infinite OR under seven fits)
e9 <- exact_cond_or(tab8, conf.level = 0.90); ft9 <- c(unname(f8$estimate), f8$conf.int); rd9 <- relmax(c(e9$estimate, e9$conf.int), ft9)
rec("M9", "S6.d.matches_fishertest (< 1%)", "interval solved at 90% instead of 95%", rd9 >= 0.01, sprintf("max relative difference %.3g", rd9))
# 10-12. Stage 2 step 2: the +Inf allowance of results.no_nonfinite_values (09_collect.R) is narrow
R6 <- fread(file.path(OUT, "results.tsv"), colClasses = "character", na.strings = NULL)
rec("M10", "results.no_nonfinite_values (allowance)", "unmutated results.tsv: only the declared +Inf values pass", nrow(nonfinite_bad(R6)) == 0L && sum(R6$value == "Inf") == 2L,
    sprintf("%d non-finite flagged; %d Inf values in results.tsv", nrow(nonfinite_bad(R6)), sum(R6$value == "Inf")))
m11 <- copy(R6); m11[quantity_id == "S6.d.p", value := "Inf"]
rec("M11", "results.no_nonfinite_values (allowance)", "+Inf in S6.d.p (not an OR or an upper bound)", nrow(nonfinite_bad(m11)) == 1L, "S6.d.p flagged")
m12 <- copy(R6); m12[quantity_id == "S6.d.b", value := "1"]
rec("M12", "results.no_nonfinite_values (allowance)", "S6.d.b = 1 with S6.d.or still +Inf", nrow(nonfinite_bad(m12)) == 2L, "S6.d.or and S6.d.hi flagged when the table has no zero cell")
# 13-22. Stage 2 step 2 FIX ROUND (2026-10-10): the id-limited gate R/step2_diff.R must fail on a wrong value. Each mutation is applied to a
# COPY of results.tsv / branches.tsv (and, for the legitimate-change case, of the expected-changes table) in a temp folder; the real gate is run
# on the copy as a separate Rscript process; the mutation counts as TRIPPED only if the exit status is non-zero AND the output names the
# mutated id and column (and the kind of problem). Controls (unmutated; a legitimate change recorded in the table) must exit 0.
source(file.path(OUT, "R", "step2_gate_lib.R"))
GATE <- file.path(OUT, "R", "step2_diff.R"); EXPT <- file.path(OUT, "step2_expected_changes.tsv"); PREDIR <- file.path(OUT, "_pre_step2_20261010")
gate_copy <- function(tag, mut = function(R, B, E) list(R = R, B = B, E = E)) {
  td <- file.path(tempdir(), paste0("gate_", tag)); dir.create(td, showWarnings = FALSE, recursive = TRUE)
  R <- g_rd(file.path(OUT, "results.tsv")); B <- g_rd(file.path(OUT, "branches.tsv")); E <- g_rd(EXPT)
  m <- mut(R, B, E)
  fwrite(m$R, file.path(td, "results.tsv"), sep = "\t", quote = TRUE); fwrite(m$B, file.path(td, "branches.tsv"), sep = "\t", quote = TRUE); fwrite(m$E, file.path(td, "exp.tsv"), sep = "\t", quote = TRUE)
  Sys.setenv(W6B_NEW = td, W6B_EXPECTED = file.path(td, "exp.tsv"), W6B_DIFFOUT = file.path(td, "diff.tsv"))   # (system2(env=) is unreliable on Windows)
  o <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript.exe"), shQuote(GATE), stdout = TRUE, stderr = TRUE))
  Sys.unsetenv(c("W6B_NEW", "W6B_EXPECTED", "W6B_DIFFOUT"))
  st <- attr(o, "status"); list(rc = if (is.null(st)) 0L else st, out = paste(o, collapse = "\n"))
}
only1 <- function(DT, id) { n <- sum(DT$quantity_id == id); stopifnot(n == 1L); which(DT$quantity_id == id) }   # a mutation target must be unique
plus <- function(x) sprintf("%.10g", as.numeric(x) + 1e-3)
gate_rec <- function(id, mutation, mut, must_have) {
  g <- gate_copy(id, mut); ok <- g$rc != 0L && all(vapply(must_have, function(x) grepl(x, g$out, fixed = TRUE), TRUE))
  rec(id, "step2_diff.R id-limited gate (exit non-zero, message names id and column)", mutation, ok, sprintf("exit %d; message must contain: %s", g$rc, paste(must_have, collapse = " | ")))
}
ctl <- g_run(OUT, PREDIR, EXPT)
g0 <- gate_copy("G0")
rec("G0", "step2_diff.R gate (control: unmutated output must pass)", "none (copy of the real results.tsv/branches.tsv/table)", g0$rc == 0L && length(ctl) == 0L && grepl("0 problems", g0$out, fixed = TRUE),
    sprintf("exit %d; %d problems (a control is TRIPPED when it passes as required)", g0$rc, length(ctl)))
gate_rec("G1", "+1e-3 on the unrelated cell S6.umod.Q.aware_2nd.p value", function(R, B, E) { i <- only1(R, "S6.umod.Q.aware_2nd.p"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S6.umod.Q.aware_2nd.p column=value", "UNEXPECTED CHANGE"))
gate_rec("G2a", "+1e-3 on an s_max whose value is not on the allowed list (S6.d.ITIH4_eGFR.s_max: only its estimand is listed)", function(R, B, E) { i <- only1(R, "S6.d.ITIH4_eGFR.s_max"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S6.d.ITIH4_eGFR.s_max column=value", "UNEXPECTED CHANGE"))
gate_rec("G2b", "+1e-3 on a listed s_max (S6.d.UMOD_eGFR.s_max): listed id, wrong content", function(R, B, E) { i <- only1(R, "S6.d.UMOD_eGFR.s_max"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S6.d.UMOD_eGFR.s_max column=value", "differs from the recorded expected content"))
gate_rec("G3", "+1e-3 on a per-fit s (new row S6.d.UMOD_eGFR.s.protein.susie.489, value)", function(R, B, E) { i <- only1(R, "S6.d.UMOD_eGFR.s.protein.susie.489"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S6.d.UMOD_eGFR.s.protein.susie.489 column=*", "differs from the recorded expected row"))
gate_rec("G3b", "+1e-3 on a pre-existing per-fit s (S6.d.UMOD_eGFR.s_protein_4fit, new sensitivity row)", function(R, B, E) { i <- only1(R, "S6.d.UMOD_eGFR.s_protein_4fit"); R[i, value := plus(value)]; list(R = R, B = B, E = E) },
         c("id=S6.d.UMOD_eGFR.s_protein_4fit column=*", "differs from the recorded expected row"))
gate_rec("G5", "an unlisted new row appended", function(R, B, E) { x <- R[1]; x[, quantity_id := "S6.zz.unlisted_new_row"]; list(R = rbind(R, x), B = B, E = E) },
         c("id=S6.zz.unlisted_new_row column=*", "UNEXPECTED NEW ROW"))
gate_rec("G6", "an unrelated branches.tsv cell changed (row 20, qualifiers)", function(R, B, E) { stopifnot(!(B$rule_ref[20] %in% character()), !any(E$file == "branches.tsv" & E$key == sprintf("row 20 (%s)", B$rule_ref[20]))); B[20, qualifiers := paste0(qualifiers, " x")]; list(R = R, B = B, E = E) },
         c("branches.tsv id=row 20", "column=qualifiers", "UNEXPECTED CHANGE"))
gate_rec("G7", "a results.tsv row removed (S6.umod.Q.aware_2nd.p)", function(R, B, E) { i <- only1(R, "S6.umod.Q.aware_2nd.p"); list(R = R[-i], B = B, E = E) },
         c("id=S6.umod.Q.aware_2nd.p column=*", "row REMOVED"))
gate_rec("G8", "a listed change reverted to the pre-step-2 text (S6.d.GUSB_eGFR.s_protein value): EXPECTED CHANGE ABSENT", function(R, B, E) {
  i <- only1(R, "S6.d.GUSB_eGFR.s_protein"); O <- g_rd(file.path(PREDIR, "results.tsv")); R[i, value := O[quantity_id == "S6.d.GUSB_eGFR.s_protein"]$value]; list(R = R, B = B, E = E) },
  c("id=S6.d.GUSB_eGFR.s_protein column=value", "EXPECTED CHANGE ABSENT"))
# G4: a LEGITIMATE allowed change (the cell is listed; the new content is recorded in the table) must pass; the same change without the table update must fail
legit <- function(R, B, E, upd) { i <- only1(R, "S6.d.GUSB_eGFR.s_max"); R[i, estimand := paste0(estimand, " [legit]")]
  j <- which(E$file == "results.tsv" & E$key == "S6.d.GUSB_eGFR.s_max" & E$column == "estimand"); stopifnot(length(j) == 1L)
  if (upd) E[j, new_md5 := g_md5(R$estimand[i])]; list(R = R, B = B, E = E) }
g4 <- gate_copy("G4", function(R, B, E) legit(R, B, E, TRUE))
rec("G4", "step2_diff.R gate (control: a legitimate allowed change recorded in the table must pass)", "S6.d.GUSB_eGFR.s_max estimand text changed AND its table md5 updated", g4$rc == 0L, sprintf("exit %d (must be 0)", g4$rc))
gate_rec("G4b", "the same legitimate-looking change WITHOUT updating the table", function(R, B, E) legit(R, B, E, FALSE), c("id=S6.d.GUSB_eGFR.s_max column=estimand", "differs from the recorded expected content"))
MT <- rbindlist(M)
fwrite(MT, file.path(OUT, "mutation_tests.tsv"), sep = "\t")
a <- fread(ASSERTS); mu <- a[assert_id == "MUT.ld_every_pair_returned"]
say("mutation tests: %d of %d tripped; PLINK LD mutation (08_umod): %s", sum(MT$result == "TRIPPED"), nrow(MT), if (nrow(mu)) tail(mu$status, 1) else "not yet run")
