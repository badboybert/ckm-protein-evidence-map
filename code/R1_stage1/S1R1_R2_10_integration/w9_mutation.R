# W9 mutation tests: every new W9 check is fed a mutated input and must trip (brief section 0 rule 3 spirit; global CLAUDE.md s9).
# Writes mutation_tests.tsv in the W9 folder; scratch files go to mutations/ inside the folder.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_10_integration/w9_common.R"))
LOGF <- file.path(OUT, Sys.getenv("W9_LOG", "run.log"))
logline("================ W9 w9_mutation.R start ================")
MD <- file.path(OUT, "mutations"); dir.create(MD, showWarnings = FALSE)
MT <- list()
trips <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")
rec <- function(id, check, mutation, tripped, detail) {
  MT[[length(MT) + 1L]] <<- data.table(mutation_id = id, check = check, mutation = mutation, tripped = tripped, detail = detail)
  logline(sprintf("mutation %-4s %-42s tripped=%s  | %s", id, check, tripped, detail))
}

# M1 md5 assertion
f <- file.path(S1, "S1R1_R2_04a_arg1grid/results.tsv")
rec("M1", "assert_md5", "expected md5 with its last hex digit changed", trips(assert_md5(f, "f4c5104dce1f1822448a6a7bdced1630", "mut")), "W4a results.tsv")

# M2 lookups: missing id, duplicated id, placeholder token, non-numeric
dt <- data.table(quantity_id = c("a", "b", "b", "c", "d"), value = c("0.5", "0.1", "0.2", "TBD", "abc"))
rec("M2a", "getq (join yield 1/1)", "quantity id absent", trips(getq(dt, "mut", "zz")), "0 rows")
rec("M2b", "getq (join yield 1/1)", "quantity id duplicated", trips(getq(dt, "mut", "b")), "2 rows")
rec("M2c", "getq (placeholder token)", "value 'TBD'", trips(getq(dt, "mut", "c")), "TBD")
rec("M2d", "getq (finite number)", "value 'abc'", trips(getq(dt, "mut", "d")), "abc")
rec("M2e", "getq control (must NOT trip)", "none", trips(getq(dt, "mut", "a")), "a = 0.5 returns cleanly; tripped must be FALSE")

# M3 reader: empty table
ef <- file.path(MD, "empty.tsv"); writeLines("quantity_id\tvalue", ef)
rec("M3", "read_tsv_chk (non-empty)", "header-only file", trips(read_tsv_chk(ef, "mut-empty")), "0 rows")
# M3b reader control: a well-formed 2-row file must read cleanly (the check must not trip on good input)
bf <- file.path(MD, "rowcount.tsv"); writeLines(c("quantity_id\tvalue", "x\t1", "y\t2"), bf)
rec("M3b", "read_tsv_chk control (must NOT trip)", "none", trips(read_tsv_chk(bf, "mut-ok")), "2 rows; tripped must be FALSE")

# M4 Holm + nominal: the observed family with ARG1-T2D moved from 0.0365 to 0.004
p <- c(0.9456049, 0.03651846, 0.06200886, 8.088989e-05, 3.797569e-12, 0.97966496, 0.11063492, 0.4472125938)
h0 <- holm(p); h1 <- holm(replace(p, 2, 0.004))
c0 <- nominal_literal(p[2], h0[2])$code; c1 <- nominal_literal(0.004, h1[2])$code
rec("M4", "nominal (literal, default) follows Holm", "ARG1-T2D p 0.0365 -> 0.004", c0 != c1, sprintf("observed %s (Holm %.4f) -> mutated %s (Holm %.4f)", c0, h0[2], c1, h1[2]))
c0i <- nominal_I1(p[2], h0[2])$code; c1i <- nominal_I1(0.004, h1[2])$code
rec("M4c", "nominal (applied I1) follows Holm", "ARG1-T2D p 0.0365 -> 0.004", c0i != c1i, sprintf("observed %s -> mutated %s", c0i, c1i))
rec("M4b", "holm() p.adjust vs manual", "manual implementation off by one rank", {
  bad <- function(p) { m <- length(p); o <- order(p); a <- numeric(m); a[o] <- pmin(1, cummax((m - seq_len(m)) * p[o])); a }
  max(abs(p.adjust(p, "holm") - bad(p))) > 1e-12 }, "a rank-shifted manual Holm disagrees with p.adjust, so the cross-check can fail")

# M5 non-rejection
rec("M5", "nominal I1 (non-rejection)", "p 0.03 -> 0.06 at Holm 0.31", nominal_I1(0.03, 0.31)$code != nominal_I1(0.06, 0.31)$code,
    sprintf("%s -> %s", nominal_I1(0.03, 0.31)$code, nominal_I1(0.06, 0.31)$code))
rec("M5b", "nominal literal (Holm side of 0.05)", "Holm p 0.31 -> 0.04", nominal_literal(0.03, 0.31)$attach != nominal_literal(0.03, 0.04)$attach,
    sprintf("attach %s -> %s", nominal_literal(0.03, 0.31)$attach, nominal_literal(0.03, 0.04)$attach))

# M6 section 3 branch: PAM PP.H3 at 1e-6 0.2437 -> 0.6
b0 <- s3_branch(0.968787, 0.0312126, 0.756325, 0.243675); b1 <- s3_branch(0.968787, 0.0312126, 0.756325, 0.6)
rec("M6", "s3_branch == W3 branch", "PAM PP.H3(1e-6) 0.244 -> 0.6", b0 != b1, sprintf("%s -> %s", b0, b1))
rec("M6b", "corroborated (condition 2)", "PAM branch with the add-on", corroborated(TRUE, b0)$value != corroborated(TRUE, b1)$value,
    sprintf("%s -> %s", corroborated(TRUE, b0)$value, corroborated(TRUE, b1)$value))
rec("M6c", "s3_branch (branch 3)", "RAD23B PP.H3(1e-5) 0.293 -> 0.6", s3_branch(0.0167325, 0.292674, 0.00169884, 0.297149) != s3_branch(0.0167325, 0.6, 0.00169884, 0.297149),
    sprintf("%s -> %s", s3_branch(0.0167325, 0.292674, 0.00169884, 0.297149), s3_branch(0.0167325, 0.6, 0.00169884, 0.297149)))

# M7 corroboration conditions 1 and 3
rec("M7a", "corroborated (condition 1)", "check not concordant", corroborated(TRUE, NA)$value != corroborated(FALSE, NA)$value, "ARG1-type pair, no section 3")
rec("M7b", "corroborated (condition 3)", "TGFB1 section 2 (a) -> (b)", corroborated(TRUE, "1", "(a)", TRUE)$value != corroborated(TRUE, "1", "(b)", TRUE)$value,
    sprintf("%s -> %s", corroborated(TRUE, "1", "(a)", TRUE)$value, corroborated(TRUE, "1", "(b)", TRUE)$value))

# M8 section 2 branch: flip the sign of FinnGen rho_B
s0 <- s2_branch(0.245, -0.332, 3.8e-12, 0.395, -0.380, 2.4e-13, 1); s1 <- s2_branch(0.245, 0.332, 3.8e-12, 0.395, -0.380, 2.4e-13, 1)
rec("M8", "s2_branch == W2 branch", "FinnGen rho_B sign flipped", s0 != s1, sprintf("%s -> %s", s0, s1))

# M9 A1 detection: a colocalization clause injected into the T2D qualifiers must take the front slot; 'does not fire' must not
q <- "nominal [W9]; sensitive to the reference LD; the sign of the correlation between the signals is not identified"
d0 <- detect_caveats("4 T2D", "primary p < 0.05", "signal B's T2D association fell short", q, TRUE)
d1 <- detect_caveats("4 T2D", "primary p < 0.05", "signal B's T2D association fell short",
                     paste(q, "| the T2D association at signal A did not colocalize with the protein (PP.H4 0.5), so the ratio there is not anchored to the protein signal"), TRUE)
rec("M9", "A1 front caveat order", "colocalization clause injected", d0$front != d1$front, sprintf("%s -> %s", d0$front, d1$front))
d2 <- detect_caveats("Rule 3(iv) x", "does not fire", "", "did not colocalize with the protein", FALSE)
rec("M9b", "A1 COLOC guard ('does not fire')", "branch reads 'does not fire'", !("COLOC" %in% d2$caveats), "COLOC not detected when the clause does not fire")
d3 <- detect_caveats("4 T2D", "x", "y", "sensitive to the reference LD", FALSE); d4 <- detect_caveats("4 T2D", "x", "y", "sensitive to the reference LD", TRUE)
rec("M9c", "A1 order (LD before nominal)", "nominal added", identical(d3$front, d4$front) && identical(d4$caveats, c("LDREF", "NOMINAL")),
    sprintf("front stays %s; caveats %s", d4$front, paste(d4$caveats, collapse = " > ")))

# M10 placeholder residual check
qr <- c("no 'nominal' [W9 ...]", "nominal: PENDING_HOLM (W9 resolves)")
rec("M10", "no_unresolved_placeholder", "one PENDING_HOLM left in place", any(grepl("PENDING_HOLM|left to W9", qr)), "check expression is TRUE on the mutated vector, so the assert would fail")

# M11 first-mention flags: Rule 2 only when the class differs
f0 <- firstmention_flags("'proxy' label; class under second-order SE: Bonferroni (first-order Bonferroni)")
f1 <- firstmention_flags("'proxy' label; class under second-order SE: FDR (first-order Bonferroni)")
rec("M11", "firstmention_flags (Rule 2)", "second-order class Bonferroni -> FDR", !identical(f0, f1), sprintf("{%s} -> {%s}", paste(f0, collapse = ","), paste(f1, collapse = ",")))

# M12 renal D6 / list (renal_d6_list in w9_common.R, the one implementation): (d) exact lower bound 0.140 -> 1.2
rd0 <- renal_d6_list(1.002789613, 0.1399375022, 164.4, 0.863, 12.48, 1.028, 12.46, 0.874, 11.70)
rd1 <- renal_d6_list(1.002789613, 1.2, 164.4, 0.863, 12.48, 1.028, 12.46, 0.874, 11.70)
rec("M12", "D6 == W6b", "(d) lower 0.140 -> 1.2", rd0$d6 != rd1$d6, sprintf("%s -> %s", rd0$d6, rd1$d6))
rec("M12b", "renal list ('the LD guard')", "(d) CI 0.140-164.4 -> 1.2-164.4", ("the LD guard" %in% rd0$list) != ("the LD guard" %in% rd1$list),
    sprintf("{%s} -> {%s}", paste(rd0$list, collapse = ", "), paste(rd1$list, collapse = ", ")))

# M13 section 12(a) between-resource Q: flip the sign of the DIAMANTE marginal difference
qf <- function(d, s) { w <- 1 / s^2; e <- sum(w * d) / sum(w); pchisq(sum(w * (d - e)^2), 1, lower.tail = FALSE) }
q0 <- qf(c(-0.1546862, -0.1532843), c(0.07870827, 0.06876239)); q1 <- qf(c(-0.1546862, 0.1532843), c(0.07870827, 0.06876239))
rec("M13", "s12(a) Cochran Q p", "DIAMANTE marginal difference sign flipped", (q0 >= 0.05) != (q1 >= 0.05), sprintf("Q p %.3f -> %.4f", q0, q1))

# ---------------------------------------------------------------- W9F (2026-10-09): every assertion W9F added
# Each mutation below runs w9_main.R end to end with W9_MUTATE=<id>: the real inputs are read and md5-asserted, ONE cell
# of one loaded table (or of an intermediate table) is changed in memory, and the run writes its asserts.tsv to
# mutations/run_<id>/. The mutation passes only if the assert it NAMES is ASSERT_FAILED there (other failures are listed).
# Q1 is a unit test of quote_ok(), which stops rather than records.
main_asserts <- fread(file.path(OUT, "asserts.tsv"), sep = "\t", colClasses = "character")
rec("W0", "control: the unmutated main run (must NOT trip)", "none", any(main_asserts$status != "PASS"),
    sprintf("%d asserts in asserts.tsv, %d ASSERT_FAILED; tripped must be FALSE", nrow(main_asserts), sum(main_asserts$status != "PASS")))
W9F_MUT <- list(
  c1_status    = c("S9.check.C1_class_verified.RAD23B_T2D", "C1_SURVIVAL: C1H.RAD23B_T2D.class status SURVIVES -> L1L2_AGREE_IMPL_DIFFERS"),
  c1_L2        = c("S9.check.C1_class_verified.RAD23B_T2D", "C1_SURVIVAL: C1H.RAD23B_T2D.class L2 concordant -> inconclusive"),
  c1_value     = c("S9.check.C1_class_verified.ARG1_T2D", "C1 results.tsv: C1.ARG1_T2D.class concordant -> inconclusive"),
  c1_vs_w1     = c("S9.check.C1_class_agrees_R1_and_W1.PAM_T2D", "W1 ext_checks.tsv: PAM class_resource_p concordant -> discordant"),
  rad_interval = c("S9.check.RAD23B_clause_fires_on_adjudicator_numbers", "adj_C1.tsv: RAD23B first-order rounding interval 0.0361-0.05086 -> 0.0361-0.0490"),
  rad_bonf     = c("S9.check.RAD23B_bonf6_agrees", "adj_C1.tsv: RAD23B bonf6_W False -> True"),
  a2_share     = c("S9.check.A2_failing_resources_rederived", "W4a results: FinnGen conditional bootstrap share 0.8384 -> 0.95"),
  a2_frag      = c("S9.check.A2_caveat_matches_W4a_fragment", "W4a branches: A2 row caveat 'in FinnGen' -> 'in DIAMANTE'"),
  cond2        = c("S9.check.cond2_invariant_to_S3_conventions.PAM_T2D", "adj_S3_coloc_pkg.tsv: PAM B_full_exposure branch 2 -> 2_addon"),
  n2_files     = c("S9.check.s6d_notch2_two_files_agree", "a3_notch2_fits.tsv: keep_shift_rows p1fr@489 sP -> 0.2"),
  w6b_smax     = c("S9.check.s6d_seven_fit_primary_is_W6b", "W6b results: S6.d.GSTA1_eGFR.s_max (seven-fit) 0.1875 -> 0.2100"),
  recount      = c("S9.check.s6d_counts_rederived_from_s", "a3_d_table.tsv: four/keep/(d) n_reclassified 13 -> 12"),
  invariant    = c("S9.check.s6d_branch_levels_invariant", "a3_d_table.tsv: seven(eight)/keep/(d) lower bound 0.255 -> 1.2"),
  impl_out     = c("S9.check.s6d_seven_keep_is_W6b_output", "W6b results: S6.d.a 2 -> 3"),
  n2_cross     = c("S9.check.s6d_notch2_only_seven_keep_crosses", "a3_d_table.tsv: four/keep/(d) NOTCH2_s_max 0.0927 -> 0.1001"),
  w6b_smax4    = c("S9.check.s6d_four_fit_sensitivity_is_W6b", "W6b results: S6.d.GSTA1_eGFR.s_max_4fit (four-fit sensitivity) 0.1795 -> 0.2100"),
  omni_p       = c("S9.check.omnibus_six_signal_p_match_W4b_row", "W4b results: six-signal DIAMANTE p 0.102 -> 0.202"),
  drop_ver     = c("S9.check.RAD23B_drop_version_matches_S3_adjudicator", "W3 results: RAD23B pmafdrop PP.H3(1e-5) 0.183189 -> 0.193189"),
  l1_smax      = c("S9.check.s6d_smax_L1L2_hold_seven_fit", "S6_SURVIVAL: S6H.d.UMOD_eGFR.s_max L1 0.2971914636 -> 0.2972914636 (+1e-4)"),
  fam_key      = c("S9.check.family_input_is_stage1b_key.dpep1_egfr", "S5_SURVIVAL: S5H.Q.egfr.p impl 0.11063492 -> 0.12"),
  F_nom_default = c("S9.check.nominal_literal_is_default", "fam: UMOD literal 'nominal' attach TRUE -> FALSE"),
  C_cond1      = c("S9.check.w9f.c1_condition1", "corroboration: RAD23B check_class_used concordant -> R1-only"),
  C_rad_clause = c("S9.check.w9f.rad23b_clause", "corroboration: RAD23B travel_w9f -> empty"),
  G_arg1       = c("S9.check.w9f.arg1_framing", "branches_integrated: ARG1 framing disclosure removed from the W4b ARG1-arginine gene-naming row"),
  G_omni       = c("S9.check.w9f.omnibus_disclosure", "branches_integrated: W4b T2D omnibus front caveat -> none declared"),
  G_a2         = c("S9.check.w9f.a2_merged", "branches_integrated: W4a T2D caveat 'in FinnGen' removed"),
  G_eas        = c("S9.check.w9f.eas_derived_route", "branches_integrated: W7 East Asian eGFR front caveat -> none declared"),
  G_s8         = c("S9.check.w9f.s8_front_note", "branches_integrated: W7 section 8 front_surface_note -> empty"),
  G_rad_s3     = c("S9.check.w9f.rad23b_s3_disclosure", "branches_integrated: W3 RAD23B pmafdrop travel -> empty"),
  G_s6d        = c("S9.check.w9f.s6d_travel", "branches_integrated: W6b (d) travel loses 'seven fits: n_reclassified'"),
  G_nominal    = c("S9.check.w9f.nominal_literal_default", "branches_integrated: W5 metabolite-Q qualifier 'nominal [W9' -> 'no 'nominal' [W9'"),
  # W9F step 2b (2026-10-09; defects D1-D4 of the independent W9F verifier). 'exact': the named assert must be the ONLY one that fails.
  Q_cq03       = c("S9.check.CQ03_reports_every_full_window_reading", "convention_questions: CQ03 effect_on_branches -> the pre-fix text 'none under Rule 1; ...' (D1)", "exact"),
  s3_pp        = c("S9.check.s3_adj_branches_rederived", "adj_S3_scenarios.tsv: PAM sub cap_both_fullwindows PP.H3(1e-6) 0.1538 -> 0.18, so its Rule 1 branch label 1 no longer follows from its PP values (D1)", "exact"),
  S_scen       = c("S9.check.cond2_invariant_to_S3_conventions.PAM_T2D", "S3 adjudicator scenario table: the adj_S3_scenarios.tsv PAM sub cap_both_fullwindows Rule 1 branch 1 -> 2_addon (a row the pre-fix scope never read) (D2)", "exact"),
  C_arg1_framing = c("S9.check.w9f.arg1_framing_corrob", "corroboration: ARG1-T2D travel_w9f (the ARG1 framing disclosure) -> empty (D3)", "exact"),
  C_scen_empty = c("S9.check.cond2_invariant_to_S3_conventions.PCSK9_CAD", "corroboration: PCSK9 s3_branch_all_adjudicator_scenarios -> empty (the assert must not pass vacuously; step 2c)", "exact"),
  # Stage 2 opening fixes (2026-10-09). 'exact': the named assert must be the ONLY one that fails.
  V_vis        = c("S9.check.visible_label_corrected", "S4_SURVIVAL: S4H.arg.p comparator status SURVIVES_VISIBLE_TARGET -> SURVIVES (the comparator's flagged set shrinks from 67 to 66)", "exact"),
  X_cond2na    = c("S9.check.cond2_not_applicable_where_s3_not_run", "corroboration: FGFR4-T2D (EAS) cond2_s3 'not applicable (section 3 not run)' -> TRUE", "exact"),
  R7_unclass   = c("S9.check.rule7_every_travelling_number_classified", "NOT_KEYS list: S5.s.protein.max removed, so a quantity id on W5 row 3 is neither key nor not-key", "exact"),
  R7_label     = c("S9.check.rule7_label_on_every_notkey_row", "branches_integrated: the Rule 7 label removed from the W4a test (iii) row, whose beside values are not-keys", "exact"),
  R7_br        = c("S9.check.rule7_label_on_W9_branches_rows", "branches (W9): the Rule 7 label removed from the East Asian eGFR scale row", "exact"),
  R7_s12       = c("S9.check.rule7_s12_labelled", "results: the Rule 7 label removed from the first section 12 row", "exact"),
  D_count      = c("S9.check.dispositions_count_equals_defect_lines", "dispositions: the first row dropped, so 59 dispositions face 60 defect lines", "exact"),
  V_rep        = c("S9.check.visible_keys_match_STAGE1B_REPORT", "STAGE1B_REPORT section 3 parse: S6H.primary.p dropped from the genuinely-visible keys", "exact"),
  R7_ms        = c("S9.check.rule7_merge_rows_labelled", "results: the Rule 7 label removed from the first S9.merge.dpep1 s row", "exact"),
  R7_travel    = c("S9.check.rule7_travel_text_ids_classified", "branches_integrated: W6b row 18 (the 6(d) fit-set travel text) notkey_quantity_ids emptied", "exact"),
  R7_adj       = c("S9.check.adjudicator_label_on_s6d_values", "branches_integrated: W6b row 18 adjudicator label emptied", "exact"),
  RES_adj6d    = c("S9.check.adjudicator_label_on_s6d_values", "results: the adjudicator label removed from the first S9.s6d row", "exact"),
  BR_adj       = c("S9.check.rule7_label_on_W9_branches_rows", "branches (W9): the adjudicator label removed from the W9F item 2a (6(d)) row", "exact"),
  R7_cor       = c("S9.check.rule7_corroboration_label", "corroboration: PAM-T2D outcome_check_label loses the r0 READ_FROM_SEEN Rule 7 label", "exact"),
  BR_c1        = c("S9.check.rule7_label_on_W9_branches_rows", "branches (W9): the Rule 7 label removed from the PAM-T2D corroboration row, which prints the C1 outcome-check label (r0 = 0.061)", "exact"),
  D_cq         = c("S9.check.dispositions_reference_existing_CQ", "dispositions: row 5 mapped to CQ99, which is not in convention_questions.tsv", "exact"),
  N_omni_val   = c("S9.check.w9f.omnibus_disclosure", "branches_integrated: W4b T2D omnibus qualifiers 'T2D DIAMANTE p 0.10' -> 'p 0.90' (D4; the verifier's V_omni_num)", "exact"))
rs <- file.path(R.home("bin"), "Rscript")
for (id in names(W9F_MUT)) {
  target <- W9F_MUT[[id]][1]
  rd <- file.path(MD, paste0("run_", id)); dir.create(rd, showWarnings = FALSE)
  af <- file.path(rd, "asserts.tsv")
  if (file.exists(af)) file.remove(af)   # never read a previous run's file
  old <- Sys.getenv("W9_MUTATE"); Sys.setenv(W9_MUTATE = id)
  rc <- system2(rs, shQuote(file.path(OUT, "w9_main.R")), stdout = FALSE, stderr = file.path(rd, "stderr.txt"))
  Sys.setenv(W9_MUTATE = old)
  if (rc != 0L || !file.exists(af)) { rec(paste0("W.", id), target, W9F_MUT[[id]][2], FALSE, sprintf("mutation run failed (exit %d); see mutations/run_%s/", rc, id)); next }
  a <- fread(af, sep = "\t", colClasses = "character")
  hit <- a$status[a$id == target]
  if (length(hit) != 1L) { rec(paste0("W.", id), target, W9F_MUT[[id]][2], FALSE, sprintf("named assert appears %d times in the mutation run", length(hit))); next }
  others <- setdiff(a$id[a$status != "PASS"], target)
  exact <- identical(W9F_MUT[[id]][3], "exact")
  rec(paste0("W.", id), target, W9F_MUT[[id]][2], identical(hit, "ASSERT_FAILED") && (!exact || !length(others)),
      sprintf("named assert %s; %d other assert(s) also failed%s%s", hit, length(others), if (length(others)) paste0(": ", paste(others, collapse = ", ")) else "",
              if (exact) " (exact: no other assert may fail)" else ""))
}
rec("Q1", "quote_ok (quoted fragment present in its pinned source)", "fragment 'returned WEAKENS' -> 'returned STRENGTHENS'",
    trips(quote_ok(GOV$prespec_red$path, "this framing was chosen after the pre-stated test returned STRENGTHENS", "mut")), "the guard stops the run")
rec("Q1b", "quote_ok control (must NOT trip)", "none", trips(quote_ok(GOV$prespec_red$path, "this framing was chosen after the pre-stated test returned WEAKENS", "ctl")),
    "the real fragment is found; tripped must be FALSE")

# Stage 2 opening fix 1: the deviations section is extracted verbatim from STAGE1B_REPORT; D-1 .. D-4 must each be present (unit mutation + control)
s1b_txt <- paste(readLines(STAGE1B_REPORT$path, encoding = "UTF-8", warn = FALSE), collapse = "
")
rec("DV1", "deviations_section (D-1 .. D-4 each present once)", "the D-3 bullet renamed D-X", trips(deviations_section(sub("- **D-3.", "- **D-X.", s1b_txt, fixed = TRUE))), "the extraction stops")
rec("DV1b", "deviations_section control (must NOT trip)", "none", trips(deviations_section(s1b_txt)), "the real report section 1 extracts; tripped must be FALSE")
MT <- rbindlist(MT)
expect_false <- MT$mutation_id %in% c("M2e", "M3b", "W0", "Q1b", "DV1b")
MT[, pass := ifelse(expect_false, !tripped, tripped)]
fwrite(MT, file.path(OUT, "mutation_tests.tsv"), sep = "\t", quote = TRUE)
logline(sprintf("mutation tests: %d of %d behaved as required (%d mutations tripped their check; %d controls did not trip)",
                sum(MT$pass), nrow(MT), sum(MT$tripped & !expect_false), sum(!MT$tripped & expect_false)))
stopifnot(all(MT$pass))
