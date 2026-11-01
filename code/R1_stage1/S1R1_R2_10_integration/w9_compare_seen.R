# W9 COMPARE_SEEN: written AFTER results.tsv, branches.tsv and branches_integrated.tsv were pinned (run.log "RESULTS PINNED"),
# and after the full pre-declaration and amendment were opened (run.log "FULL TEXT OPENED").
# Seen values are transcribed from the full texts with their line numbers. Nothing in results.tsv is changed:
# its md5 is asserted equal to the pinned value before and after. Diagnostics computed here explain differences only.
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
logline("================ W9 w9_compare_seen.R start ================")
PINNED <- c(results = "8f5829e16b47d5cedea67bce11ab5ed0", branches = "0ae8fee0b8bc6f03faec0cce2c4163cb")   # fix round 2 re-pin (2026-10-10: one new row, the assert S9.check.s6d_smax_L1L2_hold_seven_fit; was eae4c6c9...); Stage 2 step 2 re-pin (seven-fit s_max; was c78e770c... / a797508b...); Stage 2 opening fixes verifier-round re-pin (2026-10-09; was 14612cd9... / 1d55988a...; first re-pin was aad91915... / 530c41fa...; W9F step 2b re-pin: was 0b32b78e... / d1bbba48...; before W9F 5d39aef9... / 9adcc941...)
assert_md5(GOV$prespec_full$path, GOV$prespec_full$md5, "prespec_full")
assert_md5(GOV$amend1_full$path, GOV$amend1_full$md5, "amend1_full")
assert_md5(file.path(OUT, "results.tsv"), PINNED[["results"]], "W9 results.tsv (pinned)")
assert_md5(file.path(OUT, "branches.tsv"), PINNED[["branches"]], "W9 branches.tsv (pinned)")
assert_md5(file.path(S1, WS$W4a$dir, "results.tsv"), WS$W4a$results, "W4a/results.tsv")
RS <- read_tsv_chk(file.path(OUT, "results.tsv"), "W9 results.tsv")
W4A <- read_tsv_chk(file.path(S1, WS$W4a$dir, "results.tsv"), "W4a/results.tsv")
full <- readLines(GOV$prespec_full$path, encoding = "UTF-8", warn = FALSE)
amd  <- readLines(GOV$amend1_full$path, encoding = "UTF-8", warn = FALSE)
# guard: the quoted line must contain the seen text we transcribe (a moved or edited line trips this)
seen_line <- function(txt, ln, must) { if (!grepl(must, txt[ln], fixed = TRUE)) stop(sprintf("full-text line %d does not contain '%s'", ln, must)); invisible(TRUE) }
v <- function(q) getq(RS, "W9", q)
vs <- function(q) getq(RS, "W9", q, numeric = FALSE)

CS <- list()
cmp <- function(id, seen, seen_label, s1, s1_label, tol, expl_ok, expl_diff) {
  ad <- if (is.numeric(seen) && is.numeric(s1)) abs(seen - s1) else NA_real_
  differs <- if (is.na(ad)) NA else ad > tol
  CS[[length(CS) + 1L]] <<- data.table(quantity_id = id, seen_value = as.character(seen), seen_label = seen_label,
    stage1_value = if (is.numeric(s1)) format(s1, digits = 7) else as.character(s1), stage1_label = s1_label,
    abs_diff = if (is.na(ad)) "" else format(signif(ad, 3)),
    explanation = if (isTRUE(differs)) paste0("DIFFERS beyond printed precision (tol ", format(tol), "): ", expl_diff)
                  else if (isFALSE(differs)) paste0("within printed precision (tol ", format(tol), "); ", expl_ok) else expl_ok)
}

# ---- section 12(a): pooled two-resource difference (full text l.586)
seen_line(full, 586, "marginal −0.154, 95% CI −0.259 to −0.049, p 0.004; conditional −0.173, p 0.003; between-resource Q p 0.92–0.99")
# diagnostic only: which SE source reproduces the seen interval? (W9's declared estimate uses W4a's aware 2nd-order SEs)
pool <- function(d, s) { w <- 1 / s^2; e <- sum(w * d) / sum(w); se <- sqrt(1 / sum(w)); c(e, e - qnorm(0.975) * se, e + qnorm(0.975) * se, 2 * pnorm(-abs(e / se))) }
gq <- function(q) getq(W4A, "W4a", q)
diag_rows <- list()
for (est in c("marg", "cond")) {
  d <- c(gq(sprintf("S4.grid.FinnGen.%s.diff", est)), gq(sprintf("S4.grid.DIAMANTE.%s.diff", est)))
  se_aw <- c(gq(sprintf("S4.grid.FinnGen.%s.diff_se", est)), gq(sprintf("S4.grid.DIAMANTE.%s.diff_se", est)))
  se_ign <- abs(d) / qnorm(1 - c(gq(sprintf("S4.grid.FinnGen.%s.p_ign", est)), gq(sprintf("S4.grid.DIAMANTE.%s.p_ign", est))) / 2)
  se_1st <- abs(d) / qnorm(1 - c(gq(sprintf("S4.grid.FinnGen.%s.p_aware_1st", est)), gq(sprintf("S4.grid.DIAMANTE.%s.p_aware_1st", est))) / 2)
  d_sp <- c(d[1], gq(sprintf("S4.grid.DIAMANTE_SEfromP.%s.diff", est))); se_sp <- c(se_aw[1], gq(sprintf("S4.grid.DIAMANTE_SEfromP.%s.diff_se", est)))
  for (nm in c("aware_2nd (declared)", "ign (SE backed out of p_ign)", "aware_1st (SE backed out of p_aware_1st)", "aware_2nd, DIAMANTE SE-from-p")) {
    r <- switch(nm, "aware_2nd (declared)" = pool(d, se_aw), "ign (SE backed out of p_ign)" = pool(d, se_ign),
                "aware_1st (SE backed out of p_aware_1st)" = pool(d, se_1st), "aware_2nd, DIAMANTE SE-from-p" = pool(d_sp, se_sp))
    diag_rows[[length(diag_rows) + 1L]] <- data.table(estimand = est, se_source = nm, est = r[1], lo = r[2], hi = r[3], p = r[4])
  }
}
DG <- rbindlist(diag_rows)
fwrite(DG, file.path(OUT, "compare_seen_diag_s12a.tsv"), sep = "\t", quote = TRUE)
for (i in seq_len(nrow(DG))) logline(sprintf("diag s12(a) %s %-42s est %+.4f CI %+.4f to %+.4f p %.4f", DG$estimand[i], DG$se_source[i], DG$est[i], DG$lo[i], DG$hi[i], DG$p[i]))
dm <- DG[DG$estimand == "marg"]
best <- dm$se_source[which.min(abs(dm$lo - (-0.259)) + abs(dm$hi - (-0.049)))]
dexp <- sprintf("unexplained in origin, but a post-comparison diagnostic (compare_seen_diag_s12a.tsv) shows the seen interval is closest to pooling with '%s' SEs (marginal: %s); W9 pooled W4a's aware 2nd-order SEs as the brief declares, and W4a's per-cell aware CIs agree with the seen per-cell CIs (full text l.349-350). Stage 1 value unchanged",
                best, paste(sprintf("%s CI %.3f to %.3f p %.4f", dm$se_source, dm$lo, dm$hi, dm$p), collapse = "; "))
cmp("S9.s12a.marg.est", -0.154, "pooled marginal A-B difference [full l.586]", v("S9.s12a.marg.est"), "IVW FE pool, W4a aware 2nd-order SEs", 5e-4, "same estimate", dexp)
cmp("S9.s12a.marg.ci_lo", -0.259, "pooled marginal 95% CI lower [full l.586]", as.numeric(RS$ci_lo[RS$quantity_id == "S9.s12a.marg.est"]), "IVW FE pool", 5e-4, "", dexp)
cmp("S9.s12a.marg.ci_hi", -0.049, "pooled marginal 95% CI upper [full l.586]", as.numeric(RS$ci_hi[RS$quantity_id == "S9.s12a.marg.est"]), "IVW FE pool", 5e-4, "", dexp)
cmp("S9.s12a.marg.p", 0.004, "pooled marginal p [full l.586]", as.numeric(RS$p[RS$quantity_id == "S9.s12a.marg.est"]), "IVW FE pool", 5e-4, "", dexp)
cmp("S9.s12a.cond.est", -0.173, "pooled conditional A-B difference [full l.586]", v("S9.s12a.cond.est"), "IVW FE pool, W4a aware 2nd-order SEs", 5e-4, "same estimate", dexp)
cmp("S9.s12a.cond.p", 0.003, "pooled conditional p [full l.586]", as.numeric(RS$p[RS$quantity_id == "S9.s12a.cond.est"]), "IVW FE pool", 5e-4, "",
    paste0(dexp, " (conditional diagnostic rows likewise in compare_seen_diag_s12a.tsv)"))
cmp("S9.s12a.marg.Q_p", 0.99, "between-resource Q p, upper end of the printed range 0.92-0.99 [full l.586]", as.numeric(RS$p[RS$quantity_id == "S9.s12a.marg.Q"]), "Cochran Q, df 1", 0.005, "the range endpoint matches the marginal Q p", "")
cmp("S9.s12a.cond.Q_p", 0.92, "between-resource Q p, lower end of the printed range 0.92-0.99 [full l.586]", as.numeric(RS$p[RS$quantity_id == "S9.s12a.cond.Q"]), "Cochran Q, df 1", 0.005, "the range endpoint matches the conditional Q p", "")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.s12a.fate", seen_value = "supplement only, labelled post hoc; never in the Abstract, Research Insights, a figure or the response letter",
  seen_label = "section 12 fate [full l.586]", stage1_value = "results.tsv notes: POST HOC, SUPPLEMENT-ONLY; never in a branch", stage1_label = "W9 label",
  abs_diff = "", explanation = "consistent; the full-text fate also bars figures and the response letter, which W9's label does not name: Stage 2 applies the full-text fate")

# ---- section 12(b): pre-reclassification Holm (full text l.376-378)
seen_line(full, 377, "four cells: 0.0921 ign / 0.1032 aware"); seen_line(full, 378, "two originally reported cells: 0.0585 ign / 0.0516 aware")
for (cv in c("ign", "aware")) {
  s4 <- if (cv == "ign") 0.0921 else 0.1032; s2 <- if (cv == "ign") 0.0585 else 0.0516
  for (cell in c("FinnGen.marg", "FinnGen.cond", "DIAMANTE.marg", "DIAMANTE.cond"))
    cmp(sprintf("S9.s12b.grid4.%s.%s", cv, cell), s4, sprintf("Holm over the four grid cells, %s [full l.377]", cv), v(sprintf("S9.s12b.grid4.%s.%s", cv, cell)),
        "Holm over W4a's four cells", 5e-5, "deterministic from W4a's grid p values, which match the seen grid p (full l.336-343): weak evidence", "")
  for (cell in c("FinnGen.cond", "DIAMANTE.marg"))
    cmp(sprintf("S9.s12b.orig2.%s.%s", cv, cell), s2, sprintf("Holm over the two originally reported cells, %s [full l.378]", cv), v(sprintf("S9.s12b.orig2.%s.%s", cv, cell)),
        "Holm over FinnGen conditional and DIAMANTE marginal", 5e-5, "deterministic from W4a's p values: weak evidence", "")
}
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.s12b.fate", seen_value = "supplement, as a conservative bound, labelled as computed before the grid was reclassified as a stability display",
  seen_label = "section 12 fate [full l.588]", stage1_value = "results.tsv notes: POST HOC, SUPPLEMENT-ONLY; the grid is a stability display, not a family (Rule 4)", stage1_label = "W9 label",
  abs_diff = "", explanation = "consistent; Stage 2 uses the full-text label wording ('a conservative bound ... computed before the grid was reclassified')")

# ---- Rule 4 family members: unadjusted p entered into Holm
seen_line(full, 333, "marginal p 0.947 ign / 0.946 aware"); seen_line(full, 336, "ign p 0.0336, aware 0.0365"); seen_line(full, 359, "conditional p 0.0462; marginal 0.0620")
seen_line(full, 386, "5.3 × 10⁻⁵ / 5.7 × 10⁻⁵ (aware)"); seen_line(full, 444, "No eGFR-per-DPEP1 ratio has been computed"); seen_line(full, 504, "no all-signal Q has been computed")
cmp("S9.holm.arg1_arginine.p_unadj", 0.946, "arginine A vs B marginal aware p [full l.333]", v("S9.holm.arg1_arginine.p_unadj"), "W4a S4.arg.marg.p_aware_2nd", 5e-4, "weak evidence (deterministic, visible definition)", "")
cmp("S9.holm.arg1_t2d.p_unadj", 0.0365, "T2D FinnGen conditional aware p [full l.336]", v("S9.holm.arg1_t2d.p_unadj"), "W4a S4.grid.FinnGen.cond.p_aware_2nd", 5e-5, "weak evidence (deterministic)", "")
cmp("S9.holm.arg1_test3.p_cond", 0.0462, "test (iii) conditional aware p [full l.359]", v("S9.holm.arg1_test3.p_cond"), "W4a S4.t3.FinnGen.cond.503.2nd.p_aware", 5e-5, "weak evidence", "")
cmp("S9.holm.arg1_test3.p_unadj", 0.0620, "test (iii) marginal aware p = the larger, entered in the family [full l.359]", v("S9.holm.arg1_test3.p_unadj"), "max(cond, marg)", 5e-5, "the larger estimand p is the marginal one in both", "")
cmp("S9.holm.arg1_ornithine.p_unadj", 5.3e-5, "ornithine aware p, 503, scoping two_signal.tsv dev_sd test (conditional) [full l.386]", v("S9.holm.arg1_ornithine.p_unadj"),
    "W4a S4.orn.503.p_aware_2nd (declared primary: marginal, Rule 5 delta-method, aware 2nd order)", 5e-6, "",
    "different statistic and estimand: the seen value is the engine's conditional dev_sd z-space test; the family entry is the declared marginal Rule 5 test (W4a COMPARE_SEEN rows S4.orn.* explain the scale). Same side of 0.05 and the same Holm outcome (Holm p < 0.05 either way: 5.3e-5 x 7 = 3.7e-4)")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.holm.tgfb1_cad.p_unadj", seen_value = "naive within-European heterogeneity z about 6.2 (marginal, not the declared conditional estimand)",
  seen_label = "full l.212", stage1_value = format(v("S9.holm.tgfb1_cad.p_unadj"), digits = 4), stage1_label = "W2 S2.het.finngen_chd.p_aware_2nd (conditional, aware)",
  abs_diff = "", explanation = "no seen value for the declared estimand (FinnGen CHD conditional); the seen naive marginal z is a different quantity; first computed by Stage 1")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.holm.dpep1_cysgly.p_unadj|S9.holm.dpep1_egfr.p_unadj", seen_value = "none: no all-signal Q seen; 'No eGFR-per-DPEP1 ratio has been computed'",
  seen_label = "full l.430-444", stage1_value = sprintf("%s | %s", format(v("S9.holm.dpep1_cysgly.p_unadj"), digits = 4), format(v("S9.holm.dpep1_egfr.p_unadj"), digits = 4)),
  stage1_label = "W5 all-signal GLS Q, aware 2nd, 503", abs_diff = "", explanation = "first computed by Stage 1 (new, not seen)")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.holm.umod_egfr.p_unadj", seen_value = "none: 'no all-signal Q has been computed'", seen_label = "full l.504",
  stage1_value = format(v("S9.holm.umod_egfr.p_unadj"), digits = 4), stage1_label = "W6b all-signal GLS Q, aware 2nd, 503", abs_diff = "", explanation = "first computed by Stage 1 (new, not seen)")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.holm.*.p_holm|S9.holm.*.nominal", seen_value = "none: no Holm value over the 8-member family appears in the full text",
  seen_label = "full text searched for the family Holm (only the grid Holm, l.376-378, is seen)", stage1_value = paste(sprintf("%s %s", sub("^S9\\.holm\\.", "", sub("\\.p_holm$", "", RS$quantity_id[grepl("\\.p_holm$", RS$quantity_id)])),
  format(as.numeric(RS$value[grepl("\\.p_holm$", RS$quantity_id)]), digits = 3)), collapse = "; "), stage1_label = "Holm over the 8 members", abs_diff = "",
  explanation = "first computed by Stage 1 (new). Three of the eight inputs (arginine, T2D, test (iii)) equal seen values, so the T2D 'nominal' outcome was foreseeable from seen numbers")

# ---- Amendment 1 A2 (full amendment l.37-41)
seen_line(amd, 38, "FinnGen conditional: aware p 0.0365"); seen_line(amd, 39, "DIAMANTE marginal: aware p 0.0258 (ign 0.0293)"); seen_line(amd, 41, "On these values the condition would be met")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.check.A2_matches_W4a", seen_value = "FinnGen cond aware p 0.0365; DIAMANTE marg aware p 0.0258; both A-B negative; condition would be met",
  seen_label = "full amendment l.37-41", stage1_value = "met (FinnGen 0.0365, DIAMANTE 0.0258, both negative)", stage1_label = "W9 re-derivation from W4a", abs_diff = "",
  explanation = "same; weak evidence: deterministic values whose inputs were seen")

# ---- corroboration
seen_line(full, 234, "branch (a) is the arithmetic expectation"); seen_line(full, 262, "Nothing for these four regional tests")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.corrob.TGFB1_CAD.corroborated", seen_value = "section 2 'branch (a) is the arithmetic expectation' (implies not corroborated)",
  seen_label = "full l.234", stage1_value = vs("S9.corrob.TGFB1_CAD.corroborated") , stage1_label = "Principle 0 condition 3 fails (section 2 branch (a))", abs_diff = "",
  explanation = "agreement with an expectation recorded before Stage 1: weak evidence")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.corrob.{PAM_T2D,RAD23B_T2D,PCSK9_CAD,TGFB1_CAD}.s3_branch", seen_value = "none: 'Nothing for these four regional tests'",
  seen_label = "full l.262", stage1_value = paste(sapply(c("PAM_T2D", "RAD23B_T2D", "PCSK9_CAD", "TGFB1_CAD"), function(k) paste0(k, " ", vs(sprintf("S9.corrob.%s.s3_branch", k)))), collapse = "; "),
  stage1_label = "W3 branches re-derived by W9", abs_diff = "", explanation = "first computed by Stage 1 (new); so the corroboration verdicts for PAM, RAD23B and PCSK9 are new")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.corrob.*.check_class", seen_value = "R1_outcome_checks_tsv class column (seen_output; md5 7d7f1431...)", seen_label = "R1 declared record, printed beside (READ_FROM_SEEN)",
  stage1_value = paste(sapply(c("PAM_T2D", "RAD23B_T2D", "ARG1_T2D", "PCSK9_CAD", "TGFB1_CAD", "FGFR4_T2D_EAS"), function(k) paste0(k, " ", vs(sprintf("S9.corrob.%s.check_class", k)))), collapse = "; "),
  stage1_label = "W9F: condition 1 from C1's verified classes (S1R1_R2_11_extchecks; Stage 1b, from raw; comparator SURVIVES)", abs_diff = "",
  explanation = "same class for all six pairs (asserts S9.check.C1_class_agrees_R1_and_W1.*); agreement with R1 is weak evidence (deterministic quantities from the same files), but the classes now rest on C1, which both Stage 1b lenses reproduced. W1's ext_checks.tsv also agrees for the five European pairs. ARG1 signal B (not a C1 pair) still rests on the R1 record")
seen_line(full, 505, "(d), illustrative:")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.s6d.*", seen_value = "(d), illustrative: reclassifying UMOD-eGFR, UMOD-CKD and DPEP1-eGFR (three candidates) gives p 0.241 (OR 2.31, 0.63-9.42)",
  seen_label = "full l.505", stage1_value = paste(vapply(c("seven", "four"), function(fs) sprintf("%s fits: OR %s (lower %s)", fs, vs(sprintf("S9.s6d.%s.keep_shift.d_estimate.or", fs)),
    format(as.numeric(RS$ci_lo[RS$quantity_id == sprintf("S9.s6d.%s.keep_shift.d_estimate.or", fs)]), digits = 4)), ""), collapse = "; "),
  stage1_label = "W9F item 2a: seven-fit (primary, user decision 2026-10-09) and four-fit (sensitivity) readings, read from the S6 adjudicator's files (A_third)", abs_diff = "",
  explanation = "different quantity: the seen value reclassifies an illustrative three candidates, the Stage 1 (d) applies Rule 3(i) to all 16; W6b's COMPARE_SEEN holds the implementer comparison. The seven-fit (primary) reading equals W6b's output (assert S9.check.s6d_seven_keep_is_W6b_output); the four-fit reading is the sensitivity")

# ---- DPEP1 merged s (full text l.451; W5 COMPARE_SEEN holds the per-set s comparisons)
seen_line(full, 451, "inadmissible for DPEP1 (protein s 0.40)")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.merge.dpep1.admissible.cysgly_pair", seen_value = "inadmissible (protein s 0.40)", seen_label = "full l.451 (recorded scoping classification)",
  stage1_value = vs("S9.merge.dpep1.admissible.cysgly_pair"), stage1_label = "Rule 3(i) on the merged table (FALSE = inadmissible)", abs_diff = "",
  explanation = "NO CLASS_FLIP. Per-set s differences are W5's (n choice; W5 COMPARE_SEEN S5.s.*)")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "S9.merge.dpep1.admissible.egfr_pair", seen_value = "none: no DPEP1-eGFR s or class recorded ('No eGFR-per-DPEP1 ratio has been computed', l.444)",
  seen_label = "full l.444", stage1_value = vs("S9.merge.dpep1.admissible.egfr_pair"), stage1_label = "Rule 3(i) on W6b's DPEP1-eGFR fits", abs_diff = "",
  explanation = "first computed by Stage 1 (new); no recorded class, so no CLASS_FLIP is possible")
CS[[length(CS) + 1L]] <- data.table(quantity_id = "RNG", seen_value = "n/a", seen_label = "brief section 0 rule 8", stage1_value = "W9 uses no random numbers", stage1_label = "",
  abs_diff = "", explanation = "no 'same RNG and seed' statement applies")

CSD <- rbindlist(CS)
fwrite(CSD, file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = TRUE)
logline(sprintf("COMPARE_SEEN.tsv: %d rows; %d flagged DIFFERS", nrow(CSD), sum(grepl("^DIFFERS", CSD$explanation))))
assert_md5(file.path(OUT, "results.tsv"), PINNED[["results"]], "W9 results.tsv unchanged after compare")
assert_md5(file.path(OUT, "branches.tsv"), PINNED[["branches"]], "W9 branches.tsv unchanged after compare")
logline(sprintf("COMPARE_SEEN pinned md5 %s", md5(file.path(OUT, "COMPARE_SEEN.tsv"))))
