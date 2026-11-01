# W9 README tables: every table of numbers in README.md is generated here from the output files, between
# <!-- BEGIN:<name> --> and <!-- END:<name> --> markers (each marker pair must occur exactly once). Nothing in those
# blocks is typed by hand. Run last by run_all.R. Writes only README.md (the marked blocks).
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
logline("================ W9 w9_readme.R start ================")
rd <- function(f) read_tsv_chk(file.path(OUT, f), f)
RS <- rd("results.tsv"); COR <- rd("corroboration.tsv"); D6R <- rd("s6d_fitset_readings.tsv"); CQ <- rd("convention_questions.tsv")
AS <- rd("asserts.tsv"); MT <- rd("mutation_tests.tsv"); CS <- rd("COMPARE_SEEN.tsv"); PP <- rd("pin_proofs.tsv")
v <- function(q, col = "value") getq(RS, "W9", q, col, numeric = FALSE)
esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
md_table <- function(df) {
  h <- paste0("| ", paste(names(df), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  b <- apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |"))
  paste(c(h, s, b), collapse = "\n")
}
p2 <- function(x) { x <- as.numeric(x); if (x >= 0.995) "1" else formatC(x, digits = 2, format = "g", flag = "#") }

# ---- Holm (Rule 4) ----
mem <- c("tgfb1_cad", "arg1_ornithine", "arg1_t2d", "arg1_test3", "dpep1_egfr", "umod_egfr", "arg1_arginine", "dpep1_cysgly")
mem <- mem[order(as.numeric(sapply(mem, function(m) v(paste0("S9.holm.", m, ".rank")))))]
HB <- data.frame(check.names = FALSE,
  `Member (primary cell)` = sapply(mem, function(m) RS$estimand[RS$quantity_id == paste0("S9.holm.", m, ".p_unadj")]),
  `unadjusted p` = sapply(mem, function(m) p2(v(paste0("S9.holm.", m, ".p_unadj")))),
  `Holm p (m = 8)` = sapply(mem, function(m) p2(v(paste0("S9.holm.", m, ".p_holm")))),
  `"nominal": literal reading (default)` = sapply(mem, function(m) v(paste0("S9.holm.", m, ".nominal"))),
  `"nominal": applied reading I1 (beside)` = sapply(mem, function(m) v(paste0("S9.holm.", m, ".nominal_I1"))),
  `Stage 1b key status` = sapply(mem, function(m) v(paste0("S9.holm.", m, ".stage1b_status"))))
# ---- corroboration ----
CB <- data.frame(check.names = FALSE, Pair = COR$pair,
  `Condition 1: C1 class (verified)` = ifelse(COR$check_class_C1 == "not a C1 pair", "not a C1 pair", paste0(COR$check_class_C1, " (", COR$C1_class_status, ")")),
  `R1 class (beside; READ_FROM_SEEN)` = COR$check_class_R1,
  `cond. 1` = COR$cond1_concordant, `section 3 branch` = COR$s3_branch, `cond. 2` = COR$cond2_s3, `section 2 branch` = COR$s2_branch, `cond. 3` = COR$cond3_s2,
  Verdict = COR$verdict, Reason = COR$reason)
# ---- section 6(d) ----
fmtn <- function(x) { x <- as.numeric(x); ifelse(!is.finite(x), "Inf", ifelse(x >= 10, sprintf("%.1f", x), sprintf("%.3f", x))) }
DB <- data.frame(check.names = FALSE, `Fit set` = D6R$fit_set, `Build-shift rows` = D6R$build_shift_rows,
  Version = ifelse(grepl("^\\(d\\)", D6R$ING1_version), "(d) estimate", "'kept' version"),
  n_reclassified = D6R$n_reclassified, `a/b/c/d` = paste(D6R$a, D6R$b, D6R$c, D6R$d, sep = "/"),
  OR = fmtn(D6R$or), `95% CI` = paste0(fmtn(D6R$lo), "-", fmtn(D6R$hi)), p = sapply(D6R$p, p2),
  `NOTCH2 s_max` = sprintf("%.4f", as.numeric(D6R$NOTCH2_s_max)), D6 = D6R$d6_base, `sensitive to` = D6R$qualifier_list)
# ---- convention questions ----
QB <- data.frame(check.names = FALSE, ID = CQ$id, Question = CQ$question, `Main-text slot` = CQ$main_text_slot,
  `Reading A -> consequence` = paste0(CQ$reading_A, ": ", CQ$consequence_A), `Reading B -> consequence` = paste0(CQ$reading_B, ": ", CQ$consequence_B),
  `Effect on branches` = CQ$effect_on_branches, `Raised by` = CQ$raised_by)
# ---- gates ----
GB <- paste0(
  "- **Step 0 pin proofs** (`pin_proofs.tsv`): ", sum(PP$pass == "TRUE"), " of ", nrow(PP), " pass (", sum(grepl("^P-M", PP$item) & PP$pass == "TRUE"), " of them mutations that tripped).\n",
  "- **Asserts** (`asserts.tsv`): ", sum(AS$status == "PASS"), " of ", nrow(AS), " PASS; ", sum(AS$status != "PASS"), " ASSERT_FAILED.\n",
  "- **Mutation tests** (`mutation_tests.tsv`): ", sum(MT$pass == "TRUE"), " of ", nrow(MT), " behaved as required; ", sum(grepl("^W\\.", MT$mutation_id)),
  " are W9F end-to-end mutations, each tripping the assert it names.\n",
  "- **COMPARE_SEEN** (`COMPARE_SEEN.tsv`): ", nrow(CS), " rows, ", sum(grepl("^DIFFERS", CS$explanation)), " flagged DIFFERS.")
# ---- section 12 (post hoc, supplement-only) ----
f3 <- function(x) sprintf("%.3f", as.numeric(x))
SA <- data.frame(check.names = FALSE, Estimand = c("Marginal", "Conditional"),
  `Pooled A-B difference (95% CI)` = sapply(c("marg", "cond"), function(e) sprintf("%s (%s to %s)", f3(v(sprintf("S9.s12a.%s.est", e))), f3(v(sprintf("S9.s12a.%s.est", e), "ci_lo")), f3(v(sprintf("S9.s12a.%s.est", e), "ci_hi")))),
  p = sapply(c("marg", "cond"), function(e) p2(v(sprintf("S9.s12a.%s.est", e), "p"))),
  `Cochran Q p` = sapply(c("marg", "cond"), function(e) p2(v(sprintf("S9.s12a.%s.Q", e), "p"))))
SB <- data.frame(check.names = FALSE, `Set of cells` = c("Four grid cells", "Two originally reported cells"),
  aware = c(sprintf("%.4f", as.numeric(v("S9.s12b.grid4.aware.FinnGen.cond"))), sprintf("%.4f", as.numeric(v("S9.s12b.orig2.aware.FinnGen.cond")))),
  ign = c(sprintf("%.4f", as.numeric(v("S9.s12b.grid4.ign.FinnGen.cond"))), sprintf("%.4f", as.numeric(v("S9.s12b.orig2.ign.FinnGen.cond")))))
chr_df <- function(d) { d[] <- lapply(d, as.character); d }
# ---- Stage 2 opening fixes (2026-10-09): blocks generated from the new outputs ----
IRd <- rd("branches_integrated.tsv"); BRd <- rd("branches.tsv"); KL <- rd("stage1b_key_labels.tsv"); DPT <- rd("defect_dispositions.tsv")
# item 5 bullet (fix 4: scope of "travels")
ir_eas <- IRd$source_row[IRd$source == "W7" & grepl("^EAS_eGFR", IRd$branch)]; stopifnot(length(ir_eas) == 1L, grepl("derived route", IRd$travel_w9f[IRd$source == "W7" & IRd$source_row == ir_eas], fixed = TRUE),
  sum(grepl("derived route", IRd$travel_w9f, fixed = TRUE)) == 1L)
b_item5 <- sprintf('- **Item 5. East Asian eGFR scale claim**: "derived route; not new evidence" is attached to the East Asian eGFR scale claim row (W7 row %s) and is its A1 front caveat; manuscript surfaces do not exist yet.', ir_eas)
# section 6(d) scope sentence (fix 4: what W9 re-checks)
b_s6dscope <- paste0("5. **Section 6(d): fit set DECIDED 2026-10-09 (seven fits primary; four fits a sensitivity).** Every number printed in the fit-set rows is read from the S6 adjudicator's files. ",
  "W9 recounts each table's counts (n_reclassified and a/b/c/d) from the per-candidate s values and re-derives D6 and the qualifier list from each reading's CI; ",
  "the s values themselves are read from the adjudicator's files and are not recomputed by W9. ",
  "Since Stage 2 step 2, W6b's own output is the seven-fit reading (and `S6.d.<cand>.s_max_4fit` the sensitivity): the seven-fit, rows-kept rows equal W6b's own output, ",
  "and the adjudicator's four-fit and seven-fit maxima equal W6b's `s_max_4fit` and `s_max` for every computable candidate (S9.check.s6d_seven_fit_primary_is_W6b, S9.check.s6d_four_fit_sensitivity_is_W6b).")
# visibility (fix 2)
n_raw <- sum(KL$status_comparator == "SURVIVES_VISIBLE_TARGET"); n_ind <- sum(KL$status_corrected == VIS_LABEL_INDEP)
stopifnot(setequal(KL$key_id[startsWith(KL$status_corrected, "SURVIVES_VISIBLE_TARGET")], VIS_GENUINE), setequal(KL$key_id[KL$status_corrected == VIS_LABEL_DERIVABLE], VIS_DERIVABLE),
          nrow(KL) == 134L, n_raw == length(VIS_GENUINE) + n_ind)
s6d_all <- KL$key_id[nzchar(KL$qualifier) & grepl("^S6H[.]d[.].+[.]s_max$", KL$key_id)]; s6d_ind <- s6d_all[KL$status_corrected[match(s6d_all, KL$key_id)] == VIS_LABEL_INDEP]
stopifnot(length(s6d_all) == 15L)
b_visibility <- paste0(
  "   **Stage 1b status labels (corrected at Stage 2, H2 D2).** The comparator's `SURVIVES_VISIBLE_TARGET` is a lexical token match, not a semantic one: it flagged ", n_raw, " of ", nrow(KL),
  " keys. W9 prints it only for the ", length(VIS_GENUINE), " keys whose answer was printed in a file the verifiers could open (", paste0("`", VIS_GENUINE, "`", collapse = ", "), "). ",
  paste0("`", VIS_DERIVABLE, "`", collapse = ", "), " are printed '", VIS_LABEL_DERIVABLE, "', because they are closed-form functions of those visible counts. ",
  "The other ", n_ind, " keys the comparator flagged are printed '", VIS_LABEL_INDEP, "': their matches were thresholds, one-digit tokens or unrelated rows (for example `S4H.arg.p` matched '.9' in `RESOURCES.json`), so they are independent re-derivations of the value the implementer held when the comparison was made. ",
  "For the ", length(s6d_ind), " `S6H.d.<cand>.s_max` keys among them, that value was the earlier FOUR-fit s_max, whereas the printed s_max has been the seven-fit value since Stage 2 step 2. ",
  "For all ", length(s6d_all), " computable `S6H.d.<cand>.s_max` keys the verifier lenses L1 and L2 already hold seven-fit values that match the printed seven-fit values within ", S6D_SMAX_TOL_TXT, " (assert `S9.check.s6d_smax_L1L2_hold_seven_fit`). ",
  "The ", sum(KL$status_corrected[KL$key_id %in% s6d_all] == "L1L2_AGREE_IMPL_DIFFERS"), " `L1L2_AGREE_IMPL_DIFFERS` statuses among them are stale only because the implementer value they compared against was the four-fit one, and the `SURVIVES` labels on those keys were compared against the four-fit value; none is re-labelled (column `qualifier` of `stage1b_key_labels.tsv`). ",
  "The raw label stays in `stage1b_key_labels.tsv` (column `status_comparator`). The corrected list is derived in code from a constant of ", length(VIS_KEYS8), " key ids (assert `S9.check.visible_label_corrected`), and the 5 + 3 ids are checked against the lists in STAGE1B_REPORT section 3, parsed from the pinned report (assert `S9.check.visible_keys_match_STAGE1B_REPORT`).")
# condition 2 (fix 3)
na_pairs <- COR$pair[COR$cond2_s3 == COND2_NA]
b_cond2 <- paste0(
  "   - Condition 2 is asserted to hold, or to fail, identically under every section 3 scenario the S3 adjudicator ran (prior cap; DIAMANTE EAF printed 0/1).\n",
  "   - Condition 2 reads '", COND2_NA, "' for ", paste(na_pairs, collapse = ", "), ". Their verdicts rest on the other conditions (the external check; for ARG1-T2D signal B, none, because it is not a pair). ",
  "Section 3 was not run for them, so condition 2 is vacuous there; it was printed as TRUE before Stage 2 (H2 D6). The verdicts are unchanged (assert `S9.check.cond2_not_applicable_where_s3_not_run`).")
# rule 7 amendment (fix 5)
ISC <- IRd$n_key_numbers != ""   # claim or fragment rows, plus any row carrying a registered W9F travel text
R7T <- do.call(rbind, lapply(sort(unique(IRd$source[ISC])), function(w) { s <- ISC & IRd$source == w
  data.frame(check.names = FALSE, Workstream = w, `Rows in scope` = sum(s), `Rows labelled` = sum(as.integer(IRd$n_notkey_numbers[s]) > 0L),
             `Key numbers` = sum(as.integer(IRd$n_key_numbers[s])), `Not-key numbers (labelled)` = sum(as.integer(IRd$n_notkey_numbers[s]))) }))
R7T <- rbind(R7T, data.frame(check.names = FALSE, Workstream = "All", `Rows in scope` = sum(R7T[[2]]), `Rows labelled` = sum(R7T[[3]]), `Key numbers` = sum(R7T[[4]]), `Not-key numbers (labelled)` = sum(R7T[[5]])))
s12r <- RS[startsWith(RS$quantity_id, "S9.s12"), ]
b_rule7 <- paste0(
  "Rule 7 of the pre-declaration bars Stage 1 numbers from the manuscript until the section 10 panel has re-derived them. Under the middle path (deviation D-1) only the 134 headline keys were re-derived. ",
  "User decision 12 (2026-10-09) amends Rule 7: every other number that travels on an integrated claim row is admitted, labelled **", RULE7_LABEL, "**, in supplement or side positions only.\n\n",
  "- **Key vs not-key.** A key is a quantity id that is a single-id row of a `KEY_MAP.tsv` (the 104 numeric keys; the 30 categorical keys are branch or class levels, not numbers; the six `C1.<pair>.class` results rows of the C1 class keys count as keys). Every other quantity id (`NOT_KEYS.tsv`, 3,119 rows) is a not-key.\n",
  "- **Which numbers travel.** The quantity ids of each `claim` or `fragment` row of `branches_integrated.tsv`, and of any other row carrying a registered W9F travel text: the workstream's own id list where the Stage 1b normalisers wrote one (`branches_inputs_ids.tsv` for W3, W4a, W5), otherwise the ids parsed from `inputs`; plus W1's inline numbers matched by value (declared in `w9_main.R`); plus the ids that each W9F travel text is built from (registry `TRAVEL_REG` in `w9_main.R`: the 6(d) fit-set text, the RAD23B p-source clause, the six-signal omnibus p values, the RAD23B drop-version text). A row's not-key ids are written to `notkey_quantity_ids`; `S9.check.rule7_travel_text_ids_classified` asserts that every row containing a registered text carries them.\n",
  "- **Granularity.** The label is attached per row, not per printed number: a row that mixes keys and not-keys carries it, and the counts are of quantity ids, not of printed numbers. It is appended to `qualifiers_resolved` and written to `rule7_label`. It is never the A1 front caveat or the front-surface note.\n",
  "- **Named cases.** Test (iii) beside values (DIAMANTE, 489, beta-space, first-order) are in W4a row 18; the six-signal omnibus p values are in the W4b six-signal row; the section 12 pooled values are the ", nrow(s12r), " `S9.s12*` rows of `results.tsv`. All carry the label.\n",
  "- **Adjudicator values.** The section 6(d) fit-set tables and the seven-fit and NOTCH2 s values are the S6 adjudicator's third computation, not a Stage 1b key and not analyst-computed. Where they travel they carry their own label, **", ADJ_LABEL, "**: ",
  sum(IRd$rule7_adjudicator_label != ""), " integrated rows (column `rule7_adjudicator_label`), ", sum(grepl(ADJ_LABEL, BRd$qualifiers, fixed = TRUE)), " W9 `branches.tsv` row(s) and ", sum(startsWith(RS$quantity_id, "S9.s6d.")), " `S9.s6d.*` result rows (`S9.check.adjudicator_label_on_s6d_values`).\n",
  "- **Result rows.** ", nrow(s12r), " `S9.s12*` rows and ", sum(startsWith(RS$quantity_id, "S9.merge.dpep1.s") & grepl(RULE7_LABEL, RS$notes, fixed = TRUE)), " of the `S9.merge.dpep1` s and s_max rows (those whose quantity ids are not-keys) carry the analyst label in `notes` (`S9.check.rule7_s12_labelled`, `S9.check.rule7_merge_rows_labelled`).\n",
  "- **W9's own `branches.tsv`.** ", sum(grepl(RULE7_LABEL, BRd$qualifiers, fixed = TRUE)), " of ", nrow(BRd), " rows quote not-key quantity ids (from `inputs` and from the registered travel texts); the assert proves that exactly the rows whose ids include a not-key carry the label (`S9.check.rule7_label_on_W9_branches_rows`).\n",
  "- **Not labelled.** Holm-adjusted p values, corroboration verdicts and counts are computed by code from keys (deviation D-4). The side tables `corroboration.tsv` and `convention_questions.tsv` are not claim rows: they quote Stage 1b adjudicator and verifier numbers beside their sources and are not labelled.\n",
  "- **Asserts.** `S9.check.rule7_every_travelling_number_classified`, `S9.check.rule7_label_on_every_notkey_row`, `S9.check.rule7_travel_text_ids_classified`, `S9.check.adjudicator_label_on_s6d_values`, `S9.check.rule7_label_on_W9_branches_rows`, `S9.check.rule7_s12_labelled`, `S9.check.rule7_merge_rows_labelled`.\n\n",
  md_table(chr_df(R7T)))
# dispositions (fix 6)
n_cq <- sum(grepl("^CQ", DPT$disposition)); n_mo <- sum(DPT$disposition == "Methods wording only (no consequence on these inputs)"); n_fx <- sum(startsWith(DPT$disposition, "fixed ("))
DPD <- data.frame(check.names = FALSE, `#` = DPT$ord, `Section/lens` = DPT$section_lens, `Reported key` = DPT$reported_key, Defect = DPT$defect,
                  Disposition = DPT$disposition, Detail = DPT$detail)
b_disp <- paste0(
  "Every defect line of `S1R1_R2_1b_verify/VERIFIER_DEFECTS.md` (", nrow(DPT), " lines; assert `S9.check.dispositions_count_equals_defect_lines`) and its disposition. ",
  n_cq, " map to a convention question of the table above, ", n_mo, " are Methods wording only (no consequence on these inputs), ", n_fx, " are fixed in W9. ",
  "The Defect column is the first part of the verifier's text.\n\n", md_table(chr_df(DPD)))
# summary of the opening fixes
b_stage2 <- paste0(
  "Stage 2 opening fixes (STAGE1B_REPORT section 7; RESTART_STAGE2 section B). Generated files change; no numeric value changes. Full record: `FIXLOG.md`.\n\n",
  "1. **Deviations.** The section below lists D-1 to D-4 of STAGE1B_REPORT section 1, verbatim.\n",
  "2. **Visibility.** The comparator's label is corrected (see the note after the Holm table; `stage1b_key_labels.tsv`).\n",
  "3. **Corroboration condition 2.** '", COND2_NA, "' replaces TRUE for ", paste(na_pairs, collapse = ", "), ". Verdicts unchanged.\n",
  "4. **Scope.** 'travels to every surface' is now 'attached to ...' the rows it is actually on (the East Asian scale claim row, W7 row ", ir_eas, "; manuscript surfaces do not exist yet); the section 6(d) sentence states what W9 recounts and what it only reads.\n",
  "5. **Rule 7 amendment** (user decision 12). ", sum(IRd$rule7_label != ""), " of ", sum(ISC), " integrated rows in scope carry not-key quantity ids and the label; ", sum(IRd$rule7_adjudicator_label != ""), " carry the adjudicator label.\n",
  "6. **Defect dispositions.** ", nrow(DPT), " verifier-reported definition defects, each dispositioned.\n",
  "7. **Not done (deferred).** Folding cd-aware tracking into `transcript_audit` (with H2's 7 probes as tests) was NOT done in this pass; it stays an open Stage 2 opening item.")
# deviations (fix 1): STAGE1B_REPORT section 1, verbatim
assert_md5(STAGE1B_REPORT$path, STAGE1B_REPORT$md5, "STAGE1B_REPORT.md")
b_dev <- paste0(
  "The Stage 1b middle path (user, 2026-10-09) declares these deviations. They are reproduced verbatim from `S1R1_R2_1b_verify/STAGE1B_REPORT.md` section 1 (md5 ", substr(STAGE1B_REPORT$md5, 1, 8), "); that report is the record.\n\n",
  deviations_section(paste(readLines(STAGE1B_REPORT$path, encoding = "UTF-8", warn = FALSE), collapse = "\n")), "\n\n",
  "Within W9 itself, every item in the brief was run as declared and nothing is NOT RUN. Scope notes:\n",
  "- W9F reads the comparator `*_SURVIVAL.tsv` tables, as the W9F task directs. For S1, S2, S4, S5 and S7 they are read from `_restricted/compare_1b/`, because no adjudication copy exists.\n",
  "- `verify_SINGLE/verify_w9_part1.py`, the Stage 1 W9 verifier's script, still pins the pre-normaliser W6b and W7 md5s. It is the verifier's record and was not edited.\n",
  "- The D-1 text above says the not-keys travel 'unlabelled'; that was the state before Stage 2. Since the Stage 2 opening fixes they carry the Rule 7 label (user decision 12; section 'Rule 7 amendment' above).")

blocks <- list(holm = md_table(HB), corroboration = md_table(CB), s6d = md_table(DB), conventions = md_table(QB), gates = GB,
               item5 = b_item5, s6dscope = b_s6dscope, visibility = b_visibility, cond2note = b_cond2, stage2fixes = b_stage2, rule7 = b_rule7, dispositions = b_disp, deviations = b_dev,
               s12 = paste0("(a) Pooled two-resource A-B difference (inverse-variance fixed effect, W4a aware 2nd-order SEs):\n\n", md_table(SA),
                            "\n\n(b) Holm on the pre-reclassification p values:\n\n", md_table(SB),
                            "\n\nBoth tables: ", RULE7_LABEL, " (computed by W9 from W4a's not-key SEs and p values; supplement-only; Rule 7 amendment)."))
rp <- file.path(OUT, "README.md")
txt <- paste(readLines(rp, encoding = "UTF-8", warn = FALSE), collapse = "\n")
for (b in names(blocks)) {
  pat <- sprintf("(?s)<!-- BEGIN:%s -->.*?<!-- END:%s -->", b, b)
  n <- lengths(regmatches(txt, gregexpr(pat, txt, perl = TRUE)))
  if (n != 1L) stop(sprintf("README marker block '%s' occurs %d times (expected 1)", b, n))
  rep <- sprintf("<!-- BEGIN:%s -->\n%s\n<!-- END:%s -->", b, blocks[[b]], b)
  txt <- sub(pat, gsub("\\\\", "\\\\\\\\", rep), txt, perl = TRUE)
  logline(sprintf("README block %-14s written (%d chars)", b, nchar(blocks[[b]])))
}
con <- file(rp, open = "wb"); writeBin(charToRaw(enc2utf8(paste0(txt, "\n"))), con); close(con)
logline(sprintf("README.md md5 %s", md5(rp)))
