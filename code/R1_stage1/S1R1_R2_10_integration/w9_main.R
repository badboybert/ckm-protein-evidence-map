# W9 - Stage 1 integration: Holm over the Rule 4 family, the "nominal" placeholders, Principle 0 corroboration,
# Amendment 1 A1 front-surface caveats, the DPEP1 s merge, and the two section-12 items.
# Computed from the redacted pre-declaration and amendment only. The full texts are opened afterwards, by
# w9_compare_seen.R, to write COMPARE_SEEN.tsv. Nothing computed here is changed after that.
# Rerun: Rscript run_all.R   (from any directory)

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
# W9F: mutation mode. W9_MUTATE=<id> (set only by w9_mutation.R) applies one in-memory mutation to a loaded input or to
# an intermediate table (the hooks below), writes asserts.tsv and the log to mutations/run_<id>/, and writes nothing else.
MUT <- Sys.getenv("W9_MUTATE", "")
OUT_W <- if (nzchar(MUT)) file.path(OUT, "mutations", paste0("run_", MUT)) else OUT
if (nzchar(MUT)) dir.create(OUT_W, recursive = TRUE, showWarnings = FALSE)
LOGF <- file.path(OUT_W, Sys.getenv("W9_LOG", "run.log"))
T0 <- Sys.time()
logline("================ W9 w9_main.R start ================")
logline("R: ", R.version.string, " | data.table ", as.character(packageVersion("data.table")),
        " | jsonlite ", as.character(packageVersion("jsonlite")))
logline("seeds: none (W9 runs no bootstrap or simulation)")

ASSERTS <- list()
add_assert <- function(id, ok, derivation) {
  ASSERTS[[length(ASSERTS) + 1L]] <<- data.table(id = id, status = if (isTRUE(ok)) "PASS" else "ASSERT_FAILED", derivation = derivation)
  logline(sprintf("assert %-60s %s  | %s", id, if (isTRUE(ok)) "PASS" else "ASSERT_FAILED", derivation))
  invisible(ok)
}

# ---------------------------------------------------------------- 1. governing documents and inputs
logline("--- governing documents ---")
for (g in names(GOV)) assert_md5(GOV[[g]]$path, GOV[[g]]$md5, g)
resj <- jsonlite::fromJSON(GOV$resources$path, simplifyVector = FALSE)
stopifnot(length(resj) == 324L)   # W9F step 0: re-issued map (546b769a, 324 entries; 311 in the Stage 1 brief)
oc_entry <- Filter(function(e) identical(e$id, "R1_outcome_checks_tsv"), resj)
stopifnot(length(oc_entry) == 1L)
oc_entry <- oc_entry[[1]]
stopifnot(identical(oc_entry$md5, EXTRA$R1_outcome_checks_tsv$md5), identical(oc_entry$role, "seen_output"),
          identical(normalizePath(oc_entry$path, winslash = "/"), normalizePath(EXTRA$R1_outcome_checks_tsv$path, winslash = "/")))
logline("RESOURCES: R1_outcome_checks_tsv role seen_output, md5 ", oc_entry$md5, " (path and md5 agree with the brief) -> READ_FROM_SEEN")
logline("gls_q.R: md5 asserted only; W9 computes no GLS Q, so it is not sourced and its self-test is not part of this run")

logline("--- staging snapshot (start) ---")
snap <- function() {
  f <- list.files(STAGING, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  data.table(file = substring(f, nchar(STAGING) + 2L), bytes = file.size(f), md5 = unname(tools::md5sum(f)))[order(file)]
}
SNAP0 <- snap()
if (!nzchar(MUT)) fwrite(SNAP0, file.path(OUT_W, "staging_snapshot_start.tsv"), sep = "\t")
logline(sprintf("staging snapshot start: %d files, %s bytes", nrow(SNAP0), format(sum(SNAP0$bytes), big.mark = ",")))

logline("--- workstream inputs ---")
R <- list(); B <- list()
for (w in names(WS)) {
  rp <- file.path(S1, WS[[w]]$dir, "results.tsv"); bp <- file.path(S1, WS[[w]]$dir, "branches.tsv")
  assert_md5(rp, WS[[w]]$results, paste0(w, "/results.tsv"))
  assert_md5(bp, WS[[w]]$branches, paste0(w, "/branches.tsv"))
  R[[w]] <- read_tsv_chk(rp, paste0(w, "/results.tsv"))
  B[[w]] <- read_tsv_chk(bp, paste0(w, "/branches.tsv"))
  stopifnot(identical(names(B[[w]]), c("rule_ref", "inputs", "branch", "filled_template", "qualifiers")))
  stopifnot(all(c("quantity_id", "panel", "cov", "se_order", "value", "p", "estimand", "resource") %in% names(R[[w]])))
}
for (x in names(EXTRA)) assert_md5(EXTRA[[x]]$path, EXTRA[[x]]$md5, x)
EXT  <- read_tsv_chk(EXTRA$W1_ext_checks$path, "W1/ext_checks.tsv")
R1OC <- read_tsv_chk(EXTRA$R1_outcome_checks_tsv$path, "R1_outcome_checks_tsv [READ_FROM_SEEN]")
stopifnot(nrow(R1OC) == 7L, nrow(EXT) == 5L)

logline("--- W9F inputs (Stage 1b-C): C1 external checks, comparator copies, adjudication files ---")
for (x in names(W9F_IN)) assert_md5(W9F_IN[[x]]$path, W9F_IN[[x]]$md5, x)
C1R   <- read_tsv_chk(W9F_IN$C1_results$path, "C1/results.tsv")
C1S   <- read_tsv_chk(W9F_IN$C1_survival$path, "C1_SURVIVAL.tsv (comparator copy)")
S3S   <- read_tsv_chk(W9F_IN$S3_survival$path, "S3_SURVIVAL.tsv (comparator copy)")
S6S   <- read_tsv_chk(W9F_IN$S6_survival$path, "S6_SURVIVAL.tsv (comparator copy)")
C1ADJ <- read_tsv_chk(W9F_IN$C1_adj$path, "adj_C1.tsv (C1 adjudicator)")
S3SC  <- read_tsv_chk(W9F_IN$S3_scen$path, "adj_S3_scenarios.tsv (S3 adjudicator)")
S3PK  <- read_tsv_chk(W9F_IN$S3_pkg$path, "adj_S3_coloc_pkg.tsv (S3 adjudicator)")
A6D   <- read_tsv_chk(W9F_IN$S6_d_table$path, "a3_d_table.tsv (S6 adjudicator)")
A6F   <- read_tsv_chk(W9F_IN$S6_fitset$path, "a3_allcand_fitset.tsv (S6 adjudicator)")
A6N   <- read_tsv_chk(W9F_IN$S6_notch2$path, "a3_notch2_fits.tsv (S6 adjudicator)")
A6K   <- read_tsv_chk(W9F_IN$S6_checks$path, "a3_checks.tsv (S6 adjudicator)")
A6U   <- read_tsv_chk(W9F_IN$S6_umod_q$path, "a3_umod_q.tsv (S6 adjudicator)")
TPMI_J <- jsonlite::fromJSON(W9F_IN$C1_tpmi_json$path)
S3J   <- read_tsv_chk(W9F_IN$S3_join$path, "adj_S3_join_audit.tsv (S3 adjudicator)")
SURV  <- rbindlist(c(list(C1 = C1S, S3 = S3S, S6 = S6S),
                     setNames(lapply(c("S1", "S2", "S4", "S5", "S7"), function(s) read_tsv_chk(W9F_IN[[paste0(s, "_survival")]]$path, paste0(s, "_SURVIVAL.tsv (comparator)"))),
                              c("S1", "S2", "S4", "S5", "S7"))), idcol = "sec", fill = TRUE)
stopifnot(!anyDuplicated(SURV$key_id))
sv <- function(key) { r <- SURV[SURV$key_id == key]; if (nrow(r) != 1L) stop("comparator key ", key, ": ", nrow(r), " rows"); r }

# ---- Stage 2 opening fix 5: key / not-key classification of quantity ids (Rule 7 amendment) ----
# KEYQ = quantity ids that are headline keys (KEY_MAP.tsv, single-id rows: the 104 numeric keys; categorical keys are branch or class levels, not numbers); NKQ = NOT_KEYS.tsv ids.
KEYQ <- character(); NKQ <- character()
for (s_ in names(KM_MD5)) {
  km <- read_tsv_chk(W9F_IN[[paste0("keymap_", s_, "_key")]]$path, paste0("KEY_MAP ", s_))
  qc <- grep("^quantity_id", names(km), value = TRUE); stopifnot(length(qc) == 1L)
  ids_ <- km[[qc]]; KEYQ <- c(KEYQ, ids_[grepl("^(C1|S[1-7])\\.[^ ]+$", ids_)])
  NKQ <- c(NKQ, read_tsv_chk(W9F_IN[[paste0("keymap_", s_, "_notkey")]]$path, paste0("NOT_KEYS ", s_))$quantity_id)
}
stopifnot(!anyDuplicated(KEYQ), length(intersect(KEYQ, NKQ)) == 0L)
IDMAP <- lapply(c(W3 = "idmap_W3", W4a = "idmap_W4a", W5 = "idmap_W5"), function(x) read_tsv_chk(W9F_IN[[x]]$path, x))
# ---- W9F mutation hooks (w9_mutation.R sets W9_MUTATE; each mutation changes exactly one cell of a loaded table) ----
MUT_APPLIED <- FALSE
mut_cell <- function(dt, sel, col, new) {
  if (sum(sel) != 1L) stop(sprintf("MUTATION %s: target not unique (%d rows)", MUT, sum(sel)))
  if (identical(dt[[col]][which(sel)], new)) stop(sprintf("MUTATION %s: the new value equals the old one", MUT))
  logline(sprintf("MUTATION %s: %s '%s' -> '%s'", MUT, col, substr(dt[[col]][which(sel)], 1, 80), substr(new, 1, 80)))
  data.table::set(dt, which(sel), col, new); MUT_APPLIED <<- TRUE; invisible(TRUE)
}
qrow <- function(dt, q) dt$quantity_id == q
if (nzchar(MUT)) switch(MUT,
  c1_status   = mut_cell(C1S, C1S$key_id == "C1H.RAD23B_T2D.class", "status", "L1L2_AGREE_IMPL_DIFFERS"),
  c1_L2       = mut_cell(C1S, C1S$key_id == "C1H.RAD23B_T2D.class", "L2", "inconclusive"),
  c1_value    = mut_cell(C1R, qrow(C1R, "C1.ARG1_T2D.class"), "value", "inconclusive"),
  c1_vs_w1    = mut_cell(EXT, EXT$pair == "PAM_T2D", "class_resource_p", "discordant"),
  rad_interval = mut_cell(C1ADJ, C1ADJ$pair == "RAD23B_T2D", "p_wald1_rounding_interval", "0.0361-0.0490"),
  rad_bonf    = mut_cell(C1ADJ, C1ADJ$pair == "RAD23B_T2D", "bonf6_W", "True"),
  a2_share    = mut_cell(R$W4a, qrow(R$W4a, "S4.ld.boot.FinnGen.cond.aware2nd.share_same"), "value", "0.95"),
  a2_frag     = mut_cell(B$W4a, B$W4a$rule_ref == "Amendment 1 A2", "qualifiers",
                         sub("in FinnGen", "in DIAMANTE", B$W4a$qualifiers[B$W4a$rule_ref == "Amendment 1 A2"], fixed = TRUE)),
  cond2       = mut_cell(S3PK, S3PK$pair == "PAM_T2D" & S3PK$maf_mode == "sub" & S3PK$reading == "B_full_exposure", "branch", "2_addon"),
  s3_pp       = mut_cell(S3SC, S3SC$pair == "PAM_T2D" & S3SC$maf_mode == "sub" & S3SC$scenario == "cap_both_fullwindows", "H3_6", "0.18"),
  n2_files    = mut_cell(A6N, A6N$reading == "keep_shift_rows" & A6N$fit == "p1fr@489", "sP", "0.2"),
  w6b_smax    = mut_cell(R$W6b, qrow(R$W6b, "S6.d.GSTA1_eGFR.s_max"), "value", "0.2100"),
  w6b_smax4   = mut_cell(R$W6b, qrow(R$W6b, "S6.d.GSTA1_eGFR.s_max_4fit"), "value", "0.2100"),
  recount     = mut_cell(A6D, A6D$fit_set == "four" & A6D$shift_rows == "keep" & A6D$version == "reclassify_ING1", "n_reclassified", "12"),
  invariant   = mut_cell(A6D, A6D$fit_set == "eight" & A6D$shift_rows == "keep" & A6D$version == "reclassify_ING1", "lo", "1.2"),
  l1_smax     = mut_cell(S6S, S6S$key_id == "S6H.d.UMOD_eGFR.s_max", "L1", "0.2972914636"),
  impl_out    = mut_cell(R$W6b, qrow(R$W6b, "S6.d.a"), "value", "3"),
  n2_cross    = mut_cell(A6D, A6D$fit_set == "four" & A6D$shift_rows == "keep" & A6D$version == "reclassify_ING1", "NOTCH2_s_max", "0.1001"),
  omni_p      = mut_cell(R$W4b, qrow(R$W4b, "S4.omni.T2D_DM.6sig.p_aware_2nd"), "value", "0.2020789738"),
  drop_ver    = mut_cell(R$W3, qrow(R$W3, "S3.RAD23B_T2D.pmafdrop.1e-05.H3"), "value", "0.193189"),
  fam_key     = mut_cell(SURV, SURV$key_id == "S5H.Q.egfr.p", "impl", "0.12"),
  V_vis       = mut_cell(SURV, SURV$key_id == "S4H.arg.p", "status", "SURVIVES"),
  NULL)
mut_late <- function(point, expr) if (nzchar(MUT) && startsWith(MUT, point)) force(expr)   # hooks on intermediate tables
stopifnot(nrow(C1ADJ) == 6L, nrow(A6D) == 8L, nrow(A6F) == 16L * 7L, nrow(A6N) == 14L, nrow(A6U) == 4L, all(A6K$pass == "True"))

# ---- Stage 2 opening fix 2: correct the comparator's SURVIVES_VISIBLE_TARGET label (H2 D2) ----
# The comparator label is a lexical token match: it over-flags. Only the 5 keys of VIS_GENUINE had the answer printed in a file the verifiers
# could open; VIS_DERIVABLE are closed-form functions of those counts; every other key so labelled is an independent re-derivation.
SURV[, status_comparator := status]
vis_raw <- SURV$key_id[SURV$status_comparator == "SURVIVES_VISIBLE_TARGET"]
SURV[, status := ifelse(status_comparator == "SURVIVES_VISIBLE_TARGET",
                        ifelse(key_id %in% VIS_GENUINE, VIS_LABEL_GENUINE, VIS_LABEL_INDEP),
                        ifelse(key_id %in% VIS_DERIVABLE, VIS_LABEL_DERIVABLE, status_comparator))]
# the three adjudication-copy tables (C1, S3, S6) are read separately by later sections: correct their visible-target labels the same way (rows with other statuses, incl. mutated ones, are left as read)
for (tb in list(C1S, S3S, S6S)) tb[status == "SURVIVES_VISIBLE_TARGET", status := SURV$status[match(key_id, SURV$key_id)]]
add_assert("S9.check.visible_label_corrected",
           length(VIS_KEYS8) == 8L && !anyDuplicated(VIS_KEYS8) && length(vis_raw) == 67L && all(VIS_GENUINE %in% vis_raw) &&
             all(SURV$status_comparator[SURV$key_id %in% VIS_DERIVABLE] == "SURVIVES") &&
             setequal(SURV$key_id[startsWith(SURV$status, "SURVIVES_VISIBLE_TARGET")], VIS_GENUINE) &&
             setequal(SURV$key_id[SURV$status == VIS_LABEL_DERIVABLE], VIS_DERIVABLE) &&
             sum(SURV$status == VIS_LABEL_INDEP) == length(vis_raw) - length(VIS_GENUINE) && !any(c(C1S$status, S3S$status, S6S$status) == "SURVIVES_VISIBLE_TARGET"),
           sprintf("comparator labelled %d of %d keys SURVIVES_VISIBLE_TARGET; genuinely visible: %s; derivable from visible counts: %s; the other %d are independent re-derivations (token false positives, H2 D2); 8 key ids from the constant VIS_KEYS8",
                   length(vis_raw), nrow(SURV), paste(VIS_GENUINE, collapse = ", "), paste(VIS_DERIVABLE, collapse = ", "), sum(SURV$status == VIS_LABEL_INDEP)))
assert_md5(STAGE1B_REPORT$path, STAGE1B_REPORT$md5, "STAGE1B_REPORT.md")
RV <- report_vis_keys(paste(readLines(STAGE1B_REPORT$path, encoding = "UTF-8", warn = FALSE), collapse = "\n"))
mut_late("V_rep", { RV$genuine <- setdiff(RV$genuine, "S6H.primary.p"); MUT_APPLIED <- TRUE })
add_assert("S9.check.visible_keys_match_STAGE1B_REPORT", setequal(RV$genuine, VIS_GENUINE) && setequal(RV$derivable, VIS_DERIVABLE) && length(RV$genuine) == 5L && length(RV$derivable) == 3L,
           sprintf("STAGE1B_REPORT section 3 (md5 %s) lists as genuinely visible: %s; derivable: %s; both equal the VIS_GENUINE / VIS_DERIVABLE constants", substr(STAGE1B_REPORT$md5, 1, 8),
                   paste(RV$genuine, collapse = ", "), paste(RV$derivable, collapse = ", ")))
# Stage 2 step 2 fix round: the Stage 1b status of an S6H.d.<cand>.s_max key was recorded against the earlier FOUR-fit implementer value; W6b now prints the SEVEN-fit s_max.
S6D_SMAX_QUAL <- paste0("status was compared against the earlier four-fit implementer value, not the printed seven-fit s_max; the verifier lenses L1 and L2 already hold seven-fit values for this key that match the printed seven-fit s_max within ", S6D_SMAX_TOL_TXT,
  " (S6_SURVIVAL.tsv columns L1, L2; assert S9.check.s6d_smax_L1L2_hold_seven_fit); where the status is L1L2_AGREE_IMPL_DIFFERS it is stale only because the implementer value compared against was the four-fit one; a SURVIVES label was compared against the four-fit value; not re-labelled")
s6d_k <- S6S$key_id[grepl("^S6H[.]d[.].+[.]s_max$", S6S$key_id) & !is.na(suppressWarnings(as.numeric(S6S$L1)))]
s6d_dev <- vapply(s6d_k, function(k) { w <- as.numeric(R$W6b$value[R$W6b$quantity_id == sub("^S6H[.]", "S6.", k)]); stopifnot(length(w) == 1L, is.finite(w))
  max(abs(as.numeric(S6S$L1[S6S$key_id == k]) - w), abs(as.numeric(S6S$L2[S6S$key_id == k]) - w)) }, 0)
add_assert("S9.check.s6d_smax_L1L2_hold_seven_fit", length(s6d_k) == 15L && max(s6d_dev) <= S6D_SMAX_TOL && max(s6d_dev) > S6D_SMAX_TOL / 1.1 && sum(S6S$status[match(s6d_k, S6S$key_id)] == "L1L2_AGREE_IMPL_DIFFERS") == 6L,
           sprintf("%d computable S6H.d.<cand>.s_max keys: max |L1 or L2 - W6b printed seven-fit s_max| = %.3g <= %s; %d of them carry the Stage 1b status L1L2_AGREE_IMPL_DIFFERS (compared against the four-fit implementer value)",
                   length(s6d_k), max(s6d_dev), S6D_SMAX_TOL_TXT, sum(S6S$status[match(s6d_k, S6S$key_id)] == "L1L2_AGREE_IMPL_DIFFERS")))
if (!nzchar(MUT)) fwrite(SURV[, .(key_id, section = sec, status_comparator, status_corrected = status,
                                  qualifier = ifelse(grepl("^S6H[.]d[.].+[.]s_max$", key_id) & status != "MISSING", S6D_SMAX_QUAL, ""))],
                         file.path(OUT_W, "stage1b_key_labels.tsv"), sep = "\t", quote = TRUE)

RES_ROWS <- list()
addr <- function(quantity_id, prespec_ref, estimand, resource = "", panel = "", cov = "", se_order = "",
                 value, se = "", ci_lo = "", ci_hi = "", p = "", n = "", unit = "", notes = "") {
  fmt <- function(x) if (is.numeric(x)) format(x, digits = 10, scientific = FALSE, drop0trailing = TRUE, trim = TRUE) else as.character(x)
  fmtp <- function(x) if (is.numeric(x)) format(x, digits = 10, trim = TRUE) else as.character(x)
  RES_ROWS[[length(RES_ROWS) + 1L]] <<- data.table(quantity_id = quantity_id, prespec_ref = prespec_ref, estimand = estimand,
    resource = resource, panel = panel, cov = cov, se_order = se_order, value = fmtp(value), se = fmt(se),
    ci_lo = fmtp(ci_lo), ci_hi = fmtp(ci_hi), p = fmtp(p), n = fmt(n), unit = unit, notes = notes)
}
f2 <- function(x) formatC(x, digits = 2, format = "g", flag = "#")   # 2 significant figures, as printed in branches
f3 <- function(x) sprintf("%.3f", x)

# ---------------------------------------------------------------- 2. Rule 4 family and Holm
logline("--- Rule 4 family: p_aware_2nd at each named primary cell ---")
fam <- data.table(
  member = c("arg1_arginine", "arg1_t2d", "arg1_test3", "arg1_ornithine", "tgfb1_cad", "dpep1_cysgly", "dpep1_egfr", "umod_egfr"),
  label  = c("ARG1-arginine A vs B", "ARG1-T2D A vs B", "ARG1 test (iii)", "ARG1-ornithine A vs B", "TGFB1-CAD (FinnGen CHD)",
             "DPEP1-cysteinylglycine (all signals)", "DPEP1-eGFR (all signals)", "UMOD-eGFR (all signals)"),
  ws     = c("W4a", "W4a", "W4a", "W4a", "W2", "W5", "W5", "W6b"),
  qid    = c("S4.arg.marg.p_aware_2nd", "S4.grid.FinnGen.cond.p_aware_2nd", "S4.t3.FinnGen.marg.503.2nd.p_aware|S4.t3.FinnGen.cond.503.2nd.p_aware",
             "S4.orn.503.p_aware_2nd", "S2.het.finngen_chd.p_aware_2nd", "S5.Q.cysgly.all5.503.aware_2nd.p",
             "S5.Q.eGFR.all5.503.aware_2nd.p", "S6.umod.Q.aware_2nd.p"),
  cell   = c("marginal, aware, 2nd order, 503 (Rule 3(ii) named primary cell)",
             "FinnGen R10 conditional, aware, 2nd order, 503",
             "FinnGen R10, aware, 2nd order, 503; entered with the larger of the conditional and marginal p",
             "marginal, aware, 2nd order, 503",
             "FinnGen CHD conditional, aware, 2nd order, 503",
             "all 5 signals, GLS Q, aware, 2nd order, 503",
             "all 5 signals, GLS Q, aware, 2nd order, 503",
             "all 5 signals, GLS Q, aware, 2nd order, 503"),
  estimand_pat = c("marginal", "conditional", "test \\(iii\\)", "marginal", "conditional ratios", "", "", "all 5 signals|5 conditionally distinct"),
  resource_pat = c("arginine", "FinnGen", "FinnGen", "ornithine", "finngen_chd", "GCST90200364", "CKDGen", "CKDGen"))

conform <- function(w, q, estpat, respat) {
  row <- R[[w]][R[[w]]$quantity_id == q]
  ok <- nrow(row) == 1L && row$panel == "503" && row$cov == "aware" && row$se_order == "2nd" &&
        (estpat == "" || grepl(estpat, row$estimand)) && grepl(respat, row$resource)
  add_assert(paste0("S9.conform.", q), ok, sprintf("%s %s: panel=%s cov=%s se_order=%s; estimand~'%s'; resource~'%s'",
             w, q, row$panel, row$cov, row$se_order, estpat, respat))
}
fam[, p_unadj := NA_real_]; fam[, p_note := ""]
for (i in seq_len(nrow(fam))) {
  if (fam$member[i] == "arg1_test3") {
    qs <- strsplit(fam$qid[i], "|", fixed = TRUE)[[1]]
    pm <- getq(R$W4a, "W4a", qs[1]); pc <- getq(R$W4a, "W4a", qs[2])
    for (q in qs) conform("W4a", q, fam$estimand_pat[i], fam$resource_pat[i])
    fam$p_unadj[i] <- max(pm, pc)
    fam$p_note[i] <- sprintf("conditional %s, marginal %s; larger = %s", format(pc, digits = 6), format(pm, digits = 6),
                             if (pm >= pc) "marginal" else "conditional")
    T3 <- list(p_cond = pc, p_marg = pm)
  } else {
    fam$p_unadj[i] <- getq(R[[fam$ws[i]]], fam$ws[i], fam$qid[i])
    conform(fam$ws[i], fam$qid[i], fam$estimand_pat[i], fam$resource_pat[i])
  }
}
stopifnot(all(fam$p_unadj > 0 & fam$p_unadj <= 1), nrow(fam) == 8L, !anyDuplicated(fam$member))
# W9F: the Stage 1b headline key behind each member's p, and its comparator status (identity asserted by value)
fam[, key_1b := c("S4H.arg.p", "S4H.t2d.p", "S4H.t3.p_cond|S4H.t3.p_marg", "S4H.orn.p", "S2H.het_p.finngen", "S5H.Q.cysgly.p", "S5H.Q.egfr.p", "S6H.umod.Q.p")]
fam[, status_1b := ""]
for (i in seq_len(nrow(fam))) {
  ks <- strsplit(fam$key_1b[i], "|", fixed = TRUE)[[1]]
  rows <- lapply(ks, sv)
  vals <- if (fam$member[i] == "arg1_test3") c(T3$p_cond, T3$p_marg) else fam$p_unadj[i]
  same <- all(mapply(function(r, v) abs(as.numeric(r$impl) - v) <= 1e-6 * abs(v), rows, vals))
  fam$status_1b[i] <- paste(vapply(rows, function(r) sprintf("%s %s%s (L1 %s, L2 %s)", r$key_id, r$status, if (nzchar(r$alt_flag)) paste0(" ", r$alt_flag) else "",
                                                             format(as.numeric(r$L1), digits = 4), format(as.numeric(r$L2), digits = 4)), ""), collapse = "; ")
  add_assert(paste0("S9.check.family_input_is_stage1b_key.", fam$member[i]), same,
             sprintf("W9 input p %s vs comparator impl %s: %s", paste(format(vals, digits = 7), collapse = "/"), paste(vapply(rows, `[[`, "", "impl"), collapse = "/"), fam$status_1b[i]))
}
fam[, p_holm := holm(p_unadj)]
fam[, rank := rank(p_unadj, ties.method = "first")]
# W9F item 6: nominal_code / nominal_attach = the LITERAL reading (default); *_I1 = the applied reading, printed beside
fam[, c("nominal_code", "nominal_attach", "nominal_code_I1", "nominal_attach_I1") := {
  rl <- lapply(seq_len(.N), function(i) nominal_literal(p_unadj[i], p_holm[i]))
  ri <- lapply(seq_len(.N), function(i) nominal_I1(p_unadj[i], p_holm[i]))
  list(vapply(rl, `[[`, "", "code"), vapply(rl, `[[`, TRUE, "attach"), vapply(ri, `[[`, "", "code"), vapply(ri, `[[`, TRUE, "attach"))
}]
for (i in seq_len(nrow(fam))) logline(sprintf("family %-15s p %-12s Holm %-12s rank %d  literal %-24s | I1 %s", fam$member[i],
  format(fam$p_unadj[i], digits = 6), format(fam$p_holm[i], digits = 6), fam$rank[i], fam$nominal_code[i], fam$nominal_code_I1[i]))
mut_late("F_nom_default", { fam[fam$member == "umod_egfr", nominal_attach := FALSE]; MUT_APPLIED <- TRUE; logline("MUTATION F_nom_default: umod_egfr literal attach -> FALSE") })
add_assert("S9.check.nominal_literal_is_default", all(fam$nominal_attach == (fam$p_holm >= 0.05)) &&
             all(fam$nominal_attach_I1 == (fam$p_unadj < 0.05 & fam$p_holm >= 0.05)),
           sprintf("literal 'nominal' attached to %d of 8 members (every Holm p >= 0.05); I1 attaches it to %d (unadjusted p < 0.05 and Holm p >= 0.05)",
                   sum(fam$nominal_attach), sum(fam$nominal_attach_I1)))

# cross-check: the unadjusted p printed in W4a's "Rule 4 family entry" rows
fe <- B$W4a[B$W4a$rule_ref == "Rule 4 family entry"]
stopifnot(nrow(fe) == 4L)
fe_map <- c("ARG1-arginine A vs B" = "arg1_arginine", "ARG1-T2D A vs B (FinnGen conditional aware)" = "arg1_t2d",
            "ARG1 test (iii) (larger of conditional and marginal)" = "arg1_test3", "ARG1-ornithine A vs B" = "arg1_ornithine")
stopifnot(all(fe$branch %in% names(fe_map)))
for (k in seq_len(nrow(fe))) {
  pstr <- sub("^unadjusted p ", "", fe$inputs[k]); printed <- as.numeric(pstr)
  nsig <- nchar(gsub("^0+", "", gsub("[^0-9]", "", sub("e.*$", "", pstr))))   # significant digits as printed
  mine <- fam$p_unadj[fam$member == fe_map[[fe$branch[k]]]]
  add_assert(paste0("S9.check.W4a_family_entry_p.", fe_map[[fe$branch[k]]]), isTRUE(all.equal(signif(mine, nsig), printed, tolerance = 1e-12)),
             sprintf("W4a printed '%s' (%d significant digits) vs W9 lookup %s rounded to the same precision = %s", fe$inputs[k], nsig,
                     format(mine, digits = 6), format(signif(mine, nsig))))
}

for (i in seq_len(nrow(fam))) {
  m <- fam$member[i]
  base_ref <- "Rule 4 (Holm over the 8-member family); Rule 5 primary p"
  addr(paste0("S9.holm.", m, ".p_unadj"), base_ref, paste0(fam$label[i], "; ", fam$cell[i]), fam$ws[i], "503", "aware", "2nd",
       fam$p_unadj[i], p = fam$p_unadj[i], n = 8L,
       notes = paste0("from ", fam$ws[i], "/results.tsv (md5 ", WS[[fam$ws[i]]]$results, ") quantity ", fam$qid[i],
                      if (nzchar(fam$p_note[i])) paste0("; ", fam$p_note[i]) else ""))
  addr(paste0("S9.holm.", m, ".p_holm"), base_ref, paste0(fam$label[i], "; Holm-adjusted over m = 8"), fam$ws[i], "503", "aware", "2nd",
       fam$p_holm[i], p = fam$p_holm[i], n = 8L, notes = "p.adjust(method = 'holm'), cross-checked against a manual step-down implementation")
  addr(paste0("S9.holm.", m, ".rank"), base_ref, fam$label[i], fam$ws[i], value = fam$rank[i], n = 8L, notes = "rank of the unadjusted p (1 = smallest)")
  addr(paste0("S9.holm.", m, ".nominal"), "Rule 4 'nominal' (literal reading; default column, W9F item 6)", fam$label[i], fam$ws[i],
       value = fam$nominal_code[i],
       notes = switch(fam$nominal_code[i],
         NOMINAL = "Holm p >= 0.05: the statement carries 'nominal' (literal reading of Rule 4's second sentence)",
         "NOMINAL (non-rejection)" = "Holm p >= 0.05: the statement carries 'nominal' under the literal reading, although it is a non-rejection worded on the unadjusted p (Rule 4, third sentence); the applied reading I1 (S9.holm.*.nominal_I1) attaches none: the user decides",
         SURVIVES_HOLM = "Holm p < 0.05: no 'nominal'"))
  addr(paste0("S9.holm.", m, ".nominal_I1"), "Rule 4 'nominal' (applied reading I1; printed beside for the user, W9F item 6)", fam$label[i], fam$ws[i],
       value = fam$nominal_code_I1[i],
       notes = switch(fam$nominal_code_I1[i],
         NOMINAL = "unadjusted p < 0.05 and Holm p >= 0.05: 'nominal'",
         SURVIVES_HOLM = "unadjusted p < 0.05 and Holm p < 0.05: no 'nominal'",
         NON_REJECTION = "unadjusted p >= 0.05: a non-rejection worded on the unadjusted p (Rule 4); I1 attaches no 'nominal'"))
  addr(paste0("S9.holm.", m, ".stage1b_status"), "Stage 1b comparator (W9F)", paste0(fam$label[i], ": Stage 1b headline key status"), fam$ws[i],
       value = fam$status_1b[i], notes = "W9's input equals the comparator's impl value (asserted); status label corrected at Stage 2 (H2 D2): SURVIVES_VISIBLE_TARGET is kept only for S6H.primary.a, .b, .c, .d, .p (answer printed in a file the verifiers could open); S6H.primary.or/lo/hi are derivable from visible counts; every other key the comparator labelled so is an independent re-derivation (token false positive)")
}
addr("S9.holm.m", "Rule 4", "number of family members", value = 8L, n = 8L, notes = "the 8 members listed in the redacted Rule 4 and the W9 brief")
addr("S9.holm.arg1_test3.p_cond", "Rule 4; 4(iii)", "test (iii) conditional p (FinnGen aware 2nd 503)", "W4a", "503", "aware", "2nd", T3$p_cond, p = T3$p_cond,
     notes = "S4.t3.FinnGen.cond.503.2nd.p_aware")
addr("S9.holm.arg1_test3.p_marg", "Rule 4; 4(iii)", "test (iii) marginal p (FinnGen aware 2nd 503)", "W4a", "503", "aware", "2nd", T3$p_marg, p = T3$p_marg,
     notes = "S4.t3.FinnGen.marg.503.2nd.p_aware")

# ---------------------------------------------------------------- 3. re-derivations from workstream results (W9's own checks)
logline("--- re-derivation checks ---")
# section 2 branch
s2v <- list(
  rA_fg = getq(R$W2, "W2", "S2.ratio.finngen_chd.rs73045269.value"), rB_fg = getq(R$W2, "W2", "S2.ratio.finngen_chd.rs12983775.value"),
  p_fg  = getq(R$W2, "W2", "S2.het.finngen_chd.p_aware_2nd"),
  rA_ieu = getq(R$W2, "W2", "S2.ratio.ieu_a_7.rs73045269.value"), rB_ieu = getq(R$W2, "W2", "S2.ratio.ieu_a_7.rs12983775.value"),
  p_ieu = getq(R$W2, "W2", "S2.het.ieu_a_7.p_aware_2nd"), disc = getq(R$W2, "W2", "S2.discovery_sign"))
S2B <- s2_branch(s2v$rA_fg, s2v$rB_fg, s2v$p_fg, s2v$rA_ieu, s2v$rB_ieu, s2v$p_ieu, s2v$disc)
w2_main <- B$W2[grepl("^prespec s2 wording rule", B$W2$rule_ref)]
stopifnot(nrow(w2_main) == 1L)
add_assert("S9.check.s2_branch_matches_W2", identical(S2B, w2_main$branch),
           sprintf("re-derived %s from ratios FG %+.3f/%+.3f (p %s) and ieu-a-7 %+.3f/%+.3f (p %s), discovery sign %+d; W2 branch %s",
                   S2B, s2v$rA_fg, s2v$rB_fg, f2(s2v$p_fg), s2v$rA_ieu, s2v$rB_ieu, f2(s2v$p_ieu), as.integer(s2v$disc), w2_main$branch))

# section 3 branches (primary input p + MAF + N)
S3P <- c("PAM_T2D", "RAD23B_T2D", "TGFB1_CAD", "PCSK9_CAD")
S3 <- list()
w3_main <- B$W3[grepl("branches 1-4", B$W3$rule_ref)]
stopifnot(nrow(w3_main) == 4L)
for (pr in S3P) {
  h <- sapply(c("1e-05.H4", "1e-05.H3", "1e-06.H4", "1e-06.H3"), function(s) getq(R$W3, "W3", paste0("S3.", pr, ".pmaf.", s)))
  br <- s3_branch(h[1], h[2], h[3], h[4])
  wrow <- w3_main[grepl(paste0("^", pr, " primary input"), w3_main$inputs)]
  stopifnot(nrow(wrow) == 1L)
  S3[[pr]] <- list(branch = br, h = h, w3 = wrow$branch, template = wrow$filled_template)
  add_assert(paste0("S9.check.s3_branch_matches_W3.", pr), identical(br, wrow$branch),
             sprintf("PP.H4/PP.H3 1e-5 %.4f/%.4f, 1e-6 %.4f/%.4f -> branch %s; W3 branch %s", h[1], h[2], h[3], h[4], br, wrow$branch))
}
# W9F step 2b (D1, D2): every section 3 scenario the S3 adjudicator ran, from BOTH of its files (adj_S3_coloc_pkg.tsv
# readings; adj_S3_scenarios.tsv scenarios), every MAF mode, with the files' Rule 1 branch and PP.H4-only branch. W9
# re-derives both branches from the files' PP.H3/PP.H4 with its own s3_branch() before using them.
S3ADJ <- rbind(
  S3PK[, .(file = "adj_S3_coloc_pkg.tsv", pair, maf_mode, scenario = reading, n_desc = sprintf("n_ds1 %s, n_ds2 %s", n_ds1, n_ds2),
           p1 = as.numeric(p1_used), p2 = as.numeric(p2_used), H4_5 = as.numeric(H4_5), H3_5 = as.numeric(H3_5),
           H4_6 = as.numeric(H4_6), H3_6 = as.numeric(H3_6), branch, branch_alt_h4only)],
  S3SC[, .(file = "adj_S3_scenarios.tsv", pair, maf_mode, scenario, n_desc = sprintf("n_exp_rows %s, n_out_window %s", n_exp_rows, n_out_window),
           p1 = as.numeric(p1), p2 = as.numeric(p2), H4_5 = as.numeric(H4_5), H3_5 = as.numeric(H3_5),
           H4_6 = as.numeric(H4_6), H3_6 = as.numeric(H3_6), branch, branch_alt_h4only)])
stopifnot(all(is.finite(c(S3ADJ$p1, S3ADJ$p2, S3ADJ$H4_5, S3ADJ$H3_5, S3ADJ$H4_6, S3ADJ$H3_6))))
s3_rb <- s3_file_code(mapply(s3_branch, S3ADJ$H4_5, S3ADJ$H3_5, S3ADJ$H4_6, S3ADJ$H3_6))
s3_rh <- s3_file_code(mapply(s3_branch, S3ADJ$H4_5, S3ADJ$H3_5, S3ADJ$H4_6, S3ADJ$H3_6, MoreArgs = list(h4only_at_1e6 = TRUE)))
s3_both <- all(vapply(S3P, function(p) all(c("adj_S3_coloc_pkg.tsv", "adj_S3_scenarios.tsv") %in% S3ADJ$file[S3ADJ$pair == p]), TRUE))
add_assert("S9.check.s3_adj_branches_rederived",
           nrow(S3ADJ) == nrow(S3PK) + nrow(S3SC) && s3_both && all(s3_rb == S3ADJ$branch) && all(s3_rh == S3ADJ$branch_alt_h4only),
           sprintf("%d rows (adj_S3_coloc_pkg.tsv %d, adj_S3_scenarios.tsv %d; each of %s in both files: %s): the Rule 1 branch and the PP.H4-only branch re-derived by W9 from each row's PP.H3/PP.H4 equal the file's 'branch' in %d and 'branch_alt_h4only' in %d rows",
                   nrow(S3ADJ), nrow(S3PK), nrow(S3SC), paste(S3P, collapse = ", "), s3_both, sum(s3_rb == S3ADJ$branch), sum(s3_rh == S3ADJ$branch_alt_h4only)))
mut_late("S_scen", mut_cell(S3ADJ, S3ADJ$file == "adj_S3_scenarios.tsv" & S3ADJ$pair == "PAM_T2D" & S3ADJ$maf_mode == "sub" & S3ADJ$scenario == "cap_both_fullwindows",
                            "branch", "2_addon"))

# Amendment 1 A2
a2 <- list(p_fg = getq(R$W4a, "W4a", "S4.grid.FinnGen.cond.p_aware_2nd"), d_fg = getq(R$W4a, "W4a", "S4.grid.FinnGen.cond.diff"),
           p_dm = getq(R$W4a, "W4a", "S4.grid.DIAMANTE.marg.p_aware_2nd"), d_dm = getq(R$W4a, "W4a", "S4.grid.DIAMANTE.marg.diff"))
A2_MET <- a2$p_fg < 0.05 && a2$p_dm < 0.05 && sign(a2$d_fg) == sign(a2$d_dm)
w4a_a2 <- B$W4a[B$W4a$rule_ref == "Amendment 1 A2"]
stopifnot(nrow(w4a_a2) == 1L)
add_assert("S9.check.A2_matches_W4a", A2_MET == grepl("^met", w4a_a2$branch),
           sprintf("FinnGen cond p %s d %+.4f; DIAMANTE marg p %s d %+.4f -> A2 %s; W4a branch '%s'", f2(a2$p_fg), a2$d_fg, f2(a2$p_dm), a2$d_dm,
                   if (A2_MET) "met" else "not met", w4a_a2$branch))
# W9F item 4: which resource(s) fail Rule 3(ii) in the two A2 cells (named in the merged T2D claim's caveat)
w4a_r3 <- function(pat) { r <- B$W4a[startsWith(B$W4a$rule_ref, pat)]; stopifnot(nrow(r) == 1L); !grepl("NOT robust", r$branch, fixed = TRUE) }
A2R <- chk_a2_failing(
  w4a_robust = c(FinnGen = w4a_r3("Rule 3(ii) T2D A vs B (FinnGen R10 conditional, aware)"), DIAMANTE = w4a_r3("Rule 3(ii) A2: T2D A vs B DIAMANTE marginal aware")),
  share = c(FinnGen = getq(R$W4a, "W4a", "S4.ld.boot.FinnGen.cond.aware2nd.share_same"), DIAMANTE = getq(R$W4a, "W4a", "S4.ld.boot.DIAMANTE.marg.aware2nd.share_same")),
  p503 = c(FinnGen = a2$p_fg, DIAMANTE = a2$p_dm),
  p489 = c(FinnGen = getq(R$W4a, "W4a", "S4.grid.FinnGen.cond.489.p_aware_2nd"), DIAMANTE = getq(R$W4a, "W4a", "S4.grid.DIAMANTE.marg.489.p_aware_2nd")),
  d503 = c(FinnGen = a2$d_fg, DIAMANTE = a2$d_dm),
  d489 = c(FinnGen = getq(R$W4a, "W4a", "S4.grid.FinnGen.cond.489.diff"), DIAMANTE = getq(R$W4a, "W4a", "S4.grid.DIAMANTE.marg.489.diff")))
add_assert("S9.check.A2_failing_resources_rederived", A2R$ok, A2R$detail)
A2_FAILS <- A2R$fails
A2_CAVEAT <- if (length(A2_FAILS)) paste0("sensitive to the reference LD in ", paste(A2_FAILS, collapse = " and ")) else ""
add_assert("S9.check.A2_caveat_matches_W4a_fragment", !nzchar(A2_CAVEAT) || grepl(A2_CAVEAT, w4a_a2$qualifiers, fixed = TRUE),
           sprintf("W9 caveat '%s' vs W4a A2 row qualifiers '%s'", A2_CAVEAT, w4a_a2$qualifiers))

# W9F item 1: C1's external-check classes (Stage 1b; verified when both lenses reproduced them, comparator status SURVIVES)
C1_PAIRS <- c(PAM_T2D = "PAM_T2D", RAD23B_T2D = "RAD23B_T2D", ARG1_T2D = "ARG1_T2D", PCSK9_CAD = "PCSK9_CAD", TGFB1_CAD = "TGFB1_CAD", FGFR4_T2D_EAS = "FGFR4_T2D")
C1C <- list()
for (k in names(C1_PAIRS)) {
  cp <- C1_PAIRS[[k]]
  cls <- getq(C1R, "C1", paste0("C1.", cp, ".class"), numeric = FALSE)
  v <- chk_c1_class_verified(C1S, cp, cls)
  add_assert(paste0("S9.check.C1_class_verified.", k), v$ok, v$detail)
  C1C[[k]] <- list(class = cls, verified = v$ok, status = C1S$status[C1S$key_id == paste0("C1H.", cp, ".class")],
                   p_resource = getq(C1R, "C1", paste0("C1.", cp, ".p_resource")), p_wald1 = getq(C1R, "C1", paste0("C1.", cp, ".p_wald1")),
                   p_status = C1S$status[C1S$key_id == paste0("C1H.", cp, ".p_resource")],
                   alt_flag = C1S$alt_flag[C1S$key_id == paste0("C1H.", cp, ".p_resource")],
                   bonf6 = getq(C1R, "C1", paste0("C1.", cp, ".bonf6"), numeric = FALSE),
                   label = getq(C1R, "C1", paste0("C1.", cp, ".label"), numeric = FALSE))
}

# W9F item 2b: the RAD23B section 1 p-source/rounding clause, re-evaluated on the C1 adjudicator's numbers
adj_rad <- C1ADJ[C1ADJ$pair == "RAD23B_T2D"]; stopifnot(nrow(adj_rad) == 1L)
w1_extrow <- B$W1[B$W1$rule_ref == "s1 Wording rule: external check"]; stopifnot(nrow(w1_extrow) == 1L)
RCL <- chk_rad23b_clause(adj_rad, w1_extrow$branch)
add_assert("S9.check.RAD23B_clause_fires_on_adjudicator_numbers", RCL$ok, RCL$detail)
add_assert("S9.check.RAD23B_bonf6_agrees", identical(C1C$RAD23B_T2D$bonf6, "FALSE") && identical(adj_rad$bonf6_R, "False") && identical(adj_rad$bonf6_W, "False") &&
             grepl("0.05/6: not met", w1_extrow$branch, fixed = TRUE),
           sprintf("0.05/6 (0.0083): C1 bonf6 %s; adjudicator bonf6 on resource p %s / on Wald p %s; W1 branch '%s'", C1C$RAD23B_T2D$bonf6, adj_rad$bonf6_R, adj_rad$bonf6_W, w1_extrow$branch))
RAD23B_CLAUSE <- sprintf("the classification depends on the p source and on rounding of the printed DIAMANTE estimates (pre-declaration section 1; fires on the C1 adjudicator's numbers: first-order rounded-Wald p %s, rounding interval %s straddles 0.05; class at the interval ends %s); the check does not meet the six-check threshold 0.05/6 (0.0083)",
                         format(as.numeric(adj_rad$p_wald1), digits = 3), adj_rad$p_wald1_rounding_interval, gsub("_", " ", adj_rad$class_W_at_interval_ends))

# UMOD control criterion (uses the unadjusted p, Rule 4)
um <- list(hA5 = getq(R$W6b, "W6b", "S6.umod.H4.A.503"), hB5 = getq(R$W6b, "W6b", "S6.umod.H4.B.503"),
           hA9 = getq(R$W6b, "W6b", "S6.umod.H4.A.489"), hB9 = getq(R$W6b, "W6b", "S6.umod.H4.B.489"),
           p = getq(R$W6b, "W6b", "S6.umod.Q.aware_2nd.p"))
UMOD_OK <- all(c(um$hA5, um$hB5, um$hA9, um$hB9) >= 0.80) && um$p >= 0.05
w6b_umod <- B$W6b[B$W6b$rule_ref == "6 UMOD control"]
stopifnot(nrow(w6b_umod) == 1L)
add_assert("S9.check.umod_control_matches_W6b", UMOD_OK == identical(w6b_umod$branch, "behaved as a control"),
           sprintf("cond PP.H4 A/B 503 %.3f/%.3f 489 %.3f/%.3f; unadjusted Q p %s -> %s; W6b '%s'", um$hA5, um$hB5, um$hA9, um$hB9, f2(um$p),
                   if (UMOD_OK) "behaved as a control" else "did not behave", w6b_umod$branch))

# renal: D6 and the 'sensitive to {list}' qualifiers (closes W6a's PENDING_W6b rows)
rn <- list(or_lo = getq(R$W6a, "W6a", "S6.primary.or_lo"), lo_lo = getq(R$W6a, "W6a", "S6.b.locus_outcome.or_lo"),
           lo_hi = getq(R$W6a, "W6a", "S6.b.locus_outcome.or_hi"), c_lo = getq(R$W6a, "W6a", "S6.c.pooled_lo"),
           c_hi = getq(R$W6a, "W6a", "S6.c.pooled_hi"), a_lo = getq(R$W6a, "W6a", "S6.b.a_pair.or_lo"),
           a_hi = getq(R$W6a, "W6a", "S6.b.a_pair.or_hi"), d_lo = getq(R$W6b, "W6b", "S6.d.lo"), d_hi = getq(R$W6b, "W6b", "S6.d.hi", allow_inf = TRUE),
           or = getq(R$W6a, "W6a", "S6.primary.or_mle"))
RDL <- renal_d6_list(rn$or_lo, rn$d_lo, rn$d_hi, rn$lo_lo, rn$lo_hi, rn$c_lo, rn$c_hi, rn$a_lo, rn$a_hi)   # w9_common.R (one implementation)
D6 <- RDL$d6; RENAL_LIST <- RDL$list
w6b_d6 <- B$W6b[B$W6b$rule_ref == "11 D6"]; w6b_ren <- B$W6b[B$W6b$rule_ref == "6 renal OR wording; 11 D6"]
stopifnot(nrow(w6b_d6) == 1L, nrow(w6b_ren) == 1L)
add_assert("S9.check.D6_matches_W6b", identical(D6, w6b_d6$filled_template),
           sprintf("primary exact lower %.3f, (d) exact lower %.3f -> '%s'; W6b '%s'", rn$or_lo, rn$d_lo, D6, w6b_d6$filled_template))
renal_phrase <- if (length(RENAL_LIST) == 0) "" else paste0("sensitive to ", paste(RENAL_LIST[-length(RENAL_LIST)], collapse = ", "),
                                                           if (length(RENAL_LIST) > 1) " and " else "", RENAL_LIST[length(RENAL_LIST)])
add_assert("S9.check.renal_list_matches_W6b", if (nzchar(renal_phrase)) grepl(renal_phrase, w6b_ren$filled_template, fixed = TRUE) else
                                                !grepl("sensitive to", w6b_ren$filled_template, fixed = TRUE),
           sprintf("locus-outcome %.3f-%.3f, (c) %.3f-%.3f, (a) %.3f-%.3f, (d) %.3f-%.3f -> {%s}", rn$lo_lo, rn$lo_hi, rn$c_lo, rn$c_hi,
                   rn$a_lo, rn$a_hi, rn$d_lo, rn$d_hi, paste(RENAL_LIST, collapse = "; ")))

# ---------------------------------------------------------------- 4. corroboration table (Principle 0)
logline("--- corroboration (Principle 0) ---")
# W9F item 3 travel text (defined here, before the corroboration table, which carries it on its ARG1 rows: W9F step 2b, D3)
quote_ok(GOV$prespec_red$path, "this framing was chosen after the pre-stated test returned WEAKENS", "redacted prespec")
quote_ok(GOV$prespec_red$path, "the pre-stated consequence of WEAKENS was a kidney-led framing", "redacted prespec")
quote_ok(GOV$prespec_red$path, "The ARG1 verdict **WEAKENS is frozen at its R1 value**", "redacted prespec")
ARG1_FRAMING <- "ARG1 framing disclosure (Principle 0; stated in the Methods or Results wherever ARG1 is introduced as the principal worked example): this framing was chosen after the pre-stated test returned WEAKENS; the pre-stated consequence of WEAKENS was a kidney-led framing; WEAKENS is frozen at its R1 value (Rule 1)"
oc_pairs <- c("PAM-T2D" = "PAM_T2D", "RAD23B-T2D" = "RAD23B_T2D", "ARG1-T2D" = "ARG1_T2D", "PCSK9-CAD" = "PCSK9_CAD",
              "TGFB1-CAD" = "TGFB1_CAD", "FGFR4-T2D (EAS)" = "FGFR4_T2D_EAS", "ARG1-T2D signal B" = "ARG1_T2D_signalB")
stopifnot(setequal(R1OC$pair, names(oc_pairs)))
w1_ext <- B$W1[B$W1$rule_ref == "s1 Wording rule: external check"]
stopifnot(nrow(w1_ext) == 1L)
COR <- list()
for (k in seq_len(nrow(R1OC))) {
  r <- R1OC[k]; key <- oc_pairs[[r$pair]]
  wald <- as.numeric(r$wald); disc <- as.numeric(r$discovery_wald); pw <- as.numeric(r$p); pf <- as.numeric(r$p_out_file)
  stopifnot(is.finite(wald), is.finite(disc), is.finite(pw), is.finite(pf))
  same <- sign(wald) == sign(disc)
  cls_re <- if (same && pw < 0.05) "concordant" else if (same) "directionally concordant, not significant" else if (pw < 0.05) "discordant" else "inconclusive"
  add_assert(paste0("S9.check.R1_class_rederived.", key), identical(cls_re, r$class),
             sprintf("[READ_FROM_SEEN] wald %+.4f vs discovery %+.4f, Wald p %s -> %s; R1 class %s", wald, disc, f2(pw), cls_re, r$class))
  st1 <- NA_character_; st1_note <- "no Stage 1 re-derivation (not among W1's five European checks)"
  if (key %in% EXT$pair) {
    e <- EXT[EXT$pair == key]
    st1 <- e$class_resource_p
    st1_note <- sprintf("W1 ext_checks.tsv (Stage 1, from raw): resource p %s, first-order Wald p %s, second-order %s, same sign %s; class on resource p %s / 1st %s / 2nd %s",
                        e$p_resource, f2(as.numeric(e$wald_p1)), f2(as.numeric(e$wald_p2)), e$same_sign, e$class_resource_p, e$class_1st, e$class_2nd)
    add_assert(paste0("S9.check.W1_ext_class_agrees_R1.", key), identical(e$class_resource_p, r$class) && identical(e$class_1st, r$class),
               st1_note)
  }
  # W9F item 1: condition 1 comes from C1's verified class (Stage 1b, from raw, classified on the resource's own p;
  # comparator status SURVIVES). The R1 class (READ_FROM_SEEN) is printed beside it. ARG1 signal B is not a C1 pair.
  s3 <- if (key %in% names(S3)) S3[[key]]$branch else NA_character_
  s2 <- if (key == "TGFB1_CAD") S2B else NA_character_
  # W9F step 2b (D2): every scenario in BOTH S3 adjudicator files, every MAF mode, the Rule 1 and the PP.H4-only branch
  s3_adj <- if (key %in% names(S3)) paste(sort(unique(c(S3ADJ$branch[S3ADJ$pair == key], S3ADJ$branch_alt_h4only[S3ADJ$pair == key]))), collapse = "|") else "not run for this pair"
  travel <- character(0)
  if (key == "ARG1_T2D_signalB") {
    declared_class <- r$class; c1 <- NULL
    conc_basis <- "not a C1 pair (outcome-check pre-specification section 6 row); R1 declared record (READ_FROM_SEEN)"
    cv <- list(value = NA, c1 = identical(declared_class, "concordant"), c2 = TRUE, c3 = TRUE)
    verdict <- "NOT_A_PAIR"
    why <- "row of the outcome-check pre-specification section 6 (ARG1 signal B, the triangle test), not a pre-declared pair; its check is 'inconclusive', so condition 1 would fail in any case"
    label <- "UKB-free external consortium check, not independent replication (R1 record; r0 not carried for this row)"
  } else {
    c1 <- C1C[[key]]
    # an unverified class is never used silently: its assert row is ASSERT_FAILED and the reason says so (asserts never stop the run)
    declared_class <- c1$class
    conc_basis <- sprintf("C1 (S1R1_R2_11_extchecks, Stage 1b, from raw): class on the resource's own p %s, first-order Wald p %s beside; comparator C1H.%s.class %s",
                          format(c1$p_resource, digits = 4), format(c1$p_wald1, digits = 4), C1_PAIRS[[key]], c1$status)
    add_assert(paste0("S9.check.C1_class_agrees_R1_and_W1.", key), identical(c1$class, r$class) && (!(key %in% EXT$pair) || identical(EXT$class_resource_p[EXT$pair == key], c1$class)),
               sprintf("C1 %s; R1 [READ_FROM_SEEN] %s; W1 ext_checks %s", c1$class, r$class, if (key %in% EXT$pair) EXT$class_resource_p[EXT$pair == key] else "n/a"))
    cv <- corroborated(identical(declared_class, "concordant"), s3, s2, key == "TGFB1_CAD")
    verdict <- if (cv$value) "CORROBORATED" else "NOT_CORROBORATED"
    fails <- c(if (!cv$c1) "external check not concordant", if (!cv$c2) sprintf("section 3 branch %s (needs 1, or 2 without the add-on)", s3),
               if (!cv$c3) sprintf("section 2 branch %s (needs (b))", s2))
    why <- if (cv$value) "all Principle 0 conditions hold" else paste(fails, collapse = "; ")
    if (!isTRUE(c1$verified)) why <- paste0(why, "; CONDITION 1 RESTS ON A C1 CLASS THAT IS NOT VERIFIED (see S9.check.C1_class_verified.", key, ")")
    # outcome-check label from C1; its r0 is READ_FROM_SEEN (S1R1_overlap_ldsc.tsv), read by C1, not computed
    label <- if (key == "FGFR4_T2D_EAS") c1$label else paste0(c1$label, " [r0 READ_FROM_SEEN by C1 from S1R1_overlap_ldsc.tsv]; not independent replication")
    label <- paste0(label, "; ", R0_LABEL)   # Stage 2 verifier round: C1.<pair>.label is a not-key
    if (key == "RAD23B_T2D") travel <- c(travel, RAD23B_CLAUSE)
  }
  # W9F item 3 on the ARG1 corroboration rows of corroboration.tsv (W9F step 2b, D3)
  if (key %in% c("ARG1_T2D", "ARG1_T2D_signalB")) travel <- c(travel, ARG1_FRAMING)
  COR[[k]] <- data.table(pair = r$pair, pkey = key, resource = r$resource, instrument = r$instrument,
    check_class_C1 = if (is.null(c1)) "not a C1 pair" else c1$class, C1_class_status = if (is.null(c1)) "" else c1$status,
    C1_p_resource = if (is.null(c1)) NA_real_ else c1$p_resource, C1_p_wald1 = if (is.null(c1)) NA_real_ else c1$p_wald1,
    C1_p_resource_status = if (is.null(c1)) "" else paste0(c1$p_status, if (nzchar(c1$alt_flag)) paste0(" (", c1$alt_flag, ")") else ""),
    check_class_R1 = r$class, check_wald = wald, check_p_wald1 = pw, check_p_file = pf, discovery_wald = disc, same_sign = same,
    check_class_used = declared_class, concordance_basis = conc_basis, stage1_recheck = st1_note,
    s3_branch = ifelse(is.na(s3), "not run for this pair", s3), s3_branch_all_adjudicator_scenarios = s3_adj, s2_branch = ifelse(is.na(s2), "not applicable", s2),
    cond1_concordant = cv$c1, cond2_s3 = if (is.na(s3)) COND2_NA else as.character(cv$c2), cond3_s2 = cv$c3, verdict = verdict, reason = why,
    outcome_check_label = label, travel_w9f = paste(travel, collapse = " | "),
    source_flag = if (is.null(c1)) "READ_FROM_SEEN (R1_outcome_checks_tsv, md5 7d7f1431...)" else "condition 1: C1 results.tsv (md5 b43f1d94...), verified (C1_SURVIVAL); R1 class beside: READ_FROM_SEEN")
  logline(sprintf("corroboration %-18s %-17s %s", r$pair, verdict, why))
}
COR <- rbindlist(COR)
mut_late("C_scen_empty", mut_cell(COR, COR$pkey == "PCSK9_CAD", "s3_branch_all_adjudicator_scenarios", ""))
# ---- Stage 2 opening fix 3: condition 2 is 'not applicable' where section 3 was not run; the verdicts must not change ----
mut_late("R7_cor", mut_cell(COR, COR$pkey == "PAM_T2D", "outcome_check_label", sub(R0_LABEL, "", COR$outcome_check_label[COR$pkey == "PAM_T2D"], fixed = TRUE)))
add_assert("S9.check.rule7_corroboration_label",
           all(grepl(R0_LABEL, COR$outcome_check_label[COR$check_class_C1 != "not a C1 pair"], fixed = TRUE)) && sum(COR$check_class_C1 != "not a C1 pair") == 6L,
           sprintf("the outcome_check_label column of corroboration.tsv (C1.<pair>.label, not a key; r0 READ_FROM_SEEN) carries '%s' for all %d C1 pairs", R0_LABEL, sum(COR$check_class_C1 != "not a C1 pair")))
mut_late("X_cond2na", mut_cell(COR, COR$pkey == "FGFR4_T2D_EAS", "cond2_s3", "TRUE"))
add_assert("S9.check.cond2_not_applicable_where_s3_not_run",
           all((COR$cond2_s3 == COND2_NA) == (COR$s3_branch == "not run for this pair")) &&
             setequal(COR$pair[COR$cond2_s3 == COND2_NA], c("ARG1-T2D", "FGFR4-T2D (EAS)", "ARG1-T2D signal B")) &&
             identical(setNames(COR$verdict, COR$pair)[names(COR_VERDICTS_BEFORE)], COR_VERDICTS_BEFORE),
           sprintf("condition 2 prints '%s' for %s (section 3 not run); verdicts identical to the pre-fix corroboration.tsv: %s",
                   COND2_NA, paste(COR$pair[COR$cond2_s3 == COND2_NA], collapse = ", "), paste(sprintf("%s %s", names(COR_VERDICTS_BEFORE), COR_VERDICTS_BEFORE), collapse = "; ")))
CORR_LIST <- COR$pair[COR$verdict == "CORROBORATED"]
# condition 2 under every S3 adjudicator scenario (coloc prior cap; DIAMANTE EAF printed 0/1): it must not depend on them.
# W9F step 2b (D2): the scope is every row of BOTH adjudicator files (S3ADJ), every MAF mode, Rule 1 and PP.H4-only branches.
for (k in which(COR$pkey %in% names(S3))) {
  sc <- strsplit(COR$s3_branch_all_adjudicator_scenarios[k], "|", fixed = TRUE)[[1]]
  c2_all <- vapply(sc, function(b) corroborated(TRUE, s3_label(b))$c2, TRUE)
  nr <- table(S3ADJ$file[S3ADJ$pair == COR$pkey[k]])
  # W9F step 2c: never vacuous - at least one scenario branch, and rows from BOTH adjudicator files, for this pair
  nonempty <- length(sc) > 0L && all(sc != "") && all(c("adj_S3_coloc_pkg.tsv", "adj_S3_scenarios.tsv") %in% names(nr)) && all(nr > 0L)
  add_assert(paste0("S9.check.cond2_invariant_to_S3_conventions.", COR$pkey[k]), nonempty && all(c2_all == COR$cond2_s3[k]),
             sprintf("section 3 branch across every S3 adjudicator scenario (adj_S3_coloc_pkg.tsv %d rows, adj_S3_scenarios.tsv %d rows; every MAF mode; Rule 1 and PP.H4-only branches): {%s} -> condition 2 %s; W3 primary input %s -> %s",
                     as.integer(nr["adj_S3_coloc_pkg.tsv"]), as.integer(nr["adj_S3_scenarios.tsv"]),
                     paste(sc, collapse = ", "), paste(unique(c2_all), collapse = "/"), COR$s3_branch[k], COR$cond2_s3[k]))
}
if (!nzchar(MUT)) fwrite(COR, file.path(OUT_W, "corroboration.tsv"), sep = "\t", quote = TRUE)
mut_late("C_cond1", mut_cell(COR, COR$pkey == "RAD23B_T2D", "check_class_used", "R1-only"))
mut_late("C_rad_clause", mut_cell(COR, COR$pkey == "RAD23B_T2D", "travel_w9f", ""))
mut_late("C_arg1_framing", mut_cell(COR, COR$pkey == "ARG1_T2D", "travel_w9f", ""))
for (g in names(COR_GATES)) { rg <- COR_GATES[[g]](COR); add_assert(paste0("S9.check.w9f.", g), rg$ok, rg$detail) }
for (k in seq_len(nrow(COR))) {
  cr <- COR[k]
  pre <- paste0("S9.corrob.", cr$pkey, ".")
  addr(paste0(pre, "check_class"), "Principle 0 'corroborated' (condition 1); outcome-check pre-specification s5", paste0(cr$pair, " external check class used"),
       cr$resource, value = cr$check_class_used, p = if (is.na(cr$C1_p_resource)) cr$check_p_wald1 else cr$C1_p_resource,
       notes = paste0("basis: ", cr$concordance_basis, "; R1 class beside (READ_FROM_SEEN) ", cr$check_class_R1, ", Wald ", format(cr$check_wald, digits = 6),
                      ", Wald p ", format(cr$check_p_wald1, digits = 6), ", file p ", format(cr$check_p_file, digits = 6), "; ", cr$stage1_recheck,
                      if (nzchar(cr$travel_w9f)) paste0("; travels: ", cr$travel_w9f) else ""))
  addr(paste0(pre, "s3_branch"), "Principle 0 condition 2; section 3", cr$pair, "W3", value = cr$s3_branch,
       notes = paste0("W3 primary input (p + MAF + N), re-derived by W9 from PP.H3/PP.H4; branches under the S3 adjudicator's scenarios: ", cr$s3_branch_all_adjudicator_scenarios))
  addr(paste0(pre, "s2_branch"), "Principle 0 condition 3; section 2", cr$pair, "W2", value = cr$s2_branch, notes = "TGFB1-CAD only")
  addr(paste0(pre, "corroborated"), "Principle 0 'corroborated'", cr$pair, value = cr$verdict,
       notes = paste0(cr$reason, "; ", if (cr$pkey == "ARG1_T2D_signalB") "condition 1 from the R1 record (READ_FROM_SEEN; never counts as reproduced)" else "condition 1 from C1's verified class"))
}
addr("S9.corrob.R1_opening_list", "Principle 0: the R1 opening names every pair that meets the rule, and no other", "pairs meeting the corroboration rule",
     value = paste(CORR_LIST, collapse = "; "), n = length(CORR_LIST), notes = "condition 1 from C1's verified classes (Stage 1b; comparator SURVIVES for all six); R1 classes beside (READ_FROM_SEEN)")

# ---------------------------------------------------------------- 5. DPEP1 s merge (W6b's DPEP1-eGFR s into the DPEP1 s table)
logline("--- DPEP1 s merge ---")
MS <- list()
for (tr in c("protein", "cysgly")) for (st in c("susie", "p1fr", "p1np", "p2np")) for (pn in c("503", "489")) {
  q <- sprintf("S5.s.%s.%s.%s", tr, st, pn)
  MS[[length(MS) + 1L]] <- data.table(pair = "DPEP1-cysteinylglycine", trait = if (tr == "protein") "DPEP1 protein" else "cysteinylglycine",
    set = st, panel = pn, s = getq(R$W5, "W5", q), source = "W5", quantity_id = q,
    fit = "W5 DPEP1 x cysteinylglycine harmonised set")
}
for (st in c("susie.503", "p1fr.503", "p1np.503", "p2np.489")) for (tr in c("s", "s_protein")) {
  q <- sprintf("S6.dpep1egfr.%s.%s", st, tr)
  sp <- strsplit(st, ".", fixed = TRUE)[[1]]
  MS[[length(MS) + 1L]] <- data.table(pair = "DPEP1-eGFR", trait = if (tr == "s") "eGFR (CKDGen 2019 EA)" else "DPEP1 protein",
    set = sp[1], panel = sp[2], s = getq(R$W6b, "W6b", q), source = "W6b", quantity_id = q,
    fit = "W6b DPEP1 x CKDGen eGFR harmonised set (DPEP1-eGFR l3-finemap-style run)")
}
# Stage 2 step 2 (user decision 2026-10-09, seven fits at every locus): W6b now also fits p2np/503, susie/489 and p1fr/489 for DPEP1-eGFR
# (s only; S6.d.DPEP1_eGFR.s.<trait>.<set>.<panel>), so the W6b DPEP1-eGFR s table has the seven distinct fits (p1np/489 = p2np/489).
for (st in c("p2np.503", "susie.489", "p1fr.489")) for (tr in c("outcome", "protein")) {
  q <- sprintf("S6.d.DPEP1_eGFR.s.%s.%s", tr, st)
  sp <- strsplit(st, ".", fixed = TRUE)[[1]]
  MS[[length(MS) + 1L]] <- data.table(pair = "DPEP1-eGFR", trait = if (tr == "outcome") "eGFR (CKDGen 2019 EA)" else "DPEP1 protein",
    set = sp[1], panel = sp[2], s = getq(R$W6b, "W6b", q), source = "W6b", quantity_id = q,
    fit = "W6b DPEP1 x CKDGen eGFR harmonised set (additional fit for the seven-fit s_max, Stage 2 step 2)")
}
MS <- rbindlist(MS)
stopifnot(nrow(MS) == 30L, all(is.finite(MS$s)), all(MS$s > 0))
if (!nzchar(MUT)) fwrite(MS, file.path(OUT_W, "dpep1_s_merged.tsv"), sep = "\t", quote = TRUE)
smax <- function(pr, tr) max(MS$s[MS$pair == pr & MS$trait == tr])
DP <- list(cys_prot = smax("DPEP1-cysteinylglycine", "DPEP1 protein"), cys_met = smax("DPEP1-cysteinylglycine", "cysteinylglycine"),
           eg_prot = smax("DPEP1-eGFR", "DPEP1 protein"), eg_out = smax("DPEP1-eGFR", "eGFR (CKDGen 2019 EA)"))
DP$cys_adm <- DP$cys_prot < 0.10 && DP$cys_met < 0.10
DP$eg_adm  <- DP$eg_prot < 0.10 && DP$eg_out < 0.10
w5_s_prot <- getq(R$W5, "W5", "S5.s.protein.max"); w5_s_met <- getq(R$W5, "W5", "S5.s.cysgly.max")
w6b_smax <- getq(R$W6b, "W6b", "S6.d.DPEP1_eGFR.s_max"); w6b_sout <- getq(R$W6b, "W6b", "S6.d.DPEP1_eGFR.s_outcome")
add_assert("S9.check.dpep1_cysgly_smax_matches_W5", abs(DP$cys_prot - w5_s_prot) < 1e-9 && abs(DP$cys_met - w5_s_met) < 1e-9,
           sprintf("merged max protein %.4f / cys-gly %.4f vs W5 S5.s.*.max %.4f / %.4f", DP$cys_prot, DP$cys_met, w5_s_prot, w5_s_met))
add_assert("S9.check.dpep1_egfr_smax_matches_W6b", abs(DP$eg_out - w6b_sout) < 1e-9 && abs(max(DP$eg_prot, DP$eg_out) - w6b_smax) < 1e-9,
           sprintf("merged max eGFR %.4f / protein %.4f vs W6b S6.d.DPEP1_eGFR s_outcome %.4f, s_max %.4f", DP$eg_out, DP$eg_prot, w6b_sout, w6b_smax))
add_assert("S9.check.dpep1_w5_admissible_flag", identical(getq(R$W5, "W5", "S5.s.pair.admissible", numeric = FALSE), if (DP$cys_adm) "TRUE" else "FALSE"),
           "W5 S5.s.pair.admissible vs the merged-table rule")
for (k in seq_len(nrow(MS))) addr(sprintf("S9.merge.dpep1.s.%s.%s.%s.%s", ifelse(MS$pair[k] == "DPEP1-eGFR", "egfr_pair", "cysgly_pair"),
  ifelse(grepl("protein", MS$trait[k]), "protein", ifelse(MS$pair[k] == "DPEP1-eGFR", "eGFR", "cysgly")), MS$set[k], MS$panel[k]),
  "Rule 3(i); W9 merge (W6b DPEP1-eGFR s into the DPEP1 s table)", paste0(MS$trait[k], " s (susieR estimate_s_rss) on ", MS$fit[k]),
  MS$source[k], MS$panel[k], value = MS$s[k], notes = paste0(MS$source[k], " ", MS$quantity_id[k], "; locally computed s is a proxy for the pipeline fit"))
addr("S9.merge.dpep1.smax.cysgly_pair.protein", "Rule 3(i)", "max s, DPEP1 protein, over W5's 8 fits", "W5", "both", value = DP$cys_prot)
addr("S9.merge.dpep1.smax.cysgly_pair.cysgly", "Rule 3(i)", "max s, cysteinylglycine, over W5's 8 fits", "W5", "both", value = DP$cys_met)
addr("S9.merge.dpep1.smax.egfr_pair.protein", "Rule 3(i)", "max s, DPEP1 protein, over W6b's 7 distinct DPEP1-eGFR fits", "W6b", "both", value = DP$eg_prot,
     notes = sprintf("W5's protein fits (DPEP1 x cys-gly set) give max %.4f beside; either way >= 0.10", DP$cys_prot))
addr("S9.merge.dpep1.smax.egfr_pair.eGFR", "Rule 3(i)", "max s, eGFR, over W6b's 7 distinct DPEP1-eGFR fits", "W6b", "both", value = DP$eg_out)
addr("S9.merge.dpep1.admissible.cysgly_pair", "Rule 3(i)", "SuSiE with reference LD admissible for DPEP1-cysteinylglycine", value = if (DP$cys_adm) "TRUE" else "FALSE")
addr("S9.merge.dpep1.admissible.egfr_pair", "Rule 3(i)", "SuSiE with reference LD admissible for DPEP1-eGFR", value = if (DP$eg_adm) "TRUE" else "FALSE",
     notes = sprintf("W6b now fits all seven distinct fits (Stage 2 step 2; earlier four: susie/p1fr/p1np on 503, p2np on 489), as W5's 4 x 2 cross does; every DPEP1-eGFR eGFR s is >= %.3f", min(MS$s[MS$pair == "DPEP1-eGFR" & MS$trait == "eGFR (CKDGen 2019 EA)"])))

# ---------------------------------------------------------------- 5b. section 6(d): fit-set readings (W9F item 2a; DECIDED 2026-10-09)
# USER DECISION 2026-10-09 (Stage 2): the s_max that enters Rule 3(i) / section 6(d) uses ALL SEVEN distinct SuSiE fits at every locus
# (one rule everywhere); the earlier four-fit value is a labelled SENSITIVITY. W2, W4b and W6b were re-run with the seven-fit set
# (Stage 2 step 2), so W6b's own S6.d.* output is now the seven-fit reading.
# Reading "seven" (PRIMARY): every set on both panels, seven distinct fits (p1np@489 = p2np@489). Both Stage 1b lenses used seven.
# Reading "four" (SENSITIVITY): one panel per variant set, the earlier four fits (susie, p1fr, p1np on 503; p2np on 489).
# Every number printed here is still read from the S6 adjudicator's files (A_third) and cross-checked against W6b's own output
# (S6.d.<cand>.s_max, S6.d.<cand>.s_max_4fit, S6.d.*). W9 re-derives the counts from the per-candidate s values and the D6 /
# qualifier levels from each reading's CI.
logline("--- section 6(d) fit-set readings (adjudicator files; seven fits primary, four fits sensitivity) ---")
num <- function(x) suppressWarnings(as.numeric(x))
A6F$sP_n <- num(A6F$sP); A6F$sO_n <- num(A6F$sO)
smax_set <- function(four) {
  s <- vapply(c(S6D_KIDNEY, S6D_NONKID), function(cc) {   # (not 'cid': data.table would resolve it to the column)
    x <- A6F[A6F$cid == cc & (!four | A6F$in_impl4 == "True")]
    stopifnot(nrow(x) == (if (four) 4L else 7L))
    if (cc == "ING1_CAD") { stopifnot(all(is.na(x$sP_n))); return(NA_real_) }   # no UKB-PPP protein window: not computable
    max(c(x$sP_n, x$sO_n)) }, 0)
  s }
SMX <- list(four = smax_set(TRUE), seven = smax_set(FALSE))
n2_drop <- function(four) { x <- A6N[A6N$reading == "drop_shift_rows" & (!four | A6N$in4 == "True")]; max(num(x$sP), num(x$sO)) }
n2_keep <- function(four) { x <- A6N[A6N$reading == "keep_shift_rows" & (!four | A6N$in4 == "True")]; max(num(x$sP), num(x$sO)) }
# the fit-set file and the NOTCH2 rebuild agree on NOTCH2 with the shift rows kept (adjudicator check K6/K8 territory)
add_assert("S9.check.s6d_notch2_two_files_agree", abs(n2_keep(TRUE) - SMX$four[["NOTCH2_T2D"]]) < 1e-6 && abs(n2_keep(FALSE) - SMX$seven[["NOTCH2_T2D"]]) < 1e-6,
           sprintf("NOTCH2-T2D s_max, shift rows kept: four fits %.6f (fitset) / %.6f (notch2 rebuild); seven fits %.6f / %.6f",
                   SMX$four[["NOTCH2_T2D"]], n2_keep(TRUE), SMX$seven[["NOTCH2_T2D"]], n2_keep(FALSE)))
# Stage 2 step 2 identity: W6b's own seven-fit s_max (primary) reproduces the adjudicator's seven-fit maxima, and W6b's s_max_4fit (sensitivity)
# reproduces the adjudicator's four-fit maxima, for every computable candidate (tolerance 1e-6; the observed difference is logged).
w6b_smax <- vapply(setdiff(c(S6D_KIDNEY, S6D_NONKID), "ING1_CAD"), function(cid) getq(R$W6b, "W6b", sprintf("S6.d.%s.s_max", cid)), 0)
w6b_smax4 <- vapply(setdiff(c(S6D_KIDNEY, S6D_NONKID), "ING1_CAD"), function(cid) getq(R$W6b, "W6b", sprintf("S6.d.%s.s_max_4fit", cid)), 0)
add_assert("S9.check.s6d_seven_fit_primary_is_W6b", max(abs(w6b_smax - SMX$seven[names(w6b_smax)])) < 1e-6,
           sprintf("max |seven-fit s_max (adjudicator) - W6b S6.d.<cand>.s_max| = %.2g over %d computable candidates", max(abs(w6b_smax - SMX$seven[names(w6b_smax)])), length(w6b_smax)))
add_assert("S9.check.s6d_four_fit_sensitivity_is_W6b", max(abs(w6b_smax4 - SMX$four[names(w6b_smax4)])) < 1e-6,
           sprintf("max |four-fit s_max (adjudicator) - W6b S6.d.<cand>.s_max_4fit| = %.2g over %d computable candidates", max(abs(w6b_smax4 - SMX$four[names(w6b_smax4)])), length(w6b_smax4)))
D6R <- list()
for (fs in c("four", "seven")) for (shift in c("keep", "drop")) for (ver in c("reclassify_ING1", "keep_ING1")) {
  s <- SMX[[fs]]
  if (shift == "drop") s[["NOTCH2_T2D"]] <- if (fs == "four") n2_drop(TRUE) else n2_drop(FALSE)
  tb <- s6d_table(s, keep_not_computable = (ver == "keep_ING1"))
  a <- A6D[A6D$fit_set == (if (fs == "seven") "eight" else "four") & A6D$shift_rows == shift & A6D$version == ver]
  stopifnot(nrow(a) == 1L)
  lo <- num(a$lo); hi <- num(a$hi)
  rd <- renal_d6_list(rn$or_lo, lo, hi, rn$lo_lo, rn$lo_hi, rn$c_lo, rn$c_hi, rn$a_lo, rn$a_hi)
  D6R[[length(D6R) + 1L]] <- data.table(fit_set = fs, fits = if (fs == "four") "one panel per variant set: susie, p1fr, p1np on 503; p2np on 489 (the earlier implementer set; sensitivity)" else
                                          "every set on both panels: seven distinct fits (p1np@489 = p2np@489; both Stage 1b lenses; the decided primary reading, user decision 2026-10-09)",
    build_shift_rows = shift, ING1_version = if (ver == "reclassify_ING1") "(d) estimate: not-computable ING1-CAD reclassified" else "'kept' version printed beside: ING1-CAD kept",
    n_reclassified = as.integer(a$n_reclassified), a = as.integer(a$a), b = as.integer(a$b), c = as.integer(a$c), d = as.integer(a$d),
    or = num(a$or), lo = lo, hi = hi, p = num(a$p), NOTCH2_s_max = num(a$NOTCH2_s_max), NOTCH2_reclassified = num(a$NOTCH2_s_max) >= 0.10,
    d6_base = rd$d6, qualifier_list = paste(rd$list, collapse = ", "), reclassified = a$reclassified,
    w9_recount_agrees = tb$n == as.integer(a$n_reclassified) && tb$a == as.integer(a$a) && tb$b == as.integer(a$b) && tb$c == as.integer(a$c) &&
                        tb$d == as.integer(a$d) && identical(paste(tb$reclassified, collapse = ";"), a$reclassified) && abs(num(a$NOTCH2_s_max) - s[["NOTCH2_T2D"]]) < 1e-6,
    adjudicator_levels_agree = identical(gsub("_", " ", a$d6_base), rd$d6) && identical(a$LD_guard_named == "TRUE", "the LD guard" %in% rd$list),
    source = "S1R1_R2_1b_verify/S6/adjudication/A_third/a3_d_table.tsv (md5 f438c8af); W9 recount from a3_allcand_fitset.tsv / a3_notch2_fits.tsv")
}
D6R <- rbindlist(D6R)
if (!nzchar(MUT)) fwrite(D6R, file.path(OUT_W, "s6d_fitset_readings.tsv"), sep = "\t", quote = TRUE)
add_assert("S9.check.s6d_counts_rederived_from_s", all(D6R$w9_recount_agrees),
           sprintf("W9 recount of n, a, b, c, d, the reclassified list and NOTCH2 s_max from the per-candidate s values agrees with a3_d_table.tsv in %d of %d rows", sum(D6R$w9_recount_agrees), nrow(D6R)))
add_assert("S9.check.s6d_branch_levels_invariant", length(unique(paste(D6R$d6_base, D6R$qualifier_list))) == 1L && all(D6R$adjudicator_levels_agree) &&
             identical(unique(D6R$d6_base), D6) && identical(unique(D6R$qualifier_list), paste(RENAL_LIST, collapse = ", ")),
           sprintf("over %d combinations (fit set x build-shift rows x ING1 version): D6 {%s}, qualifier list {%s}; W6b: '%s' / {%s}; adjudicator columns agree in %d of %d",
                   nrow(D6R), paste(unique(D6R$d6_base), collapse = "|"), paste(unique(D6R$qualifier_list), collapse = "|"), D6, paste(RENAL_LIST, collapse = ", "),
                   sum(D6R$adjudicator_levels_agree), nrow(D6R)))
impl_row <- D6R[D6R$fit_set == "seven" & D6R$build_shift_rows == "keep" & grepl("^\\(d\\)", D6R$ING1_version)]
kept_row <- D6R[D6R$fit_set == "seven" & D6R$build_shift_rows == "keep" & grepl("^'kept'", D6R$ING1_version)]
hi_same <- function(a, b) identical(is.finite(a), is.finite(b)) && (!is.finite(a) || abs(a - b) < 1e-3)   # +Inf (a zero cell) equals +Inf
add_assert("S9.check.s6d_seven_keep_is_W6b_output",
           impl_row$n_reclassified == as.integer(getq(R$W6b, "W6b", "S6.d.reclassified", "n")) &&
           all(c(impl_row$a, impl_row$b, impl_row$c, impl_row$d) == sapply(c("a", "b", "c", "d"), function(x) getq(R$W6b, "W6b", paste0("S6.d.", x)))) &&
           abs(impl_row$lo - rn$d_lo) < 1e-6 && hi_same(impl_row$hi, rn$d_hi) &&
           all(c(kept_row$a, kept_row$b, kept_row$c, kept_row$d) == sapply(c("a", "b", "c", "d"), function(x) getq(R$W6b, "W6b", paste0("S6.d.kept.", x)))),
           sprintf("seven fits, shift rows kept: n %d, table %d/%d/%d/%d, lo %.6f, hi %s (W6b S6.d: n %s, lo %.6f, hi %s); 'kept' table %d/%d/%d/%d",
                   impl_row$n_reclassified, impl_row$a, impl_row$b, impl_row$c, impl_row$d, impl_row$lo, format(impl_row$hi),
                   getq(R$W6b, "W6b", "S6.d.reclassified", "n", numeric = FALSE), rn$d_lo, format(rn$d_hi), kept_row$a, kept_row$b, kept_row$c, kept_row$d))
S6D_SUMMARY <- function(fs, shift = "keep") {
  r <- D6R[D6R$fit_set == fs & D6R$build_shift_rows == shift & grepl("^\\(d\\)", D6R$ING1_version)]
  k <- D6R[D6R$fit_set == fs & D6R$build_shift_rows == shift & grepl("^'kept'", D6R$ING1_version)]
  sprintf("%s fits: n_reclassified %d; (d) table %d/%d vs %d/%d; OR %s (exact 95%% CI %.3f-%s), p %s; 'kept' version OR %s (lower bound %.3f)",
          fs, r$n_reclassified, r$a, r$a + r$b, r$c, r$c + r$d, if (is.finite(r$or)) sprintf("%.2f", r$or) else "Inf", r$lo,
          if (is.finite(r$hi)) sprintf("%.1f", r$hi) else "Inf", f2(r$p), if (is.finite(k$or)) sprintf("%.2f", k$or) else "Inf", k$lo)
}
N2_FLAG <- sprintf("NOTCH2-T2D borderline (Rule 3(i) cut 0.10): s_max with the build-shift rows kept: seven fits (primary) %.4f, four fits (sensitivity) %.4f; dropped: seven fits %.4f, four fits %.4f; reclassified under the seven-fit rule with the rows kept (user decision 2026-10-09: rows kept, disclosed), not under the four-fit sensitivity",
                   SMX$seven[["NOTCH2_T2D"]], SMX$four[["NOTCH2_T2D"]], n2_drop(FALSE), n2_drop(TRUE))
add_assert("S9.check.s6d_notch2_only_seven_keep_crosses", identical(D6R$NOTCH2_reclassified, D6R$fit_set == "seven" & D6R$build_shift_rows == "keep"),
           N2_FLAG)
S6D_TRAVEL <- sprintf("section 6(d) fit set (USER DECISION 2026-10-09: seven fits primary, one rule at every locus; the four-fit value is a sensitivity): %s | %s | %s | no section 6 branch level changes under any combination of fit set, build-shift rows and ING1 version (D6 '%s'; sensitive to %s)",
                      S6D_SUMMARY("seven"), S6D_SUMMARY("four"), N2_FLAG, D6, paste(RENAL_LIST, collapse = ", "))
logline(S6D_TRAVEL)
for (k in seq_len(nrow(D6R))) {
  id <- sprintf("S9.s6d.%s.%s_shift.%s", D6R$fit_set[k], D6R$build_shift_rows[k], if (grepl("^\\(d\\)", D6R$ING1_version[k])) "d_estimate" else "kept_version")
  addr(paste0(id, ".or"), "6(d); Rule 3(i); W9F item 2a (user decision 2026-10-09: seven fits primary)", paste0("(d) exact conditional-MLE OR, ", D6R$fits[k], "; build-shift rows ", D6R$build_shift_rows[k], "; ", D6R$ING1_version[k]),
       "S6 adjudicator (A_third)", value = D6R$or[k], ci_lo = D6R$lo[k], ci_hi = D6R$hi[k], p = D6R$p[k], n = D6R$n_reclassified[k],
       notes = sprintf("table %d/%d kidney vs %d/%d (a/b/c/d %d/%d/%d/%d); n = n_reclassified; NOTCH2 s_max %.4f; D6 '%s'; list {%s}; reclassified %s",
                       D6R$a[k], D6R$a[k] + D6R$b[k], D6R$c[k], D6R$c[k] + D6R$d[k], D6R$a[k], D6R$b[k], D6R$c[k], D6R$d[k], D6R$NOTCH2_s_max[k],
                       D6R$d6_base[k], D6R$qualifier_list[k], D6R$reclassified[k]))
}
for (fs in names(SMX)) for (cid in names(SMX[[fs]]))
  addr(sprintf("S9.s6d.smax.%s.%s", fs, cid), "6(d); Rule 3(i); W9F item 2a", sprintf("s_max, %s fits (shift rows kept)", fs), "S6 adjudicator (A_third)",
       value = if (is.na(SMX[[fs]][[cid]])) "not computable (no UKB-PPP protein window)" else SMX[[fs]][[cid]],
       notes = if (is.na(SMX[[fs]][[cid]])) "reclassified in the (d) estimate; kept in the 'kept' version" else
         sprintf("%s (>= 0.10 reclassifies)", if (SMX[[fs]][[cid]] >= 0.10) "reclassified" else "kept"))

# ---------------------------------------------------------------- 6. section 12 items (post hoc, supplement-only; never in a branch)
logline("--- section 12 items (post hoc, supplement-only) ---")
cells <- list(FinnGen.marg = "S4.grid.FinnGen.marg", FinnGen.cond = "S4.grid.FinnGen.cond", DIAMANTE.marg = "S4.grid.DIAMANTE.marg", DIAMANTE.cond = "S4.grid.DIAMANTE.cond")
G <- rbindlist(lapply(names(cells), function(cn) {
  q <- cells[[cn]]
  data.table(cell = cn, d = getq(R$W4a, "W4a", paste0(q, ".diff")), se = getq(R$W4a, "W4a", paste0(q, ".diff_se")),
             p_aware = getq(R$W4a, "W4a", paste0(q, ".p_aware_2nd")), p_ign = getq(R$W4a, "W4a", paste0(q, ".p_ign")))
}))
for (q in paste0(unlist(cells), ".diff_se")) conform("W4a", q, "", "")
# (a) inverse-variance fixed-effect pool of the two resources' A-B differences, per estimand
for (est in c("marg", "cond")) {
  g <- G[G$cell %in% paste0(c("FinnGen.", "DIAMANTE."), est)]
  stopifnot(nrow(g) == 2L)
  w <- 1 / g$se^2; e <- sum(w * g$d) / sum(w); s <- sqrt(1 / sum(w)); z <- e / s
  pp <- 2 * pnorm(-abs(z)); Q <- sum(w * (g$d - e)^2); Qp <- pchisq(Q, df = 1, lower.tail = FALSE)
  lab <- if (est == "marg") "marginal" else "conditional"
  nt <- sprintf("POST HOC, SUPPLEMENT-ONLY (section 12(a)); never in a branch; inputs W4a %s: FinnGen %+.5f (SE %.5f), DIAMANTE %+.5f (SE %.5f), aware 2nd-order SEs; fixed effect assumes independent resources, but DIAMANTE includes Finnish cohorts that FinnGen also holds (Amendment 1 A2), so the pooled SE is anticonservative to an unmeasured degree",
                lab, g$d[1], g$se[1], g$d[2], g$se[2])
  ref <- "section 12(a) pooled two-resource difference"
  esn <- paste0("T2D log-odds per SD ARG1, A-B difference, ", lab, "; IVW fixed-effect pool of FinnGen R10 and DIAMANTE noUKB")
  addr(paste0("S9.s12a.", est, ".est"), ref, esn, "FinnGen R10 + DIAMANTE noUKB", "503", "aware", "2nd", e, s, e - qnorm(0.975) * s, e + qnorm(0.975) * s, pp, 2L, "log-odds per SD", nt)
  addr(paste0("S9.s12a.", est, ".Q"), ref, paste0("Cochran Q between resources, ", lab), "FinnGen R10 + DIAMANTE noUKB", "503", "aware", "2nd", Q, p = Qp, n = 2L, notes = "df = 1; POST HOC, SUPPLEMENT-ONLY")
  logline(sprintf("s12(a) %-11s pooled %+.4f (SE %.4f; 95%% CI %+.4f to %+.4f) p %s; Cochran Q %.4f p %s", lab, e, s,
                  e - qnorm(0.975) * s, e + qnorm(0.975) * s, f2(pp), Q, f2(Qp)))
}
# (b) pre-reclassification Holm values
for (cv in c("aware", "ign")) {
  pv <- if (cv == "aware") G$p_aware else G$p_ign
  h4 <- holm(pv)
  for (k in seq_len(nrow(G))) addr(sprintf("S9.s12b.grid4.%s.%s", cv, G$cell[k]), "section 12(b) pre-reclassification Holm (four grid cells)",
    paste0("ARG1 T2D A vs B, ", G$cell[k], "; Holm over the 4 grid cells"), "W4a", "503", cv, "2nd", h4[k], p = h4[k], n = 4L,
    notes = sprintf("unadjusted %s p %s; POST HOC, SUPPLEMENT-ONLY; the grid is a stability display, not a family (Rule 4)", cv, format(pv[k], digits = 7)))
  sel <- G$cell %in% c("FinnGen.cond", "DIAMANTE.marg")
  h2 <- holm(pv[sel])
  for (k in which(sel)) addr(sprintf("S9.s12b.orig2.%s.%s", cv, G$cell[k]), "section 12(b) pre-reclassification Holm (two originally reported cells)",
    paste0("ARG1 T2D A vs B, ", G$cell[k], "; Holm over FinnGen conditional and DIAMANTE marginal"), "W4a", "503", cv, "2nd",
    h2[match(k, which(sel))], p = h2[match(k, which(sel))], n = 2L,
    notes = sprintf("unadjusted %s p %s; POST HOC, SUPPLEMENT-ONLY", cv, format(pv[k], digits = 7)))
  logline(sprintf("s12(b) %-5s grid4 Holm %s | orig2 Holm %s", cv, paste(sprintf("%s=%s", G$cell, format(h4, digits = 4)), collapse = " "),
                  paste(sprintf("%s=%s", G$cell[sel], format(h2, digits = 4)), collapse = " ")))
}


# ---- Stage 2 opening fix 5: Rule 7 amendment on the section 12 pooled values (computed by W9 from W4a's not-key SEs) ----
s12i <- which(startsWith(vapply(RES_ROWS, function(r) r$quantity_id, ""), "S9.s12"))
for (j in s12i) data.table::set(RES_ROWS[[j]], 1L, "notes", paste0(RES_ROWS[[j]]$notes, " | ", RULE7_LABEL,
  " (inputs W4a S4.grid.*.diff_se and p values; Rule 7 amendment, user decision 12: supplement-only position)"))
mut_late("R7_s12", { data.table::set(RES_ROWS[[s12i[1]]], 1L, "notes", sub(RULE7_LABEL, "", RES_ROWS[[s12i[1]]]$notes, fixed = TRUE)); MUT_APPLIED <- TRUE })
s12n <- vapply(RES_ROWS[s12i], function(r) r$notes, "")
add_assert("S9.check.rule7_s12_labelled", length(s12i) == 16L && all(grepl(RULE7_LABEL, s12n, fixed = TRUE)) && all(grepl("SUPPLEMENT-ONLY", s12n, fixed = TRUE)),
           sprintf("%d section 12 result rows (S9.s12a/b.*) all carry '%s' and stay supplement-only", length(s12i), RULE7_LABEL))

# ---- Stage 2 opening fix 5 (verifier round): the DPEP1 merge rows (analyst s values) and the S9.s6d rows (adjudicator third computation) ----
UALL <- unique(c(unlist(lapply(R, function(x) x$quantity_id)), C1R$quantity_id))
# Stage 2 step 2: the seven-fit re-run of W2, W4b and W6b ADDED quantity ids that no Stage 1b key list knows (the three extra s fits per trait and the
# four-fit sensitivity values, all analyst-computed s values). They are classified here as not-keys; the key maps (_restricted/) are not edited.
# The existing key ids S6.d.<cand>.s_max (and the not-key S2/S4 s_max ids) now hold the seven-fit values; the W6b seven-fit s_max reproduces the
# S6 adjudicator's third computation (S9.check.s6d_seven_fit_primary_is_W6b), but the stored Stage 1b comparator statuses are of the four-fit values.
STEP2_NK_RE <- paste0("(_4fit([.]|$)|[.]s[.](protein|outcome)[.](p2np[.]503|susie[.]489|p1fr[.]489)$|^S2[.]susie[.][a-z_0-9]+[.](p2np_503panel|susie_489panel|p1fr_489panel)[.]s$|",
                      "^S4[.]s[.][A-Za-z0-9_]+[.](p2np[.]503|susie[.]489|p1fr[.]489)$)")
step2_nk <- function(uall) { x <- grep(STEP2_NK_RE, uall, value = TRUE); x[!(x %in% KEYQ) & !(x %in% NKQ)] }
STEP2_NEW <- step2_nk(UALL)
NKQ <- c(NKQ, STEP2_NEW)
add_assert("S9.check.step2_new_ids_are_seven_fit_additions", length(STEP2_NEW) > 0L && !any(STEP2_NEW %in% KEYQ) &&
             all(grepl("^(S2|S4|S6)[.]", STEP2_NEW)),
           sprintf("%d quantity ids added by the Stage 2 step 2 seven-fit re-run (W2 %d, W4b %d, W6b %d) are classified as not-keys (analyst-computed s values); the Stage 1b key maps are untouched",
                   length(STEP2_NEW), sum(startsWith(STEP2_NEW, "S2.")), sum(startsWith(STEP2_NEW, "S4.")), sum(startsWith(STEP2_NEW, "S6."))))
# the six C1 class ids (C1.<pair>.class) are the results rows of the categorical C1 class keys: class levels, not numbers; they count as keys
KEYQ <- unique(c(KEYQ, grep("^C1[.][A-Za-z0-9_]+[.]class$", C1R$quantity_id, value = TRUE)))
rid <- vapply(RES_ROWS, function(r) r$quantity_id, "")
mg <- which(startsWith(rid, "S9.merge.dpep1.s.") | startsWith(rid, "S9.merge.dpep1.smax."))
mg_ids <- lapply(mg, function(j) { q <- rid[j]
  if (startsWith(q, "S9.merge.dpep1.smax.cysgly_pair")) MS$quantity_id[MS$pair == "DPEP1-cysteinylglycine"]
  else if (startsWith(q, "S9.merge.dpep1.smax.egfr_pair")) MS$quantity_id[MS$pair == "DPEP1-eGFR"]
  else r7_ids(RES_ROWS[[j]]$notes, UALL)$ids[1] })
mg_nk <- lapply(mg_ids, function(x) x[(x %in% NKQ) & !(x %in% KEYQ)])
for (a in seq_along(mg)) if (length(mg_nk[[a]]))
  data.table::set(RES_ROWS[[mg[a]]], 1L, "notes", paste0(RES_ROWS[[mg[a]]]$notes, " | ", RULE7_LABEL, sprintf(" (%d of %d quantity ids; Rule 7 amendment, user decision 12: supplement or side positions only)", length(mg_nk[[a]]), length(mg_ids[[a]]))))
s6r <- which(startsWith(rid, "S9.s6d."))
for (j in s6r) data.table::set(RES_ROWS[[j]], 1L, "notes", paste0(RES_ROWS[[j]]$notes, if (nzchar(RES_ROWS[[j]]$notes)) " | " else "", ADJ_LABEL, " (S6 adjudicator A_third; Rule 7 amendment, user decision 12)"))
mut_late("R7_ms", { data.table::set(RES_ROWS[[mg[1]]], 1L, "notes", sub(RULE7_LABEL, "", RES_ROWS[[mg[1]]]$notes, fixed = TRUE)); MUT_APPLIED <- TRUE })
mut_late("RES_adj6d", { data.table::set(RES_ROWS[[s6r[1]]], 1L, "notes", sub(ADJ_LABEL, "", RES_ROWS[[s6r[1]]]$notes, fixed = TRUE)); MUT_APPLIED <- TRUE })
mgn <- vapply(RES_ROWS[mg], function(r) r$notes, "")
add_assert("S9.check.rule7_merge_rows_labelled",
           length(mg) == nrow(MS) + 4L && length(unlist(lapply(mg_ids, function(x) x[!(x %in% KEYQ) & !(x %in% NKQ)]))) == 0L &&
             all(grepl(RULE7_LABEL, mgn, fixed = TRUE) == (lengths(mg_nk) > 0L)) && any(lengths(mg_nk) > 0L),
           sprintf("%d S9.merge.dpep1 s and s_max rows (DPEP1 s values from W5 and W6b results): %d quote not-key quantity ids and exactly those carry '%s'; none is unclassified", length(mg), sum(lengths(mg_nk) > 0L), RULE7_LABEL))

# ---------------------------------------------------------------- 7. branches_integrated.tsv
logline("--- branches_integrated ---")
# Rule 4 text: the LITERAL reading is the default (W9F item 6); the applied reading I1 is printed beside it
holm_txt <- function(m) {
  f <- fam[fam$member == m]
  switch(f$nominal_code,
    NOMINAL = sprintf("nominal [W9: Holm p %s >= 0.05 over the 8-member Rule 4 family; unadjusted p %s; the applied reading I1 also attaches 'nominal']", f2(f$p_holm), f2(f$p_unadj)),
    "NOMINAL (non-rejection)" = sprintf("nominal [W9, literal reading of Rule 4 (default): Holm p %s >= 0.05 over the 8-member family; this is a non-rejection worded on the unadjusted p %s (Rule 4); applied reading I1, printed beside: no 'nominal' on a non-rejection; the user decides]",
                                        f2(f$p_holm), f2(f$p_unadj)),
    SURVIVES_HOLM = sprintf("no 'nominal' [W9: Holm p %s < 0.05 over the 8-member Rule 4 family; unadjusted p %s; the same under the applied reading I1]", f2(f$p_holm), f2(f$p_unadj)))
}
member_for <- function(ws, rr, br) {
  if (ws == "W2"  && grepl("^prespec s2 wording rule", rr)) return("tgfb1_cad")
  if (ws == "W4a" && rr == "4 Arginine A vs B") return("arg1_arginine")
  if (ws == "W4a" && rr == "4 T2D") return("arg1_t2d")
  if (ws == "W4a" && rr == "4 Combined (test (iii))") return("arg1_test3")
  if (ws == "W4a" && rr == "4 Ornithine") return("arg1_ornithine")
  if (ws == "W4a" && rr == "Rule 4 family entry") return(fe_map[[br]])
  if (ws == "W5"  && grepl("wording rule: metabolite Q$", rr)) return("dpep1_cysgly")
  if (ws == "W5"  && grepl("wording rule: eGFR Q$", rr)) return("dpep1_egfr")
  if (ws == "W6b" && rr == "6 UMOD control") return("umod_egfr")
  NA_character_
}
corr_txt <- function(key) {
  cr <- COR[COR$pkey == key]
  sprintf("corroboration (W9, Principle 0): %s %s (%s; condition 1 from C1's verified class '%s', R1 class '%s' beside)%s", cr$pair,
          if (cr$verdict == "CORROBORATED") "IS externally corroborated" else "is NOT called externally corroborated", cr$reason,
          cr$check_class_used, cr$check_class_R1, if (nzchar(cr$travel_w9f)) paste0("; ", cr$travel_w9f) else "")
}

# ---- W9F travel texts (each quoted fragment is asserted present in its pinned source document) ----
# (ARG1_FRAMING and its three quote checks are in section 4, before the corroboration table)
quote_ok(GOV$prespec_red$path, "The five-signal set was chosen after seeing the five-signal Q values.", "redacted prespec")
quote_ok(GOV$prespec_red$path, "those signals were derived with the same reference LD that the GLS uses", "redacted prespec")
OMNI6 <- list(arg = getq(R$W4b, "W4b", "S4.omni.arginine.6sig.p_aware_2nd"), fg = getq(R$W4b, "W4b", "S4.omni.T2D_FG.6sigproxy.p_aware_2nd"),
              fg5 = getq(R$W4b, "W4b", "S4.omni.T2D_FG.6sigdrop5.p_aware_2nd"), dm = getq(R$W4b, "W4b", "S4.omni.T2D_DM.6sig.p_aware_2nd"))
w4b_six <- B$W4b[B$W4b$rule_ref == "4 Omnibus (six-signal sensitivity)"]; stopifnot(nrow(w4b_six) == 1L)
printed_p <- as.numeric(sub("^p ", "", regmatches(w4b_six$filled_template, gregexpr("p [0-9.]+(e-?[0-9]+)?", w4b_six$filled_template))[[1]]))
add_assert("S9.check.omnibus_six_signal_p_match_W4b_row", length(printed_p) == 4L &&
             all(vapply(unlist(OMNI6), function(v) any(abs(printed_p - signif(v, 2)) <= 1e-9 * abs(v)), TRUE)),
           sprintf("W4b results six-signal p: arginine %s, T2D FinnGen (proxy) %s, FinnGen (five available) %s, DIAMANTE %s; printed in W4b's six-signal row: %s",
                   f2(OMNI6$arg), f2(OMNI6$fg), f2(OMNI6$fg5), f2(OMNI6$dm), paste(printed_p, collapse = ", ")))
OMNI_DISC <- "post hoc set choice (section 4(iv)): the five-signal set was chosen after seeing the five-signal Q values; the stated reason: those signals were derived with the same reference LD that the GLS uses"
OMNI_TRAVEL <- list(
  "4 Omnibus (arginine)" = sprintf("%s; six-signal sensitivity beside (UKB-PPP in-sample set, 503 LD; analyst-computed, not a Stage 1b headline key): arginine p %s", OMNI_DISC, f2(OMNI6$arg)),
  "4 Omnibus (T2D)" = sprintf("%s; six-signal sensitivity beside (UKB-PPP in-sample set, 503 LD; analyst-computed, not a Stage 1b headline key): T2D FinnGen R10 p %s (rs2698107 proxied: DEVIATION, absent from FinnGen R10; the five available signals p %s), T2D DIAMANTE p %s",
                              OMNI_DISC, f2(OMNI6$fg), f2(OMNI6$fg5), f2(OMNI6$dm)))
# RAD23B section 3 PP.H3: the drop version was seen before substitution was chosen (W3 README deviation 1)
quote_ok(W9F_IN$W3_readme$path, "Those posteriors were therefore seen before substitution was chosen as primary.", "W3 README")
RDROP <- list(h3_5 = getq(R$W3, "W3", "S3.RAD23B_T2D.pmafdrop.1e-05.H3"), h3_6 = getq(R$W3, "W3", "S3.RAD23B_T2D.pmafdrop.1e-06.H3"))
adj_drop <- S3SC[S3SC$pair == "RAD23B_T2D" & S3SC$maf_mode == "drop" & S3SC$scenario == "declared"]; stopifnot(nrow(adj_drop) == 1L)
add_assert("S9.check.RAD23B_drop_version_matches_S3_adjudicator", abs(RDROP$h3_5 - num(adj_drop$H3_5)) < 1e-6 && abs(RDROP$h3_6 - num(adj_drop$H3_6)) < 1e-6,
           sprintf("drop-version PP.H3: W3 %s (1e-5) / %s (1e-6); S3 adjudicator %s / %s", format(RDROP$h3_5), format(RDROP$h3_6), adj_drop$H3_5, adj_drop$H3_6))
RAD23B_S3_DISC <- sprintf("DEVIATION disclosure (W3 README deviation 1): the drop version's posteriors (SNPs whose DIAMANTE EAF prints as 0 or 1 dropped: PP.H3 %.2f at 1e-5, %.2f at 1e-6; branch 4) were seen in the first development run before MAF substitution was chosen as primary; the choice rests on keeping the declared window; both versions are printed and their branches agree",
                          RDROP$h3_5, RDROP$h3_6)
EAS_DERIVED <- "derived route; not new evidence"
quote_ok(W9F_IN$posthoc_spec$path, "**Status: POST HOC.**", "post hoc spec")   # the section 8 declared caveat (W9F item 7)

IR <- list(); n_ph <- 0L; ph_members <- character(0)
for (w in names(B)) {
  b <- B[[w]]
  for (k in seq_len(nrow(b))) {
    rr <- b$rule_ref[k]; q0 <- b$qualifiers[k]; qr <- q0; note <- character(0); travel <- character(0)
    tpl_int <- b$filled_template[k]
    m <- NA_character_
    if (grepl("PENDING_HOLM", paste(q0, b$filled_template[k], b$branch[k]), fixed = TRUE)) {
      m <- member_for(w, rr, b$branch[k])
      if (is.na(m)) stop(sprintf("UNMAPPED nominal placeholder in %s row %d (%s)", w, k, rr))
      if (!grepl("nominal: PENDING_HOLM", q0, fixed = TRUE)) stop(sprintf("placeholder outside the qualifiers column in %s row %d", w, k))
      qr <- gsub("nominal: PENDING_HOLM", holm_txt(m), qr, fixed = TRUE)
      n_ph <- n_ph + 1L; ph_members <- c(ph_members, m)
      note <- c(note, sprintf("Rule 4 placeholder resolved: member %s", m))
    }
    if (grepl("corroboration: left to W9", qr, fixed = TRUE)) {
      key <- sub(" primary input.*$", "", b$inputs[k])
      if (!key %in% COR$pkey) stop(sprintf("UNMAPPED corroboration placeholder in %s row %d", w, k))
      qr <- gsub("corroboration: left to W9", corr_txt(key), qr, fixed = TRUE)
      note <- c(note, "corroboration placeholder resolved")
    }
    if (w == "W6a" && grepl("PENDING_W6b", paste(b$branch[k], b$filled_template[k]), fixed = TRUE)) {
      note <- c(note, sprintf("PENDING_W6b resolved from W6b (results md5 %s): (d) exact 95%% CI %.3f-%.3f %s 1 -> 'the LD guard' %s; D6 -> '%s'; renal list {%s} (re-derived by W9, asserts S9.check.D6_matches_W6b / renal_list_matches_W6b)",
        WS$W6b$results, rn$d_lo, rn$d_hi, if (incl1(rn$d_lo, rn$d_hi)) "includes" else "excludes", if (incl1(rn$d_lo, rn$d_hi)) "named" else "not named", D6,
        paste(RENAL_LIST, collapse = ", ")))
    }
    rule_row <- is_rule_row(w, rr, b$branch[k], b$filled_template[k])
    # ---- W9F item 4: Amendment 1 A2 merged into the section 4 T2D claim row
    if (w == "W4a" && rr == "Amendment 1 A2") {
      rule_row <- "fragment merged into the section 4 T2D claim row (Amendment 1 A2; W9F item 4)"
    }
    if (w == "W4a" && rr == "4 T2D" && A2_MET) {
      frag <- sprintf("in both outcome resources (FinnGen p %s; DIAMANTE p %s)", f2(a2$p_fg), f2(a2$p_dm))
      if (!identical(frag, w4a_a2$filled_template)) stop("A2 fragment rebuilt by W9 differs from W4a's A2 row: '", frag, "' vs '", w4a_a2$filled_template, "'")
      anchor <- "; the sign of the correlation between the signals is not identified"
      if (lengths(regmatches(tpl_int, gregexpr(anchor, tpl_int, fixed = TRUE))) != 1L) stop("A2 merge anchor not unique in the T2D template")
      tpl_int <- sub(anchor, paste0(", a difference seen ", frag, anchor), tpl_int, fixed = TRUE)
      n_ld <- lengths(regmatches(qr, gregexpr("sensitive to the reference LD(?! in)", qr, perl = TRUE)))
      if (n_ld != 1L) stop("A2: the generic LD qualifier occurs ", n_ld, " times in the T2D qualifiers (expected 1)")
      qr <- sub("sensitive to the reference LD(?! in)", A2_CAVEAT, qr, perl = TRUE)
      qr <- paste0(qr, "; never 'replicated'/'replication' (DIAMANTE includes Finnish cohorts that FinnGen also holds); the family entry and its Holm adjustment are unchanged (FinnGen primary)")
      note <- c(note, "Amendment 1 A2 merged (W9F item 4): fragment from W4a row 'Amendment 1 A2'; caveat names the resource(s) failing Rule 3(ii) (S9.check.A2_failing_resources_rederived); placement in the sentence is Stage 2's")
    }
    # ---- W9F item 3: the ARG1 framing disclosure travels with every ARG1 slot claim (section 4: W4a, W4b)
    if (w %in% c("W4a", "W4b") && !nzchar(rule_row)) travel <- c(travel, ARG1_FRAMING)
    # ---- W9F item 2: the five-signal omnibus rows carry the set-choice disclosure with the six-signal p beside
    if (w == "W4b" && rr %in% names(OMNI_TRAVEL)) travel <- c(travel, OMNI_TRAVEL[[rr]])
    # ---- W9F item 2c: RAD23B section 3 PP.H3 disclosure
    if (w == "W3" && grepl("^RAD23B_T2D (primary input|pmafdrop)", b$inputs[k])) travel <- c(travel, RAD23B_S3_DISC)
    # ---- W9F item 5: the East Asian eGFR scale claim
    if (w == "W7" && grepl("^EAS_eGFR", b$branch[k])) travel <- c(travel, EAS_DERIVED)
    # ---- W9F item 2a: section 6(d) fit-set readings (user decision 2026-10-09: seven fits primary, four fits sensitivity)
    if (w == "W6b" && rr == "6(d); Rule 3(i)") {
      cid <- gsub("-", "_", sub("^([A-Z0-9]+-[A-Za-z0-9]+):.*$", "\\1", b$filled_template[k]))
      if (!cid %in% names(SMX$four)) stop("unmapped 6(d) candidate row: ", b$filled_template[k])
      fmt_s <- function(x) if (is.na(x)) "not computable" else sprintf("%.4f", x)
      travel <- c(travel, sprintf("section 6(d) fit set USER DECISION 2026-10-09 (seven fits primary): s_max seven fits %s, four fits (sensitivity) %s -> %s under the seven-fit rule%s",
                                  fmt_s(SMX$seven[[cid]]), fmt_s(SMX$four[[cid]]),
                                  if (is.na(SMX$seven[[cid]]) || SMX$seven[[cid]] >= 0.10) "reclassified" else "kept",
                                  if (cid == "NOTCH2_T2D") paste0("; ", N2_FLAG) else ""))
    }
    if (w == "W6b" && rr %in% c("6(d)", "11 D6", "6 renal OR wording; 11 D6")) travel <- c(travel, S6D_TRAVEL)
    if (length(travel)) qr <- paste0(qr, if (nzchar(qr)) " | " else "", "W9F travel: ", paste(travel, collapse = " | "))

    nom_attach <- if (!is.na(m)) fam$nominal_attach[fam$member == m] else FALSE
    nom_attach_I1 <- if (!is.na(m)) fam$nominal_attach_I1[fam$member == m] else FALSE
    front_lab <- function(dc) {
      front <- if (nzchar(rule_row)) NA_character_ else dc$front
      ft <- if (nzchar(rule_row)) paste0("n/a: ", rule_row) else if (is.na(front)) "none declared for this claim" else A1_LABEL[[front]]
      if (!nzchar(rule_row) && identical(front, "LDREF") && w == "W4a" && rr == "4 T2D" && A2_MET)
        ft <- sprintf("%s (Rule 3(ii); Amendment 1 A2 form)", A2_CAVEAT)
      if (!nzchar(rule_row) && identical(front, "COLOC")) {
        cl <- regmatches(paste(b$filled_template[k], qr), regexpr("the [a-zA-Z0-9 -]+ association at signal [^,;|]+ did not colocalize with the protein \\([^)]*\\),? so the ratio there is not anchored to the protein signal", paste(b$filled_template[k], qr)))
        ft <- paste0("colocalization clause (Rule 3(iv)): ", if (length(cl)) cl else "[clause text in the row]")
      }
      if (!nzchar(rule_row) && identical(front, "FIRSTMENTION")) {
        ft <- paste0("first-mention flag: ", paste(dc$fm, collapse = " + "),
          if (length(dc$fm) > 1) sprintf(" (Principle 0 bundles these in one item; if one must be chosen, the written order Rule 2, Rule 3(i), 'proxy' gives '%s': Stage 2 flag I2)", dc$fm[1]) else "")
      }
      ft
    }
    dc  <- detect_caveats(rr, b$branch[k], b$filled_template[k], qr, nom_attach)      # literal reading (default)
    dcI <- detect_caveats(rr, b$branch[k], b$filled_template[k], qr, nom_attach_I1)   # applied reading I1 (beside)
    IR[[length(IR) + 1L]] <- data.table(source = w, source_row = k, rule_ref = rr, inputs = b$inputs[k], branch = b$branch[k],
      filled_template = b$filled_template[k], filled_template_integrated = tpl_int, qualifiers = q0, qualifiers_resolved = qr,
      travel_w9f = paste(travel, collapse = " | "),
      holm_member = ifelse(is.na(m), "", m),
      p_unadj = if (is.na(m)) "" else format(fam$p_unadj[fam$member == m], digits = 6),
      p_holm  = if (is.na(m)) "" else format(fam$p_holm[fam$member == m], digits = 6),
      nominal = if (is.na(m)) "" else fam$nominal_code[fam$member == m],
      nominal_I1 = if (is.na(m)) "" else fam$nominal_code_I1[fam$member == m],
      row_kind = if (nzchar(rule_row)) rule_row else "claim",
      a1_caveats_present = if (nzchar(rule_row) || !length(dc$caveats)) "" else paste(dc$caveats, collapse = " > "),
      a1_front_caveat = front_lab(dc),
      a1_front_caveat_I1 = if (nzchar(rule_row)) "" else front_lab(dcI),
      a1_firstmention_flags = paste(dc$fm, collapse = "; "),
      a1_outcome_check_labels = paste(dc$oclabels, collapse = "; "),
      front_surface_note = if (nzchar(rule_row)) "" else front_surface_note(b$filled_template[k], qr, rr),
      integration_note = paste(note, collapse = " | "))
  }
}
IR <- rbindlist(IR)
# ---------------------------------------------------------------- 7b. Stage 2 opening fix 5: Rule 7 amendment (user decision 12, 2026-10-09)
# Rule 7 (pre-declaration) bars unverified Stage 1 numbers from the manuscript. Under the middle path (declared deviation D-1) only the 134 headline
# keys were re-derived. The amendment admits every other number, labelled RULE7_LABEL, in supplement or side positions only. A number is a KEY if its
# quantity id is a single-id row of a KEY_MAP.tsv; every other quantity id in NOT_KEYS.tsv is a NOT-KEY. The label is attached per row: to qualifiers_resolved
# (never to the A1 front caveat) and to the new columns n_key_numbers, n_notkey_numbers, notkey_quantity_ids, rule7_label.
UALL <- unique(c(unlist(lapply(R, function(x) x$quantity_id)), C1R$quantity_id))
NKQ <- c(NKQ, step2_nk(UALL))   # Stage 2 step 2 additions (see above); idempotent
mut_late("R7_unclass", { NKQ <- setdiff(NKQ, "S5.s.protein.max"); MUT_APPLIED <- TRUE })
# W1 rows 2-8 print numbers inline (no id list in the workstream); their ids were matched by value to W1/results.tsv (p, value, CI columns) and are declared here
W1_INLINE <- list(
  "2" = "S1.rad23b.cut.second_order",
  "3" = paste0("S1.rad23b.ext.", c("resource_p_lead", "wald1_p", "wald1_lo", "wald1_hi", "wald2_p", "wald2_lo", "wald2_hi")),
  "4" = c("S1.rad23b.r2_plink", "S1.ledger.RAD23B_T2D.eaf_exposure_lead", "S1.ledger.RAD23B_T2D.eaf_outcome_lead", "S1.directlead.RAD23B_T2D.wald", "S1.directlead.RAD23B_T2D.p"),
  "5" = c("S1.ledger.ADM_AF.r2_plink", "S1.ledger.ADM_AF.eaf_exposure_lead", "S1.ledger.ADM_AF.eaf_outcome_lead", "S1.directlead.ADM_AF.wald", "S1.directlead.ADM_AF.p"),
  "6" = c("S1.ledger.PECAM1_CAD.r2_plink", "S1.ledger.PECAM1_CAD.eaf_exposure_lead", "S1.ledger.PECAM1_CAD.eaf_outcome_lead", "S1.directlead.PECAM1_CAD.wald", "S1.directlead.PECAM1_CAD.p"),
  "7" = "S1.ledger.ABO_Stroke.r2_plink",
  "8" = c("S1.ledger.ABO_T2D.r2_plink", "S1.directlead.ABO_T2D.wald", "S1.directlead.ABO_T2D.p"))
stopifnot(all(unlist(W1_INLINE) %in% UALL))
R7_NO_NUMBER_ROWS <- data.table(source = "W1", source_row = 1L)   # the evidence-tier row prints no number
r7_row_ids <- function(w, k, rr, inputs) {
  txt <- inputs; map_ok <- TRUE; extra <- character(0)
  if (w %in% names(IDMAP)) {
    m <- IDMAP[[w]][as.integer(IDMAP[[w]]$row) == k]
    if (nrow(m) == 1L) { map_ok <- identical(m$rule_ref, rr); if (grepl("^(C1|S[1-7])\\.", m$inputs)) txt <- m$inputs }
    else map_ok <- FALSE
  }
  if (w == "W3" && k == 11L) { m9 <- IDMAP$W3[as.integer(IDMAP$W3$row) == 9L]; txt <- m9$inputs }   # PCSK9 control: a function of the PCSK9 primary-input row (branches.tsv row 9)
  if (w == "W1" && as.character(k) %in% names(W1_INLINE)) extra <- W1_INLINE[[as.character(k)]]
  x <- r7_ids(txt, UALL)
  list(ids = unique(c(x$ids, extra)), unresolved = x$unresolved, map_ok = map_ok)
}
IR[, c("n_key_numbers", "n_notkey_numbers", "notkey_quantity_ids", "rule7_label", "rule7_adjudicator_label") := ""]
# Verifier round (MAJOR 1/2): numbers added by the W9F travel texts are built from quantity ids that are not in `inputs`. Each travel text is registered with the ids
# it is built from; any row (IR or BR) whose text contains a registered text carries those ids in its classification. adj = TRUE: the numbers come from an adjudicator's
# third computation (S6 A_third), which has its own label ADJ_LABEL.
S6D_FOUR_IDS <- c(paste0("S6.d.", c("a", "b", "c", "d", "p", "or", "lo", "hi", "reclassified")), paste0("S6.d.kept.", c("a", "b", "c", "d", "p", "or", "lo", "hi")))
S6D_FOUR_IDS <- S6D_FOUR_IDS[S6D_FOUR_IDS %in% UALL]
OMNI6_IDS <- c("S4.omni.arginine.6sig.p_aware_2nd", "S4.omni.T2D_FG.6sigproxy.p_aware_2nd", "S4.omni.T2D_FG.6sigdrop5.p_aware_2nd", "S4.omni.T2D_DM.6sig.p_aware_2nd")
TRAVEL_REG <- list(
  s6d_summary   = list(match = S6D_TRAVEL, ids = S6D_FOUR_IDS, adj = TRUE),
  s6d_candidate = list(match = "section 6(d) fit set USER DECISION 2026-10-09 (seven fits primary): s_max seven fits", ids = character(0), adj = TRUE),
  rad23b_clause = list(match = RAD23B_CLAUSE, ids = paste0("S1.rad23b.ext.", c("resource_p_lead", "wald1_p", "wald1_lo", "wald1_hi", "wald2_p")), adj = FALSE),
  omni_six      = list(match = "six-signal sensitivity beside", ids = OMNI6_IDS, adj = FALSE),
  rad23b_s3     = list(match = "DEVIATION disclosure (W3 README deviation 1)", ids = c("S3.RAD23B_T2D.pmafdrop.1e-05.H3", "S3.RAD23B_T2D.pmafdrop.1e-06.H3"), adj = FALSE),
  w6b_resolution = list(match = "PENDING_W6b resolved from W6b", ids = c("S6.d.lo", "S6.d.hi"), adj = FALSE))
# the C1 outcome-check label text (FinnGen overlap r0 etc.) is the not-key C1.<pair>.label; one entry per C1 pair, found only on rows that name the pair
for (key_ in names(C1C)) TRAVEL_REG[[paste0("c1_label_", key_)]] <- list(match = C1C[[key_]]$label, with = paste0("S9.corrob.", key_, "."), ids = paste0("C1.", C1_PAIRS[[key_]], ".label"), adj = FALSE)
NUMFREE_TRAVEL <- list(ARG1_FRAMING = ARG1_FRAMING, EAS_DERIVED = EAS_DERIVED, OMNI_DISC = OMNI_DISC)   # travel texts without a decimal number
reg_hit <- function(txt) names(TRAVEL_REG)[vapply(TRAVEL_REG, function(e) grepl(e$match, txt, fixed = TRUE) && (is.null(e$with) || grepl(e$with, txt, fixed = TRUE)), TRUE)]
ir_txt <- paste(IR$qualifiers_resolved, IR$filled_template_integrated, IR$travel_w9f, IR$integration_note)
ir_hit <- lapply(ir_txt, reg_hit)
R7I <- IR$row_kind == "claim" | startsWith(IR$row_kind, "fragment") | lengths(ir_hit) > 0L
r7_all <- list(); r7_nk <- list(); r7_uncl <- list(); r7_unres <- list(); r7_mapok <- logical(0)
for (i in which(R7I)) {
  x <- r7_row_ids(IR$source[i], IR$source_row[i], IR$rule_ref[i], IR$inputs[i])
  x$ids <- unique(c(x$ids, unlist(lapply(TRAVEL_REG[ir_hit[[i]]], function(e) e$ids))))
  adj_i <- any(vapply(TRAVEL_REG[ir_hit[[i]]], function(e) e$adj, TRUE))
  isk <- x$ids %in% KEYQ; isn <- (x$ids %in% NKQ) & !isk
  r7_all[[as.character(i)]] <- x$ids; r7_nk[[as.character(i)]] <- x$ids[isn]
  r7_uncl[[as.character(i)]] <- x$ids[!isk & !isn]; r7_unres[[as.character(i)]] <- x$unresolved; r7_mapok <- c(r7_mapok, x$map_ok)
  data.table::set(IR, i, "n_key_numbers", as.character(sum(isk))); data.table::set(IR, i, "n_notkey_numbers", as.character(sum(isn)))
  data.table::set(IR, i, "notkey_quantity_ids", paste(x$ids[isn], collapse = "; "))
  if (any(isn)) {
    data.table::set(IR, i, "rule7_label", RULE7_LABEL)
    data.table::set(IR, i, "qualifiers_resolved", paste0(IR$qualifiers_resolved[i], if (nzchar(IR$qualifiers_resolved[i])) " | " else "",
      sprintf("Rule 7 amendment (user decision 12): %d of %d quantity ids on this row are %s (ids in column notkey_quantity_ids); supplement or side positions only, never the front caveat",
              sum(isn), length(x$ids), RULE7_LABEL)))
  }
  if (adj_i) {
    data.table::set(IR, i, "rule7_adjudicator_label", ADJ_LABEL)
    data.table::set(IR, i, "qualifiers_resolved", paste0(IR$qualifiers_resolved[i], if (nzchar(IR$qualifiers_resolved[i])) " | " else "",
      "Rule 7 amendment (user decision 12): the section 6(d) fit-set values on this row come from the S6 adjudicator's third computation: ", ADJ_LABEL))
  }
}
is_nonum <- paste(IR$source, IR$source_row) %in% paste(R7_NO_NUMBER_ROWS$source, R7_NO_NUMBER_ROWS$source_row)
add_assert("S9.check.rule7_every_travelling_number_classified",
           all(lengths(r7_uncl) == 0L) && all(lengths(r7_unres) == 0L) && all(r7_mapok) && all(lengths(r7_all) > 0L | is_nonum[R7I]),
           sprintf("%d integrated claim/fragment rows: %d distinct quantity ids resolved against the 9 workstream results files; %d are keys (KEY_MAP.tsv), %d are not-keys (NOT_KEYS.tsv), %d are neither; %d unresolved id tokens; %d branch-input id maps disagree with the row's rule_ref; %d rows have no id and are declared no-number rows (W1 evidence tier)",
                   sum(R7I), length(unique(unlist(r7_all))), sum(unique(unlist(r7_all)) %in% KEYQ), sum(unique(unlist(r7_all)) %in% NKQ & !(unique(unlist(r7_all)) %in% KEYQ)),
                   length(unlist(r7_uncl)), length(unlist(r7_unres)), sum(!r7_mapok), sum(lengths(r7_all) == 0L)))
i_t3 <- which(IR$source == "W4a" & IR$rule_ref == "4 Combined (test (iii))"); i_six <- which(IR$source == "W4b" & IR$rule_ref == "4 Omnibus (six-signal sensitivity)")
T3_BESIDE <- c("S4.t3.DIAMANTE.cond.503.2nd.p_aware", "S4.t3.FinnGen.cond.489.2nd.p_aware", "S4.t3.FinnGen.cond.503.2nd.p_cojo", "S4.t3.FinnGen.cond.503.1st.p_aware")
mut_late("R7_label", mut_cell(IR, IR$source == "W4a" & IR$rule_ref == "4 Combined (test (iii))", "qualifiers_resolved",
                        sub(RULE7_LABEL, "", IR$qualifiers_resolved[i_t3], fixed = TRUE)))
has_lab <- grepl(RULE7_LABEL, IR$qualifiers_resolved, fixed = TRUE)
nkn <- suppressWarnings(as.integer(IR$n_notkey_numbers))
add_assert("S9.check.rule7_label_on_every_notkey_row",
           all(has_lab[R7I] == (nkn[R7I] > 0L)) && all(IR$rule7_label[R7I] == ifelse(nkn[R7I] > 0L, RULE7_LABEL, "")) && !any(has_lab[!R7I]) &&
             !any(grepl(RULE7_LABEL, IR$a1_front_caveat, fixed = TRUE)) && !any(grepl(RULE7_LABEL, IR$a1_front_caveat_I1, fixed = TRUE)) && !any(grepl(RULE7_LABEL, IR$front_surface_note, fixed = TRUE)) &&
             length(i_t3) == 1L && has_lab[i_t3] && all(T3_BESIDE %in% r7_nk[[as.character(i_t3)]]) &&
             length(i_six) == 1L && has_lab[i_six] && "S4.omni.T2D_DM.6sig.p_aware_2nd" %in% r7_nk[[as.character(i_six)]],
           sprintf("%d of %d integrated claim/fragment rows carry not-key numbers and exactly those %d carry '%s' in qualifiers_resolved and rule7_label (never in a front caveat); named cases: test (iii) beside values (DIAMANTE, 489, beta-space, first-order: %s) and the six-signal omnibus p (S4.omni.*.6sig.*) are labelled",
                   sum(nkn[R7I] > 0L, na.rm = TRUE), sum(R7I), sum(has_lab), RULE7_LABEL, paste(T3_BESIDE, collapse = ", ")))
# travel-text ids (verifier MAJOR 1/2): every registered text is found, its ids are in the row's id column (not-key ones) and classified; number-free texts have no decimal
mut_late("R7_travel", mut_cell(IR, IR$source == "W6b" & IR$source_row == 18L, "notkey_quantity_ids", ""))
reg_ids_all <- unique(unlist(lapply(TRAVEL_REG, function(e) e$ids)))
ir_hit2 <- lapply(paste(IR$qualifiers_resolved, IR$filled_template_integrated, IR$travel_w9f, IR$integration_note), reg_hit)
tr_ok <- vapply(seq_len(nrow(IR)), function(i) { e <- unique(unlist(lapply(TRAVEL_REG[ir_hit2[[i]]], function(z) z$ids))); nk_e <- e[(e %in% NKQ) & !(e %in% KEYQ)]
  all(nk_e %in% strsplit(IR$notkey_quantity_ids[i], "; ", fixed = TRUE)[[1]]) && (length(nk_e) == 0L || nzchar(IR$rule7_label[i])) }, TRUE)
add_assert("S9.check.rule7_travel_text_ids_classified",
           all(tr_ok) && all(reg_ids_all %in% UALL) && all(reg_ids_all %in% c(KEYQ, NKQ)) && all(vapply(setdiff(names(TRAVEL_REG), grep("^c1_label_", names(TRAVEL_REG), value = TRUE)), function(n) any(vapply(ir_hit2, function(h) n %in% h, TRUE)), TRUE)) &&
             !any(vapply(NUMFREE_TRAVEL, function(t) grepl("[0-9]*\\.[0-9]+|[0-9]e-[0-9]", t), TRUE)) && IR$notkey_quantity_ids[IR$source == "W6b" & IR$source_row == 18L] != "",
           sprintf("%d registered W9F travel texts (%s) built from %d quantity ids, all classified; every row containing one carries its not-key ids and the label (%d rows contain a registered text, e.g. W6b row 18, the 6(d) fit-set text); %d travel texts are declared free of decimal numbers (%s)",
                   length(TRAVEL_REG), paste(names(TRAVEL_REG), collapse = ", "), length(reg_ids_all), sum(lengths(ir_hit2) > 0L), length(NUMFREE_TRAVEL), paste(names(NUMFREE_TRAVEL), collapse = ", ")))
mut_late("R7_adj", mut_cell(IR, IR$source == "W6b" & IR$source_row == 18L, "rule7_adjudicator_label", ""))
adj_rows <- vapply(ir_hit2, function(h) any(vapply(TRAVEL_REG[h], function(z) z$adj, TRUE)), TRUE)
adj_res <- vapply(RES_ROWS[s6r], function(r) grepl(ADJ_LABEL, r$notes, fixed = TRUE), TRUE)
add_assert("S9.check.adjudicator_label_on_s6d_values",
           any(adj_rows) && all(IR$rule7_adjudicator_label[adj_rows] == ADJ_LABEL) && all(grepl(ADJ_LABEL, IR$qualifiers_resolved[adj_rows], fixed = TRUE)) &&
             all(IR$rule7_adjudicator_label[!adj_rows] == "") && length(s6r) > 0L && all(adj_res),
           sprintf("%d integrated rows print S6 fit-set values (travel text from the adjudicator's third computation) and carry '%s'; no other integrated row does; all %d S9.s6d.* result rows carry it", sum(adj_rows), ADJ_LABEL, length(s6r)))
add_assert("S9.check.no_unresolved_placeholder", !any(grepl("PENDING_HOLM|left to W9", IR$qualifiers_resolved)),
           sprintf("%d nominal placeholders resolved; none remain in qualifiers_resolved", n_ph))
add_assert("S9.check.every_member_has_a_placeholder_row", setequal(unique(ph_members), fam$member),
           sprintf("members with >= 1 resolved placeholder row: %s", paste(sort(unique(ph_members)), collapse = ", ")))
add_assert("S9.check.W6a_PENDING_W6b_closed", sum(grepl("PENDING_W6b resolved", IR$integration_note)) == sum(grepl("PENDING_W6b", paste(B$W6a$branch, B$W6a$filled_template))),
           "every W6a PENDING_W6b row carries the W6b resolution")
# W9F travel and front-surface gates (functions in w9_common.R; each is mutation-tested on a copy of this output)
irow <- function(src, sel) IR$source == src & sel
mut_late("G_arg1", mut_cell(IR, irow("W4b", IR$rule_ref == "Rule 3(iii) gene naming (ARG1-arginine)"), "travel_w9f", ""))
mut_late("G_omni", mut_cell(IR, irow("W4b", IR$rule_ref == "4 Omnibus (T2D)"), "a1_front_caveat", "none declared for this claim"))
mut_late("G_a2", mut_cell(IR, irow("W4a", IR$rule_ref == "4 T2D"), "qualifiers_resolved",
                          sub("sensitive to the reference LD in FinnGen", "sensitive to the reference LD", IR$qualifiers_resolved[irow("W4a", IR$rule_ref == "4 T2D")], fixed = TRUE)))
mut_late("G_eas", mut_cell(IR, irow("W7", grepl("^EAS_eGFR", IR$branch)), "a1_front_caveat", "none declared for this claim"))
mut_late("G_s8", mut_cell(IR, irow("W7", grepl("^s8 ", IR$rule_ref)), "front_surface_note", ""))
mut_late("G_rad_s3", mut_cell(IR, irow("W3", grepl("^RAD23B_T2D pmafdrop", IR$inputs)), "travel_w9f", ""))
mut_late("G_s6d", mut_cell(IR, irow("W6b", IR$rule_ref == "6(d)"), "travel_w9f",
                           sub("seven fits: n_reclassified", "seven fits: n", IR$travel_w9f[irow("W6b", IR$rule_ref == "6(d)")], fixed = TRUE)))
mut_late("G_nominal", mut_cell(IR, irow("W5", IR$source_row == 1L), "qualifiers_resolved",
                               sub("nominal [W9", "no 'nominal' [W9", IR$qualifiers_resolved[irow("W5", IR$source_row == 1L)], fixed = TRUE)))
# W9F step 2b (D4): a printed six-signal p value changed on the integrated T2D omnibus row (the gate checks values, not labels)
mut_late("N_omni_val", mut_cell(IR, irow("W4b", IR$rule_ref == "4 Omnibus (T2D)"), "qualifiers_resolved",
                                sub(sprintf("T2D DIAMANTE p %s", f2(OMNI6$dm)), "T2D DIAMANTE p 0.90", IR$qualifiers_resolved[irow("W4b", IR$rule_ref == "4 Omnibus (T2D)")], fixed = TRUE)))
# (the 'mutation not applied' guard is at the end of section 8b, after the last hook)
for (g in names(W9F_GATES)) { r <- W9F_GATES[[g]](IR, w4b = R$W4b); add_assert(paste0("S9.check.w9f.", g), r$ok, r$detail) }

# ---------------------------------------------------------------- 8. W9's own branches.tsv
BR <- list()
addb <- function(rule_ref, inputs, branch, filled, quals) BR[[length(BR) + 1L]] <<- data.table(rule_ref = rule_ref, inputs = inputs, branch = branch, filled_template = filled, qualifiers = quals)
ARG1_MEMBERS <- c("arg1_arginine", "arg1_t2d", "arg1_test3", "arg1_ornithine")
for (i in order(fam$rank)) {
  f <- fam[i]
  addb("Rule 4 (Holm, 8-member family)", sprintf("S9.holm.%s.p_unadj %s (%s %s); S9.holm.%s.p_holm %s; rank %d of 8", f$member, format(f$p_unadj, digits = 4), f$ws,
                                                  sub("\\|.*$", "", f$qid), f$member, format(f$p_holm, digits = 4), f$rank),
       sprintf("literal (default): %s | applied I1 (beside): %s", f$nominal_code, f$nominal_code_I1),
       sprintf("%s: %s", f$label, holm_txt(f$member)),
       paste(c(if (f$p_unadj >= 0.05) "Rule 4: a criterion met by not rejecting uses the unadjusted p; Rule 4b wording; the literal reading attaches 'nominal' (default), the applied reading I1 does not: the user decides",
               if (f$member == "arg1_t2d") "Amendment 1 A2: the family entry and its Holm adjustment are unchanged (FinnGen primary); 'nominal' follows the Holm result",
               if (f$member %in% ARG1_MEMBERS) paste0("W9F travel: ", ARG1_FRAMING)), collapse = "; "))
}
for (k in seq_len(nrow(COR))) {
  cr <- COR[k]
  addb("Principle 0 'externally corroborated'",
       sprintf("S9.corrob.%s.check_class %s [%s]; S9.corrob.%s.s3_branch %s; S9.corrob.%s.s2_branch %s", cr$pkey, cr$check_class_used,
               if (cr$pkey == "ARG1_T2D_signalB") "READ_FROM_SEEN" else sprintf("C1 verified; R1 class beside: %s, READ_FROM_SEEN", cr$check_class_R1),
               cr$pkey, cr$s3_branch, cr$pkey, cr$s2_branch),
       cr$verdict,
       if (cr$verdict == "CORROBORATED") sprintf("%s: externally corroborated", cr$pair) else
         if (cr$verdict == "NOT_A_PAIR") sprintf("%s: not a pair-level check; no corroboration statement", cr$pair) else
           sprintf("%s: not called externally corroborated; described by its own %s wording", cr$pair,
                   if (cr$pkey == "TGFB1_CAD") "section 2 and section 3" else if (cr$pkey %in% names(S3)) "section 3" else "check"),
       paste(c(cr$reason, paste0("label: ", cr$outcome_check_label),
               if (cr$pkey == "ARG1_T2D_signalB") "condition 1 rests on the R1 declared record (READ_FROM_SEEN), which never counts as reproduced" else
                 "condition 1 from C1's verified class (Stage 1b, from raw, on the resource's own p); the R1 class is printed beside it",
               if (nzchar(cr$travel_w9f)) paste0("W9F travel: ", cr$travel_w9f)), collapse = "; "))   # ARG1 rows: ARG1_FRAMING is in travel_w9f (D3)
}
addb("Principle 0: the R1 opening names every pair that meets the rule, and no other", "S9.corrob.R1_opening_list",
     sprintf("%d pairs meet the rule", length(CORR_LIST)), paste(CORR_LIST, collapse = "; "),
     "each with its outcome-check label (outcome-check pre-specification s5); the title uses neither 'corroboration' nor 'corroborated'; condition 1 from C1's verified classes")
addb("Rule 3(i) on the merged DPEP1 s table (DPEP1-eGFR)", "S9.merge.dpep1.smax.egfr_pair.protein; S9.merge.dpep1.smax.egfr_pair.eGFR",
     if (DP$eg_adm) "admissible" else "inadmissible (s >= 0.10 for at least one trait)",
     sprintf("SuSiE-based evidence (reference-LD PP.H4, credible sets) for DPEP1-eGFR is reported in the supplement, labelled inadmissible under Rule 3(i), and never cited as support (s_max protein %.3f, eGFR %.3f)", DP$eg_prot, DP$eg_out),
     "locally computed s is a proxy for the pipeline fit; W6b covered susie/p1fr/p1np on 503 and p2np on 489 (W5 covered the full 4 x 2 cross for cys-gly)")
addb("Rule 3(i) on the merged DPEP1 s table (DPEP1-cysteinylglycine; restated from W5)", "S9.merge.dpep1.smax.cysgly_pair.protein; S9.merge.dpep1.smax.cysgly_pair.cysgly",
     if (DP$cys_adm) "admissible" else "inadmissible (s >= 0.10 for at least one trait)",
     sprintf("unchanged from W5 (s_max protein %.3f, cysteinylglycine %.3f)", DP$cys_prot, DP$cys_met), "")
addb("Amendment 1 A1 (front surfaces: one caveat per claim, fixed order)", "branches_integrated.tsv columns a1_caveats_present, a1_front_caveat (literal 'nominal'), a1_front_caveat_I1",
     "applied to every claim row of W1-W7", sprintf("%d claim rows; %d carry at least one A1-ordered caveat; the front caveat differs between the literal and the I1 reading of 'nominal' in %d rows",
     sum(IR$row_kind == "claim"), sum(IR$row_kind == "claim" & nzchar(IR$a1_caveats_present)),
     sum(IR$row_kind == "claim" & IR$a1_front_caveat != IR$a1_front_caveat_I1)),
     "outcome-check pre-specification labels are governed by that document and travel in addition (column a1_outcome_check_labels); every other surface carries the full set (qualifiers_resolved)")
# ---- W9F rows (Stage 1b-C)
addb("W9F item 3: Principle 0 disclosure of the ARG1 framing (travel row)", "branches_integrated.tsv column travel_w9f (W4a and W4b claim rows); corroboration.tsv column travel_w9f (ARG1-T2D, ARG1-T2D signal B); this file's ARG1 Holm and corroboration rows",
     "travels with every ARG1 slot claim", ARG1_FRAMING,
     sprintf("%d W4a/W4b claim rows carry it (S9.check.w9f.arg1_framing) and %d corroboration.tsv rows (S9.check.w9f.arg1_framing_corrob); it is a Methods/Results disclosure, not an A1 front-surface caveat",
             sum(grepl("ARG1 framing disclosure", IR$travel_w9f, fixed = TRUE)), sum(grepl("ARG1 framing disclosure", COR$travel_w9f, fixed = TRUE))))
addb("W9F item 2: section 4(iv) five-signal omnibus set-choice disclosure (travel row)", "S4.omni.*.5sig (W4b); S4.omni.arginine.6sig.p_aware_2nd; S4.omni.T2D_FG.6sigproxy.p_aware_2nd; S4.omni.T2D_FG.6sigdrop5.p_aware_2nd; S4.omni.T2D_DM.6sig.p_aware_2nd",
     "travels with both five-signal omnibus rows", paste(unlist(OMNI_TRAVEL), collapse = " || "),
     "A1: the disclosure is a declared caveat outside Principle 0's named list (A1 step 4), so it is the rows' front caveat; the omnibus is secondary and never in the Abstract")
addb("W9F item 4: Amendment 1 A2 merged into the section 4 T2D claim", "S4.grid.FinnGen.cond.p_aware_2nd; S4.grid.DIAMANTE.marg.p_aware_2nd; S4.ld.boot.FinnGen.cond.aware2nd.share_same; S4.ld.boot.DIAMANTE.marg.aware2nd.share_same; S4.grid.*.489.p_aware_2nd",
     if (A2_MET) "met: merged" else "not met: rev-4.2 template; the other resource's p printed beside",
     IR$filled_template_integrated[IR$source == "W4a" & IR$rule_ref == "4 T2D"],
     sprintf("caveat names the resource(s) failing Rule 3(ii): '%s' (%s); never 'replicated'/'replication'", A2_CAVEAT, A2R$detail))
addb("W9F item 2a: section 6(d) fit-set reading (user decision 2026-10-09: seven fits primary)", "S9.s6d.* (from S1R1_R2_1b_verify/S6/adjudication/A_third/a3_d_table.tsv; s6d_fitset_readings.tsv); W6b S6.d.* (seven-fit output) and S6.d.<cand>.s_max_4fit (sensitivity)",
     "seven fits primary; four fits printed beside as the sensitivity", S6D_TRAVEL,
     "travels with the W6b (d), D6 and renal OR rows and with every 6(d) candidate row (branches_integrated.tsv travel_w9f)")
addb("W9F item 2b: RAD23B section 1 p-source/rounding clause (travel row)", "C1.RAD23B_T2D.class; adj_C1.tsv (p_resource, p_wald1, p_wald1_rounding_interval, p_wald2, bonf6_R, bonf6_W)",
     if (RCL$fires) "fires on the adjudicator's numbers" else "does not fire on the adjudicator's numbers", RAD23B_CLAUSE,
     sprintf("travels with every RAD23B external-check class: W1 external-check row (already carried), the RAD23B corroboration rows (corroboration.tsv, this file) and W3's RAD23B corroboration text; %s", RCL$detail))
addb("W9F item 2c: RAD23B section 3 PP.H3 drop-version disclosure (travel row)", "S3.RAD23B_T2D.pmaf.*; S3.RAD23B_T2D.pmafdrop.1e-05.H3; S3.RAD23B_T2D.pmafdrop.1e-06.H3",
     "travels with the printed PP.H3", RAD23B_S3_DISC, "W3 RAD23B primary-input and pmafdrop rows (branches_integrated.tsv travel_w9f)")
addb("W9F item 5: East Asian eGFR scale claim", "S7.scale.EAS_eGFR_2P2M.median (W7)", sprintf("attached to the East Asian eGFR scale claim row (W7 row %s) and is its A1 front caveat; manuscript surfaces do not exist yet", paste(IR$source_row[IR$source == "W7" & grepl("^EAS_eGFR", IR$branch)], collapse = ", ")),
     EAS_DERIVED, "W7 'EAS_eGFR_2P2M' row: travel_w9f and a1_front_caveat (S9.check.w9f.eas_derived_route)")
addb("W9F item 7: section 8 front-surface note", "W7 's8' row", "front_surface_note filled for its declared caveat ('post hoc')",
     IR$front_surface_note[IR$source == "W7" & grepl("^s8 ", IR$rule_ref)], "P31R_POSTHOC_SPEC.md (md5 573fd781) 'Status: POST HOC'")
BR <- rbindlist(BR)
# Stage 2 opening fix 5: the same Rule 7 amendment on W9's own branches rows that quote workstream numbers (ids parsed from `inputs`; the merged DPEP1 s rows
# quote the 24 merged s values, whose ids are MS$quantity_id)
br_txt <- paste(BR$inputs, BR$branch, BR$filled_template, BR$qualifiers)
br_hit <- lapply(br_txt, reg_hit)
br_ids <- lapply(seq_len(nrow(BR)), function(i) {
  x <- unique(c(r7_ids(BR$inputs[i], UALL)$ids, unlist(lapply(TRAVEL_REG[br_hit[[i]]], function(e) e$ids))))
  if (startsWith(BR$rule_ref[i], "Rule 3(i) on the merged DPEP1 s table (DPEP1-eGFR")) x <- unique(c(x, MS$quantity_id[MS$pair == "DPEP1-eGFR"]))
  if (startsWith(BR$rule_ref[i], "Rule 3(i) on the merged DPEP1 s table (DPEP1-cysteinylglycine")) x <- unique(c(x, MS$quantity_id[MS$pair == "DPEP1-cysteinylglycine"]))
  x })
br_nk <- lapply(br_ids, function(x) x[(x %in% NKQ) & !(x %in% KEYQ)])
br_uncl <- lapply(br_ids, function(x) x[!(x %in% KEYQ) & !(x %in% NKQ)])
for (i in which(lengths(br_nk) > 0L))
  data.table::set(BR, i, "qualifiers", paste0(BR$qualifiers[i], if (nzchar(BR$qualifiers[i])) "; " else "",
    sprintf("Rule 7 amendment (user decision 12): %d of %d quantity ids quoted on this row are %s (not-key ids: %s)", length(br_nk[[i]]), length(br_ids[[i]]), RULE7_LABEL, paste(br_nk[[i]], collapse = ", "))))
br_adj <- vapply(br_hit, function(h) any(vapply(TRAVEL_REG[h], function(e) e$adj, TRUE)), TRUE)
for (i in which(br_adj)) data.table::set(BR, i, "qualifiers", paste0(BR$qualifiers[i], "; Rule 7 amendment (user decision 12): the section 6(d) fit-set values quoted on this row are ", ADJ_LABEL))
mut_late("R7_br", mut_cell(BR, startsWith(BR$rule_ref, "W9F item 5: East Asian eGFR scale claim"), "qualifiers", sub(RULE7_LABEL, "", BR$qualifiers[startsWith(BR$rule_ref, "W9F item 5: East Asian eGFR scale claim")], fixed = TRUE)))
mut_late("BR_c1", mut_cell(BR, grepl("S9.corrob.PAM_T2D", BR$inputs, fixed = TRUE), "qualifiers", gsub(RULE7_LABEL, "", BR$qualifiers[grepl("S9.corrob.PAM_T2D", BR$inputs, fixed = TRUE)], fixed = TRUE)))
mut_late("BR_adj", mut_cell(BR, startsWith(BR$rule_ref, "W9F item 2a: section 6(d)"), "qualifiers", sub(ADJ_LABEL, "", BR$qualifiers[startsWith(BR$rule_ref, "W9F item 2a: section 6(d)")], fixed = TRUE)))
br_lab <- grepl(RULE7_LABEL, BR$qualifiers, fixed = TRUE)
add_assert("S9.check.rule7_label_on_W9_branches_rows",
           length(unlist(br_uncl)) == 0L && all(vapply(grep("^c1_label_", names(TRAVEL_REG), value = TRUE), function(n) any(vapply(br_hit, function(h) n %in% h, TRUE)), TRUE)) && all(br_lab == (lengths(br_nk) > 0L)) && sum(br_lab) >= 1L && all(grepl(ADJ_LABEL, BR$qualifiers, fixed = TRUE) == br_adj) && sum(br_adj) >= 1L &&
             br_lab[startsWith(BR$rule_ref, "W9F item 5: East Asian eGFR scale claim")],
           sprintf("%d of %d W9 branches.tsv rows quote not-key quantity ids (from `inputs` and from the registered travel texts) and exactly those carry the label; %d quoted ids are neither key nor not-key; %d rows quote adjudicator-computed 6(d) values and carry '%s'; rows with the label: %s", sum(br_lab), nrow(BR), length(unlist(br_uncl)),
                   sum(br_adj), ADJ_LABEL, paste(substr(BR$rule_ref[br_lab], 1, 40), collapse = " / ")))

# ---------------------------------------------------------------- 8b. convention questions for the user (W9F item 2d)
# One row per open definition/convention question raised by the Stage 1b adjudicators or VERIFIER_DEFECTS.md that touches
# a main-text slot. Every consequence is computed here from a pinned file; each question's source text is asserted present.
logline("--- convention questions for the user ---")
VD <- W9F_IN$verifier_defects$path
CQ <- list()
cq <- function(id, question, raised_by, slot, readA, consA, readB, consB, effect) {
  CQ[[length(CQ) + 1L]] <<- data.table(id = id, question = question, raised_by = raised_by, main_text_slot = slot,
    reading_A = readA, consequence_A = consA, reading_B = readB, consequence_B = consB, effect_on_branches = effect)
}
st <- function(key) { r <- sv(key); paste0(r$status, if (nzchar(r$alt_flag)) paste0(" ", r$alt_flag) else "") }
d6r <- function(fs, shift) D6R[D6R$fit_set == fs & D6R$build_shift_rows == shift & grepl("^\\(d\\)", D6R$ING1_version)]
# CQ01-02 section 6(d)
cq("CQ01", "Section 6(d) fit set: s_max over one panel per variant set (four fits) or every set on both panels (seven fits)? DECIDED 2026-10-09: seven fits primary, four fits a sensitivity",
   "S6 adjudicator (A_third: a3_allcand_fitset.tsv, a3_d_table.tsv); both lenses used seven fits", "renal OR sentence (D6 base phrase; 'the LD guard')",
   "seven fits (both Stage 1b lenses; the user's decision; primary)", S6D_SUMMARY("seven"), "four fits (the earlier implementer set; sensitivity)", S6D_SUMMARY("four"),
   sprintf("none: D6 '%s' and 'sensitive to %s' under every combination (S9.check.s6d_branch_levels_invariant); decided 2026-10-09", D6, paste(RENAL_LIST, collapse = ", ")))
quote_ok(VD, "the choice decides which side of 0.10 the value falls on", "VERIFIER_DEFECTS")
cq("CQ02", "NOTCH2-T2D window: keep or drop the UKB-PPP rows whose GRCh37 IDs lie in 1q21, about 25 Mb from the rest (build shift)?",
   "VERIFIER_DEFECTS.md S6/L1 and S6/L2 (S6H.d.NOTCH2_T2D.s_max); S6 adjudicator (a3_notch2_fits.tsv)", "renal OR sentence (via n_reclassified and the (d) CI)",
   "keep (implementer; the literal 'min/max ID pos37 +/- 1000' rule; user decision 2026-10-09: keep, disclosed)", sprintf("NOTCH2 s_max four fits %.4f, seven fits %.4f; n_reclassified %d / %d",
     d6r("four", "keep")$NOTCH2_s_max, d6r("seven", "keep")$NOTCH2_s_max, d6r("four", "keep")$n_reclassified, d6r("seven", "keep")$n_reclassified),
   "drop", sprintf("NOTCH2 s_max four fits %.4f, seven fits %.4f; n_reclassified %d / %d",
     d6r("four", "drop")$NOTCH2_s_max, d6r("seven", "drop")$NOTCH2_s_max, d6r("four", "drop")$n_reclassified, d6r("seven", "drop")$n_reclassified),
   "none on the branch levels; NOTCH2 crosses 0.10 only with seven fits and the rows kept, which is the decided reading (S9.check.s6d_notch2_only_seven_keep_crosses)")
# CQ03 coloc prior cap
quote_ok(VD, "does not say whether coloc.abf's internal adjust_prior cap applies", "VERIFIER_DEFECTS")
pk <- function(pr, rd) { r <- S3PK[S3PK$pair == pr & S3PK$maf_mode == "sub" & S3PK$reading == rd]; stopifnot(nrow(r) == 1L); r }
# W9F step 2b (D1): reading B reports EVERY 'full window' definition in BOTH S3 adjudicator files (S3ADJ; branches
# re-derived by W9, S9.check.s3_adj_branches_rederived), each with its PP.H4/PP.H3 and Rule 1 branch, and the effect
# column states, per definition, the change (or not) to PAM's section 3 sentence and to corroboration condition 2.
for (ph in S3_PHRASE_FRAGS) quote_ok(GOV$prespec_red$path, ph, "redacted prespec")
PAM_B0 <- s3_file_code(S3$PAM_T2D$branch)
FW <- cq03_full_windows(S3ADJ)
fw_nm <- s3_def_name(FW)
fw_mm <- lapply(seq_len(nrow(FW)), function(i) S3ADJ[S3ADJ$pair == "PAM_T2D" & S3ADJ$file == FW$file[i] & S3ADJ$scenario == FW$scenario[i]])
mm_txt <- if (all(vapply(fw_mm, function(x) length(unique(x$branch)) == 1L, TRUE)))
  # W9F step 2c: name the modes each definition was actually run under (the union over definitions overstated scope)
  paste0("MAF modes: each definition gives the same Rule 1 branch under every MAF mode the adjudicator ran for it (",
         paste(sprintf("%s: %s", fw_nm, vapply(fw_mm, function(x) paste(sort(unique(x$maf_mode)), collapse = ", "), "")), collapse = "; "), ")") else
  paste0("MAF modes: ", paste(sprintf("%s: %s", fw_nm, vapply(fw_mm, function(x) paste(sprintf("%s %s", x$maf_mode, x$branch), collapse = ", "), "")), collapse = "; "))
s3_union <- function(p) sort(unique(c(S3ADJ$branch[S3ADJ$pair == p], S3ADJ$branch_alt_h4only[S3ADJ$pair == p])))
oth <- c("RAD23B_T2D", "TGFB1_CAD", "PCSK9_CAD")
oth_same <- vapply(oth, function(p) identical(s3_union(p), s3_file_code(S3[[p]]$branch)), TRUE)
br_other <- vapply(oth, function(p) paste0(p, " {", paste(s3_union(p), collapse = ", "), "}"), "")
chg <- FW$branch != PAM_B0
grp <- split(fw_nm[chg], FW$branch[chg])
h4c <- FW$branch_alt_h4only != PAM_B0
grp_h4 <- split(fw_nm[h4c], FW$branch_alt_h4only[h4c])
c2_fw <- vapply(FW$branch, function(b) corroborated(TRUE, s3_label(b))$c2, TRUE); c2_0 <- corroborated(TRUE, s3_label(PAM_B0))$c2
cq03_eff <- paste0(
  if (!any(chg)) "none: " else "",
  sprintf("PAM section 3 sentence (W3 primary: branch %s, '%s'): ", PAM_B0, S3_PHRASE[[PAM_B0]]),
  paste(c(if (any(!chg)) sprintf("unchanged under %s", paste(fw_nm[!chg], collapse = ", ")),
          vapply(names(grp), function(b) sprintf("changes to branch %s ('%s') under %s", b, S3_PHRASE[[b]], paste(grp[[b]], collapse = ", ")), "")), collapse = "; "),
  if (any(h4c)) paste0("; under the PP.H4-only reading (PP.H4 >= 0.80 alone at 1e-6): ",
                       paste(vapply(names(grp_h4), function(b) sprintf("branch %s under %s", b, paste(grp_h4[[b]], collapse = ", ")), ""), collapse = "; ")) else "",
  if (all(c2_fw == c2_0)) sprintf("; corroboration: condition 2 is unchanged under every definition (%s under %s), so PAM-T2D stays %s and the R1 opening list does not change (S9.check.cond2_invariant_to_S3_conventions.PAM_T2D)",
                                  c2_0, { bb <- sort(unique(c(PAM_B0, FW$branch))); paste0(if (length(bb) > 1L) "branches " else "branch ", paste(bb, collapse = " and ")) },
                                  COR$verdict[COR$pkey == "PAM_T2D"]) else
    sprintf("; corroboration: condition 2 changes under %s (%s -> %s), so PAM-T2D's verdict and the R1 opening list would change", paste(fw_nm[c2_fw != c2_0], collapse = ", "), c2_0, !c2_0),
  if (all(oth_same)) "; RAD23B, TGFB1 and PCSK9: section 3 branch unchanged in every reading of both files" else
    sprintf("; section 3 branch changes in some reading for %s", paste(oth[!oth_same], collapse = ", ")))
cq("CQ03", "Section 3 priors: does coloc.abf's adjust_prior cap (p = 1/(nsnps + 1) when nsnps x p >= 1) apply, and is dataset 1 the analysed set or the whole exposure window?",
   "VERIFIER_DEFECTS.md S3/L2; S3 adjudicator (adj_S3_coloc_pkg.tsv, adj_S3_scenarios.tsv)", "section 3 sentences (PAM, RAD23B, TGFB1, PCSK9; PCSK9 control); corroboration condition 2",
   "declared p1 = p2 = 1e-4 on the analysed set (implementer; no cap triggered)",
   sprintf("PAM PP.H4 at 1e-6 %.3f (PP.H3 %.3f): branch %s", num(pk("PAM_T2D", "A_analysed")$H4_6), num(pk("PAM_T2D", "A_analysed")$H3_6), pk("PAM_T2D", "A_analysed")$branch),
   sprintf("cap applied with a full window: %d definitions in the two S3 adjudicator files", nrow(FW)),
   paste(c("PAM, primary MAF mode (substitute); every full-window definition, Rule 1 = PP.H4 >= 0.80 and >= 5 x PP.H3 at both priors",
           vapply(seq_len(nrow(FW)), function(i) s3_def_segment(FW[i]), ""), mm_txt,
           sprintf("other pairs, every row of both files (Rule 1 and PP.H4-only branches): %s", paste(br_other, collapse = "; "))), collapse = " || "),
   cq03_eff)
# CQ04 DIAMANTE EAF printed 0/1
quote_ok(VD, "Neither reading (substitute vs drop) is declared", "VERIFIER_DEFECTS")
sc <- function(pr, md) { r <- S3SC[S3SC$pair == pr & S3SC$maf_mode == md & S3SC$scenario == "declared"]; stopifnot(nrow(r) == 1L)
  list(h4_5 = num(r$H4_5), h3_5 = num(r$H3_5), h4_6 = num(r$H4_6), h3_6 = num(r$H3_6), nsnp = as.integer(r$nsnp),
       br = s3_branch(num(r$H4_5), num(r$H3_5), num(r$H4_6), num(r$H3_6))) }
eaf_txt <- function(mode) sprintf("PAM PP.H4 %.4f / %.4f (1e-5 / 1e-6), branch %s; RAD23B PP.H3 %.3f / %.3f, branch %s (nsnp PAM %d / RAD23B %d)",
  sc("PAM_T2D", mode)$h4_5, sc("PAM_T2D", mode)$h4_6, sc("PAM_T2D", mode)$br, sc("RAD23B_T2D", mode)$h3_5, sc("RAD23B_T2D", mode)$h3_6, sc("RAD23B_T2D", mode)$br,
  sc("PAM_T2D", mode)$nsnp, sc("RAD23B_T2D", mode)$nsnp)
cq("CQ04", "Section 3 outcome MAF where DIAMANTE prints EAF as 0 or 1: substitute the UKB-PPP MAF, drop the SNPs, or clamp?",
   "VERIFIER_DEFECTS.md S3/L1 (__alt rows); W3 README deviation 1; S3 adjudicator (adj_S3_scenarios.tsv)", "section 3 sentences for PAM and RAD23B",
   "substitute (implementer; primary)", eaf_txt("sub"), "drop (L1 __alt) / clamp", paste0("drop: ", eaf_txt("drop"), " | clamp: ", eaf_txt("clamp")),
   "none (branches identical); the drop version was seen before substitution was chosen: that disclosure travels with RAD23B's printed PP.H3 (W9F item 2c)")
# CQ05 FGFR4 TPMI SE convention
fg <- C1ADJ[C1ADJ$pair == "FGFR4_T2D"]; stopifnot(nrow(fg) == 1L)
cq("CQ05", "FGFR4-T2D (TPMI, SAIGE): which p classifies the check, and what SE convention do the SPA-adjusted rows use?",
   "C1 adjudicator (adj_C1.tsv, adj_C1_tpmi_consistency.json); VERIFIER_DEFECTS.md C1/L1, C1/L2", "R1 opening (corroboration of FGFR4-T2D); East Asian check sentence",
   "the resource's own p (declared: Principle 0, OCP section 5; C1)", sprintf("p %s: class %s; 0.05/6 %s", fg$p_resource, gsub("_", " ", fg$class_R), if (fg$bonf6_R == "True") "met" else "not met"),
   "first-order Wald p from beta/SE",
   sprintf("p %s (%.2f times the resource p; in SPA-adjusted rows the SE matches the one-sided p in a share %.2f of %d sampled rows): class %s; 0.05/6 %s",
           format(num(fg$p_wald1), digits = 3), num(fg$p_wald1) / num(fg$p_resource), TPMI_J$spa$frac_SE_matches_one_sided, TPMI_J$spa$n,
           gsub("_", " ", fg$class_W), if (fg$bonf6_W == "True") "met" else "not met"),
   "none (class and 0.05/6 result identical; FGFR4-T2D corroboration unchanged)")
# CQ06 UMOD palindromic signal 5; CQ13 optimize tolerance
uq <- function(pat, ref) { r <- A6U[startsWith(A6U$reading, pat) & startsWith(A6U$refine, ref)]; stopifnot(nrow(r) == 1L); r }
p_alt <- num(uq("__alt", "exact")$p)
p_fam_alt <- replace(fam$p_unadj, fam$member == "umod_egfr", p_alt); h_alt <- holm(p_fam_alt)
same_codes <- all(mapply(function(pu, ph, i) identical(nominal_literal(pu, ph)$code, fam$nominal_code[i]) && identical(nominal_I1(pu, ph)$code, fam$nominal_code_I1[i]),
                         p_fam_alt, h_alt, seq_len(nrow(fam))))
cq("CQ06", "UMOD all-signal Q: the declared five signals include rs62033284 (A/T, palindromic), but the '503' label means palindromes dropped. Keep signal 5?",
   "S6 adjudicator (a3_umod_q.tsv: __alt reading)", "UMOD control sentence; Rule 4 family (UMOD-eGFR member)",
   "all five signals (declared set; implementer)", sprintf("Q %.3f, df %s, p %s (gls_q.R refine) / %s (exact optimum)", num(uq("all", "coarse")$Q), uq("all", "coarse")$df,
                                                             format(num(uq("all", "coarse")$p), digits = 4), format(num(uq("all", "exact")$p), digits = 4)),
   "signal 5 dropped (the '503 = p1np' label)", sprintf("Q %.3f, df %s, p %s (exact); UMOD Holm p with this input %s", num(uq("__alt", "exact")$Q), uq("__alt", "exact")$df,
                                                        format(p_alt, digits = 4), format(h_alt[fam$member == "umod_egfr"], digits = 3)),
   sprintf("none: 'behaved as a control' under both (p >= 0.05); every member's literal and I1 'nominal' code unchanged: %s", same_codes))
quote_ok(VD, "does not give the interval optimize searches", "VERIFIER_DEFECTS")
cq("CQ13", "UMOD Q minimisation: refine at R optimize's default tolerance (gls_q.R) or to the exact optimum; over what interval?",
   "VERIFIER_DEFECTS.md S6/L1, S6/L2; comparator S6H.umod.Q.p; S6 adjudicator (a3_umod_q.tsv)", "UMOD control sentence; Rule 4 family",
   "gls_q.R refine (implementer)", sprintf("p %s", format(num(uq("all", "coarse")$p), digits = 6)), "exact optimum (both lenses)",
   sprintf("p %s; comparator status %s", format(num(uq("all", "exact")$p), digits = 6), st("S6H.umod.Q.p")),
   "none (p >= 0.05 either way; Holm p 1 either way)")
# CQ07 the '503' label
quote_ok(VD, "two of which (rs1126464, rs72811595) are C/G palindromes", "VERIFIER_DEFECTS")
quote_ok(VD, "which that panel excludes", "VERIFIER_DEFECTS")
cq("CQ07", "Panel label '503' = 1000G EUR with palindromes dropped (p1np), but declared quantities involve palindromic SNPs. Which governs?",
   "VERIFIER_DEFECTS.md S1/L1, S1/L2 (S1H.r2), S5/L1, S5/L2 (S5H.Q.*, S5H.sep_maxr2)", "RAD23B first mention (r2); DPEP1 Q sentences; UMOD (CQ06)",
   "the declared variant (implementer override: palindromes kept)",
   sprintf("RAD23B lead-proxy r2 %.4f (%s); DPEP1 five-signal Q p cys-gly %s (%s), eGFR %s (%s)", getq(R$W1, "W1", "S1.rad23b.r2_plink"), st("S1H.r2"),
           format(fam$p_unadj[fam$member == "dpep1_cysgly"], digits = 3), st("S5H.Q.cysgly.p"), format(fam$p_unadj[fam$member == "dpep1_egfr"], digits = 3), st("S5H.Q.egfr.p")),
   "the label (palindromes dropped)", "RAD23B r2 undefined (its lead rs7035725 is palindromic); the DPEP1 Q without rs1126464 and rs72811595 was not computed by any lens; UMOD: CQ06",
   "unknown for DPEP1 (not computed); RAD23B's printed r2 exists only under the override")
# CQ08 the cut; CQ09 the proxy rule
quote_ok(VD, "gives no alpha, no inequality direction", "VERIFIER_DEFECTS")
quote_ok(W9F_IN$W1_readme$path, "BH over 15,584 cells, with the 1,063 untested cells set to p = 1", "W1 README")
w1_mov <- B$W1[startsWith(B$W1$rule_ref, "Rule 2: a pair whose class differs")]; stopifnot(nrow(w1_mov) == 1L)
mv <- regmatches(w1_mov$filled_template, regexec("RAD23B-T2D: ([A-Za-z]+) under first-order SE, ([A-Za-z]+) under second-order SE", w1_mov$filled_template))[[1]]
stopifnot(length(mv) == 3L)
cutp <- function(q) getq(R$W1, "W1", q, "p")
yn <- function(q) if (getq(R$W1, "W1", q, numeric = FALSE) == "TRUE") "meets" else "does not meet"
cq("CQ08", "Rule 1 'the cut' for the RAD23B first mention: which alpha, which inequality, Bonferroni or BH over 15,584 tests?",
   "VERIFIER_DEFECTS.md S1/L1, S1/L2 (S1H.cut.*)", "RAD23B first-mention sentence (meets the cut at the lead / proxy / under second-order SE)",
   "Bonferroni 0.05/15,584, strict '<' (implementer; both lenses)",
   sprintf("lead %s (p %s), proxy %s (p %s), second-order %s (p %s); comparator %s / %s / %s",
           yn("S1.rad23b.cut.lead"), format(cutp("S1.rad23b.cut.lead"), digits = 3), yn("S1.rad23b.cut.proxy"), format(cutp("S1.rad23b.cut.proxy"), digits = 3),
           yn("S1.rad23b.cut.second_order"), format(cutp("S1.rad23b.cut.second_order"), digits = 3), st("S1H.cut.lead"), st("S1H.cut.proxy"), st("S1H.cut.second")),
   "BH over the family",
   sprintf("W1's SE-convention table classes RAD23B-T2D (proxy) as %s under first-order SE and %s under second-order SE ('FDR' = BH-significant over the 15,584 cells but not Bonferroni; W1 README), so under a BH reading the pair meets the cut at the proxy under both SE orders; the lead is not a family member (not evaluable); no lens computed BH (the 15,584-p vector is not in the keys' resources)",
           mv[2], mv[3]),
   "the clause 'does not meet it under second-order SE' holds under the Bonferroni reading only")
quote_ok(VD, "do not name the proxy or define the rule", "VERIFIER_DEFECTS")
quote_ok(VD, "whether the second-order Wald ratio is taken at the lead or at the proxy", "VERIFIER_DEFECTS")
cq("CQ09", "The proxy rule is a declared primary but is not defined; is the second-order cut taken at the proxy or at the lead?",
   "VERIFIER_DEFECTS.md S1/L1, S1/L2 (S1H.proxy_p, S1H.p_second_order)", "RAD23B first mention ('proxy' label; second-order cut)",
   "proxy (implementer; outcome at the proxy, exposure at the lead)", sprintf("second-order p %s (%s)", format(cutp("S1.rad23b.cut.second_order"), digits = 3), st("S1H.p_second_order")),
   "lead", "not computed by any lens (no __alt)", "first-order p unaffected; the second-order cut at the lead is unknown")
# CQ10 Rule 4 nominal
dif_rows <- IR[IR$row_kind == "claim" & IR$a1_front_caveat != IR$a1_front_caveat_I1]
cq("CQ10", "Rule 4 'nominal': literal (every statement with Holm p >= 0.05) or applied (I1: only where the unadjusted p < 0.05)?",
   "W9 interpretation I1; STAGE1B_BRIEFS_v3 W9F item 6", "every Rule 4 family sentence (ARG1 x4, TGFB1, DPEP1 x2, UMOD)",
   "literal (default column)", sprintf("'nominal' on %d members: %s", sum(fam$nominal_attach), paste(fam$label[fam$nominal_attach], collapse = "; ")),
   "applied I1 (printed beside)", sprintf("'nominal' on %d member(s): %s", sum(fam$nominal_attach_I1), paste(fam$label[fam$nominal_attach_I1], collapse = "; ")),
   sprintf("wording only; the A1 front caveat differs in %d claim rows: %s", nrow(dif_rows), paste(sprintf("%s row %d (%s)", dif_rows$source, dif_rows$source_row, dif_rows$rule_ref), collapse = "; ")))
# CQ11-12, CQ23 s computation details
k2 <- A6K$detail[A6K$check == "K2_four_fit_max_equals_impl"]; k3 <- A6K$detail[A6K$check == "K3_seven_fit_max_equals_lens_primary"]
quote_ok(VD, "I dropped all copies", "VERIFIER_DEFECTS"); quote_ok(VD, "accepted only identical or swapped", "VERIFIER_DEFECTS")
cq("CQ11", "Section 6(d) s computation, steps (i) and (iii): drop every copy or only later copies of a repeated position + allele pair; how are indel alleles matched?",
   "VERIFIER_DEFECTS.md S6/L1, S6/L2", "renal OR sentence (via s_max)", "implementer's harmonised tables", "s_max as printed by W6b",
   "both lenses (every copy dropped; indels identical or swapped only)", "s_max as printed by the lenses",
   sprintf("none found: the S6 adjudicator shows the fit set alone explains every s_max difference (K2 four-fit max vs implementer: %s; K3 seven-fit max vs lenses: %s)", k2, k3))
quote_ok(VD, "the other basis changes protein n by <0.1%", "VERIFIER_DEFECTS")
cq("CQ12", "Section 6(d) step (viii): is the protein-n MAF > 0.05 filter on the UKB-PPP A1FREQ or on the reference-panel frequency?",
   "VERIFIER_DEFECTS.md S6/L1, S6/L2", "renal OR sentence (via s_max)", "UKB-PPP A1FREQ (both lenses)", "s_max as printed",
   "reference-panel frequency", "protein n changes by < 0.1% (L1's report; not computed by W9)", "none found (K2/K3 as CQ11)")
quote_ok(VD, "There is no rule for missing genotypes in the 503 panel", "VERIFIER_DEFECTS")
cq("CQ23", "Section 6(d) step (vi): how are missing genotypes in the 503 panel handled in the Pearson r?",
   "VERIFIER_DEFECTS.md S6/L1, S6/L2", "renal OR sentence (via s_max)", "pairwise-complete (L2)", "no fit's SNPs had a missing call (both lenses)",
   "mean-imputed or dropped", "identical on these inputs", "none (latent)")
# CQ14-16 section 4
quote_ok(VD, "undeclared 489-panel minor-allele-frequency filter", "VERIFIER_DEFECTS")
cq("CQ14", "Section 4: the 489 SNP set (LOC_ARG1_6q23_snpset_p2np) has no stated construction rule (consistent with an undeclared MAF filter)",
   "VERIFIER_DEFECTS.md S4/L1 (S4H.orn.PPH4_A.489)", "ornithine sentence (Rule 3(iv) colocalization clause)", "the file as given (implementer; both lenses)",
   sprintf("ornithine PP.H4 at signal A: 503 %.3f, 489 %.3f (%s)", getq(R$W4a, "W4a", "S4.coloc_clause.ornithine.A.503.PPH4"), getq(R$W4a, "W4a", "S4.coloc_clause.ornithine.A.489.PPH4"), st("S4H.orn.PPH4_A.489")),
   "a rebuilt set", "not computed", "none on the clause: it fires on the 503 panel alone (PP.H4 < 0.80)")
quote_ok(VD, "which L1 does in R", "VERIFIER_DEFECTS"); quote_ok(VD, "set.seed(20261008)", "VERIFIER_DEFECTS"); quote_ok(VD, "set.seed(20261009)", "VERIFIER_DEFECTS")
fgs <- sv("S4H.t2d.share.fg_cond")
cq("CQ15", "Rule 3(ii) bootstrap: which random-number generator draws r* (LD in R, or general computation in Python)?",
   "VERIFIER_DEFECTS.md S4/L1", "T2D sentence ('sensitive to the reference LD in FinnGen')", "R Mersenne-Twister, seed 20261007 (implementer)",
   sprintf("FinnGen conditional share %s", fgs$impl), "the lenses' draws (L1 seed 20261008, L2 seed 20261009; both in R)", sprintf("shares %s / %s (%s)", fgs$L1, fgs$L2, fgs$status),
   "none: every share is below 0.90, so FinnGen fails Rule 3(ii) under each")
quote_ok(VD, "the window is not defined", "VERIFIER_DEFECTS")
s4n <- SURV[SURV$sec == "S4"]
cq("CQ16", "Section 4 join rule: over what window is 'a source key that occurs more than once' judged for the genome-wide outcome files?",
   "VERIFIER_DEFECTS.md S4/L2", "ARG1 T2D and test (iii) sentences", "the implementer's implicit window", "values as printed by W4a",
   "the lenses' windows", sprintf("%d of %d S4 headline keys reproduced by both lenses", sum(s4n$L1_ok == "True" & s4n$L2_ok == "True"), nrow(s4n)), "none found")
# CQ17 S3 join; CQ18-19 S5; CQ20 S7; CQ21 S2; CQ22 C1 precision
quote_ok(VD, "so the ambiguity had no effect here", "VERIFIER_DEFECTS")
cq("CQ17", "Section 3 join: may an outcome row be matched twice (exact and complement), and how are multi-base alleles complemented?",
   "VERIFIER_DEFECTS.md S3/L1, S3/L2; S3 adjudicator (adj_S3_join_audit.tsv)", "section 3 sentences", "one-to-one (implementer)", "PP as printed",
   "many-to-one / reverse complement", sprintf("complement candidates %d (base-wise) and %d (reverse), outcome rows used twice %d, over the four windows",
                                               sum(as.integer(S3J$complement_basewise_candidates)), sum(as.integer(S3J$complement_revcomp_candidates)), sum(as.integer(S3J$outcome_rows_used_twice))),
   "none (the case does not arise on these inputs)")
quote_ok(VD, "majority' without saying how an exact half counts", "VERIFIER_DEFECTS")
cq("CQ18", "DPEP1 attribution (Rule 3(i)/(iii)): which n enters s; does an exact half count as a majority; how are s-gated settings treated?",
   "VERIFIER_DEFECTS.md S5/L1, S5/L2", "DPEP1 attribution sentence ('unresolved')", "implementer's readings", "unresolved",
   "the lenses' readings (L1 wrote an __alt)", sprintf("comparator S5H.branch.attribution: %s (L1 %s, L2 %s)", st("S5H.branch.attribution"), sv("S5H.branch.attribution")$L1, sv("S5H.branch.attribution")$L2),
   "none ('unresolved' under every reading computed)")
quote_ok(VD, "I joined on chr16", "VERIFIER_DEFECTS")
cq("CQ19", "DPEP1 eGFR join: Pos_b37 + unordered alleles without a chromosome?", "VERIFIER_DEFECTS.md S5/L1", "DPEP1 eGFR Q sentence",
   "literal (no chromosome)", "could match another chromosome", "chromosome 16", sprintf("no other-chromosome rows at the five positions (L1); S5H.Q.egfr.p %s", st("S5H.Q.egfr.p")), "none")
cq("CQ20", "East Asian stratified exact test: how is a stratum with fewer than 2 rows handled (mantelhaen.test stops; the convolution passes it through)?",
   "VERIFIER_DEFECTS.md S7/L1, S7/L2", "D12-A sentence (stratified p)", "convolution", "p as printed", "mantelhaen.test", sprintf("identical here (every stratum has >= 2 rows); S7H.strat_p %s", st("S7H.strat_p")), "none (latent)")
quote_ok(VD, "I took B's pair from", "VERIFIER_DEFECTS")
s2n <- SURV[SURV$sec == "S2"]
cq("CQ21", "TGFB1 signal B (rs12983775): its allele pair is not declared; where does it come from?", "VERIFIER_DEFECTS.md S2/L1", "TGFB1 section 2 sentence",
   "UKB-PPP ID at 19:41857439 (L1)", sprintf("%d of %d S2 headline keys reproduced by both lenses", sum(s2n$L1_ok == "True" & s2n$L2_ok == "True"), nrow(s2n)), "other source", "not computed", "none found")
quote_ok(VD, "stores beta (and sebeta and af) to two decimals", "VERIFIER_DEFECTS")
cq("CQ22", "FGFR4 exposure beta: the definition asks for full precision, but the named source stores beta, SE and AF to 2 significant figures",
   "VERIFIER_DEFECTS.md C1/L1, C1/L2", "East Asian check sentence; FGFR4-T2D corroboration", "the source as stored", sprintf("beta_exp %s; Wald %s (%s)", format(getq(C1R, "C1", "C1.FGFR4_T2D.beta_exp")),
   format(getq(C1R, "C1", "C1.FGFR4_T2D.wald"), digits = 4), st("C1H.FGFR4_T2D.wald")), "a more precise exposure value", "no such value in the declared source",
   "none on the class (concordant under both p readings, CQ05); the ratio inherits the rounding")
CQ <- rbindlist(CQ)[order(id)]
stopifnot(!anyDuplicated(CQ$id))
# W9F step 2b (D1): CQ03 must report every full-window definition of both S3 adjudicator files and state its effect
mut_late("Q_cq03", mut_cell(CQ, CQ$id == "CQ03", "effect_on_branches",
                            "none under Rule 1; corroboration condition 2 unchanged for all four pairs (S9.check.cond2_invariant_to_S3_conventions.*); PAM PP.H4 at 1e-6 crosses 0.80 only with both windows full"))
# ---------------------------------------------------------------- 8c. Stage 2 opening fix 6: disposition of every verifier-reported definition defect (H2 D8)
# One row per defect line of VERIFIER_DEFECTS.md (a line starting '- **'), in file order. `key_prefix` is what the line's bold key text must start with, so a
# shifted mapping fails the count assert. Disposition = a convention_questions.tsv id, METHODS_ONLY, or 'fixed (where)'.
METHODS_ONLY <- "Methods wording only (no consequence on these inputs)"
vdl <- readLines(VD, encoding = "UTF-8", warn = FALSE)
cur <- ""; VDT <- list()
for (l_ in vdl) {
  h_ <- regmatches(l_, regexec("^## ((?:S[1-7]|C1)/L[12])\\b", l_, perl = TRUE))[[1]]
  if (length(h_) == 2L) cur <- h_[2]
  d_ <- regmatches(l_, regexec("^- \\*\\*(.+?)\\*\\*: (.*)$", l_, perl = TRUE))[[1]]
  if (length(d_) == 3L) VDT[[length(VDT) + 1L]] <- data.table(section_lens = cur, reported_key = d_[2], defect = sub("^(.{1,170}[.;:]|.{170}).*$", "\\1", d_[3], perl = TRUE))
}
VDT <- rbindlist(VDT); VDT[, ord := .I]
SCH <- "comparator schema gap: compare_col 'branch' has no column in the mine.tsv schema (level written to both value and category); a Stage 1b procedure item, no manuscript text"
RCN <- "the RESOURCES note labels the DIAMANTE row count '(wc -l)' but it is the data-row count (off by the header line); the data-row reading was asserted; correct the RESOURCES note at its next re-issue (RESOURCES.json is pinned here, not edited)"
DISP <- rbindlist(list(
  list("S1H.r2", "CQ07", ""), list("S1H.p_second_order", "CQ09", ""), list("S1H.cut.lead", "CQ08", ""), list("S1H.proxy_p", "CQ09", ""),
  list("S1H.lead_p", METHODS_ONLY, "exposure-table provenance (CSV TAB_instruments vs the .rds 24f2b98f) is inferred from the file name; partly closed because C1's beta_exp equal the OCP's printed values; state the provenance in Methods"),
  list("S1H.cut.lead", "CQ08", ""), list("S1H.cut.proxy", "CQ08", ""), list("S1H.cut.second", "CQ08; CQ09", ""), list("S1H.p_second_order", "CQ09", ""),
  list("S1H.r2", "CQ07", ""), list("S1H.lead_p", METHODS_ONLY, "same provenance statement as S1/L1 (exposure table vs the .rds md5); the p value does not depend on the exposure row"),
  list("S2H.ratio.finngen.rs73045269", "CQ21", ""), list("S2H.branch", METHODS_ONLY, SCH),
  list("S3H.PAM_T2D.PPH4.1e-5", "CQ04", ""), list("S3H.PAM_T2D.branch / S3H.RAD23B_T2D.branch", METHODS_ONLY, SCH),
  list("S3H.<pair>.branch (all four)", "CQ04", "one __alt slot cannot hold the drop-reading branch; CQ04 reports the branch under drop and clamp beside the substitute reading"),
  list("S3H.<pair>.PPH4/PPH3 (all 16 numeric keys)", "CQ17", ""),
  list("S3H.PAM_T2D.branch; S3H.RAD23B_T2D.branch", METHODS_ONLY, SCH),
  list("all 16 S3H.*.PPH3/PPH4 keys", "CQ03", ""), list("all 16 S3H.*.PPH3/PPH4 keys (join rule)", "CQ17", ""),
  list("S3H.PAM_T2D.*, S3H.RAD23B_T2D.* (row-count assertion)", METHODS_ONLY, RCN),
  list("S4H.t2d.cells_same_direction", METHODS_ONLY, SCH),
  list("S4H.orn.PPH4_A.489", "CQ14", ""), list("S4H.t2d.share.fg_cond", "CQ15", ""),
  list("S4H.t2d.p489", METHODS_ONLY, "the key's resource list omits LOC_ARG1_6q23_snpset_p2np, so the .i8 route cannot be followed from the list alone; list the file in the resource map"),
  list("S4H.t2d.diff (and every FinnGen/DIAMANTE key", "CQ16", ""),
  list("S5H.Q.cysgly.p", "CQ07", ""), list("S5H.Q.egfr.p", "CQ07; CQ19", ""), list("S5H.sep_maxr2", "CQ07", ""), list("S5H.branch.attribution", "CQ18", ""),
  list("S5H.Q.cysgly.p", "CQ07", ""), list("S5H.Q.egfr.p", "CQ07", ""),
  list("S5H.branch.attribution", METHODS_ONLY, "the p1fr raw-to-harmonised join is spelled out only for DPEP1 protein and cysteinylglycine; the same convention was applied to the other eight raw files; state it in Methods"),
  list("S5H.branch.attribution", METHODS_ONLY, "no grid s value is NA in the DPEP1-cysteinylglycine rows, so the missing-s rule did not bind; state the rule in Methods"),
  list("S6H.d.NOTCH2_T2D.s_max", "CQ02", ""), list("S6H.d.CELSR2_eGFR.s_max", "CQ12", ""), list("S6H.d.GSTA1_eGFR.s_max", "CQ11", ""),
  list("S6H.d.UMOD_eGFR.s_max", "CQ23", ""), list("S6H.umod.Q.p", "CQ13", ""),
  list("S6H.d.ING1_CAD.s_max", METHODS_ONLY, "only the optional outcome-only s (notes only) is affected, by comma-separated rsID cells; ING1-CAD s_max is not computable in any computation"),
  list("S6H.d.<candidate>.s_max (all 15 computable; step i)", "CQ11", ""), list("S6H.d.<candidate>.s_max (step viii)", "CQ12", ""),
  list("S6H.d.<candidate>.s_max (step vi)", "CQ23", ""), list("S6H.d.NOTCH2_T2D.s_max (step i)", "CQ02", ""),
  list("S6H.d.<candidate>.s_max (step iii)", "CQ11", ""), list("S6H.umod.Q.p", "CQ13", ""),
  list("S7H.strat_p", "CQ20", ""),
  list("S7H.strat_p", METHODS_ONLY, "the resource map has two entries for this file that disagree on its format (tab vs comma); the file is comma-separated; correct the RESOURCES note at its next re-issue"),
  list("S7H.ckb_median_sd", METHODS_ONLY, "the RESOURCES note on the folder-digest format does not reproduce the stored digest; the key text's corrected form does; correct the RESOURCES note at its next re-issue"),
  list("S7H.strat_p", "CQ20", ""),
  list("S7H.strat_p", METHODS_ONLY, "encoding of eas_direction_same ('True'/'False' text, empty for non-evaluable rows) assumed equal to qc_evaluable's; state the encoding in Methods"),
  list("C1H.FGFR4_T2D.wald", "CQ22", ""), list("C1H.PAM_T2D.class", METHODS_ONLY, SCH),
  list("C1H.FGFR4_T2D.p_resource", "CQ05", "the two alternative readings (first-order Wald p; SAIGE p.value.NA) collide on one __alt slot; CQ05 reports both p readings"),
  list("C1H.PAM_T2D.p_resource", METHODS_ONLY, RCN),
  list("C1H.PAM_T2D.class; C1H.RAD23B_T2D.class", METHODS_ONLY, SCH),
  list("C1H.FGFR4_T2D.wald", "CQ22", ""), list("C1H.FGFR4_T2D.p_resource", "CQ05", "the two alternative readings collide on one __alt slot; CQ05 reports both p readings"),
  list("C1H.PAM_T2D.p_resource; C1H.RAD23B_T2D.p_resource", METHODS_ONLY, RCN),
  list("C1H.PAM_T2D.wald; C1H.RAD23B_T2D.wald", METHODS_ONLY, "cosmetic: a corrupted beta character in a quoted pre-declaration fragment; the meaning is recoverable")),
  use.names = FALSE)
setnames(DISP, c("key_prefix", "disposition", "detail"))
mut_late("D_count", { DISP <- DISP[-1L]; MUT_APPLIED <- TRUE })
n_bullets <- sum(startsWith(vdl, "- "))
add_assert("S9.check.dispositions_count_equals_defect_lines",
           nrow(DISP) == nrow(VDT) && nrow(VDT) == n_bullets && n_bullets > 0L && all(startsWith(VDT$reported_key, DISP$key_prefix)),
           sprintf("VERIFIER_DEFECTS.md has %d defect lines (%d bullet lines), %d dispositions; each line's reported key starts with the key its disposition was written for", nrow(VDT), n_bullets, nrow(DISP)))
mut_late("D_cq", { DISP$disposition[5L] <- "CQ99"; MUT_APPLIED <- TRUE })
disp_ids <- strsplit(DISP$disposition, "; ", fixed = TRUE)
cq_used <- sort(unique(unlist(disp_ids))[unique(unlist(disp_ids)) %in% CQ$id])
add_assert("S9.check.dispositions_reference_existing_CQ",
           all(vapply(disp_ids, function(x) all(x %in% CQ$id) || identical(x, METHODS_ONLY) || (length(x) == 1L && startsWith(x, "fixed (")), TRUE)),
           sprintf("every disposition is a convention_questions.tsv id (used: %s), '%s' or 'fixed (where)'; %d map to CQ ids, %d to Methods wording only, %d fixed",
                   paste(cq_used, collapse = ", "), METHODS_ONLY, sum(vapply(disp_ids, function(x) all(x %in% CQ$id), TRUE)),
                   sum(DISP$disposition == METHODS_ONLY), sum(startsWith(DISP$disposition, "fixed ("))))
if (!nzchar(MUT)) fwrite(cbind(VDT[, .(ord, section_lens, reported_key, defect)], DISP[, .(disposition, detail)]), file.path(OUT_W, "defect_dispositions.tsv"), sep = "\t", quote = TRUE)
if (nzchar(MUT) && !MUT_APPLIED) stop(sprintf("MUTATION %s was not applied (unknown id or hook not reached)", MUT))
rq <- chk_cq03(CQ[CQ$id == "CQ03"], S3ADJ, PAM_B0)
add_assert("S9.check.CQ03_reports_every_full_window_reading", rq$ok, rq$detail)
if (!nzchar(MUT)) fwrite(CQ, file.path(OUT_W, "convention_questions.tsv"), sep = "\t", quote = TRUE)
logline(sprintf("convention_questions.tsv: %d questions", nrow(CQ)))
# ---------------------------------------------------------------- 9. assemble results, write, pin
for (a in ASSERTS) addr(a$id, "W9 re-derivation / conformance assert", "assert", value = a$status, notes = a$derivation)
RESULTS <- rbindlist(RES_ROWS)
stopifnot(!anyDuplicated(RESULTS$quantity_id))
bad <- RESULTS[RESULTS$value %in% c("", "NA", "NaN")]
if (nrow(bad)) stop("results rows with an empty or NA value: ", paste(bad$quantity_id, collapse = ", "))
if (!nzchar(MUT)) {
  fwrite(RESULTS, file.path(OUT_W, "results.tsv"), sep = "\t", quote = TRUE)
  fwrite(BR, file.path(OUT_W, "branches.tsv"), sep = "\t", quote = TRUE)
  fwrite(IR, file.path(OUT_W, "branches_integrated.tsv"), sep = "\t", quote = TRUE)
  fwrite(rbindlist(LOOKUPS$rows), file.path(OUT_W, "lookups.tsv"), sep = "\t", quote = TRUE)
}
fwrite(rbindlist(ASSERTS), file.path(OUT_W, "asserts.tsv"), sep = "\t", quote = TRUE)
logline(sprintf("lookups: %d quantity lookups, each returned exactly 1 row (join yield 100%%; 0 null values)", length(LOOKUPS$rows)))
logline(sprintf("asserts: %d PASS, %d ASSERT_FAILED", sum(sapply(ASSERTS, function(a) a$status == "PASS")), sum(sapply(ASSERTS, function(a) a$status != "PASS"))))

# end checks: inputs unchanged, staging unchanged
for (w in names(WS)) { assert_md5(file.path(S1, WS[[w]]$dir, "results.tsv"), WS[[w]]$results, paste0("END ", w, "/results.tsv"))
                       assert_md5(file.path(S1, WS[[w]]$dir, "branches.tsv"), WS[[w]]$branches, paste0("END ", w, "/branches.tsv")) }
for (x in names(W9F_IN)) assert_md5(W9F_IN[[x]]$path, W9F_IN[[x]]$md5, paste0("END ", x))
SNAP1 <- snap()
if (!nzchar(MUT)) fwrite(SNAP1, file.path(OUT_W, "staging_snapshot_end.tsv"), sep = "\t")
stopifnot(identical(SNAP0, SNAP1))
logline(sprintf("staging snapshot end: %d files, identical to start", nrow(SNAP1)))
if (nzchar(MUT)) { logline(sprintf("MUTATION RUN %s done: %d ASSERT_FAILED", MUT, sum(sapply(ASSERTS, function(a) a$status != "PASS")))); quit(save = "no", status = 0) }
for (f in c("results.tsv", "branches.tsv", "branches_integrated.tsv", "corroboration.tsv", "dpep1_s_merged.tsv", "s6d_fitset_readings.tsv", "convention_questions.tsv", "lookups.tsv", "asserts.tsv"))
  logline(sprintf("RESULTS PINNED  %-26s md5 %s  rows %d", f, md5(file.path(OUT_W, f)), nrow(fread(file.path(OUT_W, f), sep = "\t", colClasses = "character"))))
logline(sprintf("W9 w9_main.R done in %.1f s", as.numeric(difftime(Sys.time(), T0, units = "secs"))))
