# W9 (Stage 1 integration) - shared helpers, pinned md5s, readers.
# Sourced by w9_main.R, w9_mutation.R and w9_compare_seen.R. Writes nothing by itself.
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
suppressPackageStartupMessages(library(data.table))

V3   <- .P31R_V3
S1   <- file.path(V3, "R1_stage1")
OUT  <- file.path(S1, "S1R1_R2_10_integration")
RES  <- file.path(S1, "S1R1_R2_00_resources")
STAGING <- file.path(V3, "V4_local/staging/l3-finemap")

# ---- governing documents (md5s passed in the workflow prompt) ----
GOV <- list(
  briefs        = list(path = file.path(RES, "STAGE1_BRIEFS.md"),                          md5 = "40d0d52f7055542d774c03f1edf04e00"),
  prespec_red   = list(path = file.path(RES, "P31R_R2_STAGE1_PRESPEC_REDACTED.md"),         md5 = "a281c3464619b1b67c28e13e7392c41c"),
  amend1_red    = list(path = file.path(RES, "P31R_R2_STAGE1_PRESPEC_AMEND1_REDACTED.md"),  md5 = "f265d600973c43239bb8377e0cd39d33"),
  prespec_full  = list(path = file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"),        md5 = "3319dc39740d9f5b3b8822133d9108b7"),
  amend1_full   = list(path = file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), md5 = "377508832ea9bd384133e869c456c909"),
  gls_q         = list(path = file.path(RES, "gls_q.R"),                                    md5 = "b7390fdf5b6f2f9c1d30c298d6afdafc"),
  # W9F step 0 (2026-10-09): RESOURCES.json was re-issued during Stage 1b (6d212381 -> 968b814e -> 546b769a; 311 -> 324
  # entries). W9 reads one entry (R1_outcome_checks_tsv); its path, md5 and role are asserted unchanged in w9_main.R.
  resources     = list(path = file.path(RES, "RESOURCES.json"),                             md5 = "546b769ae93a8f566eec2b92f4be1575"),
  outcome_prespec = list(path = file.path(S1, "S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), md5 = "c9079f2de1ffe75ce31f8a1d597d1fda")
)

# ---- workstream outputs read by W9 (md5 pinned at W9 start, 2026-10-08 18:3x, after the S4/S6 fixer rounds) ----
# W9F step 0 (2026-10-09): W6b and W7 branches.tsv were rewritten in place by the Stage 1b key-id normalisers (inputs
# column only). Re-pinned after w9f_step0_pins.R proved, by byte reconstruction, that no other cell changed
# (pin_proofs.tsv): W6b 8af0cfe4db091fccc04879366a52a9fe -> d8b9beb10ec535f5a49d9e43bf59a095;
#                   W7  ff6a4d7ea6987c571a1a375abbe1a5c0 -> bb128be6fa4da3c2724c315b81529f2a.
WS <- list(
  W1  = list(dir = "S1R1_R2_01_rad23b",        results = "1fe7a0f4dcd2210160fae96b72586024", branches = "42b0780d1824d77e3f719f8c0c3dbe95"),
  W2  = list(dir = "S1R1_R2_02_tgfb1",         results = "306bf3a7821495f1b3f7369d73df81e7", branches = "de44008ce761bc3b65fdc3db1d1863af"),   # Stage 2 step 2 (seven-fit s_max); was a860054f / 566257a2
  W3  = list(dir = "S1R1_R2_03_extcoloc",      results = "d1f25136b28f38cc4620b46700d036d8", branches = "82c3f175cd9f3955adfd4ba232f64071"),
  W4a = list(dir = "S1R1_R2_04a_arg1grid",     results = "f4c5104dce1f1822448a6a7bdced1636", branches = "cb24233b571bd9922eebb6a83a7a04a4"),
  W4b = list(dir = "S1R1_R2_04b_arg1sets",     results = "d77fea1ac7a811492dcf1140dae068c9", branches = "68cc87393485c037608686df77c1bf47"),   # Stage 2 step 2 (seven-fit s_max); was e995af3e / 172bcc96; branches re-pinned in the fix round (was 769ab8f9: row 13 text now names the seven-fit argmax)
  W5  = list(dir = "S1R1_R2_05_dpep1",         results = "ef11ffa688693bd6ef4528fcb8a5f412", branches = "fc5c3873758545f40db6525cf8ff3192"),
  W6a = list(dir = "S1R1_R2_06a_renal",        results = "2a85c748f20495b74227cb4780a5f759", branches = "bfe38e2f73bf8f430547bf698f951132"),
  W6b = list(dir = "S1R1_R2_06b_renal_d_umod", results = "1f32cd299a3b1b7e9d6d82dc77863226", branches = "c85c9159cffa2ca197257394a62f403d"),   # Stage 2 step 2 (seven-fit s_max, key-id normaliser re-applied by R/10b); was eb8081b6 / d8b9beb1
  W7  = list(dir = "S1R1_R2_07_ancestry",      results = "d1d3e66459f79fa030357e3280ba0346", branches = "bb128be6fa4da3c2724c315b81529f2a")
)
# ---- W9F (Stage 1b-C) inputs: C1's external checks, the Stage 1b comparator copies and the adjudication files ----
V1B <- file.path(S1, "S1R1_R2_1b_verify")
A6  <- file.path(V1B, "S6/adjudication/A_third")
W9F_IN <- list(
  C1_results   = list(path = file.path(S1, "S1R1_R2_11_extchecks/results.tsv"),  md5 = "b43f1d947b6d3de05df2f4bba5dc80c2"),
  C1_branches  = list(path = file.path(S1, "S1R1_R2_11_extchecks/branches.tsv"), md5 = "d1e063e66d8eca599f464318f694deeb"),
  # comparator tables, as copied by the orchestrator beside each adjudication (byte-identical to _restricted/compare_1b)
  C1_survival  = list(path = file.path(V1B, "C1/adjudication/C1_SURVIVAL.tsv"),   md5 = "3e5999b3c1ac23b8bdbb867983300ec7"),
  S3_survival  = list(path = file.path(V1B, "S3/adjudication/S3_SURVIVAL.tsv"),   md5 = "776cb823394cb54c8e1b68a44590d2e5"),
  S6_survival  = list(path = file.path(V1B, "S6/adjudication/S6_SURVIVAL.tsv"),   md5 = "9b8df7945ca2c24c07ea6c3a6d0fbbfd"),
  # sections without an adjudication copy: read from the comparator folder (the W9F task names it as Stage 1b evidence;
  # only the *_SURVIVAL.tsv tables are read, never the orchestrator's HOLM_CORROB_1b files)
  S1_survival  = list(path = file.path(RES, "_restricted/compare_1b/S1_SURVIVAL.tsv"), md5 = "538a356c5c70b6a499180d642b82f2c5"),
  S2_survival  = list(path = file.path(RES, "_restricted/compare_1b/S2_SURVIVAL.tsv"), md5 = "2a19d25f1ebf83291cfe4a28682ac0cd"),
  S4_survival  = list(path = file.path(RES, "_restricted/compare_1b/S4_SURVIVAL.tsv"), md5 = "5ff0e68cf2c0345f40277ec087d08881"),
  S5_survival  = list(path = file.path(RES, "_restricted/compare_1b/S5_SURVIVAL.tsv"), md5 = "cc5e6f2a9d653094d2c6710d9976f8cc"),
  S7_survival  = list(path = file.path(RES, "_restricted/compare_1b/S7_SURVIVAL.tsv"), md5 = "86782df4d9ee3448fb48f39f840877a3"),
  S3_join      = list(path = file.path(V1B, "S3/adjudication/adj_work/adj_S3_join_audit.tsv"), md5 = "2c28b5e167b01d190ee7727532d5764b"),
  C1_adj       = list(path = file.path(V1B, "C1/adjudication/adj_work/adj_C1.tsv"), md5 = "e2df1246f3eb7595d49c0fb1f9e05edd"),
  C1_tpmi_json = list(path = file.path(V1B, "C1/adjudication/adj_work/adj_C1_tpmi_consistency.json"), md5 = "77b3ec8adc9a95e0cf835af230d5b7db"),
  S3_scen      = list(path = file.path(V1B, "S3/adjudication/adj_work/adj_S3_scenarios.tsv"), md5 = "953b41ad182581de44072164250ef8e6"),
  S3_pkg       = list(path = file.path(V1B, "S3/adjudication/adj_work/adj_S3_coloc_pkg.tsv"), md5 = "64d2b670e88c3a98656d3e6b2200a455"),
  S6_d_table   = list(path = file.path(A6, "a3_d_table.tsv"),        md5 = "f438c8af378a882e55a496772eb3438c"),
  S6_fitset    = list(path = file.path(A6, "a3_allcand_fitset.tsv"), md5 = "d77aa19e1e6b44314538c0ba0e3169f0"),
  S6_notch2    = list(path = file.path(A6, "a3_notch2_fits.tsv"),    md5 = "377dc82580cbc28c8e1519082360aa37"),
  S6_checks    = list(path = file.path(A6, "a3_checks.tsv"),         md5 = "fa03d7e6f7b9a0c019b5d530c81f812c"),
  S6_umod_q    = list(path = file.path(A6, "a3_umod_q.tsv"),         md5 = "d9dbb9b82aabbf9c15f4396f29a28012"),
  # documents quoted by W9F travel rows / the convention table; each quoted fragment is asserted present
  W3_readme    = list(path = file.path(S1, "S1R1_R2_03_extcoloc/README.md"), md5 = "6ded05630e9482b2d42df3d4169c3ad4"),
  W1_readme    = list(path = file.path(S1, "S1R1_R2_01_rad23b/README.md"), md5 = "c0befefe8be54de470ff0db3159da73b"),
  verifier_defects = list(path = file.path(V1B, "VERIFIER_DEFECTS.md"), md5 = "ffbb5c61724815dab6ee7479d6bb2319"),
  posthoc_spec = list(path = file.path(V3, "R2_rewrite/posthoc/P31R_POSTHOC_SPEC.md"), md5 = "573fd781c6d25cc2a912634e79eb2d80"),
  # Stage 2 opening fixes (2026-10-09): the workstreams' branch-input id maps (Stage 1b normalisers) and the Stage 1b report
  idmap_W3     = list(path = file.path(S1, "S1R1_R2_03_extcoloc/branches_inputs_ids.tsv"),  md5 = "080b039e9f1125ccd2e0b65d406f491e"),
  idmap_W4a    = list(path = file.path(S1, "S1R1_R2_04a_arg1grid/branches_inputs_ids.tsv"), md5 = "72a633758112c7c2a0ea74b97bf06543"),
  idmap_W5     = list(path = file.path(S1, "S1R1_R2_05_dpep1/branches_inputs_ids.tsv"),     md5 = "4f883a1a710ad257ecad5ccb74799e63")
)
# Stage 2 opening fix 5 (Rule 7 amendment): the key lists (KEY_MAP.tsv) and not-key lists (NOT_KEYS.tsv) of Stage 1b, read only to
# classify the quantity ids that travel on integrated claim rows (the SURVIVAL comparator tables remain the only verification status read)
KM_DIR <- file.path(RES, "_restricted/keymap")
KM_MD5 <- list(
  C1 = c(key = "f516a57dc1721cfb9dab582e0d34fb49", notkey = "ad557cdde3fb2f6c44021ebff29c1eb7"),
  S1 = c(key = "e9964228392c1bed135cfda557136b90", notkey = "24c6adc8b6a37d776558394d6e47c310"),
  S2 = c(key = "a9ef9764babb68f76f64403b6e884b9c", notkey = "3374cbcaeb2ed739c351b3364a6c5865"),
  S3 = c(key = "24650455ffd471d936879fb11f9735b9", notkey = "570e675867843974db8648ad2b774eb7"),
  S4 = c(key = "9531dd553c96a1ef6e0de6b8d9d3d281", notkey = "ada5f47c880842b5fc993ba3037aa509"),
  S5 = c(key = "02c7ca1a5b9edfa6915add8e9a57cd8e", notkey = "51310312a8cfba56788294938bcb887e"),
  S6 = c(key = "5121b259671ee6748666bc87a7414095", notkey = "45fff103d666bbabd126eb31ce7190bc"),
  S7 = c(key = "701ef2ec95d8c1ad6c552ab881439c9d", notkey = "3abcb1db2326a078db5fcf25b1a2a22b"))
for (s in names(KM_MD5)) {
  W9F_IN[[paste0("keymap_", s, "_key")]]    <- list(path = file.path(KM_DIR, s, "KEY_MAP.tsv"),  md5 = KM_MD5[[s]][["key"]])
  W9F_IN[[paste0("keymap_", s, "_notkey")]] <- list(path = file.path(KM_DIR, s, "NOT_KEYS.tsv"), md5 = KM_MD5[[s]][["notkey"]])
}
STAGE1B_REPORT <- list(path = file.path(V1B, "STAGE1B_REPORT.md"), md5 = "c4f3272dcab644ce3595128141086601")   # read by w9_readme.R (section 1, verbatim)

# ---- Stage 2 opening fixes: constants shared by w9_main.R, w9_readme.R and w9_mutation.R ----
# Fix 2: the comparator's SURVIVES_VISIBLE_TARGET is a lexical token match (H2 D2). Only these 8 key ids are corrected as follows:
# five were genuinely visible (the answer printed in a file the verifiers could open), three are closed-form functions of visible counts.
VIS_GENUINE   <- paste0("S6H.primary.", c("a", "b", "c", "d", "p"))
VIS_DERIVABLE <- paste0("S6H.primary.", c("or", "lo", "hi"))
VIS_KEYS8     <- c(VIS_GENUINE, VIS_DERIVABLE)
VIS_LABEL_GENUINE   <- "SURVIVES_VISIBLE_TARGET (answer printed in a file the verifiers could open)"
VIS_LABEL_DERIVABLE <- "SURVIVES (derivable from visible counts)"
VIS_LABEL_INDEP     <- "SURVIVES (independent re-derivation)"
# Stage 2 step 2 fix round 2: the S6H.d.<cand>.s_max statuses. Max |L1 or L2 - printed seven-fit s_max| over the 15 computable keys, rounded up to 2 significant digits
# (asserted in w9_main.R: S9.check.s6d_smax_L1L2_hold_seven_fit); the text is generated from this constant.
S6D_SMAX_TOL <- 2.1e-8; S6D_SMAX_TOL_TXT <- "2.1e-8"
# Fix 3: corroboration condition 2 where section 3 was not run
COND2_NA <- "not applicable (section 3 not run)"
COR_VERDICTS_BEFORE <- c("PAM-T2D" = "CORROBORATED", "RAD23B-T2D" = "NOT_CORROBORATED", "ARG1-T2D" = "CORROBORATED", "PCSK9-CAD" = "NOT_CORROBORATED",
                         "TGFB1-CAD" = "NOT_CORROBORATED", "FGFR4-T2D (EAS)" = "CORROBORATED", "ARG1-T2D signal B" = "NOT_A_PAIR")   # as in the pre-fix corroboration.tsv
# Fix 5: Rule 7 amendment (user decision 12, 2026-10-09)
RULE7_LABEL <- "analyst-computed; not independently re-derived (middle-path DEVIATION)"
R0_LABEL <- paste0(RULE7_LABEL, "; r0 READ_FROM_SEEN")   # the C1 outcome-check label text (C1.<pair>.label, with its r0)
ADJ_LABEL <- "adjudicator third computation; not a Stage 1b key"   # values computed by an adjudicator (S6 third computation), not in any key list
# the genuine / derivable visible keys, parsed from STAGE1B_REPORT section 3 (the sentence "N had the answer genuinely printed ... M are closed-form functions ...")
report_vis_keys <- function(txt) {
  ln <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  l <- ln[grepl("genuinely printed in a file the verifiers could open", ln, fixed = TRUE) & grepl("closed-form", ln, fixed = TRUE)]
  if (length(l) != 1L) stop("STAGE1B_REPORT: the section 3 sentence on visible keys was not found exactly once")
  cut <- regexpr("closed-form", l, fixed = TRUE)
  ex <- function(s) { tk <- regmatches(s, gregexpr("S6H\\.primary\\.[a-z]+(?:/[a-z]+)*", s, perl = TRUE))[[1]]
    unlist(lapply(sub("^S6H\\.primary\\.", "", tk), function(z) paste0("S6H.primary.", strsplit(z, "/", fixed = TRUE)[[1]]))) }
  list(genuine = ex(substr(l, 1, cut - 1L)), derivable = ex(substr(l, cut, nchar(l))))
}
R7_TOKRE <- "(?:C1|S[1-8])\\.(?:[^\\s;,{}]|\\{[^}]*\\})+"   # a quantity id, a brace group, a wildcard or a <placeholder>
r7_expand <- function(tok, uni) {
  tok <- sub("[:.)]+$", "", tok)
  rx <- gsub(".", "\\.", tok, fixed = TRUE)
  rx <- gsub("<[^>]+>", "[^.]+", rx)
  rx <- gsub("\\{([^}]*)\\}", "(?:\\1)", rx); rx <- gsub(",", "|", rx, fixed = TRUE)
  rx <- gsub("*", ".*", rx, fixed = TRUE)
  uni[grepl(paste0("^", rx, "$"), uni, perl = TRUE)]
}
r7_ids <- function(txt, uni) {
  toks <- regmatches(txt, gregexpr(R7_TOKRE, txt, perl = TRUE))[[1]]
  ids <- character(); unres <- character()
  for (t in toks) { e <- r7_expand(t, uni); if (length(e)) ids <- c(ids, e) else unres <- c(unres, t) }
  list(ids = unique(ids), unresolved = unres)
}
# Fix 1: the deviations section is reproduced verbatim from STAGE1B_REPORT.md section 1; D-1 .. D-4 must all be present
deviations_section <- function(txt) {
  ln <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  a <- grep("^## 1\\. Declared deviations", ln); b <- grep("^## 2\\.", ln)
  if (length(a) != 1L || length(b) != 1L || b <= a) stop("STAGE1B_REPORT section 1 / section 2 headings not found exactly once")
  sec <- ln[(a + 1L):(b - 1L)]
  for (d in paste0("D-", 1:4, ".")) if (sum(grepl(paste0("^- \\*\\*", gsub(".", "\\.", d, fixed = TRUE)), sec)) != 1L) stop("deviation ", d, " not found exactly once in STAGE1B_REPORT section 1")
  paste(sec, collapse = "\n")
}
EXTRA <- list(
  W1_ext_checks = list(path = file.path(S1, "S1R1_R2_01_rad23b/ext_checks.tsv"), md5 = "7f68f93129a8606acbea9ee5043c1f5f"),
  # seen_output, read under the brief's exception 3 (READ_FROM_SEEN); md5 also asserted against RESOURCES.json below
  R1_outcome_checks_tsv = list(path = file.path(S1, "S1R1_05b_outcome_checks/S1R1_outcome_checks.tsv"), md5 = "7d7f14314e04d9609400dbedfabac1f2")
)

md5 <- function(p) unname(tools::md5sum(p))

LOGF <- NULL
logline <- function(...) {
  s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  cat(s, "\n")
  if (!is.null(LOGF)) cat(s, "\n", file = LOGF, append = TRUE)
}

assert_md5 <- function(path, expect, label) {
  if (!file.exists(path)) stop(sprintf("MISSING input %s: %s", label, path))
  got <- md5(path)
  if (!identical(got, expect)) stop(sprintf("MD5 MISMATCH %s: expected %s got %s (%s)", label, expect, got, path))
  logline(sprintf("md5 OK  %-28s %s  %s", label, got, path))
  invisible(got)
}

# Read a TSV as all-character, from either a quoted or an unquoted file, and cross-check the row count with an
# independent reader (fread can silently return an empty table on some paths with spaces; hygiene rule 7).
read_tsv_chk <- function(path, label) {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  quoted <- startsWith(first, "\"")
  dt <- fread(path, sep = "\t", quote = if (quoted) "\"" else "", colClasses = "character",
              na.strings = NULL, encoding = "UTF-8", showProgress = FALSE)
  rd <- utils::read.delim(path, quote = if (quoted) "\"" else "", colClasses = "character",
                          comment.char = "", na.strings = character(0), check.names = FALSE,
                          encoding = "UTF-8", fileEncoding = "UTF-8")
  if (nrow(dt) == 0L) stop(sprintf("EMPTY read of %s", label))
  if (nrow(dt) != nrow(rd)) stop(sprintf("ROW COUNT DISAGREES for %s: fread %d vs read.delim %d", label, nrow(dt), nrow(rd)))
  if (!identical(dt[[1]], rd[[1]])) stop(sprintf("FIRST COLUMN DISAGREES between readers for %s", label))
  logline(sprintf("read    %-28s %d rows x %d cols (quoted=%s; fread == read.delim)", label, nrow(dt), ncol(dt), quoted))
  dt
}

# a quoted fragment must be present verbatim in the pinned document it is quoted from (W9F travel texts, the convention table)
quote_ok <- function(path, frag, label) {
  txt <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  if (!grepl(frag, txt, fixed = TRUE)) stop(sprintf("QUOTE NOT FOUND in %s: '%s'", label, frag))
  logline(sprintf("quote OK  %-18s '%s'", label, substr(frag, 1, 90)))
}

# ---- quantity lookups: every lookup must return exactly one row (join yield 1/1) ----
LOOKUPS <- new.env()
LOOKUPS$rows <- list()
PLACEHOLDER_TOKENS <- c("", "NA", "NaN", "TBD", ".", "-", "NULL", "#N/A", "PENDING", "NOT RUN")
# Stage 2 step 2: allow_inf = TRUE admits +Inf only (the conditional-MLE OR and its exact upper bound of a table with a zero cell,
# S6.d.hi under the seven-fit reading); NaN, -Inf and placeholder tokens still stop.
getq <- function(dt, ws, qid, col = "value", numeric = TRUE, allow_inf = FALSE) {
  hit <- dt[dt$quantity_id == qid]
  if (nrow(hit) != 1L) stop(sprintf("LOOKUP %s %s returned %d rows (expected exactly 1)", ws, qid, nrow(hit)))
  raw <- hit[[col]]
  if (numeric) {
    if (raw %in% PLACEHOLDER_TOKENS) stop(sprintf("LOOKUP %s %s col %s holds placeholder token '%s'", ws, qid, col, raw))
    v <- suppressWarnings(as.numeric(raw))
    if (!is.finite(v) && !(allow_inf && identical(v, Inf))) stop(sprintf("LOOKUP %s %s col %s is not a finite number: '%s'", ws, qid, col, raw))
  } else v <- raw
  LOOKUPS$rows[[length(LOOKUPS$rows) + 1L]] <- data.table(ws = ws, quantity_id = qid, col = col, value = as.character(raw))
  v
}

# ---- Holm, two ways (p.adjust and by hand) ----
holm_manual <- function(p) {
  m <- length(p); o <- order(p); adj <- numeric(m)
  adj[o] <- pmin(1, cummax((m - seq_len(m) + 1) * p[o]))
  adj
}
holm <- function(p) {
  a <- p.adjust(p, method = "holm"); b <- holm_manual(p)
  if (max(abs(a - b)) > 1e-12) stop("Holm: p.adjust and the manual implementation disagree")
  a
}

# Rule 4 "nominal" resolution, two readings (W9F item 6, 2026-10-09: the LITERAL reading is the default column; the
# applied reading I1 is printed beside it for the user).
#   literal (default): "a statement whose Holm p is >= 0.05 carries 'nominal'" -> NOMINAL whenever Holm p >= 0.05,
#                      SURVIVES_HOLM otherwise. A non-rejection (unadjusted p >= 0.05) therefore carries "nominal" too.
#   applied I1       : unadjusted p < 0.05 and Holm p >= 0.05 -> NOMINAL; both < 0.05 -> SURVIVES_HOLM;
#                      unadjusted p >= 0.05 -> NON_REJECTION, worded on the unadjusted p (Rule 4, third sentence), no "nominal".
nominal_literal <- function(p_unadj, p_holm) {
  if (p_holm >= 0.05) return(list(code = if (p_unadj >= 0.05) "NOMINAL (non-rejection)" else "NOMINAL", attach = TRUE))
  list(code = "SURVIVES_HOLM", attach = FALSE)
}
nominal_I1 <- function(p_unadj, p_holm) {
  if (p_unadj < 0.05 && p_holm >= 0.05) return(list(code = "NOMINAL", attach = TRUE))
  if (p_unadj < 0.05 && p_holm < 0.05)  return(list(code = "SURVIVES_HOLM", attach = FALSE))
  list(code = "NON_REJECTION", attach = FALSE)
}

# ---- W9F check functions. Each returns list(ok, detail); w9_main.R records it with add_assert(), and w9_mutation.R
#      feeds the SAME function a mutated copy of its real input and requires ok == FALSE. ----

# W9F item 1: a C1 class is "verified" when the comparator row for C1H.<pair>.class has status SURVIVES and the
# implementer, L1 and L2 levels all equal the class C1 printed in its results.tsv.
chk_c1_class_verified <- function(surv, pair, c1_class) {
  r <- surv[surv$key_id == paste0("C1H.", pair, ".class")]
  if (nrow(r) != 1L) return(list(ok = FALSE, detail = sprintf("C1H.%s.class: %d comparator rows", pair, nrow(r))))
  ok <- identical(r$status, "SURVIVES") && identical(r$impl, c1_class) && identical(r$L1, c1_class) && identical(r$L2, c1_class) &&
        identical(r$L1_ok, "True") && identical(r$L2_ok, "True")
  list(ok = ok, detail = sprintf("C1H.%s.class: status %s; impl %s, L1 %s, L2 %s; C1 results.tsv class %s", pair, r$status, r$impl, r$L1, r$L2, c1_class))
}

# W9F item 2b: pre-declaration section 1 RAD23B clause: add "the classification depends on the p source and on rounding
# of the printed DIAMANTE estimates" if any printed p source is >= 0.05, or any rounding interval straddles 0.05.
rad23b_clause_fires <- function(p, lo, hi) {
  stopifnot(length(lo) == length(hi), all(is.finite(c(p, lo, hi))), all(lo <= hi))
  any(p >= 0.05) || any(lo < 0.05 & hi >= 0.05)
}
parse_interval <- function(s) { x <- as.numeric(strsplit(s, "-", fixed = TRUE)[[1]]); if (length(x) != 2L || any(!is.finite(x))) stop("bad interval '", s, "'"); x }
chk_rad23b_clause <- function(adj_row, w1_branch) {
  iv <- parse_interval(adj_row$p_wald1_rounding_interval)
  p <- as.numeric(c(adj_row$p_resource, adj_row$p_wald1, adj_row$p_wald2))
  fires <- rad23b_clause_fires(p, iv[1], iv[2])
  w1_added <- grepl("p-source/rounding clause ADDED", w1_branch, fixed = TRUE)
  list(ok = isTRUE(fires) && w1_added, fires = fires,
       detail = sprintf("adjudicator (adj_C1.tsv): resource p %s, first-order rounded-Wald p %s (rounding interval %s-%s), second-order p %s -> clause %s; W1 branch: clause %s",
                        format(p[1], digits = 4), format(p[2], digits = 4), format(iv[1]), format(iv[2]), format(p[3], digits = 4),
                        if (fires) "FIRES" else "does not fire", if (w1_added) "ADDED" else "not added"))
}

# W9F item 4: the Amendment 1 A2 caveat names the resource(s) failing Rule 3(ii). Robustness is read from W4a's Rule 3(ii)
# rows and cross-checked against the bootstrap share and the 489-panel side of 0.05 in W4a's results.tsv.
chk_a2_failing <- function(w4a_robust, share, p503, p489, d503, d489) {
  # w4a_robust, share, p503, p489, d503, d489: named by resource (FinnGen, DIAMANTE)
  mine <- share >= 0.90 & ((p503 < 0.05) == (p489 < 0.05)) & (sign(d503) == sign(d489))
  fails <- names(mine)[!mine]
  list(ok = identical(unname(mine), unname(w4a_robust)), fails = fails,
       detail = sprintf("%s; W4a Rule 3(ii) rows: %s -> resources failing Rule 3(ii): %s",
                        paste(sprintf("%s share %.4f, 503 p %.4g / 489 p %.4g, A-B sign %+d/%+d -> robust %s", names(mine), share, p503, p489,
                                      as.integer(sign(d503)), as.integer(sign(d489)), mine), collapse = "; "),
                        paste(sprintf("%s %s", names(w4a_robust), ifelse(w4a_robust, "robust", "NOT robust")), collapse = ", "),
                        if (length(fails)) paste(fails, collapse = " and ") else "none"))
}

# W9F item 2a: the section 6(d) table from per-candidate s_max values (Rule 3(i): s_max >= 0.10, or not computable,
# reclassifies the candidate; the 'kept' version keeps a not-computable candidate). Primary table 11/16 vs 57/144.
S6D_KIDNEY <- c("CELSR2_eGFR", "DPEP1_eGFR", "GSTA1_eGFR", "GUSB_eGFR", "INHBB_eGFR", "ITIH4_eGFR", "NT5C3A_eGFR", "PBLD_eGFR",
                "SORT1_eGFR", "UMOD_eGFR", "UMOD_CKD")
S6D_NONKID <- c("APOBR_HF", "APOE_CAD", "ING1_CAD", "NOTCH2_T2D", "NPDC1_T2D")
s6d_table <- function(smax, keep_not_computable = FALSE, a0 = 11L, b0 = 5L, c0 = 57L, d0 = 87L) {
  stopifnot(setequal(names(smax), c(S6D_KIDNEY, S6D_NONKID)))
  recl <- names(smax)[is.na(smax) | smax >= 0.10]
  if (keep_not_computable) recl <- setdiff(recl, names(smax)[is.na(smax)])
  list(n = length(recl), a = a0 - sum(recl %in% S6D_KIDNEY), b = b0 - sum(recl %in% S6D_NONKID),
       c = c0 + sum(recl %in% S6D_KIDNEY), d = d0 + sum(recl %in% S6D_NONKID), reclassified = sort(recl))
}

# ---- W9F gates over branches_integrated (IR) and corroboration (COR). Each returns list(ok, detail) and requires at
#      least one target row, so a selector that matches nothing fails rather than passing vacuously.
has_all <- function(x, needles) all(vapply(needles, function(n) all(grepl(n, x, fixed = TRUE)), TRUE))
gate <- function(sel, needles_travel = character(0), needles_q = character(0), extra = TRUE, n_expect = NA, label = "") {
  n <- sum(sel)
  ok <- n >= 1L && (is.na(n_expect) || n == n_expect) && isTRUE(extra)
  list(ok = ok, n = n, detail = sprintf("%s: %d target row(s)%s", label, n, if (!is.na(n_expect)) sprintf(" (expected %d)", n_expect) else ""))
}
# W9F step 2b (D4): the six-signal p values printed on the two five-signal omnibus rows are checked BY VALUE against W4b
# results.tsv (as printed: 2 significant figures). Each label must occur exactly once in travel_w9f and once in
# qualifiers_resolved, and each printed number must equal signif(W4b value, 2).
OMNI_P_KEYS <- list(
  arginine = list(row = "4 Omnibus (arginine)", label = "arginine p ",                  qid = "S4.omni.arginine.6sig.p_aware_2nd"),
  T2D_FG   = list(row = "4 Omnibus (T2D)",      label = "T2D FinnGen R10 p ",           qid = "S4.omni.T2D_FG.6sigproxy.p_aware_2nd"),
  T2D_FG5  = list(row = "4 Omnibus (T2D)",      label = "the five available signals p ", qid = "S4.omni.T2D_FG.6sigdrop5.p_aware_2nd"),
  T2D_DM   = list(row = "4 Omnibus (T2D)",      label = "T2D DIAMANTE p ",              qid = "S4.omni.T2D_DM.6sig.p_aware_2nd"))
omni_printed_p_check <- function(IR, w4b) {
  if (is.null(w4b)) return(list(ok = FALSE, detail = "W4b results.tsv not supplied to the gate"))
  one <- function(k) {
    sel <- IR$source == "W4b" & IR$rule_ref == k$row
    hit <- w4b$value[w4b$quantity_id == k$qid]
    if (sum(sel) != 1L || length(hit) != 1L) return(list(ok = FALSE, txt = sprintf("%s: %d rows / %d W4b values", trimws(k$label), sum(sel), length(hit))))
    v <- as.numeric(hit)
    pat <- paste0(k$label, "([0-9.]+(e-?[0-9]+)?)")
    got <- lapply(c(travel_w9f = IR$travel_w9f[sel], qualifiers_resolved = IR$qualifiers_resolved[sel]),
                  function(x) as.numeric(sub(pat, "\\1", regmatches(x, gregexpr(pat, x))[[1]])))
    ok <- is.finite(v) && all(lengths(got) == 1L) && all(vapply(got, function(g) length(g) == 1L && is.finite(g) && abs(g - signif(v, 2)) <= 1e-9 * abs(v), TRUE))
    list(ok = ok, txt = sprintf("%s printed %s / %s vs W4b %s %s", trimws(k$label), paste(got$travel_w9f, collapse = ","), paste(got$qualifiers_resolved, collapse = ","),
                                k$qid, format(v, digits = 4)))
  }
  r <- lapply(OMNI_P_KEYS, one)
  list(ok = all(vapply(r, `[[`, TRUE, "ok")), detail = paste(vapply(r, `[[`, "", "txt"), collapse = "; "))
}
W9F_GATES <- list(
  arg1_framing = function(IR, ...) {
    s <- IR$source %in% c("W4a", "W4b") & IR$row_kind == "claim"
    g <- gate(s, label = "W4a/W4b claim rows carry the ARG1 framing disclosure in travel_w9f and qualifiers_resolved",
              extra = has_all(IR$travel_w9f[s], c("ARG1 framing disclosure", "returned WEAKENS", "kidney-led framing", "WEAKENS is frozen")) &&
                      has_all(IR$qualifiers_resolved[s], "ARG1 framing disclosure"))
    g },
  omnibus_disclosure = function(IR, w4b = NULL, ...) {
    s <- IR$source == "W4b" & IR$rule_ref %in% c("4 Omnibus (arginine)", "4 Omnibus (T2D)")
    t2d <- IR$source == "W4b" & IR$rule_ref == "4 Omnibus (T2D)"
    vchk <- omni_printed_p_check(IR, w4b)
    g <- gate(s, n_expect = 2L, label = "five-signal omnibus rows carry the set-choice disclosure with the six-signal p beside (values = W4b results.tsv), and it is their A1 front caveat",
         extra = has_all(IR$qualifiers_resolved[s], c("the five-signal set was chosen after seeing the five-signal Q values", "six-signal sensitivity beside")) &&
                 has_all(IR$qualifiers_resolved[t2d], c("T2D FinnGen R10 p", "T2D DIAMANTE p")) &&
                 all(IR$a1_front_caveat[s] == A1_LABEL[["OMNI_SETCHOICE"]]) && vchk$ok)
    g$detail <- paste0(g$detail, "; ", vchk$detail)
    g },
  a2_merged = function(IR, ...) {
    s <- IR$source == "W4a" & IR$rule_ref == "4 T2D"
    f <- IR$source == "W4a" & IR$rule_ref == "Amendment 1 A2"
    gate(s, n_expect = 1L, label = "the section 4 T2D claim carries the A2 fragment and the resource-named LD caveat; the A2 row is marked merged",
         extra = has_all(IR$filled_template_integrated[s], "in both outcome resources (FinnGen p") &&
                 has_all(IR$qualifiers_resolved[s], c("sensitive to the reference LD in FinnGen", "never 'replicated'/'replication'")) &&
                 !grepl("sensitive to the reference LD(?! in)", IR$qualifiers_resolved[s], perl = TRUE) &&
                 sum(f) == 1L && startsWith(IR$row_kind[f], "fragment merged")) },
  eas_derived_route = function(IR, ...) {
    s <- IR$source == "W7" & grepl("^EAS_eGFR", IR$branch)
    gate(s, n_expect = 1L, label = "East Asian eGFR scale claim: 'derived route; not new evidence' travels and is the A1 front caveat under both readings",
         extra = has_all(IR$qualifiers_resolved[s], "derived route; not new evidence") &&
                 all(IR$a1_front_caveat[s] == A1_LABEL[["DERIVED_ROUTE"]]) && all(IR$a1_front_caveat_I1[s] == A1_LABEL[["DERIVED_ROUTE"]])) },
  s8_front_note = function(IR, ...) {
    s <- IR$source == "W7" & grepl("^s8 ", IR$rule_ref)
    gate(s, n_expect = 1L, label = "section 8 row: front_surface_note for its declared caveat ('post hoc') and A1 front caveat",
         extra = has_all(IR$front_surface_note[s], "post hoc (P31R_POSTHOC_SPEC.md 'Status: POST HOC')") && all(IR$a1_front_caveat[s] == A1_LABEL[["S8_POSTHOC"]])) },
  rad23b_s3_disclosure = function(IR, ...) {
    s <- IR$source == "W3" & grepl("^RAD23B_T2D (primary input|pmafdrop)", IR$inputs)
    gate(s, n_expect = 2L, label = "RAD23B section 3 PP.H3 rows carry the drop-version disclosure",
         extra = has_all(IR$travel_w9f[s], c("were seen in the first development run before MAF substitution was chosen as primary", "both versions are printed"))) },
  s6d_travel = function(IR, ...) {
    s <- IR$source == "W6b" & IR$rule_ref %in% c("6(d)", "11 D6", "6 renal OR wording; 11 D6")
    cand <- IR$source == "W6b" & IR$rule_ref == "6(d); Rule 3(i)"
    n2 <- cand & grepl("^NOTCH2-T2D", IR$filled_template)
    gate(s, n_expect = 3L, label = "section 6(d) rows carry the decided fit set (seven fits primary) with the four-fit sensitivity; every candidate row carries its seven-fit and four-fit s_max",
         extra = has_all(IR$travel_w9f[s], c("USER DECISION 2026-10-09: seven fits primary", "four fits: n_reclassified", "seven fits: n_reclassified", "NOTCH2-T2D borderline", "no section 6 branch level changes")) &&
                 sum(cand) == 16L && has_all(IR$travel_w9f[cand], c("s_max seven fits", "four fits (sensitivity")) &&
                 sum(n2) == 1L && has_all(IR$travel_w9f[n2], "NOTCH2-T2D borderline")) },
  nominal_literal_default = function(IR, ...) {
    s <- nzchar(IR$holm_member)
    lit <- startsWith(IR$nominal[s], "NOMINAL")
    gate(s, label = "rows carrying a Rule 4 member: 'nominal' column is the literal reading; the qualifier text carries it, with I1 beside",
         extra = all(nzchar(IR$nominal_I1[s])) && all(grepl("nominal [W9", IR$qualifiers_resolved[s][lit], fixed = TRUE)) &&
                 all(grepl("no 'nominal' [W9", IR$qualifiers_resolved[s][!lit], fixed = TRUE)) &&
                 all(grepl("I1", IR$qualifiers_resolved[s], fixed = TRUE))) }
)
COR_GATES <- list(
  c1_condition1 = function(COR) {
    s <- COR$pkey != "ARG1_T2D_signalB"
    gate(s, n_expect = 6L, label = "condition 1 uses C1's verified class for all six pairs, with the R1 class beside",
         extra = all(COR$check_class_used[s] == COR$check_class_C1[s]) && all(COR$C1_class_status[s] == "SURVIVES") &&
                 all(nzchar(COR$check_class_R1[s])) && all(COR$cond1_concordant[s] == (COR$check_class_C1[s] == "concordant"))) },
  rad23b_clause = function(COR) {
    s <- COR$pkey == "RAD23B_T2D"
    gate(s, n_expect = 1L, label = "the RAD23B class carries the section 1 p-source/rounding clause and its 0.05/6 result",
         extra = has_all(COR$travel_w9f[s], c("the classification depends on the p source and on rounding of the printed DIAMANTE estimates", "0.05/6"))) },
  # W9F step 2b (D3): W9F item 3 ('a travel row on every ARG1 slot claim') reaches corroboration.tsv's two ARG1 rows
  arg1_framing_corrob = function(COR) {
    s <- COR$pkey %in% c("ARG1_T2D", "ARG1_T2D_signalB")
    gate(s, n_expect = 2L, label = "the ARG1 corroboration rows (ARG1-T2D; ARG1-T2D signal B) carry the ARG1 framing disclosure in travel_w9f, and no other row does",
         extra = has_all(COR$travel_w9f[s], c("ARG1 framing disclosure", "returned WEAKENS", "kidney-led framing", "WEAKENS is frozen")) &&
                 !any(grepl("ARG1 framing disclosure", COR$travel_w9f[!s], fixed = TRUE))) }
)

# Renal D6 base phrase and the 'sensitive to {list}' qualifier (pre-declaration section 11 D6 and section 6 wording rule).
# One implementation, used for W9's original W6b cross-check and for every section 6(d) reading.
incl1 <- function(lo, hi) lo <= 1 && hi >= 1
renal_d6_list <- function(primary_lo, d_lo, d_hi, lo_lo, lo_hi, c_lo, c_hi, a_lo, a_hi) {
  d6 <- if (round(primary_lo, 3) > 1.000 && round(d_lo, 3) > 1.000) "over-represented at the pair level" else "proportionally more frequent"
  lst <- c(if (incl1(lo_lo, lo_hi) || incl1(c_lo, c_hi)) "grouping", if (incl1(a_lo, a_hi)) "the acceptance rule", if (incl1(d_lo, d_hi)) "the LD guard")
  list(d6 = d6, list = lst)
}

# ---- section 3 wording branch from PP.H3/PP.H4 (Rule 1 colocalization rule: PP.H4 >= 0.80 and PP.H4 >= 5 x PP.H3) ----
coloc_rule_met <- function(h4, h3) (h4 >= 0.80) && (h4 >= 5 * h3)
# h4only_at_1e6 = TRUE: the PP.H4-only alternative the S3 adjudicator also ran (column branch_alt_h4only): PP.H4 >= 0.80
# alone in place of the Rule 1 test at 1e-6; 1e-5 keeps Rule 1.
s3_branch <- function(h4_5, h3_5, h4_6, h3_6, h4only_at_1e6 = FALSE) {
  m5 <- coloc_rule_met(h4_5, h3_5); m6 <- if (h4only_at_1e6) h4_6 >= 0.80 else coloc_rule_met(h4_6, h3_6)
  if (m5 && m6) return("1")
  if (m5 && !m6) return(if (h3_6 >= 0.5) "2 (with add-on)" else "2 (no add-on)")
  if (!m5 && (h3_5 >= 0.5 || h3_6 >= 0.5)) return("3")
  "4"
}
# W9 labels <-> the S3 adjudicator files' branch codes ("2" / "2_addon")
S3_CODE <- c("2 (no add-on)" = "2", "2 (with add-on)" = "2_addon")
s3_file_code <- function(b) unname(ifelse(b %in% names(S3_CODE), S3_CODE[b], b))
s3_label <- function(code) unname(ifelse(code %in% S3_CODE, names(S3_CODE)[match(code, S3_CODE)], code))
# section 3 wording-rule phrases (redacted pre-declaration section 3; each fragment is asserted present by w9_main.R)
S3_PHRASE <- c("1" = "supported in external outcome data, robust to the prior", "2" = "supported at the default prior",
               "2_addon" = "supported at the default prior + 'at the conservative prior the data favour distinct variants'",
               "3" = "the external data favour distinct variants", "4" = "the external data do not resolve a shared variant")
S3_PHRASE_FRAGS <- c("supported in external outcome data, robust to the prior", "supported at the default prior",
                     "at the conservative prior the data favour distinct variants", "the external data favour distinct variants",
                     "the external data do not resolve a shared variant")
# W9F step 2b (D1): the S3 adjudicator's 'full window' definitions, in BOTH of its files (adj_S3_coloc_pkg.tsv readings,
# adj_S3_scenarios.tsv scenarios). Labels paraphrase the adjudicator's own definitions (adj_s3_coloc_pkg.R header;
# adj_s3_coloc.py); an unlabelled non-declared scenario stops the run, so a definition can never be dropped silently.
S3_DECLARED_SCEN <- c("A_analysed", "declared")
S3_DEF_LABEL <- c(
  B_full_exposure         = "coloc.abf given the whole exposure window as dataset 1 (it caps p1)",
  C_full_both             = "coloc.abf given the whole exposure window and the whole valid outcome window, rows coloc.abf accepts: 0 < MAF < 1, unique SNP (it caps p1 and p2)",
  cap_dataset1_fullwindow = "p1 capped with the exposure-window row count",
  cap_both_fullwindows    = "p1 and p2 capped with the exposure-window and the raw outcome-window row counts")
cq03_full_windows <- function(s3adj, pr = "PAM_T2D", mm = "sub") {   # argument names must not be column names (data.table scoping)
  fw <- s3adj[s3adj$pair == pr & s3adj$maf_mode == mm & !(s3adj$scenario %in% S3_DECLARED_SCEN)]
  if (nrow(fw) == 0L || any(fw$pair != pr) || any(fw$maf_mode != mm)) stop("cq03_full_windows: selection failed")
  if (!all(fw$scenario %in% names(S3_DEF_LABEL))) stop("unlabelled S3 adjudicator scenario(s): ", paste(setdiff(fw$scenario, names(S3_DEF_LABEL)), collapse = ", "))
  fw
}
s3_def_name <- function(fw) paste(fw$file, fw$scenario)
s3_def_segment <- function(r) sprintf("%s %s (%s; %s; p1 %s, p2 %s): PP.H4/PP.H3 1e-5 %.3f/%.3f, 1e-6 %.3f/%.3f (ratio %.2f); Rule 1 branch %s (PP.H4-only reading %s)",
  r$file, r$scenario, S3_DEF_LABEL[[r$scenario]], r$n_desc, formatC(r$p1, digits = 3, format = "g"), formatC(r$p2, digits = 3, format = "g"),
  r$H4_5, r$H3_5, r$H4_6, r$H3_6, r$H4_6 / r$H3_6, r$branch, r$branch_alt_h4only)
# CQ03 must report every full-window definition (its own segment with its PP.H4/PP.H3 and Rule 1 branch), and its
# effect text must name, per definition, whether PAM's section 3 sentence changes and what happens to condition 2.
chk_cq03 <- function(cq, s3adj, pam_b0) {
  fw <- cq03_full_windows(s3adj)
  if (nrow(cq) != 1L) return(list(ok = FALSE, detail = sprintf("CQ03: %d rows", nrow(cq))))
  files_ok <- all(c("adj_S3_coloc_pkg.tsv", "adj_S3_scenarios.tsv") %in% fw$file)
  segs <- strsplit(cq$consequence_B, " || ", fixed = TRUE)[[1]]
  dsegs <- segs[startsWith(segs, "adj_S3_")]
  seg_ok <- vapply(seq_len(nrow(fw)), function(i) {
    s <- dsegs[startsWith(dsegs, paste0(s3_def_name(fw[i]), " ("))]
    length(s) == 1L && grepl(sprintf("1e-5 %.3f/%.3f, 1e-6 %.3f/%.3f", fw$H4_5[i], fw$H3_5[i], fw$H4_6[i], fw$H3_6[i]), s, fixed = TRUE) &&
      grepl(sprintf("Rule 1 branch %s (PP.H4-only reading %s)", fw$branch[i], fw$branch_alt_h4only[i]), s, fixed = TRUE)
  }, TRUE)
  eff <- cq$effect_on_branches
  chg <- fw$branch != pam_b0
  clause <- function(lead) { m <- regmatches(eff, regexpr(paste0(lead, "[^;]*"), eff)); if (length(m)) m else "" }
  eff_ok <- if (any(chg)) !startsWith(eff, "none") &&
      all(vapply(which(chg), function(i) grepl(s3_def_name(fw[i]), clause(sprintf("changes to branch %s \\(", fw$branch[i])), fixed = TRUE), TRUE)) else startsWith(eff, "none")
  eff_ok <- eff_ok && all(vapply(which(!chg), function(i) grepl(s3_def_name(fw[i]), clause("unchanged under "), fixed = TRUE), TRUE))
  c2 <- vapply(fw$branch, function(b) corroborated(TRUE, s3_label(b))$c2, TRUE); c2_0 <- corroborated(TRUE, s3_label(pam_b0))$c2
  c2_ok <- if (all(c2 == c2_0)) grepl("condition 2 is unchanged under every definition", eff, fixed = TRUE) else grepl("condition 2 changes under", eff, fixed = TRUE)
  list(ok = files_ok && length(dsegs) == nrow(fw) && all(seg_ok) && eff_ok && c2_ok,
       detail = sprintf("PAM full-window definitions in the two S3 adjudicator files: %d (%s); CQ03 segments %d, each with its PP.H4/PP.H3 and Rule 1 branch: %s; Rule 1 branch differs from W3 primary %s for %d (%s): effect text names it %s; condition 2 %s",
                        nrow(fw), paste(s3_def_name(fw), collapse = ", "), length(dsegs), all(seg_ok), pam_b0, sum(chg),
                        if (any(chg)) paste(sprintf("%s -> %s", s3_def_name(fw[chg]), fw$branch[chg]), collapse = ", ") else "none", eff_ok,
                        if (all(c2 == c2_0)) "unchanged under every definition, stated" else "changes, stated"))
}

# ---- section 2 wording branch ----
s2_branch <- function(rA_fg, rB_fg, p_fg, rA_ieu, rB_ieu, p_ieu, disc_sign) {
  opp <- function(a, b) sign(a) != sign(b)
  if (opp(rA_fg, rB_fg) && opp(rA_ieu, rB_ieu) && p_fg < 0.05 && p_ieu < 0.05) return("(a)")
  if (all(sign(c(rA_fg, rB_fg, rA_ieu, rB_ieu)) == disc_sign) && p_fg >= 0.05 && p_ieu >= 0.05) return("(b)")
  "(c)"
}

# ---- Principle 0 "externally corroborated" ----
corroborated <- function(check_concordant, s3, s2 = NA_character_, is_tgfb1 = FALSE) {
  c1 <- isTRUE(check_concordant)
  c2 <- is.na(s3) || s3 %in% c("1", "2 (no add-on)")
  c3 <- !is_tgfb1 || identical(s2, "(b)")
  list(value = c1 && c2 && c3, c1 = c1, c2 = c2, c3 = c3)
}

# ---- Amendment 1 A1: caveat detection and the fixed order ----
# W9F (2026-10-09) adds three declared caveats that are not in Principle 0's named list; A1 step 4 places them after the
# list's items (they never co-occur on one claim): the section 4(iv) five-signal set-choice disclosure (W9F item 2),
# the section 8 'post hoc' status (W9F item 7) and the East Asian eGFR "derived route; not new evidence" (W9F item 5).
A1_ORDER <- c("COLOC", "LDREF", "NOMINAL",
              "S3ADDON", "S2B", "SCALE", "FIRSTMENTION", "SIGNR", "RAD23B_PSOURCE",
              "EPITOPE", "RENAL_LIST", "D12_POSTHOC", "UNRESOLVED", "NOT_MEDIATION",
              "OMNI_SETCHOICE", "S8_POSTHOC", "DERIVED_ROUTE")
A1_LABEL <- c(
  COLOC = "colocalization clause (Rule 3(iv))",
  LDREF = "sensitive to the reference LD (Rule 3(ii))",
  NOMINAL = "nominal (Rule 4)",
  S3ADDON = "section 3 conservative-prior add-on",
  S2B = "the section 2(b) caveat",
  SCALE = "section 7 'scale not established'",
  FIRSTMENTION = "first-mention flag (Rules 2 and 3(i), incl. the 'proxy' label)",
  SIGNR = "the sign of the correlation between the signals is not identified",
  RAD23B_PSOURCE = "RAD23B p-source/rounding clause and its 0.05/6 result",
  EPITOPE = "rs12983775 epitope clause",
  RENAL_LIST = "renal 'sensitive to {list}'",
  D12_POSTHOC = "D12 'post hoc'",
  UNRESOLVED = "'unresolved' gene attribution",
  NOT_MEDIATION = "'not a formal mediation test'",
  OMNI_SETCHOICE = "post hoc set choice: the five-signal set was chosen after seeing the five-signal Q values (section 4(iv))",
  S8_POSTHOC = "post hoc (section 8 numbers fixed by the post hoc specification)",
  DERIVED_ROUTE = "derived route; not new evidence (East Asian eGFR scale)")

detect_caveats <- function(rule_ref, branch, filled, quals, nominal_attach) {
  low <- function(x) tolower(ifelse(is.na(x), "", x))
  tpl <- low(filled); q <- low(quals); br <- low(branch); rr <- low(rule_ref)
  # W9F: the corroboration statement W9 writes into W3's qualifiers is a separate claim (W9 branches.tsv); its own
  # travel items (e.g. the RAD23B p-source clause of the external-check class) are not caveats of the section 3 claim
  q <- sub("corroboration \\(w9, principle 0\\):.*?(?= \\| w9f travel:|$)", "", q, perl = TRUE)
  all_txt <- paste(br, tpl, q)
  hits <- character(0)
  if (grepl("did not colocalize with the protein", all_txt, fixed = TRUE) && !grepl("^does not fire", br)) hits <- c(hits, "COLOC")
  if (grepl("sensitive to the reference ld", all_txt, fixed = TRUE)) hits <- c(hits, "LDREF")
  if (isTRUE(nominal_attach)) hits <- c(hits, "NOMINAL")
  if (grepl("at the conservative prior the data favour distinct variants", tpl, fixed = TRUE)) hits <- c(hits, "S3ADDON")
  if (grepl("did not persist after mutual conditioning", tpl, fixed = TRUE)) hits <- c(hits, "S2B")
  if (grepl("scale not established", paste(br, tpl), fixed = TRUE)) hits <- c(hits, "SCALE")
  fm <- firstmention_flags(all_txt)
  if (length(fm)) hits <- c(hits, "FIRSTMENTION")
  if (grepl("sign of the correlation between the signals is not identified", all_txt, fixed = TRUE)) hits <- c(hits, "SIGNR")
  if (grepl("depends on the p source", all_txt, fixed = TRUE)) hits <- c(hits, "RAD23B_PSOURCE")
  if (grepl("epitope or assay effect is not excluded", all_txt, fixed = TRUE)) hits <- c(hits, "EPITOPE")
  if (grepl("sensitive to (grouping|the acceptance rule|the ld guard)", tpl) ) hits <- c(hits, "RENAL_LIST")
  if (grepl("d12", rr, fixed = TRUE) && grepl("post hoc", paste(tpl, q), fixed = TRUE)) hits <- c(hits, "D12_POSTHOC")
  if (grepl("unresolved", tpl, fixed = TRUE)) hits <- c(hits, "UNRESOLVED")
  if (grepl("not a formal mediation test", tpl, fixed = TRUE)) hits <- c(hits, "NOT_MEDIATION")
  if (grepl("the five-signal set was chosen after seeing the five-signal q values", all_txt, fixed = TRUE)) hits <- c(hits, "OMNI_SETCHOICE")
  if (grepl("^s8 ", rr) && grepl("post hoc", q, fixed = TRUE)) hits <- c(hits, "S8_POSTHOC")
  if (grepl("derived route; not new evidence", q, fixed = TRUE)) hits <- c(hits, "DERIVED_ROUTE")
  hits <- A1_ORDER[A1_ORDER %in% hits]
  oc <- character(0)
  if (grepl("forward-strand assumption", all_txt, fixed = TRUE)) oc <- c(oc, "forward-strand assumption")
  if (grepl("ukb-free external consortium check", all_txt, fixed = TRUE)) oc <- c(oc, "UKB-free external consortium check, not independent replication")
  if (grepl("independent east asian check", all_txt, fixed = TRUE)) oc <- c(oc, "independent East Asian check")
  list(caveats = hits, front = if (length(hits)) hits[1] else NA_character_, oclabels = oc, fm = fm)
}

# The Principle 0 item "the first-mention flags of Rules 2 and 3(i), including the 'proxy' label" bundles three flags.
# Rule 2 counts only where the pair's class actually differs under the sensitivity (W1 prints both classes).
firstmention_flags <- function(t) {
  t <- tolower(t); out <- character(0)
  r2 <- grepl("class differs between first- and second-order", t, fixed = TRUE) || grepl("rule 2: class [a-z]+ under second-order", t) ||
        grepl("class would be", t, fixed = TRUE)
  m <- regmatches(t, regexec("class under second-order se: ([a-z]+) \\(first-order ([a-z]+)\\)", t))[[1]]
  if (length(m) == 3 && m[2] != m[3]) r2 <- TRUE
  if (r2) out <- c(out, "Rule 2 class flag")
  if (grepl("flagged at first mention|first-mention flag applies|rule 3\\(i\\) first-mention flag", t)) out <- c(out, "Rule 3(i) flag")
  if (grepl("'proxy' label", t, fixed = TRUE)) out <- c(out, "'proxy' label")
  out
}

front_surface_note <- function(filled, quals, rule_ref = "") {
  t <- paste(ifelse(is.na(filled), "", filled), ifelse(is.na(quals), "", quals))
  out <- character(0)
  if (grepl("never in the Abstract", t, fixed = TRUE)) out <- c(out, "secondary (Rule 4): never in the Abstract")
  if (grepl("reported in the supplement", t, fixed = TRUE) || grepl("supplement only", t, fixed = TRUE)) out <- c(out, "supplement statement")
  if (grepl("DEVIATION", t, fixed = TRUE)) out <- c(out, "carries a DEVIATION label (not in Principle 0's list; Stage 2 places it)")
  # W9F item 7: the section 8 row's declared caveat. 'Status: POST HOC' is declared by the post hoc specification the
  # pre-declaration's section 8 names (P31R_POSTHOC_SPEC.md, md5 573fd781); its 'not independent replication' sentence is
  # the outcome-check label, which travels in addition (Amendment 1 A1, last bullet).
  if (grepl("^s8 ", rule_ref) && grepl("post hoc", t, fixed = TRUE))
    out <- c(out, "post hoc (P31R_POSTHOC_SPEC.md 'Status: POST HOC'): on a front surface the claim carries 'post hoc' as its one A1 caveat; 'not independent replication' (FinnGen overlap; outcome-check label) travels in addition")
  paste(out, collapse = "; ")
}

# Rows that evaluate a rule or carry a fragment/qualifier for another row's sentence (not a claim on their own).
# Classified by rule_ref and template shape; every classification is printed in branches_integrated.tsv.
is_rule_row <- function(ws, rule_ref, branch, filled) {
  rr <- ifelse(is.na(rule_ref), "", rule_ref); tp <- ifelse(is.na(filled), "", trimws(filled)); br <- ifelse(is.na(branch), "", branch)
  if (tp == "") return("rule-evaluation row (no sentence)")
  if (grepl("^(Rule 3\\(ii\\)|Rule 3\\(iv\\)|Rule 4 family entry|Rule 3\\(i\\) scope)", rr)) return("rule-evaluation row (its qualifier travels to the claim row)")
  if (grepl("^Amendment 1 A2$", rr)) return("fragment of the section 4 T2D sentence (A2)")
  if (grepl("descriptive only", rr) && grepl("not a wording input", tp)) return("descriptive variant (not a wording input)")
  if (grepl("^(add |no clause|Holm adjustment|computed by W6b|no LD qualifier|any sentence about|s not computed)", tp)) return("instruction row (qualifier travels to the claim row)")
  if (ws == "W6a") return("partial row superseded by W6b's completed renal rows")
  if (ws == "W6b" && grepl("^11 D6$", rr)) return("fragment (D6 base phrase) of the renal OR sentence")
  if (ws == "W6b" && grepl("^6\\(d\\); Rule 3\\(i\\)$", rr)) return("per-candidate input to section 6(d)")
  ""
}
