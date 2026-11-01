#!/usr/bin/env Rscript
## C1 -- COMPARE_SEEN (Stage 1 brief section 0 rule 10: compute first, compare second).
## Runs only after C1_extchecks.R has written results.tsv. Never modifies results.tsv or branches.tsv (asserted by md5).
## Seen sources opened here, and only here:
##   R1_outcome_checks_tsv (RESOURCES role seen_output)          -- the R1 outcome-check table
##   full pre-declaration P31R_R2_STAGE1_PRESPEC.md (3319dc39...)  -- section 1 RAD23B DIAMANTE, section 2 TGFB1 Seen lines
##   full amendment P31R_R2_STAGE1_PRESPEC_AMEND1.md (37750883...) -- searched; holds no Seen value for the six checks
##   outcome-check pre-specification section 1 table (governing; open BEFORE computing, so agreement = "answer visible")
## Every seen value is extracted from its file by a pattern asserted to match exactly once; none is typed here.

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
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
OUT <- Sys.getenv("C1_OUT", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_11_extchecks"))   # override used only by C1_compare_mutation.R (a _mut subfolder)
RES <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources")
V3  <- .P31R_V3
LOG <- file.path(OUT, "run.log")
logf <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(s, "\n", sep = "", file = LOG, append = TRUE); cat(s, "\n", sep = ""); invisible(s) }
md5 <- function(p) unname(tools::md5sum(p))
stopf <- function(...) { logf("COMPARE_SEEN STOP: ", ...); quit(status = 4) }
rd_tsv <- function(p) { x <- readLines(p, encoding = "UTF-8"); if (length(x) < 2) stopf("empty table ", p)
  fread(text = x, sep = "\t", colClasses = "character", quote = "", na.strings = NULL) }

logf("==================== C1 COMPARE_SEEN START ====================")
## ---------------------------------------------------------------- 0. order gate: results.tsv exists and is the one the main run logged
P_RES <- file.path(OUT, "results.tsv"); P_BR <- file.path(OUT, "branches.tsv")
if (!file.exists(P_RES) || !file.exists(P_BR)) stopf("results.tsv / branches.tsv missing: run C1_extchecks.R first")
lg <- readLines(LOG); ml <- tail(grep("results\\.tsv md5 [0-9a-f]{32}; branches\\.tsv md5 [0-9a-f]{32}", lg, value = TRUE), 1)
if (!length(ml)) stopf("no results md5 line from the main run in run.log")
m_res_logged <- sub(".*results\\.tsv md5 ([0-9a-f]{32});.*", "\\1", ml); m_br_logged <- sub(".*branches\\.tsv md5 ([0-9a-f]{32}).*", "\\1", ml)
m_res0 <- md5(P_RES); m_br0 <- md5(P_BR)
if (m_res0 != m_res_logged || m_br0 != m_br_logged) stopf("results/branches differ from the md5 the main run logged (", m_res0, " vs ", m_res_logged, ")")
logf("order gate OK: results.tsv md5 ", m_res0, " and branches.tsv md5 ", m_br0, " = the main run's logged md5s; seen sources are opened only now")

## ---------------------------------------------------------------- 1. seen sources, md5 asserted
RJ <- fromJSON(file.path(RES, "RESOURCES.json"), simplifyVector = FALSE)
e_r1 <- Filter(function(z) identical(z$id, "R1_outcome_checks_tsv"), RJ)[[1]]
if (!identical(e_r1$role, "seen_output")) stopf("R1_outcome_checks_tsv role is ", e_r1$role)
SRC <- rbindlist(list(
  list("R1_outcome_checks_tsv", e_r1$path, e_r1$md5),
  list("PRESPEC_FULL", file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC.md"), "3319dc39740d9f5b3b8822133d9108b7"),
  list("AMEND1_FULL", file.path(V3, "R2_rewrite/review_R2/P31R_R2_STAGE1_PRESPEC_AMEND1.md"), "377508832ea9bd384133e869c456c909"),
  list("OUTCOME_CHECK_PRESPEC", file.path(V3, "R1_stage1/S1R1_05b_outcome_checks/S1R1_PRESPEC_outcome_checks.md"), "c9079f2de1ffe75ce31f8a1d597d1fda")))
setnames(SRC, c("id", "path", "md5_expected"))
for (i in seq_len(nrow(SRC))) { m <- md5(SRC$path[i]); logf(sprintf("seen source %-22s md5 %s expected %s  %s", SRC$id[i], m, SRC$md5_expected[i], SRC$path[i]))
  if (m != SRC$md5_expected[i]) stopf("md5 mismatch for ", SRC$id[i]) }

RS <- rd_tsv(P_RES)
val <- function(qid, col = "value") { r <- RS[quantity_id == qid]; if (nrow(r) != 1L) stopf("results.tsv row not unique: ", qid); r[[col]] }

## ---------------------------------------------------------------- 2. comparison helpers
half_ulp <- function(s) {           # half a unit in the last printed digit of a numeric string
  s <- trimws(s); m <- regmatches(s, regexec("^[-+]?([0-9]*)\\.?([0-9]*)(?:[eE]([-+]?[0-9]+))?$", s, perl = TRUE))[[1]]
  if (!length(m)) return(NA_real_)
  d <- nchar(m[3]); ex <- if (nzchar(m[4])) as.numeric(m[4]) else 0
  0.5 * 10^(ex - d)
}
unminus <- function(x) gsub("\u2212", "-", x, fixed = TRUE)
CMP <- list()
add <- function(qid, seen, seen_label, s1, s1_label, abs_diff, expl) CMP[[length(CMP) + 1]] <<- data.table(
  quantity_id = qid, seen_value = seen, seen_label = seen_label, stage1_value = s1, stage1_label = s1_label, abs_diff = abs_diff, explanation = expl)
cmp_num <- function(qid, seen, seen_label, s1, s1_label, note = "", sign_flip = FALSE, why = "") {   # why = a known cause, used only if it differs
  sv <- as.numeric(unminus(seen)); tv <- as.numeric(s1); if (sign_flip) tv <- -tv
  if (is.na(sv) || is.na(tv)) { add(qid, seen, seen_label, s1, s1_label, "NA", paste0("DIFFERS: a value is missing; ", note)); return(invisible()) }
  tol <- max(half_ulp(unminus(seen)), half_ulp(s1), na.rm = TRUE) + 1e-15
  ad <- abs(sv - tv)
  ex <- if (ad <= tol) sprintf("AGREES within the coarser printed precision (tol %.3g)", tol) else
    sprintf("DIFFERS beyond printed precision (tol %.3g); %s", tol, if (nzchar(why)) paste("EXPLAINED:", why) else "UNEXPLAINED")
  add(qid, seen, seen_label, s1, s1_label, formatC(ad, digits = 4, format = "g"), paste(c(ex, note[nzchar(note)]), collapse = "; "))
}
cmp_cat <- function(qid, seen, seen_label, s1, s1_label, note = "", map = identity) {
  ok <- identical(map(seen), map(s1))
  add(qid, seen, seen_label, s1, s1_label, if (ok) "0" else "category differs",
      paste(c(if (ok) "AGREES (same category)" else "DIFFERS (category); UNEXPLAINED", note[nzchar(note)]), collapse = "; "))
}
ANSWER_VISIBLE <- "answer visible: the OCP section 1 table is a governing document that C1 read before computing, and asserted against in-run"
DETERMINISTIC_R1 <- "R1 computed this deterministically from the same files (R1 script is R by filename, not opened); agreement expected and is weak evidence"

## ---------------------------------------------------------------- 3. R1 outcome-check table (seen_output)
R1 <- rd_tsv(e_r1$path)
logf(sprintf("R1 outcome-check table: %d rows; columns %s", nrow(R1), paste(names(R1), collapse = ",")))
PMAP <- c("PAM-T2D" = "PAM_T2D", "RAD23B-T2D" = "RAD23B_T2D", "ARG1-T2D" = "ARG1_T2D", "PCSK9-CAD" = "PCSK9_CAD",
          "TGFB1-CAD" = "TGFB1_CAD", "FGFR4-T2D (EAS)" = "FGFR4_T2D")
if (!all(names(PMAP) %in% R1$pair)) stopf("R1 pair labels not all found: ", paste(setdiff(names(PMAP), R1$pair), collapse = ","))
strand_cat <- function(x) {
  if (grepl("forward-strand assumption", x, fixed = TRUE)) return("forward-strand assumption")
  if (grepl("AF concordant", x, fixed = TRUE) || grepl("frequencies agree", x, fixed = TRUE)) return("palindromic, frequencies agree")
  if (grepl("strand-flipped", x, fixed = TRUE)) return("possibly strand-flipped")
  if (x %in% c("", "\"\"", "not palindromic")) return("not palindromic")
  x
}
for (lab in names(PMAP)) {
  r <- R1[pair == lab]; if (nrow(r) != 1L) stopf("R1 row not unique: ", lab)
  pr <- PMAP[[lab]]; q <- function(k) paste0("C1.", pr, ".", k); SL <- "R1 S1R1_outcome_checks.tsv"
  ## identity of the instrument and effect allele (from the beta_exp estimand text)
  est <- val(q("beta_exp"), "estimand")
  cmp_cat(paste0(q("beta_exp"), "[instrument,EA]"), paste(r$instrument, r$ea), paste(SL, "instrument, ea"),
          sub("^.* at (rs[0-9]+) \\(.*EA ([ACGT]+) / OA.*$", "\\1 \\2", est), "C1 beta_exp estimand", DETERMINISTIC_R1)
  cmp_num(q("eaf_exp"), r$eaf_exposure, paste(SL, "eaf_exposure"), val(q("eaf_exp")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("eaf_out"), r$eaf_resource_aligned, paste(SL, "eaf_resource_aligned"), val(q("eaf_out")), "C1 results value", DETERMINISTIC_R1)
  cmp_cat(q("palindromic"), r$palindrome, paste(SL, "palindrome"), sub("^strand status: ", "", val(q("palindromic"), "notes")),
          "C1 strand status (palindromic notes)", "wording differs between R1 and C1; compared as categories", map = strand_cat)
  cmp_num(q("beta_out"), r$b_out, paste(SL, "b_out"), val(q("beta_out")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("se_out"), r$se_out, paste(SL, "se_out"), val(q("se_out")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("p_resource"), r$p_out_file, paste(SL, "p_out_file"), val(q("p_resource")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("wald"), r$wald, paste(SL, "wald"), val(q("wald")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("se1"), r$wald_se, paste(SL, "wald_se"), val(q("se1")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("se_delta"), r$wald_se_delta, paste(SL, "wald_se_delta"), val(q("se_delta")), "C1 results value", DETERMINISTIC_R1)
  cmp_num(q("p_wald1"), r$p, paste(SL, "p"), val(q("p_wald1")), "C1 results value", DETERMINISTIC_R1)
  cmp_cat(q("class"), r$class, paste(SL, "class"), val(q("class")), "C1 results value", DETERMINISTIC_R1)
  cmp_cat(q("bonf6"), r$bonferroni_6, paste(SL, "bonferroni_6"), val(q("bonf6")), "C1 results value", DETERMINISTIC_R1)
  ## n: R1 leaves DIAMANTE blank (the file has no N column)
  if (nzchar(r$n_resource)) cmp_num(paste0(q("beta_out"), "[n]"), r$n_resource, paste(SL, "n_resource"), val(q("beta_out"), "n"), "C1 beta_out n", DETERMINISTIC_R1) else
    add(paste0(q("beta_out"), "[n]"), "(blank)", paste(SL, "n_resource"), val(q("beta_out"), "n"), "C1 beta_out n", "NA",
        "R1 left n blank (DIAMANTE has no N column); C1 fills n from OCP section 2 case + control counts. Not a disagreement")
  ## extras (not C1 keys)
  dnote <- if (pr == "FGFR4_T2D") "EXTRA. R1 used the OCP-printed AGEN discovery (4 dp); C1 derives it from AGEN raw / CKB beta" else paste("EXTRA.", DETERMINISTIC_R1)
  cmp_num(paste0(q("disc_wald")), r$discovery_wald, paste(SL, "discovery_wald"), val(q("disc_wald")), "C1 results value", dnote)
  cmp_num(paste0(q("disc_wald"), "[se]"), r$discovery_se, paste(SL, "discovery_se"), val(q("disc_wald"), "se"), "C1 disc_wald se", dnote)
  hnote <- paste("EXTRA. Sign convention: R1 het_z = (check - discovery), C1 het_z = (discovery - check); C1 value negated before comparing",
                 if (pr == "FGFR4_T2D") "| FGFR4: R1 used the OCP-printed discovery 0.0482 (SE 0.0069) at 4 dp, C1 the raw-derived discovery, so het_z/het_p differ by input precision" else "")
  hwhy <- if (pr == "FGFR4_T2D") "different discovery inputs (OCP-printed 4 dp in R1; raw-derived in C1)" else ""
  cmp_num(paste0(q("het_z")), r$het_z, paste(SL, "het_z"), val(q("het_z")), "C1 results value (negated)", hnote, sign_flip = TRUE, why = hwhy)
  cmp_num(paste0(q("het_z"), "[p]"), r$het_p, paste(SL, "het_p"), val(q("het_z"), "p"), "C1 het_z p", hnote, why = hwhy)
}
## the R1 row C1 does not compute
for (lab in setdiff(R1$pair, names(PMAP)))
  add(paste0("(R1 row) ", lab), R1[pair == lab]$class, "R1 S1R1_outcome_checks.tsv class", "NOT COMPUTED", "", "NA",
      "outside the C1 brief: the OCP section 6 ARG1 signal-B test is not one of the six pre-declared checks C1 computes")

## ---------------------------------------------------------------- 4. full pre-declaration Seen lines (section 1 RAD23B, section 2 TGFB1)
PF <- unminus(readLines(SRC[id == "PRESPEC_FULL"]$path, encoding = "UTF-8"))
grab <- function(pat, what) { h <- grep(pat, PF, perl = TRUE); if (length(h) != 1L) stopf("full-text pattern for ", what, " matched ", length(h), " lines")
  m <- regmatches(PF[h], regexec(pat, PF[h], perl = TRUE))[[1]]; list(line = h, g = m[-1]) }
SLF <- "full pre-declaration"
g <- grab("resource p ([0-9.]+) at the lead \\(the declared check's instrument\\)", "RAD23B resource p")
cmp_num("C1.RAD23B_T2D.p_resource", g$g[1], sprintf("%s l.%d (section 1 Seen, DIAMANTE for RAD23B)", SLF, g$line), val("C1.RAD23B_T2D.p_resource"), "C1 results value", "seen scoping value; same file, same column")
g <- grab("first-order Wald p ([0-9.]+) \\(rounding interval ([0-9.]+)[-\u2013]([0-9.]+)\\)", "RAD23B first-order Wald p")
cmp_num("C1.RAD23B_T2D.p_wald1", g$g[1], sprintf("%s l.%d (section 1 Seen; rounding interval %s-%s)", SLF, g$line, g$g[2], g$g[3]), val("C1.RAD23B_T2D.p_wald1"), "C1 results value",
        "seen value printed at 2 s.f.; the rounding interval belongs to W1")
g <- grab("second-order ([0-9.]+) \\(interval ([0-9.]+)[-\u2013]([0-9.]+)\\)", "RAD23B second-order Wald p")
cmp_num("C1.RAD23B_T2D.p_wald2", g$g[1], sprintf("%s l.%d (section 1 Seen; interval %s-%s)", SLF, g$line, g$g[2], g$g[3]), val("C1.RAD23B_T2D.p_wald2"), "C1 results value (EXTRA)",
        "seen value printed at 2 s.f.")
g <- grab("^\\s*- (fails) 0\\.05/6 under every p source", "RAD23B 0.05/6")
cmp_cat("C1.RAD23B_T2D.bonf6", "FALSE", sprintf("%s l.%d ('fails 0.05/6 under every p source')", SLF, g$line), val("C1.RAD23B_T2D.bonf6"), "C1 results value",
        "seen statement converted to FALSE; C1 bonf6 is on the resource p (0.05/6 on p_wald1 is in the bonf6 notes)")
g <- grab("rs73045269-T:\\*\\* pQTL \\+([0-9.]+) SD; Wald \\+([0-9.]+) \\(SE ([0-9.]+)\\) in FinnGen CHD.*; \\+([0-9.]+) in CARDIoGRAMplusC4D", "TGFB1 marginal Seen")
SL2 <- sprintf("%s l.%d (section 2 Seen, rs73045269-T marginal)", SLF, g$line)
cmp_num("C1.TGFB1_CAD.beta_exp", g$g[1], SL2, val("C1.TGFB1_CAD.beta_exp"), "C1 results value", "seen at 3 dp")
cmp_num("C1.TGFB1_CAD.disc_wald", g$g[2], SL2, val("C1.TGFB1_CAD.disc_wald"), "C1 results value (EXTRA)", "seen at 3 dp")
cmp_num("C1.TGFB1_CAD.disc_wald[se]", g$g[3], SL2, val("C1.TGFB1_CAD.disc_wald", "se"), "C1 disc_wald se (EXTRA)", "seen at 3 dp")
cmp_num("C1.TGFB1_CAD.wald", g$g[4], SL2, val("C1.TGFB1_CAD.wald"), "C1 results value", "seen at 3 dp")

## the full amendment: searched, no Seen value for the six checks
AF <- readLines(SRC[id == "AMEND1_FULL"]$path, encoding = "UTF-8")
h <- grep("TPMI|CARDIoGRAM|ieu-a-7|PCSK9|PAM|RAD23B|FGFR4", AF)
logf(sprintf("full amendment: lines naming TPMI/CARDIoGRAM/ieu-a-7/PCSK9/PAM/RAD23B/FGFR4 = %d; its DIAMANTE Seen line concerns the section 4 ARG1 A-vs-B test, not a C1 quantity", length(h)))

## ---------------------------------------------------------------- 5. OCP section 1 table (governing; answer visible)
OC <- unminus(readLines(SRC[id == "OUTCOME_CHECK_PRESPEC"]$path, encoding = "UTF-8"))
RSID <- c(PAM_T2D = "rs149802978", RAD23B_T2D = "rs7035725", ARG1_T2D = "rs2781668", PCSK9_CAD = "rs11591147", TGFB1_CAD = "rs73045269", FGFR4_T2D = "rs451643")
for (pr in names(RSID)) {
  h <- grep(paste0("^\\| [^|]+\\| ", RSID[[pr]], ","), OC, perl = TRUE); if (length(h) != 1L) stopf("OCP row for ", pr, " matched ", length(h))
  cells <- trimws(strsplit(OC[h], "|", fixed = TRUE)[[1]])[-1]   # pair, instrument, EA/OA, EAF, beta (SE), discovery outcome, discovery Wald (SE), p
  num1 <- function(s) sub("^\\+", "", regmatches(s, regexpr("^[-+]?[0-9.]+", s)))
  par1 <- function(s) sub("^.*?\\(([0-9.]+)\\).*$", "\\1", s, perl = TRUE)
  SLo <- sprintf("OCP l.%d section 1 table", h); q <- function(k) paste0("C1.", pr, ".", k)
  cmp_num(q("eaf_exp"), num1(cells[4]), paste(SLo, "EAF"), val(q("eaf_exp")), "C1 results value", ANSWER_VISIBLE)
  cmp_num(q("beta_exp"), num1(cells[5]), paste(SLo, "beta_exp"), val(q("beta_exp")), "C1 results value", ANSWER_VISIBLE)
  cmp_num(q("se_exp"), par1(cells[5]), paste(SLo, "beta_exp SE"), val(q("se_exp")), "C1 results value", ANSWER_VISIBLE)
  cmp_num(q("disc_wald"), num1(cells[7]), paste(SLo, "discovery Wald"), val(q("disc_wald")), "C1 results value (EXTRA)", ANSWER_VISIBLE)
  cmp_num(paste0(q("disc_wald"), "[se]"), par1(cells[7]), paste(SLo, "discovery SE"), val(q("disc_wald"), "se"), "C1 disc_wald se (EXTRA)", ANSWER_VISIBLE)
}

## ---------------------------------------------------------------- 6. r0 label (READ_FROM_SEEN by brief)
for (pr in names(RSID)[1:5]) add(paste0("C1.", pr, ".label"), "(the r0 value in S1R1_overlap_ldsc.tsv)", "READ_FROM_SEEN source", val(paste0("C1.", pr, ".label")),
                                 "C1 label", "NA", "READ_FROM_SEEN per the C1 brief: identical by construction; never counts as reproduced")

## ---------------------------------------------------------------- 7. write; assert results/branches untouched
CS <- rbindlist(CMP)
CS[, status := fifelse(startsWith(explanation, "AGREES"), "agrees", fifelse(startsWith(explanation, "DIFFERS"), "differs", "other"))]
fwrite(CS[, .(quantity_id, seen_value, seen_label, stage1_value, stage1_label, abs_diff, explanation)], file.path(OUT, "COMPARE_SEEN.tsv"), sep = "\t", quote = FALSE)
m_res1 <- md5(P_RES); m_br1 <- md5(P_BR)
if (m_res1 != m_res0 || m_br1 != m_br0) stopf("results.tsv or branches.tsv changed during COMPARE_SEEN")
logf(sprintf("COMPARE_SEEN.tsv: %d rows; agrees %d; differs %d (unexplained %d); other %d. md5 %s", nrow(CS), sum(CS$status == "agrees"), sum(CS$status == "differs"),
             sum(grepl("UNEXPLAINED", CS$explanation)), sum(CS$status == "other"), md5(file.path(OUT, "COMPARE_SEEN.tsv"))))
for (k in which(CS$status != "agrees")) logf(sprintf("  COMPARE %-36s seen %s | C1 %s | %s", CS$quantity_id[k], CS$seen_value[k], CS$stage1_value[k], substr(CS$explanation[k], 1, 160)))
logf("results.tsv md5 unchanged ", m_res1, "; branches.tsv md5 unchanged ", m_br1)
