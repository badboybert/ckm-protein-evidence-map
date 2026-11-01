#!/usr/bin/env Rscript
# w4a_staging_beside.R -- Stage 1 W4a, step 3 (brief W4a: "Fresh R coloc at p12 = 1e-5 governs; print the staging value
# beside it (READ_FROM_SEEN, after COMPARE_SEEN only); report any change in firing").
# Runs only after COMPARE_SEEN.tsv exists. APPENDS READ_FROM_SEEN rows to results.tsv; the computed rows are not touched
# (asserted: the frozen pre-compare file is a byte-identical prefix of the new file). These rows never count as "reproduced".
#   Rscript w4a_staging_beside.R
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
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04a_arg1grid")
WORK <- file.path(OUT, "work")
source(file.path(OUT, "w4a_funcs.R"))
log_init(file.path(OUT, "run.log"))
logf("==== W4a staging-beside step (READ_FROM_SEEN) start ====")
if (!file.exists(file.path(OUT, "COMPARE_SEEN.tsv"))) stop("STOP: COMPARE_SEEN.tsv must exist before any seen output is read")
FROZEN_MD5 <- readLines(file.path(WORK, "FREEZE.txt"))[1]   # results.tsv md5 written by run_all.R before the full text was opened
rp <- file.path(OUT, "results.tsv"); fz <- file.path(WORK, "results.precompare_frozen.tsv")
lfr <- readLines(rp)
if (md5f(rp) == FROZEN_MD5) file.copy(rp, fz, overwrite = TRUE) else if (!(file.exists(fz) && md5f(fz) == FROZEN_MD5 && identical(lfr[seq_along(readLines(fz))], readLines(fz))))
  stop("STOP: results.tsv is neither the frozen version nor frozen + appended rows")
assert_md5(fz, FROZEN_MD5, "frozen pre-compare results.tsv")
RES <- fromJSON(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_00_resources/RESOURCES.json"), simplifyVector = FALSE)
names(RES) <- vapply(RES, function(e) e$id, "")
rows <- list()
for (pn in c("503", "489")) {
  id <- if (pn == "503") "LOC_ARG1_6q23_two_signal_p1np" else "LOC_ARG1_6q23_two_signal_p2np"
  e <- RES[[id]]
  assert_md5(e$path, e$md5, sprintf("%s (%s; READ_FROM_SEEN after COMPARE_SEEN)", id, e$role))
  ts <- read.delim(e$path, stringsAsFactors = FALSE, colClasses = "character")
  logf("%s: %d rows; rows used: protein P__ARG1, A rs2781668, B rs17788484", id, nrow(ts))
  for (tr in list(c("arginine", "M__GCST90200372"), c("ornithine", "M__GCST90200404"))) {
    x <- ts[ts$protein == "P__ARG1" & ts$other == tr[2] & ts$sig_A == "rs2781668" & ts$sig_B == "rs17788484", ]
    if (!nrow(x)) { logf("no staging row for %s on %s", tr[1], pn); next }
    x <- x[1, ]
    for (sg in c("A", "B")) {
      h4 <- as.numeric(x[[paste0("cond_H4_", sg)]]); z <- as.numeric(x[[paste0("M", sg)]])
      fresh_h4 <- NA; fresh_z <- NA
      cur <- fread(file = fz, sep = "\t", colClasses = "character")
      fresh_h4 <- as.numeric(cur$value[cur$quantity_id == sprintf("S4.coloc_clause.%s.%s.%s.PPH4", tr[1], sg, pn)])
      fresh_z <- as.numeric(cur$value[cur$quantity_id == sprintf("S4.coloc_clause.%s.%s.%s.z", tr[1], sg, pn)])
      nt <- sprintf("READ_FROM_SEEN (staging two_signal.tsv, signal_rule %s; engine coloc.abf p12 1e-5); fresh Stage 1 R coloc %.4f (z %.3f) governs; never counts as reproduced",
                    x$signal_rule, fresh_h4, fresh_z)
      rows[[length(rows) + 1]] <- data.frame(quantity_id = sprintf("S4.coloc_clause.%s.%s.%s.PPH4_staging", tr[1], sg, pn), prespec_ref = "Rule 3(iv) (staging value beside)",
        estimand = sprintf("conditional PP.H4 ARG1 x %s at signal %s (staging)", tr[1], sg), resource = id, panel = pn, cov = "", se_order = "",
        value = sprintf("%.7g", h4), se = "NA", ci_lo = "NA", ci_hi = "NA", p = "NA", n = x$nsnp, unit = "PP.H4", notes = nt)
      rows[[length(rows) + 1]] <- data.frame(quantity_id = sprintf("S4.coloc_clause.%s.%s.%s.z_staging", tr[1], sg, pn), prespec_ref = "Rule 3(iv) (staging value beside)",
        estimand = sprintf("conditional z of %s at signal %s (staging)", tr[1], sg), resource = id, panel = pn, cov = "", se_order = "",
        value = sprintf("%.7g", z), se = "NA", ci_lo = "NA", ci_hi = "NA", p = "NA", n = x$nsnp, unit = "z", notes = nt)
      logf("staging beside: %s signal %s panel %s: staging PP.H4 %.4f z %.3f | fresh PP.H4 %.4f z %.3f", tr[1], sg, pn, h4, z, fresh_h4, fresh_z)
    }
  }
}
SB <- rbindlist(rows)
# firing with staging values vs fresh values (per trait x signal; associated = |z| >= 1.96 on 503; fires if PP.H4 < 0.80 on either panel)
cur <- fread(file = fz, sep = "\t", colClasses = "character")
for (tn in c("arginine", "ornithine")) for (sg in c("A", "B")) {
  g <- function(pn, kind, stg) { q <- sprintf("S4.coloc_clause.%s.%s.%s.%s%s", tn, sg, pn, kind, if (stg) "_staging" else "")
    as.numeric(if (stg) SB$value[SB$quantity_id == q] else cur$value[cur$quantity_id == q]) }
  fire <- function(stg) abs(g("503", "z", stg)) >= 1.96 && (g("503", "PPH4", stg) < 0.8 || g("489", "PPH4", stg) < 0.8)
  logf("Rule 3(iv) firing %s signal %s: fresh %s, staging %s -> %s", tn, sg, fire(FALSE), fire(TRUE), if (fire(FALSE) == fire(TRUE)) "no change in firing" else "CHANGE IN FIRING")
}
logf("Rule 3(iv) T2D: no staging two_signal value exists (T2D is not in the l3-finemap harmonised table); fresh values only")
lf <- readLines(fz)
if (!identical(strsplit(lf[1], "\t")[[1]], names(SB))) stop("STOP: column mismatch")
writeLines(c(lf, do.call(paste, c(as.list(SB), sep = "\t"))), rp)   # append as text: the frozen lines are copied verbatim
ln <- readLines(rp)
if (!identical(ln[seq_along(lf)], lf)) stop("STOP: computed rows changed")
logf("results.tsv: %d computed rows (frozen prefix byte-identical) + %d READ_FROM_SEEN rows appended; new md5 %s", length(lf) - 1, nrow(SB), md5f(rp))
logf("==== W4a staging-beside step end ====")
