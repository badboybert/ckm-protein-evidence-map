# W6a mutation tests: each mutation is injected into a COPY of w6a_compute.R (inputs mutated in memory right after
# they are read; nothing on disk is changed), the copy writes into mutation_tests/<id>/, and the test passes only if
# the assertion that NAMES the defect goes red (and was green in the unmutated run).
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
suppressMessages(library(data.table))
OUT <- Sys.getenv("W6A_OUT", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))   # override to re-run into another folder (e.g. verify_L3/rerun)
SRC <- Sys.getenv("W6A_SRC", paste0(.P31R_V3, "/R1_stage1/S1R1_R2_06a_renal"))   # where the W6a scripts live
src <- readLines(file.path(SRC, "w6a_compute.R"), encoding = "UTF-8")
base <- fread(file.path(OUT, "w6a_asserts.tsv"))
MUT <- list(
  M1_drop_candidate = list(code = 'CA[Gene == "GUSB" & Endpoint == "eGFR", Status := "Dropped"]',
                           expect = c("cand.routes_agree", "pool.candidates", "primary.table", "c.haldane_or")),
  M2_gene_into_mhc  = list(code = 'INS[gene == "GUSB", `:=`(gene_chrom = 6L, gene_start = 30000000L, gene_end = 30100000L)]',
                           expect = c("pool.mhc", "pool.rows")),
  M3_coloc_rule     = list(code = 'CO[gene == "GUSB" & outcome_name == "eGFR", PP.H3 := 0.1]',
                           expect = c("coloc.core_pass", "pool.colocfail")),
  M4_first_order_se = list(code = 'INS[, se := 1e-12]',
                           expect = c("fdr.rows")),
  M5_cluster_split  = list(code = 'INS[gene == "ITIH1", `:=`(gene_start = 260000000L, gene_end = 260100000L)]',
                           expect = c("pool.clusters")),
  M6_failure_mode   = list(code = 'CO[gene == "APOBR" & outcome_name == "HF", `:=`(PP.H4 = 0.1, PP.H3 = 0.8)]',
                           expect = c("fail.961")),
  # FIXER F1: re-inject the original defect (MLE and CI from fisher.test's default-tolerance uniroot) -> the root check must trip
  M7_fishertest_default_tol = list(code = 'exact_cond_or <- function(tab, conf.level = 0.95, tol = 1e-13) { f <- fisher.test(tab, conf.level = conf.level); list(estimate = unname(f$estimate), conf.int = f$conf.int, resid = exact_cond_or_resid(tab, unname(f$estimate), f$conf.int[1], f$conf.int[2], conf.level)) }',
                           expect = c("fisher.roots.S6.primary", "fisher.roots.S6.b.locus_outcome", "fisher.roots.S6.a")),
  # FIXER F1: solving for the wrong interval (90% instead of 95%) -> the same-estimand check against fisher.test must trip
  M8_wrong_conf_level = list(code = 'EXACT_CONF <- 0.90',
                           expect = c("fisher.matches_fishertest.S6.primary")))
inj <- grep("^GE  <- fread\\(", src); stopifnot(length(inj) == 1)
outl <- grep("^OUT <- ", src); stopifnot(length(outl) == 1)
stl <- grep("^st <- tryCatch\\(system2", src); stopifnot(length(stl) == 1)
rows <- list()
for (m in names(MUT)) {
  d <- file.path(OUT, "mutation_tests", m); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  s <- src
  s[outl] <- sprintf('OUT <- "%s"', d)
  s[stl] <- 'st <- "gls_q.R self-test skipped in mutation run"'
  stopifnot(grepl('^[[:space:]]+error = function', s[stl + 1])); s[stl + 1] <- ''   # the system2 call spans two lines
  s <- append(s, c(paste0("# >>> MUTATION ", m), MUT[[m]]$code, 'lg("MUTATION APPLIED: ', m, '")'), after = inj)
  f <- file.path(d, "w6a_compute_mutated.R"); writeLines(s, f, useBytes = TRUE)
  rc <- system2(file.path(R.home("bin"), "Rscript.exe"), shQuote(f), stdout = file.path(d, "stdout.txt"), stderr = file.path(d, "stdout.txt"))
  a <- if (file.exists(file.path(d, "w6a_asserts.tsv"))) fread(file.path(d, "w6a_asserts.tsv")) else data.table(assert_id = character(), ok = logical())
  for (e in MUT[[m]]$expect) {
    was_ok <- isTRUE(base[assert_id == e, ok]); now <- a[assert_id == e, ok]
    rows[[length(rows) + 1]] <- data.table(mutation = m, injected = MUT[[m]]$code, named_assert = e, baseline_ok = was_ok,
      mutated_ok = if (length(now)) now else NA, unique_target = sum(a$assert_id == e) == 1L,
      tripped = was_ok && length(now) == 1 && isFALSE(now), exit_code = rc,
      other_failed = paste(setdiff(a[ok == FALSE, assert_id], e), collapse = ","))
  }
  cat(m, ": ", paste(sprintf("%s tripped=%s", MUT[[m]]$expect, vapply(MUT[[m]]$expect, function(e) isFALSE(a[assert_id == e, ok]), TRUE)), collapse = "; "), "\n")
}
r <- rbindlist(rows)
fwrite(r, file.path(OUT, "w6a_mutation_tests.tsv"), sep = "\t")
cat("ALL NAMED ASSERTS TRIPPED:", all(r$tripped), "\n")
cat(sprintf("[%s] mutation tests: %d/%d named asserts tripped (w6a_mutation_tests.tsv)\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sum(r$tripped), nrow(r)),
    file = file.path(OUT, "run.log"), append = TRUE)
