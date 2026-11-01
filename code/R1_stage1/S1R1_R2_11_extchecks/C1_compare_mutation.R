#!/usr/bin/env Rscript
## Mutation tests for C1_compare_seen.R: each mutation is applied to a COPY of results.tsv in OUT/_mut/<id>/ and must trip the
## comparator check it names. The real results.tsv, branches.tsv, COMPARE_SEEN.tsv and run.log are never touched (asserted by md5).
## Appends its rows to OUT/mutation_tests.txt.
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
OUT <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_11_extchecks")
CMP <- file.path(OUT, "C1_compare_seen.R")
md5 <- function(p) unname(tools::md5sum(p))
real <- file.path(OUT, c("results.tsv", "branches.tsv", "COMPARE_SEEN.tsv", "run.log")); m0 <- md5(real)
rs <- readLines(file.path(OUT, "results.tsv"))
hdr <- strsplit(rs[1], "\t", fixed = TRUE)[[1]]; iv <- match("value", hdr)

## system2(env=) is ignored on Windows and shQuote() gives POSIX quotes: set the variable in this process (inherited by the child)
## and quote for cmd. The real outputs' md5s are re-checked after every child run.
run_cmp <- function(d) {
  Sys.setenv(C1_OUT = d); on.exit(Sys.unsetenv("C1_OUT"))
  st <- system2("Rscript", shQuote(CMP, type = "cmd"), stdout = FALSE, stderr = FALSE)
  if (!identical(md5(real), m0)) stop("a child run touched the real outputs")
  st
}
run_mut <- function(id, qid, newval, drop_log = FALSE) {
  d <- file.path(OUT, "_mut", id); unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
  k <- which(startsWith(rs, paste0(qid, "\t")))
  if (length(k) != 1L) stop("mutation target not unique: ", qid, " (", length(k), ")")   # refuse loudly
  x <- rs; f <- strsplit(x[k], "\t", fixed = TRUE)[[1]]; old <- f[iv]; f[iv] <- newval; x[k] <- paste(f, collapse = "\t")
  if (identical(old, newval)) stop("mutation is a no-op for ", qid)
  writeLines(x, file.path(d, "results.tsv")); file.copy(file.path(OUT, "branches.tsv"), file.path(d, "branches.tsv"))
  ## the order gate needs the main run's md5 line; for the order-gate mutation the line names the ORIGINAL results md5
  mres <- if (drop_log) md5(file.path(OUT, "results.tsv")) else md5(file.path(d, "results.tsv"))
  writeLines(sprintf("x  results.tsv md5 %s; branches.tsv md5 %s", mres, md5(file.path(d, "branches.tsv"))), file.path(d, "run.log"))
  st <- run_cmp(d)
  cs <- file.path(d, "COMPARE_SEEN.tsv")
  out <- if (file.exists(cs)) fread(text = readLines(cs), sep = "\t", colClasses = "character", quote = "", na.strings = NULL) else NULL
  list(status = st, cs = out, old = old)
}
res <- list()
add <- function(id, check, mutation, tripped) res[[length(res) + 1]] <<- data.table(mutation_id = id, check_named = check, mutation = mutation,
                                                                                   tripped = tripped, verdict = if (isTRUE(tripped)) "PASS (tripped)" else "FAIL (did not trip)")
## C1: a numeric value perturbed beyond printed precision must be reported DIFFERS ... UNEXPLAINED on that id (and only there)
m <- run_mut("C1", "C1.ARG1_T2D.wald", "-0.112829432919")
r <- if (!is.null(m$cs)) m$cs[quantity_id == "C1.ARG1_T2D.wald"] else NULL
add("C1", "COMPARE_SEEN cmp_num (numeric differs -> UNEXPLAINED)", sprintf("C1.ARG1_T2D.wald %s -> -0.112829432919 (copy)", m$old),
    m$status == 0 && nrow(r) == 1L && grepl("^DIFFERS.*UNEXPLAINED", r$explanation) &&
      sum(grepl("UNEXPLAINED", m$cs$explanation)) == 1L)
## C2: a category flipped must be reported DIFFERS (category) UNEXPLAINED
m <- run_mut("C2", "C1.PCSK9_CAD.class", "discordant")
r <- if (!is.null(m$cs)) m$cs[quantity_id == "C1.PCSK9_CAD.class"] else NULL
add("C2", "COMPARE_SEEN cmp_cat (category differs)", "C1.PCSK9_CAD.class concordant -> discordant (copy)",
    m$status == 0 && nrow(r) == 1L && grepl("^DIFFERS \\(category\\); UNEXPLAINED", r$explanation))
## C3: order gate -- results.tsv differs from the md5 the main run logged -> the comparator must refuse (exit 4) and write nothing
m <- run_mut("C3", "C1.TGFB1_CAD.p_resource", "0.5", drop_log = TRUE)
add("C3", "COMPARE_SEEN order gate (results md5 = main run's logged md5)", "results.tsv copy altered after the logged md5",
    m$status == 4 && is.null(m$cs))
## C4: control -- an unmutated copy must give zero UNEXPLAINED (otherwise C1/C2 prove nothing)
d <- file.path(OUT, "_mut", "C4"); unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
file.copy(file.path(OUT, c("results.tsv", "branches.tsv")), d)
writeLines(sprintf("x  results.tsv md5 %s; branches.tsv md5 %s", md5(file.path(d, "results.tsv")), md5(file.path(d, "branches.tsv"))), file.path(d, "run.log"))
st <- run_cmp(d)
cs <- fread(text = readLines(file.path(d, "COMPARE_SEEN.tsv")), sep = "\t", colClasses = "character", quote = "", na.strings = NULL)
add("C4", "control (unmutated copy)", "none", st == 0 && sum(grepl("UNEXPLAINED", cs$explanation)) == 0L &&
      identical(md5(file.path(d, "COMPARE_SEEN.tsv")), md5(file.path(OUT, "COMPARE_SEEN.tsv"))))

R <- rbindlist(res)
stopifnot(identical(md5(real), m0))   # the real outputs were not touched
unlink(file.path(OUT, "_mut"), recursive = TRUE)
cat(c("", "# COMPARE_SEEN mutation tests (C1_compare_mutation.R): mutations on copies in _mut/, removed afterwards; C4 is the unmutated control.",
      apply(R, 1, paste, collapse = "\t")), sep = "\n", file = file.path(OUT, "mutation_tests.txt"), append = TRUE)
cat(sprintf("%s  COMPARE_SEEN mutation tests: %d of %d PASS; real outputs unchanged (md5)\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sum(R$tripped), nrow(R)),
    file = file.path(OUT, "run.log"), append = TRUE)
print(R)
