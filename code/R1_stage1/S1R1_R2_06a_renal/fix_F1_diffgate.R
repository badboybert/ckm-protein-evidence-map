# fix_F1_diffgate.R -- scope gate for FIXER F1 (FIXLOG.md). Used by W6a and W6b (one implementation).
# It proves that a regenerated results.tsv / branches.tsv differs from its preserved pre-fix copy only where defect F1 (the
# root-finding tolerance of the exact CI) allows:
#   results.tsv  - same rows in the same order; a cell may differ only in a row whose quantity_id matches <results id regex>,
#                  and only in the columns value, ci_lo, ci_hi and notes;
#   branches.tsv - same rows in the same order; rule_ref and branch identical everywhere (no branch decision may change);
#                  an inputs / filled_template / qualifiers cell may differ only if its text is identical once every number is
#                  masked, or if it is listed as an allowed cell <row>:<column>.
# Usage:  Rscript fix_F1_diffgate.R <folder> <prefix folder> <results id regex> [<row>:<column> ...]
# Writes <folder>/fix_F1_diffgate.tsv (every differing cell, with its verdict) and exits non-zero on any change outside scope.
# Mutation test: set FIX_F1_MUTATE=1 to also run the gate on in-memory copies carrying an out-of-scope change (a p value; a
# branch decision; a table count); each must trip, and the result is written to <folder>/fix_F1_diffgate_mutation.tsv.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
DIR <- args[1]; PRE <- args[2]; RE <- args[3]; ALLOW_CELLS <- args[-(1:3)]
rd <- function(p) fread(p, sep = "\t", colClasses = "character", na.strings = NULL, quote = "\"")
mask <- function(x) gsub("[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?", "#", x)

gate <- function(new_r, old_r, new_b, old_b) {
  out <- list(); add <- function(file, row, col, old, new, ok, why) out[[length(out) + 1]] <<- data.table(file, row, column = col, old, new, allowed = ok, why)
  if (!identical(names(new_r), names(old_r)) || nrow(new_r) != nrow(old_r) || !identical(new_r$quantity_id, old_r$quantity_id))
    add("results.tsv", "(structure)", "", sprintf("%d rows", nrow(old_r)), sprintf("%d rows", nrow(new_r)), FALSE, "row set, row order or columns differ")
  else for (cc in names(old_r)) {
    d <- which(old_r[[cc]] != new_r[[cc]])
    for (i in d) add("results.tsv", old_r$quantity_id[i], cc, old_r[[cc]][i], new_r[[cc]][i],
                     grepl(RE, old_r$quantity_id[i]) && cc %in% c("value", "ci_lo", "ci_hi", "notes"),
                     if (grepl(RE, old_r$quantity_id[i])) "F1 row" else "row outside F1 scope")
  }
  if (!identical(names(new_b), names(old_b)) || nrow(new_b) != nrow(old_b))
    add("branches.tsv", "(structure)", "", sprintf("%d rows", nrow(old_b)), sprintf("%d rows", nrow(new_b)), FALSE, "row count or columns differ")
  else for (cc in names(old_b)) {
    d <- which(old_b[[cc]] != new_b[[cc]])
    for (i in d) {
      cell <- sprintf("%d:%s", i, cc)
      ok <- !(cc %in% c("rule_ref", "branch")) && (mask(old_b[[cc]][i]) == mask(new_b[[cc]][i]) || cell %in% ALLOW_CELLS)
      add("branches.tsv", sprintf("row %d (%s)", i, old_b$rule_ref[i]), cc, old_b[[cc]][i], new_b[[cc]][i], ok,
          if (cc %in% c("rule_ref", "branch")) "branch decision changed" else if (mask(old_b[[cc]][i]) == mask(new_b[[cc]][i])) "numbers only" else if (cell %in% ALLOW_CELLS) "allowed cell (label text)" else "text changed")
    }
  }
  if (length(out)) rbindlist(out) else data.table(file = character(), row = character(), column = character(), old = character(), new = character(), allowed = logical(), why = character())
}

new_r <- rd(file.path(DIR, "results.tsv")); old_r <- rd(file.path(PRE, "results.tsv"))
new_b <- rd(file.path(DIR, "branches.tsv")); old_b <- rd(file.path(PRE, "branches.tsv"))
G <- gate(new_r, old_r, new_b, old_b)
fwrite(G, file.path(DIR, "fix_F1_diffgate.tsv"), sep = "\t", quote = TRUE)
md5f <- function(p) unname(as.character(tools::md5sum(p)))
cat(sprintf("fix_F1_diffgate: %s | results %d rows (old md5 %s, new %s) | branches %d rows (old %s, new %s) | differing cells %d (results %d, branches %d) | OUT OF SCOPE %d\n",
            DIR, nrow(new_r), md5f(file.path(PRE, "results.tsv")), md5f(file.path(DIR, "results.tsv")), nrow(new_b), md5f(file.path(PRE, "branches.tsv")),
            md5f(file.path(DIR, "branches.tsv")), nrow(G), sum(G$file == "results.tsv"), sum(G$file == "branches.tsv"), sum(!G$allowed)))

if (Sys.getenv("FIX_F1_MUTATE") == "1") {
  M <- list()
  m1 <- copy(new_r); i <- which(grepl("\\.p$", m1$quantity_id))[1]; m1$value[i] <- paste0(m1$value[i], "1")
  M[[1]] <- data.table(mutation = sprintf("results value of %s altered (a p value: outside F1 scope)", m1$quantity_id[i]), tripped = any(!gate(m1, old_r, new_b, old_b)$allowed))
  m2 <- copy(new_b); j <- which(nzchar(m2$branch))[1]; m2$branch[j] <- paste(m2$branch[j], "X")
  M[[2]] <- data.table(mutation = sprintf("branches row %d branch text altered", j), tripped = any(!gate(new_r, old_r, m2, old_b)$allowed))
  m3 <- copy(new_r); k <- which(grepl("\\.a$", m3$quantity_id))[1]; m3$value[k] <- "999"
  M[[3]] <- data.table(mutation = sprintf("results count %s altered", m3$quantity_id[k]), tripped = any(!gate(m3, old_r, new_b, old_b)$allowed))
  m4 <- copy(new_b); j4 <- which(nzchar(m4$filled_template))[1]; m4$filled_template[j4] <- paste(m4$filled_template[j4], "and more")
  M[[4]] <- data.table(mutation = sprintf("branches row %d filled_template words altered", j4), tripped = any(!gate(new_r, old_r, m4, old_b)$allowed))
  MT <- rbindlist(M); fwrite(MT, file.path(DIR, "fix_F1_diffgate_mutation.tsv"), sep = "\t")
  cat(sprintf("fix_F1_diffgate mutation tests: %d/%d tripped\n", sum(MT$tripped), nrow(MT)))
  if (!all(MT$tripped)) quit(status = 3)
}
if (any(!G$allowed)) { print(G[allowed == FALSE]); quit(status = 2) }
if (!any(G$file == "results.tsv")) { cat("fix_F1_diffgate: no results.tsv cell differs from the pre-fix copy -- the fix is not in this file\n"); quit(status = 4) }
