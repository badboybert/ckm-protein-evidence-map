# step2_gate_lib.R -- ID-LIMITED expected-change gate (Stage 2 step 2 fix round, 2026-10-10). Identical copy in W6b R/ and W4b root.
# A regenerated results.tsv / branches.tsv may differ from its baseline snapshot ONLY in the (file, key, column) cells and added rows
# listed in the expected-changes table, AND each listed cell must carry exactly the recorded content (md5 of the cell text; md5 of the
# tab-joined row for an added row). Anything else -- a changed value elsewhere, a listed cell with different content, a removed row, an
# unlisted new row, a listed change that did not happen -- is a PROBLEM line naming the id and the column; the caller stops (exit != 0).
# Table columns: file, key, column, new_md5, kind (changed|added), reason. key = quantity_id (results.tsv) or "row N (rule_ref)" (branches.tsv).
suppressMessages({library(data.table); library(digest)})
g_md5 <- function(x) vapply(enc2utf8(x), function(s) digest(s, algo = "md5", serialize = FALSE), "", USE.NAMES = FALSE)
g_rd <- function(p) fread(p, sep = "\t", colClasses = "character", na.strings = NULL, quote = "\"")
g_keys <- function(file, DT) if (file == "results.tsv") DT$quantity_id else sprintf("row %d (%s)", seq_len(nrow(DT)), DT$rule_ref)
g_rowmd5 <- function(DT) g_md5(do.call(paste, c(unname(as.list(DT)), sep = "\t")))
g_clip <- function(x) { x <- gsub("[\r\n\t]+", " ", x); ifelse(nchar(x) > 70, paste0(substr(x, 1, 70), "..."), x) }
g_check <- function(file, nw, ol, EXP) {
  P <- character(); ff <- file; E <- EXP[EXP[["file"]] == ff]
  if (!identical(names(nw), names(ol))) return(sprintf("%s: column names differ from the baseline", file))
  kn <- g_keys(file, nw); ko <- g_keys(file, ol)
  if (anyDuplicated(kn)) P <- c(P, sprintf("%s: duplicate key %s in the new file", file, kn[anyDuplicated(kn)]))
  if (anyDuplicated(ko)) P <- c(P, sprintf("%s: duplicate key %s in the baseline", file, ko[anyDuplicated(ko)]))
  Em <- setNames(E$new_md5, paste(E$key, E$column, sep = "\t")); seen <- character()
  for (k in setdiff(ko, kn)) P <- c(P, sprintf("%s id=%s column=*: row REMOVED", file, k))
  add <- setdiff(kn, ko)
  if (length(add)) {
    rm5 <- g_rowmd5(nw[match(add, kn)])
    for (j in seq_along(add)) {
      ix <- paste(add[j], "*", sep = "\t"); seen <- c(seen, ix)
      if (is.na(Em[ix])) P <- c(P, sprintf("%s id=%s column=*: UNEXPECTED NEW ROW (not in the expected-changes table)", file, add[j]))
      else if (Em[ix] != rm5[j]) P <- c(P, sprintf("%s id=%s column=*: new row differs from the recorded expected row (md5 %s, recorded %s)", file, add[j], rm5[j], Em[ix]))
    }
  }
  common <- intersect(ko, kn); io <- match(common, ko); inn <- match(common, kn)
  for (cc in names(ol)) {
    if (file == "results.tsv" && cc == "quantity_id") next
    a <- ol[[cc]][io]; b <- nw[[cc]][inn]; ch <- which(a != b)
    if (!length(ch)) next
    b5 <- g_md5(b[ch])
    for (j in seq_along(ch)) {
      ix <- paste(common[ch[j]], cc, sep = "\t"); seen <- c(seen, ix)
      if (is.na(Em[ix])) P <- c(P, sprintf("%s id=%s column=%s: UNEXPECTED CHANGE (old '%s' -> new '%s')", file, common[ch[j]], cc, g_clip(a[ch[j]]), g_clip(b[ch[j]])))
      else if (Em[ix] != b5[j]) P <- c(P, sprintf("%s id=%s column=%s: changed cell differs from the recorded expected content (new '%s', md5 %s, recorded %s)", file, common[ch[j]], cc, g_clip(b[ch[j]]), b5[j], Em[ix]))
    }
  }
  for (ix in setdiff(names(Em), seen)) { p <- strsplit(ix, "\t", fixed = TRUE)[[1]]; P <- c(P, sprintf("%s id=%s column=%s: EXPECTED CHANGE ABSENT (listed in the table, not present in the output)", file, p[1], p[2])) }
  P
}
# new_dir / pre_dir hold results.tsv and branches.tsv; exp_path is the expected-changes table. Returns the PROBLEM lines.
g_run <- function(new_dir, pre_dir, exp_path) {
  EXP <- g_rd(exp_path)
  stopifnot(identical(names(EXP), c("file", "key", "column", "new_md5", "kind", "reason")), !anyDuplicated(EXP[, .(file, key, column)]), all(nzchar(EXP$reason)))
  c(g_check("results.tsv", g_rd(file.path(new_dir, "results.tsv")), g_rd(file.path(pre_dir, "results.tsv")), EXP),
    g_check("branches.tsv", g_rd(file.path(new_dir, "branches.tsv")), g_rd(file.path(pre_dir, "branches.tsv")), EXP))
}
