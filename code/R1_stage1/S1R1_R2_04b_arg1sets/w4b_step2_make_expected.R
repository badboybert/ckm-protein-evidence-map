# HONESTY NOTE (2026-10-10): These tables pin the md5 of what the pipeline emitted at step 2; they are REGRESSION guards (reject any later
# change) and do NOT establish that the pinned values are correct. Independent anchors: W6b s_max - A_third + L1/L2 within 2.1e-8; W2/W4b
# s_max, per-fit s and text cells - none.
# w4b_step2_make_expected.R -- ONE-OFF generator of step2_expected_changes.tsv (2026-10-10). NOT part of w4b_run_all.R: the table is the
# recorded, reviewed expected-change list; regenerating it from a new output would defeat the gate. Baseline = fix1_prefix/ (the Stage 1
# FREEZE copy), so the list holds the FIXLOG F1 changes (notes of the six S4.genename.ARG1_ornithine rows; the branch rows F1 touched)
# and the Stage 2 step 2 seven-fit changes. Every cell must get a reason or the script stops.
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
source(paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets/w4b_step2_gate_lib.R"))
D <- paste0(.P31R_V3, "/R1_stage1/S1R1_R2_04b_arg1sets"); PRE <- file.path(D, "fix1_prefix")
reason_res <- function(id, col, kind) {
  if (kind == "added") {
    if (grepl("^S4[.]s[.]max_4fit[.]", id)) return("new row: the four-fit maximum kept as the labelled sensitivity")
    if (grepl("^S4[.]s[.][A-Za-z0-9_]+[.](p2np[.]503|susie[.]489|p1fr[.]489)$", id)) return("new row: s of one of the three added SuSiE fits (s-only, no credible set)")
    return(NA_character_) }
  if (col == "notes" && startsWith(id, "S4.genename.ARG1_ornithine.")) return("FIXLOG F1: notes of the ARG1-ornithine gene-naming rows (grid premise conflict)")
  if (col == "value" && grepl("^S4[.]s[.]max[.]ARG1_", id)) return("s_max is now the maximum over all seven distinct fits (was four)")
  if (col %in% c("notes", "estimand") && grepl("^S4[.](s[.]max[.]ARG1_|cs[.]|genename[.])", id)) return("text quotes the seven-fit s_max / its argmax set (Rule 3(i) restriction text)")
  NA_character_ }
reason_br <- function(col) switch(col, inputs = "ids of the seven-fit / sensitivity rows in inputs", filled_template = "template prints the seven-fit s_max", qualifiers = "qualifier text names the seven-fit and four-fit maxima",
                                  branch = "FIXLOG F1: gene-naming branch text for ARG1-ornithine", NA_character_)
rows <- list()
for (f in c("results.tsv", "branches.tsv")) {
  nw <- g_rd(file.path(D, f)); ol <- g_rd(file.path(PRE, f)); kn <- g_keys(f, nw); ko <- g_keys(f, ol)
  add <- setdiff(kn, ko)
  if (length(add)) { r <- vapply(add, reason_res, "", col = "*", kind = "added"); if (f != "results.tsv" || anyNA(r)) stop("unreasoned added rows: ", paste(add[is.na(r)], collapse = ","))
    rows[[length(rows) + 1]] <- as.data.table(list(file = f, key = add, column = "*", new_md5 = g_rowmd5(nw[match(add, kn)]), kind = "added", reason = unname(r))) }
  common <- intersect(ko, kn); io <- match(common, ko); inn <- match(common, kn)
  for (cc in names(ol)) { if (f == "results.tsv" && cc == "quantity_id") next
    ch <- which(ol[[cc]][io] != nw[[cc]][inn]); if (!length(ch)) next
    rs <- if (f == "results.tsv") mapply(reason_res, common[ch], MoreArgs = list(col = cc, kind = "changed"), USE.NAMES = FALSE) else ifelse(startsWith(common[ch], "row 5 "), "FIXLOG F1: the ARG1-ornithine gene-naming branch (branch, template, qualifiers)", reason_br(cc))
    if (anyNA(rs)) stop("unreasoned changed cells in ", f, " ", cc, ": ", paste(common[ch][is.na(rs)], collapse = ","))
    rows[[length(rows) + 1]] <- as.data.table(list(file = f, key = common[ch], column = cc, new_md5 = g_md5(nw[[cc]][inn][ch]), kind = "changed", reason = unname(rs))) } }
X <- rbindlist(rows)
X[file == "branches.tsv" & startsWith(key, "row 13 ") & column == "filled_template", reason := "template prints the seven-fit s_max and (fix round 2026-10-10) names the seven-fit argmax set instead of calling the maximum a 503-reference value"]
fwrite(X, file.path(D, "step2_expected_changes.tsv"), sep = "\t", quote = TRUE)
cat(sprintf("step2_expected_changes.tsv: %d entries (%d added rows, %d changed cells)\n", nrow(X), sum(X$kind == "added"), sum(X$kind == "changed")))
